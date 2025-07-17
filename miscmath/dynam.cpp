
//    --------------------------------------------------------------------
//
//    This file is part of Luna.
//
//    LUNA is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    Luna is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with Luna. If not, see <http://www.gnu.org/licenses/>.
//
//    Please see LICENSE.txt for more details.
//
//    --------------------------------------------------------------------

#include "miscmath/dynam.h"
#include "miscmath/miscmath.h"
#include "helper/helper.h"
#include "dsp/tv.h"
#include "dsp/spline.h"
#include "db/db.h"
#include "eval.h"
#include "edf/edf.h"
#include "stats/gpa.h" // linmod_t
#include "stats/Eigen/Dense"
#include "stats/eigen_ops.h"
#include <cstddef>

extern writer_t writer;


dynam_t::dynam_t( const std::vector<double> & y ) : y(y) 
{
  t.resize( y.size() );
  for (int i=0;i<t.size();i++) t[i] = i;
}

dynam_t::dynam_t( const std::vector<double> & y , const std::vector<double> & t ) : y(y), t(t) 
{
  if ( y.size() != t.size() ) Helper::halt( "dynam_t given unequal y and t lengths" );
}
  
dynam_t::dynam_t( const std::vector<double> & y , const std::vector<int> & ti ) : y(y)
{
  if ( y.size() != ti.size() ) Helper::halt( "dynam_t given unequal y and t lengths" );
  t.resize(y.size());
  for (int i=0;i<t.size();i++) t[i] = ti[i];
}
 
   
void dynam_t::denoise( double lambda )
{ 
  dsptools::TV1D_denoise( y , lambda );
}
  

bool dynam_t::mean_variance( double * mean , double * var )
{
  const int n = y.size();

  double my = 0;

  // mean only?
  if ( var == NULL )
    {
      if ( n < 1 ) return false;
      for (int i=0;i<n;i++) 
	my += y[i];      
      my /= n;
      *mean = my;
      return true;
    }

  // mean & variance
  if ( n < 2 ) return false;

  double myy = 0;  
  for (int i=0;i<n;i++) {
    my += y[i];
    myy += y[i] * y[i];
  }
  
  my /= n;
  myy /= n;

  *mean = my;
  *var = myy - my*my;
  
  return true;
}


bool dynam_t::linear_trend( double * beta , double * rsq , double * intercept )
{

  const int n = y.size();
  double mx = 0;
  double my = 0;
  double mxy = 0;
  double mxx = 0;
  double myy = 0;

  for (int i=0;i<n;i++) {
    my += y[i];
    mx += t[i];
    mxy += y[i] * t[i];
    mxx += t[i] * t[i];
    myy += y[i] * y[i];
  }

  mx /= n;
  my /= n;
  mxy /= n;
  mxx /= n;
  myy /= n;

  double varx = mxx - mx*mx;
  double vary = myy - my*my;

  if ( varx == 0 ) return false;

  *beta = ( mxy - ( mx *  my ) ) / varx ; 

  if ( intercept != NULL ) 
    *intercept = my - *beta * mx;
  
  if ( rsq != NULL )
    {
      if ( vary == 0 ) rsq = 0;
      else  
	{
	  *rsq = ( mxy - mx * my ) / sqrt( varx * vary );
	  *rsq *= *rsq;
	}
    }
  return true;
}
    
void dynam_t::hjorth( double * h1 , double *h2 , double *h3 )
{  
  MiscMath::hjorth( &y , h1 , h2 , h3 , ! globals::legacy_hjorth );
}


//
// group-stratified dynamics
//

gdynam_t::gdynam_t( const std::vector<int> & g , const std::vector<double> & y ) : g(g) , y(y) 
{
  if ( g.size() != y.size() ) Helper::halt( "problem in gdynam_t" );
  t.resize( y.size() );
  for (int i=0;i<t.size();i++) t[i] = i;
}

gdynam_t::gdynam_t( const std::vector<int> & g , const std::vector<double> & y , const std::vector<double> & t ) : g(g), y(y), t(t)
{
  if ( g.size() != y.size() ) Helper::halt( "problem in gdynam_t" );
  if ( g.size() != t.size() ) Helper::halt( "problem in gdynam_t" );
}

gdynam_t::gdynam_t( const std::vector<int> & g , const std::vector<double> & y , const std::vector<int> & ti ) : g(g) , y(y) 
{
  if ( g.size() != y.size() ) Helper::halt( "problem in gdynam_t" );
  if ( g.size() != ti.size() ) Helper::halt( "problem in gdynam_t" ); 
  t.resize(y.size());
  for (int i=0;i<t.size();i++) t[i] = ti[i];  
}





int gdynam_t::stratify()
{

  // maps strata to epochs
  gmap.clear();

  for (int i=0;i<g.size();i++) 
    {
      if ( g[i] > 0 )  // only strata 1, 2, 3 etc are included
	gmap[ g[i] ].insert( i );
    }
  
  // number of strata
  
  const int ng = gmap.size();
  
  if ( ng < 2 ) return ng;


  // get means for a between-strata model
  
  between.clear();
  
  // create a dynam_t for each strata (which are numbered 1..ng)
  
  for (int i=1;i<=ng;i++)
    {
      std::vector<double> yy, tt;
      std::set<int>::const_iterator ii = gmap[i].begin();

      while ( ii != gmap[i].end() ) 
	{
	  yy.push_back( y[ *ii ] );
	  tt.push_back( t[ *ii ] );
	  ++ii;
	}
      
      within[ i ].y = yy;
      within[ i ].t = tt;      
      
      // means for between-group model
      between.y.push_back( MiscMath::mean( yy ) );
      between.t.push_back( MiscMath::mean( tt ) );
      
    }

  return ng;
}

  


dissipation_t::dissipation_t( const std::vector<double> & x , 
			      const int mx ,  
			      const double winsor  )
{
  // assumptions::: positive values
  // treated as contiguous, e.g. all spliced NREM epochs
  
  std::vector<double> y = x;

  if ( mx ) y.resize( mx , 0 );
      
  const int np = y.size();
  
  if ( winsor > 0 ) 
    MiscMath::winsorize( &y, winsor );
  
  // make cumulative sum
  s.resize( np, 0 );
  double sum = 0;
  for (int p=0; p<np; p++)
    {
      if ( y[p] < 0 ) Helper::halt( "dissipation_t() expects only positive inputs" );
      sum += y[p];
      s[p] = sum ;
    }
  
  // as a proportion
  for (int p=0; p<np; p++)
    s[p] /= sum;
  
  
}

  
std::vector<double> dissipation_t::plife( const std::vector<double> & ps )
{

  // make spline
  // note t is the DV
  // i.e. to predict t (number of epochs) given a percentile value

  const int np = s.size();
  
  std::vector<double> t0( np );
  for (int p=0; p<np; p++) t0[p] = p;
  
  tk::spline spline;
  spline.set_points( s, t0 );
  const int n = ps.size();
  
  std::vector<double> res( n );

  for (int i=0; i<n; i++)
    {      
      if ( ps[i] < 0 || ps[i] > 1  ) 
	Helper::halt( "invalid spline call" );
      res[i] = spline( ps[i] );
    }
  
  return res;
}




