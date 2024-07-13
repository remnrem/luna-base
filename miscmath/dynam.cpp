
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

#include <cstddef>

extern writer_t writer;

// helper(s)

//
// helper functions
//

bool dynam_compile_cycles( edf_t & edf , const std::vector<int> & epochs, std::vector<std::string> * cycles )
{

  // data must be epoched already (i.e. will be if cycles done)
  if ( ! edf.timeline.epoched() ) return false;
  if ( cycles == NULL ) return false;
  
  // do we have any cycles? 
  const bool has_cycles = edf.timeline.epoch_annotation( "_NREMC_1" ) 
    || edf.timeline.epoch_annotation( "_NREMC_2" )
    || edf.timeline.epoch_annotation( "_NREMC_3" )
    || edf.timeline.epoch_annotation( "_NREMC_4" )
    || edf.timeline.epoch_annotation( "_NREMC_5" )
    || edf.timeline.epoch_annotation( "_NREMC_6" )
    || edf.timeline.epoch_annotation( "_NREMC_7" )
    || edf.timeline.epoch_annotation( "_NREMC_8" );
  
  if ( ! has_cycles ) return false; 
  
  // construct by iterating over current epoch set
  
  cycles->clear();
  
  for (int e=0;e<epochs.size(); e++)
    {
      
      std::string c = "."; // null                                                                                                                     

      // nb. uses legacy epoch-annotation encoding
      // take up to 8 cycles
      // nb. epoch_annot() takes 0..ne current epoch encoding (will map as needed)
      if      ( edf.timeline.epoch_annotation( "_NREMC_1" , e ) ) c = "C1";
      else if ( edf.timeline.epoch_annotation( "_NREMC_2" , e ) ) c = "C2";
      else if ( edf.timeline.epoch_annotation( "_NREMC_3" , e ) ) c = "C3";
      else if ( edf.timeline.epoch_annotation( "_NREMC_4" , e ) ) c = "C4";
      else if ( edf.timeline.epoch_annotation( "_NREMC_5" , e ) ) c = "C5";
      else if ( edf.timeline.epoch_annotation( "_NREMC_6" , e ) ) c = "C6";
      else if ( edf.timeline.epoch_annotation( "_NREMC_7" , e ) ) c = "C7";
      else if ( edf.timeline.epoch_annotation( "_NREMC_8" , e ) ) c = "C8";
      
      cycles->push_back( c );
      
      // next epoch
    }
  
  return true;
}


// wrapper(s)

void dynam_report_with_log( param_t & param,
                            const std::vector<double> & y ,
                            const std::vector<int> & t ,
                            const std::vector<std::string> * g )
{
  std::vector<double> tl( t.size() );
  for (int i=0;i<t.size();i++) tl[i] = t[i];
  dynam_report_with_log( param, y , tl , g );
}



void dynam_report( param_t & param,
                   const std::vector<double> & y ,
                   const std::vector<int> & t ,
                   const std::vector<std::string> * g )
{
  std::vector<double> tl( t.size() );
  for (int i=0;i<t.size();i++) tl[i] = t[i];
  dynam_report( param, y , tl , g ); 
}

  
void dynam_report_with_log( param_t & param,
			    const std::vector<double> & y , 
			    const std::vector<double> & t , 
			    const std::vector<std::string> * g )
{
  std::vector<double> yl( y.size() );
  for (int i=0;i<y.size();i++) yl[i] = 10 * log10( y[i] );
  dynam_report( param, yl , t , g );
}



void dynam_report( param_t & param,
		   const std::vector<double> & y_ , 
		   const std::vector<double> & t_ , 
		   const std::vector<std::string> * g_ )
{

  //
  // method 1 (assume called by dynam_report_with_log() and so from welch at least, this will always
  // be log-transformed already;  do this before outlier removal (do winsor) so call this before the
  // older steps
  //

  //
  // do we have NREM cycle information?
  //
  
  const bool has_cycles = g_ != NULL ;

  qdynam_t qd( y_.size() , g_ );

  //
  // pass epoch-level information (for times) [ assuming passed display_epoch() - 1 as per welch ] 
  //

  std::vector<int> e_( t_.size() );
  for (int i=0;i<t_.size();i++) e_[i] = t_[i];
  qd.set_epochs( e_ );

  //
  // options
  //

  const bool verbose = param.has( "dynam-verbose" );
  const bool epoch_output = param.has( "dynam-epoch" );

  if ( param.has( "dynam-min-ne" ) ) qd.set_min_ne( param.requires_int( "dynam-min-ne" ) );
  
  const double qd_winsor = param.has( "dynam-winsor" ) ? param.requires_dbl( "dynam-winsor" ) : 0.05 ;  
  qd.winsorize( qd_winsor );

  if ( param.has( "dynam-median-window" ) )
    qd.set_smoothing_median_window( param.requires_int( "dynam-median-window" ) );
  
  if ( param.has( "dynam-mean-window" ) )
    qd.set_smoothing_mean_window( param.requires_int( "dynam-mean-window" ) );
  
  if ( param.has( "dynam-norm-mean" ) ) qd.set_norm_mean( param.yesno( "dynam-norm-mean" ) );
  else if ( param.has( "dynam-norm-max" ) ) qd.set_norm_max( param.yesno( "dynam-norm-max" ) );
  
  if ( param.has( "dynam-norm-cycles" ) ) qd.set_norm_cycles( param.yesno( "dynam-norm-cycles" ) );
  
  if ( param.has( "dynam-max-cycle" ) ) qd.set_max_cycles( param.requires_int( "dynam-max-cycle" ) );
  else if ( param.has( "dynam-cycles" ) ) qd.set_cycles( param.intvector( "dynam-cycles" ) );

  // default false
  if ( param.has( "dynam-weight-cycles" ) )
    qd.set_weight_cycles( param.yesno( "dynam-weight-cycles" ) );
  
  //
  // process
  //
  
  qd.proc( y_ );
  
  //
  // output
  //

  writer.level( "TOT" , "QD" );

  qdynam_t::output_helper( qd.r1 , verbose );

  if ( has_cycles )
    {
      // between (only if not norming within each cycle)
      if ( ! qd.norm_cycles() ) 
	{
	  writer.level( "BETWEEN" , "QD" );
	  qdynam_t::output_helper( qd.rb , verbose , true );
	}

      // average within
      writer.level( "WITHIN" , "QD" );
      qdynam_t::output_helper( qd.rwa , verbose );

      // each cycle
      std::map<std::string,qdynam_results_t> & cycs = qd.rw;
      std::map<std::string,qdynam_results_t>::const_iterator cc = cycs.begin();
      while ( cc != cycs.end() )
	{
	  writer.level( "W_" + cc->first , "QD" );
	  qdynam_t::output_helper( cc->second , verbose );
	  ++cc;
	}
    }
  
  writer.unlevel( "QD" );

  //
  // q10 outputs  
  //


  // overall
  const std::vector<double> & ss = qd.r1_q10;
  const std::vector<double> & os = qd.r1_os_q10;
  
  // epochs contains nq(=10) quantiles
  
  writer.level( "TOT" , "QD" );
  for (int i=0; i< ss.size(); i++)
    {
      writer.level(i+1,"Q");
      writer.value( "SS" , ss[i] );
      writer.value( "OS" , os[i] );
    }
  writer.unlevel( "Q" );
  
  
  // cycles
  if ( has_cycles )
    {
      
      std::map<std::string,std::vector<double> > & cycs = qd.rw_q10;
      std::map<std::string,std::vector<double> >::const_iterator cc = cycs.begin();
      
      while ( cc != cycs.end() )
	{
	  writer.level( "W_" + cc->first , "QD" );
	  
	  const std::vector<double> & ss = cc->second;
	  const std::vector<double> & os = qd.rw_os_q10[ cc->first ];
		    
	  for (int i=0; i<ss.size(); i++)
	    {
	      writer.level( i + 1 , "Q" ); 
	      writer.value( "SS" , ss[i] );

	      // do not show if not norming each section
	      // (just the way things are calculated internally, we don't get this)
	      if ( qd.norm_cycles() )
		writer.value( "OS" , os[i] );
	    }
	  writer.unlevel( "Q" );
	  
	  ++cc;
	}
    }
  
  writer.unlevel( "QD" );

      
  //
  // optional outputs
  //

  if ( epoch_output )
    {

      // overall
      const std::vector<double> & ss = qd.r1_smoothed_series;
       
      if ( e_.size() == ss.size() )
	{
	  // epochs contains display_epoch() - 1

	  writer.level( "TOT" , "QD" );
	  for (int i=0; i<ss.size(); i++)
	    {
	      writer.epoch( e_[i] + 1 ); // 1-based outputs
	      writer.value( "SS" , ss[i] );
	    }
	  writer.unepoch();
	}


      // cycles
      if ( has_cycles )
	{

	  std::map<std::string,std::vector<double> > & cycs = qd.rw_smoothed_series;
	  std::map<std::string,std::vector<double> >::const_iterator cc = cycs.begin();

	  while ( cc != cycs.end() )
	    {
	      writer.level( "W_" + cc->first , "QD" );
	      
	      const std::vector<double> & ss = cc->second;
	      const std::vector<int> & ee = qd.rw_epochs[ cc->first ];
	      
	      for (int i=0; i<ss.size(); i++)
		{
		  writer.epoch( ee[i] + 1 ); // 1-based outputs
		  writer.value( "SS" , ss[i] );
		}
	      writer.unepoch();
	      
	      ++cc;
	    }
	}
      
      writer.unlevel( "QD" );
      
    }
      
  
  

  //
  // (original) method 2 (now requires dynam-ols flag)
  //


  if ( ! param.has( "dynam-ols" ) ) return;
  if ( ! param.yesno( "dynam-ols" ) ) return;
  

  //
  // Remove 'y' outliers
  //
  const double th = 3;

  std::vector<double> z0 = MiscMath::Z( y_ );

  std::vector<double> y, t;
  std::vector<std::string> g;

  for (int i=0;i<y_.size();i++)
    {
      if ( z0[i] >= -th && z0[i] <= th ) 
	{
	  y.push_back( y_[i] );
	  t.push_back( t_[i] );
	  if ( g_ != NULL ) g.push_back( (*g_)[i] ); 
	}
    }

  writer.value( "NOUT" , (int)(y_.size() - y.size()) );


  if ( t.size() == 0 ) return;

  //
  // scale 't' to be 0..1
  //

  double mnt = t[0] , mxt = t[0];

  for (int i=1;i<t.size();i++)
    {
      if ( t[i] < mnt ) mnt = t[i];
      if ( t[i] > mxt ) mxt = t[i];
    }

  std::vector<double> t01( t.size() );
  for (int i=0;i<t.size();i++)
    t01[i] = ( t[i] - mnt ) / ( mxt - mnt );



  //
  // Scale 'y' to be N(0,1)
  //

  std::vector<double> z = MiscMath::Z( y );
  

  // create actual dynam_t object

  dynam_t d( z , t01 );

  // requires at least 10 epochs 

  if ( d.size() < 10 ) return;
  
  writer.level( "UNSTRAT" , "EDYNAM" );
  
  double slope, rsq, h1, h2, h3;

  d.linear_trend( &slope, &rsq );
  d.hjorth( &h1, &h2, &h3 );
  
  writer.value( "N" , d.size() );
  writer.value( "SLOPE" , slope );
  writer.value( "RSQ" , rsq );
  writer.value( "H1" , h1 );
  writer.value( "H2" , h2 );
  writer.value( "H3" , h3 );
  
  writer.unlevel( "EDYNAM" );
  
  //
  // all done if no additional group-specification (e.g. sleep cycle)
  //

  if ( g_ == NULL ) return;


  //
  // make integer representaiton of group label
  //

  std::map<std::string,int> glabel;

  // 0 means skip; should be encoded as "" or "." by the requesting function
  
  glabel[ "" ] = 0;
  glabel[ "." ] = 0;

  int cnt = 1;
  
  std::vector<int> gint( g.size() );
  
  for (int i=0; i < g.size(); i++) 
    {
      std::map<std::string,int>::iterator ii =  glabel.find( g[i] );
      
      if ( ii == glabel.end() )
	{
	  glabel[ g[i] ] = cnt;
	  ++cnt;
	}
      
      gint[i] = glabel[ g[i] ];
      
    }

			       
  //
  // now we want to get within and between group (i.e. group-mean) based results
  //
  
  gdynam_t gdynam( gint , z , t );


  // assess strata and create the within-group dynam_t objects

  if ( gdynam.stratify() < 2 ) return;
  
  // between-group report : requires at least 3 groups

  if ( gdynam.stratify() >= 3 ) 
    {
      
      writer.level( "BETWEEN" , "EDYNAM" );
      
      double slope, rsq, h1, h2, h3;
      
      gdynam.between.linear_trend( &slope, &rsq );
      //    gdynam.between.hjorth( &h1, &h2, &h3 );
      
      writer.value( "N" , gdynam.between.size() );
      writer.value( "SLOPE" , slope );
      writer.value( "RSQ" , rsq );
      
      // these won't be meaningful typically...
      //     writer.value( "H1" , h1 );
      //     writer.value( "H2" , h2 );
      //     writer.value( "H3" , h3 );
      
      writer.unlevel( "BETWEEN" );
      
    }
  
  
  std::map<std::string,int>::const_iterator gg = glabel.begin();
  while ( gg != glabel.end() )
    {
      
      if ( gg->first == "." || gg->first == "" ) { ++gg; continue; }
      
      dynam_t & d = gdynam.within[ gg->second ];
      
      if ( d.size() >= 10 ) 
        {
	  
	  writer.level( gg->first , "EDYNAM" );
	  
	  double mean, var, slope, rsq, h1, h2, h3;
	  
	  d.mean_variance( &mean , &var );
	  d.linear_trend( &slope, &rsq );
	  d.hjorth( &h1, &h2, &h3 );
	  
	  writer.value( "N" , d.size() );
	  writer.value( "MEAN" , mean );
	  writer.value( "VAR" , var );
	  writer.value( "SLOPE" , slope );
	  writer.value( "RSQ" , rsq );
	  writer.value( "H1" , h1 );
	  writer.value( "H2" , h2 );
	  writer.value( "H3" , h3 );
	  
	  writer.unlevel( "EDYNAM" );
	}
      ++gg;
    }
  
}


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
  MiscMath::hjorth( &y , h1 , h2 , h3 );
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



//
// --------------------------------------------------------------------------------
// qdynam_t
//


qdynam_t::qdynam_t( const int ne , const std::vector<std::string> * pcycles )
  : ne(ne) 
{
  has_cycles = pcycles != NULL;
  
  if ( has_cycles )
    {      
      cycles = *pcycles;
      if ( cycles.size() != ne ) Helper::halt( "internal error in qdynam_t::qdynam_t()" );
    }

  winsor = -1;
  logscale = false;
  min_ne = 10;  // default limit (10 epochs, 5 mins)
  
  
  // default norm: only by min (set to 0) 
  norm01 = false;   // norm by min + max (dynam-norm-max=T)
  norm_mean = true; // default           (dynam-norm-mean=T) 
  // to only norm by min set dynam-norm-mean=F
  
  norm_each_section = true;  // for each within-section, do a separate norm
  
  // smoothing
  median_window = 19; // ~10 mins
  mean_window = 9;
  
  wcycles = false; // weight cycles by # epochs for BETWEEN  
  
  nq = 10; // default to 10 quantiles in Q strata
}


void qdynam_t::include( const std::vector<int> & e )
{
  // asumption: e contains elements from 0 to ne-1, when cycles if size ne
  //            and so e should be <= cycles.size() if cycles if set

  incl.resize( ne, false );

  for (int i=0;i<e.size();i++)
    {
      if ( e[i] < 0 || e[i] >= ne ) Helper::halt( "invalid epoch in qdynam_t" );
      incl[ e[i] ] = true;
    }
}

void qdynam_t::include( const std::vector<bool> & x )
{
  if ( x.size() != ne ) Helper::halt( "invalid epoch range in qdynam_t" );
  incl = x ; 
}

void qdynam_t::set_epochs( const std::vector<int> & e )
{  
  if ( ne != e.size() ) Helper::halt( "invalid epoch in qdynam_t" );
  epochs = e;
  // std::cerr << " setting " << epochs.size() << " epochs\n";
  // std::cerr << "e1 = " << epochs[0] << "\n";
}


void qdynam_t::set_max_cycles( const int n )
{
  if ( n < 1 ) return;
  // have max of 8 cycles for now
  incl_cycles.clear();
  for (int i=1; i<= (n > 8 ? 8 : n ) ; i++)
    incl_cycles.insert( "C" + Helper::int2str( i ) );
}

void qdynam_t::set_cycles( const std::vector<int> n )
{
  incl_cycles.clear();
  std::vector<int>::const_iterator ii = n.begin();
  while ( ii != n.end() )
    {
      if ( *ii >= 1 && *ii <= 8 ) 
	incl_cycles.insert( "C" + Helper::int2str( *ii ) );
      ++ii;
    }
}

void qdynam_t::winsorize( const double p )
{
  winsor = p;
}

void qdynam_t::proc( const std::vector<double> & x )
{

  // if not otherwise specified, include all epochs in analysis
  if ( incl.size() == 0 ) incl.resize( ne , true );

  bool has_mask = false;
  for (int i=0; i<ne; i++)
    {
      if ( ! incl[i] )
	{
	  has_mask = true;
	  break;
	}
    }
  
  // do we have epochs supplied (ne)
  const bool epochs_passed = epochs.size() == ne; 
  
  // calculate stats both assuming original epoch time and collapsed time
  //  (these will be the same if no include-mask has been set)
  
  std::vector<double> x1;
  std::vector<int> e1;
  std::vector<std::string> c1;

  if ( has_mask )
    {
      for (int i=0; i<ne; i++)
	{
	  if ( incl[i] )
	    {
	      // signals
	      x1.push_back( x[i] );

	      // cycles
	      if ( has_cycles ) c1.push_back( cycles[i] );

	      // epochs	      
	      e1.push_back( epochs_passed ? epochs[i] : i );
	    }
	}
    }
  else // if unmasked
    {
      // signal
      x1 = x;

      // cycles
      if ( has_cycles ) c1 = cycles;

      // epochs
      e1.resize( ne );
      for (int i=0; i<ne; i++) e1[i] = epochs_passed ? epochs[i] : i ;
    }

  // uniq cycles in this included set 
  // (will determine within-cycle, stratified outputs)
  const int nie = x1.size();
  
  // transform?
  if ( logscale )
    for (int i=0;i<nie; i++)
      x1[i] = log1p( x1[i] );
  
  // delineate cycles
  std::set<std::string> uniq_cycles; 
  if ( has_cycles ) 
    for (int i=0;i<nie; i++)
      if ( c1[i] != "." && c1[i] != "" )
	uniq_cycles.insert( c1[i] ) ;
  
  // winsorize?
  if ( winsor > 0 )
    MiscMath::winsorize( &x1 , winsor );

  // copy before norming (for original nq) 
  std::vector<double> ox1 = x1;
  
  
  //
  // stats
  //  

  //
  // 1) overall (QD = TOT)
  //

  const bool DO_SMOOTHING = true;
  const bool DO_NORMING = true;

  r1 = calc( x1 , e1 , DO_SMOOTHING, DO_NORMING );

  // store here for TOT, in case we want to output the
  // total smoothed series later
  r1_smoothed_series = ss;  
  r1_q10 = qnt( ss , nq );
  r1_os_q10 = qnt( os , nq ); 
  
  
  if ( ! has_cycles ) return;

  //
  // 2) stratified by 'cycle' 
  //

  std::vector<double> xc;
  std::vector<int> ec;

  int wtote = 0;
  
  std::set<std::string>::const_iterator cc = uniq_cycles.begin();
  while ( cc != uniq_cycles.end() )
    {

      // not including this cycle?

      if ( incl_cycles.size() != 0 )
	{
	  if ( incl_cycles.find( *cc ) == incl_cycles.end() )
	    {
	      ++cc;
	      continue;
	    }
	}

      // process
      
      std::vector<double> x2;
      std::vector<int> e2;

      for (int i=0; i<nie; i++)
	{
	  if ( c1[i] == *cc )
	    {
	      // use previously normed and smoothed values?
	      x2.push_back( norm_each_section ? ox1[i] : r1_smoothed_series[i] );
	      e2.push_back( e1[i] );
	    }
	}

      
      // only do if big enough
      if ( x2.size() >= min_ne )
	{
	  	  
	  // do calcs (sets 'os' and 'ss')
	  if ( norm_each_section ) 
	    rw[ *cc ] = calc( x2 , e2 , DO_SMOOTHING, DO_NORMING); 
	  else
	    rw[ *cc ] = calc( x2 , e2 , false , false );

	  // store
	  rw_smoothed_series[ *cc ] = ss;
	  rw_epochs[ *cc ] = e2;

	  // quantiles
	  rw_q10[ *cc ] = qnt( ss , nq );	  
	  rw_os_q10[ *cc ] = qnt( os , nq );
	  
	  // save means (for between-cycle stats)
	  // for between -cycle stats (based on means)
	  xc.push_back( rw[ *cc ].mean );	  
	  ec.push_back( MiscMath::mean( e2 ) );
	  
	  // for calculating average of within-cycle effects
	  rwa.ne++;
	  
	  int w = wcycles ? x2.size() : 1 ; // i.e. 1= no weighting

	  wtote += x2.size();
	  
	  rwa.sd     += w * rw[ *cc ].sd  ;
	  rwa.omean  += w * rw[ *cc ].omean ;
	  rwa.mean   += w * rw[ *cc ].mean ;
	  rwa.cv     += w * rw[ *cc ].cv;
	  rwa.tstat1 += w * rw[ *cc ].tstat1;
	  rwa.tstat2 += w * rw[ *cc ].tstat2;

	  rwa.corr1 += w * rw[ *cc ].corr1;
          rwa.corr2 += w * rw[ *cc ].corr2;

	  rwa.tmax += w * rw[ *cc ].tmax;
	  rwa.amax += w * rw[ *cc ].amax;
	  rwa.rmax += w * rw[ *cc ].rmax;

	  rwa.tmin += w * rw[ *cc ].tmin;
	  rwa.amin += w * rw[ *cc ].amin;
	  rwa.rmin += w * rw[ *cc ].rmin;

	  rwa.tminmax += w * rw[ *cc ].tminmax;
	  rwa.aminmax += w * rw[ *cc ].aminmax;
	  rwa.rminmax += w * rw[ *cc ].rminmax;

	}
      ++cc;
    }
  
  // between cycles (only makes sense if not norming within )
  if ( xc.size() > 1 && ! norm_each_section )
    {
      // no re-smoothing/norming needed here: (e.g. only a few data points, one per cycle)
      const bool NO_SMOOTHING = false; // i.e. F means do not do
      const bool NO_NORMING = false;
      rb = calc( xc , ec , NO_SMOOTHING , NO_SMOOTHING); 
    }

  //
  // average within cycle
  //
  
  if ( rwa.ne > 1 ) 
    {

      // either total number of epochs ( in weighted case)
      // or # of cycles (in unweighted case)
      
      double denom = wcycles ? wtote : rwa.ne ; 

      rwa.sd /= denom ;
      rwa.mean /= denom ;
      rwa.omean /= denom ;
      rwa.cv /= denom ;
      rwa.tstat1 /= denom ;
      rwa.tstat2 /= denom ;

      rwa.corr1 /= denom ;
      rwa.corr2 /= denom ;

      rwa.tmax /= denom ;
      rwa.amax /= denom ;      
      rwa.rmax /= denom ;

      rwa.tmin /= denom ;
      rwa.amin /= denom ;      
      rwa.rmin /= denom ;

      rwa.tminmax /= denom ;
      rwa.aminmax /= denom ;
      rwa.rminmax /= denom ;

      // set to -ve so we know it is # cycles, not # epochs
      // in output (otherwise we ignore some stats)
      // given we have 2+ cycles (i.e. something to avg over)
      rwa.ne = - rwa.ne; 
      
    }
  
}


qdynam_results_t qdynam_t::calc( const std::vector<double> & xx ,
				 const std::vector<int> & ee ,
				 const bool do_smoothing , 
				 const bool do_norming )
{
  
  // original (will be left as is)
  os = xx;

  // copy to smooth/norm and calculate all stats for
  ss = xx;
  
  
  // smooth? (if we're passing in the between-cycle series (e.g. may
  // only have 5-6 elements) we naturally don't want to smooth again,
  // thus the option to skip)
  
  if ( do_smoothing )
    ss = qdynam_t::smooth( ss , ee, median_window , mean_window );

  
  // norm?

  if ( do_norming )
    qdynam_t::norm( &ss , norm01 , norm_mean );


  // calculate T stat
  
  const int nn = ss.size();
  
  // grand total of time series
  double s_tot = 0; 

  // mean
  double s_mean = MiscMath::mean( ss );

  double sct = 0 , set = 0;

   
  double sct1 = 0 ;
  
  for (int i=0; i<nn; i++)
    {
      
      // sum signal, weighted by epoch numbera
      sct += ss[i] * ee[i] ;
      
      // same statistic, but if signal were completely uniform/flat 
      sct1 += s_mean * ee[i] ;      

      // same statistic, but weight by epoch order/rank rather than clock position
      set += ss[i] * i;
	          
      // get total 
      s_tot  += ss[i];

    }

  // 'clock-time' statistic (ct)
  const double sct_max = ee[nn-1] * s_tot;
  const double sct_min = ee[0] * s_tot ;
  const double stat_ct = ( sct - sct_min ) / ( sct_max - sct_min );

  // 'flat clock-time' statistic (ct1)
  const double sct1_max = ee[nn-1] * s_tot;
  const double sct1_min = ee[0] * s_tot ;
  const double stat_ct1 = ( sct1 - sct1_min ) / ( sct1_max - sct1_min );

  // epoch order/rank statistic (et)
  const double set_max = (nn-1) * s_tot;
  const double set_min = 0 * s_tot ;
  const double stat_et = ( set - set_min ) / ( set_max - set_min );

  // all above statistics scaled between min/max and so [ 0 , 1 ] range
  // scale to [ -100 , +100 ] when returning

  // clock-time statistic is adjusted by the 'expectation' under a completely flat
  // set of data-points

  // for rank-based statistic, we don't need to do this, as we know that would be 0
  // by definition
  
  // return   100 * ( (2S)-1 ) 

  qdynam_results_t r;
  r.ne = nn;

  // clock time statistic
  r.tstat1 = 100 * ( ( stat_ct * 2.0 ) - 1.0 ) ;

  // collapsed value
  r.tstat2 = 100 * (  ( stat_et * 2.0 ) - 1.0 ) ; 

  // adjust clock-time stat by expectation under flatness 
  double tstat11 = 100 * ( ( stat_ct1 * 2.0 ) - 1.0 ) ;

  // adust clock-time by expectation under flat data
  r.tstat1 -= tstat11;

  //
  // simple corrs (duh...)
  //
  std::vector<double> e1( nn ), e2( nn );;
  for (int i=0; i<nn; i++) { e1[i] = i; e2[i] = ee[i]; }  
  r.corr1 = Statistics::correlation( ss , e1 );
  r.corr2 = Statistics::correlation( ss , e2 );
  
  
  //
  // basics
  //

  // mean of original (unnormed) time-series
  r.omean = MiscMath::mean( os );
  
  // stats for smoothed, normed series
  r.sd = MiscMath::sdev( ss );
  r.mean = MiscMath::mean( ss );
  r.cv = r.sd / r.mean;
  
  // max/min slope stats
 
  double ss_min = ss[0];
  double ss_max = ss[0];
  int ss_max_i = 0;
  int ss_min_i = 0;

  for (int i=1; i<nn; i++)
    {
      if ( ss[i] < ss_min ) // takes first even if tied (from winsorizing)
	{
	  ss_min = ss[i];
	  ss_min_i = i;
	}
      
      if ( ss[i] > ss_max )
	{
	  ss_max = ss[i];
	  ss_max_i = i;
	}
    }

  
  r.tmax = ee[ ss_max_i ] - ee[ 0 ] ; // use real epoch counts
  r.amax = ss_max - ss[ 0 ] ; 
  r.rmax = r.amax / ( r.tmax + 1 ); // i.e. if max is epoch 0

  r.tmin = ee[ ss_min_i ] - ee[ 0 ] ; // use real epoch counts
  r.amin = ss_min - ss[ 0 ] ; // as above, just make -ve
  r.rmin = r.amin / ( r.tmin + 1 ); // i.e. if max is epoch 0

  r.aminmax = ss_max - ss_min ;
  r.tminmax = ee[ ss_max_i ] - ee[ ss_min_i ];  // define as max-to-min (+ve, max last)
  r.rminmax = r.aminmax / ( r.tminmax == 0 ? 1 : r.tminmax ) ; // in case flat sig
  
  return r;
}    



void qdynam_t::output_helper( const qdynam_results_t & res , const bool verbose , const bool between )
{

  // handle WITHIN case where -ve means # cycles
  writer.value( "N" , res.ne < 0 ? - res.ne : res.ne );
  
  if ( ! between ) 
    writer.value( "OMEAN" , res.omean ); 
  
  writer.value( "MEAN" , res.mean );
  writer.value( "SD" , res.sd );
  writer.value( "T" , res.tstat1 );
  writer.value( "R" , res.corr1 );
    
  if ( verbose )
    {
      writer.value( "CV" , res.cv );
      writer.value( "TR" , res.tstat2 );
      writer.value( "RR" , res.corr2 );
    }

  if ( res.ne > 10 || res.ne < 0 ) // if -ve means WITHIN, # cycles
    {
      
      writer.value( "T_P2P" , res.tminmax );
      writer.value( "A_P2P" , res.aminmax );
      
      if ( verbose )
	{
	  
	  writer.value( "AT_P2P" , res.rminmax );	  
	  
	  writer.value( "T_MX" , res.tmax );
	  writer.value( "A_MX" , res.amax );	  
	  writer.value( "AT_MX" , res.rmax );
	  
	  writer.value( "T_MN" , res.tmin );
	  writer.value( "A_MN" , res.amin );
	  writer.value( "AT_MN" , res.rmin );
      
	}
    }
}

std::vector<double> qdynam_t::qnt( const std::vector<double> & x , const int nq )
{
  std::vector<double> q( nq );
  const int n = x.size();
  const double s = n / (double)nq;

  double w = 0;

  for (int i=0;i<nq; i++)
    {

      // span from 'w' to 'w+s'
      //  if 'w' is fractional, include 1 - f of w
      double w2 = w + s;
      double t = 0 ;

      //std::cout << " going " << w << " to " << w2 << "\n";
      
      while ( 1 )
	{
	  int w1 = w;
	  double f = w - w1;
	  
	  // (fractional) part of first block 
	  t += x[w1] * ( 1 - f );

	  //std::cout << " adding " << w1 << " wgt " << (1-f) << "\n";

	  // next 
	  // (fractional, potentially zero, part of next)
	  w += ( 1 - f );
	  	  
	  if ( w2 - w > 1 ) f = 1 - f;
	  else f = w2 - w; 

	  if ( f > 0 && w1  < n )
	    {
	      t += x[ w1 ] * f;	      
	      //std::cout << " a(2nd) " << w1+1 << " wgt " << f << "\n";
	    }
	  
	  // slide one however much chunked off
	  w += f;
	  
	  // done
	  if ( w >= w2 || fabs( w - w2 ) < 1e-4 ) break;
	}
      //std::cout << " done\n";
      
      q[ i ] = t / (double) s;
      
    }
  
  return q;
}


std::vector<double> qdynam_t::smooth( const std::vector<double> & x ,
				      const std::vector<int> & e ,
				      const int w1, const int w2 )
{

  // nothing to do?
  if ( ! ( w1 > 1 || w2 > 1 ) ) return x;


  
  // to avoid bad smoothing over gaps, here we take the epoch count
  // as well and expand the series first, linear interpolation between
  // first and last (N) points, then smooth, then splice out the
  // desired components... should help to reduce edge effects
  
  if ( x.size() != e.size() ) Helper::halt( "internal logic error (1) in smooth()" );  
  const int n = x.size();

  // too small
  if ( n < w1 || n < w2 ) return x;
  
  const bool debug = false;

  if ( debug )
    for (int i=0; i<n; i++)
      std::cout << " ---> " << i << "\t" << e[i] << "\t" << x[i] << "\n";
  
  // *assume* sorted; +1 as e[] is 0-based  
  // all epochs from first observed to last observed
  int n2 = e[n-1] - e[0] + 1;

  int e0 = e[0]; // first element
  
  if ( debug ) 
    std::cout << " n = " << n << "\n"
	      << " e0 = " << e0 << "\n"
	      << " en = " << e[n-1] << "\n"
	      << " n2 = " << n2 << "\n";
  
  
  // create full time-series : which goes from e0 to e[n-1] (incl) w/ n2 elements
  std::vector<double> x2( n2 , 0 );
  std::vector<double> e2( n2 , 0 );
  std::vector<bool> fill( n2 , true );
  
  // nb. here we need to adjust for e0
  for (int i=0; i<n; i++)
    {
      x2[ e[i] - e0 ] = x[i] ;
      e2[ e[i] - e0 ] = e[i] ; 
      fill[ e[i] - e0 ] = false ;
    }

  // track original gaps
  std::vector<bool> fill_orig = fill;
  
  if ( debug )
    {
      for (int i=0; i<n2; i++)
	{
	  std::cout << "orig" << i << "\t"
		    << fill[i] << "\t"
		    << e2[i] << "\t"
		    << x2[i] << "\n";
	}
      std::cout << "\n\n";
    }
  
  // linear interpolate over gaps  
  for (int i=0; i<n2; i++)
    {
      
      // have we come across a new gap?
      // nb. here want to use fill[] which gets modified as we go
      if ( fill[i] )
	{
	  
	  if ( debug )
	    std::cout << " found gap starting " << i << "\n";
	  
	  std::vector<double> earlier;
	  std::vector<double> later;

	  // go back until data
	  int p = i;
	  while ( 1 )
	    {
	      // unlikely, but check in case
	      if ( p == 0 ) break; 
	      
	      // go back
	      --p;
	      
	      // if encounter another gap, stop
	      // n.b. need to look up original gap status here
	      // only want to fill w/ observed vals
	      if ( fill_orig[p] && earlier.size() != 0 ) break;
	      
	      // if non-gap, take values, up to mx of 3
	      if ( ! fill_orig[p] ) {

		if ( debug )
		  std::cout << "  adding earlier point " << p << "\n";

		earlier.push_back( x2[p] );
	      }
	      if ( earlier.size() == 3 ) break;	  
	    }
	  
	  // go forward 
	  p = i;	  
	  while ( 1 )
	    {
	      ++p;
	      
	      // unlikely, but check in case
	      if ( p == n2 ) break; 
	      
	      // if encounter another gap, stop
	      if ( fill_orig[p] && later.size() != 0 ) break;
	      
	      // if non-gap, take values, up to mx of 3
	      if ( ! fill_orig[p] )
		{
		  if ( debug ) std::cout << "  adding later point " << p << "\n";
		  later.push_back( x2[p] );
		}
	      if ( later.size() == 3 ) break; 
	    }

	  // should always have /something/ in earlier[] and later[]
	  
	  if ( earlier.size() == 0 || later.size() == 0 )
	    Helper::halt( "internal logic error (2) in smooth() interpolation" );
	  
	  const double emean = MiscMath::mean( earlier );
	  const double lmean = MiscMath::mean( later );

	  if ( debug )
	    std::cout << " el = " << earlier.size() << " " << later.size() << " " << emean << " " << lmean << "\n";
	  
	  // we'll be hitting this the first epoch of a gap;
	  // fill all here
	  // i.e. should always be the case the i-1 is a non-gap
	  if ( i == 0 || fill[i-1] ) Helper::halt( "internal logic error (3) in smooth()" );

	  // find size of gap
	  p = i;
	  while ( 1 )
	    {
	      ++p;
	      if ( p == n2 ) break;
	      if ( ! fill[p] ) break;
	    }
	  // i is start of gap
	  // p is now one past
	  int gap_size = p - i; 
	  
	  if ( debug ) std::cout << " gap size = " << gap_size << "\n";
	  
	  const double gradient = ( lmean - emean ) / (double)(gap_size+1);
	  
	  p = i;
	  for (int j=1; j<=gap_size; j++)
	    {

	      if ( debug ) std::cout << " p = " << p << " " << fill[p] << "\n";

	      // should not already be filled
	      if ( ! fill[p] ) Helper::halt( "internal logic error (4) in smooth()");

	      // fill the gap now
	      x2[p] = emean + j * gradient;

	      // note that we've filled this, so won't be considered
	      fill[p] = false; 
	      
	      ++p;
	    }	  
	}
    }

  std::vector<double> r2;
  
  // initial median filter
  if ( w1 > 1 )
    {
      if ( w2 > 1 )
	r2 = MiscMath::moving_average( MiscMath::median_filter( x2 , w1 ) , w2 );
      else
	r2 = MiscMath::median_filter( x2 , w1 );
    }
  else
    {
      if ( w2 > 1 )
	r2 = MiscMath::moving_average( x2 , w2 );
      else
        r2 = x2 ; // should have been dealt w/ above 
    }


  //
  // unsplice
  //
  
  std::vector<double> r( n );

  int p = 0;
  for (int i=0; i<n2; i++)
    {
      if ( ! fill_orig[i] ) // based on the original fill tracking
	{
	  if ( p == n ) Helper::halt( "internal logic error (5) in smooth()" );
	  r[p] = r2[i];
	  ++p;	  
	}
    }


  if ( debug ) 
    {
      std::cout << " FINAL\n";
      for (int i=0; i<n2; i++)
	{
	  std::cout << i << "\t"
		    << fill_orig[i] << "\t"
		    << fill[i] << "\t"
	     		<< x2[i] << "\t"
		    << r2[i] << "\n";
	}
      std::cout << "DONE\n\n";
    }
  
  
  //
  // all done, return;
  //

  return r;
  
}

void qdynam_t::norm( std::vector<double> * x , const bool do_max , const bool do_mean )
{
  
  const int n = x->size();
  
  // 1) always constrain x to be positive
  // 2)  do_max = F , do_mean = F  : leave rest of signal as is
  //     do_max = T , do_mean = F  : scale so max = 1 
  //     do_max = F , do_mean = T  : scale so mean = 1 

  double xmin = (*x)[0];
  double xmax = (*x)[0];

  // get min/max
  for (int i=1;i<n;i++)
    {
      if ( (*x)[i] < xmin ) xmin = (*x)[i];
      if ( do_max && (*x)[i] > xmax ) xmax = (*x)[i];
    }
  
  // set min to 0.0
  for (int i=0;i<n; i++)
    (*x)[i] -= xmin;

  // set 1.0 == max
  if ( do_max ) 
    {
      xmax -= xmin;
      for (int i=0;i<n;i++) (*x)[i] /= xmax;	
    }

  // set 1.0 == mean
  else if ( do_mean )
    {
      const double xmean = MiscMath::mean( *x );
      for (int i=0;i<n;i++) (*x)[i] /= xmean;
    }

}
