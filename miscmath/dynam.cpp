
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

bool dynam_compile_cycles( edf_t & edf , std::vector<std::string> * cycles, std::vector<int> * epochs  )
{

  // data must be epoched already (i.e. will be if cycles done)
  if ( !edf.timeline.epoched() ) return false;
  
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
  
  std::vector<int> _epochs;
  std::vector<std::string> _cycles;
  
  edf.timeline.first_epoch();

  while ( 1 ) 
    {
      
      int epoch = edf.timeline.next_epoch();      
      
      if ( epoch == -1 ) break;
      
      // cycle annot
      
      std::string c = "."; // null
      
      // nb. uses legacy epoch-annotation encoding
      // take up to 10 cycles
      
      if      ( edf.timeline.epoch_annotation( "_NREMC_1" , epoch ) ) c = "C1";
      else if ( edf.timeline.epoch_annotation( "_NREMC_2" , epoch ) ) c = "C2";
      else if ( edf.timeline.epoch_annotation( "_NREMC_3" , epoch ) ) c = "C3";
      else if ( edf.timeline.epoch_annotation( "_NREMC_4" , epoch ) ) c = "C4";
      else if ( edf.timeline.epoch_annotation( "_NREMC_5" , epoch ) ) c = "C5";
      else if ( edf.timeline.epoch_annotation( "_NREMC_6" , epoch ) ) c = "C6";
      else if ( edf.timeline.epoch_annotation( "_NREMC_7" , epoch ) ) c = "C7";
      else if ( edf.timeline.epoch_annotation( "_NREMC_8" , epoch ) ) c = "C8";
      else if ( edf.timeline.epoch_annotation( "_NREMC_9" , epoch ) ) c = "C9";
      else if ( edf.timeline.epoch_annotation( "_NREMC_10" , epoch ) ) c = "C10";
      
      _cycles.push_back( c );
      _epochs.push_back( epoch );

      // next epoch
    }

  if ( cycles ) *cycles = _cycles;
  if ( epochs ) *epochs = _epochs;

  return true;
}


// wrapper(s)

void dynam_report_with_log( param_t & param,
			    const std::vector<double> & y , 
			    const std::vector<double> & t , 
			    const std::vector<std::string> * g )
{
  std::vector<double> yl( y.size() );
  for (int i=0;i<y.size();i++) yl[i] = log( y[i] );
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
  
  const double qd_winsor = param.has( "dynam-winsor" ) ? param.requires_dbl( "dynam-winsor" ) : 0.02 ;  
  qd.winsorize( qd_winsor );

  if ( param.has( "dynam-median-window" ) )
    qd.set_smoothing_median_window( param.requires_int( "dynam-median-window" ) );
  
  if ( param.has( "dynam-mean-window" ) )
    qd.set_smoothing_mean_window( param.requires_int( "dynam-mean-window" ) );

  
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
      // between 
      writer.level( "BETWEEN" , "QD" );
      qdynam_t::output_helper( qd.rb , verbose );
      
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
	      
	      for (int i=0; i<ss.size(); i++)
		{
		  writer.epoch( i + 1 ); // 1-based outputs     
		  writer.value( "SS" , ss[i] );
		}
	      writer.unepoch();
	      
	      ++cc;
	    }
	}
      
      writer.unlevel( "QD" );
      
    }
      
  
  
  //
  // (original) method 2
  //
  
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
  median_window = 19; // ~10 mins
  mean_window = 9;
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
      uniq_cycles.insert( c1[i] ) ;
  
  // winsorize?
  if ( winsor > 0 )
    MiscMath::winsorize( &x1 , winsor );

  // normalize x1 to be positive
  double xmin = x1[0];

  for (int i=1;i<nie;i++)
    if ( x1[i] < xmin ) xmin = x1[i];
  
  for (int i=0;i<nie; i++)
    x1[i] -= xmin;
  
  // std::cout << "\n\n----\n";

  // for (int i=0;i<nie; i++)
  //   std::cout << "det\t" 
  //  	      << x1[i] << "\t"
  //  	      << c1[i] << "\t"
  //  	      << e1[i] << "\n";
  
  
  //
  // stats
  //  

  // 1) overall

  r1 = calc( x1 , e1 );

  // store here for TOT, in case we want to output the
  // total smoothed series later
  r1_smoothed_series = ss;
  
  if ( ! has_cycles ) return;
    
  // 2) stratified by 'cycle' 

  std::vector<double> xc;
  std::vector<int> ec;

  std::set<std::string>::const_iterator cc = uniq_cycles.begin();
  while ( cc != uniq_cycles.end() )
    {
      std::vector<double> x2;
      std::vector<int> e2;

      for (int i=0; i<nie; i++)
	{
	  if ( c1[i] == *cc )
	    {
	      x2.push_back( x1[i] );
	      e2.push_back( e1[i] );
	    }
	}

      // do calces
      rw[ *cc ] = calc( x2 , e2 ); 

      // store
      rw_smoothed_series[ *cc ] = ss;
      
      // save means (for between-cycle stats)
      // if big enough
      if ( x2.size() >= min_ne )
	{
	  // for between -cycle stats (based on means)
	  xc.push_back( rw[ *cc ].mean );
	  ec.push_back( MiscMath::mean( e2 ) );
	  
	  // for calculating average of within-cycle effects
	  rwa.ne++;
	  rwa.sd += rw[ *cc ].sd;
	  rwa.mean += rw[ *cc ].mean;
	  rwa.cv += rw[ *cc ].cv;
	  rwa.tstat1 += rw[ *cc ].tstat1;
	  rwa.tstat2 += rw[ *cc ].tstat2;

	  rwa.tmax += rw[ *cc ].tmax;
	  rwa.amax += rw[ *cc ].amax;
	  rwa.lmax += rw[ *cc ].lmax;
	  rwa.rmax += rw[ *cc ].rmax;

	  rwa.tmin += rw[ *cc ].tmin;
	  rwa.amin += rw[ *cc ].amin;
	  rwa.lmin += rw[ *cc ].lmin;
	  rwa.rmin += rw[ *cc ].rmin;

	  rwa.tminmax += rw[ *cc ].tminmax;
	  rwa.aminmax += rw[ *cc ].aminmax;
	  rwa.lminmax += rw[ *cc ].lminmax;
	  rwa.rminmax += rw[ *cc ].rminmax;

	}
      ++cc;
    }
  
  // between cycles
  if ( xc.size() > 1 )
    {
      const bool SKIP_SMOOTHING = true;
      rb = calc( xc , ec , SKIP_SMOOTHING );
    }

  // average within cycle
  if ( rwa.ne > 1 ) 
    {
      rwa.sd /= (double) rwa.ne ;
      rwa.mean /= (double) rwa.ne ;
      rwa.cv /= (double) rwa.ne ;
      rwa.tstat1 /= (double) rwa.ne ;
      rwa.tstat2 /= (double) rwa.ne ;

      rwa.tmax /= (double) rwa.ne ;
      rwa.amax /= (double) rwa.ne ;
      rwa.lmax /= (double) rwa.ne ;
      rwa.rmax /= (double) rwa.ne ;

      rwa.tmin /= (double) rwa.ne ;
      rwa.amin /= (double) rwa.ne ;
      rwa.lmin /= (double) rwa.ne ;
      rwa.rmin /= (double) rwa.ne ;

      rwa.tminmax /= (double) rwa.ne ;
      rwa.aminmax /= (double) rwa.ne ;
      rwa.lminmax /= (double) rwa.ne ;
      rwa.rminmax /= (double) rwa.ne ;

    }
  
}

qdynam_results_t qdynam_t::calc( const std::vector<double> & xx ,
				 const std::vector<int> & ee ,
				 const bool skip_smoothing )
{

  // all analysis/stats are reported for the smoothed time-series
  //  note - this does not respect discontinuities, but general
  //         principals should be fine (e.g. these are done on winsorized
  //         time series, and within-cycle analyses and between-cycle analysis
  //         should be quite robust, i.e. not many large gaps) 


  // if we're passing in the between-cycle series (e.g. may only have 5-6 elements) we naturally
  // don't want to smooth again, thus the option to skip
  
  ss = xx;
  
  if ( ! skip_smoothing )
    {
      // initial median filter
      if ( median_window > 1 ) 
	ss = MiscMath::median_filter( ss , median_window );
      
      // secondary moving average filter
      ss = MiscMath::moving_average( ss , mean_window );
    }
  
  // 1) scale: min = 0, max = whatever
  // 2) t : either scale 0..(ne-1)
  //        or (first_incl)..(last_incl)

  const int nn = ss.size();
    
  double sct = 0 , set = 0;
  double sk = 0;

  double sct1 = 0 ;
  double sk1 = 0;

  double mss = MiscMath::mean( ss );
  
  for (int i=0; i<nn; i++)
    {
      sct  += ss[i]   * ee[i] ;
      set += ss[i] * i;
      sk  += ss[i];

      sct1 += mss     * ee[i] ;
      sk1 += mss;
      
    }

  const double sct_max = ee[nn-1] * sk;
  const double sct_min = ee[0] * sk ;
  const double stat_ct = ( sct - sct_min ) / ( sct_max - sct_min );
  
  const double sct1_max = ee[nn-1] * sk1;
  const double sct1_min = ee[0] * sk1 ;
  const double stat_ct1 = ( sct1 - sct1_min ) / ( sct1_max - sct1_min );
  
  const double set_max = (nn-1) * sk;
  const double set_min = 0 * sk ;
  const double stat_et = ( set - set_min ) / ( set_max - set_min );

  
  // return [0.1] --> [-100 , +100]
  
  qdynam_results_t r;
  r.ne = nn;

  r.tstat1 = 100 * ( ( stat_ct * 2.0 ) - 1.0 ) ;

  double tstat11 = 100 * ( ( stat_ct1 * 2.0 ) - 1.0 ) ;
  // adust clock-time by expectation under flat data
  r.tstat1 -= tstat11;

  // collapsed value
  r.tstat2 = 100 * (  ( stat_et * 2.0 ) - 1.0 ) ; 
  
  //
  // basics
  //

  r.sd = MiscMath::sdev( ss );
  r.mean = MiscMath::mean( ss );
  r.cv = r.sd / r.mean;

  //
  // max/min stats
  //
 
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
  r.lmax = r.amax * r.tmax;
  r.rmax = r.amax / ( r.tmax + 1 ); // i.e. if max is epoch 0

  r.tmin = ee[ ss_min_i ] - ee[ 0 ] ; // use real epoch counts
  r.amin = ss_min - ss[ 0 ] ; // as above, just make -ve
  r.lmin = r.amin * r.tmin;
  r.rmin = r.amin / ( r.tmin + 1 ); // i.e. if max is epoch 0

  r.aminmax = ss_max - ss_min ;
  r.tminmax = ee[ ss_min_i ] - ee[ ss_max_i ];
  r.lminmax = r.aminmax * r.tminmax ;
  r.rminmax = r.aminmax / r.tminmax ;

  // std::cout << " ss_min_i = " << ss_min_i << "\n"
  // 	    << " ee[ ss_min_i ] = " << ee[ ss_min_i ] << "\n"
  // 	    << " ee[ 0 ] = " << ee[ 0 ] << "\n"
  // 	    << " ss_min = " << ss_min << "\n"
  // 	    << " ss[ 0 ] = " << ss[ 0 ] <<  "\n"
  // 	    << " r.amin = " << r.amin << "\n"
  // 	    << " r.tmin = " << r.tmin << "\n";
  // for (int i=0;i<nn;i++)
  //   std::cout << "det " << i << " " << ss[i] << "\n";
  // std::cout << "\n\n";
  
  return r;
}    


void qdynam_t::output_helper( const qdynam_results_t & res , const bool verbose )
{

  writer.value( "N" , res.ne );
  writer.value( "MEAN" , res.mean );
  writer.value( "SD" , res.sd );
  writer.value( "T" , res.tstat2 );
  
  if ( verbose )
    {
      writer.value( "CV" , res.cv );
      writer.value( "T1" , res.tstat1 );
    }

  if ( res.ne > 10 )
    {
      writer.value( "MAX_T" , res.tmax );
      writer.value( "MAX_A" , res.amax );
      
      writer.value( "MIN_T" , res.tmin );
      writer.value( "MIN_A" , res.amin );

      writer.value( "MINMAX_T" , res.tminmax );
      writer.value( "MINMAX_A" , res.aminmax );

      if ( verbose )
	{
	  writer.value( "MAX_PROD" , res.lmax );
	  writer.value( "MAX_RATIO" , res.rmax );
      
	  writer.value( "MIN_PROD" , res.lmin );
	  writer.value( "MIN_RATIO" , res.rmin );
      
	  writer.value( "MINMAX_PROD" , res.lminmax );
	  writer.value( "MINMAX_RATIO" , res.rminmax );
	}
    }
}
