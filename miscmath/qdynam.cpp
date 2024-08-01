
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

#include "miscmath/qdynam.h"
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


//
// extract NREM cycle annotations for current epoch set
//

bool qdynam_t::dynam_compile_cycles( edf_t & edf )
{

  // data must have already been epoched (by HYPNO, so that
  // we expected NREM cycle epoch-annotations)
  if ( ! edf.timeline.epoched() )
    Helper::halt( "data not epoched: run HYPNO before dynam-submodules are invoked" );      
  
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
  
  // generate current epoch code
  epochs.clear();
  uepochs.clear();
  
  int ne = edf.timeline.first_epoch();
  while ( 1 )
    {
      int epoch = edf.timeline.next_epoch();
      if ( epoch == -1 ) break;
      
      // save display epoch number - 1 
      const int disp_epoch = edf.timeline.display_epoch( epoch ) - 1;
      epochs.push_back( disp_epoch );
      uepochs.insert( disp_epoch );
      
    }  
  
  // construct by iterating over current epoch set
  cycles.clear();
  
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
      
      cycles[ epochs[e] ] = c ;
      // next epoch
    }
  
  return true;
}


void qdynam_t::add( const std::map<std::string,std::string> & faclvl,
		    const std::string & metric,
		    const int epoch , 
		    const double value )
{

  // retain original order (noting it may be passed w/ multiple vars)
  // (var output will be alphabetical)
  if ( sequences.find( faclvl ) == sequences.end() )
    osequences.push_back( faclvl );

  // test valid epoch
  if ( uepochs.find( epoch ) == uepochs.end() )
    Helper::halt( "undefined epoch, internal error in qdynam_t::add()" );
  
  // save value
  sequences[ faclvl ][ metric ][ epoch ] = value; 

}
		    

void qdynam_t::proc_all()
{

  logger << "  running dynam submodule for " << osequences.size() << " distinct strata\n";
  
  // iterate over each fac/lvl/var in the store
  
  std::vector<std::map<std::string,std::string> >::const_iterator os = osequences.begin();
  while ( os != osequences.end() )
    {
      
      std::map<std::map<std::string,std::string>,std::map<std::string,std::map<int,double> > >::const_iterator ss = sequences.find( *os );
      if ( ss == sequences.end() ) Helper::halt( "internal error in qdynam_t::proc_all()" );
      
      const std::map<std::string,std::string> & faclvl = ss->first;
      const std::map<std::string,std::map<int,double> > & vars = ss->second;
      std::map<std::string,std::map<int,double> >::const_iterator vv = vars.begin();
      while ( vv != vars.end() )
	{
	  const std::string & var = vv->first;
	  const std::map<int,double> & data = vv->second;
	  const int ne = data.size();

	  // debug output
	  if ( 0 )
	    {
	      std::cout << " strata ne = " << ne << "\n";
	      // dump output factors
	      std::map<std::string,std::string>::const_iterator fff = faclvl.begin();
	      while ( fff != faclvl.end() )
		{
		  std::cout << "  " << fff->first << " -> " << fff->second << "\n";
		  ++fff;
		}
	    }
	  
	  // pull out relevant epochs (may be a subset)
	  std::vector<double> data1( ne ) ;
	  std::vector<int> epochs1( ne );
	  std::vector<std::string> cycles1( ne );
	  int idx = 0;
	  std::map<int,double>::const_iterator ee = data.begin();
	  while ( ee != data.end() )
	    {
	      epochs1[ idx ] = ee->first;
	      data1[ idx ] = ee->second;
	      cycles1[ idx ] = cycles[ ee->first ];
	      ++idx;
	      ++ee;
	    }

	  
	  //
	  // initiate the run
	  //

	  proc( data1 , epochs1 , cycles1 );	  
	  

	  // ------------------------------------------------------------
	  //
	  // output 
	  //

	  // set up faclvl info for tracked measures
	  writer.level( var , globals::var_strat );
	  
	  // dump output factors
	  std::map<std::string,std::string>::const_iterator fff = faclvl.begin();
	  while ( fff != faclvl.end() )
	    {
	      writer.level( fff->second , fff->first );
	      //std::cout << "  " << fff->first << " -> " << fff->second << "\n";
	      ++fff;
	    }
	  

	  // now start reporting different QD 

	  writer.level( "TOT" , "QD" );

	  qdynam_t::output_helper( r1 , verbose );

	  if ( has_cycles )
	    {
	      // between (only if not norming within each cycle)
	      if ( ! norm_cycles() ) 
		{
		  writer.level( "BETWEEN" , "QD" );
		  qdynam_t::output_helper( rb , verbose , true );
		}
	      
	      // average within
	      writer.level( "WITHIN" , "QD" );
	      qdynam_t::output_helper( rwa , verbose );
	      
	      // each cycle
	      std::map<std::string,qdynam_results_t> & cycs = rw;
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
	  const std::vector<double> & ss = r1_q10;
	  const std::vector<double> & os = r1_os_q10;
	  
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
	      
	      std::map<std::string,std::vector<double> > & cycs = rw_q10;
	      std::map<std::string,std::vector<double> >::const_iterator cc = cycs.begin();
	      
	      while ( cc != cycs.end() )
		{
		  writer.level( "W_" + cc->first , "QD" );
		  
		  const std::vector<double> & ss = cc->second;
		  const std::vector<double> & os = rw_os_q10[ cc->first ];
		  
		  for (int i=0; i<ss.size(); i++)
		    {
		      writer.level( i + 1 , "Q" ); 
		      writer.value( "SS" , ss[i] );
		      
		      // do not show if not norming each section
		      // (just the way things are calculated internally, we don't get this)
		      if ( norm_cycles() )
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
	      const std::vector<double> & ss = r1_smoothed_series;
	      
	      if ( epochs.size() == ss.size() )
		{
		  // epochs contains display_epoch() - 1
		  
		  writer.level( "TOT" , "QD" );
		  for (int i=0; i<ss.size(); i++)
		    {
		      writer.epoch( epochs[i] + 1 ); // 1-based outputs
		      writer.value( "SS" , ss[i] );
		    }
		  writer.unepoch();
		}
	      

	      // cycles
	      if ( has_cycles )
		{
		  
		  std::map<std::string,std::vector<double> > & cycs = rw_smoothed_series;
		  std::map<std::string,std::vector<double> >::const_iterator cc = cycs.begin();
		  
		  while ( cc != cycs.end() )
		    {
		      writer.level( "W_" + cc->first , "QD" );
		      
		      const std::vector<double> & ss = cc->second;
		      const std::vector<int> & ee = rw_epochs[ cc->first ];
		      
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
	  // close out factors
	  //

	  writer.unlevel( globals::var_strat );

	  fff = faclvl.begin();
          while ( fff != faclvl.end() )
            {
              writer.unlevel( fff->first );
              ++fff;
            }

	  
	  
	  //
	  // next varaiable
	  //
	  
	  ++vv;
	}

      //
      // next strata
      //

      ++os;
    }

  //
  // all done
  //
  
}






// --------------------------------------------------------------------------------
// qdynam_t
//


void qdynam_t::init( edf_t & edf , param_t & param )
{

  //
  // default values
  //

  winsor = -1;

  logscale = false;

  min_ne = 10;  // default limit (10 epochs, 5 mins)

  trim_epochs.resize(2,0);
  
  norm01 = false;   // norm by min + max (dynam-norm-max=T)
  
  norm_mean = true; // default           (dynam-norm-mean=T) 
  // to only norm by min set dynam-norm-mean=F

  // for each within-section, do a separate norm
  norm_each_section = true;  
  
  // smoothing
  median_window = 19; // ~10 mins
  mean_window = 9;

  // weight cycles by # epochs for WITHIN
  wcycles = true; 

  // default to 10 quantiles in Q strata
  nq = 10; 


  //
  // options
  //
  
  verbose = param.has( "dynam-verbose" );

  epoch_output = param.has( "dynam-epoch" );

  if ( param.has( "dynam-min-ne" ) )
    set_min_ne( param.requires_int( "dynam-min-ne" ) );

  if ( param.has( "dynam-trim-epochs" ) )
    {
      std::vector<int> x = param.intvector( "dynam-trim-epochs" );
      if ( x.size() == 1 )
	{
	  int xx = x[0];
	  x.resize( 2 );
	  x[0] = x[1] = xx;
	}
      if ( x.size() == 2 ) 
	set_trim_epochs( x );
    }
  
  const double qd_winsor = param.has( "dynam-winsor" ) ? param.requires_dbl( "dynam-winsor" ) : 0.05 ;  
  winsorize( qd_winsor );

  if ( param.has( "dynam-median-window" ) )
    set_smoothing_median_window( param.requires_int( "dynam-median-window" ) );
  
  if ( param.has( "dynam-mean-window" ) )
    set_smoothing_mean_window( param.requires_int( "dynam-mean-window" ) );
  
  if ( param.has( "dynam-norm-mean" ) )
    set_norm_mean( param.yesno( "dynam-norm-mean" ) );
  else if ( param.has( "dynam-norm-max" ) )
    set_norm_max( param.yesno( "dynam-norm-max" ) );
  
  if ( param.has( "dynam-norm-cycles" ) )
    set_norm_cycles( param.yesno( "dynam-norm-cycles" ) );
  
  if ( param.has( "dynam-max-cycle" ) )
    set_max_cycles( param.requires_int( "dynam-max-cycle" ) );
  else if ( param.has( "dynam-cycles" ) )
    set_cycles( param.intvector( "dynam-cycles" ) );

  // default false
  if ( param.has( "dynam-weight-cycles" ) )
    set_weight_cycles( param.yesno( "dynam-weight-cycles" ) );


  //
  // compile cycles
  //

  has_cycles = dynam_compile_cycles( edf );
  
  
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


void qdynam_t::reinit()
{
  r1.init();
  rb.init();
  rwa.init();
  rw.clear();

  // smoothed & normed series
  r1_smoothed_series.clear();
  rw_smoothed_series.clear();
  rw_epochs.clear();

  // quantile traces (ss)
  r1_q10.clear();
  rw_q10.clear();

  // quantile traces (os)
  // repeats, for original (smoothed) series
  r1_os_q10.clear();
  rw_os_q10.clear();

}

//
// primary runner for single epoch-level time-series
//

void qdynam_t::proc( const std::vector<double> & x ,
		     const std::vector<int> & e1 ,
		     const std::vector<std::string> & c1 )
{
  
  // ensure we wipe any previous results  
  reinit();
  
  // do we have epochs supplied (ne)
  if ( e1.size() == 0 )
    {
      logger << "  *** warning, no epochs set for qdynam_t()\n";
      return;
    }
  
  // calculate stats both assuming original epoch time and collapsed time
  //  (these will be the same if no include-mask has been set)
  
  std::vector<double> x1 = x;

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
	  
	  if ( c1[ i ] == *cc )
	    {
	      // use previously normed and smoothed values?
	      x2.push_back( norm_each_section ? ox1[i] : r1_smoothed_series[i] );
	      e2.push_back( e1[i] );
	    }
	}
      
      
      // trim?
      if ( trim_epochs[0] != 0 || trim_epochs[1] != 0 )
	{
	  std::vector<double> xx = x2;
	  std::vector<int> ee = e2;
	  x2.clear();
	  e2.clear();

	  int end = xx.size() >= trim_epochs[1] ? xx.size() - trim_epochs[1] : trim_epochs[0] ; 
	  
	  for (int i=trim_epochs[0] ; i < end ; i++)
	    {
	      x2.push_back( xx[i] );
	      e2.push_back( ee[i] );
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

	  rwa.lma1 += w * rw[ *cc ].lma1;
          rwa.lmb2 += w * rw[ *cc ].lmb2;
	  rwa.r_lma1 += w * rw[ *cc ].r_lma1;
          rwa.r_lmb2 += w * rw[ *cc ].r_lmb2;

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
      
      rwa.lma1 /= denom ;
      rwa.lmb2 /= denom ;
      rwa.r_lma1 /= denom ;
      rwa.r_lmb2 /= denom ;

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
  // linear model (w/ non-linear/interaction terms)
  //

  Eigen::VectorXd Y  = Eigen::VectorXd::Zero( nn );
  Eigen::MatrixXd X  = Eigen::MatrixXd::Zero( nn , 2 ); // ee, ee^2
  Eigen::MatrixXd Xr = Eigen::MatrixXd::Zero( nn , 2 ); // ee, ee^2 based on rank
  Eigen::MatrixXd Z  = Eigen::MatrixXd::Zero( nn , 0 ); // no covariates

  const double ee_mean = MiscMath::mean( ee );  
  const double er_mean = MiscMath::mean( e1 );

  for (int i=0; i<nn; i++)
    {
      Y[i] = ss[i];
      X(i,0) = ee[i] - ee_mean;
      X(i,1) = X(i,0) * X(i,0);
      
      Xr(i,0) = i - er_mean;
      Xr(i,1) = Xr(i,0) * Xr(i,0);
    }  
  
  eigen_ops::scale( Y , true , true );
  eigen_ops::scale( X , true , true );
  eigen_ops::scale( Xr , true , true );
    
  const std::vector<std::string> yvars = { "Y" };
  const std::vector<std::string> xvars1 = { "X1" };
  const std::vector<std::string> xvars2 = { "X1","X2" };
  
  // linear term only
  linmod_t lm1( Y, yvars, X.col(0), xvars1, Z );  
  linmod_results_t results1 = lm1.run( 0 ); // i.e. ignore perms
  r.lma1 = results1.beta[ "X1" ][ "Y" ];

  // U term
  linmod_t lm2( Y, yvars, X, xvars2, Z );
  linmod_results_t results2 = lm2.run( 0 ); // i.e. ignore perms  
  r.lmb2 = results2.beta[ "X2" ][ "Y" ];

  // repeat, but w/ rank-based
  // linear term only
  linmod_t r_lm1( Y, yvars, Xr.col(0), xvars1, Z );  
  linmod_results_t r_results1 = r_lm1.run( 0 ); // i.e. ignore perms
  r.r_lma1 = r_results1.beta[ "X1" ][ "Y" ];

  // U term
  linmod_t r_lm2( Y, yvars, X, xvars2, Z );
  linmod_results_t r_results2 = r_lm2.run( 0 ); // i.e. ignore perms  
  r.r_lmb2 = r_results2.beta[ "X2" ][ "Y" ];

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

  writer.value( "LM1" , res.lma1 );
  writer.value( "LM2" , res.lmb2 );

  writer.value( "LM1R" , res.r_lma1 );
  writer.value( "LM2R" , res.r_lmb2 );

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
