
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

#include "dsp/trim.h"

#include "param.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "db/db.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "timeline/cache.h"

#include "stats/Eigen/Dense"
#include "stats/eigen_ops.h"

extern logger_t logger;
extern writer_t writer;

void dsptools::trim_lights( edf_t & edf , param_t & param )
{
  
  // Given a recording, and 1 or more signals, determine where lights off/on occurred
  // based on the assumption of high amplitude noise during the lights on epochs
  
  // by default, trim both start and stop
  
  const bool trim_start = ! param.has( "only-end" ) ;
  const bool trim_stop = ! param.has( "only-start" ) ;

  // by default, use +/- 3 SD units as outlier
  const double th = param.has( "th" ) ? param.requires_dbl( "th" ) : 3 ; 

  // by default, do not allow more than 20 epochs of (10 mins) 'good' data at either end
  const int good_th = param.has( "allow" ) ? param.requires_int( "allow" ) : 20;  

  // by default, require equivalent of 10 epochs of bad data to be flagged - i.e. or else maxima may
  // point to a single trivial outliers
  const int req_epoch = param.has( "req" ) ? param.requires_int( "req" ) : 10 ;
  
  //  by default, smoothing window (total window, in epochs) - otherwise, 4 epochs either side
  const int smooth_win = param.has( "w" ) ? param.requires_int( "w" ) : 9 ; 
  const double smooth_taper = 0.5;
  
  // anchor on sleep stages (to get median and SD), if present, unless this is set ('all') 
  bool anchor_on_sleep = ! ( param.has( "all" ) || param.has( "wake" ) );
  // or 'wake'
  bool anchor_on_wake = param.has( "wake" );
  if ( param.has( "wake" ) && param.has( "all" ) ) 
    Helper::halt( "cannot specify both 'wake' and 'all' options" );

  // use H2 also?
  const bool use_h2 = param.has("h2") ;
  
  // outputs
  const bool verbose = param.has( "epoch" ) || param.has( "verbose" );
  
  // set mask?
  const bool set_mask = param.has( "mask" ) ;

  // which signals?
  const bool NO_ANNOTS = true;

  signal_list_t signals = edf.header.signal_list(  param.requires( "sig" ) , NO_ANNOTS );
  
  const int ns = signals.size();

  if ( ns == 0 ) return;
  
  std::vector<double> Fs = edf.header.sampling_freq( signals );

  // epoch-wise storage

  int ne = edf.timeline.first_epoch();

  Eigen::MatrixXd H1 = Eigen::MatrixXd::Zero( ne , ns );
  Eigen::MatrixXd H2 = Eigen::MatrixXd::Zero( use_h2 ? ne : 0  , use_h2 ? ns : 0 );
  Eigen::MatrixXd H3 = Eigen::MatrixXd::Zero( ne , ns );
  std::vector<int> ep;

  // get stages, if present

  std::vector<bool> use( ne , true );
    
  if ( anchor_on_sleep || anchor_on_wake )
    {
      // get staging
      edf.annotations->make_sleep_stage( edf.timeline );      
      const bool has_staging = edf.timeline.hypnogram.construct( &(edf.timeline) , param , false );
      int cnt = 0;
      if ( has_staging ) 
	{
	  if ( ne != edf.timeline.hypnogram.stages.size() )
	    Helper::halt( "internal error extracting staging" );
	  
	  for (int ss=0; ss<ne; ss++)
	    {
	      bool is_sleep = edf.timeline.hypnogram.stages[ ss ] == NREM1 
		|| edf.timeline.hypnogram.stages[ ss ] == NREM2
		|| edf.timeline.hypnogram.stages[ ss ] == NREM3
		|| edf.timeline.hypnogram.stages[ ss ] == REM ;
	      
	      bool is_wake = edf.timeline.hypnogram.stages[ ss ] == WAKE;

	      if ( ( anchor_on_sleep && ! is_sleep ) || ( anchor_on_wake && ! is_wake ) ) 
		{
		  use[ss] = false;
		  ++cnt;
		}
	    }
	  logger << "  anchoring on " << ( anchor_on_sleep ? "sleep" : "wake" ) 
		 << " epochs only for normative ranges, using " << ne - cnt << " of " << ne << " epochs\n";

	  // requires at least 10 epochs) 
	  if ( ne - cnt < 10 ) 
	    {
	      logger << "  could not find 10+ valid epochs, so not anchoring on sleep/wake epochs only\n";
	      anchor_on_sleep = anchor_on_sleep = false;
	    }
	}
      else
	{
	  logger << "  could not find any valid stages, so not anchoring on sleep epochs only\n";
	  anchor_on_sleep = false; 
	}
    }
  
  //
  // iterate over epochs
  //


  // need to reset, as may have stepped through epochs in HYPNO
  ne = edf.timeline.first_epoch();

  int ecnt = 0;
  
  while ( 1 )
    {
      
      int epoch = edf.timeline.next_epoch();	  
      //      std::cout << "ep " << epoch << "\n";
      if ( epoch == -1 ) break;
      
      interval_t interval = edf.timeline.epoch( epoch );
      
      // get all signals
      eigen_matslice_t mslice( edf , signals , interval );

      Eigen::MatrixXd & X = mslice.nonconst_data_ref();

      // process each signal
      for (int s=0; s<ns; s++)
	{
	  
	  std::vector<double> v = eigen_ops::copy_vector( X.col(s) );
	  
	  double activity = 0 , mobility = 0 , complexity = 0;

	  MiscMath::centre( &v );

	  MiscMath::hjorth( &v , &activity , &mobility , &complexity , ! globals::legacy_hjorth );
	  
	  activity = activity > 0 ? log( activity ) : log( activity + 1e-12 );
	  
	  H1( ecnt , s ) = activity ;
	  if ( use_h2 ) H2( ecnt , s ) = mobility ;
	  H3( ecnt , s ) = complexity ;
	}
      
      // track epochs
      ep.push_back( epoch );
      ++ecnt;

    }
  
  //
  // we now have all epoch level data H1, H2, H3
  //
  
  
  // track overall estimates of lights off/on
  
  int lights_off = -1, lights_on = -1;
  
  // proceed channel-wise
  
  for (int s=0; s<ns; s++)
    {
      // std::cout << " S = " << s << "\n";
      
      // std::cout << "H1 " << H1.rows() << " " << H1.cols() << "\n";
      // std::cout << "H2 " << H2.rows() << " " << H2.cols() << "\n";
      // std::cout << "H3 " << H3.rows() << " " << H3.cols() << "\n";      
      
      // get IQRs
      std::vector<double> v1 = eigen_ops::copy_vector( H1.col(s) );
      std::vector<double> v2;
      if ( use_h2 ) v2 = eigen_ops::copy_vector( H2.col(s) );
      std::vector<double> v3 = eigen_ops::copy_vector( H3.col(s) );
      
      // reduce subset for stats?
      std::vector<double> r1, r2, r3;
      if ( anchor_on_sleep )
	{
	  r1 = Helper::subset( v1 , use ) ;
	  if ( use_h2 ) r2 = Helper::subset( v2 , use ) ;
	  r3 = Helper::subset( v3 , use ) ;
	}

      // get medians
      double m1 = MiscMath::median( anchor_on_sleep ? r1 : v1 );
      double m2 = use_h2 ? MiscMath::median( anchor_on_sleep ? r2 : v2 ) : 0 ;
      double m3 = MiscMath::median( anchor_on_sleep ? r3 : v3 );
      
      //            std::cout << "yy\n";
		  
      // factor of 0.7413 to get robust estimate of SD from IQR
      // double sd1 = 0.7413 * MiscMath::iqr( v1 ) ;
      // double sd2 = 0.7413 * MiscMath::iqr( v2 ) ;
      // double sd3 = 0.7413 * MiscMath::iqr( v3 ) ;

      // factor of 0.7413 to get robust estimate of SD from IQR
      
      double sd1 = MiscMath::sdev( anchor_on_sleep ? r1 : v1 ) ;
      double sd2 = use_h2 ? MiscMath::sdev( anchor_on_sleep ? r2 : v2 )  : 0 ;
      double sd3 = MiscMath::sdev( anchor_on_sleep ? r3 : v3 ) ;

      //      std::cout << "zz\n";
	    
      double lwr1 = m1 - th * sd1;
      double upr1 = m1 + th * sd1;
      
      double lwr2 = m2 - th * sd2;
      double upr2 = m2 + th * sd2;
      
      double lwr3 = m3 - th * sd3;
      double upr3 = m3 + th * sd3;
      
      // outliers
      std::vector<bool> okay( ne , true );
      Eigen::VectorXd out = Eigen::VectorXd::Zero( ne );
      int f1 = 0, f2 = 0, f3 = 0, fn = 0;
      for (int e=0; e<ne; e++)
	{
	  //	  std::cout << " out e " << e << "\n";
	  const bool o1 = v1[e] < lwr1 || v1[e] > upr1 ;
	  const bool o2 = use_h2 ? ( v2[e] < lwr2 || v2[e] > upr2 ) : false ;
	  const bool o3 = v3[e] < lwr3 || v3[e] > upr3 ;

	  if ( o1 ) ++f1;
	  if ( o2 ) ++f2;
	  if ( o3 ) ++f3;
	  
	  //if ( o1 || o2 || o3 )
	  if ( o1 || o3 ) 
	    {
	      ++fn;
	      out[e] = 1;
	      okay[e] = false;
	    }
	}

      logger << "  for " << signals.label(s) << ", flagged " << fn << " epochs (";
      logger << "H1=" << f1 ;
      if ( use_h2 ) logger << ", H2="<<f2;
      logger << ", H3=" << f3 << ")\n";
      
      logger << "  H(1) bounds: " << lwr1 << " .. " << upr1 << "\n";
      if ( use_h2 ) logger << "  H(2) bounds: " << lwr2 << " .. " << upr2 << "\n";
      logger << "  H(3) bounds: " << lwr3 << " .. " << upr3 << "\n";

      // smooth (turned off if w=0)
      if ( smooth_win ) 
	out = eigen_ops::tri_moving_average( out , smooth_win , smooth_taper );
      
      // estimate for each possible 'e'   sum(X)^2 / n
      //   where n is # of epochs (from start, or unitl end)
      //   i.e. this equals sum(X) * mean(X)
      //   i.e. weights both the total amount of X and the density... ~equal balance

      // but also impose rule that if we've had more than T epochs of 'good' data, then
      // we stop tracking
      
      // lights off
      int lights_off1 = -1;
      double max_off = 0;
      double cum = 0;
      int good = 0;
      
      std::vector<double> trk_off( ne);
      std::vector<double> trk_on( ne);
      
      for (int e=0; e<ne; e++)
	{

	  // do not allow more than (e.g.) 20 epochs of 'good' (non-outlier) epochs in this region)
	  if ( okay[e] ) ++good;
	  if ( good > good_th ) break;
	  
	  cum += out[e] ;
	  double stat = pow( cum, 2) / (double)(e+1);
	  if ( stat > max_off )
	    {
	      max_off = stat;
	      lights_off1 = e;
	    }
	  trk_off[e] = stat;
	}
      
           
      // lights on
      int lights_on1 = -1;
      double max_on = 0;
      cum = 0;
      good = 0;
      for (int e1=0; e1<ne; e1++)
        {
	  int e = ne - 1 - e1;
	  
	  if ( okay[e] ) ++good;
          if ( good > good_th )	break;

	  cum += out[e] ;
	  double stat = pow( cum, 2 )  / (double)(e1+1); // nb e1 denom, i.e. # until end
	  if ( stat > max_on )
            {
              max_on = stat;
              lights_on1 = e;
	    }
	  trk_on[e] = stat;
        }

      //
      // at this point, lights_off1 & lights_on1 are 0-based epoch numbers,
      // that are the last epoch *before* we want to set lights off, i.e.
      // we should trim inclusive of these points
      //
      
      //
      // magnitude check: stats ~equals number of epochs (weighted) 
      //  requires this to be at least e.g. 20 (10 mins) or else what's the point of changing
      //
      
      if ( max_off < req_epoch )
	lights_off1 = -1;
	
      if ( max_on < req_epoch )
	lights_on1 = -1;

      
      //
      // channel-wise output
      //

      writer.level( signals.label(s) , globals::signal_strat );

      if ( trim_start && lights_off1 >= 0 ) 
	writer.value( "EOFF" , lights_off1 );
      if ( trim_stop && lights_on1 >= 0 ) 
	writer.value( "EON" , lights_on1 );
      
      if ( verbose )
	{	  
	  for (int e=0; e<ne; e++)
	    {
	      writer.epoch( edf.timeline.display_epoch( ep[e] ) );
	      if ( trim_start ) writer.value( "XOFF" , trk_off[e] );
	      if ( trim_stop ) writer.value( "XON" , trk_on[e] );

	      // nb. inclusive <= comparison here
	      if ( lights_off1 >= 0 && e <= lights_off1 ) 
		writer.value( "TRIM" , 1 );
	      else if ( lights_on1 >= 0 && e >= lights_on1 ) 
		writer.value( "TRIM" , 1 );
	      else
		writer.value( "TRIM" , 0 );

	      writer.value( "STAT" , out[e] );
	      writer.value( "FLAG" , okay[e] ? 0 : 1 ) ;
	      writer.value( "H1" , H1(e,s) );
	      if ( use_h2 ) 
		writer.value( "H2" , H2(e,s) );
	      writer.value( "H3" , H3(e,s) );
	      
	    }
	  writer.unepoch();
	  
	}

      //
      // take earliest lights on and latest lights off
      //

      // e.g. lights_off1 is the single-channel empirical value (based on EDGER stats)
      //      lights_off is the 'max' cut based on all channels
      //      --> do we want to update that?
      
      if ( lights_off1 >= 0 )
	{
	  if ( lights_off < 0 ) lights_off = lights_off1;
	  else if ( lights_off1 < lights_off ) lights_off = lights_off1;
	}
      
      if ( lights_on1 >= 0 )
	{
	  if ( lights_on < 0 )	lights_on = lights_on1;
	  else if ( lights_on1 > lights_on ) lights_on = lights_on1;
	}
      
      // next signal
    }
  
  writer.unlevel( globals::signal_strat );
  

  //
  // final determination
  //

  const bool set_off = trim_start && lights_off >= 0 ;
  const bool set_on  = trim_stop && lights_on >= 0 ;
  
  clocktime_t starttime( edf.header.starttime );
  double epoch_seconds = edf.timeline.epoch_length() ; 
  clocktime_t clock_lights_out = starttime;
  clocktime_t clock_lights_on = starttime;
  
  if ( set_off ) 
    {            
      // adjust lights_off : 
      // +1 --> lights_off is inclusive off the
      // +1 --> for output, use 1-based output      
      
      lights_off = lights_off < ne ? lights_off + 1 : lights_off ; 
      writer.value( "EOFF" , lights_off + 1 ); // for 1-based output
      clock_lights_out.advance_seconds( epoch_seconds * lights_off ) ;
      writer.value( "LOFF" , clock_lights_out.as_string() );

    }

  if ( set_on )
    {
      // adjust lights_on :
      // -1 --> lights_off is inclusive of this
      // +1 --> for output, use 1-based output

      lights_on = lights_on > 0 ? lights_on - 1 : lights_on ;

      writer.value( "EON" , lights_on + 1 ); // for output
      clock_lights_on.advance_seconds( epoch_seconds * lights_on ) ;
      writer.value( "LON" , clock_lights_on.as_string() );
    }
  

  //
  // set MASK
  //

  if ( set_mask )
    {

      // std::cout << " set_off " << lights_off << "\n" 
      //  		<< " set_on  " << lights_on << "\n" 
      // 		<< " ne      " << ne << "\n";      

      // 
      int cnt = 0 ;

      // lights_off == 0 means lights out starts at start of first epoch
      //  i.e. so number to be masked here is only if lights_off is 1 or more 

      if ( set_off ) cnt = lights_off > 0 ? lights_off - 1 : 0 ; 
      
      // likewise, if lights_on is ne-1, then there is no masked
      //  otherwise, is ne - 1 ) - lights_on 
      //  0 1 2 3 4 5 6 7 8 9
      //                       lights_off = 9
      //                    X               8
      //                  X X               7

      if ( set_on ) cnt += lights_on < ne - 1 ? ne - 1 - lights_on : 0 ; 
      
      if ( cnt )
	{
	  logger << "\n  masking " << cnt << " epochs\n";
	  
	  // set to EXCLUDE
	  bool include_mode = false;
	  
	  // n.b. expects 1-based terms in 'epochs'
	  std::set<int> epochs;

	  if ( set_off ) 
	    for (int e=0; e<lights_off-1; e++)
	      epochs.insert( e+1 );
	  
	  if ( set_on ) 
	    for (int e=lights_on+1; e<ne; e++)
	      epochs.insert( e+1 );
	  
	  // expects 1-based terms in 'epochs'
	  edf.timeline.select_epoch_range( epochs , include_mode );	      
	
	}
    }

  
  //
  // use cache to remember LON and LOFF values? [ will enable hypno to understand these ]
  //

  if ( set_off || set_on )
    {

      cache_t<double> * cache = param.has( "cache" ) ? edf.timeline.cache.find_num( param.value( "cache" ) ) : NULL ;

      if ( cache )
	logger << "  setting cache " << param.value( "cache" ) << " to store times\n";
      
      clocktime_t starttime( edf.header.starttime );
      double epoch_seconds = edf.timeline.epoch_length();
      
      clocktime_t clock_lights_out = starttime;
      if ( set_off ) 
	clock_lights_out.advance_seconds( epoch_seconds * lights_off ) ;
      
      clocktime_t clock_lights_on = starttime;
      if ( set_on ) 
	clock_lights_on.advance_seconds( epoch_seconds * lights_on ) ;
      
      if ( set_off )
	{
	  logger << "  lights-off=" << clock_lights_out.as_string() << " (skipping " << lights_off - 1 << " epochs from start)\n";	  
	  if ( cache ) cache->add( ckey_t( "LOFF" , writer.faclvl() ) , edf.timeline.epoch_length() * lights_off );
	}
      
      if ( set_on )
	{
	  logger << "  lights-on=" << clock_lights_on.as_string() << " (skipping " << ne - 1 - lights_on << " epochs from end)\n";
	  if ( cache ) cache->add( ckey_t( "LON" , writer.faclvl() ) , edf.timeline.epoch_length() * lights_on );
	}
      
    }
  else
    logger << "  no trimming indicated: did not alter lights-off or lights-on times\n";
  
  
}


