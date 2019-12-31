

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

#include "slow-waves.h"

#include "dsp/fir.h"
#include "dsp/hilbert.h"
#include "miscmath/miscmath.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "cwt/cwt.h"
#include "db/db.h"
#include "miscmath/crandom.h"
#include "eval.h"

#include "helper/helper.h"
#include "helper/logger.h"

extern writer_t writer;

extern logger_t logger;

double slow_waves_t::nearest( const int i , int * sw_idx  ) const
{

  *sw_idx = -1;

  // in a slow wave?
  if ( in_sw[i] != -1 ) { *sw_idx =  in_sw[i]; return 0; } 
  
  // search backward, then forward
  int tbck = i , tfor = i;
  
  while ( 1 ) 
    {
      --tbck;
      if ( tbck < 0 ) { break; } 
      if ( in_sw[tbck] != -1 ) break;
    }
  
  while ( 1 ) 
    {
      ++tfor;
      if ( tfor >= in_sw.size() ) { tfor = -1; break; } 
      if ( in_sw[tfor] != -1 ) break;
    }
  
  uint64_t tp_bck = tbck > 0 ?    tp[i] - tp[tbck]  : 0;
  uint64_t tp_for = tfor > 0 ?    tp[tfor] - tp[i] : 0;
  
  double sec_bck = - ( (double)tp_bck / (double)globals::tp_1sec );
  double sec_for = (double)tp_for / (double)globals::tp_1sec;

  

  if ( tbck > 0 && tfor > 0 ) 
    {
      if ( fabs( sec_bck ) < fabs( sec_for ) ) 
	{
	  *sw_idx = in_sw[tbck];
	  return sec_bck;
	}
      else
	{
	  *sw_idx = in_sw[tfor];
	  return sec_for;
	}
    }
  
  else if ( tbck > 0 ) { *sw_idx = in_sw[tbck]; return sec_bck; } 
  else if ( tfor > 0 ) { *sw_idx = in_sw[tfor]; return sec_for; } 
  else { *sw_idx = 0; return 0; } // whole sample duration, if not found
  
}


slow_waves_t::slow_waves_t( edf_t & edf , const param_t & param )
{

  std::string signal_label = param.requires( "sig" );   

  signal_list_t signals = edf.header.signal_list( signal_label );  
  
 
  //
  // default parameters from Latchoumane et al (2017) (nb. mouse EEG)
  //

  const double thr   = param.has( "mag" ) ? param.requires_dbl( "mag" )       : 0 ;

  const double f_lwr = param.has( "f-lwr" ) ? param.requires_dbl( "f-lwr" ) : 0.2  ; 
  const double f_upr = param.has( "f-upr" ) ? param.requires_dbl( "f-upr" ) : 4.5  ; 

  const double t_lwr = param.has( "t-lwr" ) ? param.requires_dbl( "t-lwr" ) : 0 ;
  const double t_upr = param.has( "t-upr" ) ? param.requires_dbl( "t-upr" ) : 10 ;
  
  const double uV_neg = param.has("uV-neg" ) ? param.requires_dbl( "uV-neg" ) : 0 ; 
  if ( uV_neg > 0 ) Helper::halt( "uV-neg should be negative" ) ;

  const double uV_p2p = param.has("uV-p2p" ) ? param.requires_dbl( "uV-p2p" ) : 0 ; 
  if ( uV_p2p < 0 ) Helper::halt( "uV-p2p should be positive" ) ;

  const double t_neg_lwr = param.has( "t-neg-lwr" ) ? param.requires_dbl( "t-neg-lwr" ) : 0;
  const double t_neg_upr = param.has( "t-neg-upr" ) ? param.requires_dbl( "t-neg-upr" ) : 10;
  
  // for full waves, default is to find positive to negative zero-crossings
  const bool   use_alternate_neg2pos_zc = param.has( "neg2pos" ) ;
  
  // use mean instead of median for thresholds
  const bool use_mean = param.has( "th-mean" ) ;

  // use mean instead of median for report
  report_median_stats = param.has( "stats-median" ) ;
  
  logger << " stats based on " 
	 << ( report_median_stats ? "median" : "mean" ) 
	 << " over SOs\n";

  // type of SO to detect: full, half (and negative or positive only)

  slow_wave_type type = SO_FULL;
  if      ( param.has( "half-wave" ) ) type = SO_HALF;
  else if ( param.has( "negative-half-wave" ) ) type = SO_NEGATIVE_HALF;
  else if ( param.has( "positive-half-wave" ) ) type = SO_POSITIVE_HALF;


  //
  // iterate over signals
  //

  const int ns = signals.size();
  
  interval_t interval = edf.timeline.wholetrace();

  for (int s=0;s<ns;s++)
    {

      //
      // Only consider raw signal channels
      //
      
      if ( edf.header.is_annotation_channel( signals(s) ) )
	continue;
      
      logger << " estimating SO for " << signals.label(s) << "\n";
      
      writer.level( signals.label(s) , globals::signal_strat );
      
      double sr = edf.header.sampling_freq( signals )[s];
      
      slice_t slice( edf , signals(s) , interval );
      
      const std::vector<double> * d = slice.pdata();
      
      const std::vector<uint64_t> * tp = slice.ptimepoints(); 

      
      //
      // Detect slow waves
      //
      
      detect_slow_waves( *d , *tp , sr , thr, use_mean , 
			 uV_neg , uV_p2p , 
			 f_lwr , f_upr , 
			 t_lwr , t_upr , 
			 t_neg_lwr, t_neg_upr , 
			 use_alternate_neg2pos_zc , type ); 
      
      
      //
      // spectral analysis around SOs
      //
      
      phase_slow_waves();
      
      
      //
      // verbose display
      //
      
      display_slow_waves( param.has( "verbose" ) , &edf );
      
      
      //
      // Optionally, consider another signal w.r.t SO 
      //
      
      
      if ( param.has( "tl" ) ) 
	{
	  
	  int position = -1;
	  
	  if ( param.has("onset" ) ) position = 0;
	  else if ( param.has("pos" ) ) position = +1;
	  
	  double twin = param.has("window") ? param.requires_dbl( "window" ) : 3.0 ; 
	  
	  std::string label2 = param.requires( "tl" );   
	  
	  signal_list_t signals2 = edf.header.signal_list( label2 );
	  
	  const int ns2 = signals2.size();
	  
	  logger << " averaging " << label2 << " based on time-locked averaging to SO ";
	  if ( position == -1 ) logger << "negative peak";
	  else if ( position == 0 ) logger << "onset";
	  else if ( position == 1 ) logger << "positive peak";
	  logger << ", within window of +/-" << twin << " seconds\n";
	  
	  
	  for (int i=0;i<ns2;i++)
	    {
	      
	      double sr2 = edf.header.sampling_freq( signals2 )[i];
	      
	      interval_t interval2 = edf.timeline.wholetrace();
	      
	      slice_t slice2( edf , signals2(i) , interval2 );
	      
	      const std::vector<double> * d2 = slice2.pdata();
	      
	      const std::vector<uint64_t> * tp2 = slice2.ptimepoints(); 
	      
	      std::vector<double> tl_sig = time_locked_averaging( d2 , sr2 , twin , twin , position );
	      
	      if ( tl_sig.size() > 0 ) 
		{
		  writer.var( "SOTL_SIG" , "Slow wave time-locked averages" );
		  writer.level( signals2.label(i)  , "CH2" );
		  
		  int sz = tl_sig.size();
		  int sz2 = - (sz-1)/2;
		  
		  for (int j=0;j<sz;j++)
		    {
		      writer.level( sz2 , "SP" );
		      writer.value( "SOTL" , tl_sig[j] );		  
		      ++sz2;
		    }
		  writer.unlevel( "SP" );
		}
	      
	      // next channel
	    }
	  writer.unlevel( "CH2" );
	  
	}
      
      
      //
      // Next signal
      //
      
    }

  writer.unlevel( globals::signal_strat );
  
}


void slow_waves_t::display_slow_waves( bool verbose , edf_t * edf )
{

  //
  // Variables
  //

  writer.var( "SO" , "Number of slow waves" );
  writer.var( "SO_RATE" , "Slow waves per minute" );
  
  writer.var( "SO_DUR" , "Median slow wave duration (sec)" );
  writer.var( "SO_NEG_DUR" , "Median negative half-wave duration (sec)" );
  writer.var( "SO_POS_DUR" , "Median positive half-wave duration (sec)" );
  writer.var( "SO_AMP" , "Median slow wave minimim peak amplitude" );
  writer.var( "SO_P2P" , "Median slow wave peak-to-peak amplitude" );
  writer.var( "SO_SLOPE_NEG1" , "Median slow wave down-going slope on negative peak (uV/sec)" );
  writer.var( "SO_SLOPE_NEG2" , "Median slow wave up-going slope on negative peak (uV/sec)" );
  writer.var( "SO_SLOPE_POS1" , "Median slow wave up-going slope on positive peak (uV/sec)" );
  writer.var( "SO_SLOPE_POS2" , "Median slow wave down-going slope on positive peak (uV/sec)" );

  //
  // Output 
  //

  writer.value( "SO" , num_waves() );
  writer.value( "SO_RATE" , num_waves() / ( signal_duration_sec / 60.0 ) );

  writer.value( "SO_TH_NEG" , th_x );
  writer.value( "SO_TH_P2P" , th_yminusx );

  if ( num_waves() > 0 ) 
    {

      //
      // statistics are means over SOs
      //
    
      if ( ! report_median_stats )
	{
	  writer.value( "SO_DUR" , avg_duration_sec );    
	  writer.value( "SO_NEG_DUR" , avg_negative_duration_sec );    
	  writer.value( "SO_POS_DUR" , avg_positive_duration_sec );    
	  writer.value( "SO_AMP" , avg_x );
	  writer.value( "SO_P2P" , avg_yminusx );
	  
	  // may not be calculated, i.e.  if looking at half-waves  
	  if ( avg_slope_n1 != 0 ) writer.value( "SO_SLOPE_NEG1" , avg_slope_n1 );  
	  if ( avg_slope_n2 != 0 ) writer.value( "SO_SLOPE_NEG2" , avg_slope_n2 );
	  if ( avg_slope_p1 != 0 ) writer.value( "SO_SLOPE_POS1" , avg_slope_p1 );  
	  if ( avg_slope_p2 != 0 ) writer.value( "SO_SLOPE_POS2" , avg_slope_p2 );
	}

      
      //
      // or medians
      //
      
      if ( report_median_stats )
	{
	  writer.value( "SO_DUR" , median_duration_sec );    
	  writer.value( "SO_NEG_DUR" , median_negative_duration_sec );    
	  writer.value( "SO_POS_DUR" , median_positive_duration_sec );    
	  writer.value( "SO_AMP" , median_x );
	  writer.value( "SO_P2P" , median_yminusx );
      
	  // may not be calculated, i.e.  if looking at half-waves  
	  if ( median_slope_n1 != 0 ) writer.value( "SO_SLOPE_NEG1" , median_slope_n1 );  
	  if ( median_slope_n2 != 0 ) writer.value( "SO_SLOPE_NEG2" , median_slope_n2 );
	  if ( median_slope_p1 != 0 ) writer.value( "SO_SLOPE_POS1" , median_slope_p1 );  
	  if ( median_slope_p2 != 0 ) writer.value( "SO_SLOPE_POS2" , median_slope_p2 );
	}

    }


  //
  // Per SO output
  //

  for (int i=0;i<sw.size();i++)
      {

	const slow_wave_t & w = sw[i];
	
	writer.level( i+1 , globals::count_strat );
	
	writer.value( "START_IDX" , static_cast<int>(w.interval.start) );
	writer.value( "STOP_IDX"  , static_cast<int>(w.interval.stop ) );
	
	writer.value( "START" , w.interval_tp.start * globals::tp_duration );
	writer.value( "STOP"  , w.interval_tp.stop * globals::tp_duration );
	
	writer.value( "DUR"  , w.interval_tp.duration_sec() );
	
	writer.value( "UP_AMP" , w.up_amplitude );
	writer.value( "DOWN_AMP" , w.down_amplitude );
	writer.value( "P2P_AMP" , w.amplitude() );

	writer.value( "UP_IDX" , w.up_peak_sp );
	writer.value( "DOWN_IDX" , w.down_peak_sp );
	
	if ( w.type == SO_FULL || w.type == SO_NEGATIVE_HALF )
	  {
	    writer.value( "SLOPE_NEG1" , w.slope_n1() );
	    writer.value( "SLOPE_NEG2" , w.slope_n2() );
	  }
	
	if ( w.type == SO_FULL || w.type == SO_POSITIVE_HALF )
	  {
	    writer.value( "SLOPE_POS1" , w.slope_p1() );
	    writer.value( "SLOPE_POS2" , w.slope_p2() );
	  }
	
	if ( 0 && verbose )
	  {
	    int pos = 0;
	    for (uint64_t j = w.interval.start ; j <= w.interval.stop ; j++ ) 
	      {
		writer.level( pos , globals::sample_strat );
		writer.value( "FLT" , filtered[ j ] );
		writer.value( "PH" , phase[ j ] ); 
		++pos;
	      }
	    writer.unlevel( globals::sample_strat );
	  }
      }
  writer.unlevel( globals::count_strat );
  
  
  //
  // Epoch-level counts of SW, and means of other key SW statistics
  //
  
  if ( edf != NULL )
    {
      edf->timeline.first_epoch();
      
      std::vector<int> epoch_counts;
      
      while ( 1 ) 
	{
	  
	  int epoch = edf->timeline.next_epoch();      
	  
	  if ( epoch == -1 ) break;
	  
	  interval_t interval = edf->timeline.epoch( epoch );
	  
	  const int nso = sw.size();
	  
	  int so_epoch = 0;
	  
	  std::set<int> sw_in_epoch;

	  for (int i=0 ; i<nso; i++)
	    {
	      
	      const slow_wave_t & w = sw[i];
	      
	      // dummy interval for just starting point of SO
	      // i.e. for purpose of assigning to EPOCH (so as not
	      // to double-count SO that overlap epoch boundaries
	      
	      // this search strategy is v. inefficient, but in this context
	      // doesn't really seem to matter (at least w/ large epoch sizes)
	      // but we should revisit this...
	      
	      interval_t sostart( w.interval_tp.start , w.interval_tp.start );
	      
	      if ( interval.overlaps( sostart ) )
		{
		  ++so_epoch;
		  sw_in_epoch.insert( i );
		}
	      else if ( sostart.is_after( interval ) )
		break; // SO are in order, so can skip
	    }
	  
	  // record
	  epoch_counts.push_back( so_epoch );
	  
	  //
	  // per-epoch output
	  //
	  
	  writer.epoch( edf->timeline.display_epoch( epoch ) );

	  //
	  // per-epoch SO count
	  //

	  writer.value( "N" , so_epoch );

	  //
	  // and mean statistics
	  //

	  double mean_dur = 0 , mean_up_amp = 0 , mean_down_amp = 0 , mean_p2p_amp = 0;
	  
	  double mean_slope_n1 = 0 , mean_slope_n2 = 0 , mean_slope_p1 = 0  , mean_slope_p2 = 0;
	  int n_pos = 0 , n_neg = 0;

	  std::set<int>::const_iterator jj = sw_in_epoch.begin();
	  while ( jj != sw_in_epoch.end() )
	    {
	      
	      const slow_wave_t & w = sw[ *jj ];

	      mean_dur += w.interval_tp.duration_sec() ;
	      mean_up_amp += w.up_amplitude ;
	      mean_down_amp += w.down_amplitude ;
	      mean_p2p_amp += w.amplitude() ;	      
	      
	      if ( w.type == SO_FULL || w.type == SO_NEGATIVE_HALF )
		{
		  mean_slope_n1 += w.slope_n1();
		  mean_slope_n2 += w.slope_n2();
		  ++n_neg;
		}
	      
	      if ( w.type == SO_FULL || w.type == SO_POSITIVE_HALF )
		{
		  mean_slope_p1 += w.slope_p1() ;
		  mean_slope_p2 += w.slope_p2() ;
		  ++n_pos;
		}
	      ++jj;
	    }
	  

	  // output epoch-level SW morphology stats
	  
	  if ( sw_in_epoch.size() > 0 ) 
	    {
	      writer.value( "DUR"  , mean_dur / (double)sw_in_epoch.size() );
	
	      writer.value( "UP_AMP" , mean_up_amp / (double)sw_in_epoch.size() );
	      writer.value( "DOWN_AMP" , mean_down_amp / (double)sw_in_epoch.size() );
	      writer.value( "P2P_AMP" , mean_p2p_amp / (double)sw_in_epoch.size() );
	      
	      if ( n_neg > 0 ) 
		{
		  writer.value( "SLOPE_NEG1" , mean_slope_n1 / (double)n_neg );
		  writer.value( "SLOPE_NEG2" , mean_slope_n2 / (double)n_neg );
		}
	      
	      if ( n_pos > 0 ) 
		{
		  writer.value( "SLOPE_POS1" , mean_slope_p1 / (double)n_pos );
		  writer.value( "SLOPE_POS2" , mean_slope_p2 / (double)n_pos );
		}
	    }

	}
      

      // end of epoch-level reporting 

      writer.unepoch();      
    }
  
  // all done 
  
}

slow_waves_t::slow_waves_t( const std::vector<double> & unfiltered , 
			    const std::vector<uint64_t> & tp ,
			    const int sr , 
			    const double thr ,
			    const bool   use_mean , 
			    const double uV_neg , 
			    const double uV_p2p , 
			    const double f_lwr , 
			    const double f_upr ,
			    const double t_lwr ,
			    const double t_upr , 
			    const double t_neg_lwr , 
			    const double t_neg_upr , 
			    const bool neg2pos , 
			    const slow_wave_type type )
{
  
  report_median_stats = false;
  
  detect_slow_waves( unfiltered, tp , sr , thr, use_mean , uV_neg , uV_p2p , f_lwr, f_upr, 
		     t_lwr, t_upr , t_neg_lwr , t_neg_upr , neg2pos , type );
}


int slow_waves_t::detect_slow_waves( const std::vector<double> & unfiltered , 
				     const std::vector<uint64_t> & _tp ,
				     const int sr , 
				     const double thr ,
				     const bool   use_mean , 
				     const double uV_neg , 
				     const double uV_p2p , 
				     const double f_lwr , 
				     const double f_upr ,
				     const double t_lwr ,
				     const double t_upr , 
				     const double t_neg_lwr , 
				     const double t_neg_upr , 
				     const bool use_alternate_neg2pos_zc , 
				     const slow_wave_type type ,
				     const double fir_ripple ,
				     const double fir_tw )
{
  
  // track Fs for later
  Fs = sr;
  tp = _tp;

  // track total signal duration
  signal_duration_sec = unfiltered.size() / (double)sr ; 

  logger << " detecting slow waves: " << f_lwr << "-" << f_upr << "Hz\n";
  
  if ( t_lwr > 0 ) logger << "  - duration " << t_lwr << "-" << t_upr << "s\n"; 
  if ( t_neg_lwr > 0 ) logger << "  - negative half-wave duration " << t_neg_lwr << "-" << t_neg_upr << "\n";

  if ( thr > 0 ) logger << "  - relative threshold " << thr  << "x " <<  ( use_mean ? "mean" : "median" ) << "\n";

  if ( uV_neg < 0 ) logger << "  - absolute threshold based on " 
			      << uV_neg << " uV for negative peak, " 
			      << uV_p2p << " peak-to-peak\n";


  if ( type == SO_FULL ) 
    logger << "  - full waves, based on consecutive "  
	      << ( use_alternate_neg2pos_zc ? "negative-to-positive" : "positive-to-negative" ) << " zero-crossings\n";
  else if ( type == SO_HALF ) 
    logger << "  - all half waves\n";
  else if ( type == SO_NEGATIVE_HALF ) 
    logger << "  - all negative half waves\n";
  else if ( type == SO_POSITIVE_HALF ) 
    logger << "  - all positive half waves\n";


  
    
  //
  // Band-pass filter for slow waves
  //

  
  filtered = dsptools::apply_fir( unfiltered , sr , fir_t::BAND_PASS ,fir_ripple , fir_tw , f_lwr , f_upr );
  
  //filtered = band_pass_filter( unfiltered , sr , filter_order , f_lwr , f_upr );  
  //  std::cout << MiscMath::mean( d ) << " is the mean \n";
  
  const int n = filtered.size();
  
  // get zero crossings
  std::vector<int> zc;
  
  if ( use_alternate_neg2pos_zc ) 
    {
      for (int i=1;i<n;i++) if ( filtered[i] >= 0 && filtered[i-1] < 0 ) zc.push_back(i);
    }
  else
    {
      for (int i=1;i<n;i++) if ( filtered[i] < 0 && filtered[i-1] >= 0 ) zc.push_back(i);
    }
  
  // averages/medians
  std::vector<double> tmp_x, tmp_yminusx;

  // # of putative SOs
  int cnt = 0;  

  logger << " " << zc.size() << " zero crossings";
  
  // putative waves
  std::vector<slow_wave_t> waves; 

  // get intervals that are within the time threshold

  for (int i=1;i<zc.size();i++)
    {
      
      // check that interval is not discontinuous
      
      if ( timeline_t::discontinuity( tp , sr , zc[i-1] , zc[i] ) ) continue;
 
      // get duration of interval

      interval_t swint( tp[ zc[i-1] ] , tp[ zc[i] ] - 1 ); 

      const double t = swint.duration_sec();
      
      // duration criteria on whole wave? 

      if ( t_lwr > 0 && 
	   ( t < t_lwr || t > t_upr ) ) continue;
      
      
      // find negative and positive peaks
      
      double x = 100;
      double y = -99;
      int xi = 0, yi = 0, n2pi = 0;
      
      for (int j=zc[i-1];j<zc[i];j++)
	{	  
	  if ( filtered[j] < x ) { x = filtered[j]; xi = j; } 	  
	  if ( filtered[j] > y ) { y = filtered[j]; yi = j; }	  
	}


      // build putative SO object 
      
      slow_wave_t w;
      w.type = type;
      w.interval = interval_t( zc[i-1] , zc[i] );
      w.interval_tp = interval_t( tp[zc[i-1]] , tp[zc[i]] );
      w.down_amplitude = x;
      w.down_peak = tp[xi];
      w.down_peak_sp = xi;
      w.up_amplitude = y;
      w.up_peak = tp[yi];
      w.up_peak_sp = yi;


      // find min/max and also the middle zero-crossing
      // for full waves

      int peak1 = w.down_peak_sp < w.up_peak_sp ?  w.down_peak_sp  : w.up_peak_sp ;
      int peak2 = w.down_peak_sp < w.up_peak_sp ?  w.up_peak_sp  : w.down_peak_sp ;
      
      for (int j=peak1;j<=peak2;j++)
	{	  
	  // neg-to-pos crossing
	  if ( use_alternate_neg2pos_zc ) // look for pos-to-neg
	    {
	      if ( filtered[j-1] >= 0 && filtered[j] < 0 ) n2pi = j;
	    }
	  else // else neg-to-pos for the middle ZC
	    {
	      if ( filtered[j-1] < 0 && filtered[j] >= 0 ) n2pi = j;
	    }
	}

      // should not happen
      if ( n2pi == 0 ) Helper::halt( "problem" );

      w.zero_crossing = n2pi ; 
      w.zero_crossing_tp = tp[n2pi] ; 
      
      
      // duration criteria on negative half wave? 

      if ( t_neg_lwr > 0 ) 
	{
	  
	  if ( use_alternate_neg2pos_zc ) // negative half wave is mid -- stop
	    {	      
	      interval_t hwint( w.zero_crossing_tp , tp[ zc[i] ] - 1 );
	      const double t = hwint.duration_sec();
	      // duration criteria on negative half-wave?
	      if ( t < t_neg_lwr || t > t_neg_upr ) continue;
	    }
	  else // or (DEFAULT) : negative half-wave is start to mid
	    {
	      interval_t hwint( tp[ zc[i-1] ] , w.zero_crossing_tp - 1 );
	      const double t = hwint.duration_sec();
	      // duration criteria on negative half-wave?
	      if ( t < t_neg_lwr || t > t_neg_upr ) continue;
	    }
	}
      
      // accumulate averages for filtering 
      tmp_x.push_back( x );
      tmp_yminusx.push_back( y - x );

      ++cnt;
      
      waves.push_back(w);
      
    } // next putative SW
  

  // either use mean or median (default) as baseline

  if ( use_mean ) 
    {
      avg_x = MiscMath::mean( tmp_x );
      avg_yminusx = MiscMath::mean( tmp_yminusx );
    }
  else
    {
      avg_x = MiscMath::median( tmp_x );
      avg_yminusx = MiscMath::median( tmp_yminusx );
    }

  // get threshold
  th_x       = avg_x * thr;
  th_yminusx = avg_yminusx * thr;

  
  // accumulators for final averages/medians
  
  std::vector<double> acc_yminusx, acc_x, 
    acc_duration_sec, acc_negative_duration_sec, acc_positive_duration_sec, 
    acc_slope_n1 , acc_slope_n2 , acc_slope_p1 , acc_slope_p2 ;
  
  sw.clear();
  
  for (int i = 0; i < waves.size(); i++)
    {
      slow_wave_t & w = waves[i];
      
      bool accepted = true;
      //logger << "thr " << w.down_amplitude  << " " << th_x << " " << uV_neg << " " << w.down_amplitude  << " " << uV_p2p << " " <<w.up_amplitude - w.down_amplitude << "\n";
      if ( thr > 0 && w.down_amplitude > th_x ) accepted = false; // i.e. negative value, should be more neg.
      if ( thr > 0 && w.up_amplitude - w.down_amplitude < th_yminusx ) accepted = false; // pos threshold
      if ( uV_neg < 0 && w.down_amplitude > uV_neg ) accepted = false;
      if ( uV_p2p > 0 && w.up_amplitude - w.down_amplitude < uV_p2p ) accepted = false;
      
      // save this wave?
      if ( accepted ) 
	{
	  sw.push_back( w );

	  // amplitude of negative peak
	  acc_x.push_back( w.down_amplitude );

	  // peak-to-peak amplitude
	  acc_yminusx.push_back( w.up_amplitude - w.down_amplitude );

	  // SO duration ( do not add + to interval, up until last point rather than include last point)
	  acc_duration_sec.push_back( ( w.interval_tp.stop - w.interval_tp.start ) * globals::tp_duration );
	  acc_negative_duration_sec.push_back( ( w.zero_crossing_tp - w.interval_tp.start ) * globals::tp_duration );
	  acc_positive_duration_sec.push_back( ( w.interval_tp.stop - w.zero_crossing_tp ) * globals::tp_duration );

	  // slopes ( returns 0 is not calculated, i.e. because looking at half-waves)
	  const double sn1 = w.slope_n1();
	  const double sn2 = w.slope_n2();
	  const double sp1 = w.slope_p1();
	  const double sp2 = w.slope_p2();
	  
	  if ( sn1 != 0 ) acc_slope_n1.push_back( sn1 );
	  if ( sn2 != 0 ) acc_slope_n2.push_back( sn2 );
	  if ( sp1 != 0 ) acc_slope_p1.push_back( sp1 );
	  if ( sp2 != 0 ) acc_slope_p2.push_back( sp2 );
	  
	  ++cnt;
	}
    }
  
  // means 
  avg_x = acc_x.size() > 0 ? MiscMath::mean( acc_x ) : 0 ; 
  avg_yminusx = acc_yminusx.size() > 0 ? MiscMath::mean( acc_yminusx ) : 0 ; 
  avg_duration_sec = acc_duration_sec.size() > 0 ? MiscMath::mean( acc_duration_sec ) : 0 ;
  avg_negative_duration_sec = acc_negative_duration_sec.size() > 0 ? MiscMath::mean( acc_negative_duration_sec ) : 0 ;
  avg_positive_duration_sec = acc_positive_duration_sec.size() > 0 ? MiscMath::mean( acc_positive_duration_sec ) : 0 ;
  avg_slope_p1 = acc_slope_p1.size() > 0 ? MiscMath::mean( acc_slope_p1 ) : 0 ; 
  avg_slope_p2 = acc_slope_p2.size() > 0 ? MiscMath::mean( acc_slope_p2 ) : 0 ;
  avg_slope_n1 = acc_slope_n1.size() > 0 ? MiscMath::mean( acc_slope_n1 ) : 0 ;
  avg_slope_n2 = acc_slope_n2.size() > 0 ? MiscMath::mean( acc_slope_n2 ) : 0 ; 

  // medians
  median_x = acc_x.size() > 0 ? MiscMath::median( acc_x ) : 0 ;
  median_yminusx = acc_yminusx.size() > 0 ? MiscMath::median( acc_yminusx ) : 0 ;
  median_duration_sec = acc_duration_sec.size() > 0 ? MiscMath::median( acc_duration_sec ) : 0 ;
  median_negative_duration_sec = acc_negative_duration_sec.size() > 0 ? MiscMath::median( acc_negative_duration_sec ) : 0 ;
  median_positive_duration_sec = acc_positive_duration_sec.size() > 0 ? MiscMath::median( acc_positive_duration_sec ) : 0 ;
  median_slope_p1 = acc_slope_p1.size() > 0 ? MiscMath::median( acc_slope_p1 ) : 0 ;
  median_slope_p2 = acc_slope_p2.size() > 0 ? MiscMath::median( acc_slope_p2 ) : 0 ;
  median_slope_n1 = acc_slope_n1.size() > 0 ? MiscMath::median( acc_slope_n1 ) : 0 ; 
  median_slope_n2 = acc_slope_n2.size() > 0 ? MiscMath::median( acc_slope_n2 ) : 0 ; 
    
  logger << ", of which " << sw.size() << " met criteria";
  if ( thr > 0 ) logger << " (thresholds (<x, >p2p) " << th_x << " " << th_yminusx << ")";
  logger << "\n";

  return sw.size();
}


void slow_waves_t::phase_slow_waves()
{
  
  logger << " running Hilbert transform\n";
  
  const int n = filtered.size();

  //
  // get hilbert transform (this is already BP-filtered)
  //

  hilbert_t hilbert( filtered );

  //
  // populate phase for each SW
  //
  
  phase = *hilbert.phase();

  // convert to degrees with 0 as pos-to-neg crossing
  for (int i=0;i<phase.size();i++) phase[i] = MiscMath::as_angle_0_pos2neg( phase[i] );

  // map back to sample points (in SW, Y/N?)
  in_sw.resize( n , -1 );

  for (int i = 0 ; i < sw.size(); i++)
    {
      slow_wave_t & w = sw[i];
      
      w.phase.clear();
      for (int j=w.interval.start;j<=w.interval.stop;j++)
	{
	  w.phase.push_back( phase[j] );
	  in_sw[j] = i;
	}
    }

}

int slow_waves_t::getbin( double x , const std::vector<double> & th , int last_bin , int nb ) const
{    

  if ( last_bin == 0 && x < th[last_bin] ) return 0;
  if ( last_bin >  0 && x >= th[last_bin-1] && x < th[last_bin] ) return last_bin;
  
  // else search, starting with last_bin
  if ( x >= th[last_bin] )
    {
      for (int b=last_bin+1; b<nb; b++)
	if ( x < th[b] ) return b; 
    }
  else
    {
      for (int b=0; b<nb; b++)
	if ( x < th[b] ) return b; 
    }

  // must be largest bin i.e. assuming is within range
  return nb-1;
}

std::vector<double> slow_waves_t::phase_locked_averaging( const std::vector<double> * sig , int nbins , const std::vector<bool> * subset )
{

  if ( sw.size() == 0 ) 
    {
      std::vector<double> sigmean;
      return sigmean;
    }

  std::vector<double> sigmean( nbins , 0 );
  std::vector<int> sigcnt( nbins , 0 );

  // phase is 0..360 degrees, with 0 as pos-2-neg crossing
  
  double inc = 360 / (double)nbins ; 
  double th_ = inc;
  
  std::vector<double> th( nbins );
  for (int i=0;i<nbins;i++) { th[i] = th_; th_ += inc; } 

  int nb = nbins;

  // for each slow wave

  for (int i=0;i<sw.size();i++)
    {

      uint64_t left = sw[i].interval.start;
      uint64_t right = sw[i].interval.stop;
      
      int last_bin = 0;
      
      for (uint64_t p = left ; p <= right ; p++ ) 
	{
	  if ( subset == NULL || (*subset)[p] ) // apply an optional mask
	    {
	      int b = getbin( phase[p] , th , last_bin , nbins );
	      last_bin = b;
	      sigmean[b] += (*sig)[p];
	      ++sigcnt[b];	 
	    }
	}

      // next slow wave
    }

  // get mean
  for (int j=0;j<sigmean.size();j++) sigmean[j] /= (double)sigcnt[j];

  return sigmean;
}


std::vector<double> slow_waves_t::time_locked_averaging( const std::vector<double> * sig , 
							 int sr , 
							 double left, double right , 
							 int position  )
{

  if ( sw.size() == 0 ) 
  {
    std::vector<double> sigmean;
    return sigmean;
  }
  
  // position
  //   0 onset 
  //  -1 negative peak (default)
  //  +1 positive peak
  
  // around each slow wave peak, average up to sr*left points before
  // and sr*right points after

  int nleft  = sr*left;
  int nright = sr*right;
  int np = nleft + 1 + nright;
 
  std::vector<double> sigmean( np , 0 );  
  std::vector<double> sigcnt( np , 0 );

    
  for (int i=0;i<sw.size();i++)
    {
      
      // interval in sample-points, plus window
      // centered around down-peak of each SW
      // by default, but can be onset or 
      // positive peak

      int centre;
      if      ( position == -1 ) centre = sw[i].down_peak_sp;
      else if ( position == 0 ) centre = sw[i].interval.start;
      else if ( position == +1 ) centre = sw[i].up_peak_sp;
      else Helper::halt("internal error in slow_waves_t::time_locked_averaging()" );
      
      // shuffle?
      if ( 0 ) 
	{
	  int offset = CRandom::rand( 3 - 6 * sr );
	  centre += offset;
	  logger << "offset = " << offset << "\n";
	}

      
      int lower = centre - nleft ; 
      int upper = centre + nright ;
      int pos = 0;
      
      for (int j = lower ; j <= upper ; j++ )
	{
	  if ( j < 0 ) { ++pos; continue; }
	  if ( j >= sig->size() ) { ++pos; continue; }
	  if ( pos >= np ) Helper::halt("internal error in slow_waves_t");
	  // track means
	  sigmean[pos] += (*sig)[j];	  
	  sigcnt[pos]++;
	  ++pos;
	  
	  //	  std::cout << "pos " << pos << " " << np << " " << centre << " " << nleft << " " << nright << " " << lower << " " << upper << "\n";
	}
      
      // next slow wave
    }      
  
  // get means
    
  for (int j=0;j<np;j++)
    sigmean[j] /= (double)sigcnt[j];
  
  return sigmean;
  
}

