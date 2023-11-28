
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


#ifndef __SLOW_WAVES_H__
#define __SLOW_WAVES_H__

#include "intervals/intervals.h"
#include "timeline/cache.h"
#include "annot/annot.h"

#include <vector>
#include <iostream>
#include <sstream>


struct edf_t;
struct param_t;

enum slow_wave_type { SO_FULL , SO_HALF, SO_NEGATIVE_HALF , SO_POSITIVE_HALF } ;

struct slow_wave_param_t {

  slow_wave_param_t( const param_t & param )
  {
    
    // freq.
    f_lwr = param.has( "f-lwr" ) ? param.requires_dbl( "f-lwr" ) : 0.5  ; 
    f_upr = param.has( "f-upr" ) ? param.requires_dbl( "f-upr" ) : 4.0  ; 
    
    // time - for full wave
    t_lwr = param.has( "t-lwr" ) ? param.requires_dbl( "t-lwr" ) : 0 ;
    t_upr = param.has( "t-upr" ) ? param.requires_dbl( "t-upr" ) : 2 ;

    // time - separately for neg and pos HWs
    t_neg_lwr = param.has( "t-neg-lwr" ) ? param.requires_dbl( "t-neg-lwr" ) : 0;
    t_neg_upr = param.has( "t-neg-upr" ) ? param.requires_dbl( "t-neg-upr" ) : 0;

    t_pos_lwr = param.has( "t-pos-lwr" ) ? param.requires_dbl( "t-pos-lwr" ) : 0;
    t_pos_upr = param.has( "t-pos-upr" ) ? param.requires_dbl( "t-pos-upr" ) : 0;

    // amplitude thresholds

    // relative: magnitude threshold, e.g. 2 = twice the mean for a) p2p amp. AND , b) negative peak amp
    thr = param.has( "mag" ) ? param.requires_dbl( "mag" )       : 0 ;
    using_rel = thr > 0 ;
    
    // or, base this on median of all events rather than mean
    use_mean = param.has( "th-mean" ) ;

    // only look at p2p AMP
    ignore_neg_peak = param.has( "ignore-neg-peak" ) ? Helper::yesno( param.value( "ignore-neg-peak" ) ) : false; 

    // fixed
    uV_neg = param.has("uV-neg" ) ? param.requires_dbl( "uV-neg" ) : 0 ; 
    if ( uV_neg > 0 ) Helper::halt( "uV-neg should be negative" ) ;

    uV_p2p = param.has("uV-p2p" ) ? param.requires_dbl( "uV-p2p" ) : 0 ; 
    if ( uV_p2p < 0 ) Helper::halt( "uV-p2p should be positive" ) ;


    //
    // SO/delta distinctions
    //

    // e.g. Kim et al neg bottom 15% (of -ve values),  pos = top 40% (of +ve values)
    pct_neg = param.has( "pct-neg" ) ? param.requires_dbl( "pct-neg" ) / 100.0 : -1 ;
    pct_pos = param.has( "pct-pos" ) ? param.requires_dbl( "pct-pos" ) / 100.0 : -1 ; 
    
    if ( pct_neg > 1 ) Helper::halt( "pct-neg should be between 0 and 100" );
    if ( pct_pos > 1 ) Helper::halt( "pct-pos should be between 0 and 100" );

    // percentile for neg & P2P (as per mag)
    pct = param.has( "pct" ) ? param.requires_dbl( "pct" ) / 100.0 : -1;
    
    // Kim et al. peak-to-peak time interval : 0.15  0.5 
    t_p2p_min = param.has( "t-p2p-min" ) ? param.requires_dbl( "t-p2p-min" ) : 0 ;
    t_p2p_max = param.has( "t-p2p-max" ) ? param.requires_dbl( "t-p2p-max" ) : 0 ;

    SO_delta_mode = 0;
    if ( param.has( "SO-only" ) ) SO_delta_mode = 1 ;
    if ( param.has( "delta-only" ) )
      {
	if ( SO_delta_mode == 1 ) 
	  Helper::halt( "cannot specify both SO-only and delta-only" );
	SO_delta_mode = 2 ;
      }

    // default FIR settings for filter-Hilbert
    fir_ripple = param.has( "sw-ripple" ) ? param.requires_dbl( "sw-ripple" ) : 0.01 ;
    fir_tw = param.has( "sw-tw" ) ? param.requires_dbl( "sw-tw" ) : 0.5 ; 
    
    // legacy/ignored
    pos2neg_zc = ! param.has( "neg2pos" ) ;
    
    // legacy/ignored
    type = SO_FULL;
    if      ( param.has( "half-wave" ) ) type = SO_HALF;
    else if ( param.has( "negative-half-wave" ) ) type = SO_NEGATIVE_HALF;
    else if ( param.has( "positive-half-wave" ) ) type = SO_POSITIVE_HALF;

    // annotations -- first check so-annot then annot (i.e. if called from SPINDLES)
    astr = ".";

    if ( param.has( "so-annot" ) )
      astr = param.value( "so-annot" );
    else if ( param.has( "annot" ) )
      astr = param.value( "annot" );
    
    // do not skip SO detection
    skip = false;
  }
  
  // relative threshold based on mean of 
  // all SOs (based on negative peak & P2P )
  // this might not be used (i.e. set to 0) 
  // if only using absolute criteria
  double thr; // 0.75 
  bool using_rel; // if thr > 0 

  // if using thr, then only base on P2P
  // i.e. if signal polarity is uncertain
  
  bool ignore_neg_peak ; // = false , 

  // use mean versus median 
  bool   use_mean ; // = false , 
  
  // absolute uV threshold for the negative peak (x < th ) 
  double uV_neg ; // = 0    
  
  // absolute uV threshold for peak-to-peak
  double uV_p2p ; // = 0 ,   
  
  // transition frequencies for BPF signal
  double f_lwr ; // = 0.5 , 
  double f_upr ; // = 4 , 
    
  // duration thresholds for entire SW
  double t_lwr = 0.8;
  double t_upr = 2; 
  
  // duration of negative deflection only
  // probably not used if using the above
  // total wave criteria
  double t_neg_lwr ; // = 0 ,  // 0.125
  double t_neg_upr; //  = 0 ,  // 1.5  
  
  // as above, but positive half-wave
  double t_pos_lwr; // = 0 ,  
  double t_pos_upr; // = 0 ,  

  // based on SO/delta distinction
  double pct_neg ; // = 0 ,  // DOWN percentile
  double pct_pos; //  = 0 ,  // UP percentile
  double pct; // DOWN and P2P percentile
  double t_p2p_min; //  = 0  ,  // DOWN peak - UP peak min time
  double t_p2p_max; //  = 0  ,  // DOWN peak - UP peak max time 
  int    SO_delta_mode; // = 0 , // 0 ignore, 1=SO, 2=Delta
			 			 
  // default FIR settings for filter-Hilbert
  double fir_ripple ; // = 0.01 ,
  double fir_tw ; // = 0.5 			 

  // redundant/ignored for now

  // default is to find ZC pairs that are pos2neg (i.e. DOWN, then UP)
  bool pos2neg_zc ;  // = T 
  
  // SW type (ignored)
  slow_wave_type type; // = SO_FULL ,

  // annotations
  std::string astr;
  
  // current channel
  std::string ch;

  // skip SO detection per se
  bool skip;
  
};
  



struct slow_wave_t
{
  
  slow_wave_type type;
  
  interval_t interval;     //sample-units
  interval_t interval_tp;  //time-point units
      
  uint64_t   zero_crossing;   // sample points
  uint64_t   zero_crossing_tp; // time points

  double     up_amplitude;
  double     down_amplitude;
  
  uint64_t   down_peak, up_peak;
  int        down_peak_sp, up_peak_sp;

  int        SO_delta; // 0=NA, 1=SO, 2=delta
  
  std::vector<double> phase;

  double amplitude() const { return up_amplitude + fabs( down_amplitude ) ; }

  double pos_amplitude() const { return up_amplitude ; }

  double neg_amplitude() const { return fabs( down_amplitude ) ; }
  
  double slope_n1() const { 
    if ( type == SO_POSITIVE_HALF ) return 0;
    return down_amplitude / ( (double)(down_peak - interval_tp.start + 1LL ) * globals::tp_duration ); 
  } 

  double slope_n2() const { 
    if ( type == SO_POSITIVE_HALF ) return 0;
    return -down_amplitude / ( (double)(zero_crossing_tp - down_peak + 1LL ) * globals::tp_duration ); 
  } 

  double slope_p1() const { 
    if ( type == SO_NEGATIVE_HALF ) return 0;
    return up_amplitude / ( (double)(up_peak - zero_crossing_tp + 1LL ) * globals::tp_duration );
  } 

  double slope_p2() const { 
    if ( type == SO_NEGATIVE_HALF ) return 0;
    return - up_amplitude / ( (double)(interval_tp.stop - up_peak + 1LL ) * globals::tp_duration );
  } 


  double dur() const {
    return ( interval_tp.stop - interval_tp.start ) * globals::tp_duration ;
  }
  
  double mid() const {
    return zero_crossing_tp * globals::tp_duration ;    
  }

  // start --> mid
  double dur1() const {
    return ( zero_crossing_tp - interval_tp.start ) * globals::tp_duration ;    
  }
  
  // mid --> stop
  double dur2() const {
    return ( interval_tp.stop - zero_crossing_tp ) * globals::tp_duration ;
  }  

  // neg --> pos
  double trans() const { 
    return ( (double)(up_peak - down_peak + 1LL ) * globals::tp_duration );
  }
  
  // neg -> post transition as freq
  double trans_freq() const {
    return 1.0 / ( 2 * trans() ) ; 
  }

  bool is_SO() const { return SO_delta == 1; }

  bool is_delta() const { return SO_delta == 2; }
    
  std::string print() const 
  {
    std::stringstream ss;
    ss << interval_tp << " "
       << zero_crossing << " "
       << up_amplitude << " "
       << down_amplitude << " "
       << phase.size() << " (";
    for (int p=0;p<phase.size();p++) ss << " " << phase[p] ;
    ss << " )";
    return ss.str();
  }

};

struct slow_waves_t
{

  slow_waves_t() {
    report_median_stats = false;
  } ;

  slow_waves_t( edf_t & , const param_t & );

  slow_waves_t( const std::vector<double> & d , 
		const std::vector<uint64_t> & tp ,
		const int sr , 
		const slow_wave_param_t & par , 
		const std::string * cach_name_neg = NULL ,
		const std::string * cach_name_pos = NULL ,
		edf_t * edf = NULL 
		);

		
  // actual detection

  int detect_slow_waves( const std::vector<double> & d , 
			 const std::vector<uint64_t> & tp ,

			 // sample rate
			 const int sr , 

			 // SW detection parameters
			 const slow_wave_param_t & par ,			 
			 
			 // add peaks to a cache? (if label is non-null)
			 // (requires edf-> also)
			 const std::string * cache_name_neg = NULL ,

			 const std::string * cache_name_pos = NULL , 
			 
			 edf_t * edf = NULL
			 
			 );
  
  // output

  void display_slow_waves( const bool verbose = false , edf_t * edf = NULL );

  // analysis of slow waves
  
  void phase_slow_waves();
  
  // time-locked and phase-locked averaging
  
  std::vector<double> time_locked_averaging( const std::vector<double> * sig , int sr , double left, double right , int position = -1 );

  std::vector<double> phase_locked_averaging( const std::vector<double> * sig , int bins , 
					      const std::vector<bool> * subset = NULL , 
					      std::vector<int> * psigcnt = NULL );
  

  // reporting functions

  std::vector<slow_wave_t> waves() const { return sw; }
  
  bool in_slow_wave( const int i ) const { return i < 0 || i >= in_sw.size() ? false : in_sw[i] != -1 ; } 
  
  double nearest( const int i , int * ) const; 
  
  int slow_wave_number( const int i ) const { return i < 0 || i >= in_sw.size() ? false : in_sw[i]; } 

  int num_waves() const { return sw.size(); } 

  void time_locked_spectral_power( const std::vector<bool> * included = NULL );

  const std::vector<double> * p_filtered() const { return &filtered; }  
  
  std::vector<bool> sp_in_sw_vec() { 
    std::vector<bool> r( in_sw.size() ) ;
    for (int i=0;i<in_sw.size();i++) r[i] = in_sw[i] != -1;
    return r;
  } 
  
private:

  // store all individual slow waves
  std::vector<slow_wave_t> sw;
  
  // simple vector of 0/n for whether this point is in a slow wave ('n') or no (-1)
  std::vector<int> in_sw;

  // filtered signal
  std::vector<double> filtered;
  std::vector<uint64_t> tp;

  // SW phase for all sample-points
  std::vector<double> phase;
  
  // detection thresholds
  double th_x;        // negative peak (-ve)
  double th_y;        // positive peak (+ve)
  double th_yminusx;  // p2p
  
  // percentiles
  double th_pct_x;       // neg peak (-ve value)
  double th_pct_y;       // pos peak (+ve value)
  double th_pct_yminusx; // p2p

  // annotations
  std::string astr;
  bool output_halfwave_annots;

  // current channel label
  std::string ch;
  
  // total signal time in seconds (i.e. denonimnator for SW rate)
  double signal_duration_sec;
  
  // use mean or median? 
  
  bool report_median_stats;
  
  double avg_x; // minimum amplitude  SW_NEG_AMP
  double avg_y; // maximum amplitude  SW_POS_AMP
  double avg_yminusx; // peak-to-peak SW_P2P
  double avg_duration_sec; // length  SW_DUR
  double avg_negative_duration_sec; // length, negative half-wave  SW_NEG_DUR
  double avg_positive_duration_sec; // length, positive half-wave  SW_POS_DUR
  double avg_slope_n1; // down-going for neg peak   SW_SLOPE_NEG1
  double avg_slope_n2; // up-going for neg peak     SW_SLOPE_NEG2 ** main
  double avg_slope_p1; // up-going for pos peak     SW_SLOPE_POS1 
  double avg_slope_p2; // down-going for pos peak   SW_SLOPE_POS2
  double avg_trans, avg_trans_freq;

  // median versions
  double median_x; // minimum amplitude  SW_NEG_AMP
  double median_y; // maximum amplitude  SW_POS_AMP
  double median_yminusx; // peak-to-peak SW_P2P
  double median_duration_sec; // length  SW_DUR
  double median_negative_duration_sec; // length, negative half-wave  SW_NEG_DUR
  double median_positive_duration_sec; // length, positive half-wave  SW_POS_DUR
  double median_slope_n1; // down-going for neg peak   SW_SLOPE_NEG1
  double median_slope_n2; // up-going for neg peak     SW_SLOPE_NEG2 ** main
  double median_slope_p1; // up-going for pos peak     SW_SLOPE_POS1 
  double median_slope_p2; // down-going for pos peak   SW_SLOPE_POS2
  double median_trans, median_trans_frq; 

  int Fs;

  
  // helper function
  int getbin( double x , const std::vector<double> & th , int last_bin , int nb ) const; 
};


#endif
