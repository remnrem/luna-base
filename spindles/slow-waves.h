
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
#include "dynamics/qdynam.h"

#include <vector>
#include <iostream>
#include <sstream>

struct edf_t;

struct param_t;

enum slow_wave_type { SO_FULL , SO_HALF, SO_NEGATIVE_HALF , SO_POSITIVE_HALF } ;

struct slow_wave_param_t {

  slow_wave_param_t() { } ;
  
  slow_wave_param_t( const param_t & param );
  
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

  // fast/slow switcher based on transition freq
  double fast_slow_switcher_th;
  int do_fast_slow;

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

  // output options
  bool out_idx;  // for SO-level outputs
  bool out_all_slopes;  // SLOPE -->  SLOPE_POS1  POS2 NEG1 NEG2 (def = NEG2)

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

  // SO dynamics

  void epoch_dynamics( edf_t * );
  
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

  // copy of parameters/options
  slow_wave_param_t par;
  
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


  // other key options
  bool using_rel;
  
  // dynamics (optional)
  qdynam_t qd;
  bool calc_dynamics;
  
  
  // helper function
  int getbin( double x , const std::vector<double> & th , int last_bin , int nb ) const; 
};


#endif
