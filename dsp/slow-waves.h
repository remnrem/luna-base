
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

#include <vector>
#include <iostream>
#include <sstream>

struct edf_t;
struct param_t;

enum slow_wave_type { SO_FULL , SO_HALF, SO_NEGATIVE_HALF , SO_POSITIVE_HALF } ; 

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

  //  bool       neg2pos;   // find neg-to-pos zerocrossing?   otherwise, pos-to-negative

  std::vector<double> phase;

  double amplitude() const { return fabs( up_amplitude ) + fabs( down_amplitude ) ; } 

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

  slow_waves_t() { } ;

  slow_waves_t( edf_t & , const param_t & );

  slow_waves_t( const std::vector<double> & d , 
		const std::vector<uint64_t> & tp ,
		const int sr , 
		const double thr   = 0.75 , 
		const bool   use_mean = false , 
		const double uv_neg = 0 , 
		const double uv_p2p = 0 ,
		const double f_lwr = 0.5 ,    
		const double f_upr = 4 , 
		const double t_lwr = 0.8 ,
		const double t_upr = 2 , 		
		const double t_neg_lwr = 0 , 
		const double t_neg_upr = 0 , 
		const bool   neg2pos = true , 
		const slow_wave_type type = SO_FULL );

		
  // actual detection

  int detect_slow_waves( const std::vector<double> & d , 
			 const std::vector<uint64_t> & tp ,

			 // sample rate
			 const int sr , 

			 // relative threshold based on mean of 
			 // all SOs (based on negative peak)
			 // this might not be used (i.e. set to 0) 
			 // if only using absolute criteria
			 
			 const double thr   = 0.75 , 
			 
			 // use mean versus median 
			 const bool   use_mean = false , 
			 
			 // absolute uV threshold for the negative peak (x < th ) 
			 const double uV_neg = 0 ,   

			 // absolute uV threshold for peak-to-peak
			 const double uV_p2p = 0 ,   
			 
			 // transition frequencies for BPF signal
			 const double f_lwr = 0.5 , 
			 const double f_upr = 4 , 
			 
			 // duration thresholds for entire SW
			 const double t_lwr = 0.8 ,
			 const double t_upr = 2 , 
			 
			 // duration of negative deflection only
			 // probably not used if using the above
			 // total wave criteria
			 const double t_neg_lwr = 0 ,  // 0.125
			 const double t_neg_upr = 0 ,  // 1.5  

			 // count negative-to-positive zero-crossings
			 // as opposed to positive-to-negative
			 const bool neg2pos = true , 
			 
			 const slow_wave_type type = SO_FULL );
  
  // output

  void display_slow_waves( const bool verbose = false , edf_t * edf = NULL );

  // analysis of slow waves
  
  void phase_slow_waves();
  
  void time_locked_spectral_analysis( edf_t & edf , const std::vector<double> * sig , double sr, double window_sec = 1.5  ); 


  // time-locked and phase-locked averaging
  
  std::vector<double> time_locked_averaging( const std::vector<double> * sig , int sr , double left, double right , int position = -1 );
  std::vector<double> phase_locked_averaging( const std::vector<double> * sig , int bins , const std::vector<bool> * subset = NULL );
  

  // reporting functions

  std::vector<slow_wave_t> waves() const { return sw; }
  
  bool in_slow_wave( const int i ) const { return i < 0 || i >= in_sw.size() ? false : in_sw[i] != -1 ; } 
  
  double nearest( const int i , int * ) const; 
  
  int slow_wave_number( const int i ) const { return i < 0 || i >= in_sw.size() ? false : in_sw[i]; } 

  int num_waves() const { return sw.size(); } 

  void time_locked_spectral_power( const std::vector<bool> * included = NULL );

  const std::vector<double> * p_filtered() const { return &filtered; }  
  

private:

  // store all individual slow waves
  std::vector<slow_wave_t> sw;
  
  // simple vector of 0/n for whether this point is in a slow wave ('n') or no (0)
  std::vector<int> in_sw;

  // filtered signal
  std::vector<double> filtered;
  std::vector<uint64_t> tp;

  // SW phase for all sample-points
  std::vector<double> phase;
  
  // detection thresholds
  double th_x;        // negative peak (-ve)
  double th_yminusx;  // p2p

  // total signal time in seconds (i.e. denonimnator for SW rate)
  double signal_duration_sec;
  
  double avg_x; // minimum amplitude  SW_SMP
  double avg_yminusx; // peak-to-peak SW_P2P
  double avg_duration_sec; // length  SW_DUR
  double avg_negative_duration_sec; // length, negative half-wave  SW_NEG_DUR
  double avg_positive_duration_sec; // length, positive half-wave  SW_POS_DUR
  double avg_slope_n1; // down-going for neg peak   SW_SLOPE_NEG1
  double avg_slope_n2; // up-going for neg peak     SW_SLOPE_NEG2 ** main
  double avg_slope_p1; // up-going for pos peak     SW_SLOPE_POS1 
  double avg_slope_p2; // down-going for pos peak   SW_SLOPE_POS2

  // median versions
  double median_x; // minimum amplitude  SW_SMP
  double median_yminusx; // peak-to-peak SW_P2P
  double median_duration_sec; // length  SW_DUR
  double median_negative_duration_sec; // length, negative half-wave  SW_NEG_DUR
  double median_positive_duration_sec; // length, positive half-wave  SW_POS_DUR
  double median_slope_n1; // down-going for neg peak   SW_SLOPE_NEG1
  double median_slope_n2; // up-going for neg peak     SW_SLOPE_NEG2 ** main
  double median_slope_p1; // up-going for pos peak     SW_SLOPE_POS1 
  double median_slope_p2; // down-going for pos peak   SW_SLOPE_POS2

  int Fs;

  // helper function
  int getbin( double x , const std::vector<double> & th , int last_bin , int nb ) const; 
};


#endif
