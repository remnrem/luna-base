
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


#ifndef __SPINDLES_H__
#define __SPINDLES_H__

struct edf_t;
struct param_t;
struct annot_t;
struct clocktime_t;

#include "intervals/intervals.h"
#include <vector>


struct spindle_t
{

  spindle_t( uint64_t a, uint64_t b , int asp, int bsp )
    : start_sp( asp) , stop_sp( bsp ), tp( interval_t( a , b ) ) 
  { 
    include = true;
  } 

  // have start and stop in sample-points as well as time-points
  int start_sp, stop_sp;
  interval_t tp;
  
  // spindle properties
  double amp, dur, fwhm, nosc, frq, fft, symm, symm2, trend, abstrend, chirp, isa;
  
  // relative enrichment per frequency range compared to whole trace baseline
  std::map<freq_range_t,double> enrich;

  // quality score based on enrichments
  double qual;
  
  // detection
  double max_stat, mean_stat;
  
  // spindle peak (point of lowest trough)
  int peak_sp;
  
  // flag not to be included in analyis
  bool include;

  bool operator<( const spindle_t & rhs ) const { return tp < rhs.tp; } 

};


annot_t * spindle_bandpass( edf_t & , param_t & );

annot_t * spindle_wavelet( edf_t & , param_t & );


// helper function for FFT
void do_fft( const std::vector<double> * d , const int Fs , std::map<freq_range_t,double> * fft );

// helper to get spindle stats
std::map<std::string,double> spindle_stats( const std::vector<spindle_t> & spindles ) ;

void characterize_spindles( edf_t & edf , 
			    param_t & param , 
			    const int s , 
			    bool bandpass_filtered , 
			    const double target_f, 
			    const double window_f , 
			    const std::string & label, 			    
			    const std::vector<double> * averaged , 
			    const std::vector<double> * original_signal , 
			    std::vector<spindle_t>    * spindles , 			    
			    clocktime_t * starttime , 
			    std::map<freq_range_t,double> * baseline = NULL , 
			    std::map<double,double> * locked = NULL 
			    );


#endif
