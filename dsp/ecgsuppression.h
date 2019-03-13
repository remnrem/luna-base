
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

#ifndef __ECGSUPPRESS_H__
#define __ECGSUPPRESS_H__

struct edf_t;
struct param_t;
struct interval_t;

#include <vector>
#include <stdint.h>
#include <cstddef>

struct rpeaks_t
{
  
  rpeaks_t( const int n )
  {
    R_t.resize( n );
    R_i.resize( n );
    R_amp.resize( n );
    S_amp.resize( n );
  }
  
  std::vector<uint64_t> R_t;
  std::vector<uint64_t> R_i;
  std::vector<double> R_amp;
  std::vector<double> S_amp;
  
  double bpm( interval_t & , double lwr = 0 , double upr = 0 ) ;
  
  int clean( double , double );
  
  std::vector<uint64_t> beats( interval_t & ) const;
  
};


namespace dsptools 
{  
  
  // ultimate wrapper
  void ecgsuppression( edf_t & , param_t & );

  // just print BPM per epoch
  void bpm( edf_t & , param_t & );

  // find ECG peaks, based on:
  // http://www.robots.ox.ac.uk/~gari/CODE/ECGtools/ecgBag/rpeakdetect.m
  
  rpeaks_t mpeakdetect( const edf_t & edf , 
			const std::vector<double> * ecg , 
			const std::vector<uint64_t> * tp , 
			int Fs , 
			const std::vector<double> * eeg = NULL , 
			bool * force = NULL );
  
  /* int mask_ecg( edf_t & edf ,  */
  /* 		const rpeaks_t & rpeaks ); */
  
  //  void suppress_ecg( edf_t & edf , 
		     
		  
}

#endif
