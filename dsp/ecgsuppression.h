
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
#include <map>
#include <set>

struct rpeaks_t
{
  
  std::vector<uint64_t> R_t; // tp
  std::vector<uint64_t> R_i; // matching sp-index on ECG

  double npks;
  double p_inverted;
  bool inverted;
  
  double bpm( interval_t & , double lwr = 0 , double upr = 0 ) ;
  
  int clean( double , double );
  
  std::vector<uint64_t> beats( interval_t & ) const;

  // helper for annot-stratified HRV analyses
  static rpeaks_t intersect( const std::set<uint64_t> & pks ,
			     const std::set<interval_t> & ints );
  
};


struct hrv_opt_t;

struct hrv_res_t {

  hrv_res_t()
    : IMPUTED(0.0), P_INV(0.0), INV(0.0), NP(0.0), NP_TOT(0.0),
      RR(0.0), HR(0.0),
      SDNN(0.0), SDNN_R(0.0), RMSSD(0.0), RMSSD_R(0.0), pNN50(0.0),
      LF(0.0), LF_N(0.0), LF_PK(0.0), HF(0.0), HF_N(0.0), HF_PK(0.0),
      LF2HF(0.0) {}
  
  
  double IMPUTED; // proportion of intervals imputed (out-of-range)
  double P_INV;   // proportion of intervals inverted
  double INV;     // is it inverted? Y/N
  double NP;      // # of RR intervals
  double NP_TOT;  // keep as total number of peaks (i.e. sum over epochs) 
  double RR;      // mean RR interval (msec)
  double HR;      // heart rate (BPM)

  double SDNN;    // SD of normal-to-normal RR intervals
  double SDNN_R;    // SD of normal-to-normal RR intervals
  double RMSSD;   // RMS of successive diffs
  double RMSSD_R;   // RMS of successive diffs
  double pNN50;   // prop. of successive N-N intervals differing by more than 50 msec 

  double LF;
  double LF_N;
  double LF_PK;
  double HF;
  double HF_N;
  double HF_PK;
  double LF2HF;
  
  static hrv_res_t summarize( const std::vector<hrv_res_t> & x );
  
  void write( const hrv_opt_t & , const bool reduced = false ) const;
   
};


struct rr_intervals_t {
  
  rr_intervals_t( const rpeaks_t & ,
		  const hrv_opt_t & );
  
  hrv_res_t res;
  
  std::vector<double> rr;
  std::vector<double> t;
  std::vector<uint64_t> tp;
  std::vector<bool> imputed;
  
};

struct rpeak_opt_t;

namespace dsptools 
{  
  
  // ultimate wrapper
  void ecgsuppression( edf_t & , param_t & );

  // just print BPM per epoch
  void bpm( edf_t & , param_t & );

  // Welch-based HRV implementation
  void hrv( edf_t & , param_t & );
  
  // find ECG peaks, based on:
  // http://www.robots.ox.ac.uk/~gari/CODE/ECGtools/ecgBag/rpeakdetect.m
  
  rpeaks_t mpeakdetect( const edf_t & edf , 
			const std::vector<double> * ecg , 
			const std::vector<uint64_t> * tp , 
			int Fs , 
			const std::vector<double> * eeg = NULL , 
			bool * force = NULL );

  rpeaks_t mpeakdetect2( const std::vector<double> * ecg , 
			 const std::vector<uint64_t> * tp , 
			 int Fs ,
			 const rpeak_opt_t & opt );
  
		  
}

#endif
