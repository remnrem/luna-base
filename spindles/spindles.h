
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
  
  // track mid in 'tp' units (from EDF start)
  uint64_t tp_mid;
  
  // spindle properties
  double amp, dur, fwhm, nosc, frq, fft, symm, symm2, isa;
  double norm_amp_max, norm_amp_mean; 
  double chirp, frq_h1, frq_h2, frq_range; 

  // neg/pos defined by HW
  double posf, negf, posb, negb, posv, negv, allf, allb, allv;
  double posisa, negisa;
  int possp, negsp; // track denom for ISA to get an amplitude measure too (i.e. mean amp for pos/neg)

  // neg/pos defined by slope
  double pos2f, neg2f, pos2b, neg2b, pos2v, neg2v;

  // winsorized? (count of samples)
  int winsor;
  
  // quality score based on enrichments
  double qual;
    
  // flag not to be included in analyis
  bool include;

  // max trough (for spindle temporal alignment)
  int max_trough_sp;
  double max_trough_rel; //(0..1)
  
  // max amplitude (based on CWT)
  int peak_amp_sp;
  double peak_amp_rel; // (0..1)
  
  // SO coupling metrics: anchor may be based on max amp, or a temporal anchor
  double anchor_sec;
  double so_phase_anchor;
  double so_nearest;
  int so_nearest_num;
  
  double if_spindle;

  bool operator<( const spindle_t & rhs ) const { return tp < rhs.tp; } 

};



struct spindle_qc_t {

  spindle_qc_t( param_t & param )
  {

    qc_qmin = 0;
    qc_qmax = -1;
    
    // cwt
    ofc.clear();
    
    // get params
    proc( param );
    
  }
  
  bool qc_q;
    
  double qc_qmin, qc_qmax;

  // original time-series (for verbose outputs)
  const std::vector<double> * ts;
  int Fs;
  
  // CWT : non-spindle bands
  int nob; // number of other-Fcs [ ofc.size() ] 
  std::vector<double> ofc;
  std::vector<double> ocycles;
  std::vector<std::vector<double> > ocwt;

  // spindle stats
  std::vector<double> scwt; // varies for each Fc

  // verbose output
  bool verbose_failed;
  bool verbose_all; 
  
  // functions
  void proc( param_t & );
  void init( const std::vector<double> * ts , int Fs );
  void init_spindle( const std::vector<double> & scwt );
  double qual( spindle_t * spindle ); 

  // verbose output
  void output( spindle_t * spindle , const int n );
  
};



annot_t * spindle_wavelet( edf_t & , param_t & );

// helper function for FFT
void do_fft( const std::vector<double> * d , const int Fs , spindle_qc_t * q , bool baseline );
void do_fft( const std::vector<double> * d , const int Fs , const int np, spindle_qc_t * q , bool baseline );
  
// helper to get spindle stats
void spindle_stats( const std::vector<spindle_t> & spindles , std::map<std::string,double> & ) ;

// helper
//void write_if_exists( const std::string & s , const std::map<std::string,double> & means ) ;

void characterize_spindles( edf_t & edf , 
			    param_t & param , 
			    const int s , 
			    bool bandpass_filtered , 
			    const double target_f, 
			    const std::string & label, 			    
			    const std::vector<double> * averaged ,			    
			    const std::vector<double> * original_signal , 
			    std::vector<spindle_t>    * spindles , 			    
			    clocktime_t * starttime ,
			    spindle_qc_t * q = NULL ,  
			    std::map<double,double> * locked = NULL , 
			    std::vector<bool> * in_pos_hw = NULL , 
			    std::vector<bool> * in_neg_hw = NULL , 
			    std::vector<bool> * in_pos_slope = NULL , 
			    std::vector<bool> * in_neg_slope = NULL ,
			    std::vector<bool> * p_winsp = NULL 
			    );


void per_spindle_output( std::vector<spindle_t>    * spindles ,
			 param_t & param , 
			 clocktime_t               * starttime ,
			 spindle_qc_t              & q );


//
// Spindle / SO coupling
//


void spindle_so_coupling( edf_t & , param_t & );


#endif
