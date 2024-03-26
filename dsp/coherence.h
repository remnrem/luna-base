
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

#ifndef __COH_H__
#define __COH_H__

#include <vector>
#include <complex>
#include <map>
#include "defs/defs.h"

struct edf_t;
struct param_t;
struct interval_t;
struct coherence_t;

// results of a single cross-spectral analysis
// either for 1 epoch, or for the whole signal

struct scoh_t
{
  
  scoh_t() { } 
  
  scoh_t( const int nf )
  {
    resize(nf);
  }
  
  void resize( const int n )
  {
    bad.resize( n , false );
    sxx.resize( n );
    syy.resize( n );
    sxy.resize( n );    
  }
  
  void proc_and_output( const coherence_t & , 
			const bool output , 
			const double upper_freq = -1 , 
			const double lower_freq = 0 );
  
  // cross and auto spectra (vector over frequencies)
  std::vector<bool>   bad;
  std::vector<double> sxx;
  std::vector<double> syy;
  std::vector<std::complex<double> > sxy;
  
  // band power
  std::map<frequency_band_t,double> bcoh, bicoh, blcoh;
  // band bin count
  std::map<frequency_band_t,int> bn;
   

};


struct coh_t
{

  // accumulate (over epochs)

  void add( const scoh_t & c ) { epochs.push_back( c ); } 
  
  // calculate and output final (averaged) connectivity stats, 
  void calc_stats( const coherence_t & , 
		   const double upper_freq = -1 , 
		   const double lower_freq = 0 );
    
  // data 
  std::vector<scoh_t> epochs;
  
};



namespace dsptools 
{  

  // ultimate wrapper
  void coherence( edf_t & , param_t & );

  void coherence_prepare( edf_t & edf , const int signal1 , const interval_t & interval , coherence_t * coh );
    
  // for a given pair of signals, return coh object, either using 
  //  coherence_t()  (from fttwrap/) 
  
  scoh_t coherence_do( coherence_t * , int signal1 , int signal2 );
  
  
}

#endif
