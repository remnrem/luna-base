
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


#ifndef __COHFFT_H__
#define __COHFFT_H__

#include "fftw3.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <complex>

#include "helper/helper.h"
#include "helper/logger.h"
#include "miscmath/miscmath.h"
#include "defs/defs.h"
#include "dsp/coherence.h"

extern logger_t logger;

struct edf_t;

struct coherence_t;

struct precoh_t {

  // pre-calculate and store channel-wise spectra
  // map-->freq-->segment-->FFT
  
  std::map<int,std::vector<std::vector<std::complex<double> > > > psd;

  std::vector<double> frq;
  
  double normalisation_factor;

  int cutoff;

  void clear() {
    psd.clear();
    frq.clear();
    cutoff = 0;
    normalisation_factor = 1.0 ;     
  }

  void prepare( coherence_t * ,
		const int s ,
		const std::vector<double> & x );   
  
  //   psd_x[i].push_back( ( a*a + b*b ) * normalisation_factor );
  // std::complex<double> Xx( a , b );
  
  // psd_y[i].push_back( ( a*a + b*b ) * normalisation_factor );
  // std::complex<double> Yy( a , b );
  
  // std::complex<double> Xy = Xx * conj( Yy );
  // cpsd[i].push_back( normalisation_factor * Xy );
  
  
};


class coherence_t {

  friend struct precoh_t;
  
 public:

  coherence_t ( int total_points ,
		int Fs, 
		double segment_sec ,  // segment size in seconds 
		double overlap_sec , // overlap in seconds
		window_function_t W = WINDOW_HANN , 
		bool average_adj = false , 
		bool detrend = false , 
		bool zerocenter = false )		
    :  total_points( total_points) , segment_sec(segment_sec) , overlap_sec(overlap_sec),
       Fs(Fs), window(W), detrend(detrend) , zerocenter(zerocenter) , average_adj(average_adj) 
  {
    
    // segment parameters in sample points (NFFT)
    segment_points = segment_sec * Fs;
    
    // overlap in sample points
    noverlap_points1  = overlap_sec * Fs;
    
    // calculate implied overlap in actual data-points
    noverlap_segments = floor( ( total_points - noverlap_points1) 
			       / (double)( segment_points - noverlap_points1 ) );
    
    noverlap_points2          = noverlap_segments > 1 
    ? ceil( ( noverlap_segments*segment_points - total_points  ) / double( noverlap_segments - 1 ) )
    : 0 ;

    // points1 -- will truncate end of signal if the segment length plus overlap is not an exact match
    // points2 -- will reduce overlap so that last intervals sits on end of signal, i.e. ensures that 
    //            whole signal is evenly covered;  probably this is a better one to use as a default
    
    const bool ensure_even_signal_coverage = true;
    
    if ( ensure_even_signal_coverage )
      segment_increment_points = segment_points - noverlap_points2;
    else
      segment_increment_points = segment_points - noverlap_points1;
    
/*          std::cerr << "noverlap_points1\t" << noverlap_points1 << "\n";  */
/*          std::cerr << "noverlap_points2\t" << noverlap_points2 << "\n";  */
/*          std::cerr << "segment_increment_points\t" << segment_increment_points << "\n";  */
/*          std::cerr << "noverlap_segments\t" << noverlap_segments << "\n";  */
/*          std::cerr << "segment_points\t" << segment_points << "\n";  */
        
    //process(); 
  } 

  double segment_sec, overlap_sec;
  
  int total_points;
  int noverlap_segments;
  int segment_points;
  int noverlap_points1, noverlap_points2;
  int segment_increment_points;

  int Fs;
  
  scoh_t res;

  const std::vector<double> & frq() const { return precoh.frq; } 
  
 private:

  window_function_t window;
  
  bool detrend;

  bool zerocenter;

  bool average_adj;

  int N;

  static precoh_t precoh;

  
public:


  void prepare( const int s , const std::vector<double> & x )
  {
    precoh.prepare( this , s , x );
  }

  void clear()
  {
    precoh.clear();
  }
  
  void process( const int , const int );

  
};



#endif
