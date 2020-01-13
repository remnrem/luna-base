
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


#ifndef __FFTWRAP_H__
#define __FFTWRAP_H__

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

#include "fftw/cohfft.h"

extern logger_t logger;

struct edf_t;

std::map<double,double> fft_spectrum( const std::vector<double> * d , int Fs );

class FFT
{
  
  friend class coherence_t;
  friend class precoh_t;

 public:

  FFT( int N , int Fs , fft_t type = FFT_FORWARD , window_function_t window = WINDOW_NONE );
  
  ~FFT() 
    {    
      fftw_destroy_plan(p);
      fftw_free(in);
      fftw_free(out);
    }
  
 private:

  // Size (NFFT)
  int N;
  
  // Sampling rate, so we can construct the appropriate Hz for the PSD
  int Fs;
  
  // Input signal
  fftw_complex *in;

  // Output signal
  fftw_complex *out;
  
  // FFT plan from FFTW3
  fftw_plan p;

  // Forward or inverse FFT?
  fft_t type;
  
  // Optional windowing function
  window_function_t window;
  std::vector<double> w;

  // Normalisation factor given the window
  double normalisation_factor;
  

  //
  // Power band helper functions
  //

  static bool add( frequency_band_t band , double f );
  
  double width( frequency_band_t band );

 public:
  
  int cutoff;
  std::vector<double> X;
  std::vector<double> mag;
  std::vector<double> frq;
  
 public:
  
  //static std::map<double,double> power_spectra( edf_t & , const std::string & signal );
  
  bool apply( const std::vector<double> & x );
  bool apply( const double * x , const int n );
  bool apply( const std::vector<std::complex<double> > & x );

  // Extract the raw transform
  std::vector<std::complex<double> > transform() const;
  // Extract the raw transform scaled by 1/n
  std::vector<std::complex<double> > scaled_transform() const;
  
  // Extract the ?
  std::vector<double> inverse() const;

  //
  // Misc helper function: average adjacent PSD values
  //

  void average_adjacent();

  

  std::vector<double> power_bands( const std::vector<double> & lwr , const std::vector<double> & upr ) 
    { 
      // Sum power over band intervals 
      if ( lwr.size() != upr.size() ) Helper::halt( "incorrectly specified bands()" ); 
      std::vector<double> pwr( lwr.size() ); 
      for (int j=0;j<cutoff;j++) 
	{ 
	  const double & f = frq[j]; 
	  for (int i=0;i<lwr.size();i++) 
	    if ( f >= lwr[i] && f < upr[i] ) pwr[i] += X[j];	 
	}     
      return pwr; 
    } 
    
};


//
// Helper 'bin' class
//

struct bin_t { 
  
  bin_t( double w , double mx_f , double Fs )
  : w(w) , mx_f(mx_f) , Fs(Fs) { } 

  double w, mx_f, Fs;
  
  int bin( const std::vector<double> & f , 
	   const std::vector<double> & y ) ;

  // data members
  std::vector<double> bspec, bfa, bfb; 
  
};



//
// Welch's power spectral density estimate
//

class PWELCH
{

  // Note: computes the power spectral density, not the power spectrum.   
  
  // Matlab equivalent:
  //     
  // [Pxx,F] = pwelch(X,WINDOW,NOVERLAP,NFFT,Fs) returns a PSD computed as
  // a function of physical frequency (Hz).  Fs is the sampling frequency
  // specified in Hz.  If Fs is empty, it defaults to 1 Hz.
  //
  // F is the vector of frequencies at which the PSD is estimated and has
  // units of Hz.  For real signals, F spans the interval [0,Fs/2] when NFFT
  // is even and [0,Fs/2) when NFFT is odd.  For complex signals, F always
  // spans the interval [0,Fs).
  
  
 public:
  
 PWELCH( const std::vector<double> & data , 
	 int Fs, 
	 double M , 
	 int noverlap_segments , 
	 window_function_t W = WINDOW_TUKEY50 , 
	 bool average_adj = false ) 
   : data(data) , Fs(Fs) , M(M) , noverlap_segments(noverlap_segments) , 
     window(W), average_adj(average_adj) 
  {

    // calculate implied overlap in actual data-points
    // the above specifies how many segments (of size 'M' seconds) we want in the 
    // window

    process(); 
  } 
  
  
  //
  // Derived variables
  //
    
  int N;

  std::vector<double> psd;
  
  std::vector<double> freq;
  
  double psdsum( double lwr , double upr )
  {
    // add is l <= x < y 
    double r = 0;
    for (int i=0;i<N;i++) 
      {
	if ( freq[i] >= upr ) break;
	if ( freq[i] >= lwr ) r += psd[i];
      }

    // area under the curve, multiple by bin width for each area
    double fbin = freq[1] - freq[0];

    return r * fbin;

  }
  
  double psdsum( frequency_band_t b )    
  {
    if ( globals::freq_band.find( b ) == globals::freq_band.end() ) 
      return 0;    
    freq_range_t f = globals::freq_band[ b ];
    return psdsum( f.first , f.second );
  }
 
  void psdsum( std::map<freq_range_t,double> * );

  void psdmean( std::map<freq_range_t,double> * );
  
 private:
    
  void process();
  
  // input signal  
  const std::vector<double> & data;
  
  // sampling rate (points per second)
  const int Fs; 
  
  // split segment into L intervals of length M (in seconds)
  double M;  
  
  // number of overlapping segments that should be placed in the epoch/window
  int noverlap_segments;  
  
  // window function (default Tukey 50% window)
  window_function_t window; 
  
  // option to average adjacent frequency bins (default=F)
  bool average_adj;
  
};


#endif
