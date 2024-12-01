
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

void psd_shape_metrics( const std::vector<double> & f , // frq
			const std::vector<double> & p , // log(power)
			const int w ,
			double * m1 ,
			double * m2 ,
			std::vector<double> * detrended = NULL ,
			std::vector<double> * smoothed = NULL , 
			std::vector<double> * difference = NULL );


bool spectral_slope_helper( const std::vector<double> & psd , 
			    const std::vector<double> & f , 
			    const std::vector<double> & fr , 
			    const double outlier , 
			    const bool display = true , 
			    double * beta = NULL , 
			    double * betan = NULL , 
			    double * intercept = NULL , 
			    double * rsq = NULL );

void peakedness( const std::vector<double> & p , 
		 const std::vector<double> & f0 , 
		 const int peak_median_filter_n , 
		 const std::vector<double> & pr , 
		 const bool verbose );



//
// Complex FFT
//

class FFT
{
  
  friend struct coherence_t;
  friend struct precoh_t;

 public:

  FFT() { } 

  FFT( int Ndata , int Nfft , int Fs , fft_t type = FFT_FORWARD , window_function_t window = WINDOW_NONE ) 
    {
      init( Ndata , Nfft , Fs , type , window );
    }

  void init( int Ndata , int Nfft , int Fs , fft_t type = FFT_FORWARD , window_function_t window = WINDOW_NONE );

  void reset();
  
  ~FFT();
  
 private:

  // Size of data 
  int Ndata;
  
  // Sampling rate, so we can construct the appropriate Hz for the PSD
  int Fs;

  // Forward or inverse FFT?
  fft_t type;
  
  // Optional windowing function
  window_function_t window;
  std::vector<double> w;

  // Input signal
  fftw_complex *in;

  // Output signal
  fftw_complex *out;
  
  // FFT plan from FFTW3
  fftw_plan p;

  // Size (NFFT)
  int Nfft;
  
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
  
  std::vector<double> inverse() const;
  std::vector<double> unscaled_inverse() const;

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
// Real 1D DFT
//

class real_FFT
{
  
  friend struct coherence_t;
  friend struct precoh_t;
  
 public:
  
  real_FFT() { } 

  real_FFT( int Ndata , int Nfft , int Fs , window_function_t window = WINDOW_NONE ) 
    {
      init( Ndata , Nfft , Fs , window );
    }
  
  void init( int Ndata , int Nfft , int Fs , window_function_t window = WINDOW_NONE );

  void norm_fac( const double f ) { normalisation_factor = f; } 

  void reset() ;
  
  ~real_FFT();
  
 private:

  // Size of data 
  int Ndata;
  
  // Sampling rate, so we can construct the appropriate Hz for the PSD
  int Fs;

  // Forward or inverse FFT?
  fft_t type;
  
  // Optional windowing function
  window_function_t window;
  std::vector<double> w;

  // Input signal (real)
  double * in;

  // Output signal
  fftw_complex *out;
  
  // FFT plan from FFTW3
  fftw_plan p;
  
  // Size (NFFT)
  int Nfft;
  
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
  
  bool apply( const std::vector<double> & x );
  bool apply( const double * x , const int n );
    
  // Extract the raw transform
  std::vector<std::complex<double> > transform() const;

  // Extract the raw transform scaled by 1/n
  std::vector<std::complex<double> > scaled_transform() const;
  
  std::vector<double> inverse() const;

  // Misc helper function: average adjacent PSD values

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
// 1D C->R inverse DFT
//

//
// Real 1D DFT
//

class real_iFFT
{
  
 public:
  
  real_iFFT() { } 

  real_iFFT( int Ndata , int Nfft , int Fs , window_function_t window = WINDOW_NONE ) 
  {
    init( Ndata , Nfft , Fs , window );
  }
  
  void init( int Ndata , int Nfft , int Fs , window_function_t window = WINDOW_NONE );
  
  void reset();
  
  ~real_iFFT();

 private:
  
  // Size of data 
  int Ndata;
  
  // Sampling rate, so we can construct the appropriate Hz for the PSD
  int Fs;
  
  // Optional windowing function
  window_function_t window;
  std::vector<double> w;

  // Input signal (complex)
  fftw_complex *in;

  // Output signal (real)
  double *out;
  
  // FFT plan from FFTW3
  fftw_plan p;
  
  // Size (NFFT)
  int Nfft;
  
  // Normalisation factor given the window
  double normalisation_factor;

 public:
  
  int cutoff;
  std::vector<double> X;
  std::vector<double> mag;
  std::vector<double> frq;
  
 public:
  
  bool apply( const std::vector<std::complex<double> > & x );   
  std::vector<double> inverse() const;
  std::vector<double> unscaled_inverse() const;
    
};




//
// Helper 'bin' class
//

struct bin_t { 
  
  // bin_t( double w , double mx_f , double Fs )
  //   : w(w) , mn_f(0) , mx_f(mx_f) , Fs(Fs) { } 
  
  bin_t( double mn_f , double mx_f , int fac )
    : fac(fac) , mn_f(mn_f) , mx_f(mx_f) { } 
  
  double fac, mn_f , mx_f;
  
  int bin( const std::vector<double> & f , 
	   const std::vector<double> & y ) ;
  
  // data members
  std::vector<double> bspec, bfa, bfb;
  std::vector<std::string> nominal; 
  
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
	 bool use_median = false ,
	 bool calc_seg_sd = false , 
	 bool average_adj = false ,
	 bool use_nextpow2 = false ,
	 bool do_normalization = true ) 
   : data(data) , Fs(Fs) , M(M) , noverlap_segments(noverlap_segments) , 
     window(W),
     use_median(use_median), calc_seg_sd(calc_seg_sd),
     average_adj(average_adj) , use_nextpow2( use_nextpow2)  , do_normalization(do_normalization) 
  {

    // calculate implied overlap in actual data-points the above
    // specifies how many segments (of size 'M' seconds) we want in
    // the window

    process(); 
  } 
  
  
  //
  // Derived variables
  //
    
  int N;

  std::vector<double> psd;

  // optionally, SD of segments
  std::vector<double> psdsd; 
  
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

    // area under the curve, scale by bin width for each area
    double fbin = freq[1] - freq[0];

    return r * fbin;

  }

  double psdsdsum( double lwr , double upr )
  {
    // add is l <= x < y 
    double r = 0;
    for (int i=0;i<N;i++) 
      {
	if ( freq[i] >= upr ) break;
	if ( freq[i] >= lwr ) r += psdsd[i];
      }

    // area under the curve, scale by bin width for each area
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
 
  double psdsdsum( frequency_band_t b )    
  {
    if ( globals::freq_band.find( b ) == globals::freq_band.end() ) 
      return 0;    
    freq_range_t f = globals::freq_band[ b ];
    return psdsdsum( f.first , f.second );
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

  // median-averaging instead of mean
  bool use_median;

  // get SD of segments
  bool calc_seg_sd;
  
  // option to average adjacent frequency bins (default=F)
  bool average_adj;

  // always set NFFT to the next power of 2
  bool use_nextpow2;

  // do 1/N^2 norm (incl. window)
  bool do_normalization;

};


#endif
