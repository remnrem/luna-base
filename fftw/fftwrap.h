
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

extern logger_t logger;

struct edf_t;

std::map<double,double> fft_spectrum( const std::vector<double> * d , int Fs );

class FFT
{
  
  friend class coherence_t;

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


class coherence_t {

 public:

  coherence_t ( const std::vector<double> & x , 
		const std::vector<double> & y , 
		int Fs, 
		double segment_sec ,  // segment size in seconds 
		double overlap_sec , // overlap in seconds
		window_function_t W = WINDOW_HANN , 
		bool average_adj = false , 
		bool detrend = false , 
		bool zerocenter = false )		
    :  segment_sec(segment_sec) , overlap_sec(overlap_sec),
    Fs(Fs), window(W), detrend(detrend) , zerocenter(zerocenter) , average_adj(average_adj) ,  x(x), y(y) 
  {
    
    if ( x.size() != y.size() ) 
      Helper::halt( "coherence_t() called for signals of varying length" );

    total_points = x.size();

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
        
    process(); 
  } 

  double segment_sec, overlap_sec;


  
  int total_points;
  int noverlap_segments;
  int segment_points;
  int noverlap_points1, noverlap_points2;
  int segment_increment_points;

  int Fs;
  
  coh_t res;
  
 private:

  window_function_t window;
  
  bool detrend;
  bool zerocenter;

  bool average_adj;

  int N;

  const std::vector<double> & x;
  const std::vector<double> & y;  

  void process();
  
};



#endif
