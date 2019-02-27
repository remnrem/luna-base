
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

#ifndef __LUNA_DSP_FIR_H__
#define __LUNA_DSP_FIR_H__

#include <vector>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <complex>

struct param_t;

struct edf_t;

// https://ptolemy.eecs.berkeley.edu/eecs20/week12/implementation.html

struct fir_impl_t { 
  
  int length;
  std::vector<double> delayLine;
  std::vector<double> coefs;
  int count;
  
  fir_impl_t( const std::vector<double> & coefs_ ); 
  
  std::vector<double> filter( const std::vector<double> * x );

  std::vector<double> fft_filter( const std::vector<double> * x );
  
  double getOutputSample(double inputSample) 
  {
    
    delayLine[count] = inputSample;
    
    double result = 0.0;
    
    int index = count;
    
    for (int i=0; i<length; i++) 
      {
	result += coefs[i] * delayLine[index--];
	if (index < 0) index = length-1;
      }
    
    if (++count >= length) count = 0;
    
    return result;
  }
  
};



struct fir_t
{

  // class for FIR design, using implementations from "FIR filters by
  // Windowing" by A.Greensted - Feb 2010
  // http://www.labbookpages.co.uk
  // http://www.labbookpages.co.uk/audio/firWindowing.html

  enum filterType { LOW_PASS, HIGH_PASS, BAND_PASS, BAND_STOP };
  enum windowType { RECTANGULAR, BARTLETT, HANN, HAMMING, BLACKMAN };
  
  // Prototypes
  std::vector<double> create1TransSinc( int windowLength, double transFreq, double sampFreq, enum filterType type);
  std::vector<double> create2TransSinc( int windowLength, double trans1Freq, double trans2Freq, double sampFreq, enum filterType type);
  
  std::vector<double> createWindow( int windowLength, enum windowType type)
  { std::vector<double> d( windowLength , 1 ); return createWindow( &d, type ); }

  std::vector<double> createWindow( const std::vector<double> * in, enum windowType type);

  void calculateKaiserParams(double ripple, double transWidth, double sampFreq, int *windowLength, double *beta);

  std::vector<double> createKaiserWindow( const int windowLength, double beta)
  { std::vector<double> d( windowLength, 1.0); return createKaiserWindow( &d, beta ); }

  std::vector<double> createKaiserWindow( const std::vector<double> *in, double beta);

  double modZeroBessel(double x);

  int outputFFT( const std::string & , const std::vector<double> & , double sampFreq);

  void demo();

};


namespace dsptools 
{ 
  
  //
  // using Kaiser Window
  //
  
  void design_fir( param_t & param );
  
  std::vector<double> design_bandpass_fir( double ripple , double tw , double fs , double f1 , double f2 , bool eval = false );
  std::vector<double> design_bandstop_fir( double ripple , double tw , double fs , double f1 , double f2 , bool eval = false );
  std::vector<double> design_lowpass_fir( double ripple  , double tw , double fs , double f , bool eval = false );
  std::vector<double> design_highpass_fir( double ripple , double tw , double fs , double f , bool eval = false );
    
  //
  // apply FIR
  //
  
  void apply_fir( edf_t & edf , param_t & param );
  void apply_fir( edf_t & edf , int s , fir_t::filterType , double ripple , double tw , double f1, double f2 );
  std::vector<double> apply_fir( const std::vector<double> & , int fs , fir_t::filterType ftype , double ripple , double tw , double f1, double f2 );
  
}


#endif

