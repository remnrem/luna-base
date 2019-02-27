
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


#include "fir.h"

#include "helper/helper.h"
#include "helper/logger.h"
#include "fftw/fftwrap.h"
#include "db/db.h"
#include "eval.h"

#include "edf/slice.h"
#include "edf/edf.h"

extern logger_t logger;

extern writer_t writer;


fir_impl_t::fir_impl_t( const std::vector<double> & coefs_ ) 
{
  count = 0;
  length = coefs_.size();
  coefs = coefs_;
  delayLine.resize( length );
  
  // expecting a linear-phase FIR with odd number of coefficients
  if ( coefs.size() % 2 != 1 ) Helper::halt( "expecting odd number of taps in FIR" );
  int del = ( coefs.size() - 1 ) / 2 ;
  
  double checksum = 0;
  for (int i=0;i<del;i++)
    checksum += fabs( coefs[i] - coefs[ coefs.size() - 1 - i ] );
  if ( checksum > 1e-8 ) Helper::halt( "problem in filter" );
  
}


void dsptools::design_fir( param_t & param )
{
  // assume this is called only by --fir option, and so we need to set this
  
  int fs = param.requires_int( "fs" );
  double ripple = param.requires_dbl( "ripple" );
  double tw = param.requires_dbl( "tw" );

  double f1 , f2;
  
  if ( param.has( "bandpass" ) )
    {
      std::vector<double> f = param.dblvector( "bandpass" );
      if ( f.size() != 2 ) Helper::halt( "expect bandpass=f1,f2" );
      f1 = f[0];
      f2 = f[1];
      logger << " designing bandpass filter, " << f1 << "-" << f2 << "Hz, ripple=" << ripple << ", tw=" << tw << ", fs=" << fs << "\n"; 
      design_bandpass_fir( ripple , tw , fs , f1, f2 , true );
    }
  
  if ( param.has( "bandstop" ) )
    {
      std::vector<double> f = param.dblvector( "bandstop" );
      if ( f.size() != 2 ) Helper::halt( "expect bandstop=f1,f2" );
      f1 = f[0];
      f2 = f[1];
      logger << " designing bandstop filter, " << f1 << "-" << f2 << "Hz, ripple=" << ripple << ", tw=" << tw << ", fs=" << fs << "\n"; 
      design_bandstop_fir( ripple , tw , fs , f1, f2 , true );
    }
  
  if ( param.has( "lowpass" ) )
    {
      f1 = param.requires_dbl( "lowpass" );
      logger << " designing lowpass filter, " << f1 << "Hz, ripple=" << ripple << ", tw=" << tw << ", fs=" << fs << "\n"; 
      design_lowpass_fir( ripple , tw , fs , f1 , true );
    }

  if ( param.has( "highpass" ) )
    {
      f1 = param.requires_dbl( "highpass" );
      logger << " designing highpass filter, " << f1 << "Hz, ripple=" << ripple << ", tw=" << tw << ", fs=" << fs << "\n"; 
      design_highpass_fir( ripple , tw , fs , f1 , true );
    }

}
  

std::vector<double> dsptools::design_bandpass_fir( double ripple , double tw , double fs , double f1 , double f2 , bool eval )
{

  fir_t fir;

  int kaiserWindowLength;
  
  double beta;
  fir.calculateKaiserParams( ripple , tw , fs , &kaiserWindowLength, &beta);

  if ( kaiserWindowLength % 2 == 0 ) ++kaiserWindowLength;

  std::vector<double> fc = fir.create2TransSinc( kaiserWindowLength, f1 , f2, fs , fir_t::BAND_PASS );

  fc = fir.createKaiserWindow(&fc, beta);

  if ( eval )
    {
      std::string label = "BANDPASS_" 
	+ Helper::dbl2str( f1 ) + ".." + Helper::dbl2str( f2 ) + "_" + Helper::dbl2str( ripple ) + "_" + Helper::dbl2str( tw ) ;
      
      fir.outputFFT( label , fc , fs );
    }

  return fc;
}


std::vector<double> dsptools::design_bandstop_fir( double ripple , double tw , double fs , double f1 , double f2 , bool eval )
{
  fir_t fir;

  int kaiserWindowLength;
  double beta;
  

  fir.calculateKaiserParams( ripple , tw , fs , &kaiserWindowLength, &beta);

  if ( kaiserWindowLength % 2 == 0 ) ++kaiserWindowLength;

  std::vector<double> fc = fir.create2TransSinc( kaiserWindowLength, f1 , f2, fs , fir_t::BAND_STOP );
 
  fc = fir.createKaiserWindow(&fc, beta);

  if ( eval ) 
    {
      std::string label = "BANDSTOP_" 
	+ Helper::dbl2str( f1 ) + ".." + Helper::dbl2str( f2 ) + "_" + Helper::dbl2str( ripple ) + "_" + Helper::dbl2str( tw ) ;
      
      fir.outputFFT( label , fc , fs );
    }

  return fc;

}

std::vector<double> dsptools::design_lowpass_fir( double ripple  , double tw , double fs , double f , bool eval )
{
  fir_t fir;

  int kaiserWindowLength;
  double beta;

  fir.calculateKaiserParams( ripple , tw , fs , &kaiserWindowLength, &beta);
  if ( kaiserWindowLength % 2 == 0 ) ++kaiserWindowLength;

  std::vector<double> fc = fir.create1TransSinc( kaiserWindowLength, f , fs , fir_t::LOW_PASS );

  fc = fir.createKaiserWindow(&fc, beta);

  if ( eval ) 
    {
      std::string label = "LOWPASS_" 
	+ Helper::dbl2str( f ) + "_" + Helper::dbl2str( ripple ) + "_" + Helper::dbl2str( tw ) ;
      
      fir.outputFFT( label , fc , fs );
    }
  
  return fc;

}

std::vector<double> dsptools::design_highpass_fir( double ripple , double tw , double fs , double f , bool eval )
{
  fir_t fir;

  int kaiserWindowLength;
  double beta;

  fir.calculateKaiserParams( ripple , tw , fs , &kaiserWindowLength, &beta);
  if ( kaiserWindowLength % 2 == 0 ) ++kaiserWindowLength;

  std::vector<double> fc = fir.create1TransSinc( kaiserWindowLength, f , fs , fir_t::HIGH_PASS );

  fc = fir.createKaiserWindow(&fc, beta);
  
  if ( eval )
    {
      std::string label = "HIGHPASS_" 
	+ Helper::dbl2str( f ) + "_" + Helper::dbl2str( ripple ) + "_" + Helper::dbl2str( tw ) ;
      
      fir.outputFFT( label , fc , fs );
    }

  return fc;

}
    

//
// apply FIR
//

std::vector<double> dsptools::apply_fir( const std::vector<double> & x , int fs, fir_t::filterType ftype , double ripple , double tw , double f1, double f2 )
{

  std::vector<double> fc;
  
  if ( ftype == fir_t::BAND_PASS ) 
    fc = design_bandpass_fir( ripple , tw , fs , f1, f2 );    
  else if ( ftype == fir_t::BAND_STOP )
    fc = design_bandstop_fir( ripple , tw , fs , f1, f2 );
  else if ( ftype == fir_t::LOW_PASS )
    fc = design_lowpass_fir( ripple , tw , fs , f1 );
  else if ( ftype == fir_t::HIGH_PASS )
    fc = design_highpass_fir( ripple , tw , fs , f1 );
  
  //
  // Apply FIR 
  //
  
  fir_impl_t fir_impl ( fc );
  
  return fir_impl.filter( &x );  

}

void dsptools::apply_fir( edf_t & edf , param_t & param )
{

  double ripple = param.requires_dbl( "ripple" );
  
  double tw = param.requires_dbl( "tw" );
  
  double f1 , f2 ;
  
  fir_t::filterType ftype = fir_t::BAND_PASS;
  
  if ( param.has( "bandpass" ) )
    {
      ftype = fir_t::BAND_PASS;
      std::vector<double> f = param.dblvector( "bandpass" );
      if ( f.size() != 2 ) Helper::halt( "expecting bandpass=f1,f2" );
      f1 = f[0];
      f2 = f[1];
    }
  else if ( param.has( "bandstop" ) )
    {
      ftype = fir_t::BAND_STOP;
      std::vector<double> f = param.dblvector( "bandstop" );
      if ( f.size() != 2 ) Helper::halt( "expecting bandstop=f1,f2" );
      f1 = f[0];
      f2 = f[1];
    }
  else if ( param.has( "lowpass" ) )
    {
      ftype = fir_t::LOW_PASS;
      f1 = param.requires_dbl( "lowpass" );
    }
  else if ( param.has( "highpass" ) )
    {
      ftype = fir_t::HIGH_PASS;
      f1 = param.requires_dbl( "highpass" );      
    }
  else 
    Helper::halt( "need to specify FIR type as bandpass, bandstop, lowpass or highpass" );

  //
  // Signals
  //

  std::string signal_label = param.requires( "sig" );
  
  signal_list_t signals = edf.header.signal_list( signal_label );      

  std::vector<double> Fs = edf.header.sampling_freq( signals );
  
  const int ns = signals.size();

  //
  // Process each signal
  //

  for (int s=0; s<ns; s++)
    {

      //
      // skip annotation channels
      //
      
      if ( edf.header.is_annotation_channel(s) ) continue;
      
      apply_fir( edf , signals(s) , ftype , ripple, tw , f1 , f2 );
      
    }

}

void dsptools::apply_fir( edf_t & edf , int s , fir_t::filterType ftype , double ripple , double tw , double f1, double f2 )
{
      
  
  interval_t interval = edf.timeline.wholetrace();
  
  //
  // Filter data
  //
  
  std::cerr << " filtering channel " << edf.header.label[ s ] << ", ";

  //
  // Pull entire signals out
  //
  
  slice_t slice( edf , s , interval );
  
  const std::vector<double> * d = slice.pdata();
      
  int fs = edf.header.sampling_freq( s );
      
  //
  // Design FIR
  // 

  std::vector<double> fc;
      
  if ( ftype == fir_t::BAND_PASS ) 
    {
      fc = design_bandpass_fir( ripple , tw , fs , f1, f2 );    
      std::cerr << "bandpass FIR order " << fc.size() << "\n";
    }
  else if ( ftype == fir_t::BAND_STOP )
    {
      fc = design_bandstop_fir( ripple , tw , fs , f1, f2 );
      std::cerr << "bandstop FIR order " << fc.size() << "\n";
    }
  else if ( ftype == fir_t::LOW_PASS )
    {
      fc = design_lowpass_fir( ripple , tw , fs , f1 );
      std::cerr << "lowpass FIR order " << fc.size() << "\n";
    }
  else if ( ftype == fir_t::HIGH_PASS )
    {
      fc = design_highpass_fir( ripple , tw , fs , f1 );
      std::cerr << "highpass FIR order " << fc.size() << "\n";
    }
  
  //int filter_order = num_taps == -1 ? 3 * ( Fs[s] / lwr ) : num_taps ; 
  
  //
  // Apply FIR 
  //
  
  fir_impl_t fir_impl ( fc );
  
  std::vector<double> filtered = fir_impl.filter( d );
  
  //
  // Place back
  //
  
  edf.update_signal( s , &filtered );
  
}


void fir_t::demo()
{
  
  // int windowLength = 201;
  // double sampFreq = 44100;
  
  int windowLength = 201;
  double sampFreq = 200;

  // Low and high pass filters
  // double transFreq = 10000;

  // std::vector<double> lpf = create1TransSinc(windowLength, transFreq, sampFreq, LOW_PASS);
  // std::vector<double> lpf_hamming = createWindow(&lpf, HAMMING);
  // std::vector<double> lpf_blackman = createWindow(&lpf, BLACKMAN);
  
  // std::vector<double> hpf = create1TransSinc(windowLength, transFreq, sampFreq, HIGH_PASS);
  // std::vector<double> hpf_hamming = createWindow(&hpf, HAMMING);

  // outputFFT("lpf-hamming.dat", lpf_hamming, sampFreq);
  // outputFFT("lpf-blackman.dat", lpf_blackman, sampFreq);
  // outputFFT("hpf-hamming.dat", hpf_hamming, sampFreq);
  
  // Band pass and band stop filters
  // double trans1Freq = 5000;
  // double trans2Freq = 17050;
  
  double trans1Freq = 0.3;
  double trans2Freq = 30;

  std::vector<double> bpf = create2TransSinc(windowLength, trans1Freq, trans2Freq, sampFreq, BAND_PASS);
  //  std::vector<double> bpf_hamming = createWindow(&bpf, HAMMING);
  std::vector<double> bpf_hamming = createWindow(&bpf, RECTANGULAR);
  
  // std::vector<double> bsf = create2TransSinc(windowLength, trans1Freq, trans2Freq, sampFreq, BAND_STOP);
  // std::vector<double> bsf_hamming = createWindow(&bsf, HAMMING);
  
  outputFFT("bpf-hamming.dat", bpf_hamming, sampFreq);
  //  outputFFT("bsf-hamming.dat", bsf_hamming, sampFreq);

  return;

  // Kaiser Window
  // int kaiserWindowLength;
  // double beta;
  
  // calculateKaiserParams(0.001, 800, sampFreq, &kaiserWindowLength, &beta);
  
  // lpf = create1TransSinc(kaiserWindowLength, transFreq, sampFreq, LOW_PASS);
  // std::vector<double> lpf_kaiser = createKaiserWindow(&lpf, beta);
  
  // outputFFT("lpf-kaiser.dat", lpf_kaiser, sampFreq);

}


// Create sinc function for filter with 1 transition - Low and High pass filters
std::vector<double> fir_t::create1TransSinc(int windowLength, double transFreq, double sampFreq, enum filterType type)
{
  
  // Allocate memory for the window
  std::vector<double> window( windowLength );
  
  if (type != LOW_PASS && type != HIGH_PASS) 
    Helper::halt("create1TransSinc: Bad filter type, should be either LOW_PASS of HIGH_PASS");
  
  // Calculate the normalised transistion frequency. As transFreq should be
  // less than or equal to sampFreq / 2, ft should be less than 0.5
  double ft = transFreq / sampFreq;
  
  double m_2 = 0.5 * (windowLength-1);
  int halfLength = windowLength / 2;
  
  // Set centre tap, if present
  // This avoids a divide by zero
  if (2*halfLength != windowLength) {
    double val = 2.0 * ft;
    
    // If we want a high pass filter, subtract sinc function from a dirac pulse
    if (type == HIGH_PASS) val = 1.0 - val;
    
    window[halfLength] = val;
  }
  else if (type == HIGH_PASS) 
    Helper::halt("create1TransSinc: For high pass filter, window length must be odd");
   
  // This has the effect of inverting all weight values
  if (type == HIGH_PASS) ft = -ft;
  
  // Calculate taps
  // Due to symmetry, only need to calculate half the window
  for (int n=0 ; n<halfLength ; n++) 
    {
      double val = sin(2.0 * M_PI * ft * (n-m_2)) / (M_PI * (n-m_2));      
      window[n] = val;
      window[windowLength-n-1] = val;
    }
  
  return window;
}

// Create two sinc functions for filter with 2 transitions - Band pass and band stop filters
std::vector<double> fir_t::create2TransSinc(int windowLength, double trans1Freq, double trans2Freq, double sampFreq, enum filterType type)
{

  // Allocate memory for the window
  std::vector<double> window( windowLength );
  
  if (type != BAND_PASS && type != BAND_STOP) 
    Helper::halt("create2TransSinc: Bad filter type, should be either BAND_PASS or BAND_STOP");
  
  // Calculate the normalised transistion frequencies.
  double ft1 = trans1Freq / sampFreq;
  double ft2 = trans2Freq / sampFreq;
  
  double m_2 = 0.5 * (windowLength-1);
  int halfLength = windowLength / 2;
  
  // Set centre tap, if present
  // This avoids a divide by zero
  if (2*halfLength != windowLength) {
    double val = 2.0 * (ft2 - ft1);
    
    // If we want a band stop filter, subtract sinc functions from a dirac pulse
    if (type == BAND_STOP) val = 1.0 - val;
    
    window[halfLength] = val;
  }
  else 
    Helper::halt("create1TransSinc: For band pass and band stop filters, window length must be odd");
  
  // Swap transition points if Band Stop
  if (type == BAND_STOP) {
    double tmp = ft1;
    ft1 = ft2; ft2 = tmp;
  }
  
  // Calculate taps
  // Due to symmetry, only need to calculate half the window
  for (int n=0 ; n<halfLength ; n++) {
    double val1 = sin(2.0 * M_PI * ft1 * (n-m_2)) / (M_PI * (n-m_2));
    double val2 = sin(2.0 * M_PI * ft2 * (n-m_2)) / (M_PI * (n-m_2));
    
    window[n] = val2 - val1;
    window[windowLength-n-1] = val2 - val1;
  }
  
  return window;
}

// Create a set of window weights
// in - If not null, each value will be multiplied with the window weight
// out - The output weight values, if NULL and new array will be allocated
// windowLength - the number of weights
// windowType - The window type

std::vector<double> fir_t::createWindow( const std::vector<double> * in, enum windowType type)
{
  
  int windowLength = in->size();

  std::vector<double> out( windowLength , 0 );
  
  int m = windowLength - 1;
  int halfLength = windowLength / 2;
  
  // Calculate taps
  // Due to symmetry, only need to calculate half the window
  switch (type)
    {
    case RECTANGULAR:
      for (int n=0 ; n<windowLength ; n++) {
	out[n] = 1.0;
      }
      break;
      
    case BARTLETT:
      for (int n=0 ; n<=halfLength ; n++) {
	double tmp = (double) n - (double)m / 2;
	double val = 1.0 - (2.0 * fabs(tmp))/m;
	out[n] = val;
	out[windowLength-n-1] = val;
      }
      
      break;
      
    case HANN:
      for (int n=0 ; n<=halfLength ; n++) {
	double val = 0.5 - 0.5 * cos(2.0 * M_PI * n / m);
	out[n] = val;
	out[windowLength-n-1] = val;
      }
      
      break;
      
    case HAMMING:
      for (int n=0 ; n<=halfLength ; n++) {
	double val = 0.54 - 0.46 * cos(2.0 * M_PI * n / m);
	out[n] = val;
	out[windowLength-n-1] = val;
      }
      break;
      
    case BLACKMAN:
      for (int n=0 ; n<=halfLength ; n++) {
	double val = 0.42 - 0.5 * cos(2.0 * M_PI * n / m) + 0.08 * cos(4.0 * M_PI * n / m);
	out[n] = val;
	out[windowLength-n-1] = val;
      }
      break;
    }
  

  // If input has been given, multiply with out
  if ( in != NULL ) 
    for (int n=0 ; n<windowLength ; n++) 
      out[n] *= (*in)[n];
  
  return out;
}


// Transition Width (transWidth) is given in Hz
// Sampling Frequency (sampFreq) is given in Hz
// Window Length (windowLength) will be set
void fir_t::calculateKaiserParams(double ripple, double transWidth, double sampFreq, int *windowLength, double *beta)
{
  // Calculate delta w
  double dw = 2 * M_PI * transWidth / sampFreq;
  
  // Calculate ripple dB
  double a = -20.0 * log10(ripple);
  
  // Calculate filter order
  int m;
  if (a>21) m = ceil((a-7.95) / (2.285*dw));
  else m = ceil(5.79/dw);
  
  *windowLength = m + 1;
  
  if (a<=21) *beta = 0.0;
  else if (a<=50) *beta = 0.5842 * pow(a-21, 0.4) + 0.07886 * (a-21);
  else *beta = 0.1102 * (a-8.7);
}


std::vector<double> fir_t::createKaiserWindow( const std::vector<double> * in, double beta )
{

  const int windowLength = in->size();
  std::vector<double> out( windowLength );
  
  double m_2 = (double)(windowLength-1) / 2.0;
  double denom = modZeroBessel(beta);  // Denominator of Kaiser function
    
  for (int n=0 ; n<windowLength ; n++)
    {
      double val = ((n) - m_2) / m_2;
      val = 1 - (val * val);
      out[n] = modZeroBessel(beta * sqrt(val)) / denom;
    }
  
  // If input has been given, multiply with out
  if ( in != NULL ) 
    for (int n=0 ; n<windowLength ; n++) 
      out[n] *= (*in)[n];
  
  return out;
}

double fir_t::modZeroBessel(double x)
{

  double x_2 = x/2;
  double num = 1;
  double fact = 1;
  double result = 1;
  
  for (int i=1 ; i<20 ; i++) {
    num *= x_2 * x_2;
    fact *= i;
    result += num / (fact * fact);  
  }  
  return result;
}


int fir_t::outputFFT(const std::string & label, const std::vector<double> & window, double sampFreq)
{	
  

  writer.level( label , "FIR" );

  
  
  //
  // Filter coefficients
  //
  
  writer.numeric_factor( "TAP" );
  
  for (int i=0;i<window.size();i++)
    {
      writer.level( i , "TAP" );
      writer.value( "W" , window[i] );
    }
  writer.unlevel( "TAP" );
  
  //
  // Impulse response
  //

  // 2 window second around filter size
  double sz = window.size() / (double)sampFreq + 2 ;
  fir_impl_t fir_impl ( window );
  std::vector<double> xx0( sampFreq * sz );
  xx0[sampFreq*(sz/2.0)-1] = 1;
  std::vector<double> xx = fir_impl.filter( &xx0 ); 
  double idx0 = sampFreq*(sz/2.0)-1;
  
  writer.numeric_factor( "SEC" );  
  
  for (int xi=0;xi<xx.size();xi++)
    {
      double tp = 1.0/sampFreq * (xi - idx0 ); 
      writer.level( Helper::dbl2str(tp) , "SEC" );
      writer.value( "IR" , xx[xi] ); 
    }
  writer.unlevel( "SEC" );
  
  //
  // Frequency response
  //
  
  const int windowLength = window.size();


  writer.value( "FS" , sampFreq );
  writer.value( "NTAPS" , windowLength );


  double *in;
  fftw_complex *out;
  fftw_plan plan;
  int result = 0;
  
  // If the window length is short, zero padding will be used
  int fftSize = (windowLength < 2048) ? 2048 : windowLength;
  
  // Calculate size of result data
  int resultSize = (fftSize / 2) + 1;
  
  // Allocate memory to hold input and output data
  in = (double *) fftw_malloc(fftSize * sizeof(double));
  out = (fftw_complex *) fftw_malloc(resultSize * sizeof(fftw_complex));
  if (in == NULL || out == NULL) {
    result = 1;
    Helper::halt( "fir_t: could not allocate input/output data");
    goto finalise;
  }

  // Create the plan and check for success
  plan = fftw_plan_dft_r2c_1d(fftSize, in, out, FFTW_MEASURE); 
  if (plan == NULL) {
    result = 1;
    Helper::halt( "fir_t: could not create plan" );
    goto finalise;
  }
  
  // Copy window and add zero padding (if required)
  int i;
  for (i=0 ; i<windowLength ; i++) in[i] = window[i];
  for ( ; i<fftSize ; i++) in[i] = 0;
  
  // Perform fft
  fftw_execute(plan);
  
  // Output result
  for (i=0 ; i<resultSize ; i++)
    {
      double freq = sampFreq * i / fftSize;
      double mag = sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]);
      double magdB = 20 * log10(mag);
      double phase = atan2(out[i][1], out[i][0]);

      writer.level( freq , globals::freq_strat );
      writer.value( "MAG" , mag );
      writer.value( "MAG_DB" , magdB );
      writer.value( "PHASE" , phase );
      
    }
  
  writer.unlevel( globals::freq_strat );

  // Perform any cleaning up
 finalise:
  if (plan != NULL) fftw_destroy_plan(plan);
  if (in != NULL) fftw_free(in);
  if (out != NULL) fftw_free(out);
  
  return result;

  writer.unlevel( "FIR" );

}



std::vector<double> fir_impl_t::filter( const std::vector<double> * x ) 
{
  
  if ( length % 2 == 0 ) Helper::halt("fir_impl_t requries odd # of coeffs");
  
  const int n = x->size();
  
  const int delay_idx = (length-1)/2;
  
  std::vector<double> r( n ) ;
  
  const double * p = &(*x)[0];
  
  // burn in
  for (int i=0;i<delay_idx;i++) 
    getOutputSample( *(p++) );
  
  // process
  int j = 0;
  for (int i=delay_idx;i<n;i++) 
    r[j++] = getOutputSample( *(p++) );
  
  // zero-pad end of signal
  for (int i=0;i<delay_idx;i++) 
    r[j++] = getOutputSample( 0 );
  
  return r;
  
}



std::vector<double> fir_impl_t::fft_filter( const std::vector<double> * px )
{
  
  std::vector<double> x = *px;
  std::vector<double> h = coefs;
  
  // signal length
  const int M = x.size();
  
  // filter length
  const int L = h.size();
  
  // next power of 2 greater than M+L-1
  long int Nfft = MiscMath::nextpow2( M + L - 1 );
  
  // zero-padding
  x.resize( Nfft , 0 );
  h.resize( Nfft , 0 );

  // FFT
  FFT fftx( Nfft , 1 , FFT_FORWARD );
  fftx.apply( x );
  std::vector<dcomp> rfftx = fftx.transform();
  
  FFT ffth( Nfft , 1 , FFT_FORWARD );
  ffth.apply( h );
  std::vector<dcomp> rffth = ffth.transform();
  
  // convolution in the frequency domain

  std::vector<dcomp> y( Nfft );
  for (int i=0;i<rfftx.size();i++) y[i] = rfftx[i] * rffth[i]; 
  
  // inverse FFT

  FFT ifft( Nfft , 1 , FFT_INVERSE );
  ifft.apply( y );
  std::vector<dcomp> conv_tmp = ifft.transform();
  
  dcomp denom( 1.0/(double)Nfft , 0 );
  
  // Normalize
  for (int i=0;i<Nfft;i++) conv_tmp[i] *= denom;

  // Trim and return read component
  std::vector<double> conv;
  const int delay_idx = (length-1)/2;
  int j = delay_idx;
  for (int i=0;i<M;i++)
    conv.push_back( std::real( conv_tmp[j++] ) );

  return conv;
  
}
 
