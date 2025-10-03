
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
#include "param.h"

#include "edf/slice.h"
#include "edf/edf.h"
#include "dsp/ngaus.h"
#include "dsp/conv.h"

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

  bool use_kaiser = param.has( "tw" ) && param.has( "ripple" );

  // ntaps = order + 1
  // order should be even, ntaps should be odd
  bool fixed_order = param.has( "order" );
  bool from_file = param.has( "file" );


  //
  // read filter coefficients from a file, and analyse
  //
  
  if ( from_file )
    {
      std::vector<double> fc;
      std::string fir_file = param.value( "file" );
      if ( ! Helper::fileExists( fir_file ) ) Helper::halt( "could not find " + fir_file );
      std::ifstream IN1( fir_file.c_str() , std::ios::in );
      while ( ! IN1.eof() )
	{
	  double c;
	  IN1 >> c;
	  if ( IN1.eof() ) break;
	  fc.push_back( c );
	}
      IN1.close();
      std::string label = fir_file;
      fir_t fir;
      fir.outputFFT( label , fc , fs );
      return;
    }


  //
  // Get key parameters
  //
  
  if ( use_kaiser == fixed_order )
    Helper::halt( "must specify either Kaiser window format or fixed FIR order" );

  // allow separate high+low pass filters
  std::vector<double> ripple, tw;
  if ( use_kaiser )
    {
      ripple = param.dblvector( "ripple" ) ;
      tw = param.dblvector( "tw" );
    }
  bool kaiser_high_and_low = ripple.size() == 2 && tw.size() == 2;

  
  
  int order = fixed_order ? param.requires_int( "order" ) : 0 ;  
  
  double f1 , f2;

  //fir_t::windowType window = fir_t::HAMMING;
  fir_t::windowType window = fir_t::HAMMING;
  if ( param.has( "rectangular" ) ) window = fir_t::RECTANGULAR;
  else if ( param.has( "bartlett" ) ) window = fir_t::BARTLETT;
  else if ( param.has( "hann" ) ) window = fir_t::HANN;
  else if ( param.has( "blackman" ) ) window = fir_t::BLACKMAN;

  if ( param.has( "bandpass" ) )
    {
      std::vector<double> f = param.dblvector( "bandpass" );
      if ( f.size() != 2 ) Helper::halt( "expect bandpass=f1,f2" );
      f1 = f[0];
      f2 = f[1];

      if ( use_kaiser )
	{

	  if ( kaiser_high_and_low )
	    {
	      logger << " designing bandpass filter, convolving highpass (" << f1 << "Hz, ripple="<< ripple[0]<<", tw="<< tw[0]<<")"
		     << " with lowpass ("<< f2 << "Hz, ripple="<< ripple[1]<<", tw="<< tw[1]<<")\n";
	      design_bandpass_fir( ripple[0], ripple[1] , tw[0], tw[1] , fs , f1, f2 , true );
	    }
	  else
	    {
	      logger << " designing bandpass filter, " << f1 << "-" << f2 << "Hz, ripple=" << ripple[0] << ", tw=" << tw[0] << ", fs=" << fs << "\n"; 
	      design_bandpass_fir( ripple[0] , tw[0] , fs , f1, f2 , true );
	    }
	}
      else
	{
	  logger << " designing bandpass filter, " << f1 << "-" << f2 << "Hz, order=" << order << ", fs=" << fs
		 << " with a " << fir_t::window( window ) << " window\n";
	  design_bandpass_fir( order , fs , f1, f2 , window , true ); 
	}

      return;

    }
  
      
  if ( param.has( "bandstop" ) )
    {
      std::vector<double> f = param.dblvector( "bandstop" );
      if ( f.size() != 2 ) Helper::halt( "expect bandstop=f1,f2" );
      f1 = f[0];
      f2 = f[1];

      if ( use_kaiser )
	{
	  logger << " designing bandstop filter, " << f1 << "-" << f2 << "Hz, ripple=" << ripple[0] << ", tw=" << tw[0] << ", fs=" << fs << "\n"; 
	  design_bandstop_fir( ripple[0] , tw[0] , fs , f1, f2 , true );
	}
      else
	{
	  logger << " designing bandstop filter, " << f1 << "-" << f2 << "Hz, order=" << order << ", fs=" << fs 
		 << " with a " << fir_t::window( window ) << " window\n";
	  design_bandstop_fir( order , fs , f1, f2 , window , true );
	}
      return;
    }
  
      
  if ( param.has( "lowpass" ) )
    {
      f1 = param.requires_dbl( "lowpass" );
      if ( use_kaiser )
	{
	  logger << " designing lowpass filter, " << f1 << "Hz, ripple=" << ripple[0] << ", tw=" << tw[0] << ", fs=" << fs << "\n"; 
	  design_lowpass_fir( ripple[0] , tw[0] , fs , f1 , true );
	}
      else
	{
	  logger << " designing lowpass filter, " << f1 << "Hz, order=" << order << ", fs=" << fs 
		 << " with a " << fir_t::window( window ) << " window\n";
	  design_lowpass_fir( order , fs , f1 , window , true );	  
	}
      return;
    }

  
  if ( param.has( "highpass" ) )
    {
      f1 = param.requires_dbl( "highpass" );
      if ( use_kaiser )
	{
	  logger << " designing highpass filter, " << f1 << "Hz, ripple=" << ripple[0] << ", tw=" << tw[0] << ", fs=" << fs << "\n"; 
	  design_highpass_fir( ripple[0] , tw[0] , fs , f1 , true );
	}
      else
	{
	  logger << " designing highpass filter, " << f1 << "Hz, order=" << order << ", fs=" << fs 
		 << " with a " << fir_t::window( window ) << " window\n";
	  design_highpass_fir( order , fs , f1 , window , true );
	}
      return;
    }
     
}
  


//
// Kaiser window specification
//

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


// make bandpass by convolving two separate filters
std::vector<double> dsptools::design_bandpass_fir( double ripple1 , double ripple2,
						   double tw1 , double tw2,
						   double fs , double f1 , double f2 , bool eval )
{
  
  std::vector<double> fc1 = design_highpass_fir( ripple1 , tw1, fs, f1 , false ); // F = no eval
  std::vector<double> fc2 = design_lowpass_fir( ripple2 , tw2, fs, f2 , false ); // F = no eval

  // convolve
  std::vector<double> fc = convolve( fc1 , fc2 );

  fir_t fir;

  if ( eval )
    {
      std::string label = "BANDPASS_HP_"
	+ Helper::dbl2str( f1 ) + "_" + Helper::dbl2str( ripple1 ) + "_" + Helper::dbl2str( tw1 )
	+ "_LP_" 
	+ Helper::dbl2str( f2 ) + "_" + Helper::dbl2str( ripple2 ) + "_" + Helper::dbl2str( tw2 );

      fir.outputFFT( label , fc , fs );
    }

  // for (int i=0;i<fc1.size();i++)
  //   std::cout << "f1\t" << fc1[i] << "\n"; 

  // for (int i=0;i<fc2.size();i++)
  //   std::cout << "f2\t" << fc2[i] << "\n"; 
  
  // for (int i=0;i<fc.size();i++)
  //   std::cout << "fc\t" << fc[i] << "\n"; 

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
// Fixed number of taps
//

std::vector<double> dsptools::design_bandpass_fir( int order , double fs , double f1 , double f2 , fir_t::windowType & window , bool eval )
{

  fir_t fir;

  // order must be even
  if ( order % 2 == 1 ) ++order;

  std::vector<double> fc = fir.create2TransSinc( order + 1 , f1 , f2, fs , fir_t::BAND_PASS );

  fc = fir.createWindow( &fc, window );

  if ( eval )
    {
      std::string label = "BANDPASS_" 
	+ Helper::dbl2str( f1 ) + ".." + Helper::dbl2str( f2 ) + "_" + Helper::int2str( order ) ;
      
      fir.outputFFT( label , fc , fs );
    }

  return fc;
}


std::vector<double> dsptools::design_bandstop_fir( int order , double fs , double f1 , double f2 , fir_t::windowType & window , bool eval )
{
  fir_t fir;

  if ( order % 2 == 1 ) ++order;

  std::vector<double> fc = fir.create2TransSinc( order + 1, f1 , f2, fs , fir_t::BAND_STOP );
 
  fc = fir.createWindow(&fc, window );

  if ( eval ) 
    {
      std::string label = "BANDSTOP_" 
	+ Helper::dbl2str( f1 ) + ".." + Helper::dbl2str( f2 ) + "_" + Helper::int2str( order ) ;
      
      fir.outputFFT( label , fc , fs );
    }

  return fc;

}



std::vector<double> dsptools::design_lowpass_fir( int order , double fs , double f , fir_t::windowType & window , bool eval )
{
  fir_t fir;

  if ( order % 2 == 1 ) ++order;

  std::vector<double> fc = fir.create1TransSinc( order + 1 , f , fs , fir_t::LOW_PASS );

  fc = fir.createWindow( &fc, window );

  if ( eval ) 
    {
      std::string label = "LOWPASS_" 
	+ Helper::dbl2str( f ) + "_" + Helper::int2str( order ) ;
      
      fir.outputFFT( label , fc , fs );
    }
  
  return fc;

}

std::vector<double> dsptools::design_highpass_fir( int order , double fs , double f , fir_t::windowType & window , bool eval )
{
  fir_t fir;

  if ( order % 2 == 1 ) ++order;

  std::vector<double> fc = fir.create1TransSinc( order + 1, f , fs , fir_t::HIGH_PASS );

  fc = fir.createWindow(&fc, window);
  
  if ( eval )
    {
      std::string label = "HIGHPASS_" 
	+ Helper::dbl2str( f ) + "_" + Helper::int2str( order ) ;
      
      fir.outputFFT( label , fc , fs );
    }

  return fc;

}




//
// apply FIR
//

std::vector<double> dsptools::apply_fir( const std::vector<double> & x ,
					 int fs,
					 fir_t::filterType ftype ,
					 int mode ,  // Kaiser window (1) or fixed order (2)
					 const std::vector<double> & ripple , const std::vector<double> & tw , // if using Kaiser window approach
					 double f1, double f2 ,             
					 int order , fir_t::windowType window ,  // if using fixed # of taps
					 const bool use_fft , const std::string & fir_file )
{

  //
  // Filter coefficients
  //

  std::vector<double> fc;

  
  //
  // Read from file? 
  //

  if ( ftype == fir_t::EXTERNAL )
    {
      if ( ! Helper::fileExists( fir_file ) ) Helper::halt( "could not find " + fir_file );
      std::ifstream IN1( fir_file.c_str() , std::ios::in );
      while ( ! IN1.eof() )
	{
	  double c;
	  IN1 >> c;
	  if ( IN1.eof() ) break;
	  fc.push_back( c );
	}
      IN1.close();
    }
  

  //
  // else, Kaiser window ?
  //
  
  bool kaiser = mode == 1 && ftype != fir_t::EXTERNAL; 
  bool fixed_order = mode != 1 && ftype != fir_t::EXTERNAL;
  bool kaiser_high_and_low = ripple.size() == 2 && tw.size() == 2; 
  
  if ( kaiser )
    {
      if ( ftype == fir_t::BAND_PASS ) 
	{
	  if ( kaiser_high_and_low ) 
	    fc = design_bandpass_fir( ripple[0] , ripple[1] , tw[0] , tw[1], fs , f1, f2 );
	  else
	    fc = design_bandpass_fir( ripple[0] , tw[0] , fs , f1, f2 );
	}
      else if ( ftype == fir_t::BAND_STOP )
	fc = design_bandstop_fir( ripple[0] , tw[0] , fs , f1, f2 );
      else if ( ftype == fir_t::LOW_PASS )
	fc = design_lowpass_fir( ripple[0] , tw[0] , fs , f1 );
      else if ( ftype == fir_t::HIGH_PASS )
	fc = design_highpass_fir( ripple[0] , tw[0] , fs , f1 );
    }
  else if ( fixed_order ) // fixed FIR order
    {
      if ( ftype == fir_t::BAND_PASS ) 
	fc = design_bandpass_fir( order, fs , f1, f2 , window );    
      else if ( ftype == fir_t::BAND_STOP )
	fc = design_bandstop_fir( order , fs , f1, f2 , window );
      else if ( ftype == fir_t::LOW_PASS )
	fc = design_lowpass_fir( order , fs , f1 , window );
      else if ( ftype == fir_t::HIGH_PASS )
	fc = design_highpass_fir( order , fs , f1 , window );
    }


  //
  // Apply FIR 
  //
  
  fir_impl_t fir_impl ( fc );

  return use_fft ? fir_impl.fft_filter( &x ) : fir_impl.filter( &x );
  
}


void dsptools::apply_fir( edf_t & edf , param_t & param )
{

  //
  // reading filter coefficients from an external file?
  //

  const bool from_file = param.has( "file" );

  std::string fir_file = "";

  //
  // Narrow-band filter via frequency-domain Gaussian
  // 
  
  const bool ngaus = param.has( "ngaus" );
  std::vector<double> npar;
  if ( ngaus )
    {    
      npar = param.dblvector( "ngaus" );
      if ( npar.size() !=2 ) Helper::halt( "expecting ngaus=<freq>,<fwhm>" );      
    }
  const double ngaus_f = ngaus ? npar[0] : 0 ;
  const double ngaus_fwhm = ngaus ? npar[1] : 0 ;
  
  
  //
  // FIR design method
  //

  const bool use_kaiser = param.has( "tw" ) || param.has( "ripple" );

  const bool fixed_order = ! ( use_kaiser || from_file || ngaus ) ;

  //
  // Kaiser-window specification: can be 1 or 2 values
  //
  
  std::vector<double> ripple, tw;
  if ( use_kaiser )
    {
      ripple = param.dblvector( "ripple" ) ;
      tw = param.dblvector( "tw" );
    }

  //
  // fixed-order
  //
  
  const int order = fixed_order ? param.requires_int( "order" ) : 0 ;

  //
  // Windowing if not Kaiser window
  //

  fir_t::windowType window = fir_t::HAMMING;
  if ( param.has( "rectangular" ) ) window = fir_t::RECTANGULAR;
  else if ( param.has( "bartlett" ) ) window = fir_t::BARTLETT;
  else if ( param.has( "hann" ) ) window = fir_t::HANN;
  else if ( param.has( "blackman" ) ) window = fir_t::BLACKMAN;

  
  //
  // transistition frequencues
  //
  
  double f1 , f2 ;

  
  //
  // Filter type
  //

  fir_t::filterType ftype = fir_t::BAND_PASS;
  
  
  //
  // standard convolution vs FFT implementation (default)
  //
  
  const bool use_fft = param.has( "fft" ) ? param.yesno( "fft" ) : true;
  
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
  else if ( param.has( "file" ) )
    {
      ftype = fir_t::EXTERNAL;
      fir_file = param.value( "file" );
    }	    
  else if ( ! ngaus ) 
    Helper::halt( "need to specify FIR type as bandpass, bandstop, lowpass, highpass (or file or ngaus)" );


  
  //
  // Signals
  //

  std::string signal_label = param.requires( "sig" );
  
  signal_list_t signals = edf.header.signal_list( signal_label );      

  std::vector<double> Fs = edf.header.sampling_freq( signals );
  
  const int ns = signals.size();


  logger << "  filtering channel(s):";
  
  //
  // Process each signal
  //

  for (int s=0; s<ns; s++)
    {

      //
      // skip annotation channels
      //
      
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
      
      logger << " " << signals.label(s);

      if ( ngaus )
	apply_ngaus( edf, signals(s) , ngaus_f, ngaus_fwhm );
      else	
	apply_fir( edf , signals(s) , ftype , use_kaiser ? 1 : 2 ,
		   ripple, tw ,
		   f1 , f2 ,
		   order , window , 
		   use_fft , fir_file );
      
    }
  logger << "\n";
}


void dsptools::apply_ngaus( edf_t & edf , int s , const double ngaus_f , const double ngaus_fwhm )
{

  
  interval_t interval = edf.timeline.wholetrace();

  //  std::cout <<" int " << interval.start << " " << interval.stop << "\n";
    
  slice_t slice( edf , s , interval );

  //  std::cout << " got a slice\n";
  
  const std::vector<double> * d = slice.pdata();
  
  int fs = edf.header.sampling_freq( s );

  //  std::cout << " F1\n";
  std::vector<double> filtered = narrow_gaussian_t::filter( *d , fs, ngaus_f , ngaus_fwhm ) ;
  //  std::cout << " F2\n";  
  edf.update_signal( s , &filtered );
  //  std::cout << " F3\n";    
}


void dsptools::apply_fir( edf_t & edf , int s , fir_t::filterType ftype ,
			  int mode ,
			  const std::vector<double> & ripple , const std::vector<double> & tw ,
			  double f1, double f2 ,
			  int order , fir_t::windowType window , 
			  const bool use_fft , const std::string & fir_file )
{
      
  
  interval_t interval = edf.timeline.wholetrace();
  
  //
  // Pull entire signals out
  //
  
  slice_t slice( edf , s , interval );
  
  const std::vector<double> * d = slice.pdata();
      
  int fs = edf.header.sampling_freq( s );
      
  //
  // Design/apply FIR
  // 
  
  std::vector<double> filtered = apply_fir( *d , fs , ftype ,
					    mode , 
					    ripple, tw ,
					    f1 , f2 ,
					    order , window , 
					    use_fft ,
					    fir_file );  
    
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
      {
	//	std::cout << " KW " << n << " " << out[n] << " " << (*in)[n] << " " << out[n] * (*in)[n] << "\n";
	out[n] *= (*in)[n];
      }
  
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
  const double padding_sec = 2; // default = 2 
  double sz = window.size() / (double)sampFreq + padding_sec ;
  fir_impl_t fir_impl ( window );
  std::vector<double> xx0( sampFreq * sz );
  xx0[sampFreq*(sz/2.0)-1] = 1;
  std::vector<double> xx = fir_impl.filter( &xx0 ); 
  double idx0 = sampFreq*(sz/2.0)-1;
  
  //
  // Step response
  //

  // same 2 window second around filter size
  fir_impl_t fir_impl2 ( window );
  std::vector<double> xx1( sampFreq * sz , 1 );
  for (int i=sampFreq*(sz/2.0)-1;i<sampFreq * sz;i++) xx1[i] = 0;
  std::vector<double> ss = fir_impl2.filter( &xx1 ); 
  double idx1 = sampFreq*(sz/2.0)-1;

  //
  // Output IR and SR
  //
  
  // SR
  double integ = xx[0];
  
  for (int xi=0;xi<xx.size();xi++)
    {
      double tp = 1.0/sampFreq * (xi - idx1 ); 
      writer.level( Helper::dbl2str(tp) , "SEC" );
      writer.value( "IR" , xx[xi] );

      writer.value( "SR" , integ );
      integ += xx[xi];
      
      // writer.value( "SR2" , ss[xi] ); 
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

  writer.unlevel( "FIR" );
  
  return result;
 

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
  FFT fftx( Nfft , Nfft , 1 );
  fftx.apply( x );
  std::vector<dcomp> rfftx = fftx.transform();
  
  FFT ffth( Nfft , Nfft , 1 );
  ffth.apply( h );
  std::vector<dcomp> rffth = ffth.transform();
  
  // convolution in the frequency domain
  std::vector<dcomp> y( Nfft );
  for (int i=0;i<rfftx.size();i++) y[i] = rfftx[i] * rffth[i]; 
  
  // inverse FFT

  FFT ifft( Nfft , Nfft , 1 , FFT_INVERSE );
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
 
