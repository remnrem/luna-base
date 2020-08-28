
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

#ifndef __CWT_H__
#define __CWT_H__

#include "miscmath/miscmath.h"
#include "helper/helper.h"
#include "defs/defs.h"

#include <vector>
#include <cmath>
#include <iostream>

void run_cwt();

class CWT {
  
 public:
  
  CWT() 
    { 
      init(); 
    }; 

  void set_verbose( bool b ) { verbose = b; } 

  std::vector<double> get_timeframe() const { return time; }  

  void set_timeframe( const double f ) 
  {
    if ( srate == 0 ) Helper::halt( "srate not set in cwt" );
    time.clear();
    
    // generates a nice range, based on the frequency of the wavelet
    double T = 50.0 / f;
    double start = -T/2.0;
    double stop  =  T/2.0;    
    double inc = 1.0/(double)srate;    
    for (double t = start; t <= stop-inc ; t+= inc )
      time.push_back(t);    
    if ( time.size() % 2 ) // i.e. if odd, force wavelet to be even
      time.push_back( stop );

    // now set all wavelet-specific factors
    n_wavelet            = time.size();    
    n_convolution        = n_wavelet+n_data-1;
    n_conv_pow2          = MiscMath::nextpow2(n_convolution);
    half_of_wavelet_size = n_wavelet/2;  // note integer division, n_wavelet forced to be even 

    
    /* std::cout << "time size " << time.size() << "\n" */
    /* 	      << n_wavelet << "\n" */
    /* 	      << n_convolution << "\n" */
    /* 	      << n_conv_pow2 << "\n" */
    /* 	      << half_of_wavelet_size << "\n\n"; */
    

  }


  static double pick_fwhm( double f , double m = -0.7316762 , double c = 1.1022791 )
  {
    return exp( log(f) * m + c );
  }
  
  void set_sampling_rate( const double sr ) 
  { 
    srate = sr; 
  }
  
  std::vector<dcomp> wavelet(const int);
  
  void add_wavelet(const double _fc, const int n_cycles ) 
  {

    fc.push_back( _fc );
    
    // s = n/2pif
    // Fb = 2s^2
    
    const double s = n_cycles / (double)( 2 * M_PI * _fc );    

    sig.push_back(s); // not used
    fb.push_back( 2 * s*s );
    
    num_frex = fc.size();

  }


  //
  // alternate specification of wavelets, based on FWHM
  //

  std::vector<dcomp> alt_wavelet(const int);

  double alt_empirical_fwhm( const int fi );

  void add_wavelets( double minFreq, double minFreqTempRes, double maxFreq, double maxFreqTempRes, double numFreqs, double waveLength )
  {

    alt_spec = true;
    
    fc = MiscMath::logspace( minFreq , maxFreq , numFreqs );
    fwhm = MiscMath::logspace( minFreqTempRes , maxFreqTempRes , numFreqs );
    wlen.clear();
    wlen.resize( fc.size() , waveLength );

    num_frex = fc.size();
	
  }

  void alt_add_wavelet( double Freq , double TempRes , double waveLength )
  {
    alt_spec = true;
    fc.push_back( Freq );
    fwhm.push_back( TempRes );
    wlen.push_back( waveLength );
    num_frex = fc.size();
  }
  
  void set_pnts_trials( const int p , const int t )
  {
    // trials / epochs here 
    num_pnts   = p; 
    num_trials = t; // or 'epochs'
    if ( data->size() != num_pnts * num_trials ) Helper::halt( "bad pnts/trials, does not match data[] in CWT()" );
  }
  
  void load( const std::vector<double> * d ) 
  {
    
    data = d;
    
    n_data               = data->size();  // should equal 'pnts x trials'

    // by default, do not assume any 'stacking' of similar slices
    // this can be changed by calling set_pnts_trials() subsequently
    num_trials = 1;
    num_pnts = data->size();
    
    verbose = false;

  }
  
  void run();
  
  void run_wrapped();
  
  double freq(const int fi) const { return fc[fi]; }
  int    points() const { return num_pnts; }
  int    freqs() const { return num_frex; }
  double result(const int fi, const int ti) const { return eegpower[fi][ti]; }
  double raw_result(const int fi, const int ti) const { return rawpower[fi][ti]; }
  const std::vector<double> & results(const int fi) const { return rawpower[fi]; }

  // same as above 'results' function

  std::vector<double> amplitude(const int fi) const { return rawpower[fi]; }
  std::vector<double> phase(const int fi) const { return ph[fi]; }
  std::vector<dcomp> get_complex( const int fi ) { return conv_complex[fi] ; } 

  
  //
  // get-functions
  //

  std::vector<double> get_freq() const { return fc; } 
  int get_num_freq() const { return num_frex; }
  double get_srate() const { return srate; }

  //
  // options
  //

  void store_real_imag_vectors( const int b ) { store_real_imag = b; } 

  void use_alt() { alt_spec = true; }

  void alt_timeline( double t ) { set_timeframe( 50.0 / t ); }

  
 private:
  
  //
  // Frequency steps
  //

  int    num_frex;

  //
  // Timeframe and Sampling rate
  //
  
  int                 srate;
  std::vector<double> time;

  //
  // Alternate specification: 
  //
  
  bool                alt_spec; 

  
  
  
  //
  // Wavelet parameters
  //

  std::vector<double> fc;
  std::vector<double> fb;
  std::vector<double> sig;

  std::vector<double> fwhm;
  std::vector<double> wlen;
  
  //
  // Segments
  //
  
  int num_pnts;   // number of points per EPOCH 
  int num_trials;  // number of trials (always 1)
  
  
  //
  // Convolution parameters
  //
  
  int n_wavelet;           
  int n_data;              
  int n_convolution;       
  int n_conv_pow2;         
  int half_of_wavelet_size;


  //
  // Input data
  //

  const std::vector<double> * data;


  //
  // Transformation (output)
  //

  std::vector<std::vector<double> > eegpower;
  std::vector<std::vector<double> > rawpower;

  //
  // Phase information (output)
  //

  std::vector<std::vector<double> > ph; 

  //
  // Optionally, store complex wavelet transform
  //

  bool store_real_imag;
  std::vector<std::vector<dcomp> > conv_complex; 
  
  //
  // Misc
  //

  bool verbose;

  void init()
  {
    alt_spec = false;
    fc.clear();
    fb.clear();
    srate = 0;
    num_pnts = num_trials = 1;
    eegpower.clear();
    rawpower.clear();
    conv_complex.clear();
    store_real_imag = false;
  }
  
};

#endif
