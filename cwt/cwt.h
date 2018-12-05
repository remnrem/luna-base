
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
    sig.push_back(s);
    fb.push_back( 2 * s*s );
    
    num_frex = fc.size();

  }


/*   void set_wavelet_param( const double f , const int a = 3 , const int b = 10 ) */
/*   { */
/*     //    Helper::halt( "not used currently..." ); */

/*     // here 'a' and 'b' are the number of cycles requested (that can change as a function of frequency */
/*     // i.e. here, it cycles between 3 and 10 (mapping onto the example frequencies in frex[0] .. freqx[end] */

/*     // This uses the relation s = n/(2pi.f) as given in the Cohen text, page 145 */
    
/*     // we could set a==b here (especially in the case of a single frex */
/*     // (i.e. frex.size() == 1 as in the spindle example, most likely) */
/*     // here num_frex == 1 too */

/*     s = MiscMath::logspace( a, b , num_frex); */
/*     for (int i=0;i<num_frex;i++) s[i] /= 2*M_PI*frex[i]; */
/*   } */

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
  
  double freq(const int fi) const { return fc[fi]; }
  int    points() const { return num_pnts; }
  int    freqs() const { return num_frex; }
  double result(const int fi, const int ti) const { return eegpower[fi][ti]; }
  double raw_result(const int fi, const int ti) const { return rawpower[fi][ti]; }
  const std::vector<double> & results(const int fi) const { return rawpower[fi]; }  
  std::vector<double> phase(const int fi) const { return ph[fi]; }
  
  //
  // get-functions
  //

  std::vector<double> get_freq() const { return fc; } 
  int get_num_freq() const { return num_frex; }
  double get_srate() const { return srate; }

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
  // Wavelet parameters
  //

  std::vector<double> fc;
  std::vector<double> fb;

  std::vector<double> sig;

  //
  // Segments
  //
  
  int num_pnts;   // number of points per EPOCH (640)
  int num_trials;  // number of EPOCHS (99)
  
  
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
  // Misc
  //

  bool verbose;

  void init()
  {
    fc.clear();
    fb.clear();
    srate = 256;
    num_pnts = num_trials = 1;
    eegpower.clear();
    rawpower.clear();
  }
  
};

#endif
