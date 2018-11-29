
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

#include "cwt.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>

#include "../miscmath/miscmath.h"
#include "../fftw/fftwrap.h"

std::vector<dcomp> CWT::wavelet( const int fi )
{
  
  // based on 'time' variable

  // Definition: a complex Morlet wavelet is
  // cmor(x) = (pi*Fb)^{-0.5} * exp(2*i*pi*Fc*x) * exp(-(x^2)/Fb)
  
  // depending on two parameters:
  //  Fb is a bandwidth parameter
  //  Fc is a wavelet center frequency

  const int n = time.size();

  std::vector<dcomp> w(n);
  
  // constant
  dcomp k( 1.0 / sqrt( fb[fi] * M_PI ) , 0 );  // From Matlab
  
  for (int i=0;i<n;i++)
    {
      w[i] = k 
	* exp( dcomp( 2 * M_PI * fc[fi] * time[i] , 0 ) * dcomp(0,1) ) 
	* exp( dcomp(  - (time[i]*time[i]) / fb[fi] , 0 ) );
      double rw = std::real(w[i]);
      // double cw = std::imag(w[i]);
      
      // std::cout << "wavelet " << " " << sig[fi] << " " << fb[fi] << " " << fc[fi] 
      // 		<< " " << i << " " << rw  << " " << cw << "\n";

    }
  return w;
}

void CWT::run()
{


  //
  // Any baseline normalization?
  //

  bool baseline_normalization = true;
  

  
  //
  // Initialize
  //
  
  // Frex x ( time x trials )
  
  eegpower.resize( num_frex );
  rawpower.resize( num_frex );
  ph.resize( num_frex ); 
  for (int i=0;i<num_frex;i++) 
    {
      eegpower[i].resize( num_pnts , 0 );
      rawpower[i].resize( num_pnts , 0 );
      ph[i].resize( num_pnts , 0 );
    }

  const double SQRT_PI = sqrt(M_PI);

  //
  // loop through frequencies and compute synchronization
  //

  for (int fi=0;fi<num_frex;fi++)
    {
      
      // Set timeline for this wavelet, then generate the wavelet
      
      set_timeframe( fc[fi] );      

      std::vector<dcomp> w = wavelet(fi);
      
      
      //
      // Initial FFT
      //
      
      FFT eegfft( n_conv_pow2 , srate );
      eegfft.apply( *data );
      std::vector<dcomp> eegfftX = eegfft.transform();
            
      //
      // First FFT
      //
      
      FFT fft1( n_conv_pow2 , 1 , FFT_FORWARD );
      fft1.apply( w );
      std::vector<dcomp> wt = fft1.transform();


      //
      // Convolution in the frequency domain 
      //

      std::vector<dcomp> y( n_conv_pow2 );
      for (int i=0;i<eegfftX.size();i++) y[i] = eegfftX[i] * wt[i]; 
      
      //
      // Inverse FFT back to time-domain
      //

      FFT ifft( n_conv_pow2 , 1 , FFT_INVERSE );
      ifft.apply( y );
      std::vector<dcomp> eegconv_tmp = ifft.transform();
      
      dcomp denom( 1.0/(double)n_conv_pow2 , 0 );

      //
      // Normalize
      //

      for (int i=0;i<n_conv_pow2;i++) eegconv_tmp[i] *= denom;

      //
      // Trim
      //

      eegconv_tmp.resize( n_convolution );      
      std::vector<dcomp> eegconv;
      for (int i=half_of_wavelet_size-1;
	   i < ( n_convolution - half_of_wavelet_size ); 
	   i++ ) 
	eegconv.push_back( eegconv_tmp[i] );

      //
      // extract phase from the convolution
      //

      for (int i=0; i<num_pnts*num_trials; i++)
	{
	  ph[fi][i] = atan2( eegconv[i].imag() , eegconv[i].real() );
	}
      

      //
      // Put results back into pnts x trials matrix; take power
      // abs(X)^2; average over trials to get a pnts-length vector of
      // average power
      //
      
      int cnt = 0;
      std::vector<double> temppower( num_pnts , 0 );
      for (int i=0; i<num_pnts; i++)
	{
	  double x = 0;
	  for (int t=0; t<num_trials; t++)
	    x += pow( abs( eegconv[ cnt + t*num_pnts ] ) , 2 ); 
	  ++cnt;
	  temppower[i] = num_trials > 1 ? x / (double)num_trials : x ;
	}
      
      //
      // Record in freq x time-point matrix; use the 'baseline
      // correction based on 'all' time-points, i.e. to get dB
      //

      double baseline       = 0;
      int    baseline_n     = 0;
      int    baseline_start = 0;
      int    baseline_stop  = num_pnts; // 1 past index
      
      if ( baseline_normalization )
	{

	  for (int i = baseline_start; i < baseline_stop; i++ ) { baseline += temppower[i]; baseline_n++; } 
	  baseline /= (double)baseline_n;
	  
	  // i.e. express as dB over entire night, i.e. 10log10(ratio)
	  for (int i=0; i<num_pnts; i++) eegpower[fi][i] = 10*log10( temppower[i]/baseline );
	}
      else
	{
	  for (int i=0; i<num_pnts; i++) eegpower[fi][i] = 10*log10( temppower[i] );
	}

      // save non-dB version too
      rawpower[fi] = temppower;
    }
    
}
