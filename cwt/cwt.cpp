
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

#include "miscmath/miscmath.h"
#include "fftw/fftwrap.h"


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
  dcomp k( 1.0 / sqrt( fb[fi] * M_PI ) , 0 );  
  
  for (int i=0;i<n;i++)
    {
      w[i] = k 
	* exp( dcomp( 2 * M_PI * fc[fi] * time[i] , 0 ) * dcomp(0,1) ) 
	* exp( dcomp(  - (time[i]*time[i]) / fb[fi] , 0 ) );      
    }
  return w;
}

std::vector<dcomp> CWT::alt_wavelet( const int fi )
{

  // alternate formulat/parameterization of wavelets
  // based on 'time' variable

  // depending on two parameters:
  //  FWHM is a bandwidth parameter
  //  Fc is a wavelet center frequency

  // also get empirical FWHM
  
  const int n = time.size();

  std::vector<dcomp> w(n);
  
  for (int i=0;i<n;i++)
    {
      w[i] = 
	exp( dcomp( 0 , 2 * M_PI * fc[fi] * time[i] ) ) 
	* exp( dcomp(  -4 * log(2.0) * (time[i]*time[i]) / (fwhm[fi]*fwhm[fi]) , 0 ) );      
    }

   // for (int i=0;i<n;i++)
   //   std::cout << "wavelet " << fc[fi] << "\t"
   // 	       << dcomp( 0 , 2 * M_PI * fc[fi] * time[i] ) << "\t"
   // 	       << exp( dcomp( 0 , 2 * M_PI * fc[fi] * time[i] ) ) << "\t"
   // 	       << std::exp( dcomp( 0 , 2 * M_PI * fc[fi] * time[i] ) ) << "\n";

   //       std::real(w[i]) << "\t" << std::imag(w[i] )  << "\n";

  return w;
}


double CWT::alt_empirical_fwhm( const int fi )
{

  // for the alternate parameterization of CWT,
  // get empirical time-domain FWHM (from the Gaussian)
  
  // depending on two parameters:
  //  FWHM is a bandwidth parameter
  //  Fc is a wavelet center frequency

  const int n = time.size();
  
  std::vector<double> g(n);
  
  for (int i=0;i<n;i++)
    g[i] = exp( -4 * log(2.0) * (time[i]*time[i]) / (fwhm[fi]*fwhm[fi]) );      
    
  int mid_idx = MiscMath::nearest_idx( g , 1 );
  int lwr_idx = MiscMath::nearest_idx( g , 0.5 , 0 , mid_idx );
  int upr_idx = MiscMath::nearest_idx( g , 0.5 , mid_idx , -1 );
  
  return time[ upr_idx ] - time[ lwr_idx ];
  
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


  if ( store_real_imag )
    {
      conv_complex.resize( num_frex );
    }
  
  const double SQRT_PI = sqrt(M_PI);


  //
  // Does the timeframe vary for different wavlets?  If not, precompute FFT of data just once
  //

  std::set<double> tf;
  for (int i=0;i<wlen.size();i++) tf.insert(wlen[i]);
  bool fixed_wlen = alt_spec && tf.size() == 1;
    
  //
  // Initial FFT of data 
  //

  real_FFT eegfft;
  std::vector<dcomp> eegfftX;

  if ( fixed_wlen ) 
    {
      set_timeframe( 50.0 / wlen[0] );      
      eegfft.init( data->size() , n_conv_pow2 , srate );
      eegfft.apply( *data );
      eegfftX = eegfft.transform();
    }

  //
  // loop through frequencies and compute synchronization
  //

  for (int fi=0;fi<num_frex;fi++)
    {
      
       //
      // Set timeline for this wavelet (or is fixed already under alt_spec) 
      //
	    
      if ( ! alt_spec ) 
	{
	  set_timeframe( fc[fi] );  
	}
      else
	set_timeframe( 50.0 / wlen[fi] );

      
      //
      // Do we need FFT of data?
      //

      if ( ! fixed_wlen ) 
	{
	  if ( fi != 0 ) eegfft.reset();
	  eegfft.init( data->size() ,  n_conv_pow2 , srate );
	  eegfft.apply( *data );	  
	  eegfftX = eegfft.transform();	  
	}

      
      //
      // Generate wavelet 
      //

      std::vector<dcomp> w = alt_spec ? alt_wavelet(fi) : wavelet(fi);
      
      // std::cerr << "n_conv_pow2 = " << n_conv_pow2 << "\n"
      // 		<< "n_data = " << n_data << "\n"
      // 		<< "n_convolution " << n_convolution << "\n"
      // 		<< "half_of_wavelet_size " << half_of_wavelet_size << "\n"
      // 		<< "n_wavelet " << n_wavelet << "\n";
      
      
      //
      // First FFT
      //
      
      
      FFT fft1( w.size() , n_conv_pow2 , 1 );
      
      fft1.apply( w );
      
      std::vector<dcomp> wt = fft1.transform();
      

      //
      // Scaling factor to ensure similar amplitudes of original traces and wavelet-filtered signal
      // kernelFFT = 2*kernelFFT./max(kernelFFT);
      //
      
      const bool rescale_wavelet = true;

      if ( rescale_wavelet ) 
	{
	  dcomp max = MiscMath::max( wt );
	  
	  for ( int i = 0 ; i < wt.size() ; i++ )
	    wt[i] = ( dcomp( 2, 0 ) * wt[i] ) / max;
	}


      
      //
      // Convolution in the frequency domain 
      //
      std::vector<dcomp> y( n_conv_pow2 );

      for (int i=0;i<eegfftX.size();i++) y[i] = eegfftX[i] * wt[i]; 

      //
      // Inverse FFT, normalized by 1/Nfft, back to time-domain
      //

      FFT ifft( n_conv_pow2 , n_conv_pow2 , 1 , FFT_INVERSE );

      ifft.apply( y );

      std::vector<dcomp> eegconv_tmp = ifft.scaled_transform() ;
      
      
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
	ph[fi][i] = atan2( eegconv[i].imag() , eegconv[i].real() );
      
      
      //
      // optionally, extract real/imag parts
      //

      if ( store_real_imag )
	{	  
	  // nb. num_trials == 1 always... 
	  conv_complex[fi].resize( num_pnts * num_trials ); 
	  for (int i=0; i<num_pnts*num_trials; i++)
	    conv_complex[fi][i] = eegconv[i] ;

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




void CWT::run_wrapped()
{

  //
  // Alternate parameterization, for a wrapped wavelet and fixed time-frame
  // following Cox & Fell
  //

  //
  // Initialize
  //
  
  // num_frex  : number of wavelets to apply

  eegpower.resize( num_frex );
  rawpower.resize( num_frex );
  ph.resize( num_frex ); 
  for (int i=0;i<num_frex;i++) 
    {
      eegpower[i].resize( num_pnts , 0 );
      rawpower[i].resize( num_pnts , 0 );
      ph[i].resize( num_pnts , 0 );
      conv_complex.resize( num_frex );
    }

  const double SQRT_PI = sqrt(M_PI);
 
  int Ldata = data->size();

  int Ltapr = time.size();

  int Lconv1 = Ldata + Ltapr - 1;  

  int Lconv  = MiscMath::nextpow2( Lconv1 );

  // for wrapped wavelet
  int middlePaddingLength = Lconv - Ltapr ;

  int half1_idx = floor( Ltapr / 2 ) - 1 ;
  int half2_idx = floor( Ltapr / 2 ) ;

  // std::cout << " sizes : Ldata " << Ldata << "\n"
  // 	    << " Ltapr " << Ltapr << "\n"
  // 	    << " Lconv1 " << Lconv1 << "\n"
  // 	    << " Lconv " << Lconv << "\n"
  // 	    << " middlePaddingLength " << middlePaddingLength << "\n";
  
  // std::cout << "hh = " << half1_idx << " " << half1_idx << " " << Ltapr << "\n";
  
  //
  // Initial FFT of data 
  //
      
  FFT eegfft( data->size() , Lconv , srate );
  eegfft.apply( *data );
  std::vector<dcomp> eegfftX = eegfft.transform();

  //
  // loop through frequencies and compute synchronization
  //

  for (int fi=0;fi<num_frex;fi++)
    {
      
      //
      // Generate wavelet 
      //

      set_timeframe( 50.0 / wlen[fi] );

      std::vector<dcomp> w0 = alt_wavelet(fi) ;
      
      //
      // Wrap wavelet? ( second half, padding, first half )
      //
      
      std::vector<dcomp> w( Lconv , dcomp(0,0) );
      int cc = 0;
      for (int i= half2_idx; i < Ltapr ; i++ ) w[ cc++ ] = w0[ i ];
      // skip middle zeros
      cc += middlePaddingLength;
      for (int i= 0; i <= half1_idx ; i++ ) w[ cc++ ] = w0[ i ];
      
      // padded, unwrapped
      // std::vector<dcomp> w( Lconv , dcomp(0,0) );
      // for (int i= 0; i <= w0.size() ; i++ ) w[ i ] = w0[i];
      
      
      
      //
      // FFT of wrapped wavelet
      //
      
      FFT kernelFFT( w.size() , Lconv , 1 , FFT_FORWARD );

      kernelFFT.apply( w );

      std::vector<dcomp> wt = kernelFFT.transform();

      //
      // Scaling factor to ensure similar amplitudes of original traces and wavelet-filtered signal
      // kernelFFT = 2*kernelFFT./max(kernelFFT);
      //

      dcomp max = MiscMath::max( wt );

      for ( int i = 0 ; i < wt.size() ; i++ )
       	wt[i] = ( dcomp( 2, 0 ) * wt[i] ) / max;


      //
      // Convolution in the frequency domain 
      //
      
      std::vector<dcomp> y( Lconv );
      for (int i=0; i<Lconv; i++)
	y[i] = eegfftX[i] * wt[i]; 
      
               
      //
      // Inverse FFT back to time-domain
      //

      FFT ifft( y.size() , Lconv , 1 , FFT_INVERSE );
      ifft.apply( y );

      // nb. using scaled transform here
      std::vector<dcomp> eegconv = ifft.scaled_transform();
      
      //
      // Trim
      // 

      //  m = m(1:Lconv1);
      eegconv.resize( Lconv1 );
      
      // m = m(1:end-Ltapr+1);  
      eegconv.resize( Lconv1 - Ltapr + 1 ); 
	

      //
      // extract phase from the convolution
      //

      for (int i=0; i < Ldata ; i++)
	ph[fi][i] = atan2( eegconv[i].imag() , eegconv[i].real() );


      //
      // extract real/imag parts
      //

      conv_complex[fi] = eegconv;

      //
      // wavelet power: abx(X)^2, dB/normalized and raw
      //
      
      for (int i=0;i<num_pnts;i++)
	{
	  rawpower[fi][i] = pow( abs( eegconv[i] ) , 2 );
	}

      //
      // Baseline normalized
      //

      double baseline = MiscMath::mean( rawpower[fi] );

      // i.e. express as dB over entire night, i.e. 10log10(ratio)
      for (int i=0; i<num_pnts; i++) eegpower[fi][i] = 10*log10( rawpower[fi][i]/baseline );
        

      //
      // next frequency/wavelet
      //
    }
    
}
