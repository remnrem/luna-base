
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

#include "dsp/freqsampling.h"
#include "helper/helper.h"
#include "miscmath/miscmath.h"
#include "fftw/fftwrap.h"
#include "stats/statistics.h"
#include "miscmath/crandom.h"

// simple generation of a random signal designing a filter to match the FFT 
// of an original signal, then feed Guassian white noise through that

std::vector<double> freq_sampl_t::generate( const std::vector<double> & d , int sr , const double flwr , const double fupr )
{
  
  const int n = d.size();
  std::vector<double> hz = MiscMath::linspace( 0, sr, n );
  const int nf = hz.size();
  
  // match on mean and SD
  const double smean = MiscMath::mean( d );
  const double sdev  = MiscMath::sdev( d , smean );
  
  real_FFT fft;
  fft.init( n , n , sr );
  fft.apply( d );
  std::vector<dcomp> X = fft.transform();

  std::vector<double> fx( n );
  for (int i=0; i<nf; i++) 
    fx[i] = std::real(X[i]) ;
  //fx[i] = fft.mag[i];
  
  // zero-out frequencies?
  if ( flwr >= 0 ) 
    for (int i=0; i<n; i++) 
      if ( hz[i] < flwr ) fx[i] = 0; 

  if ( fupr >= 0 ) 
    for (int i=0; i<n; i++) 
      if ( hz[i] > fupr ) fx[i] = 0; 
  
  // Gain normalized
  if ( 1 ) 
    {
      double max_fx = 0;
      for (int i=0; i<nf; i++)
	if ( fx[i] > max_fx ) max_fx = fx[i];
      for (int i=0; i<nf; i++)
	fx[i] /= max_fx;
    }

  // filtdat = 2*real( ifft( bsxfun(@times,fft(data,[],2),fx) ,[],2) );
  
  // gaussian white noise
  std::vector<double> z( n );
  for (int i=0; i<n; i++) 
    z[i] = Statistics::ltqnorm( CRandom::rand() ); 

  real_FFT fft2;
  fft2.init( n , n , sr );
  fft2.apply( z );
  std::vector<dcomp> Zt = fft2.transform();
  
  for (int i=0; i<nf; i++)
    Zt[i] = Zt[i] * fx[i];
  
  real_iFFT ifft( n, n , sr );
  ifft.apply( Zt );
  std::vector<double> rdat = ifft.inverse();
  
  // scale to same as original, mean & SD
  // doubtless super inefficient here and we could 
  // get scaling upfront but there you go...

  const double smean2 = MiscMath::mean( rdat );
  const double sdev2  = MiscMath::sdev( rdat , smean2 );

  for (int i=0; i<n; i++) 
    {
      rdat[i] = ( rdat[i] - smean2 ) / sdev2;
      rdat[i] = ( rdat[i] * sdev ) + smean ;
    }

  std::cout << " checks " << smean << " " << sdev << " --> " << smean2 << " " << sdev2 << "  fin " << MiscMath::mean( rdat ) << " " << MiscMath::sdev( rdat  ) << "\n";


  // check
  // do FFT of signal back again
  //
  real_FFT fftc;
  fftc.init( n , n , sr );
  fftc.apply( rdat );
  std::vector<dcomp> Xc = fftc.transform();
  
  // for (int i=0; i<nf; i++)
  //   std::cout << " hz[" << i << "] = " << hz[i] << "\t" << fx[i] << "\txx " << fftc.mag[i] << " " << std::real( Xc[i] ) << "\n";
  // std::cout << "\n\n";

  return rdat;
  
}
