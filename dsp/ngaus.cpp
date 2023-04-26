
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

#include "dsp/ngaus.h"
#include "helper/helper.h"
#include "miscmath/miscmath.h"
#include "fftw/fftwrap.h"

// implements a narrow-band filter via frequency-domain Gaussian
// based on a Matlab implementation by MX Cohen

Eigen::VectorXd narrow_gaussian_t::filter( const Eigen::VectorXd & d , int sr, double f , double fwhm )
{
  const int n = d.size();  
  std::vector<double> v2(n);  
  Eigen::VectorXd::Map( &v2[0], d.size()) = d;
  v2 = narrow_gaussian_t::filter( v2 , sr, f, fwhm );
  Eigen::Map<Eigen::VectorXd> r( &(v2[0]) , n );
  return r;
}

std::vector<double> narrow_gaussian_t::filter( const std::vector<double> & d , int sr, double f , double fwhm )
{

  const int n = d.size();
  std::vector<double> hz = MiscMath::linspace( 0, sr, n );
  const int nf = hz.size();
  
  // normalized width
  double s = fwhm * ( 2 * M_PI - 1.0) / ( 4 * M_PI ); 

  // shifted frequencies
  std::vector<double> x( nf );
  for (int i=0; i<nf; i++) x[i] = hz[i] - f;
  
  // Gaussain
  std::vector<double> fx( nf );
  for (int i=0; i<nf; i++) fx[i] = exp( -0.5 * pow( x[i] / s , 2 ) );

  // Gain normalized
  double max_fx = 0;
  for (int i=0; i<nf; i++)
    if ( fx[i] > max_fx ) max_fx = fx[i];
  for (int i=0; i<nf; i++)
    fx[i] /= max_fx;
  
  // filtdat = 2*real( ifft( bsxfun(@times,fft(data,[],2),fx) ,[],2) );
  real_FFT fft;
  fft.init( n , n , sr );
  fft.apply( d );
  std::vector<dcomp> X = fft.transform();

  for (int i=0; i<nf; i++)
    X[i] = X[i] * fx[i];

  real_iFFT ifft( n, n , sr );
  ifft.apply( X );
  return ifft.inverse();

}
