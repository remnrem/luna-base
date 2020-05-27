
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

#include "fftw/cohfft.h"
#include "fftw/fftwrap.h"

#include "edf/edf.h"
#include "timeline/timeline.h"
#include "miscmath/miscmath.h"

#include "defs/defs.h"


precoh_t coherence_t::precoh;


void precoh_t::prepare( coherence_t * coh ,
			const int s ,
			const std::vector<double> & x )
  
{ 
  
  //
  // Initial FFT to get frequencies
  //

  if ( frq.size() == 0 )
    {
      FFT fft0( coh->segment_points , coh->Fs , FFT_FORWARD , coh->window );
      if ( coh->average_adj ) fft0.average_adjacent();

      coh->N = fft0.cutoff;
      coh->res.resize( coh->N ); // set output bucket size
      
      frq.resize( coh->N );
      for (int f=0;f< coh->N;f++) frq[f] = fft0.frq[f];
    }
  

  //
  // Accumulate PSD for X
  //

  std::vector<std::vector<std::complex<double> > > & psd_x = psd[s];

  psd_x.resize(coh->N);

  
  //
  // Iterate over segments, performing individual FFT in each
  //

  for (int p = 0; p <= coh->total_points - coh->segment_points ; p += coh->segment_increment_points )
    {
      
      if ( p + coh->segment_points > coh->total_points )
	Helper::halt( "internal error in coherence()" );
      
      FFT fftx( coh->segment_points , coh->Fs , FFT_FORWARD , coh->window );
      
      if ( coh->detrend || coh->zerocenter )
 	{	  
 	  std::vector<double> x1( coh->segment_points );
 	  for (int j=0;j< coh->segment_points;j++) x1[j] = x[p+j];
	  if  ( coh->detrend ) MiscMath::detrend(&x1);
	  else MiscMath::centre(&x1);
	  fftx.apply( x1 );
 	}
       else
 	{
 	  fftx.apply( &(x[p]) , coh->segment_points );
 	}
    
      if ( coh->average_adj )
	{
	  fftx.average_adjacent();
	}
      
    
      cutoff = fftx.cutoff;
      
      // x2 is to get full spectrum  
      normalisation_factor = 2 * fftx.normalisation_factor;

      for (int i=0;i<cutoff;i++)
	{
	  double a = fftx.out[i][0];
	  double b = fftx.out[i][1];      	  
	  std::complex<double> Xx( a , b );
	  psd_x[i].push_back( Xx );	  
	}
      
    } // next segment

}



void coherence_t::process( const int s1 , const int s2 )
{

  //
  // Accumulate PSD for X, Y and cross-spectra, by freq, then 
  //
  
  std::vector<std::vector<std::complex<double> > > & cmp_x = precoh.psd[s1];
  std::vector<std::vector<std::complex<double> > > & cmp_y = precoh.psd[s2];
  
  std::vector<std::vector<std::complex<double> > > cpsd( cmp_x.size() );
  std::vector<std::vector<double> > psd_x( precoh.cutoff ) ;
  std::vector<std::vector<double> > psd_y( precoh.cutoff ) ;

  //
  // Iterate over frequencies
  //

  for (int i=0;i< precoh.cutoff;i++)
    {

      const int ns = cmp_x[i].size();
      
      for (int j=0;j<ns;j++)
	{
	  std::complex<double> & Xx = cmp_x[i][j];	  
	  double a = std::real(Xx);
	  double b = std::imag(Xx);
	  psd_x[i].push_back( ( a*a + b*b ) * precoh.normalisation_factor );
	  
	  std::complex<double> & Yy = cmp_y[i][j];
	  a = std::real(Yy);
	  b = std::imag(Yy);
	  psd_y[i].push_back( ( a*a + b*b ) * precoh.normalisation_factor );
 
	  std::complex<double> Xy = Xx * conj( Yy );	 
	  cpsd[i].push_back( precoh.normalisation_factor * Xy );	  
	}
      
    } // next segment
  

  
  //
  // take average over segments
  //

  const double COH_EPS = 1e-10;

  for (int i=0;i<psd_x.size();i++)
    {

      // https://mne.tools/stable/generated/mne.connectivity.spectral_connectivity.html

      double sxx = MiscMath::mean( psd_x[i] );
      
      double syy = MiscMath::mean( psd_y[i] );
      
      std::complex<double> sxy = MiscMath::mean( cpsd[i] );

      
      //
      // track these
      //

      res.sxx[i] = sxx;
      res.syy[i] = syy;
      res.sxy[i] = sxy;

      // bad?

      res.bad[i] = sxx < COH_EPS || syy < COH_EPS ;


      // // cross/auto spectra (in dB)
      
      // if ( phi2 > COH_EPS )
      // 	res.cross_spectrum[i] = 5.0*log10(phi2);
      // else
      // 	{
      // 	  res.cross_spectrum[i] =  -50.0 ;
      // 	  res.bad[i] = true;
      // 	}

      // if ( sxx > COH_EPS )
      // 	res.auto_spectrum1[i] = 10.0*log10(sxx) ;
      // else
      // 	{
      // 	  res.auto_spectrum1[i] = -100.0 ;
      // 	  res.bad[i] = true;
      // 	}
      
      // if (  syy > COH_EPS )
      // 	res.auto_spectrum2[i] = 10.0*log10(syy) ;
      // else
      // 	{
      // 	  res.auto_spectrum2[i] = -100.0 ;
      // 	  res.bad[i] = true;
      // 	}
      
      // res.cross_norm1[i] = phi2 / (sxx*sxx) ;
      // res.cross_norm2[i] = phi2 / (syy*syy) ;

      
    } // next frequency

  
}


