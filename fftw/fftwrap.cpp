
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

#include "fftwrap.h"

#include "edf/edf.h"
#include "timeline/timeline.h"
#include "miscmath/miscmath.h"

#include "defs/defs.h"

FFT::FFT( int N , int Fs , fft_t type , window_function_t window )
  : N(N) , Fs(Fs), type(type), window(window), in(NULL), out(NULL), p(NULL)
{
  
  // Allocate storage for input/output
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  if ( in == NULL ) Helper::halt( "FFT failed to allocate input buffer" );
  
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  if ( out == NULL ) Helper::halt( "FFT failed to allociate output buffer" );

  // Initialise (probably not necessary, but do anyway)
  for (int i=0;i<N;i++) { in[i][0] = in[i][1] = 0; }
  
  // Generate plan
  p = fftw_plan_dft_1d( N, in, out , type == FFT_FORWARD ? FFTW_FORWARD : FFTW_BACKWARD , FFTW_ESTIMATE );

  //
  // We want to return only the positive spectrum, so set the cut-off
  //
  
  cutoff = N % 2 == 0 ? N/2+1 : (N+1)/2 ;
  X.resize(cutoff,0);
  mag.resize(cutoff,0);
  frq.resize(cutoff,0);

  //
  // Scale frequencies appropriate (not used in calculation, just for output/convenience)
  //

  double T = N/(double)Fs;
  for (int i=0;i<cutoff;i++) frq[i] = i/T;
  
  //
  // Normalisation factor for PSD  (1/value)
  // i.e. equiv. to 1/(N.Fs) in unweighted case, otherwise
  // we take the window into account
  //

  w.resize( N , 1 ); // i.e. default of no window
  
  normalisation_factor = 0;  
  if      ( window == WINDOW_TUKEY50 ) w = MiscMath::tukey_window(N,0.5);
  else if ( window == WINDOW_HANNING ) w = MiscMath::hanning_window(N);
  else if ( window == WINDOW_HANN )    w = MiscMath::hann_window(N);
  else if ( window == WINDOW_HAMMING ) w = MiscMath::hamming_window(N);
  
  for (int i=0;i<N;i++) normalisation_factor += w[i] * w[i];
  normalisation_factor *= Fs;  
  normalisation_factor = 1.0/normalisation_factor;

} 

bool FFT::apply( const std::vector<double> & x )
{
  return apply( &(x[0]) , x.size() );
}
  

bool FFT::apply( const double * x , const int n )
{
  
  //
  // Load up (windowed) input buffer
  //
  
  if ( window == WINDOW_NONE )
    for (int i=0;i<n;i++) { in[i][0] = x[i];        in[i][1] = 0; } 
  else
    for (int i=0;i<n;i++) { in[i][0] = x[i] * w[i]; in[i][1] = 0; } 
  
  //
  // Execute actual FFT
  // 
  
  fftw_execute(p);
  

  //
  // Calculate PSD
  //

  //
  // psdx = (1/(Fs*N)) * abs(xdft).^2;
  // where abs() is complex sqrt(a^2+b^2)
  //

  for (int i=0;i<cutoff;i++)
    {
      
      double a = out[i][0];
      double b = out[i][1];
      
      X[i] =  ( a*a + b*b ) * normalisation_factor;
      mag[i] = sqrt( a*a + b*b );
 
      // not for DC and Nyquist, but otherwise
      // double all entries (i.e. to preserve
      // total power, as here we have the one-
      // sided PSD
      
      if ( i > 0 && i < cutoff-1 ) X[i] *= 2;
      
    }
  
  return true;

}


bool FFT::apply( const std::vector<std::complex<double> > & x )
{
  
  const int n = x.size();
  
  if ( n > N ) Helper::halt( "error in FFT" );
  
  for (int i=0;i<n;i++)
    {
      in[i][0] = std::real( x[i] );
      in[i][1] = std::imag( x[i] );	
    }    
  
  fftw_execute(p);
  

  //
  // Calculate PSD
  //

  //
  // psdx = (1/(Fs*N)) * abs(xdft).^2;
  // where abs() is complex sqrt(a^2+b^2)
  //

  for (int i=0;i<cutoff;i++)
    {
      
      double a = out[i][0];
      double b = out[i][1];
      
      X[i] =  ( a*a + b*b ) * normalisation_factor;
      mag[i] = sqrt( a*a + b*b );

      // not for DC and Nyquist, but otherwise
      // double all entries (i.e. to preserve
      // total power, as here we have the one-
      // sided PSD
      
      if ( i > 0 && i < cutoff-1 ) X[i] *= 2;
      
    }
  
  return true;

}


std::vector<std::complex<double> > FFT::transform() const
{
  std::vector<std::complex<double> > r(N);
  for (int i=0;i<N;i++) 
    r[i] = std::complex<double>( out[i][0] , out[i][1] );
  return r;
}

std::vector<std::complex<double> > FFT::scaled_transform() const
{
  const double fac = 1.0 / (double)N;
  std::vector<std::complex<double> > r(N);
  for (int i=0;i<N;i++) 
    r[i] = std::complex<double>( out[i][0] * fac , out[i][1] * fac );
  return r;
}

std::vector<double> FFT::inverse() const
{
  // from an IFFT, get the REAL values and divide by N, i.e. this
  // should mirror the input data when the input data are REAL  
  std::vector<double> r(N);
  for (int i=0;i<N;i++) r[i] = out[i][0] / (double)N;
  return r;
}


void FFT::average_adjacent()
{
  
  // nb. SpectralTrainFig does not average the frequency values but
  // rather gives the higer one of each pair (except DC)
  
  // need to reset 'cut-off' also 
  std::vector<double> frq2;
  std::vector<double> X2;
  
  // keep first (DC/0Hz) value
  frq2.push_back( frq[0] );
  X2.push_back( X[0] );
  
  for (int i=1;i<cutoff;i+=2)
    {
      frq2.push_back( frq[i+1] );  // Note: not average; this mirrors DD
      X2.push_back( ( X[i] + X[i+1] ) / 2.0 );
    }

  X = X2;
  frq = frq2;
  cutoff = X.size();

}


bool FFT::add( frequency_band_t band , double f )
{
  const double lwr = globals::freq_band[ band ].first;
  const double upr = globals::freq_band[ band ].second;
  return f > lwr && f <= upr; 
}


double FFT::width( frequency_band_t band )
{
  return globals::freq_band[ band ].second - globals::freq_band[ band ].first;
}


void PWELCH::process()
{
  
  // From MATLAB parameterizatopm:
  //  K = (M-NOVERLAP)/(L-NOVERLAP)
  //    M = total_points = epoch size (in data-points)
  //    L = segment_size_points = segment size (in data-points)
  //    K = noverlap_segments = desired number of (overlapping) segments of size 'L' within 'M'
  //    NOVERLAP = noverlap_points 
    
  int total_points             = data.size();
  int segment_size_points      = M * Fs;   // 'nfft' in Matlab

  int noverlap_points          = noverlap_segments > 1 
    ? ceil( ( noverlap_segments*segment_size_points - total_points  ) / double( noverlap_segments - 1 ) )
    : 0 ;

  int segment_increment_points = segment_size_points - noverlap_points;
  
//   std::cout << "noverlap_points = " << noverlap_points << "\n"
//    	    << "segment_increment_points = " << segment_increment_points << "\n";
  
  //
  // Initial FFT
  //
  
  FFT fft0( segment_size_points , Fs , FFT_FORWARD , window );

  if ( average_adj ) 
    fft0.average_adjacent();
  
  psd.resize( fft0.cutoff , 0 );
  N = fft0.cutoff;

  freq.resize(N);
  for (int f=0;f<N;f++) freq[f] = fft0.frq[f];
  

  //
  // Iterate over segments, performing individual FFT in each
  //

  int segments = 0;
  
//    std::cout << "total_points = " << total_points << "\n"
//    	    << "segment_size_points (NFFT) = " << segment_size_points << "\n"
//    	    << "segment_increment_points = " << segment_increment_points << "\n"
//     	    << "total_points - segment_size_points = " << total_points - segment_size_points << "\n";
  
  for (int p = 0; p <= total_points - segment_size_points ; p += segment_increment_points )
    {
      
      ++segments;

      
//       std::cout << "seg " << segments << "\t" 
//        		<< p << " -- " << p + segment_size_points - 1 << "\n";
      
      // note -- this assumes no zero-padding will be applied, 
      // and all segments passed must be of exactly size segment_size_points
     
      
      FFT fft( segment_size_points , Fs , FFT_FORWARD , window );
      
      if ( p + segment_size_points > data.size() ) 
	Helper::halt( "internal error in pwelch()" );
      
      const bool detrend = false;
      const bool zerocentre = false;
      if ( detrend )
	{
	  std::vector<double> y( segment_size_points );
	  for (int j=0;j<segment_size_points;j++) y[j] = data[p+j];
	  MiscMath::detrend(&y);      
	  fft.apply( y );
	}
      else if ( zerocentre )
	{
	  std::vector<double> y( segment_size_points );
	  for (int j=0;j<segment_size_points;j++) y[j] = data[p+j];
	  MiscMath::centre(&y);      
	  fft.apply( y );
	}
      else
	{
	  fft.apply( &(data[p]) , segment_size_points );
	}
      
      if ( average_adj )
	fft.average_adjacent();
      
      int cutoff = fft.cutoff;
      
      for (int i=0;i<fft.cutoff;i++)
	psd[i] += fft.X[i];
      
    } // next segment
  

  //
  // take average over segments
  //

  for (int i=0;i<psd.size();i++)
    {
      psd[i] /= (double)segments;      
      //std::cout << "PWELCH " << fft0.frq[i] << " --> " << psd[i] << "\n";     
    }
  
}



void PWELCH::psdsum( std::map<freq_range_t,double> * f )
{  
  std::map<freq_range_t,double>::iterator ii = f->begin();
  while ( ii != f->end() )
    {
      ii->second = psdsum( ii->first.first , ii->first.second ) ;      
      ++ii;
    }
}

void PWELCH::psdmean( std::map<freq_range_t,double> * f )
{

  std::map<freq_range_t,double>::iterator ii = f->begin();
  while ( ii != f->end() )
    {
      
      const double & lwr = ii->first.first;
      const double & upr = ii->first.second; 

      // add is l <= x < y 
      double r = 0;
      int c = 0;
      for (int i=0;i<N;i++) 
	{
	  if ( freq[i] >= upr ) break;
	  if ( freq[i] >= lwr ) { ++c; r += psd[i]; }	  
	}
      
      ii->second = r / (double) c;
      ++ii;
    }
}



std::map<double,double> fft_spectrum( const std::vector<double> * d , int Fs )
{

  // gives relative log-power spectrum (used in scope for quick plotting)

  std::map<double,double> results;

  // using Welch if a large number of datapoints
  const int     np = d->size();
  const double sec = np/(double)Fs;
  
  // standalone FFT 
  if ( sec <= 60 ) 
    {
      
      FFT fft( np , Fs , FFT_FORWARD , WINDOW_HANNING );
      fft.apply( *d );
      int N = fft.cutoff;
      
      for (int f=0;f<N;f++)
	{
	  if ( fft.frq[f] > 0.5 && fft.frq[f] < 30 ) 
	    {
	      int freq = fft.frq[f];
	      double x = log( fft.X[f] );
	      results[freq] += x;
	    }
	}

      
      //
      // gives relative log power 
      //

      double mx = -99999;
      double mn = 99999;
      
      std::map<double,double>::iterator ii = results.begin();
      while ( ii != results.end() ) 
	{
	  //std::cout << "orig = " << ii->first << " " << ii->second << "\n";
	  double x = ii->second;
	  if ( x < mn ) mn = x;
	  if ( x > mx ) mx = x;
	  ++ii;
	}
      
      ii = results.begin();
      while ( ii != results.end() ) 
	{
	  ii->second = ( ii->second - mn ) / ( mx - mn ) ;
	  ++ii;
	}
            
    }
  return results;
}





void coherence_t::process()
{

  //
  // Initial FFT to get frequencies
  //
  
  FFT fft0( segment_points , Fs , FFT_FORWARD , window );

  if ( average_adj ) fft0.average_adjacent();
  
  N = fft0.cutoff;

  res.resize(N);

  // freqs.
  for (int f=0;f<N;f++) res.frq[f] = fft0.frq[f];

  //
  // Accumulate PSD for X, Y and cross-spectra, by freq, then 
  //

  std::vector<std::vector<double> > psd_x(N);
  std::vector<std::vector<double> > psd_y(N);
  std::vector<std::vector<std::complex<double> > > cpsd(N);

  
  //
  // Iterate over segments, performing individual FFT in each
  //

  int segments = 0;
  
  for (int p = 0; p <= total_points - segment_points ; p += segment_increment_points )
    {
      
      ++segments;

//       std::cout << "seg " << segments << "\t" 
// 		<< p << " -- " << p + segment_points - 1 << "\n";
      
      if ( p + segment_points > total_points )
	Helper::halt( "internal error in coherence()" );
      
      FFT fftx( segment_points , Fs , FFT_FORWARD , window );
      FFT ffty( segment_points , Fs , FFT_FORWARD , window );      
      
      if ( detrend || zerocenter )
 	{	  
 	  std::vector<double> x1( segment_points );
 	  for (int j=0;j<segment_points;j++) x1[j] = x[p+j];
	  if  ( detrend ) MiscMath::detrend(&x1);
	  else MiscMath::centre(&x1);
	  fftx.apply( x1 );

 	  std::vector<double> y1( segment_points );
 	  for (int j=0;j<segment_points;j++) y1[j] = y[p+j];
 	  if ( detrend ) MiscMath::detrend(&y1);      
	  else MiscMath::centre(&y1);      
 	  ffty.apply( y1 ); 	  
 	}
       else
 	{
 	  fftx.apply( &(x[p]) , segment_points );
	  ffty.apply( &(y[p]) , segment_points );
 	}

      if ( average_adj )
	{
	  fftx.average_adjacent();
	  ffty.average_adjacent();
	}

    
      int cutoff = fftx.cutoff;
      
      // x2 is to get full spectrum  
      double normalisation_factor = 2 * fftx.normalisation_factor;

      for (int i=0;i<cutoff;i++)
	{
	  // nb, normalization factors will be the same
	  
	  double a = fftx.out[i][0];
	  double b = fftx.out[i][1];
      	  psd_x[i].push_back( ( a*a + b*b ) * normalisation_factor );
	  std::complex<double> Xx( a , b );

	  a = ffty.out[i][0];
	  b = ffty.out[i][1];
	  psd_y[i].push_back( ( a*a + b*b ) * normalisation_factor );
	  std::complex<double> Yy( a , b );

	  std::complex<double> Xy = Xx * conj( Yy );
	  
	  cpsd[i].push_back( normalisation_factor * Xy );

	}

    } // next segment
  
  
  //
  // take average over segments
  //
  const double COH_EPS = 1e-10;

  for (int i=0;i<psd_x.size();i++)
    {

      double sxx = MiscMath::mean( psd_x[i] );
      double syy = MiscMath::mean( psd_y[i] );
      std::complex<double> sxy = MiscMath::mean( cpsd[i] );

      // calculate magnitude squared coherence and cross/auto spectra (in dB)
      double phi = abs( sxy ) ;
      double phi2 = phi * phi;
      
      if ( sxx < COH_EPS || syy < COH_EPS ) 
	res.coh[i] = -9;
      else
	res.coh[i] = phi2 / ( sxx * syy ); 
      
      // truncate?
      res.cross_spectrum[i] = phi2 > COH_EPS ? 5.0*log10(phi2) : -50.0 ;	  
      res.auto_spectrum1[i] = sxx > COH_EPS ?  10.0*log10(sxx) : -100.0 ;
      res.auto_spectrum2[i] = syy > COH_EPS ?  10.0*log10(syy) : -100.0 ;
      
      res.cross_norm1[i] = phi2 / (sxx*sxx) ;
      res.cross_norm2[i] = phi2 / (syy*syy) ;
  
	//      std::cerr << "res " << res.frq[i] << " " << res.coh[i] << "\n";
    }

}

