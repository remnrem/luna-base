
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
  else if ( window == WINDOW_HANN )    w = MiscMath::hann_window(N);
  else if ( window == WINDOW_HAMMING ) w = MiscMath::hamming_window(N);
  
  for (int i=0;i<N;i++) normalisation_factor += w[i] * w[i];
  normalisation_factor *= Fs;  
  normalisation_factor = 1.0/normalisation_factor;
  //  std::cerr << "norm " << normalisation_factor << "\n";
  
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
      
      //      std::cout << "det " << mag[i] << "\t" << X[i] << "\t" << ( a*a + b*b ) << "\n";
      
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
  
  // std::cout << "segment_size_points = " << segment_size_points << "\n"
  // 	    << "noverlap_points = " << noverlap_points << "\n"
  //   	    << "segment_increment_points = " << segment_increment_points << "\n";
  
  //
  // Initial FFT
  //
  
  FFT fft0( segment_size_points , Fs , FFT_FORWARD , window );

  if ( average_adj ) 
    fft0.average_adjacent();
  
  psd.resize( fft0.cutoff , 0 );
  N = fft0.cutoff;

  freq.resize(N);
  for (int f=0;f<N;f++)
    {
      freq[f] = fft0.frq[f];
      //      std::cout << "Wf " << f << "\t" << freq[f] << "\n" ;
    }

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
      
      FFT fft( np , Fs , FFT_FORWARD , WINDOW_HANN );
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


int bin_t::bin( const std::vector<double> & f , 
		const std::vector<double> & y ) 
{

  
  if ( f.size() != y.size() ) Helper::halt( "bin_t internal error" );
  
  bfa.clear();
  bfb.clear();  
  bspec.clear();
  
  if ( f.size() < 2 ) return 0;


  //
  // no binning?
  //
  
  if ( w == 0 )
    {
      for (int i=0;i<f.size();i++)
        {
          if ( f[i] > mx_f ) break;
          bfa.push_back( f[i] );
          bfb.push_back( f[i] );
          bspec.push_back( y[i] );
        }
      return bspec.size();
    }

  // assume always from 0, DC component  
  if ( f[0] == 0 ) 
    {
      bspec.push_back( y[0] );
      bfa.push_back( 0 );
      bfb.push_back( 0 );
    }

  //for (int ii=0;ii<f.size();ii++) std::cout << "Frq " << f[ii] << "\n";
  
  double nyquist = 0.5 * Fs;
    
  int num_freqs = f.size();
  
  double df = f[1] - f[0];    

  if ( w/df  < 1.0 ) Helper::halt( "bin resolution too small: min " + Helper::dbl2str( df ) );

  int freqwin = (int) ( w / df ) ;      
  
  if ( mx_f > nyquist ) mx_f = nyquist; 
  
  for (int i = 1; i < num_freqs ; i += freqwin)
      {
	
	double tem = 0.0;
	
	int k = 0;
	
	for (int j = i ; j < i + freqwin ; j++) 
	  {
	    
	    if (j > 0 && j < num_freqs - 1) // skip DC and Nyquist
	      {	      	      
		if ( f[j] <= mx_f )
		  {
		    std::cout << "adding " << f[j] << "\n";
		    tem += y[j];
		    k++;
		  }
	      }
	  }
	
	std::cout << "scanning " << f[i] << " to " << f[ i + freqwin -1 ] <<  " " << k << "\n";	

	if ( k > 0 ) 
	  {	  
	    bspec.push_back( tem/(double)k );
	    bfa.push_back( f[i-1] ); // less than 
	    bfb.push_back( f[i+k-1] ); // greater than or equal to
	    std::cout << "from " <<  f[i-1]  << " " << f[i+k-1] << "\n";
	  }
	
      }   
    
    return bspec.size();
  }
