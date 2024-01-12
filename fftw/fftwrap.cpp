
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

#include "fftw3.h"                                                                                                                                                                       

#include "edf/edf.h"
#include "timeline/timeline.h"
#include "miscmath/miscmath.h"

#include "helper/helper.h"
#include "miscmath/dynam.h"

#include "defs/defs.h"

#include "db/db.h"

extern writer_t writer;


//
// Complex FFT
//

void FFT::reset() 
{
  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);
}

FFT::~FFT() 
{    
  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);
}


void FFT::init( int Ndata_, int Nfft_, int Fs_ , fft_t type_ , window_function_t window_ )
{
  
  Ndata = Ndata_;
  Nfft = Nfft_;
  Fs = Fs_;
  type = type_;
  window = window_;

  if ( Ndata > Nfft ) Helper::halt( "Ndata cannot be larger than Nfft" );

  // Allocate storage for input/output
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nfft);
  if ( in == NULL ) Helper::halt( "FFT failed to allocate input buffer" );
  
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nfft);
  if ( out == NULL ) Helper::halt( "FFT failed to allociate output buffer" );

  // Initialise (probably not necessary, but do anyway)
  for (int i=0;i<Nfft;i++) { in[i][0] = in[i][1] = 0; }
  
  // Generate plan
  p = fftw_plan_dft_1d( Nfft, in, out , type == FFT_FORWARD ? FFTW_FORWARD : FFTW_BACKWARD , FFTW_ESTIMATE );

  //
  // We want to return only the positive spectrum, so set the cut-off
  //
  
  cutoff = Nfft % 2 == 0 ? Nfft/2+1 : (Nfft+1)/2 ;
  X.resize(cutoff,0);
  mag.resize(cutoff,0);
  frq.resize(cutoff,0);

  //
  // Scale frequencies appropriately (not used in calculation, just for output)
  //

  double T = Nfft/(double)Fs;

  for (int i=0;i<cutoff;i++) frq[i] = i/T;
  
  //
  // Normalisation factor for PSD  (1/value)
  // i.e. equiv. to 1/(N.Fs) in unweighted case, otherwise
  // we take the window into account
  //

  w.resize( Ndata , 1 ); // i.e. default of no window
  
  normalisation_factor = 0;  
  if      ( window == WINDOW_TUKEY50 ) w = MiscMath::tukey_window(Ndata,0.5);
  else if ( window == WINDOW_HANN )    w = MiscMath::hann_window(Ndata);
  else if ( window == WINDOW_HAMMING ) w = MiscMath::hamming_window(Ndata);
  
  for (int i=0;i<Ndata;i++) normalisation_factor += w[i] * w[i];
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
    for (int i=0;i<Ndata;i++) { in[i][0] = x[i];  in[i][1] = 0; } 
  else
    for (int i=0;i<Ndata;i++) { in[i][0] = x[i] * w[i]; in[i][1] = 0; } 

  for (int i=Ndata;i<Nfft;i++) { in[i][0] = 0;  in[i][1] = 0;  } 
  
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
  
  if ( n > Nfft ) Helper::halt( "error in FFT" );
  
  for (int i=0;i<Ndata;i++)
    {
      in[i][0] = std::real( x[i] );
      in[i][1] = std::imag( x[i] );	
    }    

  // zero-pad any remainder
  for (int i=Ndata;i<Nfft;i++)
    {
      in[i][0] =  in[i][1] = 0;
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
  std::vector<std::complex<double> > r(Nfft);
  for (int i=0;i<Nfft;i++) 
    r[i] = std::complex<double>( out[i][0] , out[i][1] );
  return r;
}

std::vector<std::complex<double> > FFT::scaled_transform() const
{
  // or Ndata?  check
  const double fac = 1.0 / (double)Nfft;
  std::vector<std::complex<double> > r(Nfft);
  for (int i=0;i<Nfft;i++) 
    r[i] = std::complex<double>( out[i][0] * fac , out[i][1] * fac );
  return r;
}

std::vector<double> FFT::inverse() const
{
  // from an IFFT, get the REAL values and divide by N, i.e. this
  // should mirror the input data when the input data are REAL  
  std::vector<double> r(Nfft);
  for (int i=0;i<Nfft;i++) r[i] = out[i][0] / (double)Nfft;
  return r;
}

std::vector<double> FFT::unscaled_inverse() const
{
  std::vector<double> r(Nfft);
  for (int i=0;i<Nfft;i++) r[i] = out[i][0] ;
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



// --------------------------------------------------------------------------
//
// Real 1D DFT (real to complex) 
//
// --------------------------------------------------------------------------

void real_FFT::reset() 
{
  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);
}

real_FFT::~real_FFT() 
{    
  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);
}


void real_FFT::init( int Ndata_, int Nfft_, int Fs_ , window_function_t window_ )
{
  
  Ndata = Ndata_;
  Nfft = Nfft_;
  Fs = Fs_;
  window = window_;

  if ( Ndata > Nfft ) Helper::halt( "Ndata cannot be larger than Nfft" );

  // Allocate storage for input/output
  in = (double*) fftw_malloc(sizeof(double) * Nfft);
  if ( in == NULL ) Helper::halt( "FFT failed to allocate input buffer" );
  
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nfft);
  if ( out == NULL ) Helper::halt( "FFT failed to allociate output buffer" );

  // Initialise (probably not necessary, but do anyway)
  for (int i=0;i<Nfft;i++) { in[i] = 0; }
  
  // Generate plan: nb. r2c 1D plan
  p = fftw_plan_dft_r2c_1d( Nfft, in, out , FFTW_ESTIMATE ) ;

  // We want to return only the positive spectrum, so set the cut-off  
  cutoff = Nfft % 2 == 0 ? Nfft/2+1 : (Nfft+1)/2 ;
  X.resize(cutoff,0);
  mag.resize(cutoff,0);
  frq.resize(cutoff,0);

  //
  // Scale frequencies appropriately (not used in calculation, just for output)
  //

  double T = Nfft/(double)Fs;

  for (int i=0;i<cutoff;i++) frq[i] = i/T;
  
  //
  // Normalisation factor for PSD  (1/value)
  // i.e. equiv. to 1/(N.Fs) in unweighted case, otherwise
  // we take the window into account
  //

  w.resize( Ndata , 1 ); // i.e. default of no window
  
  normalisation_factor = 0;  
  if      ( window == WINDOW_TUKEY50 ) w = MiscMath::tukey_window(Ndata,0.5);
  else if ( window == WINDOW_HANN )    w = MiscMath::hann_window(Ndata);
  else if ( window == WINDOW_HAMMING ) w = MiscMath::hamming_window(Ndata);
  
  for (int i=0;i<Ndata;i++) normalisation_factor += w[i] * w[i];
  normalisation_factor *= Fs;  
  normalisation_factor = 1.0/normalisation_factor;
    
} 

bool real_FFT::apply( const std::vector<double> & x )
{
  return apply( &(x[0]) , x.size() );
}
  

bool real_FFT::apply( const double * x , const int n )
{

  //
  // Load up (windowed) input buffer
  //
  
  if ( window == WINDOW_NONE )
    for (int i=0;i<Ndata;i++) in[i] = x[i];  
  else
    for (int i=0;i<Ndata;i++) in[i] = x[i] * w[i];

  // zero-padding
  for (int i=Ndata;i<Nfft;i++) in[i] = 0;  
  

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


std::vector<std::complex<double> > real_FFT::transform() const
{
  std::vector<std::complex<double> > r(Nfft);
  for (int i=0;i<Nfft;i++) 
    r[i] = std::complex<double>( out[i][0] , out[i][1] );
  return r;
}

std::vector<std::complex<double> > real_FFT::scaled_transform() const
{
  // or Ndata?  check
  const double fac = 1.0 / (double)Nfft;
  std::vector<std::complex<double> > r(Nfft);
  for (int i=0;i<Nfft;i++) 
    r[i] = std::complex<double>( out[i][0] * fac , out[i][1] * fac );
  return r;
}

std::vector<double> real_FFT::inverse() const
{
  // from an IFFT, get the REAL values and divide by N, i.e. this
  // should mirror the input data when the input data are REAL  
  std::vector<double> r(Nfft);
  for (int i=0;i<Nfft;i++) r[i] = out[i][0] / (double)Nfft;
  return r;
}


void real_FFT::average_adjacent()
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


bool real_FFT::add( frequency_band_t band , double f )
{
  const double lwr = globals::freq_band[ band ].first;
  const double upr = globals::freq_band[ band ].second;
  return f > lwr && f <= upr; 
}


double real_FFT::width( frequency_band_t band )
{
  return globals::freq_band[ band ].second - globals::freq_band[ band ].first;
}



// --------------------------------------------------------------------------
//
// Real 1D inverse DFT  (complex -> real)
//
// --------------------------------------------------------------------------

void real_iFFT::reset() 
{
  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);
}

real_iFFT::~real_iFFT() 
{    
  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);
}

void real_iFFT::init( int Ndata_, int Nfft_, int Fs_ , window_function_t window_ )
{
  
  Ndata = Ndata_;
  Nfft = Nfft_;
  Fs = Fs_;
  window = window_;

  if ( Ndata > Nfft ) Helper::halt( "Ndata cannot be larger than Nfft" );

  // Allocate storage for input/output
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nfft);
  if ( in == NULL ) Helper::halt( "FFT failed to allociate output buffer" );

  out = (double*) fftw_malloc(sizeof(double) * Nfft);
  if ( out == NULL ) Helper::halt( "FFT failed to allocate input buffer" );
  
  // Initialise (probably not necessary, but do anyway)
  for (int i=0;i<Nfft;i++) { in[i][0] = in[i][1] = 0; }
  
  // Generate plan: nb. c2r 1D plan  
  p = fftw_plan_dft_c2r_1d( Nfft, in, out , FFTW_ESTIMATE );

  // We want to return only the positive spectrum, so set the cut-off  
  cutoff = Nfft % 2 == 0 ? Nfft/2+1 : (Nfft+1)/2 ;
  X.resize(cutoff,0);
  mag.resize(cutoff,0);
  frq.resize(cutoff,0);

  //
  // Scale frequencies appropriately (not used in calculation, just for output)
  //

  double T = Nfft/(double)Fs;

  for (int i=0;i<cutoff;i++) frq[i] = i/T;
  
  //
  // Normalisation factor for PSD  (1/value)
  // i.e. equiv. to 1/(N.Fs) in unweighted case, otherwise
  // we take the window into account
  //

  w.resize( Ndata , 1 ); // i.e. default of no window
  
  normalisation_factor = 0;  
  if      ( window == WINDOW_TUKEY50 ) w = MiscMath::tukey_window(Ndata,0.5);
  else if ( window == WINDOW_HANN )    w = MiscMath::hann_window(Ndata);
  else if ( window == WINDOW_HAMMING ) w = MiscMath::hamming_window(Ndata);
  
  for (int i=0;i<Ndata;i++) normalisation_factor += w[i] * w[i];
  normalisation_factor *= Fs;  
  normalisation_factor = 1.0/normalisation_factor;
    
} 


bool real_iFFT::apply( const std::vector<std::complex<double> > & x )
{

  const int n = x.size();
  
  if ( n > Nfft ) Helper::halt( "error in FFT" );
  
  for (int i=0;i<Ndata;i++)
    {
      in[i][0] = std::real( x[i] );
      in[i][1] = std::imag( x[i] );	
    }    

  // zero-pad any remainder
  for (int i=Ndata;i<Nfft;i++)
    {
      in[i][0] =  in[i][1] = 0;
    }

  fftw_execute(p);

  //
  // Calculate PSD
  //

  //
  // psdx = (1/(Fs*N)) * abs(xdft).^2;
  // where abs() is complex sqrt(a^2+b^2)
  //

  // for (int i=0;i<cutoff;i++)
  //   {
      
  //     double a = out[i][0];
  //     double b = out[i][1];
      
  //     X[i] =  ( a*a + b*b ) * normalisation_factor;
  //     mag[i] = sqrt( a*a + b*b );

  //     // not for DC and Nyquist, but otherwise
  //     // double all entries (i.e. to preserve
  //     // total power, as here we have the one-
  //     // sided PSD
      
  //     if ( i > 0 && i < cutoff-1 ) X[i] *= 2;
      
  //    }

  return true;

}




std::vector<double> real_iFFT::inverse() const
{
  // from an IFFT, get the REAL values and divide by N, i.e. this
  // should mirror the input data when the input data are REAL  
  std::vector<double> r(Nfft);
  for (int i=0;i<Nfft;i++) r[i] = out[i] / (double)Nfft;
  return r;
}

std::vector<double> real_iFFT::unscaled_inverse() const
{
  std::vector<double> r(Nfft);
  for (int i=0;i<Nfft;i++) r[i] = out[i] ;
  return r;
}




// --------------------------------------------------------------------------
//
// Welch algorithm
//
// --------------------------------------------------------------------------

void PWELCH::process()
{

  //
  // note:  for 'PSD' we specify the segment overlap (note increment)
  //   from that, we get the implied number of segments;
  //   but here, we take only the implied number of segments, and recalculate
  //   the overlap/increment;  here, we stretch the increment out if needed -
  //   it will be good as long as we can have an integer number of points
  //    i.e. 128 Hz
  //      segement-sec=5 segment-overlap=2  w/ a 30-second epoch
  //      implies a 3-second increment ( = 3 * 128 = 384 samople points)

  //      0 - 5
  //      3 - 8
  //      6 - 11
  //      9 - 14
  //      12 - 17
  //      15 - 20
  //      18 - 23
  //      21 - 26
  //      24 - 29
  //
  //  so, this would be short, i.e. needs to go to 30;  however, '9' segments is passed to PWELCH, which
  //  then figures that a 3.125 second = 400 (integer) sample point increment makes this work;  so, changes
  //  the actual extent of overlap.

  //segment 1; p = 0 .. 640
  //segment 2; p = 400 .. 1040
  //segment 3; p = 800 .. 1440
  //segment 4; p = 1200 .. 1840
  //segment 5; p = 1600 .. 2240
  //segment 6; p = 2000 .. 2640
  //segment 7; p = 2400 .. 3040
  //segment 8; p = 2800 .. 3440
  //segment 9; p = 3200 .. 3840

  // i.e. given 9 segments, changes increment to 3.125 (400 points) rather than 3 seconds, which
  // evenly takes us to the end of the 30-sec window;
  
  // in contrast, if we asked for: segement-sec=5 segment-overlap=3  w/ a 30-second epoch
  //   implies a 2-second increment: however, here it is not possible to get an integer number
  //   of sample points to span the whole region...
  //   Luna will then complain.   Does this really matter?  Probably not, i.e. as edge segments
  //   get less weight in the Welch averages than the central segments in any case... but seems
  //   awkward to completely skip regions, and so to keep things simple we insist on a fixed, integer
  //   increment value.
  
  //
  // segment 1; p = 0 .. 640
  // segment 2; p = 266 .. 906
  // segment 3; p = 532 .. 1172
  // segment 4; p = 798 .. 1438
  // segment 5; p = 1064 .. 1704
  // segment 6; p = 1330 .. 1970
  // segment 7; p = 1596 .. 2236
  // segment 8; p = 1862 .. 2502
  // segment 9; p = 2128 .. 2768
  // segment 10; p = 2394 .. 3034
  // segment 11; p = 2660 .. 3300
  // segment 12; p = 2926 .. 3566
  // segment 13; p = 3192 .. 3832
  //specified Welch segment parameters:
  //   - segment size    = 640 sample points
  //   - segment overlap = 374 sample points
  //   - implied increment = 266 sample points
  //   - last covered point = 3832 (of 3840)
  // error : Welch segment size/increment does not span epoch fully
	      

  //
  // From MATLAB parameterizatopm:
  //  K = (M-NOVERLAP)/(L-NOVERLAP)
  //    M = total_points = epoch size (in data-points)
  //    L = segment_size_points = segment size (in data-points)
  //    K = noverlap_segments = desired number of (overlapping) segments of size 'L' within 'M'
  //    NOVERLAP = noverlap_points 
  //
  
  int total_points             = data.size();
  int segment_size_points      = M * Fs;   // 'nfft' in Matlab

  // handle special case: if only one segment, which is shorter than the total,
  // then allow for one extra segment (e.g. 5-sec epoch but 4-second window) 

  if ( segment_size_points < total_points && noverlap_segments == 1 )
    ++noverlap_segments;  
  
  int noverlap_points          = noverlap_segments > 1 
    ? ceil( ( noverlap_segments*segment_size_points - total_points  ) / double( noverlap_segments - 1 ) )
    : 0 ;
  
  int segment_increment_points = segment_size_points - noverlap_points;
  
  // std::cout << "segment_size_points = " << segment_size_points << "\n"
  // 	    << "noverlap_points = " << noverlap_points << "\n"
  //   	    << "segment_increment_points = " << segment_increment_points << "\n";
  

  //
  // Check segment coverage: by default, all points must be covered
  //  
    
  int seg_cnt = 1;
  int last_point_plus_one = 0;
  for (int p = 0; p <= total_points - segment_size_points ; p += segment_increment_points )
    {
      //std::cout << "segment " << seg_cnt << "; p = " << p << " .. " << segment_size_points + p << "\n"; 
      last_point_plus_one = p + segment_size_points;
      ++seg_cnt;
    }
    

  if ( last_point_plus_one > total_points )
    //      if ( last_point_plus_one != total_points ) // i.e. allow if slightly shorter
    {
      
      logger << "  specified Welch segment parameters:\n"
	     << "     - segment size    = " << segment_size_points << " sample points\n"
	     << "     - segment overlap = " << noverlap_points << " sample points\n"
	     << "     - implied increment = " << segment_increment_points << " sample points\n"
	     << "     - last covered point = " << last_point_plus_one << " (of " << total_points << ")\n";
      
      logger << " implied segments (in sample points): nb: overlap/increment may have been altered to fit\n"
	     << " which is fine - Luna just requires an increment of an *integer* number of samples\n"
	     << " (for a fixed total signal length, number of segments and segment length) can span the\n"
	     << " whole region\n";
      
      int seg_cnt = 1;
      int last_point_plus_one = 0;
      for (int p = 0; p <= total_points - segment_size_points ; p += segment_increment_points )
	{
	  logger << "segment " << seg_cnt << "; p = " << p << " .. " << segment_size_points + p << "\n"; 
	  last_point_plus_one = p + segment_size_points;
	  ++seg_cnt;
	}
      
      Helper::halt( "Welch segment size/increment does not span epoch fully" );
      
    }


  
  //
  // Initial FFT
  //
  
  real_FFT fft0( segment_size_points , use_nextpow2 ? MiscMath::nextpow2( segment_size_points ) : segment_size_points , Fs , window );
      
  if ( average_adj ) 
    fft0.average_adjacent();

  psd.resize( fft0.cutoff , 0 );
  N = fft0.cutoff;
  
  freq.resize(N);
  for (int f=0;f<N;f++)
    freq[f] = fft0.frq[f];
  
  //
  // Median/SD 
  //
  
  std::vector<std::vector<double> > tracker;
  std::vector<std::vector<double> > lntracker;

  if ( calc_seg_sd || use_median )
    {
      int segments = 0;
      for (int p = 0; p <= total_points - segment_size_points ; p += segment_increment_points )
	++segments;
      
      // freq x segment
      if ( use_median ) 
	{
	  tracker.resize(N);
	  for (int i=0; i<N; i++)
	    tracker[i].resize( segments );
	}
      
      // for CV, tracker natural log of segment power
      if ( calc_seg_sd )
	{
	  lntracker.resize(N);
          for (int i=0; i<N; i++)
            lntracker[i].resize( segments );
	  psdsd.resize( N );
	}
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
      
      
      //      FFT fft( nfft , Fs , FFT_FORWARD , window );
      
      if ( p + segment_size_points > data.size() ) 
	Helper::halt( "internal error in pwelch()" );
      
      const bool detrend = false;

      const bool zerocentre = false;

      if ( detrend )
	{
	  std::vector<double> y( segment_size_points );
	  for (int j=0;j<segment_size_points;j++) y[j] = data[p+j];
	  MiscMath::detrend(&y);      
	  fft0.apply( y ); //wass fft.apply()
	}
      else if ( zerocentre )
	{
	  std::vector<double> y( segment_size_points );
	  for (int j=0;j<segment_size_points;j++) y[j] = data[p+j];
	  MiscMath::centre(&y);      
	  fft0.apply( y ); //wass fft.apply()
	}
      else
	{
	  fft0.apply( &(data[p]) , segment_size_points ); //wass fft.apply()
	}
      
      if ( average_adj )
	fft0.average_adjacent(); //wass fft.apply()
      
      int cutoff = fft0.cutoff;
      
      for (int i=0;i<fft0.cutoff;i++)
	psd[i] += fft0.X[i];

      if ( use_median ) 
	{
	  for (int i=0;i<fft0.cutoff;i++)
	    tracker[i][segments-1] = fft0.X[i];
	}

      if ( calc_seg_sd )
	{
	  for (int i=0;i<fft0.cutoff;i++)
            lntracker[i][segments-1] = log( fft0.X[i] ); // natural log
	}

      
    } // next segment
  

  //
  // take average (mean or median) over segments
  //

  for (int i=0;i<psd.size();i++)
    {
      const double mn = psd[i] / (double)segments;
      
      if ( calc_seg_sd )
	{	  
	  // sd of natural log scaled
	  const double sd = MiscMath::sdev( lntracker[i] );
	  
	  // CV, using formula for log-normal data	  
	  psdsd[i] = sqrt( exp( sd * sd ) -1 );
	  // psdsd[i] = mn > 0 ? sd / mn : 0 ;
	}

      
      if ( use_median )
	{
	  // true means get avg. of lower & upper medians for even lists 
	  // slower, (i.e. need to process list twice) but will be more comparable
	  // to other implementations
	  psd[i] = MiscMath::median( tracker[i] , true );	  
	}
      else	
	psd[i] = mn;
      
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
      
      real_FFT fft( np , np , Fs , WINDOW_HANN );
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
  nominal.clear();

  if ( f.size() < 2 ) return 0;

  
  //
  // no binning?
  //
  
  if ( fac == 1 )
    {
      for (int i=0;i<f.size();i++)
        {
	  if ( f[i] < mn_f ) continue;
          if ( f[i] > mx_f ) break;
          bfa.push_back( f[i] );
          bfb.push_back( f[i] );
          bspec.push_back( y[i] );
	  //nominal.push_back( Helper::dbl2str( f[i] ) );
	  nominal.push_back( "" );
        }
      return bspec.size();
    }

  //
  // integer binning, using lower value as the nominal seed
  //

  int i = 0;
  for (i=0;i<f.size();i++)
    if ( f[i] >= mn_f ) break;

  // put DC always as a separate value

  if ( mn_f == 0 ) 
    {
      bspec.push_back( y[0] );
      bfa.push_back( 0 );
      bfb.push_back( 0 );
      nominal.push_back( "0" );
      ++i;
    }

  // make bins

  for ( ;i<f.size(); i+=fac )
    {
      if ( ( i + fac - 1 ) < f.size() )
	{
	  if ( f[i+fac-1] > mx_f ) break;
	  
	  double sum = 0;
	  for (int j=i;j<i+fac;j++)
	    sum += y[i];
	  bspec.push_back( sum / (double)fac );
	  
	  bfa.push_back( f[i] );
	  bfb.push_back( f[i + fac - 1] );
	  // midpoint as nominal label
	  nominal.push_back( Helper::dbl2str( f[i] ) + "-" + Helper::dbl2str( f[i + fac - 1] ) );
	  //nominal.push_back( ( f[i + fac - 1]  + f[i] ) / 2.0  ); // use mean as nominal label
	}
    }
  
  return bspec.size();
  
  //
  // Old version below... ignore
  //

  
  // assume always from 0, DC component  
  // if ( f[0] == 0 ) 
  //   {
  //     bspec.push_back( y[0] );
  //     bfa.push_back( 0 );
  //     bfb.push_back( 0 );
  //     nominal.push_back( 0 );
  //   }

  //for (int ii=0;ii<f.size();ii++) std::cout << "Frq " << f[ii] << "\n";
  


  // double nyquist = 0.5 * Fs;
    
  // int num_freqs = f.size();
  
  // double df = f[1] - f[0];    

  // if ( w/df  < 1.0 ) Helper::halt( "bin resolution too small: min " + Helper::dbl2str( df ) );

  // int freqwin = (int) ( w / df ) ;      
  
  // if ( mx_f > nyquist ) mx_f = nyquist; 

  // double lower_f = 0;

  // for (int i = 0; i < num_freqs ; i += freqwin)
  //     {
	
  // 	double tem = 0.0;
	
  // 	int k = 0;
	
  // 	for (int j = i ; j < i + freqwin ; j++) 
  // 	  {
	    
  // 	    if (j > 0 && j < num_freqs - 1) // skip DC and Nyquist
  // 	      {	      	      
  // 		if ( f[j] >= mn_f && f[j] <= mx_f )
  // 		  {
  // 		    //std::cout << "adding " << f[j] << "\n";
  // 		    tem += y[j];
  // 		    k++;
  // 		  }
  // 	      }
  // 	  }
	
  // 	std::cout << "scanning " << f[i] << " to " << f[ i + freqwin -1 ] <<  " " << k << "\n";	
	
  // 	if ( k > 0 ) 
  // 	  {	  
  // 	    bspec.push_back( tem/(double)k );
  // 	    bfa.push_back( f[i-1] ); // less than 
  // 	    bfb.push_back( f[i+k-1] ); // greater than or equal to
  // 	    nominal.push_back( lower_f + w * 0.5 ); // intended midpoint 
  // 	    std::cout << "ADDING: k=" << k << " points in " <<  f[i-1]  << " < F <= " << f[i+k-1] << "\n";
  // 	  }

  // 	lower_f += w;
	
  //     }
  
  //   return bspec.size();

}




void psd_shape_metrics( const std::vector<double> & f , // frq
			const std::vector<double> & x , // log(power)
			const int w , 
			double * m1 ,
			double * m2 ,
			std::vector<double> * detrended ,
			std::vector<double> * smoothed , 
			std::vector<double> * difference )
{

  const int n = f.size();
  if ( x.size() != n )
    {
      std::cerr << f.size() << "\t" << x.size() << "\n";
      Helper::halt( "f and x of different sizes" );
    }

  // assume p = log power
  
  // put on 0..1 scale

  double xmin, xmax;
  MiscMath::minmax( x , &xmin, &xmax );

  std::vector<double> p(n);
  for (int i=0; i<n; i++)
    p[i] = ( x[i] - xmin ) / ( xmax - xmin );

  // overall linear detrend (nb. lin/log scale)
  // and Z-normalize

  double pa, pb;
  p = MiscMath::edge_detrend( p , &pa, &pb ) ;

  // get smoothed: ss = smoothed, pp = p - ss 
  std::vector<double> ss; 
  std::vector<double> pp = MiscMath::remove_median_filter( p , w , &ss );  
    
  double ppmin, ppmax;
  MiscMath::minmax( pp , &ppmin, &ppmax );
  
  // sum( abs(diff(pp)) ) 
  // m1 --> SPK

  *m1 = 0; 
  for (int i=1; i<n; i++)
    *m1 += abs( pp[i] - pp[i-1] ) ; 
  
  // peak-to-peak of difference
  if ( 0 ) 
    {
      *m1 = ( ppmax - ppmin ) ; 
    }  

  // kurtosis of delta (i.e. assume mean = 0 )
  // m2 --> KURT
  double numer = 0 , denom = 0;
  for (int i=0; i<n; i++)
    {
      numer += pow( pp[i] , 4 );
      denom += pow( pp[i] , 2 );
    }

  numer /= (double)n;
  denom /= (double)n;
  denom *= denom;
  *m2 = numer / denom - 3.0;
    
  // verbose reports?
  if ( detrended != NULL ) *detrended = p;
  if ( smoothed != NULL ) *smoothed = ss;
  if ( difference != NULL ) *difference = pp;
   
}





bool spectral_slope_helper( const std::vector<double> & psd , 
			    const std::vector<double> & freq , 
			    const std::vector<double> & fr , 
			    const double outlier , 
			    const bool display , 
			    double * b , double * bn , double * bi , double * rsq )
{
  
  std::vector<double> slope_y, slope_x;

  for (int f=0; f<psd.size(); f++)
    {      
      if ( freq[f] < fr[0]  ) continue;
      if ( freq[f] > fr[1] ) break;

      slope_x.push_back( log( freq[f] ) );
      
      if ( psd[f] <= 0 ) Helper::halt( "negative/zero PSD in spectral slope estimation" );
      slope_y.push_back( log( psd[f] ) );
      
    }
  
  const int n = slope_y.size();

  // 
  // remove outliers (in log-PSD space after fitting initial regression line)?
  //

  if ( outlier > 0 ) 
    {
      // first detrend, i.e. fit initial slope
      // n.b. assumes uniform spacing in frequencies
      
      std::vector<double> dt_y = MiscMath::detrend( slope_y );

      double mean_y = MiscMath::mean( dt_y );
      double sd_y = MiscMath::sdev( dt_y , mean_y );
      double lwr = mean_y - outlier * sd_y;
      double upr = mean_y + outlier * sd_y;

      std::vector<bool> exc( n );

      bool remove = false;

      for (int i=0; i<n; i++)
	{
	  exc[i] = dt_y[i] < lwr || dt_y[i] > upr  ;
	  if ( exc[i] ) remove = true; 	  
	}
      
      if ( remove )
	{
	  std::vector<double> cp_y = slope_y;
	  std::vector<double> cp_x = slope_x;
	  slope_x.clear(); 
	  slope_y.clear();
	  
	  for (int i=0; i<n; i++)
	    {
	      if ( ! exc[i] ) 
		{
		  slope_y.push_back( cp_y[i] );
		  slope_x.push_back( cp_x[i] );		  
		}
	    }	  
	  
	}
    }

  // not enough data?
  if ( slope_y.size() < 3 ) return false ;
  
  dynam_t spec_slope( slope_y , slope_x );
  double beta_spec_slope, beta_rsq, intercept;
  spec_slope.linear_trend( &beta_spec_slope , &beta_rsq , &intercept );

  if ( display ) 
    {
      writer.value( "SPEC_SLOPE" , beta_spec_slope );
      writer.value( "SPEC_INTERCEPT" , intercept );
      writer.value( "SPEC_RSQ" , beta_rsq );
      writer.value( "SPEC_SLOPE_N" , (int)slope_y.size() );
    }
  
  if ( b != NULL ) *b = beta_spec_slope;
  if ( bn != NULL ) *bn = (int)slope_y.size();
  if ( bi != NULL ) *bi = intercept;
  if ( rsq != NULL ) *rsq = beta_rsq;
  return true;
}





 void peakedness( const std::vector<double> & p , 
		  const std::vector<double> & f0 , 
		  const int peak_median_filter_n , 
		  const std::vector<double> & peak_range , 
		  const bool verbose )
 {
   
   double m1, m2;
   
   // detrended / smoothed / difference
   std::vector<double> shape1, shape2, shape3;
   
   const int n = p.size();
   
   std::vector<double> frq;
   std::vector<double> logged;
   
   for (int i=0; i<n; i++)
     {
       if ( f0[i] >= peak_range[0] && f0[i] <= peak_range[1] ) 
	 {
	   frq.push_back( f0[i] );
	   logged.push_back( 10*log10( p[i] ) );
	   
	 }
     }

   // not a wide enough region for median smoothing even?
   if ( frq.size() < peak_median_filter_n * 1.5 ) return;

   psd_shape_metrics( frq ,
		      logged , 
		      peak_median_filter_n , 
		      &m1, &m2 ,
		      &shape1, &shape2, &shape3);
   
   writer.value( "SPK" , m1 );
   writer.value( "KURT" , m2 );
   
   if ( verbose ) 
     {
       const int n2 = frq.size();
       for (int i=0; i<n2; i++)
	 {
	   writer.level( frq[i] , globals::freq_strat );
	   writer.value( "DT" , shape1[i] );
	   writer.value( "SM" , shape2[i] );
	   writer.value( "DF" , shape3[i] );		      
	 }
       writer.unlevel( globals::freq_strat );		  
       
     }
 }



