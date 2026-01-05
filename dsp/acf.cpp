
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


#include "dsp/acf.h"
#include "param.h"
#include "fftw/fftwrap.h"
#include "miscmath/miscmath.h"
#include "db/db.h"
#include "edf/edf.h"
#include "edf/slice.h"

extern logger_t logger;
extern writer_t writer;

void dsptools::autocorr_channels( edf_t & edf , param_t & param )
{

  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) );

  const int maxlag = param.requires_int( "lag" );
  
  interval_t interval = edf.timeline.wholetrace();

  const int ns = signals.size();  
  
  for (int s=0;s<ns;s++)
    {
      
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;

      const double Fs = edf.header.sampling_freq( signals(s) );
      
      writer.level( signals.label(s) , globals::signal_strat );

      logger << "  estimating ACF for " << signals.label(s) << " (up to " << maxlag*(1.0/Fs) << " seconds)\n" ;

      slice_t slice( edf , signals(s) , interval );
      
      const std::vector<double> * d = slice.pdata();

      acf_t acf( *d , maxlag );
      std::vector<double> r = acf.acf();

      for (int lag=1;lag<r.size();lag++)
	{
	  writer.level( lag , "LAG" );
	  writer.value(  "ACF" , r[lag]  );
	  writer.value(  "SEC" , lag* (1.0/Fs) );
	}
      writer.unlevel( "LAG" );

    }
  writer.unlevel( globals::signal_strat );
  
}



// using FFT to compute ACF, e.g.
// https://dsp.stackexchange.com/questions/1919/efficiently-calculating-autocorrelation-using-ffts/

#include <cmath>
#include <limits>

void acf_t::calc(const std::vector<double> & d, int maxlag)
{
  const int N = (int)d.size();
  if (N == 0) { r.clear(); return; }

  if (maxlag == 0) maxlag = (int)std::llround(10.0 * std::log10((double)N));
  if (maxlag < 0) maxlag = 0;

  std::vector<double> dz = d;

  const double mean = MiscMath::mean(d);
  for (int i = 0; i < N; ++i) dz[i] -= mean;

  dz.resize(2 * (size_t)N, 0.0);
  const int nfft = (int)dz.size();

  FFT fft(nfft, nfft, 1, FFT_FORWARD, WINDOW_NONE);
  fft.apply(dz);
  std::vector<dcomp> f = fft.transform();
  for (size_t i = 0; i < f.size(); ++i) f[i] = f[i] * conj(f[i]);

  FFT ifft(nfft, nfft, 1, FFT_INVERSE);
  ifft.apply(f);
  std::vector<double> xx = ifft.inverse();

  if (xx.empty() || xx[0] == 0.0 || std::isnan(xx[0])) {
    r.assign(1, std::numeric_limits<double>::quiet_NaN());
    return;
  }

  int max_avail = (int)xx.size() - 1;
  if (max_avail < 0) max_avail = 0;
  if (maxlag > max_avail) maxlag = max_avail;

  r.resize(1 + maxlag);
  for (int i = 0; i <= maxlag; ++i)
    r[i] = xx[i] / xx[0];
}

// original version
// void acf_t::calc( const std::vector<double> & d , int maxlag )
// {
//   const int N = d.size();
//   if ( maxlag == 0 ) maxlag = 10*log10( d.size() );
//   r.resize( maxlag );

//   std::vector<double> dz = d;
//   // remove mean
//   double mean = MiscMath::mean( d );
//   for (int i=0;i<N;i++) dz[i] -= mean;

//   // zero-pad
//   dz.resize( 2 * d.size() , 0 );
//   const int nfft = dz.size();
  
//   // take FFT
//   FFT fft( nfft , nfft , 1 , FFT_FORWARD , WINDOW_NONE );
//   fft.apply( dz );
//   std::vector<dcomp> f = fft.transform();
//   for (int i=0;i<f.size();i++) f[i] = f[i] * conj( f[i] );

//   // IFFT
//   FFT ifft( nfft , nfft , 1 , FFT_INVERSE );
//   ifft.apply( f );
//   std::vector<double> xx = ifft.inverse();

//   r.resize( 1 + maxlag );
//   for (int i=0;i <= maxlag;i++)
//     r[i] = xx[i] / xx[0];

// }

