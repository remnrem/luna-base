
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


#include "cwt-design.h"
#include "cwt/cwt.h"
#include "fftw/fftwrap.h"
#include "db/db.h"

extern writer_t writer;

void dsptools::design_cwt( param_t & param )
{
  
  double fc = param.requires_dbl( "fc" );

  int cycles = param.requires_int( "cycles" );

  int fs = param.requires_int( "fs" );

  logger << " running CWT design for fc=" << fc << ", cycles=" << cycles << " and fs=" << fs << "\n";

  writer.cmd( "CWT-DESIGN" , 1 , param.dump( "" , " " ) );

  writer.level( Helper::dbl2str( fc ) + "_" + Helper::int2str( cycles ) + "_" + Helper::int2str( fs ) , "PARAM" );
  
  CWT cwt;
  
  cwt.set_sampling_rate( fs );
  
  cwt.add_wavelet( fc , cycles );
  
  cwt.set_timeframe( fc );
  
  std::vector<dcomp> w = cwt.wavelet(0);
  std::vector<double> t = cwt.get_timeframe();
  
  const int n = w.size();
  
  writer.numeric_factor( "SEC" );

  for (int i=0;i<n;i++)
    {
      writer.level( t[i] , "SEC" );
      writer.value( "REAL" , std::real( w[i] ) );
      writer.value( "IMAG" , std::imag( w[i] ) );
    }
  writer.unlevel( "SEC" );

  //
  // FFT of CWT
  //
  

  FFT fft( n , fs );
  fft.apply( w );
  std::vector<double> mag = fft.mag;
  const int nx= mag.size();

  // standardize to peak == 1.0
  double mx = 0;
  for (int i=0;i<nx;i++) if ( mag[i] > mx ) mx = mag[i];
  if ( mx > 0 ) for (int i=0;i<nx;i++) mag[i] /= mx;
  
  
  for (int i=0;i<fft.frq.size();i++)
    {
      writer.level( fft.frq[i] , globals::freq_strat );
      writer.value( "MAG" , mag[i] );
    }
  writer.unlevel( globals::freq_strat );

  writer.unlevel( "PARAM" );

}
