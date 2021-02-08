
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

#include "dsp/siggen.h"
#include "edf/edf.h"
#include "edf/slice.h"

#include <vector>

#include <cmath>

void dsptools::siggen( edf_t & edf , param_t & param )
{

  //
  // Get signals, check sampling rates
  //
  
  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) );

  if ( signals.size() == 0 ) return;
  
  const int ns = signals.size();
  
  std::vector<double> Fs = edf.header.sampling_freq( signals );
  
  int sr = Fs[0];
  for (int s=1;s<ns;s++)
    if ( Fs[s] != sr )
      Helper::halt( "all sampling rates must be similar for SIGGEN" );


  //
  // Options
  //

  // sine, square, saw, triangular

  bool sine_wave = param.has( "sine" );
  std::vector<double> sine_param;

  if ( sine_wave )
    {
      sine_param = param.dblvector( "sine" );
      if ( sine_param.size() == 2 ) sine_param.resize(3,0);
      else if ( sine_param.size() != 3 ) Helper::halt( "expecting sine=frq,amp{,phase}" );
      if ( sine_param[0] <= 0 ) Helper::halt( "frq must be positive" );
      if ( sine_param[0] >= sr / 2.0 ) Helper::halt( "frq not under Nyquist frequency, given sample rate" );
      if ( sine_param[1] <= 0 ) Helper::halt( "amp should be positive, non-zero" );
    }

  bool first_clear = param.has( "clear" );
  
  
  //
  // pull signals
  //

  for (int s=0; s<ns; s++)
    {

      
      interval_t interval = edf.timeline.wholetrace();

      slice_t slice( edf , signals(s) , interval );
      
      std::vector<double> d = * slice.pdata();

      const std::vector<uint64_t> * tp = slice.ptimepoints();

      const int np = tp->size();
      
      for ( int p=0 ; p<np; p++ )
	{
	  // time in seconds
	  double t = (*tp)[p] * globals::tp_duration;

	  // add to original value, or start from scratch?
	  double x = first_clear ? 0 : d[p];

	  // add a sine wave?
	  if ( sine_wave )
	    x += sine_param[1] * sin( 2 * M_PI * sine_param[0] * t + sine_param[2] );

	  // replace back
	  d[p] = x;
	}      

      //
      // copy whole signal back
      //
      
      edf.update_signal( signals(s) , &d );

    }
  

  //
  // all done
  //
  
  
}

