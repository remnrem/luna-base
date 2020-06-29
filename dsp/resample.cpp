
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


#include "resample.h"
#include <iostream>

#include "libsamplerate/samplerate.h"

#include "eval.h"
#include "edf/edf.h"
#include "edf/slice.h"

#include "helper/helper.h"
#include "helper/logger.h"

extern logger_t logger;

std::vector<double> dsptools::resample( const std::vector<double> * d , 
					int sr1 , int sr2 )
{


  int n = d->size();
  std::vector<float> f( n );
  for (int i=0;i<n;i++) f[i] = (*d)[i];
  
  double ratio = sr2 / (double)sr1;
  const int n2 = n * ratio;
  std::vector<float> f2( n2 );

  // pad a little at end (probably not necessary)
  for (int i=0;i<10;i++) 
    {
      ++n;
      f.push_back(0);
    }

  SRC_DATA src;
  src.data_in = &(f[0]);
  src.input_frames = n;
  src.data_out = &(f2[0]);
  src.output_frames = n2;
  src.src_ratio = ratio;

  int r = src_simple( &src, SRC_SINC_FASTEST , 1 );
  
  // problem?
  if ( r ) 
    {
      logger << src_strerror ( r ) << "\n";
      Helper::halt( "problem in resample()" );
    }

  std::vector<double> out( n2 );
  for (int i=0;i<n2;i++) out[i] = f2[i];

  return out;
}


void dsptools::resample_channel( edf_t & edf , const int s , const int nsr )
{
  
  // s is in 0..ns space  (not 0..ns_all)  

  if ( edf.header.is_annotation_channel(s) ) return;
    
  //
  // Original sample rate 
  //
  
  const int Fs = edf.header.sampling_freq(  s ) ;
  
  // already done?
  if ( Fs == nsr ) return; 


  //
  // Ouput
  //

  logger << "  resampling channel " << edf.header.label[ s ] << " from sample rate " << Fs << " to " << nsr << "\n";
  
  //
  // Pull entire signals out
  //
  
  interval_t interval = edf.timeline.wholetrace();

  slice_t slice( edf , s , interval );

  const std::vector<double> * d = slice.pdata();
  
  //
  // Resample to new SR
  //

  std::vector<double> resampled = resample( d , Fs , nsr );
  
  // 
  // Ensure that resultant signal is the exact correct length 
  //
  
  resampled.resize( edf.header.nr * edf.header.record_duration * nsr , 0 ); // i.e. zero-pad if necessary


  //
  // Update EDF header with new sampling rate
  //
  
  // note: the EDF header also contains n_samples_all[], which is mapped against the contents of the EDF, 
  // rather than the set of selected signals.   This should stay as is, in any case, i.e. it is only used 
  // to know how to skip signals when reading the EDF, i.e. and this won't change.
  
  edf.header.n_samples[ s ] = nsr * edf.header.record_duration ;
  
  //
  // Place back
  //

  edf.update_signal( s , &resampled );

}


void dsptools::resample_channel( edf_t & edf , param_t & param )
{
  
  std::string signal_label = param.requires( "sig" );
  signal_list_t signals = edf.header.signal_list( signal_label );      
  std::vector<double> Fs = edf.header.sampling_freq( signals );

  // new sampling rate for all channels
  int sr = param.requires_int("sr");

  const int ns = signals.size();
  
  for (int s=0;s<ns;s++)
    resample_channel( edf , signals(s) , sr );

}
