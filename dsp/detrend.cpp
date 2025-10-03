
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

#include "dsp/detrend.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "param.h"

void dsptools::detrend( edf_t & edf , param_t & param )
{
  
  // simple per-epoch mean removal

  // signal(s)
  const bool no_annots = true; 
  signal_list_t signals = edf.header.signal_list( param.value( "sig" ) , no_annots );

  const int ns = signals.size();

  if ( ns == 0 ) return;

  // do by epoch?
  bool by_epoch = param.has( "epoch" );
  if ( by_epoch ) edf.timeline.ensure_epoched();

  if ( by_epoch ) logger << "  iterating over epochs\n";
  else logger << "  correcting for entire signal\n";
  logger << "  removing signal mean:";
  
  // process each signal
  for (int s=0; s<ns; s++)
    {
      
      // get sample rate
      int sr = edf.header.sampling_freq( signals(s) );
      
      // get whole signal  (although we update epoch-by-epoch)
      // we need to this edit and return back 

      slice_t slice0( edf , signals(s) , edf.timeline.wholetrace() );
      std::vector<double> orig = * slice0.pdata();
  	
      // iterate over each epoch
      int ne = by_epoch ? edf.timeline.first_epoch() : 1 ;
      
      // compile new signal (by epochs)
      std::vector<std::vector<double> > filt;
      
      while ( 1 )
	{
	  
	  // next epoch
	  int epoch = by_epoch ? edf.timeline.next_epoch() : 1 ; 

	  // all done?
	  if ( epoch == -1 ) break;
	  
	  // get data
	  interval_t interval = by_epoch ? edf.timeline.epoch( epoch ) : edf.timeline.wholetrace();

	  slice_t slice( edf , signals(s) , interval );

	  const std::vector<double> * data = slice.pdata();
	  
	  double mn = MiscMath::mean( *data );
	  
	  std::vector<double> adj( data->size() );
	  for (int i=0; i<adj.size(); i++) adj[i] = (*data)[i] - mn; 

	  // do procedure, and store
	  filt.push_back( adj );

	  // done?
	  if ( ! by_epoch ) break;

	  // next epoch
	}
      
      // update signal
      int p = 0;
      const int n = orig.size();
      for (int e=0; e<filt.size(); e++)
	for (int i=0; i<filt[e].size(); i++)
	  orig[p++] = filt[e][i];
      
      logger << " " << signals.label(s);

      // and send back to the EDF
      edf.update_signal( signals(s) , &orig );      
      
    }
  logger << "\n";
			    
  
}

