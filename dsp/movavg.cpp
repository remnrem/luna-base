
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

#include "movavg.h"

#include "stats/eigen_ops.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include <vector>

void dsptools::movavg( edf_t & edf , param_t & param )
{

  const bool median = param.has( "median" );
  const bool triangular = param.has( "tri" );
  const double hwin_sec = param.requires_dbl( "hw" );
  const double tri_lwr = param.has( "lwr" ) ? param.requires_dbl( "lwr" ) : 0 ;
   
  // signal(s)
  const bool no_annots = true; 
  signal_list_t signals = edf.header.signal_list( param.value( "sig" ) , no_annots );
  
  const int ns = signals.size();

  if ( ns == 0 ) return;

  // do by epoch?
  bool by_epoch = param.has( "epoch" );
  if ( by_epoch ) edf.timeline.ensure_epoched();

  if ( by_epoch ) logger << "  iterating over epochs\n";
  logger << "  applying moving average (hwin = " << hwin_sec << ") :";
  
  // process each signal
  for (int s=0; s<ns; s++)
    {
      
      // get sample rate
      int sr = edf.header.sampling_freq( signals(s) );
      
      int hwin = hwin_sec * sr;
      logger <<" hwin = " << hwin << " " << hwin_sec << " " << sr << "\n";
      if ( hwin == 0 )
	{
	  logger << "  skipping " << signals.label(s) << ", sample rate too low\n";
	  continue;
	}

      // make a full window
      hwin = 1 + 2 * hwin;
      
      // get whole signal  (although we update epoch-by-epoch)
      // we need to this edit and return back 

      slice_t slice0( edf , signals(s) , edf.timeline.wholetrace() );
      std::vector<double> orig = * slice0.pdata();
  	
      // iterate over each epoch
      int ne = by_epoch ? edf.timeline.first_epoch() : 1 ;
      
      // compile new signal (by epochs)
      std::vector<Eigen::VectorXd> avg;
      
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

	  Eigen::ArrayXd dt = eigen_ops::copy_array( *data );

	  Eigen::VectorXd rt;
	  
	  if ( median ) 
	    rt = eigen_ops::median_filter( dt , hwin );
	  else if ( triangular )
	    rt = eigen_ops::tri_moving_average( dt , hwin , tri_lwr );
	  else
	    rt = eigen_ops::moving_average( dt , hwin );
	  
	  // do procedure, and store
	  avg.push_back( rt );
	  
	  // done?
	  if ( ! by_epoch ) break;
	  
	  // next epoch
	}
      
      // update signal
      int p = 0;
      const int n = orig.size();
      for (int e=0; e<avg.size(); e++)
	for (int i=0; i<avg[e].size(); i++)
	  orig[p++] = avg[e][i];
      
      logger << " " << signals.label(s);

      // and send back to the EDF
      edf.update_signal( signals(s) , &orig );      
      
    }
  logger << "\n";
  
}


