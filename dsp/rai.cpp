
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

#include "dsp/rai.h"

#include "param.h"
#include "edf/edf.h"
#include "edf/slice.h"

#include "db/db.h"
#include "helper/helper.h"
#include "helper/logger.h"

extern logger_t logger;
extern writer_t writer;

void dsptools::rai( edf_t & edf , param_t & param )
{

  //
  // get signals
  //
  
  std::string signal_label = param.requires( "sig" );
  const bool no_annotations = true;
  signal_list_t signals = edf.header.signal_list( signal_label , no_annotations );
  const int ns = signals.size();

  //
  // parameters
  //

  // lower uV threshold
  const double th = param.has( "th" ) ? param.requires_dbl( "th" ) : 1 ; 

  // exclusion region
  const double th2 = param.has( "th2" ) ? param.requires_dbl( "th2" ) : 2 ; 

  // i.e. count 0-1 and compare to rest (excluding 1-2)

  const bool verbose = param.has( "verbose" );
  
  //
  // assumptions: 
  //

  //  chin-EMG
  //  restricted to REM epochs
  //  band-pass filtered (10 - 100 Hz); possible notch at 50/60 Hz
  //  now 1-second epochs
  //  uV scaling


  if ( param.has( "epoch" ) && ! edf.timeline.epoched() )
    Helper::halt( "no EPOCHs set" );

  // use 30 seconds if not otherwise specified                                                                                                 
  if ( ! edf.timeline.epoched() )
    Helper::halt( "no EPOCHs set" );
  
  if ( fabs( edf.timeline.epoch_length() - 1.0 ) > 0.0001 ) 
    Helper::halt( "require 1second epochs" );
  
  //  -->  for each 1s epoch; get mean of rectified signal
       
  for (int s=0; s<ns; s++)
    {

      writer.level( signals.label(s) , globals::signal_strat );
		   
      std::vector<double> m;

      const int ne = edf.timeline.first_epoch();
      
      while ( 1 ) 
	{
	  
	  int epoch = edf.timeline.next_epoch();
	  
	  if ( epoch == -1 ) break;
	  
	  interval_t interval = edf.timeline.epoch( epoch );
	  
	  slice_t slice( edf , signals(s) , interval );
	  
	  std::vector<double> * d = slice.nonconst_pdata();

	  const int n = d->size();

	  // avg of the recified signal
	  double x = 0; 
	  for (int i=0; i<n; i++)
	    x += fabs( (*d)[i] ) ;
	  m.push_back( x / (double)n ) ;
	  
	} // next epoch
      
      // baseline correction: subtract minimum in average 60s window
      
      std::vector<double> mins = MiscMath::moving_min( m , 60 + 1 ) ;
      
      double u = 0 , v = 0 ;
      
      for (int i=0; i<mins.size(); i++)
	{
	  const double x = m[i] - mins[i];
	  if ( x < th ) ++u;
	  else if ( x > th2 ) ++v;

	  if ( verbose )
	    {
	      writer.level( i+1 , globals::count_strat );
	      writer.value( "X", x );
	    }	  
	}

      if ( verbose )
	writer.unlevel( globals::count_strat );

      // main output 

      writer.value( "REM_AI" , u / (double)(u+v) );
      writer.value( "NE" , u+v );
      
    } // next signal

  writer.unlevel( globals::signal_strat );
  
}


