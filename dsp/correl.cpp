
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


#include "correl.h"

#include "miscmath/miscmath.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "eval.h"

#include "stats/statistics.h"

#include "db/db.h"
#include "dsp/resample.h"

#include "helper/helper.h"
#include "helper/logger.h"

#include <string>


extern logger_t logger;

extern writer_t writer;
   
void dsptools::correlate_channels( edf_t & edf , param_t & param )
{
  
  std::string signal_label = param.requires( "signal" );
  
  signal_list_t signals = edf.header.signal_list( signal_label );  
  
  const int ns = signals.size();
  
  int sr = param.has( "sr" ) ? param.requires_int( "sr" ) : 0 ;
  
  writer.var( "R" , "Channel correlation (-1..1)" );

  //
  // adjust all SRs now if needed
  //

  if ( sr )
    {
      for (int s=0;s<ns;s++)
	{
	  
	  if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
	  
	  if ( edf.header.sampling_freq( signals(s) ) != sr ) 
	    {
	      logger << "resampling channel " << signals.label(s) 
		     << " from " << edf.header.sampling_freq( signals(s) )
		     << " to " << sr << "\n";
	      resample_channel( edf, signals(s) , sr );
	    }
	}
    }
  else
    {
      std::set<int> srs;
      for (int s=0;s<ns;s++)
	{
	  if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
	  srs.insert( edf.header.sampling_freq( signals(s) ) );
	}
      if ( srs.size() > 1 ) Helper::halt( "all sampling rates must be similar, use 'sr'" );
    }

  
  //
  // Epochs or whole signal?
  //
  
  bool epoched = edf.timeline.epoched() && param.has("epoch") ;

  bool verbose = param.has("verbose");

  //
  // Ensure similar sampling rates
  //
  
  for (int i=0;i<ns-1;i++)
    {
      
      if ( edf.header.is_annotation_channel( signals(i) ) ) continue;
      
      for (int j=i+1;j<ns;j++)
	{
	  
	  if ( edf.header.is_annotation_channel( signals(j) ) ) continue;
	  
	  //
	  // Store correlations
	  //
	  
	  std::vector<double> epoch_r;
	  double overall_r;
	  double mean_r;
	  double median_r;

	  // stratify output by SIGNALS
	  writer.level( signals.label(i) + "x" + signals.label(j) , "CHS" );
	  
	  
	  if ( epoched ) 
	    {
	      
	      int epoch = edf.timeline.first_epoch();      
	      
	      while ( 1 ) 
		{

		  int epoch = edf.timeline.next_epoch();      
		  if ( epoch == -1 ) break;

		  interval_t interval = edf.timeline.epoch( epoch );

 		  slice_t slice1( edf , signals(i) , interval );
 		  slice_t slice2( edf , signals(j) , interval );
		  
 		  const std::vector<double> * d1 = slice1.pdata();
 		  const std::vector<double> * d2 = slice2.pdata();

		  double r = Statistics::correlation( *d1 , *d2 );

		  //
		  // Output
		  //

		  if ( verbose )
		    {
		      writer.epoch( edf.timeline.display_epoch( epoch ) );
		      
		      writer.value( "R" , r );
		    }

		  epoch_r.push_back( r );


		} // next epoch
	      	      
	      if ( verbose )
		writer.unepoch();
	      
	      //
	      // Get mean/median correlation over epochs
	      //
	    
	      mean_r = MiscMath::mean( epoch_r );
	      median_r = MiscMath::median( epoch_r );
	      
	      writer.value( "R_MEAN" , mean_r ); 
	      writer.value( "R_MEDIAN" , median_r ); 

	    }
	  else
	    {
	      
	      //
	      // Coherence for entire signal
	      //

	      interval_t interval = edf.timeline.wholetrace();
	      
	      slice_t slice1( edf , signals(i) , interval );
	      slice_t slice2( edf , signals(j) , interval );
	      
	      const std::vector<double> * d1 = slice1.pdata();
	      const std::vector<double> * d2 = slice2.pdata();
	      
	      overall_r = Statistics::correlation( *d1 , *d2 );
	      	      
	      writer.value( "R" , overall_r ); 
	      	      
	    }	  
	}
    }
  
  writer.unlevel( "CHS" );

}

