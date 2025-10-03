
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

#include "dsp/zpeaks.h"

#include "param.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "timeline/cache.h"
#include "db/db.h"
#include "helper/helper.h"
#include "helper/logger.h"

extern logger_t logger;
extern writer_t writer;

//
// Implmentation and minor extension of peak finding heuristic :
//
// Brakel, J.P.G. van (2014). "Robust peak detection algorithm using
// z-scores". Stack Overflow. Available at:
// https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data/22640362#22640362
// (version: 2020-11-08).
//

//  -- primary code in MiscMath::smoothedZ()


// lag: higher = more smoothing / more adaptive to long term average (--> window, 'w')
//      for stationary series, use a higher lag
//      to capture time-varying trends, use a lower lag

// influence: extent to which peaks influence baseline
//  0 = no influence; 1 = complete
//  for stationary series, use low/0 influence
//  higher numbers = better abl to capture quick changes related to spiking

// threshold: number of SD units above moving mean; set based on expected rate
//   i.e. 3.5 --> p = 0.00047 --> 1/p = 1 in 2128


// place in annot and/or cache

void dsptools::zpeaks( edf_t & edf , param_t & param )
{

  //
  // parameters
  //

  // use local peak-finding threshold method ( still uses th/min0 and th2/min, and max                                                  
  const double window_sec = param.requires_dbl( "w" );
  
  const double influence = param.has( "influence" ) ? param.requires_dbl( "influence" ) : 0.01;
  
  if ( influence < 0 || influence > 1 ) Helper::halt( "influence should be between 0 and 1" );

  // core region
  
  const double threshold = param.requires_dbl( "th" );
  
  const double min_dur_sec = param.has( "sec" ) ? param.requires_dbl( "sec" ) : 0 ;

  const double max_threshold = param.has( "max" ) ? param.requires_dbl( "max" ) : 0 ;
    
  // flanking region:

  const double threshold2 = param.has( "th2" ) ? param.requires_dbl( "th2" ) :	0 ;

  const double min_dur2_sec = param.has( "sec2" ) ? param.requires_dbl( "sec2" ) : 0 ;
  
  const bool ignore_negatives = ! param.has( "negatives" ) ;

  //
  // save annotations
  //
  
  const std::string annot = param.has( "annot" ) ? param.value( "annot" ) : "";
  
  const double add_flank_sec = param.has( "add-flanking" ) ? param.requires_dbl( "add-flanking" ) : 0 ;  

  if ( annot != "" ) 
    logger << "  writing peaks to annotation " << annot
	   << " with " << add_flank_sec
	   << " seconds added each side\n" ;
  
  const bool to_cache = param.has( "cache" );
  const std::string cname = to_cache ? param.requires( "cache" ) : "";
  cache_t<int> * cache = to_cache ? edf.timeline.cache.find_int( cname ) : NULL ;
  if ( to_cache ) logger << "  writing peaks to cache " << cname << "\n" ;
  
  //
  // signals to process
  //
  
  std::string signal_label = param.requires( "sig" );
  const bool no_annotations = true;
  signal_list_t signals = edf.header.signal_list( signal_label , no_annotations );
  const int ns = signals.size();
  

  
  //
  // process data 
  //
  
  for (int s=0; s<ns; s++)
    {
      
      slice_t slice( edf , signals(s) , edf.timeline.wholetrace() );
      
      std::vector<double> * d = slice.nonconst_pdata();
      
      const std::vector<uint64_t> * tp = slice.ptimepoints();
      
      const int n = d->size();

      const int Fs = edf.header.sampling_freq( signals(s) ) ;
      
      //
      // outputs
      //
      
      writer.level( signals.label(s) , globals::signal_strat );
      std::map<std::string,std::string> faclvl = writer.faclvl();
      
      
      //
      // find peaks
      //

      std::vector<interval_t> peaks;
      
      
      // derive parameters given SR
      const int lag_sp = Fs * window_sec ;
      const int min_dur_sp = min_dur_sec * Fs;
      const int min_dur2_sp = min_dur2_sec * Fs;

      const bool verbose = false;
      
      std::vector<int> mxpks; // for caching - these are the top points within each detected peak (sample-point values)

      std::vector<int> pks = MiscMath::smoothedZ( *d, lag_sp,
						  threshold, influence, 
						  min_dur_sp,
						  max_threshold,
						  threshold2,
						  min_dur2_sp,
						  ignore_negatives, 
						  &peaks,
						  &mxpks, 
						  verbose );
      //
      // report
      //

      const int na = peaks.size();

      if ( mxpks.size() != na ) Helper::halt( "internal prob in zpeaks, from smoothedZ()" );
      
      // ignore if a peak spans a discontinuity

      std::vector<bool> okay( na , true );
      
      uint64_t add_flank_tp = add_flank_sec * globals::tp_1sec;

      double coverage = 0; 

      int valid = 0;

      double total_dur_min = ( edf.header.nr * edf.header.record_duration ) / 60.0;
       
      
      for (int i=0; i<na; i++)
	{	      

	  
	  if ( !edf.timeline.discontinuity( *tp , Fs , peaks[i].start , peaks[i].stop ) )
	    {
	      uint64_t start_tp = (*tp)[peaks[i].start] > add_flank_tp ? (*tp)[peaks[i].start] - add_flank_tp : 0 ;
	      uint64_t stop_tp  = (*tp)[peaks[i].stop] + add_flank_tp ;
	      
	      //	      std::cout << " det " << start_tp << " " << stop_tp << " \t" << peaks[i].start << " " << peaks[i].stop << "\n";
	      
	      coverage += globals::tp_duration * (double)( stop_tp - start_tp );
	      ++valid;	      
	    }
	  else
	    okay[i] = false;
	}
      
      logger << "  detected " << valid << " peaks for "
	     << signals.label(s) 
	     << "( " << valid / total_dur_min << " per minute)"
	     << ", spanning "
	     << coverage << " seconds\n";
      if ( na > valid )
	logger << "   rejected " << na - valid << " peaks that spanned discontinuities\n";
      
      
      
      //
      // save annots
      //

      if ( annot != "" )
	{

	  annot_t * a = edf.annotations->add( annot );
	  
	  const std::string ch = signals.label( s );
	  	  
	  for (int i=0; i<na; i++)
	    {
	      // std::cout << " peaks == " << peaks[i].start << "   " << peaks[i].stop << "\t"
	      // 	    << (*tp)[peaks[i].start] << "   " << (*tp)[peaks[i].stop] << "\n";

	      if ( okay[i] )
		{
		  uint64_t start_tp = (*tp)[peaks[i].start] > add_flank_tp ? (*tp)[peaks[i].start] - add_flank_tp : 0 ;
		  uint64_t stop_tp  = (*tp)[peaks[i].stop] + add_flank_tp ; 	  
		  
		  a->add( "." , interval_t( start_tp , stop_tp ) , ch );
		}
	    }
	  
	}
  

      //
      // to cache (points --> for TLOCK)
      //

      if ( to_cache )
	{
	  
	  if ( na == valid ) 
	    cache->add( ckey_t( "points" , faclvl ) , mxpks );
	  else
	    {
	      std::vector<int> mxpks2(valid);
	      int idx = 0;
	      for (int i=0;i<na;i++) 
		if ( okay[i] ) mxpks2[ idx++ ] = mxpks[i]; 
	      cache->add( ckey_t( "points" , faclvl ) , mxpks2 );

	    }

	  // for (int i=0; i<na; i++)
	  //   cache->add( ckey_t( "points" , faclvl ) , mxpks[i] );	  
	}
      
    
      // next signal
    }

  writer.unlevel( globals::signal_strat);
}
