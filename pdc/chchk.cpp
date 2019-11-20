
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


#include "pdc.h"

#include <string>
#include <iostream>
#include <cmath>
#include <set>

#include "eval.h"

#include "edf/edf.h"
#include "edf/slice.h"

#include "miscmath/miscmath.h"
#include "dsp/resample.h"

#include "db/db.h"
#include "helper/logger.h"

extern writer_t writer;

extern logger_t logger;

void pdc_t::channel_checker( edf_t & edf , param_t & param )
{
  
  //
  // have we already attached a PDLIB, i.e. 'reference'?
  //
  
  if ( ! param.has( "pd-lib" ) ) 
    Helper::halt( "required pdlib={library-file} option missing" );
  
  std::string pdlib = param.requires( "pd-lib" );
  

  //
  // Signals to check
  //

  std::string signal_label = param.requires( "sig" );

  signal_list_t signals = edf.header.signal_list( signal_label );  

  const int ns = signals.size();

  if ( ns == 0 ) Helper::halt( "no signals specified" );

  
  //
  // if not already loaded, attach library 
  // i.e. only do once when iterating over multiple EDFs
  //
  
  if ( obs.size() == 0 ) 
    {      
      
      std::vector<std::string> references;

      references.push_back( "eog" );
      references.push_back( "eeg" );
      references.push_back( "emg" );
      references.push_back( "ecg" );
      references.push_back( "air" );
      references.push_back( "leg" );
//       references.push_back( "SaO2" );
//       references.push_back( "SpO2" );
      references.push_back( "snore" );
      references.push_back( "abdo" );
      references.push_back( "nasal" );

      // add all target channel names
      //       for ( int s=0; s<ns; s++ )
      // 	add_channel( signals.label(s) );
      
      
      read_pdlib( pdlib );
      
      if ( obs.size() == 0 ) 
	Helper::halt( "no valid PDLIB specified" );
    }
  

  
  // 'm' and 't' values that will be used
  const int encoding_m = param.requires_int( "m" );
  const int encoding_t = param.requires_int( "t" );
  
  // number of top matches to take
  const int nmatch = param.has( "top" ) ? param.requires_int( "top" ) : 100 ; 

  // desired SRs, i.e. must match PDLIB for meaningful comparisons  
  int sr = param.requires_int( "sr" );
  
  // verbose mode
  bool verbose = param.has( "verbose" );

  // actual SRs
  std::vector<double> Fs = edf.header.sampling_freq( signals );

  // resampling if neeeded
  for (int s=0;s<ns;s++)
    {
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
      
      if ( edf.header.sampling_freq( signals(s) ) != sr ) 
	dsptools::resample_channel( edf, signals(s) , sr );
    }  
  
  

//   //
//   // Consider each eppch, and split into 3 10 second components;  match each against the PDLIB
//   //

//   // epoch x 3 vector
//   std::vector<std::vector<pdc_obs_t> > targets;
  
  
//   //
//   // Iterate over each epoch
//   //
  
//   int ne = edf.timeline.first_epoch();
  
//   if ( pre_grouping ) 
//     {
//       if ( ne != epoch2group.size() ) 
// 	{
// 	  std::cout << " ne = " << ne << " " << epoch2group.size() << "\n";
// 	  Helper::halt( "not all epochs pre-grouped" );
// 	}
//     }

  
//   std::map<int,int> displayepoch2internal;
  
//   int cnt = 0;
  
//   while ( 1 ) 
//     {
      
//       int epoch = edf.timeline.next_epoch();      
      
//       if ( epoch == -1 ) break;
      
//       // track, as we may use this if doing grouping, see below
//       displayepoch2internal[ edf.timeline.display_epoch( epoch ) ] = epoch; 

//       interval_t interval = edf.timeline.epoch( epoch );
      
//       // epoch counter
//       ++cnt;
      
//       // 30 ten-second intervals per epoch
//       pdc_obs_t e1(ns), e2(ns), e3(ns);
      
//       // select first, middle and last 10 seconds
//       int start1 = 0;
//       int end1   = start1 + 10 * sr - 1 ;
      
//       int start2 = 10 * sr - 1 ;
//       int end2   = start2 + 10 * sr - 1 ;
      
//       int start3 = 20 * sr - 1 ;
//       int end3   = start3 + 10 * sr - 1 ;
      
            
//       for ( int s=0; s<ns; s++ )
// 	{
	  
// 	  // only consider data tracks
	  
// 	  if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
	  
// 	  int c = channel( chmap[ signals.label(s) ] );
	  
// 	  if ( c == -1 ) continue;

// 	  // extract signal from EDF

// 	  slice_t slice( edf , signals(s) , interval );
	  
// 	  std::vector<double> * d = slice.nonconst_pdata();
	  
// 	  const std::vector<uint64_t> * tp = slice.ptimepoints();
	  
// 	  // check epoch length is exactly 30s, otherwise skip
// 	  const int na = d->size();
// 	  if ( na != sr * 30 ) continue;
	  
	  
// 	  std::vector<double> ts1, ts2, ts3;
// 	  for (int j=start1; j<=end1; j++) ts1.push_back( (*d)[j] );
// 	  for (int j=start2; j<=end2; j++) ts2.push_back( (*d)[j] );
// 	  for (int j=start3; j<=end3; j++) ts3.push_back( (*d)[j] );
	  
// 	  // record
// 	  e1.ch[c] = e2.ch[c] = e3.ch[c] = true;

// 	  e1.ts[c] = ts1;
// 	  e2.ts[c] = ts2;
// 	  e3.ts[c] = ts3;
	  
// 	}
      
//       // compile
//       e1.encode( encoding_m, encoding_t );
//       e2.encode( encoding_m, encoding_t );
//       e3.encode( encoding_m, encoding_t );
      
//       std::vector<pdc_obs_t> te;
//       te.reserve(3);
//       te.push_back( e1 );
//       te.push_back( e2 );
//       te.push_back( e3 );
//       targets.push_back( te );
      
      
//       //
//       // Optional output of PD for the to-be-staged epochs?
//       //
      
//       if ( out_pd ) 
// 	{
	  
// 	  // e1	  

// 	  std::map<std::string,int>::const_iterator cc = channels.begin();
// 	  while ( cc != channels.end() )
// 	    {
// 	      OUT1 << edf.id << "\t" 
// 		   << cc->first << "\t"
// 		   << epoch+1 << "\t"
// 		   << 1 ;
// 	      for (int p=0;p<e1.pd[cc->second].size();p++) OUT1 << "\t" << e1.pd[ cc->second ][p];
// 	      OUT1 << "\n";
	      
// 	      OUT1 << edf.id << "\t" 
// 		   << cc->first << "\t"
// 		   << epoch+1 << "\t"
// 		   << 2 ;
// 	      for (int p=0;p<e2.pd[cc->second].size();p++) OUT1 << "\t" << e2.pd[ cc->second ][p];
// 	      OUT1 << "\n";

// 	      OUT1 << edf.id << "\t" 
// 		   << cc->first << "\t"
// 		   << epoch+1 << "\t"
// 		   << 3 ;
// 	      for (int p=0;p<e3.pd[cc->second].size();p++) OUT1 << "\t" << e3.pd[ cc->second ][p];
// 	      OUT1 << "\n";
	      
// 	      ++cc;
// 	    }
// 	}
//     }

  
//   if ( out_pd ) OUT1.close();
  
//   logger << " scanned " << cnt << " epochs, extracted " << targets.size() << " time-series\n";

    
}



