
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

      // expecting only the single 'test' channel _T
      add_channel( "_T" ); 
      
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
      
      read_pdlib( pdlib );
      
      if ( obs.size() == 0 ) 
	Helper::halt( "no valid PDLIB specified" );
      else
	logger << " read " << obs.size() << " obs\n";
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



  ///
  // Iterate over signals
  //

  
  for (int s=0;s<ns;s++)
    {      
      

      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
      
      writer.level( signals.label(s) , globals::signal_strat );
      
      //
      // resample as neeeded
      //

      if ( edf.header.sampling_freq( signals(s) ) != sr ) 
	dsptools::resample_channel( edf, signals(s) , sr );

      //   // epoch x 3 vector
      //   std::vector<std::vector<pdc_obs_t> > targets;
  
      
      //
      // Iterate over each epoch
      //
      
      int ne = edf.timeline.first_epoch();
  
  
//   std::map<int,int> displayepoch2internal;
  
//   int cnt = 0;

       //
       // Always map against _T test-channel label (i.e. all pdlib channels)
       //
       
       int c = channel( "_T" );
  




        while ( 1 ) 
        	 {
	   
        	   int epoch = edf.timeline.next_epoch();      
	   
        	   if ( epoch == -1 ) break;
	   
		   interval_t interval = edf.timeline.epoch( epoch );
		   //interval_t interval = edf.timeline.wholetrace();
	   
	   // extract signal from EDF
	   
	   slice_t slice( edf , signals(s) , interval );
	   
	   std::vector<double> * d = slice.nonconst_pdata();
	   
	   // check epoch length is exactly 30s, otherwise skip
	   const int na = d->size();
	   //if ( na != sr * 30 ) continue;
	   
	   // 30-second intervals per epoch, for 1 channel
	   pdc_obs_t e1(1);
	   pdc_obs_t e1r(1);
	   
	   // record
	   e1.ch[c] = true;
	   e1.ts[c] = *d;
	   
	   // reverse
	   e1r.ch[c] = true;
	   e1r.ts[c] = *d;
	   for (int i=0;i<na;i++) 
	     e1r.ts[c][i] *= -1;

      	   // encode
	   
	   e1.encode( encoding_m, encoding_t );
	   e1r.encode( encoding_m, encoding_t );

	   
	   std::vector<pdc_obs_t> te;
	   te.reserve(3);
	   // te.push_back( e1 );
	   // targets.push_back( te );

	   std::set<pd_dist_t> matches1 = match( e1 , nmatch );
	   std::set<pd_dist_t> matches1r = match( e1r , nmatch );
	   

	   std::string match1, match1r;
	   double conf1, conf1r;
	   
	   std::map<std::string,double> summ1 = summarize2( matches1 , &match1 , &conf1 );
	   std::map<std::string,double> summ1r = summarize2( matches1r , &match1r , &conf1r );
	   
   	   // individual matches, and confidence scores:
	   std::cout << match1 << "\t" << conf1 << "\n";
	   
	   writer.epoch( edf.timeline.display_epoch( epoch ) );
	   
	   writer.value( "MATCH" , match1 );
	   writer.value( "RMATCH" , match1r );

	   writer.value( "CONF" , conf1 );
	   writer.value( "RCONF" , conf1r );

	   //   logger << " scanned " << cnt << " epochs, extracted " << targets.size() << " time-series\n";
	   
        	 }
       
        writer.unepoch();
       
    }
  
  writer.unlevel ( globals::signal_strat );
  
}



std::map<std::string,double> pdc_t::summarize2( const std::set<pd_dist_t> & matches , std::string * cat , double * conf )
{
  
  const int nmatches = matches.size();

  std::map<std::string,double> s;
  std::map<std::string,int> scnt;

  if ( nmatches == 0 ) return s;

  double pdmin = matches.begin()->d;
  double pdmax = matches.begin()->d;

  std::vector<pd_dist_t> scaled;
  
  std::set<pd_dist_t>::const_iterator ii = matches.begin();
  while ( ii != matches.end() ) 
    {
      if ( ii->d < pdmin ) pdmin = ii->d;
      if ( ii->d > pdmax ) pdmax = ii->d;      
      scaled.push_back( *ii );
      ++ii;
    }

  // scale so that within this 'best-match' set, best is 1.0 and worst is 0.
  // then sum for each label

  // and then get sum per label in top N, as we enfore the same
  // number of templates, this is equivalent to having a uniform prior
  // for each class

  double total = 0;

  std::vector<pd_dist_t>::iterator jj = scaled.begin();
  while ( jj != scaled.end() ) 
    {

      jj->d = 1 - ( ( jj->d - pdmin ) / ( pdmax - pdmin ) ) ;

      //      std::cout << "jj->d = " << jj->d << "\n";

      // sum per stage
      s[ obs[ jj->ix ].label ] += jj->d;

      // sum overall
      total += jj->d;

      // track how many obs per 
      scnt[ obs[ jj->ix ].label ]++;
      
      ++jj;
    }
  
  
  // normalize to 1.0 over all labels and take complement
  
  //  double total2 = 0 ; 
  
  std::set<std::string>::const_iterator ll = labels.begin();
  while ( ll != labels.end() )
    {
      
      if ( s.find( *ll ) == s.end() ) 
	{
	  s[ *ll ] = 0 ; 
	}
      else
	{
	  s[ *ll ] /= scnt[ *ll ] ;
	  //	  s[ *ll ] = scnt[ *ll ] / (double)label_count[ *ll ];
	  //s[ *ll ] = s[ *ll ] / total;
	  //std::cout <<  "  " << s[ *ll ]  << " is " << *ll << "\n";
	}
      
      ++ll;
    }
  
  // select the best label
  
  *conf = 0;
  *cat = ".";
  
  ll = labels.begin();
  while ( ll != labels.end() )
    {      
      if ( s[ *ll ] > *conf ) 
	{
	  *cat = *ll;
	  *conf = s[ *ll ];
	}
      
      ++ll;
    }
  

  return s;
}



