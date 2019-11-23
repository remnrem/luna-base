
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

void pdc_t::simple_sleep_scorer( edf_t & edf , param_t & param )
{
  
  //
  // have we already attached a PDLIB, i.e. 'reference'?
  //
  
  if ( ! param.has( "pd-lib" ) ) 
    Helper::halt( "required pdlib={library-file} option missing" );
  
  std::string pdlib = param.requires( "pd-lib" );
  

  //
  // which channels specified? EEG, EOG, EMG
  //
  
  bool has_eeg = param.has( "eeg" );
  bool has_eog = param.has( "eog" );
  bool has_emg = param.has( "emg" );

  std::string eeg_label = has_eeg ? Helper::unquote( param.value( "eeg" ) ) : "" ; 
  std::string eog_label = has_eog ? Helper::unquote( param.value( "eog" ) ) : "" ; 
  std::string emg_label = has_emg ? Helper::unquote( param.value( "emg" ) ) : "" ; 

  if ( eeg_label == "." ) has_eeg = false;
  if ( eog_label == "." ) has_eog = false;
  if ( emg_label == "." ) has_emg = false;

  std::string signal_label = has_eeg ? eeg_label : ""; 

  if ( has_eog ) 
    {
      if ( signal_label != "" ) signal_label += ",";
      signal_label += eog_label;
    }

  if ( has_emg ) 
    {
      if ( signal_label != "" ) signal_label += ",";
      signal_label += emg_label;
    }
  
  signal_list_t signals = edf.header.signal_list( signal_label );  

  const int ns = signals.size();

  if ( ns == 0 ) Helper::halt( "no signals specified" );

  logger << " using " << ns << " signals for SSS\n";
  
  
  std::map<std::string,std::string> chmap;

  if ( has_eeg ) { chmap[ eeg_label ] = "EEG"; } 
  if ( has_eog ) { chmap[ eog_label ] = "EOG"; }
  if ( has_emg ) { chmap[ emg_label ] = "EMG"; }
  
  //
  // if not already loaded, attach (i.e. only do once when iterating over multiple EDFs, i.e. 'targets')
  //
  
  if ( obs.size() == 0 ) 
    {      
      
      // add channel names
      //       for ( int s=0; s<ns; s++ )
      // 	add_channel( signals.label(s) );
      
      if ( has_eeg ) { add_channel( "EEG" ); }
      if ( has_eog ) { add_channel( "EOG" ); }
      if ( has_emg ) { add_channel( "EMG" ); }
      
      read_pdlib( pdlib );
      
      if ( obs.size() == 0 ) 
	Helper::halt( "no valid PDLIB specified" );
    }
  


  //
  // Parameters
  // 

  // epoch pre-grouping, e.g. based on ExE clustering
  // these pre-grouped epochs will be staged together
  // all epochs should be account for in the pre-grouping
  
  bool pre_grouping = param.has( "grouping" ) ;
  std::map<int,std::set<int> > group2epochs; // makes assumption that epoch grouping must be similar
  std::map<int,int> epoch2group; // reverse
  if ( pre_grouping ) 
    {
      std::string filename = param.requires( "grouping" ) + "-" + edf.id + ".txt";
      if ( ! Helper::fileExists( filename ) ) Helper::halt( "cannot find " + filename );
      std::ifstream IN1( filename.c_str() , std::ios::in );
      while ( ! IN1.eof() )
	{
	  std::string id;
	  int e;
	  int c;
	  IN1 >> id >> e >> c;
	  if ( IN1.eof() ) continue;
	  group2epochs[ c ].insert( e );
	  epoch2group[ e ] = c ;
	}
      IN1.close();
    }
  

  // 'm' and 't' values that will be used
  const int encoding_m = param.requires_int( "m" );
  const int encoding_t = param.requires_int( "t" );
  
  // number of top matches to take
  const int nmatch = param.has( "top" ) ? param.requires_int( "top" ) : 100 ; 

  // desired SRs, i.e. must match PDLIB for meaningful comparisons  
  int sr = param.requires_int( "sr" );
  
  // output: XML
  std::string xml_filename = param.has( "xml" ) ? param.value( "xml" ) : "";   
  if ( xml_filename == "." ) xml_filename = "";

  // verbose mode
  bool verbose = param.has( "verbose" );

  // output pdlib for the tested individuals? (as txt)
  bool out_pd = param.has( "pd" );
  
  std::ofstream OUT1;
  if ( out_pd ) 
    {
      std::string out_pdfile = param.value( "pd" );
      if ( out_pdfile != "." )
	{
	  out_pdfile += "-" + edf.id;
	  OUT1.open( out_pdfile.c_str() , std::ios::out );
	}
      else out_pd = false;
    }

  // actual SRs
  std::vector<double> Fs = edf.header.sampling_freq( signals );

  // resampling if neeeded
  for (int s=0;s<ns;s++)
    {
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
      
      if ( edf.header.sampling_freq( signals(s) ) != sr ) 
	dsptools::resample_channel( edf, signals(s) , sr );
    }  
  
  
  //
  // Consider each eppch, and split into 3 10 second components;  match each against the PDLIB
  //

  // epoch x 3 vector
  std::vector<std::vector<pdc_obs_t> > targets;
  
  
  //
  // Iterate over each epoch
  //
  
  int ne = edf.timeline.first_epoch();
  
  if ( pre_grouping ) 
    {
      if ( ne != epoch2group.size() ) 
	{
	  std::cout << " ne = " << ne << " " << epoch2group.size() << "\n";
	  Helper::halt( "not all epochs pre-grouped" );
	}
    }

  
  std::map<int,int> displayepoch2internal;
  
  int cnt = 0;
  
  while ( 1 ) 
    {
      
      int epoch = edf.timeline.next_epoch();      
      
      if ( epoch == -1 ) break;
      
      // track, as we may use this if doing grouping, see below
      displayepoch2internal[ edf.timeline.display_epoch( epoch ) ] = epoch; 

      interval_t interval = edf.timeline.epoch( epoch );
      
      // epoch counter
      ++cnt;
      
      // 30 ten-second intervals per epoch
      pdc_obs_t e1(ns), e2(ns), e3(ns);
      
      // select first, middle and last 10 seconds
      int start1 = 0;
      int end1   = start1 + 10 * sr - 1 ;
      
      int start2 = 10 * sr - 1 ;
      int end2   = start2 + 10 * sr - 1 ;
      
      int start3 = 20 * sr - 1 ;
      int end3   = start3 + 10 * sr - 1 ;
      
            
      for ( int s=0; s<ns; s++ )
	{
	  
	  // only consider data tracks
	  
	  if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
	  
	  int c = channel( chmap[ signals.label(s) ] );
	  
	  if ( c == -1 ) continue;

	  // extract signal from EDF

	  slice_t slice( edf , signals(s) , interval );
	  
	  std::vector<double> * d = slice.nonconst_pdata();
	  
	  const std::vector<uint64_t> * tp = slice.ptimepoints();
	  
	  // check epoch length is exactly 30s, otherwise skip
	  const int na = d->size();
	  if ( na != sr * 30 ) continue;
	  
	  
	  std::vector<double> ts1, ts2, ts3;
	  for (int j=start1; j<=end1; j++) ts1.push_back( (*d)[j] );
	  for (int j=start2; j<=end2; j++) ts2.push_back( (*d)[j] );
	  for (int j=start3; j<=end3; j++) ts3.push_back( (*d)[j] );
	  
	  // record
	  e1.ch[c] = e2.ch[c] = e3.ch[c] = true;

	  e1.ts[c] = ts1;
	  e2.ts[c] = ts2;
	  e3.ts[c] = ts3;
	  
	}
      
      // compile
      e1.encode( encoding_m, encoding_t );
      e2.encode( encoding_m, encoding_t );
      e3.encode( encoding_m, encoding_t );
      
      std::vector<pdc_obs_t> te;
      te.reserve(3);
      te.push_back( e1 );
      te.push_back( e2 );
      te.push_back( e3 );
      targets.push_back( te );
      
      
      //
      // Optional output of PD for the to-be-staged epochs?
      //
      
      if ( out_pd ) 
	{
	  
	  // e1	  

	  std::map<std::string,int>::const_iterator cc = channels.begin();
	  while ( cc != channels.end() )
	    {
	      OUT1 << edf.id << "\t" 
		   << cc->first << "\t"
		   << epoch+1 << "\t"
		   << 1 ;
	      for (int p=0;p<e1.pd[cc->second].size();p++) OUT1 << "\t" << e1.pd[ cc->second ][p];
	      OUT1 << "\n";
	      
	      OUT1 << edf.id << "\t" 
		   << cc->first << "\t"
		   << epoch+1 << "\t"
		   << 2 ;
	      for (int p=0;p<e2.pd[cc->second].size();p++) OUT1 << "\t" << e2.pd[ cc->second ][p];
	      OUT1 << "\n";

	      OUT1 << edf.id << "\t" 
		   << cc->first << "\t"
		   << epoch+1 << "\t"
		   << 3 ;
	      for (int p=0;p<e3.pd[cc->second].size();p++) OUT1 << "\t" << e3.pd[ cc->second ][p];
	      OUT1 << "\n";
	      
	      ++cc;
	    }
	}
    }

  
  if ( out_pd ) OUT1.close();
  
  logger << " scanned " << cnt << " epochs, extracted " << targets.size() << " time-series\n";

  
  //
  // Optional: if pre-groups specified, also add a series of targets that are the means of the corresponding
  //  epoch level (actually 10-sec level) intervals
  //
  

  std::map<int,std::string> pre_grouped_match;
  std::map<int,double> pre_grouped_conf;
  
  if ( pre_grouping )
    {
      
      std::vector<pdc_obs_t> grouped_targets;
      
      grouped_targets.resize( group2epochs.size() , pdc_obs_t(q) );
      
      std::map<int,std::set<int> >::const_iterator ii = group2epochs.begin();
      while ( ii != group2epochs.end() )
	{
	  
	  pdc_obs_t & t = grouped_targets[ ii->first ];
	  
	  std::set<int>::const_iterator ee = ii->second.begin();
	  while ( ee != ii->second.end() )
	    {
	      
	      if ( displayepoch2internal.find( *ee ) == displayepoch2internal.end() ) 
		{
		  ++ee; continue;
		}
	      
	      // all 3 10-second epochs
	      int e0 = displayepoch2internal[ *ee ];
	      t.add( targets[ e0 ][0] ) ;
	      t.add( targets[ e0 ][1] ) ;
	      t.add( targets[ e0 ][2] ) ;
	      
	      ++ee;
	    }
	  ++ii;
	}
     
      logger << " S1 \n";

      //
      // Normalize & get group-level matches
      //

      
      for (int c=0; c<grouped_targets.size(); c++)
	{
	  std::cout << "c = " << c << "\n";

	  // group labels should be 0, 1, 2 etc and so match 'c'  
	  
	  std::map<int,std::set<int> >::const_iterator ii = group2epochs.find( c );
	  if ( ii == group2epochs.end() ) Helper::halt( "internal error in pre-grouping" );
	  const std::set<int> & inds = ii->second;
	  
	  pdc_obs_t & t = grouped_targets[ c ];
	  
	  // scale so that PD some to 1.0 (i.e. we just get the average PD across all epochs belonging
	  // to that group (noting that each epoch contributes 3 times)
	  t.norm( group2epochs.size() * 3 );
	  
	  std::set<pd_dist_t> matches1 = match( t , nmatch );
	  
	  std::string match1;
	  double conf1;
	  
	  std::map<std::string,double> summ1 = summarize( matches1 , &match1 , &conf1 );
	  
	  // paste to individual epochs
	  std::cout << " ..\n";
	  std::set<int>::const_iterator ee = inds.begin();
	  while ( ee != inds.end() ) 
	    {
	      int e0 = displayepoch2internal[ *ee ];	      
	      pre_grouped_match[ e0 ] = match1;
	      pre_grouped_conf[ e0 ] = conf1;
	      ++ee;
	    }
	  
	}
      std::cout << "Done\n";
    }
  


  //
  // For each segment, find the best match; 
  //
  
  std::vector<std::string> stages;
  
  for (int e=0; e<targets.size(); e++)
    {

      writer.epoch( edf.timeline.display_epoch( e ) );

      // each epoch has three 10-second intervals
      std::set<pd_dist_t> matches1 = match( targets[e][0] , nmatch );
      std::set<pd_dist_t> matches2 = match( targets[e][1] , nmatch );
      std::set<pd_dist_t> matches3 = match( targets[e][2] , nmatch );
      
      std::string match1, match2, match3;
      double conf1, conf2, conf3;

      std::map<std::string,double> summ1 = summarize( matches1 , &match1 , &conf1 );
      std::map<std::string,double> summ2 = summarize( matches2 , &match2 , &conf2 );
      std::map<std::string,double> summ3 = summarize( matches3 , &match3 , &conf3 );
      
      //
      // Get a 2-of-3 consensus match?
      //
      std::string final_match = ".";

      if ( match1 == match2 || match1 == match3 ) final_match = match1;
      else if ( match2 == match3 ) final_match = match2;

      if ( 1 || verbose ) 
	{
	  
//  	  std::cout << "epoch " << e+1 << "\t"	    
//  		    << final_match << "\t";
	  
//  	  // individual matches, and confidence scores:
//  	  std::cout << match1 << "\t" << conf1 << "\t"
//  		    << match2 << "\t" << conf2 << "\t"
//  		    << match3 << "\t" << conf3 << "\n";

	   std::map<std::string,double>::const_iterator ss1 = summ1.begin();
	   while ( ss1 != summ1.end() )
	     {
	       writer.value ( "SS1_" + ss1->first , ss1->second );
	       //	       std::cout << " " << ss1->first << "=" << ss1->second ; 
	       ++ss1;
	     }
	   //std::cout << "\n";

	   //std::cout << "match2 " << match2 << "\t" << conf2 << "\n";
	   std::map<std::string,double>::const_iterator ss2 = summ2.begin();
	   while ( ss2 != summ2.end() )
	     {
	       writer.value ( "SS2_" + ss2->first , ss2->second );
	       //std::cout << " " << ss2->first << "=" << ss2->second ; 
	       ++ss2;
	     }
// 	   std::cout << "\n";
	   
// 	   std::cout << "match3 " << match3 << "\t" << conf3 << "\n";
	   std::map<std::string,double>::const_iterator ss3 = summ3.begin();
	   while ( ss3 != summ3.end() )
	     {
	       writer.value ( "SS3_" + ss3->first , ss3->second );
	       //std::cout << " " << ss3->first << "=" << ss3->second ; 
	       ++ss3;
	     }
// 	   std::cout << "\n";
	  
// 	   std::cout << "\n";
	  
	}

      
      if      ( final_match == "N1" ) stages.push_back( "Stage 1 sleep|1" );
      else if ( final_match == "N2" ) stages.push_back( "Stage 2 sleep|2" );
      else if ( final_match == "N3" ) stages.push_back( "Stage 3 sleep|3" );
      else if ( final_match == "R" ) stages.push_back( "REM sleep|5" );
      else if ( final_match == "W" ) stages.push_back( "Wake|0" );
      else stages.push_back( "Unsure|Unsure" );

      //
      // Output
      //
      
      
      writer.value( "SS" , final_match );
      
      writer.value( "SS1" , match1 );
      writer.value( "SS2" , match2 );
      writer.value( "SS3" , match3 );

      writer.value( "CONF1" , conf1 );
      writer.value( "CONF2" , conf2 );
      writer.value( "CONF3" , conf3 );

      if ( pre_grouping ) 
	{
	  writer.value( "SS_G" , pre_grouped_match[ e ] );
	  writer.value( "CONF_G" , pre_grouped_conf[ e ] );	  
	}


    }
  
  writer.unepoch();
  
  
  //
  // Write to file
  //

  if ( xml_filename != "" ) 
    write_xml( xml_filename + "-" + edf.id + ".xml" , stages );
  
  
}

void pdc_t::write_xml( const std::string & filename , const std::vector<std::string> & stages )
{

  std::ofstream XML( filename.c_str() , std::ios::out );

  XML << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n";
  XML << "<PSGAnnotation>"
      << "<ScoredEvents>\n";

  int timer = 0;

  for (int e = 0 ; e < stages.size() ; e++ ) 
    {

      XML << "<ScoredEvent>"
	  << "<EventType>Stages|Stages</EventType>"
	  << "<EventConcept>" << stages[e] << "</EventConcept>"
	  << "<Start>" << timer << "</Start>"
	  << "<Duration>30.0</Duration>"
	  << "</ScoredEvent>\n";

      // advance to next epoch
      timer += 30;

    }

  
  // all done
  
  XML << "</ScoredEvents>"
      << "</PSGAnnotation>";
  
  XML.close();
}


std::set<pd_dist_t> pdc_t::match( const pdc_obs_t & target , const int nbest )
{

  const int N = obs.size();
  
  std::set<pd_dist_t> dist, final;
  
  for (int i = 0 ; i < N ; i++ ) 
    {
      //      std::cout << "MM " << distance( target , obs[i] ) << " for " << i << " " << obs[i].label << "\n";
      dist.insert( pd_dist_t( distance( target , obs[i] ) , i ) );
    }

  
  int cnt = 0 ; 
  std::set<pd_dist_t>::const_iterator ii = dist.begin();
  while ( ii != dist.end() )
    {
      //std::cout << "FF " << cnt+1 << " " << ii->ix << " " << obs[ ii->ix ].label << " " << ii->d << "\n";
      final.insert( *ii );
      ++cnt;
      if ( cnt == nbest ) break;
      ++ii;
    }
  
  return final;
}

std::map<std::string,double> pdc_t::summarize( const std::set<pd_dist_t> & matches , std::string * cat , double * conf )
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
	  s[ *ll ] = scnt[ *ll ] / (double)label_count[ *ll ];
	  //s[ *ll ] = s[ *ll ] / total;
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



