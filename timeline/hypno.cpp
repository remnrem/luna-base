
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

#include "hypno.h"

#include "edf/edf.h"
#include "defs/defs.h"
#include "db/db.h"
#include "miscmath/crandom.h"
#include "eval.h"
#include "intervals/intervals.h"
#include "annot/annot.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "helper/token-eval.h"
#include "dsp/lzw.h"
#include "dsp/mse.h"

extern writer_t writer;

extern logger_t logger;

sleep_stage_t hypnogram_t::gap_treatment = UNKNOWN;


bool is_rem( sleep_stage_t s ) { return s == REM; } 
bool is_nrem( sleep_stage_t s ) { return s == NREM1 || s == NREM2 || s == NREM3 || s == NREM4; } 
bool is_nrem1( sleep_stage_t s ) { return s == NREM1; } 
bool is_nrem2( sleep_stage_t s ) { return s == NREM2; } 
bool is_nrem23( sleep_stage_t s ) { return s == NREM2 || s == NREM3; } 
bool is_nrem34( sleep_stage_t s ) { return s == NREM3 || s == NREM4; } 
bool is_nrem234( sleep_stage_t s ) { return s == NREM2 || s == NREM3 || s == NREM4; } 
bool is_wake( sleep_stage_t s ) { return s == WAKE || ( s == GAP && hypnogram_t::gap_treatment == WAKE ) ; }
bool is_wake_or_lights( sleep_stage_t s ) { return is_wake(s) || s == LIGHTS_ON; } 
bool is_sleep( sleep_stage_t s ) { return s == NREM1 || s == NREM2 || s == NREM3 || s == NREM4 || s == REM ; } 
bool is_absent( sleep_stage_t s ) { return s == UNSCORED || s == UNKNOWN || s == MOVEMENT || s == LIGHTS_ON || s == ARTIFACT || ( s == GAP && hypnogram_t::gap_treatment == UNKNOWN ); } 
bool is_gap( sleep_stage_t s ) { return s == GAP; }
bool is_obersved( sleep_stage_t s ) { return s != GAP; }

bool is_same_3class( sleep_stage_t s1 , sleep_stage_t s2 ) 
{ 
  if ( s1 == s2 ) return true;
  if ( ( s1 == NREM1 || s1 == NREM2 || s1 == NREM3 || s1 == NREM4 ) && 
       ( s2 == NREM1 || s2 == NREM2 || s2 == NREM3 || s2 == NREM4 ) ) return true;
  return false;
}

bool hypnogram_t::construct( timeline_t * t , param_t & param , const bool verbose , const std::vector<std::string> & s )
{ 
  timeline = t;
  req_pre_post_epochs = param.has( "req-pre-post" ) ? param.requires_int( "req-pre-post" ) : 4;
  flanking_3class = param.has( "flanking-collapse-nrem" ) ? Helper::yesno( param.value( "flanking-collapse-nrem") ) : true;

  if ( s.size() != timeline->num_total_epochs() ) 
    Helper::halt( "bad number of stages, " 
		  + Helper::int2str( (int)s.size() ) 
		  + " but expecting " 
		  + Helper::int2str( timeline->num_total_epochs() ) );    
  stages.resize( s.size() );
  for (int e=0;e<s.size();e++) stages[e] = globals::stage( s[e] );
  original_stages = stages;
  edit( t , param );
  if ( empty() ) 
    {
      logger << " ** warning, no valid stage-annotated epochs, bailing...\n"; 
      return false;
    }
  calc_stats( verbose );
  return true;
}

bool hypnogram_t::construct( timeline_t * t , param_t & param , const bool verbose , const std::string sslabel ) 
{

  // point to 'parent' timeline
  timeline = t ;

  // set any params
  req_pre_post_epochs = param.has( "req-pre-post" ) ? param.requires_int( "req-pre-post" ) : 4;
  flanking_3class = param.has( "flanking-collapse-nrem" ) ? Helper::yesno( param.value( "flanking-collapse-nrem") ) : true;

  // get handle
  annot_t * annot = timeline->annotations->find( sslabel );
  if ( annot == NULL ) 
    {
      logger << "  did not find any existing, valid sleep stage annotations...\n";
      return false;
    }

  //
  // set internal, epoch-level annotations used by timeline
  //
  
  std::set<std::string> values;
  values.clear(); values.insert( "W" );
  timeline->annotate_epochs(  globals::stage( WAKE ) , "SleepStage" , values );
  
  values.clear(); values.insert( "N1" );
  timeline->annotate_epochs(  globals::stage( NREM1  )  , "SleepStage" , values );

  values.clear(); values.insert( "N2" );
  timeline->annotate_epochs(  globals::stage( NREM2  )  , "SleepStage" , values );

  values.clear(); values.insert( "N3" );
  if ( collapse_nrem34 ) { values.insert( "NREM4" ); values.insert( "N4" ); }
  timeline->annotate_epochs(  globals::stage( NREM3  )  , "SleepStage" , values );
  
  if ( ! collapse_nrem34 ) 
    {
      values.clear(); values.insert( "NREM4" ); values.insert( "N4" );
      timeline->annotate_epochs(  globals::stage( NREM4 )  , "SleepStage" , values );
    }

  values.clear(); values.insert( "R" );
  timeline->annotate_epochs(  globals::stage( REM ) , "SleepStage" , values );
  
  values.clear(); values.insert( "L" );
  timeline->annotate_epochs(  globals::stage( LIGHTS_ON ) , "SleepStage" , values );

  //
  // Preliminary step to edit out any 'trailing' weird sleep epochs...
  //  i.e. if two hours of WAKE, then a single N1 epoch, then all wake...
  //  this appears to happen in, e.g. some MrOS studies
  //
  
  

  //
  // If we've masked the data, epoch count may not start at 0...
  // Although it probably doesn't make sense to use HYPNO then, 
  // we still want to be able to dump STAGE information easily w/ the
  // STAGE command...
  //
  

  // in VERBOSE (HYPNO) mode, we require the FULL epoch set
  //  -- although, note that we are adding support for HYPNO in EDF+D conttexts
  //     and so this is no different... 
  
  if ( verbose ) 
    {
      if ( timeline->num_total_epochs() !=  timeline->num_epochs() ) 
	Helper::halt( "cannot run HYPNO on masked data" );
      
      int eprev = -1;
      timeline->first_epoch();
      while ( 1 ) 
	{	  
	  int e = timeline->next_epoch();      
	  if ( e == -1 ) break;
	  if ( eprev >= 0 && timeline->display_epoch( e ) - eprev != 1 ) 
	    Helper::halt( "cannot run HYPNO on masked data" );
	  eprev =  timeline->display_epoch( e ) ;
	}
    }

  // this is number of observed epochs
  ne = timeline->num_total_epochs();
  
  timeline->first_epoch();

  // key measures populated here
  stages.clear();
  epoch_n.clear();
  // explicitly track length and start of each epoch
  //  (to allow for gaps)
  epoch_dur.clear();
  epoch_start.clear();
  epoch_gap.clear();

  // how ot handle gaps -- treat as "WAKE" or just as unknown? - default== MISSING
  gap_treatment = param.has( "gaps" ) && param.value( "gaps" ) == "W" ? WAKE : UNKNOWN; 
  
  //
  // canonical epoch sizes (for observed epochs, not gaps)
  //
  
  epoch_mins = timeline->epoch_length() / 60.0 ; 
  epoch_hrs = epoch_mins / 60.0;
  epoch_sec = timeline->epoch_length();

  //
  // We should be able to use current 0..ne epoch naming as epoch-annotations
  // still work after a restructure
  //
  
  n_conflicts = 0;

  uint64_t end_prior = 0;
  
  while ( 1 ) 
    {
      
      int e = timeline->next_epoch();

      if ( e == -1 ) break;
      
      writer.epoch( timeline->display_epoch( e ) );
      
      //
      // was there a gap prior to this epoch?
      //  - if contiguous, end of last == start of current
      //
      
      interval_t interval = timeline->epoch( e );
      
      if ( end_prior != 0 && end_prior != interval.start )
	{
	  
	  //std::cout << " found gap before epoch " << e << "( gap ends " << interval.start_sec() << "\n";
	  uint64_t gap_dur = interval.start - end_prior;
	  //std::cout << "   gap dur was " << gap_dur << " tp\n";
	  
	  // add a fake 'gap' epoch
	  // before this real one
	  stages.push_back( GAP );
	  epoch_gap.push_back( true );
	  epoch_start.push_back( interval.start_sec() );
	  epoch_dur.push_back( gap_dur * globals::tp_duration );
	  epoch_n.push_back( -1 );     // not used i.e. lookup to display epoch code 
	}

      //
      // update last prior stop point for the next epoch
      //

      end_prior = interval.stop; 
      

      //
      // for output of STAGES or HYPNO, use original EDF annotations
      //

      int e2 = timeline->original_epoch(e) ;
      
      bool wake = timeline->epoch_annotation( "W"  , e );
      bool n1   = timeline->epoch_annotation( "N1" , e );
      bool n2   = timeline->epoch_annotation( "N2" , e );
      bool n3   = timeline->epoch_annotation( "N3" , e );
      bool n4   = timeline->epoch_annotation( "NREM4" , e );
      bool rem  = timeline->epoch_annotation( "R"   , e );
      bool lights  = timeline->epoch_annotation( "L"   , e );      
      
      bool other = ! ( wake || n1 || n2 || n3 || n4 || rem || lights );
      bool conflict = ( (int)wake + (int)n1 + (int)n2 + (int)n3 + (int)n4 + (int)rem + (int)lights ) > 1;

      if ( conflict ) 
	{
	  std::cout << " conflict (>1 stage overlaps this epoch): " << interval.start << " " << interval.stop << "\n";
	}
      
      //
      // track any conflicts (i.e. if epochs not aligned to staging annotations)
      //
      
      if ( conflict )
	{
	  other = true;
	  ++n_conflicts;

	  std::stringstream ss;
	  bool delim = false;
	  if ( n1 ) { ss << "N1"; delim = true; }
	  if ( n2 ) { ss << ( delim ? "," : "" ) << "N2"; delim = true; }
	  if ( n3 ) { ss << ( delim ? "," : "" ) << "N3"; delim = true; }
	  if ( n4 ) { ss << ( delim ? "," : "" ) << ( collapse_nrem34 ? "N3" : "N4" ); delim = true; }
	  if ( rem ) { ss << ( delim ? "," : "" ) << "R"; delim = true; }
	  if ( wake ) { ss << ( delim ? "," : "" ) << "W"; delim = true; }	  
	  if ( lights ) { ss << ( delim ? "," : "" ) << "L"; delim = true; }	  
	  writer.value( "CONFLICT" , ss.str() );
	}
      
      // here (internally in hypno_t) we use UNKNOWN for all cases 
      if      ( conflict ) stages.push_back( UNKNOWN );
      else if ( other ) stages.push_back( UNKNOWN );
      else if ( wake ) stages.push_back( WAKE );
      else if ( n1 ) stages.push_back( NREM1 );
      else if ( n2 ) stages.push_back( NREM2 );
      else if ( n3 ) stages.push_back( NREM3 );
      else if ( n4 ) stages.push_back( collapse_nrem34 ? NREM3 : NREM4 );
      else if ( rem ) stages.push_back( REM );
      else if ( lights ) stages.push_back( LIGHTS_ON );
      else stages.push_back( UNKNOWN );
      
      // store original EDF 0-based encoding, to be passed to calc_stats()
      epoch_n.push_back( e2 );
      
      epoch_gap.push_back( false ); // is not a gap      
      epoch_start.push_back( interval.start_sec() ); // start of epoch in seconds (elapsed from EDF start) 
      epoch_dur.push_back( epoch_sec );              // epoch duration will be fixed for non-gaps
      
      // track times for basic epochs
    }

  writer.unepoch();

  // std::cout << " lengths\n"
  // 	    << stages.size() << " " << epoch_n.size() << " " << epoch_gap.size() << " " << epoch_start.size() << " " << epoch_dur.size() << "\n";

  //
  // make a copy of the stages
  //
  
  original_stages = stages;

  
  //
  // track total number of epochs + gaps
  //

  ne_gaps = stages.size();
  

  //
  // edit hypnogram as needed (e.g. for lights-off, excessive WASO, etc)
  //

  edit( timeline , param ); 
  
  //
  // If empty, just bail
  //
  
  if ( empty() ) 
    {
      logger << " ** warning, no valid stage-annotated epochs, bailing...\n"; 
      return false;
    }

  //
  // Report any conflicts
  //
  
  if ( n_conflicts )
    logger << "  *** found " << n_conflicts << " epoch(s) of " << ne << " with conflicting spanning annotations\n"
	   << "  *** check that epochs and annotations align as intended\n"
	   << "  *** see EPOCH 'align' or 'offset' options\n"; 
  
  
  
  //
  // finally, calculate hypno stats
  //
  
  calc_stats( verbose );

  //
  // If all missing, return false
  //
  
  return true;
}   

bool hypnogram_t::empty() const
{

  // do we have any of the five main stages (include NREM4 here too, so six)
  // i.e. if a hypnogram is all UNKNOWN, or LIGHTS, or GAP, etc, then doesn't count
  // as a valid hypnogram
  
  const int n = stages.size();

  for (int i=0; i<n; i++)
    {
      const sleep_stage_t & stg = stages[i];
      if ( stg == WAKE || stg == NREM1 || stg == NREM2
	   || stg == NREM3 || stg == NREM4 || stg == REM )
	return false;
    }
  return true;
  
}

void hypnogram_t::edit( timeline_t * timeline , param_t & param )
{

  // 
  // special case; nothing to do if empty
  //

  if ( empty() ) return;
   
  //
  // 1) Do we have any lights_on or lights_off annotations?

  // 2) Or have these been passed 
  // via the command line (e.g. as a variable, which may be individual-specific, lights-off=${LOFF} 
  // where LOFF is defined in a vars=data.txt file )
  //
  // 3) Or, have we been given a cache, which implies LON and LOFF (from trim)
  
  // lights-off=hh::mm::ss
  //  OR 
  // if numeric, assume this is seconds past EDF start
  //  explicitly check for ":" to denote hh:mm:ss

  double lights_off = -1;
  double lights_on = -1;

  clocktime_t st( timeline->edf->header.starttime );

  if ( param.has( "cache" ) )
    {
      std::string cname = param.value( "cache" );
      
      if ( ! timeline->cache.has_num( cname ) )
	logger << "  skipping... cache " << cname << " not found for this individual...\n";
      else
	{
	  
	  cache_t<double> * cache = timeline->cache.find_num( cname );

	  // expecting a single stratum... if not, complain (e.g. could reflect TAGs)
	  
	  std::set<ckey_t> ckeys_off = cache->keys( "LOFF" );
	  std::set<ckey_t> ckeys_on = cache->keys( "LON" );

	  if ( ckeys_off.size() > 1 ) Helper::halt( "expecting single stratum for cache" );
	  if ( ckeys_on.size() > 1 ) Helper::halt( "expecting single stratum for cache" );

	  // set lights off
	  if ( ckeys_off.size() )
	    {
	      std::vector<double> cx = cache->fetch( *ckeys_off.begin() );
	      if ( cx.size() != 1 ) Helper::halt( "internal error - expecting single cache entry for LOFF" );
 	      
	      lights_off = cx[0];

	      logger << "  from TRIM cache, setting lights_off = "
		     << lights_off << " secs, "
		     << lights_off/60.0 << " mins from start\n";
	      
	    }

	  // set lights on
	  if ( ckeys_on.size() )
	    {
	      std::vector<double> cx = cache->fetch( *ckeys_on.begin() );
	      if ( cx.size() != 1 ) Helper::halt( "internal error - expecting single cache entry for LON" );
	      
	      lights_on = cx[0];
	      
	      logger << "  from TRIM cache, setting lights_on = "
		     << lights_on << " secs, "
		     << lights_on/60.0 << " mins from start\n";
	      	      
	    }	  
	}
    }
  

  // else, if not already set, see about command line args
  if ( lights_off < 0 && param.has( "lights-off" ) && param.value( "lights-off" ) != "." && param.value( "lights-off" ) != "" )
    {
      const std::string loffstr = param.value( "lights-off" );
      const bool hms_mode = loffstr.find( ":" ) != std::string::npos;
      
      // argument in seconds:
      double x;
      if ( (!hms_mode) && Helper::str2dbl( loffstr , &x ) ) 
	{
	  if ( x < 0 ) 
	    {
	      logger << "  lights-off time less than 0 -- setting to 0 (EDF start)\n";
	      x = 0;
	    }

	  lights_off = x;
	  
	  logger << "  setting lights_off = " 
		 << lights_off << " secs, "
		 << lights_off/60.0 << " mins from start\n";
	}
      else if ( hms_mode )
	{
	  
	  // else assume is hh:mm:ss

	  if ( ! st.valid ) Helper::halt( "EDF does not have a valid start time - cannot use lights-off=hh:mm:ss" );
	  clocktime_t et( loffstr );
	  if ( et.valid ) 
	    {	  
	      int earlier = clocktime_t::earlier( st , et );
	      
	      if ( earlier == 2 )
		lights_off = 0;  // set to start of EDF
	      else
		lights_off = clocktime_t::ordered_difference_seconds( st , et );
	      
	      logger << "  setting lights_off = " << et.as_string() 
		     << " (" << lights_off << " secs, " 
		     << lights_off/60.0 << " mins from start)\n";
	    }
	  else
	    logger << "  invalid time for lights-off=" 
		   << loffstr << "  -- will ignore this\n";	  
	}
      else
	logger << "  invalid time for lights-off="
	       << loffstr << "  -- will ignore this\n";

    }
  
  //
  // Lights-On time 
  //

  if ( lights_on < 0 && param.has( "lights-on" ) && param.value( "lights-on" ) != "." && param.value( "lights-on" ) != "" )
    {
      
      const std::string lonstr = param.value( "lights-on" );
      const bool hms_mode = lonstr.find( ":" ) != std::string::npos;
      
      // argument in seconds:
      double x;
      if ( (!hms_mode) && Helper::str2dbl( lonstr , &x ) ) 
	{
	  if ( x < 0 ) 
	    {
	      logger << "  lights-on time less than 0 -- setting to 0 (EDF start)\n";
	      x = 0;
	    }

	  lights_on = x;
	  
	  logger << "  setting lights_on = " 
		 << lights_on << " secs, "
		 << lights_on/60.0 << " mins from start\n";
	}
      else if ( hms_mode )
	{
	  
	  // else assume is hh:mm:ss
	  
	  if ( ! st.valid ) 
	    Helper::halt( "EDF does not have a valid start time - cannot use lights-on=hh:mm:ss" );
	  
	  clocktime_t et( lonstr );
	  if ( et.valid ) 
	    {
	      // assume that lights-on is always *after* EDF start
	      lights_on = clocktime_t::ordered_difference_seconds( st , et );
	      logger << "  setting lights_on = " << et.as_string() 
		     << " (" << lights_on << " secs, " << lights_on/60.0 << " mins from start)\n";
	    }
	  else
	    logger << "  invalid time for lights-on=" << lonstr << "  -- will ignore this\n";	  
	}
      else
	logger << "  invalid time for lights-on=" << lonstr << "  -- will ignore this\n";	  
    }
  

  //
  // If not already set, see if there are lights_on and/or lights_off annotations present
  //
  
  annot_t * lights_on_annot = timeline->annotations->find( "lights_on" );
  annot_t * lights_off_annot = timeline->annotations->find( "lights_off" );

  //
  // valid combinations:                        | Off                    | On
  // ------------------------------------------------------------------------------------------
  //  a)  lights_off interval only (>epoch dur) : start of lights_off    | end of lights_off
  //  b)  lights_off interval only (<epoch dur) : start of lights_off    | end of recoding (i.e. no explicit lights_on)
  //  c)  two lights_on intervals               : end of first lights_on | start of second light on
  //  d)  one lights_off + 1 or 2 lights_on     : start of lights off    | start of last lights_on

  
  //  nb.  the final case works if two 'change-point' (i.e. 0-duration
  //  intervals) are specified as lights_off and lights_on as well
    
  // lights_off
  //  --> all epochs that end (using exact (not +1) time) before this time will be set to L
  
  // lights_on change-point:
  //  --> all epochs that start on or after this time will be set to L
  
  int n_annot_lights_on = 0;
  int n_annot_lights_off = 0;
  
  if ( lights_off_annot && lights_off < 0 ) // i.e. if not already set above
    {
      annot_map_t & loff = lights_off_annot->interval_events;
      n_annot_lights_off = loff.size();
    }
  
  if ( lights_on_annot && lights_on < 0 ) // i.e. if not already set above
    {
      annot_map_t & lon = lights_on_annot->interval_events;
      n_annot_lights_on = lon.size();
    }

  //
  // extract change-points depending on the available annotations
  //

  // condition d)
  if ( n_annot_lights_off == 1 && n_annot_lights_on > 0 )
    {
      annot_map_t & lon = lights_on_annot->interval_events;
      annot_map_t & loff = lights_off_annot->interval_events;

      annot_map_t::const_iterator aa = loff.begin();
      annot_map_t::const_iterator bb = lon.begin();
      
      // i.e. if leading and trailing LightsOn intervals given, take last
      // otherwise, the assumption is that lights_on is a change-point, so take start of first
      if ( lon.size() == 2 ) 
	++bb;

      lights_off = aa->first.interval.start_sec();
      lights_on = bb->first.interval.start_sec();
    }
  
  // condition a) and b)
  if ( n_annot_lights_off == 1 && n_annot_lights_on == 0 )
    {
      annot_map_t & loff = lights_off_annot->interval_events;
      annot_map_t::const_iterator aa = loff.begin();

      // a) this defines an interval of Lights Off
      lights_off = aa->first.interval.start_sec();
      lights_on = aa->first.interval.stop_sec(); 
      
      // b) if lights_off was of very short duration (i.e. 0-point, or < epoch)
      //    then we assume that it is a change-point marker (at start) but that 
      //    lights_on dees not occur until the end of the recording
      if ( lights_on - lights_off < timeline->epoch_length() )
	lights_on = timeline->last_time_point_tp * globals::tp_duration ; 
      
    }
  
  
  // condition c)
  if ( n_annot_lights_off == 0 && n_annot_lights_on == 2 )
    {
      annot_map_t & lon = lights_on_annot->interval_events;
      annot_map_t::const_iterator aa = lon.begin();
      lights_off = aa->first.interval.stop_sec();
      ++aa;
      lights_on = aa->first.interval.start_sec();
    }

  // check that we did not have any bad combos
  if ( n_annot_lights_off > 1 || n_annot_lights_on > 2 )
    logger <<  "  *** warning - multiple 'lights_off' and 'lights_on' annotations... ignoring\n";

  // check that, if both specified, then lights off happens before lights on

  if ( lights_off > 0 && lights_on > 0 )
    {
      if ( lights_on <= lights_off )
	{
	  logger << "  using lights-off = " << lights_off << " seconds\n"
		 << "        lights-on  = " << lights_on << " seconds\n";
	  Helper::halt( "lights_on must occur after lights_off" );
	}
    }

  
  //
  // Set any epochs to L if they occur before lights off or after lights on
  //
  
  n_lights_fixed = 0;
  n_lights_fixed_was_sleep = 0;
  
  int loff_n = 0 , lon_n = 0;

  if ( lights_off > 0 || lights_on > 0 )
    {
      
      for (int e=0; e < ne_gaps; e++)
	{

	  // ignore GAPS here
	  if ( stages[e] == GAP ) continue;
	  
	  // this epoch /ends/ before lights off?
	  // i.e. will /include/ the epoch when lights_out == epoch start
	  //  this also means we include partial epochs, if lights_out is midway, but so be it
	  
	  if ( lights_off > 0 )
	    {
	      
	      // n.b. fudge to avoid precision issues
	      //const double s = 60 * (e+1) * epoch_mins - 0.0001;
	      const double s = epoch_start[e] + epoch_dur[e] - 0.0001;
	      
	      if ( s < lights_off )
		{
		  if ( is_sleep( stages[e] ) ) ++n_lights_fixed_was_sleep;
		  stages[e] = LIGHTS_ON;
		  ++n_lights_fixed;
		  ++loff_n;
		}
	    }
	  
	  // or, is this epoch /starting/ at or after lights on?
	  if ( lights_on > 0 )
            {

	      // n.b. fudge to avoid precision issues
              //const double s = 60 * e * epoch_mins + 0.0001;
	      const double s = epoch_start[e] + 0.0001;	      

	      if ( s >= lights_on )
                {
                  if ( is_sleep( stages[e] ) )
		    ++n_lights_fixed_was_sleep;
		  stages[e] = LIGHTS_ON;
                  ++n_lights_fixed;
		  ++lon_n;
                }
            }
	}

      
      if ( lights_off > 0 )
	logger << "  set " << loff_n << " leading epochs to L based on a lights_off time of " << lights_off << " seconds from EDF start\n";

      if ( lights_on > 0 )
	logger << "  set " << lon_n << " final epochs to L based on a lights_on time of " << lights_on << " seconds from EDF start\n";

    }
  
  

  //
  // Further hypnogram edits:
  //
  //   1) end-wake/end-sleep : if few sleep epochs but a large gap until
  //   1b) alternative trimming, see below
  //
  //   2) trim-wake / trim-leading-wake / trim-trailing-wake
  //       set any W to ?
  //
  //   3) recode any leading/trailing ? as L (unless whole record is ?) 
  //
  //   4) contrained analyses (e.g. first) 
  
  //
  // Clean up edge cases - if we have a long W period and then only a
  // few sleep epochs set those to missing
  //

  //
  // trim method 1b
  //   -- trie to allow for scenario of   wwwwwwSwwwwwwwwwwwSSwSSwwwwwwwwSSSSSSSSwSSSwSSSSSSSSSSwww
  //   method 1 a would only trim to      wwwwww?wwwwwwwwwwwSSwSSwwwwwwwwSSSSSSSSwSSSwSSSSSSSSSSwww
  //    i.e. if sewcond period of W is long (as may happen in 24 recordings w/ a nap etc) 
  //   we want                            wwwwww?wwwwwwwwwww??w??wwwwwwwwSSSSSSSSwSSSwSSSSSSSSSSwww
  //  method 1b: 
  //    - make contigous blobs of sleep, connecting all periods w/in (e.g.) 10 epochs 
  //    - add border to each blob, e.g. of 10 epochs; these are the 'main sleep periods'
  //          wwwwwwSwwwwwwwwwwwSSwSSwwwwwwwwSSSSSSSSwSSSwSSSSSSSSSSwww
  //              XXXX        XXXXXXXXX    XXXXXXXXXXXXXXXXXXXXXXXXXXX
  //              XXXX        XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX   (drop W < 1hr between S) 
  // but do not allow cut past median sleep epoch
  
  const bool trim1b = param.has( "cut" );
  const bool verbose = param.has( "verbose" );
  
  if ( trim1b ) 
    {
      
      std::vector<double> p;
      if ( ! param.empty( "cut" ) ) p = param.dblvector( "cut" );

      if ( p.size() > 4 )
	Helper::halt( "expecting 0 to 4 params: cut=th,fac,gap,flank" );
      
      const double th      = p.size() > 0 ? p[0] : 50;  // diff (W vs weighted-S)
      const double fac     = p.size() > 1 ? p[1] : 3;    // S vs W weighting  
      const int egap       = p.size() > 2 ? p[2] : 10;   // 5 mins
      const int eflnk      = p.size() > 3 ? p[3] : 10;   // 5 mins
      
      // get median sleep epoch
      int n50 = 0;
      // major sleep period = msp
      std::vector<bool> msp( ne_gaps , false );      
      for (int e=0; e<ne_gaps; e++ )
	{
	  msp[e] = is_sleep( stages[e] );
	  if ( msp[e] ) ++n50;	    
	}
      
      n50 /= 2;
      int ns = 0;
      for (int e=0; e<ne_gaps; e++ )	
	{
	  if ( msp[e] ) ++ns;
	  if ( ns == n50 ) { n50 = e; break; } 
	}



      // extend blobs (from second epoch) based on egap     
      for (int e=1; e<ne_gaps; e++ )
	{
	  if ( msp[e-1] && ! msp[e] )  // S->W transition? 
	    {
	      // do we find 
	      int q = e;
	      bool fill = false;
	      while ( 1 ) { 
		++q;
		if ( q == ne_gaps ) break;
		if ( q - e > egap ) break;
		if ( msp[q] ) { fill = true; break; }
	      }
	      if ( fill ) 
		for (int j=e; j<q; j++) msp[j] = true;
	      // skip ahead
	      e = q;
	    }
	}

      // make a copy so we can ignore flankers when getting the prop W/S
      std::vector<bool> msp0 = msp;
      
      // add left flankers? (based on eflnk)
      for (int e=1; e<ne_gaps; e++ )
	{
	  if ( (!msp[e-1] ) && msp[e] )  // W->S transition? 
	    {
	      int q=e;
	      while (1)  {
		--q;
		if ( q == -1 ) break;
		if ( e - q > eflnk ) break;
		msp[q] = true;
	      }
	    }
	}

      // add right flankers? 
      for (int e=0; e<ne_gaps-1; e++ )
	{	  
	  if ( msp[e] && ! msp[e+1] )  // S->W transition? 
	    {	      
	      int q=e;
	      while (1)  {
		++q;
		if ( q == ne_gaps ) break;
		if ( q - e > eflnk ) break;		
		msp[q] = true;
	      }
	      e = q; // skip ahead	      
	    }
	}
      
      // forward scoring
      std::vector<double> scr( ne_gaps , 0 );
      int s = 0;
      for (int e=0; e<ne_gaps; e++)
	{
	  if ( msp0[e] ) s -= fac; else s++;
	  scr[e] = s;
	}
      
      // reverse scoring
      std::vector<double> rscr( ne_gaps , 0 );
      s = 0;
      for (int e=ne_gaps-1; e>=0; e--)
	{
	  if ( msp0[e] ) s -= fac; else s++;
	  rscr[e] = s;
	}

      // get cut-points (not allow in msp[] or past median sleep epoch #)
      int cut = -1 , rcut = -1;
      double mx = 0;
      double rmx = 0;
      double mx0 = 0 , rmx0 = 0;
      
      for (int e=0; e<ne_gaps; e++)
	{
	  if ( e > n50 ) break;
	  if ( msp[e] ) scr[e] = 0;
	  else if ( scr[e] > mx && scr[e] > 0 && scr[e] >= th ) { mx = scr[e] ; cut = e; }
	  if ( scr[e] > mx0 ) mx0 = scr[e];
	}
      
      for (int e=ne_gaps-1; e>0; e--)
	{
	  if ( e < n50 ) break;
	  if ( msp[e] ) rscr[e] = 0;
	  if ( rscr[e] > rmx && rscr[e] > 0 && rscr[e] >= th ) { rmx = rscr[e] ; rcut = e; }
	  if ( rscr[e] > rmx0 ) rmx0 = rscr[e];
	}
      
      // verbose dump
      if ( verbose )
	{
	  for (int e=0; e<ne_gaps; e++)
	    {
	      std::cout << "tb\t" << e+1 << "\t"
			<< is_sleep( stages[e] ) << "\t"
			<< msp[e] << "\t" 
			<< ( scr[e] > 0 ? scr[e] : 0 ) << "\t"
			<< ( rscr[e] > 0 ? rscr[e] : 0 )<< "\n";	  
	    }
	}

      logger << "  maximum left- and right-scores are " << mx0 << " and " << rmx0 << " (th = " << th << ")\n";
      
      if ( cut != -1 && rcut != -1 )   
	logger << "  placing left-cut at epoch " << cut + 1 << " and right-cut at " << rcut + 1 << "\n";
      else if ( cut != -1  )
	logger << "  placing left-cut at epoch " << cut + 1 << "\n";
      else if ( rcut != - 1 )
	logger << "  placing right-cut at epoch " << rcut + 1 << "\n";
      else
	logger << "  no cuts placed\n";
      
      // check we haven't spliced out everything
      if ( cut != 1 && rcut != -1 )
	if ( cut >= rcut ) logger << "  *** cut isn't leaving any sleep epochs\n";

      if ( cut != -1 ) 
	for (int e=0; e<=cut; e++) stages[e] = UNKNOWN;
      if ( rcut != -1 )
	for (int e=ne_gaps-1; e>= rcut ; e-- ) stages[e] = UNKNOWN;
      			   
    }

  // trimming method 1a
  
  
  const double end_wake = param.has( "end-wake" ) ? param.requires_dbl( "end-wake" ) : ( trim1b ? 0 : 120 ) ; 
  const double end_sleep =  param.has( "end-sleep" ) ? param.requires_dbl( "end-sleep" ) : ( trim1b ? 0 : 5 ) ; 
  
  n_fixed = 0;

  if ( end_wake > 0 )
    {      
      
      // count sleep backwards
      double s = 0;
      std::vector<double> rev_sleep( ne_gaps );
      for (int e=ne_gaps-1; e>= 0; e-- )
	{
	  if ( is_sleep( stages[e] ) ) s += epoch_dur[e];
	  rev_sleep[e] = s;	  
	}
      
      // go forwards counting wake
      n_fixed = 0;
      double cumul_wake = 0; // or missing
      for (int e=0; e<ne_gaps; e++)
	{
	  sleep_stage_t E1 = stages[e];
	  
	  if ( is_sleep( stages[e] ) ) 
	    {
	      if ( cumul_wake > end_wake && rev_sleep[e] < end_sleep )
		{
		  stages[e] = UNKNOWN;
		  ++n_fixed;
		}
	      else
		cumul_wake = 0; // reset the counter
	    }
	  else if ( is_wake( stages[e] ) )
	    cumul_wake += epoch_dur[e] ;
	  
	  // std::cout << " e = " << e << " " << cumul_wake << " " << rev_sleep[e] << "\t" 
	  //  	    << E1 << "\t" << stages[e] << "\t" << n_fixed << "\n";
	}
      
      
      // now do the reverse (i.e. to get rid of spurious leading S epochs)
      //  (with a long wait until 'real' sleep onset)
      
      // count sleep forwards
      s = 0;
      std::vector<double> fwd_sleep( ne );
      for (int e=0; e < ne; e++ )
	{
	  if ( is_sleep( stages[e] ) ) s += epoch_mins;
	  fwd_sleep[e] = s;
	  //std::cout << "fwd_sleep[e] = " << e << " " << fwd_sleep[e] << "\n";
	}
      
      // go backwards counting wake
      
      // reset wake/missing flags
      cumul_wake = 0; 
      for (int e=ne-1; e>=0; e--)
	{
	  sleep_stage_t E1 = stages[e];

	  if ( is_sleep( stages[e] ) ) 
	    {
	      if ( cumul_wake > end_wake && fwd_sleep[e] < end_sleep )
		{
		  // std::cout << " fixing " << cumul_wake << " " << end_wake << " / "
		  // 	    <<  fwd_sleep[e] << " " << end_sleep << "\n";
		  stages[e] = UNKNOWN;
		  ++n_fixed;
		}
	      else
		cumul_wake = 0; // reset the counter
	    }
	  else if ( is_wake( stages[e] ) )
	    cumul_wake += epoch_mins; 
	}
      
      logger << "  set " << n_fixed << " leading/trailing sleep epochs to '?' (given end-wake=" << end_wake 
	     << " and end-sleep=" << end_sleep << ")\n";
    }
  


  //
  // Set all leading and/or trailing wake to ?
  //
  
  // takes optional param in mins
  const bool trim_lead_wake = param.has( "trim-wake" ) || param.has( "trim-leading-wake" );
  const bool trim_trail_wake = param.has( "trim-wake" ) || param.has( "trim-trailing-wake" );
  
  double mins_lead_wake = 0 , mins_trail_wake = 0; 
  
  if ( param.has( "trim-wake" ) )
    {
      if ( param.empty( "trim-wake" ) ) 
	mins_lead_wake = mins_trail_wake = 0;
      else
	mins_lead_wake = mins_trail_wake = param.requires_dbl( "trim-wake" );      
    }
  else
    {
      if ( param.empty("trim-leading-wake" ) )
        mins_lead_wake = 0;
      else
	mins_lead_wake = param.requires_dbl( "trim-leading-wake" );

      if ( param.empty("trim-trailing-wake" ) )
        mins_trail_wake = 0;
      else
	mins_trail_wake = param.requires_dbl( "trim-trailing-wake" );

    }
  
  // allow this much 
  int epoch_lead_wake = mins_lead_wake / epoch_mins ; 
  int epoch_trail_wake = mins_trail_wake / epoch_mins ;

  n_ignore_wake = 0;
	
  if ( trim_lead_wake || trim_trail_wake ) 
    {
      const int ne = stages.size();
      
      n_ignore_wake = 0;
      
      if ( trim_lead_wake ) 
	{
	  int first_sleep = 0;
	  for (int e=0; e<ne; e++)
	    if ( is_sleep( stages[e] ) ) 
	      {  first_sleep = e; break; }
	  
	  first_sleep -= epoch_lead_wake;
	  for (int e=0; e<first_sleep; e++)
	    {
	      stages[e] = UNKNOWN;
	      ++n_ignore_wake;
	    }
	}
      

      if ( trim_trail_wake )
        {
	  int last_sleep = ne-1;
	  for (int e=ne-1; e!=0; e--)
	    if ( is_sleep( stages[e] ) ) 
	      { last_sleep = e; break; } 
	  
	  last_sleep += epoch_trail_wake;
	  
	  for (int e=last_sleep+1; e<ne; e++)
            {
	      stages[e] = UNKNOWN;
              ++n_ignore_wake;
            }
	}
      
      //
      // track & report
      //
      
      logger << "  set " << n_ignore_wake << " ";
      if ( trim_lead_wake && trim_trail_wake ) logger << "leading/trailing";
      else if ( trim_lead_wake ) logger << "leading";
      else logger << "trailing";
      logger << " wake epochs to '?'\n";
      
    }
  

  //
  // Recode any leading/trailing "?" as "L"; an exception to this is if 
  // there are not any wake/sleep epochs found (i.e. all is '?' or 'GAP')
  // which is consistent w/ not having *any* annotations supplied (i.e. just
  // leave as all null in that case, versus change to L (up to GAP or end)
  //
  
  const int ne = stages.size();
  
  bool has_nonnull = false;
  
  for (int e =0; e < ne ; e++)
    {
      if ( stages[e] == WAKE || stages[e] == NREM1 || stages[e] == NREM2 || stages[e] == NREM3 || stages[e] == NREM4 || stages[e] == REM )
	{
	  has_nonnull = true;
	  break;
	}
    }
  
  if ( has_nonnull )
    {
      for (int e =0; e < ne ; e++)
	{
	  if ( stages[e] == UNKNOWN ) stages[e] = LIGHTS_ON;
	  if ( stages[e] != UNKNOWN && stages[e] != LIGHTS_ON ) break;
	}
      
      for (int e = ne - 1 ; e != 0 ; e--)
	{
	  if ( stages[e] == UNKNOWN ) stages[e] = LIGHTS_ON;
	  if ( stages[e] != UNKNOWN && stages[e] != LIGHTS_ON ) break;
	}
    }
  
  
  //
  // constrain to only analyse the first (or last) N minutes after (or before) X ?
  //
  
  constrain_too_short = -1; // N/A

  const bool constrain_first = param.has( "first" );
  const bool constrain_last  = param.has( "last" );
  const bool constrain_clock = param.has( "clock" );
  

  // some checks
  if ( constrain_first && constrain_last )
    Helper::halt( "cannot specify both first and last" );
  if ( constrain_last && constrain_clock )
    Helper::halt( "cannot specifiy both clock and last" );
  if ( param.has( "anchor" ) && constrain_clock )
    Helper::halt( "cannot specifiy both clock and anchor" );  
  if ( constrain_clock && ! constrain_first )
    Helper::halt( "must specify first=<mins> if using clock=<hh:mm:ss> as anchor" );

  constrain_mins = -1 ; 
  if ( constrain_first )
    constrain_mins = param.requires_dbl( "first" );
  else if ( constrain_last )
    constrain_mins = param.requires_dbl( "last" );
  
  
  // T0 : ED start, 
  // T1 : lights out
  // T2 : sleep onset (default)
  //  make everything after N minutes missing '?'

  // or... for T4,T5 T6 do the same, but in reverse
  
  // default to T2 (sleep onset) or final wake (T4)
  constrain_anchor = constrain_last ? "T4" : "T2";
  
  if ( param.has( "anchor" ) )
    constrain_anchor = param.value( "anchor" );


  //
  // Constrain to first N mins 
  //
  
  if ( constrain_first && ! constrain_clock )
    {
      if ( constrain_anchor != "T0" && constrain_anchor != "T1" && constrain_anchor != "T2" ) 
	Helper::halt( "anchor should be T2 (sleep onset, default), T0 (EDF start) or T1 (lights out)" );
      logger << "  restricting statistics to the first " << constrain_mins << " minutes past " ;
      if ( constrain_anchor == "T2" ) logger << "sleep onset";
      else if ( constrain_anchor == "T0" ) logger << "EDF start";
      else logger << "lights out";
      logger << "\n";
      
      const int first_epochs = constrain_mins / (double)epoch_mins ; 
      
      int start = 0; // T0
      
      // T1 ( lights out ) 
      if ( constrain_anchor == "T1" )
	{
	  for (int e=0; e<ne; e++)
	    if ( stages[e] != LIGHTS_ON )
	      {
		start = e;
		break;
	      }
	} // T2 (sleep onset)
      else if ( constrain_anchor == "T2" )
	{
	  for (int e=0; e<ne; e++)
            if ( is_sleep( stages[e] ) )
              {
		start = e;
		break;
              }
	}
      
      // e.g. 1 min = 2 epoch = 0 1 
      // last defined as one past 
      int last = start + first_epochs ;
      
      if ( last > ne ) 
	{	  
	  last = ne;
	  // track actual time allowed
	  constrain_mins = ( last - start ) * epoch_mins;
	  logger << "  *** reducing first period, which is longer than available staging: " 
		 << constrain_mins << " minutes\n";
	  constrain_too_short = 1;
	}
      else 
	constrain_too_short = 0;
          
      logger << "  retaining only epochs " << start+1 << " to " << last << ";"
	     << " setting other epochs to L\n";

      for (int e=0; e<start; e++)
        stages[e] = LIGHTS_ON;

      for (int e=last; e<ne; e++)
	stages[e] = LIGHTS_ON;
      
    }


  //
  // Constrain to first N mins 
  //
  
  if ( constrain_last )
    {
      if ( constrain_anchor != "T4" && constrain_anchor != "T5" && constrain_anchor != "T6" ) 
	Helper::halt( "anchor should be T4 (final wake, default), T5 (lights on) or T6 (EDF stop)" );
      logger << "  restricting statistics to the last " << constrain_mins << " minutes prior to " ;
      if ( constrain_anchor == "T4" ) logger << "final wake";
      else if ( constrain_anchor == "T6" ) logger << "EDF stop";
      else logger << "lights on";
      logger << "\n";
      
      const int last_epochs = constrain_mins / (double)epoch_mins ; 
      
      int last = ne; // T6, 1-past encoding
      
      // T5 ( lights on ) 
      if ( constrain_anchor == "T5" )
	{
	  for (int e=ne-1; e>=0; e--)
	    if ( stages[e] != LIGHTS_ON )
	      {
		last = e+1;
		break;
	      }
	} // T4 
      else if ( constrain_anchor == "T4" )
	{
	  for (int e=ne-1; e>=0; e--)
            if ( is_sleep( stages[e] ) )
              {
		last = e+1;
		break;
              }
	}
      
      // last defined as one past
      int start = last - last_epochs ; 
      
      if ( start < 0 )
	{	  
	  start = 0;
	  // track actual time allowed
	  constrain_mins = ( last - start ) * epoch_mins;
	  logger << "  *** reducing last period, which is longer than available staging: " 
		 << constrain_mins << " minutes\n";
	  constrain_too_short = 1;
	}
      else 
	constrain_too_short = 0;
      
      logger << "  retaining only epochs " << start+1 << " to " << last << ";"
	     << " setting other epochs to L\n";
            
      for (int e=0; e<start; e++)
	stages[e] = LIGHTS_ON;
      
      for (int e=last; e<ne; e++)
        stages[e] = LIGHTS_ON;
      
      
    }

  
  //
  // Constrain to clock start time 
  //   NOTE - this allows for windows that do not overlap the EDF... they will effectively be L so window
  //          may be smaller... this will impact MINS estimates... i.e. careful to interpret W (are L effectively = W) 
  //

  if ( constrain_clock )
    {

      constrain_too_short = 0;
      
      if ( param.empty( "clock" ) )
	Helper::halt( "clock requires a hh:mm:ss value" );
      
      std::string ctime = param.value( "clock" );
          
      clocktime_t starttime( timeline->edf->header.starttime );
      bool invalid_hms = ! starttime.valid;
      if ( invalid_hms )
	Helper::halt( "EDF does not have valid start-time in header" );

      // note::: ctime might come before or after EDF start time
      clocktime_t t1( ctime );
      
      // 1: starttime is first, 2: ctime is first; 0: both same
      // if days not specified, then assume shortest distance (12 hr) rule
      // to determine which comes 'first'
      
      const int e = clocktime_t::earlier( starttime , t1 );

      const bool before_start = e == 2 ; 
      
      // if EDF starts *after* this time, flag that we'll have a short duration
      if ( before_start ) constrain_too_short = 1;
      
      // get time diff either way::
      const double secs = clocktime_t::ordered_difference_seconds( starttime , t1 ) ; 
      
      // note, start/stop_sec may be -ve if they come before the EDF start
      // we don't need to do anything there, as it is like having L anyway

      const double start_sec = before_start ?
	-1 * clocktime_t::ordered_difference_seconds( t1 , starttime ) : 
	+1 * clocktime_t::ordered_difference_seconds( starttime , t1 ) ;
      
      // get start epoch
      int start_epoch = 0;
      
      double esecs = 0; 
      
      for (int e=0;e<ne;e++)
	{
	  const double diff = fabs( esecs - start_sec );
	  
	  if ( esecs > start_sec || diff < 1e-4 )
	    {
	      start_epoch = e;
	      break;
	    }
	  
	  // advance to start of next epoch
	  esecs += epoch_dur[e];
	}

      // we now have the start epoch
      //   if the clock time was before the EDF start (0)
      //   then we need to adjust for that - i.e. will have a smaller
      //   window
      
      const double adjusted_mins = start_sec >= 0 ? constrain_mins : constrain_mins + ( start_sec / 60.0 ) ;  
      
      // get number of (whole) epochs to include
      // i.e. always try to include same # of epochs
      // even if epoch starts are not aligned w/ clock time
      // specified 
      
      const int included_epochs = adjusted_mins / (double) epoch_mins ; 
      
      const int last_epoch = start_epoch + included_epochs ; // 1-past encoding

      // std::cout << start_sec << " start_sec \n";
      // std::cout << constrain_mins << " constrain_mins \n";
      // std::cout << adjusted_mins << " adjusted_mins \n";
      // std::cout << included_epochs << " included_epochs (win size) \n";
      // std::cout << last_epoch << " last epoch (1 past)\n";
	    
      // special case: wipe all, if clock period starts after last
      // epoch, or if clock period ends before first epoch
      
      if ( included_epochs <= 0 || start_epoch >= ne )
	{
	  constrain_too_short = 1;
	  logger << "  clock period is outside observed epoch range: setting all epochs to L\n";
	  for (int e=0; e<ne; e++)
	    stages[e] = LIGHTS_ON;	  
	}
      
      // wipe all epochs before the start epoch 
      if ( start_epoch > 0 )
	{
          logger << "  setting " << start_epoch << " epochs before " << t1.as_string() << " to L\n";
          for (int e=0; e<start_epoch; e++)
            stages[e] = LIGHTS_ON;
	}

      // wipe all epochs after the end
      if ( included_epochs > 0 )
	{
	  int actual_stop_epoch = last_epoch;
	  
	  if ( last_epoch > ne )
	    {
	      actual_stop_epoch = ne;
	      constrain_too_short = 1;
	    }

	  if ( ne - actual_stop_epoch > 0 )
	    {
	      logger << "  setting " << ne - actual_stop_epoch << " end epochs to L\n";
	      for (int e=actual_stop_epoch; e<ne; e++)
		stages[e] = LIGHTS_ON;
	    }
	}
    }  


  //
  // Specify sliding window on epoch-level data?
  //

  sliding_window = false;

  if ( param.has( "slide" ) )
    {
    
      sliding_window = true;
      
      std::vector<double> s = param.dblvector( "slide" ); // mins
      
      const int nslide = s.size() / 2 ;
      if ( s.size() != nslide * 2 ) 
	Helper::halt( "slide=<width1>,<inc1>,<width2>,<inc2>,... must have even number of args (width,inc pairs in mins)" );
      
      window_width_epochs.clear();
      window_increment_epochs.clear();
      
      for (int i=0; i<nslide; i+=2)
	{
	  window_width_epochs.push_back( s[i] / epoch_mins );
	  window_increment_epochs.push_back( s[i+1] / epoch_mins );
	}


      // by default, anchor T2 unless told otherwise
      std::vector<int> slide_anchor( nslide , 2 );
      
      if ( param.has( "slide-anchor" ) )
	{
	  slide_anchor = param.intvector( "slide-anchor" );
	  if ( slide_anchor.size() != nslide ) 
	    Helper::halt( "slide-anchor arg len does not match given slide" );
	  for (int i=0; i<nslide; i++)
	    if ( slide_anchor[i] < 0 || slide_anchor[i] > 6 ) 
	      Helper::halt( "slide-anchor args should be between 0 and 6" );	  
	}
      
      window_anchor = slide_anchor;

    }

  //
  // Output ascending/descending N2? 
  //

  output_n2_asc_desc = param.has( "n2-asc-desc" ) ? param.yesno( "n2-asc-desc" ) : false ;
  n2_asc_desc_th = param.has( "n2-asc-desc-th" ) ? param.requires_dbl( "n2-asc-desc-th" ) : 0.25;

}

void hypnogram_t::calc_stats( const bool verbose )
{

  //
  // epoch size (in minutes) and number
  //

  const double epoch_mins = timeline->epoch_length() / 60.0 ; 

  // this includes gaps
  const int ne = stages.size();

  //
  // Basic summary statistics per-individual/night
  //

  // clear, in case this is run twice
  mins[ "W" ] = mins[ "N1" ] = mins[ "N2" ] = mins[ "N3" ] = mins[ "N4" ] = mins[ "R" ] = mins[ "?" ] = mins["L"] = 0;
  mins[ "N2A" ] = mins[ "N2D" ] = mins[ "N2F" ] = 0;
  
  // implicitly, this will only count in the TRT (i.e. ignore pre
  // lights-out, and post lights-on)

  // nb. this includes gaps, so count below
  for (int e = 0 ; e < ne ; e++ )
    {
      if ( epoch_gap[e] ) mins["GAP"] += epoch_dur[e]; // actual GAP size
      else if ( stages[e] == WAKE  ) mins[ "W" ] += epoch_mins;
      else if ( stages[e] == NREM1 ) mins[ "N1" ] += epoch_mins;
      else if ( stages[e] == NREM2 ) mins[ "N2" ] += epoch_mins;
      else if ( stages[e] == NREM3 ) mins[ "N3" ] += epoch_mins;
      else if ( stages[e] == NREM4 ) mins[ "N4" ] += epoch_mins;
      else if ( stages[e] == REM   ) mins[ "R" ] += epoch_mins;
      else if ( stages[e] == LIGHTS_ON ) mins[ "L" ] += epoch_mins;
           
      else mins[ "?" ] += epoch_mins; // movement, artifact, unscored or GAP 
    }

  // did we observe /any/ sleep?

  any_sleep = ( mins[ "N1" ] + mins["N2"] + mins["N3"] + mins["N4"] + mins["R"] ) > 0 ;


  // lights out/on check
  // i.e. L can only be at start and end of recording
  // so illegal to have L flanked by non-L on both sides : 010
  bool lights_back_on = false;
  for (int e=1;e<ne-1;e++) 
    {
      // went from lights off to on again 
      if ( (!epoch_gap[e-1])
	   && stages[e-1] != LIGHTS_ON
	   && stages[e]   == LIGHTS_ON ) lights_back_on = true;

      // if we previously have gone from OFF->ON, we cannot go back
      // OFF again--- i.e. only allow a single LIGHTS OFF interval
      
      if ( lights_back_on
	   && stages[e] == LIGHTS_ON
	   && stages[e+1] != LIGHTS_ON
	   && (!epoch_gap[e+1] ) )
	Helper::halt( "LIGHTS_ON periods can only be at start and end of recording" );
    }

  // lights out/on
  lights_out_epoch = 0;
  for (int e=0;e<ne;e++) 
    if ( stages[e] != LIGHTS_ON ) { lights_out_epoch = e; break; } 

  lights_on_epoch = ne; // by default, one past the end  
  for (int e=ne-1;e>=0;e--) 
    if ( stages[e] != LIGHTS_ON ) { lights_on_epoch = e+1 ;  break; }

  //  std::cout << "LON/OFF " << lights_out_epoch << " " << lights_on_epoch << "\n";
  
  any_dark = false; // i.e. otherwise implies all L all night
  for (int e=0;e<ne;e++)
    {
      if ( stages[e] != LIGHTS_ON ) 
	{
	  any_dark = true;
	  break;
	}
    }
  
  //
  // For cycle calculations below, etc,  use TRT which is [ lights_out_epoch ,  lights_on_epoch )  
  // rather than [ 0 , ne ) 
  //

  //
  // First wake epoch of final bout of wake (i.e. so this can be subtracted off WASO)
  //

  // can't occur after lights on
  final_wake_epoch = lights_on_epoch; // i.e. is defined as one past end
  for (int e = lights_on_epoch-1 ; e >= 0 ; e-- )
    if ( stages[e] != WAKE ) { final_wake_epoch = e+1; break; }
  
  //
  // first REM epoch
  //

  int first_rem_epoch = ne;
  for (int e = 0 ; e < ne ; e++ )
    if ( stages[e] == REM ) { first_rem_epoch = e; break; }


  // requires non-missing SLEEP
  // persistent sleep defined as 10 mins
  int lps_required = 10.0 / (double)epoch_mins;
  int cnt_sleep = 0;

  bool found_first_sleep = false;
  first_sleep_epoch = ne;
  first_persistent_sleep_epoch = ne;

  for (int e = 0 ; e < ne ; e++ )
    {
      if ( is_sleep( stages[e] ) )
	{
	  if ( ! found_first_sleep )
	    {
	      first_sleep_epoch = e; 
	      found_first_sleep = true;
	    }
	  
	  // LPS start?
	  bool lps = true;  
	  for (int e2=e; e2 < e + lps_required; e2++ )
	    {	      
	      if ( e2 >= ne ) { lps=false; break; }
	      if ( ! is_sleep( stages[e2] ) ) 
		{ lps = false; break; }
	    }
	  
	  if ( lps )
	    {
	      first_persistent_sleep_epoch = e;
	      break;
	    }
	}
    }
  
  // last epoch of sleep
  int last_sleep_epoch = 0;
  for (int e = ne - 1; e != 0 ; e-- )
    if ( is_sleep( stages[e] ) ) { last_sleep_epoch = e; break; } 
  
  // Store midpoint as an epoch code  
  sleep_midpoint_epoch = first_sleep_epoch + ( final_wake_epoch - first_sleep_epoch ) / 2;

  // total time in bed
  TIB = ne * epoch_mins;
  
  // total recording time (i.e. only from lights out, lights on)
  // note; lights_out_epoch is defined as 1 past end, so no +1
  int TRT_total_epochs = any_dark ? lights_on_epoch - lights_out_epoch : 0 ; 
  TRT =  TRT_total_epochs * epoch_mins;
  
  // total wake time (ignores pre lights out, post lights off)
  TWT = mins["W"];
  
  // final wake time (until lights on)
  FWT = ( lights_on_epoch - final_wake_epoch ) * epoch_mins; 

  // post-sleep time (until recording end)    
  t_postsleep = ( ne - final_wake_epoch ) * epoch_mins;
  
  // REM latency
  rem_lat_mins = ( first_rem_epoch - first_sleep_epoch ) * epoch_mins;
  
  // REM latency, excluding W
  rem_lat_nowake_mins = 0;
  for (int e = first_sleep_epoch ; e <= first_rem_epoch; e++)
    if ( is_nrem( stages[e] ) ) rem_lat_nowake_mins += epoch_mins;
  
  // Total sleep time (excludes 'other')
  TST = TRT - TWT - mins[ "?" ] ;
  
  // Total GAP time
  TGT = mins[ "GAP" ];
  
  // study starts/ends in sleep?
  starts_in_sleep = is_sleep( stages[0] );
  ends_in_sleep = is_sleep( stages[ne-1] );
  
  // std::cout << "TST TRT " << TST << " " << TRT << "\n";
  // std::cout << "dets" << lights_out_epoch << " " << first_sleep_epoch << " " << last_sleep_epoch << " " << lights_on_epoch << "\n";
  
  // sleep latency
  slp_lat = ( first_sleep_epoch - lights_out_epoch ) * epoch_mins;

  // pre-sleep time (includes lights off )
  t_presleep = first_sleep_epoch * epoch_mins ; 
  
  // latency to persistent sleep
  per_slp_lat = ( first_persistent_sleep_epoch - lights_out_epoch ) * epoch_mins;
  
  // Sleep period time : note, diff. from Luna output
  // here SPT = sleep onset to lights On (i.e. includes final wake)
  SPT = TRT - slp_lat;
  
  //  std::cout << "TWT , slp_lat , FWT " << TWT << " " << slp_lat << " " << FWT << "\n";

  // WASO (ignores leading and also trailing wake)
  //  TWT = SLP_LAT + WASO + FWT
  //  WASO = TWT - slp_lat - FWT;
  // BUT, the above might count UNDEF/OTHER in the leading/trailing wake, and so remove too much
  // easier to just figure it out by iteration
  int w = 0;
  for (int e=first_sleep_epoch;e<=last_sleep_epoch;e++)
    if ( stages[e] == WAKE ) ++w;
  WASO = w * epoch_mins;
  mins[ "WASO" ] = WASO;
  
  // SE - sleep efficiency (includes sleep latency as W) include OTHER in denom
  slp_eff_pct = ( TST / TRT ) * 100;
    
  // IGNORED *** now; sleep maintainence (ignores initial sleep latency as W) includes OTHER in denom
  slp_main_pct = ( TST / SPT ) * 100;

  // note - actual (reported) SPT is actually SPT - FWT (i.e. given to user);
  // this should equal SPT_actual as calculated here... used below just in the PCT2 calcs
  const double SPT_actual = epoch_mins * ( last_sleep_epoch - first_sleep_epoch + 1 ) ;    

  // SME -- alternate: sleep maintainence/efficiency 2 (denom is from initial sleep to final sleep)
  // i.e. ignores *both* leading and trailing W; includes OTHER in denom
  slp_eff2_pct = ( TST / SPT_actual ) * 100;


  //
  // Ascending/descending N2 
  //
  
  // Assign each N2 epoch a value to describe its relative position in
  // the hypnogram: consider 'k' epochs before and after, calculating 
  // the average of non-N2 epochs,  either W/R/N1 or N3/N4  

  // Calculate a single score: 
  // Defined only for N2 epochs
  //   Left epochs    +1  N3              Right epochs    -1  N3
  //                  -1  N1/W/R                          +1  N1/W/R
  
  const int n2_ascdesc_k = 10;  // 5 minutes
  
  // Can select extreme epochs as 'ascending' and 'descending' based on this score 
  // (e.g. >+0.5 and < -0.5)
  
  n2_ascdesc.resize( ne , 0 );
  
  for (int e=0;e<ne;e++)
    {
      
      if ( stages[e] != NREM2 ) continue;
      
      double left_wgt = 0;
      int left_n = 0; 
      int k = e-1;

      while ( k >= 0 )
	{

	  if ( stages[k] == NREM3 || stages[k] == NREM4 )
	    {
	      left_wgt += 1;
	      ++left_n;
	    }
	  
	  if ( stages[k] == NREM1 || stages[k] == REM || stages[k] == WAKE )
	    {
	      left_wgt += -1;
	      ++left_n;
	    }

	  // counted enough?
	  if ( left_n > n2_ascdesc_k ) break;
	  
	  // next left epoch
	  --k;
	}

      //
      // Right-most
      //

      double right_wgt = 0;
      int right_n = 0; 
      k = e+1;

      while ( k < ne )
	{

	  if ( stages[k] == NREM3 || stages[k] == NREM4 )
	    {
	      right_wgt += -1;
	      ++right_n;
	    }
	  
	  if ( stages[k] == NREM1 || stages[k] == REM || stages[k] == WAKE )
	    {
	      right_wgt += +1;
	      ++right_n;
	    }
	  
	  // counted enough?
	  if ( right_n > n2_ascdesc_k ) break;
	  
	  // next right epoch
	  ++k;
	}

      // std::cout << " left_wgt " << left_wgt << " " << left_n << "\n";
      // std::cout << " right_wgt " << right_wgt << " " << right_n << "\n\n";
      
      if ( left_n  > 0 ) left_wgt /= (double)left_n;
      if ( right_n > 0 ) right_wgt /= (double)right_n;

      // simple average of left/right averages
      // if no data, wgt will be 0, which is fine
      n2_ascdesc[e] = ( left_wgt + right_wgt ) / 2.0;
      //std::cout << " n2 AD = " << e << "  " << stages[e] << " " << n2_ascdesc[e] << "\n";
    }



  //
  // Track duration of N2 class
  //

  mins[ "N2A" ] = 0 ;
  mins[ "N2D" ] = 0 ;
  mins[ "N2F" ] = 0 ;

  std::vector<bool> n2asc( ne , false );
  std::vector<bool> n2dsc( ne , false );
  std::vector<bool> n2flt( ne , false );

  if ( output_n2_asc_desc )
    {
      for (int e=0;e<ne;e++)
	{
	  if ( stages[e] != NREM2 ) continue;
	  
	  // if      ( n2_ascdesc[e] >= 0.25 ) mins[ "N2_ASC" ] += epoch_mins;
	  // else if ( n2_ascdesc[e] <= -0.25 ) mins[ "N2_DSC" ] += epoch_mins;
	  // else mins[ "N2_FLT" ] += epoch_mins;      
	  
	  // default = 0.25
	  if      ( n2_ascdesc[e] >= n2_asc_desc_th ) 
	    {
	      mins[ "N2A" ] += epoch_mins;
	      n2asc[e] = true;
	    }
	  else if ( n2_ascdesc[e] <= -n2_asc_desc_th ) 
	    {
	      mins[ "N2D" ] += epoch_mins;
	      n2dsc[e] = true;
	    }
	  else 
	    {
	      mins[ "N2F" ] += epoch_mins;      
	      n2flt[e] = true;
	    }
	}
    }



  //
  // relative (PCT and DENS) calcs
  //

  if ( TST > 0 ) 
    {
      pct["N1"]  = mins[ "N1" ] / TST;
      pct["N2"]  = mins[ "N2" ] / TST;
      pct["N3"]  = mins[ "N3" ] / TST;
      pct["N4"]  = mins[ "N4" ] / TST;
      pct["R"]   = mins[ "R" ] / TST;
      
      // special cases
      pct["N2A"]  = mins[ "N2" ] > 0 ? mins["N2A" ] /  mins[ "N2" ]  : -1 ;
      pct["N2D"]  = mins[ "N2" ] > 0 ? mins["N2D" ] /  mins[ "N2" ]  : -1 ;
      pct["N2F"]  = mins[ "N2" ] > 0 ? mins["N2F" ] /  mins[ "N2" ]  : -1 ;
    }
  else
    {
      pct[ "N1" ] = 0;
      pct[ "N2" ] = 0;
      pct[ "N3" ] = 0;
      pct[ "N4" ] = 0;
      pct[ "R" ] = 0;      

      pct[ "N2A" ] = 0;
      pct[ "N2D" ] = 0;
      pct[ "N2F" ] = 0;

    }

  if ( SPT_actual > 0 ) 
    {
      pct2["N1"]  = mins[ "N1" ] / SPT_actual;
      pct2["N2"]  = mins[ "N2" ] / SPT_actual;
      pct2["N3"]  = mins[ "N3" ] / SPT_actual;
      pct2["N4"]  = mins[ "N4" ] / SPT_actual;
      pct2["R"]   = mins[ "R" ] / SPT_actual;

      pct2["N2A"]  = mins[ "N2A" ] / SPT_actual;
      pct2["N2D"]  = mins[ "N2D" ] / SPT_actual;
      pct2["N2F"]  = mins[ "N2F" ] / SPT_actual;
      
    } 
  else
    {
      pct2[ "N1" ] = 0;
      pct2[ "N2" ] = 0;
      pct2[ "N3" ] = 0;
      pct2[ "N4" ] = 0;
      pct2[ "R" ] = 0;

      pct2[ "N2A" ] = 0;
      pct2[ "N2D" ] = 0;
      pct2[ "N2F" ] = 0;

    }


  //
  // Runs test on stages;;; collapse 
  //

  if ( 0 ) 
    {

      std::vector<std::string> runs_stage5, runs_stage3;
      
      // ignore L
      for (int e=0;e<ne;e++)
	{
	  if ( is_rem( stages[e] ) )
	    {
	      runs_stage5.push_back( "R" );
	      runs_stage3.push_back( "R" );
	    }
	  else if ( is_wake( stages[e] ) ) 
	    {
	      runs_stage5.push_back( "W" );
	      runs_stage3.push_back( "W" );
	    }
	  else if ( is_nrem1( stages[e] ) ) 
	    {
	      runs_stage5.push_back( "N1" );
	      runs_stage3.push_back( "NR" );
	    }
	  else if ( is_nrem2( stages[e] ) ) 
	    {
	      runs_stage5.push_back( "N2" );
	      runs_stage3.push_back( "NR" );
	    }
	  else if ( is_nrem34( stages[e] ) ) 
	    {
	      runs_stage5.push_back( "N3" );
	      runs_stage3.push_back( "NR" );
	    }
	}
      
      
      runs_pv5 = Statistics::runs_test( runs_stage5 );
      
      runs_pv3 = Statistics::runs_test( runs_stage3 );
    }

  //
  // Sliding window for stage density?
  //

  if ( sliding_window ) do_slide();

  //
  // Bout count and duration (ignore N4 here.)
  // Also, make an output of each bout
  //

  const std::vector<std::string> ar1 = { "N1", "N2", "N2A","N2D","N2F", "N3", "N4", "NR", "R", "S", "W" , "?" , "L" , "WASO" } ;
  const std::vector<std::string> ar2 = { "N1", "N2",                    "N3", "N4", "NR", "R", "S", "W" , "?" , "L" , "WASO" } ;
  
  const std::vector<std::string> these_stages = output_n2_asc_desc ? ar1 : ar2 ; 

  // initialize these two maps, which are +='ed below
  // i.e. (if HYPNO is run twice)
  bout_5.clear();
  bout_10.clear();


  //
  // Get bout stats
  //



  std::vector<std::string>::const_iterator qq = these_stages.begin();
  while ( qq != these_stages.end() )
    {
      sleep_stage_t stage = WAKE;
      if ( *qq == "N1" ) stage = NREM1;
      if ( *qq == "N2" ) stage = NREM2;
      if ( *qq == "N3" ) stage = NREM3; 
      if ( *qq == "N4" ) stage = NREM4; 
      if ( *qq == "R" ) stage = REM;
      if ( *qq == "?" ) stage = UNKNOWN;
      if ( *qq == "L" ) stage = LIGHTS_ON;
      
      // special cases
      bool all_nrem = *qq == "NR";      
      bool all_sleep = *qq == "S";
      bool waso = *qq == "WASO";
      bool n2a = *qq == "N2A";
      bool n2d = *qq == "N2D";
      bool n2f = *qq == "N2F";


      std::vector<double> b;
      for (int e=0; e<ne; e++)
	{
	  // start of a bout?
	  bool bout_start = false;
	  if ( all_nrem )
	    bout_start = stages[e] == NREM1 || stages[e] == NREM2 || stages[e] == NREM3 || stages[e] == NREM4 ;	      	  
	  else if ( all_sleep ) 
	    bout_start = stages[e] == NREM1 || stages[e] == NREM2 || stages[e] == NREM3 || stages[e] == NREM4 || stages[e] == REM ;
	  else if ( waso )
	    bout_start = stages[e] == WAKE && e >= first_sleep_epoch && e <= last_sleep_epoch ;
	  else if ( n2a ) 
	    bout_start = n2asc[e];
	  else if ( n2d ) 
	    bout_start = n2dsc[e];
	  else if ( n2f ) 
	    bout_start = n2flt[e];
	  else
	    bout_start = stages[e] == stage;

	  // skip to next epoch, not a match
	  if ( ! bout_start ) continue;
	  
	  // assess this bout -- look ahead as many matching epochs as possible (or to end)
	  double l = epoch_mins;
	  while ( 1 )
	    {
	      ++e;

	      // end of recording?... add this final bout
	      if ( e == ne )
		{
		  b.push_back( l );
		  break;
		}

	      // end of bout?
	      if ( all_nrem )
		{
		  if ( stages[e] != NREM1 && stages[e] != NREM2 && stages[e] != NREM3 && stages[e] != NREM4 )
		    {
		      b.push_back( l );
                      break;
		    }
		}
	      else if ( all_sleep ) 
		{
		  if ( stages[e] != NREM1 && stages[e] != NREM2 && stages[e] != NREM3 && stages[e] != NREM4 && stages[e] != REM )
                    {
                      b.push_back( l );
                      break;
                    }
		}
	      else if ( waso ) 
		{
		  if ( stages[e] != WAKE || e > last_sleep_epoch )
		    {
		      b.push_back( l );
                      break;
		    }		  
		}
	      else if ( n2a ) 
		{
		  if ( ! n2asc[e] )
                    {
                      b.push_back( l );
                      break;
                    }
		}
	      else if ( n2d ) 
		{
		  if ( ! n2dsc[e] )
                    {
                      b.push_back( l );
                      break;
                    }
		}
	      else if ( n2f ) 
		{
		  if ( ! n2flt[e] )
                    {
                      b.push_back( l );
                      break;
                    }
		}
	      else
		{
		  if ( stages[e] != stage )
		    {
		      b.push_back( l );
		      break;
		    }
		}

	      // else increase bout length
	      l += epoch_mins;
	    }
	  
	  // now continue to next epoch
	}

      //
      // record stats
      //
      
      for (int bb=0; bb<b.size(); bb++)
	{
	  if ( b[bb] >= 5 ) bout_5[ *qq ] += b[bb];
	  if ( b[bb] >= 10 ) bout_10[ *qq ] += b[bb];	  
	}
      
      bout_n[ *qq ] = b.size();
      if ( b.size() != 0 )
	{
	  bout_med[ *qq ] = MiscMath::median( b , true ); // T -> handle ties properly
	  bout_mean[ *qq ] = MiscMath::mean( b );
	  bout_max[ *qq ] = MiscMath::max( b );
	}
      
      //
      // next stage
      //
      
      ++qq;
    }


  //
  // Bouts
  //
  
  bouts.clear();

  // use bout_t to get NREM downcasting
  bout_t curr( 0 , 0 , stages[0] );
  int bstart = 0;
  
  for (int e=1; e<=ne; e++)
    {
      
      if ( e == ne )
	{
	  bouts.insert( bout_t( bstart , ne-1 , curr.ss ) );
	  break;
	}

      bout_t next( 0 , 0 , stages[e] );

      if ( next.ss != curr.ss )
	{
	  bouts.insert( bout_t( bstart , e-1 , curr.ss ) );
	  curr.ss = next.ss;
	  bstart = e;
	}
    }


  
  //
  // Elapsed stage duration/time values
  //   for a given stage/class: 
  //     **anchored on first sleep**, get all epochs (or all sleep epochs) and normalize to 0..1 
  //     calculate median 
  //

  if ( found_first_sleep )
    {
      
      // first_sleep_epoch = e;

      // use a reduced set of stages types here, i.e no need for L or ? etc
      // i.e. only things within the sleep period (so no W either)

      const std::vector<std::string> ar1 =  { "N1", "N2", "N2A", "N2D", "N2F", "N3", "NR", "R", "S", "WASO" };
      const std::vector<std::string> ar2 =  { "N1", "N2",                      "N3", "NR", "R", "S", "WASO" };
      
      const std::vector<std::string> these_stages = output_n2_asc_desc ? ar1 : ar2 ; 
      
      std::vector<std::string>::const_iterator qq = these_stages.begin();
      while ( qq != these_stages.end() )
	{
	  sleep_stage_t stage = WAKE;
	  if ( *qq == "N1" ) stage = NREM1;
	  if ( *qq == "N2" ) stage = NREM2;
	  if ( *qq == "N3" ) stage = NREM3; 
	  if ( *qq == "N4" ) stage = NREM4; 
	  if ( *qq == "R" ) stage = REM;
	  if ( *qq == "?" ) stage = UNKNOWN;
	  if ( *qq == "L" ) stage = LIGHTS_ON;
	  
	  // special cases
	  bool all_nrem = *qq == "NR";      
	  bool all_sleep = *qq == "S";
	  bool waso = *qq == "WASO";
	  bool n2a = *qq == "N2A";
	  bool n2d = *qq == "N2D";
	  bool n2f = *qq == "N2F";

	  // track time from EDF start
	  double t_all = 0;
	  double t_sleep = 0;
	  
	  // elapsed seconds of this stage
	  std::vector<double> elapsed_all, elapsed_sleep;
	  std::vector<int> chk;
	  

	  //
	  // scale to 0..1 ;; 0 = sleep onset, and t_all, t_sleep track max
	  //   -- always used sleep onset as 0 (min);
	  //      but first point will have min as dur of first sleep epoch dur, 
	  //      so that is the min possible value, if we want scaling to 0..1, not 0 per se.
	  
	  const double min_all = epoch_dur[ first_sleep_epoch ];
	  const double min_sleep = epoch_dur[ first_sleep_epoch ];
	  
	  double max_all = 0;
	  double max_sleep = 0;
	  
	  // iterate over epochs, *** starting at first sleep epoch
	  for (int e=first_sleep_epoch; e<ne; e++)
	    {

	      // only consider within the sleep period
	      if ( e == final_wake_epoch ) break;

	      // track from EDF start - note
	      // add this first, as makes sense that we achieve this at the 'end'
	      // of the given epoch that adds it in. i.e. if 1 min, then at end of 2nd epoch
	      // would say 'yes, we've now had 1 min of that stage, not at start
	      
	      t_all += epoch_dur[e];
	      
	      if ( stages[e] == NREM1 || stages[e] == NREM2 || stages[e] == NREM3 || stages[e] == NREM4 || stages[e] == REM )
		t_sleep += epoch_dur[e];
	      
	      // is this the last sleep epoch?
	      
	      if ( e + 1 == final_wake_epoch ) 
		{
		  max_all = t_all;
		  max_sleep = t_sleep;
		}
	      
	      //
	      // is this a match
	      //

	      bool matches = false;
	      
	      if ( all_nrem )
		matches = stages[e] == NREM1 || stages[e] == NREM2 || stages[e] == NREM3 || stages[e] == NREM4 ;
	      else if ( all_sleep ) 
		matches = stages[e] == NREM1 || stages[e] == NREM2 || stages[e] == NREM3 || stages[e] == NREM4 || stages[e] == REM ;
	      else if ( waso )
		matches = stages[e] == WAKE && e >= first_sleep_epoch && e <= last_sleep_epoch ;
	      else if ( n2a ) 
		matches = n2asc[e] ; 
	      else if ( n2d ) 
		matches = n2dsc[e] ; 
	      else if ( n2f )
		matches = n2flt[e]; 
	      else
		matches = stages[e] == stage;
	      	      
	      if ( matches )
		{
		  // add in
		  elapsed_all.push_back( t_all );
		  elapsed_sleep.push_back( t_sleep );
		  chk.push_back( e+1 );
		}

	    } // next epoch
	  
	  const int n1 = elapsed_all.size();
	  
	  // requires at least 5 epochs 
	  if ( n1 >= 5 ) 
	    {
	      
	      for (int i=0; i<n1; i++)
		{
		  // double xx = elapsed_sleep[i];
		  // double yy = elapsed_all[i] ;
		  
		  elapsed_all[i] = ( elapsed_all[i] - min_all ) / ( max_all - min_all );
		  elapsed_sleep[i] = ( elapsed_sleep[i] - min_sleep ) / ( max_sleep - min_sleep );

		  // if ( all_nrem ) 
		  //   std::cout << " NR elapsed_sleep[i] = " << i << "\t" << xx << "\t" << yy << "\t" << elapsed_sleep[i] << "\t" << min_sleep << "\t" << max_sleep << "\n";
		  // if ( stage == REM )
		  //   std::cout << " R elapsed_sleep[i] = " << i << "\t" << xx << "\t" << yy << "\t" << elapsed_sleep[i] << "\t" << min_sleep << "\t" << max_sleep << "\n";
		  
		}

	      // get medians
	      
	      const double med_all = MiscMath::median( elapsed_all );
	      const double med_sleep = MiscMath::median( elapsed_sleep );


	      
	      // const double med_all = MiscMath::mean( elapsed_all );
	      // const double med_sleep = MiscMath::mean( elapsed_sleep );
	      
	      //std::cout << " stats med = " << med_all << " " << med_sleep << "\n";

	      
	      // store results
	      stg_timing_all[ *qq ] = med_all;
	      
	      if ( stage == NREM1 || stage == NREM2 || stage == NREM3 || stage == NREM4 || stage == REM || n2a || n2d || n2f || all_nrem )
		stg_timing_sleep[ *qq ] = med_sleep;
	      
	    }
	  
	  // next stage/class
	  ++qq;
	}
    }
  
   
  //
  // Sleep cycles : based on modified Floyd & Feinberg rules
  //
  

  // Thresholds

  // Minimum duration for a NREM period
  const double def_min_nrem_duration_mins = 15;

  // Minimum duration for REM period (cycle 2 and after)
  const double def_min_rem_duration_mins = 5;

  // Maximum duration of NREM/W allowed within a single REM episode
  const double def_rem_period_interuption_mins  = 15;

  // If skipping REM period, minimum W/N1 to terminate a NREM period
  const double def_terminating_waso_duration_mins = 15;

  // Persistent sleep is defined after 10 minutes of sleep
  const double def_persistent_sleep_mins = 10;


  // convert the above to epoch counts
  const int def_persistent_sleep_epochs = def_persistent_sleep_mins / (double)epoch_mins ; 
  const int def_rem_period_interuption_epochs = def_rem_period_interuption_mins / (double)epoch_mins;
  const int def_min_nrem_duration_epochs = def_min_nrem_duration_mins / (double)epoch_mins;
  const int def_terminating_waso_duration_epochs = def_terminating_waso_duration_mins / (double)epoch_mins;
  const int def_min_rem_duration_epochs = def_min_rem_duration_mins / (double)epoch_mins;


  //
  // 0)  Handle movement: 'impute' with the following epoch value
  //     ignored for now... 
  

  //  
  // 1)  Find periods of 'persistent sleep' (default 10mins prior sleep)
  //
  
  std::vector<std::string> persistent_sleep( ne , "" );
  for (int e=0;e<ne;e++)
    {
      
      if ( stages[ e ] == WAKE || stages[ e ] == LIGHTS_ON || stages[e] == UNKNOWN )
	{
	  persistent_sleep[e] = "W";
	  continue;
	}
      
      // 'not sleep' = WAKE, L or ?
      // otherwise, assume all other annotations are consistent with sleep
      bool okay = true;
      int ec = e - def_persistent_sleep_epochs;
      
      while ( okay )
	{
	  if ( ec < 0 ) { okay = false; break; }	  
	  if ( stages[ ec ] == WAKE || stages[ ec ] == LIGHTS_ON || stages[ec ] == UNKNOWN ) { okay = false; break; }  	  
	  if ( ++ec == e ) break;
	}

      if ( okay ) persistent_sleep[e] = "S"; 
      else persistent_sleep[e] = "W";
    }
  
  
  
  //
  // 2) Find sleep onset 
  //
  
  std::vector<std::string> sleep_onset( ne , "" );
  bool found_sleep = false;
  for (int e = 0; e < ne ; e++ ) 
    {
      if ( is_sleep( stages[e] ) ) found_sleep = true;
      sleep_onset[ e ] = found_sleep ? "S" : "W";
    }

  // and likewise, same from end of sleep
  found_sleep = false;
  for (int e = ne-1; e >= 0 ; e-- ) 
    {
      if ( is_sleep( stages[e] ) ) break;
      sleep_onset[ e ] = "W";
    }
 


  //
  // 3) Cumulative count of sleep
  //

  // after LightsOn, defined as '-1', if this matters 
  std::vector<int> sleep_count( ne , 0 );
  int cum_sleep = 0;
  for (int e = 0 ; e < ne ; e++)
    {
      if ( persistent_sleep[e] == "S" ) ++cum_sleep;
      if ( stages[e] == LIGHTS_ON && cum_sleep > 0 ) sleep_count[e] = -1;
      else sleep_count[e] = cum_sleep;
    }
  

  //
  // 4) Sleep state
  //
  
  std::vector<std::string> sleep_state( ne , "" );
  for (int e = 0 ; e < ne; e++)
    {
      if ( stages[e] == LIGHTS_ON && sleep_count[e] == 0 ) sleep_state[e] = "Prior";
      else if ( sleep_count[e] == 0 ) sleep_state[e] = "LPS"; // latency before persistent sleep
      else if ( sleep_count[e] == 1 ) sleep_state[e] = "LPO"; // onset of persistent sleep
      else 
	{
	  if ( sleep_count[e] > 1 ) sleep_state[e] = "SPT"; // sleep period time
	  else sleep_state[e] = "After"; 
	}
    }
  

  //
  // 5) Final wake ('WATA')
  //

  std::vector<bool> wata( ne , false );
  for (int e = ne-1 ; e>= 0 ; e-- )
    {
      // note: is_sleep() and is_wake() functions will return F is lights on
      // or undefined
      if ( is_sleep( stages[e] ) ) break;
      if ( is_wake_or_lights( stages[e] ) ) wata[e] = true;
    }
  
  //
  // 6) Sleep period/cycle 
  //
  
  std::vector<std::string> sleep_period( ne , "" );
  std::vector<bool> cycle_ending_waso( ne , false );
    
  for (int e = 0 ; e < ne ; e++ ) 
    {

      if ( sleep_onset[e] == "W" ) continue;

      bool previous_epoch_defined = e == 0 ? false : sleep_period[ e - 1 ] != "" ;

      if ( is_rem( stages[e] ) && previous_epoch_defined )
	{
	  // continues a new REM stage
	  sleep_period[e] = "REM";	  
	}
      else
	{

	  // check subsequent 15 mins
	  bool has_another_rem = false;
	  const int elimit = ne-1 < e + def_rem_period_interuption_epochs - 1 ? 
	    ne-1 : e + def_rem_period_interuption_epochs - 1 ;
	  for (int e2=e;e2 <= elimit; e2++)
	    {
	      if ( is_rem( stages[e2] ) ) { has_another_rem = true; break; }
	    }


	  // from start (i.e. including this one), next 15 mins has to have at least 1 other REM 
	  if ( ( e>0 && sleep_period[e-1] == "REM" ) && has_another_rem )
	    sleep_period[e] = "REM";
	  else
	    {

	      // else, if previously REM and cycle ended last epoch OR W/N1
	      if ( ( ( e>0 && sleep_period[e-1] == "REM" ) || ( e>0 && cycle_ending_waso[e-1] ) ) // previous REM, or cycle end
		   && ( is_wake( stages[e] ) || is_nrem1( stages[e] ) ) ) // AND W/N1 this epoch
		sleep_period[e] = ""; 
	      else
		{
		  
		  // note, uses potentially different duration from above, but obviosuly all this code 
		  // can be streamlined.
		  
		  bool has_another_rem = false;
		  const int elimit = ne-1 < e + def_min_nrem_duration_epochs - 1 ? 
		    ne-1 : e + def_min_nrem_duration_epochs - 1 ;
		  for (int e2=e;e2 <= elimit; e2++)
		    if ( is_rem( stages[e2] ) ) { has_another_rem = true; break; }
		  
		  // else, if previously in NREM, continue; (Q6!="")
		  if ( ( e==0 || sleep_period[e-1]=="" ) &&
		       ( is_wake( stages[e] ) || is_nrem1( stages[e] )  || has_another_rem ) )
		    sleep_period[e] = ""; 
		  else	// else, initiate a new NREM cycle, only if no REM wihtin the next 15mins
		    {
		      sleep_period[e] = "NREM";		      
		    }
		}
	    }
	}
      

      //
      // Cycle-ending WASO
      //

      //  put a cycle end IF currently NREM and in the next 15mins (including this epoch),
      //  if this is NREM but the next 15 mins has only W/N1, then end cycle here
      //  ELSE, if past epoch was end flag, and is WAKE, also flag (i.e. up to end of WASO W)
      
      
      bool no_near_sleep = true;
      const int elimit = ne-1 < e + def_terminating_waso_duration_epochs - 1 ? ne-1 : e+ def_terminating_waso_duration_epochs - 1 ;
      for (int e2=e;e2<=elimit;e2++)
	if ( is_nrem234( stages[e2] ) || is_rem( stages[e2] ) ) 
	  { no_near_sleep = false; break; } 
      
      if ( sleep_period[e] == "NREM" && no_near_sleep )
	{
	  cycle_ending_waso[e] = true;
	}
      else
	{
	  if ( e > 0 && cycle_ending_waso[ e - 1 ] && is_wake( stages[e] ) ) 
	    cycle_ending_waso[e] = true;
	}
      
    }

  
  //
  // Cycle type, number
  //

  sleep_code.resize( ne , 0 ); // 0, 1, 5 for W, NREM, REM
  sleep_cycle_number.resize( ne , 0 );

  
  // get first REM/cycle-ending epoch
  
  int first_sleep_period_rem = 99999;
  int first_cycle_ending_waso = 99999;
  
  for (int e=0; e<ne; e++)
    if ( sleep_period[e] == "REM" ) { first_sleep_period_rem = e; break; }
    
  for (int e=0; e<ne; e++)
    if ( cycle_ending_waso[e] ) { first_cycle_ending_waso = e; break; }

  for (int e=0; e<ne; e++)
    {

      // skip if a cycle-ending WASO
      if ( cycle_ending_waso[e] ) continue;

      if ( sleep_period[e] == "NREM" )
	{
	  sleep_code[e] = 1;
	}
      else
	{
	  
	  if ( sleep_period[e] == "REM" )
	    {
	      if ( e > 0 && sleep_period[e-1] == "NREM" )
		{
		  
		  // first cycle?
		  if ( e <= first_sleep_period_rem 
		       && e <= first_cycle_ending_waso )
		    sleep_code[e] = 5;
		  else
		    {
		      // check ahead... requires at least 'def_min_rem_duration_epochs' of REM
		      // nb. not sure why this would ever be ">=" the # of epochs... but keep as is for now
		      
		      int count_rem = 0;
		      const int elimit = ne-1 < e + def_min_rem_duration_epochs - 1 ? ne-1 : e+ def_min_rem_duration_epochs -1 ;
		      for (int e2=e; e2<=elimit; e2++)
			if ( sleep_period[e2] == "REM" ) ++count_rem;
		      if ( count_rem >= def_min_rem_duration_epochs ) 
			sleep_code[e] = 5;
		      else
			sleep_code[e] = 1;
		    }
		}
	      else
		{
		  if ( e>0 && sleep_period[e-1] == "REM" && sleep_code[e-1] == 5 ) 
		    sleep_code[e] = 5;
		  else 
		    sleep_code[e] = 1;
		}             
	    }
	  else
	    {
	      if ( e > 0 && sleep_period[e-1] == "REM" && sleep_code[e-1] == 1 )
		sleep_code[e] = 1; 
	      else
		{
		  if ( wata[e] ) 
		    sleep_code[e] = 0;
		  else
		    {
		      if ( sleep_period[e] == "" && e>0 && sleep_code[e-1] == 1 && ! cycle_ending_waso[e] )
			 sleep_code[e] = 1; 
		      else
			sleep_code[e] = 0;
		    }
		  
		}
	      
	    }
	  
	}
      
      // next epoch
    }


        
  //
  // Define cycles
  //

  
  int cnt_cycle = 0;

  for (int e=0; e< ne; e++)
    {
      if ( sleep_code[e] == 0 ) 
	sleep_cycle_number[e] = 0;
      else
	{
	  	  
	  const int previous_code = e == 0 ? 0 : sleep_code[e-1];

	  // change in cycle?
	  // start of a new NREM?
	  if ( sleep_code[e] - previous_code == 1  // into NREM(1) from WASO/N1(0) 
	       || previous_code - sleep_code[e] == 4 ) // from REM(5) to NREM(1)
	    {

	      // requires NREM (15mins) of 
	      
	      // find next REM and WASO epoch [ 'sleep_code' ]
	      // count eppochs (F) [ stages] 
	      
	      int elimit = ne - 1;
	      for (int e2 = e ; e2 < ne ; e2++)
		if ( sleep_code[e2] == 0 || sleep_code[e2] == 5 ) 
		  { elimit = e2; break; } 
	      
	      int cnt_nrem = 0;
	      for (int e2=e; e2 <= elimit; e2++)
		if ( is_nrem( stages[e2] ) ) ++cnt_nrem;
	      
	      // enough NREM for a new cycle?
	      if ( cnt_nrem >= def_min_nrem_duration_epochs ) 
		sleep_cycle_number[e] = ++cnt_cycle;
	      
	    }
	  else
	    sleep_cycle_number[e] = e == 0 ? 0 : sleep_cycle_number[e-1];
	  
	}
    }


  
  //
  // Get cycle/period statistics
  //

  // Count number of cyles
  
  num_nremc = 0;
  nremc_mean_duration = 0;

  std::map<int,int> cmin;
  std::map<int,int> cmax;
  std::map<int,double> secs_tot;
  std::map<int,double> secs_rem;
  std::map<int,double> secs_nrem;
  std::map<int,double> secs_other;
  
  for (int e=0;e<ne;e++)
    {
      const int & sn = sleep_cycle_number[e];
      if ( sn == 0 ) continue;
      
      // non-gapped epoch count
      const int eidx = epoch_n[e];
      const bool is_gap = epoch_gap[e];

      if ( ! is_gap )
	{
	  if ( sn > num_nremc ) num_nremc = sn;	  
	  // track epochs based on non-happed version	  
	  if ( cmin.find( sn ) == cmin.end() ) cmin[ sn ] = cmax[sn] = eidx;
	  cmax[sn] = eidx; // track max
	}
      
      // accumulate times (including any gaps)
      secs_tot[sn] += epoch_dur[e];
      if      ( is_rem( stages[e] ) ) secs_rem[sn] += epoch_dur[e];
      else if ( is_nrem( stages[e] ) ) secs_nrem[sn] += epoch_dur[e];
      else secs_other[sn] += epoch_dur[e]; // includes gaps
    }


  
  std::map<int,int>::iterator ii = cmin.begin();
  while ( ii != cmin.end())
    {
      const int & sn = ii->first ;
      
      // total cycle duration
      double dur = secs_tot[sn];
      double dur_mins = dur / 60.0;
      
      nremc_mean_duration += dur_mins;
      
      nremc_duration[ sn ] = ( secs_rem[sn] + secs_nrem[sn] + secs_other[sn] ) / 60.0;
      nremc_nrem_duration[ sn ] = secs_nrem[sn] / 60.0;
      nremc_rem_duration[ sn ] = secs_rem[sn] / 60.0;
      
      nremc_start_epoch[ sn ] = ii->second + 1 ;  // output 1-based coding
      
      nremc_epoch_duration[ sn ] = cmax[sn] - ii->second + 1; 
            
      ++ii;
    }

  if ( num_nremc > 0 ) nremc_mean_duration /= (double)num_nremc;


  // cycle positions for each epoch
  cycle_pos_relative.resize( ne , -1 );
  cycle_pos_absolute.resize( ne , -1 );
  
  std::map<int,double> elapsed_in_cycle;

  for (int e=0; e<ne; e++)
    {
      const int & sn = sleep_cycle_number[e];

      if ( sn == 0 ) continue;      
      
      int cycle_start = cmin[sn];
      
      // position within each cycle.xxxxxx
      //cycle_pos_absolute[e] = ( e - cycle_start ) * epoch_mins ;
      // cycle_pos_absolute[e] = ( e - cycle_start ) * epoch_dur[e] mins ; 
      // cycle_pos_relative[e] = cycle_pos_absolute[e] / (double)nremc_duration[sn];

      cycle_pos_absolute[e] = elapsed_in_cycle[sn] / 60.0 ; 
      cycle_pos_relative[e] = cycle_pos_absolute[e] / (double)nremc_duration[sn];

      // increment
      elapsed_in_cycle[sn] += epoch_dur[e];
      
    }


  // after the fact, track epoch-level stats
  in_persistent_sleep.resize( ne , false );
  for (int e=0; e<ne; e++)
    in_persistent_sleep[e] = persistent_sleep[e] == "S" ;
  
  // do not alter original persistent sleep definition, as that was used in NREM cycle
  // construct;  but here we need to add in the fact that we've skipped the e.g. 10 mins
  // of sleep prior to start of persistent sleep... this is included in PER_SLP_LAT but
  // original not in TPST.  so add in here...

  for (int e=1; e<ne; e++)
    {
      
      // std::cout << " e = " << e << "\t"
      //  		<< globals::stage( stages[e] ) << "\t"
      //  		<< in_persistent_sleep[e] << "\n";
      
      // start of a persistent sleep bout?
      if ( in_persistent_sleep[e] && ! in_persistent_sleep[e-1] )
	{
	  int ec = e;
	  for (int i=0; i<def_persistent_sleep_epochs; i++)
	    {
	      --ec;
	      if ( ec < 0 || in_persistent_sleep[ec] )
		{
		  logger << "  first epoch of persistent sleep e = " << e << "\n";
		  logger << "  tracking back " << def_persistent_sleep_epochs << " epochs of sleep, to also mark as persistent\n";
		  if ( ec < 0 ) 
		    logger << "  however, epoch count is less than 0 (ec = " << ec << ")\n";
		  else
		    logger << "  however, encountering epochs already marked as persistent.. this should not happen\n";
		  
		  Helper::halt( "error defining persistent sleep bouts... check stage/epoch alignment (EPOCH align)" );
		}
	      in_persistent_sleep[ec] = true;
	    }
	}
    }
  
  TpST = 0;
  for (int e=0; e<ne; e++)
    if ( in_persistent_sleep[e] ) TpST += epoch_mins;
	
  
  //
  // Sleep cycle definitions
  //
  
  // Based primarily on Feinberg & Floyd (1979) definitions
  // with some modifications
  
  // Cycle = "NREM phase of at least X mins terminated by end of REM
  // phase OR a preset duration of a bout of wake/N1
  
  // "Wake bout" = periods of wake of any duration after latency to
  // persistent sleep onset
  
  // If REM phase is skipped (i.e. cycle is terminated by a wake bout
  // of a modifiable duration in minutes), bout of wake/N1 is excluded
  // from cycle new NREM phase onset is with first stage of sleep;
  // wake between cycles is excluded from a cycle and is denoted cycle
  // 0, but wake within a cycle is still included


  // REM phase is period of REM sleep of any duration terminated by a
  // preset contiguous duration of any other stage (wake or NREM).
  // REM phase may contain NREM or wake as long as it is less than the
  // preset duration of any combination of NREM or wake.  the user may
  // define whether stage N1 can define the onset of a new NREM phase.
 
  
  //
  // Flanking epochs 
  //

  is_waso.resize( ne , false );
  for (int e = 0 ; e < ne ; e++)
    {
      if ( stages[e] == WAKE 
	   && e > first_sleep_epoch 
	   && e < final_wake_epoch ) is_waso[e] = true;
    }

  
  flanking.resize( ne , 0 );
  flanking_tot.resize( ne , 0 );
  nearest_wake.resize( ne , 0 );

  // transitions in/out of NREM, REM or W sleep
  //  (3class analysis)
  // note: if not flanking-collapse-nrem=T (default) then 
  // this will collapse NR when finding transitions, but 
  //  if any req-pre-post=X is supplied, these values will
  // be based on a 5-class (or 6-class) NREM value... a 
  // little inconsistent perhaps, but the default will be good

  nrem2rem.resize( ne , 0 ); nrem2rem_total.resize( ne , 0 );
  nrem2wake.resize( ne , 0 ); nrem2wake_total.resize( ne , 0 );

  rem2nrem.resize( ne , 0 ); rem2nrem_total.resize( ne , 0 );
  rem2wake.resize( ne , 0 ); rem2wake_total.resize( ne , 0 );
  
  wake2nrem.resize( ne , 0 ); wake2nrem_total.resize( ne , 0 );
  wake2rem.resize( ne , 0 ); wake2rem_total.resize( ne , 0 );

  transitions.clear();
  transitions5.clear();

  for (int e = 0 ; e < ne ; e++)
    {
      
      //
      // calculate the number of similar epochs 
      // (FLANKING and FLANKING_ALL)
      //
      
      int sim = 0;  
      
      for (int j=1;j<ne;j++)
	{
	  const int eleft  = e - j;
	  const int eright = e + j;
	  // too much
	  if ( eleft < 0 || eright >= ne ) { sim = j-1; break; }
	  
	  // 3 class or full comparison?
	  if ( flanking_3class )
	    {
	      if (  ( ! is_same_3class( stages[eleft] , stages[e] ) )
		    || ( ! is_same_3class( stages[eright] , stages[e] ) ) )
		{ sim = j-1; break; }	  
	    }
	  else
	    if ( stages[eleft] != stages[e] || stages[eright] != stages[e] ) 
	      { sim = j-1; break; }	  
	}
      
      
      int sim_all = 1; // self
      
      // forward
      for (int ee=e+1 ; ee<ne; ee++)
	{
	  if ( flanking_3class )
	    if ( is_same_3class( stages[ee] , stages[e] ) ) ++sim_all; else break;
	  else
	    if ( stages[ee] == stages[e] ) ++sim_all; else break;
	}

      // backward      
      for (int ee=e-1 ; ee != -1; ee--)
	{
	  if ( flanking_3class )
	    if ( is_same_3class( stages[ee] , stages[e] ) )  ++sim_all; else break;    
	  else
	    if ( stages[ee] == stages[e] ) ++sim_all; else break;    
	}

      int nw = 0;

      if ( stages[e] != WAKE )
	{
	  for (int j=1;j<ne;j++)
	    {
	      const int eleft  = e - j;
	      const int eright = e + j;
	      // too much
	      if ( eleft < 0 || eright >= ne ) { nw = j; break; }
	      if ( stages[eleft] == WAKE || stages[eright] == WAKE ) { nw = j; break; }	  
	    }
	}

      flanking[e] = sim;
      flanking_tot[e] = sim_all;
      nearest_wake[e] = nw;

      //
      // Generic transition matrix counts
      //
            
      if ( e != 0 )	
	{
	  if ( flanking_3class )
	    {
	      // make all NREM --> NREM2 (and switch output) 
	      sleep_stage_t ss1 = is_nrem( stages[ e - 1 ] ) ? NREM2 : stages[ e - 1 ] ;
	      sleep_stage_t ss2 = is_nrem( stages[ e ] ) ? NREM2 : stages[ e ] ;
	      ++transitions[ ss1 ][ ss2 ];
	    }
	  else	    
	    ++transitions[ stages[ e - 1 ] ][ stages[e] ];

	  // return 5-class transitions in either case (for STI etc)
	  ++transitions5[ stages[ e - 1 ] ][ stages[e] ];
	  
	}


    } // next epoch
      

  //
  // Loop again over epochs (as we need to calc the flanking_tot[] value for the /next/
  // epoch here
  //
  
  for (int e=0; e<ne; e++)
    {
      //
      // transitions FROM NREM? 
      //
      
      if ( is_nrem( stages[e] ) ) 
	{
	  
	  // nr to rem
	  int ei = 1;
	  while ( 1 ) 
	    {
	      if ( e + ei == ne ) { ei=0; break; }
	      if ( is_nrem( stages[ e + ei ] ) ) { ++ei; continue; }
	      if ( is_rem( stages[ e + ei ] ) && flanking_tot[e+ei] >= req_pre_post_epochs ) break;
	      ei = 0; break;
	    }
	  nrem2rem[e] = ei;
	  
	  // nr to wake
	  ei = 1;
	  while ( 1 )
	    {
	      if ( e + ei == ne ) { ei=0; break; }
	      if ( is_nrem( stages[ e + ei ] ) ) { ++ei; continue; }
	      if ( is_wake( stages[ e + ei ] ) && flanking_tot[e+ei] >= req_pre_post_epochs ) break;
	      ei = 0; break;
	    }
	  nrem2wake[e] = ei;
	}


      //
      // transitions FROM REM?
      //
      
      if ( is_rem( stages[e] ) ) 
	{
	  
	  // rem to nr
	  int ei = 1;
	  while ( 1 ) 
	    {
	      if ( e + ei == ne ) { ei=0; break; }
	      if ( is_rem( stages[ e + ei ] ) ) { ++ei; continue; }
	      if ( is_nrem( stages[ e + ei ] ) && flanking_tot[e+ei] >= req_pre_post_epochs ) break;
	      ei = 0; break;
	    }
	  rem2nrem[e] = ei;
	  
	  // rem to wake
	  ei = 1;
	  while ( 1 )
	    {
	      if ( e + ei == ne ) { ei=0; break; }
	      if ( is_rem( stages[ e + ei ] ) ) { ++ei; continue; }
	      if ( is_wake( stages[ e + ei ] ) && flanking_tot[e+ei] >= req_pre_post_epochs ) break;
	      ei = 0; break;
	    }
	  rem2wake[e] = ei;
	}

      //
      // transitions FROM wake?
      //
      
      if ( is_wake( stages[e] ) ) 
	{
	  
	  // w to nr
	  int ei = 1;
	  while ( 1 ) 
	    {
	      if ( e + ei == ne ) { ei=0; break; }
	      if ( is_wake( stages[ e + ei ] ) ) { ++ei; continue; }
	      if ( is_nrem( stages[ e + ei ] ) && flanking_tot[e+ei] >= req_pre_post_epochs ) break;
	      ei = 0; break;
	    }
	  wake2nrem[e] = ei;
	  
	  // w to rem
	  ei = 1;
	  while ( 1 )
	    {
	      if ( e + ei == ne ) { ei=0; break; }
	      if ( is_wake( stages[ e + ei ] ) ) { ++ei; continue; }
	      if ( is_rem( stages[ e + ei ] ) && flanking_tot[e+ei] >= req_pre_post_epochs ) break;
	      ei = 0; break;
	    }
	  wake2rem[e] = ei;
	}
           
    
    } // next epoch
  
  // now figure out the _total values 
  // i.e. move forward and copy largest number until we hit 0      

  int e_nrem2rem  = nrem2rem[0];
  int e_nrem2wake = nrem2wake[0]; 

  int e_rem2nrem  = rem2nrem[0];
  int e_rem2wake  = rem2wake[0]; 
  
  int e_wake2nrem = wake2nrem[0];
  int e_wake2rem  = wake2rem[0]; 

  for (int e = 1 ; e < ne ; e++ )
    {
      // NR --> 
      if ( nrem2rem[e] == 0 ) e_nrem2rem = 0;
      else if ( nrem2rem[e] > e_nrem2rem ) e_nrem2rem = nrem2rem[e];
      nrem2rem_total[e] = e_nrem2rem;
      
      if ( nrem2wake[e] == 0 ) e_nrem2wake = 0;
      else if ( nrem2wake[e] > e_nrem2wake ) e_nrem2wake = nrem2wake[e];
      nrem2wake_total[e] = e_nrem2wake;

      // REM -> 
      if ( rem2nrem[e] == 0 ) e_rem2nrem = 0;
      else if ( rem2nrem[e] > e_rem2nrem ) e_rem2nrem = rem2nrem[e];
      rem2nrem_total[e] = e_rem2nrem;

      if ( rem2wake[e] == 0 ) e_rem2wake = 0;
      else if ( rem2wake[e] > e_rem2wake ) e_rem2wake = rem2wake[e];
      rem2wake_total[e] = e_rem2wake;

      // wake ->
      if ( wake2nrem[e] == 0 ) e_wake2nrem = 0;
      else if ( wake2nrem[e] > e_wake2nrem ) e_wake2nrem = wake2nrem[e];
      wake2nrem_total[e] = e_wake2nrem;

      if ( wake2rem[e] == 0 ) e_wake2rem = 0;
      else if ( wake2rem[e] > e_wake2rem ) e_wake2rem = wake2rem[e];
      wake2rem_total[e] = e_wake2rem;
    }

  
  int first_lights_out_epoch = 0;
  for (int e=0; e<ne; e++)
    {
      if ( stages[e] != LIGHTS_ON ) 
	{
	  first_lights_out_epoch = e;
	  break;
	}
    }

  // +1 is okay even if lights on is at end	
  // as this will never be used other than to give a
  // time below (i.e. start of epoch 'ne' = last point of recording) +1LLU
  int first_lights_on_epoch = ne;
  for (int e=ne-1; e>=0; e--)
    {
      if ( stages[e] != LIGHTS_ON )
	{
	  first_lights_on_epoch = e+1;
	  break;
	}
    }


  //
  // Clocktime-based measures
  //
  
  clocktime_t starttime( timeline->edf->header.starttime );
  if ( ! starttime.valid ) 
    {
      clock_start.valid = clock_lights_out.valid = clock_sleep_onset.valid 
	= clock_sleep_midpoint.valid = clock_wake_time.valid 
	= clock_lights_on.valid = clock_stop.valid = false;
    }
  else
    {
      // T0
      clock_start = starttime;
      
      double epoch_hrs = epoch_mins / 60.0;      

      // T1
      clock_lights_out = starttime;
      clock_lights_out.advance_hrs( epoch_hrs * first_lights_out_epoch );

      // T2
      clock_sleep_onset = starttime;
      clock_sleep_onset.advance_hrs( epoch_hrs * first_sleep_epoch );
      
      // T4
      clock_wake_time = starttime;
      clock_wake_time.advance_hrs( epoch_hrs * final_wake_epoch );
      
      // T5 
      clock_lights_on = starttime;
      clock_lights_on.advance_hrs( epoch_hrs * first_lights_on_epoch );
    
      // T6 end
      clock_stop = starttime;
      clock_stop.advance_hrs( epoch_hrs * ne );

      // T3 (midpoint)
      clock_sleep_midpoint.midpoint( clock_sleep_onset , clock_wake_time );      
      
    }



}


void hypnogram_t::annotate( const std::string & annot_prefix , const std::string & suffix )
{
  
  logger << "  creating hypnogram-derived annotations, with prefix " << annot_prefix << "\n";

  // either nothing, or insert underscore if prefix specified
  const std::string prefix = annot_prefix != "" ? annot_prefix + "_" : "" ;

  // use annots in class:value format, so that one can split into class/instance IDs downstream as needed
  // if make suffix == ':' then annotations will be split into class : instance IDs when loaded downstream
  // by default, suffix is '_'
 
  // all annotations start <prefix>_
  
  // NREMC cycle
  //    number: hyp_cycle_1 hyp_cycle_2 ... 
  //    mins (abs)   : hyp_cycle_m0_10 hyp_cycle_m10_20 ...
  //    pct position (quinites) : hyp_cycle_q1 hyp_cycle_q2 ... hyp_cycle_q5
  
  // Bouts
  //    if in a 5-min / 10-min bout
  //    for stages N1, N2, N3, R, W, NR, S, WASO
  //    annotaton labels:   hyp_bout5_N1   hyp_bout10_S  etc

  // 0-duration mark-points
  
  // Transitions 
  //    hyp_tr_W_NR   hyp_tr_W_NR  etc

  // Time-points 
  //   hyp_t0_start, hyp_t1_lights_off, hyp_t2_sleep_onset,
  //   hyp_t3_sleep_midpoint, hyp_t4_final_wake, hyp_t5_lights_on, hyp_t6_stop

  // Epoch-level annotations
  //  WASO   waso
  //  Pre-sleep wake     pre_sleep_wake
  //  Post-sleep wake    post_sleep_wake

  //  Pre-sleep          pre_sleep (any epoch whether wake or not, before SO)
  //  Post-sleep         post_sleep
  
  // Cumulative elapsed durations
  //   Clock-time (by 24-hour) 
  //     hyp_clock:20 hyp_clock:23 hyp_clock:00 etc
  //   Elapsed time (hours) from sleep onset (1-baed counting)
  //     hyp_elapsed:1hr hyp_elapsed:2hr ...
  //   Elapsed stage : hours or quintiles
  //     hyp_n1:1hr hyp_n1:2hr ...
  //     hyp_n1:q1  hyp_n1:q2 .... 
  //     for : E_N1 E_N2 E_N3 E_REM E_SLEEP E_WAKE E_WASO

  // Misc

  // Ascending/descending N2 ( Y / N /not sure)
  //   hyp_N2_asc hyp_N2_dsc

  // basic 'lights on' (i.e. will reflect any truncation too)
  
  //
  // NREM cycles
  //

  std::map<int,double>::iterator cc = nremc_duration.begin();
  while ( cc != nremc_duration.end() )
    {      

      // cc->first  : 1-based cycle #
      // nremc_start_epoch[ cc->first ] : 1-based epoch start
      // nremc_epoch_duration[ cc->first ] : cycle length ( in epochs )
      
      // NREMC number
      std::string cn = Helper::int2str( cc->first );

      // annot class
      annot_t * a = timeline->annotations->add( prefix + "cycle" + suffix + "n" + cn );
      
      // epoch number start (adjust for 1-base encoding)
      int start_epoch = nremc_start_epoch[ cc->first ] - 1;
      
      // length of cycle (in epochs) , minus 1 as we'll add this to stop of first
      int length = nremc_epoch_duration[ cc->first ] - 1;
      
      // get interval
      interval_t interval = timeline->epoch( start_epoch );
      
      // adjust end-point -> convert to time-points
      interval.stop += (uint64_t)(timeline->epoch_length_tp * length ) ;
      
      // add annotation
      instance_t * instance = a->add( "." , interval , "." );

      // quintiles of position
      uint64_t len = interval.duration(); 
      uint64_t q5  = len / 5LLU ;
      for (int q=0; q<5; q++)
	{
	  annot_t * a = timeline->annotations->add( prefix + "cycle" + suffix + "q" + Helper::int2str( q+1 ) );
	  interval_t qinterval( interval.start + q5 * q , interval.start + q5 * (q+1) );
	  instance_t * instance = a->add( "." , qinterval , "." );
	}

      // mins (10-min blocks) 
      uint64_t min10 = 10 * 60 * globals::tp_1sec ; 
      int mc = 0;
      for (uint64_t m=0; m < len; m += min10 )
	{
	  if ( m + min10 <= len )
	    {	      
	      annot_t * a = timeline->annotations->add( prefix + "cycle" + suffix + "m" + Helper::int2str( mc ) + "_" + Helper::int2str( mc+10)  );
	      interval_t minterval( interval.start + m , interval.start + m + min10 );
	      instance_t * instance = a->add( "." , minterval , "." );
	      mc += 10;
	    }
	}
            
      ++cc;
    }

  //
  // Bouts
  //
  //    annotaton labels:   hyp_bout5:N1   hyp_bout10:S  etc

  std::set<bout_t>::const_iterator bb = bouts.begin();
  while ( bb != bouts.end() )
    {
      const int e1 = epoch_n[ bb->start ] ;
      const int e2 = epoch_n[ bb->stop ] ;

      interval_t interval1 = timeline->epoch( e1 );
      interval_t interval2 = timeline->epoch( e2 );
      interval_t interval( interval1.start , interval2.stop );

      // nb. NREM downcasting
      std::string stg;
      if ( bb->ss == NREM2 ) 
	stg = "NR";
      else
	stg = globals::stage( bb->ss );

      
      double len = ( ( e2 - e1 + 1 ) * timeline->epoch_length() ) / 60.0 ;

      if ( len >= 10 )
	{
	  annot_t * a = timeline->annotations->add( prefix + "bout10" + suffix + stg ) ;	  
	  a->add( "." , interval , "." );
	}
      else if ( len >= 5 )
	{
	  annot_t * a = timeline->annotations->add( prefix + "bout05" + suffix + stg ) ;	  
	  a->add( "." , interval , "." );	  
	}
      
      ++bb;
    }


  //
  // Clock-time & elapsed/cumulative stage times
  //
  
  // 24-hr clocktime: 
  //   hyp_clock:20 hyp_clock:23 hyp_clock:00 etc                                                                                                                
  
  // get initial 
  int hr = clock_start.h;

  // seconds until next hour (will always be integer s)
  double secs2next = 60 * ( 60 - clock_start.m ) - clock_start.s;
  uint64_t tp = secs2next * globals::tp_1sec ;  

  // fixed 1 hour
  const uint64_t hour_tp = 3600 * globals::tp_1sec ;
      
  // first interval (start to hour), may be fractional, as always starts at zero 
  std::string hrstr = ( hr < 10 ? "0"	: "" ) + Helper::int2str( hr );
  annot_t * ah1 = timeline->annotations->add( prefix + "clock" + suffix + hrstr );
  ah1->add( "." , interval_t( 0 , tp ) , "." );

  // move to next hour
  ++hr;
  
  // continue w/ full hours
  while ( tp < timeline->last_time_point_tp )
    {      
      std::string hrstr = ( hr < 10 ? "0" : "" ) + Helper::int2str( hr );      
      annot_t * ah1 = timeline->annotations->add( prefix + "clock" + suffix + hrstr );
      ah1->add( "." , interval_t( tp , tp + hour_tp ) , "." );
      
      // move to the next hour
      ++hr;
      if ( hr == 24 ) hr = 0;
      tp += hour_tp;
    }
  
  
  //   Elapsed time (hours) from sleep onset (1-based counting)
  //     hyp_elapsed_h1 hyp_elapsed_h2 ...                                                                                                                       
  // reset hr to first (1-based counting)

  hr = 1;
  tp = 0; // start of record
  
  while ( tp < timeline->last_time_point_tp )
    {
      std::string hrstr = "h" + Helper::int2str( hr ) ;
      annot_t * ah1 = timeline->annotations->add( prefix + "elapsed" + suffix + "T_" + hrstr );
      ah1->add( "." , interval_t( tp , tp + hour_tp ) , "." );
      // move to the next hour                                                                                                                                    
      ++hr;
      tp += hour_tp;
    }


  //   Elapsed stage : hours or quintiles
  //     hyp_n1_h1  hyp_n1_h2 ...
  //     hyp_n1_q1  hyp_n1_q2 ... hyp_n1_q5
  //     for : N1 N2 N3 NR R S W WASO

  // add annotation at epoch levels:
  
  for (int e=0; e<ne; e++)
    {
      
      // or epoch_n[] ?      
      interval_t interval = timeline->epoch( e );
      
      bool is_wake = stages[e] == WAKE ;
      bool is_waso = stages[e] == WAKE && e > first_sleep_epoch && e < final_wake_epoch ;
      bool is_n1   = stages[e] == NREM1;
      bool is_n2   = stages[e] == NREM2;
      bool is_n3   = stages[e] == NREM3 || stages[e] == NREM4;
      bool is_nr   = is_n1 || is_n2 || is_n3 ;
      bool is_rem  = stages[e] == REM ;
      bool is_sleep = is_nr || is_rem ;
      bool is_lights = stages[e] == LIGHTS_ON ;

      if ( is_lights )
	{
	  std::string alabel = prefix + "lights";
          annot_t * ah1 = timeline->annotations->add( alabel );
          ah1->add( "." , interval , "." );
	}
      
      if ( is_wake )
	{
	  // abs
	  double a = elapsed_stg_sec[ "W" ][ e ];
	  int h = 1 + floor( a / 60.0 ) ;
	  std::string alabel = prefix + "elapsed" + suffix + "W_h" + Helper::int2str( h ) ;
	  annot_t * ah1 = timeline->annotations->add( alabel );
	  ah1->add( "." , interval , "." );

	  // rel
	  double r = elapsed_stg_rel[ "W" ][ e ];
	  int q = r == 1 ? 5 : floor( r / 0.2 ) + 1 ;
          alabel = prefix + "elapsed" + suffix + "W_q" + Helper::int2str( q ) ;
          annot_t * ah2 = timeline->annotations->add( alabel );
          ah2->add( "." , interval , "." );	 	  
	}

      if ( is_waso )
	{
	  // abs
	  double a = elapsed_stg_sec[ "WASO" ][ e ];
	  int h = 1 + floor( a / 60.0 ) ;
	  std::string alabel = prefix + "elapsed" + suffix + "waso_h" + Helper::int2str( h ) ;
	  annot_t * ah1 = timeline->annotations->add( alabel );
	  ah1->add( "." , interval , "." );

	  // rel
	  double r = elapsed_stg_rel[ "WASO" ][ e ];
	  int q = r == 1 ? 5 : floor( r / 0.2 ) + 1 ;
          alabel = prefix + "elapsed" + suffix + "waso_q" + Helper::int2str( q ) ;
          annot_t * ah2 = timeline->annotations->add( alabel );
          ah2->add( "." , interval , "." );	 	  
	}

      if ( is_sleep )
	{
	  // abs
	  double a = elapsed_stg_sec[ "S" ][ e ];
	  int h = 1 + floor( a / 60.0 ) ;
	  std::string alabel = prefix + "elapsed" + suffix + "S_h" + Helper::int2str( h ) ;
	  annot_t * ah1 = timeline->annotations->add( alabel );
	  ah1->add( "." , interval , "." );

	  // rel
	  double r = elapsed_stg_rel[ "S" ][ e ];
	  int q = r == 1 ? 5 : floor( r / 0.2 ) + 1 ;
          alabel = prefix + "elapsed" + suffix + "S_q" + Helper::int2str( q ) ;
          annot_t * ah2 = timeline->annotations->add( alabel );
          ah2->add( "." , interval , "." );	 	  
	}

      if ( is_n1 )
	{
	  // abs
	  double a = elapsed_stg_sec[ "N1" ][ e ];
	  int h = 1 + floor( a / 60.0 ) ;
	  std::string alabel = prefix + "elapsed" + suffix + "N1_h" + Helper::int2str( h ) ;
	  annot_t * ah1 = timeline->annotations->add( alabel );
	  ah1->add( "." , interval , "." );

	  // rel
	  double r = elapsed_stg_rel[ "N1" ][ e ];
	  int q = r == 1 ? 5 : floor( r / 0.2 ) + 1 ;
          alabel = prefix + "elapsed" + suffix + "N1_q" + Helper::int2str( q ) ;
          annot_t * ah2 = timeline->annotations->add( alabel );
          ah2->add( "." , interval , "." );	 	  
	}

      if ( is_n2 )
	{
	  // abs
	  double a = elapsed_stg_sec[ "N2" ][ e ];
	  int h = 1 + floor( a / 60.0 ) ;
	  std::string alabel = prefix + "elapsed" + suffix + "N2_h" + Helper::int2str( h ) ;
	  annot_t * ah1 = timeline->annotations->add( alabel );
	  ah1->add( "." , interval , "." );

	  // rel
	  double r = elapsed_stg_rel[ "N2" ][ e ];
	  int q = r == 1 ? 5 : floor( r / 0.2 ) + 1 ;
          alabel = prefix + "elapsed" + suffix + "N2_q" + Helper::int2str( q ) ;
          annot_t * ah2 = timeline->annotations->add( alabel );
          ah2->add( "." , interval , "." );	 	  
	}

      if ( is_n3 )
	{
	  // abs
	  double a = elapsed_stg_sec[ "N3" ][ e ];
	  int h = 1 + floor( a / 60.0 ) ;
	  std::string alabel = prefix + "elapsed" + suffix + "N3_h" + Helper::int2str( h ) ;
	  annot_t * ah1 = timeline->annotations->add( alabel );
	  ah1->add( "." , interval , "." );

	  // rel
	  double r = elapsed_stg_rel[ "N3" ][ e ];
	  int q = r == 1 ? 5 : floor( r / 0.2 ) + 1 ;
          alabel = prefix + "elapsed" + suffix + "N3_q" + Helper::int2str( q ) ;
          annot_t * ah2 = timeline->annotations->add( alabel );
          ah2->add( "." , interval , "." );	 	  
	}
      
      if ( is_nr )
	{
	  // abs
	  double a = elapsed_stg_sec[ "NR" ][ e ];
	  int h = 1 + floor( a / 60.0 ) ;
	  std::string alabel = prefix + "elapsed" + suffix + "NR_h" + Helper::int2str( h ) ;
	  annot_t * ah1 = timeline->annotations->add( alabel );
	  ah1->add( "." , interval , "." );

	  // rel
	  double r = elapsed_stg_rel[ "NR" ][ e ];
	  int q = r == 1 ? 5 : floor( r / 0.2 ) + 1 ;
          alabel = prefix + "elapsed" + suffix + "NR_q" + Helper::int2str( q ) ;
          annot_t * ah2 = timeline->annotations->add( alabel );
          ah2->add( "." , interval , "." );	 	  
	}

      if ( is_rem )
	{
	  // abs
	  double a = elapsed_stg_sec[ "R" ][ e ];
	  int h = 1 + floor( a / 60.0 ) ;
	  std::string alabel = prefix + "elapsed" + suffix + "R_h" + Helper::int2str( h ) ;
	  annot_t * ah1 = timeline->annotations->add( alabel );
	  ah1->add( "." , interval , "." );

	  // rel
	  double r = elapsed_stg_rel[ "R" ][ e ];
	  int q = r == 1 ? 5 : floor( r / 0.2 ) + 1 ;
          alabel = prefix + "elapsed" + suffix + "R_q" + Helper::int2str( q ) ;
          annot_t * ah2 = timeline->annotations->add( alabel );
          ah2->add( "." , interval , "." );	 	  
	}


      
      
      // next epoch
    }
  

  
  //
  // Time points (0-dur events)
  //
  //   hyp_t0_start, hyp_t1_lights_off, hyp_t2_sleep_onset,
  //   hyp_t3_sleep_midpoint, hyp_t4_final_wake, hyp_t5_lights_on, hyp_t6_stop

  annot_t * a0 = timeline->annotations->add( prefix + "t0_start" );
  a0->add( "." , interval_t( tp_0_start , tp_0_start ) , "." );

  // std::cout << " tp_0_start = " << tp_0_start << "\n"
  // 	    << " t1_lights_out = " << tp_1_lights_out << "\n"
  // 	    << " t2_sleep_onset = " << tp_2_sleep_onset << "\n"
  // 	    << " t3_sleep_midpoint = " << tp_3_sleep_midpoint << "\n"
  // 	    << " t4_final_wake = " << tp_4_final_wake << "\n"
  // 	    << " t5_lights_on = " << tp_5_lights_on << "\n"
  // 	    << " t6_stop = " << tp_6_stop << "\n";
  
  annot_t * a1 = timeline->annotations->add( prefix + "t1_lights_out" );
  a1->add( "." , interval_t( tp_1_lights_out , tp_1_lights_out ) , "." );
  
  
  //   hyp_t3_sleep_midpoint, hyp_t4_final_wake, hyp_t5_lights_on, hyp_t6_stop                                                                 
  
  if ( tp_2_sleep_onset != tp_4_final_wake ) 
    {

      annot_t * a2 = timeline->annotations->add( prefix + "t2_sleep_onset" );
      a2->add( "." , interval_t( tp_2_sleep_onset , tp_2_sleep_onset ) , "." );

      annot_t * a3 = timeline->annotations->add( prefix + "t3_sleep_midpoint" );
      a3->add( "." , interval_t( tp_3_sleep_midpoint , tp_3_sleep_midpoint ) , "." );

      annot_t * a4 = timeline->annotations->add( prefix + "t4_final_wake" );
      a4->add( "." , interval_t( tp_4_final_wake , tp_4_final_wake ) , "." );
    }
  
  annot_t * a5 = timeline->annotations->add( prefix + "t5_lights_on" );
  a5->add( "." , interval_t( tp_5_lights_on , tp_5_lights_on ) , "." );

  annot_t * a6 = timeline->annotations->add( prefix + "t6_stop" );
  a6->add( "." , interval_t( tp_6_stop , tp_6_stop ) , "." );
  
  
  //
  // Epoch level annotations:
  //
  
  //
  // Transitions: hyp_tr:W_NR   hyp_tr:W_NR  etc
  // Ascending/descending N2: hyp_N2:asc hyp_N2:dsc
  // WASO : waso, pre_sleep_wake, post_sleep_wake, pre_sleep, post_sleep
  //

  // get final wake
  int final_sleep_epoch = ne-1;

  for (int e=0; e<ne; e++)
    {

      interval_t interval = timeline->epoch( e );
      
      const bool is_wake = stages[e] == WAKE ;

      const bool is_pre_sleep = e < first_sleep_epoch ;
      const bool is_post_sleep = e >= final_wake_epoch ;

      const bool is_pre_sleep_wake = is_wake && is_pre_sleep ; 
      const bool is_post_sleep_wake = is_wake && is_post_sleep ;      
      
      // transition var of '1' means at *end* of that
      // epoch there is a transition

      const bool is_nr2r = nrem2rem[e] == 1;
      const bool is_nr2w = nrem2wake[e] == 1;
      const bool is_r2nr = rem2nrem[e] == 1;
      const bool is_r2w  = rem2wake[e] == 1;
      const bool is_w2nr = wake2nrem[e] == 1;
      const bool is_w2r  = wake2rem[e] == 1;
      
      if ( is_nr2r || is_nr2w || is_r2nr || is_r2w || is_w2nr || is_w2r )
	{
	  std::string tr = "NR";
	  if ( is_r2nr || is_r2w )
	    tr = "R";
	  else if ( is_w2nr || is_w2r )
	    tr = "W";
	  
	  if ( is_r2nr || is_w2nr )
	    tr += "_NR";
	  else if ( is_nr2r || is_w2r )
	    tr += "_R";
	  else
	    tr += "_W";
	  
	  annot_t * a = timeline->annotations->add( prefix + "tr" + suffix + tr );
	  interval_t tinterval( interval.stop , interval.stop );
	  instance_t * instance = a->add( "." , tinterval , "." );
	  
	}

      // WASO?
      if ( is_waso[e] )
	{
	  annot_t * a = timeline->annotations->add( prefix + "waso" );
          instance_t * instance = a->add( "." , interval , "." );
	}

      // pre-sleep wake?
      if ( is_pre_sleep_wake )
	{
	  annot_t * a = timeline->annotations->add( prefix + "pre_sleep_wake" );
          instance_t * instance = a->add( "." , interval , "." );
	}
      
      // post-sleep wake?
      if ( is_post_sleep_wake )
	{
	  annot_t * a = timeline->annotations->add( prefix + "post_sleep_wake" );
          instance_t * instance = a->add( "." , interval , "." );
	}

      // pre-sleep?
      if ( is_pre_sleep )
	{
	  annot_t * a = timeline->annotations->add( prefix + "pre_sleep" );
          instance_t * instance = a->add( "." , interval , "." );
	}
      
      // post-sleep?
      if ( is_post_sleep )
	{
	  annot_t * a = timeline->annotations->add( prefix + "post_sleep" );
          instance_t * instance = a->add( "." , interval , "." );
	}

      // Ascending/descending N2
      if ( n2_ascdesc[e] >= 0.25 )
	{
	  annot_t * a = timeline->annotations->add( prefix + "N2" + suffix + "asc" );
	  instance_t * instance = a->add( "." , interval , "." );
	}
      
      if ( n2_ascdesc[e] <= -0.25 )
	{
	  annot_t * a = timeline->annotations->add( prefix + "N2" + suffix + "dsc" );
	  instance_t * instance = a->add( "." , interval , "." );
	}

      // persistent sleep
      if ( in_persistent_sleep[e] )
	{
	  annot_t * a = timeline->annotations->add( prefix + "persistent_sleep" );
          instance_t * instance = a->add( "." , interval , "." );
	}
      
    }

}


void hypnogram_t::output( const bool verbose ,
			  const bool epoch_lvl_output ,
			  const bool extra_output ,
			  const std::string & eannot ,
			  const std::string & annot_prefix ,
			  const std::string & annot_suffix )
{
  
  
  //
  // Minimal output (stages .eannot style) to std::cout only?
  //

  const bool minimal = eannot == ".";



  //
  // Epoch-level annotation of NREM cycles also (both STAGE and HYPNO)
  //   (really, special legacy case, used by PSD dynamics...
  //    can likely retire over time)
  
  for (int e=0;e<timeline->num_epochs() ;e++)
    if ( sleep_cycle_number[e] ) 
      {	
	const std::string cycle = "_NREMC_" + Helper::int2str( sleep_cycle_number[e] );	
	timeline->annotate_epoch( cycle , e );	
      }

  
  
  //
  // currently, this routine is hard-coded to assume 30-second epochs,
  // so for now flag if this is not the case (we can fix downstream)
  //
  
  if ( verbose )
    {
      if ( ! Helper::similar( timeline->epoch_length() , 30 , 0.001 ) ) 
	Helper::halt( "requires 30-second epochs to be set currently" );
      if ( ! Helper::similar( timeline->epoch_inc() , 30 , 0.001 ) ) 
	Helper::halt( "requires non-overlapping 30-second epochs to be set currently" );
    }

  // also saved as tp (for annot output)
  tp_0_start = 0LLU;
  tp_1_lights_out = 0LLU;
  tp_2_sleep_onset = 0LLU;
  tp_3_sleep_midpoint = 0LLU;
  tp_4_final_wake = 0LLU;
  tp_5_lights_on = 0LLU;
  tp_6_stop = 0LLU;

  
  //
  // Per individual level output (VERBOSE MODE ONLY)
  //
  
  if ( verbose )
    {
      
      // HMS :: hh:mm:ss clocktime
      // E   :: elapsed time, minutes
      // T   :: clocktime, hours past previous midnight
      
      if ( clock_lights_out.valid )
	{
	  
	  double t0 = clock_start.hours();
	  
	  // ensure all are yoked to the same midnight as T0	  
	  double t1 = clock_lights_out.hours();
	  if ( t1 < t0 ) t1 += 24.0;

	  double t2 = any_sleep ? clock_sleep_onset.hours() : 0 ;
	  if ( t2 < t0 ) t2 += 24.0;
	  
	  double t3 = any_sleep ? clock_sleep_midpoint.hours() : 0;
	  if ( t3 < t0 ) t3 += 24.0;

	  double t4 = any_sleep ? clock_wake_time.hours() : 0;
	  if ( t4 < t0 ) t4 += 24.0;

	  double t5 = clock_lights_on.hours();
	  if ( t5 < t0 ) t5 += 24.0;

	  double t6 = clock_stop.hours();
	  if ( t6 < t0 ) t6 += 24.0;
	  		  
	  // finally, if t0 is at mignight, or just after, we need to
	  // make it align to the *previous* midnight. 
	  // i.e. if t0 < 12pm, then shift everything by 24 hours

	  if ( t0 < 12 )
	    {
	      t0 += 24.0;
	      t1 += 24.0;
	      t2 += 24.0; // t2/t3/t4 may not be 
	      t3 += 24.0; // defined if no sleep,	      
	      t4 += 24.0; // but no probs, they will not be output
	      t5 += 24.0;
	      t6 += 24.0;	      
	    }
	  

	  // for annots
	  tp_0_start = 0LLU;
	  tp_1_lights_out = ( t1 - t0 ) * 60.0 * 60.0 * (double)globals::tp_1sec ; 

	  if ( any_sleep )
	    {
	      tp_2_sleep_onset = ( t2 - t0 ) * 60.0 * 60.0 * (double)globals::tp_1sec ;
	      tp_3_sleep_midpoint = ( t3 - t0 ) * 60.0 * 60.0 * (double)globals::tp_1sec ;
	      tp_4_final_wake = ( t4 - t0 ) * 60.0 * 60.0 * (double)globals::tp_1sec ;
	    }
	  tp_5_lights_on = ( t5 - t0 ) * 60.0 * 60.0 * (double)globals::tp_1sec ;
	  tp_6_stop = ( t6 - t0 ) * 60.0 * 60.0 * (double)globals::tp_1sec ;
	  
	  // outputs
	  writer.value(  "T0_START" , t0 );
	  writer.value(  "E0_START" , 0 );
	  
	  if ( any_dark ) 
	    {
	      writer.value(  "T1_LIGHTS_OFF" , t1 );
	      writer.value(  "E1_LIGHTS_OFF" , ( t1 - t0 ) * 60.0 );
	    }

	  if ( any_sleep ) 
	    {
	      writer.value(  "T2_SLEEP_ONSET" , t2 );
	      writer.value(  "E2_SLEEP_ONSET" , Helper::dbl2str( ( t2 - t0 ) * 60.0 , 3 ) );
	      
	      writer.value(  "T3_SLEEP_MIDPOINT" , t3  );
	      writer.value(  "E3_SLEEP_MIDPOINT" , Helper::dbl2str( ( t3 - t0 ) * 60.0, 3 ) );
	      
	      writer.value(  "T4_FINAL_WAKE" , t4 );
	      writer.value(  "E4_FINAL_WAKE" , Helper::dbl2str( ( t4 - t0 ) * 60.0 , 3 ) );
	    }

	  if ( any_dark )
	    {
	      writer.value(  "T5_LIGHTS_ON" , t5 );
	      writer.value(  "E5_LIGHTS_ON" , Helper::dbl2str( ( t5 - t0 ) * 60.0 , 3 ) );
	    }

	  writer.value(  "T6_STOP" , t6 );
	  writer.value(  "E6_STOP" , Helper::dbl2str( ( t6 - t0 ) * 60.0 , 3 ) );
	  
	  // same in HMS
	  writer.value(  "HMS0_START" , clock_start.as_string(':') );
	  
	  if ( any_dark )
	    writer.value(  "HMS1_LIGHTS_OFF" , clock_lights_out.as_string(':') );
	  
	  if ( any_sleep )
	    {
	      writer.value(  "HMS2_SLEEP_ONSET" , clock_sleep_onset.as_string(':') );
	      writer.value(  "HMS3_SLEEP_MIDPOINT" , clock_sleep_midpoint.as_string(':') );
	      writer.value(  "HMS4_FINAL_WAKE" , clock_wake_time.as_string(':') );
	    }
	  
	  if ( any_dark )
	    writer.value(  "HMS5_LIGHTS_ON" , clock_lights_on.as_string(':') );
	  
	  writer.value(  "HMS6_STOP" , clock_stop.as_string(':') );
	  
	}
      
      // NREM cycles
      
      if ( any_sleep )
	{
	  writer.value(  "NREMC" , num_nremc );
	  writer.value(  "NREMC_MINS" , nremc_mean_duration );
	}

      const double epoch_mins = timeline->epoch_length() / 60.0 ;

      // note: here we swap definitions of "TIB" and "TRT"
      // i.e. in the code, they are reversed.... 
      // should go back and make consistent some point soon!
      writer.value( "TRT" , TIB ); 
      writer.value( "TIB" , TRT );

      writer.value( "TGT" , TGT );
      writer.value( "TST" , TST );
      writer.value( "TST_PER" , TpST );
      writer.value( "TWT" , TWT );
      writer.value( "LOT" , mins[ "L" ] );
      writer.value( "OTHR" , mins[ "?" ] );
      if ( constrain_too_short != -1 ) 
	writer.value( "SHORT" , constrain_too_short );
      writer.value( "CONF" , n_conflicts );
      writer.value( "FIXED_SLEEP" , n_fixed ); // --> ?
      writer.value( "FIXED_WAKE" , n_ignore_wake ); // --> ?
      writer.value( "FIXED_LIGHTS" , n_lights_fixed ); // --> ?
      writer.value( "LOST" , n_lights_fixed_was_sleep * epoch_mins ); 
      writer.value( "SINS" , (int)starts_in_sleep );
      writer.value( "EINS" , (int)ends_in_sleep );

      if ( any_sleep )
	{
	  writer.value( "WASO" , WASO );

	  // nb. different definition used internally
	  writer.value( "SPT" , SPT - FWT );
	  
	  writer.value( "FWT" , FWT );
	  writer.value( "SOL" , slp_lat );

	  // SOL/FWT versions that include LightsOff
	  writer.value( "PRE" , t_presleep );
	  writer.value( "POST" , t_postsleep );
	  	  
	  // was SLP_EFF
	  writer.value( "SE" , slp_eff_pct );	  

	  // was SLP_EFF2 --> this is the new default SE (i.e. denom SPT) 
	  writer.value( "SME" , slp_eff2_pct );

	  // ignore
	  //writer.value( "SLP_MAIN_EFF" , slp_main_pct );

	  // only defined if there is at least some persistent sleep
	  if ( TpST > 0 ) 
	    {
	      writer.value( "SOL_PER" , per_slp_lat );      
	      // adjust for increased time for onset of persistent sleep veresus first sleep
	      writer.value( "SPT_PER" , SPT - FWT - ( per_slp_lat - slp_lat ) );
	    }
	  	  
	  //
	  // Sleep Fragmentation Index /  Stage Transition Index
	  //
	  
	  // SFI = # of transitions into W / TST
	  // STI = # of sleep-sleep transitions / TST
	  
	  int trans_to_w = transitions5[ NREM1 ][ WAKE ]
	    + transitions5[ NREM2 ][ WAKE ]
	    + transitions5[ NREM3 ][ WAKE ]
	    + transitions5[ REM ][ WAKE ];
	  writer.value( "SFI" , trans_to_w / (double)TST );
	  
	  int trans_within_sleep =
	    transitions5[ NREM1 ][ NREM2 ]
	    + transitions5[ NREM1 ][ NREM3 ]
	    + transitions5[ NREM1 ][ REM ]
	    + transitions5[ NREM2 ][ NREM1 ]
	    + transitions5[ NREM2 ][ NREM3 ]
	    + transitions5[ NREM2 ][ REM ]
	    + transitions5[ NREM3 ][ NREM1 ]
	    + transitions5[ NREM3 ][ NREM2 ]
	    + transitions5[ NREM3 ][ REM ]
	    + transitions5[ REM ][ NREM1 ]
	    + transitions5[ REM ][ NREM2 ]
	    + transitions5[ REM ][ NREM3 ];
	  writer.value( "TI_S" , trans_within_sleep / (double)TST );
	  
	  // REM - NREM transition only
	  int rem_nrem_trans = 
	    + transitions5[ NREM1 ][ REM ]
	    + transitions5[ NREM2 ][ REM ]
	    + transitions5[ NREM3 ][ REM ]
            + transitions5[ REM ][ NREM1 ]
            + transitions5[ REM ][ NREM2 ]
            + transitions5[ REM ][ NREM3 ];
          writer.value( "TI_RNR" , rem_nrem_trans / (double)TST );

	  // 3-class (NREM / REM / WAKE) - denom = SPT 
	  int s3_trans =
            + transitions5[ NREM1 ][ REM ] // NREM - REM
            + transitions5[ NREM2 ][ REM ]
            + transitions5[ NREM3 ][ REM ]
            + transitions5[ REM ][ NREM1 ] 
            + transitions5[ REM ][ NREM2 ]
            + transitions5[ REM ][ NREM3 ]	    
	    + transitions5[ REM ][ WAKE ]  // REM - WAKE
	    + transitions5[ WAKE ][ REM ] 
            + transitions5[ NREM1 ][ WAKE ] // NREM - REM
            + transitions5[ NREM2 ][ WAKE ]
            + transitions5[ NREM3 ][ WAKE ]
            + transitions5[ WAKE ][ NREM1 ] 
            + transitions5[ WAKE ][ NREM2 ]
            + transitions5[ WAKE ][ NREM3 ];	    
          writer.value( "TI_S3" , s3_trans / (double)(SPT - FWT) );
	  
	  
	  //
	  // REM latency
	  //
	  
	  if ( mins[ "R" ] > 0 )
	    {
	      writer.value( "REM_LAT" , rem_lat_mins );
	      writer.value( "REM_LAT2" , rem_lat_nowake_mins );
	    }
	}
                 

      // ignore for now... 
      /// metrics need normalization by sequence length
      if ( 0 && any_sleep ) 
	{
	  if ( runs_pv5 >= 0 ) writer.value( "RUNS" , runs_pv5 );
	  if ( runs_pv3 >= 0 ) writer.value( "RUNS3" , runs_pv3 );	  
	}
      
    }

  //
  // Stage density sliding window
  //

  if ( sliding_window )
    {
      
      // only add W, INC & ANCHOR strata if needed
      // otherwise just STAGE x MINS 
      const int nslide = window_width_epochs.size();
      
      std::map<int,std::map<std::string,std::vector<double> > >::const_iterator ww = window_stats.begin();
      while ( ww != window_stats.end() )
	{
	  
	  const int p_width = window_width_epochs[ ww->first ]; 
	  const int p_inc = window_increment_epochs[ ww->first ]; 
	  const int p_anchor = window_anchor[ ww->first ]; 

	  std::stringstream wstr;
	  wstr << "W" << p_width * epoch_mins << "_" << p_inc * epoch_mins << "_T" << p_anchor;

	  //if ( nslide > 1 ) 
	    writer.level( wstr.str() , "W" );

	  // for each stage
	  std::map<std::string,std::vector<double> >::const_iterator zz =ww->second.begin();
	  while ( zz != ww->second.end() )
	    {
	      const std::string & stg = zz->first;
	      writer.level( stg , globals::stage_strat );
	      
	      // for each time-point
	      int tidx = 0;
	      std::vector<double>::const_iterator tt = zz->second.begin();
	      while ( tt != zz->second.end() )
		{
		  writer.level( window_timer[ ww->first ][ tidx ] , "MINS" );
		  writer.value( "DENS" , *tt );
		  ++tidx;
		  ++tt;
		}
	      writer.unlevel( "MINS" );
	      ++zz;
	    }
	  writer.unlevel( globals::stage_strat );
	  ++ww;
	}
      //      if ( nslide > 1 ) 
	writer.unlevel( "W" );
    }

  
  //
  // LZW compression index, and sample entropy
  //

  if ( verbose && any_sleep )
    {

      // nb - this will include lights on/out, but those should not add much to complexity values
      
      std::vector<char> sc(stages.size(),'?');
      std::vector<char> sc3(stages.size(),'?'); // 3-class NR/R/W
      
      for (int e=0; e<stages.size(); e++)
	{
	  if      ( stages[e] == NREM1 ) { sc[e] = 'A'; sc3[e] = 'A'; } 
	  else if ( stages[e] == NREM2 ) { sc[e] = 'B'; sc3[e] = 'A'; } 
	  else if ( stages[e] == NREM3 || stages[e] == NREM4 ) { sc[e] = 'C'; sc3[e] = 'A'; } 
	  else if ( stages[e] == REM ) { sc[e] = 'D'; sc3[e] = 'B'; }
	  else if ( stages[e] == WAKE ) { sc[e] = 'E'; sc3[e] = 'C'; }       
	}
           
      
      // 5-class
      std::string seq( sc.begin() , sc.end() );
      double r5 = 0;
      lzw_t lzw( seq , &r5 );
      
      // 3-class
      std::string seq3( sc3.begin() , sc3.end() );      
      double r3 = 0;
      lzw_t lzw3( seq3 , &r3 );
      
      // also get smallest possible value for fixed length 
      const int n1 = seq.size();
      double r0 = 0;
      std::string seq0( n1 , 'A' );
      lzw_t lzw0( seq0 , &r0 );
      
      //  r0    = min possile ratio
      //  r3/r5 = obs
      //  1.0   = max 

      double scaled5 = ( r5 - r0 ) / ( 1.0 - r0 );
      double scaled3 = ( r3 - r0 ) / ( 1.0 - r0 );
      writer.value( "LZW" , scaled5 );
      writer.value( "LZW3" , scaled3 );
      
    }
  

  // mse_t se;
  // for (int m=1; m<=7; m++)                                                                                      
  //   {
  //     double mse = se.sampen( seq , m );
  //     std::cout << " Mse = " << m << "\t" << mse << "\n";
  //   }	
  
  //
  // NREM cycle summary stats
  //
  
  if ( verbose && any_sleep && extra_output )
    {
      writer.var( "NREMC_START" , "NREM cycle start epoch" );
      writer.var( "NREMC_NREM_MINS" , "NREM cycle NREM duration (mins)" );
      writer.var( "NREMC_REM_MINS" , "NREM cycle REM duration (mins)" );
      writer.var( "NREMC_OTHER_MINS" , "NREM cycle other duration (mins)" );
      writer.var( "NREMC_MINS" , "NREM cycle total duration (mins)" );
      writer.var( "NREMC_N" , "NREM cycle total duration (epochs)" );
    }
    

  //
  // Stage-stratified outputs
  //

  if ( verbose && any_sleep )
    {

      const std::vector<std::string> these_stages_without_n4 = { "N1", "N2", "N3", "NR", "R", "S", "W" , "?" , "L" , "WASO" } ;
      const std::vector<std::string> these_stages_with_n4 =    { "N1", "N2", "N3", "N4" , "NR", "R", "S", "W" , "?" , "L" , "WASO" } ;
      // with N2_ASC/DSC/FLT
      const std::vector<std::string> these_stages_without_n4_n2ad = { "N1", "N2", "N2A", "N2D", "N2F", "N3", "NR", "R", "S", "W" , "?" , "L" , "WASO" } ;
      const std::vector<std::string> these_stages_with_n4_n2ad =    { "N1", "N2", "N2A", "N2D", "N2F", "N3", "N4" , "NR", "R", "S", "W" , "?" , "L" , "WASO" } ;
      
      const std::vector<std::string> & these_stages = collapse_nrem34 ? 
	( output_n2_asc_desc ? these_stages_without_n4_n2ad : these_stages_with_n4_n2ad ) :
	( output_n2_asc_desc ? these_stages_with_n4_n2ad : these_stages_with_n4 );
      
      // get total NR stats
      mins[ "NR" ] = mins[ "N1" ] + mins[ "N2" ] + mins[ "N3" ] + mins[ "N4" ];
      pct[ "NR" ] = pct[ "N1" ] + pct[ "N2" ] + pct[ "N3" ] + pct[ "N4" ];
      pct2[ "NR" ] = pct2[ "N1" ] + pct2[ "N2" ] + pct2[ "N3" ] + pct2[ "N4" ];
      
      mins[ "S" ] = mins[ "N1" ] + mins[ "N2" ] + mins[ "N3" ] + mins[ "N4" ] + mins[ "R" ] ;
      pct[ "S" ] = pct[ "N1" ] + pct[ "N2" ] + pct[ "N3" ] + pct[ "N4" ] + pct[ "R" ];  // should be 1.0
      pct2[ "S" ] = pct2[ "N1" ] + pct2[ "N2" ] + pct2[ "N3" ] + pct2[ "N4" ] + pct2[ "R" ];  // should be 1.0 - SME 
 
      std::vector<std::string>::const_iterator ss = these_stages.begin();
      while ( ss != these_stages.end() )
	{
	  
	  // skip?
	  if ( ! output_n2_asc_desc ) 
	    {
	      if ( *ss == "N2A" || *ss == "N2D" || *ss == "N2F"  ) 
		{
		  ++ss;
		  continue;
		}
	    }

	  
	  writer.level( *ss , globals::stage_strat );
	  writer.value( "MINS" , mins[ *ss] );
	  
	  // sleep stage as % of TST
	  if ( *ss == "N1" || *ss == "N2" || *ss == "N3" || *ss == "N4" || *ss == "NR" || *ss == "R" 
	       || *ss == "S" || *ss == "N2A" || *ss == "N2D" || *ss == "N2F" )
	    {
	      writer.value( "PCT" , pct[ *ss] );
	      writer.value( "DENS" , pct2[ *ss] );
	    }
	  
	  writer.value( "BOUT_N" , bout_n[ *ss] );
	  writer.value( "BOUT_MX" , bout_max[ *ss] );
	  writer.value( "BOUT_MN" , bout_mean[ *ss] );
	  writer.value( "BOUT_MD" , bout_med[ *ss] );
	  writer.value( "BOUT_05" , bout_5[ *ss] );
	  writer.value( "BOUT_10" , bout_10[ *ss] );

	  //
	  // dur/time trackers (if this stage class was included)
	  //
	  
	  if ( stg_timing_all.find( *ss ) != stg_timing_all.end() )	    
	    writer.value( "TA" , stg_timing_all[ *ss ] );
	  if ( stg_timing_sleep.find( *ss ) != stg_timing_sleep.end() )	    
	    writer.value( "TS" , stg_timing_sleep[ *ss ] );
	  
	  // next stage/class type
	  ++ss;
	}
           
      writer.unlevel( globals::stage_strat );      


      //
      // Bouts
      //
      
      if ( extra_output )
	{
	  int bn = 0;
	  std::set<bout_t>::const_iterator bb = bouts.begin();
	  while ( bb != bouts.end() )
	    {
	      writer.level( ++bn ,"N" );
	      
	      const int e1 = epoch_n[ bb->start ] ;
	      const int e2 = epoch_n[ bb->stop ] ;
	      
	      clocktime_t ct1 = timeline->edf->header.starttime;
	      ct1.advance_seconds( timeline->epoch_length() * e1 );
	      
	      // nb. +1 to get to /end/ of the last epoch
	      clocktime_t ct2 = timeline->edf->header.starttime;
	      ct2.advance_seconds( timeline->epoch_length() * ( e2 + 1 ) ) ;
	      
	      if ( bb->ss == NREM2 ) 
		writer.value( "STAGE" , "NR" );	  
	      else
		writer.value( "STAGE" , globals::stage( bb->ss ) );	    
	      
	      writer.value( "FIRST_EPOCH" , e1 + 1 );
	      writer.value( "LAST_EPOCH" , e2 + 1 );
	      
	      writer.value( "START" , ct1.as_string() );
	      writer.value( "STOP"  , ct2.as_string() );
	      
	      writer.value( "MINS" , ( ( e2 - e1 + 1 ) * timeline->epoch_length() ) / 60.0 );
	      
	      ++bb;
	    }
	  writer.unlevel( "N" );
	}
      
    }
  
  
  //
  // Cycle-stratified outputs (verbose mode only), and transitions
  //

  if ( verbose && extra_output ) 
    {
      
      if ( any_sleep ) 
	{
	  
	  std::map<int,double>::iterator cc = nremc_duration.begin();
	  while ( cc != nremc_duration.end() )
	    {
	      writer.level( cc->first , globals::cycle_strat );
	      
	      writer.value(  "NREMC_START" , nremc_start_epoch[ cc->first ] );
	      writer.value(  "NREMC_NREM_MINS" , nremc_nrem_duration[ cc->first ] );
	      writer.value(  "NREMC_REM_MINS" , nremc_rem_duration[ cc->first ] );
	      writer.value(  "NREMC_OTHER_MINS" , cc->second - nremc_nrem_duration[ cc->first ] - nremc_rem_duration[ cc->first ] );
	      writer.value( "NREMC_MINS" , cc->second );
	      writer.value( "NREMC_N" , nremc_epoch_duration[ cc->first ] );
	      
	      ++cc;
	    }
	  
	  writer.unlevel( globals::cycle_strat );
	  
	  //
	  // Transitions
	  //
	  
	  std::vector<sleep_stage_t> ss = { NREM1 , NREM2 , NREM3 , REM , WAKE };
	  std::vector<std::string> ss_str = { "N1" , "N2" , "N3" , "R" , "W" };

	  std::vector<sleep_stage_t> ss3 = { NREM2 , REM , WAKE };
	  std::vector<std::string> ss3_str = { "NR" , "R" , "W" };
	  
	  if ( flanking_3class ) 
	    {
	      ss = ss3; 
	      ss_str = ss3_str; 
	    }

	  std::map<sleep_stage_t,int> marg_pre, marg_post;
	  int tot = 0;
          for (int ss1=0; ss1<ss.size(); ss1++)	    
	    for (int ss2=0; ss2<ss.size(); ss2++)
	      {
		tot += transitions[ ss[ss1] ][ ss[ss2] ];
		marg_pre[ ss[ss1] ] += transitions[ ss[ss1] ][ ss[ss2] ];
		marg_post[ ss[ss2] ] += transitions[ ss[ss1] ][ ss[ss2] ];
	      }

	  for (int ss1=0; ss1<ss.size(); ss1++)
	    {
	      writer.level( ss_str[ss1] , "PRE" );

	      for (int ss2=0; ss2<ss.size(); ss2++)
		{
		  writer.level(  ss_str[ss2] , "POST" );
		  writer.value( "N" , transitions[ ss[ss1] ][ ss[ss2] ] );

		  // probabilities: joint 
		  if ( tot > 0 ) 
		    writer.value( "P" , transitions[ ss[ss1] ][ ss[ss2] ] / (double)tot );
		  
		  // P( post | pre ) 
		  if ( marg_pre[ ss[ss1] ] > 0 )
		    writer.value( "P_POST_COND_PRE" , transitions[ ss[ss1] ][ ss[ss2] ] / (double)marg_pre[ ss[ss1] ] );

		  // P( pre | post )
		  if ( marg_post[ ss[ss2] ] > 0 ) 
		    writer.value( "P_PRE_COND_POST" , transitions[ ss[ss1] ][ ss[ss2] ] / (double)marg_post[ ss[ss2] ] );

		}
	      writer.unlevel( "POST" );
	    }
	  writer.unlevel( "PRE" );
	  
	}

    }

  
  //
  // Per epoch level output
  //


  // stage information and time only in non-verbose mode

  std::map<sleep_stage_t,int> stagen;
  stagen[ WAKE ] = 1;
  stagen[ REM ] = 0;
  stagen[ NREM1 ] = -1;
  stagen[ NREM2 ] = -2;
  stagen[ NREM3 ] = -3;
  stagen[ NREM4 ] = collapse_nrem34 ? -3 : -4;
  
  // all 'bad' here -- treat as 'UNKNOWN'
  stagen[ UNKNOWN ] = 2; 
  stagen[ UNSCORED ] = 2; // these others should not happen, but in case...
  stagen[ MOVEMENT ] = 2;
  stagen[ ARTIFACT ] = 2;
  stagen[ LIGHTS_ON ] = 3;


  const int ne = timeline->num_epochs();
  
  clocktime_t starttime( clock_start );
  
  
  //
  // output in non-verbsoe mode (STAGES command)
  //
  
  if ( ! verbose )
    {

      if ( eannot == "." )
	{
	  logger << "  writing epoch-level sleep stages to standard out\n";
	  for (int e=0;e<ne_gaps;e++)
	    if ( ! epoch_gap[e] )
	      std::cout << globals::stage( stages[ e ] ) << "\n";	
	  return;
	}      
      else if ( eannot != "" )
	{
	  logger << "  writing epoch-level sleep stages to " << eannot << "\n";
	  std::ofstream EOUT( Helper::expand( eannot ).c_str() , std::ios::out );
	  for (int e=0;e<ne_gaps;e++)
	    if ( ! epoch_gap[e] )		      
	      EOUT << globals::stage( stages[ e ] ) << "\n";
	  EOUT.close();
	  return;
	}
      
      // Typical STAGE command

      // actual existing epoch count
      int ecnt = 0;
      
      // elapsed time from start of first stage epoch
      double mins = 0;
      
      for (int e=0;e<ne_gaps;e++)
	{
	  
	  // skip gaps in the output
	  const bool is_gap = epoch_gap[e];
	  if ( is_gap )
	    {
	      mins += epoch_dur[e];
	      continue;	      
	    }
	  
	  // get actual epoch number
	  const int eidx = epoch_n[e];
	  
	  // epoch-level stratification
	  // epoch_n is the original 0-based epoch encoding
	  // so, for diplay +1 , but for other calculations
	  // we want to keep this original encoding
	  
	  //std::cout << " e " << e << " / " << ne << " --> " << epoch_n.size() << " " << epoch_n[e] << "\n";
	  
	  writer.epoch( eidx + 1 );
	  
	  // new - use actual epoch encoding (it's what it's there for!)
	  interval_t interval = timeline->epoch( ecnt );
	  ++ecnt;
	  
	  const double sec0 = interval.start * globals::tp_duration;

	  // clock time based on EDF header	  
	  if ( starttime.valid ) 
	    {

	      // old - assume epoch 1 starts at 0 / EDF start
	      //clocktime_t current_clock_time = starttime;	      
	      //current_clock_time.advance_seconds( epoch_sec * epoch_n[e] );
	      
              clocktime_t present = starttime;
              present.advance_seconds( sec0 );
	      
	      //std::cout << sec0 << " is sec0 " << present.as_string( ':' ) << "\n";
	      
	      writer.value( "CLOCK_TIME" , present.as_string( ':' ) );
	      
	      if ( verbose ) 
		writer.value( "CLOCK_HOURS" ,  present.as_numeric_string() );
	      
	    }
	  
	  // time in minutes (from start of stage-aligned epochs)	  
	  writer.value( "MINS" ,  mins / 60.0 ); // eidx * epoch_mins );
	  
	  // time from EDF start (seconds)
	  writer.value( "START_SEC" , sec0 ); 
	  
	  // stages	  
	  writer.value( "STAGE" , globals::stage( stages[e] ) );
	  
	  // i.e. prior to anything being set to L or ?
	  writer.value( "OSTAGE" , globals::stage( original_stages[e] ) );
	  
	  writer.value( "STAGE_N" , stagen[ stages[e] ] );

	  // track mins
	  mins += epoch_dur[e];
	  
	}
      
      writer.unepoch();
      
        
      return;

    }

  

  //
  // ... otherwise, the rest of this function is verbose mode only
  //

  
  // Outputs
  // Per epoch, we have
  //   a) stage (done above)
  //   b) elapsed time
  //   c) elapsed sleep
  //   d) period number
  //   e) N2 measure of direction

  
  //
  // output epoch level data ?
  //   if not, quit - unless we also need to make annotations afterwards
  //

  const bool annotate_features = annot_prefix != "";
  
  if ( ! epoch_lvl_output )
    if ( ! annotate_features ) 
      return;

  
  double elapsed_n1 = 0 , elapsed_n2 = 0 , elapsed_n34 = 0 , elapsed_rem = 0;
    
  double elapsed_sleep = 0 , elapsed_wake = 0 , elapsed_waso = 0 ;

  //  std::cout << " ne ne_gaps " << ne << " " << ne_gaps << "\n";
  
  double elapsed_mins_from_epoch1 = 0;
  
  for (int e=0;e<ne_gaps;e++)
    {
      
      // skip gaps in the output
      const bool is_gap = epoch_gap[e];
      
      if ( is_gap )
	{
	  // for MINS output
	  elapsed_mins_from_epoch1 += epoch_dur[e];

	  // here, implies GAP counts of ? unknown
	  // i.e. do not increment any elapsed_* counters
	  continue;
	}
      
      // get actual epoch number
      const int eidx = epoch_n[e];
      
      // epoch-level stratification
      
      writer.epoch( timeline->display_epoch( eidx ) );
      
      // new - use actual epoch encoding (it's what it's there for!)                                              
      interval_t interval = timeline->epoch( eidx );

      const double sec0 = interval.start * globals::tp_duration;
      
      if ( starttime.valid ) 
	{
	  
	  clocktime_t present = starttime;
	  present.advance_seconds( sec0 );
	  
	  writer.value( "CLOCK_TIME" , present.as_string( ':' ) );
	  
	  if ( verbose )
	    writer.value( "CLOCK_HOURS" ,  present.as_numeric_string() );
	  
	}

      // time in minutes (from EPOCH 1, not EDF start, i.e. if EPOCH align)
      writer.value( "MINS" ,  elapsed_mins_from_epoch1 / 60.0 ); // e * epoch_mins );
      elapsed_mins_from_epoch1 += epoch_dur[e];

      // flag if comes after a GAP
      writer.value( "AFTER_GAP" , e != 0 && epoch_gap[e-1] ? 1 : 0 );
      
      // time from EDF start (seconds)
      writer.value( "START_SEC" , sec0 );

      // stages      
      writer.value( "STAGE" , globals::stage( stages[e] ) );    
      writer.value( "OSTAGE" , globals::stage( original_stages[e] ) );    
      writer.value( "STAGE_N" , stagen[ stages[e] ] );
      

      // stage stats
      writer.value( "E_WAKE" , elapsed_wake );
      writer.value( "E_WASO" , elapsed_waso );
      writer.value( "E_SLEEP" , elapsed_sleep );
      writer.value( "E_N1" , elapsed_n1 );
      writer.value( "E_N2" , elapsed_n2 );
      writer.value( "E_N3" , elapsed_n34 );
      writer.value( "E_REM" , elapsed_rem );		  
      
      // and as percentages
      writer.value( "PCT_E_SLEEP" , TST>0 ? elapsed_sleep / TST : 0 );

      writer.value( "PCT_E_N1" , mins[ "N1" ] > 0 ? elapsed_n1 / mins[ "N1" ] : 0 );
      writer.value( "PCT_E_N2" , mins[ "N2" ] > 0 ? elapsed_n2 / mins[ "N2" ] : 0 );
      writer.value( "PCT_E_N3" , (mins["N3"] + mins["N4"]) > 0 ? elapsed_n34 / (mins["N3"]+mins["N4"]) : 0 );
      writer.value( "PCT_E_REM" , mins["R"] > 0 ? elapsed_rem / mins["R"] : 0 );

      // track if making annots?
      if ( annotate_features )
	{
	  elapsed_stg_sec[ "N1" ].push_back( elapsed_n1 );
	  elapsed_stg_sec[ "N2" ].push_back( elapsed_n2 );
	  elapsed_stg_sec[ "N3" ].push_back( elapsed_n34 );
	  elapsed_stg_sec[ "NR" ].push_back( elapsed_n1 + elapsed_n2 + elapsed_n34 );
	  elapsed_stg_sec[ "R" ].push_back( elapsed_rem );
	  elapsed_stg_sec[ "S" ].push_back( elapsed_sleep );
	  elapsed_stg_sec[ "W" ].push_back( elapsed_wake );
	  elapsed_stg_sec[ "WASO" ].push_back( elapsed_waso );
	  
	  elapsed_stg_rel[ "N1" ].push_back( mins[ "N1" ] > 0 ? elapsed_n1 / mins[ "N1" ] : 0 );
          elapsed_stg_rel[ "N2" ].push_back( mins[ "N2" ] > 0 ? elapsed_n2 / mins[ "N2" ] : 0 );
          elapsed_stg_rel[ "N3" ].push_back( (mins["N3"] + mins["N4"]) > 0 ? elapsed_n34 / (mins["N3"]+mins["N4"]) : 0 ) ;
	  elapsed_stg_rel[ "NR" ].push_back( (mins["N1"] + mins["N2"] + mins["N3"] + mins["N4"]) > 0 ? (elapsed_n1+elapsed_n2+elapsed_n34) / (mins["N1"] + mins["N2"]+mins["N3"]+mins["N4"]) : 0 ) ;
	  elapsed_stg_rel[ "R" ].push_back( mins["R"] > 0 ? elapsed_rem / mins["R"] : 0 );
          elapsed_stg_rel[ "S" ].push_back( mins[ "S" ] > 0 ? elapsed_sleep / mins[ "S" ] : 0  );
          elapsed_stg_rel[ "W" ].push_back( mins[ "W" ] > 0 ? elapsed_wake / mins[ "W" ] : 0  );
          elapsed_stg_rel[ "WASO" ].push_back( mins[ "WASO" ] > 0 ? elapsed_waso / mins[ "WASO" ] : 0  );
	}
      
	  
      
      // track elapsed time
      if ( stages[e] == WAKE ) 
	{
	  elapsed_wake += epoch_mins;
	  if ( e > first_sleep_epoch && e < final_wake_epoch ) elapsed_waso += epoch_mins;
	}
      else if ( stages[e] == NREM1 ) { elapsed_sleep += epoch_mins; elapsed_n1  += epoch_mins; }
      else if ( stages[e] == NREM2 ) { elapsed_sleep += epoch_mins; elapsed_n2  += epoch_mins; }
      else if ( stages[e] == NREM3 || stages[e] == NREM4 ) { elapsed_sleep += epoch_mins; elapsed_n34 += epoch_mins; }
      else if ( stages[e] == REM ) { elapsed_sleep += epoch_mins; elapsed_rem += epoch_mins; }
      
      // persistent sleep
      
      writer.value( "PERSISTENT_SLEEP" , in_persistent_sleep[e] );

      // in/out of sleep period
      
      writer.value( "SPT" , (int)(e >= first_sleep_epoch && e < final_wake_epoch) );
      writer.value( "POST" , (int)(e >= final_wake_epoch) );
      writer.value( "PRE" , (int)(e < first_sleep_epoch ) );
      
      // cycles
      
      if ( sleep_cycle_number[e] )
	{
	  writer.value( "CYCLE" , sleep_cycle_number[e] );
	  writer.value( "PERIOD" , sleep_code[e] == 5 ? "REMP" : sleep_code[e]==1 ? "NREMP" : "." ) ;
	  writer.value( "CYCLE_POS_REL" , cycle_pos_relative[e] );
	  writer.value( "CYCLE_POS_ABS" , cycle_pos_absolute[e] );
	}
      
      // flanking epochs

      writer.value( "FLANKING" , flanking[e] );
      writer.value( "FLANKING_ALL" , flanking_tot[e] );
      writer.value( "NEAREST_WAKE" , nearest_wake[e] );
      writer.value( "WASO" , is_waso[e] );

      writer.value( "TR_NR2R" , nrem2rem[e] );
      writer.value( "TOT_NR2R" , nrem2rem_total[e] );
      writer.value( "TR_NR2W" , nrem2wake[e] );
      writer.value( "TOT_NR2W" , nrem2wake_total[e] );

      writer.value( "TR_R2NR" , rem2nrem[e] );
      writer.value( "TOT_R2NR" , rem2nrem_total[e] );
      writer.value( "TR_R2W" , rem2wake[e] );
      writer.value( "TOT_R2W" , rem2wake_total[e] );

      writer.value( "TR_W2NR" , wake2nrem[e] );
      writer.value( "TOT_W2NR" , wake2nrem_total[e] );
      writer.value( "TR_W2R" , wake2rem[e] );
      writer.value( "TOT_W2R" , wake2rem_total[e] );

      // N2 ascending/descending status
      
      if ( stages[e] == NREM2 ) 
	writer.value( "N2_WGT" , n2_ascdesc[e] );
    
           
    } // next epoch

  writer.unepoch();



  //
  // Add annotation to denote multiple hypnogram features, e.g. including NREM cycle?
  // Note - this is done *after* output (meaning that annot will only work w/ HYPNO, not stages)
  //
          
  if ( annotate_features )
    {
      annotate( annot_prefix , annot_suffix );      
    }

   

}

void dummy_hypno()
{
  annotation_set_t annotations;
  edf_t edf( &annotations );
  param_t param;

  // dummy values
    
  hypnogram_t h;
  h.timeline = &edf.timeline;
  
  while( ! std::cin.eof() )
    {
      std::string s;
      std::cin >> s;
      if ( std::cin.eof() ) break;
      if ( s == "W" ) h.stages.push_back( WAKE );
      else if ( s == "N1" ) h.stages.push_back( NREM1 );
      else if ( s == "N2" ) h.stages.push_back( NREM2 );
      else if ( s == "N3" ) h.stages.push_back( NREM3 );
      else if ( s == "N4" ) h.stages.push_back( NREM4 );
      else if ( s == "R"  ) h.stages.push_back( REM );
      else if ( s == "L"  ) h.stages.push_back( LIGHTS_ON );
      else if ( s == "?"  ) h.stages.push_back( UNKNOWN );
      else logger << "did not recognize " << s << "\n";
    }

  logger << "read " << h.stages.size() << "\n";

  edf.header.starttime = "10:00:00";
  
  // fudge so it works in this non-standard case...
  edf.id = "_DUMMY_";
  h.fudge( 30 , h.stages.size() );

  // make a copy of stages 
  h.original_stages = h.stages;
  
  h.edit( h.timeline , param );
  if ( h.empty() ) 
    {
      logger << " ** warning, no valid stage-annotated epochs, bailing...\n"; 
      return; // bail if empty
    }
  h.calc_stats( true );
  h.output( true , true , true ); // verbose mode == T 

}

void hypnogram_t::fudge( double es, int ne )
{
  timeline->epoch_length_tp = es * globals::tp_1sec;
  timeline->epochs.resize( ne );
}



void hypnogram_t::do_slide()
{

  window_stats.clear();
  window_timer.clear();

  // requires at least some sleep:
  
  bool sleep = false;
  for (int e=0; e<ne; e++)
    if ( is_sleep( stages[e] ) )
      {
	sleep = true;
	break;
      }
  if ( ! sleep ) return;
  
  // consider each window specification

  const int nslide = window_width_epochs.size();

  for (int i=0; i<nslide; i++)
    {
      const int wid = window_width_epochs[i];
      const int inc = window_increment_epochs[i];
      const int anc = window_anchor[i];
           
      // std::map<int,std::map<std::string,std::vector<double> > > > > window_stats;
      // std::map<int,std::vector<double> > > > window_timer;

      // build window_timer, and start/stop epoch codes
      std::vector<double> wt;
      std::vector<int> e1, e2;

      int align = +1; // forwards
      if ( anc == 3 ) align = 0; // outwards both ways
      if ( anc > 3 ) align = -1; // backwards
      
      if ( wid == 0 || inc == 0 ) 
	Helper::halt( "sliding window width and increment must be >0" );
            
      // nb, below b is +1 past end 
      if ( anc < 3 ) 
	{	  
	  int a0 = 0; // EDF start
	  if ( anc == 1 ) a0 = lights_out_epoch;
	  else if ( anc == 2 ) a0 = first_sleep_epoch;
	  
	  // end point
	  int end = ne;
	  if ( anc == 1 ) end = lights_on_epoch;
	  if ( anc == 2 ) end = final_wake_epoch;
	  
	  for (int a=a0; a < ne; a += inc)
	    {
	      int b = a + wid;
	      if ( b <= end )
		{
		  e1.push_back( a );
		  e2.push_back( b );
		  wt.push_back( epoch_mins * ( ( a + b )/2.0 - a0 ) ); 
		}
	    }	  
	}

      // reverse
      if ( anc > 3 )
        {
          // "start" is always +1 past 
	  int a0 = ne; // EDF end (+1)
          if ( anc == 5 ) a0 = lights_on_epoch;
          else if ( anc == 4 ) a0 = final_wake_epoch;

          // end point (left most), 0-based
	  int end = 0;
          if ( anc == 5) end = lights_out_epoch;
          if ( anc == 4) end = first_sleep_epoch;

	  // b - a 
          for (int a=a0; a >= 0; a -= inc)
            {
              int b = a - wid; // left
              if ( b >= end )   // exact check
                {
		  //std::cout << " adding " << b << " to " << a << "\n";
                  e1.push_back( b ); // order flipped
                  e2.push_back( a ); 
                  wt.push_back( epoch_mins * ( ( a + b )/2.0 - a0 ) );
                }
            }
        }

      // midpoint based, delimited by SPT
      if ( anc == 3 ) 
	{
	  //std::cout << " epochs = " << first_sleep_epoch << " " << sleep_midpoint_epoch << " " << final_wake_epoch << "\n";
	  
	  // shift by half window, to cover (as here, right is aligned to MID)
	  int a0 = sleep_midpoint_epoch + wid/2 ;
	  int end = first_sleep_epoch;
	  
	  // nb. MINS is based on center, so wt calc is changed...  starts at 0 not w/2

	  // first back
	  for (int a=a0; a >= 0; a -= inc)
            {
              int b = a - wid; // left

	      if ( b >= end )   // exact check
		{                     		  
		  e1.push_back( b ); 
                  e2.push_back( a );
                  wt.push_back( epoch_mins * ( ( a + b )/2.0 - sleep_midpoint_epoch ) );		  
                }
            }
	  
	  // flip to get right order
	  std::vector<int> r1 = e1;
	  std::vector<int> r2 = e2;
	  std::vector<double> rt = wt;
	  const int nn = e1.size();
	  for (int i=0;i<nn;i++)
	    {
	      e1[nn-1-i] = r1[i];
	      e2[nn-1-i] = r2[i];
	      wt[nn-1-i] = rt[i];
	    }

	  // then forward, note, the same shift - 
	  a0 = sleep_midpoint_epoch - wid/2;
	  // but we don't need to repeat the zero-aligned one, so shift now
	  a0 += inc;

	  end = final_wake_epoch;
	  for (int a=a0; a < ne; a += inc)
            {
              int b = a + wid;
              if ( b <= end )
                {
		  e1.push_back( a );
                  e2.push_back( b );
                  wt.push_back( epoch_mins * ( ( a + b )/2.0 - sleep_midpoint_epoch ) );		  
                }
            }
	}


      // store times
      window_timer[i] = wt;

      // get densities      
      const std::vector<std::string> these_stages
	= { "N1", "N2", "N3", "NR", "R", "S", "WASO","W" } ;
      
      const int nw = wt.size();

      std::vector<std::string>::const_iterator qq = these_stages.begin();
      while ( qq != these_stages.end() )
	{
	  
	  sleep_stage_t stage = WAKE;
	  if ( *qq == "N1" ) stage = NREM1;
	  if ( *qq == "N2" ) stage = NREM2;
	  if ( *qq == "N3" ) stage = NREM3;
	  if ( *qq == "N4" ) stage = NREM4;
	  if ( *qq == "R" ) stage = REM;
	  // special cases
	  bool all_nrem = *qq == "NR";
	  bool all_sleep = *qq == "S";
	  bool waso = *qq == "WASO";

	  std::vector<double> dens( nw );
	  
	  for (int j=0; j<nw; j++)
	    {
	      double d = 0;
	      for (int e=e1[j];e<e2[j];e++)
		{
		  bool match = false;
		  if ( all_nrem ) 
		    match = stages[e] == NREM1 || stages[e] == NREM2 || stages[e] == NREM3 || stages[e] == NREM4 ;
		  else if ( all_sleep ) 
		    match = stages[e] == NREM1 || stages[e] == NREM2 || stages[e] == NREM3 || stages[e] == NREM4 || stages[e] == REM ;
		  else if ( waso ) 
		    match = stages[e] == WAKE && e >= first_sleep_epoch && e < final_wake_epoch;
		  else 
		    match = stages[e] == stage;
		  
		  if ( match ) d++;
		}
	      
	      // normalize to 0..1 rate
	      d /= (double)( e2[j] - e1[j] );
	      
	      dens[ j ] = d ; 
	      
	    }
	  
	  // store
	  window_stats[ i ][ *qq ] = dens;
	  ++qq;
	}
    }
}

