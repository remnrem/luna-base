

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

extern writer_t writer;

extern logger_t logger;


bool is_rem( sleep_stage_t s ) { return s == REM; } 
bool is_nrem( sleep_stage_t s ) { return s == NREM1 || s == NREM2 || s == NREM3 || s == NREM4; } 
bool is_nrem1( sleep_stage_t s ) { return s == NREM1; } 
bool is_nrem23( sleep_stage_t s ) { return s == NREM2 || s == NREM3; } 
bool is_nrem234( sleep_stage_t s ) { return s == NREM2 || s == NREM3 || s == NREM4; } 
bool is_wake( sleep_stage_t s ) { return s == WAKE; }
bool is_wake_or_lights( sleep_stage_t s ) { return s == WAKE || s == LIGHTS_ON; } 
bool is_sleep( sleep_stage_t s ) { return s == NREM1 || s == NREM2 || s == NREM3 || s == NREM4 || s == REM ; } 


bool hypnogram_t::construct( timeline_t * t , const bool verbose , const std::vector<std::string> & s )
{ 
  timeline = t;
  if ( s.size() != timeline->num_total_epochs() ) 
    Helper::halt( "bad number of stages, " 
		  + Helper::int2str( (int)s.size() ) 
		  + " but expecting " 
		  + Helper::int2str( timeline->num_total_epochs() ) );    
  stages.resize( s.size() );
  for (int e=0;e<s.size();e++) stages[e] = globals::stage( s[e] );
  calc_stats( verbose );
  return true;
} 

bool hypnogram_t::construct( timeline_t * t , const bool verbose , const std::string sslabel ) 
{

  // point to 'parent' timeline
  timeline = t ;
  
  // get handle
  annot_t * annot = timeline->annotations( sslabel );
  if ( annot == NULL ) 
    {
      logger << " could not find any valid sleep stage annotations... quitting\n";
      return false;
      //Helper::halt( "[" + sslabel + "] not set" );
    }

  //
  // set internal, epoch-level annotations used by timeline
  //
  
  std::set<std::string> values;
  values.clear(); values.insert( "wake" );
  timeline->annotate_epochs(  globals::stage( WAKE ) , "SleepStage" , values );
  
  values.clear(); values.insert( "NREM1" );
  timeline->annotate_epochs(  globals::stage( NREM1  )  , "SleepStage" , values );

  values.clear(); values.insert( "NREM2" );
  timeline->annotate_epochs(  globals::stage( NREM2  )  , "SleepStage" , values );

  values.clear(); values.insert( "NREM3" );
  timeline->annotate_epochs(  globals::stage( NREM3  )  , "SleepStage" , values );

  values.clear(); values.insert( "NREM4" );
  timeline->annotate_epochs(  globals::stage( NREM4 )  , "SleepStage" , values );

  values.clear(); values.insert( "REM" );
  timeline->annotate_epochs(  globals::stage( REM ) , "SleepStage" , values );


  //
  // If we've masked the data, epoch count may not start at 0...
  // Although it probably doesn't make sense to use HYPNO then, 
  // we still want to be able to dump STAGE information easily w/ the
  // STAGE command...
  //
  

  // in VERBOSE (HYPNO) mode, we require the FULL epoch set
  
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

  const int ne = timeline->num_total_epochs();
  
  timeline->first_epoch();

  stages.clear();
  epoch_n.clear();

  //
  // We should be able to use current 0..ne epoch naming as epoch-annotations
  // still work after a restructure
  //
  
  while ( 1 ) 
    {

      int e = timeline->next_epoch();
      
      if ( e == -1 ) break;

      // for output of STAGES or HYPNO, use original EDF annotations though
      int e2 = timeline->original_epoch(e) ;
      
      bool wake = timeline->epoch_annotation( "wake"  , e );
      bool n1   = timeline->epoch_annotation( "NREM1" , e );
      bool n2   = timeline->epoch_annotation( "NREM2" , e );
      bool n3   = timeline->epoch_annotation( "NREM3" , e );
      bool n4   = timeline->epoch_annotation( "NREM4" , e );
      bool rem  = timeline->epoch_annotation( "REM"   , e );      

      bool other = ! ( wake || n1 || n2 || n3 || n4 || rem );
      bool conflict = ( (int)wake + (int)n1 + (int)n2 + (int)n3 + (int)n4 + (int)rem ) > 1;
      if ( conflict ) other = true;

      if      ( conflict ) stages.push_back( UNSCORED );
      else if ( other ) stages.push_back( UNSCORED );
      else if ( wake ) stages.push_back( WAKE );
      else if ( n1 ) stages.push_back( NREM1 );
      else if ( n2 ) stages.push_back( NREM2 );
      else if ( n3 ) stages.push_back( NREM3 );
      else if ( n4 ) stages.push_back( NREM4 );
      else if ( rem ) stages.push_back( REM );
      else stages.push_back( UNSCORED );
      
      // store original EDF 0-based encoding, to be passed to calc_stats()
      epoch_n.push_back( e2 );
      
    }

   calc_stats( verbose );

   return true;
}   


void hypnogram_t::calc_stats( const bool verbose )
{

  //
  // epoch size (in minutes) and number
  //

  const double epoch_mins = timeline->epoch_length() / 60.0 ; 
  
  const int ne = stages.size();
  

  //
  // Recode any leading/trailing "?" as "L"
  //
  
  for (int e =0; e < ne ; e++)
    {
      if ( stages[e] == UNSCORED ) stages[e] = LIGHTS_ON;
      if ( stages[e] != UNSCORED && stages[e] != LIGHTS_ON ) break;
    }
  
  for (int e = ne - 1 ; e != 0 ; e--)
    {
      if ( stages[e] == UNSCORED ) stages[e] = LIGHTS_ON;
      if ( stages[e] != UNSCORED && stages[e] != LIGHTS_ON ) break;
    }

  
  //
  // Basic summary statistics per-individual/night
  //

  mins_wake = mins_n1 = mins_n2 = mins_n3 = mins_n4 = mins_rem = mins_other = 0;
  
  // implicitly, this will only count in the TRT (i.e. ignore pre
  // lights-out, and post lights-on)

  for (int e = 0 ; e < ne ; e++ )
    {
      if      ( stages[e] == WAKE  ) mins_wake += epoch_mins;
      else if ( stages[e] == NREM1 ) mins_n1 += epoch_mins;
      else if ( stages[e] == NREM2 ) mins_n2 += epoch_mins;
      else if ( stages[e] == NREM3 ) mins_n3 += epoch_mins;
      else if ( stages[e] == NREM4 ) mins_n4 += epoch_mins;
      else if ( stages[e] == REM   ) mins_rem += epoch_mins;
      else if ( stages[e] != LIGHTS_ON ) mins_other += epoch_mins; // movement, artifact, unscored, but only 
    }

  // did we observe /any/ sleep?

  any_sleep = ( mins_n1 + mins_n2 + mins_n3 + mins_n4 + mins_rem ) > 0 ;

  // lights out/on check
  // i.e. L can only be at start and end of recording
  // so illegal to have L flanked by non-L on both sides : 010
  bool lights_back_on = false;
  for (int e=1;e<ne-1;e++) 
    {
      // went from lights off to on again 
      if ( stages[e-1] != LIGHTS_ON && stages[e]   == LIGHTS_ON ) lights_back_on = true;

      // if we previously have gone from OFF->ON, we cannot go back OFF again--- i.e. only allow 
      // a single LIGHTS OFF interval
      
      if ( lights_back_on && stages[e] == LIGHTS_ON && stages[e+1] != LIGHTS_ON ) 
	Helper::halt( "LIGHTS_ON periods can only be at start and end of recording" );
    }

  
  // lights out/on
  int lights_out_epoch = 0;
  for (int e=0;e<ne-1;e++) 
    {
      if ( stages[e] != LIGHTS_ON ) { lights_out_epoch = e; break; } 
    }

  int lights_on_epoch = ne; // by default, one past the end  
  for (int e=ne-1;e>0;e--) 
    if ( stages[e] != LIGHTS_ON ) { lights_on_epoch = e+1 ;  break; }

  
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
  
  // total time in bed
  TIB = ne * epoch_mins;
  
  // total recording time (i.e. only from lights out, lights on)
  // note; lights_out_epoch is defined as 1 past end, so no +1
  int TRT_total_epochs = lights_on_epoch - lights_out_epoch ;
  TRT =  TRT_total_epochs * epoch_mins;
  
  // total wake time (ignores pre lights out, post lights off)
  TWT = mins_wake;
  
  // final wake time 
  FWT = ( lights_on_epoch - final_wake_epoch ) * epoch_mins; 

  // REM latency
  rem_lat_mins = ( first_rem_epoch - first_sleep_epoch ) * epoch_mins;
  
  // Total sleep time (excludes 'other')
  TST = TRT - TWT - mins_other ; 

  // std::cout << "TST TRT " << TST << " " << TRT << "\n";
  // std::cout << "dets" << lights_out_epoch << " " << first_sleep_epoch << " " << last_sleep_epoch << " " << lights_on_epoch << "\n";
  
  // sleep latency
  slp_lat = ( first_sleep_epoch - lights_out_epoch ) * epoch_mins;
  
  // latency to persistent sleep
  per_slp_lat = ( first_persistent_sleep_epoch - lights_out_epoch ) * epoch_mins;

  // Sleep period time
  SPT = TRT - slp_lat;

  //  std::cout << "TWT , slp_lat , FWT " << TWT << " " << slp_lat << " " << FWT << "\n";

  // WASO (ignores leading and also trailing wake)
  //  WASO = TWT - slp_lat - FWT;
  // BUT, the above might count UNDEF/OTHER in the leading/trailing wake, and so remove too much
  // easier to just figure it out by iteration
  int w = 0;
  for (int e=first_sleep_epoch;e<=last_sleep_epoch;e++)
    if ( stages[e] == WAKE ) ++w;
  WASO = w * epoch_mins;

  // sleep efficiency (includes sleep latency as W) include OTHER in denom
  slp_eff_pct = ( TST / TRT ) * 100;
    
  // sleep maintainence (ignores initial sleep latency as W) includes OTHER in denom
  slp_main_pct = ( TST / SPT ) * 100;

  // alternate: sleep maintainence/efficiency 2 (denom is from initial sleep to final sleep)
  // i.e. ignores both leading and trailing W; includes OTHER in denom
  slp_eff2_pct = ( TST / ( epoch_mins * ( last_sleep_epoch - first_sleep_epoch + 1 ) ) ) * 100 ; 

  if ( TST > 0 ) 
    {
      pct_n1  = mins_n1  / TST;
      pct_n2  = mins_n2  / TST;
      pct_n3  = mins_n3  / TST;
      pct_n4  = mins_n4  / TST;
      pct_rem = mins_rem / TST;
    }
  else
    {
      pct_n1  = 0;
      pct_n2  = 0;
      pct_n3  = 0;
      pct_n4  = 0;
      pct_rem = 0;
      
    }

   


  //
  // Sleep cycles : based on modified Floyd & Feinberg rules
  //
  

  // Thresolds

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
  
  TpST = 0;

  std::vector<std::string> persistent_sleep( ne , "" );
  for (int e=0;e<ne;e++)
    {
      
      if ( stages[ e ] == WAKE || stages[ e ] == LIGHTS_ON || stages[e] == UNSCORED )
	{
	  persistent_sleep[e] = "W";
	  continue;
	}
      
      // otherwise, assume all other annotations are consistent with sleep
      bool okay = true;
      int ec = e - def_persistent_sleep_epochs;
      
      while ( okay )
	{
	  if ( ec < 0 ) { okay = false; break; }	  
	  if ( stages[ ec ] == WAKE || stages[ ec ] == LIGHTS_ON ) { okay = false; break; }  	  
	  if ( ++ec == e ) break;
	}

      if ( okay ) 
	{
	  persistent_sleep[e] = "S"; 
	  TpST += epoch_mins;
	}
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


  // for (int e=0;e<ne;e++)
  //   {
  //     std::cout << "C\t" 
  // 		<< e+1 << "\t"
  // 		<< stages[e] << "\t"
  // 		<< sleep_period[e] << "\t"
  // 		<< sleep_code[e] <<"\t"
  // 		<< wata[e] << "\t"
  // 		<< cycle_ending_waso[e] << "\t"
  // 		<< sleep_cycle_number[e] << "\n";
		
  //   }
  
  //
  // Get cycle/period statistics
  //

  // Count number of cyles
  
  num_nremc = 0;
  nremc_mean_duration = 0;

  std::map<int,int> cmin;
  std::map<int,int> cmax;
  std::map<int,int> counts_rem;
  std::map<int,int> counts_nrem;
  std::map<int,int> counts_other;

  for (int e=0;e<ne;e++)
    {
      const int & sn = sleep_cycle_number[e];
      if ( sn == 0 ) continue;
      if ( sn > num_nremc ) num_nremc = sn;
      if ( cmin.find( sn ) == cmin.end() ) cmin[ sn ] = cmax[sn] = e;
      cmax[sn] = e; // track max
      if ( is_rem( stages[e] ) ) counts_rem[sn]++;
      else if ( is_nrem( stages[e] ) ) counts_nrem[sn]++;
      else counts_other[sn]++;
    }


  
  std::map<int,int>::iterator ii = cmin.begin();
  while ( ii != cmin.end())
    {
      const int & sn = ii->first ;
      
      // total cycle duration
      double dur = cmax[ sn ] - ii->second + 1;
      double dur_mins = dur * epoch_mins ; 
      
      nremc_mean_duration += dur_mins;
 
      nremc_duration[ sn ] = ( counts_rem[sn] + counts_nrem[sn] + counts_other[sn] ) * epoch_mins ; 
      nremc_nrem_duration[ sn ] = counts_nrem[sn] * epoch_mins ; 
      nremc_rem_duration[ sn ] = counts_rem[sn]  * epoch_mins ; 
      
      nremc_start_epoch[ sn ] = ii->second + 1 ;  // output 1-based coding

      nremc_epoch_duration[ sn ] = counts_rem[sn] + counts_nrem[sn] + counts_other[sn] ;
            
      ++ii;
    }

  if ( num_nremc > 0 ) nremc_mean_duration /= (double)num_nremc;

  // cycle positions
  cycle_pos_relative.resize( ne , -1 );
  cycle_pos_absolute.resize( ne , -1 );
  for (int e=0; e<ne; e++)
    {
      const int & sn = sleep_cycle_number[e];
      if ( sn == 0 ) continue;      
      int cycle_start = cmin[sn];
      
      // position within each cycle.
      cycle_pos_absolute[e] = ( e - cycle_start ) * epoch_mins ; 
      cycle_pos_relative[e] = cycle_pos_absolute[e] / (double)nremc_duration[sn];
    }


  // after the fact, track epoch-level stats
  in_persistent_sleep.resize( ne , false );
  for (int e=0; e<ne; e++)
    if ( persistent_sleep[e] == "S" ) in_persistent_sleep[e] = true;


  
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
      
      if ( left_n  > 0 ) left_wgt /= (double)left_n;
      if ( right_n > 0 ) right_wgt /= (double)right_n;

      // simple average of left/right averages
      // if no data, wgt will be 0, which is fine
      n2_ascdesc[e] = ( left_wgt + right_wgt ) / 2.0;

    }
  

  
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
  nearest_wake.resize( ne , 0 );

  nrem2rem.resize( ne , 0 ); nrem2rem_total.resize( ne , 0 );
  nrem2wake.resize( ne , 0 ); nrem2wake_total.resize( ne , 0 );

  for (int e = 0 ; e < ne ; e++)
    {
      
      //
      // calculate the number of similar epochs 
      // (FLANKING_SIM)
      //
      
      int sim = 0;  
      
      for (int j=1;j<ne;j++)
	{
	  const int eleft  = e - j;
	  const int eright = e + j;
	  // too much
	  if ( eleft < 0 || eright >= ne ) { sim = j-1; break; }
	  if ( stages[eleft] != stages[e] || stages[eright] != stages[e] ) { sim = j-1; break; }	  
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
      nearest_wake[e] = nw;

      //
      // transitions FROM N2?
      //
            
      if ( stages[e] == NREM2 )
	{
	  
	  // n2 to rem
	  int ei = 1;
	  while ( 1 ) 
	    {
	      if ( e + ei == ne ) { ei=0; break; }
	      if ( stages[ e + ei ] == NREM2 ) { ++ei; continue; }
	      if ( stages[ e + ei ] == REM   ) break;
	      ei = 0; break;
	    }
	  nrem2rem[e] = ei;
	  
	  // n2 to wake
	  ei = 1;
	  while ( 1 )
	    {
	      if ( e + ei == ne ) { ei=0; break; }
	      if ( stages[ e + ei ] == NREM2 ) { ++ei; continue; }
	      if ( stages[ e + ei ] == WAKE   ) break;
	      ei = 0; break;
	    }
	  nrem2wake[e] = ei;
	}
    
    } // next epoch
  
  // now figure out the _total values 
  // i.e. move forward and copy largest number until we hit 0      
  int e_rem  = nrem2rem[0];
  int e_wake = nrem2wake[0]; 
  
  for (int e = 1 ; e < ne ; e++ )
    {
      if ( nrem2rem[e] == 0 ) e_rem = 0;
      else if ( nrem2rem[e] > e_rem ) e_rem = nrem2rem[e];
      nrem2rem_total[e] = e_rem;

      if ( nrem2wake[e] == 0 ) e_wake = 0;
      else if ( nrem2wake[e] > e_wake ) e_wake = nrem2wake[e];
      nrem2wake_total[e] = e_wake;
    }
    
  //
  // Clocktime-based measures
  //
  
  clocktime_t starttime( timeline->edf->header.starttime );
  if ( ! starttime.valid ) 
    {
      clock_lights_out.valid = clock_sleep_onset.valid 
	= clock_sleep_midpoint.valid = clock_wake_time.valid 
	= clock_lights_on.valid = false;
    }
  else
    {
      clock_lights_out     = starttime;
        
      double epoch_hrs = epoch_mins / 60.0;

      clock_sleep_onset    = starttime;
      clock_sleep_onset.advance( epoch_hrs * first_sleep_epoch );
      
      clock_wake_time      = starttime;
      clock_wake_time.advance( epoch_hrs * final_wake_epoch );
      
      clock_lights_on      = starttime;
      clock_lights_on.advance( epoch_hrs * ne );
    
      clock_sleep_midpoint.midpoint( clock_sleep_onset , clock_wake_time );      
      
    }



}

void hypnogram_t::output( const bool verbose )
{

  // currently, this routine is hard-coded to assume 30-second epochs,
  // so for now flag if this is not the case (we can fix downstream)

  if ( verbose )
    {
      if ( ! Helper::similar( timeline->epoch_length() , 30 , 0.001 ) ) 
	Helper::halt( "requires 30-second epochs to be set currently" );
    }

  
  //
  // Per individual level output (VERBOSE MODE ONLY)
  //

  if ( verbose )
    {
      writer.var( "T1_LIGHTS_OFF"     , "Lights out time [0,24)" );
      writer.var( "T2_SLEEP_ONSET"    , "Sleep onset time [0,24)" );
      writer.var( "T3_SLEEP_MIDPOINT" , "Sleep mid-point time [0,24)" );
      writer.var( "T4_FINAL_WAKE"     , "Final wake time [0,24)" );
      writer.var( "T5_LIGHTS_ON"      , "Lights on time [0,24)" );
      
      writer.var( "NREMC" , "Number of NREM cycles" );
      writer.var( "NREMC_MINS" , "Average NREM cycle duration (mins)" );
      
      writer.var( "TIB" , "Time in Bed (hours): EDF start --> EDF stop" );
      writer.var( "TRT" , "Total Recording Time (hours): LIGHTS_OUT --> LIGHTS_ON" );
      writer.var( "TST" , "Total Sleep Time (hours): SLEEP_ONSET --> FINAL_WAKE" );
      writer.var( "TPST" , "Total persistent Sleep Time (hours): PERSISTENT_SLEEP_ONSET --> FINAL_WAKE" );
      
      writer.var( "TWT" , "Total Wake Time (hours): all WAKE" );
      writer.var( "WASO" , "Wake After Sleep Onset (hours)" );
      
      writer.var( "SLP_LAT" , "Sleep latency" );
      writer.var( "PER_SLP_LAT" , "Persistent sleep latency" );
      
      writer.var( "SLP_EFF"      , "Sleep efficiency: LIGHTS_OUT --> LIGHTS_ON" );
      writer.var( "SLP_MAIN_EFF" , "Sleep maintenance efficiency: LIGHTS_OUT --> LIGHTS_ON" );
      writer.var( "SLP_EFF2"     , "Sleep efficiency: SLEEP_ONSET --> FINAL_WAKE" );

      writer.var( "REM_LAT" , "REM latency (from SLEEP_ONSET)" );
      
      writer.var( "PCT_N1" , "Proportion of sleep that is N1" );
      writer.var( "PCT_N2" , "Proportion of sleep that is N2" );
      writer.var( "PCT_N3" , "Proportion of sleep that is N3" );
      writer.var( "PCT_N4" , "Proportion of sleep that is N4" );
      writer.var( "PCT_REM" , "Proportion of sleep that is REM" );
      
      writer.var( "MINS_N1" , "Proportion of sleep that is N1" );
      writer.var( "MINS_N2" , "Proportion of sleep that is N2" );
      writer.var( "MINS_N3" , "Proportion of sleep that is N3" );
      writer.var( "MINS_N4" , "Proportion of sleep that is N4" );
      writer.var( "MINS_REM" , "Proportion of sleep that is REM" );
      
      // values
      writer.value(  "T1_LIGHTS_OFF" , clock_lights_out.as_numeric_string() );

      if ( any_sleep ) 
	{
	  writer.value(  "T2_SLEEP_ONSET" , clock_sleep_onset.as_numeric_string() );
	  writer.value(  "T3_SLEEP_MIDPOINT" , clock_sleep_midpoint.as_numeric_string() );
	  writer.value(  "T4_FINAL_WAKE" , clock_wake_time.as_numeric_string() );
	}
      writer.value(  "T5_LIGHTS_ON" , clock_lights_on.as_numeric_string() );
      
      if ( any_sleep )
	{
	  writer.value(  "NREMC" , num_nremc );
	  writer.value(  "NREMC_MINS" , nremc_mean_duration );
	}

      writer.value( "TIB" , TIB );
      writer.value( "TRT" , TRT );
      writer.value( "TST" , TST );
      writer.value( "TPST" , TpST );
      writer.value( "TWT" , TWT );
      
      if ( any_sleep )
	{
	  writer.value( "WASO" , WASO );
	  writer.value( "SLP_LAT" , slp_lat );
	  writer.value( "PER_SLP_LAT" , per_slp_lat );      
	  writer.value( "SLP_EFF" , slp_eff_pct );
	  writer.value( "SLP_MAIN_EFF" , slp_main_pct );
	  writer.value( "SLP_EFF2" , slp_eff2_pct );
	}

      if ( mins_rem > 0 )
	writer.value( "REM_LAT" , rem_lat_mins );

      if ( any_sleep )
	{
	  writer.value( "PCT_N1" , pct_n1 );
	  writer.value( "PCT_N2" , pct_n2 );
	  writer.value( "PCT_N3" , pct_n3 );
	  writer.value( "PCT_N4" , pct_n4 );
	  writer.value( "PCT_REM" , pct_rem);
	}

      writer.value( "MINS_N1" , mins_n1 );
      writer.value( "MINS_N2" , mins_n2 );
      writer.value( "MINS_N3" , mins_n3 );
      writer.value( "MINS_N4" , mins_n4 );
      writer.value( "MINS_REM" , mins_rem);

    }
  
  //
  // Cycle-specific output (verbose mode only)
  //

  if ( verbose ) 
    {

      if ( any_sleep ) 
	{
	  writer.var( "NREMC_START" , "NREM cycle start epoch" );
	  writer.var( "NREMC_NREM_MINS" , "NREM cycle NREM duration (mins)" );
	  writer.var( "NREMC_REM_MINS" , "NREM cycle REM duration (mins)" );
	  writer.var( "NREMC_OTHER_MINS" , "NREM cycle other duration (mins)" );
	  writer.var( "NREMC_MINS" , "NREM cycle total duration (mins)" );
	  writer.var( "NREMC_N" , "NREM cycle total duration (epochs)" );
	  
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
  stagen[ NREM4 ] = -4;
  stagen[ UNSCORED ] = 2;
  stagen[ UNKNOWN ] = 2; // this should not happen
  stagen[ MOVEMENT ] = 2;
  stagen[ ARTIFACT ] = 2;
  stagen[ LIGHTS_ON ] = 2;

  writer.var( "MINS" , "Elapsed time since start of recording (minutes)" );
  writer.var( "CLOCK_TIME" , "Clock time (hh:mm:ss)" );

  writer.var( "STAGE" , "Sleep stage, string label" );
  writer.var( "STAGE_N" , "Sleep stage, numeric encoding" );

  // epoch size (in minutes)
  const double epoch_mins = timeline->epoch_length() / 60.0 ; 
  const double epoch_hrs = epoch_mins / 60.0;
  
  const int ne = timeline->num_epochs();
  
  clocktime_t starttime( clock_lights_out );


  
  // output in non-verbsoe mode (STAGES command)

  if ( ! verbose )
    {
      
      for (int e=0;e<ne;e++)
	{
	  
	  // epoch-level stratification
	  // epoch_n is the original 0-based epoch encoding
	  // so, for diplay +1 , but for other calculations
	  // we want to keep this original encoding
	  
	  writer.epoch( epoch_n[e] + 1 );
	  
	  
	  // clock time based on EDF header
	  
	  if ( starttime.valid ) 
	    {
	      clocktime_t current_clock_time = starttime;
	      
	      current_clock_time.advance( epoch_hrs * epoch_n[e] );
	      
	      writer.value( "CLOCK_TIME" , current_clock_time.as_string() );      
	      
	      if ( verbose ) 
		writer.value( "CLOCK_HOURS" ,  current_clock_time.as_numeric_string() );
	      
	    }
	  
	  // time in minutes
	  
	  writer.value( "MINS" ,  epoch_n[e] * epoch_mins );
	  
	  // stages
	  
	  writer.value( "STAGE" , globals::stage( stages[e] ) );
	  
	  writer.value( "STAGE_N" , stagen[ stages[e] ] );
	  
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
  // Add cycle epoch-annotation
  //
  
  for (int e=0;e<ne;e++)
    if ( sleep_cycle_number[e] ) 
      {	
	const std::string cycle = "_NREMC_" + Helper::int2str( sleep_cycle_number[e] );	
	timeline->annotate_epoch( cycle , e );	
      }


  double elapsed_n1 = 0 , elapsed_n2 = 0 , elapsed_n34 = 0 , elapsed_rem = 0;
    
  double elapsed_sleep = 0 , elapsed_wake = 0 , elapsed_waso = 0 ;

  
  // header

  writer.var( "CLOCK_HOURS" , "Clock time [0,24) hours" );
        
  writer.var( "E_WAKE" , "Elapsed wake (mins)" );
  writer.var( "E_WASO" , "Elapsed WASO (mins)" );
  writer.var( "E_SLEEP" , "Elapsed sleep (mins)" );

  writer.var( "E_N1" , "Elapsed N1 (mins)" );
  writer.var( "E_N2" , "Elapsed N2 (mins)" );
  writer.var( "E_N3" , "Elapsed N3 (mins)" );
  writer.var( "E_REM" , "Elapsed REM (mins)" );

  writer.var( "PCT_E_SLEEP" , "Elapsed sleep (percent of all sleep)" );
  writer.var( "PCT_E_N1" , "Elapsed N1 (percent of all N1)" );
  writer.var( "PCT_E_N2" , "Elapsed N2 (percent of all N2)" );
  writer.var( "PCT_E_N3" , "Elapsed N3 (percent of all N3)" );
  writer.var( "PCT_E_REM" , "Elapsed REM (percent of all REM)" );

  writer.var( "PERSISTENT_SLEEP" , "Persistent sleep yes/no? (1=Y)" );

  writer.var( "CYCLE" , "NREMC number" );
  writer.var( "PERIOD" , "NREMC period (NREM/REM)" );
  
  writer.var( "CYCLE_POS_REL" , "Position within NREMC, relative" );
  writer.var( "CYCLE_POS_ABS" , "Position within NREMC, absolute (mins)" );

  writer.var( "FLANKING_SIM" , "Number of similar epochs w.r.t. stage" );
  writer.var( "NEAREST_WAKE" , "Number of epochs until the nearest wake" );

  writer.var( "WASO" , "Epoch is WASO (1=Y)" );
    
  
  // these next four are all reported for the  NREM epoch
  writer.var( "NREM2REM" , "If NREM epoch, number of NREM if next non-NREM is REM" );
  writer.var( "NREM2REM_TOTAL" , "If NREM epoch, total number of contiguous NREM if next non-NREM is REM" );

  writer.var( "NREM2WAKE" , "If NREM epoch, number of NREM if next non-NREM is WAKE" );
  writer.var( "NREM2WAKE_TOTAL" , "If NREM epoch, total number of contiguous NREM if next non-NREM is WAKE" );

  writer.var ("N2_WGT" , "Score for descending/ascending N2 epochs (-1 to +1)" );

  
  
  // output
  for (int e=0;e<ne;e++)
    {
      
      // epoch-level stratification

      // nb. we can use display_epoch() here, not epoch_n[], as 
      // HYPNO (verbose=T) should never be called on a masked/noncontiguous 
      // set of epochs...

      writer.epoch( timeline->display_epoch( e ) );
      

      // clock time based on EDF header
      if ( starttime.valid ) 
	{
	  clocktime_t current_clock_time = starttime;
	  
	  current_clock_time.advance( epoch_hrs * e );
	  
	  writer.value( "CLOCK_TIME" , current_clock_time.as_string() );      
	  
	  if ( verbose ) 
	    writer.value( "CLOCK_HOURS" ,  current_clock_time.as_numeric_string() );

	}

      // time in minutes
      writer.value( "MINS" ,  e * epoch_mins );
      
      // stages      
      writer.value( "STAGE" , globals::stage( stages[e] ) );    
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

      writer.value( "PCT_E_N1" , mins_n1 > 0 ? elapsed_n1 / mins_n1 : 0 );
      writer.value( "PCT_E_N2" , mins_n2 > 0 ? elapsed_n2 / mins_n2 : 0 );
      writer.value( "PCT_E_N3" , (mins_n3+mins_n4) > 0 ? elapsed_n34 / (mins_n3+mins_n4) : 0 );
      writer.value( "PCT_E_REM" , mins_rem > 0 ? elapsed_rem / mins_rem : 0 );

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
      
      // cycles
      
      if ( sleep_cycle_number[e] )
	{
	  writer.value( "CYCLE" , sleep_cycle_number[e] );
	  writer.value( "PERIOD" , sleep_code[e] == 5 ? "REMP" : sleep_code[e]==1 ? "NREMP" : "." ) ;
	  writer.value( "CYCLE_POS_REL" , cycle_pos_relative[e] );
	  writer.value( "CYCLE_POS_ABS" , cycle_pos_absolute[e] );
	}
      
      // flanking epochs

      writer.value( "FLANKING_SIM" , flanking[e] );
      writer.value( "NEAREST_WAKE" , nearest_wake[e] );
      writer.value( "WASO" , is_waso[e] );

      writer.value( "NREM2REM" , nrem2rem[e] );
      writer.value( "NREM2REM_TOTAL" , nrem2rem_total[e] );

      writer.value( "NREM2WAKE" , nrem2wake[e] );
      writer.value( "NREM2WAKE_TOTAL" , nrem2wake_total[e] );

      // N2 ascending/descending status

      if ( stages[e] == NREM2 ) 
	writer.value( "N2_WGT" , n2_ascdesc[e] );
           
    } // next epoch

  writer.unepoch();


}

void dummy_hypno()
{
  edf_t edf;
  
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
      else if ( s == "?"  ) h.stages.push_back( UNSCORED );
      else logger << "did not recognize " << s << "\n";
    }

  logger << "read " << h.stages.size() << "\n";

  edf.header.starttime = "10:00:00";
  
  // fudge so it works in this non-standard case...
  edf.id = "_DUMMY_";
  h.fudge( 30 , h.stages.size() );

  h.calc_stats( true );
  h.output( true ); // verbose mode == T 

}

void hypnogram_t::fudge( double es, int ne )
{
  timeline->epoch_length_tp = es * globals::tp_1sec;
  timeline->epochs.resize( ne );
}




