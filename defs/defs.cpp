
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

#include "defs.h"
#include "miscmath/crandom.h"
#include "eval.h"
#include "db/db.h"

#include <iostream>

extern writer_t writer;

std::string globals::version;
std::string globals::date;

std::string globals::annot_folder;

std::map<frequency_band_t,freq_range_t> globals::freq_band;

sleep_stage_label_t globals::sleep_stage;
sleep_stage_label_lookup_t globals::sleep_stage_labels;

char globals::folder_delimiter;
std::string globals::current_tag;
std::string globals::indiv_wildcard;
bool globals::skip_edf_annots;

std::set<std::string> globals::excludes;

int globals::sample_list_min;
int globals::sample_list_max;
std::string globals::sample_list_id;

std::string globals::edf_timetrack_label;
int globals::edf_timetrack_size;

uint64_t globals::tp_1sec;
double globals::tp_duration;

bool globals::problem;

std::map<std::string,sleep_stage_t> sleep_stage_labels;

param_t globals::param;

void (*globals::bail_function) ( const std::string & );

bool globals::silent; 
bool globals::Rmode;

std::string globals::epoch_strat;
std::string globals::time_strat;
std::string globals::freq_strat;
std::string globals::signal_strat;
std::string globals::stage_strat;
std::string globals::cycle_strat;
std::string globals::band_strat;
std::string globals::annot_strat;
std::string globals::count_strat;
std::string globals::sample_strat;

std::string & globals::SQLITE_SCRATCH_FOLDER() { static std::string s = ""; return s; }

void globals::api()
{
  silent = true;
  writer.nodb();
}

void globals::R()
{
  Rmode = true;
  api();
}

void globals::init_defs()
{
  
  //
  // Version
  //
  
  version = "v0.2";
  date    = "12-Dec-2018";
  
  //
  // Set up RNG
  //

  CRandom::srand(time(0));

  //
  // OPtional bail function after halt() is called
  //
  
  bail_function = NULL;

  silent = false;
  
  //
  // Annotations
  //
  
  annot_folder = "";
  
  //
  // Frequency bands req. bins as defined Manoach et al. (2014) (Table 4)
  //
  
  freq_band[ SLOW  ] = freq_range_t( 0.5 , 1    );

  freq_band[ DELTA ] = freq_range_t( 1   , 4    );
  freq_band[ THETA ] = freq_range_t( 4   , 8    );
  freq_band[ ALPHA ] = freq_range_t( 8   , 12   );

  freq_band[ SIGMA      ] = freq_range_t( 12  , 15   );
  freq_band[ LOW_SIGMA  ] = freq_range_t( 12   ,  13.5 );
  freq_band[ HIGH_SIGMA ] = freq_range_t( 13.5 ,  15   );

  freq_band[ BETA  ] = freq_range_t( 15  , 30   );
  freq_band[ GAMMA ] = freq_range_t( 30  , 1000 );

  freq_band[ TOTAL ] = freq_range_t( 0  , 1000 );
  

  //
  // Sleep stage encoding (as per SHHS)
  //

  sleep_stage[ WAKE  ]    = "Wake";
  sleep_stage[ LIGHTS_ON ] = "LightsOn";
  sleep_stage[ NREM1 ]    = "NREM1";
  sleep_stage[ NREM2 ]    = "NREM2";
  sleep_stage[ NREM3 ]    = "NREM3";
  sleep_stage[ NREM4 ]    = "NREM4";
  sleep_stage[ REM   ]    = "REM";
  sleep_stage[ MOVEMENT ] = "Movement";
  sleep_stage[ UNKNOWN ]  = "Unknown";


  //
  // Common/NSRR labels
  //

  // e.g. SOF study
  sleep_stage_labels[ "SRO:Wake" ]             = WAKE;
  sleep_stage_labels[ "SRO:Stage1Sleep" ]      = NREM1;
  sleep_stage_labels[ "SRO:Stage2Sleep" ]      = NREM2;
  sleep_stage_labels[ "SRO:Stage3Sleep" ]      = NREM3;
  sleep_stage_labels[ "SRO:Stage4Sleep" ]      = NREM4;
  sleep_stage_labels[ "SRO:Stage34Sleep" ]     = NREM3; // if this exists?
  sleep_stage_labels[ "SRO:RapidEyeMovement" ] = REM;

  // e.g. SHHS
  sleep_stage_labels[ "SDO:WakeState" ]                   = WAKE;
  sleep_stage_labels[ "SDO:NonRapidEyeMovementSleep-N1" ] = NREM1;
  sleep_stage_labels[ "SDO:NonRapidEyeMovementSleep-N2" ] = NREM2;
  sleep_stage_labels[ "SDO:NonRapidEyeMovementSleep-N3" ] = NREM3;
  sleep_stage_labels[ "SDO:NonRapidEyeMovementSleep-N4" ] = NREM4;
  sleep_stage_labels[ "SDO:RapidEyeMovementSleep" ]       = REM;

  // other NSRR
  sleep_stage_labels[ "Wake|0" ]          = WAKE;
  sleep_stage_labels[ "Stage 1 sleep|1" ] = NREM1;
  sleep_stage_labels[ "Stage 2 sleep|2" ] = NREM2;
  sleep_stage_labels[ "Stage 3 sleep|3" ] = NREM3;
  sleep_stage_labels[ "Stage 4 sleep|4" ] = NREM4;
  sleep_stage_labels[ "REM sleep|5" ]     = REM;
  sleep_stage_labels[ "Unsure|Unsure" ]   = UNKNOWN;

  // Basic

  sleep_stage_labels[ "Wake" ]     = WAKE;  
  sleep_stage_labels[ "NREM1" ]    = NREM1;  
  sleep_stage_labels[ "NREM2" ]    = NREM2;  
  sleep_stage_labels[ "NREM3" ]    = NREM3;  
  sleep_stage_labels[ "NREM4" ]    = NREM4;  
  sleep_stage_labels[ "REM" ]      = REM;
  sleep_stage_labels[ "Movement" ] = MOVEMENT;
  sleep_stage_labels[ "Unscored" ] = UNKNOWN; 
  sleep_stage_labels[ "L" ] = LIGHTS_ON;

  // minimal 

  sleep_stage_labels[ "W" ]     = WAKE;  
  sleep_stage_labels[ "N1" ]    = NREM1;  
  sleep_stage_labels[ "N2" ]    = NREM2;  
  sleep_stage_labels[ "N3" ]    = NREM3;  
  sleep_stage_labels[ "N4" ]    = NREM4;  
  sleep_stage_labels[ "R" ]     = REM;
  sleep_stage_labels[ "?" ]     = UNKNOWN;  
  sleep_stage_labels[ "M" ]     = MOVEMENT;  
  
  // mouse: make generic 'NR' -> NREM2
  sleep_stage_labels[ "W" ]     = WAKE;  
  sleep_stage_labels[ "NR" ]    = NREM2;  
  sleep_stage_labels[ "R" ]     = REM;
  sleep_stage_labels[ "?" ]     = UNKNOWN;  
  
  //
  // Time-units
  //
  
  tp_1sec  = 1000000000000LLU; // 1e-12 sec resolution
  tp_duration    = 1.0 / (double)tp_1sec;
  
  
  //
  // Common output stratifiers
  //

  freq_strat   = "F";
  signal_strat = "CH";
  stage_strat  = "S";
  cycle_strat  = "C";
  band_strat   = "B";
  annot_strat  = "ANN";
  count_strat  = "N";
  epoch_strat  = "E";
  time_strat   = "T";
  sample_strat = "SP"; // sample-point

  //
  // Misc.
  //

  folder_delimiter = '/';
  
  skip_edf_annots = false;

  current_tag = "";

  indiv_wildcard = "^";

  sample_list_min = -1;
  sample_list_max = -1;
  sample_list_id = "";

  problem = false;
  
  edf_timetrack_label = "_TT";
  edf_timetrack_size = 15; // i.e. up to 30 chars



}

std::string globals::band( frequency_band_t b )
{
  switch( b )
    {
    case SLOW  : return "SLOW";
    case ALPHA : return "ALPHA";
    case BETA  : return "BETA";
    case GAMMA : return "GAMMA";
    case THETA : return "THETA";
    case DELTA : return "DELTA";
    case SIGMA : return "SIGMA";
    case HIGH_SIGMA : return "FAST_SIGMA";
    case LOW_SIGMA  : return "SLOW_SIGMA";    
    case TOTAL : return "TOTAL";
    default : return "UNKNOWN";
    } 
}


std::string globals::stage( sleep_stage_t s )
{
  if ( sleep_stage.find(s) == sleep_stage.end() ) return "?";
  return sleep_stage[s];
}

std::string globals::stage( int s )
{
  if      ( s == 0 ) return stage( WAKE );
  else if ( s == 1 ) return stage( NREM1 );
  else if ( s == 2 ) return stage( NREM2 );
  else if ( s == 3 ) return stage( NREM3 );
  else if ( s == 4 ) return stage( NREM4 );
  else if ( s == 5 ) return stage( REM );
  else return stage( UNKNOWN );
}


sleep_stage_t globals::stage( const std::string & s )
{
  sleep_stage_label_lookup_t::const_iterator ii = sleep_stage_labels.find( s );
  if ( ii == sleep_stage_labels.end() ) return UNKNOWN;
  return ii->second;
}

double globals::band_width( frequency_band_t b )
{
  if ( freq_band.find(b) == freq_band.end() ) return 0;
  freq_range_t r = freq_band[b];
  return r.second - r.first;
}

std::string globals::print( const freq_range_t & r )
{
  std::stringstream ss;
  ss << r.first << ".." << r.second ;
  return ss.str();
}
