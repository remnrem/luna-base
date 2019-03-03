
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
#include "helper/logger.h"
#include "annot/annot.h"
#include "annot/nsrr-remap.h"

#include <iostream>

extern writer_t writer;

extern logger_t logger;

std::string globals::version;
std::string globals::date;

std::string globals::annot_folder;
std::vector<std::string> globals::annot_files;

bool globals::read_ftr;
std::set<std::string> globals::specified_annots;
bool globals::remap_nsrr_annots;

std::map<globals::atype_t,std::string> globals::type_name;
std::map<std::string,globals::atype_t> globals::name_type;

bool globals::enforce_epoch_check;
int globals::default_epoch_len;

std::map<frequency_band_t,freq_range_t> globals::freq_band;

sleep_stage_label_t globals::sleep_stage;
sleep_stage_label_lookup_t globals::sleep_stage_labels;

char globals::folder_delimiter;
std::string globals::project_path;

bool globals::assume_pm_starttime;

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
std::string globals::annot_instance_strat;
std::string globals::annot_meta_strat;
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
  

  version = "v0.9";
  date    = "12-Feb-2019";

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
  // Annotation folder
  //
  
  annot_folder = "";

  annot_files.clear();
  
  
  //
  // Requested to load specific annotations only?
  //

  specified_annots.clear();
  
  //
  // Automatically remap NSRR annotations (nsrr-remap=Y)
  //

  remap_nsrr_annots = true;

  nsrr_t::init();

  //
  // By default, read and extract all FTR; if this is a pain, can be turned off (ftr=0)
  //
  
  read_ftr = true;
  
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
  // Primary sleep stage encoding 
  //
  
  sleep_stage[ WAKE  ]     = "wake";
  sleep_stage[ LIGHTS_ON ] = "L";
  sleep_stage[ NREM1 ]     = "NREM1";
  sleep_stage[ NREM2 ]     = "NREM2";
  sleep_stage[ NREM3 ]     = "NREM3";
  sleep_stage[ NREM4 ]     = "NREM4";
  sleep_stage[ REM   ]     = "REM";
  sleep_stage[ MOVEMENT ]  = "M";
  sleep_stage[ UNSCORED ]  = "?";
  sleep_stage[ UNKNOWN ]   = ".";


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
  sleep_stage_labels[ "Unsure|Unsure" ]   = UNSCORED;

  // Basic
  sleep_stage_labels[ "wake" ]     = WAKE;  
  sleep_stage_labels[ "NREM1" ]    = NREM1;  
  sleep_stage_labels[ "NREM2" ]    = NREM2;  
  sleep_stage_labels[ "NREM3" ]    = NREM3;  
  sleep_stage_labels[ "NREM4" ]    = NREM4;  
  sleep_stage_labels[ "REM" ]      = REM;
  sleep_stage_labels[ "Movement" ] = MOVEMENT;
  sleep_stage_labels[ "Unscored" ] = UNSCORED;
  sleep_stage_labels[ "L" ] = LIGHTS_ON;

  // minimal 
  sleep_stage_labels[ "W" ]     = WAKE;  
  sleep_stage_labels[ "N1" ]    = NREM1;  
  sleep_stage_labels[ "N2" ]    = NREM2;  
  sleep_stage_labels[ "N3" ]    = NREM3;  
  sleep_stage_labels[ "N4" ]    = NREM4;  
  sleep_stage_labels[ "R" ]     = REM;
  sleep_stage_labels[ "?" ]     = UNSCORED;
  sleep_stage_labels[ "M" ]     = MOVEMENT;  
  
  // mouse: make generic 'NR' -> NREM2
  sleep_stage_labels[ "W" ]     = WAKE;  
  sleep_stage_labels[ "NR" ]    = NREM2;  
  sleep_stage_labels[ "R" ]     = REM;
  sleep_stage_labels[ "?" ]     = UNSCORED;
  
  //
  // Time-units
  //
  
  // 1e-12 sec resolution
//   tp_1sec  = 1000000000000LLU; 
//   tp_duration    = 1.0 / (double)tp_1sec;

  // 1e-9 sec resolution
  tp_1sec  = 1000000000LLU;
  tp_duration    = 1.0 / (double)tp_1sec;
  
  //
  // Common output stratifiers
  //

  freq_strat   = "F";
  signal_strat = "CH";
  stage_strat  = "S";
  cycle_strat  = "C";
  band_strat   = "B";
  annot_strat  = "ANNOT";
  annot_instance_strat  = "INST";
  annot_meta_strat  = "META";
  count_strat  = "N";
  epoch_strat  = "E";
  time_strat   = "T";
  sample_strat = "SP"; // sample-point

  //
  // Misc.
  //

  project_path = "";
  
  folder_delimiter = '/';
  
  skip_edf_annots = false;

  current_tag = "";

  indiv_wildcard = "^";

  sample_list_min = -1;
  sample_list_max = -1;
  sample_list_id = "";

  assume_pm_starttime = true;
  
  problem = false;
  
  edf_timetrack_label = "_TT";
  edf_timetrack_size = 15; // i.e. up to 30 chars

  
  // whether to assume 30-sec and enfore epoch check when first attaching 
  // an .eannot file

  enforce_epoch_check = true;

  default_epoch_len = 30;
    
  //
  // Annot types
  //

  type_name[ A_NULL_T ] = "null";
  type_name[ A_FLAG_T ] = "flag";
  type_name[ A_MASK_T ] = "mask";

  type_name[ A_TXT_T ] = "txt";
  type_name[ A_INT_T ] = "int";
  type_name[ A_DBL_T ] = "num";
  type_name[ A_BOOL_T ] = "bool";

  type_name[ A_TXTVEC_T ] = "txtvec";
  type_name[ A_INTVEC_T ] = "intvec";
  type_name[ A_DBLVEC_T ] = "numvec";
  type_name[ A_BOOLVEC_T ] = "boolvec";
  
  // flags (i.e. no value)

  name_type[ "FLAG" ] = A_FLAG_T; 
  name_type[ "flag" ] = A_FLAG_T; 

  name_type[ "MASK" ] = A_MASK_T;
  name_type[ "mask" ] = A_MASK_T;
  
  // scalars 

  name_type[ "TXT" ] = A_TXT_T;
  name_type[ "txt" ] = A_TXT_T;

  name_type[ "INT" ] = A_INT_T; 
  name_type[ "int" ] = A_INT_T; 

  name_type[ "NUM" ] = A_DBL_T; 
  name_type[ "num" ] = A_DBL_T; 

  name_type[ "BOOL" ] = A_BOOL_T; 
  name_type[ "bool" ] = A_BOOL_T; 
  name_type[ "YN" ] = A_BOOL_T; 
  name_type[ "yn" ] = A_BOOL_T; 


  // vectors

  name_type[ "TXTVEC" ] = A_TXTVEC_T;
  name_type[ "txtvec" ] = A_TXTVEC_T;

  name_type[ "INTVEC" ] = A_INTVEC_T; 
  name_type[ "intvec" ] = A_INTVEC_T; 

  name_type[ "NUMVEC" ] = A_DBLVEC_T; 
  name_type[ "numvec" ] = A_DBLVEC_T; 

  name_type[ "BOOLVEC" ] = A_BOOLVEC_T; 
  name_type[ "boolvec" ] = A_BOOLVEC_T; 
  name_type[ "YNVEC" ] = A_BOOLVEC_T; 
  name_type[ "ynvec" ] = A_BOOLVEC_T; 

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
  else if ( s == 6 ) return stage( UNSCORED );
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
