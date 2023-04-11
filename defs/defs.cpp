
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
#include "edf/signal-list.h"

#include <iostream>

extern writer_t writer;

extern logger_t logger;

std::string globals::version;
std::string globals::date;

int globals::retcode;

cmddefs_t & globals::cmddefs()
{
  static cmddefs_t * ans = new cmddefs_t();
  return *ans;
}

//std::string globals::annot_folder;
std::vector<std::string> globals::annot_files;
bool globals::allow_space_delim = false;
char globals::annot_class_inst_combiner = '_';
bool globals::combine_annot_class_inst = false;
char globals::class_inst_delimiter = ':';
char globals::annot_keyval_delim = '=' ; 
std::string globals::annot_disc_segment = "segment";
std::string globals::annot_disc_gap = "gap";
bool globals::annot_disc_drop_spanning = true;

bool globals::sanitize_everything = true;

std::set<std::string> globals::annot_alignment;
bool globals::autofix_edf;

int globals::time_format_dp;

bool globals::read_ftr;
std::set<std::string> globals::specified_annots;

std::map<globals::atype_t,std::string> globals::type_name;
std::map<std::string,globals::atype_t> globals::name_type;

channel_map_t globals::chmap1; 
channel_map_t globals::chmap2; 
std::map<std::string,channel_type_t> globals::label2ch;
std::map<channel_type_t,std::string> globals::ch2label;
std::map<std::string,channel_type_t> globals::sig2type; 

int globals::enforce_epoch_check;
int globals::default_epoch_len;

std::map<frequency_band_t,freq_range_t> globals::freq_band;

std::string globals::sleep_stage_prefix;
sleep_stage_label_t globals::sleep_stage;
sleep_stage_label_lookup_t globals::sleep_stage_labels;
bool globals::sleep_stage_assume_epoch_duration;

bool globals::replace_channel_spaces;
bool globals::uppercase_channels;
bool globals::replace_annot_spaces;
char globals::space_replacement;
bool globals::set_0dur_as_ellipsis;

char globals::folder_delimiter;
std::string globals::project_path;
char globals::file_list_delimiter;

std::string globals::mkdir_command;

std::string globals::txt_table_prepend;
std::string globals::txt_table_append;

bool globals::assume_pm_starttime;
int globals::assume_pm_starttime_hour;

std::string globals::current_tag;
std::string globals::indiv_wildcard;

int globals::anon_idroot_cnt;

bool globals::force_edf;
bool globals::skip_edf_annots;
bool globals::skip_nonedf_annots;
bool globals::set_annot_inst2hms;
bool globals::set_annot_inst2hms_force;
bool globals::skip_sl_annots;

std::set<std::string> globals::id_excludes;
std::set<std::string> globals::id_includes;

std::set<std::string> globals::sl_annot_extensions;
std::map<std::string,sample_list_t> globals::sl_data;
bool globals::sl_visit_edf;
bool globals::sl_link_across_folders;

int globals::sample_list_min;
int globals::sample_list_max;
std::string globals::sample_list_id;

bool globals::anon; 
std::string globals::force_starttime;
std::string globals::force_startdate;

bool globals::write_naughty_list;
std::string globals::naughty_list;

std::string globals::edf_timetrack_label;
int globals::edf_timetrack_size;
std::string globals::edf_annot_label;

uint64_t globals::tp_1sec;
uint64_t globals::tp_1000thsec;
double globals::tp_duration;

bool globals::problem;
bool globals::empty;

bool globals::bail_on_fail;

std::map<std::string,sleep_stage_t> sleep_stage_labels;

param_t globals::param;

void (*globals::bail_function) ( const std::string & );

void (*globals::logger_function) ( const std::string & );

bool globals::silent;
bool globals::verbose; 
bool globals::Rmode;
bool globals::Rdisp;
bool globals::devel;

std::string globals::epoch_strat;
std::string globals::time_strat;
std::string globals::freq_strat;
std::string globals::sec_strat;
std::string globals::signal_strat;
std::string globals::signal1_strat;
std::string globals::signal2_strat;
std::string globals::stage_strat;
std::string globals::cycle_strat;
std::string globals::band_strat;
std::string globals::annot_strat;
std::string globals::annot_instance_strat;
std::string globals::annot_meta_strat;
std::string globals::count_strat;
std::string globals::sample_strat;
std::string globals::anchor_strat;
std::string globals::cluster_strat;
std::string globals::var_strat;
std::string globals::value_strat;
std::string globals::feature_strat;

std::string & globals::SQLITE_SCRATCH_FOLDER() { static std::string s = ""; return s; }

void globals::api()
{
  silent = true;
  writer.nodb();
}

void globals::R( bool disp )
{
  Rmode = true;
  Rdisp = disp;
  api();
}


void globals::init_defs()
{

  
  //
  // Version
  //
  
  version = "v0.28.0";
  
  date    = "10-Apr-2023";

  //
  // Return code
  //

  retcode = 0;

  //
  // initialize cmddefs_t (note, using pattern to avoid static initialization issues)
  //
  
  static cmddefs_t & x = cmddefs();
  x.init();
  
  
  //
  // Set up RNG
  //

  CRandom::srand(time(0));

  //
  // Optional bail function after halt() is called
  //
  
  bail_function = NULL;

  //
  // Optional redirect of logger?
  //

  logger_function = NULL; 
  
  //
  // Output
  //

  silent = false;
  
  Rmode = false;

  Rdisp = false;

  verbose = false; 

  devel = false;
		 
  time_format_dp = 3;

  //
  // Annotation folder
  //
  
  //  annot_folder = "";

  annot_files.clear();
  

  //
  // Spaces in channel names
  //

  replace_channel_spaces = true;

  uppercase_channels = false;
  
  replace_annot_spaces = true;
   
  space_replacement = '_';

  //
  // On writing .annot (only), set 0 duration events to '...' 
  //

  set_0dur_as_ellipsis = false;

  
  //
  // Which annots, if any, should be aligned to the starts of EDF records?
  //
    
  annot_alignment.clear();

  
  //
  // Requested to load specific annotations only?
  //

  specified_annots.clear();

  //
  // try to fix EDF if 'corrupt'
  //

  autofix_edf = false;

  
  //
  // Automatically remap stage annotations; 
  // nb. if annot-remap=F then this is subsequently wiped
  //     if nssr-remap=T then we also later call nsrr_t::init_nsrr_mappings()
  //
  
  nsrr_t::init();

  
  //
  // --build
  //

  sl_visit_edf = true;
  
  
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
  freq_band[ GAMMA ] = freq_range_t( 30  , 50 );

  freq_band[ TOTAL ] = freq_range_t( 0.5  , 50 );
  

  //
  // Primary sleep stage encoding 
  //

  // i.e. to track different schemes, SUDS_N1,  sleep_stage_prefix == "SUDS"
  // SOAP 's'
  // POPS 'p' (predicted)
  
  sleep_stage_prefix       = ""; // by default not  

  // if we find 0-duration annots, assume epoch length
  sleep_stage_assume_epoch_duration = true;
  
  sleep_stage[ WAKE  ]     = "W";
  sleep_stage[ LIGHTS_ON ] = "L";
  sleep_stage[ NREM1 ]     = "N1";
  sleep_stage[ NREM2 ]     = "N2";
  sleep_stage[ NREM3 ]     = "N3";
  sleep_stage[ NREM4 ]     = "NREM4"; // redundant
  sleep_stage[ REM   ]     = "R";
  sleep_stage[ MOVEMENT ]  = "M";
  sleep_stage[ UNSCORED ]  = "U"; // ambiguous - 'unscorable'
  sleep_stage[ UNKNOWN ]   = "?"; // missing
  sleep_stage[ GAP ]       = "GAP"; // GAP

  // minimal: default
  sleep_stage_labels[ "W" ]     = WAKE;  
  sleep_stage_labels[ "N1" ]    = NREM1;  
  sleep_stage_labels[ "N2" ]    = NREM2;  
  sleep_stage_labels[ "N3" ]    = NREM3;  
  sleep_stage_labels[ "N4" ]    = NREM4; // should not encounter
  sleep_stage_labels[ "R" ]     = REM;
  sleep_stage_labels[ "U" ]     = UNSCORED;
  sleep_stage_labels[ "?" ]     = UNKNOWN;
  sleep_stage_labels[ "M" ]     = MOVEMENT;  
  sleep_stage_labels[ "L" ]     = LIGHTS_ON;
  sleep_stage_labels[ "G" ]     = GAP ; 

  
  //
  // POPS predictions
  //
  
  // sleep_stage_labels[ "pW" ]  = WAKE;
  // sleep_stage_labels[ "pN1" ] = NREM1;
  // sleep_stage_labels[ "pN2" ] = NREM2;
  // sleep_stage_labels[ "pN3" ] = NREM3;
  // sleep_stage_labels[ "pR" ]  = REM;
  // sleep_stage_labels[ "p?" ]  = UNKNOWN;

  
  // //
  // // SOAP predictions
  // //
  
  // sleep_stage_labels[ "sW" ]  = WAKE;
  // sleep_stage_labels[ "sN1" ] = NREM1;
  // sleep_stage_labels[ "sN2" ] = NREM2;
  // sleep_stage_labels[ "sN3" ] = NREM3;
  // sleep_stage_labels[ "sR" ]  = REM;
  // sleep_stage_labels[ "s?" ]  = UNKNOWN;

  
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

  // other
  sleep_stage_labels[ "Stage1" ] = NREM1;
  sleep_stage_labels[ "Stage2" ] = NREM2;
  sleep_stage_labels[ "Stage3" ] = NREM3;
  sleep_stage_labels[ "Stage4" ] = NREM4;
  sleep_stage_labels[ "S1" ] = NREM1;
  sleep_stage_labels[ "S2" ] = NREM2;
  sleep_stage_labels[ "S3" ] = NREM3;
  sleep_stage_labels[ "S4" ] = NREM4;

  sleep_stage_labels[ "lights" ] = LIGHTS_ON;
  sleep_stage_labels[ "unknown" ] = UNKNOWN;
  sleep_stage_labels[ "missing" ] = UNKNOWN;
  sleep_stage_labels[ "A" ] = UNKNOWN;
  sleep_stage_labels[ "artifact" ] = UNKNOWN;

  sleep_stage_labels[ "G" ] = GAP;
  sleep_stage_labels[ "-" ] = GAP;

  
  // other NSRR
  sleep_stage_labels[ "Wake|0" ]          = WAKE;
  sleep_stage_labels[ "Stage 1 sleep|1" ] = NREM1;
  sleep_stage_labels[ "Stage 2 sleep|2" ] = NREM2;
  sleep_stage_labels[ "Stage 3 sleep|3" ] = NREM3;
  sleep_stage_labels[ "Stage 4 sleep|4" ] = NREM4;
  sleep_stage_labels[ "REM sleep|5" ]     = REM;
  //  sleep_stage_labels[ "Unsure|Unsure" ]   = UNSCORED;

  // Basic
  sleep_stage_labels[ "wake" ]     = WAKE;  
  sleep_stage_labels[ "NREM1" ]    = NREM1;  
  sleep_stage_labels[ "NREM2" ]    = NREM2;  
  sleep_stage_labels[ "NREM3" ]    = NREM3;  
  sleep_stage_labels[ "NREM4" ]    = NREM4;  
  sleep_stage_labels[ "REM" ]      = REM;
  sleep_stage_labels[ "Movement" ] = MOVEMENT;
  sleep_stage_labels[ "Unscored" ] = UNSCORED;

  
  // mouse: make generic 'NR' -> NREM2
  sleep_stage_labels[ "W" ]     = WAKE;  
  sleep_stage_labels[ "NR" ]    = NREM2;  
  sleep_stage_labels[ "R" ]     = REM;
  sleep_stage_labels[ "?" ]     = UNSCORED;
  

  //
  // Channel types
  //

  init_channel_types();

  
  //
  // Time-units
  //
  
  // 1e-9 sec resolution
  tp_1sec  = 1000000000LLU;

  // 1e-6 sec resolution in 1/1000th of a second; used when reading in floating point
  // times in seconds, Helper::sec2tp(), to avoid floating point issues beyond 1/1000th sec 
  tp_1000thsec  = 1000000LLU;

  tp_duration    = 1.0 / (double)tp_1sec;
  
  //
  // Common output stratifiers
  //

  freq_strat   = "F";
  signal_strat = "CH";
  signal1_strat = "CH1";
  signal2_strat = "CH2";
  stage_strat  = "SS";
  cycle_strat  = "C";
  band_strat   = "B";
  annot_strat  = "ANNOT";
  annot_instance_strat  = "INST";
  annot_meta_strat  = "META";
  count_strat  = "N";
  epoch_strat  = "E";
  time_strat   = "T";
  sample_strat = "SP"; // sample-point
  anchor_strat = "ANCHOR";
  cluster_strat = "K";
  var_strat     = "VAR";
  value_strat   = "VAL";
  feature_strat     = "FTR";
  sec_strat = "SEC";
  
  //
  // Misc.
  //

  project_path = "";

#ifdef WINDOWS
  folder_delimiter = '\\';
  mkdir_command = "mkdir";
#else
  folder_delimiter = '/';
  mkdir_command = "mkdir -p";
#endif  

  file_list_delimiter = ',';
  
  txt_table_prepend = "";
  txt_table_append = "";

  force_edf = false;
  skip_edf_annots = false;
  skip_nonedf_annots = false;
  skip_sl_annots = false;
  
  set_annot_inst2hms = false;
  set_annot_inst2hms_force = false;
   
  current_tag = "";

  indiv_wildcard = "^";

  anon_idroot_cnt = 0; // count incremented each time set 
    
  sample_list_min = -1;
  sample_list_max = -1;
  sample_list_id = "";

  anon = false; 

  force_starttime = "";
  force_startdate = "";
  
  write_naughty_list = false;
  naughty_list = "";
  
  assume_pm_starttime = false;
  assume_pm_starttime_hour = 4; // i.e. if 04:00 given or later (up to & including 12:00) assume +12 hours
  
  // otherwise, leave as is
  //  00:00 -> 00:00 (as is)
  //  01:00 -> 01:00 (as is)
  //  03:00 -> 03:00 (as is)
  
  //  04:00 -> 16:00 (*shift)
  //  06:00 -> 18:00 (*shift)
  //  12:00 -> 00:00 (*shift)

  //  13:00 -> 13:00 (as is)
  //  18:00 -> 18:00 (as is)

  problem = false;
  empty = false;
  
  bail_on_fail = true;
  
  edf_timetrack_label = "_TT";
  edf_timetrack_size = 15; // i.e. up to 30 chars
  edf_annot_label = "edf_annot";


  // if running 'collapse' EDF+D --> EDF
  // add annot to show gaps 
  annot_disc_segment = "segment";
  annot_disc_gap = "gap";
  annot_disc_drop_spanning = true;
  // S1 S1 S1 G1 G1 S2 S2 G2 S3 S3

  // S1 S1 S1 S2 S2 S3 S3
  // -------- ----- -----
  // segment1   2     3 
  // gap     1     2   (2-tp duration marker, so hits both epochs)
  
  // whether to assume 30-sec and enfore epoch check when first attaching 
  // an .eannot file;   default is allow up to 5 epochs difference

  enforce_epoch_check = 5;

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
    case DENOM : return "TOTAL";      
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
  if ( sleep_stage_prefix == "" ) 
    {
      sleep_stage_label_lookup_t::const_iterator ii = sleep_stage_labels.find( s );
      if ( ii == sleep_stage_labels.end() ) return UNKNOWN;
      return ii->second;
    }
  
  // otherwise, we are expecting to see the prefix at the string start
  if ( sleep_stage_prefix != s.substr( 0 , sleep_stage_prefix.size() ) ) return UNKNOWN;
  sleep_stage_label_lookup_t::const_iterator ii = sleep_stage_labels.find( s.substr( sleep_stage_prefix.size() ) );
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


std::string globals::map_channel_label( const std::string & s )
{
  return ch2label[ map_channel( s ) ] ;
}

channel_type_t globals::map_channel( const std::string & s )
{

  // for wildcards/case-insensitive matches
  const std::string s2 = Helper::toupper( s );
  
  // special case: for IGNORE first (i.e. for OFF and Status channels)  
  if ( chmap2[ IGNORE ].find( s ) != chmap2[ IGNORE ].end() ) return IGNORE;
  // wildcard match for IGNORE
  const std::set<std::string> & ss = chmap1[ IGNORE ];
  std::set<std::string>::const_iterator jj = ss.begin();
  while ( jj != ss.end() )
    {
      if ( s2.find( *jj ) != std::string::npos ) // any match?
	return IGNORE;
      ++jj;
    }

  // otherwise, try all other TYPES
  // given chmap1 and chmap2, figure out what type this channel 's' is
  // try for an exact match first
  channel_map_t::const_iterator cc = chmap2.begin();
  while ( cc != chmap2.end() )
    {
      if ( cc->second.find( s ) != cc->second.end() )
	return cc->first;
      ++cc;
    }
  
  cc = chmap1.begin();
  while ( cc != chmap1.end() )
    {
      const std::set<std::string> & ss = cc->second;
      std::set<std::string>::const_iterator jj = ss.begin();
      while ( jj != ss.end() )
	{
 	  if ( s2.find( *jj ) != std::string::npos ) // any match?
	    return cc->first;
	  ++jj;
	}
      // next class
      ++cc;
    }
    
  // otherwise, not found...
  return GENERIC;
}


void globals::clear_channel_map()
{
  chmap1.clear();
  chmap2.clear();
}

void globals::add_channel_map( const std::string & s , const std::string & lab )
{
  if ( label2ch.find( lab ) == label2ch.end() ) Helper::halt( "bad channel type: " + lab );
  add_channel_map( s , label2ch[ lab ] );
}

void globals::add_channel_map_exact( const std::string & s , const std::string & lab )
{
  if ( label2ch.find( lab ) == label2ch.end() ) Helper::halt( "bad channel type: " + lab );
  add_channel_map_exact( s , label2ch[ lab ] );
}

void globals::add_channel_map( const std::string & s , channel_type_t ch )
{
  // check not already specified: if it is, erase from other type
  if ( sig2type.find( s ) != sig2type.end() )
    {
      channel_type_t ch = sig2type[ s ];
      // case-insensitive match
      if ( chmap1[ ch ].find( Helper::toupper( s ) ) != chmap1[ ch ].end() )
	chmap1[ ch ].erase( chmap1[ ch ].find( Helper::toupper( s ) ) );
      // exact match
      if ( chmap2[ ch ].find( s ) != chmap2[ ch ].end() )
	chmap2[ ch ].erase( chmap2[ ch ].find( s ) );
    }

  // now added: case-insensitive, so all UPPER
  chmap1[ ch ].insert( Helper::toupper( s ) );
  
}

void globals::add_channel_map_exact( const std::string & s , channel_type_t ch )
{
  // check not already specified: if it is, erase from other type
  if ( sig2type.find( s ) != sig2type.end() )
    {
      channel_type_t ch = sig2type[ s ];
      // case-insensitive match
      if ( chmap1[ ch ].find( Helper::toupper( s ) ) != chmap1[ ch ].end() )
	chmap1[ ch ].erase( chmap1[ ch ].find( Helper::toupper( s ) ) );
      // exact match
      if ( chmap2[ ch ].find( s ) != chmap2[ ch ].end() )
	chmap2[ ch ].erase( chmap2[ ch ].find( s ) );
    }

  chmap2[ ch ].insert( s );
}

void globals::channel_type( const std::string & label , channel_type_t ch )
{
  label2ch[ label ] = ch;
  ch2label[ ch ] = label;
}


std::string globals::list_channels( channel_type_t ch ,
				    const std::vector<std::string> & signals , 
				    const std::string & delim )
{

  std::stringstream ss;
  bool first = true;
  // consider each channel, and extract ones matching type 'ch'
  for (int s=0; s<signals.size(); s++) 
    {
      // match this type?
      if ( map_channel( signals[s] ) == ch )
	{
	  if  ( ! first ) ss << delim;
	  ss << signals[s];
	  first = false;
	}
    }
  return ss.str();
}


void globals::init_channel_types()
{
  channel_type( "EEG" , EEG );
  channel_type( "REF" , REF );
  channel_type( "IC"  , IC );
  channel_type( "IMF" , IMF );
  channel_type( "EOG" , EOG );
  channel_type( "GENERIC" , GENERIC );
  channel_type( "ECG" , ECG );
  channel_type( "EMG" , EMG );
  channel_type( "LEG" , LEG );
  channel_type( "AIRFLOW" , AIRFLOW );
  channel_type( "EFFORT" , EFFORT );
  channel_type( "OXYGEN" , OXYGEN );
  channel_type( "POSITION" , POSITION );
  channel_type( "LIGHT" , LIGHT );
  channel_type( "SNORE" , SNORE );
  channel_type( "HR" , HR );
  channel_type( "IGNORE" , IGNORE );

  // wild-cards (any partial match/case-insensitive)

  add_channel_map( "OFF" , IGNORE );
  add_channel_map( "STATUS" , IGNORE );

  // canonical/base signals

  add_channel_map_exact( "csEEG" , EEG );
  add_channel_map_exact( "csCEN" , EEG );
  add_channel_map_exact( "csFRT" , EEG );
  add_channel_map_exact( "csC3" , EEG );
  add_channel_map_exact( "csC4" , EEG );
  add_channel_map_exact( "csF3" , EEG );
  add_channel_map_exact( "csF4" , EEG );
  add_channel_map_exact( "csO1" , EEG );
  add_channel_map_exact( "csO2" , EEG );
  
  add_channel_map_exact( "csEOG" , EOG );
  add_channel_map_exact( "csLOC" , EOG );
  add_channel_map_exact( "csROC" , EOG );

  add_channel_map_exact( "csEMG" , EMG );

  add_channel_map_exact( "csECG" , ECG );

  add_channel_map_exact( "csCAN" , AIRFLOW );
  add_channel_map_exact( "csTHM" , AIRFLOW );

  add_channel_map_exact( "csTHX" , EFFORT );
  add_channel_map_exact( "csABD" , EFFORT );

  add_channel_map_exact( "csOXY" , OXYGEN );
    
  
  // EEG --- full 64-EEG montage labels to be added in too
  add_channel_map( "EEG" , EEG );
  add_channel_map( "C3" , EEG );
  add_channel_map( "C4" , EEG );
  add_channel_map( "F3" , EEG );
  add_channel_map( "F4" , EEG );
  add_channel_map( "T3" , EEG );
  add_channel_map( "T5" , EEG );
  add_channel_map( "T6" , EEG );
  add_channel_map( "T4" , EEG );
  add_channel_map( "O1" , EEG );
  add_channel_map( "O2" , EEG );
  add_channel_map( "CZ" , EEG );
  add_channel_map( "FZ" , EEG );
  add_channel_map( "PZ" , EEG );
  add_channel_map( "OZ" , EEG );
  add_channel_map( "FPZ" , EEG );

  add_channel_map( "FP2" , EEG );
  add_channel_map( "FP1" , EEG );
  add_channel_map( "AF8" , EEG );
  add_channel_map( "AF7" , EEG );
  add_channel_map( "F8" , EEG );
  add_channel_map( "F6" , EEG );
  add_channel_map( "F2" , EEG );
  add_channel_map( "F1" , EEG );
  add_channel_map( "F5" , EEG );
  add_channel_map( "F7" , EEG );
  add_channel_map( "FC6" , EEG );
  add_channel_map( "FC2" , EEG );
  add_channel_map( "FC1" , EEG );
  add_channel_map( "FC5" , EEG );

  add_channel_map( "T8" , EEG );
  add_channel_map( "C6" , EEG );
  add_channel_map( "C2" , EEG );
  add_channel_map( "C1" , EEG );
  add_channel_map( "C5" , EEG );
  add_channel_map( "T7" , EEG );
  
  add_channel_map( "TP8" , EEG );
  add_channel_map( "CP6" , EEG );
  add_channel_map( "CP4" , EEG );
  add_channel_map( "CP2" , EEG );
  add_channel_map( "CP1" , EEG );
  add_channel_map( "CP3" , EEG );
  add_channel_map( "CP5" , EEG );
  add_channel_map( "TP7" , EEG );
  
  add_channel_map( "P8" , EEG );
  add_channel_map( "P6" , EEG );
  add_channel_map( "P4" , EEG );
  add_channel_map( "P2" , EEG );
  add_channel_map( "P1" , EEG );
  add_channel_map( "P3" , EEG );
  add_channel_map( "P5" , EEG );
  add_channel_map( "P7" , EEG );
  
  add_channel_map( "PO8" , EEG );
  add_channel_map( "PO4" , EEG );
  add_channel_map( "POZ" , EEG );
  add_channel_map( "PO3" , EEG );
  add_channel_map( "PO7" , EEG );

  // REF : note, use exact encoding, as C3-M1 is an EEG channel
  add_channel_map_exact( "M1" , REF );
  add_channel_map_exact( "A1" , REF );
  add_channel_map_exact( "M2" , REF );
  add_channel_map_exact( "A2" , REF );
  
  // IC
  add_channel_map( "IC_" , IC );

  // IMF
  add_channel_map( "IMF_" , IMF );
  
  // EOG
  add_channel_map( "EOG" , EOG );
  add_channel_map( "LOC" , EOG );
  add_channel_map( "ROC" , EOG );
  add_channel_map( "E1" , EOG );
  add_channel_map( "E2" , EOG );
  
  // ECG
  add_channel_map( "ECG" , ECG );
  add_channel_map( "EKG" , ECG );
  add_channel_map_exact( "LA" , ECG );
  add_channel_map_exact( "RA" , ECG );
  add_channel_map_exact( "LL" , ECG );

  // EMG
  add_channel_map( "EMG" , EMG );
  add_channel_map( "CHIN" , EMG );

  // LEG
  add_channel_map( "LEG" , LEG );
  add_channel_map( "LAT" , LEG );
  add_channel_map( "RAT" , LEG );
  
  // AIRFLOW
  add_channel_map( "FLOW" , AIRFLOW );
  add_channel_map( "NASAL" , AIRFLOW );
  add_channel_map( "THERM" , AIRFLOW );


  // EFFORT
  add_channel_map( "ABD" , EFFORT );
  add_channel_map( "CHEST" , EFFORT );
  add_channel_map( "THOR" , EFFORT );
  add_channel_map( "SUM" , EFFORT );

  // OXYGEN
  add_channel_map( "SPO2" , OXYGEN );
  add_channel_map( "SAO2" , OXYGEN );
  add_channel_map( "SP02" , OXYGEN );
  add_channel_map( "SA02" , OXYGEN );
  add_channel_map( "OX"   , OXYGEN ); 
  
  // HR
  add_channel_map( "HR" , HR );
  add_channel_map_exact( "HRate" , HR );
  add_channel_map( "PULSE" , HR );
  add_channel_map_exact( "PR" , HR );
      
  // POSITION
  add_channel_map( "POS" , POSITION );

  // LIGHT
  add_channel_map( "LIGHT" , LIGHT );

  // SNORE
  add_channel_map( "SNORE" , SNORE );

  // GENERIC (add to avoid clash w/ EEG because of O2)
  add_channel_map( "etco2" , GENERIC );
  add_channel_map( "etc02" , GENERIC );

  add_channel_map( "DIF5" , GENERIC );
  add_channel_map( "DIF6" , GENERIC );

  add_channel_map( "DC1" , GENERIC );
  add_channel_map( "DC2" , GENERIC );
  add_channel_map( "DC3" , GENERIC );
  add_channel_map( "DC4" , GENERIC );
  add_channel_map( "DC5" , GENERIC );
  add_channel_map( "DC6" , GENERIC );
  add_channel_map( "DC7" , GENERIC );
  add_channel_map( "DC8" , GENERIC );
  add_channel_map( "DC9" , GENERIC );
  add_channel_map( "DC10" , GENERIC );
  
}


std::string globals::dump_channel_map()
{
  std::stringstream ss;

  channel_map_t::const_iterator cc = chmap2.begin();
  while ( cc != chmap2.end() )
    {
      const std::set<std::string> & kk = cc->second;
      std::set<std::string>::const_iterator jj = kk.begin();
      while ( jj != kk.end() )
	{
	  ss << "EXACT\t" << *jj << "\t" << ch2label[ cc->first ] << "\n";
	  ++jj;
	}
      ++cc;
    }

  cc = chmap1.begin();
  while ( cc != chmap1.end() )
    {
      const std::set<std::string> & kk = cc->second;
      std::set<std::string>::const_iterator jj = kk.begin();
      while ( jj != kk.end() )
	{
	  ss << "PARTIAL\t" << *jj << "\t" << ch2label[ cc->first ] << "\n";
	  ++jj;
	}
      ++cc;
    }
    
    return ss.str();
}
