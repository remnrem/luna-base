
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

#include "luna.h"
#include "db/db.h"
#include <iomanip>

#include "optdefs.h"

extern globals global;

optdefs_t::optdefs_t()
{
  init();
}

void optdefs_t::init()
{
  domains.clear();
  domain2opt.clear();
  odesc.clear();
  otype.clear();

  // register all known options
  
  // ID inclusion/exclusion
  globals::optdefs().add( "inputs", "id" , OPT_STRVEC_T , "Select these IDs from the sample list" );
  globals::optdefs().add( "inputs", "skip" , OPT_STRVEC_T , "Skip these IDs from the sample list" );  
  globals::optdefs().add( "inputs", "vars" , OPT_FILE_T , "Specify file(s) of individual-specific variables" );
  globals::optdefs().add( "inputs", "ids" , OPT_STRVEC_T , "Select only one or more IDs from the sample list" );
  globals::optdefs().add( "inputs", "exclude" , OPT_FILE_T , "File of IDs to exclude" );
  globals::optdefs().add( "inputs", "include" , OPT_FILE_T , "File of IDs to include" );
  globals::optdefs().add( "inputs", "path" , OPT_PATH_T , "Set project/sample-list 'current path' (i.e. for relative sample list file)");
  globals::optdefs().add( "inputs", "preload" , OPT_BOOL_T , "Read all EDF(+) records on firat attaching" );

  // logging
  globals::optdefs().add( "logging" , "verbose" , OPT_BOOL_T , "Set verbose logging" );
  globals::optdefs().add( "logging" , "silent" , OPT_BOOL_T , "Suppress console logging" );
  globals::optdefs().add( "logging" , "log" , OPT_FILE_T , "Write console output to this file" );
  globals::optdefs().add( "logging" , "write-log" , OPT_BOOL_T , "Turn off log-saving (i.e. used via API)" );
  globals::optdefs().add( "logging" , "mirror" , OPT_BOOL_T , "Mirror inputs in console log" );


  // reading EDFs
  globals::optdefs().add( "signals", "force-edf" , OPT_BOOL_T , "Read EDF+ as EDF (i.e. ignore gaps, annotations)");
  globals::optdefs().add( "signals", "sig" , OPT_STRVEC_T , "One or more signals to import from the EDF" );
  globals::optdefs().add( "signals", "order-signals" , OPT_BOOL_T , "Order signals alphabetically" );
  globals::optdefs().add( "signals", "anon" , OPT_BOOL_T , "Do not read IDs in" );
  globals::optdefs().add( "signals", "fix-edf" , OPT_BOOL_T , "Attempt to correct truncated/extended EDFs" );  
  globals::optdefs().add( "signals", "digital" , OPT_BOOL_T , "(debug-mode) read digital values, do not map to physical values" );
  globals::optdefs().add( "signals", "force-digital-minmax" , OPT_BOOL_T , "(debug-mode) Force digitial min/max (if it is invalid)" );
  globals::optdefs().add( "signals", "force-digital-min" , OPT_INT_T , "(debug-mode) Force digitial min value" );
  globals::optdefs().add( "signals", "force-digital-max" , OPT_INT_T , "(debug-mode) Force digitial max value" );
  
  // annotations
  globals::optdefs().add( "annotations", "keep-annot-spaces" , OPT_BOOL_T , "Keep spaces as is for annotation labels" );
  globals::optdefs().add( "annotations", "add-ellipsis" , OPT_BOOL_T , "Mark 0-duration intervals as '...' (WRITE-ANNOTS --> .annot)" );  
  globals::optdefs().add( "annotations", "class-instance-delimiter" , OPT_CHAR_T , "Annotation class/instance delimiter, default = :" );
  globals::optdefs().add( "annotations", "combine-annots" , OPT_CHAR_T , "Combine class/instance delimiter, default = _" );
  globals::optdefs().add( "annotations", "annot-segment" , OPT_STR_T , "Annotation label to mark EDF+D segments; default = segment" );
  globals::optdefs().add( "annotations", "annot-gap" , OPT_STR_T , "Annotation label to mark EDF+D gaps; default = gap" );
  globals::optdefs().add( "annotations", "annot-whitelist" , OPT_BOOL_T , "Skip annotations not whitelisted (i.e. with an explicit remap)" );
  globals::optdefs().add( "annotations", "annot-unmapped" , OPT_BOOL_T , "Skip whitelisted annotations (i.e. without an explicit remap)" );  
  globals::optdefs().add( "annotations", "edf-annot-class" , OPT_STRVEC_T , "For EDF+ annotations, treat these are full classes" );
  globals::optdefs().add( "annotations", "edf-annot-class-all" , OPT_STRVEC_T , "Treat all EDF+ annotations as full classes (edf-annot-class=*)" );
  globals::optdefs().add( "annotations", "tab-only" , OPT_BOOL_T , "Set to F to allow space-delimiters in .annot files" );
  globals::optdefs().add( "annotations", "inst-hms" , OPT_BOOL_T , "If T, set blank annotation instances to hh:mm:ss" );
  globals::optdefs().add( "annotations", "force-inst-hms" , OPT_BOOL_T , "If T, force all annotation instances to hh:mm:ss" );
  globals::optdefs().add( "annotations", "skip-edf-annots" , OPT_BOOL_T , "Skip any EDF+ annotations" );
  globals::optdefs().add( "annotations", "skip-sl-annots" , OPT_BOOL_T , "Skip any sample-list annotations" );
  globals::optdefs().add( "annotations", "skip-annots" , OPT_BOOL_T , "Skip all sample-list annotations [also skip-all-annots]" );
  globals::optdefs().add( "annotations", "annot-file" , OPT_FILE_T , "One or more additional annotation files" );
  globals::optdefs().add( "annotations", "annot" , OPT_STRVEC_T , "Only load these annotations(s) based on class ID ('.' for all)" );
  globals::optdefs().add( "annotations", "raw-annot" , OPT_STRVEC_T , "As annot option, but without label sanitization ('.' for all)" );
  globals::optdefs().add( "annotations", "ignore-annot" , OPT_STRVEC_T , "Exclude these annotations ('.' for all)" );
  globals::optdefs().add( "annotations", "ignore-raw-annot" , OPT_STRVEC_T , "As ignore-annot, but w/out label sanitization" );
  globals::optdefs().add( "annotations", "annot-remap" , OPT_BOOL_T , "If false, wipe all stage/preloaded annotation remappings" );
  globals::optdefs().add( "annotations", "nsrr-remap" , OPT_BOOL_T , "If true, add in extra NSRR-centric annotation remappings" );

  // annotation meta-date
  globals::optdefs().add( "metadata", "annot-keyval" , OPT_CHAR_T , "Assignment for key=value annotation meta-data (default '=')");
  globals::optdefs().add( "metadata", "annot-meta-delim1" , OPT_CHAR_T , "Delimiter for annotation meta-data (default '|')");
  globals::optdefs().add( "metadata", "annot-meta-delim2" , OPT_CHAR_T , "Alternate delimiter for annotation meta-data (default ';')");
  globals::optdefs().add( "metadata", "annot-meta-default-num" , OPT_BOOL_T , "Assume annotation meta-data is numeric unless told otherwise (default T)" );
  globals::optdefs().add( "metadata", "num-atype" , OPT_STRVEC_T , "Specify numeric metadata type(s)" );
  globals::optdefs().add( "metadata", "str-atype" , OPT_STRVEC_T , "Specify string metadata type(s) [or txt-atype]" );
  globals::optdefs().add( "metadata", "int-atype" , OPT_STRVEC_T , "Specify integer metadata type(s)" );
  globals::optdefs().add( "metadata", "bool-atype" , OPT_STRVEC_T , "Specify boolean metadata type(s)" );

  // aliasing/remapping
  globals::optdefs().add( "aliasing", "alias" , OPT_STR_T , "Channel alias primary|alias1|alias2|..." );
  globals::optdefs().add( "aliasing", "remap" , OPT_STR_T , "Specifying annotation remapping: primary|alias1|alias2|..." );
  globals::optdefs().add( "aliasing", "retain-case" , OPT_BOOL_T , "If aliasing a primary, retain input case (default = T)" );  
  globals::optdefs().add( "aliasing", "sanitize" , OPT_BOOL_T , "Sanitize labels (signals & annots)" );
  globals::optdefs().add( "aliasing", "spaces" , OPT_CHAR_T , "Character to replace spaces with; default = _" );
  globals::optdefs().add( "aliasing", "upper" , OPT_BOOL_T , "Set signal labels to uppercase" );
  globals::optdefs().add( "aliasing", "keep-spaces" , OPT_BOOL_T , "Keep spaces as is for channel & annotation labels" );
  globals::optdefs().add( "aliasing", "keep-channel-spaces" , OPT_BOOL_T , "Keep spaces as is for channel labels" );

  // epochs
  globals::optdefs().add( "epochs", "epoch-check" , OPT_INT_T , "Tolerance of EDF/.eannot epoch check (default 5)" );
  globals::optdefs().add( "epochs", "epoch-len" , OPT_INT_T , "Set default epoch length (seconds, default 30)" );

  // times/dates
  globals::optdefs().add( "datetimes", "date-format" , OPT_SPECIAL_T , "Set input date format (MDY, DMY or YMD)" );
  globals::optdefs().add( "datetimes", "write-date-format" , OPT_SPECIAL_T , "Set output date format (MDY, DMY or YMD) - currently not used" );
  globals::optdefs().add( "datetimes", "edf-date-format" , OPT_SPECIAL_T , "Set input EDF header date format (MDY, DMY or YMD)" );
  globals::optdefs().add( "datetimes", "starttime" , OPT_TIME_T , "Force EDF start time" );
  globals::optdefs().add( "datetimes", "startdate" , OPT_DATE_T , "Force EDF start date" );
  globals::optdefs().add( "datetimes", "assume-pm-start" , OPT_INT_T , "Shift +12 hours if EDF starts after in AM at or after this hour (set 'n' to turn off)" );
  globals::optdefs().add( "datetimes", "default-starttime" , OPT_TIME_T , "Set default start time for EDFs" );
  globals::optdefs().add( "datetimes", "no-default-starttime" , OPT_BOOL_T , "Do not apply a default start time" );
  globals::optdefs().add( "datetimes", "default-startdate" , OPT_DATE_T , "Set default start date for EDFs" );
  globals::optdefs().add( "datetimes", "no-default-startdate" , OPT_BOOL_T , "Do not apply a default start date" );
  globals::optdefs().add( "datetimes", "sec-dp" , OPT_INT_T , "Set decimal places for certain time outputs" );

  
  // scripts
  globals::optdefs().add( "scripts", "wildcard" , OPT_CHAR_T , "Set ID wildcard character; default = ^" );
  globals::optdefs().add( "scripts", "show-assignments" , OPT_BOOL_T , "Log all variable assignments" );

  // staging
  globals::optdefs().add( "staging", "ss-prefix" , OPT_STR_T , "Set sleep-stage prefix (e.g. pN1, pN2, etc)" );
  globals::optdefs().add( "staging", "ss-pops" , OPT_BOOL_T , "Implies ss-prefix=p" );
  globals::optdefs().add( "staging", "ss-soap" , OPT_BOOL_T , "Implies ss-prefix=s" );
  globals::optdefs().add( "staging", "assume-stage-duration" , OPT_BOOL_T , "Assume 0-dur sleep stages are of epoch-duration" );

  // channel types
  globals::optdefs().add( "types", "ch-match" , OPT_STRVEC_T , "Specify partial channel-type match(es) (type|label1|label2)" );
  globals::optdefs().add( "types", "ch-exact" , OPT_STRVEC_T , "Specify exact channel-type match(es) (type|label1|label2)" );
  globals::optdefs().add( "types", "ch-clear" , OPT_BOOL_T , "Wipe channel-type information" );
  
  // general output options
  globals::optdefs().add( "outputs", "tt-prefix" , OPT_STR_T , "Tag to add to -t output filenames [also tt-prepend]" );
  globals::optdefs().add( "outputs", "tt-suffix" , OPT_STR_T , "Tag to add at end of -t output filenames [also tt-append]" );    
  globals::optdefs().add( "outputs", "compressed" , OPT_BOOL_T , "Compress (gzip) all -t output" );

  // stats/numeric
  globals::optdefs().add( "numeric", "srand" , OPT_INT_T , "Set random seed (long unsigned int)" );
  globals::optdefs().add( "numeric", "legacy-hjorth" , OPT_BOOL_T , "Use legacy Hjorth complexity calculation" );
  globals::optdefs().add( "numeric", "slow" , OPT_NUM_INTERVAL_T , "Set SLOW [lwr,upr) band (default 0.5-1)" );
  globals::optdefs().add( "numeric", "delta" , OPT_NUM_INTERVAL_T , "Set DELTA [lwr,upr) band (default 1-4)" );
  globals::optdefs().add( "numeric", "theta" , OPT_NUM_INTERVAL_T , "Set THETA [lwr,upr) band (default 4-8)" );
  globals::optdefs().add( "numeric", "alpha" , OPT_NUM_INTERVAL_T , "Set ALPHA [lwr,upr) band (default 8-11)" );
  globals::optdefs().add( "numeric", "sigma" , OPT_NUM_INTERVAL_T , "Set SIGMA [lwr,upr) band (default 11-15)" );
  globals::optdefs().add( "numeric", "beta" , OPT_NUM_INTERVAL_T , "Set BETA [lwr,upr) band (default 15-30)" );
  globals::optdefs().add( "numeric", "gamma" , OPT_NUM_INTERVAL_T , "Set GAMMA [lwr,upr) band (default 30-50)" );
  globals::optdefs().add( "numeric", "total" , OPT_NUM_INTERVAL_T , "Set TOTAL [lwr,upr) band (default 0.5-50)" );

  // misc / possibly legacy
  globals::optdefs().add( "misc", "fail-list" , OPT_FILE_T , "Write failing IDs to this file" );
  globals::optdefs().add( "misc", "bail-on-fail" , OPT_BOOL_T , "Behavior when an internal error flag is raised" );
  globals::optdefs().add( "misc", "align-annots" , OPT_STRVEC_T , "Align annotations [check still in use]" );
  

  
}

void optdefs_t::add( const std::string & domain, const std::string & opt , const opt_type_t & type , const std::string & desc )
{
  // retain domain-ordering based on first exposure
  if ( domain2opt.find( domain ) == domain2opt.end() )
    domains.push_back( domain );
  domain2opt[ domain ].push_back( opt );
  odesc[ opt ] = desc;
  otype[ opt ] = type;
}



opt_type_t optdefs_t::get_type( const std::string & t ) const
{
  std::map<std::string,opt_type_t>::const_iterator ii = otype.find( t );
  if ( ii == otype.end() ) return OPT_UNDEFINED_T;
  return ii->second;
}

std::string optdefs_t::get_desc( const std::string & t ) const
{
  std::map<std::string,std::string>::const_iterator ii = odesc.find( t );
  if ( ii == odesc.end() ) return ".";
  return ii->second;
}

std::vector<std::string> optdefs_t::get_opts( const std::string & domain ) const
{
  std::vector<std::string> res;
  std::map<std::string,std::vector<std::string> >::const_iterator ii = domain2opt.find( domain );
  if ( ii == domain2opt.end() ) return res;
  return ii->second;
}





