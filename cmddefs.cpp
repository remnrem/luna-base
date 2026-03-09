
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
#include <cctype>
#include <iomanip>

//#pragma GCC push_options
//#pragma GCC optimize ("O0")

#include "cmddefs.h"

extern globals global;

namespace {

const int verbose_help_wrap = 80;

std::string wrap_text( const std::string & text , const int width )
{
  std::stringstream in( text );
  std::stringstream out;
  std::string line;
  bool first_paragraph = true;

  while ( std::getline( in , line ) )
    {
      bool blank = true;
      for (int i=0; i<line.size(); i++)
        if ( ! std::isspace( static_cast<unsigned char>( line[i] ) ) )
          {
            blank = false;
            break;
          }

      if ( blank )
        {
          out << "\n";
          first_paragraph = false;
          continue;
        }

      std::stringstream words( line );
      std::string word;
      int col = 0;

      if ( ! first_paragraph ) out << "\n";
      first_paragraph = false;

      while ( words >> word )
        {
          if ( col == 0 )
            {
              out << word;
              col = word.size();
            }
          else if ( col + 1 + word.size() <= width )
            {
              out << " " << word;
              col += 1 + word.size();
            }
          else
            {
              out << "\n" << word;
              col = word.size();
            }
        }
    }

  return out.str();
}

}


cmddefs_t::cmddefs_t()
{
  //init();
}


void cmddefs_t::init()
{

  //
  // parameters
  //
  
  allz = false;
  nonez = false;
  
  /////////////////////////////////////////////////////////////////////////
  //
  // Document domains, commands, parameters, output tables and variables 
  //
  /////////////////////////////////////////////////////////////////////////


  //
  // base URL
  //

  url_root = "https://zzz.nyspi.org/luna/ref/";


  //
  // Domains
  //

  add_domain( "summ"       , "Summaries"        , "Basic summary commands" );
  add_domain( "annot"      , "Annotations"      , "Adding and displaying annotations" );
  add_domain( "expr"       , "Expressions"      , "Evaluating more advanced annotation-based expressions" );
  add_domain( "epoch"      , "Epochs"           , "Epoching signals and epoch-level annotations" );
  add_domain( "mask"       , "Masks"            , "Masking epochs based on annotations and other criteria" );
  add_domain( "freeze"     , "Freezes & caches" , "EDF freezes and cache mechanisms" );
  add_domain( "canon"      , "Canonoical signals" , "Canonical signal mapping" );
  add_domain( "manip"      , "Manipulations"    , "Manipulating signal data" );
  // add_domain( "align"      , "Record alignment" , "Signal/annotation alignment" );
  add_domain( "output"     , "Outputs"          , "Commands to output signals in different formats" );
  add_domain( "filter"     , "FIR filters"      , "FIR filter design and application" );
  add_domain( "artifact"   , "Artifacts"        , "Artifacts detection/correction routines" );
  add_domain( "physio"     , "Physiology"       , "ECG, EMG and related physiological analyses" );
  add_domain( "hypno"      , "Hypnograms"       , "Characterizations of hypnograms" );
  add_domain( "stage"      , "Staging"          , "Automated staging/stage evaluation" );
  add_domain( "actig"      , "Actigraphy"       , "Actigraphy: circadian metrics and wake/sleep scoring" );
  add_domain( "power"      , "Time/frequency analysis" , "TF including power spectral density estimation" );
  add_domain( "trans"      , "NREM transients (spindle/SO)"  , "Spindles and slow oscillations" );
  add_domain( "cc"         , "Coupling/connectvitiy" , "Coherence and other topographical analyses" );
  add_domain( "interval"   , "Interval-based analysis" , "Analyses and summaries based on time-domain intervals" );
  add_domain( "psc"        , "Principal spectral components" , "PSC command" );
  add_domain( "spatial"    , "Topographical analysis" , "EEG channel locations, interpolation and surface Laplacian" );
  add_domain( "multi"      , "Multi-channel analysis" , "ICA and PCA" );
  add_domain( "ms"         , "EEG microstate analysis" , "Segmentation, backfitting and sequence analysis" );
  add_domain( "cluster"    , "Clustering"  , "PDC-based epoch/channel clustering " ); 
  add_domain( "assoc"      , "Association"      , "Association models" );
  add_domain( "pred"       , "Prediction"       , "Prediction models" );
  add_domain( "simul"      , "Simulation"       , "Basic signal simulation" );
  add_domain( "helpers"    , "Helper utilities" , "Misc. utility functions" );
  add_domain( "exp"        , "Experimental"     , "Experimental features: under heavy development, for internal use only" );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // SUMMARIES
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // DESC
  //

  add_cmd( "summ" , "DESC" , "Simple description of an EDF, sent to the console" );
  add_verb( "DESC" ,
	    "DESC prints a brief human-readable summary of the EDF to the console.\n"
	    "\n"
	    "Use it for a quick sanity check of the file contents. With channels specified,\n"
	    "it can also act as a simple channel-name dump." );
  add_param( "DESC" , "channels" , "" , "Only write channel names, one-per-line" );

  //
  // SUMMARY
  //

  add_cmd( "summ" , "SUMMARY" , "More verbose description, sent to the console" );
  add_verb( "SUMMARY" ,
	    "SUMMARY prints a longer human-readable report on the EDF to the console.\n"
	    "\n"
	    "Compared to DESC, it is intended as a more narrative overview of the file\n"
	    "contents rather than structured tabular output." );

  //
  // HEADERS
  //

  add_cmd( "summ" , "HEADERS" , "Tabulate (channel-specific) EDF header information" );
  add_url( "HEADERS" , "summaries/#headers" );
  add_verb( "HEADERS" ,
	    "HEADERS emits structured EDF header information as Luna output tables.\n"
	    "\n"
	    "It reports both file-level metadata, such as start/stop times and duration,\n"
	    "and channel-level metadata, such as sample rate, scaling and channel type." );
  add_param( "HEADERS" , "sig" , "C3,C4" , "Restrict output to these signals" );
  add_param( "HEADERS" , "signals" , "" , "Include a comma-delimited list of selected signals" );
  
  add_table( "HEADERS" , "" , "Basic EDF header information" );
  add_var( "HEADERS" , "" , "NR" , "Number of records" );
  add_var( "HEADERS" , "" , "NS" , "Number of signals/channels" );
  add_var( "HEADERS" , "" , "EDF_ID" , "ID in the EDF header" );
  add_var( "HEADERS" , "" , "START_TIME" , "Start time in the EDF header" );
  add_var( "HEADERS" , "" , "STOP_TIME" , "Stop time" );
  add_var( "HEADERS" , "" , "START_DATE" , "Start date in the EDF header" );
  add_var( "HEADERS" , "" , "REC_DUR" , "Duration of each record (seconds)" );
  add_var( "HEADERS" , "" , "TOT_DUR_SEC" , "Current EDF duration (seconds)" );
  add_var( "HEADERS" , "" , "TOT_DUR_HMS" , "Current EDF duration (hh:mm:ss)" );
  add_var( "HEADERS" , "" , "EDF_TYPE" , "EDF, EDF+C or EDF+D" );
  add_var( "HEADERS" , "" , "NS_ALL" , "Number of signals in original EDF" );
  add_var( "HEADERS" , "" , "REC_DUR_HMS" , "Original recording duration (hh:mm:ss)" );
  add_var( "HEADERS" , "" , "REC_DUR_SEC" , "Original recording duration (seconds)" );
  add_var( "HEADERS" , "" , "STOP_DATE" , "Stop date" );
  add_var( "HEADERS" , "" , "SIGNALS" , "Comma-delimited list of selected signals [signals]" );
  
  add_table( "HEADERS" , "CH" , "Per-channel header information" );
  add_var( "HEADERS" , "CH" , "DMAX" , "Digital max" );
  add_var( "HEADERS" , "CH" , "DMIN" , "Digital min" );
  add_var( "HEADERS" , "CH" , "PDIM", "Physical dimension" );
  add_var( "HEADERS" , "CH" , "PMAX", "Physical max" );
  add_var( "HEADERS" , "CH" , "PMIN", "Physical min" );
  add_var( "HEADERS" , "CH" , "SR", "Sample rate (Hz)" );
  add_var( "HEADERS" , "CH" , "SENS", "Sensitivity (unit/bit)" );
  add_var( "HEADERS" , "CH" , "TRANS", "Transducer type" );
  add_var( "HEADERS" , "CH" , "POS", "Position in EDF" );
  add_var( "HEADERS" , "CH" , "TYPE", "Channel type (from Luna TYPES)" );

  //
  // TYPES
  //

  add_cmd( "summ" , "TYPES" , "Show channel type mappings" );
  add_url( "TYPES" , "summaries/#types" );
  add_verb( "TYPES" ,
	    "TYPES reports Luna's current mapping from channel labels to channel types.\n"
	    "\n"
	    "This is useful for checking how canonicalization and type inference have\n"
	    "classified the signals in the current EDF." );

  //
  // VARS
  //

  add_cmd( "summ" , "VARS" , "Dump currently defined variables" );
  add_url( "VARS" , "summaries/#vars" );
  add_verb( "VARS" ,
	    "VARS lists currently defined Luna variables for the current individual.\n"
	    "\n"
	    "This includes variables created earlier in the command stream, whether they\n"
	    "are individual-specific, and their current values." );
  add_table( "VARS" , "VAR" , "Defined variables" );
  add_var( "VARS" , "VAR" , "INDIV" , "Indicator for an individual-specific variable (0/1)" );
  add_var( "VARS" , "VAR" , "VAL" , "Variable value" );

  //
  // DUPES
  //

  add_cmd( "summ" , "DUPES" , "Find near-duplicate signals" );
  add_url( "DUPES" , "summaries/#dupes" );
  add_verb( "DUPES" ,
	    "DUPES looks for channels that are effectively duplicates, as well as flat or\n"
	    "otherwise invalid signals.\n"
	    "\n"
	    "It can compare signals using either digital or physical values, and is mainly\n"
	    "intended as a QC step for detecting redundant or problematic channels." );
  add_param( "DUPES" , "sig" , "C3,C4" , "Restrict analysis to these channels" );
  add_param( "DUPES" , "physical" , "" , "Compare physical values instead of digital values" );
  add_param( "DUPES" , "eps" , "0.01" , "Epsilon used for physical-value comparisons" );
  add_param( "DUPES" , "prop" , "0.1" , "Minimum discordant proportion required to call channels different" );
  add_table( "DUPES" , "" , "Duplicate and flat-signal summary" );
  add_var( "DUPES" , "" , "INVALID" , "Number of signals with invalid ranges" );
  add_var( "DUPES" , "" , "FLAT" , "Number of flat signals" );
  add_var( "DUPES" , "" , "DUPES" , "Number of duplicated channel pairs" );
  add_table( "DUPES" , "CH" , "Per-channel duplicate and flat-signal flags" );
  add_var( "DUPES" , "CH" , "INVALID" , "Signal has an invalid range" );
  add_var( "DUPES" , "CH" , "FLAT" , "Signal is flat" );
  add_var( "DUPES" , "CH" , "DUPE" , "Signal is part of at least one duplicated pair" );
  add_table( "DUPES" , "CHS" , "Duplicated channel pairs" );
  add_var( "DUPES" , "CHS" , "DUPE" , "Channel pair is duplicated" );

  //
  // CONTAINS
  //

  add_cmd( "summ" , "CONTAINS" , "Tests for particular signals/annotations/staging being present" );
  add_url( "CONTAINS" , "summaries/#contains" );
  add_verb( "CONTAINS" ,
	    "CONTAINS checks whether required channels, annotations or valid sleep staging\n"
	    "are present in the current EDF.\n"
	    "\n"
	    "It can be used as a gatekeeper in command streams, either by setting a variable\n"
	    "or by causing the current EDF to be skipped when requirements are not met." );
  add_param( "CONTAINS" , "sig" , "EMG,ECG" , "Test for these signals" );
  add_param( "CONTAINS" , "annot" , "apnea,hypopnea" , "Test for these annotations" );
  add_param( "CONTAINS" , "annots" , "apnea,hypopnea" , "Test for these annotations" );
  add_param( "CONTAINS" , "stages" , "" , "Test for valid staging" );
  add_param( "CONTAINS" , "var" , "HAS_SIGS" , "Set this variable to T/F instead of using the return code" );
  add_param( "CONTAINS" , "skip" , "" , "Skip to next EDF on failure" );
  add_param( "CONTAINS" , "skip-if-none" , "" , "Skip only if none of the requested signals are present" );

  add_table( "CONTAINS" , "" , "Base" );
  add_var( "CONTAINS" , "" , "OVERLAP" , "Sleep stage annotations overlapped one another (0/1)" );
  add_var( "CONTAINS" , "" , "STAGES" , "Valid staging present (0/1)" );
  add_var( "CONTAINS" , "" , "STAGE_COUNTS" , "Sleep stage counts" );
  add_var( "CONTAINS" , "" , "UNIQ_STAGES" , "Number of unique stage labels" );
  add_var( "CONTAINS" , "" , "NA_REQ" , "Number of required annots" );
  add_var( "CONTAINS" , "" , "NA_OBS" , "Number of observed annots" );

  add_var( "CONTAINS" , "" , "NS_OBS" , "Number of observed channels" );
  add_var( "CONTAINS" , "" , "NS_REQ" , "Number of required channels" );
  add_var( "CONTAINS" , "" , "NS_TOT" , "Total number of channels" );

  add_table( "CONTAINS" , "ANNOT" , "Annotation information" );
  add_var( "CONTAINS" , "ANNOT" , "PRESENT" , "Annotation present" );

  add_table( "CONTAINS" , "CH" , "Channel information" );
  add_var( "CONTAINS" , "CH" , "PRESENT" , "Channel present" );

  //
  // ALIASES
  //

  add_cmd( "summ" , "ALIASES" , "Tabulate channel and annotation alias replacements" );
  add_url( "ALIASES" , "summaries/#aliases" );
  add_verb( "ALIASES" ,
	    "ALIASES reports channel and annotation label substitutions that Luna has applied.\n"
	    "\n"
	    "This is mainly useful for auditing remapping and harmonization steps, so you can\n"
	    "see the original labels that were replaced by the current canonical labels." );
  
  add_table( "ALIASES" , "CH" , "Channel aliasing" );
  add_var( "ALIASES" , "CH" , "ORIG" , "Original channel label in EDF" );

  add_table( "ALIASES" , "ANNOT" , "Annotation aliasing" );
  add_var( "ALIASES" , "ANNOT" , "ORIG" , "Original annotation label" );

  //
  // TAG
  //

  add_cmd( "summ" , "TAG" , "Generic command to add a tag (level/factor) to the output" );  
  add_url( "TAG" , "summaries/#tag" );
  add_verb( "TAG" ,
	    "TAG adds a user-defined factor/level pair to downstream Luna output.\n"
	    "\n"
	    "It does not analyze the EDF itself; instead it annotates subsequent output so\n"
	    "results can be stratified or grouped by run, condition or any other label." );
  add_param( "TAG" , ""    , "RUN/L1" , "Add tag with level L1 to factor RUN in output" );
  add_param( "TAG" , "tag" , "RUN/L1" , "Identical to the above, but explicitly using the tag option" );

  //
  // STATS
  //

  add_cmd( "summ"   , "STATS" , "Basic signal statistics (min/max, mean, RMS, etc)" );
  add_url( "STATS" , "summaries/#stats" );
  add_verb( "STATS" ,
	    "STATS computes basic descriptive statistics for one or more signals.\n"
	    "\n"
	    "It can summarize whole-signal distributions, optionally emit per-epoch values,\n"
	    "tabulate digital encodings, and derive simple dynamics summaries from epoch-level\n"
	    "statistics." );
  add_param( "STATS" , "sig" , "C3,C4" , "Restrict analysis to these channels" );
  add_param( "STATS" , "epoch" , "" , "Calculate per-epoch statistics" );
  add_param( "STATS" , "encoding" , "" , "Tabulate the digital encoding distribution" );
  add_param( "STATS" , "min" , "" , "Minimal output with mean only" );
  add_param( "STATS" , "minimal" , "" , "Minimal output with mean only" );
  add_param( "STATS" , "pct" , "F" , "Toggle percentile output (default T)" );
  add_param( "STATS" , "kurt3" , "" , "Report classical kurtosis instead of excess kurtosis" );
  add_param( "STATS" , "sr-under" , "128" , "Restrict analysis to channels at or below this sampling rate" );
  add_param( "STATS" , "dynam" , "" , "Run quantile dynamics on epoch-level statistics" );
  
  add_table( "STATS" , "CH" , "Whole-night, per-channel statistics, based on all epochs" );
  add_var( "STATS" , "CH" , "MEAN_MD" , "Median of epoch-level means [epoch]" );
  add_var( "STATS" , "CH" , "MEDIAN_MD" , "Median of epoch-level medians [epoch]" );
  add_var( "STATS" , "CH" , "RMS_MN" , "Mean of epoch-level RMS values [epoch]" );
  add_var( "STATS" , "CH" , "SKEW_MN" , "Mean of epoch-level skewness values [epoch]" );
  add_var( "STATS" , "CH" , "KURT_MN" , "Mean of epoch-level kurtosis values [epoch]" );
  add_var( "STATS" , "CH" , "RMS_MD" , "Median of epoch-level RMS values [epoch]" );
  add_var( "STATS" , "CH" , "SKEW_MD" , "Median of epoch-level skewness values [epoch]" );
  add_var( "STATS" , "CH" , "KURT_MD" , "Median of epoch-level kurtosis values [epoch]" );
  add_var( "STATS" , "CH" , "MIN" , "Signal minimum (from data, not EDF header)" );
  add_var( "STATS" , "CH" , "MAX" , "Signal maximum (from data, not EDF header)" );
  add_var( "STATS" , "CH" , "MEAN" , "Signal mean" );
  add_var( "STATS" , "CH" , "RMS" , "Signal root mean square" );
  add_var( "STATS" , "CH" , "SKEW" , "Signal skewness" );
  add_var( "STATS" , "CH" , "KURT" , "Signal kurtosis" );
  add_var( "STATS" , "CH" , "SD" , "Signal SD" );

  add_var( "STATS" , "CH" , "P01" , "1st percentile" );
  add_var( "STATS" , "CH" , "P02" , "2nd percentile" );
  add_var( "STATS" , "CH" , "P05" , "5th percentile" );
  add_var( "STATS" , "CH" , "P10" , "10th percentile" );
  add_var( "STATS" , "CH" , "P20" , "20th percentile" );
  add_var( "STATS" , "CH" , "P30" , "30th percentile" );
  add_var( "STATS" , "CH" , "P40" , "40th percentile" );
  add_var( "STATS" , "CH" , "P50" , "50th percentile" );
  add_var( "STATS" , "CH" , "P60" , "60th percentile" );
  add_var( "STATS" , "CH" , "P70" , "70th percentile" );
  add_var( "STATS" , "CH" , "P80" , "80th percentile" );
  add_var( "STATS" , "CH" , "P90" , "90th percentile" );
  add_var( "STATS" , "CH" , "P95" , "95th percentile" );
  add_var( "STATS" , "CH" , "P98" , "98th percentile" );
  add_var( "STATS" , "CH" , "P99" , "99th percentile" );


  add_var( "STATS" , "CH" , "MAX_ENCODING" , "Possible # of unique values" );
  add_var( "STATS" , "CH" , "OBS_ENCODING" , "Observed # of unique values" );
  add_var( "STATS" , "CH" , "PCT_ENCODING" , "Obs/possible unique values" );

  add_table( "STATS" , "CH,VAL" , "Encoding value distribution [encoding]" );
  add_var( "STATS" , "CH,VAL" , "CNT" , "Number of observations" );

  add_var( "STATS" , "CH" , "NE" ,  "Total number of epochs in record [epoch]" );
  add_var( "STATS" , "CH" , "NE1" , "Number of unmasked epochs actually used in calculations [epoch]" );

  add_table( "STATS" , "CH,E" , "Per-epoch, per-channel statistics for unmasked epochs only" );
  add_var( "STATS" , "CH,E" , "N" , "Number of samples in the epoch" );
  add_var( "STATS" , "CH,E" , "MIN" , "Signal minimum (from data, not EDF header)" );
  add_var( "STATS" , "CH,E" , "MAX" , "Signal maximum (from data, not EDF header)" );
  add_var( "STATS" , "CH,E" , "MEAN" , "Signal mean" );
  add_var( "STATS" , "CH,E" , "MEDIAN" , "Signal median" );
  add_var( "STATS" , "CH,E" , "RMS" , "Signal root mean square" );
  add_var( "STATS" , "CH,E" , "SD" , "Signal SD" );
  add_var( "STATS" , "CH,E" , "SKEW" , "Signal skewness" );
  add_var( "STATS" , "CH,E" , "KURT" , "Signal kurtosis" );

  add_var( "STATS" , "CH,E" , "P01" , "1st percentile" );
  add_var( "STATS" , "CH,E" , "P02" , "2nd percentile" );
  add_var( "STATS" , "CH,E" , "P05" , "5th percentile" );
  add_var( "STATS" , "CH,E" , "P10" , "10th percentile" );
  add_var( "STATS" , "CH,E" , "P20" , "20th percentile" );
  add_var( "STATS" , "CH,E" , "P30" , "30th percentile" );
  add_var( "STATS" , "CH,E" , "P40" , "40th percentile" );
  add_var( "STATS" , "CH,E" , "P50" , "50th percentile" );
  add_var( "STATS" , "CH,E" , "P60" , "60th percentile" );
  add_var( "STATS" , "CH,E" , "P70" , "70th percentile" );
  add_var( "STATS" , "CH,E" , "P80" , "80th percentile" );
  add_var( "STATS" , "CH,E" , "P90" , "90th percentile" );
  add_var( "STATS" , "CH,E" , "P95" , "95th percentile" );
  add_var( "STATS" , "CH,E" , "P98" , "98th percentile" );
  add_var( "STATS" , "CH,E" , "P99" , "99th percentile" );

  //
  // SIGSTATS
  //

  add_cmd( "summ" , "SIGSTATS" , "Epoch-level signal statistics (Hjorth, RMS, clipping, catch22)" );
  add_url( "SIGSTATS" , "summaries/#sigstats" );
  add_verb( "SIGSTATS" ,
	    "SIGSTATS computes signal-quality and complexity summaries, centered on Hjorth\n"
	    "statistics and related epoch-level metrics.\n"
	    "\n"
	    "Depending on options, it can report clipping, flatness, RMS, max-threshold counts,\n"
	    "Petrosian fractal dimension, permutation entropy, catch22/catch24 features, and\n"
	    "second-order or epoch-series summaries of those metrics." );
  add_param( "SIGSTATS" , "sig" , "C3,C4" , "Restrict analysis to these channels" );
  add_param( "SIGSTATS" , "verbose" , "" , "Report epoch-level statistics" );
  add_param( "SIGSTATS" , "epoch" , "" , "Report epoch-level statistics (same as verbose)" );
  add_param( "SIGSTATS" , "rms" , "" , "Report RMS values" );
  add_param( "SIGSTATS" , "clipped" , "0.05" , "Report proportion of clipped sample points" );
  add_param( "SIGSTATS" , "flat" , "1e-6" , "Report proportion of flat sample points using this epsilon" );
  add_param( "SIGSTATS" , "max" , "200" , "Report proportion of sample points above this absolute value" );
  add_param( "SIGSTATS" , "sr-over" , "128" , "Restrict analysis to channels at or above this sampling rate" );
  add_param( "SIGSTATS" , "catch22" , "" , "Calculate catch22 features" );
  add_param( "SIGSTATS" , "catch24" , "" , "Calculate catch22 plus mean and SD" );
  add_param( "SIGSTATS" , "epoch-catch22" , "" , "Calculate catch22 on epoch-to-epoch statistic series" );
  add_param( "SIGSTATS" , "epoch-catch24" , "" , "Calculate catch24 on epoch-to-epoch statistic series" );
  add_param( "SIGSTATS" , "pfd" , "" , "Calculate Petrosian fractal dimension" );
  add_param( "SIGSTATS" , "dynam" , "" , "Run quantile dynamics on epoch-level Hjorth statistics" );
  add_param( "SIGSTATS" , "pe" , "" , "Calculate permutation entropy for default embedding dimensions" );
  add_param( "SIGSTATS" , "pe-m" , "3,4,5,6,7" , "Embedding dimensions for permutation entropy" );
  add_param( "SIGSTATS" , "pe-t" , "1" , "Lag for permutation entropy" );
  add_param( "SIGSTATS" , "hjorth2" , "" , "Calculate second-order Hjorth statistics" );
  add_param( "SIGSTATS" , "hjorth2-win" , "1" , "Window size in seconds for second-order Hjorth statistics" );
  add_param( "SIGSTATS" , "hjorth2-inc" , "1" , "Window increment in seconds for second-order Hjorth statistics" );

  add_table( "SIGSTATS" , "CH" , "Per-channel whole-signal statistics" );
  add_var( "SIGSTATS" , "CH" , "H1" , "First Hjorth parameter (activity)" );
  add_var( "SIGSTATS" , "CH" , "H2" , "Second Hjorth parameter (mobility)" );
  add_var( "SIGSTATS" , "CH" , "H3" , "Third Hjorth parameter (complexity)" );
  add_var( "SIGSTATS" , "CH" , "CLIP" , "Proportion of clipped sample points" );
  add_var( "SIGSTATS" , "CH" , "FLAT" , "Proportion of flat sample points" );
  add_var( "SIGSTATS" , "CH" , "MAX" , "Proportion of sample points exceeding the max threshold" );
  add_var( "SIGSTATS" , "CH" , "RMS" , "Signal root mean square" );
  add_var( "SIGSTATS" , "CH" , "EC22_NC" , "Number of contiguous epoch blocks used for epoch-catch22 summaries" );
  add_var( "SIGSTATS" , "CH" , "EC22_NE" , "Total number of epochs contributing to epoch-catch22 summaries" );

  add_var( "SIGSTATS" , "CH" , "acf_first_min" , "First minimum of the autocorrelation function" );
  add_var( "SIGSTATS" , "CH" , "acf_timescale" , "First 1/e crossing of the autocorrelation function" );
  add_var( "SIGSTATS" , "CH" , "ami2" , "Histogram-based automutual information at lag 2" );
  add_var( "SIGSTATS" , "CH" , "ami_timescale" , "First minimum of the automutual information function" );
  add_var( "SIGSTATS" , "CH" , "centroid_freq" , "Welch spectral centroid frequency" );
  add_var( "SIGSTATS" , "CH" , "dfa" , "Detrended fluctuation analysis scaling statistic" );
  add_var( "SIGSTATS" , "CH" , "embedding_dist" , "Goodness of exponential fit to embedding distance distribution" );
  add_var( "SIGSTATS" , "CH" , "entropy_pairs" , "Entropy of successive pairs in a symbolized series" );
  add_var( "SIGSTATS" , "CH" , "forecast_error" , "Error of a 3-point rolling mean forecast" );
  add_var( "SIGSTATS" , "CH" , "high_fluctuation" , "Proportion of high incremental changes in the series" );
  add_var( "SIGSTATS" , "CH" , "low_freq_power" , "Power in the lowest 20% of frequencies" );
  add_var( "SIGSTATS" , "CH" , "mean" , "Series mean" );
  add_var( "SIGSTATS" , "CH" , "mode_10" , "10-bin histogram mode" );
  add_var( "SIGSTATS" , "CH" , "mode_5" , "5-bin histogram mode" );
  add_var( "SIGSTATS" , "CH" , "outlier_timing_neg" , "Timing statistic for negative outliers" );
  add_var( "SIGSTATS" , "CH" , "outlier_timing_pos" , "Timing statistic for positive outliers" );
  add_var( "SIGSTATS" , "CH" , "periodicity" , "Wang periodicity metric" );
  add_var( "SIGSTATS" , "CH" , "rs_range" , "Rescaled-range fluctuation analysis scaling statistic" );
  add_var( "SIGSTATS" , "CH" , "SD" , "Series standard deviation" );
  add_var( "SIGSTATS" , "CH" , "stretch_decreasing" , "Longest stretch of consecutive decreases" );
  add_var( "SIGSTATS" , "CH" , "stretch_high" , "Longest stretch of above-mean values" );
  add_var( "SIGSTATS" , "CH" , "transition_matrix" , "Transition matrix column variance in a symbolized series" );
  add_var( "SIGSTATS" , "CH" , "trev" , "Time-reversibility statistic" );
  add_var( "SIGSTATS" , "CH" , "whiten_timescale" , "Change in autocorrelation timescale after differencing" );

  add_table( "SIGSTATS" , "CH,EC22" , "Epoch-series catch22 summaries" );
  add_var( "SIGSTATS" , "CH,EC22" , "H1" , "Catch22 summary for the epoch-wise H1 series" );
  add_var( "SIGSTATS" , "CH,EC22" , "H2" , "Catch22 summary for the epoch-wise H2 series" );
  add_var( "SIGSTATS" , "CH,EC22" , "H3" , "Catch22 summary for the epoch-wise H3 series" );
  add_var( "SIGSTATS" , "CH,EC22" , "acf_first_min" , "Catch22 summary for the epoch-wise acf_first_min series" );
  add_var( "SIGSTATS" , "CH,EC22" , "acf_timescale" , "Catch22 summary for the epoch-wise acf_timescale series" );
  add_var( "SIGSTATS" , "CH,EC22" , "ami2" , "Catch22 summary for the epoch-wise ami2 series" );
  add_var( "SIGSTATS" , "CH,EC22" , "ami_timescale" , "Catch22 summary for the epoch-wise ami_timescale series" );
  add_var( "SIGSTATS" , "CH,EC22" , "centroid_freq" , "Catch22 summary for the epoch-wise centroid_freq series" );
  add_var( "SIGSTATS" , "CH,EC22" , "dfa" , "Catch22 summary for the epoch-wise dfa series" );
  add_var( "SIGSTATS" , "CH,EC22" , "embedding_dist" , "Catch22 summary for the epoch-wise embedding_dist series" );
  add_var( "SIGSTATS" , "CH,EC22" , "entropy_pairs" , "Catch22 summary for the epoch-wise entropy_pairs series" );
  add_var( "SIGSTATS" , "CH,EC22" , "forecast_error" , "Catch22 summary for the epoch-wise forecast_error series" );
  add_var( "SIGSTATS" , "CH,EC22" , "high_fluctuation" , "Catch22 summary for the epoch-wise high_fluctuation series" );
  add_var( "SIGSTATS" , "CH,EC22" , "low_freq_power" , "Catch22 summary for the epoch-wise low_freq_power series" );
  add_var( "SIGSTATS" , "CH,EC22" , "mean" , "Catch22 summary for the epoch-wise mean series" );
  add_var( "SIGSTATS" , "CH,EC22" , "mode_10" , "Catch22 summary for the epoch-wise mode_10 series" );
  add_var( "SIGSTATS" , "CH,EC22" , "mode_5" , "Catch22 summary for the epoch-wise mode_5 series" );
  add_var( "SIGSTATS" , "CH,EC22" , "outlier_timing_neg" , "Catch22 summary for the epoch-wise outlier_timing_neg series" );
  add_var( "SIGSTATS" , "CH,EC22" , "outlier_timing_pos" , "Catch22 summary for the epoch-wise outlier_timing_pos series" );
  add_var( "SIGSTATS" , "CH,EC22" , "periodicity" , "Catch22 summary for the epoch-wise periodicity series" );
  add_var( "SIGSTATS" , "CH,EC22" , "rs_range" , "Catch22 summary for the epoch-wise rs_range series" );
  add_var( "SIGSTATS" , "CH,EC22" , "SD" , "Catch22 summary for the epoch-wise SD series" );
  add_var( "SIGSTATS" , "CH,EC22" , "stretch_decreasing" , "Catch22 summary for the epoch-wise stretch_decreasing series" );
  add_var( "SIGSTATS" , "CH,EC22" , "stretch_high" , "Catch22 summary for the epoch-wise stretch_high series" );
  add_var( "SIGSTATS" , "CH,EC22" , "transition_matrix" , "Catch22 summary for the epoch-wise transition_matrix series" );
  add_var( "SIGSTATS" , "CH,EC22" , "trev" , "Catch22 summary for the epoch-wise trev series" );
  add_var( "SIGSTATS" , "CH,EC22" , "whiten_timescale" , "Catch22 summary for the epoch-wise whiten_timescale series" );

  add_table( "SIGSTATS" , "CH,E" , "Per-channel per-epoch statistics [epoch]" );
  add_var( "SIGSTATS" , "CH,E" , "H1" , "First Hjorth parameter (activity)" );
  add_var( "SIGSTATS" , "CH,E" , "H2" , "Second Hjorth parameter (mobility)" );
  add_var( "SIGSTATS" , "CH,E" , "H3" , "Third Hjorth parameter (complexity)" );
  add_var( "SIGSTATS" , "CH,E" , "H1H1" , "Second-order Hjorth summary for H1 activity" );
  add_var( "SIGSTATS" , "CH,E" , "H1H2" , "Second-order Hjorth summary for H1 mobility" );
  add_var( "SIGSTATS" , "CH,E" , "H1H3" , "Second-order Hjorth summary for H1 complexity" );
  add_var( "SIGSTATS" , "CH,E" , "H2H1" , "Second-order Hjorth summary for H2 activity" );
  add_var( "SIGSTATS" , "CH,E" , "H2H2" , "Second-order Hjorth summary for H2 mobility" );
  add_var( "SIGSTATS" , "CH,E" , "H2H3" , "Second-order Hjorth summary for H2 complexity" );
  add_var( "SIGSTATS" , "CH,E" , "H3H1" , "Second-order Hjorth summary for H3 activity" );
  add_var( "SIGSTATS" , "CH,E" , "H3H2" , "Second-order Hjorth summary for H3 mobility" );
  add_var( "SIGSTATS" , "CH,E" , "H3H3" , "Second-order Hjorth summary for H3 complexity" );
  add_var( "SIGSTATS" , "CH,E" , "CLIP" , "Proportion of clipped sample points" );
  add_var( "SIGSTATS" , "CH,E" , "FLAT" , "Proportion of flat sample points" );
  add_var( "SIGSTATS" , "CH,E" , "MAX" , "Proportion of sample points exceeding the max threshold" );
  add_var( "SIGSTATS" , "CH,E" , "RMS" , "Signal root mean square" );
  add_var( "SIGSTATS" , "CH,E" , "PFD" , "Petrosian fractal dimension" );
  add_var( "SIGSTATS" , "CH,E" , "PE3" , "Permutation entropy with embedding dimension 3" );
  add_var( "SIGSTATS" , "CH,E" , "PE4" , "Permutation entropy with embedding dimension 4" );
  add_var( "SIGSTATS" , "CH,E" , "PE5" , "Permutation entropy with embedding dimension 5" );
  add_var( "SIGSTATS" , "CH,E" , "PE6" , "Permutation entropy with embedding dimension 6" );
  add_var( "SIGSTATS" , "CH,E" , "PE7" , "Permutation entropy with embedding dimension 7" );
  add_var( "SIGSTATS" , "CH,E" , "acf_first_min" , "First minimum of the autocorrelation function" );
  add_var( "SIGSTATS" , "CH,E" , "acf_timescale" , "First 1/e crossing of the autocorrelation function" );
  add_var( "SIGSTATS" , "CH,E" , "ami2" , "Histogram-based automutual information at lag 2" );
  add_var( "SIGSTATS" , "CH,E" , "ami_timescale" , "First minimum of the automutual information function" );
  add_var( "SIGSTATS" , "CH,E" , "centroid_freq" , "Welch spectral centroid frequency" );
  add_var( "SIGSTATS" , "CH,E" , "dfa" , "Detrended fluctuation analysis scaling statistic" );
  add_var( "SIGSTATS" , "CH,E" , "embedding_dist" , "Goodness of exponential fit to embedding distance distribution" );
  add_var( "SIGSTATS" , "CH,E" , "entropy_pairs" , "Entropy of successive pairs in a symbolized series" );
  add_var( "SIGSTATS" , "CH,E" , "forecast_error" , "Error of a 3-point rolling mean forecast" );
  add_var( "SIGSTATS" , "CH,E" , "high_fluctuation" , "Proportion of high incremental changes in the series" );
  add_var( "SIGSTATS" , "CH,E" , "low_freq_power" , "Power in the lowest 20% of frequencies" );
  add_var( "SIGSTATS" , "CH,E" , "mean" , "Series mean" );
  add_var( "SIGSTATS" , "CH,E" , "mode_10" , "10-bin histogram mode" );
  add_var( "SIGSTATS" , "CH,E" , "mode_5" , "5-bin histogram mode" );
  add_var( "SIGSTATS" , "CH,E" , "outlier_timing_neg" , "Timing statistic for negative outliers" );
  add_var( "SIGSTATS" , "CH,E" , "outlier_timing_pos" , "Timing statistic for positive outliers" );
  add_var( "SIGSTATS" , "CH,E" , "periodicity" , "Wang periodicity metric" );
  add_var( "SIGSTATS" , "CH,E" , "rs_range" , "Rescaled-range fluctuation analysis scaling statistic" );
  add_var( "SIGSTATS" , "CH,E" , "SD" , "Series standard deviation" );
  add_var( "SIGSTATS" , "CH,E" , "stretch_decreasing" , "Longest stretch of consecutive decreases" );
  add_var( "SIGSTATS" , "CH,E" , "stretch_high" , "Longest stretch of above-mean values" );
  add_var( "SIGSTATS" , "CH,E" , "transition_matrix" , "Transition matrix column variance in a symbolized series" );
  add_var( "SIGSTATS" , "CH,E" , "trev" , "Time-reversibility statistic" );
  add_var( "SIGSTATS" , "CH,E" , "whiten_timescale" , "Change in autocorrelation timescale after differencing" );

  //
  // TABULATE
  //

  add_cmd( "summ" , "TABULATE" , "Tabulates discrete values in a signal" );
  add_url( "TABULATE" , "summaries/#tabulate" );
  add_verb( "TABULATE" ,
	    "TABULATE counts distinct sample values in one or more signals.\n"
	    "\n"
	    "It is mainly intended for discrete or quantized signals, where you want to know\n"
	    "how many values occur and how often each value is observed." );
  add_param( "TABULATE" , "sig" , "C3,C4" , "Restrict analysis to these channels" );
  add_param( "TABULATE" , "req" , "1000", "Count distinct values w/ at least this many observations" );
  add_param( "TABULATE" , "prec" , "2" , "Round values to this number of decimal places before tabulation" );

  add_table( "TABULATE" , "CH" , "Per-channel tabulation" );
  add_var( "TABULATE" , "CH" , "NV" , "Number of discrete values observed for this channel" );

  add_table( "TABULATE" , "CH,REQ" , "Per-channel required-count summaries" );
  add_var( "TABULATE" , "CH,REQ" , "NV" , "Number of values observed at least REQ times" );

  add_table( "TABULATE" , "CH,VALUE" , "Per-channel/value tabulation statistics" );
  add_var( "TABULATE" , "CH,VALUE" , "N" , "Number of sample points for this value/channel" );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // ANNOTATIONS
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // --xml
  //

  add_cmd( "helpers" , "--xml" , "View NSRR XML annotations" );
  add_verb( "--xml" ,
            "--xml reads a single NSRR-style XML annotation file and prints a "
            "compact console summary of the parsed annotations. It is a simple "
            "inspection utility and does not require an EDF or sample list." );
  add_param( "--xml" , "{xml}" , "f1.xml" , "A single XML file" );

  //
  // --xml2
  //

  add_cmd( "helpers" , "--xml2" , "View NSRR XML annotations verbosely" );
  add_verb( "--xml2" ,
            "--xml2 is the verbose counterpart to --xml. It dumps the full parsed "
            "XML tree to the console so the raw structure and nested fields can be "
            "inspected directly." );
  add_param( "--xml2" , "{xml}" , "f1.xml" , "A single XML file" );

  //
  // S2A
  //

  add_cmd( "annot" , "S2A" , "Signal-to-annotation" );
  add_url( "S2A" , "annotations/#s2a" );
  add_verb( "S2A" ,
	    "S2A creates annotations from signal values.\n"
	    "\n"
	    "Depending on options, it can threshold or otherwise convert a signal into\n"
	    "interval annotations, and in waveform mode it also reports per-channel QC\n"
	    "counts for detected waves." );
  add_param( "S2A" , "sig" , "C3,C4" , "Signal(s) to analyze" );
  add_param( "S2A" , "encoding" , "N2,0.5,+0.1" , "Encoding mode: label,value,window tuples" );
  add_param( "S2A" , "encoding2" , "N2,0.5,1.0" , "Encoding mode: label,lower,upper tuples" );
  add_param( "S2A" , "bins" , "0,360,12" , "For generic mode: bins=min,max,n; for wave mode: request phase-bin annotations" );
  add_param( "S2A" , "q" , "5" , "Create q quantile-based annotation bins" );
  add_param( "S2A" , "pos" , "75" , "Create annotations for values above this threshold" );
  add_param( "S2A" , "neg" , "-75" , "Create annotations for values below this threshold" );
  add_param( "S2A" , "pos-pct" , "0.95" , "Create annotations for values above this percentile" );
  add_param( "S2A" , "neg-pct" , "0.05" , "Create annotations for values below this percentile" );
  add_param( "S2A" , "class" , "STATE" , "Use one annotation class, with labels stored as instance IDs" );
  add_param( "S2A" , "add-channel-label" , "" , "Append channel labels to created annotation class names" );
  add_param( "S2A" , "bin-label" , "B" , "Base label for generated bins" );
  add_param( "S2A" , "no-bin-label" , "" , "Do not prepend an automatic bin label" );
  add_param( "S2A" , "span-gaps" , "" , "Allow annotations to span EDF discontinuities" );
  add_param( "S2A" , "waves" , "WAVE" , "Wave mode: create wave annotations from a phase-angle signal" );
  add_param( "S2A" , "phase-ext" , "_ht_ang" , "Suffix used to locate the phase-angle signal in wave mode" );
  add_param( "S2A" , "mag-sig" , "_ht_mag" , "Suffix used to locate the magnitude signal in wave mode" );
  add_param( "S2A" , "add-channel-inst-label" , "" , "Add channel names to wave-mode instance IDs" );
  add_param( "S2A" , "add-channel-class-label" , "" , "Add channel names to wave-mode class labels" );
  add_param( "S2A" , "pos2neg" , "" , "Wave mode: segment positive-to-negative crossings" );
  add_param( "S2A" , "t-min" , "0.3" , "Wave mode: minimum duration in seconds" );
  add_param( "S2A" , "t-max" , "2.0" , "Wave mode: maximum duration in seconds" );
  add_param( "S2A" , "t-min-phbin" , "0.01" , "Wave mode: minimum per-phase-bin duration in seconds" );
  add_param( "S2A" , "t-max-phbin" , "0.25" , "Wave mode: maximum per-phase-bin duration in seconds" );
  add_param( "S2A" , "mag-percentile" , "0.8" , "Wave mode: keep only waves above this magnitude percentile" );
  add_param( "S2A" , "mag-z" , "0.5" , "Wave mode: keep only waves above this magnitude z-score" );
  add_param( "S2A" , "monotonic" , "" , "Wave mode: require monotonic phase-bin progression" );
  add_param( "S2A" , "slope" , "" , "Wave mode: add slope annotations for monotonic waves" );
  add_param( "S2A" , "state" , "" , "Wave mode: add state annotations for monotonic waves" );

  add_table( "S2A" , "CH" , "Channel-level metrics (for waves option)" );
  add_var( "S2A" , "CH" , "N"  , "Final number of included waves" );
  add_var( "S2A" , "CH" , "N0" , "Original (pre-QC) number of included waves" );
  add_var( "S2A" , "CH" , "EXC1_DUR"  , "Number of exclusions due to duration criteria" );
  add_var( "S2A" , "CH" , "EXC2_MONO" , "Number of exclusions due to monotonic phase constraint" );
  add_var( "S2A" , "CH" , "EXC3_MAG" , "Number of exclusions due to magnitude criteria" );
  add_var( "S2A" , "CH" , "EXC4_PDUR"  , "Number of exclusions due to phase-bin duration criteria" );

  //
  // S2C
  //

  add_cmd( "annot" , "S2C" , "Signal-to-cycle" );
  add_url( "S2C" , "annotations/#s2c" );
  add_verb( "S2C" ,
	    "S2C segments oscillatory activity into cycle annotations and optional cycle-locked summaries.\n"
	    "\n"
	    "It is designed for extracting cycles, half-waves, peak points and aligned signal summaries,\n"
	    "with outputs ranging from seed-level summaries to per-cycle and phase/time-binned metrics." );

  add_table( "S2C" , "SEED" , "Seed-level summary statistics" );
  add_var( "S2C" , "SEED" , "N_PRE" , "Number of detected cycles before magnitude filtering" );
  add_var( "S2C" , "SEED" , "N_POST" , "Number of detected cycles after magnitude filtering" );
  add_var( "S2C" , "SEED" , "DUR_POS" , "Mean positive half-wave duration (seconds)" );
  add_var( "S2C" , "SEED" , "DUR_POS_MD" , "Median positive half-wave duration (seconds)" );
  add_var( "S2C" , "SEED" , "DUR_POS_MAD" , "MAD of positive half-wave duration (seconds) [emit-mad]" );
  add_var( "S2C" , "SEED" , "DUR_NEG" , "Mean negative half-wave duration (seconds)" );
  add_var( "S2C" , "SEED" , "DUR_NEG_MD" , "Median negative half-wave duration (seconds)" );
  add_var( "S2C" , "SEED" , "DUR_NEG_MAD" , "MAD of negative half-wave duration (seconds) [emit-mad]" );
  add_var( "S2C" , "SEED" , "DUR" , "Mean total cycle duration (seconds)" );
  add_var( "S2C" , "SEED" , "DUR_MD" , "Median total cycle duration (seconds)" );
  add_var( "S2C" , "SEED" , "DUR_MAD" , "MAD of total cycle duration (seconds) [emit-mad]" );
  add_var( "S2C" , "SEED" , "DUR_ZC_NEG" , "Mean duration from zero-crossing to negative peak (seconds)" );
  add_var( "S2C" , "SEED" , "DUR_ZC_NEG_MD" , "Median duration from zero-crossing to negative peak (seconds)" );
  add_var( "S2C" , "SEED" , "DUR_ZC_NEG_MAD" , "MAD of zero-crossing to negative peak duration [emit-mad]" );
  add_var( "S2C" , "SEED" , "DUR_NEG_ZC" , "Mean duration from negative peak to next zero-crossing (seconds)" );
  add_var( "S2C" , "SEED" , "DUR_NEG_ZC_MD" , "Median duration from negative peak to next zero-crossing (seconds)" );
  add_var( "S2C" , "SEED" , "DUR_NEG_ZC_MAD" , "MAD of negative peak to zero-crossing duration [emit-mad]" );
  add_var( "S2C" , "SEED" , "DUR_ZC_POS" , "Mean duration from zero-crossing to positive peak (seconds)" );
  add_var( "S2C" , "SEED" , "DUR_ZC_POS_MD" , "Median duration from zero-crossing to positive peak (seconds)" );
  add_var( "S2C" , "SEED" , "DUR_ZC_POS_MAD" , "MAD of zero-crossing to positive peak duration [emit-mad]" );
  add_var( "S2C" , "SEED" , "DUR_POS_ZC" , "Mean duration from positive peak to next zero-crossing (seconds)" );
  add_var( "S2C" , "SEED" , "DUR_POS_ZC_MD" , "Median duration from positive peak to next zero-crossing (seconds)" );
  add_var( "S2C" , "SEED" , "DUR_POS_ZC_MAD" , "MAD of positive peak to zero-crossing duration [emit-mad]" );
  add_var( "S2C" , "SEED" , "AMP_POS" , "Mean positive-peak amplitude (from baseline)" );
  add_var( "S2C" , "SEED" , "AMP_POS_MD" , "Median positive-peak amplitude" );
  add_var( "S2C" , "SEED" , "AMP_POS_MAD" , "MAD of positive-peak amplitude [emit-mad]" );
  add_var( "S2C" , "SEED" , "AMP_NEG" , "Mean absolute negative-peak amplitude (from baseline)" );
  add_var( "S2C" , "SEED" , "AMP_NEG_MD" , "Median absolute negative-peak amplitude" );
  add_var( "S2C" , "SEED" , "AMP_NEG_MAD" , "MAD of absolute negative-peak amplitude [emit-mad]" );
  add_var( "S2C" , "SEED" , "AMP_P2P" , "Mean per-cycle peak-to-peak amplitude (v_pos - v_neg)" );
  add_var( "S2C" , "SEED" , "AMP_P2P_MD" , "Median per-cycle peak-to-peak amplitude" );
  add_var( "S2C" , "SEED" , "AMP_P2P_MAD" , "MAD of per-cycle peak-to-peak amplitude [emit-mad]" );
  add_var( "S2C" , "SEED" , "SLOPE_POS" , "Mean absolute positive slope" );
  add_var( "S2C" , "SEED" , "SLOPE_POS_MD" , "Median absolute positive slope" );
  add_var( "S2C" , "SEED" , "SLOPE_POS_MAD" , "MAD of absolute positive slope [emit-mad]" );
  add_var( "S2C" , "SEED" , "SLOPE_NEG" , "Mean absolute negative slope" );
  add_var( "S2C" , "SEED" , "SLOPE_NEG_MD" , "Median absolute negative slope" );
  add_var( "S2C" , "SEED" , "SLOPE_NEG_MAD" , "MAD of absolute negative slope [emit-mad]" );
  add_var( "S2C" , "SEED" , "SHARP_POS" , "Mean positive-peak sharpness" );
  add_var( "S2C" , "SEED" , "SHARP_POS_MD" , "Median positive-peak sharpness" );
  add_var( "S2C" , "SEED" , "SHARP_POS_MAD" , "MAD of positive-peak sharpness [emit-mad]" );
  add_var( "S2C" , "SEED" , "SHARP_NEG" , "Mean negative-peak sharpness" );
  add_var( "S2C" , "SEED" , "SHARP_NEG_MD" , "Median negative-peak sharpness" );
  add_var( "S2C" , "SEED" , "SHARP_NEG_MAD" , "MAD of negative-peak sharpness [emit-mad]" );
  add_var( "S2C" , "SEED" , "DUR_RATIO" , "Mean ratio of pos/neg half-wave duration [emit-asym]" );
  add_var( "S2C" , "SEED" , "DUR_RATIO_MD" , "Median ratio of pos/neg half-wave duration [emit-asym]" );
  add_var( "S2C" , "SEED" , "DUR_RATIO_MAD" , "MAD of pos/neg half-wave duration ratio [emit-asym,emit-mad]" );
  add_var( "S2C" , "SEED" , "AMP_ASYM" , "Mean amplitude asymmetry [emit-asym]" );
  add_var( "S2C" , "SEED" , "AMP_ASYM_MD" , "Median amplitude asymmetry [emit-asym]" );
  add_var( "S2C" , "SEED" , "AMP_ASYM_MAD" , "MAD of amplitude asymmetry [emit-asym,emit-mad]" );
  add_var( "S2C" , "SEED" , "SLOPE_ASYM" , "Mean slope asymmetry [emit-asym]" );
  add_var( "S2C" , "SEED" , "SLOPE_ASYM_MD" , "Median slope asymmetry [emit-asym]" );
  add_var( "S2C" , "SEED" , "SLOPE_ASYM_MAD" , "MAD of slope asymmetry [emit-asym,emit-mad]" );
  add_var( "S2C" , "SEED" , "SHARP_ASYM" , "Mean sharpness asymmetry [emit-asym]" );
  add_var( "S2C" , "SEED" , "SHARP_ASYM_MD" , "Median sharpness asymmetry [emit-asym]" );
  add_var( "S2C" , "SEED" , "SHARP_ASYM_MAD" , "MAD of sharpness asymmetry [emit-asym,emit-mad]" );

  add_table( "S2C" , "SEED,BIN" , "12-bin cycle-locked seed summary (B01..B12)" );
  add_var( "S2C" , "SEED,BIN" , "MEAN" , "Mean seed value in cycle bin" );

  add_table( "S2C" , "SEED,CYCLE" , "Per-cycle metrics (optional)" );
  add_var( "S2C" , "SEED,CYCLE" , "T0_S" , "Cycle start time (seconds)" );
  add_var( "S2C" , "SEED,CYCLE" , "T1_S" , "Cycle end time (seconds)" );
  add_var( "S2C" , "SEED,CYCLE" , "DUR_POS" , "Positive half-wave duration (seconds)" );
  add_var( "S2C" , "SEED,CYCLE" , "DUR_NEG" , "Negative half-wave duration (seconds)" );
  add_var( "S2C" , "SEED,CYCLE" , "DUR_ZC_NEG" , "Duration from zero-crossing to negative peak (seconds)" );
  add_var( "S2C" , "SEED,CYCLE" , "DUR_NEG_ZC" , "Duration from negative peak to next zero-crossing (seconds)" );
  add_var( "S2C" , "SEED,CYCLE" , "DUR_ZC_POS" , "Duration from zero-crossing to positive peak (seconds)" );
  add_var( "S2C" , "SEED,CYCLE" , "DUR_POS_ZC" , "Duration from positive peak to next zero-crossing (seconds)" );
  add_var( "S2C" , "SEED,CYCLE" , "AMP_POS" , "Positive-peak amplitude (from baseline)" );
  add_var( "S2C" , "SEED,CYCLE" , "AMP_NEG" , "Absolute negative-peak amplitude (from baseline)" );
  add_var( "S2C" , "SEED,CYCLE" , "AMP_P2P" , "Peak-to-peak amplitude (v_pos - v_neg)" );
  add_var( "S2C" , "SEED,CYCLE" , "SLOPE_POS" , "Absolute positive slope" );
  add_var( "S2C" , "SEED,CYCLE" , "SLOPE_NEG" , "Absolute negative slope" );
  add_var( "S2C" , "SEED,CYCLE" , "SHARP_POS" , "Positive-peak sharpness" );
  add_var( "S2C" , "SEED,CYCLE" , "SHARP_NEG" , "Negative-peak sharpness" );
  add_var( "S2C" , "SEED,CYCLE" , "DUR_RATIO" , "Ratio of pos/neg half-wave duration [emit-asym]" );
  add_var( "S2C" , "SEED,CYCLE" , "AMP_ASYM" , "Amplitude asymmetry [emit-asym]" );
  add_var( "S2C" , "SEED,CYCLE" , "SLOPE_ASYM" , "Slope asymmetry [emit-asym]" );
  add_var( "S2C" , "SEED,CYCLE" , "SHARP_ASYM" , "Sharpness asymmetry [emit-asym]" );

  add_table( "S2C" , "SEED,SIG" , "Seed-by-signal summary statistics" );
  add_var( "S2C" , "SEED,SIG" , "D" , "Mean waveform peak-to-peak amplitude" );
  add_var( "S2C" , "SEED,SIG" , "TAU_MAX_DEG" , "Phase (deg) of peak in mean waveform" );
  add_var( "S2C" , "SEED,SIG" , "SINFIT_OVER_P2P" , "Sin/cos fit amplitude divided by peak-to-peak amplitude" );
  add_var( "S2C" , "SEED,SIG" , "D_SD" , "SD of per-cycle peak-to-peak amplitude" );
  add_var( "S2C" , "SEED,SIG" , "CIRC_MEAN_DEG" , "Circular mean phase of per-cycle peaks (deg)" );
  add_var( "S2C" , "SEED,SIG" , "R" , "Resultant vector length of per-cycle peak phases" );
  add_var( "S2C" , "SEED,SIG" , "CC_LAG" , "Cross-correlation lag (seconds)" );
  add_var( "S2C" , "SEED,SIG" , "CC_R" , "Cross-correlation at lag" );

  add_var( "S2C" , "SEED,SIG" , "DT" , "Mean lag to nearest local extremum (seconds)" );
  add_var( "S2C" , "SEED,SIG" , "DT_MD" , "Median lag to nearest local extremum (seconds)" );
  add_var( "S2C" , "SEED,SIG" , "DT_MAD" , "MAD of lag to nearest local extremum [emit-mad]" );
  add_var( "S2C" , "SEED,SIG" , "TD_USED_PCT" , "Proportion of cycles with valid time-domain windows" );
  add_var( "S2C" , "SEED,SIG" , "TD_P2P" , "Mean window peak-to-peak amplitude" );
  add_var( "S2C" , "SEED,SIG" , "TD_P2P_MD" , "Median window peak-to-peak amplitude" );
  add_var( "S2C" , "SEED,SIG" , "TD_P2P_MAD" , "MAD of window peak-to-peak amplitude [emit-mad]" );

  add_table( "S2C" , "SEED,SIG,BIN" , "12-bin cycle-locked means (B01..B12)" );
  add_var( "S2C" , "SEED,SIG,BIN" , "MEAN" , "Mean signal value in cycle bin" );

  add_table( "S2C" , "SEED,SIG,PH" , "Phase-bin mean waveform" );
  add_var( "S2C" , "SEED,SIG,PH" , "MEAN" , "Mean amplitude in phase bin" );
  add_var( "S2C" , "SEED,SIG,PH" , "SE" , "Bootstrap SE of mean [emit-se]" );
  add_var( "S2C" , "SEED,SIG,PH" , "CI_LO" , "Bootstrap CI lower bound" );
  add_var( "S2C" , "SEED,SIG,PH" , "CI_HI" , "Bootstrap CI upper bound" );

  add_table( "S2C" , "SEED,SIG,PH,AMP" , "Phase-by-amplitude density" );
  add_var( "S2C" , "SEED,SIG,PH,AMP" , "DENS" , "Density in phase/amplitude bin" );

  add_table( "S2C" , "SEED,SIG,SEC" , "Time-domain grid (time-locked)" );
  add_var( "S2C" , "SEED,SIG,SEC" , "MEAN" , "Mean amplitude at time bin" );
  add_var( "S2C" , "SEED,SIG,SEC" , "SE" , "Bootstrap SE of mean [emit-se]" );
  add_var( "S2C" , "SEED,SIG,SEC" , "CI_LO" , "Bootstrap CI lower bound" );
  add_var( "S2C" , "SEED,SIG,SEC" , "CI_HI" , "Bootstrap CI upper bound" );

  add_table( "S2C" , "SEED,SIG,SEC,AMP" , "Time-by-amplitude density" );
  add_var( "S2C" , "SEED,SIG,SEC,AMP" , "DENS" , "Density in time/amplitude bin" );

  //
  // DROP-ANNOTS
  //

  add_cmd( "annot" , "DROP-ANNOTS" , "Drop annotations" );
  add_url( "DROP-ANNOTS" , "annotations/#drop-annots" );
  add_verb( "DROP-ANNOTS" ,
	    "DROP-ANNOTS removes one or more loaded annotation classes from the current EDF.\n"
	    "\n"
	    "Use it to discard unwanted annotations before downstream analysis or export." );
  add_param( "DROP-ANNOTS" , "annot" , "N4,M" , "Drop annotations 'N4' and 'M'" );

  //
  // ANNOTS
  //

  add_cmd( "annot" , "ANNOTS" , "List annotations" );
  add_url( "ANNOTS" , "annotations/#annots" );
  add_verb( "ANNOTS" ,
	    "ANNOTS lists loaded annotations either as class summaries, instance summaries or per-epoch overlaps.\n"
	    "\n"
	    "It can restrict output to annotations overlapping unmasked regions, or include masked annotations\n"
	    "with explicit mask-status fields for each instance." );
  //  add_note( "ANNOTS" , "Any formatted free text goes here,\n shown at end of verbose help link\n");

  add_param( "ANNOTS" , "annot" , "N2,N3" , "Restrict output to these annotations" );
  add_param( "ANNOTS" , "epoch" , "" , "Show epoch-level summaries" );
  add_param( "ANNOTS" , "show-masked" , "" , "Show masked annotations (default is not to do so)" );
  add_param( "ANNOTS" , "any" , "" , "Keep annotations that have any overlap with one or more unmasked epochs (default)" );
  add_param( "ANNOTS" , "all" , "" , "Only keep annotations that are completely within unmasked epochs" );
  add_param( "ANNOTS" , "start" , "" , "Keep annotations that start in an unmasked epoch" );

  add_table( "ANNOTS" , "ANNOT" , "Class-level annotation summary" );
  add_var( "ANNOTS" , "ANNOT" , "COUNT" , "Number of instances of that annotation class" );
  add_var( "ANNOTS" , "ANNOT" , "DUR" , "Combined duration (seconds) of all instances of that annotation class" );

  add_table( "ANNOTS" , "ANNOT,INST" , "Instance-level annotation summary" );
  add_var( "ANNOTS" , "ANNOT,INST" , "COUNT" , "Number of instances of that annotation class and instance ID" );
  add_var( "ANNOTS" , "ANNOT,INST" , "DUR" , "Combined duration (seconds) of all instances of that annotation class and instance ID" );

  add_table( "ANNOTS" , "ANNOT,INST,T" , "Instance-level annotation tabulation" );
  add_var( "ANNOTS" , "ANNOT,INST,T" , "START" , "Start time (secs) of this instance" );
  add_var( "ANNOTS" , "ANNOT,INST,T" , "STOP" , "Stop time (secs) of this instance" );
  add_var( "ANNOTS" , "ANNOT,INST,T" , "DUR" , "Annotation duration (secs)" );
  add_var( "ANNOTS" , "ANNOT,INST,T" , "VAL" , "The meta-data for this instance, if any exists (otherwise missing NA)" );
  add_var( "ANNOTS" , "ANNOT,INST,T" , "ALL_MASKED" , "Annotation interval is completely masked (1/0) [show-masked]" );
  add_var( "ANNOTS" , "ANNOT,INST,T" , "ALL_UNMASKED" , "Annotation interval is completely unmasked (1/0) [show-masked]" );
  add_var( "ANNOTS" , "ANNOT,INST,T" , "SOME_MASKED" , "Annotation interval overlaps any masked region (1/0) [show-masked]" );
  add_var( "ANNOTS" , "ANNOT,INST,T" , "SOME_UNMASKED" , "Annotation interval overlaps any unmasked region (1/0) [show-masked]" );
  add_var( "ANNOTS" , "ANNOT,INST,T" , "START_MASKED" , "Annotation start lies in a masked region (1/0) [show-masked]" );

  add_var( "ANNOTS" , "ANNOT,INST,T" , "CH" , "Any associated channel(s)" );
  add_var( "ANNOTS" , "ANNOT,INST,T" , "START_ELAPSED_HMS" , "Annotation start (elapsed hh:mm:ss)" );
  add_var( "ANNOTS" , "ANNOT,INST,T" , "START_HMS" , "Annotation start" );
  add_var( "ANNOTS" , "ANNOT,INST,T" , "STOP_ELAPSED_HMS" , "Annotation stop (elapsed hh:mm:ss)" );
  add_var( "ANNOTS" , "ANNOT,INST,T" , "STOP_HMS" , "Annotation stop" );

  add_table( "ANNOTS" , "E,INTERVAL,INST,CH" , "Per-epoch instance-level annotation tabulation" );
  add_var( "ANNOTS" , "E,INTERVAL,INST,CH" , "AMASK" , "Annotation instance is masked/excluded under the current overlap rule (1/0) [epoch]" );
  add_var( "ANNOTS" , "E,INTERVAL,INST,CH" , "EMASK" , "Epoch is masked/excluded (1/0) [epoch]" );

  //
  // AXA
  //

  add_cmd( "annot" , "AXA" , "Annotation cross-tabs" );
  add_url( "AXA" , "annotations/#axa" );
  add_verb( "AXA" ,
	    "AXA summarizes pairwise relationships between seed annotations and other annotations.\n"
	    "\n"
	    "It reports overlap, distance-to-nearest, counts and total spanned duration, optionally\n"
	    "stratified within channels or by annotation instance." );
  add_param( "AXA" , "annot" , "N2,N3,arousal,apnea" , "Annotations to query" );
  add_param( "AXA" , "by-instance" , "" , "Define distinct annotations using class/instance ID" );
  add_param( "AXA" , "match-instance" , "" , "Only compare annotations with matching instance IDs [by-instance]" );
  add_param( "AXA" , "flatten" , "" , "Flatten events within each annotation before comparison" );
  add_param( "AXA" , "within-channel" , "" , "Restrict comparisons to annotations within the same channel" );
  add_param( "AXA" , "start" , "" , "Use annotation starts as the anchor for nearest-distance comparisons" );
  add_param( "AXA" , "stop" , "" , "Use annotation stops as the anchor for nearest-distance comparisons" );
  add_param( "AXA" , "w" , "5" , "Time window in seconds for nearest-distance comparisons" );
  add_param( "AXA" , "truncate" , "" , "Truncate nearest distances to the search window instead of excluding them" );
  add_param( "AXA" , "trunc" , "" , "Alias for truncate" );
  add_param( "AXA" , "verbose" , "" , "Write additional event-level information to the console" );
  
  add_table( "AXA" , "SEED,ANNOT" , "AXA pairwise metrics" );
  add_var( "AXA" , "SEED,ANNOT" , "P" , "Mean proportion of each seed spanned by annotation" );
  //add_var( "AXA" , "SEED,ANNOT" , "P_MD" , "Median proportion of each seed spanned by annotation" );
  add_var( "AXA" , "SEED,ANNOT" , "T" , "Mean time (secs) spanned by annotation per seed event" );
  //add_var( "AXA" , "SEED,ANNOT" , "T_MD" , "Median time (secs) spanned by annotation per seed event" );
  add_var( "AXA" , "SEED,ANNOT" , "N" , "Mean number of spanning annotations per seed event" );
  //add_var( "AXA" , "SEED,ANNOT" , "N_MD" , "Median number of spanning annotations per seed event" );
  add_var( "AXA" , "SEED,ANNOT" , "A" , "Proportion of seeds with any spanning annotation" );
  add_var( "AXA" , "SEED,ANNOT" , "D" , "Time (distance) to nearest" );
  add_var( "AXA" , "SEED,ANNOT" , "DABS" , "|Time| (distance) to nearest" );
  add_var( "AXA" , "SEED,ANNOT" , "D_N" , "Number of seeds with a nearest event" );
  add_var( "AXA" , "SEED,ANNOT" , "TOT_N" , "Total number of (flattened) annotations spanning all seeds" );
  add_var( "AXA" , "SEED,ANNOT" , "TOT_T" , "Total duration (secs) of (flattened) annotations spanning all seeds" );

  add_table( "AXA" , "CH,SEED,ANNOT" , "AXA pairwise metrics, within-channel" );
  add_var( "AXA" , "CH,SEED,ANNOT" , "P" , "Mean proportion of each seed spanned by annotation" );
  //add_var( "AXA" , "CH,SEED,ANNOT" , "P_MD" , "Median proportion of each seed spanned by annotation" );
  add_var( "AXA" , "CH,SEED,ANNOT" , "T" , "Mean time (secs) spanned by annotation per seed event" );
  //add_var( "AXA" , "CH,SEED,ANNOT" , "T_MD" , "Median time (secs) spanned by annotation per seed event" );
  add_var( "AXA" , "CH,SEED,ANNOT" , "N" , "Mean number of spanning annotations per seed event" );
  //add_var( "AXA" , "CH,SEED,ANNOT" , "N_MD" , "Median number of spanning annotations per seed event" );
  add_var( "AXA" , "CH,SEED,ANNOT" , "A" , "Proportion of seeds with any spanning annotation" );
  add_var( "AXA" , "CH,SEED,ANNOT" , "D" , "Time (distance) to nearest" );
  add_var( "AXA" , "CH,SEED,ANNOT" , "DABS" , "|Time| (distance) to nearest" );
  add_var( "AXA" , "CH,SEED,ANNOT" , "D_N" , "Number of seeds with a nearest event" );
  add_var( "AXA" , "CH,SEED,ANNOT" , "TOT_N" , "Total number of (flattened) annotations spanning all seeds" );
  add_var( "AXA" , "CH,SEED,ANNOT" , "TOT_T" , "Total duration (secs) of (flattened) annotations spanning all seeds" );

  //
  // SPANNING
  //

  add_cmd( "annot" , "SPANNING" , "Annotation coverage" );
  add_url( "SPANNING" , "annotations/#spanning" );
  add_verb( "SPANNING" ,
	    "SPANNING measures how much of the recording is covered by one or more annotation classes.\n"
	    "\n"
	    "It reports total annotation duration, collapsed covered duration, uncovered duration and any\n"
	    "invalid annotations that extend beyond the EDF timeline." );
  add_param( "SPANNING" , "annot" , "N1,N2,N3,R,W" , "Spanning annotation group" );

  add_table( "SPANNING" , "N" , "Invalid annotations" );
  add_var( "SPANNING" , "N" , "ANNOT" , "Annotation class" );
  add_var( "SPANNING" , "N" , "INST" , "Annotation instance" );
  add_var( "SPANNING" , "N" , "START" , "Start (seconds)" );
  add_var( "SPANNING" , "N" , "STOP" , "Stop (seconds)" );
  add_var( "SPANNING" , "N" , "DUR" , "Duration beyond the EDF timeline (seconds)" );

  add_table( "SPANNING" , "" , "Spanning summary report" );
  add_var( "SPANNING" , "" , "REC_SEC" , "EDF recording duration (seconds)" );
  add_var( "SPANNING" , "" , "REC_HMS" , "EDF recording duration (hh:mm:ss)" );

  add_var( "SPANNING" , "" , "ANNOT_N" , "Number of annotations in group" );
  add_var( "SPANNING" , "" , "ANNOT_SEC" , "Total (potentially overlapping) annotation duration (secs)" );
  add_var( "SPANNING" , "" , "ANNOT_HMS" , "Total (potentially overlapping) annotation duration (hh:mm:ss)" );

  add_var( "SPANNING" , "" , "ANNOT_OVERLAP" , "Do any annotations in group overlap w/ one another (0/1)?" );

  add_var( "SPANNING" , "" , "INVALID_N" , "Number of annotations that over-extend EDF duration" );
  add_var( "SPANNING" , "" , "VALID_N" , "Number of valid annotations, ANNOT_N - INVALID_N" );
  add_var( "SPANNING" , "" , "NSEGS" , "Number of contiguous annotation blocks" );

  add_var( "SPANNING" , "" , "INVALID_SEC" , "Total duration of all annotation beyond EDF end" );

  add_var( "SPANNING" , "" , "SPANNED_PCT" , "% of EDF spanned by 1+ of these annotations" );
  add_var( "SPANNING" , "" , "SPANNED_SEC" , "Duration of EDF spanned by 1+ of these annotations (secs)" );
  add_var( "SPANNING" , "" , "SPANNED_HMS" , "Duration of EDF spanned by 1+ of these annotations (hh:mm:ss)" );

  add_var( "SPANNING" , "" , "UNSPANNED_PCT" , "% of EDF unspanned by 1+ of these annotations" );
  add_var( "SPANNING" , "" , "UNSPANNED_SEC" , "Duration of EDF unspanned by 1+ of these annotations (secs)" );
  add_var( "SPANNING" , "" , "UNSPANNED_HMS" , "Duration of EDF unspanned by 1+ of these annotations (hh:mm:ss)" );

  //
  // WRITE-ANNOTS
  //

  add_cmd( "annot" , "WRITE-ANNOTS" , "Write annotations" );
  add_url( "WRITE-ANNOTS" , "outputs/#write-annots" );
  add_verb( "WRITE-ANNOTS" ,
	    "WRITE-ANNOTS exports the current annotation set to disk.\n"
	    "\n"
	    "It can write XML output or Luna's .annot format, making it useful after creating,\n"
	    "remapping or otherwise editing annotations inside a command stream." );
  
  add_param( "WRITE-ANNOTS" , "file" , "out.annot" , "Write annotations to this file" );
  add_param( "WRITE-ANNOTS" , "annot-dir" , "annots/" , "Write one annotation file per EDF into this directory" );
  add_param( "WRITE-ANNOTS" , "annot" , "N2,N3" , "Restrict output to these annotations" );
  add_param( "WRITE-ANNOTS" , "prefix" , "sp_" , "Also include annotations matching these prefixes" );
  add_param( "WRITE-ANNOTS" , "xml" , "" , "Force XML output format" );
  add_param( "WRITE-ANNOTS" , "hms" , "" , "Use clock times rather than elapsed seconds in .annot output" );
  add_param( "WRITE-ANNOTS" , "dhms" , "" , "Use date+time stamps in .annot output" );
  add_param( "WRITE-ANNOTS" , "collapse" , "" , "Collapse internal EDF+D offsets back to standard EDF timing" );
  add_param( "WRITE-ANNOTS" , "min-dur" , "30" , "Enforce a minimum duration for output annotations" );
  add_param( "WRITE-ANNOTS" , "meta" , "F" , "Toggle writing metadata fields" );
  add_param( "WRITE-ANNOTS" , "tab-meta" , "" , "Write metadata in tabular columns instead of key=value form" );
  add_param( "WRITE-ANNOTS" , "headers" , "" , "Include .annot header lines" );
  add_param( "WRITE-ANNOTS" , "specials" , "" , "Include Luna special annotations in output" );
  add_param( "WRITE-ANNOTS" , "minimal" , "" , "Minimal output format" );
  add_param( "WRITE-ANNOTS" , "min" , "" , "Alias for minimal" );
  add_param( "WRITE-ANNOTS" , "set-id" , "" , "Write the EDF ID as annotation metadata" );
  add_param( "WRITE-ANNOTS" , "set-id-key" , "ID" , "Metadata key name used with set-id" );
  add_param( "WRITE-ANNOTS" , "offset" , "2.5" , "Apply this time offset in seconds to all annotations on write" );
  add_param( "WRITE-ANNOTS" , "remap" , "N2|Stage2" , "On-the-fly class label remapping for output" );

  //
  // META
  //

  add_cmd( "annot" , "META" , "Annotation meta-data" );
  add_verb( "META" ,
	    "META adds or updates per-instance meta-data fields on existing annotations.\n"
	    "\n"
	    "Use it when annotations need extra attached values for later tabulation, filtering or export." );
  add_param( "META" , "annot" , "spindle" , "Annotations to update" );
  add_param( "META" , "md" , "DUR" , "Metadata field name to set" );
  add_param( "META" , "dur" , "" , "Store the annotation duration" );
  add_param( "META" , "sig" , "C3" , "Use signal-based summaries over each annotation interval" );
  add_param( "META" , "mean" , "" , "Store the mean value of the signal/window" );
  add_param( "META" , "min" , "" , "Store the minimum value of the signal/window" );
  add_param( "META" , "max" , "" , "Store the maximum value of the signal/window" );
  add_param( "META" , "range" , "" , "Store the signal range over the interval/window" );
  add_param( "META" , "other" , "arousal" , "Use other annotations to derive metadata" );
  add_param( "META" , "flatten" , "" , "Flatten other annotations before overlap-based calculations" );
  add_param( "META" , "overlap" , "" , "Store whether any other annotation overlaps the interval" );
  add_param( "META" , "complete-overlap" , "" , "Store whether the interval is completely spanned by another annotation" );
  add_param( "META" , "whole-other" , "" , "Store whether the interval completely spans another annotation" );
  add_param( "META" , "count" , "" , "Store the number of overlapping other annotations" );
  add_param( "META" , "nearest" , "60" , "Store signed time to the nearest other annotation" );
  add_param( "META" , "nearest-start" , "60" , "Store signed time to the nearest annotation start" );
  add_param( "META" , "nearest-stop" , "60" , "Store signed time to the nearest annotation stop" );
  add_param( "META" , "nearest-midpoint" , "60" , "Store signed time to the nearest annotation midpoint" );
  add_param( "META" , "w" , "0.5" , "Expand the interval symmetrically by this many seconds before evaluating" );
  add_param( "META" , "w-left" , "0.5" , "Expand the interval left by this many seconds" );
  add_param( "META" , "w-right" , "0.5" , "Expand the interval right by this many seconds" );

  //
  // ESPAN
  //

  add_cmd( "annot" , "ESPAN" , "Epoch annotation coverage" );
  add_verb( "ESPAN" ,
	    "ESPAN summarizes, epoch by epoch, how much time is covered by selected annotations.\n"
	    "\n"
	    "It is useful when you want annotation coverage expressed on the same epoch grid used by\n"
	    "other Luna analyses." );
  add_param( "ESPAN" , "annot" , "spindle" , "Annotations to summarize by epoch" );
  add_param( "ESPAN" , "sec" , "T" , "Emit spanned seconds per epoch" );
  add_param( "ESPAN" , "pct" , "T" , "Emit proportion of the epoch spanned" );
  add_param( "ESPAN" , "cnt" , "T" , "Emit the count of overlapping annotation instances" );
  add_param( "ESPAN" , "has" , "T" , "Emit a 0/1 flag for any overlap in the epoch" );
  add_param( "ESPAN" , "verbose" , "" , "Emit original and clipped per-instance overlaps" );

  add_table( "ESPAN" , "E" , "Per-epoch summary across all selected annotations" );
  add_var( "ESPAN" , "E" , "SEC" , "Seconds spanned by any selected annotation in the epoch [sec]" );
  add_var( "ESPAN" , "E" , "PCT" , "Proportion of the epoch spanned [pct]" );
  add_var( "ESPAN" , "E" , "HAS" , "Epoch has any selected annotation overlap (1/0) [has]" );
  add_var( "ESPAN" , "E" , "CNT" , "Count of overlapping annotation instances across all selected annotations [cnt]" );

  add_table( "ESPAN" , "E,ANNOT" , "Per-epoch summary within annotation class" );
  add_var( "ESPAN" , "E,ANNOT" , "SEC" , "Seconds spanned by this annotation in the epoch [sec]" );
  add_var( "ESPAN" , "E,ANNOT" , "PCT" , "Proportion of the epoch spanned by this annotation [pct]" );
  add_var( "ESPAN" , "E,ANNOT" , "HAS" , "This annotation overlaps the epoch (1/0) [has]" );
  add_var( "ESPAN" , "E,ANNOT" , "CNT" , "Count of overlapping instances of this annotation [cnt]" );

  add_table( "ESPAN" , "E,ANNOT,INST" , "Verbose per-instance overlaps within epoch [verbose]" );
  add_var( "ESPAN" , "E,ANNOT,INST" , "START" , "Original annotation start time (seconds)" );
  add_var( "ESPAN" , "E,ANNOT,INST" , "STOP" , "Original annotation stop time (seconds)" );
  add_var( "ESPAN" , "E,ANNOT,INST" , "DUR" , "Original annotation duration (seconds)" );
  add_var( "ESPAN" , "E,ANNOT,INST" , "XSTART" , "Clipped overlap start within the epoch (seconds)" );
  add_var( "ESPAN" , "E,ANNOT,INST" , "XSTOP" , "Clipped overlap stop within the epoch (seconds)" );
  add_var( "ESPAN" , "E,ANNOT,INST" , "XDUR" , "Clipped overlap duration within the epoch (seconds)" );

  //
  // A2S
  //

  add_cmd( "annot" , "A2S" , "Annotation-to-signal" );
  add_verb( "A2S" ,
	    "A2S creates a derived signal from annotations.\n"
	    "\n"
	    "Typically this means generating a 0/1 or otherwise annotation-driven time series that can be\n"
	    "used by later signal-based commands." );
  add_param( "A2S" , "annot" , "spindle" , "Annotations to convert to signals" );
  add_param( "A2S" , "sr" , "100" , "Sample rate of the derived signal" );
  add_param( "A2S" , "label" , "SPINDLE" , "Output signal label(s); defaults to the annotation names" );
  add_param( "A2S" , "numeric-inst" , "" , "Use numeric instance IDs as signal values instead of 0/1" );

  //
  // REMAP
  //

  add_cmd( "annot" , "REMAP" , "Remap annotations" );
  add_url( "REMAP" , "annotations/#remap" );
  add_verb( "REMAP" ,
	    "REMAP applies annotation label remappings after annotations have already been loaded.\n"
	    "\n"
	    "It reads one or more remapping files, updates annotation class names in memory and tracks the\n"
	    "old-to-new mapping so the changes can be inspected later via ALIASES." );
  add_param( "REMAP" , "file" , "remap.txt" , "One or more remapping definition files" );
  add_param( "REMAP" , "remap-col" , "" , "Require a leading 'remap' column in the input file" );
  add_param( "REMAP" , "optional-remap-col" , "" , "Allow, but do not require, a leading 'remap' column" );
  add_param( "REMAP" , "allow-spaces" , "T" , "Allow quoted space-delimited as well as tab-delimited fields" );
  add_param( "REMAP" , "verbose" , "" , "Log parsed remappings and remapped annotations" );

  //
  // MAKE-ANNOTS
  //

  add_cmd( "annot" , "MAKE-ANNOTS" , "Make annotations" );
  add_url( "MAKE-ANNOTS" , "annotations/#make-annots" );
  add_verb( "MAKE-ANNOTS" ,
	    "MAKE-ANNOTS is a multi-mode utility for creating new annotations from existing annotations or epochs.\n"
	    "\n"
	    "The main modes are:\n"
	    "  expr        Combine two annotations using A|B, A*B, A+B or A-B\n"
	    "  epoch       Create one flattened background annotation across unmasked epochs\n"
	    "  epoch-num   Create one distinct annotation class per unmasked epoch\n"
	    "  split       Break an annotation into epoch-bounded fragments\n"
	    "  flatten     Collapse overlapping intervals within one annotation\n"
	    "  complement  Create the inverse of one or more annotations across the record\n"
	    "  pool        Pool intervals across multiple annotations, optionally flattening overlaps\n"
	    "  w / w-left / w-right\n"
	    "              Expand annotations into windows around existing intervals\n"
	    "  midpoint / start / stop\n"
	    "              Reduce intervals to 0-duration markers\n"
	    "\n"
	    "In general, annot names the new annotation class to create. Some modes operate on\n"
	    "existing annotations via orig-annot, split, flatten, complement, pool or pool-flatten;\n"
	    "the expr mode instead uses expr to define a binary annotation set operation." );
  add_param( "MAKE-ANNOTS" , "annot" , "NEW" , "Name of the annotation class to create" );
  add_param( "MAKE-ANNOTS" , "expr" , "A|B" , "Binary set expression: A|B, A*B, A+B or A-B" );
  add_param( "MAKE-ANNOTS" , "ch" , "C3" , "Assign this channel label to created annotations in expr mode" );
  add_param( "MAKE-ANNOTS" , "epoch" , "EPOCH" , "Create annotations spanning each unmasked epoch" );
  add_param( "MAKE-ANNOTS" , "epoch-num" , "EPOCH" , "Create one distinct annotation class per unmasked epoch" );
  add_param( "MAKE-ANNOTS" , "w" , "1,5,10" , "For epoch mode, create symmetric edge windows of these durations (seconds)" );
  add_param( "MAKE-ANNOTS" , "collapse-edges" , "" , "In epoch mode, collapse left/right edge windows into one class" );
  add_param( "MAKE-ANNOTS" , "edge" , "EDGE" , "Base label for edge-window annotations in epoch mode" );
  add_param( "MAKE-ANNOTS" , "s" , "1" , "Sequence increment for edge-window instance labels" );
  add_param( "MAKE-ANNOTS" , "split" , "OLD" , "Split an annotation into epoch-bounded pieces" );
  add_param( "MAKE-ANNOTS" , "flatten" , "OLD" , "Flatten overlapping intervals from an existing annotation" );
  add_param( "MAKE-ANNOTS" , "complement" , "A,B" , "Create the complement of one or more annotations across the record" );
  add_param( "MAKE-ANNOTS" , "pool" , "A,B" , "Pool intervals from multiple annotations without flattening" );
  add_param( "MAKE-ANNOTS" , "pool-flatten" , "A,B" , "Pool intervals from multiple annotations and flatten overlaps" );
  add_param( "MAKE-ANNOTS" , "orig-annot" , "OLD" , "Source annotation class for windowing or point-reduction modes" );
  add_param( "MAKE-ANNOTS" , "w-left" , "2" , "Expand source annotations left by this many seconds" );
  add_param( "MAKE-ANNOTS" , "w-right" , "2" , "Expand source annotations right by this many seconds" );
  add_param( "MAKE-ANNOTS" , "midpoint" , "" , "Reduce each source annotation to a 0-duration midpoint marker" );
  add_param( "MAKE-ANNOTS" , "start" , "" , "Reduce each source annotation to a 0-duration start marker" );
  add_param( "MAKE-ANNOTS" , "stop" , "" , "Reduce each source annotation to a 0-duration stop marker" );


  //
  // DAYS
  //

  add_cmd( "actig" , "DAYS" , "Create day and clock-time annotations for multi-day recordings" );
  add_url( "DAYS" , "annotations/#days" );
  add_verb( "DAYS" ,
	    "DAYS creates annotations marking calendar days and optionally hourly clock-time\n"
	    "periods across the recording. Designed for multi-day recordings such as actigraphy.\n"
	    "\n"
	    "Day annotations (day01, day02, ...) span from the recording start to the first\n"
	    "day boundary, then each subsequent 24-hour period. By default days start at\n"
	    "noon; use anchor to change (e.g. anchor=0 for midnight-to-midnight).\n"
	    "\n"
	    "With the hours option, hourly annotations (00h..23h) are added, each spanning\n"
	    "the corresponding clock hour across all days. With halves, AM and PM annotations\n"
	    "are created spanning midnight-to-noon and noon-to-midnight respectively.\n"
	    "With weekend, each anchored day window that begins on Saturday or Sunday also\n"
	    "receives a weekend annotation." );

  add_param( "DAYS" , "anchor" , "12" , "Hour (0-23) at which a new day begins (default: 12 = noon)" );
  add_param( "DAYS" , "prefix" , "day" , "Prefix for day annotation labels (default: day)" );
  add_param( "DAYS" , "hours" , "" , "Also create hourly annotations (00h, 01h, ..., 23h)" );
  add_param( "DAYS" , "hour-prefix" , "" , "Prefix for hourly annotation labels" );
  add_param( "DAYS" , "halves" , "" , "Also create AM and PM annotations" );
  add_param( "DAYS" , "weekend" , "" , "Also create weekend annotations for anchored day windows starting on Saturday/Sunday" );
  add_param( "DAYS" , "weekend-label" , "WEEKEND" , "Label for weekend annotations (default: WEEKEND)" );
  add_param( "DAYS" , "verbose" , "" , "Output day-level details to the output database" );


  //
  // ACTIG
  //

  add_cmd( "actig" , "ACTIG" , "Actigraphy nonparametric metrics and wake/sleep scoring" );
  add_url( "ACTIG" , "annotations/#actig" );
  add_verb( "ACTIG" ,
	    "ACTIG computes standard nonparametric circadian rhythm metrics\n"
	    "and optionally performs wake/sleep scoring from an activity signal.\n"
	    "Gaps (non-wear, bad data) are handled automatically: run MASK and\n"
	    "RE before ACTIG to mark bad intervals. Epochs with insufficient\n"
	    "samples (< gap-min-pct of expected) are flagged as gaps and\n"
	    "excluded from all calculations. Gap runs are annotated as ACTIG_G.\n"
	    "\n"
	    "Nonparametric metrics (always computed):\n"
	    "  IS   Interdaily stability (0=random, 1=perfectly regular)\n"
	    "  IV   Intradaily variability (~0=smooth, >1=fragmented)\n"
	    "  L5   Mean activity in the least active 5h window\n"
	    "  M10  Mean activity in the most active 10h window\n"
	    "  RA   Relative amplitude: (M10-L5)/(M10+L5)\n"
	    "\n"
	    "Wake/sleep scoring methods (with score option):\n"
	    "  luna       Generic robust scorer (default): device- and epoch-size-\n"
	    "             agnostic; uses log-transform, robust normalization,\n"
	    "             local smoothing (smooth), burst density (burst), and\n"
	    "             run-length persistence rules (min-sleep, max-gap, min-wake).\n"
	    "             Works with any activity metric (ENMO, counts, MESA, etc.).\n"
	    "             Optional diagnostic channels (channels): adds _Z, _L, _D,\n"
	    "             _CS, _SW signals to the EDF for visual inspection.\n"
	    "  cole       Cole-Kripke (1992) weighted moving average; requires\n"
	    "             epoch=60 and counts-per-minute activity signal.\n"
	    "  threshold  Simple threshold on raw activity level.\n"
	    "\n"
	    "Fragmentation metrics (with score option, or existing S/W annotations):\n"
	    "  SFI      Sleep fragmentation index: sleep->wake transitions per hour of sleep\n"
	    "  SFI_ACT  Actigraphy-style FI = MI + IMM1;\n"
	    "           MI = % wake epochs in sleep period,\n"
	    "           IMM1 = % sleep bouts <= 1 epoch\n"
	    "\n"
	    "With prescored, ACTIG uses existing sleep/wake annotations instead\n"
	    "of activity-based scoring. prescored with no value defaults to S,W;\n"
	    "or specify prescored=sleep_label,wake_label.\n"
	    "With score, ACTIG creates\n"
	    "W (wake), S (sleep), and ACTIG_G (gap) annotations by default\n"
	    "(override with wake= / sleep= / gap-out=). Outputs TST/WASO/GAP\n"
	    "overall and per-day DAY strata.\n"
	    "Per-day summaries include VALID_MIN, VALID_PCT, INCLUDED, and an\n"
	    "optional day-level QC layer. INCLUDED remains coverage-only\n"
	    "(valid minutes >= day-min-valid). Cross-day averages and debt\n"
	    "windows use INCLUDED days that are not QC_DAY_EXCLUDED." );

  add_param( "ACTIG" , "sig" , "Activity" , "Activity signal to analyze (required unless prescored)" );
  add_param( "ACTIG" , "epoch" , "60" , "Epoch length in seconds for binning (default: 60)" );
  add_param( "ACTIG" , "bin" , "60" , "NP metric bin size in minutes (default: 60)" );
  add_param( "ACTIG" , "sum" , "" , "Use sum (not mean) when binning epochs" );
  add_param( "ACTIG" , "l" , "5" , "Window size in hours for least-active metric (default: 5)" );
  add_param( "ACTIG" , "m" , "10" , "Window size in hours for most-active metric (default: 10)" );
  add_param( "ACTIG" , "score" , "" , "Also perform wake/sleep scoring" );
  add_param( "ACTIG" , "prescored" , "" , "Use existing annotations as sleep/wake instead of scoring. Empty => S,W; else sleep,wake labels." );
  add_param( "ACTIG" , "method" , "luna" , "Scoring method: luna (default), cole, or threshold" );
  add_param( "ACTIG" , "thresh" , "" , "Activity threshold for threshold method (default: median of valid epochs)" );
  add_param( "ACTIG" , "cole-thresh" , "1.0" , "Wake/sleep cutoff for Cole-Kripke score (default: 1.0)" );
  add_param( "ACTIG" , "smooth" , "10" , "Luna method: smoothing window in minutes (default: 10)" );
  add_param( "ACTIG" , "burst" , "10" , "Luna method: burst-density window in minutes (default: 10)" );
  add_param( "ACTIG" , "burst-z" , "0.5" , "Luna method: z-score threshold to count an epoch as a burst (default: 0.5)" );
  add_param( "ACTIG" , "quiet-z" , "-0.5" , "Luna method: smoothed z-score below which candidate sleep is triggered (default: -0.5)" );
  add_param( "ACTIG" , "active-frac" , "0.20" , "Luna method: max burst-density fraction allowed for candidate sleep (default: 0.20)" );
  add_param( "ACTIG" , "min-sleep" , "15" , "Luna method: minimum sleep run duration in minutes (default: 15)" );
  add_param( "ACTIG" , "max-gap" , "2" , "Luna method: max within-sleep non-sleep gap to fill in minutes (default: 2)" );
  add_param( "ACTIG" , "min-wake" , "5" , "Luna method: consecutive wake minutes required to end a sleep bout (default: 5)" );
  add_param( "ACTIG" , "channels" , "" , "Luna method: add diagnostic EDF channels (_Z, _L, _D, _CS, _SW) to the EDF" );
  add_param( "ACTIG" , "wake" , "W" , "Label for wake annotations when score is used (default: W)" );
  add_param( "ACTIG" , "sleep" , "S" , "Label for sleep annotations when score is used (default: S)" );
  add_param( "ACTIG" , "gap-out" , "ACTIG_G" , "Label for gap annotations (default: ACTIG_G)" );
  add_param( "ACTIG" , "day-anchor" , "12" , "Hour (0-23) at which per-day boundaries occur (default: 12 = noon)" );
  add_param( "ACTIG" , "gap-min-pct" , "50" , "Min % of expected samples per epoch to be considered valid (default: 50)" );
  add_param( "ACTIG" , "day-min-valid" , "960" , "Min valid minutes per day for inclusion in cross-day averages (default: 960 = 16h)" );
  add_param( "ACTIG" , "qc-day" , "T" , "Enable day-level QC pass after scoring (default: T)" );
  add_param( "ACTIG" , "qc-exclude-flat" , "T" , "Enable flatline-style technical exclusion flag (default: T)" );
  add_param( "ACTIG" , "qc-exclude-lowvar" , "T" , "Enable collapsed-variability technical exclusion flag (default: T)" );
  add_param( "ACTIG" , "qc-exclude-nearfloor" , "T" , "Enable near-floor technical exclusion flag (default: T)" );
  add_param( "ACTIG" , "qc-warn-implausible" , "T" , "Enable combined plausibility warning (default: T)" );
  add_param( "ACTIG" , "qc-warn-longsleep" , "T" , "Enable long-sleep warning (default: T)" );
  add_param( "ACTIG" , "qc-warn-highsleep" , "T" , "Enable high sleep-fraction warning (default: T)" );
  add_param( "ACTIG" , "qc-exclude-out" , "ACTIG_QC_EXCLUDED" , "Annotation label for QC-excluded day windows; empty disables day-QC annotations" );
  add_param( "ACTIG" , "qc-flat-frac-th" , "0.80" , "Flatline flag threshold: fraction of adjacent valid epoch pairs with abs(delta) <= qc-flat-delta-th" );
  add_param( "ACTIG" , "qc-flat-delta-th" , "0.0" , "Flatline flag delta threshold on adjacent valid epochs" );
  add_param( "ACTIG" , "qc-lowvar-frac-th" , "0.80" , "Low-variability flag threshold on fraction of valid epochs in a rolling low-variability state" );
  add_param( "ACTIG" , "qc-lowvar-cv-th" , "0.05" , "Low-variability flag threshold on day-level coefficient of variation" );
  add_param( "ACTIG" , "qc-nearfloor-frac-th" , "0.85" , "Near-floor flag threshold on fraction of valid epochs near the day-specific lower tail" );
  add_param( "ACTIG" , "qc-nearfloor-q-th" , "0.05" , "Quantile used to define the day-specific near-floor reference (default: P05)" );
  add_param( "ACTIG" , "qc-min-active-epochs" , "24" , "Minimum active epochs required to avoid low-structure technical flags" );
  add_param( "ACTIG" , "qc-warn-sleep-pct" , "0.85" , "Warning threshold on SCORE_SLEEP_PCT expressed as a 0..1 fraction" );
  add_param( "ACTIG" , "qc-warn-longsleep-h" , "16" , "Warning threshold on longest sleep bout in hours" );
  add_param( "ACTIG" , "qc-warn-max-sleep-h" , "20" , "Higher warning threshold on extremely long sleep bout in hours" );
  add_param( "ACTIG" , "qc-warn-low-wakeruns" , "1" , "Warning threshold on daily wake-run count" );
  add_param( "ACTIG" , "np-step" , "1" ,
	     "Sliding-window step size in minutes for L5/M10 onset search (default: 1). "
	     "Must divide 60 evenly (e.g. 1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 30, 60). "
	     "Finer steps give sub-hour onset precision from the 1-hour-binned profile." );
  add_param( "ACTIG" , "np-traditional" , "" ,
	     "Revert to traditional 60-minute onset search step (integer-hour onsets, HH:00:00). "
	     "Equivalent to np-step=60." );

  add_table( "ACTIG" , "" , "Individual-level output" );
  add_var( "ACTIG" , "" , "NP_IS" , "Interdaily stability" );
  add_var( "ACTIG" , "" , "NP_IV" , "Intradaily variability" );
  add_var( "ACTIG" , "" , "NP_RA" , "Relative amplitude: (M-L)/(M+L)" );
  add_var( "ACTIG" , "" , "NP_L5" , "Mean activity in the least-active 5h window" );
  add_var( "ACTIG" , "" , "NP_L5_ONSET" , "Clock time of L5 window onset (HH:MM:SS; sub-hour precision with default np-step=1)" );
  add_var( "ACTIG" , "" , "NP_L5_ONSET_MIN" , "L5 onset in decimal minutes from midnight (0=00:00, 120=02:00, 1380=23:00); use for regression/correlation" );
  add_var( "ACTIG" , "" , "NP_M10" , "Mean activity in the most-active 10h window" );
  add_var( "ACTIG" , "" , "NP_M10_ONSET" , "Clock time of M10 window onset (HH:MM:SS)" );
  add_var( "ACTIG" , "" , "NP_M10_ONSET_MIN" , "M10 onset in decimal minutes from midnight" );
  add_var( "ACTIG" , "" , "NP_NE" , "Total epoch bin count (including gaps)" );
  add_var( "ACTIG" , "" , "NP_NE_VALID" , "Valid (non-gap) epoch bin count" );
  add_var( "ACTIG" , "" , "NP_NE_GAP" , "Gap epoch bin count" );
  add_var( "ACTIG" , "" , "NP_NBINS" , "Total NP bin count (including gaps)" );
  add_var( "ACTIG" , "" , "NP_NBINS_VALID" , "Valid NP bin count" );
  add_var( "ACTIG" , "" , "NP_NDAYS" , "Number of complete recording days" );
  add_var( "ACTIG" , "" , "SCORE_TST_MIN" , "Total sleep time (minutes, valid epochs only)" );
  add_var( "ACTIG" , "" , "SCORE_WASO_MIN" , "Total wake time (minutes, valid epochs only)" );
  add_var( "ACTIG" , "" , "SCORE_GAP_MIN" , "Total gap time (minutes)" );
  add_var( "ACTIG" , "" , "SCORE_SLEEP_PCT" , "Sleep % of valid (non-gap) scored epochs" );
  add_var( "ACTIG" , "" , "FRAG_SFI" , "Luna-style SFI: sleep->wake transitions per hour sleep" );
  add_var( "ACTIG" , "" , "FRAG_SFI_N" , "Number of sleep->wake transitions" );
  add_var( "ACTIG" , "" , "FRAG_SFI_ACT" , "Actigraphy FI: FRAG_MI_PCT + FRAG_IMM1_PCT" );
  add_var( "ACTIG" , "" , "FRAG_MI_PCT" , "Movement index (% wake in sleep period, valid epochs only)" );
  add_var( "ACTIG" , "" , "FRAG_IMM1_PCT" , "Immobility fragmentation (% sleep bouts <=1 minute)" );
  add_var( "ACTIG" , "" , "FRAG_IMM_BOUT_N" , "Total sleep bout count in sleep period" );
  add_var( "ACTIG" , "" , "FRAG_IMM1_BOUT_N" , "Sleep bout count <=1 minute in sleep period" );
  add_var( "ACTIG" , "" , "SCORE_WAKE_RUN_N" , "Number of contiguous wake annotation runs" );
  add_var( "ACTIG" , "" , "SCORE_SLEEP_RUN_N" , "Number of contiguous sleep annotation runs" );
  add_var( "ACTIG" , "" , "SCORE_GAP_RUN_N" , "Number of contiguous gap annotation runs" );
  add_var( "ACTIG" , "" , "DAY_N" , "Total number of days" );
  add_var( "ACTIG" , "" , "DAY_N_INCLUDED" , "Days meeting day-min-valid threshold (coverage-only)" );
  add_var( "ACTIG" , "" , "DAY_N_EXCLUDED" , "Days below day-min-valid threshold (coverage-only)" );
  add_var( "ACTIG" , "" , "DAY_N_QC_OK" , "Days included after QC with no warnings" );
  add_var( "ACTIG" , "" , "DAY_N_QC_EXCLUDED" , "Days excluded after applying coverage + day QC" );
  add_var( "ACTIG" , "" , "DAY_N_QC_WARN" , "Days included after QC but flagged with warnings" );
  add_var( "ACTIG" , "" , "DAY_N_INCLUDED_POSTQC" , "Days contributing to cross-day averages after coverage + QC filtering" );
  add_var( "ACTIG" , "" , "SCORE_TST_DAYAVG_MIN" , "Mean TST across post-QC included days (minutes)" );
  add_var( "ACTIG" , "" , "SCORE_WASO_DAYAVG_MIN" , "Mean WASO across post-QC included days (minutes)" );
  add_var( "ACTIG" , "" , "SCORE_SLEEP_PCT_DAYAVG" , "Mean sleep % across post-QC included days" );

  add_table( "ACTIG" , "DAY" , "Per-day output (day-anchor boundaries)" );
  add_var( "ACTIG" , "DAY" , "SCORE_TST_MIN" , "Daily total sleep time (minutes)" );
  add_var( "ACTIG" , "DAY" , "SCORE_WASO_MIN" , "Daily wake time (minutes)" );
  add_var( "ACTIG" , "DAY" , "SCORE_GAP_MIN" , "Daily gap time (minutes)" );
  add_var( "ACTIG" , "DAY" , "VALID_MIN" , "Daily valid (non-gap) minutes" );
  add_var( "ACTIG" , "DAY" , "VALID_PCT" , "Daily valid % of total day duration" );
  add_var( "ACTIG" , "DAY" , "INCLUDED" , "1 if day meets day-min-valid threshold, 0 otherwise" );
  add_var( "ACTIG" , "DAY" , "QC_DAY_OK" , "1 if day is included after QC and has no warnings" );
  add_var( "ACTIG" , "DAY" , "QC_DAY_EXCLUDED" , "1 if day is excluded after applying coverage + day QC" );
  add_var( "ACTIG" , "DAY" , "QC_DAY_WARN" , "1 if day remains included but is flagged with warnings" );
  add_var( "ACTIG" , "DAY" , "QC_TECH_FLAG_N" , "Number of technical day-QC flags triggered" );
  add_var( "ACTIG" , "DAY" , "SCORE_SLEEP_PCT" , "Daily sleep % (of valid scored epochs)" );
  add_var( "ACTIG" , "DAY" , "FRAG_SFI" , "Daily sleep->wake transitions per hour sleep" );
  add_var( "ACTIG" , "DAY" , "FRAG_SFI_N" , "Daily sleep->wake transition count" );
  add_var( "ACTIG" , "DAY" , "FRAG_SFI_ACT" , "Daily actigraphy FI: MI_PCT + IMM1_PCT" );
  add_var( "ACTIG" , "DAY" , "FRAG_MI_PCT" , "Daily movement index (% wake in sleep period)" );
  add_var( "ACTIG" , "DAY" , "FRAG_IMM1_PCT" , "Daily immobility fragmentation (% bouts <=1 min)" );
  add_var( "ACTIG" , "DAY" , "FRAG_IMM_BOUT_N" , "Daily sleep bout count in sleep period" );
  add_var( "ACTIG" , "DAY" , "FRAG_IMM1_BOUT_N" , "Daily sleep bout count <=1 minute" );
  add_var( "ACTIG" , "DAY" , "SCORE_EPOCH_N" , "Daily valid (non-gap) epoch count" );
  add_var( "ACTIG" , "DAY" , "SCORE_WAKE_EPOCH_N" , "Daily wake epoch count" );
  add_var( "ACTIG" , "DAY" , "QC_FLAT_FRAC" , "Daily fraction of adjacent valid epoch pairs with minimal change" );
  add_var( "ACTIG" , "DAY" , "QC_LOWVAR_CV" , "Daily coefficient of variation of valid raw activity" );
  add_var( "ACTIG" , "DAY" , "QC_NEARFLOOR_FRAC" , "Daily fraction of valid epochs near the day-specific lower tail" );
  add_var( "ACTIG" , "DAY" , "QC_ACTIVE_EPOCH_N" , "Daily active-structure epoch count based on within-day activity reference" );
  add_var( "ACTIG" , "DAY" , "QC_LONGEST_QUIET_RUN_MIN" , "Daily longest contiguous run of near-floor epochs (minutes)" );
  add_var( "ACTIG" , "DAY" , "QC_LONGEST_LOWVAR_RUN_MIN" , "Daily longest contiguous rolling low-variability run (minutes)" );
  add_var( "ACTIG" , "DAY" , "QC_ACT_MED" , "Daily median raw activity across valid epochs" );
  add_var( "ACTIG" , "DAY" , "QC_ACT_MAD" , "Daily MAD of raw activity across valid epochs" );
  add_var( "ACTIG" , "DAY" , "QC_ACT_SD" , "Daily SD of raw activity across valid epochs" );
  add_var( "ACTIG" , "DAY" , "QC_ACT_P05" , "Daily P05 raw activity across valid epochs" );
  add_var( "ACTIG" , "DAY" , "QC_ACT_P95" , "Daily P95 raw activity across valid epochs" );
  add_var( "ACTIG" , "DAY" , "QC_LONGEST_SLEEP_BOUT_MIN" , "Daily longest contiguous sleep bout (minutes)" );
  add_var( "ACTIG" , "DAY" , "QC_LONGEST_WAKE_BOUT_MIN" , "Daily longest contiguous wake bout (minutes)" );
  add_var( "ACTIG" , "DAY" , "QC_FLAG_FLAT" , "Daily technical flag: flatline-like activity" );
  add_var( "ACTIG" , "DAY" , "QC_FLAG_LOWVAR" , "Daily technical flag: collapsed variability and low active structure" );
  add_var( "ACTIG" , "DAY" , "QC_FLAG_NEARFLOOR" , "Daily technical flag: near-floor day with low active structure" );
  add_var( "ACTIG" , "DAY" , "QC_WARN_HIGHSLEEP" , "Daily warning: high sleep fraction" );
  add_var( "ACTIG" , "DAY" , "QC_WARN_LONGSLEEP" , "Daily warning: very long sleep bout" );
  add_var( "ACTIG" , "DAY" , "QC_WARN_LOWWAKERUNS" , "Daily warning: too few wake runs" );

  add_param( "ACTIG" , "debt" , "" , "Compute local sleep debt relative to a target night" );
  add_param( "ACTIG" , "debt-target" , "" ,
	     "Target night: 1-based day number or date (YYYY-MM-DD or DD.MM.YY); required with debt" );
  add_param( "ACTIG" , "debt-recent" , "2" ,
	     "Number of nights immediately before target forming the 'recent' window (default: 2)" );
  add_param( "ACTIG" , "debt-base" , "7" ,
	     "Number of nights before the recent window forming the baseline (default: 7)" );
  add_param( "ACTIG" , "debt-min-base" , "3" ,
	     "Minimum valid baseline nights for DEBT_TST_DELTA / DEBT_TST_REL (default: 3)" );
  add_param( "ACTIG" , "debt-min-z" , "5" ,
	     "Minimum valid baseline nights for z-scores and DEBT_INDEX (default: 5)" );
  add_param( "ACTIG" , "debt-w" , "0.5" ,
	     "Weight given to TST in DEBT_INDEX (0=fragmentation only, 1=TST only; default: 0.5)" );

  add_var( "ACTIG" , "" , "DEBT_TARGET_DAY" , "Target day number (1-based) used for debt calculation" );
  add_var( "ACTIG" , "" , "DEBT_N_RECENT"   , "Valid included nights in the recent window" );
  add_var( "ACTIG" , "" , "DEBT_N_BASE"     , "Valid included nights in the baseline window" );
  add_var( "ACTIG" , "" , "DEBT_TST_RECENT" , "Mean TST (min) over recent nights" );
  add_var( "ACTIG" , "" , "DEBT_TST_BASE"   , "Median TST (min) over baseline nights" );
  add_var( "ACTIG" , "" , "DEBT_FRAG_RECENT", "Mean SFI over recent nights (sleep->wake transitions per hour)" );
  add_var( "ACTIG" , "" , "DEBT_FRAG_BASE"  , "Median SFI over baseline nights" );
  add_var( "ACTIG" , "" , "DEBT_TST_DELTA"  , "Baseline median TST minus recent mean TST (min; positive = deficit); needs debt-min-base" );
  add_var( "ACTIG" , "" , "DEBT_TST_REL"    , "Fractional TST deficit: DELTA / baseline median (0..1); needs debt-min-base" );
  add_var( "ACTIG" , "" , "DEBT_TST_Z"      , "Z-score of recent mean TST vs baseline distribution (negative = deprived); needs debt-min-z" );
  add_var( "ACTIG" , "" , "DEBT_FRAG_Z"     , "Z-score of recent mean SFI vs baseline distribution (positive = more fragmented); needs debt-min-z" );
  add_var( "ACTIG" , "" , "DEBT_INDEX"      , "Composite debt: w*(-TST_Z) + (1-w)*FRAG_Z; higher = more deprived; needs debt-min-z" );


  //
  // S2A
  //

  //
  // S2C
  //

  add_param( "S2C" , "sig" , "C3,C4" , "Signal(s) to analyze" );
  add_param( "S2C" , "seed" , "C3" , "Segmenting seed signal(s); defaults to sig" );
  add_param( "S2C" , "all-by-all" , "" , "Compare all sig channels against all seg channels" );

  add_param( "S2C" , "waves" , "waves" , "Base label for created cycle annotations" );
  add_param( "S2C" , "add-channel-inst-label" , "" , "Add channel names to instance IDs" );
  add_param( "S2C" , "add-channel-class-label" , "" , "Add channel names to class labels" );
  add_param( "S2C" , "half-waves" , "" , "Also annotate POS/NEG half-waves; optional value sets base label, e.g. half-waves=HW" );
  add_param( "S2C" , "peak-points" , "PEAK" , "Add 0-duration POS/NEG peak-point annotations with this base label" );
  add_param( "S2C" , "waves-bins" , "bins" , "Annotate 12 phase bins (B01..B12) using the cycle warp" );

  add_param( "S2C" , "pos2neg" , "" , "Segment on POS->NEG zero-crossings (default POS->NEG)" );

  add_param( "S2C" , "t-min" , "2.0" , "Minimum cycle duration (seconds)" );
  add_param( "S2C" , "t-max" , "20.0" , "Maximum cycle duration (seconds)" );
  add_param( "S2C" , "t-min-neg" , "1.0" , "Minimum negative half-wave duration (seconds)" );
  add_param( "S2C" , "t-max-neg" , "10.0" , "Maximum negative half-wave duration (seconds)" );
  add_param( "S2C" , "t-min-pos" , "1.0" , "Minimum positive half-wave duration (seconds)" );
  add_param( "S2C" , "t-max-pos" , "10.0" , "Maximum positive half-wave duration (seconds)" );
  add_param( "S2C" , "no-halfwave-t" , "" , "Disable automatic half-wave duration defaults" );

  add_param( "S2C" , "mag-percentile" , "0.8" , "Minimum percentile of cycle magnitude (0-1)" );
  add_param( "S2C" , "mag-z" , "0.5" , "Minimum z-score of cycle magnitude" );

  add_param( "S2C" , "bootstrap" , "" , "Enable bootstrap SE/CI for mean grids" );
  add_param( "S2C" , "bootstrap-n" , "1000" , "Number of bootstrap resamples" );
  add_param( "S2C" , "bootstrap-ci" , "0.95" , "Bootstrap CI level" );

  add_param( "S2C" , "lag-window" , "0" , "Lag window half-width in seconds (0=full cycle)" );
  add_param( "S2C" , "lag-abs" , "" , "Use absolute value for lag peak search" );

  add_param( "S2C" , "emit-per-cycle" , "" , "Emit per-cycle metrics" );

  add_param( "S2C" , "amp-bins" , "10" , "Number of amplitude bins for density matrices" );

  add_param( "S2C" , "time-domain" , "" , "Emit time-domain grid and window metrics" );
  add_param( "S2C" , "time-window" , "100" , "Time-window width in seconds" );
  add_param( "S2C" , "time-bin" , "1" , "Time bin size in seconds" );
  add_param( "S2C" , "time-min-n" , "1" , "Minimum N per time bin" );
  add_param( "S2C" , "time-lock" , "pos" , "Time-lock anchor: pos or neg" );
  add_param( "S2C" , "emit-td-grid" , "" , "Emit time-domain SEC and SECxAMP grids" );
  add_param( "S2C" , "emit-td-summary" , "" , "Emit time-domain TD_* summary metrics" );

  add_param( "S2C" , "emit-ph-grid" , "" , "Emit phase-bin (PH) mean waveform" );
  add_param( "S2C" , "emit-ph-amp" , "" , "Emit phase-by-amplitude density matrix" );

  add_param( "S2C" , "emit-seed" , "" , "Emit SEED-level summary statistics" );
  add_param( "S2C" , "emit-sig" , "" , "Emit SEEDxSIG summary statistics" );
  add_param( "S2C" , "emit-asym" , "" , "Emit asymmetry metrics (DUR_RATIO/AMP_ASYM/SLOPE_ASYM/SHARP_ASYM)" );

  add_param( "S2C" , "emit-se" , "" , "Emit SE variables in bootstrap outputs" );
  add_param( "S2C" , "emit-mad" , "" , "Emit MAD variables" );
  add_param( "S2C" , "emit-asymm" , "" , "Alias for emit-asym" );
  add_param( "S2C" , "emit-asymmetry" , "" , "Alias for emit-asym" );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // EXPRESSIONS
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // --eval
  //

  add_cmd( "helpers" , "--eval" , "Evaluate expressions from stdin (or --eval-verbose)" );
  add_param( "--eval" , "{expr}" , "2+2" , "Evaluate expression (via stdin)" );

  //
  // EVAL
  //

  add_cmd( "expr" , "EVAL" , "Annotation expressions" );
  add_url( "EVAL" , "evals/#eval" );
  add_verb( "EVAL" ,
	    "EVAL creates a new annotation track by evaluating a Boolean expression over existing\n"
	    "annotations. It can work epoch-by-epoch, or in interval mode where new intervals are\n"
	    "constructed from contiguous true regions across the annotation timeline.\n"
	    "\n"
	    "Use annot to name the new annotation class in epoch mode, or interval to switch into\n"
	    "interval mode and name the derived interval annotation. globals can be used to maintain\n"
	    "running accumulator variables while the expression is evaluated." );
  add_param( "EVAL" , "annot" , "N2_OR_REM" , "Name of the new annotation class in epoch mode" );
  add_param( "EVAL" , "interval" , "A3" , "Use interval mode and name the derived interval annotation class" );
  add_param( "EVAL" , "expr" , "\"A && B\"" , "Annotation expression to evaluate" );
  add_param( "EVAL" , "globals" , "N,TOTAL" , "Accumulator variables to maintain while evaluating the expression" );

  //
  // TRANS
  //

  add_cmd( "expr" , "TRANS" , "Signal expressions" );
  add_url( "TRANS" , "evals/#trans" );
  add_verb( "TRANS" ,
	    "TRANS evaluates an expression over one or more signals. In signal mode it creates or\n"
	    "updates a channel; in annotation mode it creates a new annotation track from the true\n"
	    "regions of the expression.\n"
	    "\n"
	    "Set sig to an existing or new channel label to write a signal, or set sig=* together\n"
	    "with annot to emit annotations instead. All referenced signals must share a sampling rate." );
  add_param( "TRANS" , "sig" , "C3_DENOISED" , "Target signal label, or * to switch to annotation mode" );
  add_param( "TRANS" , "annot" , "AROUSAL_LIKE" , "Name of the new annotation class when sig=*" );
  add_param( "TRANS" , "expr" , "\"C3_M2 > 75\"" , "Signal expression to evaluate" );
  add_param( "TRANS" , "verbose" , "" , "Print additional progress information" );

  //
  // DERIVE
  //

  add_cmd( "expr" , "DERIVE" , "Per-record derived variables" );
  add_url( "DERIVE" , "evals/#derive" );
  add_verb( "DERIVE" ,
	    "DERIVE evaluates an expression once over the full record and writes the requested\n"
	    "accumulator variables as per-record outputs. It is intended for deriving observation-level\n"
	    "metrics from existing annotations and their meta-data, rather than creating new annotations.\n"
	    "\n"
	    "Use var to name the accumulator variables you want to expose in the output. req can be used\n"
	    "to require specific annotation meta-data fields before instances contribute, and if an epoch\n"
	    "mask is present then any/all/start control how annotations are filtered against that mask." );
  add_param( "DERIVE" , "expr" , "\"SUM = SUM + DUR\"" , "Expression to evaluate over the whole record" );
  add_param( "DERIVE" , "var" , "SUM,N" , "Accumulator variables to output" );
  add_param( "DERIVE" , "req" , "sp.meta" , "Require these annotation meta-data fields before including an instance" );
  add_param( "DERIVE" , "any" , "" , "Keep annotations with any overlap with an unmasked region (default)" );
  add_param( "DERIVE" , "all" , "" , "Only keep annotations fully contained in unmasked regions" );
  add_param( "DERIVE" , "start" , "" , "Keep annotations that start in an unmasked region" );
  add_table( "DERIVE" , "" , "Per-record derived outputs" );
  add_var( "DERIVE" , "" , "REQN" , "Number of annotation instances that satisfied req filters" );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // EPOCHS
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // EPOCH
  //

  add_cmd( "epoch" , "EPOCH" , "Define epochs" );
  add_url ( "EPOCH" , "epochs/#epoch" );
  add_verb( "EPOCH" ,
	    "EPOCH defines the epoch grid for subsequent epoch-based analyses. It supports standard\n"
	    "fixed-length epochs, generic epochs built from annotations, contiguous-segment epochs,\n"
	    "and table/dump modes that report the current epoch structure without changing it.\n"
	    "\n"
	    "In standard mode, set the epoch duration and increment, optionally with an offset or by\n"
	    "aligning to the first matching annotation. In generic mode, annot defines the source\n"
	    "annotations and options such as fixed, midpoint/start/stop, w, shift, trunc and flatten\n"
	    "control how annotation-derived epochs are constructed." );

  add_param( "EPOCH" , "len" , "30" , "Epoch length (seconds), defaults to 30" );
  add_param( "EPOCH" , "dur" , "30" , "Same as len" );
  add_param( "EPOCH" , "inc" , "30" , "Epoch increment (seconds), defaults to len (i.e. no overlap)" );
  add_param( "EPOCH" , "epoch" , "30,15" , "Same as len=30 inc=15" );
  add_param( "EPOCH" , "offset" , "22:00:00" , "Epoch start offset in seconds or clock time" );
  add_param( "EPOCH" , "align" , "N1,N2,N3,R,W,?,L,U,M" , "Align epochs to the first matching annotation" );
  add_param( "EPOCH" , "splice" , "" , "Define standard epochs while splicing out EDF+D gaps" );
  add_param( "EPOCH" , "splice-gaps" , "" , "Same as splice" );
  add_param( "EPOCH" , "require" , "10" , "Stop processing that EDF if there are not at least N epochs" );
  add_param( "EPOCH" , "verbose" , "" , "Output epoch-level information" );
  add_param( "EPOCH" , "dump" , "" , "Report the current epoch table without changing epochs" );
  add_param( "EPOCH" , "table" , "" , "Same as dump" );
  add_param( "EPOCH" , "masked" , "" , "In dump/table mode, include masked as well as unmasked epochs" );
  add_param( "EPOCH" , "clear" , "" , "Unepoch all signals" );
  add_param( "EPOCH" , "min" , "" , "Print only the number of epochs to stdout" );
  add_param( "EPOCH" , "contig" , "" , "Create one epoch per contiguous valid data segment" );
  add_param( "EPOCH" , "annot" , "spindle" , "Create generic epochs from these source annotations" );
  add_param( "EPOCH" , "fixed" , "5" , "For generic epochs, impose a fixed epoch duration in seconds" );
  add_param( "EPOCH" , "only-one" , "" , "With fixed generic epochs, add at most one epoch per annotation" );
  add_param( "EPOCH" , "else" , "OTHER" , "Create a complementary annotation/epoch set for regions not matched by annot" );
  add_param( "EPOCH" , "midpoint" , "" , "For generic epochs, center the window on the annotation midpoint" );
  add_param( "EPOCH" , "start" , "" , "For generic epochs, center the window on the annotation start" );
  add_param( "EPOCH" , "stop" , "" , "For generic epochs, center the window on the annotation stop" );
  add_param( "EPOCH" , "w" , "0.5" , "For generic epochs, symmetric window size around the chosen anchor point" );
  add_param( "EPOCH" , "w-before" , "0.25" , "For generic epochs, left window size around the chosen anchor point" );
  add_param( "EPOCH" , "w-after" , "0.75" , "For generic epochs, right window size around the chosen anchor point" );
  add_param( "EPOCH" , "shift" , "0.25" , "For generic epochs, shift the resulting window by this many seconds" );
  add_param( "EPOCH" , "trunc" , "0.1" , "For generic epochs, truncate this many seconds from the end" );
  add_param( "EPOCH" , "flatten" , "F" , "For generic epochs, keep overlapping/contiguous epochs separate" );
  add_param( "EPOCH" , "debug" , "" , "Verbose debug logging for generic epoch construction" );

  add_table( "EPOCH" , "" , "Epoch-level summaries" );
  add_var( "EPOCH" , "" , "DUR" , "Epoch duration (seconds)" );
  add_var( "EPOCH" , "" , "INC" , "Epoch increment (seconds)" );
  add_var( "EPOCH" , "" , "NE" , "Number of epochs" );
  add_var( "EPOCH" , "" , "NE_MASKED" , "Number of masked epochs when masked epochs are included" );
  add_var( "EPOCH" , "" , "FIXED_DUR" , "Fixed epoch duration in seconds, or 0 for variable-length generic epochs" );
  add_var( "EPOCH" , "" , "GENERIC" , "0/1 indicator that generic rather than standard epochs are in use" );
  add_var( "EPOCH" , "" , "CONTIG" , "0/1 indicator that epochs were defined from contiguous valid-data segments" );
  add_var( "EPOCH" , "" , "OFFSET" , "Epoch offset from EDF start (seconds)" );
  add_var( "EPOCH" , "" , "TOT_DUR" , "Total epoch duration" );
  add_var( "EPOCH" , "" , "TOT_PCT" , "Percent of record epoched" );
  add_var( "EPOCH" , "" , "TOT_REC" , "Total record duration" );
  add_var( "EPOCH" , "" , "TOT_SPANNED" , "Total duration spanned by epoch" );
  add_var( "EPOCH" , "" , "TOT_UNSPANNED" , "Total record duration not spanned by any epoch" );

  add_table( "EPOCH" , "E" , "Per-epoch interval information [verbose]" );
  add_var( "EPOCH" , "E" , "E1" , "Current epoch number (which may differ from E if the EDF has been restructured)" );
  add_var( "EPOCH" , "E" , "EMASK" , "0/1 indicator that this epoch is currently masked [masked dump mode]" );
  add_var( "EPOCH" , "E" , "HMS" , "Clock-time for epoch start (hh:mm:ss)" );
  add_var( "EPOCH" , "E" , "LABEL" , "Epoch label" );
  add_var( "EPOCH" , "E" , "INTERVAL" , "String label of epoch interval (seconds)" );
  add_var( "EPOCH" , "E" , "DUR" , "Epoch duration (seconds)" );
  add_var( "EPOCH" , "E" , "MID" , "Midpoint of epoch (seconds elapsed from EDF start)" );
  add_var( "EPOCH" , "E" , "START" , "Start of epoch (seconds elapsed from EDF start)" );
  add_var( "EPOCH" , "E" , "STOP" , "Stop of epoch (seconds elapsed from EDF start)" );
  add_var( "EPOCH" , "E" , "TP" , "Interval in time-points" );

  //
  // EPOCH-ANNOT
  //

  add_cmd( "epoch" , "EPOCH-ANNOT" , "Load epoch annotations" );
  add_url( "EPOCH-ANNOT" , "epochs/#epoch-annot" );  
  add_verb( "EPOCH-ANNOT" ,
	    "EPOCH-ANNOT reads one label per epoch from a simple text file and maps those values onto\n"
	    "the current epoch grid as epoch-level annotations. The file must contain exactly one row per\n"
	    "total epoch in the current EDF, and optional recodes can rename labels as they are loaded.\n"
	    "\n"
	    "This is useful for attaching staging or other epoch-wise labels after epochs have already\n"
	    "been defined." );
  add_param( "EPOCH-ANNOT" , "file" , "annots/id1.epochs" , "File path/name to read annotations from [required]" );
  add_param( "EPOCH-ANNOT" , "recode" , "NREM1=N1,NREM2=N2" , "Comma-delimited list of recodings (from=to)");

  /////////////////////////////////////////////////////////////////////////////////
  //
  // MASKS
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // MASK
  //

  add_cmd( "mask" , "MASK" , "Standard epoch mask" );
  add_url( "MASK" , "masks/#mask" );
  add_verb( "MASK" ,
	    "MASK modifies the standard epoch mask used by Luna. It supports annotation-based and\n"
	    "expression-based masking, direct epoch/time selections, random and stability-based\n"
	    "selection, and utilities that trim or reshape the current mask.\n"
	    "\n"
	    "Most forms operate in one of three modes: mask matching epochs, unmask matching epochs,\n"
	    "or force the whole mask from the match result. Annotation selectors also support OR/AND\n"
	    "logic, instance-value filters, and complete-overlap matching with a leading '+'." );

  add_param( "MASK" , "if"    , "NREM2" , "Mask NREM2 epochs, unmask all others" );
  add_param( "MASK" , "if-any" , "N2,REM" , "Mask matching epochs and unmask all others, using OR logic" );
  add_param( "MASK" , "if-all" , "N2,REM" , "Mask matching epochs and unmask all others, using AND logic" );
  add_param( "MASK" , "ifnot" , "NREM2" , "Unmask NREM2 epochs, mask all others" );
  add_param( "MASK" , "ifnot-any" , "N2,REM" , "Unmask epochs lacking any of these annotations, mask all others" );
  add_param( "MASK" , "ifnot-all" , "N2,REM" , "Unmask epochs lacking all of these annotations, mask all others" );
  add_param( "MASK" , "expr" , "A>2" , "Mask epochs with A>2, unmask all others" );
  add_param( "MASK" , "not-expr" , "A>2" , "Unmask epochs with A>2, mask all others" );

  add_param( "MASK" , "mask-if" , "NREM2" , "Mask NREM2 epochs" );
  add_param( "MASK" , "mask-if-any" , "N2,REM" , "Mask epochs matching any of these annotations" );
  add_param( "MASK" , "mask-if-all" , "N2,REM" , "Mask epochs matching all of these annotations" );
  add_param( "MASK" , "mask-ifnot" , "NREM2" , "Mask non-NREM2 epochs" );
  add_param( "MASK" , "mask-ifnot-any" , "N2,REM" , "Mask epochs lacking any of these annotations" );
  add_param( "MASK" , "mask-ifnot-all" , "N2,REM" , "Mask epochs lacking all of these annotations" );
  add_param( "MASK" , "mask-expr" , "A>2" , "Mask epochs with A>2" );

  add_param( "MASK" , "unmask-if" , "NREM2" , "Unmask NREM2 epochs" );
  add_param( "MASK" , "unmask-if-any" , "N2,REM" , "Unmask epochs matching any of these annotations" );
  add_param( "MASK" , "unmask-if-all" , "N2,REM" , "Unmask epochs matching all of these annotations" );
  add_param( "MASK" , "unmask-ifnot" , "NREM2" , "Unmask non-NREM2 epochs" );
  add_param( "MASK" , "unmask-ifnot-any" , "N2,REM" , "Unmask epochs lacking any of these annotations" );
  add_param( "MASK" , "unmask-ifnot-all" , "N2,REM" , "Unmask epochs lacking all of these annotations" );
  add_param( "MASK" , "unmask-expr" , "A>2" , "Unmask epochs with A>2" );

  add_param( "MASK" , "none" , "" ,  "Clear mask (i.e. unmask all)" );
  add_param( "MASK" , "clear" , "" , "Clear mask (i.e. unmask all)" );
  add_param( "MASK" , "include-all" , "" , "Clear mask (i.e. unmask all)" );

  add_param( "MASK" , "all" , "" , "Mask all epochs" );
  add_param( "MASK" , "total" , "" , "Mask all epochs" );
  add_param( "MASK" , "exclude-all" , "" , "Mask all epochs" );

  add_param( "MASK" , "epoch" , "1-10" , "Select epochs 1 to 10" );
  add_param( "MASK" , "mask-epoch" , "1-10" , "Mask epochs 1 to 10, leaving all others unchanged" );
  add_param( "MASK" , "sec" , "60-120" , "Select epochs overlapping this interval" );
  add_param( "MASK" , "hms" , "8:00-9:00" , "Select epochs overlapping this interval" );
  add_param( "MASK" , "dhms" , "d1-22:00:00-d1-23:00:00" , "Select epochs overlapping this date/time interval" );

  add_param( "MASK" , "random" , "20" , "Select 20 random (currently unmasked) epochs" );
  add_param( "MASK" , "random-from" , "20,N2,N3" , "Select random epochs from epochs containing one or more annotations" );

  add_param( "MASK" , "flip" , "" , "Reverse all masks" );
  add_param( "MASK" , "leading" , "W" , "Remove all leading epochs matching W" );
  add_param( "MASK" , "trailing" , "W" , "Remove all trailing epochs matching W" );
  add_param( "MASK" , "mask-leading" , "W" , "Mask leading epochs matching one or more annotations" );
  add_param( "MASK" , "mask-trailing" , "W" , "Mask trailing epochs matching one or more annotations" );
  add_param( "MASK" , "unmask-leading" , "W" , "Unmask epochs after trimming leading matches" );
  add_param( "MASK" , "unmask-trailing" , "W" , "Unmask epochs before trimming trailing matches" );
  add_param( "MASK" , "leading-trailing" , "W" , "Mask leading and trailing epochs matching one or more annotations" );
  add_param( "MASK" , "mask-leading-trailing" , "W" , "Same as leading-trailing" );
  add_param( "MASK" , "unmask-leading-trailing" , "W" , "Unmask the interior bounded by leading/trailing matches" );
  add_param( "MASK" , "flanked" , "REM,2" , "Select only REM epochs flanked by 2+ REM epochs before/after" );
  add_param( "MASK" , "regional" , "1,2" , "Retain epochs only if enough flanking epochs are unmasked" );
  add_param( "MASK" , "unmask-interior" , "" , "Unmask the contiguous interior between masked leading/trailing regions" );
  add_param( "MASK" , "first" , "100" , "Retain only the first N epochs" );
  add_param( "MASK" , "trim" , "W,10" , "Trim leading/trailing epochs with this annotation, optionally allowing N edge epochs" );
  add_param( "MASK" , "retain" , "N2,N3,REM" , "Retain only epochs containing one or more of these epoch annotations" );
  add_param( "MASK" , "stable-unique" , "10,N2,N3" , "Retain runs of at least X epochs with a stable unique annotation set" );
  add_param( "MASK" , "stable-any" , "10,N2,N3" , "Retain runs of at least X epochs containing any of these annotations" );
  
  // MASK / EMASK 
  
  add_table( "MASK" , "" , "Mask-operation summaries" );
  add_table( "MASK" , "EMASK" , "Mask-operation summaries by mask expression/label" );
  add_var( "MASK" , "", "N_MATCHES" , "Number of epochs that match the condition (e.g. having annotation A)");
  add_var( "MASK" , "", "N_MASK_SET" , "Number of previously unmasked epochs that were masked by this operation");
  add_var( "MASK" , "", "N_MASK_UNSET" , "Number of previously masked epochs that were unmasked by this operation");
  add_var( "MASK" , "", "N_UNCHANGED" , "Number of epochs whose mask status was not changed by this operation");
  add_var( "MASK" , "", "N_RETAINED" , "Number of epochs retained after this operation");
  add_var( "MASK" , "", "N_TOTAL" , "Total number of epochs");
  add_var( "MASK" , "", "MASK_MODE" , "Mask mode" );
  add_var( "MASK" , "", "MATCH_LOGIC" , "Match logic" );
  add_var( "MASK" , "", "MATCH_TYPE" , "Match type" );
  add_var( "MASK" , "EMASK", "N_MATCHES" , "Number of epochs that match the condition (e.g. having annotation A)");
  add_var( "MASK" , "EMASK", "N_MASK_SET" , "Number of previously unmasked epochs that were masked by this operation");
  add_var( "MASK" , "EMASK", "N_MASK_UNSET" , "Number of previously masked epochs that were unmasked by this operation");
  add_var( "MASK" , "EMASK", "N_UNCHANGED" , "Number of epochs whose mask status was not changed by this operation");
  add_var( "MASK" , "EMASK", "N_RETAINED" , "Number of epochs retained after this operation");
  add_var( "MASK" , "EMASK", "N_TOTAL" , "Total number of epochs");
  add_var( "MASK" , "EMASK", "MASK_MODE" , "Mask mode" );
  add_var( "MASK" , "EMASK", "MATCH_LOGIC" , "Match logic" );
  add_var( "MASK" , "EMASK", "MATCH_TYPE" , "Match type" );

  //
  // DUMP-MASK
  //

  add_cmd( "mask" , "DUMP-MASK" , "Dump the standard mask" );
  add_url( "DUMP-MASK" , "masks/#dump-mask" );
  add_verb( "DUMP-MASK" ,
	    "DUMP-MASK writes the current standard epoch mask at the epoch level. It can also create\n"
	    "a new annotation track that marks masked or unmasked epochs for later inspection.\n"
	    "\n"
	    "By default it writes epoch-level output; output can be suppressed if you only want the\n"
	    "derived annotation." );
  add_param( "DUMP-MASK" , "annot" , "MASK" , "Create an annotation marking masked or unmasked epochs" );
  add_param( "DUMP-MASK" , "annot-unmasked" , "T" , "When creating annot, mark unmasked rather than masked epochs" );
  add_param( "DUMP-MASK" , "output" , "F" , "Suppress standard output and only create the annotation" );

  add_table( "DUMP-MASK" , "E" , "Epoch-level mask tabulation" );
  add_var( "DUMP-MASK" , "E" , "EMASK" , "Mask status: 0=unmasked (included), 1=masked (excluded)" );

  //
  // RE
  //

  add_cmd( "mask"   , "RE" , "Restructure on standard mask" );
  add_url( "RE" , "masks/#restructure" );
  add_verb( "RE" ,
	    "RE physically restructures the in-memory EDF by removing data excluded by the current\n"
	    "standard epoch mask, then resetting the mask on the retained data.\n"
	    "\n"
	    "This converts logical masking into a smaller working dataset and can optionally preserve\n"
	    "caches and require a minimum number of epochs after restructuring." );
  add_param( "RE" , "verbose" , "" , "Verbose restructure logging" );
  add_param( "RE" , "preserve-cache" , "T" , "Preserve cached data across the restructure" );
  add_param( "RE" , "require" , "100" , "Set the empty flag if fewer than this many epochs remain" );

  add_table( "RE" , "" , "Restructured data duration" );
  add_var( "RE" , "" , "DUR1" , "Duration pre-restructuring (secs)");
  add_var( "RE" , "" , "DUR2" , "Duration post-restructuring (secs)");
  add_var( "RE" , "" , "NR1" , "Number of records pre-restructuring");
  add_var( "RE" , "" , "NR2" , "Number of records post-restructuring");
  add_var( "RE" , "" , "NA" , "Number of annotation channels after restructuring" );
  add_var( "RE" , "" , "NS" , "Number of data channels after restructuring" );

  //
  // RESTRUCTURE
  //

  add_cmd( "mask"   , "RESTRUCTURE" , "Restructure on standard mask" );
  add_url( "RESTRUCTURE" , "masks/#restructure" );
  add_verb( "RESTRUCTURE" ,
	    "RESTRUCTURE is the long-form alias for RE. It removes masked data from the in-memory EDF,\n"
	    "resets the standard mask on the retained data, and reports the pre/post record durations and\n"
	    "channel counts.\n"
	    "\n"
	    "Use this when downstream commands should operate only on retained data rather than consulting\n"
	    "the standard mask dynamically." );
  add_param( "RESTRUCTURE" , "verbose" , "" , "Verbose restructure logging" );
  add_param( "RESTRUCTURE" , "preserve-cache" , "T" , "Preserve cached data across the restructure" );
  add_param( "RESTRUCTURE" , "require" , "100" , "Set the empty flag if fewer than this many epochs remain" );

  add_table( "RESTRUCTURE" , "" , "Restructured data duration" );
  add_var( "RESTRUCTURE" , "" , "DUR1" , "Duration pre-restructuring (secs)");
  add_var( "RESTRUCTURE" , "" , "DUR2" , "Duration post-restructuring (secs)");
  add_var( "RESTRUCTURE" , "" , "NR1" , "Number of records pre-restructuring");
  add_var( "RESTRUCTURE" , "" , "NR2" , "Number of records post-restructuring");
  add_var( "RESTRUCTURE" , "" , "NA" , "Number of annotation channels after restructuring" );
  add_var( "RESTRUCTURE" , "" , "NS" , "Number of data channels after restructuring" );

  //
  // CHEP
  //

  add_cmd( "mask" , "CHEP" , "Channel-epoch masks" );
  add_url( "CHEP" , "masks/#chep" );
  add_verb( "CHEP" ,
	    "CHEP manages the channel-by-epoch mask, which tracks bad data at finer resolution than the\n"
	    "standard epoch mask. It can load or save CHEP files, mark whole channels as bad, collapse the\n"
	    "channel-epoch mask into the standard epoch mask, or promote heavily affected channels to whole-\n"
	    "channel masks or drops.\n"
	    "\n"
	    "Use sig to restrict the channels considered for collapse and reporting. dump writes the CHEP\n"
	    "grid to the standard output tables; otherwise CHEP still prints a console summary." );
  add_param( "CHEP" , "sig" , "C3,C4" , "Restrict CHEP operations and reporting to these signals" );
  add_param( "CHEP" , "clear" , "" , "Clear CHEP mask" );
  add_param( "CHEP" , "load" , "file.txt" , "Load CHEP from file.txt" );
  add_param( "CHEP" , "bad-channels" , "C3,C5" , "Manually specify bad channels" );
  add_param( "CHEP" , "epochs" , "0.1,2" , "Promote epochs to the standard mask if enough channels are CHEP-masked" );
  add_param( "CHEP" , "channels" , "0.5,10" , "Promote heavily masked channels to whole-channel CHEP masks" );
  add_param( "CHEP" , "drop-channels" , "0.5,10" , "Drop heavily masked channels from the EDF" );
  add_param( "CHEP" , "black-and-white" , "" , "When collapsing channels, also clear CHEP masks from channels that pass" );
  add_param( "CHEP" , "dump" , "" , "Write current CHEP mask to output" );
  add_param( "CHEP" , "save" , "file.txt" , "Write CHEP mask to file.txt" );

  // CHEP output....  
  add_table( "CHEP" , "CH" , "CHEP mask channel-wise summaries" );
  add_var( "CHEP" , "CH" , "CHEP" , "Masked epochs" );

  add_table( "CHEP" , "E" , "CHEP mask epoch-wise summaries" );
  add_var( "CHEP" , "E" , "CHEP" , "Masked channels" );

  add_table( "CHEP" , "CH,E" , "CHEP mask" );
  add_var( "CHEP" , "CH,E" , "CHEP" , "CHannel/EPoch mask" );

  //
  // CHEP-MASK
  //

  add_cmd( "artifact" , "CHEP-MASK" , "CHannel/EPoch outlier detection" );
  add_url( "CHEP-MASK" , "artifacts/#chep-mask" );
  add_verb( "CHEP-MASK" ,
            "Detect channel-by-epoch outliers and add them to the CHEP mask.\n"
            "CHEP-MASK computes Hjorth features per channel and epoch, then applies one or "
            "more iterative thresholds within-channel across epochs (ep-th), within-epoch "
            "across channels (ch-th), and/or globally across all channel-epoch pairs "
            "(chep-th). The *0 forms ignore any existing CHEP mask during that stage." );
  add_param( "CHEP-MASK" , "sig" , "C3,C4,F3,F4" , "Signals to include in CHEP outlier detection [required]" );
  add_param( "CHEP-MASK" , "ep-th" , "3,2.5" , "Iterative within-channel across-epoch Hjorth threshold(s)" );
  add_param( "CHEP-MASK" , "ep-th0" , "3,2.5" , "As ep-th, but ignore any existing CHEP mask during this stage" );
  add_param( "CHEP-MASK" , "ch-th" , "3" , "Iterative within-epoch across-channel Hjorth threshold(s)" );
  add_param( "CHEP-MASK" , "ch-th0" , "3" , "As ch-th, but ignore any existing CHEP mask during this stage" );
  add_param( "CHEP-MASK" , "chep-th" , "3" , "Iterative global channel-epoch Hjorth threshold(s)" );
  add_param( "CHEP-MASK" , "chep-th0" , "3" , "As chep-th, but ignore any existing CHEP mask during this stage" );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // FREEZES & CACHES
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // FREEZE
  //

  add_cmd( "freeze" , "FREEZE" , "Create a named freeze" );
  add_verb( "FREEZE" ,
	    "FREEZE stores a named snapshot of the current in-memory EDF so that later commands can\n"
	    "branch, experiment, and return to that state with THAW.\n"
	    "\n"
	    "The freeze name can be provided as a bare single argument or via tag=name." );
  add_param( "FREEZE" , "{tag}" , "baseline" , "Name of the freeze snapshot" );
  add_param( "FREEZE" , "tag" , "baseline" , "Name of the freeze snapshot" );

  //
  // THAW
  //

  add_cmd( "freeze" , "THAW" , "Restore a named freeze" );
  add_verb( "THAW" ,
	    "THAW replaces the current in-memory EDF with a previously frozen snapshot.\n"
	    "\n"
	    "The freeze name can be provided as a bare single argument or via tag=name. Options control\n"
	    "whether the snapshot is removed after thawing, whether caches are preserved, and whether a\n"
	    "missing freeze should halt immediately under strict mode." );
  add_param( "THAW" , "{tag}" , "baseline" , "Name of the freeze snapshot to restore" );
  add_param( "THAW" , "tag" , "baseline" , "Name of the freeze snapshot to restore" );
  add_param( "THAW" , "remove" , "T" , "Remove the freeze after restoring it" );
  add_param( "THAW" , "preserve-cache" , "T" , "Preserve caches when restoring the frozen EDF" );
  add_param( "THAW" , "strict" , "T" , "Halt if the requested freeze cannot be restored" );

  //
  // CLEAN-FREEZER
  //

  add_cmd( "freeze" , "CLEAN-FREEZER" , "Empty the freezer" );
  add_verb( "CLEAN-FREEZER" ,
	    "CLEAN-FREEZER removes all currently stored freeze snapshots.\n"
	    "\n"
	    "This only affects the internal freezer; it does not change the current in-memory EDF." );

  //
  // CACHE
  //

  add_cmd( "freeze" , "CACHE" , "Cache operations" );
  add_verb( "CACHE" ,
	    "CACHE manages Luna's internal caches. It can record selected command output into a named\n"
	    "cache, clear all caches, load or import caches from disk, or dump a chosen cache for\n"
	    "inspection.\n"
	    "\n"
	    "record configures writer-based caching of future output into a named cache. dump inspects an\n"
	    "existing cache of a specified type, while load/import restore cache contents from saved data." );
  add_param( "CACHE" , "cache" , "PSD1" , "Name of the cache to create, use, or import into" );
  add_param( "CACHE" , "record" , "PSD,PSD,CH,B" , "Record output from command,variable,{strata} into cache" );
  add_param( "CACHE" , "text" , "" , "When recording, use a text/string cache" );
  add_param( "CACHE" , "str" , "" , "Alias for text" );
  add_param( "CACHE" , "uppercase-keys" , "" , "Uppercase cached keys when recording" );
  add_param( "CACHE" , "clear" , "" , "Clear all caches" );
  add_param( "CACHE" , "load" , "cache.txt" , "Load caches from a previously saved cache file" );
  add_param( "CACHE" , "import" , "out.txt" , "Import cache values from long-format output" );
  add_param( "CACHE" , "factors" , "CH,B" , "When importing, specify which factors are present" );
  add_param( "CACHE" , "v" , "PSD" , "When importing, limit to these variables" );
  add_param( "CACHE" , "dump" , "" , "Dump a cache to stdout" );
  add_param( "CACHE" , "int" , "C1" , "Dump this integer cache" );
  add_param( "CACHE" , "num" , "C1" , "Dump this numeric cache" );
  add_param( "CACHE" , "tp" , "C1" , "Dump this time-point cache" );
  add_param( "CACHE" , "test" , "" , "Run internal cache test code" );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // CANONICAL SIGNALS
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // CANONICAL
  //

  add_cmd( "canon" , "CANONICAL" , "Canonicalize signals" );
  add_url( "CANONICAL" , "canonical/" );
  add_verb( "CANONICAL" ,
	    "CANONICAL maps EDF channels onto canonical signal definitions and standardizes the resulting\n"
	    "signals according to the current canonical rule engine.\n"
	    "\n"
	    "It reads one or more rule files from file, optionally restricts which canonical targets are\n"
	    "applied via inc/exc, and can run in dry-run mode with check. Options also control grouping,\n"
	    "whether non-canonical originals are dropped, whether prefiltering fields are retained, and\n"
	    "whether to dump the parsed rule templates without applying them." );
  add_param( "CANONICAL" , "file" , "csfile.txt" , "Canonical signal definition file" );
  add_param( "CANONICAL" , "inc" , "csEEG,csEOG" , "Only apply these canonical targets" );
  add_param( "CANONICAL" , "exc" , "csEMG" , "Exclude these canonical targets" );
  add_param( "CANONICAL" , "group" , "GRP1" , "Optional group filter applied within the rule file(s)" );
  add_param( "CANONICAL" , "prefix" , "defs/" , "Prefix to prepend to relative rule-file paths" );
  add_param( "CANONICAL" , "dump" , "" , "Dump parsed rule templates without applying them" );
  add_param( "CANONICAL" , "drop-originals" , "T" , "Drop non-canonical original signals after processing" );
  add_param( "CANONICAL" , "check" , "T" , "Dry-run mode: validate mappings without modifying the EDF" );
  add_param( "CANONICAL" , "mapper-util-mode" , "" , "Only check label mappings without full canonical processing" );
  add_param( "CANONICAL" , "verbose" , "" , "Verbose canonicalization logging" );
  add_param( "CANONICAL" , "prefiltering" , "T" , "Retain prefiltering header fields on canonical signals" );

  add_table( "CANONICAL" , "" , "Canonical signal summaries" );
  add_var( "CANONICAL" , "" , "CS_SET" , "Number of canonical signals set" );
  add_var( "CANONICAL" , "" , "CS_NOT" , "Number of canonical signals not set" );
  add_var( "CANONICAL" , "" , "USED_CH" , "Number of used EDF channels" );
  add_var( "CANONICAL" , "" , "UNUSED_CH" , "Number of unused EDF channels" );

  add_table( "CANONICAL" , "CS" , "Canonical signal information" );
  add_var( "CANONICAL" , "CS" , "DEFINED" , "Is canonical signal present/defined?" );
  add_var( "CANONICAL" , "CS" , "SIG" , "Primary signal" );
  add_var( "CANONICAL" , "CS" , "REF" , "Reference signal" );

  add_table( "CANONICAL" , "CH" , "EDF channel information" );
  add_var( "CANONICAL" , "CH" , "DROPPED" , "Original channel dropped" );
  add_var( "CANONICAL" , "CH" , "USED" , "0/1 indicator that this channel was used in constructing canonical signals" );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // MANIPULATIONS
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // SIGNALS
  //

  add_cmd( "manip" , "SIGNALS" , "Retain/remove specific EDF channels" );
  add_url( "SIGNALS" , "manipulations/#signals" );
  add_verb( "SIGNALS" ,
	    "SIGNALS keeps, drops, or picks channels from the current EDF. It can also require that\n"
	    "certain channels exist, or rename a picked channel after resolving alternatives.\n"
	    "\n"
	    "Use keep/drop for explicit inclusion or exclusion, pick to choose the first available\n"
	    "channel from an ordered list, and req when missing requested channels should flag the EDF." );
  add_param( "SIGNALS" , "drop" , "EMG,ECG" , "Drop channels EMG and ECG" );
  add_param( "SIGNALS" , "keep" , "C3,C4" , "Drop all channels except C3 and C4" );
  add_param( "SIGNALS" , "pick" , "C3,C3-M2,C4" , "Pick the first available channel from this ordered list" );
  add_param( "SIGNALS" , "req" , "C3,C4" , "Require these channels to exist, otherwise flag the EDF" );
  add_param( "SIGNALS" , "rename" , "EEG" , "Rename the picked channel to this label" );

  //
  // COPY
  //

  add_cmd( "manip" , "COPY" , "Duplicate one or more EDF channels" );
  add_url( "COPY" , "manipulations/#copy" );
  add_verb( "COPY" ,
	    "COPY duplicates one or more signals to new channel labels.\n"
	    "\n"
	    "Use tag or pretag to derive new labels from each source channel, or new to provide one\n"
	    "explicit destination label when copying a single input signal." );
  add_param( "COPY" , "sig" , "C3,C4" , "List of channels to duplicate" );
  add_param( "COPY" , "tag" , "_COPY" , "Append this tag to each copied channel label" );
  add_param( "COPY" , "pretag" , "cp_" , "Prepend this tag to each copied channel label" );
  add_param( "COPY" , "new" , "C3_COPY" , "Explicit new label, only valid when copying one signal" );
  add_param( "COPY" , "silent" , "" , "Suppress copy-operation logging" );

  //
  // RESAMPLE
  //

  add_cmd( "manip" , "RESAMPLE" , "Resample signal(s)" );
  add_url( "RESAMPLE" , "manipulations/#resample" );
  add_verb( "RESAMPLE" ,
	    "RESAMPLE changes the sampling rate of one or more signals.\n"
	    "\n"
	    "It supports multiple conversion methods, can be restricted to downsampling only, and can\n"
	    "optionally allow upsampling only above a given original sample-rate threshold." );
  add_param( "RESAMPLE" , "sig" , "C3,C4" , "List of channels to resample" );
  add_param( "RESAMPLE" , "sr" , "200" , "New sampling rate (Hz) [required]" );
  add_param( "RESAMPLE" , "downsample" , "" , "Only downsample; do not upsample lower-rate signals" );
  add_param( "RESAMPLE" , "upsample-if" , "100" , "Allow upsampling only for signals already at or above this rate" );
  add_param( "RESAMPLE" , "best" , "" , "Use the best-quality sinc resampler" );
  add_param( "RESAMPLE" , "medium" , "" , "Use the medium-quality sinc resampler" );
  add_param( "RESAMPLE" , "fastest" , "" , "Use the fastest sinc resampler" );
  add_param( "RESAMPLE" , "linear" , "" , "Use linear interpolation" );
  add_param( "RESAMPLE" , "zoh" , "" , "Use zero-order hold interpolation" );
  add_param( "RESAMPLE" , "method" , "medium" , "Resampling method: best, medium, fastest, linear or zoh" );

  //
  // REFERENCE
  //

  add_cmd( "manip" , "REFERENCE" , "Re-reference signals" );
  add_url( "REFERENCE" , "manipulations/#reference" );
  add_verb( "REFERENCE" ,
	    "REFERENCE subtracts one or more reference channels from the target signals.\n"
	    "\n"
	    "It can either update signals in place or create new referenced channels, and pairwise mode\n"
	    "supports one-to-one referencing across matching signal/reference lists." );
  add_param( "REFERENCE" , "sig" , "C3,C4" , "List of signals to re-reference" );
  add_param( "REFERENCE" , "ref" , "A1,A2" , "Signal(s) providing the reference [required]" );
  add_param( "REFERENCE" , "new" , "C3_A2" , "Create new referenced channel(s) with these label(s)" );
  add_param( "REFERENCE" , "pairwise" , "" , "Apply references pairwise across sig and ref lists" );
  add_param( "REFERENCE" , "sr" , "200" , "Sampling rate for newly created referenced channels" );

  //
  // DEREFERENCE
  //

  add_cmd( "manip" , "DEREFERENCE" , "De-reference signals" );
  add_url( "DEREFERENCE" , "manipulations/#dereference" );
  add_verb( "DEREFERENCE" ,
	    "DEREFERENCE removes a previously applied reference by adding the reference signal back.\n"
	    "\n"
	    "The parameter surface mirrors REFERENCE: operate in place or create new channels, and use\n"
	    "pairwise mode when sig and ref should be matched one-to-one." );
  add_param( "DEREFERENCE" , "sig" , "C3-A2" , "List of de-referenced signals" );
  add_param( "DEREFERENCE" , "ref" , "A2" , "Reference signal(s) to add back [required]" );
  add_param( "DEREFERENCE" , "new" , "C3" , "Create new de-referenced channel(s) with these label(s)" );
  add_param( "DEREFERENCE" , "pairwise" , "" , "Apply references pairwise across sig and ref lists" );
  add_param( "DEREFERENCE" , "sr" , "200" , "Sampling rate for newly created de-referenced channels" );

  //
  // uV
  //

  add_cmd( "manip" , "uV" , "Rescale units to uV" );
  add_url( "uV" , "manipulations/#uv" );
  add_verb( "uV" ,
	    "uV rescales one or more signals so their physical units are microvolts.\n"
	    "\n"
	    "Only standard voltage conversions supported by the EDF header units are applied." );
  add_param( "uV" , "sig" , "C3,C4" , "List of signals to convert" );

  //
  // mV
  //

  add_cmd( "manip" , "mV" , "Rescale units to mV" );
  add_url( "mV" , "manipulations/#mv" );
  add_verb( "mV" ,
	    "mV rescales one or more signals so their physical units are millivolts.\n"
	    "\n"
	    "Only standard voltage conversions supported by the EDF header units are applied." );
  add_param( "mV" , "sig" , "C3,C4" , "List of signals to convert" );

  //
  // FLIP
  //

  add_cmd( "manip" , "FLIP" , "Flip signal polarity" );
  add_url( "FLIP" , "manipulations/#flip" );
  add_verb( "FLIP" ,
	    "FLIP multiplies one or more signals by -1, reversing their polarity.\n"
	    "\n"
	    "The command reports which channels were flipped." );
  add_param( "FLIP" , "sig" , "C3,C4" , "List of signals to flip" );

  add_table( "FLIP" , "CH" , "Tracking flipped channels" );
  add_var( "FLIP" , "CH" , "FLIP" , "Channel flipped" );

  //
  // ZC
  //

  add_cmd( "manip" , "ZC" , "Mean-center signals" );
  add_url( "ZC" , "manipulations/#zc" );
  add_verb( "ZC" ,
	    "ZC removes the mean from one or more signals.\n"
	    "\n"
	    "By default it works across the full signal; add epoch to mean-center each epoch separately." );
  add_param( "ZC" , "sig" , "C3,C4" , "List of signals to mean-center" );
  add_param( "ZC" , "epoch" , "" , "Mean-center each epoch separately" );

  //
  // ROBUST-NORM
  //

  add_cmd( "manip" , "ROBUST-NORM" , "Robust normalization" );
  add_url( "ROBUST-NORM" , "manipulations/#robust-norm" );
  add_verb( "ROBUST-NORM" ,
	    "ROBUST-NORM robustly centers and scales signals, optionally epoch-wise.\n"
	    "\n"
	    "It supports simple de-meaning, IQR normalization, winsorization, a second standardization\n"
	    "pass, and independent control over centering and scaling." );
  add_param( "ROBUST-NORM" , "sig" , "C3,C4" , "List of signals to normalize" );
  add_param( "ROBUST-NORM" , "epoch" , "" , "Normalize each epoch separately" );
  add_param( "ROBUST-NORM" , "simple-demean" , "" , "Only subtract the mean, without scaling" );
  add_param( "ROBUST-NORM" , "IQR" , "" , "Use IQR-based normalization" );
  add_param( "ROBUST-NORM" , "center" , "F" , "Control robust centering" );
  add_param( "ROBUST-NORM" , "scale" , "F" , "Control robust scaling" );
  add_param( "ROBUST-NORM" , "winsor" , "3" , "Winsorize at this threshold before optional second normalization" );
  add_param( "ROBUST-NORM" , "second-norm" , "T" , "Apply a second non-robust normalization after winsorization" );
  add_param( "ROBUST-NORM" , "silent" , "T" , "Suppress normalization logging" );

  //
  // COMBINE
  //

  add_cmd( "manip" , "COMBINE" , "Combine channels" );
  add_url( "COMBINE" , "manipulations/#combine" );
  add_verb( "COMBINE" ,
	    "COMBINE creates a new signal from two or more input channels.\n"
	    "\n"
	    "Exactly one of sum, mean, or median names the new output channel. allow-missing relaxes the\n"
	    "requirement that all requested input channels be present." );
  add_param( "COMBINE" , "sig" , "C3,C4" , "Input channels to combine" );
  add_param( "COMBINE" , "sum" , "EEG_SUM" , "Create this new channel as the sum of sig" );
  add_param( "COMBINE" , "mean" , "EEG_MEAN" , "Create this new channel as the mean of sig" );
  add_param( "COMBINE" , "median" , "EEG_MED" , "Create this new channel as the median of sig" );
  add_param( "COMBINE" , "allow-missing" , "" , "Allow only a subset of requested channels to be present" );

  //
  // SCALE
  //

  add_cmd( "manip" , "SCALE" , "Rescale a channel" );
  add_url( "SCALE" , "manipulations/#scale" );
  add_verb( "SCALE" ,
	    "SCALE updates one or more signals using explicit min/max scaling and/or clipping.\n"
	    "\n"
	    "Use min-max to set a new scale range, clip-min and clip-max to clamp values, or combine\n"
	    "those operations in one pass." );
  add_param( "SCALE" , "sig" , "C3,C4" , "Signals to rescale" );
  add_param( "SCALE" , "min-max" , "-100,100" , "Set a new min/max scale for the signal" );
  add_param( "SCALE" , "clip-min" , "-500" , "Clip values below this threshold" );
  add_param( "SCALE" , "clip-max" , "500" , "Clip values above this threshold" );

  //
  // SHIFT
  //

  add_cmd( "manip" , "SHIFT" , "Shift a signal" );
  add_url( "SHIFT" , "manipulations/#shift" );
  add_verb( "SHIFT" ,
	    "SHIFT circularly shifts one or more signals by a given number of sample points.\n"
	    "\n"
	    "Use no-wrap to shift with zero filling instead of wrapping values around the ends." );
  add_param( "SHIFT" , "sig" , "C3,C4" , "Signals to shift" );
  add_param( "SHIFT" , "sp" , "100" , "Shift amount in sample points" );
  add_param( "SHIFT" , "no-wrap" , "" , "Do not wrap shifted samples around the signal edges" );

  //
  // SCRAMBLE
  //

  add_cmd( "manip" , "SCRAMBLE" , "Scramble a signal" );
  add_url( "SCRAMBLE" , "manipulations/#scramble" );
  add_verb( "SCRAMBLE" ,
	    "SCRAMBLE randomly permutes all sample values within a signal.\n"
	    "\n"
	    "The operation is applied independently to each selected signal across the whole trace." );
  add_param( "SCRAMBLE" , "sig" , "C3,C4" , "Signals to scramble" );

  //
  // RECORD-SIZE
  //

  add_cmd( "manip" , "RECORD-SIZE" , "Change EDF record size" );
  add_url( "RECORD-SIZE" , "manipulations/#record-size" );
  add_verb( "RECORD-SIZE" ,
	    "RECORD-SIZE rewrites the in-memory EDF to a new record duration and then writes the result\n"
	    "to disk via the normal WRITE machinery.\n"
	    "\n"
	    "This is primarily a file-writing operation, and unless no-problem is set it marks the EDF as\n"
	    "done so processing continues with the next record in a sample list." );
  add_param( "RECORD-SIZE" , "dur" , "1" , "New EDF record/block size" );
  add_param( "RECORD-SIZE" , "edf-dir" , "edfs/" , "Folder for writing new EDFs" );
  add_param( "RECORD-SIZE" , "edf-tag" , "rec1" , "Tag added to new EDFs" );
  add_param( "RECORD-SIZE" , "sample-list" , "s2.lst" , "Generate a sample-list pointing to the new EDFs" );
  add_param( "RECORD-SIZE" , "no-problem" , "" , "Do not set the problem flag after writing the new EDF" );

  add_table( "RECORD-SIZE" , "" , "Restructured data duration" );
  add_var( "RECORD-SIZE", "" , "NR1" , "Pre-restructure number of records" );
  add_var( "RECORD-SIZE", "" , "NR2" , "Post-restructure number of records" );
  add_var( "RECORD-SIZE", "" , "DUR1" , "Pre-restructure duration (seconds)" );
  add_var( "RECORD-SIZE", "" , "DUR2" , "Post-restructure duration (seconds)" );

  //
  // EDF-MINUS
  //

  add_cmd( "manip" , "EDF-MINUS" , "Realign EDF records" );
  add_url( "EDF-MINUS" , "manipulations/#edf-minus" );
  add_verb( "EDF-MINUS" ,
	    "EDF-MINUS writes a new EDF after retaining and optionally realigning selected segments,\n"
	    "signals, annotations, and epochs.\n"
	    "\n"
	    "It supports annotation-based requirements and alignments, segment retention policies, signal\n"
	    "subsetting, output filename control, and optional clock-time formatting for written annotations." );
  add_param( "EDF-MINUS" , "out" , "subj1-minus" , "Output filename root [required]" );
  add_param( "EDF-MINUS" , "sig" , "C3,C4" , "Restrict output to these signals" );
  add_param( "EDF-MINUS" , "max-sr" , "1024" , "Maximum allowed sample rate in the output" );
  add_param( "EDF-MINUS" , "align" , "N1,N2" , "Realign retained segments to these annotations" );
  add_param( "EDF-MINUS" , "unaligned" , "" , "Do not perform annotation-based alignment" );
  add_param( "EDF-MINUS" , "dur" , "30" , "Alignment unit in seconds" );
  add_param( "EDF-MINUS" , "require" , "N2,N3,REM" , "Require retained segments to overlap these annotations" );
  add_param( "EDF-MINUS" , "require-whole" , "N2,N3" , "Require whole-annotation overlap for retained segments" );
  add_param( "EDF-MINUS" , "require-dur" , "30" , "Require at least this much overlap in seconds" );
  add_param( "EDF-MINUS" , "prefix" , "A_" , "Prefix applied to written annotation labels" );
  add_param( "EDF-MINUS" , "policy" , "all" , "Segment inclusion policy" );
  add_param( "EDF-MINUS" , "segments" , "all" , "Which retained segments to keep in the output" );
  add_param( "EDF-MINUS" , "hms" , "" , "Write annotations using clock-time formatting" );
  add_param( "EDF-MINUS" , "dhms" , "" , "Write annotations using date/time formatting" );
  add_table( "EDF-MINUS" , "" , "EDF-MINUS summary" );
  add_var( "EDF-MINUS" , "" , "OKAY" , "0/1 indicator that EDF-MINUS completed successfully" );
  add_var( "EDF-MINUS" , "" , "MSG" , "Failure message if EDF-MINUS could not proceed" );
  add_table( "EDF-MINUS" , "SEG" , "Segment-level EDF-MINUS output" );
  add_var( "EDF-MINUS" , "SEG" , "ORIG" , "Original segment interval" );
  add_var( "EDF-MINUS" , "SEG" , "EDIT" , "Edited segment interval" );
  add_var( "EDF-MINUS" , "SEG" , "INCLUDED" , "0/1 indicator that this segment was retained" );
  add_var( "EDF-MINUS" , "SEG" , "DUR_ORIG" , "Original segment duration (seconds)" );
  add_var( "EDF-MINUS" , "SEG" , "DUR_EDIT" , "Edited segment duration (seconds)" );
  add_table( "EDF-MINUS" , "SEG,ANNOT" , "Per-segment annotation counts" );
  add_var( "EDF-MINUS" , "SEG,ANNOT" , "N_ALIGN" , "Number of alignment annotations in this segment" );
  add_var( "EDF-MINUS" , "SEG,ANNOT" , "N_REQ" , "Number of required annotations in this segment" );
  add_var( "EDF-MINUS" , "SEG,ANNOT" , "N_ALL" , "Total number of this annotation in this segment" );

  //
  // ANON
  //

  add_cmd( "manip" , "ANON" , "Strip ID information" );
  add_url( "ANON" , "manipulations/#anon" );
  add_verb( "ANON" ,
	    "ANON anonymizes the EDF header by clearing or replacing identifying fields.\n"
	    "\n"
	    "By default it nulls the EDF ID and resets the start date. insert-id uses the sample-list ID,\n"
	    "and root generates sequential anonymized IDs under a supplied prefix." );
  add_param( "ANON" , "insert-id" , "" , "Set the EDF ID from the sample-list ID" );
  add_param( "ANON" , "root" , "anon" , "Generate sequential anonymized IDs using this root" );

  //
  // SET-HEADERS
  //

  add_cmd( "manip" , "SET-HEADERS" , "Set EDF headers" );
  add_url( "SET-HEADERS" , "manipulations/#set-headers" );
  add_verb( "SET-HEADERS" ,
	    "SET-HEADERS directly assigns selected EDF header values, either globally or for a selected\n"
	    "set of signals.\n"
	    "\n"
	    "Use sig to target per-channel header fields such as transducer, unit, and prefiltering;\n"
	    "record-level fields such as id and start-date apply to the EDF as a whole." );
  add_param( "SET-HEADERS" , "id" , "subj1" , "Set the EDF patient/recording ID" );
  add_param( "SET-HEADERS" , "recording-info" , "Startdate X X X X" , "Set the EDF recording-info field" );
  add_param( "SET-HEADERS" , "start-date" , "01.01.85" , "Set the EDF start date" );
  add_param( "SET-HEADERS" , "start-time" , "22.00.00" , "Set the EDF start time" );
  add_param( "SET-HEADERS" , "sig" , "C3,C4" , "Signals whose per-channel header fields should be updated" );
  add_param( "SET-HEADERS" , "transducer" , "EEG" , "Set the transducer type for sig" );
  add_param( "SET-HEADERS" , "physical-dimension" , "uV" , "Set the physical dimension for sig" );
  add_param( "SET-HEADERS" , "unit" , "uV" , "Alias for physical-dimension" );
  add_param( "SET-HEADERS" , "prefiltering" , "0.3-35Hz" , "Set the prefiltering field for sig" );

  //
  // SET-VAR
  //

  add_cmd( "manip" , "SET-VAR" , "Set Luna variables" );
  add_url( "SET-VAR" , "manipulations/#set-var" );
  add_verb( "SET-VAR" ,
	    "SET-VAR evaluates a single expression and stores the result as an individual-level Luna\n"
	    "variable for the current EDF.\n"
	    "\n"
	    "The command takes a single key=value assignment, where the right-hand side is parsed as an\n"
	    "expression and the left-hand side is the variable name to store." );
  add_param( "SET-VAR" , "{var=expr}" , "AGE=45" , "Set an individual-level variable from an expression" );

  //
  // SET-TIMESTAMPS
  //

  add_cmd( "manip" , "SET-TIMESTAMPS" , "Set EDF time-stamps" );
  add_url( "SET-TIMESTAMPS" , "manipulations/#set-timestamps" );
  add_verb( "SET-TIMESTAMPS" ,
	    "SET-TIMESTAMPS replaces the EDF record time-stamps from a file containing one elapsed-time\n"
	    "value per EDF record.\n"
	    "\n"
	    "This converts the dataset to EDF+D timing internally, rebuilds the discontinuous timeline,\n"
	    "and regenerates the in-memory time-track annotations." );
  add_param( "SET-TIMESTAMPS" , "file" , "timestamps.txt" , "File containing one elapsed-time value per EDF record" );

  //
  // RENAME
  //

  add_cmd( "manip" , "RENAME" , "Rename channels" );
  add_url( "RENAME" , "manipulations/#rename" );
  add_verb( "RENAME" ,
	    "RENAME changes one or more signal labels.\n"
	    "\n"
	    "Mappings can come from sig/new on the command line or from a two-column tab-delimited file." );
  add_param( "RENAME" , "sig" , "C3,C4" , "Signals to rename" );
  add_param( "RENAME" , "new" , "EEG1,EEG2" , "New labels corresponding to sig" );
  add_param( "RENAME" , "file" , "rename.tsv" , "Two-column tab-delimited old/new rename file" );

  //
  // ENFORCE-SR
  //

  add_cmd( "manip" , "ENFORCE-SR" , "Require a sample rate" );
  add_url( "ENFORCE-SR" , "manipulations/#enforce-sr" );
  add_verb( "ENFORCE-SR" ,
	    "ENFORCE-SR drops signals that do not satisfy sample-rate constraints.\n"
	    "\n"
	    "It can enforce compatibility with a requested EDF record duration and/or keep only channels\n"
	    "whose sample rates fall within a specified range." );
  add_param( "ENFORCE-SR" , "sig" , "C3,C4" , "Signals to test against the sample-rate constraints" );
  add_param( "ENFORCE-SR" , "dur" , "1" , "Require signals to fit an integer number of samples into this record duration" );
  add_param( "ENFORCE-SR" , "sr" , "100,256" , "Only retain signals with sample rates inside this range" );

  //
  // MINMAX
  //

  add_cmd( "manip" , "MINMAX" , "Set digital/physical min/max" );
  add_url( "MINMAX" , "manipulations/#minmax" );
  add_verb( "MINMAX" ,
	    "MINMAX harmonizes physical min/max ranges across signals, with optional clipping.\n"
	    "\n"
	    "With no min/max values it standardizes PMIN/PMAX across the selected channels; with min and/or\n"
	    "max it clips and updates the header ranges, with force extending ranges when needed." );
  add_param( "MINMAX" , "sig" , "C3,C4" , "Signals whose min/max ranges should be harmonized" );
  add_param( "MINMAX" , "min" , "-500" , "Clip or set the physical minimum" );
  add_param( "MINMAX" , "max" , "500" , "Clip or set the physical maximum" );
  add_param( "MINMAX" , "force" , "" , "Force expansion of ranges rather than only clipping inward" );

  //
  // TIME-TRACK
  //

  add_cmd( "manip" , "TIME-TRACK" , "Add a time-track signal" );
  add_url( "TIME-TRACK" , "manipulations/#time-track" );
  add_verb( "TIME-TRACK" ,
	    "TIME-TRACK adds an EDF+ time-track annotation channel to the current EDF.\n"
	    "\n"
	    "This is useful before writing EDF+ output when an explicit time-track is required." );

  //
  // SCALE
  //

  //
  // RECTIFY
  //

  add_cmd( "manip" , "RECTIFY" , "Rectify signals" );
  add_url( "RECTIFY" , "manipulations/#rectify" );
  add_verb( "RECTIFY" ,
	    "RECTIFY replaces each sample with its absolute value.\n"
	    "\n"
	    "The operation is applied independently to each selected signal." );
  add_param( "RECTIFY" , "sig" , "C3,C4" , "Signals to rectify" );

  //
  // REVERSE
  //

  add_cmd( "manip" , "REVERSE" , "Reverse a signal" );
  add_url( "REVERSE" , "manipulations/#reverse" );
  add_verb( "REVERSE" ,
	    "REVERSE flips one or more signals in the time domain.\n"
	    "\n"
	    "The command reports which channels were reversed." );
  add_param( "REVERSE" , "sig" , "C3,C4" , "Signals to reverse in time" );
  add_table( "REVERSE" , "CH" , "Tracking reversed channels" );
  add_var( "REVERSE" , "CH" , "REVERSE" , "Channel reversed" );

  //
  // MOVING-AVERAGE
  //

  add_cmd( "manip" , "MOVING-AVERAGE" , "Moving average or median" );
  add_url( "MOVING-AVERAGE" , "manipulations/#moving-average" );
  add_verb( "MOVING-AVERAGE" ,
	    "MOVING-AVERAGE smooths signals using a centered moving average, moving median, or triangular\n"
	    "kernel.\n"
	    "\n"
	    "The half-window is specified in seconds. Add epoch to smooth each epoch independently rather\n"
	    "than across the full signal." );
  add_param( "MOVING-AVERAGE" , "sig" , "C3,C4" , "Signals to smooth" );
  add_param( "MOVING-AVERAGE" , "hw" , "0.5" , "Half-window size in seconds" );
  add_param( "MOVING-AVERAGE" , "median" , "" , "Use a moving median instead of a mean" );
  add_param( "MOVING-AVERAGE" , "tri" , "" , "Use a triangular moving average" );
  add_param( "MOVING-AVERAGE" , "lwr" , "0" , "Triangular-kernel lower-bound parameter" );
  add_param( "MOVING-AVERAGE" , "epoch" , "" , "Apply smoothing separately within each epoch" );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // ALIGNMENT
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // EDF
  //

  // add_cmd( "align" , "EDF" , "Convert EDF+D data to EDF/EDF+C by dropping discontinuities" );

  //
  // ALIGN-EPOCHS
  //

  // add_cmd( "align" , "ALIGN-EPOCHS" , "Align epochs to an external epoch scheme" );

  //
  // ALIGN-ANNOTS
  //

  // add_cmd( "align" , "ALIGN-ANNOTS" , "Align annotations to the current EDF timeline" );

  //
  // INSERT
  //

  // add_cmd( "align" , "INSERT" , "Insert records or gaps into an EDF timeline" );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // OUTPUTS
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // WRITE
  //

  add_cmd( "output" , "WRITE" , "Write a new EDF file" );
  add_url( "WRITE" , "outputs/#write" );
  add_verb( "WRITE" ,
            "Write the current in-memory EDF to disk after any prior masking, restructuring, "
            "channel changes, or annotation edits.\n"
            "This command can write a renamed copy, optionally append it to a sample list, "
            "and force EDF/EDF+C style output when needed." );
  add_param( "WRITE" , "edf" , "subj1_v2" , "Explicit output EDF root name (written as subj1_v2.edf)" );
  add_param( "WRITE" , "edf-dir" , "edfs/" , "Output folder for the written EDF" );
  add_param( "WRITE" , "edf-tag" , "v2" , "Append this tag to the existing EDF root name" );
  add_param( "WRITE" , "edfz" , "T" , "Write a compressed .edf.gz / EDFZ-style file" );
  add_param( "WRITE" , "sample-list" , "v2.lst" , "Append the written EDF to this sample list" );
  add_param( "WRITE" , "with-annots" , "" , "When appending to sample-list, also append linked annotation files" );
  add_param( "WRITE" , "force-edf" , "" , "Force EDF/EDF+C style output after restructuring" );
  add_param( "WRITE" , "null-starttime" , "" , "With force-edf, zero the start time in the written EDF" );
  add_param( "WRITE" , "EDF+D" , "" , "Preserve EDF+D style output even for quasi-continuous data" );
  add_param( "WRITE" , "channels" , "C3,C4,EOG-L" , "Write only these channels, in this order" );
  
  add_table( "WRITE", "" , "Summary from the implicit pre-WRITE restructure step" );
  add_var( "WRITE", "" , "NR1" , "Number of EDF records before restructuring" );
  add_var( "WRITE", "" , "NR2" , "Number of EDF records after restructuring" );
  add_var( "WRITE", "" , "DUR1" , "Duration before restructuring (seconds)" );
  add_var( "WRITE", "" , "DUR2" , "Duration after restructuring (seconds)" );

  //
  // MATRIX
  //

  add_cmd( "output" , "MATRIX" , "Dump signals to a file" );
  add_url( "MATRIX" , "outputs/#matrix" );
  add_verb( "MATRIX" ,
            "Write one or more uniformly sampled signals to a tabular matrix file, optionally "
            "with epoch-level or interval annotation indicators.\n"
            "This is intended for exporting aligned signal matrices for downstream analysis rather "
            "than for writing a new EDF." );
  add_param( "MATRIX" , "file" , "signals.txt" , "Output filename [required]" );
  add_param( "MATRIX" , "sig" , "C3,C4" , "Signals to export; all selected signals must share a sampling rate" );
  add_param( "MATRIX" , "annot" , "N2,SPINDLE" , "Add requested interval or epoch annotations to the matrix output" );
  add_param( "MATRIX" , "head" , "10" , "Only write the first N unmasked epochs" );
  add_param( "MATRIX" , "min" , "" , "Minimal output with signals only" );
  add_param( "MATRIX" , "minimal" , "" , "Alias for min" );
  add_param( "MATRIX" , "format2" , "" , "Write the alternate multi-row matrix format" );

  //
  // DUMP-RECORDS
  //

  add_cmd( "output" , "DUMP-RECORDS" , "Dump records to standard output" );
  add_url( "DUMP-RECORDS" , "outputs/#dump-records" );
  add_verb( "DUMP-RECORDS" ,
            "Print a verbose record-by-record dump of annotations and sample values.\n"
            "This is primarily a diagnostic view of the current EDF timeline rather than a "
            "structured writer output command." );
  add_param( "DUMP-RECORDS" , "no-signals" , "" , "Do not show signal data" );
  add_param( "DUMP-RECORDS" , "no-annots" , "" , "Do not show annotation information" );

  //
  // RECS
  //

  add_cmd( "output" , "RECS" , "Dump EDF record structure" );
  add_url( "RECS" , "outputs/#recs" );
  add_verb( "RECS" ,
            "Print a compact record-level view of the EDF timeline, including retained record "
            "indices, intervals, and any epochs spanned by each record." );

  //
  // SEGMENTS
  //

  add_cmd( "output" , "SEGMENTS" , "Report on contiguous segments in an EDF/EDF+" );
  add_url( "SEGMENTS" , "outputs/#segments" );
  add_verb( "SEGMENTS" ,
            "Summarize contiguous retained data segments and gaps in the current EDF timeline.\n"
            "For EDF+D files this reports each continuous run of records, and it can optionally "
            "add segment/gap annotations or flag the largest segments." );
  add_param( "SEGMENTS" , "annot" , "" , "Create DISC-SEGMENT and DISC-GAP annotations" );
  add_param( "SEGMENTS" , "largest" , "LARGEST_SEG" , "Create an annotation marking the largest segment(s)" );
  add_param( "SEGMENTS" , "requires-min" , "10" , "With largest, only annotate segments at least this many minutes long" );

  add_table( "SEGMENTS" , "" , "Segment and gap summary" );
  add_var( "SEGMENTS" , "" , "NSEGS" , "Number of contiguous segments" );
  add_var( "SEGMENTS" , "" , "NGAPS" , "Number of discontinuity gaps" );
  add_var( "SEGMENTS" , "" , "LARGEST" , "Duration of the largest segment (seconds)" );

  add_table( "SEGMENTS" , "SEG" , "Information on each segment" );
  add_var( "SEGMENTS" , "SEG" , "DUR_HR" , "Segment duration (hours)" );
  add_var( "SEGMENTS" , "SEG" , "DUR_MIN" , "Segment duration (minutes)" );
  add_var( "SEGMENTS" , "SEG" , "DUR_SEC" , "Segment duration (seconds)" );

  add_var( "SEGMENTS" , "SEG" , "START" , "Segment start (seconds)" );
  add_var( "SEGMENTS" , "SEG" , "START_HMS" , "Segment start (hh:mm:ss)" );
    
  add_var( "SEGMENTS" , "SEG" , "STOP" , "Segment stop (seconds)" );
  add_var( "SEGMENTS" , "SEG" , "STOP_HMS" , "Segment stop (hh:mm:ss)" );

  add_table( "SEGMENTS" , "GAP" , "Information on each discontinuity gap" );
  add_var( "SEGMENTS" , "GAP" , "DUR_HR" , "Gap duration (hours)" );
  add_var( "SEGMENTS" , "GAP" , "DUR_MIN" , "Gap duration (minutes)" );
  add_var( "SEGMENTS" , "GAP" , "DUR_SEC" , "Gap duration (seconds)" );
  add_var( "SEGMENTS" , "GAP" , "START" , "Gap start (seconds)" );
  add_var( "SEGMENTS" , "GAP" , "START_HMS" , "Gap start (hh:mm:ss)" );
  add_var( "SEGMENTS" , "GAP" , "STOP" , "Gap stop (seconds)" );
  add_var( "SEGMENTS" , "GAP" , "STOP_HMS" , "Gap stop (hh:mm:ss)" );

  //
  // HEAD
  //

  add_cmd( "output" , "HEAD" , "Show a short interval of signal data" );
  add_url( "HEAD" , "outputs/#head" );
  add_verb( "HEAD" ,
            "Print a small sample-by-sample window from one or more signals, usually from a "
            "single epoch.\n"
            "This is a quick inspection tool for checking signal values and timing, not a "
            "structured export format." );
  add_param( "HEAD" , "sig" , "C3,C4" , "Signals to display; all selected signals must share a sampling rate" );
  add_param( "HEAD" , "epoch" , "1" , "Epoch number to display (1-based)" );
  add_param( "HEAD" , "sec" , "2" , "Only show this many seconds from the start of the epoch" );

  //
  // SEDF
  //

  add_cmd( "output" , "SEDF" , "Generate a summary EDF" );
  add_url( "SEDF" , "outputs/#sedf" );
  add_verb( "SEDF" ,
            "Create a per-epoch summary EDF from one or more input signals.\n"
            "For EEG-like channels this writes Hjorth activity, mobility, and complexity as new "
            "summary channels; for other channel types it writes per-epoch mean, minimum, and "
            "maximum. The output is a continuous .sedf file with one sample per epoch." );
  add_param( "SEDF" , "sig" , "C3,C4,EMG" , "Signals to summarize [required]" );
  add_param( "SEDF" , "sedf-dir" , "sedf/" , "Output folder for written .sedf files" );
  add_param( "SEDF" , "sample-list" , "sedf.lst" , "Append the written .sedf file to this sample list" );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // FIR FILTERS
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // --fir
  //

  add_cmd( "helpers" , "--fir" , "Design FIR (--fir-design)" );
  add_param( "--fir" , "fs" , "fs=256" , "Sampling rate" );
  add_param( "--fir" , "bandpass", "bandpass=0.3,35", "Band-pass filter between 0.3 and 35 Hz" );
  add_param( "--fir" , "lowpass" , "lowpass=35", "Low-pass filter with cutoff of 35 Hz" );
  add_param( "--fir" , "highpass" , "highpass=0.3", "High-pass filter with cutoff of 0.3 Hz" );
  add_param( "--fir" , "bandstop", "bandstop=55,65", "Band-stop filter between 0.3 and 35 Hz" );
  add_param( "--fir" , "ripple", "0.01", "Ripple (proportion); can be two values for bandpass filters)" );
  add_param( "--fir" , "tw" , "1" , "Transition width (in Hz)" );
  add_param( "--fir" , "file" , "coef.txt" , "Read FIR coefficients from a file" );
  add_param( "--fir" , "order" , "10" , "Fix FIR order" );
  add_param( "--fir" , "rectangular" , "" , "Specify a rectangular window" );
  add_param( "--fir" , "bartlett" , "" , "Specify a Bartlett window" );
  add_param( "--fir" , "hann" , "" , "Specify a Hann window" );
  add_param( "--fir" , "blackman" , "" , "Specify a Blackman window" );

  add_table( "--fir" , "" , "FIR design parameters" );
  add_var( "--fir" , "" , "FIR" , "Label for FIR filter (constructed from input parameters)" );
  add_var( "--fir" , "" , "FS" , "Sampling rate (from fs input parameter)" );
  add_var( "--fir" , "" , "NTAPS" , "Filter order (number of taps)" );

  add_table( "--fir" , "F,FIR" , "Frequency response characteristics" );
  add_var( "--fir" , "F,FIR" , "F" , "Frequency (Hz)" );
  add_var( "--fir" , "F,FIR" , "FIR" , "FIR filter label" );
  add_var( "--fir" , "F,FIR" , "MAG" , "Magnitude" );
  add_var( "--fir" , "F,FIR" , "MAG_DB" , "Magnitude (dB)" );
  add_var( "--fir" , "F,FIR" , "PHASE" , "Phase" );

  add_table( "--fir" , "FIR,SEC" , "Impulse response" );
  add_var( "--fir" , "FIR,SEC" , "SEC" , "Time (seconds)" );
  add_var( "--fir" , "FIR,SEC" , "FIR" , "FIR filter label" );
  add_var( "--fir" , "FIR,SEC" , "IR" , "Impulse response" );

  //
  // FILTER
  //

  add_cmd( "filter" , "FILTER" , "Filter one or more signals" );
  add_url( "FILTER" , "fir-filters/#filter" );
  add_verb( "FILTER" ,
            "Apply one filtering mode in place to one or more signals.\n"
            "FILTER chooses among: a Kaiser-window FIR design "
            "(ripple + tw + one of bandpass/lowpass/highpass/bandstop), a fixed-order "
            "windowed FIR design (order + optional window + one band specification), an "
            "external FIR coefficient file (file), a narrow-band Gaussian filter (ngaus), "
            "or an IIR Butterworth/Chebyshev filter (butterworth or chebyshev + one band "
            "specification). The filtered signal replaces the original channel." );
  add_param( "FILTER" , "sig" , "C3,C4" , "Signals to filter [required]" );
  add_param( "FILTER" , "bandpass" , "0.3,35" , "Band-pass filter with lower and upper cutoff frequencies" );
  add_param( "FILTER" , "lowpass"  , "35"  , "Low-pass filter cutoff frequency" );
  add_param( "FILTER" , "highpass" , "0.3" , "High-pass filter cutoff frequency" );
  add_param( "FILTER" , "bandstop" , "55,65" , "Band-stop filter with lower and upper cutoff frequencies" );
  add_param( "FILTER" , "ripple" , "0.02" , "Kaiser-window ripple specification" );
  add_param( "FILTER" , "tw" , "1" , "Kaiser-window transition width (Hz)" );
  add_param( "FILTER" , "order" , "100" , "Fixed FIR order when not using Kaiser-window design" );
  add_param( "FILTER" , "file" , "coef.txt" , "Read FIR coefficients from this external file" );
  add_param( "FILTER" , "rectangular" , "" , "Use a rectangular FIR window" );
  add_param( "FILTER" , "bartlett" , "" , "Use a Bartlett FIR window" );
  add_param( "FILTER" , "hann" , "" , "Use a Hann FIR window" );
  add_param( "FILTER" , "blackman" , "" , "Use a Blackman FIR window" );
  add_param( "FILTER" , "fft" , "T" , "Use FFT-based FIR convolution (default T)" );
  add_param( "FILTER" , "fix-nyquist" , "0.5" , "Clamp upper transition frequencies below Nyquist by this amount (Hz)" );
  add_param( "FILTER" , "ngaus" , "13,2" , "Use a narrow-band Gaussian filter centered at freq with this FWHM" );
  add_param( "FILTER" , "butterworth" , "4" , "Use an IIR Butterworth filter of this order" );
  add_param( "FILTER" , "chebyshev" , "4,0.5" , "Use an IIR Chebyshev filter with order and epsilon" );
  add_param( "FILTER" , "silent" , "" , "Suppress filtering progress messages" );

  //
  // FILTER-DESIGN
  //

  add_cmd( "filter" , "FILTER-DESIGN" , "Design and summarize a FIR filter" );
  add_url( "FILTER-DESIGN" , "fir-filters/#filter-design" );
  add_verb( "FILTER-DESIGN" ,
            "Design one FIR filter and emit its impulse and frequency response summaries.\n"
            "FILTER-DESIGN only covers FIR modes: either a Kaiser-window design "
            "(fs + ripple + tw + one band specification), a fixed-order windowed design "
            "(fs + order + optional window + one band specification), or summarizing an "
            "external coefficient file (fs + file). It does not run the IIR or narrow-Gaussian "
            "modes that FILTER supports." );
  add_param( "FILTER-DESIGN" , "fs" , "200" , "Sampling rate used to design and evaluate the filter [required]" );
  add_param( "FILTER-DESIGN" , "bandpass" , "0.3,35" , "Band-pass filter with lower and upper cutoff frequencies" );
  add_param( "FILTER-DESIGN" , "lowpass"  , "35" , "Low-pass filter cutoff frequency" );
  add_param( "FILTER-DESIGN" , "highpass" , "0.3" , "High-pass filter cutoff frequency" );
  add_param( "FILTER-DESIGN" , "bandstop" , "55,65" , "Band-stop filter with lower and upper cutoff frequencies" );
  add_param( "FILTER-DESIGN" , "ripple" , "0.02" , "Kaiser-window ripple specification" );
  add_param( "FILTER-DESIGN" , "tw" , "1" , "Kaiser-window transition width (Hz)" );
  add_param( "FILTER-DESIGN" , "order" , "100" , "Fixed FIR order when not using Kaiser-window design" );
  add_param( "FILTER-DESIGN" , "file" , "coef.txt" , "Read FIR coefficients from this external file and summarize them" );
  add_param( "FILTER-DESIGN" , "rectangular" , "" , "Use a rectangular FIR window" );
  add_param( "FILTER-DESIGN" , "bartlett" , "" , "Use a Bartlett FIR window" );
  add_param( "FILTER-DESIGN" , "hann" , "" , "Use a Hann FIR window" );
  add_param( "FILTER-DESIGN" , "blackman" , "" , "Use a Blackman FIR window" );
  add_param( "FILTER-DESIGN" , "fix-nyquist" , "0.5" , "Clamp upper transition frequencies below Nyquist by this amount (Hz)" );

  add_table( "FILTER-DESIGN" , "" , "FIR design parameters" );
  add_var( "FILTER-DESIGN" , "" , "FIR" , "Filter label constructed from input parameters" );
  add_var( "FILTER-DESIGN" , "" , "FS" , "Sampling rate used to evaluate the filter" );
  add_var( "FILTER-DESIGN" , "" , "NTAPS" , "Number of FIR taps" );

  add_table( "FILTER-DESIGN" , "F,FIR" , "Frequency-response characteristics" );
  add_var( "FILTER-DESIGN" , "F,FIR" , "F" , "Frequency (Hz)" );
  add_var( "FILTER-DESIGN" , "F,FIR" , "FIR" , "Filter label" );
  add_var( "FILTER-DESIGN" , "F,FIR" , "MAG" , "Magnitude response" );
  add_var( "FILTER-DESIGN" , "F,FIR" , "MAG_DB" , "Magnitude response (dB)" );
  add_var( "FILTER-DESIGN" , "F,FIR" , "PHASE" , "Phase response" );

  add_table( "FILTER-DESIGN" , "FIR,SEC" , "Impulse response" );
  add_var( "FILTER-DESIGN" , "FIR,SEC" , "SEC" , "Time (seconds)" );
  add_var( "FILTER-DESIGN" , "FIR,SEC" , "FIR" , "Filter label" );
  add_var( "FILTER-DESIGN" , "FIR,SEC" , "IR" , "Impulse-response coefficient" );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // ARTIFACTS
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // EDGER
  //

  add_cmd( "artifact" , "EDGER" , "Flag likely leading/trailing artifact" );
  add_url( "EDGER" , "artifacts/#edger" );
  add_verb( "EDGER" ,
            "Detect likely artifactual leading and trailing portions of a recording.\n"
            "EDGER summarizes epoch-level Hjorth features, finds sustained outlier-heavy "
            "regions near the start or end of the night, and can either report the inferred "
            "lights-off/on boundaries or convert them into the standard epoch mask." );

  add_param( "EDGER" , "sig" , "C3,C4" , "Signals to use" );
  add_param( "EDGER" , "cache" , "c1" , "Cache lights-off/on for subsequent hypnograms" );
  add_param( "EDGER" , "only-start" , "" , "Only trim the start" );
  add_param( "EDGER" , "only-end" , "" , "Only trim the end" );
  add_param( "EDGER" , "th" , "3" , "Hjorth outlier threshold in SD units" );
  add_param( "EDGER" , "allow" , "20" , "Do not allow more than 20 epochs (10 mins) 'good' data at either end (default)" );
  add_param( "EDGER" , "req" , "10" , "By default, require equivalent of 10 epochs of bad data to be flagged" );
  add_param( "EDGER" , "w" , "9" , "By default, smoothing window (total window, in epochs) - otherwise, 4 epochs either side (w=9)" );
  add_param( "EDGER" , "all" , "" , "Use all epochs to norm metrics (default: sleep only)" );
  add_param( "EDGER" , "wake" , "" , "Use only wake epochs to norm metrics (default: sleep only)" );
  add_param( "EDGER" , "h2" , "" , "Include second Hjorth parameter (default is not to)" );
  add_param( "EDGER" , "epoch" , "" , "Emit epoch-level output (same as verbose)" );
  add_param( "EDGER" , "verbose" , "" , "Emit epoch-level output" );
  add_param( "EDGER" , "mask" , "" , "Set epoch mask based, for subsequent RE" );

  add_table( "EDGER" , "" , "Individual-level summaries" );
  add_var( "EDGER" , "" , "EOFF", "Epoch for lights off, if set (i.e. leading period trimming)" );
  add_var( "EDGER" , "" , "LOFF", "Clock-time for lights off, if set (i.e. leading period trimming)" );
  add_var( "EDGER" , "" , "EON" , "Epoch for lights on, if set (i.e. trailing period trimming)" );
  add_var( "EDGER" , "" , "LON" , "Clock-time for lights on, if set (i.e. leading period trimming)" );

  add_table( "EDGER" , "CH" , "Individual-level summaries" );
  add_var( "EDGER" , "CH" , "EOFF", "Epoch for lights off, if set based on this channel" );
  add_var( "EDGER" , "CH" , "EON" , "Epoch for lights on, if set, based on this channel" );

  add_table( "EDGER" , "CH,E" , "Individual-level summaries" );
  add_var( "EDGER" , "CH,E" , "H1", "First Hjorth statistic" );
  add_var( "EDGER" , "CH,E" , "H2", "Second Hjorth statistic (if h2 set)" );
  add_var( "EDGER" , "CH,E" , "H3", "Third Hjorth statistic" );
  add_var( "EDGER" , "CH,E" , "FLAG", "Indicator of whether epoch is bad (0=good/1=bad)" );
  add_var( "EDGER" , "CH,E" , "STAT", "Smoothed flag (0/1) indicator" );
  add_var( "EDGER" , "CH,E" , "XON", "Primary statistic to determine lights on (trailing artifact)" );
  add_var( "EDGER" , "CH,E" , "XOFF", "Primary statistic to determine lights off (leading artifact)" );
  add_var( "EDGER" , "CH,E" , "TRIM", "Indicator (0/1) for whether this epoch should be trimmed" );

  //
  // ARTIFACTS
  //

  add_cmd( "artifact" , "ARTIFACTS" , "Detect EEG artifacts following Buckelmueller et al." );
  add_url( "ARTIFACTS" , "artifacts/#artifacts" );
  add_verb( "ARTIFACTS" ,
            "Apply the epoch-level Buckelmuller et al. artifact detector to one or more EEG "
            "signals.\n"
            "This compares each epoch's delta and beta power to local moving averages, flags "
            "outlier epochs, and optionally converts those detections into the standard mask." );
  add_param( "ARTIFACTS" , "sig" , "C3,C4" , "Restrict analysis to these channels" );
  add_param( "ARTIFACTS" , "verbose" , "" , "Report epoch-level statistics" );
  add_param( "ARTIFACTS" , "epoch" , "" , "Alias for verbose epoch-level statistics" );
  add_param( "ARTIFACTS" , "no-mask" , "" , "Do not set mask for outlier epochs" );

  add_table( "ARTIFACTS" , "CH" , "Per-channel output" );
  add_var( "ARTIFACTS" , "CH" , "FLAGGED_EPOCHS" , "Number of epochs failing" );
  add_var( "ARTIFACTS" , "CH" , "ALTERED_EPOCHS" , "Number of epochs actually masked" );
  add_var( "ARTIFACTS" , "CH" , "TOTAL_EPOCHS" , "Number of epochs tested" );


  add_table( "ARTIFACTS" , "CH,E" , "Per-channel per-epoch output [verbose]" );
  add_var( "ARTIFACTS" , "CH,E" , "DELTA" , "Delta power" );
  add_var( "ARTIFACTS" , "CH,E" , "DELTA_AVG" , "Local average delta power" );
  add_var( "ARTIFACTS" , "CH,E" , "DELTA_FAC" , "Relative delta factor" );

  add_var( "ARTIFACTS" , "CH,E" , "BETA" , "Beta power" );
  add_var( "ARTIFACTS" , "CH,E" , "BETA_AVG" , "Local average beta power" );
  add_var( "ARTIFACTS" , "CH,E" , "BETA_FAC" , "Relative beta factor" );

  add_var( "ARTIFACTS" , "CH,E" , "DELTA_MASK" , "Masked based on delta power?" );
  add_var( "ARTIFACTS" , "CH,E" , "BETA_MASK" , "Masked based on beta power?" );
  add_var( "ARTIFACTS" , "CH,E" , "MASK" , "Is this epoch masked?" );

  //
  // SUPPRESS-ECG
  //

  add_cmd( "artifact" , "SUPPRESS-ECG" , "Detect/remove cardiac-contamination from the EEG" );
  add_url( "SUPPRESS-ECG" , "artifacts/#suppress-ecg" );
  add_verb( "SUPPRESS-ECG" ,
            "Estimate and suppress ECG contamination in EEG channels using an ECG reference.\n"
            "The command detects R peaks in the ECG channel, builds an average artifact "
            "template around those peaks, optionally masks epochs with unreliable heart-rate "
            "estimates, and subtracts the estimated cardiac artifact from each target signal." );
  add_param( "SUPPRESS-ECG" , "sig" , "C3,C4" , "Target EEG signals to correct" );
  add_param( "SUPPRESS-ECG" , "ecg" , "ECG" , "ECG reference channel used to detect cardiac peaks [required]" );
  add_param( "SUPPRESS-ECG" , "sr" , "125" , "Working sample rate for the ECG and target EEG channels" );
  add_param( "SUPPRESS-ECG" , "no-suppress" , "" , "Do not alter any EEG channels" );
  add_param( "SUPPRESS-ECG" , "mask-bad-epochs" , "" , "Mask epochs with invalid ECG-derived heart-rate estimates" );
  
  add_table( "SUPPRESS-ECG" , "" , "Individual-level summaries" );
  add_var( "SUPPRESS-ECG" , "" , "BPM" , "Mean heart rate (bpm)" );
  add_var( "SUPPRESS-ECG" , "" , "BPM_L95" , "Lower 95% confidence interval for mean HR" );
  add_var( "SUPPRESS-ECG" , "" , "BPM_U95" , "Upper 95% confidence interval for mean HR" );
  add_var( "SUPPRESS-ECG" , "" , "BPM_N_REMOVED" , "Number of epochs flagged as having invalid HR estimates" );
  add_var( "SUPPRESS-ECG" , "" , "BPM_PCT_REMOVED" , "Proportion of epochs flagged as having invalid HR estimates" );

  add_table( "SUPPRESS-ECG" , "E" , "Epoch-level metrics" );
  add_var( "SUPPRESS-ECG" , "E" , "BPM" , "HR for this epoch" );
  add_var( "SUPPRESS-ECG" , "E" , "BPM_MASK" , "Was this epoch invalid?" );
  
  add_table( "SUPPRESS-ECG" , "CH" , "Channel-level metrics" );
  add_var( "SUPPRESS-ECG" , "CH" , "ART_RMS" , "Root mean square of correction signature" );
  
  add_table( "SUPPRESS-ECG" , "CH,SP" , "Details of artifact signature" );
  add_var( "SUPPRESS-ECG" , "CH,SP" , "ART" , "Estimated artifact template value at each sample point in the correction window" );

  //
  // ALTER
  //

  add_cmd( "artifact" , "ALTER" , "Regression- or EMD-based artifact correction" );
  add_url( "ALTER" , "artifacts/#alter" ); 
  add_verb( "ALTER" ,
            "Correct artifact-contaminated signals using one or more reference channels.\n"
            "ALTER either fits a segment-wise regression model to remove variance explained by "
            "the corrector channels, or runs an EMD-based correction mode that removes selected "
            "components before reconstructing the target signal." );
  add_param( "ALTER" , "sig" , "C3,C4" , "Signals for analysis" );
  add_param( "ALTER" , "corr" , "EOG-R,EOG-L" , "Template signal(s)" );
  add_param( "ALTER" , "emd" , "1" , "Use EMD-based correction mode (non-zero integer)" );
  add_param( "ALTER" , "th" , "0.9" , "EMD correlation threshold for removing components" );
  add_param( "ALTER" , "emd-corr" , "" , "Run EMD of corrector channels" );

  add_param( "ALTER" , "segment-sec" , "4" , "Segment size" );
  add_table( "ALTER" , "NC" , "EMD correction summaries" );
  add_var( "ALTER" , "NC" , "REMOVED" , "Number of segments where this many components were removed" );

  //
  // LINE-DENOISE
  //

  add_cmd( "artifact" , "LINE-DENOISE" , "Remove line noise from signals" );
  add_url( "LINE-DENOISE" , "artifacts/#line-denoise" );
  add_verb( "LINE-DENOISE" ,
            "Remove line noise by interpolating over narrow spectral bands around one or more "
            "target frequencies.\n"
            "The denoiser can run on whole signals or epoch by epoch, replacing the amplitude "
            "in the noise band with the mean amplitude of neighboring frequency bins while "
            "retaining the original phase." );
  add_param( "LINE-DENOISE" , "sig" , "C3,C4" , "Signals to denoise [required]" );
  add_param( "LINE-DENOISE" , "f" , "50,100" , "Target line-noise frequencies to suppress [required]" );
  add_param( "LINE-DENOISE" , "w" , "1,1" , "Noise-band and neighboring-band widths in Hz" );
  add_param( "LINE-DENOISE" , "epoch" , "" , "Run the denoiser separately within each epoch" );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // HYPNOGRAMS
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // STAGE
  //

  add_cmd( "hypno" , "STAGE" , "Output sleep stages per epoch" );
  add_url( "STAGE" , "hypnograms/#stage" );
  add_verb( "STAGE" ,
            "Resolve sleep staging annotations and emit a simple per-epoch stage table.\n"
            "STAGE is the lightweight hypnogram view: it can derive stages from loaded "
            "annotations or from an external file, optionally force regeneration, and either "
            "print epoch rows or write a simple .eannot-style stage series." );
  add_param( "STAGE" , "file" , "stages.txt" , "Read stage labels from an external file instead of annotations" );
  add_param( "STAGE" , "N1" , "NREM1" , "Annotation label used for N1 sleep" );
  add_param( "STAGE" , "N2" , "NREM2" , "Annotation label used for N2 sleep" );
  add_param( "STAGE" , "N3" , "NREM3" , "Annotation label used for N3 sleep" );
  add_param( "STAGE" , "N4" , "NREM4" , "Annotation label used for N4 sleep" );
  add_param( "STAGE" , "R" , "REM" , "Annotation label used for REM sleep" );
  add_param( "STAGE" , "W" , "W" ,  "Annotation label used for wake" );
  add_param( "STAGE" , "L" , "LIGHTS_ON" ,  "Annotation label used for lights-on / lights markers" );
  add_param( "STAGE" , "?" , "-9" , "Annotation label used for unknown or other epochs" );
  add_param( "STAGE" , "force" , "" , "Rebuild sleep staging annotations even if they already exist" );
  add_param( "STAGE" , "eannot" , "stages.eannot" , "Write stages to a Luna .eannot-style file" );
  add_param( "STAGE" , "min" , "" , "Minimal stage output to stdout in .eannot style" );

  add_table( "STAGE" , "E" , "Stage annotations per-epoch" );
  add_var( "STAGE" , "E" , "CLOCK_TIME" , "Clock time (hh:mm:ss)" );
  add_var( "STAGE" , "E" , "MINS" , "Elapsed time from start of EDF (minutes)" );
  add_var( "STAGE" , "E" , "STAGE" , "Sleep stage (text value)" );
  add_var( "STAGE" , "E" , "STAGE_N" , "Numeric encoding of sleep stage" );

  //
  // HYPNO
  //

  add_cmd( "hypno" , "HYPNO" , "Sleep macro-architecture summaries" );
  add_url( "HYPNO" , "hypnograms/#hypno" );
  add_verb( "HYPNO" ,
            "Construct a full hypnogram summary from epoch-level sleep stages.\n"
            "HYPNO reports whole-night macro-architecture metrics such as total sleep, sleep "
            "latency, REM latency, wake after sleep onset, sleep efficiency, cycle counts, and "
            "timing landmarks such as lights off, sleep onset, midpoint, final wake, and "
            "lights on. It also provides stage-level summaries, optional epoch-level state "
            "annotations, NREM-cycle summaries, bout tables, and transition matrices when those "
            "views are requested." );

  add_param( "HYPNO" , "file" , "stages.txt" , "Optionally, read stages from file" );
  add_param( "HYPNO" , "N1" , "NREM1"  , "Annotation label used for N1 sleep" );
  add_param( "HYPNO" , "N2" , "NREM2" , "Annotation label used for N2 sleep" );
  add_param( "HYPNO" , "N3" , "NREM3" , "Annotation label used for N3 sleep" );
  add_param( "HYPNO" , "N4" , "NREM4" , "Annotation label used for N4 sleep" );
  add_param( "HYPNO" , "R" , "REM" , "Annotation label used for REM sleep" );
  add_param( "HYPNO" , "W" , "W" , "Annotation label used for wake" );
  add_param( "HYPNO" , "L" , "LIGHTS_ON" , "Annotation label used for lights markers" );
  add_param( "HYPNO" , "?" , "-9" , "Annotation label used for unknown/other" );
  add_param( "HYPNO" , "force" , "" , "Rebuild sleep staging annotations even if they already exist" );
  add_param( "HYPNO" , "epoch" , "" , "Emit epoch-level hypnogram outputs" );
  add_param( "HYPNO" , "verbose" , "T" , "Control cycle and transition outputs (default T)" );
  add_param( "HYPNO" , "annot" , "hyp_" , "Add hypnogram-derived annotations with this optional prefix" );
  add_param( "HYPNO" , "annot-cycles" , "hyp_" , "Legacy alias for annot" );


  add_table( "HYPNO" , "" , "Individual-level output" );
  add_var( "HYPNO" , "" , "TRT" , "Total recording time (minutes)" );
  add_var( "HYPNO" , "" , "TIB" , "Time in bed (minutes)" );
  add_var( "HYPNO" , "" , "TGT" , "Total gap time in discontinuous recordings (minutes)" );
  add_var( "HYPNO" , "" , "TST" , "Total sleep time (minutes)" );
  add_var( "HYPNO" , "" , "TST_PER" , "Persistent total sleep time (minutes)" );
  add_var( "HYPNO" , "" , "TWT" , "Total wake time (minutes)" );
  add_var( "HYPNO" , "" , "WASO" , "Wake after sleep onset (minutes)" );
  add_var( "HYPNO" , "" , "FWT" , "Final wake time before lights on (minutes)" );
  add_var( "HYPNO" , "" , "LOT" , "Lights-on stage duration (minutes)" );
  add_var( "HYPNO" , "" , "LOST" , "Minutes of sleep relabeled as lights-on" );
  add_var( "HYPNO" , "" , "SPT" , "Sleep period time excluding final wake (minutes)" );
  add_var( "HYPNO" , "" , "SPT_PER" , "Persistent sleep period time (minutes)" );
  add_var( "HYPNO" , "" , "PRE" , "Pre-sleep wake duration before first sleep epoch (minutes)" );
  add_var( "HYPNO" , "" , "POST" , "Post-sleep wake duration after final wake (minutes)" );
  add_var( "HYPNO" , "" , "OTHR" , "Unknown or other stage duration (minutes)" );
  add_var( "HYPNO" , "" , "CONF" , "Number of epochs with conflicting stage assignments" );
  add_var( "HYPNO" , "" , "SHORT" , "Indicator that the sleep period was constrained as too short" );
  add_var( "HYPNO" , "" , "FIXED_SLEEP" , "Epochs fixed while enforcing persistent sleep rules" );
  add_var( "HYPNO" , "" , "FIXED_WAKE" , "Epochs fixed due to excessive WASO" );
  add_var( "HYPNO" , "" , "FIXED_LIGHTS" , "Epochs fixed due to lights-on handling" );
  add_var( "HYPNO" , "" , "SINS" , "Indicator that the study starts in sleep" );
  add_var( "HYPNO" , "" , "EINS" , "Indicator that the study ends in sleep" );

  add_var( "HYPNO" , "" , "MINS_ASC_N2" , "Ascending N2 duration (minutes)" );
  add_var( "HYPNO" , "" , "MINS_DSC_N2" , "Descending N2 duration (minutes)" );
  add_var( "HYPNO" , "" , "MINS_FLT_N2" , "Flat N2 duration (minutes)" );
  add_var( "HYPNO" , "" , "PCT_ASC_N2" , "Proportion of N2 that is ascending" );
  add_var( "HYPNO" , "" , "PCT_DSC_N2" , "Proportion of N2 that is descending" );
  add_var( "HYPNO" , "" , "PCT_FLT_N2" , "Proportion of N2 that is flat" );

  add_var( "HYPNO" , "" , "T0_START" , "Recording start, hrs since prior midnight " );
  add_var( "HYPNO" , "" , "T1_LIGHTS_OFF" , "Lights off, hrs since prior midnight" );
  add_var( "HYPNO" , "" , "T2_SLEEP_ONSET" , "Sleep onset, hrs since prior midnight" );
  add_var( "HYPNO" , "" , "T3_SLEEP_MIDPOINT" , "Sleep midpoint, hrs since prior midnight" );
  add_var( "HYPNO" , "" , "T4_FINAL_WAKE" , "Final wake, hrs since prior midnight" );
  add_var( "HYPNO" , "" , "T5_LIGHTS_ON" , "Lights on, hrs since prior midnight" );
  add_var( "HYPNO" , "" , "T6_STOP" , "Study stop, hrs since prior midnight" );

  add_var( "HYPNO" , "" , "E0_START" , "Recording start, elapsed time" );
  add_var( "HYPNO" , "" , "E1_LIGHTS_OFF" , "Lights off, elapsed time" );
  add_var( "HYPNO" , "" , "E2_SLEEP_ONSET" , "Sleep onset, elapsed time" );
  add_var( "HYPNO" , "" , "E3_SLEEP_MIDPOINT" , "Sleep midpoint, elapsed time" );
  add_var( "HYPNO" , "" , "E4_FINAL_WAKE" , "Final wake, elapsed time" );
  add_var( "HYPNO" , "" , "E5_LIGHTS_ON" , "Lights on, elapsed time" );
  add_var( "HYPNO" , "" , "E6_STOP" , "Study stop, elapsed time" );


  add_var( "HYPNO" , "" , "HMS0_START" , "Recording start, clock time" );
  add_var( "HYPNO" , "" , "HMS1_LIGHTS_OFF" , "Lights off, clock time" );
  add_var( "HYPNO" , "" , "HMS2_SLEEP_ONSET" , "Sleep onset, clock time" );
  add_var( "HYPNO" , "" , "HMS3_SLEEP_MIDPOINT" , "Sleep midpoint, clock time" );
  add_var( "HYPNO" , "" , "HMS4_FINAL_WAKE" , "Final wake, clock time" );
  add_var( "HYPNO" , "" , "HMS5_LIGHTS_ON" , "Lights on, clock time" );
  add_var( "HYPNO" , "" , "HMS6_STOP" , "Study stop, clock time" );
  
  
  add_var( "HYPNO" , "" , "SE" , "Sleep efficiency" );
  add_var( "HYPNO" , "" , "SME" , "Sleep efficiency (alternate defn.)" );
  add_var( "HYPNO" , "" , "SOL" , "Sleep latency (minutes from lights off)" );
  add_var( "HYPNO" , "" , "SOL_PER" , "Persistent sleep latency (mins from lights off)" );
  add_var( "HYPNO" , "" , "REM_LAT" , "REM latency (minutes from onset of sleep)" );
  add_var( "HYPNO" , "" , "REM_LAT2" , "REM latency (excluding wake)" );
  add_var( "HYPNO" , "" , "NREMC" , "Number of sleep cycles" );
  add_var( "HYPNO" , "" , "NREMC_MINS" , "Mean duration of each sleep cycle" );

  add_var( "HYPNO" , "" , "SFI" , "Sleep Fragmentation Index" );
  add_var( "HYPNO" , "" , "TI_S" , "Sleep Transition Index" );
  add_var( "HYPNO" , "" , "TI_S3" , "Sleep Transition Index, 3-stage classification" );
  add_var( "HYPNO" , "" , "TI_RNR" , "Sleep Transition Index: REM-NREM only" );
  add_var( "HYPNO" , "" , "RUNS" , "Lempel-Ziv style run complexity metric on 5-stage scoring" );
  add_var( "HYPNO" , "" , "RUNS3" , "Lempel-Ziv style run complexity metric on 3-stage scoring" );
  add_var( "HYPNO" , "" , "LZW" , "Normalized LZW complexity on 5-stage scoring" );
  add_var( "HYPNO" , "" , "LZW3" , "Normalized LZW complexity on 3-stage scoring" );
  
  add_table( "HYPNO" , "SS" , "Stage-stratified output" );
  add_var( "HYPNO" , "SS" , "MINS" , "Stage duration (mins)" );
  add_var( "HYPNO" , "SS" , "PCT" , "Stage duration (% of TST)" );
  add_var( "HYPNO" , "SS" , "DENS" , "Stage duration as a proportion of sleep period time" );
  add_var( "HYPNO" , "SS" , "BOUT_N" , "Number of bouts" );
  add_var( "HYPNO" , "SS" , "BOUT_MX" , "Maximum bout duration (minutes)" );
  add_var( "HYPNO" , "SS" , "BOUT_MN" , "Mean bout duration" );
  add_var( "HYPNO" , "SS" , "BOUT_MD" , "Median bout duration" );
  add_var( "HYPNO" , "SS" , "BOUT_05" , "Stage duration in bouts at least 5 minutes long" );
  add_var( "HYPNO" , "SS" , "BOUT_10" , "Stage duration (only bouts 10+ mins)" );
  add_var( "HYPNO" , "SS" , "TA" , "Median epoch timing (vs all from sleep onset-offset)");
  add_var( "HYPNO" , "SS" , "TS" , "Median epoch timing (vs elapsed sleep from sleep onset-offset)");
  

  add_table( "HYPNO" , "C" , "NREM cycle-level output" );
  add_var( "HYPNO" , "C" , "NREMC_START" , "First epoch number of this NREM cycle" );
  add_var( "HYPNO" , "C" , "NREMC_MINS" , "Total duration of this cycle (mins)" );
  add_var( "HYPNO" , "C" , "NREMC_NREM_MINS" , "Duration of NREM in this cycle (mins)" );
  add_var( "HYPNO" , "C" , "NREMC_REM_MINS" , "Duration of REM in this cycle (mins)" );
  add_var( "HYPNO" , "C" , "NREMC_OTHER_MINS" , "Minutes of wake and unscored epochs" );
  add_var( "HYPNO" , "C" , "NREMC_N" , "Total number of epochs in this cycle" );


  add_table( "HYPNO", "N" , "Bouts" ); 
  add_var( "HYPNO" , "N" , "STAGE" , "Stage label for this bout" );
  add_var( "HYPNO" , "N" , "FIRST_EPOCH" , "First epoch" );
  add_var( "HYPNO" , "N" , "LAST_EPOCH" , "Last epoch" );
  add_var( "HYPNO" , "N" , "START" , "Start (clocktime)" );
  add_var( "HYPNO" , "N" , "STOP" , "Stop (clocktime) [ end of last epoch ]" );
  add_var( "HYPNO" , "N" , "MINS" , "Bout duration (minutes)" );
  

  add_table( "HYPNO" , "E" , "Epoch-level output" );
  add_var( "HYPNO" , "E" , "CLOCK_HOURS" , "Start time of epoch (hours since midnight)" );
  add_var( "HYPNO" , "E" , "CLOCK_TIME" , "Start time of epoch (hh:mm:ss)" );
  add_var( "HYPNO" , "E" , "MINS" , "Elapsed minutes" );
  add_var( "HYPNO" , "E" , "AFTER_GAP" , "Indicator that this epoch immediately follows a discontinuity gap" );
  add_var( "HYPNO" , "E" , "START_SEC" , "Start time (seconds since start of EDF)" );
  add_var( "HYPNO" , "E" , "STAGE" , "Text description of sleep stage" );
  add_var( "HYPNO" , "E" , "OSTAGE" , "Original stage label (pre any modifications)" );
  add_var( "HYPNO" , "E" , "STAGE_N" , "Numeric encoding of sleep stage" );
  add_var( "HYPNO" , "E" , "PERSISTENT_SLEEP" , "Flag to indicate persistent sleep" );
  add_var( "HYPNO" , "E" , "SPT" , "Indicator that this epoch lies within the sleep period time window" );
  add_var( "HYPNO" , "E" , "PRE" , "Indicator that this epoch is before sleep onset" );
  add_var( "HYPNO" , "E" , "POST" , "Indicator that this epoch is after final wake" );
  add_var( "HYPNO" , "E" , "WASO" , "Flag to indicate wake after sleep onset" );
  add_var( "HYPNO" , "E" , "E_N1" , "Cumulative elapsed N1 sleep (minutes)" );
  add_var( "HYPNO" , "E" , "E_N2" , "Cumulative elapsed N2 sleep (minutes)" );
  add_var( "HYPNO" , "E" , "E_N3" , "Cumulative elapsed N3 sleep (minutes)" );
  add_var( "HYPNO" , "E" , "E_REM" , "Cumulative elapsed REM (minutes)" );
  add_var( "HYPNO" , "E" , "E_SLEEP" , "Cumulative elapsed sleep (minutes)" );
  add_var( "HYPNO" , "E" , "E_WAKE" , "Cumulative elapsed wake (minutes)" );
  add_var( "HYPNO" , "E" , "E_WASO" , "Cumulative elapsed WASO (minutes)" );
  add_var( "HYPNO" , "E" , "PCT_E_N1" , "Cumulative elapsed N1 as proportion of total N1 sleep" );
  add_var( "HYPNO" , "E" , "PCT_E_N2" , "Cumulative elapsed N2 as proportion of total N2 sleep" );
  add_var( "HYPNO" , "E" , "PCT_E_N3" , "Cumulative elapsed N3 as proportion of total N3 sleep" );
  add_var( "HYPNO" , "E" , "PCT_E_REM" , "Cumulative elapsed REM as proportion of total REM sleep" );
  add_var( "HYPNO" , "E" , "PCT_E_SLEEP" , "Cumulative elapsed sleep as proportion of total sleep" );
  add_var( "HYPNO" , "E" , "FLANKING" , "Minimum number of similarly staged epochs forward or backward" );
  add_var( "HYPNO" , "E" , "FLANKING_ALL" , "Total number of similarly staged epochs in this contiguous run" );

  add_var( "HYPNO" , "E" , "N2_WGT" , "Score to indicate ascending versus descending N2 sleep" );
  add_var( "HYPNO" , "E" , "NEAREST_WAKE" , "Number of epochs (forward or backwards) since nearest wake epoch" );


  add_var( "HYPNO" , "E" , "CYCLE" , "Cycle number, if this epoch is in a sleep cycle" );
  add_var( "HYPNO" , "E" , "CYCLE_POS_ABS" , "Absolute position of this epoch in the current NREM cycle (mins)" );
  add_var( "HYPNO" , "E" , "CYCLE_POS_REL" , "Relative position of this epoch in the current NREM cycle (0-1)" );
  add_var( "HYPNO" , "E" , "PERIOD" , "Cycle period: NREMP or REMP, or missing if not in a cycle" );

  add_var( "HYPNO" , "E" , "TR_NR2R" , "Number of epochs from this NREM epoch to a REM transition" );
  add_var( "HYPNO" , "E" , "TOT_NR2R" , "Total number of contiguous NREM epochs followed by REM" );
  add_var( "HYPNO" , "E" , "TR_NR2W" , "Number of epochs from this NREM epoch to a wake transition" );
  add_var( "HYPNO" , "E" , "TOT_NR2W" , "Total number of contiguous NREM epochs followed by wake" );
  add_var( "HYPNO" , "E" , "TR_R2W" , "Number of epochs from this REM epoch to a wake transition" );
  add_var( "HYPNO" , "E" , "TOT_R2W" , "Total number of contiguous REM epochs followed by wake" );
  add_var( "HYPNO" , "E" , "TR_R2NR" , "Number of epochs from this REM epoch to a NREM transition" );
  add_var( "HYPNO" , "E" , "TOT_R2NR" , "Total number of contiguous REM epochs followed by NREM" );
  add_var( "HYPNO" , "E" , "TR_W2R" , "Number of epochs from this wake epoch to a REM transition" );
  add_var( "HYPNO" , "E" , "TOT_W2R" , "Total number of contiguous wake epochs followed by REM" );
  add_var( "HYPNO" , "E" , "TR_W2NR" , "Number of epochs from this wake epoch to a NREM transition" );
  add_var( "HYPNO" , "E" , "TOT_W2NR" , "Total number of contiguous wake epochs followed by NREM" );

  
  add_table( "HYPNO" , "PRE,POST" , "Stage transitions" );
  add_var( "HYPNO", "PRE,POST" , "N" , "Number of transitions" );
  add_var( "HYPNO", "PRE,POST" , "P" , "Overall proportion of transitions in this cell" );
  add_var( "HYPNO" , "PRE,POST" , "P_POST_COND_PRE" , "P( S+1 | S )" );
  add_var( "HYPNO" , "PRE,POST" , "P_PRE_COND_POST" , "P( S | S+1 )" );

  //
  // DYNAM
  //

  add_cmd( "hypno" , "DYNAM" , "Summarize epoch-level outputs by NREM cycles" );
  add_url( "DYNAM" , "hypnograms/#dynam" );
  add_verb( "DYNAM" ,
            "Run the generic quantile-dynamics analysis on epoch-level variables, usually to "
            "summarize how a measure changes across sleep cycles or normalized sleep time.\n"
            "DYNAM reads one or more epoch-level tabular outputs, groups rows by epoch and any "
            "requested factors, and then derives quantile-profile summaries across the night." );
  add_param( "DYNAM" , "inputs" , "psd1.db" , "One or more input epoch-level tables" );
  add_param( "DYNAM" , "vars" , "PSD" , "Variables to include from the input files" );
  add_param( "DYNAM" , "facs" , "CH,B" , "Factor columns that define separate dynamics profiles" );
  add_param( "DYNAM" , "no-id" , "T" , "Ignore the ID column when matching rows to the current EDF" );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // STAGING
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // SOAP
  //

  add_cmd( "stage" , "SOAP" , "Single-observation accuracies and probabilities" );
  add_url( "SOAP" , "suds/#soap" );
  add_verb( "SOAP" ,
            "Fit a single-channel SOAP model against existing staging and report "
            "agreement, posterior probabilities, and stage summaries.\n\n"
            "SOAP builds an epoch-level feature set from the signal, reduces that "
            "feature space to a smaller set of components, and then uses an LDA/QDA "
            "classifier to predict stage probabilities from the observed staging. "
            "It is primarily a within-record evaluation tool: how well can one "
            "signal recover the manual scoring already present in this EDF?" );

  add_param( "SOAP" , "sig" , "C3,C4" , "Run SOAP separately for these signal(s)" );
  add_param( "SOAP" , "model" , "_1" , "SOAP model to load" );
  add_param( "SOAP" , "read-weights" , "weights.txt" , "Read pre-computed model weights" );
  add_param( "SOAP" , "write-weights" , "weights.txt" , "Write model weights for reuse" );
  add_param( "SOAP" , "force-reload" , "" , "Force model reload before fitting" );
  add_param( "SOAP" , "epoch" , "" , "Emit epoch-level posterior output" );
  add_param( "SOAP" , "3-class" , "" , "Pool N1/N2/N3 into a single NREM class" );
  add_param( "SOAP" , "trim" , "10" , "Trim leading/trailing wake before fitting" );
  add_param( "SOAP" , "segment-sec" , "4" , "Welch segment length (seconds)" );
  add_param( "SOAP" , "segment-overlap" , "2" , "Welch segment overlap (seconds)" );
  add_param( "SOAP" , "pc" , "0.01" , "ANOVA p-value threshold for retaining components" );
  add_param( "SOAP" , "all-c" , "" , "Retain all components" );
  add_param( "SOAP" , "within" , "1.5" , "Within/total variance threshold for components" );
  add_param( "SOAP" , "th" , "5,3" , "Trainer epoch outlier thresholds" );
  add_param( "SOAP" , "robust" , "0.1" , "Use robust scaling, optionally with winsorization" );
  add_param( "SOAP" , "self-prob" , "0.9" , "Require this posterior for self-classification" );
  add_param( "SOAP" , "self-kappa" , "0.5" , "Require this self-classification kappa" );
  add_param( "SOAP" , "th-hjorth" , "5" , "Hjorth outlier threshold for target epochs" );
  add_param( "SOAP" , "norm-X" , "T" , "Standardize predictor features" );
  add_param( "SOAP" , "norm-U" , "F" , "Standardize singular vectors" );
  add_param( "SOAP" , "best-guess" , "F" , "Use hard calls rather than weighted posteriors" );
  add_param( "SOAP" , "ignore-prior" , "" , "Ignore existing target-stage priors" );
  add_param( "SOAP" , "flat-priors" , "" , "Use flat priors" );
  add_param( "SOAP" , "fixed-priors" , "0.2,0.2,0.2,0.2,0.2" , "Specify fixed stage priors" );
  add_param( "SOAP" , "save" , "" , "Cache the fitted target for subsequent RESOAP use" );
  add_param( "SOAP" , "feature-matrix" , "" , "Dump the predictor matrix to stdout" );
  add_param( "SOAP" , "dump-features" , "soap.features" , "Write the predictor matrix to a file" );
  add_param( "SOAP" , "dump-stage-assocs" , "soap.assocs" , "Write feature/stage associations" );
  add_param( "SOAP" , "dump-svd" , "soap.svd" , "Write SVD components to a file" );
  add_param( "SOAP" , "trans" , "30" , "Summarize transition-centered posteriors in new epochs" );
  add_param( "SOAP" , "req-left" , "2" , "Required stable epochs before a transition" );
  add_param( "SOAP" , "req-right" , "2" , "Required stable epochs after a transition" );
  add_param( "SOAP" , "left" , "6" , "Epochs to show before each transition" );
  add_param( "SOAP" , "right" , "6" , "Epochs to show after each transition" );
  
  add_table( "SOAP" , "" , "Overall accuracies" );
  add_var( "SOAP" , "" , "NSS" , "Number of distinct observed stages" );
  add_var( "SOAP" , "" , "ACC" , "Accuracy" );
  add_var( "SOAP" , "" , "ACC3" , "Accuracy, 3-class" );
  add_var( "SOAP" , "" , "F1" , "F1 metric" );
  add_var( "SOAP" , "" , "F13" , "F1 metric, 3-class" );
  add_var( "SOAP" , "" , "F1_WGT" , "F1 metric, weighted" );
  add_var( "SOAP" , "" , "K" , "Kappa" );
  add_var( "SOAP" , "" , "K3" , "Kappa, 3-class" );
  add_var( "SOAP" , "" , "MCC" , "Matthews correlation coef" );
  add_var( "SOAP" , "" , "MCC3" , "Matthews correlation coef, 3-class" );
  add_var( "SOAP" , "" , "PREC" , "Precision" );
  add_var( "SOAP" , "" , "PREC3" , "Precision, 3-class" );  
  add_var( "SOAP" , "" , "RECALL" , "Recall" );
  add_var( "SOAP" , "" , "RECALL3" , "Recall, 3-class" );
  add_var( "SOAP" , "" , "RECALL_WGT" , "Recall, weighted" );
  
  add_table( "SOAP" , "E" , "Epoch-level output" );
  add_var( "SOAP" , "E" , "DISC" , "Discordant epoch" );
  add_var( "SOAP" , "E" , "DISC3" , "Discordant epoch, 3-class" );
  add_var( "SOAP" , "E" , "INC" , "Epoch included" );
  add_var( "SOAP" , "E" , "PP_N1" , "N1 posterior probability" );
  add_var( "SOAP" , "E" , "PP_N2" , "N2 posterior probability" );
  add_var( "SOAP" , "E" , "PP_N3" , "N3 posterior probability" );
  add_var( "SOAP" , "E" , "PP_NR" , "NR posterior probability" );
  add_var( "SOAP" , "E" , "PP_R" , "REM posterior probability" );
  add_var( "SOAP" , "E" , "PP_W" , "Wake posterior probability" );
  add_var( "SOAP" , "E" , "PRED" , "Predicted stage" );
  add_var( "SOAP" , "E" , "PRIOR" , "Original stage" );
  
  add_table( "SOAP" , "SS" , "Stage-level output" );
  add_var( "SOAP" , "SS" , "DUR_OBS" , "Observed stage duration (for included epochs)" );
  add_var( "SOAP" , "SS" , "DUR_PRD" , "Predicted stage duration (for included epochs)" );
  add_var( "SOAP" , "SS" , "F1" , "Stage-specific F1" );
  add_var( "SOAP" , "SS" , "PREC" , "Stage-specific precision" );
  add_var( "SOAP" , "SS" , "RECALL" , "Stage-specific recall" );

  add_table( "SOAP" , "VAR" , "PSC info" );
  add_var( "SOAP" , "VAR" , "INC" , "Component included" );
  add_var( "SOAP" , "VAR" , "PV" , "1-way ANOVA p-value for association w/ observed stage" );

  add_table( "SOAP" , "ETYPE" , "Epoch-type accuracy" );
  add_var( "SOAP" , "ETYPE" , "ACC" , "Accuracy" );
  add_var( "SOAP" , "ETYPE" , "N" , "Epoch count" );

  add_table( "SOAP" , "ETYPE,SS" , "Stage-specific epoch-type accuracy" );
  add_var( "SOAP" , "ETYPE,SS" , "ACC" , "Accuracy" );
  add_var( "SOAP" , "ETYPE,SS" , "N" , "Epoch count" );

  add_table( "SOAP" , "NSS,PRED,OBS" , "Confusion matrix" );
  add_var( "SOAP" , "NSS,PRED,OBS" , "N" , "Number" );
  add_var( "SOAP" , "NSS,PRED,OBS" , "P" , "Proportion" );

  //
  // POPS
  //

  add_cmd( "stage" , "POPS" , "Population-based staging" );
  add_url( "POPS" , "staging/#pops" );
  add_verb( "POPS" ,
            "Run the lower-level POPS staging engine on the current EDF using a "
            "feature specification and model library.\n\n"
            "POPS is the population model. It derives a larger multi-feature "
            "representation from the EDF, applies a pre-specified feature set and "
            "trained LightGBM model, and can optionally layer in elapsed-sleep "
            "priors and SOAP-style post-processing. In addition to stage calls, it "
            "can emit model diagnostics such as confusion tables, feature "
            "definitions, SHAP values, and channel-equivalent summaries." );

  add_param( "POPS" , "train" , "" , "Build POPS training datasets" );
  add_param( "POPS" , "features" , "m1.ftr" , "Feature specification file" );
  add_param( "POPS" , "data" , "pops/lib/^" , "Filename for binary training files" );
  add_param( "POPS" , "model" , "m1.model" , "LGBM model file to write to or read from" );
  add_param( "POPS" , "config" , "m1.config" , "LGBM configuration file" );
  add_param( "POPS" , "path" , "." , "Base path for POPS resources" );
  add_param( "POPS" , "lib" , "s2" , "POPS library root" );
  add_param( "POPS" , "force-reload" , "" , "Force re-reading the POPS specs" );
  add_param( "POPS" , "ignore-obs-staging" , "" , "Ignore existing observed staging" );
  add_param( "POPS" , "apply-ranges" , "T" , "Apply library-specific feature ranges" );
  add_param( "POPS" , "apply-priors" , "F" , "Apply elapsed-sleep priors from the library" );
  add_param( "POPS" , "priors-c" , "0.001" , "Intercept term for elapsed-sleep priors" );
  add_param( "POPS" , "priors-rolling" , "T" , "Use rolling elapsed-sleep priors" );
  add_param( "POPS" , "priors-weighted" , "T" , "Use weighted elapsed-sleep prior counts" );
  add_param( "POPS" , "priors-es-min" , "20" , "Elapsed-sleep prior bin width (minutes)" );
  add_param( "POPS" , "priors-es-max" , "380" , "Maximum elapsed-sleep prior time (minutes)" );
  add_param( "POPS" , "priors-nr-min" , "10" , "NREM prior bin width (minutes)" );
  add_param( "POPS" , "priors-nr-max" , "60" , "Maximum NREM prior time (minutes)" );
  add_param( "POPS" , "priors-nr-allow" , "5" , "Allowed non-NREM duration within NREM priors" );
  add_param( "POPS" , "inc-vars" , "SIG_STATS,SPINDLE" , "Restrict POPS predictors to these variables" );
  add_param( "POPS" , "exc-vars" , "ARTIFACT" , "Exclude these predictors" );
  add_param( "POPS" , "soap" , "0.5" , "Enable SOAP-style probability updates" );
  add_param( "POPS" , "soap-nc" , "10" , "Number of components in SOAP updates" );
  add_param( "POPS" , "soap-lwr" , "1" , "Lower SOAP likelihood bound" );
  add_param( "POPS" , "soap-upr" , "100" , "Upper SOAP likelihood bound" );
  add_param( "POPS" , "soap-steps" , "100" , "Number of SOAP likelihood steps" );
  add_param( "POPS" , "soap-grid" , "0.8" , "Run a SOAP confidence grid search" );
  add_param( "POPS" , "verbose" , "" , "Verbose model and feature diagnostics" );
  add_param( "POPS" , "stage-assoc" , "F" , "Run stage-association summaries" );
  add_param( "POPS" , "epoch-SHAP" , "" , "Emit epoch-level SHAP values" );
  add_param( "POPS" , "SHAP-epoch" , "" , "Alias for epoch-SHAP" );
  add_param( "POPS" , "3-class" , "" , "Pool N1/N2/N3 into a single NREM class" );
  add_param( "POPS" , "trim" , "10" , "Trim leading/trailing wake" );
  add_param( "POPS" , "fft-median" , "T" , "Use median Welch aggregation" );
  add_param( "POPS" , "lwr" , "0.5" , "Lower frequency bound" );
  add_param( "POPS" , "upr" , "45" , "Upper frequency bound" );
  add_param( "POPS" , "segment-sec" , "4" , "Welch segment length (seconds)" );
  add_param( "POPS" , "segment-overlap" , "2" , "Welch segment overlap (seconds)" );
  add_param( "POPS" , "iid-weights" , "w1,w2" , "Individual-level trainer weights" );
  add_param( "POPS" , "dump-weights" , "weights.txt" , "Write fitted model weights" );
  add_param( "POPS" , "fix" , "W,R,N1,N2,N3" , "Sample a fixed number of epochs per class" );
  add_param( "POPS" , "alias" , "CEN|C3" , "Specify channel aliases" );
  add_param( "POPS" , "equiv" , "CEN|C3,C3_N" , "Specify channel-equivalent mappings" );

  add_table( "POPS" , "" , "POPS metrics" );
  add_var( "POPS", "" , "ACC" , "Accuracy" );
  add_var( "POPS", "" , "ACC3" , "Accuracy for 3-class model" );
  add_var( "POPS", "" , "K" , "Kappa statistic" );
  add_var( "POPS", "" , "K3" , "Kappa for 3-class model" );
  add_var( "POPS", "" , "F1" , "F1 statistic" );
  add_var( "POPS", "" , "F13" , "F1 for 3-class model" );
  add_var( "POPS", "" , "F1_WGT" , "F1 weighted" );
  add_var( "POPS", "" , "CONF" , "Mean confidence (max. posterior)" );
  add_var( "POPS", "" , "MCC" , "Matthews correlation coefficient" );
  add_var( "POPS", "" , "MCC3" , "Matthews correlation coefficient, 3-class" );
  add_var( "POPS", "" , "PREC" , "Precision" );
  add_var( "POPS", "" , "PREC_WGT" , "Precision, weighted" );
  add_var( "POPS", "" , "PREC3" , "Precision, 3-class" );
  add_var( "POPS", "" , "RECALL" , "Recall" );
  add_var( "POPS", "" , "RECALL3" , "Recall, 3-class" );
  add_var( "POPS", "" , "RECALL_WGT" , "Recall,weighted" );
  add_var( "POPS", "" , "SLP_LAT_OBS" , "Observed sleep latency" );
  add_var( "POPS", "" , "SLP_LAT_PRD" , "Predicted sleep latency" );
  add_var( "POPS", "" , "REM_LAT_OBS" , "Observed REM latency" );
  add_var( "POPS", "" , "REM_LAT_PRD" , "Predicted REM latency" );
  
  
 
  add_table( "POPS" , "E" , "POPS predictions" );
  add_var( "POPS" , "E" , "FLAG" , "-1/0/1/2 excluded/match/disc5/disc3" );
  add_var( "POPS" , "E" , "CONF" , "Confidence score" );
  add_var( "POPS" , "E" , "PP_N1" , "Posterior probability of N1" );
  add_var( "POPS" , "E" , "PP_N2" , "Posterior probability of N2" );
  add_var( "POPS" , "E" , "PP_N3" , "Posterior probability of N3" );
  add_var( "POPS" , "E" , "PP_R" , "Posterior probability of REM" );
  add_var( "POPS" , "E" , "PP_W" , "Posterior probability of wake" );
  add_var( "POPS" , "E" , "PRED", "Predicted stage" );
  add_var( "POPS" , "E" , "PRIOR" , "Observed stage (if known)" );
  add_var( "POPS" , "E" , "START" , "Start time (hh:mm:ss)" );
  add_var( "POPS" , "E" , "STOP" , "Stop time (hh:mm:ss)");
  
  add_table( "POPS" , "SS" , "Sleep-stage summaries" );
  add_var( "POPS" , "SS" , "OBS" , "Observed stage duration (for included epochs)");
  add_var( "POPS" , "SS" , "ORIG" , "Observed stage duration (all epochs)");
  add_var( "POPS" , "SS" , "PRF" , "Predicted stage duration, weighted" );
  add_var( "POPS" , "SS" , "PR1" , "Predicted stage duration, based on most likely" );
  add_var( "POPS" , "SS" , "F1" , "F1 statistic" );
  add_var( "POPS" , "SS" , "RECALL" , "Recall" );
  add_var( "POPS" , "SS" , "PREC" , "Precision" );

  add_table( "POPS" , "FTR" , "Feature definitions" );
  add_var( "POPS" , "FTR" , "BLOCK" , "Block label" );
  add_var( "POPS" , "FTR" , "FINAL" , "Column, if included" );
  add_var( "POPS" , "FTR" , "INC" , "Included?" );
  add_var( "POPS" , "FTR" , "LABEL" , "Feature label" );
  add_var( "POPS" , "FTR" , "LABEL_ORIG" , "Feature label" );
  add_var( "POPS" , "FTR" , "LEVEL" , "Level (1/2)" );
  add_var( "POPS" , "FTR" , "ROOT" , "Root label" );

  add_table( "POPS" , "ETYPE" , "Error type" );
  add_var( "POPS" , "ETYPE" , "ACC" , "Accuracy" );
  add_var( "POPS" , "ETYPE" , "N" , "Count" );

  add_table( "POPS" , "SS,ETYPE" , "Stage-specific error type" );
  add_var( "POPS" , "SS,ETYPE" , "ACC" , "Accuracy" );
  add_var( "POPS" , "SS,ETYPE" , "N" , "Count" );

  add_table( "POPS" , "SS,FTR" , "SHAP values" );
  add_var( "POPS" , "SS,FTR" , "SHAP" , "SHAP values" );
 
  add_table( "POPS" , "E,SS,FTR" , "Epoch-level SHAP values" );
  add_var( "POPS" , "E,SS,FTR" , "SHAP" , "SHAP values" );
  set_compressed( "POPS" , tfac_t( "E,SS,FTR" ) );
  
  add_table( "POPS" , "PRED,OBS" , "Confusion matrix" );
  add_var( "POPS" , "PRED,OBS" , "N" , "Count" );
  add_var( "POPS" , "PRED,OBS" , "P" , "Proportion" );


  add_table( "POPS" , "CHEQ" , "Channel-equivalent stats" );
  add_var( "POPS", "CHEQ" , "ACC" , "Accuracy" );
  add_var( "POPS", "CHEQ" , "ACC3" , "Accuracy for 3-class model" );
  add_var( "POPS", "CHEQ" , "K" , "Kappa statistic" );
  add_var( "POPS", "CHEQ" , "K3" , "Kappa for 3-class model" );
  add_var( "POPS", "CHEQ" , "F1" , "F1 statistic" );
  add_var( "POPS", "CHEQ" , "F13" , "F1 for 3-class model" );
  add_var( "POPS", "CHEQ" , "F1_WGT" , "F1 weighted" );
  add_var( "POPS", "CHEQ" , "CONF" , "Mean confidence (max. posterior)" );
  add_var( "POPS", "CHEQ" , "MCC" , "Matthews correlation coefficient" );
  add_var( "POPS", "CHEQ" , "MCC3" , "Matthews correlation coefficient, 3-class" );
  add_var( "POPS", "CHEQ" , "PREC" , "Precision" );
  add_var( "POPS", "CHEQ" , "PREC_WGT" , "Precision, weighted" );
  add_var( "POPS", "CHEQ" , "PREC3" , "Precision, 3-class" );
  add_var( "POPS", "CHEQ" , "RECALL" , "Recall" );
  add_var( "POPS", "CHEQ" , "RECALL3" , "Recall, 3-class" );
  add_var( "POPS", "CHEQ" , "RECALL_WGT" , "Recall,weighted" );
  add_var( "POPS", "CHEQ" , "REM_LAT_OBS" , "Observed REM latency" );
  add_var( "POPS", "CHEQ" , "REM_LAT_PRD" , "Predicted REM latency" );
  add_var( "POPS", "CHEQ" , "SLP_LAT_OBS" , "Observed sleep latency" );
  add_var( "POPS", "CHEQ" , "SLP_LAT_PRD" , "Predicted sleep latency" );
  
  add_table( "POPS" , "E,CHEQ" , "POPS predictions" );
  add_var( "POPS" , "E,CHEQ" , "FLAG" , "-1/0/1/2 excluded/match/disc5/disc3" );
  add_var( "POPS" , "E,CHEQ" , "CONF" , "Confidence score" );
  add_var( "POPS" , "E,CHEQ" , "PP_N1" , "Posterior probability of N1" );
  add_var( "POPS" , "E,CHEQ" , "PP_N2" , "Posterior probability of N2" );
  add_var( "POPS" , "E,CHEQ" , "PP_N3" , "Posterior probability of N3" );
  add_var( "POPS" , "E,CHEQ" , "PP_R" , "Posterior probability of REM" );
  add_var( "POPS" , "E,CHEQ" , "PP_W" , "Posterior probability of wake" );
  add_var( "POPS" , "E,CHEQ" , "PRED", "Predicted stage" );
  add_var( "POPS" , "E,CHEQ" , "PRIOR" , "Observed stage (if known)" );
  add_var( "POPS" , "E,CHEQ" , "START" , "Start time (hh:mm:ss)" );
  add_var( "POPS" , "E,CHEQ" , "STOP" , "Stop time (hh:mm:ss)");
  
  add_table( "POPS" , "SS,CHEQ" , "Sleep-stage summaries" );
  add_var( "POPS" , "SS,CHEQ" , "OBS" , "Observed stage duration (for included epochs)");
  add_var( "POPS" , "SS,CHEQ" , "ORIG" , "Observed stage duration (all epochs)");
  add_var( "POPS" , "SS,CHEQ" , "PRF" , "Predicted stage duration, weighted" );
  add_var( "POPS" , "SS,CHEQ" , "PR1" , "Predicted stage duration, based on most likely" );
  add_var( "POPS" , "SS,CHEQ" , "F1" , "F1 statistic" );
  add_var( "POPS" , "SS,CHEQ" , "RECALL" , "Recall" );
  add_var( "POPS" , "SS,CHEQ" , "PREC" , "Precision" );

  add_table( "POPS" , "FTR,CHEQ" , "Feature stats" );
  add_var( "POPS" , "FTR,CHEQ" , "BAD" , "Number of bad epochs" );
  add_var( "POPS" , "FTR,CHEQ" , "DROPPED" , "Feature completely dropped" );
  add_var( "POPS" , "FTR,CHEQ" , "PROP" , "Proportion of bad epochs" );


  add_table( "POPS" , "ETYPE,CHEQ" , "Error type" );
  add_var( "POPS" , "ETYPE,CHEQ" , "ACC" , "Accuracy" );
  add_var( "POPS" , "ETYPE,CHEQ" , "N" , "Count" );

  add_table( "POPS" , "ETYPE,CHEQ,SS" , "Error type" );
  add_var( "POPS" , "ETYPE,CHEQ,SS" , "ACC" , "Accuracy" );
  add_var( "POPS" , "ETYPE,CHEQ,SS" , "N" , "Count" );

  add_table( "POPS" , "PRED,OBS,CHEQ" , "Confusion matrix" );
  add_var( "POPS" , "PRED,OBS,CHEQ" , "N" , "Count" );
  add_var( "POPS" , "PRED,OBS,CHEQ" , "P" , "Proportion" );

  //
  // REBASE
  //

  add_cmd( "stage" , "REBASE" , "Translate staging between epoch lengths" );
  add_url( "REBASE" , "staging/#rebase" );
  add_verb( "REBASE" ,
            "Use SOAP posteriors to translate existing staging from the current "
            "epoch definition to a new epoch duration.\n\n"
            "REBASE assumes an explicit prior EPOCH command has already set the "
            "source epoching. It fits a SOAP model in that space, projects "
            "posteriors onto the requested duration, and reports stage summaries "
            "for the rebased epochs." );
  add_param( "REBASE" , "dur" , "30" , "Target epoch duration in seconds" );
  add_param( "REBASE" , "model" , "_1" , "SOAP model to load" );
  add_param( "REBASE" , "sig" , "C4_M1" , "Signal to use for SOAP re-basing" );
  add_param( "REBASE" , "segment-sec" , "4" , "Welch segment length (seconds)" );
  add_param( "REBASE" , "segment-overlap" , "2" , "Welch segment overlap (seconds)" );
  add_table( "REBASE" , "SS" , "Rebased stage summaries" );
  add_var( "REBASE" , "SS" , "PRF" , "Predicted stage duration, weighted" );
  add_var( "REBASE" , "SS" , "PR1" , "Predicted stage duration, hard calls" );

  //
  // PLACE
  //

  add_cmd( "stage" , "PLACE" , "Localize staged epochs on the EDF timeline" );
  add_url( "PLACE" , "staging/#place" );
  add_verb( "PLACE" ,
            "Use SOAP to align an external stage sequence to the EDF timeline.\n\n"
            "PLACE searches over integer epoch offsets, and optionally smaller "
            "micro-shifts, to find where a supplied stage list best fits the EDF. "
            "It reports overlap statistics and kappa-based fit summaries for each "
            "candidate alignment plus the best overall match." );
  add_param( "PLACE" , "stages" , "stages.eannot" , "External stage sequence to place" );
  add_param( "PLACE" , "model" , "_1" , "SOAP model to load" );
  add_param( "PLACE" , "sig" , "C4_M1" , "Signal to use for SOAP placement" );
  add_param( "PLACE" , "3-class" , "" , "Pool N1/N2/N3 into a single NREM class" );
  add_param( "PLACE" , "nr" , "" , "Alias for 3-class placement" );
  add_param( "PLACE" , "out" , "placed.eannot" , "Write the best-fit placed stages" );
  add_param( "PLACE" , "offset" , "0" , "Fix the placement offset in epochs" );
  add_param( "PLACE" , "micro-shift" , "1,30" , "Search smaller temporal shifts: unit,max" );
  add_param( "PLACE" , "edf-overlap" , "0.1" , "Required minimum EDF overlap fraction" );
  add_param( "PLACE" , "stg-overlap" , "0.5" , "Required minimum stage-file overlap fraction" );
  add_param( "PLACE" , "prefix" , "PLACED" , "Prefix for generated annotations" );
  add_table( "PLACE" , "" , "Best placement summary" );
  add_var( "PLACE" , "" , "OFFSET" , "Best offset in epochs" );
  add_var( "PLACE" , "" , "K" , "Best 5-stage kappa" );
  add_var( "PLACE" , "" , "OLAP_N" , "Number of overlapping epochs" );
  add_var( "PLACE" , "" , "OLAP_EDF" , "Fraction of EDF epochs overlapped" );
  add_var( "PLACE" , "" , "OLAP_STG" , "Fraction of supplied stages overlapped" );
  add_var( "PLACE" , "" , "START_EDF" , "First matched EDF epoch" );
  add_var( "PLACE" , "" , "STOP_EDF" , "Last matched EDF epoch" );
  add_var( "PLACE" , "" , "START_STG" , "First matched staged epoch" );
  add_var( "PLACE" , "" , "STOP_STG" , "Last matched staged epoch" );
  add_table( "PLACE" , "OFFSET" , "Placement scan across offsets" );
  add_var( "PLACE" , "OFFSET" , "FIT" , "Whether a valid fit was obtained" );
  add_var( "PLACE" , "OFFSET" , "NS" , "Number of supplied staged epochs used" );
  add_var( "PLACE" , "OFFSET" , "NE" , "Number of EDF epochs used" );
  add_var( "PLACE" , "OFFSET" , "SS" , "Observed stage set at that offset" );
  add_var( "PLACE" , "OFFSET" , "OLAP_N" , "Number of overlapping epochs" );
  add_var( "PLACE" , "OFFSET" , "OLAP_EDF" , "Fraction of EDF epochs overlapped" );
  add_var( "PLACE" , "OFFSET" , "OLAP_STG" , "Fraction of supplied stages overlapped" );
  add_var( "PLACE" , "OFFSET" , "K" , "5-stage kappa at that offset" );
  add_var( "PLACE" , "OFFSET" , "S" , "5-stage kappa scaled to the best fit" );
  add_var( "PLACE" , "OFFSET" , "K3" , "3-stage kappa at that offset" );
  add_var( "PLACE" , "OFFSET" , "S3" , "3-stage kappa scaled to the best fit" );
  add_table( "PLACE" , "MS,OFFSET" , "Placement scan with micro-shifts" );
  add_var( "PLACE" , "MS,OFFSET" , "FIT" , "Whether a valid fit was obtained" );
  add_var( "PLACE" , "MS,OFFSET" , "NS" , "Number of supplied staged epochs used" );
  add_var( "PLACE" , "MS,OFFSET" , "NE" , "Number of EDF epochs used" );
  add_var( "PLACE" , "MS,OFFSET" , "SS" , "Observed stage set at that offset" );
  add_var( "PLACE" , "MS,OFFSET" , "OLAP_N" , "Number of overlapping epochs" );
  add_var( "PLACE" , "MS,OFFSET" , "OLAP_EDF" , "Fraction of EDF epochs overlapped" );
  add_var( "PLACE" , "MS,OFFSET" , "OLAP_STG" , "Fraction of supplied stages overlapped" );
  add_var( "PLACE" , "MS,OFFSET" , "K" , "5-stage kappa at that shift" );
  add_var( "PLACE" , "MS,OFFSET" , "S" , "5-stage kappa scaled to the best fit" );
  add_var( "PLACE" , "MS,OFFSET" , "K3" , "3-stage kappa at that shift" );
  add_var( "PLACE" , "MS,OFFSET" , "S3" , "3-stage kappa scaled to the best fit" );

  //
  // RUN-POPS
  //

  add_cmd( "stage" , "RUN-POPS" , "High-level POPS wrapper" );
  add_url( "RUN-POPS" , "staging/#run-pops" );
  add_verb( "RUN-POPS" ,
            "Run the standard Luna POPS preprocessing and scoring pipeline in one "
            "command.\n\n"
            "RUN-POPS copies the requested signals, optionally re-references them, "
            "resamples to 128 Hz, band-pass filters, normalizes, optionally runs "
            "EDGER, and then invokes POPS with the assembled temporary signals." );
  add_param( "RUN-POPS" , "sig" , "C3,C4" , "Primary EEG signal(s)" );
  add_param( "RUN-POPS" , "ref" , "M2,M1" , "Reference signal(s), matching sig length" );
  add_param( "RUN-POPS" , "args" , "trim=10 3-class" , "Additional arguments passed to POPS" );
  add_param( "RUN-POPS" , "ignore-obs" , "F" , "Ignore existing observed staging" );
  add_param( "RUN-POPS" , "lib" , "s2" , "POPS library root" );
  add_param( "RUN-POPS" , "path" , "." , "Base path for POPS resources" );
  add_param( "RUN-POPS" , "filter" , "T" , "Band-pass filter copied signals before POPS" );
  add_param( "RUN-POPS" , "edger" , "T" , "Run EDGER on the copied signals" );

  //
  // EVAL-STAGES
  //

  add_cmd( "stage" , "EVAL-STAGES" , "Evaluate staging against a reference" );
  add_url( "EVAL-STAGES" , "staging/#eval-stages" );
  add_verb( "EVAL-STAGES" ,
            "Evaluate the current staging against an external stage file using the "
            "same summary machinery as POPS.\n\n"
            "This command is useful when you already have a stage sequence and want "
            "agreement metrics, stage summaries, confusion tables, and related "
            "diagnostics without rerunning a full staging model." );
  add_param( "EVAL-STAGES" , "file" , "stages.eannot" , "Reference stage file to compare against" );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // TIME/FREQUENCY ANALYSIS
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // --cwt
  //

  add_cmd( "helpers" , "--cwt" , "Or --cwt-design" );

  //
  // EMD
  //

  add_cmd( "power" , "EMD" , "Empirical mode decomposition" );
  add_url( "EMD" , "power-spectra/#emd" );
  add_verb( "EMD" ,
            "Decompose a signal into intrinsic mode functions using empirical mode "
            "decomposition.\n\n"
            "EMD is a data-adaptive method that iteratively sifts the signal into "
            "oscillatory components without imposing fixed frequency bands. The "
            "resulting IMFs can be interpreted as progressively slower modes plus a "
            "residual trend." );
  add_param( "EMD" , "sig" , "C3,C4" , "Select signals for EMD" );
  add_param( "EMD" , "tag" , "_C_" , "IMF channel tag, if not _IMF_" );
  add_param( "EMD" , "sift" , "20" , "Maximum number of sifting operations" );
  add_param( "EMD" , "imf" , "10" , "Maximum number of IMF to extract" );

  //
  // ACF
  //

  add_cmd( "power" , "ACF" , "Autocorrelation function" );
  add_url( "ACF" , "power-spectra/#acf" );
  add_verb( "ACF" ,
            "Estimate the autocorrelation structure of a signal.\n\n"
            "This method summarizes how strongly the signal resembles lagged "
            "versions of itself across time offsets, which can be useful for "
            "detecting rhythmic structure, regularity, and dominant timescales." );
  add_param( "ACF" , "sig" , "C3,C4" , "Restrict analysis to these channels" );
  add_param( "ACF" , "lag" , "100" , "Maximum lag in samples" );
  add_table( "ACF" , "CH,LAG" , "Autocorrelation by channel and lag" );
  add_var( "ACF" , "CH,LAG" , "ACF" , "Autocorrelation at this lag" );
  add_var( "ACF" , "CH,LAG" , "SEC" , "Lag in seconds" );

  //
  // DFA
  //

  add_cmd( "power" , "DFA" , "Detrended fluctuation analysis" );
  add_verb( "DFA" ,
            "DFA runs a Fourier-domain detrended fluctuation analysis following "
            "the implementation of Nolte et al. (2019). Rather than returning a "
            "single exponent, it evaluates the fluctuation function over a grid of "
            "time scales and reports the local slope at each scale.\n\n"
            "The command can be applied to the whole trace or epoch by epoch. It "
            "can also first apply a narrowband filter-Hilbert transform and then "
            "run DFA either on the filtered waveform or on its amplitude envelope, "
            "making it useful for both broadband and band-limited scaling "
            "analyses." );

  add_param( "DFA" , "sig" , "C3,C4" , "Signals to analyze" );
  add_param( "DFA" , "n" , "100" , "Number of analysis windows/scales" );
  add_param( "DFA" , "min" , "0.1" , "Minimum window size in seconds" );
  add_param( "DFA" , "m" , "2" , "Log-scale span exponent for the window grid" );
  add_param( "DFA" , "epoch" , "" , "Run DFA separately for each epoch" );
  add_param( "DFA" , "alpha-lwr" , "0.5" , "Lower time scale, in seconds, for the summary alpha fit" );
  add_param( "DFA" , "alpha-upr" , "10" , "Upper time scale, in seconds, for the summary alpha fit" );
  add_param( "DFA" , "f-lwr" , "0.5" , "Lower frequency for optional narrowband filtering" );
  add_param( "DFA" , "f-upr" , "4" , "Upper frequency for optional narrowband filtering" );
  add_param( "DFA" , "ripple" , "0.01" , "FIR ripple for narrowband filtering" );
  add_param( "DFA" , "tw" , "1" , "FIR transition width for narrowband filtering" );
  add_param( "DFA" , "envelope" , "F" , "Use the filtered signal rather than the Hilbert envelope if set to F" );

  add_table( "DFA" , "CH,SEC" , "Whole-trace DFA results by time scale" );
  add_var( "DFA" , "CH,SEC" , "FLUCT" , "Fluctuation magnitude at this time scale" );
  add_var( "DFA" , "CH,SEC" , "SLOPE" , "Local DFA slope at this time scale" );

  add_table( "DFA" , "CH" , "Whole-trace summary DFA fit" );
  add_var( "DFA" , "CH" , "ALPHA" , "Summary log-log slope of fluctuation versus scale" );
  add_var( "DFA" , "CH" , "R2" , "R-squared of the summary log-log fit" );

  add_table( "DFA" , "CH,E,SEC" , "Epoch-level DFA results by time scale" );
  add_var( "DFA" , "CH,E,SEC" , "FLUCT" , "Fluctuation magnitude at this time scale" );
  add_var( "DFA" , "CH,E,SEC" , "SLOPE" , "Local DFA slope at this time scale" );

  add_table( "DFA" , "CH,E" , "Epoch-level summary DFA fit" );
  add_var( "DFA" , "CH,E" , "ALPHA" , "Summary log-log slope of fluctuation versus scale" );
  add_var( "DFA" , "CH,E" , "R2" , "R-squared of the summary log-log fit" );

  //
  // PSD
  //

  add_cmd( "power"   , "PSD" , "Power spectral density estimation (Welch)" );
  add_url( "PSD" , "power-spectra/#psd" );
  add_verb( "PSD" ,
            "Estimate power spectra and band power using Welch's method.\n\n"
            "The signal is divided into overlapping segments, each segment is "
            "windowed and Fourier transformed, and the resulting periodograms are "
            "averaged to obtain a more stable spectral estimate than a single FFT. "
            "Luna can then summarize spectra as absolute power, relative power, "
            "spectral slope, per-epoch spectra, and band-level dynamics." );

  add_param( "PSD" , "sig" , "C3,C4" , "Restrict analysis to these channels" );
  add_param( "PSD" , "epoch" , "" , "Calculate per-epoch band power" );
  add_param( "PSD" , "max" , "100" , "Specify max frequency for power spectra" );
  add_param( "PSD" , "bin" , "1" , "Specify bin-size for power spectra" );
  add_param( "PSD" , "spectrum" , "" , "Calculate power spectra" );
  add_param( "PSD" , "epoch-spectrum" , "" , "Calculate per-epoch power spectra" );
  add_param( "PSD" , "dB" , "" , "Report power in decibel units" );
  add_param( "PSD" , "peaks" , "" , "Estimate of spectral peaks/artifacts" );

  add_param( "PSD" , "no-window" , "" , "No windowing on FFT segments" );
  add_param( "PSD" , "hann"    , "" , "Use Hann window" );
  add_param( "PSD" , "hamming" , "" , "Use Hamming window" );
  add_param( "PSD" , "tukey50" , "" , "Use Tukey(50%) window (default)" );

  add_param( "PSD" , "average-adj" , "" , "Average adjacent frequency bins" );
  
  add_param( "PSD" , "dynamics" , "" , "Power dynamics (experimental/undocumented)" );

  add_param( "PSD" , "kurtosis" , "" , "Output kurtosis for band-power (dB-scaled)" );
    
  
  add_table( "PSD" , "CH" , "Channel-level output" );
  add_var( "PSD" , "CH" , "NE" , "Number of epochs" );
  add_var( "PSD" , "CH" , "KURT" , "Peak (PSD kurtosis)" );
  add_var( "PSD" , "CH" , "SPK" , "Sum PSD peakedness" );
  add_var( "PSD" , "CH" , "SPEC_SLOPE" , "Spectral slope" );
  add_var( "PSD" , "CH" , "SPEC_SLOPE_N" , "Spectral slope number of points" );
  add_var( "PSD" , "CH" , "SPEC_SLOPE_MD" , "Spectral slope (median)" );
  add_var( "PSD" , "CH" , "SPEC_SLOPE_MN" , "Spectral slope (mean over epochs)" );
  add_var( "PSD" , "CH" , "SPEC_SLOPE_SD" , "Spectral slope (SD over epochs)" );

  add_table( "PSD" , "CH,B" , "Whole-night, per-channel band power" );
  add_var( "PSD" , "CH,B" , "PSD" , "Power" );
  add_var( "PSD" , "CH,B" , "RELPSD" , "Relative power" );
  add_var( "PSD" , "CH,B" , "KURT" , "Kurtosis" );

  add_table( "PSD" , "CH,B1,B2" , "Whole-night, per-channel band power ratios" );
  add_var( "PSD" , "CH,B1,B2" , "RATIO" , "Band power ratios" );
  
  add_table( "PSD" , "CH,F" , "Whole-night, per-channel power" );
  add_var( "PSD" , "CH,F" , "PSD" , "Power (mean over epochs)" );
  add_var( "PSD" , "CH,F" , "PSD_MD" , "Power (median over epochs)" );
  add_var( "PSD" , "CH,F" , "PSD_SD" , "Power (SD over epochs)" );
  add_var( "PSD" , "CH,F" , "PSD_CV" , "Power (CV over epochs)" );
  add_var( "PSD" , "CH,F" , "SEGCV_MN" , "Segment CV (mean)" );
  add_var( "PSD" , "CH,F" , "SEGCV_MD" , "Segment CV (median)" );
  add_var( "PSD" , "CH,F" , "SEGCV_SD" , "Segment CV (SD)" );

  add_table( "PSD" , "CH,B,E" , "Whole-night, per-channel per-epoch band power" );
  add_var( "PSD" , "CH,B,E" , "PSD" , "Power" );
  add_var( "PSD" , "CH,B,E" , "RELPSD" , "Relative power" );

  add_table( "PSD" , "CH,F,E" , "Whole-night, per-channel per-epoch power" );
  add_var( "PSD" , "CH,F,E" , "PSD" , "Power" );
  add_var( "PSD" , "CH,F,E" , "CV" , "CV" );
  set_compressed( "PSD" , tfac_t( "CH,F,E" ) );

  add_table( "PSD" , "CH,E", "Epoch/channel level stats" );
  add_var( "PSD" , "CH,E" , "KURT" , "Peak (PSD kurtosis)" );
  add_var( "PSD" , "CH,E" , "SPK" , "Sum PSD peakedness" );
  add_var( "PSD" , "CH,E" , "SPEC_SLOPE" , "Spectral slope" );
  add_var( "PSD" , "CH,E" , "SPEC_SLOPE_N" , "Spectral slope number of points" );
  
  // dynamics
  add_table( "PSD", "CH,F,VAR,QD" , "PSD spectra dynamics" );
  add_var( "PSD" , "CH,F,VAR,QD" , "A_P2P" , "Peak-to-peak amplitude of the dynamics profile" );
  add_var( "PSD" , "CH,F,VAR,QD" , "T_P2P" , "Epoch distance from trough to peak in the dynamics profile" );
  add_var( "PSD" , "CH,F,VAR,QD" , "U" , "Linear trend correlation with epoch order" );
  add_var( "PSD" , "CH,F,VAR,QD" , "U2" , "Quadratic trend coefficient with epoch order" );
  add_var( "PSD" , "CH,F,VAR,QD" , "MEAN" , "Mean of the normalized smoothed dynamics profile" );
  add_var( "PSD" , "CH,F,VAR,QD" , "OMEAN" , "Mean of the original smoothed dynamics profile" );
  add_var( "PSD" , "CH,F,VAR,QD" , "SD" , "Standard deviation of the dynamics profile" );
  add_var( "PSD" , "CH,F,VAR,QD" , "N" , "Number of epochs or sections contributing to the profile" );

  add_table( "PSD", "CH,B,VAR,QD" , "PSD band dynamics" );
  add_var( "PSD" , "CH,B,VAR,QD" , "A_P2P" , "Peak-to-peak amplitude of the band dynamics profile" );
  add_var( "PSD" , "CH,B,VAR,QD" , "T_P2P" , "Epoch distance from trough to peak in band dynamics" );
  add_var( "PSD" , "CH,B,VAR,QD" , "U" , "Linear trend correlation with epoch order" );
  add_var( "PSD" , "CH,B,VAR,QD" , "U2" , "Quadratic trend coefficient with epoch order" );
  add_var( "PSD" , "CH,B,VAR,QD" , "MEAN" , "Mean of the normalized smoothed band-dynamics profile" );
  add_var( "PSD" , "CH,B,VAR,QD" , "OMEAN" , "Mean of the original smoothed band-dynamics profile" );
  add_var( "PSD" , "CH,B,VAR,QD" , "SD" , "Standard deviation of the band-dynamics profile" );
  add_var( "PSD" , "CH,B,VAR,QD" , "N" , "Number of epochs or sections contributing to the profile" );

  add_table( "PSD", "CH,F,VAR,QD,Q" , "PSD spectra dynamics quantiles" );
  add_var( "PSD" , "CH,F,VAR,QD,Q" , "OS" , "Original smoothed dynamics value at this quantile" );
  add_var( "PSD" , "CH,F,VAR,QD,Q" , "SS" , "Normalized smoothed dynamics value at this quantile" );

  add_table( "PSD", "CH,B,VAR,QD,Q" , "PSD band dynamics quantiles" );
  add_var( "PSD" , "CH,B,VAR,QD,Q" , "OS" , "Original smoothed band-dynamics value at this quantile" );
  add_var( "PSD" , "CH,B,VAR,QD,Q" , "SS" , "Normalized smoothed band-dynamics value at this quantile" );

  //
  // PCOUPL
  //

  add_cmd( "power" , "PCOUPL" , "Generic phase/event coupling analysis" );
  add_url( "PCOUPL" , "power-spectra/#pcoupl" );
  add_verb( "PCOUPL" ,
            "Quantify how discrete events align to the phase of an oscillatory "
            "signal.\n\n"
            "PCOUPL band-pass filters the signal, extracts instantaneous phase via "
            "a Hilbert transform, and then evaluates whether event anchors occur at "
            "preferred phases more often than expected under permutation. It is a "
            "general framework for phase locking between annotations and ongoing "
            "rhythm." );

  add_param( "PCOUPL" , "sig" , "C3,C4" , "Signals" );
  add_param( "PCOUPL" , "events" , "arousal" , "One or more annotation classes" );  
  add_param( "PCOUPL" , "lwr" , "3" , "Lower frequency for filter-Hilbert" );
  add_param( "PCOUPL" , "upr" , "8" , "Upper frequency for filter-Hilbert" );
  add_param( "PCOUPL" , "anchor" , "start" , "Optional anchor (start/middle/stop)" );
  add_param( "PCOUPL" , "nreps" , "1000" , "Number of permutations" );
  add_param( "PCOUPL" , "tw" , "0.5" , "Transition width for Kaiser window FIR" );
  add_param( "PCOUPL" , "ripple" , "0.01" , "Ripple for Kaiser window FIR" );
  add_param( "PCOUPL" , "perm-whole-trace" , "" , "Permute signals across whole recording (not within epoch)" );  
  add_param( "PCOUPL" , "fixed-epoch-dur" , "20" , "If using generic epochs, set a fixed epoch size for permutation" );
  
  add_table( "PCOUPL" , "ANNOT,CH" , "Phase coupling statistics" );
  add_var( "PCOUPL" , "ANNOT,CH" , "ANGLE" , "Mean phase angle (degrees)" );
  add_var( "PCOUPL" , "ANNOT,CH" , "MAG", "Coupling magnitude (observed statistic)" );
  add_var( "PCOUPL" , "ANNOT,CH" , "MAG_Z", "Permutation-based Z-score for coupling magnitude" );
  add_var( "PCOUPL" , "ANNOT,CH" , "MAG_NULL", "Mean coupling statistic under the null" );
  add_var( "PCOUPL" , "ANNOT,CH" , "MAG_EMP", "Empirical p-value" );
  add_var( "PCOUPL" , "ANNOT,CH" , "PV", "Asymptotic p-value" );
  add_var( "PCOUPL" , "ANNOT,CH" , "SIGPV_NULL", "Proportion of asymptotic p<0.05 under the null" );
  add_var( "PCOUPL" , "ANNOT,CH" , "N" , "Number of events" );

  add_table( "PCOUPL" , "ANNOT,CH,PHASE" , "Phase-bin overlap statistics" );
  add_var( "PCOUPL" , "ANNOT,CH,PHASE" , "OVERLAP", "Observed count of event anchors per signal phase bin" );
  add_var( "PCOUPL" , "ANNOT,CH,PHASE" , "OVERLAP_EXP", "Expected count based on permutations" );
  add_var( "PCOUPL" , "ANNOT,CH,PHASE" , "OVERLAP_EMP", "Empirical p-value based on permutations" );
  add_var( "PCOUPL" , "ANNOT,CH,PHASE" , "OVERLAP_Z", "Z-score based on permutations" );

  //
  // ASYMM
  //

  add_cmd( "exp" , "ASYMM" , "EEG asymmetry" , true );
  add_url( "ASYMM" , "power-spectra/#asymm" );
  add_verb( "ASYMM" ,
            "Compare left-versus-right spectral power across sleep states, cycles, "
            "and transitions.\n\n"
            "ASYMM relies on spectral estimates, typically from cached PSD output, "
            "and summarizes asymmetry as log-ratios between homologous channels. It "
            "can then contrast those asymmetries across wake, NREM, REM, NREM "
            "cycles, and transition-centered windows." );

  add_param( "ASYMM" , "left" , "C3" , "Left channel(s)" );
  add_param( "ASYMM" , "right" , "C4" , "Right channel(s)" );
  add_param( "ASYMM" , "nreps" , "500" , "Replicates for transition randomisation test" );
  add_param( "ASYMM" , "cache-var" , "PER" , "Cached variable (if not PSD)" );
  add_param( "ASYMM" , "epoch" , "" , "Epoch level output" );
  add_param( "ASYMM" , "trans" , "" , "Transition-centric output" );
  
  add_table( "ASYMM" , "B,CHS" , "Band-based primary asymmetry stats" );
  add_var( "ASYMM", "B,CHS" , "L_SLEEP" , "Left power during sleep" );
  add_var( "ASYMM", "B,CHS" , "R_SLEEP" , "Right power during sleep" );
  add_var( "ASYMM", "B,CHS" , "LR_SLEEP" , "log2(L/R) during sleep" );
  add_var( "ASYMM", "B,CHS" , "LR_WAKE" , "log2(L/R) during wake" );
  add_var( "ASYMM", "B,CHS" , "LR_NREM" , "log2(L/R) during NREM" );
  add_var( "ASYMM", "B,CHS" , "LR_REM" , "log2(L/R) during REM" );
  add_var( "ASYMM", "B,CHS" , "Z_REM" , "NREM-normalized REM log2(L/R)");
  add_var( "ASYMM", "B,CHS" , "LOGP" , "NREM-normalized REM log2(L/R) -log10(p)");
  add_var( "ASYMM", "B,CHS" , "ABS_Z_REM" , "Absolute NREM-normalized REM log2(L/R)");
  add_var( "ASYMM", "B,CHS" , "ABS_LOGP" , "Absolute NREM-normalized REM log2(L/R) -log10(p)");
  add_var( "ASYMM", "B,CHS" , "NC" , "Number of included NREM cycles" );
  add_var( "ASYMM", "B,CHS" , "TR_NR2R_N" , "Number of NR-to-R transitions" );
  add_var( "ASYMM", "B,CHS" , "TR_R2NR_N" , "Number of R-to-NR transitions" );

  add_table( "ASYMM" , "F,CHS" , "Frequency-bin-based primary asymmetry stats" );
  add_var( "ASYMM", "F,CHS" , "L_SLEEP" , "Left power during sleep" );
  add_var( "ASYMM", "F,CHS" , "R_SLEEP" , "Right power during sleep" );
  add_var( "ASYMM", "F,CHS" , "LR_SLEEP" , "log2(L/R) during sleep" );
  add_var( "ASYMM", "F,CHS" , "LR_WAKE" , "log2(L/R) during wake" );
  add_var( "ASYMM", "F,CHS" , "LR_NREM" , "log2(L/R) during NREM" );
  add_var( "ASYMM", "F,CHS" , "LR_REM" , "log2(L/R) during REM" );
  add_var( "ASYMM", "F,CHS" , "Z_REM" , "NREM-normalized REM log2(L/R)");
  add_var( "ASYMM", "F,CHS" , "LOGP" , "NREM-normalized REM log2(L/R) -log10(p)");
  add_var( "ASYMM", "F,CHS" , "ABS_Z_REM" , "Absolute NREM-normalized REM log2(L/R)");
  add_var( "ASYMM", "F,CHS" , "ABS_LOGP" , "Absolute NREM-normalized REM log2(L/R) -log10(p)");
  add_var( "ASYMM", "F,CHS" , "NC" , "Number of included NREM cycles" );
  add_var( "ASYMM", "F,CHS" , "TR_NR2R_N" , "Number of NR-to-R transitions" );
  add_var( "ASYMM", "F,CHS" , "TR_R2NR_N" , "Number of R-to-NR transitions" );

  add_table( "ASYMM" , "E,B,CHS" , "Epoch-level frequency-band output" );
  add_var( "ASYMM", "E,B,CHS" , "C" , "NREM cycle number" );
  add_var( "ASYMM", "E,B,CHS" , "INC" , "Included in analysis, 1=Y");
  add_var( "ASYMM", "E,B,CHS" , "CONSIDER" , "Considered this epoch" );
  add_var( "ASYMM", "E,B,CHS" , "L" , "Left power" );
  add_var( "ASYMM", "E,B,CHS" , "R" , "Right power" );
  add_var( "ASYMM", "E,B,CHS" , "LR" , "log2(L/R)" );
  add_var( "ASYMM", "E,B,CHS" , "OUT" , "Flagged as outlier, 1=Y" );
  add_var( "ASYMM", "E,B,CHS" , "SS" , "Sleep stage (W/R/NR)" );

  add_table( "ASYMM" , "E,F,CHS" , "Epoch-level frequency-bin output" );
  add_var( "ASYMM", "E,F,CHS" , "C" , "NREM cycle number" );
  add_var( "ASYMM", "E,F,CHS" , "INC" , "Included in analysis, 1=Y");
  add_var( "ASYMM", "E,F,CHS" , "CONSIDER" , "Considered this epoch" );
  add_var( "ASYMM", "E,F,CHS" , "L" , "Left power" );
  add_var( "ASYMM", "E,F,CHS" , "R" , "Right power" );
  add_var( "ASYMM", "E,F,CHS" , "LR" , "log2(L/R)" );
  add_var( "ASYMM", "E,F,CHS" , "OUT" , "Flagged as outlier, 1=Y" );
  add_var( "ASYMM", "E,F,CHS" , "SS" , "Sleep stage (W/R/NR)" );

  add_table( "ASYMM" , "C,B,CHS" , "Cycle-level frequency-band output" );
  add_var( "ASYMM", "C,B,CHS" , "LR_REM" , "log2(L/R) in REM" );
  add_var( "ASYMM", "C,B,CHS" , "LR_NREM" , "log2(L/R) in NREM" );
  add_var( "ASYMM", "C,B,CHS" , "LR_LEADING_NREM" , "log2(L/R) in leading NREM" );
  add_var( "ASYMM", "C,B,CHS" , "LR_TRAILING_NREM" , "log2(L/R) in trailing NREM" );
  add_var( "ASYMM", "C,B,CHS" , "Z_REM" , "Normalized REM log2(L/R)" );
  add_var( "ASYMM", "C,B,CHS" , "P" , "REM-vs-NREM p-value" );
  add_var( "ASYMM", "C,B,CHS" , "P_NREM" , "Leading-vs-trailing NREM p-value" );
  add_var( "ASYMM", "C,B,CHS" , "LOGP" , "REM-vs-NREM p-value, log-scaled" );
  add_var( "ASYMM", "C,B,CHS" , "LOGP_NREM" , "Leading-vs-trailing NREM p-value, log-scaled" );
  add_var( "ASYMM", "C,B,CHS" , "N_NREM" , "Number of NREM epochs" );
  add_var( "ASYMM", "C,B,CHS" , "N_REM" , "Number of REM epochs" );
  add_var( "ASYMM", "C,B,CHS" , "INC" , "Includede this cycle" );
  
  add_table( "ASYMM" , "C,F,CHS" , "Cycle-level frequency-bin output" );
  add_var( "ASYMM", "C,F,CHS" , "LR_REM" , "log2(L/R) in REM" );
  add_var( "ASYMM", "C,F,CHS" , "LR_NREM" , "log2(L/R) in NREM" );
  add_var( "ASYMM", "C,F,CHS" , "LR_LEADING_NREM" , "log2(L/R) in leading NREM" );
  add_var( "ASYMM", "C,F,CHS" , "LR_TRAILING_NREM" , "log2(L/R) in trailing NREM" );
  add_var( "ASYMM", "C,F,CHS" , "Z_REM" , "Normalized REM log2(L/R)" );
  add_var( "ASYMM", "C,F,CHS" , "P" , "REM-vs-NREM p-value" );
  add_var( "ASYMM", "C,F,CHS" , "P_NREM" , "Leading-vs-trailing NREM p-value" );
  add_var( "ASYMM", "C,F,CHS" , "LOGP" , "REM-vs-NREM p-value, log-scaled" );
  add_var( "ASYMM", "C,F,CHS" , "LOGP_NREM" , "Leading-vs-trailing NREM p-value, log-scaled" );
  add_var( "ASYMM", "C,F,CHS" , "N_NREM" , "Number of NREM epochs" );
  add_var( "ASYMM", "C,F,CHS" , "N_REM" , "Number of REM epochs" );
  add_var( "ASYMM", "C,F,CHS" , "INC" , "Includede this cycle" );

  add_table( "ASYMM" , "B,CHS,TR" , "Transition-based frequency-bin output" );
  add_var( "ASYMM", "B,CHS,TR" , "NR2R" , "NREM-to-REM transition means" );
  add_var( "ASYMM", "B,CHS,TR" , "R2NR" , "REM-to-NREM transition means" );
  add_var( "ASYMM", "B,CHS,TR" , "NR2R_Z" , "NREM-to-REM transition means, Z-score" );
  add_var( "ASYMM", "B,CHS,TR" , "R2NR_Z" , "REM-to-NREM transition means, Z-score" );
  
  add_table( "ASYMM" , "F,CHS,TR" , "Transition-based frequency-bin output" );
  add_var( "ASYMM", "F,CHS,TR" , "NR2R" , "NREM-to-REM transition means" );
  add_var( "ASYMM", "F,CHS,TR" , "R2NR" , "REM-to-NREM transition means" );
  add_var( "ASYMM", "F,CHS,TR" , "NR2R_EMP" , "NREM-to-REM transition means, null empirical expectation" );
  add_var( "ASYMM", "F,CHS,TR" , "R2NR_EMP" , "REM-to-NREM transition means, null empirical expectation" );

  add_table( "ASYMM" , "E" , "Epoch-level hypnogram statistic" );
  add_var( "ASYMM", "E" , "STAGE" , "Sleep stage" );
  add_var( "ASYMM", "E" , "STAGE" , "Sleep stage (numeric encoding)" );
  add_var( "ASYMM", "E" , "STAGE" , "Original sleep stage" );
  add_var( "ASYMM", "E" , "STAGE" , "Elapsed minutes" );
  add_var( "ASYMM", "E" , "CLOCK_TIME" , "Clock time hh:mm:ss" );

  //
  // FIP
  //

  add_cmd( "power" , "FIP" , "Frequency/interval transformation" );
  add_url( "FIP" , "psc/#fip" );
  add_verb( "FIP" ,
            "Map signal power into a joint frequency-by-interval representation.\n\n"
            "FIP uses a continuous-wavelet style representation to summarize how "
            "oscillatory structure varies jointly over frequency and characteristic "
            "time interval or cycle count. It is useful when a simple band-power "
            "summary is too coarse and the target phenomenon spans both time and "
            "frequency scales." );

  add_param( "FIP" , "sig" , "C3,C4" , "Channels to analyse" );
  add_param( "FIP" , "t-lwr" , "0.1" , "Lower time bound" );
  add_param( "FIP" , "t-upr" , "4" , "Upper time bound" );
  add_param( "FIP" , "t-inc" , "0.1" , "Time increment" );

  add_param( "FIP" , "f-lwr" , "0.1" , "Lower frequency bound" );
  add_param( "FIP" , "f-upr" , "0.1" , "Lower frequency bound" );
  add_param( "FIP" , "f-inc" , "0.1" , "Frequecy increment (Hz), linear" );
  add_param( "FIP" , "f-log" , "20" , "Instead of f-inc, uniform on log scale, e.g. 20 steps" );


  add_param( "FIP" , "by-cycles" , "" , "Use cycles instead of time" );
  add_param( "FIP" , "c-lwr" , "1" , "Lower cycle value" );
  add_param( "FIP" , "c-upr" , "7" , "Upper cycle value" );
  add_param( "FIP" , "c-inc" , "1" , "Cycle increment" );
  add_param( "FIP" , "cycles" , "12" , "Set CWT cycles" );

  add_param( "FIP" , "th" , "2" , "Set Z-scale (CWT) threshold" );
  add_param( "FIP" , "log" , "" , "Log-scale Z (CWT)" );
  add_param( "FIP" , "norm" , "" , "Set CWT cycles" );

  add_table( "FIP" , "CH,TBIN,F" , "F/I plot" );
  add_var( "FIP", "CH,TBIN,F" , "FIP" , " FIP value" );
  add_var( "FIP", "CH,TBIN,F" , "ZIP" , " ZIP value" );

  //
  // FFT
  //

  add_cmd( "power"   , "FFT" , "Discrete Fourier Transform" );
  add_url( "FFT" , "power-spectra/#fft" );
  add_verb( "FFT" ,
            "Compute the direct discrete Fourier transform of the selected signals.\n\n"
            "Unlike Welch-style PSD estimation, FFT here exposes the raw spectral "
            "decomposition of the current trace, including real and imaginary "
            "coefficients as well as amplitude and PSD-derived summaries. It is most "
            "useful when you want the transform itself rather than a smoothed "
            "spectral estimate." );
  add_param( "FFT" , "sig" , "C3,C4" , "Channels to analyse" );
  add_param( "FFT" , "verbose" , "" , "Additional output variables" );

  add_table( "FFT" , "CH,F", "Channel-wise power spectra" );
  add_var( "FFT" , "CH,F" , "PSD" , "Power spectral density" );
  add_var( "FFT" , "CH,F" , "DB" , "10log10(PSD)" );
  add_var( "FFT" , "CH,F" , "RE" , "Real value of transform" );
  add_var( "FFT" , "CH,F" , "IM" , "Imaginary value of transform" );
  add_var( "FFT" , "CH,F" , "UNNORM_AMP" , "Unnormalized amplitude" );
  add_var( "FFT" , "CH,F" , "NORM_AMP" , "Normalized amplitude" );

  //
  // MTM
  //

  add_cmd( "power"   , "MTM" , "Power spectral density estimation (Welch)" );
  add_url( "MTM" , "power-spectra/#mtm" );
  add_verb( "MTM" ,
            "Estimate spectra with the multitaper method.\n\n"
            "MTM replaces a single taper with a set of orthogonal Slepian tapers, "
            "computes a spectrum under each taper, and combines them for a lower-"
            "variance estimate with good leakage control. Luna can summarize the "
            "result as whole-record spectra, band power, segment-level output, and "
            "spectral-shape statistics." );

  add_param( "MTM" , "sig" , "C3,C4" , "Restrict analysis to these channels" );
  add_param( "MTM" , "min" , "0.5" ,   "Lower frequency range" );
  add_param( "MTM" , "max" , "100" ,   "Upper frequency range" );
  add_param( "MTM" , "segment-sec" , "30" , "Segment size, seconds" );
  add_param( "MTM" , "segment-inc" , "30" ,   "Segment step, seconds" );
  add_param( "MTM" , "dB" , "" , "Decibel scale output" );
  add_param( "MTM" , "epoch" , "" , "Report per-epoch statistics" );
  
  add_table( "MTM" , "CH" , "Whole-night, per-channel stats" );
  add_var( "MTM" , "CH" , "SPEC_SLOPE" , "Spectral slope" );
  add_var( "MTM" , "CH" , "SPEC_SLOPE_N" , "Spectral slope number of points" );
  add_var( "MTM" , "CH" , "SPEC_SLOPE_MD" , "Spectral slope (median)" );
  add_var( "MTM" , "CH" , "SPEC_SLOPE_MN" , "Spectral slope (mean over epochs)" );
  add_var( "MTM" , "CH" , "SPEC_SLOPE_SD" , "Spectral slope (SD over epochs)" );
  add_var( "MTM" , "CH" , "WREL_PK_FREQ" , "WREL peak (frequency)" );
  add_var( "MTM" , "CH" , "WREL_PK_AMPL" , "WREL peak (amplitude)" );
  add_var( "MTM" , "CH" , "WMTM_PK_FREQ" , "WMTM peak (frequency)" );
  add_var( "MTM" , "CH" , "WMTM_PK_AMPL" , "WMTM peak (amplitude)" );
  add_var( "MTM" , "CH" , "MTM_PK_FREQ" , "MTM peak (frequency)" );
  add_var( "MTM" , "CH" , "MTM_PK_AMPL" , "MTM peak (amplitude)" );


  add_table( "MTM", "CH,SEG", "Segment timing details" );
  add_var( "MTM" , "CH,SEG" , "START" , "Start time (seconds)" );
  add_var( "MTM" , "CH,SEG" , "STOP" , "Stop time (seconds)");
  add_var( "MTM" , "CH,SEG" , "DISC" , "Spans a discontinuity (0/1=N/Y)");
  
  add_table( "MTM" , "CH,F" , "Whole-night, per-channel power" );
  add_var( "MTM" , "CH,F" , "MTM" , "Power" );
  add_var( "MTM" , "CH,F" , "MTM_MD" , "Median power" );
  add_var( "MTM" , "CH,F" , "MTM_SD" , "Power SD" );
  add_var( "MTM" , "CH,F" , "WMTM" , "Weighted power (variable epoch size)" );
  add_var( "MTM" , "CH,F" , "WREL" , "Weighted relative power (variable epoch size)" );


  add_table( "MTM" , "CH,B1,B2" , "Whole-night, per-channel bandpower ratios" );
  add_var( "MTM" , "CH,B1,B2" , "RATIO" , "Band power ratio" );
  add_var( "MTM" , "CH,B1,B2" , "RATIO_MD" , "Median band power ratio" );
  add_var( "MTM" , "CH,B1,B2" , "RATIO_SD" , "Band power ratio SD" );

  add_table( "MTM" , "CH,B" , "Whole-night, per-channel bandpower" );
  add_var( "MTM" , "CH,B" , "MTM" , "Band power" ); 
  add_var( "MTM" , "CH,B" , "MTM_MD" , "Median band power" ); 
  add_var( "MTM" , "CH,B" , "MTM_SD" , "Band power SD" ); 
  add_var( "MTM" , "CH,B" , "REL" , "Relative band power" ); 
  add_var( "MTM" , "CH,B" , "REL_MD" , "Median relative band power" ); 
  add_var( "MTM" , "CH,B" , "REL_SD" , "Relative band power SD" ); 
  add_var( "MTM" , "CH,B" , "SPECCV" , "Spectral CV" );
  add_var( "MTM" , "CH,B" , "SPECCV_MD" , "Median spectral CV" );
  add_var( "MTM" , "CH,B" , "SPECKURT" , "Spectral kurtosis" );
  add_var( "MTM" , "CH,B" , "SPECKURT_MD" , "Median spectral kurtosis" );
  add_var( "MTM" , "CH,B" , "SPECSKEW" , "Spectral skewness" );
  add_var( "MTM" , "CH,B" , "SPECSKEW_MD" , "Median spectral skewness" );

  add_table( "MTM" , "B" , "Whole-night, misc channel-averaged band metrics" );
  add_var( "MTM" , "B" , "SPECCV" , "Spectral CV" );
  add_var( "MTM" , "B" , "SPECCV_MD" , "Median spectral CV" );
  add_var( "MTM" , "B" , "SPECKURT" , "Spectral kurtosis" );
  add_var( "MTM" , "B" , "SPECKURT_MD" , "Median spectral kurtosis" );
  add_var( "MTM" , "B" , "SPECSKEW" , "Spectral skewness" );
  add_var( "MTM" , "B" , "SPECSKEW_MD" , "Median spectral skewness" );

  add_table( "MTM" , "CH,SP,TAPER" , "Taper coefficients" );
  add_var( "MTM" , "CH,SP,TAPER" , "W" , "Weight" );

  add_table( "MTM" , "CH,TAPER" , "Taper lambdas" );
  add_var( "MTM" , "CH,TAPER" , "LAMBDA" , "Lambda" );

  add_table( "MTM" , "CH,F,SEG" , "Whole-night, per-channel per-epoch power" );
  add_var( "MTM" , "CH,F,SEG" , "MTM" , "Power" );
  set_compressed( "MTM" , tfac_t( "CH,F,SEG" ) );

  //
  // IRASA
  //

  add_cmd( "power"   , "IRASA" , "Irregular-Resampling Auto-Spectral Analysis (IRASA)" );
  add_url( "IRASA" , "power-spectra/#irasa" );
  add_verb( "IRASA" ,
            "Separate periodic oscillatory peaks from the aperiodic 1/f-like "
            "background using IRASA.\n\n"
            "IRASA repeatedly resamples the signal by non-integer factors, which "
            "shifts narrowband oscillations while leaving the fractal background more "
            "stable. By combining these spectra, Luna estimates distinct periodic and "
            "aperiodic components and their associated spectral slopes." );

  add_param( "IRASA" , "sig" , "C3,C4" , "Restrict analysis to these channels" );
  add_param( "IRASA" , "lwr" , "1" ,   "Lower frequency range" );
  add_param( "IRASA" , "upr" , "20" ,   "Upper frequency range" );
  add_param( "IRASA" , "h-min" , "1.05" , "Minimum h" );
  add_param( "IRASA" , "h-max" , "1.95" , "Maximum h" );
  add_param( "IRASA" , "h-cnt" , "17" , "Number of h steps (min-max)" );
  add_param( "IRASA" , "dB" , "" , "Decibel scale output" );
  add_param( "IRASA" , "epoch" , "" , "Report per-epoch statistics" );
  
  add_table( "IRASA" , "CH" , "Whole-night, per-channel stats" );
  add_var( "IRASA" , "CH" , "SPEC_SLOPE" , "Spectral slope" );
  add_var( "IRASA" , "CH" , "SPEC_SLOPE_N" , "Spectral slope number of points" );
  add_var( "IRASA" , "CH" , "SPEC_SLOPE_RSQ" , "Spectral slope R-sq" );

  add_table( "IRASA" , "CH,E" , "Per-epoch, per-channel stats" );
  add_var( "IRASA" , "CH,E" , "SPEC_SLOPE" , "Spectral slope" );
  add_var( "IRASA" , "CH,E" , "SPEC_SLOPE_N" , "Spectral slope number of points" );
  add_var( "IRASA" , "CH,E" , "SPEC_SLOPE_RSQ" , "Spectral slope R-sq" );

  add_table( "IRASA" , "CH,F" , "Whole-night, per-channel stats" );
  add_var( "IRASA" , "CH,F" , "APER" , "Aperiodic PSD component" );
  add_var( "IRASA" , "CH,F" , "PER" , "Periodic PSD component" );
  add_var( "IRASA" , "CH,F" , "LOGF" , "Log-transformed frequency" );

  add_table( "IRASA" , "CH,E,F" , "Epoch-level, per-channel stats" );
  add_var( "IRASA" , "CH,E,F" , "APER" , "Aperiodic PSD component" );
  add_var( "IRASA" , "CH,E,F" , "PER" , "Periodic PSD component" );
  add_var( "IRASA" , "CH,E,F" , "LOGF" , "Log-transformed frequency" );
  set_compressed( "IRASA" , tfac_t( "CH,E,F" ) );

  //
  // MSE
  //

  add_cmd( "power" , "MSE" , "Multi-scale entropy statistics" );
  add_url( "MSE" , "power-spectra/#mse" );
  add_verb( "MSE" ,
            "Estimate sample-entropy style complexity across multiple temporal "
            "scales.\n\n"
            "MSE progressively coarse-grains the signal and then computes entropy at "
            "each scale. This gives a scale-dependent view of irregularity, helping "
            "distinguish short-timescale noise-like variability from more structured "
            "longer-timescale dynamics." );

  add_param( "MSE" , "sig" , "C3,C4" , "Restrict analysis to these channels" );
  add_param( "MSE" , "m" , "3" , "Embedding dimension (default 2)" );
  add_param( "MSE" , "r" , "0.2" , "Matching tolerance in standard deviation units (default 0.15)" );
  add_param( "MSE" , "s" , "1,15,2" , "Consider scales 1 to 15, in steps of 2 (default 1 to 10 in steps of 1)" );
  add_param( "MSE" , "verbose" , "" , "Emit epoch-level MSE statistics" );
  
  add_table( "MSE" , "CH,SCALE" , "MSE per channel and scale" );
  add_var( "MSE" , "CH,SCALE" , "MSE" , "Multi-scale entropy" );

  add_table( "MSE" , "CH,E,SCALE" , "MSE per epoch, channel and scale" );
  add_var( "MSE" , "CH,E,SCALE" , "MSE" , "Multi-scale entropy" );

  //
  // LZW
  //

  add_cmd( "power" , "LZW" , "LZW compression index" );
  add_url( "LZW" , "power-spectra/#lzw" );
  add_verb( "LZW" ,
            "Estimate signal complexity using an Lempel-Ziv-Welch compression "
            "measure.\n\n"
            "The signal is coarse-grained and discretized into bins, then encoded "
            "with an LZW-style parsing scheme. More compressible sequences yield "
            "lower complexity, whereas less predictable sequences yield higher "
            "complexity." );

  add_param( "LZW" , "nsmooth" , "2" , "Coarse-graining parameter (similar to scale s in MSE)" );
  add_param( "LZW" , "nbins" , "5" , "Matching tolerance in standard deviation units (default 10)" );
  add_param( "LZW" , "epoch" , "" , "Emit epoch-level LZW statistics" );

  add_table( "LZW" , "CH" , "LZW per channel" );
  add_var( "LZW" , "CH" , "LZW" , "Compression index" );

  add_table( "LZW" , "CH,E" , "LZW per channel, per epoch" );
  add_var( "LZW" , "CH,E" , "LZW" , "Compression index" );

  //
  // HILBERT
  //

  add_cmd( "power" , "HILBERT" , "Applies filter-Hilbert transform" );
  add_url( "HILBERT" , "power-spectra/#hilbert" );
  add_verb( "HILBERT" ,
            "Construct the analytic signal in a chosen band using a band-pass filter "
            "followed by a Hilbert transform.\n\n"
            "The method first isolates the target band with an FIR filter and then "
            "computes the analytic representation, from which Luna derives the "
            "instantaneous magnitude envelope and, optionally, instantaneous phase." );

  add_param( "HILBERT" , "sig" , "C3,C4" , "Restrict analysis to these channels" );
  add_param( "HILBERT" , "f" , "0.5,4" , "Lower and upper transition frequencies" );
  add_param( "HILBERT" , "ripple" , "0.02" , "FIR filter ripple (as proportion)" );
  add_param( "HILBERT" , "tw" , "0.5" , "Transition width (in Hz)" );
  add_param( "HILBERT" , "tag" , "v1" , "Optional tag to be added to new signals" );
  add_param( "HILBERT" , "phase" , "" , "As well as magnitude, generate signal with instantaneous phase" );

  //
  // CWT
  //

  add_cmd( "power" , "CWT" , "Applies a continuous wavelet transform (convolution with a complex Morlet wavelet)" );
  add_url( "CWT" , "power-spectra/#cwt" );
  add_verb( "CWT" ,
            "Apply a complex Morlet continuous wavelet transform to the signal.\n\n"
            "CWT measures oscillatory activity by convolving the signal with a "
            "localized wavelet centered at the requested frequency. Compared with a "
            "fixed-window Fourier approach, this gives a more explicitly "
            "time-localized analytic representation and can emit both magnitude and "
            "phase channels." );

  add_param( "CWT" , "sig" , "C3,C4" , "Restrict analysis to these channels" );
  add_param( "CWT" , "fc" , "15" , "Wavelet center frequency" );
  add_param( "CWT" , "cycles" , "12" , "Bandwidth of the wavelet (number of cycles, default 7)" );
  add_param( "CWT" , "tag" , "v1" , "Additional tag to be added to the new signal" );
  add_param ( "CWT" , "phase" , "" , "Generate a second new signal with wavelet's phase" );

  //
  // CWT-DESIGN
  //

  add_cmd( "power" , "CWT-DESIGN" , "Display the properties of a complex Morlet wavelet transform" );
  add_url( "CWT-DESIGN" , "power-spectra/#cwt-design" );
  add_verb( "CWT-DESIGN" ,
            "Inspect the effective time- and frequency-domain properties of a "
            "Morlet wavelet before applying it.\n\n"
            "This command reports the wavelet coefficients and frequency response for "
            "a chosen center frequency and cycle count, helping you understand the "
            "trade-off between temporal precision and spectral precision." );

  add_param( "CWT-DESIGN" , "fs" , "200" , "Sampling rate" );
  add_param( "CWT-DESIGN" , "fc" , "15" , "Wavelet center frequency" );
  add_param( "CWT-DESIGN" , "cycles" , "7" , "Bandwidth of the wavelet (number of cycles)" );
  add_param( "CWT-DESIGN" , "fwhm" , "2" , "Alternative frequency-domain full-width at half maximum" );
  add_param( "CWT-DESIGN" , "len" , "20" , "Wavelet length in seconds for fwhm mode" );

  add_table( "CWT-DESIGN" , "PARAM" , "Wavelet summary statistics" );
  add_var( "CWT-DESIGN" , "PARAM" , "FWHM" , "Empirical time-domain full-width at half maximum [fwhm]" );
  add_var( "CWT-DESIGN" , "PARAM" , "FWHM_F" , "Frequency-domain full-width at half maximum" );
  add_var( "CWT-DESIGN" , "PARAM" , "FWHM_LWR" , "Lower half-maximum frequency" );
  add_var( "CWT-DESIGN" , "PARAM" , "FWHM_UPR" , "Upper half-maximum frequency" );

  add_table( "CWT-DESIGN" , "PARAM,F" , "Frequency response for wavelet" );
  add_var( "CWT-DESIGN" , "PARAM,F" , "MAG" , "Magnitude of response (arbitrary units)" );

  add_table( "CWT-DESIGN" , "PARAM,SEC" , "Wavelet coefficients" );
  add_var( "CWT-DESIGN" , "PARAM,SEC" , "REAL" , "Real part of wavelet" );
  add_var( "CWT-DESIGN" , "PARAM,SEC" , "IMAG" , "Imaginary part of wavelet" );

  //
  // 1FNORM
  //

  add_cmd( "power" , "1FNORM" , "Applies a differentiator filter to remove 1/f trends in signals" );
  add_url( "1FNORM" , "power-spectra/#1fnorm" );
  add_verb( "1FNORM" ,
            "Apply a differentiator-style normalization intended to attenuate broad "
            "1/f structure in the signal.\n\n"
            "This operation emphasizes faster fluctuations relative to the slower "
            "background trend, making subsequent oscillatory analyses less dominated "
            "by low-frequency power-law structure." );

  add_param( "1FNORM" , "sig" , "C3,C4" , "Restrict analysis to these channels" );

  //
  // TV
  //

  add_cmd( "power" , "TV" , "Applies of fast algorithm for 1D total variation denoising" );

  add_url( "TV" , "power-spectra/#tv" );
  add_verb( "TV" ,
            "Denoise a signal with one-dimensional total-variation regularization.\n\n"
            "TV denoising finds a piecewise-smooth approximation that balances data "
            "fidelity against the total amount of variation. It is especially useful "
            "for suppressing high-frequency noise while preserving sharper edges than "
            "a conventional linear smoother." );

  add_param( "TV" , "sig" , "C3,C4" , "Restrict analysis to these channels" );
  add_param( "TV" , "lambda" , "10" , "Smoothing parameter (0 to infinity)" );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // NREM TRANSIENTS
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // SPINDLES
  //

  add_cmd( "trans" , "SPINDLES" , "Sleep spindle detection" );
  add_url( "SPINDLES" , "spindles-so/#spindles" );
  add_verb( "SPINDLES" ,
            "Detect sleep spindles from one or more EEG channels using a "
            "wavelet-based detector centered on one or more target spindle "
            "frequencies.\n\n"
            "In the default mode, SPINDLES computes a band-limited CWT signal, "
            "smooths it in time, and applies primary and boundary thresholds to "
            "define candidate spindle cores and flanking intervals. Final events "
            "must then satisfy duration and quality-control rules, and Luna "
            "reports summary, epoch-level, and optionally per-spindle metrics "
            "such as amplitude, duration, frequency, symmetry, chirp, and "
            "integrated spindle activity.\n\n"
            "Variants include scanning a frequency grid, using an alternate "
            "wavelet width specification, estimating empirical thresholds, "
            "starting from precomputed spindle annotations, merging detections "
            "into m-spindles, caching peaks or midpoints, and adding spindle "
            "annotations or PDF/FTR exports.\n\n"
            "If the so option is given, SPINDLES first calls the SO detector on "
            "the same record, then adds spindle/SO overlap and phase-coupling "
            "analyses. Those outputs summarize whether spindles occur within "
            "detected slow oscillations, at what SO phase they occur, and how "
            "that coupling compares to permutation-based null distributions." );

  add_param( "SPINDLES" , "sig" , "C3,C4" , "Restrict analysis to these channels" );
  add_param( "SPINDLES" , "fc" , "11,15" , "Restrict analysis to these channels (otherwise, all channels are included)" );
  add_param( "SPINDLES" , "fwhm" , "1.5" , "Specify the wavelet width as FWHM instead of number of cycles" );
  add_param( "SPINDLES" , "cycles" , "12" , "Number of cycles (default 7)" );
  add_param( "SPINDLES" , "precomputed" , "SP" , "Use precomputed spindle annotations instead of detecting new spindles" );
  add_param( "SPINDLES" , "th" , "6" , "Multiplicative threshold for core spindle detection (default 4.5)" );
  add_param( "SPINDLES" , "th2" , "3" , "Multiplicative threshold for non-core spindle detection (default=2)" );
  add_param( "SPINDLES" , "median" , "" , "Flag to indicate that the median, not mean, is used for thresholding" );
  add_param( "SPINDLES" , "winsor" , "4" , "Winsorize extreme coefficients before thresholding" );
  add_param( "SPINDLES" , "zpks" , "30" , "Use smoothed z-score peaks rather than averaged CWT as the detection statistic" );
  add_param( "SPINDLES" , "influence" , "0.01" , "Influence parameter for z-score peak detection" );
  add_param( "SPINDLES" , "q" , "0.3" , "Quality metric criterion for individual spindles (default 0)" );

  add_param( "SPINDLES" , "fc-lower" , "9" , "Lower limit if iterating over multiple F_C values" );
  add_param( "SPINDLES" , "fc-upper" , "16" , "Upper limit if iterating over multiple F_C values" );
  add_param( "SPINDLES" , "fc-step" , "2" , "Increment step if iterating over multiple F_C values" );
  add_param( "SPINDLES" , "th-max" , "10" , "Maximum threshold for spindle core (default: none)" );
  add_param( "SPINDLES" , "min" , "1" , "Minimum duration for an entire spindle (default 0.5 seconds)" );
  add_param( "SPINDLES" , "min0" , "0.3" , "Minimum duration for a spindle core (default 0.3 seconds)" );
  add_param( "SPINDLES" , "max" , "2" , "Maximum duration for an entire spindle (default 3 seconds)" );
  add_param( "SPINDLES" , "win" , "0.2" , "Smoothing window for wavelet coefficients (default 0.1 seconds)" );
  add_param( "SPINDLES" , "local" , "120" , "Use local window (in seconds) to define baseline for spindle detection" );

  add_param( "SPINDLES" , "epoch" , "" , "Show epoch-level counts" );
  add_param( "SPINDLES" , "per-spindle" , "" , "Show per-spindle output" );
  add_param( "SPINDLES" , "dynam" , "" , "Add per-epoch summaries suitable for ultradian dynamics analyses" );
  
  add_param( "SPINDLES" , "empirical" , "" , "Empirically determine thresholds" );
  hidden_param( "SPINDLES" , "set-empirical" , "" , "Use empirically determined thresholds for spindle detection" );
  hidden_param( "SPINDLES" , "verbose-empirical" , "" , "Output extensive information on threshold estimation" );

  add_param( "SPINDLES" , "merge" , "0.2", "Merge two putative spindles if within this interval (default 0.5 seconds)" );
  add_param( "SPINDLES" , "collate" , "" , "Within each channel, collate overlapping spindles of similar frequencies" );
  add_param( "SPINDLES" , "collate-channels" , "" , "As above, except merge across channels also" );
  add_param( "SPINDLES" , "th-frq" , "1" , "Frequency criterion for merging spindles (default 2 Hz)" );
  add_param( "SPINDLES" , "list-all-spindles" , "" , "List all spindles that comprise each m-spindle" );

  add_param( "SPINDLES" , "th-interval" , "0.5" , "Merge if the ratio of intersection to union is at least this (default 0, i.e. any overlap)" );
  hidden_param( "SPINDLES" , "th-interval-cross-channel" , "" , "not currently used" );
  hidden_param( "SPINDLES" , "th-interval-within-channel" , "" , "not currently used" );
  add_param( "SPINDLES" , "window" , "0.5" , "Set window around each spindle when defining temporal overlap" );
  add_param( "SPINDLES" , "hms" , "" , "Show clock-time of each m-spindle" );

  add_param( "SPINDLES" , "ftr" , "tag" , "Produce FTR files for all spindles, with the tag in the filename" );
  add_param( "SPINDLES" , "ftr-dir" , "/path/to/folder" , "Folder for FTR files" );
  add_param( "SPINDLES" , "annot" , "SP" , "Write detected spindles as annotations, optionally with this label" );
  add_param( "SPINDLES" , "out" , "spindles.annot" , "Write spindle annotations to an output file" );
  add_param( "SPINDLES" , "cache" , "spindles" , "Cache spindle-level summaries" );
  add_param( "SPINDLES" , "cache-peaks" , "spindles-peaks" , "Cache spindle peak sample-points" );
  add_param( "SPINDLES" , "cache-mids" , "spindles-mids" , "Cache spindle midpoints" );
  add_param( "SPINDLES" , "cache-peaks-sec" , "spindles-sec" , "Cache spindle peaks in elapsed seconds" );
  add_param( "SPINDLES" , "prop" , "" , "Estimate spindle propagation across channels" );
  add_param( "SPINDLES" , "collate-within-channel" , "" , "Merge overlapping spindle detections within channel only" );
  add_param( "SPINDLES" , "collate-annot" , "MSP" , "Write merged spindle events as annotations" );
  add_param( "SPINDLES" , "pdf" , "tag" , "Write per-channel spindle PDF summaries" );
  hidden_param( "SPINDLES" , "show-coef" , "" , "Request (very verbose) coefficient output (to stdout)" );

  // output

  add_table( "SPINDLES" , "CH,F" , "Individual-level output" );
  add_var( "SPINDLES" , "CH,F" , "DENS" , "Spindle density (count per minute)" );
  add_var( "SPINDLES" , "CH,F" , "CDENS" , "Coupled spindle density (count per minute)" );
  add_var( "SPINDLES" , "CH,F" , "UDENS" , "Uncoupled spindle density (count per minute)" );
  add_var( "SPINDLES" , "CH,F" , "AMP" , "Mean spindle amplitude (uV or mV units)" );
  add_var( "SPINDLES" , "CH,F" , "ACT_MX" , "Mean max spindle activity (normed CWT)" );
  add_var( "SPINDLES" , "CH,F" , "ACT_MN" , "Mean average spindle activity (normed CWT)" );
  add_var( "SPINDLES" , "CH,F" , "DUR" , "Mean spindle duration (core+flanking region)" );
  add_var( "SPINDLES" , "CH,F" , "NOSC" , "Mean number of oscillations per spindle" );
  add_var( "SPINDLES" , "CH,F" , "FWHM" , "Mean spindle FWHM (full width at half maximum)" );
  add_var( "SPINDLES" , "CH,F" , "ISA_S" , "Mean integrated spindle activity (ISA) per spindle" );
  add_var( "SPINDLES" , "CH,F" , "ISA_M" , "Mean integrated spindle activity (ISA) per minute" );
  add_var( "SPINDLES" , "CH,F" , "ISA_T" , "Total integrated spindle activity (ISA)" );
  add_var( "SPINDLES" , "CH,F" , "FRQ" , "Mean spindle frequency (from counting zero-crossings)" );
  add_var( "SPINDLES" , "CH,F" , "FFT" , "Mean spindle frequency (from FFT)" );
  add_var( "SPINDLES" , "CH,F" , "CHIRP" , "Mean chirp metric per spindle" );
  add_var( "SPINDLES" , "CH,F" , "SYMM" , "Mean spindle symmetry metric" );
  add_var( "SPINDLES" , "CH,F" , "SYMM2" , "Mean spindle folded-symmetry metric" );
  add_var( "SPINDLES" , "CH,F" , "Q" , "Mean spindle quality metric" );
  add_var( "SPINDLES" , "CH,F" , "DISPERSION" , "Mean dispersion index of epoch spindle count" );
  add_var( "SPINDLES" , "CH,F" , "DISPERSION_P" , "P-value for test of over-dispersion" );
  add_var( "SPINDLES" , "CH,F" , "MINS" , "Total duration of signal entered into the analysis (minutes)" );
  add_var( "SPINDLES" , "CH,F" , "NE" , "Number of epochs" );
  add_var( "SPINDLES" , "CH,F" , "N01" , "Number of spindles prior to merging" );
  add_var( "SPINDLES" , "CH,F" , "N02" , "Number of spindles post merging, prior to QC" );

  add_var( "SPINDLES" , "CH,F" , "EMPTH" , "Empirically-determined threshold" );
  add_var( "SPINDLES" , "CH,F" , "EMPF" , "Relative frequency of above-thresholds points based on EMPTH" );
  add_var( "SPINDLES" , "CH,F" , "MEAN_OVER_MEDIAN" , "Ratio of mean to median, to index skewness of the wavelet coefficients" );

  add_var( "SPINDLES" , "CH,F" , "CWT_TH" , "CWT threshold" );
  add_var( "SPINDLES" , "CH,F" , "FRNG2" , "Range of spindle frequencies" );
  add_var( "SPINDLES" , "CH,F" , "FRQ1" , "Frequency in spindle first half" );
  add_var( "SPINDLES" , "CH,F" , "FRQ2" , "Frequency in spindle second half" );
  add_var( "SPINDLES" , "CH,F" , "FVAR2" , "Variation in spindle frequency" );
  add_var( "SPINDLES" , "CH,F" , "N" , "Number of spindles" );
  add_var( "SPINDLES" , "CH,F" , "P01" , "Pre-QC" );
  add_var( "SPINDLES" , "CH,F" , "P02" , "Mid-QC" );
  add_var( "SPINDLES" , "CH,F" , "SEC_AMP" , "Midpoint based on CWT" );
  add_var( "SPINDLES" , "CH,F" , "SEC_P2P" , "Midpoint based on peak-to-peak" );
  add_var( "SPINDLES" , "CH,F" , "SEC_TROUGH" , "Midpoint based on trough" );
  add_var( "SPINDLES" , "CH,F" , "SYMM_AMP" , "Mean spindle symmetry metric (based on CWT)" );
  add_var( "SPINDLES" , "CH,F" , "SYMM_TROUGH" , "Mean spindle symmetry metric (based on trough)" );
  
  
  add_table( "SPINDLES" , "CH,F,TH" , "Between-class variance over range of thresholds" );
  add_var( "SPINDLES" , "CH,F,TH" , "SIGMAB" , "Between-class variance for given threshold" );

  add_table( "SPINDLES" , "CH,E,F" , "Epoch-level output [epoch]" ); 
  add_var( "SPINDLES" , "CH,E,F" , "N" , "Number of spindles observed in that epoch (for that target frequency/channel)" );

  add_table( "SPINDLES" , "CH,F,SPINDLE" , "Spindle-level output [per-spindle]" ); 
  add_var( "SPINDLES" , "CH,F,SPINDLE" , "AMP" , "Spindle amplitude (uV or mV units)" );
  add_var( "SPINDLES" , "CH,F,SPINDLE" , "ACT_MX" , "Max spindle activity (normed CWT)" );
  add_var( "SPINDLES" , "CH,F,SPINDLE" , "ACT_MN" , "Average spindle activity (normed CWT)" );
  add_var( "SPINDLES" , "CH,F,SPINDLE" , "CHIRP" , "Spindle chirp (-1 to +1)" );
  add_var( "SPINDLES" , "CH,F,SPINDLE" , "DUR" , "Spindle duration (seconds)" );
  add_var( "SPINDLES" , "CH,F,SPINDLE" , "FWHM" , "Spindle FWHM (seconds)" );
  add_var( "SPINDLES" , "CH,F,SPINDLE" , "NOSC" , "Number of oscillations" );
  add_var( "SPINDLES" , "CH,F,SPINDLE" , "FRQ" , "Spindle frequency based on counting zero-crossings in bandpass filtered signal" );
  add_var( "SPINDLES" , "CH,F,SPINDLE" , "FFT" , "Spindle frequency based on FFT" );
  add_var( "SPINDLES" , "CH,F,SPINDLE" , "ISA" , "Integrated spindle activity" );
  add_var( "SPINDLES" , "CH,F,SPINDLE" , "MAXSTAT" , "Maximum wavelet statistic" );
  add_var( "SPINDLES" , "CH,F,SPINDLE" , "MEANSTAT" , "Mean wavelet statistic" );
  add_var( "SPINDLES" , "CH,F,SPINDLE" , "Q" , "Quality metric" );
  add_var( "SPINDLES" , "CH,F,SPINDLE" , "PASS" , "Flag (0/1) for whether this spindle passes the quality metric criterion" );
  add_var( "SPINDLES" , "CH,F,SPINDLE" , "START" , "Start position of the spindle (seconds elapsed since start of EDF)" );
  add_var( "SPINDLES" , "CH,F,SPINDLE" , "PEAK" , "Peak/mid position of the spindle (seconds elapsed since start of EDF)" );
  add_var( "SPINDLES" , "CH,F,SPINDLE" , "STOP" , "Stop position of the spindle (seconds elapsed since start of EDF)" );
  add_var( "SPINDLES" , "CH,F,SPINDLE" , "START_SP" , "Start position of the spindle (in sample-units relative to current in-memory EDF)" );
  add_var( "SPINDLES" , "CH,F,SPINDLE" , "PEAK_SP" , "Peak/mid position of the spindle (in sample-units relative to the current in-memory EDF)" );
  add_var( "SPINDLES" , "CH,F,SPINDLE" , "STOP_SP" , "Stop position of the spindle (in sample-units relative to the current in-memory EDF)" );
  add_var( "SPINDLES" , "CH,F,SPINDLE" , "SYMM" , "Symmetry index (relative position of peak)" );
  add_var( "SPINDLES" , "CH,F,SPINDLE" , "SYMM2" , "Folded symmetry index (0=symmetrical, 1=asymmetrical)" );
  hidden_var( "SPINDLES" , "CH,F,SPINDLE" , "IF" , "Mean frequency per spindle over duration [if]" );

  add_table( "SPINDLES" , "CH,F,B,SPINDLE" , "Band enrichment (per-spindle)" );
  add_var( "SPINDLES" , "CH,F,B,SPINDLE" , "ENRICH" , "Spindle enrichment" );
  
  hidden_table( "SPINDLES" , "CH,F,RELLOC" , "Mean IF stratified by relative location in spindle [if]" );
  hidden_var( "SPINDLES" , "CH,F,RELLOC" , "IF" , "Mean frequency of all spindles, per relative position within the spindle (five bins)" );
  
  hidden_table( "SPINDLES", "F,CH,PHASE,RELLOC" , "Mean IF stratified by phase and relative location in spindle [if]" );
  hidden_var( "SPINDLES" , "F,CH,PHASE,RELLOC" , "SOPL_CHIRP" , "Spindle chirp" );
  
  add_table( "SPINDLES" , "" , "Individual-level summaries of m-spindles [collate]" );
  add_var( "SPINDLES" , "" , "MSP_DENS" , "m-spindle density" );
  add_var( "SPINDLES" , "" , "MSP_N" , "m-spindle count" );
  add_var( "SPINDLES" , "" , "MSP_MINS" , "Denominator for density, i.e. minutes of signal analyzed" );

  add_table( "SPINDLES" , "F" , "m-spindle density stratified by m-spindle frequency [collate]" );
  add_var( "SPINDLES" , "F" , "MSP_DENS" , "m-spindle density conditional on m-spindle frequency" );
  

  add_table( "SPINDLES" , "F,SEED" , "Spindle propagation seed summaries" );
  add_var( "SPINDLES" , "F,SEED" , "R" , "Relative SEED position among overlapping CHs" );
  add_var( "SPINDLES" , "F,SEED" , "T" , "Relative SEED time among overlapping CHs" );

  add_table( "SPINDLES" , "F,CH,SEED" , "Spindle propagation seed-channel stats" );
  add_var( "SPINDLES" , "F,CH,SEED" , "A" , "Channel amplitude relative to SEED" );
  add_var( "SPINDLES" , "F,CH,SEED" , "A_PRESEED" , "Channel amplitude relative to SEED (CH<SEED)" );
  add_var( "SPINDLES" , "F,CH,SEED" , "A_POSTSEED" , "Channel amplitude relative to SEED (SEED<CH)" );
  add_var( "SPINDLES" , "F,CH,SEED" , "N" , "Count above threshold CH-peaks" );
  add_var( "SPINDLES" , "F,CH,SEED" , "N_PRESEED" , "Count above threshold pre-SEED CH-peaks" );
  add_var( "SPINDLES" , "F,CH,SEED" , "N_POSTSEED" , "Count above threshold post-SEED CH-peaks" );
  add_var( "SPINDLES" , "F,CH,SEED" , "T" , "CH-peak time relative to SEED" );
  add_var( "SPINDLES" , "F,CH,SEED" , "T_PRESEED" , "CH-peak time relative to SEED (CH<SEED)" );
  add_var( "SPINDLES" , "F,CH,SEED" , "T_POSTSEED" , "CH-peak time relative to SEED (SEED<CH)" );
  add_var( "SPINDLES" , "F,CH,SEED" , "P" , "Proportion above threshold CH-peaks" );
  add_var( "SPINDLES" , "F,CH,SEED" , "P_PRESEED" , "Proportion above threshold pre-SEED CH-peaks" );
  add_var( "SPINDLES" , "F,CH,SEED" , "P_POSTSEED" , "Proportion above threshold post-SEED CH-peaks" );
  add_var( "SPINDLES" , "F,CH,SEED" , "PP" , "CH-SEED pre/post metric" );

  
  add_table( "SPINDLES" , "SPINDLE,F,CH,SEED" , "Spindle propagation seed-channel stats" );
  add_var( "SPINDLES" , "SPINDLE,F,CH,SEED" , "REL" , "Relative position" );
  add_var( "SPINDLES" , "SPINDLE,F,CH,SEED" , "T" , "Time" );


  add_table( "SPINDLES" , "MSPINDLE" , "Merged-spindle output [collate]" );
  add_var( "SPINDLES" , "MSPINDLE" , "MSP_DUR","Duration of this m-spindle" );
  add_var( "SPINDLES" , "MSPINDLE" , "MSP_F","Estimated frequency of this m-spindle" );
  add_var( "SPINDLES" , "MSPINDLE" , "MSP_FL","Lower frequency of this m-spindle" );
  add_var( "SPINDLES" , "MSPINDLE" , "MSP_FU","Upper frequency of this m-spindle" );
  add_var( "SPINDLES" , "MSPINDLE" , "MSP_SIZE","Number of spindles in this m-spindle" );
  add_var( "SPINDLES" , "MSPINDLE" , "MSP_STAT","Statistic for m-spindle" );
  add_var( "SPINDLES" , "MSPINDLE" , "MSP_START","Start time (seconds elapsed from EDF start) of m-spindle" );
  add_var( "SPINDLES" , "MSPINDLE" , "MSP_STOP","Stop time (seconds elapsed from EDF start) of m-spindle" );
  add_var( "SPINDLES" , "MSPINDLE" , "MSP_START_HMS" , "Merged spindle start clock-time (if 'hms')" );
  add_var( "SPINDLES" , "MSPINDLE" , "MSP_STOP_HMS" , "Merged spindle stop clock-time (if 'hms')" );
  
  add_table( "SPINDLES" , "CH,MSPINDLE" , "Within-channel merged-spindle output [collate]" );
  add_var( "SPINDLES" , "CH,MSPINDLE" , "MSP_DUR","Duration of this m-spindle" );
  add_var( "SPINDLES" , "CH,MSPINDLE" , "MSP_F","Estimated frequency of this m-spindle" );
  add_var( "SPINDLES" , "CH,MSPINDLE" , "MSP_FL","Lower frequency of this m-spindle" );
  add_var( "SPINDLES" , "CH,MSPINDLE" , "MSP_FU","Upper frequency of this m-spindle" );
  add_var( "SPINDLES" , "CH,MSPINDLE" , "MSP_SIZE","Number of spindles in this m-spindle" );
  add_var( "SPINDLES" , "CH,MSPINDLE" , "MSP_STAT","Statistic for m-spindle" );
  add_var( "SPINDLES" , "CH,MSPINDLE" , "MSP_START","Start time (seconds elapsed from EDF start) of m-spindle" );
  add_var( "SPINDLES" , "CH,MSPINDLE" , "MSP_STOP","Stop time (seconds elapsed from EDF start) of m-spindle" );
  add_var( "SPINDLES" , "CH,MSPINDLE" , "MSP_START_HMS" , "Merged spindle start clock-time (if 'hms')" );
  add_var( "SPINDLES" , "CH,MSPINDLE" , "MSP_STOP_HMS" , "Merged spindle stop clock-time (if 'hms')" );
 
  
  add_table( "SPINDLES" , "SPINDLE,MSPINDLE" , "Spindle to m-spindle mappings (from 'list-all-spindles') [collate]" );
  add_var( "SPINDLES" , "SPINDLE,MSPINDLE" , "SCH", "Spindle label (channel:target frequency)" );
  add_var( "SPINDLES" , "SPINDLE,MSPINDLE" , "FFT", "Spindle estimated frequency (via FFT)" );
  add_var( "SPINDLES" , "SPINDLE,MSPINDLE" , "START" , "Spindle start time (elapsed seconds from EDF start)" );
  add_var( "SPINDLES" , "SPINDLE,MSPINDLE" , "STOP" , "Spindle stop time (elapsed seconds from EDF start)" );
  
  // experimental
  hidden_param( "SPINDLES" , "if" , "" , "Estimate instantaneous frequency of spindles" );
  hidden_param( "SPINDLES" , "if-frq" , "1" , "Window around target frequency (default 2 hz)" );
  hidden_param( "SPINDLES" , "tlock" , "" , "Flag to request (verbose) average, peak-locked waveforms" );
  hidden_param( "SPINDLES" , "verbose-coupling" , "" , "Add extra tables of EEG/CWT phase/time-locked to SO" );


  // show-coef verbose output

  add_table( "SPINDLES" , "F,CH,T" , "Verbose threshold/coefficient output [show-coeff]" );
  add_var( "SPINDLES" , "F,CH,T" , "SEC" , "Time (sec)" );
  add_var( "SPINDLES" , "F,CH,T" , "RAWCWT" , "Raw CWT coefficient" );
  add_var( "SPINDLES" , "F,CH,T" ,"CWT" , "CWT coefficient" );
  add_var( "SPINDLES" , "F,CH,T" ,"AVG" , "Averaged CWT coefficient" );
  add_var( "SPINDLES" , "F,CH,T" ,"AVG_CORR" , "Averaged baseline-corrected CWT coefficient" );
  add_var( "SPINDLES" , "F,CH,T" ,"CWT_TH" , "CWT primary threshold" );
  add_var( "SPINDLES" , "F,CH,T" ,"CWT_TH2" , "CWT secondary threshold" );
  add_var( "SPINDLES" , "F,CH,T" ,"CWT_THMAX" , "CWT maximum threshold" );
  add_var( "SPINDLES" , "F,CH,T" ,"PUTATIVE" , "Pre-QC spindle" );
  add_var( "SPINDLES" , "F,CH,T" ,"SPINDLE" , "Post-QC spindle" );
  

  //
  // SO (duplicated from SO command below)
  //


  add_param( "SPINDLES" , "so" , "C3" , "Detect SOs, optionally on this channel, and evaluate spindle/SO coupling" );
  add_param( "SPINDLES" , "phase" , "C3" , "Use the phase of another channel/band for coupling analyses without discrete SO detection" );
  
  add_param( "SPINDLES" , "mag" , "2" , "SO, relative mangitude threshold (times mean/median)" );
  add_param( "SPINDLES" , "uV-neg" , "-40" , "SO, absolute negative peak uV amplitude threshold" );
  add_param( "SPINDLES" , "uV-p2p" , "80" , "SO, absolute peak-to-peak uV amplitude threshold" );
    
  add_param( "SPINDLES" , "f-lwr" , "0.2" , "SO filter, lower transition frequency" );
  add_param( "SPINDLES" , "f-upr" , "4.5" , "SO filter, upper transition frequency" );
  
  add_param( "SPINDLES" , "t-lwr" , "0" , "SO, lower duration (secs)" );
  add_param( "SPINDLES" , "t-upr" , "3" , "SO, upper duration (secs)" );
  
  add_param( "SPINDLES" , "t-neg-lwr" , "0" , "SO, lower duration for negative peak (secs)" );
  add_param( "SPINDLES" , "t-neg-upr" , "1" , "SO, upper duration for negative peak (secs)" );
  
  hidden_param( "SPINDLES" , "neg2pos" , "" , "SO, Use negative-to-positive zero crossings" );
  add_param( "SPINDLES" , "th-mean" , "" , "SO, use mean not median" );
  add_param( "SPINDLES" , "stats-median" , "" , "SO, use median (not mean) when reporting stats over SOs" );  
 
  add_table( "SPINDLES" , "CH" , "SO channel-level statistics" );
  add_var( "SPINDLES" , "CH" , "SO" , "Number of SO detected" );
  add_var( "SPINDLES" , "CH" , "SO_RATE" , "SO per minute" );
  add_var( "SPINDLES" , "CH" , "SO_AMP_NEG" , "SO amplitude (negative peak)" );
  add_var( "SPINDLES" , "CH" , "SO_AMP_POS" , "SO amplitude (positive peak)" );
  add_var( "SPINDLES" , "CH" , "SO_AMP_P2P" , "SO peak-to-peak amplitude" );
  add_var( "SPINDLES" , "CH" , "SO_DUR" , "SO duration (secs)" );
  add_var( "SPINDLES" , "CH" , "SO_DUR_NEG" , "Negative peak SO duration (secs)" );
  add_var( "SPINDLES" , "CH" , "SO_DUR_POS", "Positive peak SO duration (secs)" );
  add_var( "SPINDLES" , "CH" , "SO_TRANS" , "SO transition (secs)" );
  add_var( "SPINDLES" , "CH" , "SO_TRANS_FREQ" , "SO transition freq (Hz)" );
  add_var( "SPINDLES" , "CH" , "SO_SLOPE_POS1" , "Positive peak rising slope" );
  add_var( "SPINDLES" , "CH" , "SO_SLOPE_POS2" , "Positive peak falling slope" );
  add_var( "SPINDLES" , "CH" , "SO_SLOPE_NEG1" , "Negative peak falling slope" );
  add_var( "SPINDLES" , "CH" , "SO_SLOPE_NEG2" , "Negative peak rising slope" );
  add_var( "SPINDLES" , "CH" , "SO_SLOPE" , "Negative peak rising slope" );
  add_var( "SPINDLES" , "CH" , "SO_TH_NEG" , "Negative peak threshold [mag]" );
  add_var( "SPINDLES" , "CH" , "SO_TH_P2P" , "Peak-to-peak threshold [mag]" );

  add_table( "SPINDLES" , "CH,E" , "Epoch-level SO statistics" );
  add_var( "SPINDLES" , "CH,E" , "N" , "Number of SO detected" );
  add_var( "SPINDLES" , "CH,E" , "DOWN_AMP" , "Number of SO detected" );
  add_var( "SPINDLES" , "CH,E" , "UP_AMP" , "Number of SO detected" );
  add_var( "SPINDLES" , "CH,E" , "P2P_AMP" , "Number of SO detected" );
  add_var( "SPINDLES" , "CH,E" , "SLOPE_NEG1" , "Negative peak falling slope" );
  add_var( "SPINDLES" , "CH,E" , "SLOPE_NEG2" , "Negative peak rising slope" );
  add_var( "SPINDLES" , "CH,E" , "SLOPE_POS1" , "Positive peak rising slope" );
  add_var( "SPINDLES" , "CH,E" , "SLOPE_POS2" , "Positive peak falling slope" );
  
  add_table( "SPINDLES" , "CH,N" , "per-SO statistics" );
  add_var( "SPINDLES" , "CH,N" , "DOWN_AMP" , "Negative peak SO amplitude" );
  add_var( "SPINDLES" , "CH,N" , "DOWN_IDX" , "Negative peak sample index" ); 
  add_var( "SPINDLES" , "CH,N" , "UP_AMP" , "Positive peak SO ampltiude" ); 
  add_var( "SPINDLES" , "CH,N" , "UP_IDX" , "Positive peak sample index" ); 
  add_var( "SPINDLES" , "CH,N" , "START" , "Start of SO (in seconds elapsed from start of EDF)" ); 
  add_var( "SPINDLES" , "CH,N" , "START_IDX" , "Start of SO (in sample-point units)" ); 
  add_var( "SPINDLES" , "CH,N" , "STOP" , "Stop of SO (in seconds elapsed from start of EDF)" ); 
  add_var( "SPINDLES" , "CH,N" , "STOP_IDX" , "Stop of SO (in sample-point units)" ); 
  add_var( "SPINDLES" , "CH,N" , "DUR" , "SO duration (sec)" ); 
  add_var( "SPINDLES" , "CH,N" , "DUR1" , "SO HW1 duration (sec)" ); 
  add_var( "SPINDLES" , "CH,N" , "DUR2" , "SO HW2 duration (sec)" ); 
  add_var( "SPINDLES" , "CH,N" , "TRANS" , "SO transition (sec)" );
 
  add_var( "SPINDLES" , "CH,N" , "P2P_AMP" , "SO peak-to-peak amplitude" ); 
  add_var( "SPINDLES" , "CH,N" , "SLOPE_POS1" , "Positive peak rising slope" ); 
  add_var( "SPINDLES" , "CH,N" , "SLOPE_POS2" , "Positive peak falling slope" ); 
  add_var( "SPINDLES" , "CH,N" , "SLOPE_NEG1" , "Negative peak falling slope" );
  add_var( "SPINDLES" , "CH,N" , "SLOPE_NEG2" , "Negative peak rising slope" );
  
  //
  // SP/SO coupling (w/out ANCHOR)
  //

  add_table( "SPINDLES" , "CH,F" , "SP/SO coupling stats" );
  add_var( "SPINDLES" , "CH,F" , "COUPL_MAG" , "SO/SP coupling: magnitude (original statistic)" );
  add_var( "SPINDLES" , "CH,F" , "COUPL_MAG_NULL" , "SO/SP coupling: meanmagnitude under null" );
  add_var( "SPINDLES" , "CH,F" , "COUPL_MAG_Z" , "SO/SP coupling: magnitude (empirical Z)" );
  add_var( "SPINDLES" , "CH,F" , "COUPL_MAG_EMP" , "SO/SP coupling: magnitude (empirical P)" );  
  add_var( "SPINDLES" , "CH,F" , "COUPL_OVERLAP" , "SO/SP coupling: overlap (original statistic)" );
  add_var( "SPINDLES" , "CH,F" , "COUPL_OVERLAP_NULL" , "SO/SP coupling: mean overlap under null" );
  add_var( "SPINDLES" , "CH,F" , "COUPL_OVERLAP_Z" , "SO/SP coupling: overlap (empirical Z)" );
  add_var( "SPINDLES" , "CH,F" , "COUPL_OVERLAP_EMP" , "SO/SP coupling: overlap (empirical P)" );  
  add_var( "SPINDLES" , "CH,F" , "COUPL_ANGLE" , "SO/SP coupling: mean SO phase angle at spindle peak" );
  add_var( "SPINDLES" , "CH,F" , "COUPL_PV" , "SO/SP coupling: asymptotic ITPC p-value" );
  add_var( "SPINDLES" , "CH,F" , "COUPL_SIGPV_NULL" , "SO/SP coupling: null rate of asymptotic ITPC p-value < 0.05" );


  //
  // SP/SO coupling options
  //
  
  add_param( "SPINDLES" , "nreps" , "1000" , "SO/SP coupling: number of replications for SP/SO coupling" );
  add_param( "SPINDLES" , "perm-whole-trace" , "" , "SO/SP coupling: Do not use within-epoch shuffling" );
  add_param( "SPINDLES" , "all-spindles" , "" , "SO/SP coupling: Sonsider all spindles, whether ot not they overlap a SO" );
  add_param( "SPINDLES" , "stratify-by-phase" , "" , "SO/SP coupling: Overlap statistics per SO phase bin" );

  add_table( "SPINDLES" , "ANCHOR,CH,F" , "SP/SO coupling stats" );
  add_var( "SPINDLES" , "ANCHOR,CH,F" , "COUPL_MAG" , "SO/SP coupling: magnitude (original statistic)" );
  add_var( "SPINDLES" , "ANCHOR,CH,F" , "COUPL_MAG_NULL" , "SO/SP coupling: meanmagnitude under null" );
  add_var( "SPINDLES" , "ANCHOR,CH,F" , "COUPL_MAG_Z" , "SO/SP coupling: magnitude (empirical Z)" );
  add_var( "SPINDLES" , "ANCHOR,CH,F" , "COUPL_MAG_EMP" , "SO/SP coupling: magnitude (empirical P)" );  
  add_var( "SPINDLES" , "ANCHOR,CH,F" , "COUPL_OVERLAP" , "SO/SP coupling: overlap (original statistic)" );
  add_var( "SPINDLES" , "ANCHOR,CH,F" , "COUPL_OVERLAP_NULL" , "SO/SP coupling: mean overlap under null" );
  add_var( "SPINDLES" , "ANCHOR,CH,F" , "COUPL_OVERLAP_Z" , "SO/SP coupling: overlap (empirical Z)" );
  add_var( "SPINDLES" , "ANCHOR,CH,F" , "COUPL_OVERLAP_EMP" , "SO/SP coupling: overlap (empirical P)" );  
  add_var( "SPINDLES" , "ANCHOR,CH,F" , "COUPL_ANGLE" , "SO/SP coupling: mean SO phase angle at spindle peak" );
  add_var( "SPINDLES" , "ANCHOR,CH,F" , "COUPL_PV" , "SO/SP coupling: asymptotic ITPC p-value" );
  add_var( "SPINDLES" , "ANCHOR,CH,F" , "COUPL_SIGPV_NULL" , "SO/SP coupling: null rate of asymptotic ITPC p-value < 0.05" );

  add_table( "SPINDLES" , "CH,F,PHASE" , "SO-phase stratified spindle overlap" );
  add_var( "SPINDLES" , "CH,F,PHASE" , "COUPL_OVERLAP" , "SO/SP coupling: overlap (original statistic)" );
  add_var( "SPINDLES" , "CH,F,PHASE" , "COUPL_OVERLAP_EMP" , "SO/SP coupling: overlap (empirical P)" );
  add_var( "SPINDLES" , "CH,F,PHASE" , "COUPL_OVERLAP_Z" , "SO/SP coupling: overlap (Z statistic)" );
  add_var( "SPINDLES" , "CH,F,PHASE" , "SOPL_CHIRP" , "Spindle frequency | SO phase" );
   
  // spindle-level SO-coupling output
  add_var( "SPINDLES" , "CH,F,SPINDLE" , "PEAK" , "Spindle peak (seconds)" );
  add_var( "SPINDLES" , "CH,F,SPINDLE" , "SO_NEAREST" , "SO/SP coupling: time to nearest SO (0 if in one)" );
  add_var( "SPINDLES" , "CH,F,SPINDLE" , "SO_NEAREST_NUM" , "SO/SP coupling: number of nearest SO" );
  add_var( "SPINDLES" , "CH,F,SPINDLE" , "SO_PHASE_PEAK" , "SO/SP coupling: SO phase at spindle peak, if in SO" );


  add_table( "SPINDLES", "CH,PHASE" , "Raw EEG by SO phase" );
  add_var( "SPINDLES" , "CH,PHASE", "SOPL_EEG" , "Average EEG" );
  
  add_table( "SPINDLES", "CH,SP" , "Raw EEG by time from SO negative peak" );
  add_var( "SPINDLES" , "CH,SP", "SOTL_EEG" , "Average EEG" );

  add_table( "SPINDLES", "CH,F,PHASE" , "Spindle CWT by SO phase" );
  add_var( "SPINDLES" , "CH,F,PHASE", "SOPL_CWT" , "Spindle CWT" );

  add_table( "SPINDLES", "CH,F,SP" , "Spindle CWT by time from SO negative peak" );
  add_var( "SPINDLES" , "CH,F,SP", "SOTL_CWT" , "Spindle CWT" );

  //
  // SO
  //

  add_cmd( "trans" , "SO" , "Slow oscillation detection" );
  add_url( "SO" , "spindles-so/#so" );
  add_verb( "SO" ,
            "Detect slow oscillations and related half-wave events from one or "
            "more EEG channels after low-frequency filtering.\n\n"
            "SO identifies zero-crossing-defined half-waves, measures their "
            "durations, amplitudes, slopes, and transition frequencies, and "
            "retains events that satisfy the requested amplitude and timing "
            "criteria. The default use is to detect full slow oscillations, but "
            "the command can also operate in half-wave mode, positive-only or "
            "negative-only mode, and optional SO-versus-delta classification "
            "based on percentile and duration rules.\n\n"
            "Outputs include channel-level macro summaries, epoch-level counts "
            "and means, and optionally per-event details such as start/stop, "
            "peak indices, amplitudes, durations, and slopes. Additional modes "
            "can write annotations, cache positive or negative peaks, emit "
            "per-event dynamics-friendly output, and perform time-locked "
            "averaging of another signal around SO onset or peak." );
  
  add_param( "SO" , "sig" , "C3,C4" , "Restrict analysis to these channels" );
  
  add_param( "SO" , "mag" , "2" , "Relative mangitude threshold (times mean/median)" );
  add_param( "SO" , "uV-neg" , "-40" , "Absolute negative peak uV amplitude threshold" );
  add_param( "SO" , "uV-p2p" , "80" , "Absolute peak-to-peak uV amplitude threshold" );
  add_param( "SO" , "f-lwr" , "0.2" , "Lower transition frequency" );
  add_param( "SO" , "f-upr" , "4.5" , "Upper transition frequency" );
  add_param( "SO" , "t-lwr" , "0" , "Lower duration (secs)" );
  add_param( "SO" , "t-upr" , "3" , "Upper duration (secs)" );
  add_param( "SO" , "t-neg-lwr" , "0" , "Lower duration for negative peak (secs)" );
  add_param( "SO" , "t-neg-upr" , "1" , "Upper duration for negative peak (secs)" );
  add_param( "SO" , "t-pos-lwr" , "0" , "Lower duration for positive peak (secs)" );
  add_param( "SO" , "t-pos-upr" , "1" , "Upper duration for positive peak (secs)" );
  add_param( "SO" , "t-p2p-min" , "0.8" , "Minimum peak-to-peak duration when classifying SO versus delta" );
  add_param( "SO" , "t-p2p-max" , "2" , "Maximum peak-to-peak duration when classifying SO versus delta" );
  add_param( "SO" , "pct" , "75" , "Percentile-based threshold for SO/delta classification" );
  add_param( "SO" , "pct-neg" , "75" , "Percentile threshold for negative peak amplitude" );
  add_param( "SO" , "pct-pos" , "75" , "Percentile threshold for positive peak amplitude" );
  add_param( "SO" , "neg2pos" , "" , "Use negative-to-positive zero crossings" );
  add_param( "SO" , "th-mean" , "" , "Use mean not median" );
  add_param( "SO" , "stats-median" , "" , "Use median (not mean) when reporting stats over SOs" );  
  add_param( "SO" , "SO-only" , "" , "Restrict retained events to those classified as SOs" );
  add_param( "SO" , "delta-only" , "" , "Restrict retained events to those classified as delta waves" );
  add_param( "SO" , "half-wave" , "" , "Detect generic half-wave events rather than full SO cycles" );
  add_param( "SO" , "negative-half-wave" , "" , "Restrict to negative half-wave events" );
  add_param( "SO" , "positive-half-wave" , "" , "Restrict to positive half-wave events" );
  add_param( "SO" , "so-annot" , "SO" , "Write detected events as annotations with this label" );
  add_param( "SO" , "cache-pos" , "so-pos" , "Cache positive-peak sample-points" );
  add_param( "SO" , "cache-neg" , "so-neg" , "Cache negative-peak sample-points" );
  add_param( "SO" , "dynam" , "" , "Emit per-event features suitable for downstream dynamics analyses" );
  
  add_param( "SO" , "tl" , "C3" , "Output signal time-locked to detected SOs" );
  add_param( "SO" , "onset" , "" , "Sync to SO onset for tl option" );
  add_param( "SO" , "pos" , "" , "Sync to positive peak for tl option" );
  add_param( "SO" , "window" , "2" , "Specify window size (seconds) for tl option" );
  
  add_table( "SO" , "CH" , "Channel-level statistics" );
  add_var( "SO" , "CH" , "SO" , "Number of SO detected" );
  add_var( "SO" , "CH" , "SO_RATE" , "SO per minute" );
  add_var( "SO" , "CH" , "SO_AMP_NEG" , "SO amplitude (negative peak)" );
  add_var( "SO" , "CH" , "SO_AMP_POS" , "SO amplitude (positive peak)" );
  add_var( "SO" , "CH" , "SO_AMP_P2P" , "SO peak-to-peak amplitude" );
  add_var( "SO" , "CH" , "SO_DUR" , "SO duration (secs)" );      
  add_var( "SO" , "CH" , "SO_DUR_NEG" , "Negative peak duration (secs)" );
  add_var( "SO" , "CH" , "SO_DUR_POS", "Positive peak duration (secs)" );
  add_var( "SO" , "CH" , "SO_TRANS" , "SO transition (secs)" );      
  add_var( "SO" , "CH" , "SO_TRANS_FREQ" , "SO transition freq (Hz)" );      
  add_var( "SO" , "CH" , "SO_SLOPE_POS1" , "Positive peak rising slope" );
  add_var( "SO" , "CH" , "SO_SLOPE_POS2" , "Positive peak falling slope" );
  add_var( "SO" , "CH" , "SO_SLOPE_NEG1" , "Negative peak falling slope" );
  add_var( "SO" , "CH" , "SO_SLOPE_NEG2" , "Negative peak rising slope" );
  add_var( "SO" , "CH" , "SO_SLOPE"      , "Negative peak rising slope" );
  add_var( "SO" , "CH" , "SO_TH_NEG" , "Negative peak threshold [mag]" );
  add_var( "SO" , "CH" , "SO_TH_P2P" , "Peak-to-peak threshold [mag]" );

  add_table( "SO" , "CH,E" , "Epoch-level statistics" );
  add_var( "SO" , "CH,E" , "N" , "Number of SO detected" );
  add_var( "SO" , "CH,E" , "AMP_NEG" , "Mean negative peak amplitude" );
  add_var( "SO" , "CH,E" , "AMP_POS" , "Mean positive peak amplitude" );
  add_var( "SO" , "CH,E" , "AMP_P2P" , "Mean peak-to-peak SO amplitude" );
  add_var( "SO" , "CH,E" , "DUR" , "Mean SO duration" );
  add_var( "SO" , "CH,E" , "DUR_POS" , "Mean positive HW duration" );
  add_var( "SO" , "CH,E" , "DUR_NEG" , "Mean negative HW duration" );
  add_var( "SO" , "CH,E" , "TRANS" , "Mean SO transition time (sec)" );
  add_var( "SO" , "CH,E" , "TRANS_FREQ" , "Mean SO transition freq (Hz)" );
  add_var( "SO" , "CH,E" , "SLOPE_NEG1" , "Negative peak falling slope" );
  add_var( "SO" , "CH,E" , "SLOPE_NEG2" , "Negative peak rising slope" );
  add_var( "SO" , "CH,E" , "SLOPE_POS1" , "Positive peak rising slope" );
  add_var( "SO" , "CH,E" , "SLOPE_POS2" , "Positive peak falling slope" );
  add_var( "SO" , "CH,E" , "SLOPE" , "Negative peak rising slope" );
  
  add_table( "SO" , "CH,N" , "Per-SO statistics" );
  add_var( "SO" , "CH,N" , "AMP_NEG" , "Negative peak amplitude" );  
  add_var( "SO" , "CH,N" , "AMP_POS" , "Positive peak ampltiude" );
  add_var( "SO" , "CH,N" , "AMP_P2P" , "Peak-to-peak ampltiude" ); 
  add_var( "SO" , "CH,N" , "DUR" , "SO duration" ); 
  add_var( "SO" , "CH,N" , "DUR_NEG" , "SO HW1 duration" ); 
  add_var( "SO" , "CH,N" , "DUR_POS" , "SO HW2 duration" ); 
  add_var( "SO" , "CH,N" , "TRANS" , "SO transition (sec)" );
  add_var( "SO" , "CH,N" , "SLOPE_POS1" , "Positive peak rising slope" ); 
  add_var( "SO" , "CH,N" , "SLOPE_POS2" , "Positive peak falling slope" ); 
  add_var( "SO" , "CH,N" , "SLOPE_NEG1" , "Negative peak falling slope" );
  add_var( "SO" , "CH,N" , "SLOPE_NEG2" , "Negative peak rising slope" );
  add_var( "SO" , "CH,N" , "SLOPE" , "Negative peak rising slope" );
  add_var( "SO" , "CH,N" , "IDX_NEG" , "Negative peak sample index" );
  add_var( "SO" , "CH,N" , "IDX_POS" , "Positive peak sample index" ); 
  add_var( "SO" , "CH,N" , "START" , "Start of SO (in seconds elapsed from start of EDF)" ); 
  add_var( "SO" , "CH,N" , "START_IDX" , "Start of SO (in sample-point units)" ); 
  add_var( "SO" , "CH,N" , "STOP" , "Stop of SO (in seconds elapsed from start of EDF)" ); 
  add_var( "SO" , "CH,N" , "STOP_IDX" , "Stop of SO (in sample-point units)" ); 
  
  add_table( "SO" , "CH,CH2,SP" , "SO time-locked signal averaging [tl]" );
  add_var( "SO" , "CH,CH2,SP" , "SOTL" , "SO time-locked signal average" );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // CROSS-SIGNAL ANALYSES
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // COH
  //

  add_cmd( "cc" , "COH" , "Pairwise channel coherence" );
  add_url( "COH" , "cross-signal-analysis/#coh" );
  add_verb( "COH" ,
            "Estimate pairwise linear coupling in the frequency domain via "
            "cross-spectral coherence.\n\n"
            "COH computes auto- and cross-spectra for pairs of channels and then "
            "summarizes their normalized frequency-specific dependence as "
            "magnitude-squared coherence. Luna also reports imaginary and lagged "
            "coherence variants, which are less sensitive to instantaneous shared "
            "components and volume-conduction-like effects." );
  
  add_param( "COH" , "sig" , "C3,C4" , "Restrict analysis to these channels (all-by-all pairs)" );
  add_param( "COH" , "sig1" , "C3,C4" , "Restrict analysis to sig1 x sig2 channel pairs only" );
  add_param( "COH" , "sig2" , "F3,F4" , "Restrict analysis to sig1 x sig2 channel pairs only" );

  add_param( "COH" , "sr" , "125" , "Set sample rate (i.e. if different for some channels)" );
  add_param( "COH" , "spectrum" , "" , "Show full coherence spectra as well as bands" );
  add_param( "COH" , "max" , "50" , "Upper frequency for spectra" );
  add_param( "COH" , "epoch" , "" , "Show per-epoch coherence" );
  add_param( "COH" , "epoch-spectrum" , "" , "Show per-epoch full coherence spectra" );

  add_table( "COH" , "B,CH1,CH2" , "Coherence for power bands" );
  add_var( "COH" , "B,CH1,CH2" , "COH" , "Magnitude-squared coherence" );
  add_var( "COH" , "B,CH1,CH2" , "ICOH" , "Imaginary coherence" );
  add_var( "COH" , "B,CH1,CH2" , "LCOH" , "Lagged coherence" );
  
  add_table( "COH" , "F,CH1,CH2" , "Full cross-spectra coherence [spectrum]" );
  add_var( "COH" , "F,CH1,CH2" , "COH" , "Magnitude-squared coherence" );
  add_var( "COH" , "F,CH1,CH2" , "ICOH" , "Imaginary coherence" );
  add_var( "COH" , "F,CH1,CH2" , "LCOH" , "Lagged coherence" );
  add_var( "COH" , "F,CH1,CH2" , "CSPEC" , "Cross-spectrum in dB" );

  add_table( "COH" , "B,CH1,CH2,E" , "Epoch-level band coherence" );
  add_var( "COH" , "B,CH1,CH2,E" , "COH" , "Magnitude-squared coherence" );
  add_var( "COH" , "B,CH1,CH2,E" , "ICOH" , "Imaginary coherence" );
  add_var( "COH" , "B,CH1,CH2,E" , "LCOH" , "Lagged coherence" );
  
  add_table( "COH" , "CH1,CH2,E,F" , "Epoch-level coherence" );
  add_var( "COH" , "CH1,CH2,E,F" , "COH" , "Magnitude-squared coherence" );
  add_var( "COH" , "CH1,CH2,E,F" , "ICOH" , "Imaginary coherence" );
  add_var( "COH" , "CH1,CH2,E,F" , "LCOH" , "Lagged coherence" );
  add_var( "COH" , "CH1,CH2,E,F" , "CSPEC" , "Cross-spectrum in dB" );

  // as these files can get large...
  set_compressed( "COH" , tfac_t( "CH1,CH2,B,E" ) );
  set_compressed( "COH" , tfac_t( "CH1,CH2,F,E" ) );

  //
  // PSI
  //

  add_cmd( "cc" , "PSI" , "Phase slope index" );
  add_url( "PSI" , "cross-signal-analysis/#psi" );
  add_verb( "PSI" ,
            "Estimate directed connectivity with the phase slope index.\n\n"
            "PSI summarizes whether the phase difference between two channels tends "
            "to increase or decrease with frequency across a local band. A "
            "consistent slope implies a preferred lead-lag relationship, making PSI "
            "useful for directionality analyses that are more robust than raw phase "
            "differences at a single frequency." );

  add_param( "PSI" , "sig" , "C3,C4" , "Restrict analysis to these channels (all-by-all pairs)" );
  add_param( "PSI" , "epoch" , "" , "Epoch level analysis" );

  add_param( "PSI" , "f" , "3" , "Frequency center(s)" );
  add_param( "PSI" , "f-lwr" , "3" , "Lower frequency range" );
  add_param( "PSI" , "f-upr" , "25" , "Upper frequency range" );
  add_param( "PSI" , "w" , "5" , "Window width (Hz)" );
  add_param( "PSI" , "r" , "1" , "Window increment (Hz)" );

  add_table( "PSI" , "F" , "Phase-slope index parameters" );
  add_var( "PSI" , "F" , "F1" , "Lower frequency bound" );
  add_var( "PSI" , "F" , "F2" , "Upper frequency bound" );
  add_var( "PSI" , "F" , "NF" , "NUmber of frequency bins" );

  add_table( "PSI" , "F,CH" , "Net (single-channel) Phase-slope index parameters" );
  add_var( "PSI" , "F,CH" , "PSI" , "Net Phase-slope index, normalized" );
  add_var( "PSI" , "F,CH" , "PSI_RAW" , "Net Phase-slope index, raw" );
  add_var( "PSI" , "F,CH" , "STD" , "Net Phase-slope index, SD" );
  add_var( "PSI" , "F,CH" , "APSI" , "Absolute net phase-slope index, normalized" );
  add_var( "PSI" , "F,CH" , "APSI_RAW" , "Absolute net phase-slope index, raw" );
  add_var( "PSI" , "F,CH" , "ASTD" , "Absolute net phase-slope index, SD" );

  add_table( "PSI" , "F,CH1,CH2" , "Pairwise Phase-slope index parameters" );
  add_var( "PSI" , "F,CH1,CH2" , "PSI" , "Phase-slope index, normalized" );
  add_var( "PSI" , "F,CH1,CH2" , "PSI_RAW" , "Phase-slope index, raw" );
  add_var( "PSI" , "F,CH1,CH2" , "STD" , "Phase-slope index, SD" );

  add_table( "PSI" , "E,F,CH" , "Epoch-level net (single-channel) Phase-slope index parameters" );
  add_var( "PSI" , "E,F,CH" , "PSI" , "Net Phase-slope index, normalized" );
  add_var( "PSI" , "E,F,CH" , "PSI_RAW" , "Net Phase-slope index, raw" );
  add_var( "PSI" , "E,F,CH" , "STD" , "Net Phase-slope index, SD" );

  add_table( "PSI" , "E,F,CH1,CH2" , "Epoch-level pairwise Phase-slope index parameters" );
  add_var( "PSI" , "E,F,CH1,CH2" , "PSI" , "Phase-slope index, normalized" );
  add_var( "PSI" , "E,F,CH1,CH2" , "PSI_RAW" , "Phase-slope index, raw" );
  add_var( "PSI" , "E,F,CH1,CH2" , "STD" , "Phase-slope index, SD" );

  set_compressed( "PSI" , tfac_t( "E,F,CH" ) );
  set_compressed( "PSI" , tfac_t( "E,F,CH1,CH2" ) );

  //
  // SYNC
  //

  add_cmd( "exp" , "SYNC" , "Global phase synchrony" , true );
  add_url( "SYNC" , "cross-signal-analysis/#sync" );
  add_verb( "SYNC" ,
            "SYNC estimates multichannel phase synchrony using the Kuramoto order "
            "parameter (KOP). Unlike pairwise methods such as IPC, XCORR, or "
            "TSYNC, it collapses across a set of channels to summarize how tightly "
            "their phases align as a group.\n\n"
            "Two methods are available. In fft mode, SYNC computes a per-epoch FFT "
            "for each channel and derives KOP at each frequency bin from the phase "
            "angles. In the filter-Hilbert mode, it band-pass filters the full "
            "signals in one or more frequency bands, extracts analytic phase, and "
            "then summarizes the mean KOP within each epoch. The Hilbert path can "
            "also add new KOP signals back into the EDF." );
  
  add_param( "SYNC" , "sig" , "C3,C4" , "Restrict analysis to these channels" );
  add_param( "SYNC" , "fft" , "" , "Use per-epoch FFT phase instead of filter-Hilbert phase" );
  add_param( "SYNC" , "min" , "0.5" , "Lower frequency bound for fft mode" );
  add_param( "SYNC" , "max" , "20" , "Upper frequency bound for fft mode" );
  add_param( "SYNC" , "f-lwr" , "0.5" , "Lower bound(s) of Hilbert frequency bands" );
  add_param( "SYNC" , "f-upr" , "4" , "Upper bound(s) of Hilbert frequency bands" );
  add_param( "SYNC" , "f" , "10" , "Center frequency or frequencies for Hilbert mode" );
  add_param( "SYNC" , "w" , "3" , "Band half-width or sliding-band width for Hilbert mode" );
  add_param( "SYNC" , "r" , "1" , "Step size when constructing sliding Hilbert bands" );
  add_param( "SYNC" , "tw" , "1" , "FIR transition width for Hilbert filtering" );
  add_param( "SYNC" , "ripple" , "0.025" , "FIR ripple for Hilbert filtering" );
  add_param( "SYNC" , "tag" , "KOP" , "Prefix for new KOP channels added to the EDF" );
  add_param( "SYNC" , "no-new-channels" , "" , "Do not add KOP channels back into the EDF" );

  add_table( "SYNC" , "E,F" , "Epoch-wise analysis" );
  add_var( "SYNC" , "E,F" , "KOP" , "Kuramoto order parameter for that epoch and frequency/frequency-band" );
  set_compressed( "SYNC" , tfac_t( "E,F" ) );

  //
  // XCORR
  //

  add_cmd( "cc" , "XCORR" , "Cross-correlation" );
  add_url( "XCORR" , "cross-signal-analysis/#xcorr" );
  add_verb( "XCORR" ,
            "Estimate pairwise temporal delay by maximizing cross-correlation "
            "between channels.\n\n"
            "XCORR operates directly in the time domain, scanning over positive and "
            "negative lags to find the offset at which two signals align most "
            "strongly. It can report delay estimates from whole-record averages, "
            "epoch-level analyses, and the full lag profile." );

  add_param( "XCORR" , "sig" , "C3,C4" , "Restrict analysis to these channel" );
  add_param( "XCORR" , "w" , "10" , "Restrict to window of +/- 10 seconds" );
  add_param( "XCORR" , "verbose" , "10" , "Restrict to window of +/- 10 seconds" );
  add_param( "XCORR" , "epoch" , "10" , "Epoch-level outputs" );

  add_table( "XCORR" , "CH1,CH2" , "Pairwise outputs" );
  add_var( "XCORR" , "CH1,CH2" , "D" , "Delay in seconds (+ve means CH1 lags CH2)" );
  add_var( "XCORR" , "CH1,CH2" , "S" , "Delay in samples" );
  add_var( "XCORR" , "CH1,CH2" , "D_MN" , "Delay in seconds, mean over epochs" );
  add_var( "XCORR" , "CH1,CH2" , "S_MN" , "Delay in samples, mean over epochs" );
  add_var( "XCORR" , "CH1,CH2" , "D_MD" , "Delay in seconds, median over epochs" );
  add_var( "XCORR" , "CH1,CH2" , "S_MD" , "Delay in samples, median over epochs" );
  
  add_table( "XCORR" , "CH1,CH2,D" , "Lag-wise x-corrs" );
  add_var( "XCORR" , "CH1,CH2,D" , "T" , "Lag time (seconds)" );
  add_var( "XCORR" , "CH1,CH2,D" , "XCORR" , "Cross correlation" );

  add_table( "XCORR" , "CH1,CH2,E" , "Epoch-wise outputs" );
  add_var( "XCORR" , "CH1,CH2,E" , "D" , "Delay in seconds" );

  //
  // TSYNC
  //

  add_cmd( "cc" , "TSYNC" , "Hilbert-based phase-delay synchrony" );
  add_url( "TSYNC" , "cross-signal-analysis/#tsync" );
  add_verb( "TSYNC" ,
            "Estimate pairwise synchrony from Hilbert-derived phase and amplitude "
            "signals.\n\n"
            "Unlike XCORR, which works directly on the raw time series and searches "
            "for the lag that maximizes cross-correlation, TSYNC assumes analytic "
            "signals already exist (`*_ht_ph` and `*_ht_mag`) and then evaluates how "
            "phase locking changes as one signal is shifted relative to the other. "
            "Unlike IPC, which summarizes zero-lag phase coherence only, TSYNC "
            "returns a lag profile of phase locking across a delay window." );
  add_param( "TSYNC" , "sig" , "C3,C4" , "Explicit list of signal roots to analyze" );
  add_param( "TSYNC" , "w" , "0.25" , "Half-window in seconds for delay estimation" );
  add_param( "TSYNC" , "verbose" , "" , "Emit full delay profiles" );
  add_param( "TSYNC" , "ht" , "" , "Use associated _ht_ph and _ht_mag signals (recommended mode)" );
  add_table( "TSYNC" , "CH1,CH2,D" , "Lag-wise synchrony profiles" );
  add_var( "TSYNC" , "CH1,CH2,D" , "PHL" , "Average phase-locking value at this lag" );

  //
  // GP
  //

  add_cmd( "cc" , "GP" , "Granger prediction analysis" );
  add_url( "GP" , "cross-signal-analysis/#gp" );
  add_verb( "GP" ,
            "Estimate directed predictability between channels using Granger-style "
            "autoregressive modeling.\n\n"
            "GP fits univariate and bivariate autoregressive models within each "
            "epoch, compares prediction errors to quantify directed influence in the "
            "time domain, and can optionally decompose that influence over "
            "frequency. The result is an asymmetric X->Y versus Y->X measure of "
            "predictive information flow." );
  add_param( "GP" , "sig" , "C3,C4" , "Signals to analyze, all pairwise" );
  add_param( "GP" , "w" , "500" , "Analysis window in milliseconds" );
  add_param( "GP" , "order" , "50" , "Autoregressive model order in milliseconds" );
  add_param( "GP" , "bic" , "20" , "Evaluate BIC up to this AR order" );
  add_param( "GP" , "f" , "1,25,20" , "Linear frequency grid: lower,upper,count" );
  add_param( "GP" , "f-log" , "1,25,20" , "Log-spaced frequency grid: lower,upper,count" );
  add_table( "GP" , "CH1,CH2" , "Time-domain Granger summaries" );
  add_var( "GP" , "CH1,CH2" , "Y2X" , "Directed predictability from CH2 to CH1" );
  add_var( "GP" , "CH1,CH2" , "X2Y" , "Directed predictability from CH1 to CH2" );
  add_table( "GP" , "CH1,CH2,F" , "Frequency-resolved Granger summaries" );
  add_var( "GP" , "CH1,CH2,F" , "Y2X" , "Frequency-domain predictability from CH2 to CH1" );
  add_var( "GP" , "CH1,CH2,F" , "X2Y" , "Frequency-domain predictability from CH1 to CH2" );
  add_table( "GP" , "CH1,CH2,E" , "Epoch-level Granger summaries" );
  add_var( "GP" , "CH1,CH2,E" , "Y2X" , "Directed predictability from CH2 to CH1" );
  add_var( "GP" , "CH1,CH2,E" , "X2Y" , "Directed predictability from CH1 to CH2" );
  add_var( "GP" , "CH1,CH2,E" , "BIC" , "Best BIC score [bic]" );
  add_table( "GP" , "CH1,CH2,E,F" , "Epoch-level frequency-resolved Granger summaries" );
  add_var( "GP" , "CH1,CH2,E,F" , "Y2X" , "Frequency-domain predictability from CH2 to CH1" );
  add_var( "GP" , "CH1,CH2,E,F" , "X2Y" , "Frequency-domain predictability from CH1 to CH2" );

  //
  // CORREL
  //

  add_cmd( "cc" , "CORREL" , "Pairwise channel correlation" );
  add_url( "CORREL" , "cross-signal-analysis/#correl" );
  add_verb( "CORREL" ,
            "Estimate pairwise linear association between channels using Pearson "
            "correlation.\n\n"
            "CORREL compares channels directly in the time domain, either over the "
            "whole available trace or epoch by epoch. In addition to pairwise "
            "coefficients, Luna can summarize each channel's correlation profile "
            "across the rest of the montage." );
  
  add_param( "CORREL" , "sig" , "C3,C4" , "Restrict analysis to these channels (all-by-all pairs)" );
  add_param( "CORREL" , "sig1" , "C3,C4" , "Restrict analysis to sig1 x sig2 channel pairs only" );
  add_param( "CORREL" , "sig2" , "F3,F4" , "Restrict analysis to sig1 x sig2 channel pairs only" );
  
  add_param( "CORREL" , "sr" , "128" , "Resample channels to this sample rate if needed" );
  add_param( "CORREL" , "epoch" , "" , "Display per-epoch, and estimate mean and median correlation across epochs" );

  add_param( "CORREL" , "ch-low" , "0.1" , "Number of correlations below threshold for this channel" );
  add_param( "CORREL" , "ch-high" , "0.98" , "Number of correlations above threshold for this channel" );
  
  add_table( "CORREL" , "CH1,CH2" , "Whole-signal correlations for pairs of channels" );
  add_var( "CORREL" , "CH1,CH2" , "R", "Pearson product moment correlation" );
  add_var( "CORREL" , "CH1,CH2" , "R_MEAN", "(If epoch is specified) the mean of epoch-level correlations" );
  add_var( "CORREL" , "CH1,CH2" , "R_MEDIAN" ,  "(If epoch is specified) the median of epoch-level correlations" );

  add_table( "CORREL" , "CH" , "Channel-level summaries of correlations" );
  add_var( "CORREL" , "CH" , "SUMM_LOW", "Number of correlations below ch-low threshold" );
  add_var( "CORREL" , "CH" , "SUMM_HIGH", "Number of correlations above ch-high threshold" );
  
  add_var( "CORREL" , "CH" , "SUMM_MEAN", "Mean correlation for this channel" );
  add_var( "CORREL" , "CH" , "SUMM_N", "Number of channel pairs contributing to summaries" );
  add_var( "CORREL" , "CH" , "SUMM_MIN", "Min correlation for this channel" );
  add_var( "CORREL" , "CH" , "SUMM_MAX", "Max correlation for this channel" );

  add_table( "CORREL" , "CHS" , "Sets of highly correlated channels" );
  add_var( "CORREL" , "CHS" , "SET", "Comma-delimited channel set" );
  add_var( "CORREL" , "CHS" , "N", "Number of channels in this set" );
  
  add_table( "CORREL" , "CH1,CH2,E" , "Whole-signal correlations for pairs of channels" );
  add_var( "CORREL" , "CH1,CH2,E" , "R", "Pearson product moment correlation" );
  set_compressed( "CORREL" , tfac_t( "CH1,CH2,E" ) );

  add_table( "CORREL" , "" , "Whole-record correlation summaries" );
  add_var( "CORREL" , "" , "SUMM_HIGH_N", "Number of channels with correlations above ch-high" );
  add_var( "CORREL" , "" , "SUMM_HIGH_CHS", "Comma-delimited channels exceeding ch-high" );

  //
  // MI
  //

  add_cmd( "cc" , "MI" , "Mutual information" );
  add_url( "MI" , "cross-signal-analysis/#mi" );
  add_verb( "MI" ,
            "Estimate pairwise dependence with information-theoretic measures "
            "instead of linear correlation alone.\n\n"
            "MI discretizes paired signals, estimates the joint and marginal "
            "distributions, and then derives mutual information, total correlation, "
            "dual total correlation, and entropy summaries. This can capture "
            "nonlinear dependence that Pearson correlation would miss." );

  add_param( "MI" , "sig" , "C3,C4,F3,F4" , "Optionally specify channels (defaults to all)" );
  add_param( "MI" , "epoch" , "" , "Report MI and other measures per epoch" );
  add_param( "MI" , "scott" , "" , "Use Scott's rule to determine bin number" );
  add_param( "MI" , "sturges" , "" , "Use Sturges' rule to determine bin number" );
  add_param( "MI" , "permute" , "1000" , "Estimate empirical significance via permutation, with N replicates" );

  add_table( "MI" , "CH1,CH2" , "Output for the whole signal pairs" );
  add_var( "MI" , "CH1,CH2" , "MI" , "Mutual information" );
  add_var( "MI" , "CH1,CH2" , "TOTCORR" , "Total correlation" );
  add_var( "MI" , "CH1,CH2" , "DTOTCORR" , "Dual total correlation" );
  add_var( "MI" , "CH1,CH2" , "JINF" , "Joint entropy" );
  add_var( "MI" , "CH1,CH2" , "INFA" , "Marginal entropy of first signal" );
  add_var( "MI" , "CH1,CH2" , "INFB" , "Marginal entropy of second signal" );
  add_var( "MI" , "CH1,CH2" , "NBINS" , "Number of bins" );
  add_var( "MI" , "CH1,CH2" , "EMP" , "Empirical significance [permute]" );
  add_var( "MI" , "CH1,CH2" , "Z" , "Z statistic [permute]" );

  add_table( "MI" , "CH1,CH2,E" , "Output per epoch" );
  add_var( "MI" , "CH1,CH2,E" , "MI" , "Mutual information" );
  add_var( "MI" , "CH1,CH2,E" , "TOTCORR" , "Total correlation" );
  add_var( "MI" , "CH1,CH2,E" , "DTOTCORR" , "Dual total correlation" );
  add_var( "MI" , "CH1,CH2,E" , "JINF" , "Joint entropy" );
  add_var( "MI" , "CH1,CH2,E" , "INFA" , "Marginal entropy of first signal" );
  add_var( "MI" , "CH1,CH2,E" , "INFB" , "Marginal entropy of second signal" );
  set_compressed( "MI" , tfac_t( "CH1,CH2,E" ) );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // INTERVAL-BASED ANALYSIS
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // TLOCK
  //

  add_cmd( "interval" , "TLOCK" , "Time-locked averaging" );
  add_url( "TLOCK" , "intervals/#tlock" );
  add_verb( "TLOCK" ,
            "Builds fixed-width windows around either cached sample-point events\n"
            "or existing epochs, then summarizes the aligned traces across time.\n"
            "This is a generic event-locked averaging command rather than a\n"
            "spindle-specific tool.\n"
            "\n"
            "In cache mode, TLOCK reads sample indices from a cache such as one\n"
            "created by PEAKS or Z-PEAKS and extracts a symmetric window around\n"
            "each point. In epoch mode, each epoch becomes one aligned window,\n"
            "with an optional redefinition of where time zero sits via mid.\n"
            "\n"
            "Current structured output is the number of accepted intervals plus,\n"
            "for each time offset, the mean, SD and median. If the cache carries\n"
            "extra strata, those are preserved as additional output factors." );
  add_param( "TLOCK" , "sig" , "C3,C4" , "One or more signals to summarize" );
  add_param( "TLOCK" , "cache" , "c1" , "Cache of sample-point events to align to" );
  add_param( "TLOCK" , "epoch" , "" , "Use epochs rather than cached sample-point events" );
  add_param( "TLOCK" , "w" , "2" , "Half-window size in seconds around each cached event" );
  add_param( "TLOCK" , "mid" , "0" , "In epoch mode, redefine time zero within the epoch" );
  add_param( "TLOCK" , "same-channel" , "" , "In cache mode, only use seed strata that match the readout channel" );
  add_param( "TLOCK" , "channel-postfix" , "_SIGMA" , "Allow same-channel matching after adding a postfix to readout channels" );
  add_param( "TLOCK" , "seed-postfix" , "_SIGMA" , "Allow same-channel matching after adding a postfix to seed channels" );
  add_param( "TLOCK" , "tolog" , "" , "Take logs of the input signal before averaging" );
  add_param( "TLOCK" , "np" , "0.1" , "Normalize each interval using edge regions spanning this fraction of the window" );
  add_param( "TLOCK" , "zero" , "" , "After edge normalization, center each interval at zero" );
  add_param( "TLOCK" , "th" , "4" , "Remove interval outliers beyond this SD threshold before averaging" );
  add_param( "TLOCK" , "win" , "4" , "Winsorize interval outliers beyond this SD threshold before averaging" );

  add_table( "TLOCK" , "CH" , "Time-locked summary counts" );
  add_var( "TLOCK" , "CH" , "N" , "Number of accepted aligned intervals" );
  add_table( "TLOCK" , "CH,SEC" , "Time-locked signal summaries by offset" );
  add_var( "TLOCK" , "CH,SEC" , "M" , "Mean aligned signal value" );
  add_var( "TLOCK" , "CH,SEC" , "SD" , "SD of aligned signal values" );
  add_var( "TLOCK" , "CH,SEC" , "MD" , "Median aligned signal value" );

  //
  // PEAKS
  //

  add_cmd( "interval" , "PEAKS" , "Local extrema finder" );
  add_url( "PEAKS" , "intervals/#peaks" );
  add_verb( "PEAKS" ,
            "Finds local maxima by default, or minima if requested, within one\n"
            "or more signals. The command can work on the whole trace or within\n"
            "epochs, and it can ignore saturated plateaus by masking clipped\n"
            "runs before searching for extrema.\n"
            "\n"
            "PEAKS does not emit tabular writer output. Instead it is mainly used\n"
            "to populate a cache of sample indices for downstream commands such\n"
            "as TLOCK, or to create annotations around detected peaks." );
  add_param( "PEAKS" , "sig" , "C3,C4" , "One or more signals in which to find extrema" );
  add_param( "PEAKS" , "cache" , "c1" , "Optionally write detected peak sample points to a cache" );
  add_param( "PEAKS" , "annot" , "PKS" , "Optionally write detected peaks to an annotation class" );
  add_param( "PEAKS" , "w" , "0.1" , "If writing annotations, expand each peak by +/- this many seconds" );
  add_param( "PEAKS" , "epoch" , "" , "Find extrema separately within each epoch" );
  add_param( "PEAKS" , "clipped" , "3" , "Ignore clipped runs of at least this many identical extreme samples" );
  add_param( "PEAKS" , "min" , "" , "Also find minima" );
  add_param( "PEAKS" , "min-only" , "" , "Find minima only" );
  add_param( "PEAKS" , "percentile" , "20" , "Optionally retain only the top percentile of detected extrema" );

  //
  // Z-PEAKS
  //

  add_cmd( "interval" , "Z-PEAKS" , "Smoothed z-score peak detector" );
  add_url( "Z-PEAKS" , "intervals/#z-peaks" );
  add_verb( "Z-PEAKS" ,
            "Detects peaks using a moving-baseline z-score heuristic. A local\n"
            "running mean and SD are estimated over a lag window, candidate peak\n"
            "regions are called when the signal exceeds a z threshold, and those\n"
            "regions can optionally be expanded with a secondary flank threshold.\n"
            "\n"
            "Like PEAKS, this command is mainly used to create annotations or a\n"
            "cache of peak sample points for later interval-based analyses." );
  add_param( "Z-PEAKS" , "sig" , "C3,C4" , "One or more signals in which to detect peaks" );
  add_param( "Z-PEAKS" , "w" , "5" , "Lag window in seconds for the moving baseline" );
  add_param( "Z-PEAKS" , "th" , "3.5" , "Primary z-score threshold for the core peak region" );
  add_param( "Z-PEAKS" , "sec" , "0.2" , "Minimum duration in seconds for the core region" );
  add_param( "Z-PEAKS" , "influence" , "0.01" , "How much detected peaks influence the moving baseline" );
  add_param( "Z-PEAKS" , "max" , "10" , "Optional upper z-score limit for accepted peaks" );
  add_param( "Z-PEAKS" , "th2" , "2" , "Optional secondary threshold for flanking regions" );
  add_param( "Z-PEAKS" , "sec2" , "0.1" , "Minimum duration in seconds for flanking regions" );
  add_param( "Z-PEAKS" , "negatives" , "" , "Also allow negative-going excursions" );
  add_param( "Z-PEAKS" , "annot" , "ZPKS" , "Optionally write detected peaks to an annotation class" );
  add_param( "Z-PEAKS" , "add-flanking" , "0.1" , "If writing annotations, extend each detected interval on both sides" );
  add_param( "Z-PEAKS" , "cache" , "c1" , "Optionally write peak sample points to a cache" );

  //
  // OVERLAP
  //

  add_cmd( "interval" , "OVERLAP" , "Permutation overlap statistics" );
  add_url( "OVERLAP" , "intervals/#overlap" );
  add_verb( "OVERLAP" ,
            "Compares seed annotations against other annotations using overlap,\n"
            "nearest-neighbour distance, and optional offset-based enrichment\n"
            "profiles. Significance is assessed by permutation within background\n"
            "intervals or, optionally, by local event-based permutation.\n"
            "\n"
            "The command can pool or preserve channels, reduce annotations to\n"
            "starts, stops, midpoints or relative-position markers, constrain the\n"
            "permutation space with background intervals, and report pairwise,\n"
            "seed-wise and pileup summaries." );
  add_param( "OVERLAP" , "seed" , "SP,SO" , "Seed annotation class or classes" );
  add_param( "OVERLAP" , "other" , "AROUSAL" , "Optional non-seed annotation classes to compare against" );
  add_param( "OVERLAP" , "bg" , "NREM" , "Background intervals that define where events can be evaluated or permuted" );
  add_param( "OVERLAP" , "xbg" , "ART" , "Exclusionary background intervals to remove from bg" );
  add_param( "OVERLAP" , "edges" , "0.5" , "Trim this many seconds from both edges of each background interval" );
  add_param( "OVERLAP" , "nreps" , "1000" , "Number of permutations" );
  add_param( "OVERLAP" , "w" , "10" , "Window in seconds for nearest-neighbour distance summaries" );
  add_param( "OVERLAP" , "mw" , "60" , "Window in seconds for marker-based distances" );
  add_param( "OVERLAP" , "overlap" , "0" , "Minimum proportional overlap required to count as a match" );
  add_param( "OVERLAP" , "offset" , "0.5,0.5,5" , "Offset-window grid as size,increment,max in seconds" );
  add_param( "OVERLAP" , "offset-seed" , "SP" , "Restrict offset analyses to these seed annotations" );
  add_param( "OVERLAP" , "offset-other" , "SO" , "Restrict offset analyses to these other annotations" );
  add_param( "OVERLAP" , "pool-channels" , "" , "Pool annotation channels into shared annotation classes" );
  add_param( "OVERLAP" , "within-channel" , "" , "Only compare annotations within the same channel" );
  add_param( "OVERLAP" , "midpoint" , "SP" , "Reduce all or selected annotations to their midpoints" );
  add_param( "OVERLAP" , "start" , "SP" , "Reduce all or selected annotations to their start times" );
  add_param( "OVERLAP" , "stop" , "SP" , "Reduce all or selected annotations to their stop times" );
  add_param( "OVERLAP" , "rp" , "SO|peak" , "Reduce selected annotations to a relative-position meta-data marker" );
  add_param( "OVERLAP" , "f" , "0.5" , "Add flanking time to seeds before analysis" );
  add_param( "OVERLAP" , "d2-quant" , "" , "Use lead-vs-lag balance rather than signed seconds for D2" );
  add_param( "OVERLAP" , "dist-excludes-overlapping" , "" , "Exclude overlapping events from nearest-distance calculations" );
  add_param( "OVERLAP" , "seed-seed" , "" , "Include seed-vs-seed comparisons" );
  add_param( "OVERLAP" , "pileup" , "" , "Report seed-group pileup summaries" );
  add_param( "OVERLAP" , "shuffle-others" , "" , "Permute non-seed annotations as well as seeds" );
  add_param( "OVERLAP" , "fixed" , "SP" , "Do not permute these annotations" );
  add_param( "OVERLAP" , "max-shuffle" , "30" , "Constrain shuffles to within this many seconds" );
  add_param( "OVERLAP" , "event-perm" , "5" , "Use local event-based permutation within this neighbourhood" );
  add_param( "OVERLAP" , "marker" , "LIGHTS_OFF" , "Marker annotations for event-permutation summaries" );
  add_param( "OVERLAP" , "align" , "SP_C3|SP_C4" , "Permute aligned annotation groups together" );
  add_param( "OVERLAP" , "ordered" , "" , "Preserve order in seed-seed group summaries" );
  add_param( "OVERLAP" , "contrasts" , "A|B-C|D" , "Evaluate user-defined contrasts between overlap terms" );
  add_param( "OVERLAP" , "matched" , "MATCHED" , "Create a new annotation containing matched seeds" );
  add_param( "OVERLAP" , "unmatched" , "UNMATCHED" , "Create a new annotation containing unmatched seeds" );
  add_param( "OVERLAP" , "m-count" , "1" , "Require at least this many matches to count as matched" );
  add_param( "OVERLAP" , "chs-inc" , "C3,C4" , "Restrict to these annotation channels" );
  add_param( "OVERLAP" , "chs-exc" , "EMG" , "Exclude these annotation channels" );
  add_param( "OVERLAP" , "flt" , "wgt,10,." , "Filter annotations by numeric meta-data limits" );
  add_param( "OVERLAP" , "add-shuffled-annots" , "SP" , "After the first permutation, write shuffled annotations back into the EDF" );
  add_param( "OVERLAP" , "shuffled-annots-tag" , "s_" , "Prefix for shuffled annotation classes" );
  add_param( "OVERLAP" , "a-list" , "annots.lst" , "In multi-sample mode, list of IDs and annotation files" );
  add_param( "OVERLAP" , "merged" , "merged" , "In multi-sample mode, basename for the temporary merged annotation file" );
  add_param( "OVERLAP" , "set-id" , "" , "In multi-sample mode, add source individual IDs as annotation meta-data" );
  add_param( "OVERLAP" , "set-id-key" , "ID" , "Meta-data key to use for multi-sample IDs" );
  add_param( "OVERLAP" , "keep-metadata" , "" , "Retain original annotation meta-data when merging multi-sample input" );

  add_table( "OVERLAP" , "SEEDS" , "Seed-group pileup summaries" );
  add_var( "OVERLAP" , "SEEDS" , "OBS" , "Observed pileup count for this seed combination" );
  add_var( "OVERLAP" , "SEEDS" , "EXP" , "Expected pileup count under permutation" );
  add_var( "OVERLAP" , "SEEDS" , "P" , "Empirical p-value for pileup enrichment" );
  add_var( "OVERLAP" , "SEEDS" , "Z" , "Z-score for pileup enrichment" );

  add_table( "OVERLAP" , "SEED" , "Seed-level overlap summaries" );
  add_var( "OVERLAP" , "SEED" , "PROP" , "Observed proportion of seed time overlapped by any other event" );
  add_var( "OVERLAP" , "SEED" , "PROP_EXP" , "Expected proportion of seed time overlapped under permutation" );
  add_var( "OVERLAP" , "SEED" , "PROP_P" , "Empirical p-value for proportional overlap" );
  add_var( "OVERLAP" , "SEED" , "PROP_Z" , "Z-score for proportional overlap" );

  add_table( "OVERLAP" , "SEED,OTHERS" , "Seed versus pooled-other overlap counts" );
  add_var( "OVERLAP" , "SEED,OTHERS" , "N_OBS" , "Observed count of pooled-other overlaps" );
  add_var( "OVERLAP" , "SEED,OTHERS" , "N_EXP" , "Expected count of pooled-other overlaps under permutation" );
  add_var( "OVERLAP" , "SEED,OTHERS" , "N_Z" , "Z-score for pooled-other overlap count" );

  add_table( "OVERLAP" , "SEED,OTHER" , "Seed-other pairwise overlap and distance summaries" );
  add_var( "OVERLAP" , "SEED,OTHER" , "SEED_ANNOT" , "Seed annotation class label" );
  add_var( "OVERLAP" , "SEED,OTHER" , "SEED_CH" , "Seed annotation channel label" );
  add_var( "OVERLAP" , "SEED,OTHER" , "OTHER_ANNOT" , "Other annotation class label" );
  add_var( "OVERLAP" , "SEED,OTHER" , "OTHER_CH" , "Other annotation channel label" );
  add_var( "OVERLAP" , "SEED,OTHER" , "N_OBS" , "Observed number of overlapping events" );
  add_var( "OVERLAP" , "SEED,OTHER" , "N_EXP" , "Expected number of overlapping events under permutation" );
  add_var( "OVERLAP" , "SEED,OTHER" , "N_P" , "Empirical p-value for overlap count" );
  add_var( "OVERLAP" , "SEED,OTHER" , "N_Z" , "Z-score for overlap count" );
  add_var( "OVERLAP" , "SEED,OTHER" , "D1_OBS" , "Observed mean absolute nearest distance" );
  add_var( "OVERLAP" , "SEED,OTHER" , "D1_EXP" , "Expected mean absolute nearest distance" );
  add_var( "OVERLAP" , "SEED,OTHER" , "D1_P" , "Empirical p-value for absolute nearest distance" );
  add_var( "OVERLAP" , "SEED,OTHER" , "D1_Z" , "Z-score for absolute nearest distance" );
  add_var( "OVERLAP" , "SEED,OTHER" , "D_N" , "Number of seed events contributing to nearest-distance summaries" );
  add_var( "OVERLAP" , "SEED,OTHER" , "D_N_EXP" , "Expected number of seed events contributing to nearest-distance summaries" );
  add_var( "OVERLAP" , "SEED,OTHER" , "D2_OBS" , "Observed signed distance summary" );
  add_var( "OVERLAP" , "SEED,OTHER" , "D2_SEC" , "Observed mean signed nearest distance in seconds" );
  add_var( "OVERLAP" , "SEED,OTHER" , "D2_LEAD" , "Proportion of seeds for which the other event leads" );
  add_var( "OVERLAP" , "SEED,OTHER" , "D2_LAG" , "Proportion of seeds for which the other event lags" );
  add_var( "OVERLAP" , "SEED,OTHER" , "D2_EXP" , "Expected signed distance summary under permutation" );
  add_var( "OVERLAP" , "SEED,OTHER" , "D2_P" , "Empirical p-value for the signed distance summary" );
  add_var( "OVERLAP" , "SEED,OTHER" , "D2_Z" , "Z-score for the signed distance summary" );

  add_table( "OVERLAP" , "SEED,OTHER,OFFSET" , "Offset-window overlap enrichment profiles" );
  add_var( "OVERLAP" , "SEED,OTHER,OFFSET" , "INT" , "Human-readable description of the offset window" );
  add_var( "OVERLAP" , "SEED,OTHER,OFFSET" , "N_OBS" , "Observed overlap count in this offset window" );
  add_var( "OVERLAP" , "SEED,OTHER,OFFSET" , "N_EXP" , "Expected overlap count in this offset window" );
  add_var( "OVERLAP" , "SEED,OTHER,OFFSET" , "N_P" , "Empirical p-value for offset-window overlap" );
  add_var( "OVERLAP" , "SEED,OTHER,OFFSET" , "N_Z" , "Z-score for offset-window overlap" );

  //
  // MEANS
  //

  add_cmd( "interval" , "MEANS" , "Signal means by annotation" );
  add_url( "MEANS" , "intervals/#means" );
  add_verb( "MEANS" ,
            "Computes signal means over the spans defined by one or more\n"
            "annotations. Means are reported by annotation class and, if\n"
            "requested, by annotation instance as well.\n"
            "\n"
            "The command can also compute left and right flanking means around\n"
            "each annotation, and for instance-level output it can optionally\n"
            "rescale means within an annotation class to a 0-1 range." );
  add_param( "MEANS" , "sig" , "C3,C4" , "One or more signals to summarize" );
  add_param( "MEANS" , "annot" , "FS,SS" , "One or more annotations" );
  add_param( "MEANS" , "w" , "5" , "Optional flanking window before/after each annotation" );
  add_param( "MEANS" , "by-instance" , "" , "Output means by annotation instance as well as by class" );
 
  add_table( "MEANS" , "CH,ANNOT" , "Annotation class means, by channel" );
  add_var( "MEANS" , "CH,ANNOT" , "M" , "Mean signal value across all instances of that annotation" );
  add_var( "MEANS" , "CH,ANNOT" , "S" , "Total annotation span in seconds" );
  add_var( "MEANS" , "CH,ANNOT" , "L" , "Left-flanking mean (if 'w' set)" );
  add_var( "MEANS" , "CH,ANNOT" , "R" , "Right-flanking mean (if 'w' set)" );

  add_table( "MEANS" , "CH,ANNOT,INST" , "Annotation class/instance means, by channel" );
  add_var( "MEANS" , "CH,ANNOT,INST" , "M" , "Mean signal value for that annotation instance" );  
  add_var( "MEANS" , "CH,ANNOT,INST" , "S" , "Annotation span in seconds" );
  add_var( "MEANS" , "CH,ANNOT,INST" , "L" , "Left-flanking mean (if 'w' set)" );
  add_var( "MEANS" , "CH,ANNOT,INST" , "R" , "Right-flanking mean (if 'w' set)" );
  add_var( "MEANS" , "CH,ANNOT,INST" , "M1" , "Within-class 0-1 normalized mean" );  

  /////////////////////////////////////////////////////////////////////////////////
  //
  // CROSS-FREQUENCY
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // IPC
  //

  add_cmd( "cc" , "IPC" , "Instantaneous phase coherence" );
  add_url( "IPC" , "cc/#ipc" );
  add_verb( "IPC" ,
	    "Computes pairwise zero-lag phase coupling from Hilbert-based\n"
	    "analytic signals. For each sample, IPC forms the phase\n"
	    "difference between two channels and summarizes it as signed\n"
	    "instantaneous phase coherence, a weighted variant, phase-locking\n"
	    "value, mean phase offset, and the fraction of samples that are\n"
	    "in phase.\n"
	    "\n"
	    "Unlike TSYNC, IPC does not scan over lags: it summarizes coupling\n"
	    "at the observed alignment of the two signals. Unlike SYNC, it is\n"
	    "a pairwise measure rather than a global synchrony index across\n"
	    "many channels." );

  add_param( "IPC" , "sig" , "C3,C4,F3,F4" , "Optionally specify channels (defaults to all)" );
  add_param( "IPC" , "sig1" , "C3" , "Alternatively specify seed channels" );
  add_param( "IPC" , "sig2" , "C4,F3,F4" , "Alternatively specify non-seed channels" );

  add_param( "IPC" , "add-channels" , "" , "Add channels (default pairwise, unless set to = 'seed'" );
  add_param( "IPC" , "prefix" , "" , "Label for new channels, defaults to IPC_" );
  
  add_table( "IPC" , "CH1,CH2" , "Primary IPC output" );
  add_var( "IPC" , "CH1,CH2" , "N_TOT" , "Total number of sample points considered" );
  add_var( "IPC" , "CH1,CH2" , "N_USED" , "Number of valid sample points used" );
  add_var( "IPC" , "CH1,CH2" , "IPC" , "Mean signed instantaneous phase coherence" );
  add_var( "IPC" , "CH1,CH2" , "WIPC" , "Weighted mean instantaneous phase coherence" );
  add_var( "IPC" , "CH1,CH2" , "PLV" , "Phase-locking value" );
  add_var( "IPC" , "CH1,CH2" , "PHASE" , "Circular mean phase difference" );
  add_var( "IPC" , "CH1,CH2" , "P_INPHASE" , "Fraction of samples with positive phase coherence" );

  //
  // CC
  //

  add_cmd( "cc" , "CC" , "Coupling and connectivity" );
  add_url( "CC" , "cc/#cc" );
  add_verb( "CC" ,
            "Estimate wavelet- or Hilbert-based coupling and connectivity metrics, "
            "including dPAC and wPLI.\n\n"
            "CC constructs analytic representations in one or two frequency bands "
            "and then evaluates either within-channel phase-amplitude coupling "
            "(dPAC) or cross-channel phase-lag connectivity (wPLI). It can scan "
            "across frequency grids, run permutation-based normalization, and emit "
            "both summary and epoch-level results." );

  add_param( "CC" , "sig" , "C3,C4,F3,F4" , "Optionally specify channels (defaults to all)" );
  add_param( "CC" , "pac" , "" , "Estimate within-channel phase-amplitude coupling metrics" );
  add_param( "CC" , "xch" , "" , "Estimate between-channel connectivity metrics" );
  add_param( "CC" , "nreps" , "1000" , "Number of replications" );
  add_param( "CC" , "fc" , "11,15" , "Wavelet center frequency/frequencies (phase)" );
  add_param( "CC" , "fwhm" , "1,1" , "Wavelet FWHM value(s) (phase)" );
  add_param( "CC" , "fc2" , "11,15" , "For PAC: as fc for amplitude" );
  add_param( "CC" , "fwhm2" , "1,1" , "For PAC: as fwhm for amplitude" );
  add_param( "CC" , "fc-range" , "1,20" , "Range of fc values" );
  add_param( "CC" , "num" , "10" , "Number of steps for fc-range" );
  add_param( "CC" , "linear" , "" , "Uniform ranged fc in linear space (versus log)" );
  add_param( "CC" , "fc2-range" , "20,40" , "Range of fc2 values" );
  add_param( "CC" , "fwhm-range" , "5,0.25" , "Range of fwhm values" );
  add_param( "CC" , "fwhm2-range" , "5,0.25" , "Range of fwhm2 values" );
  add_param( "CC" , "no-epoch-output" , "" , "Do not output epoch-level results" );

  add_table( "CC" , "CH1,CH2,F1,F2" , "Primary CC output" );
  add_var( "CC" , "CH1,CH2,F1,F2" , "CFC" , "Cross-frequency coupling 0/1" );
  add_var( "CC" , "CH1,CH2,F1,F2" , "XCH" , "Cross-channel coupling 0/1" );
  add_var( "CC" , "CH1,CH2,F1,F2" , "wPLI" , "Weighted phase lag index" );
  add_var( "CC" , "CH1,CH2,F1,F2" , "wPLI_Z" , "Z-normalized weighted phase lag index" );
  add_var( "CC" , "CH1,CH2,F1,F2" , "dPAC" , "dPAC metric" );
  add_var( "CC" , "CH1,CH2,F1,F2" , "dPAC_Z" , "Z-normalized dPAC metric" );

  add_table( "CC" , "E,CH1,CH2,F1,F2" , "Epoch-level CC output" );
  add_var( "CC" , "E,CH1,CH2,F1,F2" , "CFC" , "Cross-frequency coupling 0/1" );
  add_var( "CC" , "E,CH1,CH2,F1,F2" , "XCH" , "Cross-channel coupling 0/1" );
  add_var( "CC" , "E,CH1,CH2,F1,F2" , "wPLI" , "Weighted phase lag index" );
  add_var( "CC" , "E,CH1,CH2,F1,F2" , "wPLI_Z" , "Z-normalized weighted phase lag index" );
  add_var( "CC" , "E,CH1,CH2,F1,F2" , "dPAC" , "dPAC metric" );
  add_var( "CC" , "E,CH1,CH2,F1,F2" , "dPAC_Z" , "Z-normalized dPAC metric" );
  set_compressed( "CC" , tfac_t( "E,CH1,CH2,F1,F2" ) );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // PSC
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // --psc
  //

  add_cmd( "psc" , "--psc" , "Build a PSC projection model" );
  add_url( "--psc" , "psc/#build-psc" );
  add_verb( "--psc" ,
            "Constructs a principal spectral component model from one or more\n"
            "tabular spectral output files spanning many observations. The input\n"
            "matrix is built from selected variables, channels and frequencies,\n"
            "optionally filtered by ID lists, transformed by absolute-value or dB\n"
            "operations, restricted by frequency range, and standardized.\n"
            "\n"
            "Method-wise, --psc performs a PCA/SVD-style decomposition of the\n"
            "assembled feature matrix, optionally after iterative outlier removal.\n"
            "It emits component loadings and variance-explained summaries, and it\n"
            "can save a projection file that the runtime PSC command later applies\n"
            "to individual EDFs." );
  add_param( "--psc" , "spectra" , "psd1.txt,psd2.txt" , "Input spectral summary files" );
  add_param( "--psc" , "ch" , "C3,C4" , "Optionally restrict to these channels" );
  add_param( "--psc" , "v" , "PSD,COH" , "Variables to extract into the feature matrix" );
  add_param( "--psc" , "abs" , "PSI" , "Take absolute values of these variables before decomposition" );
  add_param( "--psc" , "dB" , "PSD" , "Log-transform these variables before decomposition" );
  add_param( "--psc" , "f-lwr" , "0.5" , "Lower frequency bound for included features" );
  add_param( "--psc" , "f-upr" , "30" , "Upper frequency bound for included features" );
  add_param( "--psc" , "epoch" , "" , "Treat rows as epoch-level observations using ID:E" );
  add_param( "--psc" , "signed-pairwise" , "" , "Treat paired CH1-CH2 features as signed" );
  add_param( "--psc" , "inc-ids" , "good.ids" , "Include only these IDs" );
  add_param( "--psc" , "ex-ids" , "bad.ids" , "Exclude these IDs" );
  add_param( "--psc" , "drop-incomplete-rows" , "" , "Drop observations with missing requested features" );
  add_param( "--psc" , "cmd-var" , "SPINDLES,DENS" , "Map variables back to command prefixes when naming features" );
  add_param( "--psc" , "proj" , "proj1.txt" , "Write the projection file to disk" );
  add_param( "--psc" , "dump" , "mat1.txt" , "Write component scores and original input features to disk" );
  add_param( "--psc" , "q" , "5" , "Optionally stratify dumped features by PSC quantiles" );
  add_param( "--psc" , "norm" , "" , "Standardize input features before decomposition" );
  add_param( "--psc" , "nc", "10" , "Number of PSCs to extract" );
  add_param( "--psc" , "th", "5,3" , "Iterative SD thresholds for outlier removal" );

  add_table( "--psc" , "I" , "Component-level PSC model summaries" );
  add_var( "--psc" , "I" , "W" , "Singular value for that component" );
  add_var( "--psc" , "I" , "VE" , "Variance explained by that component" );
  add_var( "--psc" , "I" , "CVE" , "Cumulative variance explained" );
  add_var( "--psc" , "I" , "INC" , "Whether that component was retained in the requested model" );

  add_table( "--psc" , "J" , "Feature metadata for PSC model variables" );
  add_var( "--psc" , "J" , "CH" , "Channel label for that feature" );
  add_var( "--psc" , "J" , "F" , "Frequency for that feature" );
  add_var( "--psc" , "J" , "VAR" , "Base variable name for that feature" );
  
  add_table( "--psc" , "I,J" , "PSC loading matrix" );
  add_var( "--psc" , "I,J" , "V" , "Loading of feature J on component I" );

  //
  // PSC
  //

  add_cmd( "psc" , "PSC" , "Project cached spectra into PSC space" );
  add_url( "PSC" , "psc/#project-psc" );
  add_verb( "PSC" ,
            "Applies a previously trained principal spectral component model to\n"
            "the current EDF. The model itself is created by the helper command\n"
            "--psc from a reference set of spectra; PSC then reads that saved\n"
            "projection, pulls the corresponding features from a cache, and\n"
            "projects the observation into the learned low-dimensional space.\n"
            "\n"
            "Method-wise, PSC is a PCA/SVD-style projection on a spectral feature\n"
            "matrix whose columns are typically channel-by-frequency summaries\n"
            "such as PSD, MTM, PSI or coherence-derived measures. At application\n"
            "time, the cached features are matched to the saved variable list,\n"
            "mean-centered using the reference-population means, optionally\n"
            "scaled by the reference SDs, multiplied by any saved per-feature\n"
            "scales, and then projected through the stored loading matrix.\n"
            "\n"
            "The runtime PSC command does not refit components and does not emit\n"
            "loadings or variance-explained tables. It only reports projected PSC\n"
            "scores for the current observation. The model-definition outputs such\n"
            "as component weights and variable loadings belong to --psc." );
  add_param( "PSC" , "cache" , "c1" , "Numeric cache holding the spectral features to project" );
  add_param( "PSC" , "proj" , "proj1.txt" , "Projection file created previously by --psc" );
  add_param( "PSC" , "nc" , "5" , "Optionally retain only the first k components from the projection file" );
  add_param( "PSC" , "drop" , "3,4" , "Zero out the listed 1-based components before projection" );
  add_param( "PSC" , "keep" , "1,2" , "Retain only the listed 1-based components before projection" );
  add_param( "PSC" , "norm" , "" , "Standardize cached features by the reference-population SDs before projection" );
  add_param( "PSC" , "signed" , "" , "For paired CH1-CH2 features, flip the sign if only the reversed pair is found" );
  add_param( "PSC" , "cmd-var" , "SPINDLES,DENS" , "Map feature names to CACHE command prefixes when resolving variables" );
  add_param( "PSC" , "clear" , "" , "Clear any attached projection before loading a new one" );

  add_table( "PSC" , "PSC" , "Projected principal spectral component scores" );
  add_var( "PSC" , "PSC" , "U" , "Projected PSC score for this observation" );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // TOPOGRAPHICAL ANALYSIS
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // CLOCS
  //

  add_cmd( "spatial" , "CLOCS" , "Attach channel location metadata" );
  add_url( "CLOCS" , "topo/#clocs" );
  add_verb( "CLOCS" ,
            "Attaches Cartesian channel-location metadata to the current EDF.\n"
            "If no file is provided, Luna loads its built-in default electrode\n"
            "location set. These coordinates are then used by downstream spatial\n"
            "commands such as SL and INTERPOLATE.\n"
            "\n"
            "CLOCS does not alter signal values itself and does not emit tabular\n"
            "writer output. It simply loads or initializes the montage geometry\n"
            "that other commands depend on." );
  add_param( "CLOCS" , "file" , "clocs.txt" , "Cartesian channel-location file" );
  add_param( "CLOCS" , "verbose" , "" , "Show additional diagnostics when loading a location file" );

  //
  // SL
  //

  add_cmd( "spatial" , "SL" , "Surface Laplacian transform" );
  add_url( "SL" , "topo/#sl" );
  add_verb( "SL" ,
            "Applies a spherical-spline surface Laplacian to a set of channels.\n"
            "This is a spatial high-pass transform that emphasizes local scalp\n"
            "gradients and attenuates broadly shared activity by combining all\n"
            "selected channels through their geometric arrangement on the head.\n"
            "\n"
            "Method-wise, Luna converts channel locations to a unit sphere,\n"
            "builds inter-electrode distance relationships, evaluates Legendre\n"
            "polynomials, and forms the standard G and H spline matrices. The\n"
            "resulting transform is then applied to the selected signals and\n"
            "written back in place." );
  add_param( "SL" , "sig" , "C3,C4,F3,F4" , "Signals to include in the spatial transform" );
  add_param( "SL" , "m" , "4" , "Spline flexibility parameter" );
  add_param( "SL" , "order" , "10" , "Maximum Legendre polynomial order" );
  add_param( "SL" , "lambda" , "1e-5" , "Regularization term added to the spline system" );

  //
  // INTERPOLATE
  //

  add_cmd( "spatial" , "INTERPOLATE" , "Interpolate bad channels using channel locations" );
  add_url( "INTERPOLATE" , "topo/#interpolate" );
  add_verb( "INTERPOLATE" ,
            "Interpolates bad channel-epoch pairs using attached channel locations\n"
            "and the current CHEP mask. For each epoch, channels marked bad by\n"
            "CHEP are reconstructed from the remaining good channels using the\n"
            "same spherical-spline interpolation framework that underlies the\n"
            "surface-Laplacian machinery.\n"
            "\n"
            "This command requires epoched data and an active CHEP mask. If too\n"
            "few good channels remain in an epoch, Luna marks the entire epoch\n"
            "instead of attempting interpolation. Successful interpolations are\n"
            "written back into the EDF and the corresponding CHEP flags are\n"
            "cleared." );
  add_param( "INTERPOLATE" , "sig" , "C3,C4,F3,F4" , "Signals eligible for interpolation" );
  add_param( "INTERPOLATE" , "m" , "4" , "Spline flexibility parameter" );
  add_param( "INTERPOLATE" , "order" , "10" , "Maximum Legendre polynomial order" );
  add_param( "INTERPOLATE" , "lambda" , "1e-5" , "Regularization term added to the interpolation system" );
  add_param( "INTERPOLATE" , "req-chs" , "8" , "Require at least this many good channels to attempt interpolation" );
  add_param( "INTERPOLATE" , "req-chs-prop" , "0.5" , "Require at least this proportion of channels to be good" );

  add_table( "INTERPOLATE" , "" , "Whole-run interpolation summary" );
  add_var( "INTERPOLATE" , "" , "NE_MASKED" , "Number of epochs masked because too few good channels were available" );
  add_var( "INTERPOLATE" , "" , "NE_NONE" , "Number of epochs with no bad channels requiring interpolation" );
  add_var( "INTERPOLATE" , "" , "NE_INTERPOLATED" , "Number of epochs in which at least one channel was interpolated" );
  add_var( "INTERPOLATE" , "" , "NCHEP_INTERPOLATED" , "Total number of bad channel-epoch pairs interpolated" );

  add_table( "INTERPOLATE" , "CH" , "Channel-level interpolation counts" );
  add_var( "INTERPOLATE" , "CH" , "NE_INTERPOLATED" , "Number of epochs in which that channel was interpolated" );
  add_var( "INTERPOLATE" , "CH" , "PCT_INTERPOLATED" , "Proportion of analyzed epochs in which that channel was interpolated" );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // MULTI-CHANNEL ANALYSIS
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // ICA
  //

  add_cmd( "multi" , "ICA" , "Independent component analysis" );
  add_url( "ICA" , "multi/#ica" );
  add_verb( "ICA" ,
            "Runs a whole-record independent component decomposition across a set\n"
            "of channels using a FastICA-style model. The selected signals are\n"
            "treated as a multivariate sample-by-channel matrix, whitened and\n"
            "decomposed into statistically independent source time courses plus\n"
            "their associated mixing and unmixing matrices.\n"
            "\n"
            "In practice, Luna can add the recovered IC time courses back into the\n"
            "EDF as new channels, write the decomposition matrices to disk, and\n"
            "export the A matrix in a simple text format for later use by ADJUST." );
  add_param( "ICA" , "sig" , "C3,C4,F3,F4" , "Signals to include in the ICA decomposition" );
  add_param( "ICA" , "A" , "A.txt" , "Required output file for the mixing matrix in text form" );
  add_param( "ICA" , "nc" , "10" , "Number of independent components to estimate" );
  add_param( "ICA" , "tag" , "IC_" , "Prefix for any newly added IC channels" );
  add_param( "ICA" , "file" , "ica_" , "Optional fileroot for matrix dumps" );
  add_param( "ICA" , "S" , "" , "With file=, also write the IC source matrix S" );
  add_param( "ICA" , "original-signals" , "" , "With file=, also write the original signal matrix X" );
  add_param( "ICA" , "no-new-channels" , "" , "Do not add recovered IC time courses as new EDF channels" );

  add_table( "ICA" , "IC1,IC2" , "ICA unmixing matrix" );
  add_var( "ICA" , "IC1,IC2" , "W" , "Unmixing matrix element" );

  add_table( "ICA" , "IC,CH" , "ICA mixing matrix" );
  add_var( "ICA" , "IC,CH" , "A" , "Mixing weight of a channel on that independent component" );

  add_table( "ICA" , "KCH,KIC" , "ICA whitening matrix" );
  add_var( "ICA" , "KCH,KIC" , "K" , "Whitening-matrix element" );

  //
  // ADJUST
  //

  add_cmd( "multi" , "ADJUST" , "Adjust signals using a multi-channel decomposition" );
  add_url( "ADJUST" , "multi/#adjust" );
  add_verb( "ADJUST" ,
            "Uses previously derived independent-component time courses and a saved\n"
            "A matrix to subtract selected components from one or more signals.\n"
            "This is a targeted back-projection workflow: component contributions\n"
            "are reconstructed for each channel and then removed from the original\n"
            "signals.\n"
            "\n"
            "Selection of components can be explicit, or driven by simple criteria\n"
            "such as extreme spatial weights in the A matrix or strong correlation\n"
            "with reference channels. The adjusted signals are written back in\n"
            "place, and Luna reports how similar each adjusted signal remains to\n"
            "its original version." );
  add_param( "ADJUST" , "A" , "A.txt" , "Required ICA mixing-matrix file produced by ICA" );
  add_param( "ADJUST" , "sig" , "C3,C4,F3,F4" , "Signals to adjust" );
  add_param( "ADJUST" , "adj" , "IC_1,IC_2" , "Independent-component channels to subtract from the signals" );
  add_param( "ADJUST" , "tag" , "IC_" , "Expected prefix for IC labels in the A matrix" );
  add_param( "ADJUST" , "force" , "IC_1,IC_2" , "Always include these components in the adjustment set" );
  add_param( "ADJUST" , "spatial" , "2.5" , "Select components whose maximum absolute spatial Z-score exceeds this threshold" );
  add_param( "ADJUST" , "corr-sig" , "ECG,EMG" , "Reference signals used to select components by correlation" );
  add_param( "ADJUST" , "corr-th" , "0.3,0.4" , "Absolute correlation thresholds paired with corr-sig" );

  add_table( "ADJUST" , "CH" , "Per-signal adjustment summaries" );
  add_var( "ADJUST" , "CH" , "R" , "Correlation between the original and adjusted signal" );

  //
  // SVD
  //

  add_cmd( "multi" , "SVD" , "Singular value decomposition" );
  add_url( "SVD" , "multi/#svd" );
  add_verb( "SVD" ,
            "Runs a whole-record singular value decomposition on a multichannel\n"
            "signal matrix. Signals are assembled into a sample-by-channel matrix,\n"
            "robustly centered, optionally variance-normalized and winsorized, and\n"
            "then decomposed into orthogonal components ordered by singular value.\n"
            "\n"
            "This is effectively a PCA-style channel decomposition of the current\n"
            "recording. Luna reports the singular spectrum and channel loadings,\n"
            "and it can add the leading component time courses back into the EDF\n"
            "as new channels." );
  add_param( "SVD" , "sig" , "C3,C4,F3,F4" , "Signals to include in the decomposition" );
  add_param( "SVD" , "nc" , "10" , "Number of leading components to retain as new channels" );
  add_param( "SVD" , "tag" , "U_" , "Prefix for any newly added component channels" );
  add_param( "SVD" , "norm" , "" , "Standardize each channel before decomposition" );
  add_param( "SVD" , "winsor" , "0.05" , "Winsorize each channel at this tail proportion before decomposition" );
  add_param( "SVD" , "no-new-channels" , "" , "Do not add component time courses as new EDF channels" );

  add_table( "SVD" , "C" , "Component-level SVD summaries" );
  add_var( "SVD" , "C" , "W" , "Singular value" );
  add_var( "SVD" , "C" , "VE" , "Variance explained by that component" );
  add_var( "SVD" , "C" , "CVE" , "Cumulative variance explained" );
  add_var( "SVD" , "C" , "INC" , "Whether that component was retained in the requested set" );

  add_table( "SVD" , "FTR,C" , "Channel loadings on SVD components" );
  add_var( "SVD" , "FTR,C" , "V" , "Loading of that channel on that component" );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // EEG MICROSTATE ANALYSIS
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // MS
  //

  add_cmd( "ms" , "MS" , "EEG microstate analysis" );
  add_verb( "MS" ,
            "Performs EEG microstate analysis from multichannel scalp topographies.\n"
            "In the basic single-record workflow, Luna computes global field power\n"
            "(GFP), selects GFP peaks, segments the peak maps with modified\n"
            "k-means across one or more candidate values of K, chooses an optimal\n"
            "microstate solution, backfits those prototype maps to all time\n"
            "points, optionally rejects ambiguous or very short segments, and then\n"
            "summarizes state coverage, duration, transitions and sequence\n"
            "complexity.\n"
            "\n"
            "The same command also supports a multi-sample workflow split into\n"
            "three stages: aggregate peaks from many EDFs, segment those pooled\n"
            "peaks into prototype maps, and then backfit saved prototypes to new\n"
            "records. It can also add annotations or state signals, cache state\n"
            "transition points, and compute k-mer enrichment statistics on the\n"
            "resulting microstate sequences." );
  add_param( "MS" , "sig" , "C3,C4,F3,F4" , "Channels used for microstate analysis" );
  add_param( "MS" , "k" , "4,5,6" , "Candidate numbers of microstate classes for segmentation" );
  add_param( "MS" , "peaks" , "peaks.edf" , "Multi-sample mode: aggregate GFP peaks into a pooled EDF" );
  add_param( "MS" , "segment" , "prototypes.txt" , "Segment a pooled peak EDF and write prototype maps" );
  add_param( "MS" , "backfit" , "prototypes.txt" , "Backfit previously saved prototype maps to the current record" );
  add_param( "MS" , "canonical" , "canon.txt" , "Map discovered states onto a canonical prototype set" );
  add_param( "MS" , "all-points" , "" , "Skip GFP peak selection and use all time points for segmentation" );
  add_param( "MS" , "standardize" , "" , "Standardize maps before modified k-means clustering" );
  add_param( "MS" , "gfp-max" , "3" , "Discard GFP peaks above mean + T*SD" );
  add_param( "MS" , "gfp-min" , "1" , "Discard GFP peaks below mean - T*SD" );
  add_param( "MS" , "gfp-kurt" , "3" , "Discard GFP peaks with unusually high map kurtosis" );
  add_param( "MS" , "npeaks" , "1000" , "Randomly retain at most this many GFP peaks" );
  add_param( "MS" , "lambda" , "0.1" , "Backfitting regularization parameter" );
  add_param( "MS" , "ambig" , "1.1,0.5" , "Mark samples ambiguous using relative and absolute SPC thresholds" );
  add_param( "MS" , "min-msec" , "20" , "Reject state segments shorter than this duration during smoothing" );
  add_param( "MS" , "epoch" , "" , "In backfit mode, analyze each epoch separately" );
  add_param( "MS" , "grouped" , "" , "Use grouped upper/lower-case state labels" );
  add_param( "MS" , "add-annot" , "MS_" , "Add annotation intervals for inferred microstates" );
  add_param( "MS" , "add-sig" , "MS_" , "Add binary state-indicator signals" );
  add_param( "MS" , "add-spc-sig" , "SPC_" , "Add state-wise spatial-correlation signals" );
  add_param( "MS" , "add-conf" , "CONF" , "With add-spc-sig, also add a confidence signal" );
  add_param( "MS" , "cache" , "c1" , "Cache sample points at which microstate transitions occur" );
  add_param( "MS" , "kmers" , "2,6,1000" , "Compute k-mer enrichment for lengths min..max with n permutations" );
  add_param( "MS" , "dump-gfp" , "gfp.txt" , "Dump the matrix used for GFP-peak segmentation" );
  add_param( "MS" , "gfp" , "states.txt" , "Dump final per-sample states and GFP values" );
  add_param( "MS" , "write-states" , "states.txt" , "Write the final state sequence to a file" );
  add_param( "MS" , "verbose" , "" , "Emit detailed peak-level output" );

  add_table( "MS" , "" , "Whole-run microstate summaries" );
  add_var( "MS" , "" , "GFP_MEAN" , "Mean GFP over candidate GFP peaks before thresholding" );
  add_var( "MS" , "" , "GFP_SD" , "SD of GFP over candidate GFP peaks before thresholding" );
  add_var( "MS" , "" , "NP0" , "Initial number of GFP peaks before threshold-based exclusions" );
  add_var( "MS" , "" , "NP" , "Final number of GFP peaks retained" );
  add_var( "MS" , "" , "OPT_K" , "Selected optimal number of microstate classes" );
  add_var( "MS" , "" , "AMBIG" , "Proportion of samples marked ambiguous" );
  add_var( "MS" , "" , "LZW" , "Lempel-Ziv complexity of the final state sequence" );
  add_var( "MS" , "" , "GEV" , "Total global explained variance of the final solution" );

  add_table( "MS" , "SP" , "Verbose GFP peak output" );
  add_var( "MS" , "SP" , "GFP" , "GFP at that retained peak sample" );

  add_table( "MS" , "CH,K" , "Optimal prototype maps" );
  add_var( "MS" , "CH,K" , "A" , "Prototype map weight for that channel and state" );

  add_table( "MS" , "K1,K2" , "Spatial correlations between optimal prototype maps" );
  add_var( "MS" , "K1,K2" , "SPC" , "Spatial correlation between two prototype maps" );

  add_table( "MS" , "KN,CH,K" , "Prototype maps for all tested K values" );
  add_var( "MS" , "KN,CH,K" , "A" , "Prototype map weight for that tested K solution" );

  add_table( "MS" , "NK" , "Fit summaries for all tested K values" );
  add_var( "MS" , "NK" , "R2" , "Explained variance for that K solution" );
  add_var( "MS" , "NK" , "MSE" , "Reconstruction MSE for that K solution" );
  add_var( "MS" , "NK" , "SIG2" , "Residual variance estimate for that K solution" );
  add_var( "MS" , "NK" , "SIG2_MCV" , "Modified cross-validation residual variance for that K solution" );

  add_table( "MS" , "K" , "Per-state microstate summaries" );
  add_var( "MS" , "K" , "GFP" , "Mean GFP for samples assigned to that state" );
  add_var( "MS" , "K" , "COV" , "Coverage of that state over all samples" );
  add_var( "MS" , "K" , "COV2" , "Coverage of that state over unambiguous samples" );
  add_var( "MS" , "K" , "OCC" , "Occurrence rate of that state over all samples" );
  add_var( "MS" , "K" , "OCC2" , "Occurrence rate of that state over unambiguous samples" );
  add_var( "MS" , "K" , "DUR" , "Mean duration of that state's contiguous runs" );
  add_var( "MS" , "K" , "WGT" , "Mean probabilistic coverage of that state across all samples" );
  add_var( "MS" , "K" , "SPC" , "Mean spatial correlation for samples assigned to that state" );
  add_var( "MS" , "K" , "GEV" , "Global explained variance attributable to that state" );
  add_var( "MS" , "K" , "N" , "Number of runs of that state" );

  add_table( "MS" , "PRE,POST" , "State transition probabilities" );
  add_var( "MS" , "PRE,POST" , "P" , "Probability of transitioning from PRE to POST" );

  add_table( "MS" , "L,S" , "Observed-sequence k-mer enrichment" );
  add_var( "MS" , "L,S" , "NG" , "Number of equivalent sequence-group members" );
  add_var( "MS" , "L,S" , "SG" , "Equivalent sequence-group label" );
  add_var( "MS" , "L,S" , "OBS" , "Observed count of that sequence" );
  add_var( "MS" , "L,S" , "EXP" , "Expected count of that sequence under permutation" );
  add_var( "MS" , "L,S" , "P" , "Empirical p-value for that sequence" );
  add_var( "MS" , "L,S" , "Z" , "Z-score for that sequence" );
  add_var( "MS" , "L,S" , "W_OBS" , "Observed count pooled over equivalent sequences" );
  add_var( "MS" , "L,S" , "W_EXP" , "Expected count pooled over equivalent sequences" );
  add_var( "MS" , "L,S" , "W_P" , "Empirical p-value pooled over equivalent sequences" );
  add_var( "MS" , "L,S" , "W_Z" , "Z-score pooled over equivalent sequences" );

  add_table( "MS" , "L,SG" , "Equivalent-group k-mer enrichment" );
  add_var( "MS" , "L,SG" , "NG" , "Number of sequences in that equivalent group" );
  add_var( "MS" , "L,SG" , "OBS" , "Observed count of that equivalent group" );
  add_var( "MS" , "L,SG" , "EXP" , "Expected count of that equivalent group" );
  add_var( "MS" , "L,SG" , "P" , "Empirical p-value for that equivalent group" );
  add_var( "MS" , "L,SG" , "Z" , "Z-score for that equivalent group" );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // CLUSTERING
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // EXE
  //

  add_cmd( "cluster" , "EXE" , "PDC-based epoch and channel clustering" );
  add_url( "EXE" , "exp/#exe" );
  add_verb( "EXE" ,
            "EXE builds an all-by-all distance matrix between observations using "
            "a permutation distribution complexity (PDC) encoding of the input "
            "signals. The default mode pools multiple channels to cluster epochs "
            "using a joint multichannel representation.\n\n"
            "With uni, EXE instead builds separate per-channel epoch-by-epoch "
            "matrices and reports channel-stratified solutions. With cat, it "
            "concatenates channel and epoch identity to form a larger observation "
            "set, allowing joint clustering of channel-epoch pairs; if the record "
            "is not epoched, cat falls back to clustering channels over the whole "
            "trace.\n\n"
            "The command can also write the raw distance matrix to disk, derive a "
            "set of representative exemplar epochs by recursive distance-based "
            "splitting, and optionally run hierarchical clustering with either a "
            "fixed number of clusters or a maximum cluster-size constraint." );

  add_param( "EXE" , "sig" , "C3,C4,F3,F4" , "Signals to include in the PDC distance matrix" );
  add_param( "EXE" , "uni" , "" , "Run separate univariate channel-wise analyses instead of one pooled multichannel analysis" );
  add_param( "EXE" , "cat" , "" , "Concatenate channels and epochs into one joint observation set" );
  add_param( "EXE" , "mat" , "dist.txt" , "Write the pairwise distance matrix to this file" );
  add_param( "EXE" , "sr" , "128" , "Resample signals to this sample rate before PDC encoding" );
  add_param( "EXE" , "entropy" , "" , "Use the entropy heuristic to choose the PDC embedding parameters" );
  add_param( "EXE" , "representative" , "4" , "Extract this many representative exemplar splits" );
    
  add_param( "EXE" , "m" , "5" , "PDC embedding dimension" );
  add_param( "EXE" , "t" , "1" , "PDC span" );

  add_param( "EXE" , "k" , "10" , "Target number of clusters for the final clustering solution" );
  add_param( "EXE" , "mx" , "20" , "Maximum allowed cluster size for the final clustering solution" );
    
  add_table( "EXE" , "E" , "Epoch-level outputs for pooled multichannel mode" );
  add_var( "EXE" , "E" , "CL" , "Assigned cluster label for that epoch [cluster]" );
  add_var( "EXE" , "E" , "K" , "Representative split label for that epoch [representative]" );
  add_var( "EXE" , "E" , "KE" , "Representative exemplar epoch for that epoch's split [representative]" );

  add_table( "EXE" , "CH,E" , "Channel-stratified outputs for uni or cat modes" );
  add_var( "EXE" , "CH,E" , "CL" , "Assigned cluster label for that channel/epoch observation [cluster]" );
  add_var( "EXE" , "CH,E" , "K" , "Representative split label for that channel/epoch observation [representative]" );
  add_var( "EXE" , "CH,E" , "KE" , "Representative exemplar epoch for that observation's split [representative]" );

  add_table( "EXE" , "K" , "Representative split summaries in pooled multichannel mode" );
  add_var( "EXE" , "K" , "E" , "Representative exemplar epoch for split K" );
  add_var( "EXE" , "K" , "N" , "Number of epochs assigned to split K" );

  add_table( "EXE" , "CH,K" , "Channel-stratified representative split summaries" );
  add_var( "EXE" , "CH,K" , "E" , "Representative exemplar epoch for split K" );
  add_var( "EXE" , "CH,K" , "N" , "Number of epochs assigned to split K" );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // ASSOCIATION
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // --gpa-prep
  //

  add_cmd( "assoc" , "--gpa-prep" , "Generic permutation-based association prep (--gpa-prep)" );
  add_url( "--gpa-prep" , "assoc/#gpa" );

  add_param( "--gpa-prep" , "dat" , "b.1" , "Write to this binary data file" );
  add_param( "--gpa-prep" , "specs" , "specs.json" , "Read from this JSON specification file" );
  add_param( "--gpa-prep" , "inputs" , "b.1" , "Read these input text files" );

  //
  // GPA
  //

  add_cmd( "assoc" , "GPA" , "Generic permutation-based association" );
  add_url( "GPA" , "assoc/#gpa" );
  add_verb( "GPA" ,
            "GPA fits linear-model associations between predictors and outcomes "
            "drawn from a manifest-driven binary data matrix, typically created by "
            "--gpa-prep. The main mode tests one or more Y variables against one "
            "or more X variables, with optional covariates Z, and can report "
            "either asymptotic p-values alone or permutation-based empirical "
            "significance.\n\n"
            "The command also supports several related workflows: filtering the "
            "loaded variable set by manifest properties, reporting descriptive "
            "statistics instead of regression, stratifying summaries by selected "
            "variables, and running comparison-style enrichment tests for selected "
            "predictors. Output therefore depends on mode: X,Y holds the main "
            "association results, X holds comparison summaries, X,N the null "
            "comparison distribution, and VAR or VAR,STRATUM the descriptive "
            "statistics." );

  add_param( "GPA" , "dat" , "b.1" , "Read from this binary data file (created by --gpa-prep)" );
  add_param( "GPA" , "vars" , "A,B,C" , "Include these variables" );
  add_param( "GPA" , "xvars" , "A,B,C" , "Exclude these variables" );
  add_param( "GPA" , "X" , "TST" , "Predictor variable(s)" );
  add_param( "GPA" , "Z" , "AGE,SEX" , "Covariates (nuisance variables)" );
  add_param( "GPA" , "Y" , "V1,V2" , "Explicitly set dependent variables" );
  add_param( "GPA" , "Xg" , "spindles,slow" , "Set predictor variables by group" );
  add_param( "GPA" , "Zg" , "demo,clinical" , "Set covariates by group" );
  add_param( "GPA" , "Yg" , "spindles,slow" , "Set dependent variables by group" );
  add_param( "GPA" , "nreps" , "1000" , "Number of permutations" );
  add_param( "GPA" , "adj" , "" , "Run all adjusted p-value corrections" );
  add_param( "GPA" , "bonf" , "" , "Add Bonferroni-adjusted p-values" );
  add_param( "GPA" , "holm" , "" , "Add Holm-adjusted p-values" );
  add_param( "GPA" , "fdr-by" , "" , "Add FDR(B&Y)-adjusted p-values" );
  add_param( "GPA" , "fdr" , "F" , "Turn off default FDR(B&H) adjusted p-values" );
  add_param( "GPA" , "dump" , "" , "Dump the data matrix to stdout" );
  add_param( "GPA" , "manifest" , "" , "Dump the variable manifest to stdout" );
  add_param( "GPA" , "summary" , "" , "Show a summary of the loaded data" );
  add_param( "GPA" , "summarize" , "" , "Alias for summary" );
  add_param( "GPA" , "desc" , "" , "Show a summary of the loaded data" );
  add_param( "GPA" , "stats" , "" , "Report descriptive statistics instead of association tests" );
  add_param( "GPA" , "strata" , "SEX" , "When used with stats, stratify summaries by these variables" );
  add_param( "GPA" , "verbose" , "" , "Verbose output" );
  add_param( "GPA" , "X-factors" , "" , "Add X-variable manifest details to the X,Y output" );
  add_param( "GPA" , "facs" , "F,CH" , "Include only variables stratified by this set of factors" );
  add_param( "GPA" , "xfacs" , "F,CH" , "Exclude variables stratified by this set of factors" );
  add_param( "GPA" , "grps" , "spindles,slow" , "Include variables assigned to these groups" );
  add_param( "GPA" , "xgrps" , "spindles,slow" , "Exclude variables assigned to these groups" );
  add_param( "GPA" , "n-req" , "10", "Drop columns with fewer than this many non-missing values" );
  add_param( "GPA" , "n-prop" , "0.1" , "Drop columns with more than this proportion of missing values" );
  add_param( "GPA" , "retain-cols" , "" , "Retain columns with missing values or invariant values" );
  add_param( "GPA" , "faclvls" , "B/SIGMA|BETA,CH/CZ" , "Include only variables with these factor levels" );
  add_param( "GPA" , "xfaclvls" , "B/SIGMA|BETA,CH/CZ" , "Exclude variables with these factor levels" );
  add_param( "GPA" , "lvars" , "10,100-200,200" , "Include only these variables by manifest list order" );
  add_param( "GPA" , "xlvars" , "10,100-200,200" , "Exclude these variables by manifest list order" );
  add_param( "GPA" , "nvars" , "10,100-200,200" , "Include these variables by manifest number" );
  add_param( "GPA" , "xnvars" , "10,100-200,200" , "Exclude these variables by manifest number" );
  add_param( "GPA" , "retain-rows" , "" , "Retain rows with missing values" );
  add_param( "GPA" , "subset" , "+MALE" , "Include only individuals positive (>0) for these variables" );
  add_param( "GPA" , "inc-ids" , "id1,id2" , "Include only these individuals" );
  add_param( "GPA" , "ex-ids" , "id1,id2", "Exclude these individuals" );
  add_param( "GPA" , "all-by-all" , "" , "Set all X to be all Y" );
  add_param( "GPA" , "winsor" , "0.05" , "Winsorize all variables at this threshold" );
  add_param( "GPA" , "knn" , "10" , "Impute missing values with k-nearest neighbors" );
  add_param( "GPA" , "force-zeros" , "" , "Replace remaining missing values with zero" );
  add_param( "GPA" , "qc" , "F" , "Turn off QC checks (if set to F)" );
  add_param( "GPA" , "p" , "0.01" , "Only output results below this nominal significance threshold" );
  add_param( "GPA" , "padj" , "0.05" , "Only output results below this adjusted significance threshold" );
  add_param( "GPA" , "adj-all-X" , "" , "Adjust PADJ across all X variables rather than within each X" );
  add_param( "GPA" , "progress" , "F" , "Turn off permutation progress reporting (if set to F)" );
  add_param( "GPA" , "comp" , "" , "Run comparison-style enrichment tests for the selected X variables" );
  add_param( "GPA" , "comp-verbose" , "" , "Dump the null comparison statistics for each X variable" );
    
  add_table( "GPA" , "X,Y" , "Association results per predictor (X) and outcome (Y)" );
  add_var( "GPA" , "X,Y" , "B" , "Regression coefficient" );
  add_var( "GPA" , "X,Y" , "T" , "t-statistic" );
  add_var( "GPA" , "X,Y" , "N" , "Number of observations" );
  add_var( "GPA" , "X,Y" , "P" , "Asymptotic p-value" );
  add_var( "GPA" , "X,Y" , "P_FDR" , "FDR(B&H)-adjusted p-value" );
  add_var( "GPA" , "X,Y" , "P_FDR_BY" , "FDR(B&Y)-adjusted p-value" );
  add_var( "GPA" , "X,Y" , "P_HOLM" , "Holm-adjusted p-value" );
  add_var( "GPA" , "X,Y" , "P_BONF" , "Bonferroni-adjusted p-value" );
  add_var( "GPA" , "X,Y" , "EMP" , "Empirical p-value" );
  add_var( "GPA" , "X,Y" , "EMPADJ" , "Adjusted empirical p-value" );
  add_var( "GPA" , "X,Y" , "GROUP" , "Dependent-variable group from the manifest" );
  add_var( "GPA" , "X,Y" , "BASE" , "Dependent-variable base name from the manifest" );
  add_var( "GPA" , "X,Y" , "STRAT" , "Dependent-variable factor-level string" );
  add_var( "GPA" , "X,Y" , "XGROUP" , "Predictor-variable group from the manifest [X-factors]" );
  add_var( "GPA" , "X,Y" , "XBASE" , "Predictor-variable base name from the manifest [X-factors]" );
  add_var( "GPA" , "X,Y" , "XSTRAT" , "Predictor-variable factor-level string [X-factors]" );

  add_table( "GPA" , "X" , "Comparison-test outputs per predictor [comp]" );
  add_var( "GPA" , "X" , "MATCH_OBS" , "Observed comparison statistic" );
  add_var( "GPA" , "X" , "MATCH_EXP" , "Expected comparison statistic under the null" );
  add_var( "GPA" , "X" , "Z" , "Standardized comparison statistic" );
  add_var( "GPA" , "X" , "P" , "One-sided empirical p-value for the comparison test" );
  add_var( "GPA" , "X" , "P2" , "Two-sided empirical p-value for the comparison test" );

  add_table( "GPA" , "X,N" , "Null comparison statistics per predictor [comp-verbose]" );
  add_var( "GPA" , "X,N" , "S" , "Null comparison statistic for permutation N" );

  add_table( "GPA" , "VAR" , "Descriptive statistics per variable [stats]" );
  add_var( "GPA" , "VAR" , "NOBS" , "Number of non-missing observations" );
  add_var( "GPA" , "VAR" , "MEAN" , "Sample mean" );
  add_var( "GPA" , "VAR" , "SD" , "Sample standard deviation" );

  add_table( "GPA" , "VAR,STRATUM" , "Stratified descriptive statistics per variable [stats,strata]" );
  add_var( "GPA" , "VAR,STRATUM" , "NOBS" , "Number of non-missing observations in that stratum" );
  add_var( "GPA" , "VAR,STRATUM" , "MEAN" , "Sample mean in that stratum" );
  add_var( "GPA" , "VAR,STRATUM" , "SD" , "Sample standard deviation in that stratum" );

  //
  // CPT
  //

  add_cmd( "assoc" , "CPT" , "Cluster permutation testing" );
  add_url( "CPT" , "assoc/#cpt" );
  add_verb( "CPT" ,
            "CPT performs mass-univariate linear-model testing on a set of "
            "dependent variables and then applies cluster-based permutation "
            "correction over an adjacency graph defined in space, frequency, "
            "and/or time.\n\n"
            "The command loads an IV/covariate table and one or more DV tables, "
            "optionally filters IDs, channels, and numeric ranges, applies "
            "outlier removal or winsorization, and computes variable-level "
            "statistics. It can then form supra-threshold clusters using adjacent "
            "features and estimate corrected significance from permutations. "
            "Outputs include per-variable results, cluster-level summaries, and "
            "explicit cluster memberships." );

  add_param( "CPT" , "iv-file" , "demo.txt" , "Tab-delimited file containing the primary independent variable and covariates" );
  add_param( "CPT" , "iv" , "DIS" , "Primary independent variable (a column in iv-file)" );
  add_param( "CPT" , "covar" , "AGE,SEX" , "Covariates, coded numerically, as columns in iv-file" );
  add_param( "CPT" , "dv-file" , "spec.txt,psd.txt" , "One or more dependent-variable files in long format" );  
  add_param( "CPT" , "dv" , "DENS,AMP" , "Dependent variables to analyze" );
  add_param( "CPT" , "all-dvs" , "" , "Use all dependent variables from the DV files" );  
  add_param( "CPT" , "ch" , "C3,C4" , "Restrict analysis to these channels" );
  add_param( "CPT" , "th" , "5" , "SD units for individual-level DV outlier removal" );
  add_param( "CPT" , "winsor" , "0.02" , "Threshold for winsorization of DVs" );
  add_param( "CPT" , "clocs" , "clocs.txt" , "File containing channel location information" );
  add_param( "CPT" , "nreps" , "1000" , "Number of permutations to perform" );
  add_param( "CPT" , "th-spatial" , "0.5" , "Threshold for defining adjacent channels (Euclidean distance, 0 to 2)" );
  add_param( "CPT" , "th-freq" , "1" , "Threshold for defining adjacent frequencies (Hz)" );
  add_param( "CPT" , "th-time" , "0.5" , "Threshold for defining adjacent time points (seconds)" );
  add_param( "CPT" , "th-cluster" , "2" , "Absolute t-statistic threshold for cluster inclusion" );
  add_param( "CPT" , "no-clustering" , "" , "Skip cluster formation and report only variable-level results" );
  add_param( "CPT" , "no-clocs" , "" , "Ignore channel locations and do not use spatial adjacency" );
  add_param( "CPT" , "dB", "", "Take the log of selected DVs" );
  add_param( "CPT" , "abs", "", "Take the absolute value of selected DVs" );
  add_param( "CPT" , "f-lwr", "0.5", "Ignore values for frequencies below this threshold" );
  add_param( "CPT" , "f-upr", "25", "Ignore values for frequencies above this threshold" );
  add_param( "CPT" , "t-lwr", "0.0", "Ignore values for times below this threshold" );
  add_param( "CPT" , "t-upr", "1.0", "Ignore values for times above this threshold" );
  add_param( "CPT" , "complete-obs", "", "Flag an error instead of dropping incomplete rows" );
  add_param( "CPT" , "dv-mat", "dv.txt", "Write the dependent-variable analysis matrix to this file" );
  add_param( "CPT" , "all-mat", "all.txt", "Write the full design matrix to this file" );
  add_param( "CPT" , "ex-ids", "id001,id002", "Individual IDs to exclude" );
  add_param( "CPT" , "inc-ids", "@{include.txt}", "Individual IDs to include" );
  add_param( "CPT" , "1-sided", "", "Assume a one-sided test (B > 0)" );
  add_param( "CPT" , "one-sided", "", "Alias for 1-sided" );
  add_param( "CPT" , "dump-adj", "", "Verbose dump of the adjacency structure" );
  add_param( "CPT" , "verbose", "", "Verbose output" );
  
  add_table( "CPT" , "VAR" , "Variable-level association output" );
  add_var( "CPT" , "VAR", "CH" , "Channel name" );
  add_var( "CPT" , "VAR", "CH1" , "First channel for variables stratified by channel pairs" );
  add_var( "CPT" , "VAR", "CH2" , "Second channel for variables stratified by channel pairs" );  
  add_var( "CPT" , "VAR", "F" , "Frequency (Hz) for variables stratified by frequency" );
  add_var( "CPT" , "VAR", "T" , "Time for variables stratified by time" );    
  add_var( "CPT" , "VAR", "B" , "Regression coefficient" );
  add_var( "CPT" , "VAR", "STAT" , "t-statistic" );
  add_var( "CPT" , "VAR", "PU" , "Uncorrected empirical p-value" );
  add_var( "CPT" , "VAR", "PC" , "Family-wise corrected empirical p-value" );
  add_var( "CPT" , "VAR", "CLST" , "Cluster number K for clustered variables, else 0" );

  add_table( "CPT" , "K" , "Cluster-level permutation results" );
  add_var( "CPT" , "K" , "N" , "Number of variables in this cluster" );
  add_var( "CPT" , "K" , "P" , "Empirical significance value for this cluster" );
  add_var( "CPT" , "K" , "SEED" , "Seed variable for this cluster" );
  
  add_table( "CPT" , "K,M" , "Cluster membership output" );
  add_var( "CPT", "K,M" , "VAR" , "Variable name, i.e. member M of cluster K" );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // PREDICTION
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // PREDICT
  //

  add_cmd( "pred" , "PREDICT" , "Apply a prediction model" );
  add_verb( "PREDICT" ,
            "PREDICT applies a pre-specified linear prediction model to a single "
            "EDF. The model file defines the required features, how they are "
            "normalized, the regression coefficients, and optional extras such as "
            "k-nearest-neighbor imputation, bias correction, and softplus "
            "post-processing.\n\n"
            "Features are usually read from a named numeric cache, although the "
            "command can also run cache-free if matching Luna variables already "
            "exist. Missing values can be imputed from an attached reference "
            "dataset, and atypical features can optionally be re-imputed if they "
            "fall beyond a z-distance threshold. Output includes the final "
            "prediction, optional bias-corrected or observed values, and per-"
            "feature diagnostics showing raw values, normalized values, model "
            "coefficients, and any imputation decisions." );

  add_param( "PREDICT" , "model" , "model.txt" , "Prediction model definition file" );
  add_param( "PREDICT" , "cache" , "ftrs" , "Numeric cache holding the input features" );
  add_param( "PREDICT" , "dump-model" , "" , "Dump the parsed model definition to the console" );
  add_param( "PREDICT" , "data" , "train.txt" , "Override the model's attached reference dataset for kNN imputation" );
  add_param( "PREDICT" , "knn" , "10" , "Override the model's k-nearest-neighbor imputation setting" );
  add_param( "PREDICT" , "retain-data" , "" , "Retain previously loaded kNN reference data" );
  add_param( "PREDICT" , "th" , "3" , "Re-impute features whose kNN distance exceeds this threshold" );
  add_param( "PREDICT" , "drop" , "FTR1,FTR2" , "Force these model terms to be treated as missing" );

  add_table( "PREDICT" , "" , "Overall prediction summary" );
  add_var( "PREDICT" , "" , "NF" , "Total number of model features" );
  add_var( "PREDICT" , "" , "NF_OBS" , "Number of non-missing features observed before imputation" );
  add_var( "PREDICT" , "" , "OKAY" , "1 if the prediction could be made, else 0" );
  add_var( "PREDICT" , "" , "Y" , "Predicted value" );
  add_var( "PREDICT" , "" , "Y1" , "Bias-corrected predicted value, if defined by the model" );
  add_var( "PREDICT" , "" , "YOBS" , "Observed value, if supplied by the model" );
  add_var( "PREDICT" , "" , "DIFF" , "Predicted minus observed difference" );

  add_table( "PREDICT" , "FTR" , "Per-feature prediction diagnostics" );
  add_var( "PREDICT" , "FTR" , "X" , "Observed raw feature value, if present" );
  add_var( "PREDICT" , "FTR" , "Z" , "Normalized feature value after any imputation" );
  add_var( "PREDICT" , "FTR" , "D" , "kNN distance for this feature [kNN]" );
  add_var( "PREDICT" , "FTR" , "IMP" , "1 if this feature was initially imputed [kNN]" );
  add_var( "PREDICT" , "FTR" , "REIMP" , "1 if this feature was re-imputed for atypicality [kNN]" );
  add_var( "PREDICT" , "FTR" , "M" , "Model mean used for normalization" );
  add_var( "PREDICT" , "FTR" , "SD" , "Model standard deviation used for normalization" );
  add_var( "PREDICT" , "FTR" , "B" , "Model coefficient for this feature" );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // SIMULATION
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // SIMUL
  //

  add_cmd( "simul" , "SIMUL" , "Simulate a signal from spectral constraints" );
  add_verb( "SIMUL" ,
            "SIMUL generates synthetic signals by constructing a target power "
            "spectrum, assigning random phases, and inverting back to the time "
            "domain. The basic mode creates a new or replacement channel from a "
            "PSD defined either explicitly, by a 1/f functional form, or from an "
            "input file.\n\n"
            "Two alternate modes reuse existing data more directly. With cache, "
            "SIMUL rebuilds a new signal epoch by epoch from cached PSD values, "
            "typically from prior PSD output. With fft, it uses an existing signal "
            "as the template for FFT-based reseeding. Additional options can zero "
            "selected frequency ranges, inject pulse windows, add the simulated "
            "signal back into an existing channel, or emit the PSD used for "
            "simulation." );

  add_param( "SIMUL" , "sig" , "C3" , "Signal to create or update" );
  add_param( "SIMUL" , "sr" , "128" , "Sample rate for a newly created signal" );
  add_param( "SIMUL" , "add" , "" , "Add the simulated signal to an existing channel instead of replacing it" );
  add_param( "SIMUL" , "white" , "" , "Generate white noise" );
  add_param( "SIMUL" , "file" , "psd.txt" , "Read the target PSD from a file with F and PSD columns" );
  add_param( "SIMUL" , "alpha" , "1" , "Use a 1/f^alpha spectral slope" );
  add_param( "SIMUL" , "intercept" , "1" , "Intercept term for the 1/f spectral model" );
  add_param( "SIMUL" , "frq" , "10,20" , "Center frequencies for explicitly specified peaks" );
  add_param( "SIMUL" , "psd" , "1,0.5" , "Power values for explicitly specified peaks" );
  add_param( "SIMUL" , "w" , "1" , "Gaussian peak width for frq/psd mode" );
  add_param( "SIMUL" , "pulses" , "5,1.0" , "Keep only N non-overlapping pulse windows of duration T seconds" );
  add_param( "SIMUL" , "impulse" , "10,1,20" , "Add one or more impulses as T,A,D triplets" );
  add_param( "SIMUL" , "zero" , "0.5,30" , "Zero out spectral power outside the specified frequency range" );
  add_param( "SIMUL" , "verbose" , "" , "Output the PSD used for simulation" );
  add_param( "SIMUL" , "cache" , "psd" , "Reconstruct a new signal epoch by epoch from cached PSD values" );
  add_param( "SIMUL" , "new" , "SIM" , "Output channel name for cache or fft reseeding modes" );
  add_param( "SIMUL" , "fft" , "" , "Use FFT reseeding mode based on an existing signal" );

  add_table( "SIMUL" , "F" , "Spectrum used for direct simulation [verbose]" );
  add_var( "SIMUL" , "F" , "LF" , "Natural log of frequency" );
  add_var( "SIMUL" , "F" , "P" , "Power spectral density at this frequency" );
  add_var( "SIMUL" , "F" , "LP" , "Natural log of power spectral density" );
  add_var( "SIMUL" , "F" , "DB" , "Power spectral density in dB" );

  add_table( "SIMUL" , "CH,E,F" , "Epoch-specific spectra used in cache mode [verbose]" );
  add_var( "SIMUL" , "CH,E,F" , "LF" , "Natural log of frequency" );
  add_var( "SIMUL" , "CH,E,F" , "P" , "Power spectral density at this frequency" );
  add_var( "SIMUL" , "CH,E,F" , "LP" , "Natural log of power spectral density" );
  add_var( "SIMUL" , "CH,E,F" , "DB" , "Power spectral density in dB" );

  //
  // SIGGEN
  //

  add_cmd( "simul" , "SIGGEN" , "Generate simple synthetic waveforms" );
  add_verb( "SIGGEN" ,
            "SIGGEN creates or updates a single signal channel using a small set "
            "of simple waveform generators. At present it supports sine waves and "
            "step-like impulses, either writing a new channel or adding the "
            "generated waveform to an existing one.\n\n"
            "This is a lighter-weight utility than SIMUL. It does not construct a "
            "signal from a target PSD or cached spectra; instead it directly draws "
            "the requested waveform in the time domain over the full record." );

  add_param( "SIGGEN" , "sig" , "SIM" , "Single signal to create or update" );
  add_param( "SIGGEN" , "sr" , "128" , "Sample rate for a newly created signal" );
  add_param( "SIGGEN" , "add" , "" , "Add the generated waveform to an existing channel instead of replacing it" );
  add_param( "SIGGEN" , "sine" , "10,1,0" , "Generate a sine wave as frequency, amplitude, and optional phase" );
  add_param( "SIGGEN" , "impulse" , "0.5,1,20" , "Add one or more impulses as T,A,D triplets" );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // HELPER UTILITIES
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // -h
  //

  add_cmd( "helpers" , "-h" , "Help functions" );

  //
  // --version
  //

  add_cmd( "helpers" , "--version" , "Show version (or -v)" );

  //
  // --build
  //

  add_cmd( "helpers" , "--build" , "Generate a sample list automatically" );
  add_verb( "--build" ,
            "--build scans one or more folders recursively and constructs a Luna "
            "sample list by pairing EDFs with associated annotation files. It is "
            "intended as a project bootstrap utility for assembling a dataset "
            "without manually writing the sample list.\n\n"
            "Options control how IDs are derived, which annotation extensions are "
            "recognized, whether similarly named files are matched across folders, "
            "and whether NSRR annotation naming conventions are assumed." );
  add_param( "--build" , "-nsrr" , "" , "Use NSRR annotation extension, i.e. `-ext=-nsrr.xml`" );
  add_param( "--build" , "-edfid" , "" , "Use filename as ID, instead of looking in each EDF header" );
  add_param( "--build" , "-nospan" , "" , "Do not match similarly-named files across folders" );
  add_param( "--build" , "-ext" , "-ext=txt,eannot,annot" , "Consider these extensions as annotation files" );

  //
  // --validate
  //

  add_cmd( "helpers" , "--validate" , "Validate all EDFs and annotation files in a project" );
  add_verb( "--validate" ,
            "--validate checks every EDF and attached annotation file referenced by "
            "a sample list. It is primarily a project hygiene tool for finding bad "
            "headers, missing files, unreadable annotations, and other dataset "
            "integrity problems before running larger analyses.\n\n"
            "Depending on options, it can also validate slice definitions and "
            "write lists of files or records that should be excluded." );
  add_param( "--validate" , "slist" , "slist=s.lst" , "Specficy the sample list" );
  add_param( "--validate" , "slice" , "" , "Also validate any slice definitions in the sample list" );
  add_param( "--validate" , "exclude-list" , "bad.lst" , "Write excluded records to this file" );
  add_param( "--validate" , "edf-exclude-list" , "bad-edf.lst" , "Write EDF-level exclusions to this file" );
  add_param( "--validate" , "annot-exclude-list" , "bad-annot.lst" , "Write annotation-level exclusions to this file" );

  add_table( "--validate" , "" , "Primary VALIDATE output" );
  add_var( "--validate" , "" , "EDF" , "Valid/invalid EDF (1=valid)" );
  add_var( "--validate" , "" , "ANNOTS" , "Valid/invalid annotations (1=valid)" );

  //
  // --repath
  //

  add_cmd( "helpers" , "--repath" , "Change root filepaths in a sample list" );
  add_verb( "--repath" ,
            "--repath rewrites paths inside an existing sample list by replacing "
            "one path prefix with another. This is useful when moving a project "
            "between machines or storage locations while keeping IDs and file "
            "names the same." );
  add_param( "--repath" , "{1st arg}" , "/home/john/" , "First argument: match string" );
  add_param( "--repath" , "{2nd arg}" , "/home/mary/" , "Second argument: replacement string" );

  //
  // --merge
  //

  add_cmd( "helpers" , "--merge" , "Merge (concatenate) multiple EDFs" );  
  add_verb( "--merge" ,
            "--merge concatenates two or more EDFs end to end to produce a single "
            "new EDF. It is intended for record-wise binding of compatible files, "
            "such as adjoining sessions from the same individual.\n\n"
            "The helper can take explicit EDF filenames or sample-list driven "
            "inputs, and lets you set the output EDF filename and the ID written "
            "into the new EDF header." );
  add_param( "--merge" , "slist" , "slist=s.lst" , "Specficy the sample list" );
  add_param( "--merge" , "id" , "id=id001" , "Specficy the ID in the new EDF header" );
  add_param( "--merge" , "edf" , "edf=m.edf" , "Filename for the resulting EDF (instead of merged.edf)" );
  add_param( "--merge" , "{edfs}" , "f1.edf f2.edf" , "Two or more EDFs" );

  //
  // --bind
  //

  add_cmd( "helpers" , "--bind" , "Merge (column/channel bind) multiple EDFs" );
  add_verb( "--bind" ,
            "--bind performs a column-wise merge of multiple EDFs, combining their "
            "channels into one output EDF rather than concatenating them in time. "
            "It is useful when different channel sets were recorded separately but "
            "need to be assembled into one file.\n\n"
            "As with --merge, the helper accepts either explicit EDF filenames or "
            "sample-list based inputs and lets you set the output filename and EDF "
            "ID." );
  add_param( "--bind" , "slist" , "slist=s.lst" , "Specficy the sample list" );
  add_param( "--bind" , "id" , "id=id001" , "Specficy the ID in the new EDF header" );
  add_param( "--bind" , "edf" , "edf=m.edf" , "Filename for the resulting EDF (instead of merged.edf)" );

  //
  // EXIT
  //

  add_cmd( "helpers" , "EXIT" , "Exit command processing" );

  //
  // REQUIRES
  //

  add_cmd( "helpers" , "REQUIRES" , "Require a minimum Luna version" );
  add_verb( "REQUIRES" ,
            "REQUIRES is a simple guard command that stops processing if the "
            "running Luna version is older than the version requested.\n\n"
            "It is mainly intended for scripts and pipelines that depend on "
            "newer command behavior and want to fail early with a clear message "
            "rather than run under an older binary." );
  add_param( "REQUIRES" , "version" , "1.3.4" , "Minimum required Luna version" );

  //
  // READ
  //

  add_cmd( "helpers" , "READ" , "Read signal data without further processing" );

  //
  // HRV
  //

  add_cmd( "physio" , "HRV" , "Heart-rate variability metrics" );
  add_verb( "HRV" ,
            "HRV derives beat-to-beat intervals from an ECG channel, computes "
            "time-domain and optionally frequency-domain heart-rate variability "
            "metrics, and can summarize them over the whole trace, by epoch, or "
            "within annotation-defined regions.\n\n"
            "The command first detects R-peaks, converts them to RR intervals, "
            "and then reports RR, heart-rate, SDNN, RMSSD, pNN50, and optional "
            "LF/HF spectral metrics. It can also emit R-peak and RR annotations, "
            "stratify by annotation class or instance, and tune the underlying "
            "R-peak detector and Welch spectrum settings." );
  add_param( "HRV" , "sig" , "ECG" , "ECG channel or channels to analyze" );
  add_param( "HRV" , "epoch" , "" , "Emit epoch-level HRV output when epochs are defined" );
  add_param( "HRV" , "annot" , "REM" , "Also summarize HRV within these annotation classes" );
  add_param( "HRV" , "by-instance" , "" , "When using annot, stratify by annotation instance as well" );
  add_param( "HRV" , "freq-domain" , "T" , "Enable frequency-domain HRV metrics" );
  add_param( "HRV" , "time-domain" , "T" , "Enable time-domain HRV metrics" );
  add_param( "HRV" , "lwr" , "0.3" , "Lower valid RR interval bound in seconds" );
  add_param( "HRV" , "upr" , "2" , "Upper valid RR interval bound in seconds" );
  add_param( "HRV" , "w" , "5" , "Median-filter width for RR interval cleaning" );
  add_param( "HRV" , "ns" , "512" , "Welch segment length for frequency-domain HRV" );
  add_param( "HRV" , "add-annot" , "Rpk" , "Add generic R-peak annotations" );
  add_param( "HRV" , "add-annot-rr" , "RRint" , "Add generic RR-interval annotations" );
  add_param( "HRV" , "add-annot-ch" , "Rpk" , "Add channel-specific R-peak annotations" );
  add_param( "HRV" , "add-annot-rr-ch" , "RRint" , "Add channel-specific RR-interval annotations" );
  add_param( "HRV" , "rp-lag" , "0.2" , "R-peak detector lag window in seconds" );
  add_param( "HRV" , "rp-infl" , "0.01" , "R-peak detector influence parameter" );
  add_param( "HRV" , "rp-th" , "3.5" , "Primary threshold for R-peak detection" );
  add_param( "HRV" , "rp-th2" , "1.5" , "Secondary threshold for R-peak detection" );
  add_param( "HRV" , "rp-max" , "2" , "Maximum allowed R-peak width" );
  add_param( "HRV" , "rp-dur" , "0.02" , "Minimum R-peak duration in seconds" );
  add_param( "HRV" , "rp-dur2" , "0.04" , "Secondary minimum R-peak duration in seconds" );
  add_param( "HRV" , "rp-ripple" , "0.02" , "Kaiser ripple parameter for ECG band-pass filtering" );
  add_param( "HRV" , "rp-tw" , "1" , "Transition width for ECG band-pass filtering" );
  add_param( "HRV" , "rp-flwr" , "5" , "Lower ECG band-pass edge for R-peak detection" );
  add_param( "HRV" , "rp-fupr" , "25" , "Upper ECG band-pass edge for R-peak detection" );
  add_param( "HRV" , "rp-w" , "0.02" , "Median-filter window for the R-peak detector" );

  add_table( "HRV" , "CH" , "Whole-trace or epoch-summarized HRV metrics" );
  add_var( "HRV" , "CH" , "IMPUTED" , "Number of imputed RR intervals" );
  add_var( "HRV" , "CH" , "P_INV" , "Estimated probability that the ECG is inverted" );
  add_var( "HRV" , "CH" , "INV" , "0/1 flag indicating ECG inversion" );
  add_var( "HRV" , "CH" , "NP" , "Number of retained R-peaks" );
  add_var( "HRV" , "CH" , "NP_TOT" , "Total number of R-peaks before exclusions or summarization" );
  add_var( "HRV" , "CH" , "RR" , "Mean RR interval" );
  add_var( "HRV" , "CH" , "HR" , "Mean heart rate" );
  add_var( "HRV" , "CH" , "SDNN" , "Standard deviation of NN intervals" );
  add_var( "HRV" , "CH" , "SDNN_R" , "Relative SDNN" );
  add_var( "HRV" , "CH" , "RMSSD" , "Root mean square successive RR difference" );
  add_var( "HRV" , "CH" , "RMSSD_R" , "Relative RMSSD" );
  add_var( "HRV" , "CH" , "pNN50" , "Proportion of successive RR differences exceeding 50 ms" );
  add_var( "HRV" , "CH" , "LF" , "Low-frequency HRV power" );
  add_var( "HRV" , "CH" , "HF" , "High-frequency HRV power" );
  add_var( "HRV" , "CH" , "LF_N" , "Normalized low-frequency power" );
  add_var( "HRV" , "CH" , "HF_N" , "Normalized high-frequency power" );
  add_var( "HRV" , "CH" , "LF_PK" , "Peak low-frequency HRV frequency" );
  add_var( "HRV" , "CH" , "HF_PK" , "Peak high-frequency HRV frequency" );
  add_var( "HRV" , "CH" , "LF2HF" , "Low-to-high frequency power ratio" );

  add_table( "HRV" , "CH,E" , "Epoch-level HRV metrics [epoch]" );
  add_var( "HRV" , "CH,E" , "IMPUTED" , "Number of imputed RR intervals" );
  add_var( "HRV" , "CH,E" , "P_INV" , "Estimated probability that the ECG is inverted" );
  add_var( "HRV" , "CH,E" , "INV" , "0/1 flag indicating ECG inversion" );
  add_var( "HRV" , "CH,E" , "NP" , "Number of retained R-peaks" );
  add_var( "HRV" , "CH,E" , "RR" , "Mean RR interval" );
  add_var( "HRV" , "CH,E" , "HR" , "Mean heart rate" );
  add_var( "HRV" , "CH,E" , "SDNN" , "Standard deviation of NN intervals" );
  add_var( "HRV" , "CH,E" , "SDNN_R" , "Relative SDNN" );
  add_var( "HRV" , "CH,E" , "RMSSD" , "Root mean square successive RR difference" );
  add_var( "HRV" , "CH,E" , "RMSSD_R" , "Relative RMSSD" );
  add_var( "HRV" , "CH,E" , "pNN50" , "Proportion of successive RR differences exceeding 50 ms" );
  add_var( "HRV" , "CH,E" , "LF" , "Low-frequency HRV power" );
  add_var( "HRV" , "CH,E" , "HF" , "High-frequency HRV power" );
  add_var( "HRV" , "CH,E" , "LF_N" , "Normalized low-frequency power" );
  add_var( "HRV" , "CH,E" , "HF_N" , "Normalized high-frequency power" );
  add_var( "HRV" , "CH,E" , "LF_PK" , "Peak low-frequency HRV frequency" );
  add_var( "HRV" , "CH,E" , "HF_PK" , "Peak high-frequency HRV frequency" );
  add_var( "HRV" , "CH,E" , "LF2HF" , "Low-to-high frequency power ratio" );

  add_table( "HRV" , "CH,ANNOT" , "Annotation-stratified HRV metrics [annot]" );
  add_var( "HRV" , "CH,ANNOT" , "NP" , "Number of retained R-peaks" );
  add_var( "HRV" , "CH,ANNOT" , "RR" , "Mean RR interval" );
  add_var( "HRV" , "CH,ANNOT" , "HR" , "Mean heart rate" );
  add_var( "HRV" , "CH,ANNOT" , "SDNN" , "Standard deviation of NN intervals" );
  add_var( "HRV" , "CH,ANNOT" , "SDNN_R" , "Relative SDNN" );
  add_var( "HRV" , "CH,ANNOT" , "RMSSD" , "Root mean square successive RR difference" );
  add_var( "HRV" , "CH,ANNOT" , "RMSSD_R" , "Relative RMSSD" );
  add_var( "HRV" , "CH,ANNOT" , "pNN50" , "Proportion of successive RR differences exceeding 50 ms" );

  add_table( "HRV" , "CH,ANNOT,INST" , "Annotation-instance-stratified HRV metrics [by-instance]" );
  add_var( "HRV" , "CH,ANNOT,INST" , "NP" , "Number of retained R-peaks" );
  add_var( "HRV" , "CH,ANNOT,INST" , "RR" , "Mean RR interval" );
  add_var( "HRV" , "CH,ANNOT,INST" , "HR" , "Mean heart rate" );
  add_var( "HRV" , "CH,ANNOT,INST" , "SDNN" , "Standard deviation of NN intervals" );
  add_var( "HRV" , "CH,ANNOT,INST" , "SDNN_R" , "Relative SDNN" );
  add_var( "HRV" , "CH,ANNOT,INST" , "RMSSD" , "Root mean square successive RR difference" );
  add_var( "HRV" , "CH,ANNOT,INST" , "RMSSD_R" , "Relative RMSSD" );
  add_var( "HRV" , "CH,ANNOT,INST" , "pNN50" , "Proportion of successive RR differences exceeding 50 ms" );

  //
  // RAI
  //

  add_cmd( "physio" , "RAI" , "REM atonia index" );
  add_verb( "RAI" ,
            "RAI computes the REM atonia index from a chin-EMG-like channel "
            "using 1-second epochs.\n\n"
            "For each epoch, Luna averages the rectified signal, subtracts a "
            "moving-minimum baseline, and then compares the residual amplitude "
            "to lower and upper thresholds. The REM atonia index is the "
            "proportion of epochs below the lower threshold relative to epochs "
            "below the lower threshold plus epochs above the upper threshold, "
            "excluding the intermediate region." );
  add_param( "RAI" , "sig" , "EMG" , "EMG channel or channels to analyze" );
  add_param( "RAI" , "th" , "1" , "Lower atonia threshold" );
  add_param( "RAI" , "th2" , "2" , "Upper exclusion threshold" );
  add_param( "RAI" , "verbose" , "" , "Emit per-epoch baseline-corrected amplitudes" );

  add_table( "RAI" , "CH" , "Channel-level REM atonia summaries" );
  add_var( "RAI" , "CH" , "REM_AI" , "REM atonia index" );
  add_var( "RAI" , "CH" , "NE" , "Number of epochs contributing to the index" );

  add_table( "RAI" , "CH,N" , "Per-epoch corrected amplitudes [verbose]" );
  add_var( "RAI" , "CH,N" , "X" , "Baseline-corrected mean rectified amplitude for that epoch" );

  //
  // AROUSALS
  //

  add_cmd( "physio" , "AROUSALS" , "Sleep arousal detection" );
  add_verb( "AROUSALS" ,
            "AROUSALS detects candidate arousals from EEG and optional EMG "
            "channels using short overlapping windows and a small set of "
            "derived spectral and complexity features.\n\n"
            "Luna first re-epochs the record into short overlapping windows, "
            "derives EEG and EMG feature summaries, and then classifies windows "
            "into baseline, artifact, micro-arousal, and arousal-like patterns. "
            "The command writes summary counts and feature means by class, adds "
            "annotation tracks for detected events, and can optionally create "
            "derived channels for the underlying features." );
  add_param( "AROUSALS" , "eeg" , "C3,C4" , "EEG channels used to derive arousal features" );
  add_param( "AROUSALS" , "emg" , "EMG" , "Optional EMG channels used to augment arousal detection" );
  add_param( "AROUSALS" , "win" , "1.0" , "Epoch window length in seconds for feature extraction" );
  add_param( "AROUSALS" , "inc" , "0.5" , "Epoch increment in seconds for feature extraction" );
  add_param( "AROUSALS" , "winsor" , "0.005" , "Winsorization fraction for outlier handling" );
  add_param( "AROUSALS" , "no-winsor" , "" , "Disable winsorization of feature values" );
  add_param( "AROUSALS" , "annot" , "l" , "Prefix for added arousal annotations" );
  add_param( "AROUSALS" , "add" , "a_" , "Add derived channels with this prefix" );
  add_param( "AROUSALS" , "prefix" , "ar_" , "Prefix for newly derived feature channels" );
  add_param( "AROUSALS" , "per-channel" , "T" , "Retain channel-specific feature metrics when adding channels" );

  add_table( "AROUSALS" , "" , "Overall arousal summary" );
  add_var( "AROUSALS" , "" , "MINS" , "Total analyzed duration in minutes" );
  add_var( "AROUSALS" , "" , "N" , "Number of detected arousals" );
  add_var( "AROUSALS" , "" , "AI" , "Arousal index per hour" );
  add_var( "AROUSALS" , "" , "DUR" , "Mean duration of detected arousals" );
  add_var( "AROUSALS" , "" , "N_MICRO" , "Number of detected micro-arousals" );
  add_var( "AROUSALS" , "" , "AI_MICRO" , "Micro-arousal index per hour" );
  add_var( "AROUSALS" , "" , "DUR_MICRO" , "Mean duration of detected micro-arousals" );
  add_var( "AROUSALS" , "" , "N_ART" , "Number of detected artifact events" );
  add_var( "AROUSALS" , "" , "AI_ART" , "Artifact-event index per hour" );

  add_table( "AROUSALS" , "CLS" , "Class-specific feature summaries" );
  add_var( "AROUSALS" , "CLS" , "NE" , "Number of windows assigned to this class" );
  add_var( "AROUSALS" , "CLS" , "PWR" , "Mean total-power feature for this class" );
  add_var( "AROUSALS" , "CLS" , "BETA" , "Mean relative beta feature for this class" );
  add_var( "AROUSALS" , "CLS" , "EMG" , "Mean EMG feature for this class" );
  add_var( "AROUSALS" , "CLS" , "SIGMA" , "Mean sigma-band feature for this class" );
  add_var( "AROUSALS" , "CLS" , "CMPLX" , "Mean complexity feature for this class" );

  //
  // OTSU
  //

  add_cmd( "helpers" , "OTSU" , "Calculate thresholds based on Otsu's method (internal channel)" );
  add_verb( "OTSU" ,
            "OTSU estimates one or more binary thresholds directly from one or "
            "more EDF channels using Otsu's method. For each signal, it evaluates "
            "candidate cut-points over the observed data range and reports the "
            "threshold that maximizes the between-class variance.\n\n"
            "This is the in-EDF version of the helper: it operates on channel data "
            "from the current record and reports both the empirical threshold and "
            "the full threshold scan if requested through the standard output "
            "tables." );
  add_param( "OTSU" , "sig" , "C3,C4" , "Signals on which to estimate Otsu thresholds" );
  add_param( "OTSU" , "k" , "100" , "Number of candidate threshold bins" );
  add_param( "OTSU" , "verbose" , "" , "Verbose flag accepted by the implementation" );

  add_table( "OTSU" , "CH" , "Otsu threshold summary by channel" );
  add_var( "OTSU" , "CH" , "EMPTH" , "Estimated Otsu threshold" );
  add_var( "OTSU" , "CH" , "EMPF" , "Empirical percentile at the estimated threshold" );

  add_table( "OTSU" , "CH,TH" , "Threshold scan across candidate cut-points" );
  add_var( "OTSU" , "CH,TH" , "SIGMAB" , "Between-class variance at that threshold" );
  add_var( "OTSU" , "CH,TH" , "F" , "Empirical percentile at that threshold" );

  //
  // --otsu
  //

  add_cmd( "helpers" , "--otsu" , "Calculate thresholds based on Otsu's method (external data)" );
  add_verb( "--otsu" ,
            "--otsu is the command-line form of Otsu thresholding. It reads a "
            "numeric vector from external input, evaluates candidate thresholds, "
            "and reports the cut-point maximizing between-class variance.\n\n"
            "Unlike OTSU, this helper does not operate on EDF channels. It is a "
            "general-purpose threshold finder for externally supplied values." );
  add_param( "--otsu" , "slist" , "slist=s.lst" , "Specficy the sample list" );
  add_param( "--otsu" , "k" , "100" , "Number of candidate threshold bins" );

  add_table( "--otsu" , "TH" , "Threshold scan across candidate cut-points" );
  add_var( "--otsu" , "TH" , "SIGMAB" , "Between-class variance at that threshold" );
  add_var( "--otsu" , "TH" , "F" , "Empirical percentile at that threshold" );

  /////////////////////////////////////////////////////////////////////////////////
  //
  // EXPERIMENTAL
  //
  /////////////////////////////////////////////////////////////////////////////////

//   add_cmd( "exp" , "ICA" , "Implementation of fastICA" );
//   add_url( "ICA" , "exp/#mi" );
  
//   add_param( "ICA" , "tag" , "comp" , "Add tag 'comp' to new channels, defaults to 'ICA'" );

//   add_param( "ICA" , "file" , "path/to/
//   std::string component_tag = param.has( "tag" ) ? param.value( "tag" ) : "ICA";

//   bool write_S_matrix = param.has( "file" );

//   std::string S_matrix_fileroot = write_S_matrix ? param.value( "file" ) : "xxx";

//   bool original_signals = param.has( "original-signals" );

//   bool do_not_add_channels = param.has( "no-new-channels" );

//   int compc = param.has( "compc" ) ? param.requires_int( "compc" ) : ns ;

  //
  // POL
  //

  add_cmd( "artifact" , "POL" , "Signal polarity diagnostics" , true );
  add_verb( "POL" ,
            "POL evaluates whether an EEG-like signal appears to have the expected "
            "polarity by comparing upward and downward half-waves after slow-band "
            "filtering. It is a diagnostic command intended to flag channels whose "
            "sign may be inverted or whose slow-wave morphology is atypical.\n\n"
            "The default method extracts thresholded half-wave segments from a "
            "band-passed signal and compares Hjorth-style and spectral summaries "
            "between the two polarities. Several variants change how half-waves "
            "are defined: zc2zc controls whether full zero-crossing-to-zero-"
            "crossing segments are used, mirror/double alter the alignment of "
            "upward and downward segments, raw bypasses the band-passed waveform "
            "for some summaries, d-mode switches to delta-oriented up/down "
            "analysis, and ht invokes the Hilbert-based polarity check." );

  add_param( "POL" , "sig" , "C3,C4" , "Signals to analyze" );
  add_param( "POL" , "th" , "1.0" , "SD threshold for selecting candidate half-waves" );
  add_param( "POL" , "not-zc2zc" , "" , "Do not force extraction to span zero crossing to zero crossing" );
  add_param( "POL" , "flim" , "5.0" , "Upper frequency limit for PSD-based summaries" );
  add_param( "POL" , "f-lwr" , "0.5" , "Lower band-pass frequency" );
  add_param( "POL" , "f-upr" , "4.0" , "Upper band-pass frequency" );
  add_param( "POL" , "not-mirror" , "" , "Do not mirror opposite-polarity segments before comparison" );
  add_param( "POL" , "double" , "" , "Use doubled upward segments instead of mirrored segments" );
  add_param( "POL" , "raw" , "" , "Use the raw signal rather than the band-passed signal for selected summaries" );
  add_param( "POL" , "d-mode" , "" , "Use delta up/down mode rather than positive/negative polarity mode" );
  add_param( "POL" , "ht" , "" , "Use the Hilbert-transform-based polarity check" );

  add_table( "POL" , "CH" , "Channel-level polarity summary" );
  add_var( "POL" , "CH" , "UP_TIME" , "Mean duration of upward half-waves" );
  add_var( "POL" , "CH" , "DOWN_TIME" , "Mean duration of downward half-waves" );
  add_var( "POL" , "CH" , "MN" , "Mean signal level summary" );
  add_var( "POL" , "CH" , "MD" , "Median signal level summary" );
  add_var( "POL" , "CH" , "T_DIFF" , "t-statistic comparing upward versus downward signal level" );
  add_var( "POL" , "CH" , "T_H1" , "t-statistic comparing Hjorth activity" );
  add_var( "POL" , "CH" , "T_H2" , "t-statistic comparing Hjorth mobility" );
  add_var( "POL" , "CH" , "T_H3" , "t-statistic comparing Hjorth complexity" );
  add_var( "POL" , "CH" , "UP_H1" , "Mean upward Hjorth activity [d-mode]" );
  add_var( "POL" , "CH" , "UP_H2" , "Mean upward Hjorth mobility [d-mode]" );
  add_var( "POL" , "CH" , "UP_H3" , "Mean upward Hjorth complexity [d-mode]" );
  add_var( "POL" , "CH" , "DOWN_H1" , "Mean downward Hjorth activity [d-mode]" );
  add_var( "POL" , "CH" , "DOWN_H2" , "Mean downward Hjorth mobility [d-mode]" );
  add_var( "POL" , "CH" , "DOWN_H3" , "Mean downward Hjorth complexity [d-mode]" );
  add_var( "POL" , "CH" , "N_H" , "Number of half-wave comparisons used for time-domain summaries" );
  add_var( "POL" , "CH" , "N_FFT" , "Number of half-wave comparisons used for spectral summaries" );

  add_table( "POL" , "CH,F" , "Frequency-specific polarity summary" );
  add_var( "POL" , "CH,F" , "T_PSD" , "t-statistic comparing absolute PSD [if emitted]" );
  add_var( "POL" , "CH,F" , "T_RELPSD" , "t-statistic comparing relative PSD" );
  add_var( "POL" , "CH,F" , "UP_PSD" , "Mean upward PSD [d-mode]" );
  add_var( "POL" , "CH,F" , "DOWN_PSD" , "Mean downward PSD [d-mode]" );
  add_var( "POL" , "CH,F" , "UP_RELPSD" , "Mean upward relative PSD [d-mode]" );
  add_var( "POL" , "CH,F" , "DOWN_RELPSD" , "Mean downward relative PSD [d-mode]" );

  // A.txt
  // 

  // K : cols x compc                                                                                                                                       
  // A : compc x compc                                                                                                                                      
  // W : compc x compc                                                                                                                                      

  // S : as original data                                                                                                                                   


  // ICA   Independent component analysis
  // EMD   Empirical mode decomposition  
  // ED    Diagnostic for electrical bridging
  // POL   Polarity check heuristic for sleep EEG
  // FIP   Frequency-interval plots

  // TSLIB Build library for SSS
  // SSS   Simple sleep stager
  // SLICE Short-time FFT for specified intervals

}



tfac_t::tfac_t( const strata_t & s )
{
  if ( s.levels.size() == 0 ) return;
  
  std::map<factor_t,level_t>::const_iterator ff = s.levels.begin();
  while ( ff != s.levels.end() )
    {
      if ( ff->first.factor_name[0] != '_' && ! globals::cmddefs().is_tag( ff->first.factor_name ) )
   fac.insert( ff->first.factor_name );
      ++ff;
    }
}


tfac_t::tfac_t( const std::string & s , const std::string & delim ) { 
  std::vector<std::string> tok = Helper::parse( s , delim );
  for (int i=0;i<tok.size();i++) 
    {
      if ( tok[i][0] != '_' && ! globals::cmddefs().is_tag( tok[i] ) )
   fac.insert( tok[i] );
    } 
}

std::string tfac_t::as_string( const std::string & delim ) const { 
  if ( fac.size() == 0 ) return "{baseline}";
  std::stringstream ss;
  std::set<std::string>::const_iterator ii = fac.begin();
  while ( ii != fac.end() ) {
    if ( ii != fac.begin() )
      ss << delim;
    ss << *ii;
    ++ii;
  }
  return ss.str();
}


bool tfac_t::operator< ( const tfac_t & rhs ) const { 
  if ( fac.size() < rhs.fac.size() ) return true;
  if ( fac.size() > rhs.fac.size() ) return false;
  std::set<std::string>::const_iterator ii = fac.begin();
  std::set<std::string>::const_iterator jj = rhs.fac.begin();
  while ( ii != fac.end() ) {
    if ( *ii < *jj ) return true;
    if ( *ii > *jj ) return false;      
    ++ii;
    ++jj;
  }
  return false;
} 


bool tfac_t::operator== ( const tfac_t & rhs ) const { 
  if ( fac.size() != rhs.fac.size() ) return false;
  std::set<std::string>::const_iterator ii = fac.begin();
  std::set<std::string>::const_iterator jj = rhs.fac.begin();
  while ( ii != fac.end() ) {
    if ( *ii != *jj ) return false;
    ++ii;
    ++jj;
  }
  return true;
} 


//
// Help display commands
//

std::string cmddefs_t::help_domains( ) const
{
  std::stringstream ss;
  std::map<std::string,std::string>::const_iterator ii = domain_desc.begin();
  while ( ii != domain_desc.end() ) { 
    ss << std::left << std::setw( 10 ) << ii->first << " " 
       << std::setw( 28 ) << domain_label.find( ii->first )->second << "\n";
    //<< ii->second << "\n";
    ++ii;
  }
  return ss.str();
}

std::string cmddefs_t::help_domain( const std::string & d ) const
{
  std::map<std::string,std::string>::const_iterator ii = domain_desc.find( d );
  if ( ii != domain_desc.end() ) return ii->second;
  return "";
}

bool cmddefs_t::check( const std::string & cmd ) const
{
  return cmds.find( cmd ) != cmds.end() ;
}

bool cmddefs_t::check( const std::string & cmd , const std::set<std::string> & k , std::set<std::string> * unknown ) const
{

  if ( k.size() == 0 ) return true;
  
  if ( cmds.find( cmd ) == cmds.end() ) return false;

  if ( pdesc.find( cmd ) == pdesc.end() ) 
    {
      *unknown = k;
      return false;
    }

  bool okay = true;

  const std::map<std::string,std::string> & p = pdesc.find( cmd )->second ;


  // check if command allows wild-cards, i.e. TAG 
  // this is specified by having "" as a registered parameter. 
  // in this case, we let anything go
  
  if ( p.find( "" ) != p.end() ) return true;
  
  // otherwise, explicitly check, one-by-one
  
  std::set<std::string>::const_iterator kk = k.begin();
  
  while ( kk != k.end() )
    {
      if ( p.find( *kk ) == p.end() ) 
   {
     unknown->insert( *kk );
     okay = false; 
   }
      ++kk;
    }
  
  return okay;
}


// all commands
std::string cmddefs_t::help_commands() const
{
  std::stringstream ss;
  std::map<std::string,std::set<std::string> >::const_iterator ii = dcmds.begin();
  while ( ii != dcmds.end() ) { 
    const std::set<std::string> & dc = ii->second;
    std::set<std::string>::const_iterator jj = dc.begin();
    while ( jj != dc.end() ) {      
      ss << help( *jj , true , false ) ;
      ++jj;
    }
    ss << "\n";
    ++ii;
  }
  return ss.str();  
}

// all commands in a domain
std::string cmddefs_t::help_commands( const std::string & d , const bool primary ) const 
{
  std::stringstream ss;
  std::map<std::string,std::set<std::string> >::const_iterator ii = dcmds.find( d );
  if ( ii == dcmds.end() ) return "";
  const std::set<std::string> & c = ii->second;
  std::set<std::string>::const_iterator jj = c.begin();
  while ( jj != c.end() ) {
    if ( ( ! primary ) || is_primary_cmd( *jj ) ) 
      ss << help( *jj , false , false ) ;
    ++jj;
  }
  return ss.str();    
}
  
// verbose describe commands [cmd] 

std::string cmddefs_t::help( const std::string & cmd , bool show_domain_label , bool verbose , bool primary ) const
{
  if ( cmds.find( cmd ) == cmds.end() ) return "";
  
  std::stringstream ss;
  if ( ! verbose ) 
    {
      if ( show_domain_label )
     ss << std::left << std::setw( 18 ) << domain_label.find( cdomain.find( cmd )->second )->second  << " " ;
      
      ss << std::left << std::setw( 12 ) << cmd << " "
	 << cmds.find( cmd )->second << "\n";    
    }
  else
    {
      ss << "\n";
      ss << cmd << " : " << cmds.find( cmd )->second 
	 << " (" << domain_label.find( cdomain.find( cmd )->second )->second  << ")\n"; 
      
      if ( curl.find( cmd ) != curl.end() ) 
	ss << std::string( cmd.size() , ' ' ) << " : " << url_root << curl.find( cmd )->second << "\n";

      if ( cverb.find( cmd ) != cverb.end() )
	ss << "\n" << wrap_text( cverb.find( cmd )->second , verbose_help_wrap ) << "\n";
      
      //
      // params
      //

      ss << "\nParameters:"
	 << "\n===========\n\n";

      std::map<std::string,std::map<std::string,std::string> >::const_iterator ii = pdesc.find( cmd );
      if ( ii == pdesc.end() ) ss << "   none\n";
      else
	{
	  std::map<std::string,std::string>::const_iterator jj = ii->second.begin();
	  while ( jj != ii->second.end() ) {
	    
	    if ( primary && ! is_primary_par( cmd , jj->first ) )
	      {
		++jj;
		continue;
	      }
	    
	    ss << "  " << std::left << std::setw( 24 ) << jj->first
	       << jj->second;
	    
	    // requirements?
	    std::string req = preq.find( cmd )->second.find( jj->first )->second ; 
	    if ( req != "" )
	      ss << " [req. " << req << "]";

	    // example (if any)
	    std::string ex = px.find( cmd )->second.find( jj->first )->second ; 
	    if ( ex != "" )
	      ss << "  e.g. " << jj->first << "=" << ex;
	    
	    ss << "\n";
	    
	    ++jj;
	  }
	  ++ii;
	}    


      //
      // outputs
      //

      ss << "\nOutputs:"
	 << "\n========\n\n";
      
      if ( otables.find( cmd ) == otables.end() )
	ss << "   none\n";
      else
	{
	  const std::map<tfac_t,std::string> & tab = otables.find( cmd )->second;
	  std::map<tfac_t,std::string>::const_iterator ii = tab.begin();
	  while ( ii != tab.end() )
	    {
	      const tfac_t & tfac = ii->first;

	      // skip if not primary?
	      if ( primary && ! is_primary_tbl( cmd , tfac ) )
		{
		  ++ii;
		  continue;
		}
	      
	      ss << "   " << std::left << std::setw( 24 ) << tfac.as_string( " x " ) 
		 << ii->second << "\n";
	      
	      ss << "   " << std::left << std::string( 60 , '-' )
		<< "\n";

	      // dump as compressed text ?
	      bool tdump = false;
	      if ( allz ) tdump = true;
	      else if ( nonez ) tdump = false;
	      else tdump = ofacs.find( cmd ) != ofacs.end() && ofacs.find( cmd )->second.find( tfac )->second ;
	      
	      if ( tdump ) ss << "   (compressed output)\n";
	      
	      // variables?
	      if ( ovars.find( cmd ) != ovars.end() )
		{
		  const std::map<tfac_t,std::map<std::string,std::string> > & t = ovars.find( cmd )->second;
		  if ( t.find( ii->first ) != t.end() ) 
		    {
		      const std::map<std::string,std::string> & v = t.find( ii->first )->second;
		      std::map<std::string,std::string>::const_iterator vv = v.begin();
		      while ( vv != v.end() ) {

			if ( primary && ! is_primary_var( cmd , tfac , vv->first ) ) 
			  {
			    ++vv;
			    continue;
			  }
			
			ss << "     " 
			   << std::left 
			   << std::setw(21) 
			   << vv->first << " " 
			   << vv->second << "\n";
			++vv;
		      }
		       
		    }
		}
	      
	      // spacer for next table
	      ss << "\n";
	      ++ii;
	    }
	  
	}

      // any notes? 
      std::map<std::string,std::string>::const_iterator nn = cnotes.find( cmd );
      if ( nn != cnotes.end() )
	ss << "\n" << nn->second << "\n";

    }

  
  return ss.str();
  
}



//
// test if a table exists
//

bool cmddefs_t::exists( const std::string & cmd ,
			const tfac_t & tfac ) const
{
  
  if ( cmds.find( cmd ) == cmds.end() )
      return false;

  if ( ofacs.find( cmd ) == ofacs.end() ) return false; 

  bool rv = ofacs.find( cmd )->second.find( tfac ) != ofacs.find( cmd )->second.end() ;

  return rv; 
}


//
// Output?
//

bool cmddefs_t::out_compressed( const std::string & cmd , const tfac_t & tfac ) const
{
  // all or none factors supplied?
  if ( allz ) return true;
  if ( nonez ) return false;

  // cmd not found [use ofacs instead of cmd, i.e. for commands w/ no registered tables, such as DESC
  //  if ( cmds.find( cmd ) == cmds.end() ) return false;
  
  // table not found 
  if ( ofacs.find( cmd ) == ofacs.end() ) return false;
  if ( ofacs.find( cmd )->second.find( tfac ) == ofacs.find( cmd )->second.end() ) return false;
  
  return ofacs.find( cmd )->second.find( tfac )->second ;
  
}


void cmddefs_t::set_compressed( const std::string & cmd , const tfac_t & tfac , const bool b )
{
  if ( cmds.find( cmd ) == cmds.end() ) return;
  if ( ofacs[ cmd ].find( tfac ) == ofacs[ cmd ].end() ) return;
  ofacs[ cmd ][ tfac ] = b;
}


void cmddefs_t::set_compressed( const std::string & cmd , const bool b )
{
  std::map<std::string,std::map<tfac_t,bool> >::iterator ii = ofacs.find( cmd );
  if ( ii == ofacs.end() ) return;

  std::map<tfac_t,bool>::iterator jj = ii->second.begin();
  while ( jj != ii->second.end() )
    {
      jj->second = b;
      ++jj;
    }
}


std::set<std::string> cmddefs_t::variables( const std::string & cmd ,  const param_t * param , const tfac_t & tfac  )
{
  // TODO: apply param-specified restriction on the variable list...
  //  --> REVISION: can now be done via REPORT hide/show mechanism, as added below
  std::set<std::string> r;

  // cmd is hidden
  if ( is_hidden_cmd( cmd ) ) return r;

  // cmd does not exist
  std::map<std::string,std::map<tfac_t,std::map<std::string,std::string> > >::const_iterator ii = ovars.find( cmd );
  if ( ii == ovars.end() ) return r;

  // table is hidden
  if ( is_hidden_table( cmd, tfac ) ) return r;

  // table does not exist
  const std::map<tfac_t,std::map<std::string,std::string> > & v2 = ii->second;
  std::map<tfac_t,std::map<std::string,std::string> >::const_iterator jj = v2.find( tfac );
  if ( jj == v2.end() ) return r;

  // hidden variables?  
  // std::map<std::string,std::map<tfac_t,std::map<std::string,bool> > > vhide;  // variables
  const std::map<std::string,bool> & hvars = vhide[ cmd ][ tfac ];
  
  // get variables
  const std::map<std::string,std::string> & v3 = jj->second;
  std::map<std::string,std::string>::const_iterator kk = v3.begin();
  while ( kk != v3.end() )
    {
      std::map<std::string,bool>::const_iterator hh = hvars.find( kk->first );
      if ( hh != hvars.end() && ! hh->second )  // i.e. vhide says 'show'
	r.insert( kk->first );
      ++kk;
    }
  return r;
}

//
// 'primary' commands (for -h listing) 
//

bool cmddefs_t::is_primary_cmd( const std::string & cmd ) const
{
  return pri_cmd.find( cmd ) != pri_cmd.end() ;
}

bool cmddefs_t::is_primary_par( const std::string & cmd , const std::string & param ) const
{
  std::map<std::string,std::set<std::string> >::const_iterator ii = pri_par.find( cmd );
  if ( ii == pri_par.end() ) return false;
  return ii->second.find( param ) != ii->second.end();
}

bool cmddefs_t::is_primary_tbl( const std::string & cmd , const tfac_t & tfac ) const
{
  std::map<std::string,std::set<tfac_t> >::const_iterator ii = pri_tbl.find( cmd );
  if ( ii == pri_tbl.end() ) return false;
  return ii->second.find( tfac ) != ii->second.end();
}

bool cmddefs_t::is_primary_var( const std::string & cmd , const tfac_t & tfac , const std::string & var ) const
{
  std::map<std::string,std::map<tfac_t,std::set<std::string> > >::const_iterator ii = pri_var.find( cmd );
  if ( ii == pri_var.end() ) return false;
  std::map<tfac_t,std::set<std::string> >::const_iterator jj = ii->second.find( tfac );
  if ( jj == ii->second.end() ) return false;
  return jj->second.find( var ) != jj->second.end();
}


// domain description
void cmddefs_t::add_domain( const std::string & domain , const std::string & label ,  const std::string & desc )
{
  if ( domain_label.find( domain ) == domain_label.end() )
    domain_ordered.push_back( domain );  
  domain_label[ domain ] = label;
  domain_desc[ domain ] = desc;
}

bool cmddefs_t::is_domain( const std::string & d ) 
{
  return domain_label.find( d ) != domain_label.end();
}

// command description 
void cmddefs_t::add_cmd1( const std::string & domain , const std::string & cmd , const std::string & desc , const bool hide )
{
  pri_cmd.insert( cmd );
  add_cmd( domain, cmd, desc, hide );
}

void cmddefs_t::add_cmd( const std::string & domain , const std::string & cmd , const std::string & desc , const bool hide )
{
  dcmds[ domain ].insert( cmd );
  cmds[ cmd ] = desc ; 
  cdomain[ cmd ] = domain;
  chide[ cmd ] = hide ;
}

// hidden command description 
void cmddefs_t::hidden_cmd( const std::string & domain , const std::string & cmd , const std::string & desc )
{
  add_cmd( domain, cmd , desc , true );    
}

bool cmddefs_t::is_cmd( const std::string & c ) 
{
  return cmds.find( c ) != cmds.end();
}


// command URLs , e.g.  zzz.bwh.harvard.edu/luna/ref/
void cmddefs_t::add_url( const std::string & cmd , const std::string & url ) 
{
  if ( cmds.find( cmd ) == cmds.end() ) Helper::halt( cmd + " not registered" );
  curl[ cmd ] = url;
}

void cmddefs_t::add_verb( const std::string & cmd , const std::string & verb )
{
  if ( cmds.find( cmd ) == cmds.end() ) Helper::halt( cmd + " not registered" );
  cverb[ cmd ] = verb;
}

void cmddefs_t::add_note( const std::string & cmd , const std::string & note ) 
{
  if ( cmds.find( cmd ) == cmds.end() ) Helper::halt( cmd + " not registered" );
  cnotes[ cmd ] = note;
}

// parameters for this command
void cmddefs_t::add_param1( const std::string & cmd , const std::string & param , 
			    const std::string & ex ,  // "" if none
			    const std::string & desc , 
			    const std::string & requirements  ,
			    const bool hide  )
{
  pri_par[ cmd ].insert( param );
  add_param( cmd , param , ex , desc , requirements , hide );
}

void cmddefs_t::add_param( const std::string & cmd , const std::string & param , 
			   const std::string & ex ,  // "" if none
			   const std::string & desc , 
			   const std::string & requirements  ,
			   const bool hide  )
{
  pdesc[ cmd ][ param ] = desc;
  preq[ cmd ][ param ] = requirements;
  px[ cmd ][ param ] = ex;
  phide[ cmd ][ param ] = hide;
}


// hide parameter for this command
void cmddefs_t::hidden_param( const std::string & cmd , const std::string & param , 
			    const std::string & ex ,  // "" if none
			    const std::string & desc , 
			    const std::string & requirements )
{
  add_param( cmd , param , ex , desc , requirements , true );
}


// output from this command , "CMD" , "F,B,CH,E" , "desc" , is compressed Y/N
void cmddefs_t::add_table1( const std::string & cmd , const std::string & factors , const std::string & desc , bool isz , bool hide )
{
  pri_tbl[ cmd ].insert( tfac_t( factors ) );
  add_table( cmd , factors , desc , isz , hide );
}

void cmddefs_t::add_table( const std::string & cmd , const std::string & factors , const std::string & desc , bool isz , bool hide )
{
  tfac_t tfac( factors );
  otables[ cmd ][ tfac ] = desc ; 
  ofacs[ cmd ][ tfac ] = isz ; 
  ohide[ cmd ][ tfac ] = hide;
}

void cmddefs_t::hidden_table( const std::string & cmd , const std::string & factors , const std::string & desc , bool isz )
{
  add_table( cmd , factors , desc , isz , true );
}

// ensure table exists
void cmddefs_t::ensure_table( const std::string & cmd , const std::string & factors )
{

  std::map<std::string,std::map<tfac_t,std::string> >::const_iterator ff = otables.find( cmd );
  // first need to ensure command exists
  if ( ff == otables.end() )
    {
      // will add a dummy domain '.' but that should not matter      
      add_cmd( "." , cmd , "." );
    }
  
  tfac_t tfac( factors );

  // only add if doesn't already exist
  if ( ff->second.find( tfac ) == ff->second.end() )
    {
      otables[ cmd ][ tfac ] = "." ;
      ofacs[ cmd ][ tfac ] = false ; 
      ohide[ cmd ][ tfac ] = false ;
    }

}


// add variable
void cmddefs_t::add_var1( const std::string & cmd , const std::string & factors , const std::string & var , const std::string & desc , const bool hide )
{
  pri_var[ cmd ][ tfac_t( factors ) ].insert( var );
  add_var( cmd , factors , var , desc , hide );
}

void cmddefs_t::add_var( const std::string & cmd , const std::string & factors , const std::string & var , const std::string & desc , const bool hide )
{
  tfac_t tfac( factors );
  ovars[ cmd ][ tfac ][ var ] = desc;
  vhide[ cmd ][ tfac ][ var ] = hide;
  otout[ cmd ][ tfac ][ var ] = true;  
}

// add hidden variable
void cmddefs_t::hidden_var( const std::string & cmd , const std::string & factors , const std::string & var , const std::string & desc )
{
  add_var( cmd , factors , var , desc , true );
}

// register (new) var for output
void cmddefs_t::register_var( const std::string & cmd , const std::string & factors , const std::string & var , const bool value )
{
  tfac_t tfac( factors );
  otout[ cmd ][ tfac ][ var ] = value;  
}


void cmddefs_t::all_compressed( bool b ) { allz = b; } 

bool cmddefs_t::all_compressed() const { return allz; }

void cmddefs_t::none_compressed( bool b ) { nonez = b; } 

bool cmddefs_t::none_compressed() const { return nonez; }


void cmddefs_t::add_tag( const std::string & tag ) { tags.insert( tag ); } 

void cmddefs_t::clear_tags() { tags.clear(); }

bool cmddefs_t::is_tag( const std::string & tag ) const { return tags.find( tag ) != tags.end(); } 



bool cmddefs_t::is_hidden_cmd( const std::string & c ) const
{
  std::map<std::string,bool>::const_iterator cc = chide.find( c );
  if ( cc == chide.end() ) return false;
  return cc->second;
}

bool cmddefs_t::is_hidden_param( const std::string & c , const std::string & p ) const
{
  std::map<std::string,std::map<std::string,bool> >::const_iterator cc = phide.find( c );
  if ( cc == phide.end() ) return false;
  std::map<std::string,bool>::const_iterator pp = cc->second.find( p );
  if ( pp == cc->second.end() ) return false;
  return pp->second;
}

bool cmddefs_t::is_hidden_table( const std::string & c , const tfac_t & tfac ) const
{
  std::map<std::string,std::map<tfac_t,bool> >::const_iterator cc = ohide.find( c );
  if ( cc == ohide.end() ) return false;
  std::map<tfac_t,bool>::const_iterator tt = cc->second.find( tfac );
  if ( tt == cc->second.end() ) return false;
  return tt->second;
}

bool cmddefs_t::is_hidden_var( const std::string & c , const tfac_t & tfac , const std::string & v ) const
{
  std::map<std::string,std::map<tfac_t,std::map<std::string,bool> > >::const_iterator cc = vhide.find( c );
  if ( cc == vhide.end() ) return false;
  std::map<tfac_t,std::map<std::string,bool> >::const_iterator tt = cc->second.find( tfac );
  if ( tt == cc->second.end() ) return false;
  std::map<std::string,bool>::const_iterator vv = tt->second.find( v );
  if ( vv == tt->second.end() ) return false;
  return vv->second;
}


// std::map<std::string,bool> chide;  // cmds
// std::map<std::string,std::map<tfac_t,bool> > ohide;  // tables
// std::map<std::string,std::map<tfac_t,std::map<std::string,bool> > > vhide;  // variables


void cmddefs_t::show_all( const bool show )
{
  std::map<std::string,bool>::iterator cc = chide.begin();
  while ( cc != chide.end() )
    {
      show_cmd( cc->first , show );
      ++cc;
    }  
}

void cmddefs_t::show_cmd( const std::string & cmd , const bool show )
{
  // don't worry whether this command is actually valid here
  chide[ cmd ] = ! show ;

  // hide/show all tales
  std::map<std::string,std::map<tfac_t,bool> >::iterator tt = ohide.find( cmd );
  if ( tt != ohide.end() )
    {
      std::map<tfac_t,bool> & table = tt->second;
      std::map<tfac_t,bool>::iterator qq = table.begin();
      while ( qq != table.end() )
   	{
   	  show_table( cmd , qq->first , show );
   	  ++qq;
   	}
    }
}

void cmddefs_t::show_table( const std::string & cmd , const tfac_t & factors , const bool show )
{
  ohide[ cmd ][ factors ] = ! show;
  
  // if showing this table, we need to make sure the cmd is shown too
  // if hiding this table, we can leave the cmd as is
  if ( show ) chide[ cmd ] = false;
    
  // and then also hide/show all variables in this table
  std::map<std::string,std::map<tfac_t,std::map<std::string,bool> > >::iterator it1 = vhide.find( cmd );
  if ( it1 == vhide.end() ) return;
  std::map<tfac_t,std::map<std::string,bool> >::iterator it2 = it1->second.find( factors );
  if ( it2 == it1->second.end() ) return;
  
  std::map<std::string,bool> & vars = it2->second;
  std::map<std::string,bool>::iterator vv = vars.begin();
  while ( vv != vars.end() )
    {
      show_var( cmd , factors , vv->first , show );
      ++vv;
    }
}

void cmddefs_t::show_table( const std::string & cmd , const std::string & factors , const bool show )
{
  show_table( cmd , tfac_t( factors ), show );
}


void cmddefs_t::show_var( const std::string & cmd , const tfac_t & factors , const std::string & var , const bool show )
{

  std::map<std::string,std::map<tfac_t,std::map<std::string,bool> > >::iterator it1 = vhide.find( cmd );
  if ( it1 == vhide.end() ) return;
  std::map<tfac_t,std::map<std::string,bool> >::iterator it2 = it1->second.find( factors );
  if ( it2 == it1->second.end() ) return;
  
  // always add
  it2->second[ var ] = ! show ;

  // if showing this variable, we need to make sure the cmd and table
  // are shown; but we don't want to run show_cmd() as that also recursively will 
  // set all tables/vars below; so just tweak the main vars here
  
  // if hiding this variable, we can leave as is

  if ( show ) 
    {
      chide[ cmd ] = false;
      ohide[ cmd ][ factors ] = false;
    }

}
  
void cmddefs_t::show_var( const std::string & cmd , const std::string & factors , const std::string & var , const bool show )
{
  show_var( cmd , tfac_t( factors ), var , show );
}


void cmddefs_t::hide_all()
{
  show_all( false );
}

void cmddefs_t::hide_cmd( const std::string & cmd )
{
  show_cmd( cmd , false );
}

void cmddefs_t::hide_table( const std::string & cmd , const std::string & factors )
{
  show_table( cmd, factors, false );
}

void cmddefs_t::hide_table( const std::string & cmd , const tfac_t & factors )
{
  show_table( cmd , factors , false );
}
  
void cmddefs_t::hide_var( const std::string & cmd , const std::string & factors , const std::string & var )
{
  show_var( cmd , factors , var , false );
}

void cmddefs_t::hide_var( const std::string & cmd , const tfac_t & factors , const std::string & var )
{
  show_var( cmd , factors , var , false );
}


//
// fetchers
//

std::vector<std::string> cmddefs_t::fetch_doms( const bool all ) const
{
  return domain_ordered;
}

std::vector<std::string> cmddefs_t::fetch_cmds( const std::string & dom , const bool all ) const
{
  std::vector<std::string> res;
  std::map<std::string,std::set<std::string> >::const_iterator ff = dcmds.find( dom );
  if ( ff != dcmds.end() )
    {
      const std::set<std::string> & s = ff->second;
      std::set<std::string>::const_iterator ss = s.begin();
      while ( ss != s.end() )
	{
	  res.push_back( *ss );
	  ++ss;
	}
    }  
  return res;
}

std::vector<std::string> cmddefs_t::fetch_params( const std::string & cmd, const bool all ) const
{
  std::map<std::string,std::map<std::string,std::string> >::const_iterator ff = pdesc.find( cmd );
  std::vector<std::string> res;
  if ( ff != pdesc.end() )
    {
      const std::map<std::string,std::string> & s = ff->second;
      std::map<std::string,std::string>::const_iterator ss = s.begin();
      while( ss != s.end() )
	{
	  res.push_back( ss->first );
	  ++ss;
	}
    }
  return res;  
}

std::vector<std::string> cmddefs_t::fetch_tbls( const std::string & cmd, const bool all ) const
{
  std::map<std::string,std::map<tfac_t,std::string> >::const_iterator ff = otables.find( cmd );
  std::vector<std::string> res;
  if ( ff != otables.end() )
    {
      const std::map<tfac_t,std::string> & s = ff->second;
      std::map<tfac_t,std::string>::const_iterator ss = s.begin();
      while ( ss != s.end() )
	{
	  res.push_back( ss->first.as_string() );
	  ++ss;
	}
    }
  return res;
}

std::vector<std::string> cmddefs_t::fetch_vars( const std::string & cmd, const std::string & tbl, const bool all ) const
{
  std::vector<std::string> res;
  std::map<std::string,std::map<tfac_t,std::map<std::string,std::string> > >::const_iterator ff = ovars.find( cmd );
  if ( ff == ovars.end() ) return res;
  std::map<tfac_t,std::map<std::string,std::string> >::const_iterator gg = ff->second.find( tfac_t( tbl ) );
  if ( gg == ff->second.end() ) return res;
  const std::map<std::string,std::string> & s = gg->second;
  std::map<std::string,std::string>::const_iterator ss = s.begin();
  while ( ss != s.end() )
    {
      res.push_back( ss->first );
      ++ss;
    }
  return res;
}


std::string cmddefs_t::fetch_desc_dom( const std::string & dom ) const
{
  std::map<std::string,std::string>::const_iterator ff = domain_desc.find( dom );
  return ff != domain_desc.end() ? ff->second : ".";
}

std::string cmddefs_t::fetch_label_dom( const std::string & dom ) const
{
  std::map<std::string,std::string>::const_iterator ff = domain_label.find( dom );
  return ff != domain_label.end() ? ff->second : ".";
}

std::string cmddefs_t::fetch_desc_cmd( const std::string & cmd ) const
{
  std::map<std::string,std::string>::const_iterator ff = cmds.find( cmd );
  return ff != cmds.end() ? ff->second : ".";
}

std::string cmddefs_t::fetch_desc_param( const std::string & cmd, const std::string & param ) const
{
  std::map<std::string,std::map<std::string,std::string> >::const_iterator ff = pdesc.find( cmd );
  if ( ff == pdesc.end() ) return "";
  const std::map<std::string,std::string> & s = ff->second;  
  std::map<std::string,std::string>::const_iterator ss = s.find( param );
  if ( ss == s.end() ) return "";
  return ss->second;
}

std::string cmddefs_t::fetch_desc_tbl( const std::string & cmd, const std::string & tbl ) const
{
  std::map<std::string,std::map<tfac_t,std::string> >::const_iterator ff = otables.find( cmd );
  if ( ff == otables.end() ) return "";
  const std::map<tfac_t,std::string> & s = ff->second;
  std::map<tfac_t,std::string>::const_iterator ss = s.find( tfac_t( tbl ) );
  if ( ss == s.end() ) return "";
  return ss->second;
}

std::string cmddefs_t::fetch_desc_var( const std::string & cmd, const std::string & tbl, const std::string & var ) const
{
  std::map<std::string,std::map<tfac_t,std::map<std::string,std::string> > >::const_iterator ff = ovars.find( cmd );
  if ( ff == ovars.end() ) return "";
  const	std::map<tfac_t,std::map<std::string,std::string> > & s = ff->second;
  std::map<tfac_t,std::map<std::string,std::string> >::const_iterator ss = s.find( tfac_t( tbl ) );
  if ( ss == s.end() ) return "";
  const std::map<std::string,std::string> & v = ss->second;
  std::map<std::string,std::string>::const_iterator vv = v.find( var );
  if ( vv == v.end() ) return "";
  return vv->second;
}



//#pragma GCC pop_options
