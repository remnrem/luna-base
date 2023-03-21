
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
#include <iomanip>

//#pragma GCC push_options
//#pragma GCC optimize ("O0")

#include "cmddefs.h"

extern globals global;


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

  url_root = "http://zzz.bwh.harvard.edu/luna/ref/";


  //
  // Domains
  //

  add_domain( "summ"       , "Summaries"       , "Basic summary commands" );
  add_domain( "annot"      , "Annotations"     , "Adding and displaying annotations" );
  add_domain( "expr"       , "Expressions"     , "Evaluating more advanced annotation-based expressions" );
  add_domain( "epoch"      , "Epochs"          , "Epoching signals and epoch-level annotations" );
  add_domain( "mask"       , "Masks"           , "Masking epochs based on annotations and other criteria" );
  add_domain( "manip"      , "Manipulations"   , "Manipulating signal data" );
  add_domain( "output"     , "Outputs"         , "Commands to output signals in different formats" );
  add_domain( "filter"     , "FIR filters"     , "FIR filter design and application" );
  add_domain( "artifact"   , "Artifacts"       , "Artifacts detection/correction routines" );
  add_domain( "hypno"      , "Hypnograms"      , "Characterizations of hypnograms" );
  add_domain( "staging"    , "Staging"         , "Automated staging/stage evaluation" );
  add_domain( "power"      , "Power spectra"   , "Power spectral density estimation" );
  add_domain( "transients" , "Spindles and SO" , "Spindles and slow oscillations" );
  add_domain( "topo"       , "Cross-signal"    , "Coherence and other topographical analyses" );
  add_domain( "cfc"        , "Cross-frequency" , "Phase-amplitude coupling" );
  add_domain( "misc"       , "Misc"            , "Misc. commands" );
  add_domain( "exp"        , "Experimental"    , "Experimental features: under heavy development, for internal use only" );
  add_domain( "cmdline"    , "Command-line"    , "Functions that do not operate on EDFs" );  


  /////////////////////////////////////////////////////////////////////////////////
  //
  // COMMAND-LINE OPTIONS
  //
  /////////////////////////////////////////////////////////////////////////////////
    
  add_cmd( "cmdline" , "--build" , "Scan folders recursively to geneate a sample list" );
  add_param( "--build" , "-edfid" , "" , "Use filename as ID, instead of looking in each EDF header" );
  add_param( "--build" , "-nospan" , "" , "Do not match similarly-named files across folders" );
  add_param( "--build" , "-ext" , "-ext=txt,eannot,annot" , "Consider these extensions as annotation files" );
  
  add_cmd( "cmdline" , "--xml" , "Dump annotations from an XML annotation file (to console)" );
  add_cmd( "cmdline" , "--xml2" , "Dump entire XML tree (to console)" );
  add_cmd( "cmdline" , "--eval" , "" );
  add_cmd( "cmdline" , "--pdlib" , "" );
  add_cmd( "cmdline" , "--fir" , " Or --fir-design" );
  add_cmd( "cmdline" , "--cwt" , "Or --cwt-design" );
  add_cmd( "cmdline" , "--eval-verbose" , "" );
  add_cmd( "cmdline" , "-h" , "Help functions" );
  add_cmd( "cmdline" , "--version" , "Show version (or -v)" );


  // -o
  // -t 
  // -a 
  // -s
  // @parameter file
  // x=y ... including special variables
  // ID
  // 1 2 


  /////////////////////////////////////////////////////////////////////////////////
  //
  // SUMMARIES
  //
  /////////////////////////////////////////////////////////////////////////////////


  //
  // DESC
  //

  add_cmd( "summ" , "DESC" , "Simple description of an EDF, sent to the console" );
  add_param( "DESC" , "channels" , "" , "Only write channel names, one-per-line" );

  
  //
  // SUMMARY
  //

  add_cmd( "summ" , "SUMMARY" , "More verbose description, sent to the console" );


  //
  // HEADERS
  //
  
  add_cmd( "summ" , "HEADERS" , "Tabulate (channel-specific) EDF header information" );
  
  add_table( "HEADERS" , "" , "Basic EDF header information" );
  add_var( "HEADERS" , "" , "NR" , "Number of records" );
  add_var( "HEADERS" , "" , "NS" , "Number of signals/channels" );
  add_var( "HEADERS" , "" , "EDF_ID" , "ID in the EDF header" );
  add_var( "HEADERS" , "" , "START_TIME" , "Start time in the EDF header" );
  add_var( "HEADERS" , "" , "STOP_TIME" , "Stop time" );
  add_var( "HEADERS" , "" , "START_DATE" , "Start date in the EDF header" );
  add_var( "HEADERS" , "" , "REC_DUR" , "Duration of each record (seconds)" );
  add_var( "HEADERS" , "" , "TOT_DUR_SEC" , "Total duration of EDF (seconds)" );
  add_var( "HEADERS" , "" , "TOT_DUR_HMS" , "Total duration of EDF (hh:mm:ss string)" );

  add_table( "HEADERS" , "CH" , "Per-channel header information" );
  add_var( "HEADERS" , "CH" , "DMAX" , "Digital max" );
  add_var( "HEADERS" , "CH" , "DMIN" , "Digital min" );
  add_var( "HEADERS" , "CH" , "PDIM", "Physical dimension" );
  add_var( "HEADERS" , "CH" , "PMAX", "Physical min" );
  add_var( "HEADERS" , "CH" , "PMIN", "Physical max" );
  add_var( "HEADERS" , "CH" , "SR", "Sample rate (Hz)" );
  add_var( "HEADERS" , "CH" , "SENS", "Sensitivity (unit/bit)" );
  add_var( "HEADERS" , "CH" , "TRANS", "Transducer type" );

  add_var( "HEADERS" , "CH" , "SENS", "Sensitivity (unit/bit)" );
  add_var( "HEADERS" , "CH" , "SENS", "Sensitivity (unit/bit)" );
  add_var( "HEADERS" , "CH" , "SENS", "Sensitivity (unit/bit)" );

  
  //
  // ALIASES
  //

  add_cmd( "summ" , "ALIASES" , "Tabulate channel and annotation alias replacements" );
  
  add_table( "ALIASES" , "CH" , "Channel aliasing" );
  add_var( "ALIASES" , "CH" , "ORIG" , "Original channel label in EDF" );

  add_table( "ALIASES" , "ANNOT" , "Annotation aliasing" );
  add_var( "ALIASES" , "ANNOT" , "ORIG" , "Original annotation label" );

  //
  // TAG
  //

  add_cmd( "summ" , "TAG" , "Generic command to add a tag (level/factor) to the output" );  
  add_param( "TAG" , ""    , "RUN/L1" , "Add tag with level L1 to factor RUN in output" );
  add_param( "TAG" , "tag" , "RUN/L1" , "Identical to the above, but explicitly using the tag option" );


  //
  // STATS
  //

  add_cmd( "summ"   , "STATS" , "Basic signal statistics (min/max, mean, RMS, etc)" );
  add_param( "STATS" , "sig" , "C3,C4" , "Restrict analysis to these channels" );
  add_param( "STATS" , "epoch" , "" , "Calculate per-epoch statistics" );
  
  add_table( "STATS" , "CH" , "Whole-night, per-channel statistics, based on all epochs" );
  add_var( "STATS" , "CH" , "MIN" , "Signal minimum (from data, not EDF header)" );
  add_var( "STATS" , "CH" , "MAX" , "Signal maximum (from data, not EDF header)" );
  add_var( "STATS" , "CH" , "MEAN" , "Signal mean" );
  add_var( "STATS" , "CH" , "MEDIAN" , "Signal median" );
  add_var( "STATS" , "CH" , "RMS" , "Signal root mean square" );

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
  add_var( "STATS" , "CH" , "MEDIAN_MEAN" , "Median of all per-epoch means [epoch]" );
  add_var( "STATS" , "CH" , "MEDIAN_MEDIAN" , "Median of all per-epoch medians [epoch]" );
  add_var( "STATS" , "CH" , "MEDIAN_RMS" , "Median of all per-epoch RMS [epoch]" );

  add_table( "STATS" , "CH,E" , "Per-epoch, per-channel statistics for unmasked epochs only" );
  add_var( "STATS" , "CH,E" , "MIN" , "Signal minimum (from data, not EDF header)" );
  add_var( "STATS" , "CH,E" , "MAX" , "Signal maximum (from data, not EDF header)" );
  add_var( "STATS" , "CH,E" , "MEAN" , "Signal mean" );
  add_var( "STATS" , "CH,E" , "MEDIAN" , "Signal median" );
  add_var( "STATS" , "CH,E" , "RMS" , "Signal root mean square" );

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
  // TLOCK  
  //

  add_cmd( "summ"   , "TLOCK" , "Time-locked signal summaries" );
  add_param( "TLOCK" , "sig" , "C3,C4" , "Restrict analysis to these channels" );
  add_param( "TLOCK" , "cache" , "" , "Sample-point peak cache (req.)" );

  add_param( "TLOCK" , "tolog" , "" , "Take log of input signals" );
  add_param( "TLOCK" , "verbose" , "" , "Output individual intervals" );

  add_param( "TLOCK" , "w" , "2" , "Window size around peaks (seconds) (req.)");
  add_param( "TLOCK" , "phase" , "20" , "Expect phase values (radians) and summarize in, e.g. 20 bins" );
  
  // spindle peaks
  add_table( "TLOCK" , "CH,sCH,sF" , "Spindle peak time-locked counts" );
  add_var( "TLOCK" , "CH,sCH,sF" , "N" , "Number of included peaks" );
  add_var( "TLOCK" , "CH,sCH,sF" , "N_ALL" , "Total number of included peaks" );

  add_table( "TLOCK" , "CH,SEC,sCH,sF" , "Spindle peak time-locked summaries" );
  add_var( "TLOCK" , "CH,SEC,sCH,sF" , "M" , "Signal mean" );

  add_table( "TLOCK" , "N,SEC,CH,sCH,sF" , "Spindle peak time-locked counts" );
  add_var( "TLOCK" , "N,SEC,CH,sCH,sF" , "V" , "Signal value" );
  set_compressed( "TLOCK" , tfac_t( "N,SEC,CH,sCH,sF" ) );


  //
  // PEAKS
  //

  add_cmd( "misc"   , "PEAKS" , "Peak finder (maxima)" );
  add_param( "PEAKS" , "sig" , "C3,C4" , "Restrict analysis to these channels" );
  add_param( "PEAKS" , "cache" , "c1" , "Sample-point peak cache (req.)" );

  add_param( "PEAKS" , "epoch" , "" , "Also find minima" );
  add_param( "PEAKS" , "clipped" , "3" , "Clipped regions = 3 consecutive points (0 = ignore clipping)" );
  add_param( "PEAKS" , "min" , "" , "Also find minima" );
  add_param( "PEAKS" , "min-only" , "" , "Only find minima" );
  add_param( "PEAKS" , "percentile" , "20" , "Only report top 20% of peaks" );

    
  /////////////////////////////////////////////////////////////////////////////////
  //
  // ANNOTATIONS
  //
  /////////////////////////////////////////////////////////////////////////////////
  
  //
  // --xml
  //

  add_cmd( "annot" , "--xml" , "Quickly view an NSRR XML annotation file" );
  add_note( "--xml" , "Command line option\n  luna --xml file.xml" );

  //
  // --xml2
  //

  add_cmd( "annot" , "--xml2" , "Verbose dump of XML tree" );
  add_note( "--xml2" , "Command line option\n  luna --xml2 file.xml" );

  //
  // ANNOTS
  //
  
  add_cmd( "annot" , "ANNOTS" , "Tabulate all annotations" );
  add_url( "ANNOTS" , "annotations/#annots" );
  //  add_note( "ANNOTS" , "Any formatted free text goes here,\n shown at end of verbose help link\n");

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
  add_var( "ANNOTS" , "ANNOT,INST,T" , "START" , "Start time (seconds) of this instance" );
  add_var( "ANNOTS" , "ANNOT,INST,T" , "STOP" , "Stop time (seconds) of this instance" );
  add_var( "ANNOTS" , "ANNOT,INST,T" , "VAL" , "The meta-data for this instance, if any exists (otherwise missing NA)" );
  add_var( "ANNOTS" , "ANNOT,INST,T" , "ALL_MASKED" , "? [show-masked]" );
  add_var( "ANNOTS" , "ANNOT,INST,T" , "ALL_UNMASKED" , "? [show-masked]" );
  add_var( "ANNOTS" , "ANNOT,INST,T" , "SOME_MASKED" , "? [show-masked]" );
  add_var( "ANNOTS" , "ANNOT,INST,T" , "SOME_UNMASKED" , "? [show-masked]" );
  add_var( "ANNOTS" , "ANNOT,INST,T" , "START_MASKED" , "? [show-masked]" );

  add_table( "ANNOTS" , "E,INTERVAL,INST" , "Per-epoch instance-level annotation tabulation" );
  add_var( "ANNOTS" , "E,INTERVAL,INST" , "AMASK" , "Annotation instance mask status (1=masked/excluded) [epoch]" );
  add_var( "ANNOTS" , "E,INTERVAL,INST" , "EMASK" , "Epoch mask status (1=masked/excluded) [epoch]" );


  //
  // SPANNING
  //

  add_cmd( "annot" , "SPANNING" , "Report duration spanned or not by group of annotations" );
  add_url( "SPANNING" , "annotations/#spanning" );
  add_param( "SPANNING" , "annot" , "N1,N2,N3,R,W" , "Spanning annotation group" );

  add_table( "SPANNING" , "N" , "Invalid annotations" );
  add_var( "SPANNING" , "N" , "ANNOT" , "Annotation class" );
  add_var( "SPANNING" , "N" , "INST" , "Annotation instance" );
  add_var( "SPANNING" , "N" , "START" , "Start (seconds)" );
  add_var( "SPANNING" , "N" , "STOP" , "Stop (seconds)" );

  add_table( "SPANNING" , "" , "Spanning summary report" );
  add_var( "SPANNING" , "" , "REC_SEC" , "EDF recording duration (seconds)" );
  add_var( "SPANNING" , "" , "REC_HMS" , "EDF recording duration (hh:mm:ss)" );

  add_var( "SPANNING" , "" , "ANNOT_N" , "Number of annotations in group" );
  add_var( "SPANNING" , "" , "ANNOT_SEC" , "Total (potentially overlapping) annotation duration (secs)" );
  add_var( "SPANNING" , "" , "ANNOT_HMS" , "Total (potentially overlapping) annotation duration (hh:mm:ss)" );

  add_var( "SPANNING" , "" , "ANNOT_OVERLAP" , "Do any annotations in group overlap w/ one another (0/1)?" );

  add_var( "SPANNING" , "" , "INVALID_N" , "Number of annotations that over-extend EDF duration" );
  add_var( "SPANNING" , "" , "VALID_N" , "Number of valid annotations, ANNOT_N - INVALID_N" );

  add_var( "SPANNING" , "" , "INVALID_SEC" , "Total duration of all annotation beyond EDF end" );

  add_var( "SPANNING" , "" , "SPANNED_PCT" , "% of EDF spanned by 1+ of these annotations" );
  add_var( "SPANNING" , "" , "SPANNED_SEC" , "Duration of EDF spanned by 1+ of these annotations (secs)" );
  add_var( "SPANNING" , "" , "SPANNED_HMS" , "Duration of EDF spanned by 1+ of these annotations (hh:mm:ss)" );

  add_var( "SPANNING" , "" , "UNSPANNED_PCT" , "% of EDF unspanned by 1+ of these annotations" );
  add_var( "SPANNING" , "" , "UNSPANNED_SEC" , "Duration of EDF unspanned by 1+ of these annotations (secs)" );
  add_var( "SPANNING" , "" , "UNSPANNED_HMS" , "Duration of EDF unspanned by 1+ of these annotations (hh:mm:ss)" );


  /////////////////////////////////////////////////////////////////////////////////
  //
  // EPOCHS
  //
  /////////////////////////////////////////////////////////////////////////////////
  
  // EPOCH 

  add_cmd( "epoch" , "EPOCH" , "Set epochs" );
  add_url ( "EPOCH" , "epochs/#epoch" );

  add_param( "EPOCH" , "len" , "30" , "Epoch length (seconds), defaults to 30" );
  add_param( "EPOCH" , "dur" , "30" , "Same as len" );
  add_param( "EPOCH" , "inc" , "30" , "Epoch increment (seconds), defaults to len (i.e. no overlap)" );
  add_param( "EPOCH" , "epoch" , "30,15" , "Same as len=30 inc=15" );
  add_param( "EPOCH" , "require" , "10" , "Stop processing that EDF if there are not at least N epochs" );
  add_param( "EPOCH" , "verbose" , "" , "Output epoch-level information" );
  add_param( "EPOCH" , "clear" , "" , "Unepoch all signals" );

  add_table( "EPOCH" , "" , "Epoch-level summaries" );
  add_var( "EPOCH" , "" , "DUR" , "Epoch duration (seconds)" );
  add_var( "EPOCH" , "" , "INC" , "Epoch increment (seconds)" );
  add_var( "EPOCH" , "" , "NE" , "Number of epochs" );

  add_table( "EPOCH" , "E" , "Per-epoch interval information [verbose]" );
  add_var( "EPOCH" , "E" , "E1" , "Current epoch number (which may differ from E if the EDF has been restructured)" );
  add_var( "EPOCH" , "E" , "HMS" , "Clock-time for epoch start (hh:mm:ss)" );
  add_var( "EPOCH" , "E" , "INTERVAL" , "String label of epoch interval (seconds)" );
  add_var( "EPOCH" , "E" , "MID" , "Midpoint of epoch (seconds elapsed from EDF start)" );
  add_var( "EPOCH" , "E" , "START" , "Start of epoch (seconds elapsed from EDF start)" );
  add_var( "EPOCH" , "E" , "STOP" , "Stop of epoch (seconds elapsed from EDF start)" );
  add_var( "EPOCH" , "E" , "TP" , "Interval in time-points" );

  // EPOCH-ANNOT

  add_cmd( "epoch" , "EPOCH-ANNOT" , "Attach epoch-level annotations from a file, to an epoched EDF" );
  add_url( "EPOCH-ANNOT" , "epochs/#epoch-annot" );  
  add_param( "EPOCH-ANNOT" , "file" , "annots/id1.epochs" , "File path/name to read annotations from [required]" );
  add_param( "EPOCH-ANNOT" , "recode" , "NREM1=N1,NREM2=N2" , "Comma-delimited list of recodings (from=to)");
  
  
  
  /////////////////////////////////////////////////////////////////////////////////
  //
  // MASKS
  //
  /////////////////////////////////////////////////////////////////////////////////


  // MASK

  add_cmd( "mask" , "MASK" , "Mask epochs based on annotations and other features" );
  add_url( "MASK" , "masks/#mask" );

  add_param( "MASK" , "if"    , "NREM2" , "Mask NREM2 epochs, unmask all others" );
  add_param( "MASK" , "ifnot" , "NREM2" , "Unmask NREM2 epochs, mask all others" );
  add_param( "MASK" , "expr" , "A>2" , "Mask epochs with A>2, unmask all others" );
  add_param( "MASK" , "not-expr" , "A>2" , "Unmask epochs with A>2, mask all others" );

  add_param( "MASK" , "mask-if" , "NREM2" , "Mask NREM2 epochs" );
  add_param( "MASK" , "mask-ifnot" , "NREM2" , "Mask non-NREM2 epochs" );
  add_param( "MASK" , "mask-expr" , "A>2" , "Mask epochs with A>2" );

  add_param( "MASK" , "unmask-if" , "NREM2" , "Unask NREM2 epochs" );
  add_param( "MASK" , "unmask-ifnot" , "NREM2" , "Unask non-NREM2 epochs" );
  add_param( "MASK" , "unmask-expr" , "A>2" , "Unmask epochs with A>2" );

  add_param( "MASK" , "none" , "" ,  "Clear mask (i.e. unmask all)" );
  add_param( "MASK" , "clear" , "" , "Clear mask (i.e. unmask all)" );
  add_param( "MASK" , "include-all" , "" , "Clear mask (i.e. unmask all)" );

  add_param( "MASK" , "all" , "" , "Mask all epochs" );
  add_param( "MASK" , "total" , "" , "Mask all epochs" );
  add_param( "MASK" , "exclude-all" , "" , "Mask all epochs" );


  add_param( "MASK" , "epoch" , "1-10" , "Select epochs 1 to 10" );
  add_param( "MASK" , "sec" , "60-120" , "Select epochs overlapping this interval" );
  add_param( "MASK" , "hms" , "8:00-9:00" , "Select epochs overlapping this interval" );

  add_param( "MASK" , "random" , "20" , "Select 20 random (currently unmasked) epochs" );

  add_param( "MASK" , "flip" , "" , "Reverse all masks" );
  add_param( "MASK" , "leading" , "W" , "Remove all leading epochs matching W" );
  add_param( "MASK" , "flanked" , "REM,2" , "Select only REM epochs flanked by 2+ REM epochs before/after" );
  
  
  add_table( "MASK" , "EMASK" , "Output stratified by mask" );
  add_var( "MASK" , "EMASK", "N_MATCHES" , "Number of epochs that match the condition (e.g. having annotation A)");
  add_var( "MASK" , "EMASK", "N_MASK_SET" , "Number of previously unmasked epochs that were masked by this operation");
  add_var( "MASK" , "EMASK", "N_MASK_UNSET" , "Number of previously masked epochs that were unmasked by this operation");
  add_var( "MASK" , "EMASK", "N_UNCHANGED" , "Number of epochs whose mask status was not changed by this operation");
  add_var( "MASK" , "EMASK", "N_RETAINED" , "Number of epochs retained after this operation");
  add_var( "MASK" , "EMASK", "N_TOTAL" , "Total number of epochs");
  
  // DUMP-MASK
  
  add_cmd( "mask" , "DUMP-MASK" , "Output epoch-level mask information" );
  add_url( "DUMP-MASK" , "masks/#dump-mask" );

  add_table( "DUMP-MASK" , "E" , "Epoch-level mask tabulation" );
  add_var( "DUMP-MASK" , "E" , "EMASK" , "Mask status: 0=unmasked (included), 1=masked (excluded)" );
 

  // RE (or RESTRUCTURE)

  add_cmd( "mask"   , "RE" , "Restructure an EDF (drop channels/epochs)" );
  add_url( "RE" , "masks/#restructure" );

  add_table( "RE" , "" , "Restructured data duration" );
  add_var( "RE" , "" , "DUR1" , "Duration pre-restructuring (secs)");
  add_var( "RE" , "" , "DUR2" , "Duration post-restructuring (secs)");
  add_var( "RE" , "" , "NR1" , "Duration pre-restructuring (records)");
  add_var( "RE" , "" , "NR2" , "Duration post-restructuring (records)");

  add_cmd( "mask"   , "RESTRUCTURE" , "Restructure an EDF (drop channels/epochs)" );
  add_url( "RESTRUCTURE" , "masks/#restructure" );

  add_table( "RESTRUCTURE" , "" , "Restructured data duration" );
  add_var( "RESTRUCTURE" , "" , "DUR1" , "Duration pre-restructuring (secs)");
  add_var( "RESTRUCTURE" , "" , "DUR2" , "Duration post-restructuring (secs)");
  add_var( "RESTRUCTURE" , "" , "NR1" , "Duration pre-restructuring (records)");
  add_var( "RESTRUCTURE" , "" , "NR2" , "Duration post-restructuring (records)");


  //
  // CHEP
  //

  add_cmd( "mask" , "CHEP" , "CHannel/EPoch masks" );
  add_url( "CHEP" , "masks/#chep" );
  add_param( "CHEP" , "clear" , "" , "Clear CHEP mask" );
  add_param( "CHEP" , "load" , "file.txt" , "Load CHEP from file.txt" );
  add_param( "CHEP" , "bad-channels" , "C3,C5" , "Manually specify bad channels" );
  add_param( "CHEP" , "epochs" , "2,0.1" , "Mask epochs with 2 or more bad channels, or >10% bad channels" );
  add_param( "CHEP" , "channels" , "10,0.5" , "Mask channels with 10 or more bad epochs, or >50% bad epochs" );
  add_param( "CHEP" , "dump" , "" , "Write current CHEP mask to output" );
  add_param( "CHEP" , "save" , "file.txt" , "Write CHEP mask to file.txt" );

  // CHEP output....  
  add_table( "CHEP" , "CH" , "CHEP mask channel-wise summaries" );
  add_var( "CHEP" , "CH" , "CHEP" , "Masked epochs" );

  add_table( "CHEP" , "E" , "CHEP mask epoch-wise summaries" );
  add_var( "CHEP" , "E" , "CHEP" , "Masked channels" );

  add_table( "CHEP" , "CH,E" , "CHEP mask" );
  add_var( "CHEP" , "CH,E" , "CHEP" , "CHannel/EPoch mask" );
  
  
  /////////////////////////////////////////////////////////////////////////////////
  //
  // MANIPULATIONS
  //
  /////////////////////////////////////////////////////////////////////////////////
  
  // SIGNALS

  add_cmd( "manip" , "SIGNALS" , "Retain/remove specific EDF channels" );
  add_url( "SIGNALS" , "manipulatons/#signals" );
  add_param( "SIGNALS" , "drop" , "EMG,ECG" , "Drop channels EMG and ECG" );
  add_param( "SIGNALS" , "keep" , "C3,C4" , "Drop all channels except C3 and C4" );

  
  // CONTAINS

  add_cmd( "manip" , "CONTAINS" , "Tests for particular signals/annotations/staging being present" );
  add_url( "CONTAINS" , "manipulatons/#contains" );
  add_param( "CONTAINS" , "sig" , "EMG,ECG" , "Test for these signals" );
  add_param( "CONTAINS" , "annots" , "apnea,hypopnea" , "Test for these annotations" );
  add_param( "CONTAINS" , "stages" , "" , "Test for valid staging" );
  add_param( "CONTAINS" , "skip" , "" , "Skip to next EDF on failure" );

  add_table( "CONTAINS" , "" , "Base" );
  add_var( "CONTAINS" , "" , "STAGE_COUNTS" , "Sleep stage counts" );
  add_var( "CONTAINS" , "" , "UNIQ_STAGES" , "Number of unique stage labels" );
  add_var( "CONTAINS" , "" , "NA_REQ" , "Number of required annots" );
  add_var( "CONTAINS" , "" , "NA_OBS" , "Number of observed annots" );

  add_var( "CONTAINS" , "" , "NS_OBS" , "Number of required channels" );
  add_var( "CONTAINS" , "" , "NS_REQ" , "Number of observed channels" );
  add_var( "CONTAINS" , "" , "NS_TOT" , "Tot number of channels" );

  add_table( "CONTAINS" , "ANNOT" , "Annotation informationm" );
  add_var( "CONTAINS" , "CH" , "PRESENT" , "Annotation resent" );

  add_table( "CONTAINS" , "CH" , "Channel informationm" );
  add_var( "CONTAINS" , "CH" , "PRESENT" , "Channel present" );


  
  // COPY

  add_cmd( "manip" , "COPY" , "Duplicate one or more EDF channels" );
  add_url( "COPY" , "manipulations/#copy" );
  add_param( "COPY" , "sig" , "C3,C4" , "List of channels to duplicate" );
  add_param( "COPY" , "tag" , "V2" , "Tag add to new channel names, e.g. C3_V2 [required] " );
    

  // CANONICAL
  add_cmd( "manip" , "CANONICAL" , "Create canonical signals" );
  add_url( "CANONICAL" , "manipulations/#canonical" );
  add_param( "CANONICAL" , "file" , "csfile.txt" , "File with canonical signal definitions" );
  add_param( "CANONICAL" , "group" , "GRP1" , "Group (from csfile.txt)" );
  add_param( "CANONICAL" , "cs" , "EEG,LOC,ROC" , "Optional: only calculate these CS" );

  add_table( "CANONICAL" , "" , "Canonical signal summaries" );
  add_var( "CANONICAL" , "" , "CS_SET" , "Number of canonical signals set" );
  add_var( "CANONICAL" , "" , "CS_NOT" , "Number of canonical signals not set" );
  add_var( "CANONICAL" , "" , "USED_CH" , "Number of used EDF channels" );
  add_var( "CANONICAL" , "" , "UNUSED_CH" , "Number of ununsed EDF channels" );

  add_table( "CANONICAL" , "CS" , "Canonical signal information" );
  add_var( "CANONICAL" , "CS" , "DEFINED" , "Is canonical signal present/defined?" );
  add_var( "CANONICAL" , "CS" , "SIG" , "Primary signal" );
  add_var( "CANONICAL" , "CS" , "REF" , "Reference signal" );
  add_var( "CANONICAL" , "CS" , "SR" , "Sample rate" );
  add_var( "CANONICAL" , "CS" , "UNITS" , "Units for canonical signal" );
  add_var( "CANONICAL" , "CS" , "NOTES" , "Optional, notes" );

  add_table( "CANONICAL" , "CH" , "EDF channel information" );
  add_var( "CANONICAL" , "CH" , "DROPPED" , "Original channel dropped" );
  add_var( "CANONICAL" , "CH" , "USED" , "Not used in constructing canonical signals" );

  // RESAMPLE 

  add_cmd( "manip" , "RESAMPLE" , "Resample signal(s)" );
  add_url( "RESAMPLE" , "manipulations/#resample" );
  add_param( "RESAMPLE" , "sig" , "C3,C4" , "List of channels to resample" );
  add_param( "RESAMPLE" , "sr" , "200" , "New sampling rate (Hz) [required]" );
  
  // REFERENCE

  add_cmd( "manip" , "REFERENCE" , "Resample signal(s)" );
  add_url( "REFERENCE" , "manipulations/#resample" );
  add_param( "REFERENCE" , "sig" , "C3,C4" , "List of signals to re-reference" );
  add_param( "REFERENCE" , "ref" , "A1,A2" , "Signal(s) providing the reference [required]" );

  // uV

  add_cmd( "manip" , "uV" , "Converts a signal to uV units" );
  add_url( "uV" , "manipulations/#uv" );
  add_param( "uV" , "sig" , "C3,C4" , "List of signals to convert" );

  // mV

  add_cmd( "manip" , "mV" , "Converts a signal to mV units" );
  add_url( "mV" , "manipulations/#mv" );
  add_param( "mV" , "sig" , "C3,C4" , "List of signals to convert" );

  // FLIP

  add_cmd( "manip" , "FLIP" , "Flips the polarity of a signal" );
  add_url( "FLIP" , "manipulations/#flip" );
  add_param( "FLIP" , "sig" , "C3,C4" , "List of signals to flip" );

  add_table( "FLIP" , "CH" , "Tracking flipped channels" );
  add_var( "FLIP" , "CH" , "FLIP" , "Channel flipped" );
  
  
  // RECORD-SIZE

  add_cmd( "manip" , "RECORD-SIZE" , "Alters the record size of an EDF, and writes a new EDF" );
  add_url( "RECORD-SIZE" , "manipulations/#record-size" );
  add_param( "RECORD-SIZE" , "dur" , "1" , "New EDF record/block size" );
  add_param( "RECORD-SIZE" , "edf-dir" , "edfs/" , "Folder for writing new EDFs" );
  add_param( "RECORD-SIZE" , "edf-tag" , "rec1" , "Tag added to new EDFs" );
  add_param( "RECORD-SIZE" , "sample-list" , "s2.lst" , "Generate a sample-list pointing to the new EDFs" );

  add_table( "RECORD-SIZE" , "" , "Restructured data duration" );
  add_var( "RECORD-SIZE", "" , "NR1" , "Pre-restructure number of records" );
  add_var( "RECORD-SIZE", "" , "NR2" , "Post-restructure number of records" );
  add_var( "RECORD-SIZE", "" , "DUR1" , "Pre-restructure duration (seconds)" );
  add_var( "RECORD-SIZE", "" , "DUR2" , "Post-restructure duration (seconds)" );

  
  
  // ANON

  add_cmd( "manip" , "ANON" , "Strips EDF ID and and Start Date headers" );
  add_url( "ANON" , "manipulations/#anon" );


  /////////////////////////////////////////////////////////////////////////////////
  //
  // OUTPUTS
  //
  /////////////////////////////////////////////////////////////////////////////////

  // WRITE

  add_cmd( "output" , "WRITE" , "Write a new EDF file" );
  add_url( "WRITE" , "outputs/#write" );
  add_param( "WRITE" , "edf-dir" , "edfs/" , "Set folder where new EDFs should be written" );
  add_param( "WRITE" , "edf-tag" , "v2" , "Add a tag to each new EDF filename" );
  add_param( "WRITE" , "sample-list" , "v2.lst" , "Name of the new sample-list" );
  
  add_table( "WRITE", "" , "Misc output from pre-WRITE restructure" );
  add_var( "WRITE", "" , "NR1" , "Pre-restructure number of records" );
  add_var( "WRITE", "" , "NR2" , "Post-restructure number of records" );
  add_var( "WRITE", "" , "DUR1" , "Pre-restructure duration (seconds)" );
  add_var( "WRITE", "" , "DUR2" , "Post-restructure duration (seconds)" );


  // MATRIX

  add_cmd( "output" , "MATRIX" , "Dumps signal information to a file" );
  add_url( "MATRIX" , "outputs/#matrix" );
  add_param( "MATRIX" , "file" , "signals.txt" , "Required parameter, to specify the filename for the output" );
  add_param( "MATRIX" , "sig" , "C3,C4" , "Restrict output to these signal(s)" );
  add_param( "MATRIX" , "hms" , "" , "Add a clock-time column in hh:mm:ss format" );
  add_param( "MATRIX" , "hms2" , "" , "Add a clock-time column in hh:mm:ss:microsecond format" );
  add_param( "MATRIX" , "annot" , "X,Y" , "Add columns with values 1/0 to indicate the presence/absence of that annotation" );
  add_param( "MATRIX" , "min" , "" , "Minimal output to show only signal information (no headers or lead columns)" );
  

  // DUMP-RECORDS

  add_cmd( "output" , "DUMP-RECORDS" , "Writes detailed annotation and signal data to standard output" );
  add_url( "DUMP-RECORDS" , "outputs/#dump-records" );
  add_param( "DUMP-RECORDS" , "no-signals" , "" , "Do not show signal data" );
  add_param( "DUMP-RECORDS" , "no-annots" , "" , "Do not show annotation information" );


  // RECS
  
  add_cmd( "output" , "RECS" , "Dumps information on EDF record structure to standard out" );
  add_url( "RECS" , "outputs/#recs" );

  // SEGMENTS

  add_cmd( "output" , "SEGMENTS" , "Report on contiguous segments in an EDF/EDF+" );
  add_url( "SEGMENTS" , "outputs/#segments" );

  add_table( "SEGMENTS" , "" , "Number of contiguous segments" );
  add_var( "SEGMENTS" , "" , "NSEGS" , "Number of contiguous segments" );

  add_table( "SEGMENTS" , "SEG" , "Information on each segment" );
  add_var( "SEGMENTS" , "SEG" , "DUR_HR" , "Segment duration (hours)" );
  add_var( "SEGMENTS" , "SEG" , "DUR_MIN" , "Segment duration (minutes)" );
  add_var( "SEGMENTS" , "SEG" , "DUR_SEC" , "Segment duration (seconds)" );

  add_var( "SEGMENTS" , "SEG" , "START" , "Segment start (seconds)" );
  add_var( "SEGMENTS" , "SEG" , "START_HMS" , "Segment start (hh:mm:ss)" );
    
  add_var( "SEGMENTS" , "SEG" , "STOP" , "Segment stop (seconds)" );
  add_var( "SEGMENTS" , "SEG" , "STOP_HMS" , "Segment stop (hh:mm:ss)" );


  // WRITE-ANNOTS
  
  add_cmd( "output" , "WRITE-ANNOTS" , "Write all annotations to file" );
  add_url( "WRITE-ANNOTS" , "outputs/#write-annots" );
  
  add_param( "WRITE-ANNOTS" , "file" , "f1.xml" , "Required filename for output" );
  add_param( "WRITE-ANNOTS" , "luna" , "" , "Output in Luna .annot format instead of XML" );

  
  /////////////////////////////////////////////////////////////////////////////////
  //
  // FILTERS
  //
  /////////////////////////////////////////////////////////////////////////////////

  // FILTER

  add_cmd( "filter" , "FILTER" , "Apply a FIR filter to one or more signals");
  add_url( "FILTER" , "fir-filters/#filter" );
  add_param( "FILTER" , "sig" , "C3,C4" , "Restrict analysis to these channels" );

  add_param( "FILTER" , "bandpass" , "0.3,35" , "Band-pass filter between 0.3 and 35 Hz" );
  add_param( "FILTER" , "lowpass"  , "35"  , "Low-pass filter with cutoff of 35 Hz" );
  add_param( "FILTER" , "highpass" , "0.3" , "High-pass filter with cutiff of 0.3 Hz" );
  add_param( "FILTER" , "bandstop" , "55,65" , "Band-stop filter between 55 and 65 Hz" );
  add_param( "FILTER" , "ripple" , "0.02" , "Ripple (as a proportion)" );
  add_param( "FILTER" , "tw" , "1" , "Transition width (in Hz)" );


  // FILTER-DESIGN

  add_cmd( "filter" , "FILTER-DESIGN" , "Apply a FIR filter to one or more signals");
  add_url( "FILTER-DESIGN" , "fir-filters/#filter-design" );
  add_param( "FILTER-DESIGN" , "bandpass" , "0.3,35" , "Band-pass filter between 0.3 and 35 Hz" );
  add_param( "FILTER-DESIGN" , "lowpass"  , "35" , "Low-pass filter with cutoff of 35 Hz" );
  add_param( "FILTER-DESIGN" , "highpass" , "0.3" , "High-pass filter with cutiff of 0.3 Hz" );
  add_param( "FILTER-DESIGN" , "bandstop" , "55,65" , "Band-stop filter between 55 and 65 Hz" );
  add_param( "FILTER-DESIGN" , "ripple" , "0.02" , "Ripple (as a proportion)" );
  add_param( "FILTER-DESIGN" , "tw" , "1" , "Transition width (in Hz)" );
  add_param( "FILTER-DESIGN" , "fs" , "200" , "Specify sample rate (in Hz)" ); 


  /////////////////////////////////////////////////////////////////////////////////
  //
  // ARTIFACTS
  //
  /////////////////////////////////////////////////////////////////////////////////

  // SIGSTATS

  add_cmd( "artifact" , "SIGSTATS" , "Per-epoch outlier detection (RMS, Hjorth parameters, clipped signals)" );
  add_url( "SIGSTATS" , "artifacts/#sigstats" );
  add_param( "SIGSTATS" , "sig" , "C3,C4" , "Restrict analysis to these channels" );

  add_param( "SIGSTATS" , "verbose" , "" , "Report epoch-level statistics" );
  add_param( "SIGSTATS" , "epoch" , "" , "Report epoch-level statistics (same as verbose)" );
  add_param( "SIGSTATS" , "chep" , "" , "Set CHEP mask for outlier epochs" );
  add_param( "SIGSTATS" , "astats" , "3,3" , "Between-epoch, betwee-channel filtering" );
  add_param( "SIGSTATS" , "cstats" , "2" , "Within-epoch, between-channel filtering" );

  add_param( "SIGSTATS" , "rms" , "" , "Calculate/mask on RMS" );
  add_param( "SIGSTATS" , "clipped" , "0.05" , "Calculate/mask on signal clipping" );
  add_param( "SIGSTATS" , "flat" , "0.05" , "Calculate/mask on signal clipping" );
  add_param( "SIGSTATS" , "max" , "0.05" , "Calculate/mask on signal clipping" );

  add_param( "SIGSTATS" , "threshold" , "2,2" , "Set eppoch masks based on SD unit (iterative) outlier detection" );
  add_param( "SIGSTATS" , "th" , "2,2" , "Same as 'threshold'" );
  add_param( "SIGSTATS" , "cstats" , "2" , "Run channel-comparisons, with threshold in SD units" );
  add_param( "SIGSTATS" , "cstats-unmasked-only" , "" , "Channel-comparisons only for unmasked epochs" );

  add_table( "SIGSTATS" , "CH" , "Per-channel whole-signal statistics" );
  add_var( "SIGSTATS" , "CH" , "CLIP" , "Proportion of clipped sample points" );
  add_var( "SIGSTATS" , "CH" , "FLAT" , "Proportion of flat sample points" );
  add_var( "SIGSTATS" , "CH" , "MAX" , "Proportion of max sample points" );
  add_var( "SIGSTATS" , "CH" , "H1" , "First Hjorth parameter (activity)" );
  add_var( "SIGSTATS" , "CH" , "H2" , "Second Hjorth parameter (mobility)" );
  add_var( "SIGSTATS" , "CH" , "H3" , "Third Hjorth parameter (complexity)" );
  add_var( "SIGSTATS" , "CH" , "RMS" , "Signal root mean square" );

  add_var( "SIGSTATS" , "CH" , "P_H1" , "Proportion flagged epochs for H1 [cstats]" );
  add_var( "SIGSTATS" , "CH" , "P_H2" , "Proportion flagged epochs for H2 [cstats]" );
  add_var( "SIGSTATS" , "CH" , "P_H3" , "Proportion flagged epochs for H3 [cstats]" );
  add_var( "SIGSTATS" , "CH" , "P_OUT" , "Proportion flagged epochs for H1, H2 or H3 [cstats]" );

  add_var( "SIGSTATS" , "CH" , "Z_H1" , "Z score for H1 [cstats]" );
  add_var( "SIGSTATS" , "CH" , "Z_H2" , "Z score for H2 [cstats]" );
  add_var( "SIGSTATS" , "CH" , "Z_H3" , "Z score for H3 [cstats]" );

  add_var( "SIGSTATS" , "CH" , "CNT_ACT" , "Number of epochs flagged based on H1 [mask]" );
  add_var( "SIGSTATS" , "CH" , "CNT_MOB" , "Number of epochs flagged based on H2 [mask]" );
  add_var( "SIGSTATS" , "CH" , "CNT_CMP" , "Number of epochs flagged based on H3 [mask]" );
  add_var( "SIGSTATS" , "CH" , "CNT_CLP" , "Number of epochs flagged based on clipping metric" );
  add_var( "SIGSTATS" , "CH" , "CNT_RMS" , "Number of epochs flagged based on RMS" );

  add_var( "SIGSTATS" , "CH" , "FLAGGED_EPOCHS" , "Number of epochs flagged as outliers [mask]" );
  add_var( "SIGSTATS" , "CH" , "ALTERED_EPOCHS" , "Number of epochs whose mask was altered [mask]" );
  add_var( "SIGSTATS" , "CH" , "TOTAL_EPOCHS" , "Total number of masked epochs [mask]" );

  add_table( "SIGSTATS" , "CH,E" , "Per-channel per-epoch statistics [epoch]" );
  add_var( "SIGSTATS" , "CH,E" , "H1" , "First Hjorth parameter (activity)" );
  add_var( "SIGSTATS" , "CH,E" , "H2" , "Second Hjorth parameter (mobility)" );
  add_var( "SIGSTATS" , "CH,E" , "H3" , "Third Hjorth parameter (complexity)" );
  hide_var( "SIGSTATS" , "CH,E" , "CLIP" , "Proportion of clipped sample points" );
  hide_var( "SIGSTATS" , "CH,E" , "FLAT" , "Proportion of flat sample points" );
  hide_var( "SIGSTATS" , "CH,E" , "MAX" , "Proportion of max sample points" );
  hide_var( "SIGSTATS" , "CH,E" , "RMS" , "Signal root mean square" );


  // ARTIFACTS

  add_cmd( "artifact" , "ARTIFACTS" , "Detect EEG artifacts following Buckelmueller et al." );
  add_url( "ARTIFACTS" , "artifacts/#artifacst" );
  add_param( "ARTIFACTS" , "sig" , "C3,C4" , "Restrict analysis to these channels" );
  add_param( "ARTIFACTS" , "verbose" , "" , "Report epoch-level statistics" );
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


  // SUPPRESS-ECG

  add_cmd( "artifact" , "SUPPRESS-ECG" , "Detect/remove cardiac-contamination from the EEG" );
  add_url( "SUPPRESS-ECG" , "artifacts/#suppress-ecg" );
  add_param( "SUPPRESS-ECG" , "sig" , "C3,C4" , "Restrict analysis to these channels" );
  add_param( "SUPPRESS-ECG" , "sr" , "125" , "Set sample rate for ECG/EEG channels" );
  add_param( "SUPPRESS-ECG" , "no-suppress" , "" , "Do not alter any EEG channels" );
  
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
  add_var( "SUPPRESS-ECG" , "CH,SP" , "ART_RMS" , "Estimate correction factor, for each sample point in a 2-second window" );


  // EMD

  add_cmd( "power" , "EMD" , "Empirical mode decomposition" );
  add_url( "EMD" , "power-spectra/#emd" );
  add_param( "EMD" , "sig" , "C3,C4" , "Select signals for EMD" );
  add_param( "EMD" , "tag" , "_C_" , "IMF channel tag, if not _IMF_" );
  add_param( "EMD" , "sift" , "20" , "Maximum number of sifting operations" );
  add_param( "EMD" , "imf" , "10" , "Maximum number of IMF to extract" );
    

  // ALTER

  add_cmd( "artifact" , "ALTER" , "Regression- or EMD-based artifact correction" );
  add_url( "ALTER" , "artifacts/#alter" ); 
  add_param( "ALTER" , "sig" , "C3,C4" , "Signals for analysis" );
  add_param( "ALTER" , "corr" , "EOG-R,EOG-L" , "Template signal(s)" );
  add_param( "ALTER" , "emd" , "" , "Use EMD instead of raw regression" );
  add_param( "ALTER" , "th" , "0.9" , "Threshold" );
  add_param( "ALTER" , "emd-corr" , "" , "Run EMD of corrector channels" );

  add_param( "ALTER" , "segment-sec" , "4" , "Segment size" );
  add_param( "ALTER" , "segment-step" , "2" , "Segment step (half size by default)" );
  
  
  /////////////////////////////////////////////////////////////////////////////////
  //
  // HYPNOGRAMS
  //
  /////////////////////////////////////////////////////////////////////////////////


  add_cmd( "hypno" , "STAGE" , "Output sleep stage annotations, per epoch" );
  add_url( "STAGE" , "hypnograms/#stage" );
  add_param( "STAGE" , "N1" , "NREM1" , "Set the annotation used for N1 sleep" );
  add_param( "STAGE" , "N2" , "NREM2" , "Set the annotation used for N2 sleep" );
  add_param( "STAGE" , "N3" , "NREM3" , "Set the annotation used for N3 sleep" );
  add_param( "STAGE" , "REM" , "REM" , "Set the annotation used for REM sleep" );
  add_param( "STAGE" , "wake" , "W" ,  "Set the annotation used for N3 sleep" );
  add_param( "STAGE" , "?" , "-9" , "Set the annotation used for unknown/other" );

  add_table( "STAGE" , "E" , "Stage annotations per-epoch" );
  add_var( "STAGE" , "E" , "CLOCK_TIME" , "Clock time (hh:mm:ss)" );
  add_var( "STAGE" , "E" , "MINS" , "Elapsed time from start of EDF (minutes)" );
  add_var( "STAGE" , "E" , "STAGE" , "Sleep stage (text value)" );
  add_var( "STAGE" , "E" , "STAGE_N" , "Numeric encoding of sleep stage" );


  add_cmd( "hypno" , "HYPNO" , "Metrics based on sleep stage annotations" );
  add_url( "HYPNO" , "hypnograms/#hypno" );

  add_param( "HYPNO" , "file" , "stages.txt" , "Optionally, read stages from file" );
  add_param( "HYPNO" , "N1" , "NREM1"  , "Set the annotation used for N1 sleep" );
  add_param( "HYPNO" , "N2" , "NREM2" , "Set the annotation used for N2 sleep" );
  add_param( "HYPNO" , "N3" , "NREM3" , "Set the annotation used for N3 sleep" );
  add_param( "HYPNO" , "REM" , "REM" , "Set the annotation used for REM sleep" );
  add_param( "HYPNO" , "wake" , "W" , "Set the annotation used for N3 sleep" );
  add_param( "HYPNO" , "?" , "-9" , "Set the annotation used for unknown/other" );


  add_table( "HYPNO" , "" , "Individual-level output" );
  add_var( "HYPNO" , "" , "TRT" , "Total sleep time" );
  add_var( "HYPNO" , "" , "TST" , "Total sleep time" );
  add_var( "HYPNO" , "" , "TST_PER" , "Total persistent sleep time" );
  add_var( "HYPNO" , "" , "TIB" , "Time in bed" );
  add_var( "HYPNO" , "" , "SPT" , "Sleep period time" );
  add_var( "HYPNO" , "" , "SPT_PER" , "Persistent sleep period time" );
  add_var( "HYPNO" , "" , "TWT" , "Total wake time" );
  add_var( "HYPNO" , "" , "WASO" , "Wake after sleep onset" );
  add_var( "HYPNO" , "" , "FWT" , "Final wake time" );
  add_var( "HYPNO" , "" , "LOT" , "Lights On time" );
  add_var( "HYPNO" , "" , "LOST" , "Lights On sleep time" );
  add_var( "HYPNO" , "" , "SINS" , "Study Starts In Sleep" );
  add_var( "HYPNO" , "" , "EINS" , "Study Ends In Sleep" );

  add_var( "HYPNO" , "" , "OTHR" , "Unknown stage duration" );
  add_var( "HYPNO" , "" , "CONF" , "Number of epochs with conflicting stage assignments" );
    

  add_var( "HYPNO" , "" , "FIXED_WAKE" , "Epochs fixed due to excessive WASO" );
  add_var( "HYPNO" , "" , "FIXED_LIGHTS" , "Epochs fixed due to Lights On" );

  add_var( "HYPNO" , "" , "MINS_ASC_N2" , "Duration of ascending N2 (mins)" );
  add_var( "HYPNO" , "" , "MINS_DSC_N2" , "Duration of descending N2 (mins)" );
  add_var( "HYPNO" , "" , "MINS_FLT_N2" , "Duration of flat N2 (mins)" );
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
  
  add_table( "HYPNO" , "SS" , "Stage-stratified output" );
  add_var( "HYPNO" , "SS" , "MINS" , "Stage duration (mins)" );
  add_var( "HYPNO" , "SS" , "PCT" , "Stage duration (% of TST)" );
  add_var( "HYPNO" , "SS" , "BOUT_N" , "Number of bouts" );
  add_var( "HYPNO" , "SS" , "BOUT_MN" , "Mean bout duration" );
  add_var( "HYPNO" , "SS" , "BOUT_MD" , "Median bout duration" );
  add_var( "HYPNO" , "SS" , "BOUT_5" , "Stage duration (only bouts 5+ mins)" );
  add_var( "HYPNO" , "SS" , "BOUT_10" , "Stage duration (only bouts 10+ mins)" );

  add_table( "HYPNO" , "C" , "NREM cycle-level output" );
  add_var( "HYPNO" , "C" , "NREMC_START" , "First epoch number of this NREM cycle" );
  add_var( "HYPNO" , "C" , "NREMC_MINS" , "Total duration of this cycle (mins)" );
  add_var( "HYPNO" , "C" , "NREMC_NREM_MINS" , "Duration of NREM in this cycle (mins)" );
  add_var( "HYPNO" , "C" , "NREMC_REM_MINS" , "Duration of REM in this cycle (mins)" );
  add_var( "HYPNO" , "C" , "NREMC_OTHER_MINS" , "Minutes of wake and unscored epochs" );


  add_table( "HYPNO", "N" , "Bouts" ); 
  add_var( "HYPNO" , "N" , "FIRST_EPOCH" , "First epoch" );
  add_var( "HYPNO" , "N" , "LAST_EPOCH" , "Last epoch" );
  add_var( "HYPNO" , "N" , "START" , "Start (clocktime)" );
  add_var( "HYPNO" , "N" , "STOP" , "Stop (clocktime) [ end of last epoch ]" );
  add_var( "HYPNO" , "N" , "MINS" , "Bout duration (minutes)" );
  

  add_table( "HYPNO" , "E" , "Epoch-level output" );
  add_var( "HYPNO" , "E" , "CLOCK_HOURS" , "Start time of epoch (hours since midnight)" );
  add_var( "HYPNO" , "E" , "CLOCK_TIME" , "Start time of epoch (hh:mm:ss)" );
  add_var( "HYPNO" , "E" , "MINS" , "Start time of epoch (minutes since start of EDF)" );
  add_var( "HYPNO" , "E" , "STAGE" , "Text description of sleep stage" );
  add_var( "HYPNO" , "E" , "STAGE_N" , "Numeric encoding of sleep stage" );
  add_var( "HYPNO" , "E" , "PERSISTENT_SLEEP" , "Flag to indicate persistent sleep" );
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
  add_var( "HYPNO" , "E" , "FLANKING_SIM" , "Number of similarly-staged epochs,either forwards or backwards" );

  add_var( "HYPNO" , "E" , "N2_WGT" , "Score to indicate ascending versus descending N2 sleep" );
  add_var( "HYPNO" , "E" , "NEAREST_WAKE" , "Number of epochs (forward or backwards) since nearest wake epoch" );
  add_var( "HYPNO" , "E" , "NREM2REM" , "Number of epochs from this N2 epoch to the N2/REM transition" );
  add_var( "HYPNO" , "E" , "NREM2REM_TOTAL" , "Total number of contiguous N2 epochs until a REM transition" );
  add_var( "HYPNO" , "E" , "NREM2WAKE" , "Number of epochs from this N2 epoch to the N2/Wake transition" );
  add_var( "HYPNO" , "E" , "NREM2WAKE_TOTAL" , "Total number of contiguous N2 epochs until a Wake transition" );
  add_var( "HYPNO" , "E" , "CYCLE" , "Cycle number, if this epoch is in a sleep cycle" );
  add_var( "HYPNO" , "E" , "CYCLE_POS_ABS" , "Absolute position of this epoch in the current NREM cycle (mins)" );
  add_var( "HYPNO" , "E" , "CYCLE_POS_REL" , "Relative position of this epoch in the current NREM cycle (0-1)" );
  add_var( "HYPNO" , "E" , "PERIOD" , "Cycle period: NREMP or REMP, or missing if not in a cycle" );

  add_table( "HYPNO" , "C" , "NREM cycle-level output" );

  add_table( "HYPNO" , "PRE,POST" , "Stage transitions" );
  add_var( "HYPNO", "PRE,POST" , "N" , "Number of transitions" );
  add_var( "HYPNO" , "PRE,POST" , "P_POST_COND_PRE" , "P( S+1 | S )" );
  add_var( "HYPNO" , "PRE,POST" , "P_PRE_COND_POST" , "P( S | S+1 )" );



  /////////////////////////////////////////////////////////////////////////////////
  //
  // SUDS
  //
  /////////////////////////////////////////////////////////////////////////////////

  add_cmd( "staging"   , "SOAP" , "Single Observation Accuracies and Probabilities" );
  add_url( "SOAP" , "suds/#soap" );

  add_param( "SOAP" , "sig" , "C3,C4" , "Restrict analysis to these channels" );
  add_param( "SOAP" , "nc" , "10" , "Number of principal components" );
  
  add_table( "SOAP" , "" , "Overall accuracies" );
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
  // SUDS
  //

  add_cmd( "staging"   , "SUDS" , "Staging Using the Dynamics of Sleep" );
  add_url( "SUDS" , "suds/#suds" );

  add_param( "SUDS" , "sig" , "C3,C4" , "Restrict analysis to these channels" );
  add_param( "SUDS" , "nc" , "10" , "Number of principal components" );
  add_param( "SUDS" , "lambda" , "2" , "Regularization factor" );
  add_param( "SUDS" , "th" , "5,3" , "Statistical outlier removal" );
  add_param( "SUDS" , "robust" , "0.1" , "Robust standardization" );
  add_param( "SUDS" , "wgt-exp" , "4" , "Expoential weighting" );

  add_table( "SUDS" , "" , "SUDS metrics" );
  add_var( "SUDS", "" , "ACC" , "Accuracy" );
  add_var( "SUDS", "" , "ACC3" , "Accuracy for 3-class model" );  
  add_var( "SUDS", "" , "K" , "Kappa statistic" );
  add_var( "SUDS", "" , "K3" , "Kappa for 3-class model" );
  add_var( "SUDS", "" , "F1" , "F1 statistic" );
  add_var( "SUDS", "" , "F13" , "F1 for 3-class model" );
  add_var( "SUDS", "" , "F1_WGT" , "F1 weighted" );
  add_var( "SUDS", "" , "MAXPP" , "Mean maximum posterior" );
  add_var( "SUDS", "" , "MCC" , "Matthews correlation coefficient" );
  add_var( "SUDS", "" , "MCC3" , "Matthews correlation coefficient, 3-class" );
  add_var( "SUDS", "" , "PREC" , "Precision" );
  add_var( "SUDS", "" , "PREC_WGT" , "Precision, weighted" );
  add_var( "SUDS", "" , "PREC3" , "Precision, 3-class" );
  add_var( "SUDS", "" , "RECALL" , "Recall" );
  add_var( "SUDS", "" , "RECALL3" , "Recall, 3-class" );
  add_var( "SUDS", "" , "RECALL_WGT" , "Recall,weighted" );
  add_var( "SUDS", "" , "R_WGT" , "Correlation between weight and K3" );

  add_table( "SUDS" , "E" , "Epoch-level SUDS predictions" );
  add_var( "SUDS" , "E" , "DISC" , "Discordant prior/predicted w.r.t 5-classes" );
  add_var( "SUDS" , "E" , "DISC3" , "Discordant prior/predicted w.r.t 3-classes" );
  add_var( "SUDS" , "E" , "INC" , "0/1 for whether epoch was included in analysis" );
  add_var( "SUDS" , "E" , "PP_N1" , "Posterior probability of N1" );
  add_var( "SUDS" , "E" , "PP_N2" , "Posterior probability of N2" );
  add_var( "SUDS" , "E" , "PP_N3" , "Posterior probability of N3" );
  add_var( "SUDS" , "E" , "PP_R" , "Posterior probability of REM" );
  add_var( "SUDS" , "E" , "PP_W" , "Posterior probability of wake" );
  add_var( "SUDS" , "E" , "PRED", "Predicted stage" );
  add_var( "SUDS" , "E" , "PRIOR" , "Observed stage (if known)" );
  
  add_table( "SUDS" , "SS" , "Sleep-stage summaries" );
  add_var( "SUDS" , "SS" , "DUR_OBS" , "Observed stage duration (for included epochs)");
  add_var( "SUDS" , "SS" , "DUR_PRD" , "Predicted stage duration, weighted" );
  add_var( "SUDS" , "SS" , "DUR_PRD2" , "Predicted stage duration, based on most likely" );
  add_var( "SUDS" , "SS" , "F1" , "F1 statistic" );
  add_var( "SUDS" , "SS" , "RECALL" , "Recall" );
  add_var( "SUDS" , "SS" , "PREC" , "Precision" );

  add_table( "SUDS" , "TRAINER" , "Trainer-level metrics" );
  add_var( "SUDS" , "TRAINER" , "K3" , "3-class kappa" );
  add_var( "SUDS" , "TRAINER" , "NS" , "Number of unique stages in prediction" );
  add_var( "SUDS" , "TRAINER" , "N_N1" , "N1 duration" );
  add_var( "SUDS" , "TRAINER" , "N_N2" , "N2 duration" );
  add_var( "SUDS" , "TRAINER" , "N_N3" , "N3 duration" );
  add_var( "SUDS" , "TRAINER" , "N_REM" , "REM duration" );
  add_var( "SUDS" , "TRAINER" , "N_W" , "Wake duration" );
  add_var( "SUDS" , "TRAINER" , "WGT" , "Trainer weight" );
  
  add_table( "SUDS" , "WTRAINER" , "Weight-trainer metrics" );
  add_var( "SUDS" , "WTRAINER" , "K3" , "Mean weight trainer K3" );

  add_table( "SUDS" , "E,TRAINER" , "Verbose trainer metrics" );
  add_var( "SUDS" , "E,TRAINER" , "PP_N1" , "Posterior probability of N1" );
  add_var( "SUDS" , "E,TRAINER" , "PP_N2" , "Posterior probability of N2" );
  add_var( "SUDS" , "E,TRAINER" , "PP_N3" , "Posterior probability of N3" );
  add_var( "SUDS" , "E,TRAINER" , "PP_R" , "Posterior probability of REM" );
  add_var( "SUDS" , "E,TRAINER" , "PP_W" , "Posterior probability of wake" );
  add_var( "SUDS" , "E,TRAINER" , "PRED" , "Predicted (most likely) stage" );
  set_compressed( "SUDS" , tfac_t( "E,TRAINER" ) );

  add_table( "SUDS" , "NSS,PRED,OBS" , "Confusion matrix" );
  add_var( "SUDS" , "NSS,PRED,OBS" , "N" , "Number" );
  add_var( "SUDS" , "NSS,PRED,OBS" , "P" , "Proportion" );


  //
  // POPS
  //

  add_cmd( "staging"   , "POPS" , "Population-based staging" );
  add_url( "POPS" , "staging/#pops" );

  add_param( "POPS" , "train" , "" , "Build POPS training datasets" );
  add_param( "POPS" , "features" , "m1.ftr" , "Feature specification file" );
  add_param( "POPS" , "data" , "pops/lib/^" , "Filename for bimary training files" );
  add_param( "POPS" , "model" , "m1.model" , "LGBM model file to write to/read from" );
  add_param( "POPS" , "config" , "m1.config" , "LGBM configuration file" );

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

  
  /////////////////////////////////////////////////////////////////////////////////
  //
  // SPECTRAL
  //
  /////////////////////////////////////////////////////////////////////////////////


  //
  // PSD
  //

  add_cmd( "power"   , "PSD" , "Power spectral density estimation (Welch)" );
  add_url( "PSD" , "power-spectra/#psd" );

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

  add_table( "PSD" , "CH,F" , "Whole-night, per-channel power" );
  add_var( "PSD" , "CH,F" , "PSD" , "Power (mean over epochs)" );
  add_var( "PSD" , "CH,F" , "PSD_MD" , "Power (median over epochs)" );
  add_var( "PSD" , "CH,F" , "PSD_SD" , "Power (SD over epochs)" );


  add_table( "PSD" , "CH,B,E" , "Whole-night, per-channel per-epoch band power" );
  add_var( "PSD" , "CH,B,E" , "PSD" , "Power" );
  add_var( "PSD" , "CH,B,E" , "RELPSD" , "Relative power" );

  add_table( "PSD" , "CH,F,E" , "Whole-night, per-channel per-epoch power" );
  add_var( "PSD" , "CH,F,E" , "PSD" , "Power" );
  set_compressed( "PSD" , tfac_t( "CH,F,E" ) );

  add_table( "PSD" , "CH,E", "Epoch/channel level stats" );
  add_var( "PSD" , "CH,E" , "KURT" , "Peak (PSD kurtosis)" );
  add_var( "PSD" , "CH,E" , "SPK" , "Sum PSD peakedness" );
  add_var( "PSD" , "CH,E" , "SPEC_SLOPE" , "Spectral slope" );
  add_var( "PSD" , "CH,E" , "SPEC_SLOPE_N" , "Spectral slope number of points" );


  //
  // ASYMM
  //

  add_cmd( "power" , "ASYMM" , "EEG asymmetry" );
  add_url( "ASYMM" , "power-spectra/#asymm" );

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
  // PSC  (nb. --psc uses PSC command label internally)
  //

  add_cmd( "cmdline" , "--psc" , "Create PSC from power spectra from multiple individuals" );
  add_url( "--psc" , "psc/#build-psc" );
  add_param( "--psc" , "spectra" , "psd1.txt,psd2.txt" , "File list of PSD/COH output" );
  add_param( "--psc" , "v" , "PSD,COH" , "List of variables to extract" );
  add_param( "--psc" , "log" , "PSD" , "Take log of these variables" );
  add_param( "--psc" , "proj" , "proj1.txt" , "Write projection file to disk" );
  add_param( "--psc" , "output-input", "mat1.txt" , "Write constructed input matrix to disk" );
  add_param( "--psc" , "nc", "10" , "Number of PSCs to extract (default 10)" );
  add_param( "--psc" , "th", "5,3" , "Iterative SD thresholds for outlier removal" );
  
  add_cmd( "power"   , "PSC" , "Calculate/apply Power spectral density estimation (Welch)" );
  add_url( "PSC" , "psc/#project-psc" );
  add_param( "PSC" , "proj" , "proj1.txt" , "PSC projection file (from --psc) ");
  add_param( "PSC" , "nc" , "5" , "Number of components (if subset of projection desired (default all)" );

  add_table( "PSC" , "PSC" , "Principal spectral components" );
  add_var( "PSC", "PSC" , "U" , "Principal spectral component value" );

  add_table( "PSC" , "I" , "Singular values/variance explained" );
  add_var( "PSC", "I" , "W" , "Singular value" );
  add_var( "PSC", "I" , "VE" , "Variance explained" );

  add_table( "PSC" , "J" , "Variable labels" );
  add_var( "PSC" , "J" , "CH" , "Channel label" );
  add_var( "PSC" , "J" , "F" , "Frequency" );
  add_var( "PSC" , "J" , "VAR" , "Variable" );
  
  add_table( "PSC" , "I,J" , "V matrix" );
  add_var( "PSC" , "I,J" , "V" , "V" );
  

  //
  // MTM
  //

  add_cmd( "power"   , "MTM" , "Power spectral density estimation (Welch)" );
  add_url( "MTM" , "power-spectra/#mtm" );

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

  add_table( "MTM", "CH,SEG", "Segment timing details" );
  add_var( "MTM" , "CH,SEG" , "START" , "Start time (seconds)" );
  add_var( "MTM" , "CH,SEG" , "STOP" , "Stop time (seconds)");
  add_var( "MTM" , "CH,SEG" , "DISC" , "Spans a discontinuity (0/1=N/Y)");
  
  add_table( "MTM" , "CH,F" , "Whole-night, per-channel power" );
  add_var( "MTM" , "CH,F" , "MTM" , "Power" );

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

  add_param( "CWT-DESIGN" , "sr" , "200" , "Sampling rate" );
  add_param( "CWT-DESIGN" , "fc" , "15" , "Wavelet center frequency" );
  add_param( "CWT-DESIGN" , "cycles" , "7" , "Bandwidth of the wavelet (number of cycles)" );

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

  add_param( "1FNORM" , "sig" , "C3,C4" , "Restrict analysis to these channels" );

  //
  // TV
  //

  add_cmd( "power" , "TV" , "Applies of fast algorithm for 1D total variation denoising" );

  add_url( "TV" , "power-spectra/#tv" );

  add_param( "TV" , "sig" , "C3,C4" , "Restrict analysis to these channels" );
  add_param( "TV" , "lambda" , "10" , "Smoothing parameter (0 to infinity)" );


  /////////////////////////////////////////////////////////////////////////////////
  //
  // SPINDLES/SO
  //
  /////////////////////////////////////////////////////////////////////////////////

  add_cmd( "transients" , "SPINDLES" , "Wavelet-based sleep spindle detection" );
  add_url( "SPINDLES" , "spindles-so/#spindles" );

  add_param( "SPINDLES" , "sig" , "C3,C4" , "Restrict analysis to these channels" );

  add_param( "SPINDLES" , "fc" , "11,15" , "Restrict analysis to these channels (otherwise, all channels are included)" );
  add_param( "SPINDLES" , "cycles" , "12" , "Number of cycles (default 7)" );
  add_param( "SPINDLES" , "th" , "6" , "Multiplicative threshold for core spindle detection (default 4.5)" );
  add_param( "SPINDLES" , "th2" , "3" , "Multiplicative threshold for non-core spindle detection (default=2)" );
  add_param( "SPINDLES" , "median" , "" , "Flag to indicate that the median, not mean, is used for thresholding" );
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
  
  add_param( "SPINDLES" , "empirical" , "" , "Empirically determine thresholds" );
  hide_param( "SPINDLES" , "set-empirical" , "" , "Use empirically determined thresholds for spindle detection" );
  hide_param( "SPINDLES" , "verbose-empirical" , "" , "Output extensive information on threshold estimation" );

  add_param( "SPINDLES" , "merge" , "0.2", "Merge two putative spindles if within this interval (default 0.5 seconds)" );
  add_param( "SPINDLES" , "collate" , "" , "Within each channel, collate overlapping spindles of similar frequencies" );
  add_param( "SPINDLES" , "collate-channels" , "" , "As above, except merge across channels also" );
  add_param( "SPINDLES" , "th-frq" , "1" , "Frequency criterion for merging spindles (default 2 Hz)" );
  add_param( "SPINDLES" , "list-all-spindles" , "" , "List all spindles that comprise each m-spindle" );

  add_param( "SPINDLES" , "th-interval" , "0.5" , "Merge if the ratio of intersection to union is at least this (default 0, i.e. any overlap)" );
  hide_param( "SPINDLES" , "th-interval-cross-channel" , "" , "not currently used" );
  hide_param( "SPINDLES" , "th-interval-within-channel" , "" , "not currently used" );
  add_param( "SPINDLES" , "window" , "0.5" , "Set window around each spindle when defining temporal overlap" );
  add_param( "SPINDLES" , "hms" , "" , "Show clock-time of each m-spindle" );

  add_param( "SPINDLES" , "ftr" , "tag" , "Produce FTR files for all spindles, with the tag in the filename" );
  add_param( "SPINDLES" , "ftr-dir" , "/path/to/folder" , "Folder for FTR files" );
  hide_param( "SPINDLES" , "show-coef" , "" , "Request (very verbose) coefficient output (to stdout)" );

  // output

  add_table( "SPINDLES" , "CH,F" , "Individual-level output" );
  add_var( "SPINDLES" , "CH,F" , "DENS" , "Spindle density (count per minute)" );
  add_var( "SPINDLES" , "CH,F" , "AMP" , "Mean spindle amplitude (uV or mV units)" );
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

  add_table( "SPINDLES" , "CH,F,TH" , "Between-class variance over range of thresholds" );
  add_var( "SPINDLES" , "CH,F,TH" , "SIGMAB" , "Between-class variance for given threshold" );

  add_table( "SPINDLES" , "CH,E,F" , "Epoch-level output [epoch]" ); 
  add_var( "SPINDLES" , "CH,E,F" , "N" , "Number of spindles observed in that epoch (for that target frequency/channel)" );

  add_table( "SPINDLES" , "CH,F,SPINDLE" , "Spindle-level output [per-spindle]" ); 
  add_var( "SPINDLES" , "CH,F,SPINDLE" , "AMP" , "Spindle amplitude (uV or mV units)" );
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
  hide_var( "SPINDLES" , "CH,F,SPINDLE" , "IF" , "Mean frequency per spindle over duration [if]" );
  
  hide_table( "SPINDLES" , "CH,F,RELLOC" , "Mean IF stratified by relative location in spindle [if]" );
  hide_var( "SPINDLES" , "CH,F,RELLOC" , "IF" , "Mean frequency of all spindles, per relative position within the spindle (five bins)" );
  
  hide_table( "SPINDLES", "F,CH,PHASE,RELLOC" , "Mean IF stratified by phase and relative location in spindle [if]" );
  hide_var( "SPINDLES" , "F,CH,PHASE,RELLOC" , "SOPL_CHIRP" , "Spindle chirp" );
  
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
  hide_param( "SPINDLES" , "if" , "" , "Estimate instantaneous frequency of spindles" );
  hide_param( "SPINDLES" , "if-frq" , "1" , "Window around target frequency (default 2 hz)" );
  hide_param( "SPINDLES" , "tlock" , "" , "Flag to request (verbose) average, peak-locked waveforms" );
  hide_param( "SPINDLES" , "verbose-coupling" , "" , "Add extra tables of EEG/CWT phase/time-locked to SO" );


  // show-coef verbose output

  add_table( "SPINDLES" , "F,CH,T" , "Verbose threshold/coefficient output [show-coeff]" );
  add_var( "SPINDLES" , "F,CH,T" , "RAWCWT" , "Raw CWT coefficient" );
  add_var( "SPINDLES" , "F,CH,T" ,"CWT" , "CWT coefficient" );
  add_var( "SPINDLES" , "F,CH,T" ,"CWT_TH" , "CWT primary threshold" );
  add_var( "SPINDLES" , "F,CH,T" ,"CWT_TH2" , "CWT secondary threshold" );
  add_var( "SPINDLES" , "F,CH,T" ,"CWT_THMAX" , "CWT maximum threshold" );
  


  //
  // SO (duplicated from SO command below)
  //


  add_param( "SPINDLES" , "so" , "" , "Detects slow oscillations and spindle/SO coupling" );
  
  add_param( "SPINDLES" , "mag" , "2" , "SO, relative mangitude threshold (times mean/median)" );
  add_param( "SPINDLES" , "uV-neg" , "-40" , "SO, absolute negative peak uV amplitude threshold" );
  add_param( "SPINDLES" , "uV-p2p" , "80" , "SO, absolute peak-to-peak uV amplitude threshold" );
    
  add_param( "SPINDLES" , "f-lwr" , "0.2" , "SO filter, lower transition frequency" );
  add_param( "SPINDLES" , "f-upr" , "4.5" , "SO filter, upper transition frequency" );
  
  add_param( "SPINDLES" , "t-lwr" , "0" , "SO, lower duration (secs)" );
  add_param( "SPINDLES" , "t-upr" , "3" , "SO, upper duration (secs)" );
  
  add_param( "SPINDLES" , "t-neg-lwr" , "0" , "SO, lower duration for negative peak (secs)" );
  add_param( "SPINDLES" , "t-neg-upr" , "1" , "SO, upper duration for negative peak (secs)" );
  
  hide_param( "SPINDLES" , "neg2pos" , "" , "SO, Use negative-to-positive zero crossings" );
  add_param( "SPINDLES" , "th-mean" , "" , "SO, use mean not median" );
  add_param( "SPINDLES" , "stats-median" , "" , "SO, use median (not mean) when reporting stats over SOs" );  
 
  add_table( "SPINDLES" , "CH" , "SO channel-level statistics" );
  add_var( "SPINDLES" , "CH" , "SO" , "Number of SO detected" );
  add_var( "SPINDLES" , "CH" , "SO_RATE" , "SO per minute" );
  add_var( "SPINDLES" , "CH" , "SO_AMP" , "SO amplitude (negative peak)" );
  add_var( "SPINDLES" , "CH" , "SO_P2P" , "SO peak-to-peak amplitude" );
  add_var( "SPINDLES" , "CH" , "SO_DUR" , "SO duration (secs)" );	   
  add_var( "SPINDLES" , "CH" , "SO_TRANS" , "SO transition (secs)" );
  add_var( "SPINDLES" , "CH" , "SO_TRANS_FREQ" , "SO transition freq (Hz)" );
  add_var( "SPINDLES" , "CH" , "SO_NEG_DUR" , "Negative peak SO duration (secs)" );
  add_var( "SPINDLES" , "CH" , "SO_POS_DUR", "Positive peak SO duration (secs)" );
  add_var( "SPINDLES" , "CH" , "SO_P2P" , "Peak-to-peak SO amplitude" );
  add_var( "SPINDLES" , "CH" , "SO_SLOPE_POS1" , "Positive peak rising slope" );
  add_var( "SPINDLES" , "CH" , "SO_SLOPE_POS2" , "Positive peak falling slope" );
  add_var( "SPINDLES" , "CH" , "SO_SLOPE_NEG1" , "Negative peak falling slope" );
  add_var( "SPINDLES" , "CH" , "SO_SLOPE_NEG2" , "Negative peak rising slope" );
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

  add_cmd( "transients" , "SO" , "Detect slow oscillations" );
  add_url( "SO" , "spindles-so/#so" );
  
  add_param( "SO" , "sig" , "C3,C4" , "Restrict analysis to these channels" );
  
  add_param( "SO" , "mag" , "2" , "Relative mangitude threshold (times mean/median)" );
  add_param( "SO" , "uV-neg" , "-40" , "Absolute negative peak uV amplitude threshold" );
  add_param( "SO" , "uV-p2p" , "80" , "Absolute peak-to-peak uV amplitude threshold" );
  
  add_param( "SO" , "mag" , "2" , "Relative mangitude threshold (times mean/median)" );
  
  add_param( "SO" , "f-lwr" , "0.2" , "Lower transition frequency" );
  add_param( "SO" , "f-upr" , "4.5" , "Upper transition frequency" );
  
  add_param( "SO" , "t-lwr" , "0" , "Lower duration (secs)" );
  add_param( "SO" , "t-upr" , "3" , "Upper duration (secs)" );
  
  add_param( "SO" , "t-neg-lwr" , "0" , "Lower duration for negative peak (secs)" );
  add_param( "SO" , "t-neg-upr" , "1" , "Upper duration for negative peak (secs)" );
  
  add_param( "SO" , "neg2pos" , "" , "Use negative-to-positive zero crossings" );
  add_param( "SO" , "th-mean" , "" , "Use mean not median" );
  add_param( "SO" , "stats-median" , "" , "Use median (not mean) when reporting stats over SOs" );  
 
  add_param( "SO" , "tl" , "C3" , "Output signal time-locked to detected SOs" );
  add_param( "SO" , "onset" , "" , "Sync to SO onset for tl option" );
  add_param( "SO" , "pos" , "" , "Sync to positive peak for tl option" );
  add_param( "SO" , "window" , "2" , "Specify window size (seconds) for tl option" );
  
  add_table( "SO" , "CH" , "Channel-level statistics" );
  add_var( "SO" , "CH" , "SO" , "Number of SO detected" );
  add_var( "SO" , "CH" , "SO_RATE" , "SO per minute" );
  add_var( "SO" , "CH" , "SO_AMP" , "SO amplitude (negative peak)" );
  add_var( "SO" , "CH" , "SO_P2P" , "SO peak-to-peak amplitude" );
  add_var( "SO" , "CH" , "SO_DUR" , "SO duration (secs)" );	   
  add_var( "SO" , "CH" , "SO_NEG_DUR" , "Negative peak duration (secs)" );
  add_var( "SO" , "CH" , "SO_POS_DUR", "Positive peak duration (secs)" );
  add_var( "SO" , "CH" , "SO_TRANS" , "SO transition (secs)" );	   
  add_var( "SO" , "CH" , "SO_TRANS_FREQ" , "SO transition freq (Hz)" );	   

  add_var( "SO" , "CH" , "SO_P2P" , "Peak-to-peak amplitude" );
  add_var( "SO" , "CH" , "SO_SLOPE_POS1" , "Positive peak rising slope" );
  add_var( "SO" , "CH" , "SO_SLOPE_POS2" , "Positive peak falling slope" );
  add_var( "SO" , "CH" , "SO_SLOPE_NEG1" , "Negative peak falling slope" );
  add_var( "SO" , "CH" , "SO_SLOPE_NEG2" , "Negative peak rising slope" );
  add_var( "SO" , "CH" , "SO_TH_NEG" , "Negative peak threshold [mag]" );
  add_var( "SO" , "CH" , "SO_TH_P2P" , "Peak-to-peak threshold [mag]" );

  add_table( "SO" , "CH,E" , "Epoch-level statistics" );
  add_var( "SO" , "CH,E" , "N" , "Number of SO detected" );
  add_var( "SO" , "CH,E" , "DOWN_AMP" , "Number of SO detected" );
  add_var( "SO" , "CH,E" , "UP_AMP" , "Number of SO detected" );
  add_var( "SO" , "CH,E" , "P2P_AMP" , "Number of SO detected" );
  add_var( "SO" , "CH,E" , "SLOPE_NEG1" , "Negative peak falling slope" );
  add_var( "SO" , "CH,E" , "SLOPE_NEG2" , "Negative peak rising slope" );
  add_var( "SO" , "CH,E" , "SLOPE_POS1" , "Positive peak rising slope" );
  add_var( "SO" , "CH,E" , "SLOPE_POS2" , "Positive peak falling slope" );
  
  add_table( "SO" , "CH,N" , "Per-SO statistics" );
  add_var( "SO" , "CH,N" , "DOWN_AMP" , "Negative peak amplitude" );
  add_var( "SO" , "CH,N" , "DOWN_IDX" , "Negative peak sample index" ); 
  add_var( "SO" , "CH,N" , "UP_AMP" , "Positive peak ampltiude" ); 
  add_var( "SO" , "CH,N" , "UP_IDX" , "Positive peak sample index" ); 
  add_var( "SO" , "CH,N" , "START" , "Start of SO (in seconds elapsed from start of EDF)" ); 
  add_var( "SO" , "CH,N" , "START_IDX" , "Start of SO (in sample-point units)" ); 
  add_var( "SO" , "CH,N" , "STOP" , "Stop of SO (in seconds elapsed from start of EDF)" ); 
  add_var( "SO" , "CH,N" , "STOP_IDX" , "Stop of SO (in sample-point units)" ); 
  add_var( "SO" , "CH,N" , "DUR" , "SO duration" ); 
  add_var( "SO" , "CH,N" , "DUR1" , "SO HW1 duration" ); 
  add_var( "SO" , "CH,N" , "DUR2" , "SO HW2 duration" ); 
  add_var( "SO" , "CH,N" , "TRANS" , "SO transition (sec)" );
  add_var( "SO" , "CH,N" , "P2P_AMP" , "SO peak-to-peak amplitude" ); 
  add_var( "SO" , "CH,N" , "SLOPE_POS1" , "Positive peak rising slope" ); 
  add_var( "SO" , "CH,N" , "SLOPE_POS2" , "Positive peak falling slope" ); 
  add_var( "SO" , "CH,N" , "SLOPE_NEG1" , "Negative peak falling slope" );
  add_var( "SO" , "CH,N" , "SLOPE_NEG2" , "Negative peak rising slope" );
  
  add_table( "SO" , "CH,CH2,SP" , "SO time-locked signal averaging [tl]" );
  add_var( "SO" , "CH,CH2,SP" , "SOTL" , "SO time-locked signal average" );

  
  
  /////////////////////////////////////////////////////////////////////////////////
  //
  // CROSS-SIGNAL
  //
  /////////////////////////////////////////////////////////////////////////////////

  //
  // COH
  //

  add_cmd( "topo" , "COH" , "Pairwise channel coherence" );
  add_url( "COH" , "cross-signal-analysis/#coh" );
  
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

  add_table( "COH" , "B,CH1,CH2,E" , "Epoch-level band coherence" );
  add_var( "COH" , "B,CH1,CH2,E" , "COH" , "Magnitude-squared coherence" );
  add_var( "COH" , "B,CH1,CH2,E" , "ICOH" , "Imaginary coherence" );
  add_var( "COH" , "B,CH1,CH2,E" , "LCOH" , "Lagged coherence" );
  
  add_table( "COH" , "CH1,CH2,E,F" , "Epoch-level coherence" );
  add_var( "COH" , "CH1,CH2,E,F" , "COH" , "Magnitude-squared coherence" );
  add_var( "COH" , "CH1,CH2,E,F" , "ICOH" , "Imaginary coherence" );
  add_var( "COH" , "CH1,CH2,E,F" , "LCOH" , "Lagged coherence" );

  // as these files can get large...
  set_compressed( "COH" , tfac_t( "CH1,CH2,B,E" ) );
  set_compressed( "COH" , tfac_t( "CH1,CH2,F,E" ) );
      
  //
  // PSI 
  //

  add_cmd( "topo" , "PSI" , "Phase slope index" );
  add_url( "PSI" , "cross-signal-analysis/#psi" );

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

  add_cmd( "topo" , "SYNC" , "Global phase synchrony" );
  add_url( "SYNC" , "cross-signal-analysis/#sync" );
  
  add_param( "SYNC" , "sig" , "C3,C4" , "Restrict analysis to these channel" );

  add_table( "SYNC" , "E,F" , "Epoch-wise analysis" );
  add_var( "SYNC" , "E,F" , "KOP" , "Magnitude-squared coherence" );
  set_compressed( "SYNC" , tfac_t( "E,F" ) );

  
  //
  // CORREL
  //

  add_cmd( "topo" , "CORREL" , "Pairwise signal correlation coefficients" );
  add_url( "CORREL" , "cross-signal-analysis/#correl" );
  
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

  add_table( "CORREL" , "CH,CH" , "Channel-level summaries of whole-signal correlations" );
  add_var( "CORREL" , "CH" , "SUMM_LOW", "Number of correlations below ch-low threshold" );
  add_var( "CORREL" , "CH" , "SUMM_HIGH", "Number of correlations aboive ch-high threshold" );
  
  add_var( "CORREL" , "CH" , "SUMM_MEAN", "Mean correlation for this channel" );
  add_var( "CORREL" , "CH" , "SUMM_MIN", "Min correlation for this channel" );
  add_var( "CORREL" , "CH" , "SUMM_MAX", "Max correlation for this channel" );
  
  add_table( "CORREL" , "CH1,CH2,E" , "Whole-signal correlations for pairs of channels" );
  add_var( "CORREL" , "CH1,CH2,E" , "R", "Pearson product moment correlation" );
  set_compressed( "CORREL" , tfac_t( "CH1,CH2,E" ) );

  
  //
  // MI
  //
  
  add_cmd( "topo" , "MI" , "Calculates pairwise mutual information metrics across channels" );
  add_url( "MI" , "cross-signal-analysis/#mi" );

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
  

  //
  // INTERPOLATE
  //
  
  
  
  /////////////////////////////////////////////////////////////////////////////////
  //
  // CFC
  //
  /////////////////////////////////////////////////////////////////////////////////

  
  //
  // CC
  //

  add_cmd( "topo" , "CC" , "Calculates dPAC and wPLI" );
  add_url( "CC" , "cc/#cc" );

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
  add_var( "CC" , "CH1,CH2,F1,F2" , "CFC" , "Cross-frequency coupling 0/1" );
  add_var( "CC" , "CH1,CH2,F1,F2" , "XCH" , "Cross-channel coupling 0/1" );
  add_var( "CC" , "CH1,CH2,F1,F2" , "dPAC" , "dPAC metric" );
  add_var( "CC" , "CH1,CH2,F1,F2" , "dPAC_Z" , "Z-normalized dPAC metric" );

  add_table( "CC" , "E,CH1,CH2,F1,F2" , "Epoch-level CC output" );
  add_var( "CC" , "E,CH1,CH2,F1,F2" , "CFC" , "Cross-frequency coupling 0/1" );
  add_var( "CC" , "E,CH1,CH2,F1,F2" , "CFC" , "Cross-frequency coupling 0/1" );
  add_var( "CC" , "E,CH1,CH2,F1,F2" , "XCH" , "Cross-channel coupling 0/1" );
  add_var( "CC" , "E,CH1,CH2,F1,F2" , "dPAC" , "dPAC metric" );
  add_var( "CC" , "E,CH1,CH2,F1,F2" , "dPAC_Z" , "Z-normalized dPAC metric" );
  set_compressed( "CC" , tfac_t( "E,CH1,CH2,F1,F2" ) );
  
  
  //
  // GLM
  //
  
  

  /////////////////////////////////////////////////////////////////////////////////
  //
  // MISC
  //
  /////////////////////////////////////////////////////////////////////////////////


  // HR     Estimate per-epoch heart rate from ECG
  // SPIKE  Create a synthetic signal by combining part of one signal with another
  // ZR     Calculate per-epoch Z-ratio


  /////////////////////////////////////////////////////////////////////////////////
  //
  // EXPERIMENTAL
  //
  /////////////////////////////////////////////////////////////////////////////////


  //
  // EXE
  //


  add_cmd( "exp" , "EXE" , "Epoch-by-epoch PDC-based clustering" );
  add_url( "EXE" , "exp/#exe" );

  add_param( "EXE" , "sig" , "C3,C4,F3,F4" , "Optionally specify channels (defaults to all)" );
  add_param( "EXE" , "uni" , "" , "For N signals, run N univariate analyses, rather than a single multi-signal one" );
  add_param( "EXE" , "representative" , "4" , "Extract N representative epochs" );
    
  add_param( "EXE" , "m" , "5" , "PDC embedding dimension" );
  add_param( "EXE" , "t" , "1" , "PDC span" );

  add_param( "EXE" , "k" , "10" , "Number of clusters" );
    
  add_table( "EXE" , "E,CH" , "Epoch cluster assignment" );
  add_var( "EXE" , "E,CH" , "CL" , "Cluster code [cluster]" );
  add_var( "EXE" , "E,CH" , "K" , "Representative split [representative]");
  add_var( "EXE" , "E,CH" , "KE" , "Representative epoch [representative]");

  add_table( "EXE" , "CH,K" , "Representative split info [representative]" );
  add_var( "EXE" , "CH,K" , "E" , "Representative epoch for split K" );
  add_var( "EXE" , "CH,K" , "N" , "Number of epochs in split K" );

  add_table( "EXE" , "E,CH" , "Epoch cluster assignment" );
  add_var( "EXE" , "E,CH" , "CL" , "Cluster code" );

  
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
      if ( ! hidden_cmd( *jj ) ) 
	ss << help( *jj , true , false ) ;
      ++jj;
    }
    ss << "\n";
    ++ii;
  }
  return ss.str();  
}

// all commands in a domain
std::string cmddefs_t::help_commands( const std::string & d ) const 
{
  std::stringstream ss;
  std::map<std::string,std::set<std::string> >::const_iterator ii = dcmds.find( d );
  if ( ii == dcmds.end() ) return "";
  const std::set<std::string> & c = ii->second;
  std::set<std::string>::const_iterator jj = c.begin();
  while ( jj != c.end() ) { 
    if ( ! hidden_cmd( *jj ) )
      ss << help( *jj , false , false ) ;
    ++jj;
  }
  return ss.str();    
}
  
// verbose describe commands [cmd] 

std::string cmddefs_t::help( const std::string & cmd , bool show_domain_label , bool verbose ) const
{
  if ( cmds.find( cmd ) == cmds.end() ) return "";
  if ( hidden_cmd( cmd ) ) return "";
	
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

	    if ( hidden_param( cmd , jj->first ) )
	      {
		++jj;
		continue;
	      }
	    
	    ss << "  " << std::left << std::setw( 12 ) << jj->first;
	    
  	    // example (if any)
	    std::string ex = px.find( cmd )->second.find( jj->first )->second ; 
	    if ( ex != "" )
	      {
		std::string msg = jj->first + "=" + ex;
		ss << std::left << std::setw( 20 ) << msg;
	      }
 	    else
 	      ss << std::left << std::setw(20) << " ";;
	    
	    // description
	    ss << std::left << std::setw( 12 ) << jj->second ;
	    
	    // requirements?
	    std::string req = preq.find( cmd )->second.find( jj->first )->second ; 
	    if ( req != "" )
	      ss << " [req. " << req << "]";
	    
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
	      
	      if ( hidden_table( cmd , tfac ) )
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
			if ( ! hidden_var( cmd , tfac , vv->first ) ) 
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



std::set<std::string> cmddefs_t::variables( const std::string & cmd ,  const param_t * param , const tfac_t & tfac  )
{
  // add param restriction on the variable list...
  std::set<std::string> r;
  std::map<std::string,std::map<tfac_t,std::map<std::string,std::string> > >::const_iterator ii = ovars.find( cmd );
  if ( ii == ovars.end() ) return r;
  const std::map<tfac_t,std::map<std::string,std::string> > & v2 = ii->second;
  std::map<tfac_t,std::map<std::string,std::string> >::const_iterator jj = v2.find( tfac );
  if ( jj == v2.end() ) return r;
  const std::map<std::string,std::string> & v3 = jj->second;
  std::map<std::string,std::string>::const_iterator kk = v3.begin();
  while ( kk != v3.end() )
    {
      r.insert( kk->first );
      ++kk;
    }
  return r;
}


// domain description
void cmddefs_t::add_domain( const std::string & domain , const std::string & label ,  const std::string & desc )
{
  domain_label[ domain ] = label;
  domain_desc[ domain ] = desc;
}

bool cmddefs_t::is_domain( const std::string & d ) 
{
  return domain_label.find( d ) != domain_label.end();
}

// command description 
void cmddefs_t::add_cmd( const std::string & domain , const std::string & cmd , const std::string & desc , const bool hide )
{
  dcmds[ domain ].insert( cmd );
  cmds[ cmd ] = desc ; 
  cdomain[ cmd ] = domain;
  chide[ cmd ] = hide ;
}

// hidden command description 
void cmddefs_t::hide_cmd( const std::string & domain , const std::string & cmd , const std::string & desc )
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

void cmddefs_t::add_note( const std::string & cmd , const std::string & note ) 
{
  if ( cmds.find( cmd ) == cmds.end() ) Helper::halt( cmd + " not registered" );
  cnotes[ cmd ] = note;
}

// parameters for this command
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
void cmddefs_t::hide_param( const std::string & cmd , const std::string & param , 
			    const std::string & ex ,  // "" if none
			    const std::string & desc , 
			    const std::string & requirements )
{
  add_param( cmd , param , ex , desc , requirements , true );
}


// output from this command , "CMD" , "F,B,CH,E" , "desc" , is compressed Y/N
void cmddefs_t::add_table( const std::string & cmd , const std::string & factors , const std::string & desc , bool isz , bool hide )
{
  tfac_t tfac( factors );
  otables[ cmd ][ tfac ] = desc ; 
  ofacs[ cmd ][ tfac ] = isz ; 
  ohide[ cmd ][ tfac ] = hide;
}

void cmddefs_t::hide_table( const std::string & cmd , const std::string & factors , const std::string & desc , bool isz )
{
  add_table( cmd , factors , desc , isz , true );
}

// add variable
void cmddefs_t::add_var( const std::string & cmd , const std::string & factors , const std::string & var , const std::string & desc , const bool hide )
{
  tfac_t tfac( factors );
  ovars[ cmd ][ tfac ][ var ] = desc;
  vhide[ cmd ][ tfac ][ var ] = hide;
}

// add hidden variable
void cmddefs_t::hide_var( const std::string & cmd , const std::string & factors , const std::string & var , const std::string & desc )
{
  add_var( cmd , factors , var , desc , true );
}



void cmddefs_t::all_compressed( bool b ) { allz = b; } 

bool cmddefs_t::all_compressed() const { return allz; }

void cmddefs_t::none_compressed( bool b ) { nonez = b; } 

bool cmddefs_t::none_compressed() const { return nonez; }


void cmddefs_t::add_tag( const std::string & tag ) { tags.insert( tag ); } 

void cmddefs_t::clear_tags() { tags.clear(); }

bool cmddefs_t::is_tag( const std::string & tag ) const { return tags.find( tag ) != tags.end(); } 



bool cmddefs_t::hidden_cmd( const std::string & c ) const
{
  std::map<std::string,bool>::const_iterator cc = chide.find( c );
  if ( cc == chide.end() ) return false;
  return cc->second;
}

bool cmddefs_t::hidden_param( const std::string & c , const std::string & p ) const
{
  std::map<std::string,std::map<std::string,bool> >::const_iterator cc = phide.find( c );
  if ( cc == phide.end() ) return false;
  std::map<std::string,bool>::const_iterator pp = cc->second.find( p );
  if ( pp == cc->second.end() ) return false;
  return pp->second;
}

bool cmddefs_t::hidden_table( const std::string & c , const tfac_t & tfac ) const
{
  std::map<std::string,std::map<tfac_t,bool> >::const_iterator cc = ohide.find( c );
  if ( cc == ohide.end() ) return false;
  std::map<tfac_t,bool>::const_iterator tt = cc->second.find( tfac );
  if ( tt == cc->second.end() ) return false;
  return tt->second;
}

bool cmddefs_t::hidden_var( const std::string & c , const tfac_t & tfac , const std::string & v ) const
{
  std::map<std::string,std::map<tfac_t,std::map<std::string,bool> > >::const_iterator cc = vhide.find( c );
  if ( cc == vhide.end() ) return false;
  std::map<tfac_t,std::map<std::string,bool> >::const_iterator tt = cc->second.find( tfac );
  if ( tt == cc->second.end() ) return false;
  std::map<std::string,bool>::const_iterator vv = tt->second.find( v );
  if ( vv == tt->second.end() ) return false;
  return vv->second;
}

//#pragma GCC pop_options
