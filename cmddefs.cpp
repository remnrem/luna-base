

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


#include "cmddefs.h"

#include "luna.h"
#include <iomanip>

extern globals global;

void cmddefs_t::init()
{

  //
  // base URL
  //

  url_root = "http://zzz.bwh.harvard.edu/luna/ref/";


  //
  // Domains
  //

  add_domain( "summ"     , "Summaries"       , "Basic summary commands" );
  add_domain( "annot"    , "Annotations"     , "Adding and displaying annotations" );
  add_domain( "expr"     , "Expressions"     , "Evaluating more advanced annotation-based expressions" );
  add_domain( "epoch"    , "Epochs"          , "Epoching signals and epoch-level annotations" );
  add_domain( "mask"     , "Masks"           , "Masking epochs based on annotations and other criteria" );
  add_domain( "manip"    , "Manipulations"   , "Manipulating signal data" );
  add_domain( "output"   , "Outputs"         , "Commands to output signals in different formats" );
  add_domain( "filter"   , "FIR filters"     , "FIR filter design and application" );
  add_domain( "artifact" , "Artifacts"       , "Artifacts detection/correction routines" );
  add_domain( "hypno"    , "Hypnograms"      , "Characterizations of hypnograms" );
  add_domain( "power"    , "Power spectra"   , "Power spectral density estimation" );
  add_domain( "spindles" , "Spindles and SO" , "Spindles and slow oscillations" );
  add_domain( "topo"     , "Cross-signal analyses"    , "Coherence and other topographical analyses" );
  add_domain( "cfc"      , "Cross-frequency coupling" , "Phase-amplitude and spindle/SO coupling" );
  add_domain( "misc"     , "Misc"           , "Misc. commands" );
  add_domain( "exp"      , "Experimental"   , "Experimental features, under heavy development / for internal use only" );

  
  //
  // Commands, parameters, output tables and variables 
  //
  

  //
  // DESC
  //

  add_cmd( "summ" , "DESC" , "Simple description of an EDF, sent to the console" );



  //
  // SUMMARY
  //

  add_cmd( "summ" , "SUMMARY" , "More verbose description, sent to the console" );


  //
  // HEADERS
  //
  
  add_cmd( "summ" , "HEADERS" , "Tabulate (channel-specific) EDF header information" );

  add_table( "HEADERS" , "." , "Basic EDF header information" );
  add_var( "HEADERS" , "." , "NR" , "Number of records" );
  add_var( "HEADERS" , "." , "NS" , "Number of signals/channels" );
  add_var( "HEADERS" , "." , "REC.DUR" , "Duration of each record (seconds)" );
  add_var( "HEADERS" , "." , "TOT.DUR.SEC" , "Total duration of EDF (seconds)" );
  add_var( "HEADERS" , "." , "TOT.DUR.HMS" , "Total duration of EDF (hh:mm:ss string)" );

  add_table( "HEADERS" , "CH" , "Per-channel header information" );
  add_var( "HEADERS" , "CH" , "DMAX" , "Digital max" );
  add_var( "HEADERS" , "CH" , "DMIN" , "Digital min" );
  add_var( "HEADERS" , "CH" , "PDIM", "Physical dimension" );
  add_var( "HEADERS" , "CH" , "PMAX", "Physical min" );
  add_var( "HEADERS" , "CH" , "PMIN", "Physical max" );
  add_var( "HEADERS" , "CH" , "SR", "Sample rate (Hz)" );


  //
  // TAG
  //

  add_cmd( "summ" , "TAG" , "Generic command to add a tag (level/factor) to the output" );  
  add_param( "TAG" , ""    , "RUN/L1" , "" , "Add tag with level L1 to factor RUN in output" );
  add_param( "TAG" , "tag" , "RUN/L1" , "" , "Identical to the above, but explicitly using the tag option" );

  //
  // STATS
  //

  add_cmd( "summ"   , "STATS" , "Basic signal statistics (min/max, mean, RMS, etc)" );
  add_param( "STATS" , "sig" , "C3,C4" , "" , "Restrict analysis to these channels" );
  add_param( "STATS" , "epoch" , "" , "" , "Calculate per-epoch statistics" );
  
  add_table( "STATS" , "CH" , "Whole-night, per-channel statistics, based on all epochs" );
  add_var( "STATS" , "CH" , "MIN" , "Signal minimum (from data, not EDF header)" );
  add_var( "STATS" , "CH" , "MAX" , "Signal maximum (from data, not EDF header)" );
  add_var( "STATS" , "CH" , "MEAN" , "Signal mean" );
  add_var( "STATS" , "CH" , "MEDIAN" , "Signal median" );
  add_var( "STATS" , "CH" , "RMS" , "Signal root mean square" );

  add_var( "STATS" , "CH" , "NE" ,  "Total number of epochs in record [epoch]" );
  add_var( "STATS" , "CH" , "NE1" , "Number of unmasked epochs actually used in calculations [epoch]" );
  add_var( "STATS" , "CH" , "MEDIAN.MEAN" , "Median of all per-epoch means [epoch]" );
  add_var( "STATS" , "CH" , "MEDIAN.MEDIAN" , "Median of all per-epoch medians [epoch]" );
  add_var( "STATS" , "CH" , "MEDIAN.RMS" , "Median of all per-epoch RMS [epoch]" );

  add_table( "STATS" , "CH,E" , "Per-epoch, per-channel statistics for unmasked epochs only" );
  add_var( "STATS" , "CH,E" , "MIN" , "Signal minimum (from data, not EDF header)" );
  add_var( "STATS" , "CH,E" , "MAX" , "Signal maximum (from data, not EDF header)" );
  add_var( "STATS" , "CH,E" , "MEAN" , "Signal mean" );
  add_var( "STATS" , "CH,E" , "MEDIAN" , "Signal median" );
  add_var( "STATS" , "CH,E" , "RMS" , "Signal root mean square" );

    
  /////////////////////////////////////////////////////////////////////////////////
  //
  // ANNOTATIONS
  //
  /////////////////////////////////////////////////////////////////////////////////
  
  //
  // --xml
  //

  add_cmd( "annot" , "--xml" , "Quickly view an NSRR XML annotation file" );


  //
  // ANNOTS
  //
  
  add_cmd( "annot" , "ANNOTS" , "Tabulate all annotations" );
  add_url( "ANNOTS" , "annotations/#annots" );
  //  add_note( "ANNOTS" , "Any formatted free text goes here,\n shown at end of verbose help link\n");

  add_param( "ANNOTS" , "epoch" , "" , "" , "Show epoch-level summaries" );
  add_param( "ANNOTS" , "show-masked" , "" , "" , "Show masked annotations (default is not to do so)" );
  add_param( "ANNOTS" , "any" , "" , "" , "Keep annotations that have any overlap with one or more unmasked epochs (default)" );
  add_param( "ANNOTS" , "all" , "" , "" , "Only keep annotations that are completely within unmasked epochs" );
  add_param( "ANNOTS" , "start" , "" , "" , "Keep annotations that start in an unmasked epoch" );

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
  add_var( "ANNOTS" , "E,INTERVAL,INST" , "ANNOT_MASK" , "Annotation instance mask status (1=masked/excluded) [epoch]" );
  add_var( "ANNOTS" , "E,INTERVAL,INST" , "EPOCH_MASK" , "Epoch mask status (1=masked/excluded) [epoch]" );



  /////////////////////////////////////////////////////////////////////////////////
  //
  // EPOCHS
  //
  /////////////////////////////////////////////////////////////////////////////////
  

  //
  // Domain: 
  //

  add_cmd( "epoch" , "EPOCH" , "Set epochs" );
  add_url ( "EPOCH" , "epochs/#epoch" );

  add_table( "EPOCH" , "" , "Epoch-level summaries" );
  add_var( "EPOCH" , "" , "DUR" , "Epoch duration (seconds)" );
  add_var( "EPOCH" , "" , "INC" , "Epoch increment (seconds)" );
  add_var( "EPOCH" , "" , "NE" , "Number of epochs" );

  

  /////////////////////////////////////////////////////////////////////////////////
  //
  // MASKS
  //
  /////////////////////////////////////////////////////////////////////////////////


//   MASKMask epochs based on annotations and other features

//     ifif=N2Mask N2 epochs; unmask non-N2 epochs
// 			     ifnotifnot=N2Mask non-N2 epochs; unmask N2 epochs
// 								mask-ifmask-if=N2Mask N2 epochs; leave non-N2 epochs as they are
// 												   mask-ifnotmask-ifnot=N2Mask non-N2 epochs, leave N2 epochs as they are
// 												   unmask-ifunmask-if=N2Unmask N2 epochs; leave non-N2 epochs as they are
// 																	    unmask-ifnotunmask-ifnot=N2Unmask non-N2 epochs; leave N2 epochs as they are


// 																							       exprexpr="if(annot1) && annot2.v2 > 0.95"Mask/unmask epochs for which the expression is true/false
// 																																       not-exprnot-expr="if(annot1) && annot2.v2 > 0.95"Unmask/mask epochs for which the expression is true/false
// 																																										       mask-exprmask-expr="if(annot1) && annot2.v2 > 0.95"Mask epochs for which the expression is true
// 																																																				  unmask-exprunmask-expr="if(annot1) && annot2.v2 > 0.95"Unmask epochs for which the expression is true


// 																																																														   epochepoch=20-50Mask epochs outside of this range
// 																																																														   mask-epochmask-epoch=20-50As above, except mask epochs inside this range
// 																																																														   secsec=0-60Set mask to include all epochs that span the interval from 0 to 60 seconds (i.e. these are unmasked, all other epochs are masked)
// 																																																														   hmshms=8:00-9:00Set mask to include all epochs that span the interval from 8am to 9am


// 																																																														   noneMASK noneSet to include all epochs (default when first loading a new EDF)
// 																																																														   allMASK allSet to exclude (i.e. mask) all epochs
// 																																																														   flipMASK flipFlip all epoch masks
// 																																																														   randomMASK random=50Select up to 50 from currently unmasked epochs
// 																																																														   leadingMASK leading=wakeRemove all epochs up to the first epoch without this annotation
// 																																																														   flankedMASK flanked=N2,2Include only N2 epochs flanked by at least 2 other N2 epochs


// 																																																														   N_MATCHESNumber of epochs that match the condition (e.g. having annotation A)
// 																																																														   N_MASK_SETNumber of previously unmasked epochs that were masked by this operation
// 																																																														   N_MASK_UNSETNumber of previously masked epochs that were unmasked by this operation
// 																																																														   N_UNCHANGEDNumber of epochs whose mask status was not changed by this operation
// 																																																														   N_RETAINEDNumber of epochs retained after this operation
// 																																																														   N_TOTALTotal number of epochs




//     DUMP-MASKOutput epoch-level mask information

// 																																																														   Out: 
// Epoch-level tabulation (strata: E)


//       EPOCH_MASK Mask status: 0 is unmasked (included), and 1 is masked (i.e. excluded)



//     RESTRUCTURE (or RE)Remove masked out epochs (and channels)


  /////////////////////////////////////////////////////////////////////////////////
  //
  // MANIPULATIONS
  //
  /////////////////////////////////////////////////////////////////////////////////


  add_cmd( "manip"   , "RE" , "Restructure an EDF (drop channels/epochs)" );
  add_table( "RE" , "" , "Restructured data duration" );
  add_var( "RE" , "" , "DUR1" , "Duration pre-restructuring (secs)");
  add_var( "RE" , "" , "DUR2" , "Duration post-restructuring (secs)");
  add_var( "RE" , "" , "NR1" , "Duration pre-restructuring (records)");
  add_var( "RE" , "" , "NR2" , "Duration post-restructuring (records)");

  add_cmd( "manip"   , "RESTRUCTURE" , "Restructure an EDF (drop channels/epochs)" );
  add_table( "RESTRUCTURE" , "" , "Restructured data duration" );
  add_var( "RESTRUCTURE" , "" , "DUR1" , "Duration pre-restructuring (secs)");
  add_var( "RESTRUCTURE" , "" , "DUR2" , "Duration post-restructuring (secs)");
  add_var( "RESTRUCTURE" , "" , "NR1" , "Duration pre-restructuring (records)");
  add_var( "RESTRUCTURE" , "" , "NR2" , "Duration post-restructuring (records)");


// SIGNALSRetain/remove specific EDF channels

//   drop drop=EMG,ECGDrop channels EMG and ECG
//   keep keep=C3,C4Drop all channels except C3 and C4


// COPYDuplicate one or more EDF channels

//   sig=C3,C4List of channels to duplicate
//   tagtag=DELTAA required option, this is added to make the new channel name, e.g. C3 becomes C3_DELTA


// RESAMPLEResample signal(s)

//   sig sig=C3,C4Signal list
//   sr sr=100New sampling rate (Hz)


// REFERENCERe-reference signals

//   sigsig=C3,C4Signal(s) to re-reference
//   refref=A1,A2Signal(s) to provide the reference


// uVRescale units to uV

//   sigsig=C3,C4Signal(s) to convert

// mVRescale units to mV

//   sigsig=C3,C4Signal(s) to convert

// FLIPFlip polarity of signal

//   sigsig=C3,C4Signals to flip

// TIME-TRACK Add a time-track to an EDF
//   // no param

// RECORD-SIZEChange EDF record size

//   durdur=1New EDF record/block size
//   edf-diredf-dir=edfs/Folder for writing new EDFs
//   edf-tagedf-tag=rec1Tag added to new EDFs
//   sample-listsample-list=s2.lstGenerate a sample-list pointing to the new EDFs

// ANONStrip ID information from EDF header

// 					   // no param


  /////////////////////////////////////////////////////////////////////////////////
  //
  // OUTPUTS
  //
  /////////////////////////////////////////////////////////////////////////////////



  /////////////////////////////////////////////////////////////////////////////////
  //
  // FILTERS
  //
  /////////////////////////////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////////////////////////
  //
  // ARTIFACTS
  //
  /////////////////////////////////////////////////////////////////////////////////



  /////////////////////////////////////////////////////////////////////////////////
  //
  // HYPNOGRAMS
  //
  /////////////////////////////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////////////////////////
  //
  // SPECTRAL
  //
  /////////////////////////////////////////////////////////////////////////////////


  //
  // PSD
  //

  add_cmd( "power"   , "PSD" , "Power spectral density estimation (Welch)" );
  add_param( "PSD" , "sig" , "C3,C4" , "" , "Restrict analysis to these channels" );
  add_param( "PSD" , "epoch" , "" , "" , "Calculate per-epoch statistics" );
  add_param( "PSD" , "max" , "100" , "" , "Calculate per-epoch statistics" );
  add_param( "PSD" , "spectrum" , "" , "" , "Calculate per-epoch statistics" );
  add_param( "PSD" , "epoch-spectrum" , "" , "" , "Calculate per-epoch statistics" );
  
  add_table( "PSD" , "CH" , "Number of epochs" );
  add_var( "PSD" , "CH" , "NE" , "Number of epochs" );

  add_table( "PSD" , "CH,B" , "Whole-night, per-channel band power" );
  add_var( "PSD" , "CH,B" , "PSD" , "Power" );
  add_var( "PSD" , "CH,B" , "RELPSD" , "Relative power" );

  add_table( "PSD" , "CH,F" , "Whole-night, per-channel power" );
  add_var( "PSD" , "CH,F" , "PSD" , "Power" );

  add_table( "PSD" , "CH,B,E" , "Whole-night, per-channel per-epoch band power" );
  add_var( "PSD" , "CH,B,E" , "PSD" , "Power" );
  add_var( "PSD" , "CH,B,E" , "RELPSD" , "Relative power" );

  add_table( "PSD" , "CH,F,E" , "Whole-night, per-channel per-epoch power" );
  add_var( "PSD" , "CH,F,E" , "PSD" , "Power" );
  set_compressed( "PSD" , tfac_t( "CH,F,E" ) );



  //
  // MTM
  //

  add_cmd( "power"   , "MTM" , "Power spectral density estimation (Welch)" );
  add_param( "MTM" , "sig" , "C3,C4" , "" , "Restrict analysis to these channels" );
  add_param( "MTM" , "epoch" , "" , "" , "Calculate per-epoch statistics" );
  add_param( "MTM" , "max" , "100" , "" , "Calculate per-epoch statistics" );
  add_param( "MTM" , "dB" , "" , "" , "Decibel scale output" );
  add_param( "MTM" , "spectrum" , "" , "" , "Calculate per-epoch statistics" );
  add_param( "MTM" , "epoch-spectrum" , "" , "" , "Calculate per-epoch statistics" );
  
  add_table( "MTM" , "CH" , "Number of epochs" );
  add_var( "MTM" , "CH" , "NE" , "Number of epochs" );

  add_table( "MTM" , "CH,F" , "Whole-night, per-channel power" );
  add_var( "MTM" , "CH,F" , "MTM" , "Power" );

  add_table( "MTM" , "CH,F,E" , "Whole-night, per-channel per-epoch power" );
  add_var( "MTM" , "CH,F,E" , "MTM" , "Power" );
  set_compressed( "MTM" , tfac_t( "CH,F,E" ) );


  /////////////////////////////////////////////////////////////////////////////////
  //
  // SPINDLES/SO
  //
  /////////////////////////////////////////////////////////////////////////////////



  /////////////////////////////////////////////////////////////////////////////////
  //
  // CROSS-SIGNAL
  //
  /////////////////////////////////////////////////////////////////////////////////



  /////////////////////////////////////////////////////////////////////////////////
  //
  // CFC
  //
  /////////////////////////////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////////////////////////
  //
  // MISC
  //
  /////////////////////////////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////////////////////////
  //
  // EXPERIMENTAL
  //
  /////////////////////////////////////////////////////////////////////////////////



}





tfac_t::tfac_t( const std::string & s ) { 
  std::vector<std::string> tok = Helper::parse( s , "," );
  for (int i=0;i<tok.size();i++) 
    {
      if ( tok[i][0] != '_' && ! globals::cmddefs.is_tag( tok[i] ) )
	fac.insert( tok[i] );
    } 
}

cmddefs_t::cmddefs_t()
{
  init();
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
      ss << help( *jj , false ) ;
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
    ss << help( *jj , false ) ;
    ++jj;
  }
  return ss.str();    
}
  
// verbose describe commands [cmd] 

std::string cmddefs_t::help( const std::string & cmd , bool verbose ) const
{
  if ( cmds.find( cmd ) == cmds.end() ) return "";
  std::stringstream ss;
  if ( ! verbose ) 
    {
      ss << std::left << std::setw( 18 ) << domain_label.find( cdomain.find( cmd )->second )->second  << " " 
	 << std::left << std::setw( 12 ) << cmd << " "
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

	    ss << "  " << std::left << std::setw( 12 ) << jj->first;

	    // requirements?
	    std::string req = px.find( cmd )->second.find( jj->first )->second ; 
	    if ( req != "" )
	      ss << std::left << std::setw(12) << req;
	    else
	      ss << std::left << std::setw(12) << " ";;

	    // example (if any)
	    std::string ex = px.find( cmd )->second.find( jj->first )->second ; 
	    if ( ex != "" )
	      ss << " e.g. " << jj->first << "=" << ex;
	    ss << " " << jj->second << "\n";
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
	      
	      ss << "   " << std::left << std::setw( 24 ) << tfac.as_string( " x " ) 
		 << ii->second << "\n";
	      
	      ss << "   " << std::left << std::string( 60 , '-' )
		<< "\n";

	      // dump as compressed text ?
	      bool tdump = allz || ofacs.find( cmd )->second.find( tfac )->second;
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
    }

  
  return ss.str();
  
}



//
// test if a table exists
//

bool cmddefs_t::exists( const std::string & cmd ,
			const tfac_t & tfac ) const
{

  if ( cmds.find( cmd ) == cmds.end() ) return false;
  return ofacs.find( cmd )->second.find( tfac ) != ofacs.find( cmd )->second.end() ;
}


//
// Output?
//

bool cmddefs_t::out_compressed( const std::string & cmd , const tfac_t & tfac ) const
{
  // all output is in compressed plain-text
  if ( allz ) return true;
  
  // cmd not found
  if ( cmds.find( cmd ) == cmds.end() ) return false;
  
  // table not found
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

