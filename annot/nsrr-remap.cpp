
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


#include "nsrr-remap.h"
#include "defs/defs.h"
#include "eval.h"

std::string nsrr_t::remap( const std::string & s1 )
{
  
  // always trim obvious whitespace
  std::string a0 = Helper::trim( s1 );
  
  // swap internal spaces for a different character?
  std::string a1 = globals::replace_annot_spaces
    ? Helper::search_replace( a0 , ' ' , globals::space_replacement ) 
    : a0 ;
   
  // santization (but perhaps allowing for spaces?)
  // and then retrim on sanitized characters
  if ( globals::sanitize_everything )
    {
      if ( globals::replace_annot_spaces )
	a1 = Helper::trim( Helper::sanitize( a1 ) , '_' ) ;
      else // allow spaces in a sanitized version still
	a1 = Helper::trim( Helper::sanitize( a1 , ' ' ) , '_' ) ;
    }
  
  //
  // reduce multiple (internal) spaces or underscores to one
  //
  
  a1 = Helper::squash( Helper::squash( a1 , ' ' ) , '_' );
  
  //
  // do nothing
  //
  // if ( ! do_remap )
  //   {
  //     return a1;
  //   }

  //
  // do remapping
  //
  
  std::string a = a1;
  
  std::string a_uc = Helper::toupper( a );


  //
  // found as a primary? ( return preferred case in pmap )
  //
    
  if ( pmap.find( a_uc ) != pmap.end() )
    {  
      if ( unmapped ) return "";
      else return pmap[ a_uc ];
    }

  //
  // found as an alias?
  //
  
  bool found = amap.find( a_uc ) != amap.end() ;

  //
  // not found as is (or map is empty)
  //
  
  if ( ! found )
    {

      std::string ca = a;
      
      // we may still want to swap spaces
      // note - this is prob redundant now, as we
      // will always replace spaces above..
      
      if ( globals::replace_annot_spaces )
	ca = Helper::search_replace( ca , ' ' , globals::space_replacement );
      
      std::string ca_uc = Helper::toupper( ca );
      
      // recheck that space-swapped version does not have an alias
      found = amap.find( ca_uc ) != amap.end() ;

      if ( ! found )
	{
	  // if in white-list mode, we return nothing (i.e. an invalid/skipped annot)
	  if ( whitelist ) return "";

	  // else leave case of original unchanged and return 
	  return ca ; 
	}
      else
	{
	  // skipping whitelisted terms?
	  if ( unmapped ) return "";

	  return amap[ ca_uc ]; // if using as as a key, must be UC

	}
    }

  // skipping whitelisted terms?
  if ( unmapped ) return "";

  // else return the remapped term (which may contain spaces, if that
  // was so requested...)
  return amap[ a_uc ] ;
  
}


void nsrr_t::annot_remapping( const std::string & s )
{
  // same structure for annotation re-mappings as for signal channels;
  // code copied from cmd_t::signal_alias()
  
  // the primary alias can occur multiple times, and have multiple 
  // labels that are mapped to it; matches are case-insensitive, as
  // all aliases and checks stored in upper-case
  
  // however: two rules
  // 1. many-to-one mapping means the same label cannot have multiple primary aliases
  // 2. following, and on principle of no transitive properties, alias cannot have alias
  
  // X|Y|Z
  // X|A|B
  
  // W|A  bad, A already mapped
  // V|X  bad, X already mapped
  // i.e. things can only occur once in the RHS, or multiple times in the LHS
  

  // format canonical|alias1|alias2 , etc.
  std::vector<std::string> tok = Helper::quoted_parse( s , "|" );    
  if ( tok.size() < 2 ) Helper::halt( "bad format for annotation remapping:  canonical|alias 1|alias 2\n" + s );
  
  const std::string primary = Helper::unquote( tok[0] );
  const std::string uc_primary = Helper::toupper( primary );

  // check that no | characters.. will cause problems downstream, so avoid here
  if ( primary.find( "|" ) != std::string::npos )
    Helper::halt( "primary annotation labels cannot contain pipe (|) characters" );
  
  
  
  // store primary mapping (but only on first occurrence)
  // i.e. the primary might be listed w/ variable case in the remap file
  // behavior here is to take the first; error if different
  if ( pmap.find( uc_primary ) == pmap.end() )
    {
      pmap[ uc_primary ] = primary;
    }
  else
    {
      if ( pmap[ uc_primary ] != primary )
	Helper::halt( "inconsistent case in remaps for primary: " + pmap[ uc_primary ] + " & " + primary ); 
    }
    
  // primary already specified
  if ( amap.find( uc_primary ) != amap.end() )
    Helper::halt( primary + " specified as both primary annotation and mapped term" );

  // add mappings
  for (int j=1;j<tok.size();j++) 
    {
      
      // impose rules
      // optional swap spaces for another character too
	  const std::string mapped = globals::sanitize_everything 
	    ? Helper::sanitize(  Helper::toupper( Helper::unquote( tok[j] ) ) )
	    : ( globals::replace_annot_spaces
		? Helper::search_replace( Helper::toupper( Helper::unquote( tok[j] ) ) , ' ' , globals::space_replacement )
		: Helper::toupper( Helper::unquote( tok[j] ) )
		);		

      // if same, skip
      if ( mapped == uc_primary ) continue;
      
      // otherwise, check no circular definitions
      if ( bmap.find( mapped ) != bmap.end() )
	Helper::halt( mapped + " specified as both primary annotation and mapped term" );
      
      if ( amap.find( mapped ) != amap.end() )
	if ( Helper::toupper( primary ) != Helper::toupper( amap[ mapped ] ) )
	  Helper::halt( mapped + " specified twice in alias file w/ different primary remaping" );
      
      // otherwise, set 
      amap[ mapped ] = primary;
      
      bmap[ uc_primary ].push_back( mapped );
    }
    
}


// add a new annotation remap
void nsrr_t::add( const std::string & p , const std::string & a )
{

  const std::string p2 = globals::sanitize_everything
    ? Helper::sanitize(  Helper::unquote( p ) )
    : ( globals::replace_annot_spaces
	? Helper::search_replace( Helper::unquote( p )  , ' ' , globals::space_replacement )
	: Helper::unquote( p ) );
	
  const std::string a2 = globals::sanitize_everything
    ? Helper::sanitize(  Helper::unquote( a ) )
    : ( globals::replace_annot_spaces
	? Helper::search_replace( Helper::unquote( a )  , ' ' , globals::space_replacement )
	: Helper::unquote( a ) );

  // swap spaces for another char?
  // const std::string p2 = globals::replace_annot_spaces ? Helper::search_replace( p , ' ' , globals::space_replacement ) : p;
  // const std::string a2 = globals::replace_annot_spaces ? Helper::search_replace( a , ' ' , globals::space_replacement ) : a;

  //std::cout << " mapping " << a2 << " " << p2 << "\n";
  amap[ Helper::toupper( a2 ) ] = p2;
  bmap[ Helper::toupper( p2 ) ].push_back( Helper::toupper( a2 ) );
  pmap[ Helper::toupper( p2 ) ] = p2;
}


// clear existing NSRR annotations
void nsrr_t::clear()
{
  amap.clear();
  bmap.clear();
  pmap.clear();  
}


void nsrr_t::init() 
{

  // just one for now (given remapping of ANNOTS as below) 
  // ${sleep}
  
  cmd_t::vars[ "sleep" ] = "N1,N2,N3,R";
  
  add( "N1" , "NREM1" );
  add( "N1" , "NREM1 sleep" );
  add( "N1" , "N1 sleep" ); 
  add( "N1" , "Stage 1 sleep|1" ); 
  add( "N1" , "Sleep stage N1" ); 
  add( "N1" , "Stage N1" );
  add( "N1" , "Stage NREM1" ); 
  
  add( "N2" , "NREM2" );
  add( "N2" , "NREM2 sleep" );
  add( "N2" , "N2 sleep" );
  add( "N2" , "Stage 2 sleep|2" ); 
  add( "N2" , "Sleep stage N2" );
  add( "N2" , "Stage N2" );
  add( "N2" , "Stage NREM2" ); 
  
  add( "N3" , "NREM3" );
  add( "N3" , "NREM3 sleep" );
  add( "N3" , "N3 sleep" ); 
  add( "N3" , "Stage 3 sleep|3" ); 
  add( "N3" , "Sleep stage N3" ); 
  add( "N3" , "Stage N3" );
  add( "N3" , "Stage NREM3" ); 

  // NNREM4 --> N3
  add( "N3" , "N4" ); 
  add( "N3" , "NREM4" );
  add( "N3" , "NREM4 sleep" );
  add( "N3" , "N4 sleep" ); 
  add( "N3" , "Stage 4 sleep|4" ); 
  add( "N3" , "Sleep stage N4" ); 
  add( "N3" , "Stage N4" );
  add( "N3" , "Stage NREM4" ); 

  add( "R" , "REM" );
  add( "R" , "REM sleep" ); 
  add( "R" , "REM sleep|5" ); 
  add( "R" , "Sleep stage R" ); 
  add( "R" , "Stage R" );
  add( "R" , "Stage REM" );
  
  add( "W" , "Wake" ); 
  add( "W" , "Wake|0" ); 
  add( "W" , "Sleep stage W" );
  add( "W" , "Stage W" );
  add( "W" , "Stage Wake" );
  add( "W" , "Wake stage" ); 

  add( "NR" , "Sleep stage N" ); 
  add( "NR" , "Sleep stage NREM" ); 
  add( "NR" , "NREM" );
  add( "NR" , "NREM sleep" );
  add( "NR" , "NR sleep" ); 

  add( "U" , "Unscored" ); 
  add( "U" , "Unscored|9" ); 

  add( "?" , "Unknown" ); 
  add( "?" , "Sleep stage ?" ); 
  add( "?" , "Stage ?" );
  
  add( "M" , "Movement|6" ); 

  add( "L" , "Lights" ); 

  add( "lights_on" , "Lights On" ); 
  add( "lights_on" , "LightsOn" ); 

  add( "lights_off" , "Lights Off" );
  add( "lights_off" , "LightsOff" ); 
    
  //
  // EDF+ annotations
  //   - by default, take sleep stages (remapped) and arousals as class-level vars
  //   - can be over-ridden with edf-annot-class=X,Y,Z
  
  edf_annot_class( "N1,N2,N3,R,W,?,arousal,LM,NR" );
  
}


void nsrr_t::init_nsrr_mappings() 
{
  
  //
  // auto-populated variables -- note these are defined here even if nsrr-remap mode if 
  // later turned off.  not a problem, as the variables will be overwritten by anything specified on
  // the command line, this will not yet have been parsed
  //
    
  //
  // actual remapping terms
  //

  // this should be kept fairly update-to-date w/ the NSRR NAP file harm.annots
  // assuming harm.annots only has one line per row, the following will generate the text below:
  // but do not include sleep stage terms (i.e. they are added above, separaetly)
  
  // grep ^remap resources/harm.annots | cut -f2 | cut -d"|" -f1  > tmp.1
  // grep ^remap resources/harm.annots | cut -f2 | cut -d"|" -f2- | tr -d '"' > tmp.2  
  // paste tmp.1 tmp.2 | awk -F"\t" ' { printf "add( \"" $1 "\" , \""  $2  "\" ); \n" } '
      
add( "arousal" , "Arousal ()" ); 
add( "arousal" , "Arousal|Arousal ()" );
add( "arousal" , "Arousal|Arousal" ); 
add( "arousal" , "Arousal|Arousal (Standard)" ); 
add( "arousal" , "Arousal_(STANDARD)" ); 
add( "arousal" , "Arousal|Arousal_(Arousal)" ); 
add( "arousal" , "ASDA arousal|Arousal (ADSA)" ); 
add( "arousal" , "ASDA arousal|Arousal (ASDA)" ); 
add( "arousal" , "Arousal (ASDA)" ); 
add( "arousal" , "Arousal_(Asda)" ); 
add( "arousal" , "EEG arousal" );
add( "arousal:spontaneous" , "Arousal (ARO SPONT)" ); 
add( "arousal:spontaneous" , "Spontaneous arousal|Arousal (apon aro)" ); 
add( "arousal:spontaneous" , "Spontaneous arousal|Arousal (ARO SPONT)" ); 
add( "arousal:spontaneous" , "Spontaneous arousal|Arousal (SPON ARO)" ); 
add( "arousal:respiratory" , "Arousal resulting from respiratory effort|Arousal (ARO RES)" ); 
add( "arousal:respiratory" , "RERA" ); 
add( "arousal:respiratory" , "Arousal (ARO RES)" ); 
add( "arousal:respiratory" , "Arousal resulting from respiratory effort|Arousal (RESP ARO)" ); 
add( "arousal:respiratory" , "Respiratory effort related arousal|RERA" ); 
add( "arousal:external" , "External arousal|Arousal (External Arousal)" ); 
add( "arousal:external" , "Arousal_(External_Arousal)" ); 
add( "arousal:cheshire" , "Arousal resulting from Chin EMG|Arousal (Cheshire)" ); 
add( "arousal:cheshire" , "Arousal_(CHESHIRE)" ); 
add( "arousal:lm" , "arousal_lm" ); 
add( "arousal:lm" , "lml_arousal" ); 
add( "arousal:lm" , "lmr_arousal" ); 
add( "arousal:lm" , "lmb_arousal" ); 
add( "arousal:lm" , "Arousal_(ARO_Limb)" ); 
add( "arousal:plm" , "arousal_plm" ); 
add( "arousal:plm" , "Arousal_resulting_from_periodic_leg_movement|Arousal_(PLM)" ); 
add( "arousal:plm" , "Arousal_resulting_from_periodic_leg_movement|Arousal_(PLM_ARO)" ); 
add( "apnea" , "Apnea" ); 
add( "apnea:obstructive" , "Obstructive apnea|Obstructive Apnea" ); 
add( "apnea:obstructive" , "Obstructive Apnea" ); 
add( "apnea:obstructive" , "apnea_obstructive" ); 
add( "apnea:obstructive" , "Obstructive_apnea|APNEA-OBSTRUCTIVE" ); 
add( "apnea:central" , "Central Apnea" ); 
add( "apnea:central" , "apnea_central" ); 
add( "apnea:central" , "Central apnea|Central Apnea" ); 
add( "apnea:central" , "Central_apnea|APNEA-CENTRAL" ); 
add( "apnea:mixed" , "Mixed Apnea" ); 
add( "apnea:mixed" , "apnea_mixed" ); 
add( "apnea:mixed" , "Mixed apnea|Mixed Apnea" ); 
add( "apnea:mixed" , "Mixed apnea|APNEA-MIXED" ); 
add( "hypopnea" , "Hypopnea|Hypopnea" ); 
add( "hypopnea:obstructive" , "hypopnea_obstructive" ); 
add( "hypopnea:obstructive" , "Obstructive_Hypopnea" ); 
add( "hypopnea:central" , "hypopnea_central" ); 
add( "periodic_breathing" , "Periodic Breathing" ); 
add( "periodic_breathing" , "Periodic breathing|Periodic Breathing" ); 
add( "respiratory_paradox" , "Respiratory Paradox" ); 
add( "snoring" , "Snoring" ); 
add( "cheynestokes_breathing" , "cheynestokes_breathing" ); 
add( "desat" , "SpO2 desaturation" ); 
add( "desat" , "SpO2 desaturation|SpO2 desaturation" ); 
add( "desat" , "SpO2 desaturation|DESAT" ); 
add( "unsure" , "Unsure|Unsure|Unsure" ); 
add( "movement" , "Movement" ); 
add( "PLM" , "Periodic leg movement" ); 
add( "PLM" , "Periodic leg movement|PLM" ); 
add( "PLM:left" , "Periodic leg movement - left|PLM (Left)" ); 
add( "PLM:left" , "PLM (Left)" ); 
add( "PLM:right" , "Periodic leg movement - right|PLM (Right)" ); 
add( "PLM:right" , "PLM (Right)" ); 
add( "LM" , "Limb Movement" ); 
add( "LM" , "Limb movement|Limb Movement" ); 
add( "LM:left" , "Limb Movement (Left)" ); 
add( "LM:left" , "Limb movement - left|Limb Movement (Left)" ); 
add( "LM:right" , "Limb Movement (Right)" ); 
add( "LM:right" , "Limb movement - right|Limb Movement (Right)" ); 
add( "artifact" , "Signal artifact|SIGNAL-ARTIFACT" ); 
add( "artifact:respiratory" , "Respiratory artifact" ); 
add( "artifact:respiratory" , "Respiratory artifact|Respiratory artifact" ); 
add( "artifact:proximal_pH" , "Proximal pH artifact" ); 
add( "artifact:proximal_pH" , "Proximal_pH_artifact|Proximal_pH_artifact" ); 
add( "artifact:distal_pH" , "Distal pH artifact" ); 
add( "artifact:pH" , "Proximal_pH|Distal_pH_artifact" ); 
add( "artifact:blood_pressure" , "Blood pressure artifact" ); 
add( "artifact:blood_pressure" , "Blood_pressure_artifact|Blood_pressure_artifact" ); 
add( "artifact:TcCO2" , "TcCO2 artifact" ); 
add( "artifact:TcCO2" , "TcCO2 artifact|TcCO2 artifact" ); 
add( "artifact:SpO2" , "SpO2 artifact" ); 
add( "artifact:SpO2" , "SpO2 artifact|SpO2 artifact" ); 
add( "artifact:EtCO2" , "EtCO2 artifact" ); 
add( "artifact:EtCO2" , "EtCO2 artifact|EtCO2 artifact" ); 
add( "artifact:body_temperature" , "Body_temperature_artifact|Body_temperature_artifact" ); 
add( "position:left" , "Body position change to left|POSITION-LEFT" ); 
add( "position:right" , "Body position change to right|POSITION-RIGHT" ); 
add( "position:prone" , "Body position change to prone|POSITION-PRONE" ); 
add( "position:supine" , "Body position change to supine|POSITION-SUPINE" ); 
add( "position:upright" , "Body position change to upright|POSITION-UPRIGHT" ); 
add( "arrhythmia:bradycardia" , "Bradycardia" ); 
add( "arrhythmia:bradycardia" , "Bradycardia|Bradycardia" ); 
add( "arrhythmia:tachycardia" , "Tachycardia" ); 
add( "arrhythmia:tachycardia" , "Tachycardia|Tachycardia" ); 
add( "arrhythmia:narrow_complex_tachycardia" , "Narrow Complex Tachycardia" ); 
add( "arrhythmia:narrow_complex_tachycardia" , "Narrow complex tachycardia|Narrow Complex Tachycardia" ); 
add( "notes" , "Technician Notes" ); 
  
  // END of auto-generated code (from NSRR/NAP harm.annots)


}


void nsrr_t::edf_annot_class( const std::string & s )
{

  if ( s == "*" )
    {
      all_edf_class = true;
      return;
    }
  
  edf_class.clear();
  std::vector<std::string> tok = Helper::parse( s , "," );
  for (int i=0; i<tok.size(); i++)
    edf_class.insert( tok[i] );
}

bool nsrr_t::as_edf_class(  const std::string & s )
{
  return all_edf_class || edf_class.find( s ) != edf_class.end() ;
}


