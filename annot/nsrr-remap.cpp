
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

// annotation map
std::map<std::string,std::string> nsrr_t::amap;
std::map<std::string,std::vector<std::string> > nsrr_t::bmap;
std::map<std::string,std::string> nsrr_t::pmap;

// annot list acts as an annot white list
bool nsrr_t::whitelist = false;
bool nsrr_t::unmapped = false;

std::string nsrr_t::remap( const std::string & a )
{
  
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
    pmap[ uc_primary ] = primary;
  else
    {
      if ( pmap[ uc_primary ] != primary )
	Helper::halt( "inconsistent case in remaps for primary: " + pmap[ uc_primary ] + " & " + primary ); 
    }
    
  // primary already specified
  if ( amap.find( uc_primary ) != amap.end() )
    Helper::halt( primary + " specified as both primary annotation and mapped term" );
  
  for (int j=1;j<tok.size();j++) 
    {
      
      // impose rules
      // optional swap spaces for another character too
      const std::string mapped = globals::replace_annot_spaces ?
	Helper::search_replace( Helper::toupper( Helper::unquote( tok[j] ) ) , ' ' , globals::space_replacement ) :
	Helper::toupper( Helper::unquote( tok[j] ) ) ;

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
  // swap spaces for another char?
  const std::string p2 = globals::replace_annot_spaces ? Helper::search_replace( p2 , ' ' , globals::space_replacement ) : p;
  const std::string a2 = globals::replace_annot_spaces ? Helper::search_replace( a2 , ' ' , globals::space_replacement ) : a;

  amap[ Helper::toupper( a2 ) ] = p2;
  bmap[ Helper::toupper( p2 ) ].push_back( Helper::toupper( a2 ) );
  pmap[ Helper::toupper( p2 ) ] = p2;
}

// clear all existing 
void nsrr_t::clear()
{
  amap.clear();
  bmap.clear();
  pmap.clear();
}


void nsrr_t::init() 
{
  
  //
  // auto-populated variables -- note these are defined here even if nsrr-remap mode if 
  // later turned off.  not a problem, as the variables will be overwritten by anything specified on
  // the command line, this will not yet have been parsed
  //

  // just one for now (given remapping of ANNOTS as below) 
  // ${sleep}

  cmd_t::vars[ "sleep" ] = "N1,N2,N3,R";
  
  
  //
  // actual remapping terms
  //

  // this should be kept fairly update-to-date w/ the NSRR NAP file harm.annots
  // assuming harm.annots only has one line per row, the following will generate the text below:
  
  // grep ^remap resources/harm.annots | cut -f2 | cut -d"|" -f1  > tmp.1
  // grep ^remap resources/harm.annots | cut -f2 | cut -d"|" -f2- | tr -d '"' > tmp.2  
  // paste tmp.1 tmp.2 | awk -F"\t" ' { printf "add( \"" $1 "\" , \""  $2  "\" ); \n" } '
      

  add( "arousal" , "Arousal ()" ); 
  add( "arousal" , "Arousal|Arousal ()" ); 
  add( "arousal" , "Arousal|Arousal (Standard)" ); 
  add( "arousal" , "ASDA arousal|Arousal (ADSA)" ); 
  add( "arousal" , "ASDA arousal|Arousal (ASDA)" ); 
  add( "arousal" , "Arousal (ASDA)" ); 
  add( "arousal/spontaneous" , "Arousal (ARO SPONT)" ); 
  add( "arousal/spontaneous" , "Spontaneous arousal|Arousal (apon aro)" ); 
  add( "arousal/spontaneous" , "Spontaneous arousal|Arousal (ARO SPONT)" ); 
  add( "arousal/spontaneous" , "Spontaneous arousal|Arousal (SPON ARO)" ); 
  add( "arousal/external" , "External arousal|Arousal (External Arousal)" ); 
  add( "arousal/RERA" , "Arousal resulting from respiratory effort|Arousal (ARO RES)" ); 
  add( "arousal/RERA" , "RERA" ); 
  add( "arousal/RERA" , "Arousal (ARO RES)" ); 
  add( "arousal/RERA" , "Arousal resulting from respiratory effort|Arousal (RESP ARO)" ); 
  add( "arousal/RERA" , "Respiratory effort related arousal|RERA" ); 
  add( "arousal/cheshire" , "Arousal resulting from Chin EMG|Arousal (Cheshire)" ); 
  add( "apnea/obstructive" , "Obstructive apnea|Obstructive Apnea" ); 
  add( "apnea/obstructive" , "Obstructive Apnea" ); 
  add( "apnea/central" , "Central Apnea" ); 
  add( "apnea/central" , "Central apnea|Central Apnea" ); 
  add( "apnea/mixed" , "Mixed Apnea" ); 
  add( "apnea/mixed" , "Mixed apnea|Mixed Apnea" ); 
  add( "apnea/mixed" , "Mixed apnea|APNEA-MIXED" ); 
  add( "hypopnea" , "Hypopnea" );
  add( "hypopnea" , "Hypopnea|Hypopnea" ); 
  add( "periodic_breathing" , "Periodic Breathing" ); 
  add( "periodic_breathing" , "Periodic breathing|Periodic Breathing" ); 
  add( "respiratory_paradox" , "Respiratory Paradox" ); 
  add( "desat" , "SpO2 desaturation" ); 
  add( "desat" , "SpO2 desaturation|SpO2 desaturation" ); 
  add( "desat" , "SpO2 desaturation|DESAT" ); 
  add( "unsure" , "Unsure|Unsure|Unsure" ); 
  add( "N1" , "NREM1" ); 
  add( "N1" , "Stage 1 sleep|1" ); 
  add( "N2" , "NREM2" ); 
  add( "N2" , "Stage 2 sleep|2" ); 
  add( "N3" , "NREM3" ); 
  add( "N3" , "Stage 3 sleep|3" ); 
  add( "N3" , "N4" ); 
  add( "N3" , "NREM4" ); 
  add( "N3" , "Stage 4 sleep|4" ); 
  add( "R" , "REM" ); 
  add( "R" , "REM sleep|5" ); 
  add( "W" , "Wake" ); 
  add( "W" , "Wake|0" ); 
  add( "U" , "Unscored" ); 
  add( "U" , "Unscored|9" ); 
  add( "?" , "Unknown" ); 
  add( "M" , "movement" ); 
  add( "M" , "Movement|6" ); 
  add( "L" , "Lights" ); 
  add( "L" , "Lights On" ); 
  add( "L" , "LightsOn" ); 
  add( "PLM" , "Periodic leg movement" ); 
  add( "PLM" , "Periodic leg movement|PLM" ); 
  add( "PLM/left" , "Periodic leg movement - left|PLM (Left)" ); 
  add( "PLM/right" , "Periodic leg movement - right|PLM (Right)" ); 
  add( "PLM/right" , "PLM (Right)" ); 
  add( "LM" , "Limb Movement" ); 
  add( "LM" , "Limb movement|Limb Movement" ); 
  add( "LM/left" , "Limb Movement (Left)" ); 
  add( "LM/left" , "Limb movement - left|Limb Movement (Left)" ); 
  add( "LM/right" , "Limb Movement (Right)" ); 
  add( "LM/right" , "Limb movement - right|Limb Movement (Right)" ); 
  add( "artifact" , "Signal artifact|SIGNAL-ARTIFACT" ); 
  add( "artifact/respiratory" , "Respiratory artifact" ); 
  add( "artifact/respiratory" , "Respiratory artifact|Respiratory artifact" ); 
  add( "artifact/proximal_pH" , "Proximal pH artifact" ); 
  add( "artifact/distal_pH" , "Distal pH artifact" ); 
  add( "artifact/blood_pressure" , "Blood pressure artifact" ); 
  add( "artifact/TcCO2" , "TcCO2 artifact" ); 
  add( "artifact/TcCO2" , "TcCO2 artifact|TcCO2 artifact" ); 
  add( "artifact/SpO2" , "SpO2 artifact" ); 
  add( "artifact/SpO2" , "SpO2 artifact|SpO2 artifact" ); 
  add( "artifact/EtCO2" , "EtCO2 artifact" ); 
  add( "artifact/EtCO2" , "EtCO2 artifact|EtCO2 artifact" ); 
  add( "position/left" , "Body position change to left|POSITION-LEFT" ); 
  add( "position/right" , "Body position change to right|POSITION-RIGHT" ); 
  add( "position/prone" , "Body position change to prone|POSITION-PRONE" ); 
  add( "position/supine" , "Body position change to supine|POSITION-SUPINE" ); 
  add( "position/upright" , "Body position change to upright|POSITION-UPRIGHT" ); 
  add( "arrhythmia/bradycardia" , "Bradycardia" ); 
  add( "arrhythmia/tachycardia" , "Tachycardia" ); 
  add( "arrhythmia/narrow_complex_tachycardia" , "Narrow Complex Tachycardia" ); 
  add( "arrhythmia/narrow_complex_tachycardia" , "Narrow complex tachycardia|Narrow Complex Tachycardia" ); 
  add( "notes" , "Technician Notes" ); 
  
  // END of auto-generated code (from NSRR/NAP harm.annots)
  
}

