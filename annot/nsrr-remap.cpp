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

namespace {

struct remap_definition_t {
  std::string primary;
  std::vector<std::string> aliases;
  bool strict;
};

// Keep definitions independent of the label settings in effect at registration.
std::vector<remap_definition_t> definitions;
bool lookup_dirty = true;
bool lookup_sanitize = false;
bool lookup_replace_spaces = false;
char lookup_space_replacement = '_';

std::string normalize_label( const std::string & label )
{
  std::string normalized = Helper::trim( label );
  if ( globals::replace_annot_spaces )
    normalized = Helper::search_replace( normalized , ' ' , globals::space_replacement );
  if ( globals::sanitize_everything )
    normalized = Helper::trim( globals::replace_annot_spaces
                              ? Helper::sanitize( normalized )
                              : Helper::sanitize( normalized , ' ' ) , '_' );
  return Helper::squash( Helper::squash( normalized , ' ' ) , '_' );
}

std::string raw_label( const std::string & label )
{
  return Helper::repair_label_commas( Helper::unquote( label ) );
}

}

void nsrr_t::ensure_mappings()
{
  if ( ! lookup_dirty && lookup_sanitize == globals::sanitize_everything
       && lookup_replace_spaces == globals::replace_annot_spaces
       && lookup_space_replacement == globals::space_replacement ) return;

  amap.clear();
  bmap.clear();
  pmap.clear();
  std::map<std::string,std::string> alias_sources, primary_sources;

  for ( const auto & definition : definitions )
    {
      const std::string primary = normalize_label( definition.primary );
      const std::string key = Helper::toupper( primary );
      const std::string raw_primary = Helper::toupper( definition.primary );

      if ( primary_sources.count( key ) && primary_sources[key] != raw_primary )
        Helper::halt( "annotation normalization collision between primary labels: "
                      + primary_sources[key] + " & " + definition.primary );
      primary_sources[key] = raw_primary;

      if ( definition.strict )
        {
          if ( pmap.count( key ) && pmap[key] != primary )
            Helper::halt( "inconsistent case in remaps for primary: " + pmap[key] + " & " + primary );
          if ( amap.count( key ) )
            Helper::halt( primary + " specified as both primary annotation and mapped term" );
        }
      pmap[key] = primary;

      for ( const auto & alias : definition.aliases )
        {
          const std::string mapped = Helper::toupper( normalize_label( alias ) );
          const std::string raw_alias = Helper::toupper( alias );
          if ( definition.strict && mapped == key ) continue;

          if ( definition.strict && pmap.count( mapped ) )
            Helper::halt( mapped + " specified as both primary annotation and mapped term" );

          if ( amap.count( mapped ) && Helper::toupper( amap[mapped] ) != key )
            {
              // add() deliberately replaces an exact alias (e.g. Movement).
              // Different raw aliases must not silently collide after normalization.
              if ( definition.strict || alias_sources[mapped] != raw_alias )
                Helper::halt( "annotation normalization collision for alias: " + alias
                              + " (" + amap[mapped] + " & " + primary + ")" );
            }
          amap[mapped] = primary;
          alias_sources[mapped] = raw_alias;
        }
    }

  // Derive the reverse map from the final aliases so overrides leave no stale entries.
  for ( const auto & alias : amap )
    bmap[ Helper::toupper( alias.second ) ].push_back( alias.first );

  lookup_sanitize = globals::sanitize_everything;
  lookup_replace_spaces = globals::replace_annot_spaces;
  lookup_space_replacement = globals::space_replacement;
  lookup_dirty = false;
}

std::string nsrr_t::remap( const std::string & label )
{
  ensure_mappings();
  const std::string normalized = normalize_label( label );
  const std::string key = Helper::toupper( normalized );
  const auto primary = pmap.find( key );
  if ( primary != pmap.end() ) return unmapped ? "" : primary->second;
  const auto alias = amap.find( key );
  if ( alias != amap.end() ) return unmapped ? "" : alias->second;
  return whitelist ? "" : normalized;
}

void nsrr_t::annot_remapping( const std::string & s )
{
  // canonical|alias1|alias2; user rules retain the no-transitive-alias rule.
  const std::vector<std::string> tok = Helper::quoted_parse( s , "|" );
  if ( tok.size() < 2 )
    Helper::halt( "bad format for annotation remapping:  canonical|alias 1|alias 2\n" + s );
  remap_definition_t definition;
  definition.primary = raw_label( tok[0] );
  definition.strict = true;
  if ( definition.primary.find( "|" ) != std::string::npos )
    Helper::halt( "primary annotation labels cannot contain pipe (|) characters" );
  for ( size_t j = 1; j < tok.size(); ++j )
    definition.aliases.push_back( raw_label( tok[j] ) );
  definitions.push_back( definition );
  lookup_dirty = true;
}

void nsrr_t::add( const std::string & primary , const std::string & alias )
{
  definitions.push_back( { raw_label( primary ) , { raw_label( alias ) } , false } );
  lookup_dirty = true;
}

void nsrr_t::clear()
{
  definitions.clear();
  amap.clear();
  bmap.clear();
  pmap.clear();
  lookup_dirty = true;
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
  add( "N1" , "Sleep stage 1" ); 
  add( "N1" , "Stage N1" );
  add( "N1" , "StageN1" );
  add( "N1" , "Stage 1" );
  add( "N1" , "Stage1" );
  add( "N1" , "Stage NREM1" ); 
  add( "N1" , "StageNREM1" ); 
  add( "N1" , "S1" );
  
  add( "N2" , "NREM2" );
  add( "N2" , "NREM2 sleep" );
  add( "N2" , "N2 sleep" );
  add( "N2" , "Stage 2 sleep|2" ); 
  add( "N2" , "Sleep stage N2" );
  add( "N2" , "Sleep stage 2" );
  add( "N2" , "Stage N2" );
  add( "N2" , "StageN2" );
  add( "N2" , "Stage 2" );
  add( "N2" , "Stage2" );
  add( "N2" , "Stage NREM2" ); 
  add( "N2" , "StageNREM2" ); 
  add( "N2" , "S2" ); 
  
  add( "N3" , "NREM3" );
  add( "N3" , "NREM3 sleep" );
  add( "N3" , "N3 sleep" ); 
  add( "N3" , "Stage 3 sleep|3" ); 
  add( "N3" , "Sleep stage N3" );
  add( "N3" , "Sleep stage 3" ); 
  add( "N3" , "Stage N3" );
  add( "N3" , "StageN3" );
  add( "N3" , "Stage 3" );
  add( "N3" , "Stage3" );
  add( "N3" , "Stage NREM3" ); 
  add( "N3" , "StageNREM3" ); 
  add( "N3" , "S3" );

  // NNREM4 --> N3
  add( "N3" , "N4" ); 
  add( "N3" , "NREM4" );
  add( "N3" , "NREM4 sleep" );
  add( "N3" , "N4 sleep" ); 
  add( "N3" , "Stage 4 sleep|4" ); 
  add( "N3" , "Sleep stage N4" );
  add( "N3" , "Sleep stage 4" ); 
  add( "N3" , "Stage N4" );
  add( "N3" , "StageN4" );
  add( "N3" , "Stage 4" );
  add( "N3" , "Stage4" );
  add( "N3" , "Stage NREM4" ); 
  add( "N3" , "S4" );

  add( "R" , "REM" );
  add( "R" , "REM sleep" ); 
  add( "R" , "REM sleep|5" );
  add( "R" , "REM|5" );
  add( "R" , "Sleep stage REM" );
  add( "R" , "Sleep stage R" ); 
  add( "R" , "Stage R" );
  add( "R" , "StageR" );
  add( "R" , "Stage REM" );
  add( "R" , "StageREM" );
  
  add( "W" , "Wake" ); 
  add( "W" , "Wake|0" ); 
  add( "W" , "Sleep stage Wake" );
  add( "W" , "Sleep stage W" );
  add( "W" , "Stage W" );
  add( "W" , "StageW" );
  add( "W" , "Stage Wake" );
  add( "W" , "Wake stage" ); 

  add( "NR" , "Sleep stage N" ); 
  add( "NR" , "Sleep stage NR" ); 
  add( "NR" , "Sleep stage NREM" );
  add( "NR" , "Stage N" ); 
  add( "NR" , "Stage NR" ); 
  add( "NR" , "Stage NREM" ); 
  add( "NR" , "StageN" ); 
  add( "NR" , "StageNR" ); 
  add( "NR" , "StageNREM" ); 
  add( "NR" , "NREM" );
  add( "NR" , "NREM sleep" );
  add( "NR" , "NR sleep" ); 

  // non-sleep stage
  add( "U" , "Unscored" ); 
  add( "U" , "Unscored|9" ); 

  add( "?" , "Unknown" ); 
  add( "?" , "Sleep stage ?" ); 
  add( "?" , "Stage ?" );
  add( "?" , "missing" );

  add( "M" , "Movement" ); 
  add( "M" , "Movement|6" ); 

  // lights 
  add( "L" , "Lights" );
  add( "L" , "lights" ); 
  add( "L" , "Light" ); 

  add( "lights_on" , "Lights On" ); 
  add( "lights_on" , "LightsOn" ); 

  add( "lights_off" , "Lights Off" );
  add( "lights_off" , "LightsOff" ); 

  // gaps
  add( "G", "GAP" );
  add( "G", "-" );
  
  
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

 // some additional edge-case NSRR-specific stage remappings
 
 add( "W" , "SRO:Wake" );
 add( "N1" , "SRO:Stage1Sleep" );
 add( "N2" , "SRO:Stage2Sleep" );
 add( "N3" , "SRO:Stage3Sleep" );
 add( "N3" , "SRO:Stage4Sleep" );
 add( "N3" , "SRO:Stage34Sleep" );
 add( "R" , "SRO:RapidEyeMovement" );

 add( "W", "SDO:WakeState" );
 add( "N1", "SDO:NonRapidEyeMovementSleep-N1" );
 add( "N2", "SDO:NonRapidEyeMovementSleep-N2" );
 add( "N3", "SDO:NonRapidEyeMovementSleep-N3" );
 add( "N3", "SDO:NonRapidEyeMovementSleep-N4" );
 add( "R", "SDO:RapidEyeMovementSleep" );

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


