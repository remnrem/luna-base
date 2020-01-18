
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

std::string nsrr_t::remap( const std::string & a )
{

  // check that remapping is enforced
  if ( ! globals::remap_nsrr_annots ) return a;
  
  // not found in map 
  if ( amap.find( a ) == amap.end() ) return a;

  // else return the remapped term
  return amap[ a ] ;
  
}

void nsrr_t::init() 
{

  //
  // auto-populated variables -- note these are defined here even if nsrr-remap mode if 
  // later turned off.  not a problem, as the variables will be overwritten by anything specified on
  // the command line, this will not yet have been parsed
  //
  
  cmd_t::vars[ "arousal" ] = "arousal_standard,arousal_spontaneous,arousal_external,arousal_respiratory,arousal_plm,arousal_cheshire";
  
  cmd_t::vars[ "apnea" ] = "apnea_obstructive,apnea_central,apnea_mixed,hypopnea";
  
  cmd_t::vars[ "artifact" ] 
    = "artifact_respiratory,artifact_proximal_pH,artifact_distal_pH,artifact_blood_pressure,artifact_TcCO2,artifact_SpO2,artifact_EtCO2";
  
  cmd_t::vars[ "arrhythmia" ] = "bradycardia,tachycardia,tachycardia_narrowcomplex";
	       
  cmd_t::vars[ "plm" ] = "plm_left,plm_right";
  
  cmd_t::vars[ "n1" ] = "NREM1";
  
  cmd_t::vars[ "n2" ] = "NREM2";
  
  cmd_t::vars[ "n3" ] = "NREM3,NREM4";
  
  cmd_t::vars[ "rem" ] = "REM";
  
  cmd_t::vars[ "wake" ] = "wake";
  
  cmd_t::vars[ "sleep" ] = "NREM1,NREM2,NREM3,NREM4,REM";
  
  
  //
  // actual remapping terms
  //
  
  amap[ "Arousal ()" ] = "arousal_standard";
  amap[ "Arousal|Arousal ()" ] = "arousal_standard";
  amap[ "Arousal|Arousal (Arousal)" ] = "arousal_standard";
  amap[ "Arousal|Arousal (Standard)" ] = "arousal_standard";
  amap[ "Arousal|Arousal (STANDARD)" ] = "arousal_standard";
  amap[ "ASDA arousal|Arousal (ADSA)" ] = "arousal_standard";
  amap[ "ASDA arousal|Arousal (asda)" ] = "arousal_standard";
  amap[ "ASDA arousal|Arousal (ASDA)" ] = "arousal_standard";
  amap[ "Arousal (Asda)" ] = "arousal_standard";
  amap[ "Arousal (ASDA)" ] = "arousal_standard";

  amap[ "Arousal (ARO SPONT)" ] = "arousal_spontaneous";
  amap[ "Spontaneous arousal|Arousal (apon aro)" ] = "arousal_spontaneous" ;
  amap[ "Spontaneous arousal|Arousal (ARO SPONT)" ] = "arousal_spontaneous" ;
  amap[ "Spontaneous arousal|Arousal (spon aro)" ] = "arousal_spontaneous" ;
  amap[ "Spontaneous arousal|Arousal (SPON ARO)" ] = "arousal_spontaneous";

  amap[ "External arousal|Arousal (External Arousal)" ] = "arousal_external";
  
  amap[ "Arousal resulting from respiratory effort|Arousal (ARO RES)" ] = "arousal_respiratory";
  amap[ "Arousal resulting from respiratory effort|Arousal (RESP ARO)" ] = "arousal_respiratory";
  amap[ "RERA" ] = "arousal_respiratory";
  amap[ "Arousal (ARO RES)" ] = "arousal_respiratory";
  amap[ "Respiratory effort related arousal|RERA" ] = "arousal_respiratory";
      
  amap[ "Arousal resulting from Chin EMG|Arousal (Cheshire)" ] = "arousal_cheshire";
  amap[ "Arousal resulting from Chin EMG|Arousal (CHESHIRE)" ] = "arousal_cheshire";
      
  amap[ "Arousal (ARO Limb)" ] = "arousal_plm";
  amap[ "Arousal resulting from periodic leg movement|Arousal (PLM)" ] = "arousal_plm";
  amap[ "Arousal resulting from periodic leg movement|Arousal (PLM ARO)" ] = "arousal_plm";
    
  // apnea and hypopnea
  
  amap [ "Obstructive Apnea" ] = "apnea_obstructive" ;
  amap[ "Obstructive apnea|Obstructive Apnea" ] = "apnea_obstructive" ;
  
  amap[ "Central Apnea" ] = "apnea_central";
  amap[ "Central apnea|Central Apnea" ] = "apnea_central";
  
  amap[ "Mixed Apnea" ] = "apnea_mixed";
  amap[ "Mixed apnea|Mixed Apnea" ] = "apnea_mixed";
  
  amap[ "Hypopnea" ] = "hypopnea";
  amap[ "Hypopnea|Hypopnea" ] = "hypopnea";

  // note -- this is some kind of apnea event... need to update
  amap[ "Unsure" ] = "unsure";
  amap[ "Unsure|Unsure" ] = "unsure";
  
  
  // sleep staging
      
  amap[ "NREM1" ] = "NREM1";
  amap[ "Stage 1 sleep|1" ] = "NREM1";
      
  amap[ "NREM2" ] = "NREM2";
  amap[ "Stage 2 sleep|2" ] = "NREM2";
  
  amap[ "NREM3" ] = "NREM3";
  amap[ "Stage 3 sleep|3" ] = "NREM3";
  
  amap[ "NREM4" ] = "NREM4" ;
  amap[ "Stage 4 sleep|4" ] = "NREM4";
  
  amap[ "REM" ] = "REM";
  amap[ "REM sleep|5" ] = "REM";
  
  amap[ "Unscored" ] = "unscored";
  amap[ "Unscored|9" ] = "unscored";
  
  amap[ "Wake" ] = "wake";
  amap[ "Wake|0" ] = "wake";
  
  amap[ "Movement|6" ] = "movement";
  
  
  
  // limb movements
  
  amap[ "Periodic leg movement - left|PLM (Left)" ] = "plm_left";
  amap[ "Periodic leg movement - right|PLM (Right)" ] = "plm_right";
  amap[ "PLM (Left)" ] = "plm_left";
  amap[ "PLM (Right)" ] = "plm_right";
  
  amap[ "Limb Movement (Left)" ] = "plm_left";
  amap[ "Limb movement - left|Limb Movement (Left)" ] = "plm_left";
  amap[ "Limb Movement (Right)" ] = "plm_right";      
  amap[ "Limb movement - right|Limb Movement (Right)" ] = "plm_right";
      
      
  // artifacts 
  amap[ "Respiratory artifact|Respiratory artifact" ] = "artifact_respiratory";

  amap[ "Respiratory artifact" ] = "artifact_respiratory";
      
  amap[ "Proximal pH artifact" ] = "artifact_proximal_pH";

  amap[ "Distal pH artifact" ] = "artifact_distal_pH";

  amap[ "Blood pressure artifact" ] = "artifact_blood_pressure";
  
  amap[ "TcCO2 artifact" ] = "artifact_TcCO2" ;       
  
  amap[ "SpO2 artifact" ] = "artifact_SpO2";
     
  amap[ "SpO2 artifact|SpO2 artifact" ] = "artifact_SpO2";

  amap[ "EtCO2 artifact" ] = "artifact_EtCO2";
  amap[ "EtCO2 artifact|EtCO2 artifact" ] = "artifact_EtCO2";
      
  // respiratory

  amap[ "Periodic Breathing" ] = "periodic_breathing";
  amap[ "Periodic breathing|Periodic Breathing" ] = "periodic_breathing";

  amap[ "Respiratory Paradox" ] = "respiratory_paradox";

  // desats

  amap[ "SpO2 desaturation" ] = "desat";
  
  amap[ "SpO2 desaturation|SpO2 desaturation" ] = "desat";
      

  // cardiac


  amap[ "Bradycardia" ] = "bradycardia";

  amap[ "Tachycardia" ] = "tachycardia";

  amap[ "Narrow Complex Tachycardia"] = "tachycardia_narrowcomplex";      

  amap[ "Narrow complex tachycardia|Narrow Complex Tachycardia" ] = "tachycardia_narrowcomplex";
  
  // misc

  amap[ "Technician Notes" ] = "notes";
  

}

