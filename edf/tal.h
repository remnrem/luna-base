
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

#ifndef __TAL_H__
#define __TAL_H__

#include <stdint.h>
#include <vector>
#include <ostream>
#include <string>

struct edf_t;

struct interval_t;

struct tal_element_t 
{

tal_element_t( double onset = 0 , 
	       double duration = 0 , 
	       const std::string & name = "" ) 
  : onset(onset) , duration( duration ) , name(name) { }
  
  double onset;
  
  double duration;
  
  std::string name;  

  friend std::ostream & operator<<( std::ostream & out , const tal_element_t & t );

};

struct tal_t 
{
  
  tal_t( edf_t * edf , int channel , int r );
  
  void decode( const std::string & str );
  
  std::string encode() const; 

  void add( const tal_element_t & t ) { d.push_back( t ); } 

  int size() const { return d.size(); } 

  int record;
  
  std::vector<tal_element_t> d;
  
  friend std::ostream & operator<<( std::ostream & out , const tal_t & t );
};



  // "\x21"
  // "\x20"

  // Time-stamped Annotation List (TAL)
  
  // Onset\21Duration\20

  // Onset and Duration must be encoded only with: '+', '-',  '.' and '0'-'9' characters, respectively
  
  // Onset must start with a '+' or a '-' character

  // Specifies the amount of seconds by which the onset of the
  // annotated event follows ('+') or precedes ('-') the
  // startdate/time of the file, that is specified in the header.

  // Duration must not contain any '+' or '-' and specifies the 
  // duration of the annotated event in seconds. 
  
  // If such a specification is not relevant, Duration can be skipped in 
  // which case its preceding 21 must also be skipped.

  // Both Onset and Duration can contain a dot ('.') but only if the 
  // fraction of a second is specified (up to arbitrary accuracy). 

  // After the time stamp, a list of annotations all sharing the same Onset and Duration may follow. 
  // Each annotation is followed by a single 20 and may not contain any 20. 
  // A 0-byte (the unprintable ASCII character with byte value 0) follows after the last 20 of this TAL. 

  // So the TAL ends with a 20 followed by a 0.

  // In each data record, the first TAL must start at the first byte of the 'EDF Annotations signal'. 
  // Subsequent TALs in the same data record must follow immediately after the trailing 0 of the preceding TAL. 
  
  // A TAL, including its trailing 0, may not overflow into another data
  // record. Each event is annotated only once, even if its duration makes
  // it extend into the time period of other data records. Unused bytes of
  // the 'EDF Annotations' signal in the remainder of the data record are
  // also filled with 0-bytes. Additional 'EDF Annotations' signals may be
  // defined according to the same specification.

  // For example, if the technician switches off the lights and closes the door 3 minutes after startdate/time, 
  // this can be stored as the 28-bytes TAL '+180{20}Lights off{20}Close door{20}{0}' without the quotes. 

  // Alternatively, the two events can be stored as two separate shorter TALs 
  // '+180{20}Lights off{20}{0}+180{20}Close door{20}{0}', also without the quotes. 
  
  // The TAL '+1800.2{21}25.5{20}Apnea{20}{0}' codes a 25.5s apnea that 
  // begins 30 minutes and 0.2s after starttime.
  

#endif
