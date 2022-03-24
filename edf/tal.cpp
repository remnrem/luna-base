
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


#include "tal.h"
#include "edf.h"

#include "defs/defs.h"
#include "intervals/intervals.h"
#include <iostream>

tal_t edf_t::tal( const int signal , const int rec )
{
  
  // extract data from EDF record into TAL format
  tal_t t( this , signal , rec );
  
  const int np = 2 * header.n_samples[signal];
  
  std::string s( np , '\x00' ); // maximum annot. size

  //
  // need to load the record?
  //
  
  if ( ! loaded( rec ) ) 
    {
      edf_record_t record( this ); 
      record.read( rec );
      records.insert( std::map<int,edf_record_t>::value_type ( rec , record ) );
    }

  //
  // Pull data
  //
  
  // std::cout << "records " << ( records.find( rec ) != records.end() ) << "\n";
  // std::cout << "s = " << records[rec].data.size() << "\n";
  // std::cout << "signal = " <<signal << "\n";

  const std::vector<int16_t> & raw = records.find(rec)->second.data[signal];

  const int np_used = raw.size();
  
  if ( np_used > np ) 
    Helper::halt( "problem in getting TAL" );
  
  for (int i=0;i<np_used;i++)
    {      
      s[i] = (char)raw[i];      
    }

  //
  // convert raw string to TAL
  //

  t.decode(s);

  return t;  
  
}



std::ostream & operator<<( std::ostream & out , const tal_element_t & t )
{
  out << "<" << t.onset << "|";
  if ( t.duration != 0 ) out << t.duration;

  if ( t.name == globals::edf_timetrack_label ) 
    out << "|(time-stamp, secs)";
  else 
    if ( t.name != "" ) out << "|" << t.name;  
  out << ">";

  return out;
}

std::ostream & operator<<( std::ostream & out , const tal_t & t )
{
  for (int i=0;i<t.d.size();i++) 
    out << t.d[i] << "\n";
  return out;
}


tal_t::tal_t( edf_t * edf , int signal , int r )
{

  record = r;
  
  if ( record < 0 || record >= edf->header.nr_all )
    Helper::halt( "bad record # requested" );
  
  if ( ! edf->timeline.retained( record ) )
    Helper::halt( "bad record # requested" );

  if ( signal < 0 || signal >= edf->header.ns ) 
    Helper::halt( "bad signal # requested" );
  
  if ( ! edf->header.is_annotation_channel( signal ) ) 
    Helper::halt( "not an annotation channel" );

}

void tal_t::decode( const std::string & str )
{
  
  //  std::cout << "DECODE [" << str << "]\n";

  // Onset { [\x15] Duration } [\x14] {text} [\x14 ] 
  // Onset { [\x15] Duration } [\x14] {text} [\x14 ] 
  // [\x00 ]
  
  // Onset [\x14] {text} [\x14]  
  // [\x00]

	    
  // ASCII 20 --> \x14
  // ASCII 21 --> \x15
  
  d.clear();

  int p = 0;
  int cur = 0;
  
  // 1) Split by '\x00' for each TA element at a high level
  //    This should not keep empty tokens, so end of file 
  //    will be ignored effectively.
  
  const bool NO_EMPTIES = false;
 
  std::vector<std::string> toks = Helper::char_split( str , '\x00' , NO_EMPTIES );

  bool added_time_stamp = false;
  
  for (int t=0;t<toks.size();t++)
    {           
      
      std::vector<std::string> subs = Helper::char_split( toks[t] , '\x14' , NO_EMPTIES );
      
      // std::cout << "subs size = " << subs.size() << "\n";

      // this should contain at least one field (time)
      if ( subs.size() < 1 ) continue;

      // time can optionally have a duration
      std::vector<std::string> ts = Helper::char_split( subs[0] , '\x21' , NO_EMPTIES );
      
      //std::cout << "ts size = " << ts.size() << "\n";

      double onset = 0;
      double duration = 0;
      
      if ( ! ( ts.size() == 1 || ts.size() == 2 ) )
	continue;

      if ( ! Helper::str2dbl( ts[0] , &onset ) )
	Helper::halt( "problem converting time-stamp, " + ts[0] );
      
      if ( ts.size() == 2 )
	{
	  if ( ! Helper::str2dbl( ts[0] , &onset ) )
	    Helper::halt( "problem converting time-stamp, " + ts[0] );	  
	}

      // time-stamp?

      if ( ! added_time_stamp )
	{
	  tal_element_t tae( onset , duration , globals::edf_timetrack_label );
	  d.push_back( tae );
	  added_time_stamp = true;	  
	  //std::cout << "adding TS : " << tae << "\n";
	}

      // textual label(s)?
      //  (we might be here w/ 'skip EDF annots' set to true; which
      //  means we skip this -- but we'll be here if it is a EDF+D, meaning
      //  that we need the time-track from the EDF Annots

      if ( ! globals::skip_edf_annots )
	for (int j=1;j<subs.size();j++)
	  {
	    tal_element_t tae( onset , duration , subs[j] );
	    d.push_back( tae );
	    //std::cout << "adding TAE : " << tae << "\n";
	  }
      
    }
  
}

std::string tal_t::encode() const
{
  return "";
}



