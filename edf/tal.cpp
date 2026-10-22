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
#include <iomanip>
#include <iostream>
#include <sstream>

namespace {
  // render a raw TAL fragment (which contains non-printable
  // delimiters '\x00', '\x14', '\x15', and possibly garbage bytes
  // from a corrupted/mis-parsed EDF+ record) in a form that is
  // actually legible in an error message
  std::string tal_debug_string( const std::string & s )
  {
    std::stringstream ss;
    for (unsigned int i=0; i<s.size(); i++)
      {
	const unsigned char c = (unsigned char)s[i];
	if ( c == 0x00 ) ss << "<NUL>";
	else if ( c == 0x14 ) ss << "<TAL>";
	else if ( c == 0x15 ) ss << "<DUR>";
	else if ( c < 0x20 || c > 0x7e )
	  {
	    ss << "<0x" << std::hex << std::setw(2) << std::setfill('0')
	       << (int)c << std::dec << std::setfill(' ') << ">";
	  }
	else ss << (char)c;
      }
    return ss.str();
  }
}

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

std::string tal_t::export_annots() const
{
  std::stringstream ss;
  // EDFZ stores this serialized representation in its embedded index.
  // Use round-trip-safe precision so that annotation times are not
  // truncated by the stream's default six-significant-digit precision.
  ss << std::setprecision( 17 );
  bool printed = false;
  //  std::cout << "ds = " << d.size() << "\n";
  for (int i=0; i<d.size(); i++)
    {
      const tal_element_t & el = d[i];
      
      // no need to print time-track
      if ( el.name == globals::edf_timetrack_label ) continue;
      
      if ( printed ) ss << ",";
      
      ss << "\"" << el.onset << "|"
	 << el.duration << "|"
	 << ( el.name == "" ? "." : el.name )
	 << "\"";
      
      printed = true;
    }

  //  std::cout << " string = [" << ss.str() << "]\n";
  return ss.str();
}



tal_t::tal_t( edf_t * edf , int signal , int r )
{

  record = r;

  pedf = edf;

  channel = signal;

  if ( record < 0 || record >= edf->header.nr_all )
    Helper::halt( "bad record # requested" );

  if ( ! edf->timeline.retained( record ) )
    Helper::halt( "bad record # requested" );

  if ( signal < 0 || signal >= edf->header.ns )
    Helper::halt( "bad signal # requested" );

  if ( ! edf->header.is_annotation_channel( signal ) )
    Helper::halt( "not an annotation channel" );

}

std::string tal_t::tal_debug_context() const
{
  std::stringstream ss;

  ss << "record " << record
     << ", channel #" << channel
     << ( pedf != NULL && channel >= 0 && channel < pedf->header.ns ?
	  " (" + pedf->header.label[channel] + ")" : "" );

  return ss.str();
}

void tal_t::decode( const std::string & str )
{
  
  //std::cout << "DECODE [" << str << "]\n";

  // Onset { [\x15] Duration } [\x14] {text} [\x14 ] 
  // Onset { [\x15] Duration } [\x14] {text} [\x14 ] 
  // [\x00 ]
  
  // Onset [\x14] {text} [\x14]  
  // [\x00]

  //  std::cout << "\n";
  
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
      
      //std::cout << "subs size = " << subs.size() << "\n";

      // this should contain at least one field (time)
      if ( subs.size() < 1 ) continue;

      // time can optionally have a duration
      std::vector<std::string> ts = Helper::char_split( subs[0] , '\x15' , NO_EMPTIES );
      
      //std::cout << "ts size = " << ts.size() << "\n";

      double onset = 0;
      double duration = 0;
      
      if ( ! ( ts.size() == 1 || ts.size() == 2 ) )
	{
	  // this element is always skipped (not gated by skip-bad-tals): it
	  // has no usable onset candidate at all, so there's nothing to
	  // halt on in the first place
	  logger << "  *** warning: skipping malformed TAL element (expecting 1 or 2 "
		 << "'\\x15'-delimited onset/duration fields, found " << ts.size() << ") in "
		 << tal_debug_context() << "\n"
		 << "              raw TAL   : [" << tal_debug_string( toks[t] ) << "]\n"
		 << "              full record: [" << tal_debug_string( str ) << "]\n";
	  continue;
	}

      if ( ! Helper::str2dbl( ts[0] , &onset ) )
	{
	  const std::string msg = "problem converting time-stamp (expecting a +/- numeric onset) in "
	    + tal_debug_context() + "\n"
	    + "       offending field : [" + tal_debug_string( ts[0] ) + "]\n"
	    + "       raw TAL element : [" + tal_debug_string( toks[t] ) + "]\n"
	    + "       full raw record : [" + tal_debug_string( str ) + "]";

	  if ( ! globals::skip_bad_tals )
	    Helper::halt( msg + "\n"
			  + "       hint: this TAL element could instead be skipped (rather than\n"
			  + "             halting) by setting the option skip-bad-tals=T" );

	  logger << "  *** warning: skipping malformed TAL element (skip-bad-tals=T):\n" << msg << "\n";
	  continue;
	}

      if ( ts.size() == 2 )
	{
	  if ( ! Helper::str2dbl( ts[1] , &duration ) )
	    {
	      const std::string msg = "problem converting duration in "
		+ tal_debug_context() + "\n"
		+ "       offending field : [" + tal_debug_string( ts[1] ) + "]\n"
		+ "       raw TAL element : [" + tal_debug_string( toks[t] ) + "]\n"
		+ "       full raw record : [" + tal_debug_string( str ) + "]";

	      if ( ! globals::skip_bad_tals )
		Helper::halt( msg + "\n"
			      + "       hint: this TAL element could instead be skipped (rather than\n"
			      + "             halting) by setting the option skip-bad-tals=T" );

	      logger << "  *** warning: skipping malformed TAL element (skip-bad-tals=T):\n" << msg << "\n";
	      continue;
	    }
	}
      
      // time-stamp?
      //      std::cout << std::setprecision(12);
      if ( ! added_time_stamp )
	{
	  tal_element_t tae( onset , duration , globals::edf_timetrack_label );
	  d.push_back( tae );
	  added_time_stamp = true;	  
	  //  std::cout << "adding TS : " << tae << "\n";
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
