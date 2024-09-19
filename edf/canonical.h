
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

#ifndef __CANONICAL_H__
#define __CANONICAL_H__

#include <vector>
#include <string>
#include <map>
#include <set>
#include <sstream>

#include "helper/helper.h"

struct signal_list_t;
struct edf_t;
struct edf_header_t;
struct param_t;

struct cansigs_t {
  std::set<std::string> used;
  std::map<std::string,bool> okay;
  std::map<std::string,std::string> sig;
  std::map<std::string,std::string> ref;
};


struct canon_rule_t {
  
  canon_rule_t( const std::vector<std::string> & lines );
  canon_rule_t( const std::string & l );
  
  std::string canonical_label;

  // unless
  std::set<std::string> unless;

  // group
  std::set<std::string> group;

  // requires
  std::vector<std::string> req_sig;
  std::vector<std::string> req_ref;
  std::map<std::string,std::string> req_transducer;
  std::map<std::string,std::string> req_unit;
  int req_scale;
  int req_sr_min;
  int req_sr_max;
  
  // sets
  int set_sr;
  std::string set_unit;

  // "<<-" relabelling rule
  bool relabel_canonical;
  std::vector<std::string> original_canonical_label;

  // special: closed
  bool closed;
  
};

struct canon_edf_signal_t {
  
  canon_edf_signal_t( edf_header_t & hdr , const int slot );
  
  canon_edf_signal_t( const std::string & label ,
		      const int sr ,
		      const std::string & unit ,
		      const std::string & trandsucer ,
		      const int scale);
  
  canon_edf_signal_t( const std::string & label )
    : label(label)
  { }
  
  std::string label;

  int sr;
  
  std::string unit;

  std::string transducer;

  int scale;
  
  bool operator<( const canon_edf_signal_t & rhs ) const
  {
    return label < rhs.label;
  }
  
};
  
struct canonical_t 
{
  
  canonical_t( edf_t & edf , param_t & param );

  std::vector<std::string> preprocess( const std::string & filename );
  
  int read( const std::string & filename );

  bool apply_this( const std::string & );
  
  void proc( );

  // general (static) rules

  static std::vector<canon_rule_t> rules;

  static std::map<int,std::string> scale_codes;
  
  // for this EDF

  edf_t & edf;

  std::set<canon_edf_signal_t> signals;
  
  // options
  bool drop_originals;
  bool dry_run;
  bool only_check_labels; // no output (for use w/ the mapper utility)
  bool verbose;
  bool retain_prefiltering;
  bool dump_definition_file; // if using templates, show the actual thing (debug mode)
   
  // output class used in mapper-cgi-util mode
  cansigs_t retval;

  // only inc these canonica rules
  std::set<std::string> canins;
  
  // ignore these canonical rules
  std::set<std::string> canouts;
  
  // can be assigned to multiple groups
  std::set<std::string> group;

  // helpers
  bool is_in( const std::set<std::string> & a , const std::set<std::string> & b )
  {
    std::set<std::string>::const_iterator aa = a.begin();
    while ( aa != a.end() )
      {
	if ( b.find( *aa ) != b.end() ) return true;
	++aa;
      }
    return false;
  }
  
  bool is_in( const std::string & a , const std::set<std::string> & b )
  {
    return b.find( a ) != b.end() ; 
  }

  bool match( const std::vector<std::string> & a , std::set<canon_edf_signal_t> & b , std::string * match )
  {
    // find first instance in 'a' (req sig list) that is in b (EDF) 
    for (int i=0; i<a.size(); i++)
      {
	canon_edf_signal_t s1( a[i] );
	if ( b.find( s1 ) != b.end() )
	  {
	    *match = a[i];
	    return true;
	  }
      }
    return false;    
  }

  // special match for a reference -- allows for linked references
  bool ref_match( const std::vector<std::string> & a , std::set<canon_edf_signal_t> & b , std::string * match );

  static bool empty_field( const std::string & s )
  {
    if ( s == "" || s == "." ) return true;
    const std::string s2 = Helper::trim( s );
    if ( s2 == "" || s2 == "." ) return true;
    const std::string s3 = Helper::trim( s , '_' , '_' );
    if ( s3 == "" || s3 == "." ) return true;    
    return false;
  }
  
  std::string print( const std::set<std::string> & s )
  {
    std::stringstream ss;
    std::set<std::string>::const_iterator ii = s.begin();
    while ( ii != s.end() )
      {
	ss << ( ii != s.begin() ? "," : "" ) << *ii;
	++ii;
      }
    return ss.str();
  }

  std::string print( const std::vector<std::string> & s )
  {
    std::stringstream ss;
    std::vector<std::string>::const_iterator ii = s.begin();
    while ( ii != s.end() )
      {
	ss << ( ii != s.begin() ? "," : "" ) << *ii;
	++ii;
      }
    return ss.str();
  }

  std::string print( const std::map<std::string,std::string> & s )
  {
    std::stringstream ss;
    std::map<std::string,std::string>::const_iterator ii = s.begin();
    while ( ii != s.end() )
      {
	ss << ( ii != s.begin() ? "," : "" ) << ii->first << "=" << ii->second;
	++ii;
      }
    return ss.str();
  }
  
};



#endif
