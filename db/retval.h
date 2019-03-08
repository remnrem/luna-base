
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


#ifndef __LUNA_RETVAL_H__
#define __LUNA_RETVAL_H__

#include <map>
#include <set>
#include <string>
#include <vector>
#include <sstream>

#include "helper/helper.h"

struct retval_cmd_t;
struct retval_var_t;
struct retval_factor_t;
struct retval_strata_t;
struct retval_indiv_t;
struct retval_value_t;

struct strata_t;
struct timepoint_t;

typedef std::map<retval_cmd_t, 
  std::map<retval_factor_t,
  std::map<retval_var_t, 
  std::map<retval_strata_t, 
  std::map<retval_indiv_t,
  retval_value_t > > > > > retval_data_t;

// when working with Luna via an API or through R, this structure
// provides a way to return results from multiple commands for a
// single dataset (i.e. edf/individual)

// we will also track ID in the retval_value_t struct, so that databases
// containing multiple individuals can be read as a single retval_t and
// passed to R, for example

// it is designed to be a plug-in for writer, i.e. writer.var() and
// writer.value(), writer.level(), writer.unlevel(), and
// writer.epoch() all still work as expected 

// [?this may become a better, more light-weight replacement for destrat]


//   retval_t
//     retval_cmd_t  e.g. SPINDLES
//       retval_factors_t   e.g. which 'table'   "F CH"
//         retval_var_t     e.g. DENS
//           retval_indiv_t e.g. id001
//           retval_strata_t --> value   e.g. F=11, CH=C3, DENS=2.22


struct retval_cmd_t {
  retval_cmd_t( const std::string & n ) : name(n) { } 
  std::string name;
  std::set<retval_var_t> vars;  
  bool operator< ( const retval_cmd_t & rhs ) const { return name < rhs.name; }
};


struct retval_indiv_t {
  retval_indiv_t( const std::string & n ) : name(n) { } 
  std::string name;
  bool operator< ( const retval_indiv_t & rhs ) const { return name < rhs.name; }
};


struct retval_value_t {
  
  retval_value_t() { is_dbl = is_int = is_str = false; } 

  retval_value_t( double d ) :              is_dbl(true) , is_int(false), is_str(false), d(d)  { i=0; s="";  }
  retval_value_t( int64_t i ) :             is_dbl(false), is_int(true),  is_str(false), i(i)  { d=0; s=""; }
  retval_value_t( const std::string & s ) : is_dbl(false), is_int(false), is_str(true),  s(s)  { d=0; i=0; } 

  bool is_dbl, is_int, is_str;
  
  double d;
  std::string s;
  int64_t i; // to handle time-point information
  
  std::string print() const {
    if ( is_str ) return s;

    std::stringstream ss;
    if      ( is_dbl ) ss  <<d << "d";
    else if ( is_int ) ss  <<i << "i";    
    else ss << ".";
    return ss.str();    
  }
  
};

struct retval_var_t {

  retval_var_t( const std::string & n ) : name(n)
  {
    // assume int unless otherwise receive
    has_double = has_string = false;
  }
  
  std::string name;

  // track type as populated
  bool has_string;
  bool has_double;

  bool is_string() const { return has_string; }
  bool is_double() const { return has_double && ! has_string; }
  bool is_int() const { return ! ( has_double || has_string ); }

  char type() const {
    if ( is_string() ) return 'S';
    if ( is_double() ) return 'D';
    return 'I';
  }
    
  bool operator< ( const retval_var_t & rhs ) const
  {
    return name < rhs.name;
  }

};


struct retval_t { 

  // core datastore
  retval_data_t data;

  // track variable types
  std::set<std::string> var_has_strings, var_has_doubles;
  
  // to stdout
  void dump();
  
  // add a double
  void add( const retval_indiv_t & id, 
	    const retval_cmd_t & cmd ,
	    const retval_factor_t & fac ,
	    const retval_var_t & var ,
	    const retval_strata_t & stratum ,
	    const double x )
  {
    var_has_doubles.insert( var.name );
    data[cmd][fac][var][stratum][id] = retval_value_t( x );
  }
  
  void add( const retval_indiv_t & id , 
	    const retval_cmd_t & cmd ,
	    const retval_factor_t & fac ,
	    const retval_var_t & var ,
	    const retval_strata_t & stratum ,
	    int  x )
  {
    data[cmd][fac][var][stratum][id] = retval_value_t( (int64_t)x );
  }
  
  void add( const retval_indiv_t & id, 
	    const retval_cmd_t & cmd ,
	    const retval_factor_t & fac ,
	    const retval_var_t & var ,
	    const retval_strata_t & stratum ,
	    int64_t  x )
  {
    data[cmd][fac][var][stratum][id] = retval_value_t( x );
  }

  void add( const retval_indiv_t & id, 
	    const retval_cmd_t & cmd ,
	    const retval_factor_t & fac ,
	    const retval_var_t & var ,
	    const retval_strata_t & stratum ,
	    uint64_t  x )
  {
    data[cmd][fac][var][stratum][id] = retval_value_t( (int64_t)x );
  }

  void add( const retval_indiv_t & id, 
	    const retval_cmd_t & cmd ,
	    const retval_factor_t & fac ,
	    const retval_var_t & var ,
	    const retval_strata_t & stratum ,
	    const std::string & x )
  {  
    var_has_strings.insert( var.name );
    data[cmd][fac][var][stratum][id] = retval_value_t( x );
  }
    
};



struct retval_factor_t {
  
  retval_factor_t( const strata_t & s , const timepoint_t & tp );
  
  std::set<std::string> factors; // i.e. just factor names, not levels.
  // i.e. these specify which 'virtual table' we are looking at
  //  e.g.   F
  //         F CH
  //         F CH E

  void add( const std::string & f ) { factors.insert( f ); } 

  std::string print() const { return Helper::stringize( factors ); } 
  
  bool operator< ( const retval_factor_t & rhs ) const
  {
    if ( factors.size() < rhs.factors.size() ) return true;
    if ( factors.size() > rhs.factors.size() ) return false;
    
    std::set<std::string>::const_iterator ii = factors.begin();
    std::set<std::string>::const_iterator jj = rhs.factors.begin();
    while ( ii != factors.end() ) 
      {
	if ( *ii < *jj ) return true;
	if ( *jj < *ii ) return false;
	++ii;
	++jj;
      }

    return false;
  }

};

struct retval_factor_level_t {

  // null
  retval_factor_level_t() { is_str = is_int = is_dbl = false; }

  retval_factor_level_t( const std::string & f , const std::string & s )
  { factor = f ; is_str = true ; is_int = false ; is_dbl = false; str_level = s; }

  retval_factor_level_t( const std::string & f , int i )
  { factor = f ; is_int = true ; is_str = false ; is_dbl = false; int_level = i; }

  retval_factor_level_t( const std::string & f , double d )
  { factor = f ; is_dbl = true ; is_str = false ; is_int = false; dbl_level = d; }
		        
  std::string factor;

  // ensure correct sorting
  bool        is_str, is_int, is_dbl;

  std::string str_level;
  int         int_level;
  double      dbl_level;
  
  std::string print() const
  {
    std::stringstream ss;
    if      ( is_str ) ss << factor << "=" << str_level;
    else if ( is_int ) ss << factor << "=" << int_level;
    else if ( is_dbl ) ss << factor << "=" << dbl_level;
    else ss << ".";
    return ss.str();
  }
  
  bool operator< ( const retval_factor_level_t & rhs ) const {

    if ( factor < rhs.factor ) return true;
    if ( factor > rhs.factor ) return false;
    
    if ( is_str ) return str_level < rhs.str_level;
    if ( is_int ) return int_level < rhs.int_level;
    if ( is_dbl ) return dbl_level < rhs.dbl_level;

    // empty / should not occur
    return false;
  }
    
};



struct retval_strata_t { 

  // set of text key=value pairs, e.g. 
  // e.g. SS=N2, F=11,
  
  retval_strata_t( strata_t & strata , timepoint_t & tp );
  
  std::set<retval_factor_level_t> factors;

  retval_factor_level_t find( const std::string & f ) const
  {
    std::set<retval_factor_level_t>::const_iterator ff = factors.begin();
    while ( ff != factors.end() )
      {
	if ( ff->factor == f ) return *ff;
	++ff;
      }
    return retval_factor_level_t();
  }

  
  void add( const retval_factor_level_t & fl ) { factors.insert( fl ) ; } 

  std::string print() const
  {
    std::stringstream ss;
    std::set<retval_factor_level_t>::const_iterator ff = factors.begin();
    while ( ff != factors.end() )
      {
	if ( ff != factors.begin() ) ss << ";";
	ss << ff->print();
	++ff;
      }
    return ss.str();
  }
  
  
  bool operator< (const retval_strata_t & rhs ) const { 
    if ( factors.size() < rhs.factors.size() ) return true;
    if ( factors.size() > rhs.factors.size() ) return false;
    std::set<retval_factor_level_t>::const_iterator ii = factors.begin();
    std::set<retval_factor_level_t>::const_iterator jj = rhs.factors.begin();
    while ( ii != factors.end() ) 
      {
	if ( *ii < *jj ) return true;
	if ( *jj < *ii ) return false;
	++ii;
	++jj;
      }

    return false;
  }

};



#endif
