
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

#ifndef __LUNA_PARAM_H__
#define __LUNA_PARAM_H__

#include <string>
#include <map>
#include <set>
#include <vector>
#include <sstream>

struct edf_t;

//
// Helper to parse command syntax
//

struct param_t
{

 public:
  
  void add( const std::string & option , const std::string & value = "" ); 

  void add_hidden( const std::string & option , const std::string & value = "" );

  int size() const;
  
  void parse( const std::string & s );
  
  void update( const std::string & id , const std::string & wc );

  void clear();
  
  bool has(const std::string & s ) const;
  
  bool empty(const std::string & s ) const;
  
  // if ! has(X) return F
  // else return yesno(value(X))
  // default1 = returned if 's' not present
  // default2 = returned if 's' if present, but not value specified (i.e. not s=F or s=T
  bool yesno(const std::string & s , const bool default1 = false , const bool default2 = true ) const;

  std::string value( const std::string & s , const bool uppercase = false ) const;
 
  bool single() const;  

  std::string single_value() const ;

  std::string single_pair(std::string * ) const ;
  
  std::string requires( const std::string & s , const bool uppercase = false ) const;
  
  int requires_int( const std::string & s ) const;
  
  double requires_dbl( const std::string & s ) const;

  std::string dump( const std::string & indent = "  ", const std::string & delim = "\n" ) const;

  std::set<std::string> strset( const std::string & k , const std::string delim = "," , const bool uppercase = false ) const;

  std::set<std::string> strset_xsigs( const std::string & k , const std::string delim = "," , const bool uppercase = false ) const;
  
  std::vector<std::string> strvector( const std::string & k , const std::string delim = "," , const bool uppercase = false ) const;
  
  std::vector<std::string> strvector_xsigs( const std::string & k , const std::string delim = "," , const bool uppercase = false ) const;
  
  std::vector<double> dblvector( const std::string & k , const std::string delim = "," ) const;

  std::vector<int> intvector( const std::string & k , const std::string delim = "," ) const;

  std::set<std::string> keys() const;

private:

  std::map<std::string,std::string> opt;

  std::set<std::string> hidden;


};


#endif
