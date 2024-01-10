
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


#ifndef __LUNA_RTABLES_H__
#define __LUNA_RTABLES_H__

#include <map>
#include <set>
#include <string>
#include <vector>
#include <sstream>
#include <variant>

#include "helper/helper.h"
#include "db/retval.h"

// Luna 'return-tables'
//  i.e. convery interval retval_t (which is connected to the writer)
//       to a set of tables, that the lunAPI can return

struct retval_t;

// internal elements
typedef std::variant<std::string,double,int,std::monostate> rtable_elem_t;
typedef std::vector<std::vector<rtable_elem_t> > rtable_data_t;

// return elements:
// data *plus column variables* in friendly format
typedef std::tuple<std::vector<std::string>,rtable_data_t> rtable_return_t;

// multiple return-tables
typedef std::map<std::string,std::map<std::string,rtable_return_t> > rtables_return_t;


struct rtable_t {
  
  rtable_t();

  std::vector<std::string> cols;
  
  // ( col name , data matrix )-tuple
  rtable_data_t data;

  // must be of similar row count... (at input)
  int nrows;

  std::string dump() const;
  
  void checkrows( int n );
  
  void add( const std::string & v , const std::vector<std::string> & x );
  
  void add( const std::string & v , const std::vector<std::string> & x , const std::vector<bool> & m );

  // doubles
  
  void add( const std::string & v , const std::vector<double> & x );
  
  void add( const std::string & v , const std::vector<double> & x , const std::vector<bool> & m );

  // doubles
  
  void add( const std::string & v , const std::vector<int> & x );
  
  void add( const std::string & v , const std::vector<int> & x , const std::vector<bool> & m );
    
};


struct rtables_t {

  rtables_t() { };
  
  rtables_t( const retval_t & retval );

  void clear();
  
  std::vector<std::string> commands() const;
  
  std::vector<std::pair<std::string,std::string> > list() const;
    
  rtable_t table( const std::string & cmd , const std::string & strata ) const;

  rtable_return_t data( const std::string & cmd , const std::string & strata ) const; 
  
  rtables_return_t data() const;

  void dump() const;
  
  std::map<std::string,std::map<std::string,rtable_t> > tables;
  
};



#endif
