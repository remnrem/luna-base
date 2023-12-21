
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

#ifndef __LUNAPI_H__
#define __LUNAPI_H__

#include "luna.h"

struct lunapi_t {
  

  //
  // global variables 
  //

  static void init(); 
  
  static void var( const std::string & key , const std::string & value );

  static std::string var( const std::string & key );

  //
  // new instance
  //

  lunapi_t( const std::string & id ) : id(id) 
  {
    reset();  
  }

  //
  // attach data 
  //

  bool attach_edf( const std::string & filename );

  bool attach_annot( const std::string & filename );
  
  //
  // individual variables
  //
  
  void ivar( const std::string & key , const std::string & value );
  
  std::string ivar( const std::string & key ) const;


  //
  // commands
  //
  
  retval_t eval( const std::string & );


  //
  // helpers
  //

  void reset();
   
  
private:

  // 0 empty; +1 attached okay, -1 problem 
  int state;
  
  std::string id;
  
  std::string filename;
  
  edf_t edf;
  
};

#endif
