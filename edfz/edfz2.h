
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

#ifndef __LUNA_EDFZ2_H__
#define __LUNA_EDFZ2_H__

#include <iostream>
#include <cstdlib>
#include <vector>
#include <map>
#include "helper/helper.h"
#include "helper/zfstream.h"
#include <fstream>

typedef unsigned char byte_t;

// replacement to BGZF... use basic zfstream wrapper around libz and
// abandon random access.. i.e. will preread all in a single go, but
// also will use larger record sizes for better compression...

struct edfz2_t { 

  edfz2_t() 
  {
    filename = "";
  }
  
  bool open_for_reading( const std::string & fn );
  
  bool open_for_writing( const std::string & fn );
  
  void close();
  
  size_t read( byte_t * p , const int n );
  
  void write( byte_t * p , const int n );
  
  void writestring( const std::string & s , int n );

  void writestring( const int & s , int n );
  
  void writestring( const double & s , int n );
  
    
  //
  // Members
  //
  
  gzofstream zout;

  gzifstream zin;
      
  std::string filename;
    
};


#endif

