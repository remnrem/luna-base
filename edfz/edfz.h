
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

#ifndef __LUNA_EDFZ_H__
#define __LUNA_EDFZ_H__

#include <iostream>
#include "bgzf.h"
#include <cstdlib>
#include <vector>
#include <map>
#include "helper/helper.h"
#include <fstream>


struct edfz_t { 

  edfz_t() 
  {
    file = NULL;
    filename = "";
    record_size = 0;
    mode = 0;
    index.clear();
  }

  bool open_for_reading( const std::string & fn );

  bool open_for_writing( const std::string & fn );

  void close();

  size_t read( byte_t * p , const int n );

  int64_t write( byte_t * p , const int n );

  void writestring( const std::string & s , int n );

  void writestring( const int & s , int n );

  void writestring( const double & s , int n );

  bool read_record( int r, byte_t * p , const int n );

  bool read_offset( int64_t offset , byte_t * p , const int n );

  bool is_attached();
  
  int64_t tell();
  
  bool seek(int64_t offset);

  bool eof();

  void clear_index();
  
  void add_index( int r , int64_t offset , uint64_t tp , const std::string & a );
  
  int64_t get_index( int r );

  uint64_t get_tindex( int r ); 

  std::string get_annots( int r );
  
  bool read_index();

  bool write_index( const int rs );

  //
  // Members
  //
  
  BGZF * file;
  
  std::string filename;

  int mode;  // 0 closed, -1 read from , +1 write to
  
  // record index number -> index into .edfz
  std::map<int,int64_t> index;
  
  // record index number -> time-stamp
  //  (so we don't need to read whole EDF+ to get records)
  std::map<int,uint64_t> tindex;

  // also track EDF annots separately, rather than load from EDF+
  std::map<int,std::string> annots;
  
  // as specified by EDF header
  int record_size;
  
};


#endif

