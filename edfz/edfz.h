
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

  bool open_for_reading( const std::string & fn )
  {    
    filename = fn;

    // index must exist; this also sets the record size
    if ( ! read_index() ) return false;
    
    if ( ! bgzf_is_bgzf( filename.c_str() ) ) 
      return false;
    file = bgzf_open( filename.c_str() , "r" );
    mode = -1;
    return file != NULL;
  }

  bool open_for_writing( const std::string & fn )
  {    
    filename = fn;
    file = bgzf_open( filename.c_str() , "w" );
    mode = +1;
    return file != NULL;
  }

  void close()
  {
    if ( file == NULL ) return;

    if ( bgzf_close( file ) == -1 ) 
      Helper::halt( "problem closing " + filename );
  }
  
  inline size_t read( byte_t * p , const int n )
  {
    return bgzf_read( file , p , n );
  }

  // primary write, returns the index 
  inline int64_t write( byte_t * p , const int n )
  {
    int64_t offset = tell();
    if ( bgzf_write( file , p , n ) != n ) return -1;
    return offset;
  }

  inline void writestring( const std::string & s , int n )
  {
    std::string c = s;
    c.resize(n,' ');
    write( (byte_t*)c.data() ,  n );
  }
  
  inline void writestring( const int & s , int n )
  {
    std::string c = Helper::int2str(s);
    c.resize(n,' ');
    write( (byte_t*)c.data() ,  n );
  }

  inline void writestring( const double & s , int n )
  {
    std::string c = Helper::dbl2str_fixed(s,n);
    c.resize(n,' ');
    write( (byte_t*)c.data() , n  );
  }
  
  // primary read, given an index (for record)
  inline bool read_record( int r, byte_t * p , const int n )
  {
    std::map<int,int64_t>::const_iterator ii = index.find( r );
    if ( ii == index.end() ) return false;
    if ( ! seek( ii->second ) ) return false;
    return bgzf_read( file , p , n ) == n ;
  }

  // for header
  inline bool read_offset( int64_t offset , byte_t * p , const int n )
  {
    if ( ! seek( offset ) ) return false;
    return bgzf_read( file , p , n ) == n ;
  }
  

  inline bool is_attached()
  {
    return file != NULL;     
  }
  
  inline int64_t tell() 
  {
    return bgzf_tell( file );
  }
    
  inline bool seek(int64_t offset)
  {
    return bgzf_seek( file , offset , SEEK_SET ) == 0 ; 
  }

  inline bool eof() 
  {
    return bgzf_check_EOF( file );
  }

  void clear_index() 
  {
    index.clear();
  }

  void add_index( int r , int64_t offset )
  {
    index[ r ] =  offset ; 
  }
  
  int64_t get_index( int r ) 
  {
    std::map<int,int64_t>::const_iterator ii = index.find( r );
    if ( ii == index.end() ) return -1;
    return ii->second;
  }

  bool read_index()
  {
    std::string indexname = filename + ".idx";
    if ( ! Helper::fileExists( indexname ) ) return false;
    index.clear();
    std::ifstream I1( indexname.c_str() , std::ios::in );
    // record size first
    I1 >> record_size;
    int r = 0;
    while ( ! I1.eof() )
      {
	int64_t offset;
	I1 >> offset;
	if ( I1.eof() ) break;
	index[ r ] = offset ;
	++r;
      }    
    I1.close();
    return true;
  }

  bool write_index( const int rs )
  {
    record_size = rs;
    std::string indexname = filename + ".idx";
    std::ofstream O1( indexname.c_str() , std::ios::out );
    // first write record size
    O1 << record_size << "\n";
    std::map<int,int64_t>::const_iterator ii = index.begin();
    while ( ii != index.end() )
      {
	O1 << ii->second << "\n";
	++ii;
      }
    O1.close();
    return true;
  }
  
  BGZF * file;
  
  std::string filename;

  int mode;  // 0 closed, -1 read from , +1 write to

  // record index number -> index into .edfz
  std::map<int,int64_t> index;
  
  // as specified by EDF header
  int record_size;
  
};


#endif


// int main() 
// {

//   edfz_t edfz;

  
//   bool okay = edfz.open_for_writing( "my.edfz" );

//   std::map<int64_t,int> offsets;

//   // read files from STDIN
  
//   while ( !std::cin.eof() ) 
//     {
//       std::string line;
      
//       std::getline( std::cin , line );

//       if ( std::cin.eof() ) break;

//       // add back
//       line += "\n";

//       int64_t offset = edfz.tell();
      
//       offsets[ offset ] = line.size();

//       edfz.write( (byte_t*)line.c_str() , line.size() );
	
      
//     }
  
//   edfz.close();
    
  
//   //
//   // Now, track through backwards
//   //

//   edfz_t reader;
  
//   reader.open_for_reading( "my.edfz" );

//   std::cout << "now going to read\n";

//   std::map<int64_t,int>::const_iterator ii = offsets.begin();
//   while ( ii != offsets.end() )
//     {

//       int64_t offset = ii->first;

//       int n = ii->second;
      
//       std::cout << offset << " " << n << "\n";

//       std::vector<byte_t> buffer( n );
      
//       byte_t * p = &(buffer[0]);
      
//       reader.read( p , n );

//       std::cout << "[";
//       for (int i=0;i<n;i++) std::cout << (char)buffer[i];
//       std::cout << "]\n\n";

//       ii++;
//     }

//   reader.close();
    

  
//   return 0;

// }
