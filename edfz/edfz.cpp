
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

#include "edfz/edfz.h"


bool edfz_t::open_for_reading( const std::string & fn )
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

bool edfz_t::open_for_writing( const std::string & fn )
{    
  filename = fn;
  file = bgzf_open( filename.c_str() , "w" );
  mode = +1;
  return file != NULL;
}

void edfz_t::close()
{
  if ( file == NULL ) return;
  
  if ( bgzf_close( file ) == -1 ) 
    Helper::halt( "problem closing " + filename );
}

size_t edfz_t::read( byte_t * p , const int n )
{
  return bgzf_read( file , p , n );
}

// primary write, returns the index 
int64_t edfz_t::write( byte_t * p , const int n )
{
  int64_t offset = tell();
  if ( bgzf_write( file , p , n ) != n ) return -1;
  return offset;
}

void edfz_t::writestring( const std::string & s , int n )
{
  std::string c = s;
  c.resize(n,' ');
  write( (byte_t*)c.data() ,  n );
}

void edfz_t::writestring( const int & s , int n )
{
  std::string c = Helper::int2str(s);
  c.resize(n,' ');
  write( (byte_t*)c.data() ,  n );
}

void edfz_t::writestring( const double & s , int n )
{
  std::string c = Helper::dbl2str_fixed(s,n);
  c.resize(n,' ');
  write( (byte_t*)c.data() , n  );
}

// primary read, given an index (for record)
bool edfz_t::read_record( int r, byte_t * p , const int n )
{
  std::cout << "\n\n --> entering read_record \n";
  std::map<int,int64_t>::const_iterator ii = index.find( r );
  if ( ii == index.end() ) return false;
  
  if ( ! seek( ii->second ) ) return false;
  
  // tmp
  int64_t offset = tell();
  std::cout << " offset = " << offset << " for rec " << r << " " << ii->second << "\n";
  // --tmp
  
  ssize_t ss = bgzf_read( file , p , n );
  std::cout << "done read\n";
  return ss == n ; 
  //  return bgzf_read( file , p , n ) == n ;
}

// for header
bool edfz_t::read_offset( int64_t offset , byte_t * p , const int n )
{
  if ( ! seek( offset ) ) return false;
  return bgzf_read( file , p , n ) == n ;
}


bool edfz_t::is_attached()
{
  return file != NULL;     
}

int64_t edfz_t::tell() 
{
  return bgzf_tell( file );
}

bool edfz_t::seek(int64_t offset)
{
  return bgzf_seek( file , offset , SEEK_SET ) == 0 ; 
}

bool edfz_t::eof() 
{
  return bgzf_check_EOF( file );
}

void edfz_t::clear_index() 
{
  index.clear();
  tindex.clear();
  annots.clear();
}

void edfz_t::add_index( int r , int64_t offset , uint64_t tp , const std::string & a )
{
  index[ r ] = offset ;
  tindex[ r ] = tp ;
  annots[ r ] = a ;
}

int64_t edfz_t::get_index( int r ) 
{
  std::map<int,int64_t>::const_iterator ii = index.find( r );
  if ( ii == index.end() ) return -1;
  return ii->second;
}

uint64_t edfz_t::get_tindex( int r ) 
{
  std::map<int,uint64_t>::const_iterator ii = tindex.find( r );
  if ( ii == tindex.end() ) return 0;
  return ii->second;
}

std::string edfz_t::get_annots( int r )
{
  std::map<int,std::string>::const_iterator ii = annots.find( r );
  if ( ii == annots.end() ) return ".";
  return ii->second;
}


bool edfz_t::read_index()
{
  std::string indexname = filename + ".idx";
  if ( ! Helper::fileExists( indexname ) ) return false;
  index.clear();
  std::ifstream I1( indexname.c_str() , std::ios::in );
  
  // index version
  std::string line;
  Helper::safe_getline( I1 , line );
  if ( line != "EDFZv1" )
    Helper::halt( "expecting EDFZv1 format index: please remake the index" );
  
  // record size
  Helper::safe_getline( I1 , line );
  if ( ! Helper::str2int( line , &record_size ) )
    Helper::halt( "expecting EDFZv1 format index: second entry = # records" );

  // indices
  int r = 0;
  while ( ! I1.eof() )
    {
      int64_t offset;
      uint64_t tp;
      
      std::string line;
      Helper::safe_getline( I1 , line );
      
      if ( I1.eof() ) break;
			  
      std::vector<std::string> tok = Helper::parse( line , "\t" );
      
      if ( tok.size() != 3 )
	Helper::halt( "invalid .idx line:\n" + line );
      
      if ( ! Helper::str2signed_int64( tok[0] , &offset ) )
	Helper::halt( "bad .idx:\n" + line );
      
      if ( ! Helper::str2int64( tok[1] , &tp ) )
	Helper::halt( "bad .idx:\n" + line );

      index[ r ] = offset ;
      tindex[ r ] = tp;
      annots[ r ] = tok[2];
      ++r;
    }    
  I1.close();
  return true;
}

bool edfz_t::write_index( const int rs )
{
  record_size = rs;
  std::string indexname = filename + ".idx";
  std::ofstream O1( indexname.c_str() , std::ios::out );
  
  // index version
  O1 << "EDFZv1\n";    
  
  // first write record size
  O1 << record_size << "\n";
  
  // offsets, and time-stamps
  std::map<int,int64_t>::const_iterator ii = index.begin();
  while ( ii != index.end() )
    {
      O1 << ii->second << "\t"
	 << tindex[ ii->first ] << "\t"
	 << annots[ ii->first ] << "\n";
      ++ii;
    }
  O1.close();
  return true;
}
