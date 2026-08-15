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

#include "stats/dpp-io.h"

#include "helper/helper.h"
#include "helper/luna_io.h"
#include "helper/logger.h"

#include <cstdint>
#include <limits>

extern logger_t logger;

void dpp_io::bwrite( std::ofstream & O , const std::string & s )
{
  // uint16_t (not pops/io.cpp's uint8_t): a length-prefixed ID longer than
  // 255 bytes would otherwise be silently truncated, corrupting the record
  // layout for every subsequent field. This format is new/unshipped, so
  // there is no on-disk compatibility to preserve by matching POPS exactly.
  uint16_t l = (uint16_t) s.size();
  O.write( (char*)( &l ) , sizeof(uint16_t) );
  O.write( s.c_str() , l );
}

void dpp_io::bwrite( std::ofstream & O , int i )
{
  O.write( (char*)( &i ) , sizeof(int) );
}

void dpp_io::bwrite( std::ofstream & O , double d )
{
  O.write( (char*)( &d ) , sizeof(double) );
}

std::string dpp_io::bread_str( std::ifstream & I )
{
  uint16_t len;
  I.read( (char*)( &len ) , sizeof(uint16_t) );
  std::vector<char> b( len );
  if ( len > 0 ) I.read( &b[0] , len );
  return std::string( b.begin() , b.end() );
}

int dpp_io::bread_int( std::ifstream & I )
{
  int i;
  I.read( (char*)( &i ) , sizeof(int) );
  return i;
}

double dpp_io::bread_dbl( std::ifstream & I )
{
  double d;
  I.read( (char*)( &d ) , sizeof(double) );
  return d;
}

void dpp_io::save( const std::string & f , const dpp_matrix_t & m , int n_features , bool append )
{
  const int nr = (int) m.time_sec.size();

  logger << "  writing binary DPP feature matrix (" << nr << " rows, "
	 << n_features << " features) to " << f << "\n";

  std::ios_base::openmode mode = std::ios::binary | std::ios::out;
  mode |= append ? std::ios::app : std::ios::trunc;

  std::ofstream OUT1 = LunaIO::open_ofstream( Helper::expand( f ) , mode );

  bwrite( OUT1 , m.id );
  bwrite( OUT1 , nr );
  bwrite( OUT1 , n_features );

  for (int i=0; i<nr; i++)
    {
      bwrite( OUT1 , m.time_sec[i] );
      for (int j=0; j<n_features; j++)
	bwrite( OUT1 , j < (int)m.X[i].size() ? m.X[i][j] : std::numeric_limits<double>::quiet_NaN() );
    }

  OUT1.close();
}

std::vector<dpp_matrix_t> dpp_io::load( const std::string & f , int expected_n_features )
{
  std::vector<dpp_matrix_t> out;

  std::ifstream IN1 = LunaIO::open_ifstream( Helper::expand( f ) , std::ios::binary | std::ios::in );

  while ( 1 )
    {
      std::string id = bread_str( IN1 );
      if ( IN1.eof() || IN1.bad() ) break;

      dpp_matrix_t m;
      m.id = id;

      const int nr = bread_int( IN1 );
      const int nf = bread_int( IN1 );

      if ( expected_n_features >= 0 && nf != expected_n_features )
	Helper::halt( "data in " + f + " does not match the expected feature schema" );

      m.time_sec.resize( nr );
      m.X.resize( nr );
      for (int i=0; i<nr; i++)
	{
	  m.time_sec[i] = bread_dbl( IN1 );
	  m.X[i].resize( nf );
	  for (int j=0; j<nf; j++)
	    m.X[i][j] = bread_dbl( IN1 );
	}

      out.push_back( m );
    }

  IN1.close();
  return out;
}

std::vector<dpp_matrix_t> dpp_io::load_files( const std::vector<std::string> & files , int expected_n_features )
{
  std::vector<dpp_matrix_t> out;
  for (int i=0; i<(int)files.size(); i++)
    {
      std::vector<dpp_matrix_t> part = load( files[i] , expected_n_features );
      out.insert( out.end() , part.begin() , part.end() );
    }
  return out;
}
