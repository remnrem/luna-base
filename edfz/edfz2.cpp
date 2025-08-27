
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

#include "edfz/edfz2.h"

bool edfz2_t::open_for_reading( const std::string & fn )
{    
  filename = fn;
  zin.open( fn.c_str() , std::ios::in | std::ios::binary );
  return zin.good();
}

bool edfz2_t::open_for_writing( const std::string & fn )
{    
  filename = fn;
  zout.open( fn.c_str() , std::ios::out | std::ios::binary );  
  return zout.good();
}

void edfz2_t::close()
{
  if ( zin.is_open() ) zin.close();
  if ( zout.is_open() ) zout.close();
}

size_t edfz2_t::read( byte_t * p , const int n )
{
  zin.read(reinterpret_cast<char*>(p), n);
  return zin.gcount();
}

// primary write, returns the index 
void edfz2_t::write( byte_t * p , const int n )
{
  zout.write( reinterpret_cast<char*>(p), n );
}

void edfz2_t::writestring( const std::string & s , int n )
{
  std::string c = s;
  c.resize(n,' ');
  write( (byte_t*)c.data() ,  n );
}

void edfz2_t::writestring( const int & s , int n )
{
  std::string c = Helper::int2str(s);
  c.resize(n,' ');
  write( (byte_t*)c.data() ,  n );
}

void edfz2_t::writestring( const double & s , int n )
{
  std::string c = Helper::dbl2str_fixed(s,n);
  c.resize(n,' ');
  write( (byte_t*)c.data() , n  );
}

