
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

#ifndef __MERGE_HELPERS_H__
#define __MERGE_HELPERS_H__

#include <sstream>
#include <string>

bool fileExists( const std::string & f );
bool iequals(const std::string& a, const std::string& b);
bool imatch(const std::string& a, const std::string& b , unsigned int min = 0 );
std::string toupper( const std::string & s );
std::string sanitize( const std::string & s );
std::vector<std::string> parse(const std::string & item, const std::string & s , bool empty = false );
std::vector<std::string> char_split( const std::string & s , const char c , bool empty = false  );
std::vector<std::string> char_split( const std::string & s , const char c , const char c2 , bool empty = false );
std::vector<std::string> char_split( const std::string & s , const char c , const char c2 , const char c3 , bool empty = false );
void halt( const std::string & msg );
std::string expand( const std::string & f );
bool file_extension( const std::string & f, const std::string & ext );
std::string remove_extension( const std::string & f, const std::string & ext );
std::istream& safe_getline(std::istream& is, std::string& t);
std::string search_replace( const std::string & s , char a , char b );
bool str2int(const std::string & s , int * i);
bool str2dbl(const std::string & s , double * d);

template <class T>
bool from_string(T& t,
		 const std::string& s,
		 std::ios_base& (*f)(std::ios_base&))
{
  std::istringstream iss(s);
  return !(iss >> f >> t).fail();
}

#endif
