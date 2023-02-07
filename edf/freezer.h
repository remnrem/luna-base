
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

#ifndef __FREEZER_H__
#define __FREEZER_H__

#include <vector>
#include <string>
#include <map>

struct edf_t;

struct freezer_t
{
  
  freezer_t() { }
  
  void freeze( const std::string & s , edf_t & );
  
  bool thaw( const std::string & s , edf_t * , bool clean = false );
  
  void clean( const std::string & s );
  
private:
  
  std::map<std::string,edf_t*> store;

  void edf2edf( const edf_t & from , edf_t & to );

  
};

#endif
