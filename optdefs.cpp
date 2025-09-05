
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

#include "luna.h"
#include "db/db.h"
#include <iomanip>

#include "optdefs.h"

extern globals global;

optdefs_t::optdefs_t()
{
  odesc.clear();
  otype.clear();
}

void optdefs_t::add( const std::string & opt , const opt_type_t & type , const std::string & desc )
{
  odesc[ opt ] = desc;
  otype[ opt ] = type;
}

opt_type_t optdefs_t::get_type( const std::string & t ) const
{
  std::map<std::string,opt_type_t>::const_iterator ii = otype.find( t );
  if ( ii == otype.end() ) return OPT_UNDEFINED_T;
  return ii->second;
}

std::string optdefs_t::get_desc( const std::string & t ) const
{
  std::map<std::string,std::string>::const_iterator ii = odesc.find( t );
  if ( ii == odesc.end() ) return ".";
  return ii->second;
}

std::vector<std::string> optdefs_t::get_opts() const
{
  std::vector<std::string> res;
  std::map<std::string,std::string>::const_iterator ii = odesc.begin();
  while ( ii != odesc.end() )
    {
      res.push_back( ii->first );
      ++ii;
    }
  return res;    
}





