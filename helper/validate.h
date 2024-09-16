
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

#ifndef __VALIDATE_H__
#define __VALIDATE_H__

#include "defs/defs.h"
#include "edf/edf.h"


namespace Helper
{
  void validate_slist( param_t & );
  std::vector<std::tuple<std::string,std::string,bool> > validate_slist_lunapi_mode( const std::vector<std::tuple<std::string,std::string,std::set<std::string> > > & );
}

#endif
