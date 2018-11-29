
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

#include "edf.h"
#include "helper/helper.h"
#include "miscmath/miscmath.h"

void edf_t::covar( const std::string & s1 , const std::string & s2 )
{

  //
  // Attach signals
  //
  
  signal_list_t signals1 = header.signal_list( s1 );
  signal_list_t signals2 = header.signal_list( s2 );

  if ( signals1.size() == 0 || signals2.size() == 0 ) 
    Helper::halt( "covar function requires both signals1/signals2 parameters" );
    

}

