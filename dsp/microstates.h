
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

#ifndef __MICROSTATES_H__
#define __MICROSTATES_H__

struct edf_t;
struct param_t;
struct signal_list_t;

#include <vector>
#include "stats/matrix.h"

namespace dsptools
{
  void microstates( edf_t & edf , param_t & param );
}

  
struct microstates_t {
  
  microstates_t( const Data::Matrix<double> & X , 
		 const signal_list_t & signals , 
		 const std::vector<int> & ks ); 
     
  Data::Vector<int> solution( const int k );

};

#endif
