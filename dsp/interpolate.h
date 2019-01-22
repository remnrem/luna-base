
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

#ifndef __INTERPOLATE_H__
#define __INTERPOLATE_H__

struct edf_t;
struct param_t;
struct topo_t;

#include <vector>

#include "stats/matrix.h"

namespace dsptools
{
  void leave_one_out( edf_t & edf , param_t & param );

  Data::Matrix<double> interpolate2D( const std::vector<double> & x , 
				      const std::vector<double> & y , 
				      const std::vector<double> & z , // values 
				      const double xmin , 
				      const double xmax ,
				      const int    nx , 
				      const double ymin , 
				      const double ymax ,
				      const int    ny ) ;
  
  void interpolate2D( topo_t * topo , const std::vector<double> & );
  
  
}

#endif
