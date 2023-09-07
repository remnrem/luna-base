
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

#ifndef __INSERT_H__
#define __INSERT_H__

#include <vector>
#include <string>
#include <map>
#include "stats/Eigen/Dense"

struct edf_t;
struct param_t;

struct edf_inserter_t
{
  
  edf_inserter_t( edf_t & edf , param_t & param );
  
  void insert( edf_t & edf , edf_t & edf2 , const std::string & siglabel ,
	       const double offset , const double * stretch,
	       const std::string annot_label ); 
  
};

#endif
