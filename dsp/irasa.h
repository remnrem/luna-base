
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

#ifndef __IRASA_H__
#define __IRASA_H__

#include <vector>

struct edf_t;
struct param_t;

void irasa_wrapper( edf_t & edf , param_t & param );

struct irasa_t { 
  
  irasa_t( edf_t & edf ,
	   const std::vector<double> & x ,
	   const int sr ,
	   const double epoch_sec,
	   const int ne, 
	   const double h_min ,
	   const double h_max ,
	   const int h_cnt , 
 	   const double lwr,
	   const double upr,
	   const double seg_sec,
	   const double overlap_sec,
	   const int converter ,
	   const bool epoch_lvl_output ,
	   const bool logout , 	  
	   const std::vector<double> & slope_range , 
 	   const double slope_outlier ,
	   const int window , 
	   const bool segment_median ,
	   const bool epoch_median 	   
	   );

  
  int n;
  std::vector<double> frq;
  std::vector<double> periodic;
  std::vector<double> aperiodic;
  std::vector<double> aperiodic_raw;
  
};


#endif
