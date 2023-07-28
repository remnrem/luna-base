
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


#ifndef __LUNA_XCORR_H__
#define __LUNA_XCORR_H__

#include "stats/Eigen/Dense"
#include <vector>

struct edf_t;
struct param_t;

namespace dsptools {
  void xcorr( edf_t & , param_t & );
}  

struct xcorr_t {

  xcorr_t( std::vector<double> a ,
	   std::vector<double> b ,
	   const int mxlag = 0 ,
	   const int center = 0 );
    
  std::vector<double> C;
  std::vector<int> lags;
  int mx;

    
};


#endif
