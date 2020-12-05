
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

#ifndef __TLOCK_H__
#define __TLOCK_H__

struct edf_t;
struct param_t;
struct signal_list_t;

#include <vector>
#include "stats/matrix.h"

namespace dsptools
{
  void tlock( edf_t & edf , param_t & param );
}

  
struct tlock_t {
  
  tlock_t( const std::vector<double> & t , const int  );
  
  void add( const std::vector<double> * x , 
	    const int , const int , 
	    const bool take_log = false , 
	    const int angle_bin = 0 );
  
  Data::Vector<double> average() const ;

  Data::Matrix<double> angles() const ;
  
private:

  // for regular means::: track the whole matrix (just in case we want to ) 
  //   each row is a sample-point in the interval window
  //   each column is a new epoch
  
  // for angles: the matrix will be 360 / angle_bin wide
  //   each row is a sample-point in the interval window
  //   each column is an angle bin: do not track individual epochs
  //   but rather do the summatation in-place

  Data::Matrix<double> X;

  // time axis
  std::vector<double> t;

  // normalisation points
  int norm_points;
};


#endif
