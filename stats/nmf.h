
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

#ifndef __NMF_H__
#define __NMF_H__

#include "stats/Eigen/Dense"
#include "helper/helper.h"

#include <vector>
#include <map>

struct nmf_t {
  
  nmf_t( const Eigen::MatrixXd & V ,
	 const int maxiter = 500 ,
	 const double EPS = 0.00001 ); 

  void factorize( const int num_sources );
    
  // V = W H 
  Eigen::MatrixXd V;
  Eigen::MatrixXd W;
  Eigen::MatrixXd H;
  std::vector<int> rows;
  std::vector<bool> included;
  int maxiter;  
  double EPS;
  int iter;
  
};

#endif
