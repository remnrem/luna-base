
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

#ifndef __PRED_KNN_H__
#define __PRED_KNN_H__

#include "stats/Eigen/Dense"
#include <string>
#include "helper/helper.h"

struct model_knn_t {
  
  static void load( const std::string & f , int * rows = NULL , int * cols = NULL );
  
  // use k-nearest observations to fill in missing data
  
  static Eigen::VectorXd impute( const Eigen::VectorXd & f ,
				 const std::vector<bool> & missing ,
				 int k = 10 );

  static void clear();
  
private:

  static Eigen::MatrixXd X;  
  
};
  
#endif

