
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

#ifndef __PRED_H__
#define __PRED_H__

#include "stats/Eigen/Dense"
#include <map>
#include <vector>
#include <set>
#include <string>
#include "helper/helper.h"
#include "models/model.h"
#include "models/knn.h"

struct edf_t;
struct param_t;

struct prediction_t {

  prediction_t( edf_t & , param_t & );  

  double prediction();
  
private:

  void output() const;

  std::string id;
  
  prediction_model_t model;
  
  // derived feature vector (raw, normalized)
  Eigen::VectorXd X;
  Eigen::VectorXd Z;    

  // kNN for missing data / evaluate input feature typicality
  model_knn_t knn;
  Eigen::VectorXd D; 
  std::vector<bool> missing;
  std::vector<bool> missing2; // re-imputed
  
  // predicted value (raw, bias-adjusted)
  double y;
  
  double y1; 

};

#endif
