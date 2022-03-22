
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

#ifdef HAS_LGBM

#ifndef __LUNA_LGBM_ASSOC_H__
#define __LUNA_LGBM_ASSOC_H__

#include "lgbm/lgbm.h"
#include "helper/helper.h"
#include "stats/Eigen/Dense"

#include <string>
#include <limits>

struct param_t;

struct assoc_t {
  
  assoc_t( param_t & );

private:

  void attach_covariates( param_t & param );  
  
  //
  // training
  //
  
  void attach_phenotypes( param_t & param );

  void attach_test_phenotypes( param_t & param );
    
  void attach_ids( param_t & param );

  void import_training( param_t & );
  
  void save( param_t & );
  
  void load( param_t & );
  
  void train( param_t & );

  void save_model( param_t & );

  void save_varlist( const std::string & );
  
  //
  // prediction
  //

  // LGBM model and model.vars files
  void load_model( param_t & );

  void load_varlist( const std::string & );
  
  void import_testdata( param_t & );
  
  void predict( param_t & );

  void SHAP( param_t & );
  
  const double NaN_value = std::numeric_limits<double>::quiet_NaN();

  //
  // members
  //
  
  static lgbm_t lgbm;

  // variables (model_ are saved w/ model to help parse test data)
  std::vector<std::string> varlist; // PSD_CH_C3_F_22
  std::set<std::string> model_vars; // e.g. PSD
  std::set<std::string> model_strats; // e.g. F 

  // indiv by var matrix
  Eigen::MatrixXd Xtrain;
  Eigen::MatrixXd Xvalid;
  Eigen::MatrixXd X; // test/prediction

  // DV (QT or binary)
  Eigen::VectorXd Ytrain;
  Eigen::VectorXd Yvalid;
  Eigen::VectorXd Y; // test/prediction

  // variables
  std::string phenotype_label;
  
  // IDs
  std::set<std::string> valid_set;
  std::vector<double> train_phe, valid_phe, test_phe;
  std::vector<std::string> train_ids, valid_ids, test_ids;

  //
  // Misc options
  //

  bool allow_missing_values;
  
};

// struct collector_t { 
//   std::string var;
//   std::map<
  
// }
  
#endif
#endif
