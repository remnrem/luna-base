
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

#ifndef __LUNA_LGBM_MASSOC_H__
#define __LUNA_LGBM_MASSOC_H__

#include "lgbm/lgbm.h"
#include "helper/helper.h"
#include "stats/Eigen/Dense"
#include "stats/matrix.h" // legacy, for TLOCK inputs

#include <string>
#include <limits>

struct param_t;

// like assoc_t, but simplified for a matrix-import
// which can be internally passed too.  Main use-case
// to pass time-series information from TLOCK to be
// saved as a binary matrix, for subsequent training;
// we assume the feature matrix is ONLY features --
// i.e. all phenotype labels, and all ID codes (train/valid/test)
// will be passed separately via text

struct massoc_t {
  
  // when called externally, i.e. for testing/prediction
  massoc_t( param_t & );
  
  // when called internally (e.g. from TLOCK, where we just have one big feature matrix
  
  massoc_t( const std::string & id ,
	    const std::vector<std::string> & rowids ,
	    const std::vector<std::string> & eids ,
	    const std::vector<std::string> & colids ,
	    const Data::Matrix<double> & X ,
	    const std::string & filename );
  
private:


  int mode;

  //
  // inputs
  //

  // load train/valid/test status
  void attach_ids( param_t & param );
  
  // prune missing training/validation obs
  void prune();
  void prune1( const int n,
	       const std::vector<bool> & missing ,
	       std::vector<std::string> * iids,
	       std::vector<std::string> * ids,
	       std::vector<std::string> * eids,
	       Eigen::MatrixXd * X,
	       std::vector<double> * Y );
  
  // load all features
  void load( const std::string & , const int force_destin = 0 );

  // save a feature dataset (convenience wrapper for single iid case)
  void save( const std::string & ,
	     const std::vector<std::string> & rowids ,
	     const std::vector<std::string> & cntids ,
	     const std::vector<std::string> & colids , 
	     const Eigen::MatrixXd & X ,
	     const std::string & filename );

  // save a feature dataset (,ay have >1 IID, i.e. if from split)
  void save( const std::vector<std::string> & ,
	     const std::vector<std::string> & rowids ,
	     const std::vector<std::string> & cntids ,
	     const std::vector<std::string> & colids , 
	     const Eigen::MatrixXd & X ,
	     const std::string & filename );

  // split train+valid and test into two separate files
  void split( const std::string & , const std::string & ,
	      const std::string & , const std::string & ,
	      const std::set<std::string> * = NULL );

  // cbind() for two datasets (must have identical row count and IIDs and EIDs)
  // will append the IDs to each other
  void merge( const std::string & outfile );

  // load phenotypes/labels
  void attach_phenotypes( param_t & param );

  //
  // train model
  //

  void train( param_t & );

  void save_model( param_t & );


  //
  // prediction
  //

  void load_model( param_t & );
  
  void predict( param_t & );
  
  void SHAP( param_t & );
  
  const double NaN_value = std::numeric_limits<double>::quiet_NaN();

  //
  // members
  //
  
  static lgbm_t lgbm;

  /* // indiv-IDs (i.e. for splits, outputs) */
  /* std::vector<std::string> iids; */

  // col/row IDs (i.e. tracking indiv + event-level info)
  std::vector<std::string> vars;
  std::vector<std::string> ids;

  /* int nv; */
  /* int ni_train; */
  /* int ni_valid; */
  /* int ni_test; */

  // pool of potential IIDs
  std::set<std::string> training_pool;
  std::set<std::string> validation_pool;
  std::set<std::string> test_pool;


  // IIDs in data
  std::vector<std::string> training_iids;
  std::vector<std::string> validation_iids;
  std::vector<std::string> test_iids;

  // ID in data (i.e. strata/spindle type)
  std::vector<std::string> training_ids;
  std::vector<std::string> validation_ids;
  std::vector<std::string> test_ids;

  // Event count in data (i.e. 1,2,3 etc)
  std::vector<std::string> training_eids;
  std::vector<std::string> validation_eids;
  std::vector<std::string> test_eids;
    
  // indiv by var matrix
  Eigen::MatrixXd Xtrain;
  Eigen::MatrixXd Xvalid;
  Eigen::MatrixXd Xtest;

  // DV (QT or binary)
  std::vector<double> Ytrain;
  std::vector<double> Yvalid;
  std::vector<double> Ytest;

  // variable (to be pulled from ivars)
  std::string phenotype_label;
  
  // Misc options
  bool allow_missing_values;
  
};

  
#endif
#endif
