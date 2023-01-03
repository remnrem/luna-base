\
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

#ifndef __LUNA_POPS_INDIV_H__
#define __LUNA_POPS_INDIV_H__

struct edf_t;
struct param_t;

#include "pops/pops.h"


struct pops_sol_t {
  // epoch numbers
  std::vector<int> E;  
  // predicted stages
  std::vector<int> S;  
  // predictions
  Eigen::MatrixXd P;  
};

struct pops_indiv_t {
  
  // normal train/test route
  pops_indiv_t( edf_t & , param_t & );

  // special case #1: to eval against existing, external staging
  pops_indiv_t( edf_t & , param_t & , const std::string & f1 );
  // special case #2: to eval two external files:: no attached EDFs
  pops_indiv_t( param_t & , const std::string & f1 , const std::string & f2 );

  bool staging( edf_t & , param_t & );
  
  void save1( const std::string & , const std::string & );
  
  void level1( edf_t & );

  void level2( const bool quiet_mode = false );

  void apply_ranges(double,double);

  void apply_incexcvars();
  
  void ftr_summaries();
  
  void predict( const int iter = 0 );
  
  void SHAP();
  
  void apply_soap();
  Eigen::MatrixXd soap_X( bool * okay );
  double simple_soap( const Eigen::MatrixXd & X , const std::vector<int> & );
  void grid_soap();
  
  void apply_espriors( const std::string & f );

  static Eigen::VectorXd update_posteriors( const Eigen::VectorXd & posteriors ,
					    const Eigen::VectorXd & original_priors ,
					    const Eigen::VectorXd * new_priors ,
					    const Eigen::VectorXd * rescale = NULL );

  int update_predicted( std::vector<int> * cnts = NULL );  
  
  void summarize( pops_sol_t * sol = NULL );
  
  void print_confusion_matrix();
  
  void combine( std::vector<pops_sol_t> & sols ,
		int method ,
		double min_conf );

  void add_annots( edf_t & , const std::string & prefix = "p" );

  // track
  edf_t * pedf;
  
  // trainer/target?
  bool trainer;

  // has manual staging?
  bool has_staging;
  
  // number of epochs
  int ne;

  // total number of epochs
  int ne_total;
  
  // level 1 features
  Eigen::MatrixXd X1;

  // full level 1 features (i.e. copy before setting NaNs)
  //  - this is used in SOAP
  Eigen::MatrixXd X1f;
    
  // staging
  std::vector<int> S;
  std::vector<int> Sorig;

  // epoch number
  std::vector<int> E;

  // predictions
  Eigen::MatrixXd P;
  std::vector<int> PS; // (final) predicted stages
  


  //
  // Helper function -- compare this individual to an external set of calls
  //
  
  void eval_stages();
  
  //
  // binary I/O functions
  //
  
  inline static void bwrite( std::ofstream & O , const std::string & s );
  inline static void bwrite( std::ofstream & O , int i );
  inline static void bwrite( std::ofstream & O , double d );
  inline static std::string bread_str( std::ifstream & I );
  inline static int bread_int( std::ifstream & I );
  inline static double bread_dbl( std::ifstream & I );
  inline static void bskip_dbl( std::ifstream & I , const int n );
  inline static void bskip_int( std::ifstream & I , const int n );
  
  
};


#endif
#endif

