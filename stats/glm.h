
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

#ifndef __GLM_H__
#define __GLM_H__

#include "matrix.h"
#include <vector>

class GLM {
  
 public:

  enum model_t { LOGISTIC , LINEAR } ;       
  
  GLM( const model_t & m ) 
    {
      model = m;
      nind = np = 0;
      nc = 0; cluster = false;
      all_valid = false;      
      t = 1; 
      standard_beta = false;
      RSS = -1;
      ci(0.95);
      vif(10);
    }
  
  // 'data-entry'
  void set( Data::Vector<double> & y , Data::Matrix<double> & x , 
	    std::vector<int> * cl = NULL , 
	    std::vector<bool> * mask = NULL );

  void ci(double ci);
  void vif(double vif);

  // fit
  bool fit() 
  { 
    if ( ! check_VIF() ) return false;
    return model == LOGISTIC ? fit_logistic() : fit_linear() ; 
  }
  
  // also add a version that performs joint tests, etc
  
  // look at results
  bool valid() const;
  void valid(const bool b) { all_valid = b; };
  
  Data::Vector<double> beta() const;
  Data::Vector<double> se() const;

  // assumes a single "test" has been set (default = 1st parameter)

  bool   test_valid() const;
  double test_var() const;
  double test_coef() const;
  double test_se() const;
  double test_pval() const;
  double test_statistic() const;
  double test_lower_ci() const;
  double test_upper_ci() const;

  // for permutation
  double statistic() const;
  
 private:

  model_t model;
  
  int nind, np;

  // 'test parameter', e.g. for permutation
  int t;

  // Dependent, independent variables

  // Strata/clustering information
  
  bool               cluster;
  std::vector<int>   clst;
  int                nc;
  
  // references to actual data (not stored in class)

  Data::Vector<double> Y;
  Data::Vector<double> pr;
  Data::Vector<double> V;
  Data::Matrix<double> X;

  
  // State variables

  bool all_valid;
  std::vector<bool> is_valid;
  
  //
  // Derived summary/model statistics
  //

  // beta and Sigma (co-efficients, covariance matrix)
  
  Data::Vector<double> coef; 
  Data::Matrix<double> S;
  

  //
  // Misc. utility functions
  //

  double linear_hypothesis( Data::Matrix<double> & , Data::Vector<double> & );

  void HuberWhite();

  bool check_VIF();

  //
  // linear models
  //

  bool fit_linear();
  bool fit_univariate_linear();
  bool fit_logistic();

  void set_variance();
  void standardise();
  
  double meanY, varY;
  double ci_zt;
  double vif_threshold;
  double RSS;
  bool standard_beta;
  Data::Vector<double> sig;
  Data::Vector<double> w;
  Data::Matrix<double> u;
  Data::Matrix<double> v;



 public:

  Data::Vector<double> get_var();
  Data::Vector<double> get_SE();


  
  std::string summary();
  bool display( Data::Vector<double> * coef = NULL , 
		Data::Vector<double> * se = NULL , 
		Data::Vector<double> * pvalue = NULL ,
		std::vector<bool> * mask = NULL , 
		Data::Vector<double> * lowci = NULL , 
		Data::Vector<double> * uprci = NULL , 
		Data::Vector<double> * statistic = NULL );

  double calc_RSS(); 
  double calc_rsqr();
  double calc_adj_rsqr();
  double calc_MallowC( GLM *);
  double calc_FTest( GLM * );

  //
  // logistic models
  //

  double get_loglik();

};

#endif
