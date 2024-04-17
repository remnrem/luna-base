
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

#ifndef __LDA_H__
#define __LDA_H__

#include "stats/Eigen/Dense"
#include "helper/helper.h"

#include <vector>
#include <map>

typedef Eigen::Array<double, 1, Eigen::Dynamic> ArrayXd;

struct lda_model_t {
  
  lda_model_t()
  {
    valid = false;
    errmsg = "";
  }

   
  bool valid;

  std::string errmsg;
  
  ArrayXd prior;
  std::map<std::string,int> counts;
  Eigen::MatrixXd means;
  Eigen::MatrixXd scaling; // diagonal matrix
  int n;
  ArrayXd svd;
  std::vector<std::string> labels;

  ArrayXd prop_trace() const
  {
    ArrayXd t( svd.size() );

    double sum = 0;

    for (int i=0;i<svd.size();i++)
      {
	t[i] = svd[i] *	svd[i];
	sum += svd[i] * svd[i];
      }

    for	(int i=0;i<svd.size();i++)
      t[i] /= sum;

    return t;
  }
  
};


struct lda_posteriors_t {

  // cols = classes
  // rows = observations
  
  Eigen::MatrixXd pp;
  
  // most likely class label for each observation
  std::vector<std::string> cl;
  std::vector<int> cli; // as above, but index

  
};


class lda_t {
  
 public:
  
 lda_t( const std::vector<std::string> & y ,
	const Eigen::MatrixXd & X )	
   : y(y) , X(X)
  {
    
    tol = 1e-4;
    
    missing = "?";
    
  } 

  // support for adding extra columns (i.e. without caller having to remake X)
 lda_t( const std::vector<std::string> & y ,
	const Eigen::MatrixXd & X1 , 
	const Eigen::MatrixXd & X2 
	)
   : y(y) 
  {
    
    const int nr = X1.rows() ;
    if ( nr != X2.rows() ) Helper::halt( "internal error in lda_t" );
    
    const int nc1 = X1.cols();
    const int nc2 = X2.cols();
    
    X = Eigen::MatrixXd::Zero( nr , nc1 + nc2 );
    
    for (int c=0; c<nc1; c++)
      X.col( c ) = X1.col( c );
    for (int c=0; c<nc2; c++)
      X.col( nc1 + c ) = X2.col( c );
    
    tol = 1e-4;

    missing = "?";

    silent = false;
    
  }
  
  lda_model_t fit( const bool flat_priors = false , const  std::vector<std::string> * pr = NULL );

  static lda_posteriors_t predict( const lda_model_t & , const Eigen::MatrixXd & X );

  static lda_posteriors_t predict( const lda_model_t & model , const Eigen::MatrixXd & X , const Eigen::MatrixXd & X2 )
  {
    const int nr = X.rows() ;
    if ( nr != X2.rows() ) Helper::halt( "internal error in lda_t" );
    const int nc1 = X.cols();
    const int nc2 = X2.cols();
    Eigen::MatrixXd XX = Eigen::MatrixXd::Zero( nr , nc1 + nc2 );

    for (int c=0; c<nc1; c++)
      XX.col( c ) = X.col( c );
    for (int c=0; c<nc2; c++)
      XX.col( nc1 + c ) = X2.col( c );

    return predict( model , XX );
  }

  bool silent ; 
  
 private:
  
  std::vector<std::string> y;

  Eigen::MatrixXd X;

  double tol;
  
  std::string missing;
  
};

#endif
