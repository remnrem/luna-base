
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

#ifndef __QDA_H__
#define __QDA_H__

#include "stats/Eigen/Dense"
#include "helper/helper.h"

#include <vector>
#include <map>

typedef Eigen::Array<double, 1, Eigen::Dynamic> ArrayXd;

struct qda_model_t {
  
  qda_model_t()
  {
    valid = false;
    errmsg = "";
  }

   
  bool valid;  
  std::string errmsg;
  
  ArrayXd prior;
  std::map<std::string,int> counts;
  ArrayXd rows;
  Eigen::MatrixXd means;
  std::vector<Eigen::MatrixXd> scaling;
  std::vector<double> ldet;
  int n; 
  std::vector<std::string> labels;
  
};


struct qda_posteriors_t {

  // cols = classes
  // rows = observations
  
  Eigen::MatrixXd pp;
  
  // most likely class label for each observation
  std::vector<std::string> cl;
  std::vector<int> cli; // as above, but index

  
};


class qda_t {
  
 public:
  
 qda_t( const std::vector<std::string> & y ,
	const Eigen::MatrixXd & X )	
   : y(y) , X(X)
  {
    
    tol = 1e-8;
    
    missing = "?";
    
  } 
  
  // support for adding extra columns (i.e. without caller having to remake X)
  qda_t( const std::vector<std::string> & y ,
	 const Eigen::MatrixXd & X1 , 
	 const Eigen::MatrixXd & X2 
	 )
    : y(y) 
  {
    
    const int nr = X1.rows() ;
    if ( nr != X2.rows() ) Helper::halt( "internal error in qda_t" );
    
    const int nc1 = X1.cols();
    const int nc2 = X2.cols();
    
    X = Eigen::MatrixXd::Zero( nr , nc1 + nc2 );
    
    for (int c=0; c<nc1; c++)
      X.col( c ) = X1.col( c );
    for (int c=0; c<nc2; c++)
      X.col( nc1 + c ) = X2.col( c );
    
    tol = 1e-4;

    missing = "?";

  }
  
  qda_model_t fit( const bool flat_priors = false );
  
  static qda_posteriors_t predict( const qda_model_t & , const Eigen::MatrixXd & X );

  static qda_posteriors_t predict( const qda_model_t & model , const Eigen::MatrixXd & X , const Eigen::MatrixXd & X2 )
  {
    const int nr = X.rows() ;
    if ( nr != X2.rows() ) Helper::halt( "internal error in qda_t" );
    const int nc1 = X.cols();
    const int nc2 = X2.cols();
    Eigen::MatrixXd XX = Eigen::MatrixXd::Zero( nr , nc1 + nc2 );

    for (int c=0; c<nc1; c++)
      XX.col( c ) = X.col( c );
    for (int c=0; c<nc2; c++)
      XX.col( nc1 + c ) = X2.col( c );

    return predict( model , XX );
  }

 private:
  
  std::vector<std::string> y;

  Eigen::MatrixXd X;
  
  double tol;
  
  std::string missing;
  
};

#endif
