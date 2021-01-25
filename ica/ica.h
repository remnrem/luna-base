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

#ifndef __EIGEN_ICA_H__
#define __EIGEN_ICA_H__

#include "stats/Eigen/Dense"
#include "helper/helper.h"

// C/C++ & Eigen implementation of the fastICA R package
// https://cran.r-project.org/web/packages/fastICA/

struct edf_t; 
struct param_t; 

void eigen_ica_wrapper0( edf_t & , param_t & ); 

struct eigen_ica_t {
  
  eigen_ica_t( Eigen::MatrixXd & X , int compc )
  {        
    maxit = 200;
    tol = 0.0001;
    alpha = 1;
    row_norm = false;
    if ( ! proc(X,compc) ) Helper::halt( "problem in eigen_ica_t" );    
  }
  
  Eigen::MatrixXd K;
  Eigen::MatrixXd W;
  Eigen::MatrixXd A;
  Eigen::MatrixXd S;
  
  bool proc( Eigen::MatrixXd & , int compc );

  void fastICA( Eigen::MatrixXd & X , 
		const int compc , 
		Eigen::MatrixXd & W , 
		Eigen::MatrixXd & A , 
		Eigen::MatrixXd & K , 
		Eigen::MatrixXd & S ) ;

  Eigen::MatrixXd ica_parallel( const Eigen::MatrixXd & X ,
				const int nc );
				

  
  int    maxit;
  double tol;
  int    alpha;
  bool   row_norm; // standardize input 
  
};


#endif
