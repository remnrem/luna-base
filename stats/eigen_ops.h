
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

#ifndef __EIGEN_OPS_H__
#define __EIGEN_OPS_H__

#include "stats/Eigen/Dense"

namespace eigen_ops { 

  void random_normal( Eigen::MatrixXd & m );  

  void scale( Eigen::MatrixXd & m , bool );  
  
}

  // inline double min( double a, double b) { return a<b ? a : b; } 

  // inline double abs( double a, double b) { return x < 0 ? -x : x ; }

  // void mat_zeroize(Eigen::MatrixXd & , const int rows = 0 , const int cols = 0);
  
  // void vect_zeroize(Data::Vector<double> & , const int cols = 0 );
  
  // void vect_apply_fx(Eigen::ArrayXd & v, double (*fx)(double,double), double par);
  
  // void mat_apply_fx(Eigen::MatrixXd & M, double (*fx)(double,double), double par);
  
  // void mat_mean_rows(Eigen::MatrixXd & M, Eigen::RowArrayXd & v);
  
  // double mat_max_diag(Eigen::MatrixXd & M );
  
  // double mat_max_abs_diag(Eigen::MatrixXd & M );
  
  // void mat_center( Eigen::MatrixXd & M,  Eigen::RowArrayXd & means );
  
  // void mat_decenter(Eigen::MatrixXd & M, Eigen::RowArrayXd & means);
  

#endif 
