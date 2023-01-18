
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
#include <vector>
#include <map>

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

namespace eigen_ops { 
    
  void random_normal( Eigen::MatrixXd & m );  
  
  bool scale( Eigen::Ref<Eigen::MatrixXd> m , const bool,  const bool , const bool ignore_invariants = false , std::vector<int> * zeros = NULL );  

  bool robust_scale( Eigen::Ref<Eigen::MatrixXd> m , const bool , bool , double , bool second_rescale = true , const bool ignore_invariants = false , std::vector<int> * zeros = NULL );    

  Eigen::VectorXd percentile_scale( const Eigen::VectorXd & x , const double pct , const int nsegs );
  
  void accumulate( Eigen::Ref<Eigen::VectorXd> m , const int ctype );

  void deriv( Eigen::Ref<Eigen::VectorXd> m , const int hw );
  
  double sdev( const Eigen::VectorXd & x );

  double between_within_group_variance( const std::vector<std::string> & grp , const Eigen::VectorXd & x );
  
  // remove linear trend of each column by linear regression
  bool detrend( Eigen::Ref<Eigen::MatrixXd> m );

  Eigen::VectorXd median_filter( const Eigen::VectorXd & x , const int );

  Eigen::VectorXd tri_moving_average( const Eigen::VectorXd & x , int s , double w );
  
  Eigen::VectorXd moving_average( const Eigen::VectorXd & x , int s );
  
  Eigen::VectorXd unit_scale( const Eigen::VectorXd & x , double xmin , double xmax );

  Eigen::VectorXd unit_scale( const Eigen::VectorXd & x );

  std::vector<double> copy_vector( const Eigen::VectorXd & e );

  Eigen::ArrayXd copy_array( const std::vector<double> & e );

  std::vector<double> copy_array( const Eigen::ArrayXd & e );
  
  std::map<int,std::vector<double> > group_means( const Eigen::MatrixXd & x , const std::vector<int> & g );

  Eigen::VectorXd canonical_correlation( const Eigen::MatrixXd & X , const Eigen::MatrixXd & Y );

  Eigen::MatrixXd load_mat( const std::string & file ,
			    std::vector<std::string> * header = NULL ,
			    std::vector<std::string> * ids = NULL ,
			    std::vector<std::string> * labels = NULL );

  
  
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
