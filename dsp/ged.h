
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

#ifndef __GED_H__
#define __GED_H__

struct edf_t;
struct param_t;

#include "stats/Eigen/Dense"
#include <vector>

void ged_wrapper( edf_t & , param_t & param );

// narrow-band signal vs broad-band 
void ged_runmode1( edf_t & , param_t & , Eigen::MatrixXd & X , int );

// subset of time-points (e.g. peaks) versus another subset (or all)
void ged_runmode2( edf_t & , param_t & , Eigen::MatrixXd & X , const std::vector<uint64_t> * tp , int sr );

struct ged_t
{

  ged_t() { }
  
  void data( const Eigen::MatrixXd & Sd, const Eigen::MatrixXd & Rd ) ;

  void covar( const Eigen::MatrixXd & S_, const Eigen::MatrixXd & R_ )
  {
    S = S_;
    R = R_;
  }

  void calc();
  
  
  // covariance matrices
  Eigen::MatrixXd S;

  Eigen::MatrixXd R;

  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es;
  
  // (unsorted) eigenvectors
  Eigen::MatrixXd W;

  // (unsorted) eigenvalues
  Eigen::VectorXd L;

  // index of largest eigenvalue
  int largest_idx;

  Eigen::VectorXd map( const int e , const Eigen::MatrixXd & C , int * );

  Eigen::VectorXd time_series( const int e , const Eigen::MatrixXd & D , const int );
  
  
};


#endif
