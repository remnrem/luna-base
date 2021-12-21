
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

#ifndef __TCLST_H__
#define __TCLST_H__

struct edf_t;
struct param_t;
struct signal_list_t;

#include <vector>
#include "stats/matrix.h"
#include "stats/cluster.h"

namespace dsptools
{
  void tclst( edf_t & edf , param_t & param );
}


struct tclst_t {
  
  // list of samples : sp
  // window around each 
  // for N points, construct the N x N distance matrix from X (signal) and P (phase)
  //  
  tclst_t( const std::vector<Eigen::MatrixXd> * X ,
	   const std::vector<Eigen::MatrixXd> * P ,
	   const std::vector<std::string> & chs ,
	   const std::vector<double> & t , 
	   const int k1 , const int k2, 
	   const int hcK , 
	   const bool use_complex_dist ) ; 
  
  // for any one of P points, we will expect an interval for 1 or more channels
  //  -- create the P x P distance matrix based on signal amp & phase distance
  //     summed over all channels

  // number of internals (peaks)
  int n;

  // time axis
  std::vector<double> t;

  void cluster( const int k );
    
  // cluster solution: number of clustrers
  int k;

  // distance matrix
  Data::Matrix<double> D;

  // size
  // hierarchical cluster solution
  cluster_solution_t sol;

  // solution from K-means (K --> X maps) 
  std::map<int,Data::Matrix<double> > kmeans;
  std::map<int,std::vector<int> > ksol;
  std::map<int,double> varexp;

  // normalisation points
  int norm_points;

  // verbose mode
  bool verbose;

};


#endif
