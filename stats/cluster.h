
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

#ifndef __LUNA_CLUSTER_H__
#define __LUNA_CLUSTER_H__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <set>
#include <algorithm>
#include <cmath>

#include "stats/matrix.h"

// solution
struct cluster_solution_t { 

  // number of clusters in best solution (based on min silhouette)
  int k;

  // length n-obs, o-based cluster assignment
  std::vector<int> best;
  
/*   // silhouette per observation at best */
/*   std::vector<double> Si; */

/*   // average silhouette over all obs, for each K */
/*   std::vector<double> silhouette; */

};

// (naive) clustering routine
struct cluster_t {
  
  cluster_solution_t build( const Data::Matrix<double> & D );

  // Helper function: find the maximum distance between two clusters
  double cldist( const Data::Matrix<double> & , std::vector<int> &, std::vector<int> &);

  // Helper function: group average link
  double groupAvgLink( const Data::Matrix<double> &, std::vector<int> &, std::vector<int> &);
  
};


#endif
