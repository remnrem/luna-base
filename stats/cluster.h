

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

