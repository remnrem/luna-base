
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

#ifndef __KMEANS_H__
#define __KMEANS_H__

#include "stats/matrix.h"
#include "Eigen/Dense"
#include <vector>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <limits>


struct kmeans_t { 
  
  int n;

  double between;
  double within;
  Data::Vector<double> within_ss;
  
  struct point_t { 

    point_t() { } 
    
    point_t(const int n) { x.resize(n); group = 0; }
    
    point_t( const Data::Vector<double> & d ) 
    {
      group = 0;
      x.resize( d.size() );
      for (int i=0; i<d.size(); i++) x[i] = d[i];
    }
    
    void clear() { 
      group = 0; 
      for (int i=0;i<x.size();i++) x[i] = 0;
    }
    
    // 'group' used as temp n-counter
    void scale() {
      for (int i=0;i<x.size();i++) x[i] /= group;
    }

    // 'add' for centroid tracking
    void add( point_t & a )
    {
      group++;
      for (int i=0;i<x.size();i++) x[i] += a.x[i];
    }

    std::vector<double> x;
    int group;
  };
  

  double randf(double m);
  
  double dist2( const point_t & a, const point_t & b );

  int nearest( const point_t & pt, 
	       const std::vector<point_t> & cent,
	       double * d2 ,
	       int * lim = NULL 
	       );
  
  void kpp( std::vector<point_t> & pts,  
	    std::vector<point_t> & cent );


  //
  // primary function
  //
  
  Data::Matrix<double> lloyd( const Data::Matrix<double> & X , int nk , std::vector<int> * sol = NULL );
    
  std::vector<point_t> lloyd( std::vector<point_t> & pts, int nk );

  void variance_explained( const std::vector<point_t> & pts , const std::vector<point_t> & cent );

  
  //
  // Tests
  //

  Data::Matrix<double> kmeans( const Data::Matrix<double> & X , const int nk , std::vector<int> * sol = NULL );

  void test2();

};




// modified K-means for EEG, following ...

// Pascual-Marqui, R. D., Michel, C. M., & Lehmann, D. (1995).
//         Segmentation of brain electrical activity into microstates: model
//         estimation and validation. IEEE Transactions on Biomedical
//         Engineering.


struct modkmeans_out_t {

  // .Z_all    - Cell containing the microstate activations for each number of
  //             microstates defined in K_range. The dimensions of each cell
  //                 is (K x samples).
  
  Eigen::MatrixXd Z;
    
  //     .A_all    - Cell containing the spatial distribution for each number of
  //                 microstates defined in K_range. The dimensions of each cell
  //                 is (channels x K).
  
  Eigen::MatrixXd A;
  
  //     .L_all    - Cell containing the labels for each number of microstates
  //                 defined in K_range. The dimensions of each cell is
  //                 (1 x samples).
  
  std::vector<int> L;
  
  //     .R2       - Explained variance for the best solution of each K in K_range.
  double R2;

  double sig2;
  
  //     .sig2_modk     - Noise variance for the best solution of each K in K_range.  
  double sig2_modk;
    
  //     .sig2_modk_mcv - Modified predictive residual variance of best solution
  //                     of each K in K_range.  
  double sig2_modk_mcv;
  
  //     .MSE      - Mean squared error for the best solution of each K in K_range.  
  double MSE;

  // number of iterations
  int iter;
};


struct modkmeans_all_out_t {

  // final, optimal solution
  int K;

  // A_opt   - Spatial distribution of microstates (channels x K)                                                                                                                     
  Eigen::MatrixXd A;
  
  // L_opt   - Label of the most active microstate at each timepoint (1 x samples).                                                                                                   
  std::vector<int> L;
  
  //
  // Verbose info for each K
  //

  std::map<int,modkmeans_out_t> kres;
  
};



struct modkmeans_t {

  modkmeans_t( const std::vector<int> & ks ,
	       const bool normalize = false ,
	       const int nreps = 10 ,
	       const int max_iterations = 1000 ,
	       const double threshold = 1e-6 ,
	       const bool verbose = false )
    : ks(ks),
    normalize(normalize) ,
    nreps(nreps) ,
    max_iterations(max_iterations),
    threshold(threshold) ,
    verbose(verbose)
  {
    
  }

  modkmeans_all_out_t fit( const Data::Matrix<double> & );


private:

  modkmeans_out_t segmentation( const Eigen::MatrixXd & , 
				int K , 
				double const1 );
  

  double eigen_correlation( const Eigen::VectorXd & a , const Eigen::VectorXd & b )
  {
    const int n = a.size();
    if ( b.size() != n ) return 0;
    if ( n < 2 ) return 0;    
    // r = cov(1,2) / sqrt( var(1).var(2) )    
    Eigen::MatrixXd mat(n,2);
    mat.col(0) = a;
    mat.col(1) = b;
    Eigen::MatrixXd centered = mat.rowwise() - mat.colwise().mean();
    Eigen::MatrixXd cov = (centered.adjoint() * centered) / double(mat.rows() - 1);
    return cov(0,1) / sqrt( cov(0,0) * cov(1,1) ) ;    
  }

  
  //
  // Members
  //

  Eigen::MatrixXd X;
  //Data::Matrix<double> X;
  
  // K to try
  std::vector<int> ks;

  // normalize EEG with average channel STD
  bool normalize;

  // # of times to randomly restart the algorithm (default = 10)
  int nreps;

  // (default 1000)
  int max_iterations;

  // default 1e-6
  double threshold;

  // not used:: always use CV
  int fitmeas; // 0 = CV (default)
  
  bool verbose;

  // optional smoothing; # sample on each side
  int b; 

  // smoothing weightt; def = 5  
  double lambda;

  //
  // Outputs
  //

  
  
};



#endif




