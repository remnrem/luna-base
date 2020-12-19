
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
#include <vector>

// https://rosettacode.org/wiki/K-means%2B%2B_clustering
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
	       double * d2 );
  
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

  

#endif
