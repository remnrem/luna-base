
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

#ifndef __GC_H__
#define __GC_H__

#include "stats/Eigen/Dense"
#include <map>
#include <vector>

struct edf_t;
struct param_t;
struct signal_list_t;

void gc_wrapper( edf_t & edf , param_t & param );

struct gc_t { 
  
  gc_t( const Eigen::MatrixXd & X ,
	const signal_list_t & signals ,
	int sr , 
	double timewin_ms , 	
	double order_ms ,
	const std::vector<double> * frqs = NULL , 
	int compute_bic = 0 
	);
  
  void report();
  
  // time-domain
  double y2x, x2y;

  // frequency-stratified  
  std::map<double,double> tf_y2x, tf_x2y;
  
};



struct armorf_t {

  armorf_t( const Eigen::MatrixXd & X , 
	    const int Nr , 
	    const int Nl , 
	    const int p ) ;
  
  Eigen::MatrixXd coeff;
  Eigen::MatrixXd E;
  Eigen::MatrixXd K;
  
};



#endif
