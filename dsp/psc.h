
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

#ifndef __PSC_H__
#define __PSC_H__

struct edf_t;
struct param_t;

#include "stats/Eigen/Dense"
#include "helper/helper.h"

#include <vector>
#include <set>

struct psc_t {

  void construct( param_t & );

  void attach( param_t & );
  
  void project( edf_t & edf , param_t & );
  
  static void clear_proj() 
  {
    vname.clear();
    means.resize(0);
    sds.resize(0);
    W.resize(0);
    V.resize(0,0);
  }

  // members

  static std::vector<std::string> vname;  
  static Eigen::Array<double, 1, Eigen::Dynamic> means;
  static Eigen::Array<double, 1, Eigen::Dynamic> sds;

  // number of PSC
  int nc;  

  static Eigen::VectorXd W;
  static Eigen::MatrixXd V;   
  const double EPS = 1e-6;
  
};


struct psc_sort_t {

  psc_sort_t( int idx , double value ) : idx(idx) , value(value) { } 

  int idx;

  double value;

  bool operator<( const psc_sort_t & rhs ) const 
  {
    return value < rhs.value;
  }

  static std::vector<int> quantile( const std::set<psc_sort_t> & d , const int q )
  {
    
    const int n = d.size();
    const int nq = n / q ;
    // need to add +1 to this many groups?, i.e. if not a perfect division
    int ex = n - ( q * nq );
    
    std::vector<int> r(n);
    
    int curr_q = 0;
    int curr_n = 0;
    
    std::set<psc_sort_t>::const_iterator qq = d.begin();
    while ( qq != d.end() )
      {
	// sanity check
	if ( qq->idx < 0 || qq->idx >= n ) 
	  Helper::halt( "internal error in psc_t" );
	
	r[ qq->idx ] = curr_q;
	
	// updates?
	
	++curr_n;
	
	if ( curr_n == nq + ( ex > 0 ? 1 : 0 ) )
	  {
	    ++curr_q;
	    curr_n = 0;
	    --ex; // we've used up this remainder
	  }
	
	++qq;
      }
    
    return r;
  }
};



#endif
