
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


#ifndef __DFA_H__
#define __DFA_H__

#include "hilbert.h"
#include <vector>
#include <set>
#include <stats/Eigen/Dense>

namespace dsptools { 
  void dfa_wrapper( edf_t & edf , param_t & param ) ;
}

struct dfa_t {

  dfa_t();

  // default, 0.1s to 10s 
  void set_windows( double sr, double l = 0.1 , int m = 2 , int c = 100 ) ;

  void filter_hilbert( const double flwr1 , const double fupr1 ,
		       const double ripple1 , const double tw1 )    
  {
    flwr = flwr1;
    fupr = fupr1;
    ripple = ripple1;
    tw = tw1;
  }

  void proc( const std::vector<double> * d );

  double sr, flwr, fupr, ripple, tw;
  
  std::vector<double> w;
  std::vector<double> t;
  std::vector<double> fluctuations;
  std::vector<double> slopes;
  
};


#endif
