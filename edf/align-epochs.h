
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

#ifndef __ALIGN_EPOCHS_H__
#define __ALIGN_EPOCHS_H__

#include <vector>
#include <string>
#include <map>
#include "stats/Eigen/Dense"

struct edf_t;
struct param_t;

struct align_epochs_t
{
  
  align_epochs_t( edf_t & edf , param_t & param );
    
private:

  int ne;
  int ne2;
  int ns;

  std::map<int,Eigen::MatrixXd> X1, X2;
  std::vector<int> E1, E2;
  
  // track signal slots from EDF(loaded) -> EDF2(subset)
  std::vector<int> slot2;

  // final mapping: E2 --> E1 
  std::map<int,int> mapping;

  // second best choice
  std::map<int,int> mapping2;
  
  int best_match( const int e , double * , double * , int * ) const;

  double dist( const int e1 , const int e2 ) const; 

  double th;
  double th2; // second best 

  const double DEPS = 1e-20;

  bool assume_order;
  bool resolve_order;
  
  int verbose;
  int verbose2;
  
};

#endif
