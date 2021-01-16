
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

#include <vector>

struct psc_t {

  void construct( param_t & );

  void attach( param_t & );
  
  void project( edf_t & edf , param_t & );

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

#endif
