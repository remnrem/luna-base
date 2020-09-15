
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

#include "stats/matrix.h"

#include <vector>

struct psc_t {

  void construct( param_t & );

  void attach( param_t & );
  
  void project( edf_t & edf , param_t & );


  // data
  static std::vector<std::string> vname;
  static Data::Vector<double> vmean;
  static Data::Vector<double> vsd;

  // number of PSC
  int nc;
  
  static Data::Vector<double> W;

  static Data::Matrix<double> DW; // 1/W as diag matrix

  static Data::Matrix<double> V;

  const double EPS = 1e-6;

  
};

#endif
