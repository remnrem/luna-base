
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


#include "dsp/ssa.h"

#include <iostream>
#include "edf/edf.h"
#include "edf/slice.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "eval.h"
#include "db/db.h"

extern logger_t logger;
extern writer_t writer;

void dsptools::ssa_wrapper( edf_t & edf , param_t & param )
{
  
}


ssa_t::ssa_t( const std::vector<double> * x , const int l )
{
  
  const int n = x->size();
  
  // copy
  Eigen::VectorXd D = Eigen::VectorXd::Zero( n );
  for (int i=0; i<n; i++) D[i] = (*x)[i];
  
  // window length = l
  // columns in trajectory matrix
  const int k = n - l + 1;
  
  X = Eigen::MatrixXd::Zero( l , k );
  
  // stack
  for (int i=0; i<k; i++)
    X.col(i) = D.segment(i,l);

  std::cout << "X\n" << X << "\n";
  
}
