
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


#include "ica/ica.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "eval.h"
#include "db/db.h"
#include "stats/matrix.h"

extern writer_t writer;

extern logger_t logger;

// wrapper around libICA function

bool ica_t::proc( Data::Matrix<double> & X , int compc )
{

  int rows = X.dim1();

  if ( rows == 0 ) return false;

  int cols = X.dim2();
  
  // mean center
  Data::Vector<double> means( cols );
  mat_center( X , means );

  W.resize(compc, compc);
  A.resize(compc, compc);
  K.resize(cols , compc);
  S.resize(rows, cols);
  
  // compute ICA

  fastICA(X, compc, K, W, A, S);

  return true;
}
