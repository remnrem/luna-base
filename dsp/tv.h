
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

#ifndef __TV_H__
#define __TV_H__

struct edf_t;
struct param_t;

#include <vector>
#include <stdint.h>

namespace dsptools 
{  
  
  // ultimate wrapper
  void tv( edf_t & , param_t & );
  
  std::vector<double> TV1D_denoise_copy(const std::vector<double> & input, const double lambda);

  void TV1D_denoise(std::vector<double> & input, const double lambda);
  
  // alternate method, also w/ sparsity + smoothness
  void fused_lasso(double* input, double* output, const int width, const double lambda, const double mu);


  std::vector<double> tv( const std::vector<double> & );
		  
}

#endif
