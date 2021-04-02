
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

#ifndef __LINEDENOISE_H__
#define __LINEDENOISE_H__

#include <vector>

struct edf_t;
struct param_t;

// wrapper
namespace dsptools { 

  void line_denoiser( edf_t & edf , param_t & param );

}

// core function
std::vector<double> line_denosier( const std::vector<double> * x ,
				   const int Fs ,
				   const std::vector<double> & fl ,
				   const double w_noise ,
				   const double w_neighbour );
                                                                                       \

#endif
