
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

#ifndef __WRAPPERS_H__
#define __WRAPPERS_H__

struct edf_t;
struct param_t;

#include <vector>
#include <stdint.h>
#include <cstddef>

namespace dsptools 
{  
  
  void cwt( edf_t & , param_t & );

  void hilbert( edf_t & , param_t & );
  
  void run_cwt( const std::vector<double> & data , const int ,
		const double fc , const int num_cycles , 
		std::vector<double> * mag , 
		std::vector<double> * phase = NULL ); 
  
  void run_hilbert( const std::vector<double> & data , const int Fs , 
		    const double flwr , const double fupr , const double ripple , const double tw , 
		    std::vector<double> * mag , 
		    std::vector<double> * phase = NULL , 
		    std::vector<double> * frequency = NULL ); 
  
}

#endif
