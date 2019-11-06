
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

#ifndef __LUNA_SL_H__
#define __LUNA_SL_H__

struct edf_t;
struct param_t;
struct clocs_t;
struct signal_list_t;

#include <vector>

#include "stats/matrix.h"


struct sl_t { 
  
  sl_t( const clocs_t & clocs , const signal_list_t & signals , 
	int m_ = 4 , int order_ = 10 , double lambda_ = 1e-5 ); 
  
  bool apply( const Data::Matrix<double> & inp , Data::Matrix<double> & out );

 private:
  
  int m; 
  
  int order;
  
  double lambda;

  Data::Matrix<double> G;
  
  Data::Matrix<double> invG;

  Data::Matrix<double> H;

  std::vector<double> GsinvS;
  
  double sumGsinvS;
  
}; 


namespace dsptools
{
  
  void surface_laplacian_wrapper( edf_t & edf , param_t & param );
    
}

#endif
