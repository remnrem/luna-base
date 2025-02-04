
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

#ifndef __AROUSALS_H__
#define __AROUSALS_H__

#include "stats/Eigen/Dense"
#include "intervals/intervals.h"

struct edf_t;

struct param_t;

struct arousals_t {
  
  arousals_t( edf_t & edf , param_t & param );
  
  bool hjorth( const Eigen::VectorXd & x ,
	       double * , double * , double * ,
	       const bool mean_center = false ) const;

private:

  std::vector<interval_t> combine( const std::vector<int> & , const std::vector<uint64_t> & , const double , const double ) const;
  
};


#endif
