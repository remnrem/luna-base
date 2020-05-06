
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


#ifndef __ACF_H__
#define __ACF_H__

#include <vector>

struct edf_t;
struct param_t;

struct acf_t {
  acf_t( const std::vector<double> & d , int maxlag = 0 ) { calc(d,maxlag); }
  void calc( const std::vector<double> & d , int maxlag = 0 );
  std::vector<double> acf() const { return r; }
  std::vector<double> r;
};

namespace dsptools { 
  void autocorr_channels( edf_t & edf , param_t & param );
}

#endif
