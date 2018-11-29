
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


#ifndef __DSP_REDUCE_H__
#define __DSP_REDUCE_H__

#include <stdint.h>
#include <vector>

struct reduce_t 
{
  reduce_t( const std::vector<double> & , int np );
  reduce_t( const std::vector<double> * , const std::vector<uint64_t> * , uint64_t , uint64_t , int np );
  bool reduced;
  std::vector<double> hi, lo;
  std::vector<double> mean;
  std::vector<double> sd;
  std::vector<int> n;
};

#endif
