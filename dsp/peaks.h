
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

#ifndef __PEAKS_H__
#define __PEAKS_H__

struct edf_t;
struct param_t;
struct signal_list_t;

#include <vector>

namespace dsptools
{
  void peaks( edf_t & edf , param_t & param );
}

  
struct peaks_t {
  
 peaks_t()
  {
    // option defaults
    max = true;
    min = false;    
    ignore_clipped = true;
    th_clipped = 3; // 3 or more equal sample points == 'clipped'
    percentile = 0; // default --> all; expecting from 1..100  (100 implies all too)
  }

  // find peaks, and save them to the pre-specified cache  
  void detect( const std::vector<double> * x , const std::vector<int> * sp = 0 );
  
  // find max and/or min
  bool max;
  bool min;
  
  double percentile;

  // record clipped points
  bool ignore_clipped;
  // number of contiguous samples to call something clipped
  int th_clipped;  
  // to testing clipping
  const double EPS = 1e-6;

  // peak locations/values
  std::vector<int> pk;
  std::vector<double> values;
  std::vector<bool> ismin;
};

#endif
