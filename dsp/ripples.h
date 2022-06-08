
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

#ifndef __RIPPLES_H__
#define __RIPPLES_H__

struct edf_t;
struct param_t;
struct signal_list_t;

#include <vector>

#include "intervals/intervals.h"

namespace dsptools {
  
  void ripple_wrapper( edf_t & edf , param_t & param );

}


struct ripple_t {

  ripple_t( const uint64_t start , const uint64_t stop ,
	    const int start_sp , const int stop_sp )
    : pos( start , stop ) , start_sp( start_sp ) , stop_sp( stop_sp ) , wgt(1.0) , frq(0.0) 
  { } 
  
  interval_t pos;
  int start_sp, stop_sp;
  double wgt;
  double frq;
  
};

  
struct ripples_t {
  
  // detector
  ripples_t( const std::vector<double> & x ,
	     const std::vector<uint64_t> & tp ,
	     const int sr ,
	     const double flwr ,
	     const double fupr ,
	     const double kwin_ripple ,
	     const double kwin_tw , 
	     const int hfbands , 
	     const double th ,
	     const double req_msec , 
	     const int req_peaks_flt ,
	     const int req_peaks_raw );
  
  std::vector<ripple_t> ripples;
  
};

#endif

