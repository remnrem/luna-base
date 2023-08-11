
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

#ifndef __PERI_H__
#define __PERI_H__

struct edf_t;
struct param_t;
struct signal_list_t;

#include <vector>
#include <stats/Eigen/Dense>

namespace dsptools
{
  void peri( edf_t & edf , param_t & param );
}

struct peri_param_t
{
  peri_param_t( param_t & param );

  // time window: centering, censoring 
  double time0;
  double time_left, time_right;

  // CWT
  bool cwt_do;
  std::vector<double> cwt_f;
  std::vector<double> cwt_fwhm;
  double cwt_timelength;

  // XCORR
  bool xcorr_do;
  double xcorr_w_sec;
  double xcorr_c_sec;

  
  // GP
  
  
};

struct peri_t {

  // epoch-by-epoch data
  // epoch(ne) x sample-point(fs*epoch-dur) x channel(ns) 
  //   each epoch must have a similar size - or else be zero-padded 
  
  peri_t( const std::vector<Eigen::MatrixXd> & X ,
	  const peri_param_t & pp , 
	  const signal_list_t & , 
	  const int fs );
    
};


#endif
