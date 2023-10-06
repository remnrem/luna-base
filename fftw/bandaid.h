
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

#ifndef __BANDAID_H__
#define __BANDAID_H__

#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <complex>

#include "helper/helper.h"
#include "helper/logger.h"
#include "miscmath/miscmath.h"
#include "defs/defs.h"

extern logger_t logger;

struct param_t;

struct bandaid_t {

  bandaid_t();

  void init();

  int size() const;
  
  void define_bands( param_t & );

  // used w/ PSD
  void track_bands_per_epoch( double , double, double, double, double, double, double, double, double, double ); 
  double psdsum( const std::vector<double> & f , const std::vector<double> & x , const freq_range_t & b );
  
  // these next two used w/ MTM 
  void calc_bandpower( const std::vector<double> & f , const std::vector<double> & x );

  void track();
  
  double fetch( frequency_band_t b ) const;
  
  void freq_band_settings( const std::string & b , double * r0 , double * r1 );

  std::map<frequency_band_t,std::vector<double> > track_band;
  
  std::vector<frequency_band_t> bands;
  
  // holder of 'current' values
  double slow, delta, theta, alpha, sigma, beta, gamma;
  double low_sigma, high_sigma;
  double denom, total;
  
};


#endif
