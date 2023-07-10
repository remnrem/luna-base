
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

#ifndef __LUNA_DSP_IIR_H__
#define __LUNA_DSP_IIR_H__

#include <vector>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <complex>

#include "dsp/filter.h"

struct edf_t;
struct param_t;

namespace dsptools { 
  
  void apply_iir( edf_t & edf , param_t & param );
}


enum iir_type_t {
  BUTTERWORTH_LOWPASS ,
  BUTTERWORTH_HIGHPASS ,
  BUTTERWORTH_BANDPASS ,
  BUTTERWORTH_BANDSTOP , 
  CHEBYSHEV_LOWPASS ,
  CHEBYSHEV_HIGHPASS ,
  CHEBYSHEV_BANDPASS ,
  CHEBYSHEV_BANDSTOP
};

struct iir_t {

  iir_t();
  
  void init( iir_type_t , int order , double p1, double p2 , double p3 = 0 , double p4 = 0 );
  
  ~iir_t();
  
  std::vector<double> apply( const std::vector<double> & x );

private:

  BWLowPass  * bwlp;
  BWHighPass * bwhp;
  BWBandPass * bwbp;
  BWBandStop * bwbs;

  CHELowPass  * chelp;
  CHEHighPass * chehp;
  CHEBandPass * chebp;
  CHEBandStop * chebs;

  
};


#endif

