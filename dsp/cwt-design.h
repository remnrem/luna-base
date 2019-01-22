
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

#ifndef __CWT_DESIGN_H__
#define __CWT_DESIGN_H__

#include <vector>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <complex>

#include "fftw3.h"
#include "helper/helper.h"
#include "eval.h"

namespace dsptools 
{ 
  void design_cwt( param_t & param );  
} 


#endif
