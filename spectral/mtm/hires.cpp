
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

// All mtm functions are adapted code from: Lees, J. M. and J. Park
// (1995): Multiple-taper spectral analysis: A stand-alone
// C-subroutine: Computers & Geology: 21, 199-236.

#include "mtm.h"

#include <cstdio>
#include <cmath>


int mtm_t::hires(double *sqr_spec,  double *el, int nwin, int num_freq, double *ares)
{

  for (int j = 0; j < num_freq; j++)
    ares[j] = 0.;
  
  for (int i = 0; i < nwin; i++)
    {
      const int k = i * num_freq;
      double a = 1. / (el[i] * nwin);
      for ( int j = 0; j < num_freq; j++)
	{
	  const int kpoint = j + k;
	  ares[j] += a * ( sqr_spec[kpoint] );
	}
    }
  
  for (int j = 0; j < num_freq; j++)
    {
      if( ares[j]>0.0 ) 
	ares[j] = sqrt(ares[j]);
      else
	printf("sqrt problem in hires pos=%d %f\n", j, ares[j]);
    }
  
  return 1;
}      



