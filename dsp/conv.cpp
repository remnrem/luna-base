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


#include "conv.h"

#include <cstddef>
#include <cstdio>


std::vector<double> dsptools::convolve( const std::vector<double> & signal , 
					const std::vector<double> & kernel )
{

  const int nsig = signal.size();
  const int nkern = kernel.size();
  const int nconv = nsig + nkern - 1 ; 

  std::vector<double> result( nconv , 0 );

  for (int n=0; n < nconv; n++ )
    {
      
      size_t kmin, kmax;

      kmin = (n >= nkern - 1) ? n - (nkern - 1) :  0;
      kmax = (n < nsig - 1)   ? n               :  nsig - 1;

      for (size_t k = kmin; k <= kmax; k++)
	result[n] += signal[k] * kernel[n - k];       
    }
  
  return result;

}

