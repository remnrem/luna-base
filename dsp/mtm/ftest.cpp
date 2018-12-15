
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
#include <cstdlib>
#include <cmath>

#include "nrutil.h"


void mtm::get_F_values(double *sr, double *si, int nf, int nwin, double *Fvalue, double *b)
{

  // 
  // b is fft of slepian eigentapers at zero freq sr si are the
  // eigenspectra amu contains line frequency estimates and f-test
  // parameter
 
  double          sum, sumr, sumi, sum2;
  int             i, j, k;
  double         *amur, *amui;
  sum = 0.;
  amur = dvector((long) 0, (long) nf);
  amui = dvector((long) 0, (long) nf);
  
  for (i = 0; i < nwin; i++)
    sum = sum + b[i] * b[i];
  
  for (i = 0; i < nf; i++) {
    amur[i] = 0.;
    amui[i] = 0.;
    for (j = 0; j < nwin; j++) {
      k = i + j * nf;
      amur[i] = amur[i] + sr[k] * b[j];
      amui[i] = amui[i] + si[k] * b[j];
    }
    amur[i] = amur[i] / sum;
    amui[i] = amui[i] / sum;
    sum2 = 0.;
    for (j = 0; j < nwin; j++) {
      k = i + j * nf;
      sumr = sr[k] - amur[i] * b[j];
      sumi = si[k] - amui[i] * b[j];
      sum2 = sum2 + sumr * sumr + sumi * sumi;
    }
    Fvalue[i] = (double) (nwin - 1) * (SQR(amui[i]) + SQR(amur[i])) * sum / sum2;
    /* percival and walden, eq 499c, p499 */
    /* sum = Hk(0) squared  */
  }
  return;
}
