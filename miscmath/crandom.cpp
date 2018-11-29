
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

#include "miscmath/crandom.h"

#include <cstdlib>
#include <cstdio>
#include <cstddef>
#include <ctime>
#include <vector>

using namespace std;

int         CRandom::iy=0;
vector<int> CRandom::iv;

double CRandom::last = 0;

const int CRandom::IA=16807;
const int CRandom::IM=2147483647;
const int CRandom::IQ=127773;
const int CRandom::IR=2836;
const int CRandom::NTAB=32;
const int CRandom::NDIV=(1+(IM-1)/NTAB);

const double CRandom::EPS=3.0e-16;
const double CRandom::AM=1.0/IM;
const double CRandom::RNMX=(1.0-EPS);

int CRandom::idum=0;


//
// Set seed
//

void CRandom::srand ( long unsigned i )
{

  idum = -i;
    
  CRandom::iv.resize(NTAB);

  if (idum <= 0 || !iy) {
    if (-idum < 1) idum=1;
    else idum = -idum;
    for (int j=NTAB+7;j>=0;j--) {
      int k=idum/IQ;
      idum=IA*(idum-k*IQ)-IR*k;
      if (idum < 0) idum += IM;
      if (j < NTAB) iv[j] = idum;
    }
    iy=iv[0];
  }
    
}

//
// Return the next random number
//

double CRandom::rand ()
{
  int j,k;
  double temp;

  k=idum/IQ;
  idum=IA*(idum-k*IQ)-IR*k;
  if (idum < 0) idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else 
    {
      last = temp;
      return temp;
    }
}


//
// Return a random integer between 0 and fac-1
//

int CRandom::rand (int n)
{
  int r = int(rand() * n);
  if (r == n) r--;
  return r;
}

