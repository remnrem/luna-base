
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

#ifndef __CRANDOM_H__
#define __CRANDOM_H__

#include <vector>

class CRandom 
{
 public:
  
  static const int IA;
  static const int IM;
  static const int IQ;
  static const int IR;
  static const int NTAB;
  static const int NDIV;
  
  static const double EPS;
  static const double AM;
  static const double RNMX;
  
  // Current seed
  static int idum;
  
  static int iy;
  static std::vector<int> iv; 

  static double last;

  static void srand(long unsigned iseed = 0);
  static double rand();
  static int rand (int);
  
};

#endif
