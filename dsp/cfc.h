
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


#ifndef __CFC_H__
#define __CFC_H__

struct edf_t;
struct param_t;

#include <vector>

namespace dsptools 
{  
  // ultimate wrapper
  void cfc( edf_t & , param_t & );
}


struct cfc_t
{
  
  // load data
  cfc_t( const std::vector<double> & d , 
	 const double a1,
	 const double a2,
	 const double b1,
	 const double b2 , 
	 const double sr , 
	 const double tw = 1, 
	 const double ripple = 0.01 );
  
  bool glm();
  
private:

  std::vector<double> d;  
  double a1, a2, b1, b2;
  double sr;
  double ripple, tw;

  
  // output
public:  
  double r_PAC;
  double c_AMP, z_AMP;
  double r2_TOT;
  
};


#endif
