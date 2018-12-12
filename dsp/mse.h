
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


#ifndef __MSE_H__
#define __MSE_H__


#include <vector>
#include <map>

struct mse_t
{
  
  //pattern length (m) and the similarity criterion (r)

  std::vector<double> coarse_graining( const std::vector<double> & x , int j );

  double SD( const std::vector<double> & x );
  
  double sampen( const std::vector<double> y , int M , double r );

  // MSE
  
  mse_t( const int scale_min = 1 , 
	 const int scale_max = 20 , 
	 const int scale_step = 1 ,
	 const int m = 2 , 
	 const double r = 0.15)
  :     m(m) , 
        r(r) ,
        scale_min(scale_min) , 
        scale_max(scale_max) , 
        scale_step(scale_step)  
    {   }
  
  std::map<int,double> calc( const std::vector<double> & d );
  
private:

  double sample_entropy( const std::vector<double> & y , double sd = 1.0 );
  
  double proc();

  // pattern length for SampEn
  // default 2; range 1..10

  double m;
  
  // similarity criterion for SampEn
  // default 0.15; range >0 

  double r; 
  
  // scale min, max, step

  int scale_min, scale_max, scale_step;
  
};


#endif
