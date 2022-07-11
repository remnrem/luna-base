
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


#ifndef __MTSPEC_H__
#define __MTSPEC_H__

#include "stats/matrix.h"

struct mt_spectrogram_t
{

  mt_spectrogram_t( const Data::Matrix<double> & X ,
		    const int sr ,
		    const double nw ,
		    const int t ,
		    const double seg_sec ,
		    const double seg_inc ,
		    const double min_f ,
		    const double max_f , 
		    const bool dB ,
		    const bool mean_center ); 
  
  // input X = rows = observations/epochs 
  //           cols = samples
  
  
  // output Z = rows = frequencies
  //            cols = samples (time w/ within segment, defined by seg_inc)
  
  Data::Matrix<double> Z, Z_median, ZZ;
  std::vector<double> frq; // rows (freq)
  std::vector<double> t;   // cols (time)
    
    
};



#endif
