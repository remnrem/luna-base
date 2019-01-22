
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

#ifndef __POLARITY_H__
#define __POLARITY_H__

#include "edf/edf.h"
#include "intervals/intervals.h"
#include <vector>


namespace dsptools 
{ 

  void polarity( edf_t & edf , const param_t & param );
  
  std::vector<bool> make_mask( const std::vector<double> & x , double th );
  
  void polarity_check( const std::vector<double> & x0 , const std::vector<uint64_t> * tp , int fs , 
		       double th ,  // threshold for extracted filtered segments
		       bool zc2zc , // extract around peaks up to ZC's
		       double flim , // calculate PSD up to flim Hz only
		       double f_lwr , // lower BPF transition frequency
		       double f_upr , // upper 
		       bool mirror_mode ,   // mirror odd numbered up/down segments, i.e. make 'wave-like' signal
		       bool double_up , // instead of mirror-mode alternate segments, double-enter each (up and down)
		       bool analyse_bpf_signal , // analysis the BPF signal, not raw data
		       bool dmode // segment 'upward' and 'downward' rather than pos and neg
		       );


  void ht_polarity_check( const std::vector<double> & x0 , 
			  const std::vector<uint64_t> * tp , 
			  int fs , 
			  double f_lwr , // lower BPF transition frequency
			  double f_upr // upper 
			  );
    
}



#endif
