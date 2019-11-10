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

#ifndef __ICA_H__
#define __ICA_H__

#include "ica/lica_matrix.h"
#include "ica/lica_fastICA.h"

#include <vector>
#include "stats/matrix.h"
#include "helper/helper.h"


// Notes
// https://en.wikipedia.org/wiki/FastICA
// http://research.ics.aalto.fi/ica/book/

// Code here modified from C/C++ implementation of R fastICA package
// http://tumic.wz.cz/fel/online/libICA/

struct edf_t; 
struct param_t; 

void ica_wrapper0( edf_t & , param_t & ); 

struct ica_t {
  
  ica_t( Data::Matrix<double> & X , int compc )
  {        
    if ( ! proc(X,compc) ) Helper::halt( "problem in ica_t" );    
  }

  Data::Matrix<double> K;
  Data::Matrix<double> W;
  Data::Matrix<double> A;
  Data::Matrix<double> S;
  
  bool proc( Data::Matrix<double> & , int compc );


  //
  // Helper functions
  //
  
  static void cpp_svdcmp( Data::Matrix<double> & A, Data::Vector<double> & W, Data::Matrix<double> & V);
   
};


#endif
