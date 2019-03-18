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


extern "C" {
#include "ica_matrix.h"
#include "svdcmp.h"
#include "libICA.h"
}

#include <vector>

#include "helper/helper.h"

// Notes
// https://en.wikipedia.org/wiki/FastICA
// http://research.ics.aalto.fi/ica/book/

// Code adapted from:
// http://tumic.wz.cz/fel/online/libICA/
// wrapper around fastICA() 


struct edf_t; 
struct param_t; 

void ica_wrapper0( edf_t & , param_t & ); 

struct ica_t {
  
  ica_t( const std::vector<std::vector<double> > & X , int compc )
  {        
    if ( ! proc(X,compc) ) Helper::halt( "problem in ica_t" );    
  }
  
  ica_t( mat X , int rows , int cols, int compc )
  {
    if ( ! proc(X,rows, cols,compc) ) Helper::halt( "problem in ica_t" );    
  }

  std::vector<std::vector<double> > K;
  std::vector<std::vector<double> > W;
  std::vector<std::vector<double> > A;
  std::vector<std::vector<double> > S;
  
  bool proc( const std::vector<std::vector<double> > & , int compc );
  bool proc( mat , int rows, int cols, int compc );

  void copy( std::vector<std::vector<double> > & M , double ** O , int r , int c )
  {    
    M.resize( r );
    for (int i=0;i<r;i++) M[i].resize(c);    
    for (int i=0;i<r;i++)
      for (int j=0;j<c;j++)
	M[i][j] = O[i][j];    
  }
  
};


#endif
