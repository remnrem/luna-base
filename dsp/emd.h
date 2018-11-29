
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


#ifndef __EMD_H__
#define __EMD_H__

#include "hilbert.h"
#include <vector>
#include <set>



namespace dsptools { 
  void emd_wrapper( edf_t & edf , param_t & param ) ;
}


void test_emd();


struct extrema_t { 
  extrema_t( const std::vector<double> & d );
  
  std::vector<int> minindex_start, minindex_stop;
  std::vector<int> maxindex_start, maxindex_stop;
  
  // return unique list of points 
  std::vector<int> maxindex();
  std::vector<int> minindex();
  std::vector<int> cross();
  
  std::vector<int> cross_start; // index points of zero-crossings
  std::vector<int> cross_stop; // index points of zero-crossings
  
  int nmax , nmin;
  int nextrema; // number of extrema 
  int ncross;   // and zero-crossings
  bool imf() const { return nextrema == ncross || nextrema == ncross + 1 ; }
};


struct emd_t
{
  
  // Empirical mode decomposition
  emd_t( const std::vector<double> & d , double );
  
  // ensemble EMD
  int n_iter() const { return iter; } 
  void n_iter( const int i ) { iter = i; }
  
  void set_noise_sd( double d ) { noise_sd = d; }
  void set_sd_threshold( double d ) { sd_threshold = d; }
  
  double Fs;
  double tol;
  int stop_mode;
  int max_sift;
  int max_imf;
  
  std::vector<double> sift( const std::vector<double> & ); 
  std::vector<double> envelope_mean( const std::vector<double> & ); 
  
  std::vector<std::vector<double> > imf;
  std::vector<double> residual;
  
  extrema_t extrema( const std::vector<double> & ) ;
  
  int iter;              // number of iterations (>1 implies ensemble EMD)
  
  double sd_threshold;   // determine when to stop sifting
  
  double noise_sd;
  
  
};


#endif
