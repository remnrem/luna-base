
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

#ifndef __MI_H__
#define __MI_H__

struct edf_t;
struct param_t;

#include <vector>

struct mi_t
{

  mi_t() { } 

  mi_t( const std::vector<double> & a , const std::vector<double> & b );
  
  // use previous calculated bins (i.e. done on whole datasets
  // so that same grid can be used per epoch
  
  void force_thresholds( const std::vector<double> & _tha, 
			 const std::vector<double> & _thb );


  // bin rule : 1 = FD, 2 = Scott, 3 = Sturges

  void set_nbins( const int b );
  int set_nbins_sturges();
  int set_nbins_scott();
  int set_nbins_fd();

  // given bins, set thresholds
  int set_thresholds( const int bin_rule = 1 );

  // given thresholds, bin data
  void bin_data();

  // given binned data, calc MI
  void calc();

  // having calc'ed MI, get null distribution
  void  permute( const int rep , double * , double * );  

  
  // univariate information/entropy
  double infa; double pvala;
  double infb; double pvalb;
  
  // joint/mutual information
  double jointinf; double pvaljoint;
  double mutinf;   double pvalmut;
  
  double total_corr; double dual_total_corr;

  int n;

  // thresholds (depends on data and nbins)
  int nbins;
  std::vector<double> tha;
  std::vector<double> thb;
  
private:

  double eps;
    
  // copy of data
  std::vector<double> da;
  std::vector<double> db;

  // binned data
  std::vector<int> bina;
  std::vector<int> binb;
  
};


namespace dsptools 
{  
  // wrapper
  void compute_mi( edf_t & , param_t & );  
}

#endif
