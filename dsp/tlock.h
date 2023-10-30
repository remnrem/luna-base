
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

#ifndef __TLOCK_H__
#define __TLOCK_H__

struct edf_t;
struct param_t;
struct signal_list_t;


#include <vector>
#include "stats/matrix.h"
#include "timeline/timeline.h"

namespace dsptools
{
  void tlock( edf_t & edf , param_t & param );
}


// this is called per signal

struct tlock_t {
  
  tlock_t( edf_t & edf , int );
  
  // build via epochs
  void epoch_builder( const int slot );
  
  // build via cache sample-points
  void cache_builder( cache_t<int> * cache ,
		      const double hw_sec ,
		      const int slot , 
		      const std::string siglabel , 
		      const bool same_channel = false , 
		      const std::string channel_postfix = "" );
  
  
  // normalizations
  
  // normalization only interval edges? o(e.g. 0.1 is 10%+10% of edges
  void edge_normalization( Data::Vector<double> * x , const int  p ) const;
  

  // outputs  
  bool spectrogram( ) const;

  void clearX();
  
  int set_window(int hw);
  void set_window_epoch(int sp ); // know sr

  edf_t & edf;
  
  // time axis
  std::vector<double> t;  
  int sr;

  // for regular means: track the whole matrix (just in case we want to ) 
  //   but only when 'verbose' option is true;
  //   each row is a sample-point in the interval window
  //   each column is a new epoch
  
  // for angles: the matrix will be 360 / angle_bin wide
  //   each row is a sample-point in the interval window
  //   each column is an angle bin: do not track individual epochs
  //   but rather do the summatation in-place

  int ni;
  int np;
  Data::Matrix<double> X;
  
  
  //
  // main input/builder
  //
  
  void add( const std::vector<double> * x , 
	    const int , const int );

  //
  // main outputs
  //

  void outputs();
  
  Data::Vector<double> average( const double th , const double winsor ) const ;
  
  Data::Vector<double> median( const double th , const double winsor ) const ;
  
  Data::Matrix<double> remove_outliers( const Data::Matrix<double> & Y , const double th , const double winsor ) const;
  
  Data::Matrix<double> angles() const ;
  

  // normalisation points
  double norm_pct;
  bool zero_trace;
  double outlier_th;
  double outlier_winsor;
  bool take_log;
  int angle_bins;
  double emid;
  
  // means only? i.e. no need to build the matrix [ todo ] 
  bool means_only;
  Data::Vector<double> means;
  int count;

  // verbose mode
  bool verbose;
    
};


#endif
