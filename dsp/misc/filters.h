
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


#ifndef __FILTERS_H__
#define __FILTERS_H__

#include "filt.h"

#include "../helper/helper.h"

struct edf_t;
struct param_t;
struct signal_list_t;

void band_pass_filter( edf_t & , param_t & param );

std::vector<double> band_pass_filter( const std::vector<double> & input , 
				      double Fs , 
				      int num_taps , 
				      double lwr , 
				      double upr );

class MyFilter 
{
  
 public:

  MyFilter() { f = NULL; }
  
  void low_pass( int num_taps , double sampling_frequency , double transition_frequency )
  {
    f = new Filter( LPF , num_taps , sampling_frequency , transition_frequency );
    if( f->get_error_flag() != 0 ) Helper::halt("problem initializing LPF");    
  }
  
  void high_pass( int num_taps , double sampling_frequency , double transition_frequency )
  {
    f = new Filter( HPF , num_taps , sampling_frequency , transition_frequency );
    if( f->get_error_flag() != 0 ) Helper::halt("problem initializing HPF");    
  }
  
  void band_pass( int num_taps , double sampling_frequency , double lower_frequency , double upper_frequency )
  {
    f = new Filter( BPF , num_taps , sampling_frequency , lower_frequency , upper_frequency );
    if( f->get_error_flag() != 0 ) Helper::halt("problem initializing BPF: code " + Helper::int2str( f->get_error_flag() ) );    
  }

  double filter(double x) { return f->do_sample(x); }

  std::vector<double> coefs() { return f->get_taps(); } 

  // forward/reverse 0-phase distortion filter
  std::vector<double> filter( const std::vector<double> & x )
    {
      
      // forward
      int i = 0;
      std::vector<double> y( x.size() );
      std::vector<double>::const_iterator fi = x.begin();
      while ( fi != x.end() ) { y[i++] = f->do_sample( *fi ); ++fi; }

      // reverse
      i = 0;
      std::vector<double>::reverse_iterator ri = y.rbegin();
      while ( ri != y.rend() ) { *ri = f->do_sample( *ri ); ++ri; }
      
      return y;
    }
  
  ~MyFilter() { if (f) delete f; f = NULL; } 
  
 private:
  
  Filter * f;
  
};

#endif
