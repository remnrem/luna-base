

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



// This is a simple wrapper around some filter design libraries:

// 1) fidlib.h   https://uazu.net/fidlib/ (fidlib.h, Jim Peters)
// 2) fir.h      http://www.labbookpages.co.uk/audio/firWindowing.html (A.Greensted)
// 3) filt.cpp   http://www.cardinalpeak.com/blog?p=1841  (Mike Perkins)


#ifndef __FILTER_WRAP_H__
#define __FILTER_WRAP_H__

#include "fidlib.h"
#include "fir.h"
#include "filt.h"

#include "../helper/helper.h"

#include <string>
#include <vector>

struct edf_t;
struct param_t;
struct signal_list_t;


/* void band_pass_filter( edf_t & , param_t & param ); */

/* std::vector<double> band_pass_filter( const std::vector<double> & input ,  */
/* 				      double Fs ,  */
/* 				      int num_taps ,  */
/* 				      double lwr ,  */
/* 				      double upr ); */




struct filter_t 
{
  
filter_t(double Fs, int lib) : Fs(Fs), lib(lib) 
  {

    // library 1
    f1 = NULL;
    lwr_freq = -1;
    upr_freq = -1;
    num_taps = -1;

       
  }



  void design( const std::string & spec )
  {
    //
  }


  //
  // Show coeficients
  //

  std::vector<double> coef() const;
  

  //
  // Apply filter
  //

  std::vector<double> apply( const std::vector<double> & x )
  {

    if ( lib == 1 ) return filter1( x );
    
    std::vector<double> empty;
    return empty;
    
  }
  

  //
  // library specific : lib1
  //

  void design_bandpass1( double l , double u , double t = -1 )
  {
    if ( lib != 1 ) Helper::halt( "requires lib==1" );
    lwr_freq = l;
    upr_freq = u;
    if ( t > 0 ) num_taps = t;

    if ( f1 != NULL ) delete f1;    
    f1 = new Filter( BPF , num_taps , Fs, lwr_freq , upr_freq );
    if( f1->get_error_flag() != 0 ) Helper::halt( "problem initializing type 1 BPF: code " + Helper::int2str( f1->get_error_flag() ) );    
    
  }
  
  void design_lowpass1( double f , int t = -1 )
  {
    lwr_freq = f;
    upr_freq = -1;
    if ( t > 0 ) num_taps = t;

    // set num_taps to a default here

    if ( f1 != NULL ) delete f1;
    f1 = new Filter( LPF , num_taps , Fs , lwr_freq );
    if( f1->get_error_flag() != 0 ) Helper::halt("problem initializing LPF");    
  }
  
  void high_pass( double f , int t = -1 ) 
  {
    lwr_freq = -1;
    upr_freq = f;
    if ( t > 0 ) num_taps = t;
    else { } // set num_taps

    if ( f1 != NULL ) delete f1;
    f1 = new Filter( HPF , num_taps , Fs , upr_freq );
    if( f1->get_error_flag() != 0 ) Helper::halt("problem initializing HPF");    
  }

  // apply to a single sample
  double filter1( double x) { return f1->do_sample(x); }
  
  // forward/reverse 0-phase distortion filter
  std::vector<double> filter1( const std::vector<double> & x );





  //
  // Generic clean-up
  //

  ~filter_t() 
  { 
    // clean up 
    if (f1) delete f1; 
    f1 = NULL; 
  } 
  


private:

  
  // library
  //  1: filt.h
  //  2: fir.h
  //  3: fidlib.h

  int lib;

  // sample rate 

  double Fs;
  
  // transition frequencies

  double lwr_freq;
  double upr_freq;

  // (if specified, number of taps)
  int num_taps;

  
  // Library 1 

  Filter * f1;


  
};

#endif
