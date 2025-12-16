
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

#ifndef __WRAPPERS_H__
#define __WRAPPERS_H__

struct edf_t;
struct param_t;

#include "dsp/fir.h"

#include <vector>
#include <stdint.h>
#include <cstddef>

namespace dsptools 
{  

  void make_bands( edf_t & , param_t & );
  
  void cwt( edf_t & , param_t & );

  void hilbert( edf_t & , param_t & );

  void qdynam( edf_t & , param_t & );
   
  void run_cwt( const std::vector<double> & data , const int ,
		const double fc , const int num_cycles , 
		std::vector<double> * mag , 
		std::vector<double> * phase = NULL ); 

  void alt_run_cwt( const std::vector<double> & data ,
		    const int Fs, 
		    const double fc ,
		    const double FWHM ,
		    const double tlen , 
		    const bool wrapped , 
		    std::vector<double> * mag , 
		    std::vector<double> * phase );

  
  // Hilbert
  void run_hilbert( const std::vector<double> & data , const int Fs , 		    
		    std::vector<double> * mag , 
		    std::vector<double> * phase = NULL ,
		    std::vector<double> * angle = NULL , 
		    std::vector<double> * frequency = NULL ); 

  // filter-Hilbert, Kaiser
  void run_hilbert( const std::vector<double> & data , const int Fs , 
		    const double flwr , const double fupr , const double ripple , const double tw , 
		    std::vector<double> * mag , 
		    std::vector<double> * phase = NULL ,
		    std::vector<double> * angle = NULL , 
		    std::vector<double> * frequency = NULL ); 

  // filter-Hilbert, from file
  void run_hilbert( const std::vector<double> & data , const int Fs , 
		    const std::string & fir_file , 
		    std::vector<double> * mag , 
		    std::vector<double> * phase = NULL ,
		    std::vector<double> * angle = NULL , 
		    std::vector<double> * frequency = NULL ); 

  // filter-Hilbert, fixed order
  void run_hilbert( const std::vector<double> & data , const int Fs , 
		    const double flwr , const double fupr , const int order , const fir_t::windowType , 
		    std::vector<double> * mag , 
		    std::vector<double> * phase = NULL ,
		    std::vector<double> * angle = NULL , 
		    std::vector<double> * frequency = NULL ); 


  // FFT

  void fft( edf_t & , param_t & );

  void cmdline_fft( param_t & );
  
  void run_fft( const std::vector<double> & x , const int Fs , const bool verbose );
  
  // Otsu 
  
  void otsu( edf_t & , param_t & );
  
  void cmdline_otsu( param_t & );
  
  void run_otsu( const std::vector<double> & x , const int k );

    
  // Helper
  
  std::vector<double> readcin();
  
}

#endif
