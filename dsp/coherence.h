
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

#ifndef __COH_H__
#define __COH_H__

#include <vector>

struct edf_t;
struct param_t;
struct interval_t;

struct coh_t
{

  coh_t() { } 

  coh_t( const int n )
  {
    resize(n);
  }
  
  void resize( const int n )
  {
    frq.resize( n );

    coh.resize( n );
    cross_spectrum.resize( n );

    auto_spectrum1.resize( n );
    auto_spectrum2.resize( n );
    
    cross_norm1.resize( n );
    cross_norm2.resize( n );
        
  }

  
  // results of a coherence analysis
  
  std::vector<double> frq;
  
  std::vector<double> auto_spectrum1;
  std::vector<double> auto_spectrum2;
  std::vector<double> cross_spectrum;

  std::vector<double> cross_norm1;
  std::vector<double> cross_norm2;

  
  std::vector<double> coh;
  
};


namespace dsptools 
{  

  // ultimate wrapper
  void coherence( edf_t & , param_t & , bool legacy = false );

  // for a given pair of signals, return coh object, either using 
  //  coherence_t()  (from fttwrap/) or legacy_coherence (below)

  coh_t coherence( edf_t & , int signal1 , int signal2 , const interval_t & , bool legacy = false );

  // legacy code... remove in fture 
  coh_t legacy_coherence( const std::vector<double> * s1 , 
			  const std::vector<double> * s2 , 
			  double sampfreq );
    
  
  void coh_lremv( std::vector<double> & , const int );

  void coh_r2tx(int nthpo, double* cr0, double* cr1, double * ci0, double * ci1);

  void coh_r4tx(int nthpo, 
		double * cr0, 
		double * cr1, 
		double * cr2, 
		double * cr3, 
		double * ci0, 
		double * ci1, 
		double * ci2, 
		double * ci3);
  
  void coh_r8tx(int nx, int nthpo, int length, 
		double * cr0, double * cr1, double * cr2, double * cr3, double * cr4, double * cr5, 
		double * cr6, double * cr7, double * ci0, double * ci1,
		double * ci2, double * ci3, double * ci4, double * ci5, double * ci6, double * ci7);

  void coh_fft842(int in, int n, double * , double * );

  // end of legacy code...

}

#endif
