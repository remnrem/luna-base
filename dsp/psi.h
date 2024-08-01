
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

#ifndef __PSI_H__
#define __PSI_H__

struct edf_t;
struct param_t;
struct signal_list_t;

#include "stats/matrix.h"
#include "fftw/fftwrap.h"
#include <complex>
#include <vector>
#include "miscmath/qdynam.h"

namespace dsptools
{
  void psi_wrapper( edf_t & edf , param_t & param );
}

  
struct psi_t {

  psi_t( const Data::Matrix<double> * data , int eplen , int seglen , const int fs )
    : data(data) , eplen( eplen ) , seglen( seglen ) , fs(fs)
  {
    if ( seglen > eplen )
      Helper::halt( "epoch length is smaller than segment length" );

    // populate frequency index given fs and seglen
    // windowing is done manually by PSI routine
    fftseg.init( seglen , seglen , fs , WINDOW_NONE );
    frqs.clear();
    for (int i=0; i<fftseg.cutoff; i++)
      frqs.push_back( fftseg.frq[i] );
  }
  
  int max_freq_idx()
  {
    // max index
    int m = 0;
    for (int i=0;i<freqbins.size();i++)
      for (int j=0;j<freqbins[i].size();j++)
	if ( freqbins[i][j] > m ) m = freqbins[i][j];
    return m;
  }
     

  void add_freqbin()
  {
    // add all except DC term
    std::vector<int> t;
    for (int f=0;f<frqs.size(); f++)
      if ( frqs[f] > 0 ) t.push_back( f ); // store index     
    freqbins.push_back(t);
  }

  void add_freqbin( const double l , const double u )
  {
    std::vector<int> t;
    for (int f=0;f<frqs.size(); f++)
      if ( frqs[f] >= l && frqs[f] <= u )
	t.push_back( f ); // store index   
    freqbins.push_back(t);
  }

  void calc();

  void report( const signal_list_t & ,  bool by_epoch = false ,
	       qdynam_t * qd = NULL ,
	       const int e = -1 );
  
  std::vector<Data::Matrix<std::complex<double> > > data2cs_event( const Data::Matrix<double> * , int maxfreqbin );

  Data::Matrix<double> cs2ps( const std::vector<Data::Matrix<std::complex<double> > > & cs );
  
  void ps();
  
  
  //
  // Inputs
  //
  
  const Data::Matrix<double> * data;

  int eplen;
  int seglen;
  int segshift;
  int fs;
  
  std::vector<double> frqs;
  std::vector<std::vector<int> > freqbins;

  //
  // Internal 
  //

  int nchan; // number of channels

  real_FFT fftseg;
  
  qdynam_t * qd;
  int qe; // this epoch (as psi_t() is called per epoch)

  //
  // Outputs
  //

  int n_models;

  std::vector<Data::Matrix<double> > psi;
  std::vector<Data::Matrix<double> > std_psi;

  std::vector<Data::Vector<double> > psi_sum;
  std::vector<Data::Vector<double> > std_psi_sum;

  std::vector<Data::Vector<double> > apsi_sum;
  std::vector<Data::Vector<double> > std_apsi_sum;

};

#endif
