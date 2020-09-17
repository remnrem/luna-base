
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

#include "stats/matrix.h"
#include <complex>
#include <vector>

struct psi_t {
  
 psi_t( const Data::Matrix<double> * data , int seglen , int eplen )
  : data(data) , seglen( seglen ) , eplen( eplen )
  {
    calc();
  }
  
  int max_freq()
  {
    int m = 0;
    for (int i=0;i<freqbins.size();i++)
      for (int j=0;j<freqbins[i].size();j++)
	if ( freqbins[i][j] > m ) m = freqbins[i][j];
    return m;
  }
     
  void add_freqbin( const int l , const int u )
  {
    std::vector<int> t;
    for (int f=l;f<=u;f++) t.push_back(f);
    freqbins.push_back(t);
  }

  void calc();
  
  std::vector<Data::Matrix<std::complex<double> > > data2cs_event( const Data::Matrix<double> * , int maxfreqbin );

  Data::Matrix<double> cs2ps( const std::vector<Data::Matrix<std::complex<double> > > & cs );
  
  void ps();
  
  
  //
  // Inputs
  //
  
  const Data::Matrix<double> * data;

  int seglen;
  int segshift;
  int eplen;

  // nb. integer frequencies used
  std::vector<std::vector<int> > freqbins;


  //
  // Internal 
  //

  int nchan; // number of channels
  
  //
  // Outputs
  //

  Data::Matrix<double> psi;
  Data::Matrix<double> stdpsi;

  Data::Vector<double> psisum;
  Data::Vector<double> stdpsisum;

  
};

#endif
