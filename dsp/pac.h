
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


#ifndef __PAC_H__
#define __PAC_H__

struct edf_t;
struct param_t;

#include <vector>

struct pac_t
{

  // phase-amplitude coupling
  
  pac_t( const std::vector<double> * d , 
	 double a , double b ,
	 const int sr ,
	 const int nr = 1000 )
  {
    frq4phase.clear();
    frq4pow.clear();
    frq4phase.push_back(a);
    frq4pow.push_back(b);
    data = d;
    na = nb = 1;
    srate = sr;
    nreps = nr;
    size();
  }
  
  pac_t( const std::vector<double> * d ,
	 const std::vector<double> & a ,
	 const std::vector<double> & b ,
	 const int sr ,
	 const int nr = 1000 )
  {
    frq4phase = a;
    frq4pow = b;
    na = a.size();
    nb = b.size();
    data = d;
    srate = sr;
    nreps = nr;
    size();
  }

  void size()
  {
    z.resize( na );
    for (int a=0;a<na;a++) z[a].resize(nb,0);
    pval.resize( na );
    for (int a=0;a<na;a++) pval[a].resize(nb,0);    
  }
  
  void init()
  {
    na = nb = 0;
    z.clear();
    pval.clear();
  }

  double zpac( int a = 0 , int b = 0 ) const 
  {
    if ( a < 0 | a > na ) return -9;
    if ( b < 0 | b > nb ) return -9;
    return z[a][b];
  }
  
  double ppac( int a = 0 , int b = 0 ) const 
  {
    if ( a < 0 | a > na ) return -9;
    if ( b < 0 | b > nb ) return -9;
    return pval[a][b];
  }
  
  std::vector<std::vector<double> > zpac_all() const
  {
    return z;
  }

  std::vector<std::vector<double> > ppac_all() const
  {
    return pval;
  }
  
  bool calc();

  //
  // data
  //

  const std::vector<double> * data;
  std::vector<double> frq4phase;
  std::vector<double> frq4pow;
  std::vector<std::vector<double> > z;
  std::vector<std::vector<double> > pval;
  int srate;
  int na,nb;
  int nreps;
};



namespace dsptools 
{  
  // ultimate wrapper
  void pac( edf_t & , param_t & );
}




#endif
