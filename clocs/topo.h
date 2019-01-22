
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


#ifndef __TOPO_H__
#define __TOPO_H__

#include <string>
#include <map>
#include <vector>

#include "stats/matrix.h"
#include "edf/edf.h"

struct chid_t { 
  chid_t(int n) : n(n) , label("") { } 
  chid_t() { } 
  int n;
  std::string label;

  bool operator<( const chid_t & rhs ) const { return n < rhs.n; } 
};

struct topoloc_t { 
  topoloc_t() { 
    th = r = 0;
    x = y = 0;
  }

  // assume entry in degrees, then radius
  topoloc_t( const double _th , const double _r )
  {
    // degress->radians
    th = M_PI/180.0 * _th;
    r = _r;
    // cartesian co-ordinates
    x = r * cos( th );
    y = r * sin( th );
  }

  double th,r,x,y;    
};

struct topo_t {
  
  std::map<chid_t,topoloc_t> cxy;
  
  std::map<std::string,int> lab2n;
  
  int size() const { return lab2n.size(); } 
  
  std::set<std::string> channels() const { 
    std::set<std::string> c;
    std::map<std::string,int>::const_iterator ii = lab2n.begin();
    while ( ii != lab2n.end() ) { c.insert( ii->first ); ++ii; } 
    return c; }

  int label2n( const std::string & s );

  int load( const std::string & );
  
  bool add( const std::string & , const topoloc_t & );
  
  void grid( double, double, int, double, double, int);

  void grid( int, int );
  
  void max_radius( double f );

  void squeeze( double f );
  
  void pos();

  void dump();

  bool scaled_xy( const std::string & ch , double * , double * ); 
  
  Data::Matrix<double> interpolate( const std::map<std::string,double> & );
  
  // input co-ords
  int inp_n;
  std::vector<double> inp_xy; // 2n
  std::vector<bool> has_ch;

  int nx,ny;
  int out_n;
  std::vector<double> out_xy;  // 2n
  std::vector<double> out_z;   // n
  std::vector<bool>   out_inc; // nx*ny; t if point is in 'n' selected to interpolate

  
};

#endif
