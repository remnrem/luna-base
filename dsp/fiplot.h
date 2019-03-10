
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

#ifndef __FIPLOT_H__
#define __FIPLOT_H__

#include <vector>
#include <map>
#include <stdint.h>

struct edf_t;
struct param_t;

void fiplot_wrapper( edf_t & edf , const param_t & param );

struct fipair_t 
{ 
  fipair_t(double w, double n) : w(w), n(n) { } 
  fipair_t()  { w=0; n=0; } 
  double w; double n; 
};

struct fibin_t { std::map<double,fipair_t> r; };

struct fipoint_t 
{ 
  fipoint_t( int _i, int _j , double _h )
  {
    i = _i; j = _j; h = _h; t = j - i + 1 ;
  }

  int i; 
  int j;
  int t;
  double h; 

  bool operator< ( const fipoint_t & rhs ) const 
  {
    // sort by duration, longest first
    if ( t > rhs.t ) return true;
    if ( t < rhs.t ) return false;
    return i < rhs.i;
  }
}; 

struct fiplot_t
{
  fiplot_t( const std::vector<double> & x , const std::vector<uint64_t> * tp , const int _fs  ,
	    double t_lwr , double t_upr , double t_inc , bool cycles , 
	    double f_lwr , double f_upr , double f_inc , bool logspace = true )
  {
    fs = _fs;
    set_t( t_lwr, t_upr, t_inc , cycles );
    set_f( f_lwr, f_upr, f_inc , logspace );
    proc(x,tp,fs);
  }
  
  void set_t( double lwr , double upr , double inc , bool cyc ) 
  {
    t_lwr = lwr;
    t_upr = upr;
    t_inc = inc;
    cycles = cyc;
  }

  void set_f( double lwr , double upr , double inc , bool logspace ) ;

  void proc( const std::vector<double> & x , 
	     const std::vector<uint64_t> * tp , 
	     const int fs ); 
    
  int nt; // number of time points
  int nf; // number of frequencies
  
  int fs; // sample rate

  std::vector<double> frqs;
  double f_lwr, f_upr, f_inc;
  
  double t_lwr, t_upr, t_inc;
  bool cycles;

  std::vector<double> cwt( const std::vector<double> & x , const int fs, const double fc , const int num_cycles );
  
  fibin_t intervalize( const std::vector<double> & x , 
		       const std::vector<uint64_t> * tp , 
		       const int fs , 
		       const double t_lwr , 
		       const double t_upr , 
		       const double t_inc , 
		       const bool cycles , 
		       const double fc );
  
  
  
};


#endif
