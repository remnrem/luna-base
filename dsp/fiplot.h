
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
#include <cstddef>
#include <map>
#include <stdint.h>

struct edf_t;
struct param_t;

void fiplot_wrapper( edf_t & edf , const param_t & param , const std::vector<double> * raw = NULL , const int * sr = NULL );

struct fipair_t 
{ 
  fipair_t(double w, double n) : w(w), n(n) { } 
  fipair_t()  { w=0; n=0; } 
  double w; double n; 
};

struct fibin_t 
{ 

  std::vector<double> t;
  
  // map keyed of time (double) but rather than use a double as a key, 
  // store the actual time values in 't' and use the integer index for 
  // map access
  
  std::map<int,fipair_t> r; 
  
  // i.e. t[i] --> time
  //      r[i] --> frequency 
  
};

struct fipoint_t 
{ 
  fipoint_t( int _i, int _j , double _h , bool _trunc = false )
  {
    i = _i; j = _j; h = _h; t = j - i + 1 ;
    trunc = _trunc;
  }

  int i; 
  int j;
  int t;
  double h; 
  bool trunc; // if truncated, will not add to final stats
  
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
	    const int _th , const bool norm_ , const bool logit_ ,  
	    double t_lwr , double t_upr , double t_inc , bool cycles , 
	    double f_lwr , double f_upr , double f_inc , int num_cyc , bool logspace = true ,
	    bool _verb = false )
  {
    fs = _fs;
    th = _th;
    normalize = norm_;
    logit = logit_;
    set_t( t_lwr, t_upr, t_inc , cycles );
    set_f( f_lwr, f_upr, f_inc , logspace , num_cyc );
    verbose = _verb;
    proc(x,tp,fs);
  }
  
  void set_t( double lwr , double upr , double inc , bool cyc )
  {
    t_lwr = lwr;
    t_upr = upr;
    t_inc = inc;
    cycles = cyc;
  }

  void set_f( double lwr , double upr , double inc , bool logspace , int num_cycles ) ;

  void proc( const std::vector<double> & x , 
	     const std::vector<uint64_t> * tp , 
	     const int fs ); 
    
  int nt; // number of time points
  int nf; // number of frequencies
  
  int fs; // sample rate
  double th; // multiplicative threshold
  bool normalize; // 0..1 normalization?
  bool logit;  // use log-space values

  bool verbose;
  
  std::vector<double> frqs;
  double f_lwr, f_upr, f_inc;
  
  double t_lwr, t_upr, t_inc;

  // whether to show results by cycles (rather than time)
  bool cycles;

  // CWT number of cycles
  int num_cycles;
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
