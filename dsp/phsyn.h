
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

#ifndef __PHASE_SYNCHRONY_H__
#define __PHASE_SYNCHRONY_H__

#include <vector>
#include "defs/defs.h"

struct edf_t;
struct param_t;

struct phsyn_t
{
  
  phsyn_t( const std::vector<double> & x , 
	   const int sr , 
	   const std::vector<freq_range_t> & f1 , 
	   const std::vector<freq_range_t> & f2 , 
	   const int nbins = 20 , 
	   const int nreps = 1000 , 
	   const double ripple = 0.01 , 
	   const double tw = 0.5 , 
	   const int es = 0 
	   )
  : x(x) , sr(sr), f1(f2), f2(f2), nbins(nbins) , nreps(nreps) , ripple(ripple) , tw(tw) , es(es) 
  { 
  }
  
  // functions
  void calc();
  
  // inputs
  const std::vector<double> & x;
  const int sr;
  const std::vector<freq_range_t> & f1;
  const std::vector<freq_range_t> & f2;
  const int nbins;
  const int nreps;
  const double ripple;
  const double tw;
  const int es;

  // outputs
  std::vector<std::vector<double> > obs;
  std::vector<std::vector<double> > perm;

  std::vector<std::vector<int> > pv; // empirical p-value
  std::vector<std::vector<double> > z; // z-value (accumulator)
  std::vector<std::vector<double> > z2; // sum of squares


private:

  bool bin( double d , int * b , const std::vector<double> & th , const int nbins );

};



namespace dsptools 
{  
  // wrapper
  void phsyn( edf_t & , param_t & );
}

#endif
