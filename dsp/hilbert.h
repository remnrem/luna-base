
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


#ifndef __HILBERT_H__
#define __HILBERT_H__

struct edf_t;
struct param_t;

#include <vector>

struct hilbert_t
{
  
  hilbert_t() { } 
  
  // Hilbert transform
  hilbert_t( const std::vector<double> & d );
  
  // filter-Hilbert
  hilbert_t( const std::vector<double> & d , const int sr , double lwr , double upr , double ripple, double tw );
  
  
  // extract instantaneous phase, magnitude
  const std::vector<double> * phase() const;
  //  const std::vector<double> * angle() const;
  const std::vector<double> * magnitude() const;
  const std::vector<double> * signal() const;
  std::vector<double> instantaneous_frequency(double) const;

  double phase_events( const std::vector<int> & , std::vector<double> * , double *, int *, double *) const;
  double phase_events( const std::vector<int> & , std::vector<double> * , double , double *, int *, double *, std::vector<bool> * ) const;
  
private:

  void proc();
  void unwrap(std::vector<double> * ) const;

  std::vector<double> input;
  std::vector<double> ph;
  std::vector<double> mag;
  
};


#endif
