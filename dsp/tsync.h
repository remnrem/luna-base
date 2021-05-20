
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


#ifndef __LUNA_TSYNC_H__
#define __LUNA_TSYNC_H__

#include "stats/Eigen/Dense"
#include <map>

struct edf_t;
struct param_t;

namespace dsptools {
  void tsync( edf_t & , param_t & );
}  

struct tsync_t {

  // xcorr
  tsync_t( const Eigen::MatrixXd & X ,
	   double sr ,
	   int w
	   );

  std::map<int,std::map<int,std::map<int,double> > > xcorr;
  std::map<int,std::map<int,int> > delay;

  // phase-based measures (assuming HT X -> A & P )
  tsync_t( const Eigen::MatrixXd & P ,
	   const Eigen::MatrixXd & A ,
	   double sr ,
	   int w
	   );
  
  double pdiff( double a1 , double a2 )
  {
    double d = a1 - a2;
    if ( fabs(d) > M_PI )
      {
	if ( d < 0 ) d += 2 * M_PI;
	else d -= 2 * M_PI;      
      }
    return d;
  }

  std::map<int,std::map<int,std::map<int,double> > > phase_diff;
  std::map<int,std::map<int,std::map<int,double> > > phase_lock;
  
    
  
};


#endif
