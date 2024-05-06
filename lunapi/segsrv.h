
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

#ifndef __SEGSRV_H__
#define __SEGSRV_H__

#include "luna.h"
#include "stats/Eigen/Dense"

struct lunapi_inst_t;
typedef std::shared_ptr<lunapi_inst_t> lunapi_inst_ptr;

struct segsrv_t {
  
public:

  // set up
  segsrv_t( lunapi_inst_ptr ); 

  // initialize (or reset) to current internal state
  int populate( const std::vector<std::string> & chs ,
		const std::vector<std::string> & anns );
  
  // set window
  bool set_window( double a , double b );

  // get signals
  Eigen::VectorXf get_signal( const std::string & ch ) const;

  // and times
  Eigen::VectorXf get_timetrack( const std::string & ch ) const;
  
  // also, a mask of gapped regions in a given time-window   
  // e.g. if 0 .. 300 is total window and gaps between 30-60 and 270 - 350
  //  --> (30,60) , (270,300) 
  std::set<std::pair<double,double> > get_gaps() const;
  
private:

  void init();
  
  // person
  lunapi_inst_ptr p;

  // segments, gaps
  std::set<interval_t> segments, gaps;

  // to look up indices given two timepoints
  // sr: time->idx
  std::map<int,std::map<double,int> > tidx;

  // current time windows (elapsed seconds)
  double awin, bwin;
  
  // maximum (last time-point in seconds)
  double smax;

  // current window in idx points (for a given SR)
  std::map<int,int> aidx, bidx; 

  // add actual data
  bool add_channel( const std::string & );

  // given two times and a sample rate, get indices
  bool get_tidx( double a, double b , int sr , int * aidx, int *bidx ) const;
  
  // signal data
  std::map<std::string,int> srmap;
  std::map<std::string,Eigen::VectorXf> sigmap;
  std::map<int,Eigen::VectorXf> tmap;

};


#endif
