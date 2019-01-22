
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

#ifndef __SLICE_H__
#define __SLICE_H__

#include <stdint.h>
#include <vector>
#include <string>

#include "stats/matrix.h"

struct signal_list_t;
struct interval_t;
struct edf_t;
struct timeline_t;


class slice_t
{

 public:
  
  
  slice_t( edf_t & edf , 
	   int signal , 
	   const interval_t & interval , 
	   int    downsample = 1 );
  
  const std::vector<double> * pdata() const 
  { 
    return &data; 
  }
  
  std::vector<double> * nonconst_pdata() 
  { 
    return &data; 
  }

  const std::vector<uint64_t> * ptimepoints() const 
  { 
    return &time_points; 
  }
  
  const std::vector<int> * precords() const 
  { 
    return &records;
  }

  int size() const 
  { 
    return data.size(); 
  }
  
  interval_t duration() const ;

 private:

  // input
  edf_t & edf;
  const int signal;
  const interval_t & interval;
  const int downsample;
  
  // output
  std::vector<double> data;
  std::vector<uint64_t> time_points;
  std::vector<int> records;
  
  double start, stop;
  
};

class mslice_t {
  
 public:
  
  mslice_t( edf_t & edf , 
	    const signal_list_t & , 
	    const interval_t & interval , 
	    int    downsample = 1 );
  
  ~mslice_t()
    {
      for (int s=0;s<channel.size();s++)
	{
	  if ( channel[s] != NULL ) delete channel[s];
	  channel[s] = NULL;
	}
    }
  
  std::vector<slice_t*> channel;
  std::vector<std::string> labels; 

  Data::Matrix<double> extract();

  int size() const { return channel.size(); } 

  std::string label(const int s) const { return labels[s]; } 
  
};


#endif
