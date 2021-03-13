
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
#include "stats/Eigen/Dense"

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
	   int    downsample = 1 , 
	   bool digital = false );
  
  const std::vector<double> * pdata() const 
  { 
    return &data; 
  }
  
  std::vector<double> * nonconst_pdata() 
  { 
    return &data; 
  }

  // digitial (16-bit) EDF type
  const std::vector<int16_t> * ddata() const 
  { 
    return &dig_data; 
  }
  
  std::vector<int16_t> * nonconst_ddata() 
  { 
    return &dig_data; 
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

  void clear()
  {
    data.clear();
    dig_data.clear();
    time_points.clear();
    records.clear();
    start = stop = 0;
  }

 private:
  
  // input
  edf_t & edf;
  const int signal;
  const interval_t & interval;
  const int downsample;
  
  // output
  std::vector<double> data; // unless digital == T 
  std::vector<int16_t> dig_data;// optional, if digital == T
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

  void clear() // nb. leaves in an invalid state, so do not try to reuse...
  { 
    for (int s=0;s<channel.size();s++)
      {
	if ( channel[s] != NULL ) delete channel[s];
	channel[s] = NULL;
      }
  }
  
};




class matslice_t {
  
 public:
  
  matslice_t( edf_t & edf , 
	      const signal_list_t & , 
	      const interval_t & interval );
  
  ~matslice_t()
    {
      clear();      
    }

  
  const std::vector<double> * col( const int s ) const { return data.col_pointer(s)->data_pointer(); } 
  
  const Data::Matrix<double> & data_ref() const { return data; } 

  Data::Matrix<double> & nonconst_data_ref() { return data; } 

  int size() const { return labels.size(); } 
  
  std::string label(const int s) const { return labels[s]; } 

  const std::vector<uint64_t> * ptimepoints() const 
  { 
    return &time_points; 
  }

  void clear()
  { 
    data.clear();
    labels.clear();
    time_points.clear();
  }

 private:

  Data::Matrix<double> data;  

  std::vector<uint64_t> time_points;

  std::vector<std::string> labels; 

  
};



class eigen_matslice_t {
  
 public:
  
  eigen_matslice_t( edf_t & edf , 
		    const signal_list_t & , 
		    const interval_t & interval );
  
  ~eigen_matslice_t()
  {
    clear();      
  }
  
  
  const Eigen::MatrixXd & data_ref() const { return data; } 

  Eigen::MatrixXd & nonconst_data_ref() { return data; } 

  int size() const { return labels.size(); } 
  
  std::string label(const int s) const { return labels[s]; } 

  const std::vector<uint64_t> * ptimepoints() const 
  { 
    return &time_points; 
  }

  void clear()
  { 
    data.resize(0,0);
    labels.clear();
    time_points.clear();
  }

 private:

  Eigen::MatrixXd data;  

  std::vector<uint64_t> time_points;

  std::vector<std::string> labels; 

  
};






#endif
