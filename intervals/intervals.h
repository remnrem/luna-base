
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

#ifndef __INTERVALS_H__
#define __INTERVALS_H__

#include "defs/defs.h"

#include <stdint.h>

#include <sstream>

#include <ostream>
#include <set>


//
// intervals are defined from the start (inclusive) to one point past the end 
// zero-duration intervals are valid ( start == end ) however
//


struct interval_t
{
  
  friend std::ostream & operator<<( std::ostream & out , const interval_t & rhs );
  
  interval_t() { start = 0 ; stop = 0; }
  
  interval_t( uint64_t start , uint64_t stop ) : start(start) , stop(stop) { } 
  
  bool empty() const { return start == 0 && stop == 0; } 
  
  uint64_t duration() const { return stop - start ; } 
  
  double duration_sec() const { return ( stop - start ) / (double)globals::tp_1sec; } 
  
  void set_leftright( uint64_t l , uint64_t r ) 
  {
    start = l;
    stop  = r; // assumes this is 1 past the end already
    if ( start > stop ) 
      {
	uint64_t t = stop;
	stop = start;
	start = t; 
      } 
  }
  
  void set_window( uint64_t c , uint64_t w )
  {    
    uint64_t w2 = w * 0.5;
    start = w2 > c ? 0 : c - w2 ;
    stop  = c + w2 + 1LLU ; // on point past the end
  }
  
  void expand( uint64_t w )
  {
    if ( start >= w ) start -= w; 
    else start = 0;
    stop += w;
  }

  void shift_left( uint64_t w )
  {
    if ( start >= w ) start -= w;
    else start = 0;
    stop -= w;
  }

  void shift_right( uint64_t w )
  {
    start += w;
    stop += w;
  }
  
  uint64_t start;
  
  uint64_t stop;

  
  double start_sec() const { return start/(double)globals::tp_1sec; }

  double mid_sec() const { return mid()/(double)globals::tp_1sec; }

  // one time point past end
  double stop_sec() const { return stop/(double)globals::tp_1sec; }

  // for better display, gives stop exactly at last time-point
  // but handle empty and 0-duration timestamps
  double stop_sec_exact() const
  {
    return stop == 0 ? 0
      : ( start == stop ? stop/(double)globals::tp_1sec : (stop-1LLU)/(double)globals::tp_1sec ); 
  }
  
  bool operator<( const interval_t & rhs ) const 
  {
    if ( start == rhs.start ) return stop < rhs.stop;
    return start < rhs.start;
  }

  bool operator>( const interval_t & rhs ) const 
  {
    if ( start == rhs.start ) return stop > rhs.stop;
    return start > rhs.start;
  }


  bool overlaps( const interval_t & b ) const 
  {
    // note: window defined as start .. (stop-1)
    // allow for special case where 0-point duration at first point of second interval
    return ( start < b.stop && stop > b.start ) || start == b.start ;    
  }

  double prop_overlap( const interval_t & b ) const 
  {

    if ( ! overlaps( b ) ) return 0;

    uint64_t min_start = start < b.start ? start  : b.start;
    uint64_t max_start = start < b.start ? b.start : start ;  
	  
    uint64_t min_stop = stop < b.stop ? stop  : b.stop ;
    uint64_t max_stop = stop < b.stop ? b.stop : stop ; 
    
    uint64_t o_intersection = min_stop - max_start ;
    uint64_t o_union        = max_stop - min_start ;

    // two completely overlapping timestamps
    if ( o_union == 0 ) return 1.0;
    double metric = o_intersection / (double)o_union;	      
    
    return metric;
  }


  bool is_after( const interval_t & b ) const 
  {
    // as stop is 1-point after the end, use <=
    return start >= b.stop;
  }
  
  bool is_before( const interval_t & b ) const 
  {
    // as stop is 1-point after the end, use <=
    return stop <= b.start;
  }

  bool is_completely_spanned_by( const interval_t & b ) const
  {    
    return b.start <= start && b.stop >= stop;
  }
  
  bool contains( const uint64_t & tp ) const
  {
    // note: stop is 1 past the end
    // but allow for 0-duration time-stamp, i.e. explicit
    // test for tp == start as well (as here, tp will also == stop )
    
    return tp == start || ( tp > start && tp < stop );
  }
  
  uint64_t mid() const
  {
    // handle 0-dur timestamps
    return stop == start ? start : start + ( stop - start ) / (uint64_t)2;
  }

  std::string as_string( const int prec = 2 ) const 
  {
    std::stringstream ss;
    double start_sec = start / (double) globals::tp_1sec;
    double stop_sec  = stop  / (double) globals::tp_1sec;
    
    ss.precision( prec );
    ss << std::fixed << start_sec << "->" << stop_sec;
    return ss.str();
  }
  
  std::string as_tp_string() const 
  {
    std::stringstream ss;
    ss << start << "->" << stop;
    return ss.str();
  }

  static int intersect( const std::set<interval_t> & a, 
			const std::set<interval_t> & b, 
			std::set<interval_t> * botha,
			std::set<interval_t> * bothb,
			std::set<interval_t> * cons,
			std::set<interval_t> * uns,
			std::set<interval_t> * onlya,
			std::set<interval_t> * onlyb, 
			double th = 0.5 , uint64_t win = 0 );
  
  
};



struct feature_t 
{

  feature_t() 
  {
    has_value = false;
    has_colour = false;  
  } 

  feature_t( const interval_t & i , const std::string & s , const std::string & l ) 
  : feature(i),  label(l) , signal(s) 
  { 
    has_value = false;
    has_colour = false;
    colour = "";
    value = 0;
  } 

  // actual feature size
  interval_t  feature;

  // optionally, a bounding interval,
  // i.e. to define 'near overlap' 
  interval_t  window;
  
  std::string label;
  std::string signal;  

  bool has_colour;
  std::string colour;

  bool has_value;
  double value;
  
  std::map<std::string,std::string> data;

  void add_data( const std::map<std::string,std::string> & m )
  {
    std::map<std::string,std::string>::const_iterator mm = m.begin();
    while ( mm != m.end() )
      {
	data[ mm->first ] = mm->second;
	++mm;
      }
  }

  void clear_data() { data.clear(); } 

  std::string print_data( const std::string & delim1 = " = " , 
			  const std::string & delim2 = "\n" , 
			  const std::string & pre = "  " ) const
  {
    std::stringstream ss;
    std::map<std::string,std::string>::const_iterator dd = data.begin();
    while ( dd != data.end() ) 
      {
	ss << pre << dd->first << delim1 << dd->second << delim2;
	++dd;
      }
    return ss.str();
  }

  std::string as_string( const std::string & delim = "|" ) const;

  bool operator<( const feature_t & rhs ) const
  {
    if ( feature < rhs.feature ) return true;
    if ( rhs.feature < feature ) return false;
    
    if ( signal < rhs.signal ) return true;
    if ( signal > rhs.signal ) return false;
    
    return label < rhs.label;
  }
  
};


#endif
