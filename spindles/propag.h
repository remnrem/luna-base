
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


#ifndef __SPROP_H__
#define __SPROP_H__

#include "spindles/spindles.h"

#include <string>
#include <vector>

struct sp_idx_t { 
  
  sp_idx_t( uint64_t fe9 , 
  	    const std::string & ch )
  : fe9(fe9) , ch(ch) 
  { }
  
  uint64_t fe9;
  std::string ch;
  
  bool operator<( const sp_idx_t & rhs ) const 
  {
    if ( fe9 < rhs.fe9 ) return true;
    if ( fe9 > rhs.fe9 ) return false;
    return ch < rhs.ch ;
  }

};

struct sp_dat_t { 
  
  sp_dat_t() { } 

  sp_dat_t( const std::vector<spindle_t> & sp, 
  	    const std::vector<double> & coeff )
  : sp( sp ) , coeff( coeff ) 
  { }
  
  std::vector<spindle_t> sp;
  std::vector<double> coeff;

};
  

struct sp_props_t { 

  // store as, does not change
  void add_tp ( const std::vector<uint64_t> & tp );

  void add( double f , const std::string & ch , const std::vector<spindle_t> & sp, 
	    const std::vector<double> & cwt );

  double analyse( const std::set<double> & f , 
		  const std::set<std::string> & c , 
		  const std::string & seed , 
		  const double w = 1.0 ,
		  const bool verbose = false 
		  );
  
  // members
  std::vector<uint64_t> tps;
  std::map<sp_idx_t,sp_dat_t> data;
};




#endif
