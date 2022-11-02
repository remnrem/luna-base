
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

#ifndef __LATERAL_ASYMM_H__
#define __LATERAL_ASYMM_H__

#include <map>
#include <set>
#include <vector>
#include <string>

struct param_t;
struct edf_t;

struct lat_results_t {
  
};

struct lat_t
{

  enum stg_t {
    ASYMM_SS_IGNORE ,
    ASYMM_SS_WAKE ,
    ASYMM_SS_NREM ,
    ASYMM_SS_REM     
  };

  // set up
  lat_t( edf_t & , param_t & );

  // do analysis
  void proc( edf_t & , param_t & );

  // helper
  bool eval_transitions( const std::vector<double> & log2lr , 
			 const bool permute , 
			 const std::vector<bool> & outlier , 
			 const int tr_start , 
			 std::vector<double> * res_tr_R2NR_R_mean , 
			 std::vector<double> * res_tr_R2NR_NR_mean , 
			 std::vector<double> * res_tr_NR2R_NR_mean , 
			 std::vector<double> * res_tr_NR2R_R_mean , 
			 int * res_tr_R2NR_cnt , 
			 int * tr_NR2R_cnt );

  // options
  bool epoch_level_output;
  bool trans_level_output;
  int tr_start;
  int nreps;
  
  // main data stores:
  
  // freq -> epoch -> channel -> power
  std::map<double,std::map<int,std::map<std::string,double> > > f2e2ch2psd;

  // band -> epoch -> channel -> power
  std::map<std::string,std::map<int,std::map<std::string,double> > > b2e2ch2psd;
  
  // stages
  std::vector<stg_t> S;

  // transitions 
  std::vector<int> T_R2NR, T_NR2R;

  // cyles
  int num_cycles;
  std::vector<int> C;
  
  // epochs
  std::vector<int> E;
  
  // channel pairs/groups: if C3+F3+P3 means sum of all 
  std::vector<std::string> left;
  std::vector<std::string> right;
  std::map<std::string,std::set<std::string> > sum_parts; // if A+B --> { A, B } 
  
  // proc a single F/B and channel/group pair
  lat_results_t analyse( const std::vector<double> & l ,
			 const std::vector<double> & r );
  

  
};


#endif
