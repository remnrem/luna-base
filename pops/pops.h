
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


#ifdef HAS_LGBM

#ifndef __LUNA_POPS_H__
#define __LUNA_POPS_H__

#include "stats/lgbm.h"
#include "stats/Eigen/Dense"

struct param_t;
struct edf_t;

#include "pops/options.h"
#include "pops/spec.h"
#include "pops/indiv.h"

enum pops_stage_t
  {
    POPS_WAKE = 0 ,
    POPS_REM = 1 ,
    POPS_N1 = 2 , // also, NR for 3-stage model
    POPS_N2 = 3 ,
    POPS_N3 = 4 ,
    POPS_UNKNOWN = 9 
  };



struct pops_t {
  
  // two main entry points

  // 1) when predicting, we need to load the LGBM model
  void load_model( param_t & );

  
  // 2) main wrapper:: create a level 2 feature library and save
  // i.e. for trainer and/or validation library
  void make_level2_library( param_t & );

  // load level 1 data
  void load1( const std::string & f );

  // derive level 2 stats (from pops_t::specs)
  void level2();

  // fit and save a LGBM model (--> pops_t::lgbm)
  void fit_model( const std::string & f );

    

  
  static pops_opt_t opt;
  
  static lgbm_t lgbm;

  static pops_specs_t specs;

  //
  // cohort-level data 
  //

  Eigen::MatrixXd X1;
  Eigen::MatrixXd X2;
  Eigen::MatrixXd X;
  std::vector<int> S;
  std::vector<int> E;
  std::vector<int> Istart, Iend;

  //
  // helpers
  //
  
  static void outliers( const Eigen::VectorXd & x ,
			const double d ,
			const std::vector<int> & staging , 
			std::vector<int> * staging2 ); 
  
  static std::string label( pops_stage_t s )
  {
    if ( s == POPS_N1 ) return pops_opt_t::n_stages == 3 ? "NR" : "N1" ;
    if ( s == POPS_N2 ) return "N2";
    if ( s == POPS_N3 ) return "N3";
    if ( s == POPS_REM ) return "R";
    if ( s == POPS_WAKE ) return "W";
    return "?";
  }

  
};



#endif
#endif
