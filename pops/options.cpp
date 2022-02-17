
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

#include "pops/options.h"

#include "lgbm/lgbm.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "stats/eigen_ops.h"
#include "db/db.h"

extern logger_t logger;
extern writer_t writer;

#include "eval.h"


std::string pops_opt_t::pops_path;
std::string pops_opt_t::pops_root;

bool pops_opt_t::verbose;

int pops_opt_t::trim_wake_epochs;
int pops_opt_t::n_stages;

bool pops_opt_t::welch_median;
double pops_opt_t::lwr;
double pops_opt_t::upr;
double pops_opt_t::spectral_resolution;
double pops_opt_t::fft_seg_sec;
double pops_opt_t::fft_inc_sec;

double pops_opt_t::epoch_len;
double pops_opt_t::epoch_inc;

bool pops_opt_t::epoch_level_SHAP;

std::vector<double> pops_opt_t::slope_range{ 30.0 , 45.0 } ;
double pops_opt_t::slope_th  = 3;
double pops_opt_t::slope_epoch_th = 5;

std::map<std::string,std::set<std::string> > pops_opt_t::aliases;

void pops_opt_t::set_options( param_t & param )
{
  
  // i.e. prepend this to any non-absolute paths;
  pops_path = param.has( "path" ) ? param.value( "path" ) : "" ;
  
  // assume <path/root>.ftr
  // assume <path/root>.mod
  // assume <path/root>.conf
  pops_root = param.has( "root" ) ? param.value( "root" ) : "" ;  
  if ( pops_root != "" && pops_path != "" ) 
    Helper::halt( "can only specify 'root' or 'path'" );
  
  verbose = param.has( "verbose" );
  
  epoch_level_SHAP = param.has( "epoch-SHAP" );
  
  n_stages = param.has( "3-class" ) ? 3 : 5;

  trim_wake_epochs = param.has( "trim" ) ? param.requires_int( "trim" ) : -1;
    
  welch_median = param.yesno( "fft-median" );
  
  lwr = param.has( "lwr" ) ? param.requires_dbl( "lwr" ) : 0.5;
  
  upr = param.has( "upr" ) ? param.requires_dbl( "upr" ) : 45;
  
  fft_seg_sec = param.has( "segment-sec" ) ? param.requires_dbl( "segment-sec" ) : 4 ; 

  fft_inc_sec = param.has( "segment-overlap" ) ? param.requires_dbl( "segment-overlap" ) : 2 ; 
  
  spectral_resolution = 1.0 / fft_seg_sec ; 

  // these set via main epoch mechanism (edf.timeline) [ eval.cpp ] prior to calling pops_t
  epoch_len = 30;  
  epoch_inc = 30;
  
  // channel aliases
  //  (added when reading spec.) 
  if ( param.has( "alias" ) ) 
    {
      std::vector<std::string> tok = param.strvector( "alias" );
      // primary|second,primary|secondary
      for (int i=0; i<tok.size(); i++)
	{
	  std::vector<std::string> tok2 = Helper::parse( tok[i] , "|=" );
	  if ( tok2.size() < 2 ) Helper::halt( "bad format for alias=main|second,main2=second2" );
	  for (int j=1; j<tok2.size(); j++)
	    aliases[ tok2[0] ].insert( tok2[j] );
	}
    }
  

}

#endif
