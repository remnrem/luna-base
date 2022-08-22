
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
#include "pops/pops.h"

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

bool pops_opt_t::if_root_apply_ranges;
bool pops_opt_t::if_root_apply_espriors;

std::set<std::string> pops_opt_t::inc_vars;
std::set<std::string> pops_opt_t::exc_vars;

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
std::map<std::string,std::string> pops_opt_t::replacements;
std::map<std::string,std::string> pops_opt_t::replacements_rmap; // track reverse mapping

std::vector<std::string> pops_opt_t::equivs;
std::string pops_opt_t::equiv_root;
std::string pops_opt_t::equiv_swapin;

std::vector<std::string> pops_opt_t::iweights;
bool pops_opt_t::dump_model_weights;
std::string pops_opt_t::model_weights_file;

bool pops_opt_t::soap_results;
double pops_opt_t::soap_threshold;

void pops_opt_t::set_options( param_t & param )
{
  
  // i.e. prepend this to any non-absolute paths;
  
  pops_path = param.has( "path" ) ? Helper::expand( param.value( "path" ) ) : "" ;
  pops_root = param.has( "lib" ) ? Helper::expand( param.value( "lib" ) ) : "" ;
  
  // assume path/lib.ftr
  //        path/lib.mod
  //        path/lib.conf
  // optional:
  //        path/lib.ranges
  //        path/lib.espriors
  //        path/(SVD files)

  // under root-specification, able to use/not use ranges, es-priors
  if_root_apply_ranges = param.has( "apply-ranges" ) ? param.yesno( "apply-ranges" ) : true ;
  if_root_apply_espriors = param.has( "apply-es-priors" ) ? param.yesno( "apply-es-priors" ) : true;

  pops_t::ES_rolling = param.yesno( "es-rolling" );
  pops_t::ES_fractional_count = param.yesno( "es-weighted" );
  
  // vars

  if ( param.has( "inc-vars" ) ) inc_vars = param.strset( "inc-vars" );
  if ( param.has( "exc-vars" ) ) exc_vars = param.strset( "exc-vars" );

  // SOAP

  soap_results = param.has( "soap" );

  soap_threshold = param.empty( "soap" ) ? 0.5 : param.requires_dbl( "soap" ); 
    

  // misc
  
  verbose = param.has( "verbose" );
  
  epoch_level_SHAP = param.has( "epoch-SHAP" ) || param.has( "SHAP-epoch" ) ;
  
  n_stages = param.has( "3-class" ) ? 3 : 5;

  trim_wake_epochs = param.has( "trim" ) ? param.requires_int( "trim" ) : -1;
  
  welch_median = param.yesno( "fft-median" );
  
  lwr = param.has( "lwr" ) ? param.requires_dbl( "lwr" ) : 0.5;
  
  upr = param.has( "upr" ) ? param.requires_dbl( "upr" ) : 45;
  
  fft_seg_sec = param.has( "segment-sec" ) ? param.requires_dbl( "segment-sec" ) : 4 ; 

  fft_inc_sec = param.has( "segment-overlap" ) ? param.requires_dbl( "segment-overlap" ) : 2 ; 
  
  spectral_resolution = 1.0 / fft_seg_sec ; 

  // if data already epoched, these set via main epoch mechanism (edf.timeline) [ eval.cpp ] prior to calling pops_t
  epoch_len = globals::default_epoch_len;
  epoch_inc = globals::default_epoch_len;
  
  // training weights (for indiv-level vars)
  if ( param.has( "iid-weights" ) )
    iweights = param.strvector( "iid-weights" );

  // dump model weights to file? (after LGBM fitting)
  dump_model_weights = param.has( "dump-weights" );
  model_weights_file = dump_model_weights ? param.value( "dump-weights" ) : "" ; 

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

  // channel replacements : i.e. if feature has C4_M1, but we want to use C3_M2 and *not* C4_M1 (i.e. not as an 'equivalent' channel)
  if ( param.has( "replace" ) )
    {
      if ( param.empty( "replace" ) )
	Helper::halt( "no replace old,new(,old,new,...)" );
      std::vector<std::string> tok = param.strvector( "replace" );
      
      if ( tok.size() % 2 != 0 )
	Helper::halt( "expecting replace=old,new(,old,new) - i.e. an even number of args" );
      
      for (int i=0; i<tok.size(); i+=2 )
	{
	  if (  tok[i] == tok[i+1] )
	    Helper::halt( "invalid replacement (same label)" );
	  replacements[ tok[i] ] = tok[i+1];
	  replacements_rmap[ tok[i+1] ] = tok[i];	  
	}
    }
  
  // channel equivalents
  //  i.e. actually different channels; map to the preferred term in the model file
  //  currently, only allow for a single channel to be rotated (i.e. swap in multiple
  //  equivalent versions, and test for best predictions / consensus)
  //   e.g. train on C4
  //     -->    alias  C4 <- C4_M1 C4_A1 etc
  //     -->    equiv  C4,C3,C1,C3,P4,F3,F4
  //      i.e. if we have those, then try doing everything with that

  // default:no equivalence channel
  equiv_root = equiv_swapin = "";
  equivs.clear();

  if ( param.has( "equiv" ) )
    {
      std::vector<std::string> chs = param.strvector( "equiv" );
      if ( chs.size() < 2 )
	Helper::halt( "equiv requires two or more channels" );

      equiv_root = chs[0];

      // note:: includes self-equiv, i==0
      for (int i=0; i<chs.size(); i++)
	equivs.push_back( chs[i] );
      
    }

}

#endif
