
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

bool pops_opt_t::run_stage_associations;

double pops_opt_t::ES_es_tbin;
double pops_opt_t::ES_nr_tbin;

double pops_opt_t::ES_es_tmax;
double pops_opt_t::ES_nr_tmax;

double pops_opt_t::ES_non_NREM_mins;
double pops_opt_t::ES_c;

bool pops_opt_t::ES_rolling;
bool pops_opt_t::ES_fractional_count;

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

std::vector<std::map<std::string,std::string> > pops_opt_t::equivs;
std::map<std::string,std::string> pops_opt_t::equiv_swapins;
std::string pops_opt_t::equiv_label;

std::vector<std::string> pops_opt_t::iweights;
bool pops_opt_t::dump_model_weights;
std::string pops_opt_t::model_weights_file;

bool pops_opt_t::sample_fixed_n;
std::vector<int> pops_opt_t::fixed_n;

bool pops_opt_t::soap_results;
double pops_opt_t::soap_threshold;
int pops_opt_t::soap_nc;
bool pops_opt_t::soap_grid;
double pops_opt_t::lk_lwr;
double pops_opt_t::lk_upr;
double pops_opt_t::lk_steps;
double pops_opt_t::soap_grid_mean_conf;

bool pops_opt_t::eval_mode = false;


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
  //        path/lib.priors
  //        path/(SVD files)

  // under root-specification, able to use/not use ranges, es-priors
  if_root_apply_ranges = param.has( "apply-ranges" ) ? param.yesno( "apply-ranges" ) : true ;
  if_root_apply_espriors = param.has( "apply-priors" ) ? param.yesno( "apply-priors" ) : false ;

   
  // intercept (i.e. to avoid 0-weight probs for any cell)
  ES_c = param.has( "priors-c" )  ? param.requires_dbl( "priors-c" ) : 0.001 ;
  
  ES_rolling = param.yesno( "priors-rolling" );

  ES_fractional_count = param.yesno( "priors-weighted" );

  // elapsed sleep priors

  // bin size (mins) &  max time (mins)
  ES_es_tbin = param.has( "priors-es-min" ) ? param.requires_dbl( "priors-es-min" ) : 20 ;
  ES_es_tmax = param.has( "priors-es-max" ) ? param.requires_dbl( "priors-es-max" ) : 380 ;

  ES_nr_tbin = param.has( "priors-nr-min" ) ? param.requires_dbl( "priors-nr-min" ) : 10 ;
  ES_nr_tmax = param.has( "priors-nr-max" ) ? param.requires_dbl( "priors-nr-max" ) : 60 ;
  ES_non_NREM_mins = param.has( "priors-nr-allow" ) ? param.requires_dbl( "priors-nr-allow" ) : 5 ; 
  
  // NOT USED NOW: intercept (i.e. to avoid 0-weight probs for any cell)
  ES_c = param.has( "priors-c" )  ? param.requires_dbl( "priors-c" ) : 0.001 ;  

  
  // vars

  if ( param.has( "inc-vars" ) ) inc_vars = param.strset( "inc-vars" );
  if ( param.has( "exc-vars" ) ) exc_vars = param.strset( "exc-vars" );

  // SOAP

  soap_results = param.has( "soap" );

  soap_threshold = param.empty( "soap" ) ? 0.5 : param.requires_dbl( "soap" ); 

  soap_nc = param.has( "soap-nc" ) ? param.requires_int( "soap-nc" ) : 10 ;

  lk_lwr = param.has( "soap-lwr" ) ? param.requires_dbl( "soap-lwr" ) : 1;
  lk_upr = param.has( "soap-upr" ) ? param.requires_dbl( "soap-upr" ) : 100;
  lk_steps = param.has( "soap-steps" ) ? param.requires_int( "soap-steps" ) : 100;
  soap_grid = param.has( "soap-grid" );
  soap_grid_mean_conf = param.has( "soap-grid" ) && ! param.empty( "soap-grid" ) ? param.requires_dbl( "soap-grid" ) : 0.8; 


  // misc
  
  verbose = param.has( "verbose" );

  run_stage_associations = param.has( "stage-assoc" ) ? param.yesno( "stage-assoc" ) : false ;
  
  epoch_level_SHAP = param.has( "epoch-SHAP" ) || param.has( "SHAP-epoch" ) ;
  
  n_stages = param.has( "3-class" ) ? 3 : 5;

  trim_wake_epochs = param.has( "trim" ) ? param.requires_int( "trim" ) : -1;
  
  welch_median = param.has( "fft-median" ) ? param.yesno( "fft-median" ) : true;
  
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

  // randomly select exactly N of each class? fix=W,R,N1,N2,N3
  if ( param.has( "fix" ) )
    {
      fixed_n = param.intvector( "fix" );
      if ( fixed_n.size() != 5 ) Helper::halt( "expecting 5 values for fix=W,R,NR..." );
      sample_fixed_n = true;
    }
  
  // channel aliases
  //  (added when reading spec.) 
  aliases.clear();
  if ( param.has( "alias" ) ) 
    {
      
      // expecting single |, two cols
      std::vector<std::string> tok = Helper::parse( param.value( "alias" ) , "|" );
      
      if ( tok.size() != 2 ) 
	Helper::halt( "bad format for alias=main,main2,...|second,second2,..." );
      
      // primary1,primary2,...|secondary1,secondary2,...
      // i.e. both sides should have same number of elements
      std::vector<std::string> pri = Helper::parse( tok[0] , "," );
      std::vector<std::string> sec = Helper::parse( tok[1] , "," );
      if ( pri.size() != sec.size() )
	Helper::halt( "bad format for alias=main,main2,...|second,second2,..." );
      
      
      // but to map multiple aliases, can do =
      // pri1,pri2|sec1=sec1b=sec1c,sec2
      // i.e. maps sec1, sec1b and sec1c --> pri
      //      maps sec2 --> pri2
      
      for (int i=0; i<pri.size(); i++)
	{
	  std::vector<std::string> sec2 = Helper::parse( sec[i] , "=" );
	  for (int j=0; j<sec2.size(); j++)
	    aliases[ pri[i] ].insert( sec2[j] );
	}
    }
  
  // channel replacements : 
  // i.e. if feature has C4_M1, but we want to use C3_M2 and *not* C4_M1 (i.e. not as an 'equivalent' channel)
  replacements.clear();
  replacements_rmap.clear();
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
  //  we allow for multiple channels to be rotated, but all originals
  //  must have the same number of alternatives
  //   e.g.
  //      equiv=C4,F4|C3,F3|C1,F1|C2,F2
  
  //   means default channels (in the .ftr file) are 'C4' and 'F4'
  //    but three other models will be run (with 'C3' and 'F3', then with
  //    'C2' and 'F2' etc)
  //   and the final predictions will be based on the set of four predictions

  //     -->    equiv  C4|C3|C1|C3|P4|F3|F4
  //      i.e. if we have those, then try doing everything with that

  //   if an equivalence channel does not exist, just skip that step (for all)   
  
  // default:no equivalence channel
  equiv_swapins.clear();
  equivs.clear();
  
  if ( param.has( "equiv" ) )
    {
      // expect | delimited equivalence sets
      std::vector<std::string> eqs = Helper::parse( param.value( "equiv" ) , "|" );
      if ( eqs.size() < 2 )
	Helper::halt( "equiv requires two or more sets of channels" );
      
      // within each |, we expect the same number of channels (comma-delimited)
      std::vector<std::string> originals = Helper::parse( eqs[0] , "," );
      const int neq = originals.size();
      
      // add 'self' 
      std::map<std::string,std::string> eq1;
      for (int k=0; k<neq; k++)
	eq1[ originals[k] ] = originals[k];
      equivs.push_back( eq1 );
      
      // track the originals (but as a set)
      //equiv_root = chs[0];
      
      // add each equivalence set
      for (int j=1; j<eqs.size(); j++)
	{
	  std::vector<std::string> eq = Helper::parse( eqs[j] , "," );
	  if ( eq.size() != neq ) 
	    Helper::halt( "same number of equiv channels must be specified each set:\n" 
			  + eqs[0] + "\n" + eqs[j] );	  

	  std::map<std::string,std::string> eq1;
	  for (int k=0; k<neq; k++)
	    eq1[ originals[k] ] = eq[k];
	  equivs.push_back( eq1 );
	  
	}      

    }

}

#endif
