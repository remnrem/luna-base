
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

#ifndef __LUNA_POPS_OPTIONS_H__
#define __LUNA_POPS_OPTIONS_H__

#include <vector>
#include <map>
#include <set>
#include <string>

struct param_t;

struct pops_opt_t {
    
  static void set_options( param_t & );

  // channels
  static std::map<std::string,std::set<std::string> > aliases;
  static std::map<std::string,std::string> replacements;
  static std::map<std::string,std::string> replacements_rmap;
  
  // for a single channel only, can run POPS prediction swapping
  // in multiple channels instead of one 
  static std::vector<std::map<std::string,std::string> > equivs;
  static std::map<std::string,std::string> equiv_swapins;
  static std::string equiv_label;

  // files
  static std::string pops_path;
  static std::string pops_root;

  static bool if_root_apply_ranges;
  static bool if_root_apply_espriors;

  // elapsed sleep priors
  static double ES_es_tbin;
  static double ES_nr_tbin;

  static double ES_es_tmax;
  static double ES_nr_tmax;
  
  static double ES_non_NREM_mins;
  static double ES_c;
  
  static bool ES_rolling;
  static bool ES_fractional_count;
  
  
  // variables
  static std::set<std::string> inc_vars, exc_vars;

  // feature/stage associations
  static bool run_stage_associations;
  
  // misc
  static bool verbose;
  static int n_stages;
  static int trim_wake_epochs;

  static double epoch_len;
  static double epoch_inc;

  static bool welch_median;
  static double lwr;
  static double upr;
  static double fft_seg_sec;
  static double fft_inc_sec;
  static double spectral_resolution;
  static std::vector<double> slope_range;
  static double slope_th;
  static double slope_epoch_th;

  static std::vector<std::string> iweights;
  static bool dump_model_weights;
  static std::string model_weights_file;

  static bool sample_fixed_n;
  static std::vector<int> fixed_n;
  
  // post-SOAP

  static bool soap_results;
  static double soap_threshold;
  static int soap_nc;
  static bool soap_grid;
  static double soap_grid_mean_conf;
  static double lk_lwr;
  static double lk_upr;
  static double lk_steps;
  
  // outputs
  static bool epoch_level_SHAP;


  // no P (i.e. if in eval mode)
  static bool eval_mode;

};

#endif
#endif

