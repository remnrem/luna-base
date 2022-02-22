
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
  
  static std::map<std::string,std::set<std::string> > aliases;

  static std::string pops_path;
  static std::string pops_root;

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
  
  // outputs
  static bool epoch_level_SHAP;
};

#endif
#endif

