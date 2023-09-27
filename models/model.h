
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

#ifndef __PRED_MOD_H__
#define __PRED_MOD_H__

#include "stats/Eigen/Dense"
#include <map>
#include <vector>
#include <set>
#include <string>
#include "helper/helper.h"

struct edf_t;
struct param_t;


struct model_term_t {

  model_term_t()
  {
    clear();
  }

  void clear()
  {
    label = cmd = var = ".";
    // if cmd == var == ".", implies a value was set (not from cache) 
    chs.clear();   // separate out channels from other strata
    strata.clear();
    coef = mean = 0;
    sd = 1;
    required = false;
    has_value = false;
    value = ".";
  }
  
  std::string label;
  std::string cmd;
  std::string var;
  std::map<std::string,std::string> strata;
  std::vector<std::string> chs;
  
  double coef;
  double mean;
  double sd;
  bool required;
  
  // if the model file sets a VALUE (i.e. instead of specify a cache location
  // store it here (on reading model file - which can contain indiv-specific
  // variable substitutions) before passing to the feature vector; in this case,
  // CMD will be '.'
  
  bool has_value;

  std::string value; // may be missing, so treat as string here
  
  bool operator<( const model_term_t & rhs ) const { return label < rhs.label; } 

};


struct prediction_model_t {

  prediction_model_t()
  {
    terms.clear();
    specials.clear();
    specials_str.clear();
  }

  void read( const std::string & f , const std::string & id );

  int size() const { return terms.size(); }

  void populate();

  void dump() const;

  std::set<std::string> channels() const;
  
  std::set<model_term_t> terms;

  Eigen::VectorXd coef;
  
  Eigen::VectorXd mean;
  
  Eigen::VectorXd sd;

  // special variables

  std::map<std::string,double> specials;

  std::map<std::string,std::string> specials_str;

  // helper
  std::string s( const std::string & k )
  {
    std::map<std::string,std::string>::const_iterator ss = specials_str.find( k );
    if ( ss == specials_str.end() ) return ".";
    return ss->second;
  }

  // required feature?
  int min_features() const
  {
    std::map<std::string,double>::const_iterator ss = specials.find( "minf" );
    if ( ss == specials.end() ) return ss->second;
    return ss->second;
  }

  std::vector<std::string> header() const;
  
};


#endif
