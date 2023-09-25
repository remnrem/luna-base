
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
    chs.clear();   // separate out channels from other strata
    strata.clear();
    coef = mean = 0;
    sd = 1; 
  }
  
  model_term_t( const std::string & label ,
		const std::string & cmd ,
		const std::string & var ,
		const std::string & s ,
		const std::string & c , 
		double coef ,
		double mean ,
		double sd )  
		: label( label ) , cmd( cmd ) , var( var ) ,
    strata( Helper::mapize( s , ',', '/' ) ),
    chs( Helper::parse( c , ',' ) ),
    coef( coef ) , mean( mean ) , sd( sd )       
  {    
  }
  
  std::string label;
  std::string cmd;
  std::string var;
  std::map<std::string,std::string> strata;
  std::vector<std::string> chs;

  double coef;
  double mean;
  double sd;
  
  bool operator<( const model_term_t & rhs ) const { return label < rhs.label; } 
};


struct prediction_model_t {

  prediction_model_t()
  {
    title = ".";
    terms.clear();
    specials.clear();
    specials_str.clear();
  }

  std::string title; 
  
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
  
};


struct prediction_t {

  prediction_t( edf_t & , param_t & );  

  double prediction();
  
private:

  void output() const;
  

  std::string id;
  
  prediction_model_t model;
  
  // derived feature vector (raw, normalized)
  Eigen::VectorXd X;
  Eigen::VectorXd Z;    

  // predicted value (raw, bias-adjusted)
  double y;
  
  double y1; 

};

#endif
