
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


#ifndef __SIMASSOC_H__
#define __SIMASSOC_H__

// internal tool: to simulate and assess association between spectra and quantitative
// phenotypes using;  not supported or designed for external use otherwise

#include <map>
#include <vector>
#include <set>
#include <string>
#include "../stats/Eigen/Dense"

struct param_t;

struct simassoc_t {

  void load( const std::string & fileroot );

  void load_covar( const std::string & fileroot );

  void describe_cols();
  
  void generative_model();
  
  void simulate();
  
  void assoc();
  
  void output();
  
  // data (power)
  Eigen::MatrixXd X;

  // optional covariates
  Eigen::MatrixXd Z;
  
  // dictionary (type -> cols)
  std::vector<std::string> hdr;
  std::map<std::string,std::set<int> > dict;

  // true coefficients (that may be indiv-specific, thus
  // a matrix the same size as X
  Eigen::MatrixXd W;
  
  // phenotype
  Eigen::VectorXd Y;
  
  //
  // parameters
  //

  //
  // generation
  //

  double var_exp = 0.05;

  enum genmod_t {

    // simple BAND 
    BAND , 
    // i1 = 0, 1, ..., 6 = SLOW, DELTA, ... , BETA

    // shifted BAND (still 0/1 weights)
    //  but they are allowed to vary between individuals
    SHIFTED_BAND ,
    // i1 = band, as above
    // p1 = shift_left (mean) 
    // p2 = shift_left (SD) // i.e. rows in W
    // p3 = shift_right (mean)
    // p4 = shift_right (SD) // i.e. rows in W

    // Gaussian dist (e.g. mean = 8 Hz, SD = +/- 2 Hz)
    //  nb. if SD = 0, this is taken as a single bin
    SINGLE_BIN  ,
    // p1 = mean freq (scaled so W = 1 here)
    // p2 = SD 
        
  };

  

  // parameters
  int i1, i2, i3;
  double p1, p2, p3, p4;
  std::string s1, s2;

  
  //
  // association testing
  //

  // number of replicates to get empirical p-value per 
  int emp_nreps = 1000;

  // threshold for inner loop
  double emp_alpha = 0.05;
  
  //
  // power
  //

  // total number of replicates (to assess power) 
  int nreps = 100;

  // 
  double alpha = 0.05;
  
  
};


void build_param_from_cmdline( param_t * param );


#endif
