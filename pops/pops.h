
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

#include "lgbm/lgbm.h"
#include "stats/Eigen/Dense"

struct param_t;
struct edf_t;

#include "pops/options.h"
#include "pops/spec.h"
#include "pops/indiv.h"

struct pops_indiv_t;

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
  
  pops_t() { } 
  
  pops_t( param_t & param ) 
  {
    pops_opt_t::set_options( param );
  }

  // two main entry points
  
  // 1) when predicting, we need to load the LGBM model
  //  void load_model( param_t & );
  
  // 2) main wrapper:: create a level 2 feature library and save
  // i.e. for trainer and/or validation library
  void make_level2_library( param_t & );

  // load list of validation (vs training) IDs
  void load_validation_ids( const std::string & );

  // load level 1 data
  void load1( const std::string & f );
  
  // load stages/epochs only from level 1 data (for --es-priors standalone only)
  void load1_stages_only( const std::string & f );
  
  // derive level 2 stats (from pops_t::specs)
  // this is also co-opted by prediction mode 
  void level2( const bool training = true , const bool quiet = false );

  // test mean differences by stage (overall, and within person)
  void stage_association();
  void stage_association1( Eigen::VectorXd & ftr , const std::vector<std::string> & ss );
			   
  // dump feature matrix
  void dump_matrix( const std::string & f );

  // dump/read a range file (mean/SD for training features)
  void dump_ranges( const std::string & f );
  static void read_ranges( const std::string & f );
    
  // fit and save a LGBM model (--> pops_t::lgbm)
  void fit_model( const std::string & f , const lgbm_label_t & w );

  // write elapsed-sleep priors
  void write_elapsed_sleep_priors( const std::string & f );

  // as a standalone function
  void make_espriors( param_t & );
  
  // for using level2() in the context of prediction
  void from_single_target( const pops_indiv_t & );
  void copy_back( pops_indiv_t * );  
  

  //  static pops_opt_t opt;
  
  static lgbm_t lgbm;

  static std::string lgbm_model_loaded;
  
  static pops_specs_t specs;

  //
  // elapsed-sleep priors
  //
  
  static Eigen::MatrixXd ES_probs;           // P( ES, %NR, %REM | stage ) 
  static Eigen::VectorXd ES_global_priors;
  static std::vector<double> ES_mins;        // total mins elapsed sleep
  static std::vector<double> ES_prior_nrem;  // prior/recent NREM duration
  static std::map<int,std::map<int,int> > ES_rowmap; // ES-bin, prior_nrem_bin --> row of ES table
  
  
  //
  // cohort-level data 
  //

  Eigen::MatrixXd X1;
  std::map<std::string,Eigen::MatrixXd> V; // when reading in level2 SVD scoring
  std::map<std::string,Eigen::MatrixXd> W;  
  
  std::vector<int> S;
  std::vector<int> E;
  std::vector<int> Istart, Iend;
  std::vector<std::string> I; // trainer/validation IDs, yoked to Istart[]
  std::set<std::string> holdouts; // validation IDs
  int ni_validation; // holdouts actually present
  int nrows_training, nrows_validation;

  //
  // feature ranges
  //
  
  static std::map<std::string,double> range_mean;
  static std::map<std::string,double> range_sd;

  //
  // Get indiv-weights (from ivars)
  //
  
  bool attach_indiv_weights( const std::string & w , bool );
  bool dump_weights();

  void report_counts();
  
  //
  // Sample fixed number of obs per stage
  //

  void sample_fixed_n();
  
  //
  // helpers
  //

  static std::string update_filepath( const std::string & s );
  
  static Eigen::MatrixXd add_time_track( const int , const int );

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

  static std::string label5( pops_stage_t s )
  {
    if ( s == POPS_N1 ) return "N1" ;
    if ( s == POPS_N2 ) return "N2";
    if ( s == POPS_N3 ) return "N3";
    if ( s == POPS_REM ) return "R";
    if ( s == POPS_WAKE ) return "W";
    return "?";
  }

  static std::string label3( pops_stage_t s )
  {
    if ( s == POPS_N1 ) return "NR" ;
    if ( s == POPS_REM ) return "R";
    if ( s == POPS_WAKE ) return "W";
    return "?";
  }

  static std::vector<std::string> labels5; 
  static std::vector<std::string> labels3; 
  
  static std::map<int,std::map<int,int> > tabulate( const std::vector<int> & a , 
						    const std::vector<int> & b , 
						    const bool print  );
	     
	     

  static std::vector<int> NRW( const std::vector<int> & s ) 
  {
    std::vector<int> r = s;
    for (int i=0; i<s.size(); i++)
      if ( r[i] == POPS_N2 || r[i] == POPS_N3 ) 
	r[i] = POPS_N1;
    return r;  
  }
  
};


struct pops_nan_report_t {
  pops_nan_report_t( const Eigen::MatrixXd & m );
  bool any() const;
  std::map<int,int> rows; 
  std::map<int,int> cols;
};

struct pops_stats_t { 

  pops_stats_t( const std::vector<int> & obs , 
		const std::vector<int> & pred , 
		const int nstages = 5 , 
		const int type = 0 , 
		const int ostage = -1 );
  
  
  bool valid;
  int n; // 5 or 3 (stage numbers)
  int nobs; // number of valid epochs
  double kappa, acc, mcc;
  
  double macro_precision, macro_recall, macro_f1 ;
  double avg_weighted_precision, avg_weighted_recall, avg_weighted_f1;  
  std::vector<double> precision, recall, f1;
  
};

#endif
#endif
