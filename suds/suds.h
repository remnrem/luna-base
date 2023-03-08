
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


#ifndef __SUDS_H__
#define __SUDS_H__

#include <string>

#include "stats/Eigen/Dense"

#include <vector>
#include <set>
#include <map>
#include "stats/lda.h"
#include "stats/qda.h"
#include "stats/eigen_ops.h"

#include "eval.h"
#include "helper/helper.h"
#include "miscmath/miscmath.h"
#include "edf/signal-list.h"


// helper to get posteriors from either LDA or QDA

struct posteriors_t {
  
  posteriors_t() { } 
  
  posteriors_t( const lda_posteriors_t & rhs )
  {
    pp = rhs.pp;
    cl = rhs.cl;
    cli = rhs.cli;
  }
  
  posteriors_t( const qda_posteriors_t & rhs )
  {
    pp = rhs.pp;
    cl = rhs.cl;
    cli = rhs.cli;    
  }  
  Eigen::MatrixXd pp;
  std::vector<std::string> cl;
  std::vector<int> cli; // as above, but index
    
};


struct edf_t;

struct param_t;

struct suds_helper_t { 
  
  suds_helper_t(   edf_t & edf,  param_t & param ) 
  : edf(edf) , param(param) , nge(0) , ne(0), ns(0) , trimmed(0) , ambig(0) { } 
  
  edf_t & edf;
  
  param_t & param;
  
  int nge;
  int ne;
  int ns;
  std::vector<bool> retained;
  signal_list_t signals;
  std::string siglab;
  bool has_prior_staging;
  std::vector<bool> valid;

  int trimmed; // by leading/trailing wake & max-epoch reqs
  int ambig;   // by self-prob (SOAP) on trainers

};


enum suds_stage_t
  {
   SUDS_WAKE = 0 ,
   SUDS_N1 = 1 ,
   SUDS_N2 = 2 ,
   SUDS_N3 = 3 ,
   SUDS_NR = 4 , // generic NR (for 3-stage model)
   SUDS_REM = 5 ,
   SUDS_ARTIFACT = 6 , 
   SUDS_UNKNOWN = 7 ,
   SUDS_LIGHTS 
  };

enum suds_feature_t
  {
    SUDS_LOGPSD ,   // LWR , UPR (fixed 0.25 Hz intervals) 
    SUDS_RELPSD ,   // LWR , UPR , NORM-LWR , NORM-UPR ( scale PSD by sum of NORM )  
    SUDS_CVPSD ,    // LWR , UPR 
    SUDS_SLOPE , // Fixed 30-45 Hz ; fixed other param
    SUDS_SKEW , 
    SUDS_KURTOSIS ,
    SUDS_HJORTH ,
    SUDS_FD ,
    SUDS_PE ,
    SUDS_MEAN , 
    SUDS_SMOOTH , // replace
    SUDS_DENOISE ,  // replace
    SUDS_SMOOTH2 , // duplicate
    SUDS_DENOISE2 , // duplicate
    SUDS_TIME // integer: 0, 1(linear), 2(plus cubic) etc
  };


struct suds_spec_t {  
  suds_feature_t ftr;
  std::string ch; 
  std::map<std::string,double> arg;
  Eigen::ArrayXd w;
  int cols(int *) const ;
};


struct suds_channel_t {
  suds_channel_t( const std::string & ch , const int sr ) : ch(ch) , sr(sr) { }
  suds_channel_t() { }
  std::string ch;
  int sr;
};

// for sorting trainers
struct trkap_t {
  trkap_t( const std::string & id , const double k ) : id(id) , k(k) { }
  std::string id;
  double k;
  bool operator<( const trkap_t & rhs ) const {
    if ( k < rhs.k ) return true;
    if ( k > rhs.k ) return false;
    return id < rhs.id ; 			 
  }
};


struct suds_model_t {
  
  void init();
  static std::map<std::string,suds_feature_t> lab2ftr;
  static std::map<suds_feature_t,std::string> ftr2lab;

  // populate specs[]
  bool read( const std::string & ,
	     const std::string & winfile = "" , 
	     const std::string & woutfile = "" ,
	     const std::string & default_channel = "C4_M1" );

  // use default
  void default_model();
  
  bool write( const std::string & ); // why?

  bool loaded() const { return specs.size() != 0 ; } 

  // give columns for a spec/channel combo
  bool has( suds_feature_t , const std::string & ch );
  std::vector<int> cols( suds_feature_t , const std::string & ch );
  std::map<suds_feature_t,std::map<std::string,std::vector<int> > > ftr2ch2col; // track which features/cols selected
  void build_colmap(); // function to construct ftr2ch2col, called after read()
  void check_args();
  
  // channel info
  std::map<std::string,suds_channel_t> chs;
  
  // all features (in order)
  std::vector<suds_spec_t> specs;
  // but, also as a check, each feature can only be applied to each channel once
  std::map<suds_feature_t,std::map<std::string,suds_spec_t> > fcmap;
    
  // number of features
  int cols() const;
  
  // get labels
  std::vector<std::string> labels();
  
  // feature weights (prior to SVD)
  void write_weights( const std::string & weightfile );
  void read_weights( const std::string & weightfile );
  void set_weights();
  Eigen::ArrayXd W;
  
  // SVD parameters
  int nc;
  
};

struct suds_helper_t;

struct suds_indiv_t {

  suds_indiv_t() { } 

  suds_indiv_t( const std::string & id ) : id(id) { } 

  //
  // Analysis commands
  //

  // 'SOAP', i.e. fit model to self and track performance; no new prediction
  void evaluate( edf_t & edf , param_t & param );
	
  // 'REBASE' i.e. call SOAP but change epoch duration after building the model
  void rebase( edf_t & edf , param_t & param , double );
  
  // 'PLACE' i.e. figure out where existing stages should go
  void place( edf_t & edf , param_t & param , const std::string & stagefile );
  
  // 'RESOAP-FILL', i.e. change one stage and refit model (will only be 
  // called on suds_t::cached object, populated after a prior SOAP
  // alter a single epoch 
  void resoap_alter1( edf_t & , int epoch , suds_stage_t stage );

  // pick N of each stage at random
  void resoap_pickN( edf_t & , int );

  // re-fit actual model after making above types of changes
  void resoap( edf_t & , bool ); 

  // update 'ambiguous' epochs (given a confidence threshold)
  // using SOAP; assumes U matrix etc already calculated, i.e.
  // as this is only called from SUDS
  
  int resoap_update_pp( std::vector<std::string> * , 			
			Eigen::MatrixXd * pp ,
			const std::vector<std::string> & labels , 
			const bool global_mode );
  
  //
  // wrapper to add a trainer and save to the libary
  //

  void add_trainer( edf_t & edf , param_t & param );
  
  //
  // wrapper to add as a model (i.e. not w/ data, so no re-preds) 
  //
  
  void add_fit( const std::string & fit );
  
  //
  // main driver: process either trainer or target; (wrapper)
  //
  
  int proc( edf_t & edf , param_t & param , bool trainer = false );
  
  //
  // modules called by proc()
  //
  
  int proc_check_channels( suds_helper_t * );
  
  int proc_extract_observed_stages( suds_helper_t * );
  
  int proc_build_feature_matrix( suds_helper_t * );

  int proc_initial_svd_and_qc( suds_helper_t * );
  
  int proc_main_svd( suds_helper_t * );
  
  int proc_prune_cols( suds_helper_t * );
  
  int proc_prune_rows( suds_helper_t * );
  
  int proc_class_labels( suds_helper_t * );

  int proc_coda( suds_helper_t * );

  void report_epoch_counts( const std::string & l = "" );
  
  //
  // write trainers to file (text only) 
  //
  
  void write( edf_t & edf , param_t & param ) const;
  
  //
  // read/write helpers
  //
  
  inline static void bwrite( std::ofstream & O , const std::string & s );
  inline static void bwrite( std::ofstream & O , int );
  inline static void bwrite( std::ofstream & O , double );

  inline static std::string bread_str( std::ifstream & );
  inline static int bread_int( std::ifstream & );
  inline static double bread_dbl( std::ifstream & );

  // to jump over X
  inline static void bskip_dbl( std::ifstream & , const int n);
				
  //
  // fit LDA/QDA, i.e. after reloading U
  //

  void fit_qlda();


  //
  // make predictions given a different individual's signal data
  //
  
  posteriors_t predict( const suds_indiv_t & trainer , const bool use_qda , double * cancor_u = NULL , double * cancor_vw = NULL );
  
  //
  // add a prediction from one trainer
  //

  void add( const std::string & id , const posteriors_t & , double * a = NULL , double * b = NULL );  
  
  
  //
  // self-classify / run SOAP (which epochs are not self-predicted?)
  //

  int self_classify( std::vector<bool> *  , Eigen::MatrixXd *  pp = NULL );
  
  //
  // summarize stage durations (based on predictions, and note
  // discordance w/ obs, if present)
  // return number of excluded ('bad') epochs
  //

  int summarize_stage_durations( const Eigen::MatrixXd & , const std::vector<std::string> & , int , double );


  //
  // summarize stage durations (based on predictions, and note
  // discordance w/ obs, if present)
  //

  void summarize_epochs( const Eigen::MatrixXd & , 
			 const std::vector<std::string> & , 
			 int , edf_t & );


  //
  // in context of REBASE typically, given short-epoch summaries around (observed)
  // stage transitions: i.e. more detailed look around stage transitions, w/ 4 or 5 second
  // epochs
  //

  void summarize_transitions( const Eigen::MatrixXd & pp , 
			      const std::vector<std::string> & labels, 
			      const int show_left, const int req_left ,
			      const int show_right, const int req_right ,
			      const double epoch_sec , 
			      const int ne_all ,
			      edf_t & edf , param_t & param );

  
  //
  // add (internally) as an annotation
  //
  
  void add_annots( const Eigen::MatrixXd & , const std::vector<std::string> & , int , edf_t & );
		   
  
  //
  // output obs vs prd kappas (5 and 3 level)
  //

  void summarize_kappa( const std::vector<std::string> & prd , const bool to_console = false );

  //
  // stage & context specific accuracies
  //

  void summarize_acc( const std::vector<std::string> & prd );
  
  //
  // dump epoch by predictor matrix (from SOAP verbose)
  //
  
  void dump_predictor_matrix( edf_t & , const std::string & filename );

  //
  // helper to dump predictor / feature / component associations w/ stage
  //
  
  void dump_stage_associations( const std::string & filename );
  
  //
  // dump trainer matrix
  //
  
  void dump_trainer_epoch_matrix( edf_t & , std::map<trkap_t,std::vector<suds_stage_t> > & , std::map<std::string,double> & wgt , const std::string & filename );
  
  
  //
  // dump SVD components into separate files
  //
  
  void dump_svd( const std::string & froot );
  
  
  //
  // get KL weights across trainers
  //

  Eigen::ArrayXd wgt_kl() const;
  
  //
  // Member variables
  //

  // individual ID
  std::string id;

  // trainer or no?  (i.e. has manual staging?)
  bool trainer;
  
  // number of final epochs
  int nve;

  // nunmber of features
  int nf;
  
  // number of retained components (may be < suds_t::nc)
  // but U, V, W etc will have dimension nc
  int nc; 
  
  // original feature matrix
  //  - created for targets (or when making a trainer)
  //  - when predicting, do not need to load this for trainers
  //    (unless they are weight trainers also)
  Eigen::MatrixXd X;
  
  // SVD
  Eigen::MatrixXd U;  // based on own data
  //Eigen::MatrixXd U_projected; // can be projected into this space
  Eigen::ArrayXd W;
  Eigen::MatrixXd V;

  // Hjorth mean/variance (QC) by channel
  Eigen::Array<double, 1, Eigen::Dynamic> mean_h1, sd_h1;
  Eigen::Array<double, 1, Eigen::Dynamic> mean_h2, sd_h2;
  Eigen::Array<double, 1, Eigen::Dynamic> mean_h3, sd_h3;
  // for targets only, keep epoch level Hjorths (epoch x signal)
  Eigen::MatrixXd h1, h2, h3;

  
  // LDA/QDA
  std::vector<std::string> y;
  lda_model_t lda_model; // calculated on reload for trainers
  qda_model_t qda_model; // calculated on reload for trainers
  
  
  // staging
  std::vector<suds_stage_t> obs_stage;       // always all epochs
  std::vector<suds_stage_t> obs_stage_valid; // will match prd_stage
  std::vector<suds_stage_t> prd_stage;

  std::map<std::string,int> counts;

  // set priors?
  std::vector<double> get_priors( const std::vector<double> & p ) const;

  //
  // retained epochs (for targets only)
  //
  
  std::vector<int> epochs;

  //
  // target predictions/staging (generated by 'add()' ) 
  //
  
  std::map<std::string,Eigen::MatrixXd > target_posteriors;
  std::map<std::string,std::vector<suds_stage_t> > target_predictions;
  std::map<std::string,double > cancor_u;
  std::map<std::string,double > cancor_vw;
    
  bool operator<( const suds_indiv_t & rhs ) const {
    return id < rhs.id;
  }
  
};


struct suds_t { 

  //
  // SUDS model
  //

  static suds_model_t model;
  
  //
  // read trainer data from disk (binary only)
  //
  
  static std::string suds_lib_version;
  
  static void attach_db( const std::string & , bool , bool );
  
  static void attach_lib( const std::string & );
  
  static void attach_db_prefit( const std::string & fitfile );
  
  static void attach_hjorth_limits( const std::string & hjorthfile );
  
  static std::vector<suds_indiv_t*> binary_reload( const std::string & filename , bool load_rawx = false );
  
  // convert from text --> binary format for a library file [ +/- feature matrix ] 
  static void text2binary( const std::string & , const std::string & , const bool );
  
  static void combine_trainers( param_t & param );

  static void score( edf_t & edf , param_t & param );
  
  static suds_indiv_t cached;

  static Eigen::MatrixXd add_time_track( const int nr , const int tt );
  
  // final stage/elapsed sleep model
  static void read_elapsed_stages( const std::string & f );
  
  static void set_options( param_t & param );

  //
  // LDA/QDA
  //

  static bool qda;

  //
  // SUDS parameters, needed to be the same across all individuals
  //

  static int soap_mode;
  
  static bool verbose;

  static bool epoch_lvl_output;

  static bool one_by_one;

  static int nc;
  
  static double spectral_resolution;

  // trim leading/trailing wake in trainers?
  static int trim_wake_epochs;
  
  // slope parameters
  static std::vector<double> slope_range;
  static double slope_th;
  static double slope_epoch_th;

  static bool flat_priors;
  
  static std::vector<double> fixed_priors;

  static bool es_model;

  static std::string es_filename;

  static int ns;

  static int nf;
  
  static int n_stages; // 5 or 3 class problem

  static bool pick3then5;
  
  static int fake_ids;

  static std::string fake_id_root;

  static bool use_fixed_trainer_req;

  static std::vector<int> fixed_trainer_req;
  
  static bool self_classification;

  static double self_classification_prob;

  static double self_classification_kappa;

  // MTM param (not used now)
  static bool use_mtm;
  static double mt_tw;  
  static int mt_nt;

  // Welch param
  static double lwr;
  static double upr;

  static bool use_seg_median;

  // SOAP updates 
  static double soap_update_th;
  static int    soap_update_min_epochs;
  
  static double soap_global_update_th;
  static int    soap_global_update_min_epochs;
  
  // Weighting 
  static bool use_repred_weights;

  static bool use_median_repred_weights;
  
  static bool equal_wgt_in_selected;

  static bool use_mcc;
  
  static bool use_5class_repred;

  static bool use_rem_repred;

  static bool use_kl_weights;
    
  static bool use_soap_weights;

  static bool use_maxpp_weights;

  static double wgt_percentile;
  
  static int wgt_exp;

  static bool wgt_mean_normalize;

  static double wgt_mean_th;

  static bool wgt_flip;

  static bool cheat;

  static std::string single_trainer;
  static std::string single_wtrainer;

  static double denoise_fac;

  static bool standardize_X;

  static bool standardize_U;
  
  static bool robust_standardization;

  static double winsor1; // initial PSD
  static double winsor2; // final PSC
  
  static bool use_best_guess;
  
  static bool ignore_target_priors;

  static int required_epoch_n;
  static int max_epoch_n;
  static int equalize_stages;

  static double required_comp_p;
  static double betwithin_ratio;
  
  static std::vector<double> outlier_ths;
  
  // based on trainer mean/SD (averaged) Hjorth parameters
  static Eigen::ArrayXd hjorth1_lwr95;
  static Eigen::ArrayXd hjorth1_upr95;

  static Eigen::ArrayXd hjorth2_lwr95;
  static Eigen::ArrayXd hjorth2_upr95;
  
  static Eigen::ArrayXd hjorth3_lwr95;
  static Eigen::ArrayXd hjorth3_upr95;

  static double hjorth_outlier_th;
  
  static std::string eannot_file;
  static bool        eannot_ints;
  static std::string eannot_prepend;
  static bool        mem_annot;
  
  static std::string mat_dump_file;

  static bool use_bands;  // no PSC, bands instead 

  static bool cache_target; // for iterative SOAP updates
  
private: 

  // trainer library
  static std::map<std::string,suds_indiv_t*> bank;

  // weight-trainer library
  static std::map<std::string,suds_indiv_t*> wbank;
  
  // precomputed ES model weights:  ES , N1 , N2 , N3 , R , W 
  static Eigen::MatrixXd ES_probs;
  static std::vector<double> ES_mins;

  static Eigen::MatrixXd apply_es_model( const Eigen::MatrixXd & , const std::vector<std::string> & stg );
  
public:

  // clear library
  static void empty_banks()
  {
    // trainers
    std::map<std::string,suds_indiv_t*>::iterator ii = bank.begin();
    while ( ii != bank.end() )
      {
	// delete from bank
	delete ii->second;

	// if wbank also pointed to this person, make sure not to delete twice
	std::map<std::string,suds_indiv_t*>::iterator ww = wbank.find( ii->first );
	if ( ww != wbank.end() ) ww->second = NULL;

	++ii;
      }

    // weight trainers
    ii = wbank.begin();
    while ( ii != wbank.end() )
      {
	if ( ii->second != NULL ) // already deleted from above?
	  delete ii->second;
	++ii;	
      }

    bank.clear();
    
  }


  //
  // Misc helpers 
  //
  
  static void make01( Eigen::MatrixXd & r ) { 
    const int n = r.rows();
    const int ns = r.cols();
    for (int i=0;i<n;i++)
      {
	int m = 0;
	double mx = r(i,0);
	for (int j=1; j<ns ;j++) 
	  if ( r(i,j) > mx ) { mx = r(i,j) ; m = j; } 
	for (int j=0; j<ns; j++) r(i,j) = 0;
	r(i,m) = 1;
      } // next row/epoch

  }


  static double mean_maxpp( const Eigen::MatrixXd & pp ) {
    // mean maxpp 
    const int n = pp.rows();
    Eigen::VectorXd m = Eigen::VectorXd::Zero( n );
    for (int i=0; i<n; i++) m[i] = maxpp( pp.row(i).array() );
    return m.mean();
  }

  static double median_maxpp( const Eigen::MatrixXd & pp ) {
    // mean maxpp 
    const int n = pp.rows();
    Eigen::VectorXd m = Eigen::VectorXd::Zero( n );
    for (int i=0; i<n; i++) m[i] = maxpp( pp.row(i).array() );
    return MiscMath::median( eigen_ops::copy_vector( m ) ) ;
  }

  static double maxpp( const Eigen::ArrayXd & r ) { 
    double mx = r[0];    
    for (int j=1;j<suds_t::n_stages;j++) 
      if ( r[j] > mx ) mx = r[j] ;
    return mx;
  }

  static std::string max_inrow( const Eigen::ArrayXd & r , const std::vector<std::string> & labels ) { 

    const int ns = r.size();
    
    if ( ns != labels.size() )
      Helper::halt( "internal error, max()" );

    // track NR/R/W decision first
    double pp_n1 = 0 , pp_n2 = 0 , pp_n3 = 0 , pp_rem = 0 , pp_wake = 0;
      
    // any : return most likely of 5 classes
    // nr  : calc NR / R /W 
    // int m = 0;    
    // double mx = r[0];
    
    for (int j=0; j<ns; j++) 
      {
	if      ( labels[j] == "N2" ) pp_n2 = r[j];
	else if ( labels[j] == "R" ) pp_rem = r[j];
	else if ( labels[j] == "W" ) pp_wake = r[j];
	else if ( labels[j] == "N1" ) pp_n1 = r[j];
	else if ( labels[j] == "N3" ) pp_n3 = r[j];
      }
    
    // NR most likely?
    double pp_nr = pp_n1 + pp_n2 + pp_n3;
    if ( pp_nr > pp_rem && pp_nr > pp_wake )
      {
	if ( pp_n1 >= pp_n2 && pp_n1 >= pp_n3 ) return "N1";
	if ( pp_n2 >= pp_n1 && pp_n2 >= pp_n3 ) return "N2";
	return "N3";
      }

    return pp_rem > pp_wake ? "R" : "W" ; 

  }


  static std::vector<std::string> max( const Eigen::MatrixXd & r , const std::vector<std::string> & labels ) {

    const int ne = r.rows();

    std::vector<std::string> p( ne );

    for (int i=0; i<ne; i++) 
      p[i] = max_inrow( r.row(i) , labels );

    return p;
  }


  // either 5 or 3 stage fixed labels
  static std::vector<std::string> labels;
  static std::vector<std::string> labels3;
  static std::vector<std::string> labels5;
  static std::vector<std::string> labelsR;  
  
  static int num( const std::string & ss ) {
    if ( suds_t::n_stages == 5 )
      {
	if ( ss == "N1" ) return -1;
	if ( ss == "N2" ) return -2;
	if ( ss == "N3" ) return -3;
	if ( ss == "R" ) return 0;
	if ( ss == "W" ) return 1;
	return 2; // unknown/bad/missing
      }

    // 3-stage classification
    if ( ss == "NR" ) return -1;
    if ( ss == "R" ) return 0;
    if ( ss == "W" ) return 1;
    return 2; // unknown/bad/missing
    
  }


  // downcast (but handle case where ss is already downcast / 3-stage too)
  static std::string NRW( const std::string & ss ) {     
    if ( ss == "R" ) return "R";
    if ( ss == "N1" || ss == "N2" || ss == "N3" || ss == "NR" ) return "NR";
    if ( ss == "?" ) return "?";
    return "W";
  }

  static std::vector<std::string> NRW( const std::vector<std::string> & ss ) { 
    std::vector<std::string> s( ss.size() );
    for (int i=0;i<ss.size();i++) s[i] = NRW( ss[i] );
    return s;

  }

  // for repred: reduce to just REM versus not , for example
  static std::string Rnot( const std::string & ss ) {     
    if ( ss == "R" ) return "R";
    if ( ss == "?" ) return "?";
    return "NOT";
  }
  
  static std::vector<std::string> Rnot( const std::vector<std::string> & ss ) { 
    std::vector<std::string> s( ss.size() );
    for (int i=0;i<ss.size();i++) s[i] = Rnot( ss[i] );
    return s;

  }


  static std::vector<suds_stage_t> type( const std::vector<std::string> & s )
  {
    std::vector<suds_stage_t> pp( s.size() );
    for (int i=0;i<s.size();i++) pp[i] = type( s[i] );
    return pp;
  }
  
  static std::vector<std::string> str( const std::vector<suds_stage_t> & s )
  {
    std::vector<std::string> pp( s.size() );
    for (int i=0;i<s.size();i++) pp[i] = str( s[i] );
    return pp;
  }

  static std::string str( const suds_stage_t & s )
  {
    if ( s == SUDS_WAKE ) return "W";
    if ( s == SUDS_N1 ) return "N1";
    if ( s == SUDS_N2 ) return "N2";
    if ( s == SUDS_N3 ) return "N3";
    if ( s == SUDS_NR ) return "NR";
    if ( s == SUDS_REM ) return "R";
    if ( s == SUDS_ARTIFACT ) return "BAD";
    if ( s == SUDS_UNKNOWN ) return "?";
    if ( s == SUDS_UNKNOWN ) return "L";       
    return "?";
  }
  
  static suds_stage_t type( const std::string & s )
  {
    if ( s == "W" ) return SUDS_WAKE;
    if ( s == "N1" ) return SUDS_N1;
    if ( s == "N2" ) return SUDS_N2;
    if ( s == "N3" ) return SUDS_N3;
    if ( s == "NR" ) return SUDS_NR;
    if ( s == "R" ) return SUDS_REM;
    if ( s == "BAD" ) return SUDS_ARTIFACT;
    if ( s == "?" ) return SUDS_UNKNOWN;
    if ( s == "L" ) return SUDS_LIGHTS;
    return SUDS_UNKNOWN;
  }

  
  static std::map<std::string,std::map<std::string,int> > tabulate( const std::vector<std::string> & a , 
								    const std::vector<std::string> & b , 
								    const bool print = false );

  static void trainer_1x1_evals( const suds_indiv_t & , 
				 const Eigen::ArrayXd & wgt , 
				 const std::vector<std::string> & );
  
  static std::pair<double,int> context_acc_stats( const std::vector<int> & obs_ ,
						  const std::vector<int> & pred_ ,
						  const std::vector<int> & epochs_ ,
						  const int type ,
						  const int ostage );
    
  
  
};


#endif
