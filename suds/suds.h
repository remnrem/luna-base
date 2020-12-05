
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
#include "stats/matrix.h"
#include <vector>
#include <set>
#include <map>
#include "stats/lda.h"

#include "eval.h"
#include "helper/helper.h"

struct edf_t;

struct param_t;


enum suds_stage_t
  {
   SUDS_WAKE = 0 ,
   SUDS_N1 = 1 ,
   SUDS_N2 = 2 ,
   SUDS_N3 = 3 ,
   SUDS_NR = 4 , // generic NR (for 3-stage model)
   SUDS_REM = 5 ,
   SUDS_ARTIFACT = 6 , 
   SUDS_UNKNOWN = 7 
  };


struct suds_indiv_t {

  suds_indiv_t() { } 

  suds_indiv_t( const std::string & id ) : id(id) { } 

  // 'SELF-SUDS", i.e. fit model to self and track performance; no new prediction
  void evaluate( edf_t & edf , param_t & param );
		
  // wrapper to add a trainer and save to the libary
  void add_trainer( edf_t & edf , param_t & param );

  // process either trainer or target;
  int proc( edf_t & edf , param_t & param , bool trainer = false );
  
  // write trainers to file
  void write( edf_t & edf , param_t & param ) const;

  // read trainer data from disk
  void reload( const std::string & filename , bool load_rawx = false );

  // fit LDA, io.e. after reloading U
  void fit_lda();

  // make predictions given a different individuals signal data
  lda_posteriors_t predict( const suds_indiv_t & trainer );

  // add a prediction from one trainer
  void add( const std::string & id , const lda_posteriors_t & );

  // self-classify (which epochs are not self-predicted?)
  std::vector<bool> self_classify( int * count , const bool verbose = false , Data::Matrix<double> * pp = NULL );

  // summarize stage durations (based on predictions, and note
  // discordance w/ obs, if present)
  void summarize_stage_durations( const Data::Matrix<double> & , const std::vector<std::string> & , int , double );
  
  // summarize stage durations (based on predictions, and note
  // discordance w/ obs, if present)
  void summarize_epochs( const Data::Matrix<double> & , const std::vector<std::string> & , int , edf_t & );

  // output obs vs prd kappas (5 and 3 level)
  void summarize_kappa( const std::vector<std::string> & prd , const bool to_console = false );

  // get KL weights across trainers
  Data::Vector<double> wgt_kl() const;

  // individual ID
  std::string id;

  // trainer or no?  (i.e. has manual staging?)
  bool trainer;
  
  // number of final epochs
  int nve;

  // number of spectral variables
  int nbins;

  // number of retained components (may be < suds_t::nc)
  // but U, V, W etc will have dimension nc
  int nc;
  
  // spectral data: only loaded for 'weight trainers'
  // not needed to be reloaded for standard trainers
  Data::Matrix<double> PSD;
  
  // SVD
  Data::Matrix<double> U;  // based on own data
  Data::Matrix<double> U_projected; // can be projected into this space
  Data::Vector<double> W;
  Data::Matrix<double> V;

  // Hjorth (mean/variance, per signal)
  Data::Vector<double> mean_h2, sd_h2;
  Data::Vector<double> mean_h3, sd_h3;
  // for targets only, keep epoch level Hjorths (epoch x signal)
  Data::Matrix<double> h2, h3;

  // LDA
  std::vector<std::string> y;
  lda_model_t model; // calculated on reload for trainers
  
  // staging
  std::vector<suds_stage_t> obs_stage;       // always all epochs
  std::vector<suds_stage_t> obs_stage_valid; // will match prd_stage
  std::vector<suds_stage_t> prd_stage;
  
  std::map<std::string,int> counts;

  //
  // retained epochs
  //

  std::vector<int> epochs;

  //
  // target predictions/staging
  //
  
  std::map<std::string,Data::Matrix<double> > target_posteriors;
  std::map<std::string,std::vector<suds_stage_t> > target_predictions;
  
  bool operator<( const suds_indiv_t & rhs ) const {
    return id < rhs.id;
  }
  
};


struct suds_t { 

  //  friend struct suds_indiv_t;
    
  static void attach_db( const std::string & , bool );
  
  static void score( edf_t & edf , param_t & param );


  static void set_options( param_t & param )
  {

    // total/max number of PSC components
    nc = param.has( "nc" ) ? param.requires_int( "nc" ) : 20 ;

    // flat priors?
    flat_priors = param.has( "flat-priors" );

    // bands instead of PSC? 
    use_bands = param.has( "bands" );
    
    // require p<T for each component, in oneway ANOVA a/ stage ( default = 1 );
    required_comp_p = param.has( "pc" ) ? param.requires_dbl( "pc" ) : 0.01;
    if ( param.has( "all-c" ) ) required_comp_p = 99;
    
    // smoothing factor (multiple of SD)
    denoise_fac = param.has( "lambda" ) ? param.requires_dbl( "lambda" ) : 2 ; 

    // epoch-level outlier removal for trainers
    if ( param.has( "th" ) ) outlier_ths = param.dblvector( "th" );

    // self-classification threshold? only use trainer epochs where self-classification == observed score
    self_classification = param.has( "self" ); // epoch-level pruning

    // require posterior > threshold rather than most likely call
    self_classification_prob = param.has( "self-prob" ) ? param.requires_dbl( "self-prob" ) : 99;
    
    // require this individual level kappa for self-prediction to keep a trainer
    self_classification_kappa = param.has( "self-kappa" ) ? param.requires_dbl( "self-kappa" ) : 0;
    
    // must be within X SD units of trainer Hjorth distribution for target epoch to be included 
    hjorth_outlier_th = param.has( "th-hjorth" ) ? param.requires_dbl( "th-hjorth" ) : 3 ;
    
    standardize_psd  = param.has( "norm-psd" ) ? Helper::yesno( param.value( "norm-psd" ) ) : true ;
    standardize_psc  = param.has( "norm-psc" ) ? Helper::yesno( param.value( "norm-psc" ) ) : false ;

    // how to combine predictions across trainers?  by default, use weighted PP (rather than set max to 1 vs 0)
    use_best_guess = param.has( "best-guess" ) ? Helper::yesno( param.value( "best-guess" ) ) : false;

    // if target staging present, ignore ; e.g. if it is all 'UNKNOWN'
    ignore_target_priors = param.has( "ignore-prior" );
    
    // weights: take top N % (if 0 use all) based on the weighting
    wgt_percentile = param.has( "pct" ) ? param.requires_dbl( "pct" ) : 0 ;
    if ( wgt_percentile < 0 || wgt_percentile > 100 ) Helper::halt( "pct should be between 0 and 100" );

    // among selected, (pct) weight equally
    equal_wgt_in_selected = param.has( "equalize-weights" );
    
    // instead of normalizing to 0..1 range, X / mean( X ) 
    wgt_mean_normalize = param.has( "repred-mean" );    

    // threshold on the means (only take above average scorers, default = 1) 
    wgt_mean_th = ( param.has( "repred-mean" ) && param.value( "repred-mean" ) != "T" ) ? param.requires_dbl( "repred-mean" ) : 1 ;

    if ( wgt_mean_normalize && ( wgt_percentile > 0 || equal_wgt_in_selected  ) ) 
      Helper::halt( "cannot specify pct and/or equalize-weights and repred-mean together" );

    // exponential on weight
    wgt_exp = param.has( "wgt-exp" ) ? param.requires_int( "wgt-exp" ) : 0 ;
    if ( wgt_exp < 0 ) Helper::halt( "wgt-exp must be a positive integer" );

    if ( param.has( "wgt-exp" ) && ( wgt_mean_normalize || ( wgt_percentile > 0 || equal_wgt_in_selected  ) ) )
      Helper::halt( "cannot specify wgt-exp along with pct, or repred-mean" );
    
    // wgt1: (do not) use backskip re-weighting
    use_repred_weights = param.has( "repred-weights" ) ? Helper::yesno( param.value( "repred-weights" ) ) : true ;

    // use MCC instead of kappa for weighting
    use_mcc = param.has( "mcc" );
    
    // repred uses 5 classes, not NRW
    use_5class_repred = param.has( "repred5" );
    
    // repred uses only REM vs non-REM
    use_rem_repred = param.has( "repredR" ) ;

    // wgt2: use kl_weights
    use_kl_weights = param.has( "kl-weights" ) ? Helper::yesno( param.value( "kl-weights" ) ) : false ;

    // total weight is either wgt1, wgt2 or a simple average of both wgt1+wgt2

    // allow cheating (trainer is target)
    cheat = param.has( "cheat" );
    
    // NR/R/W or N1/N2/N3/R/W classification?
    n_stages = param.has( "3-stage" ) ? 3 : 5;

    labels5 = { "N1" , "N2" , "N3" , "REM" , "W" };
    labels3 = { "NR" , "R" , "W" };
    labelsR = { "R" , "NOT" }; // just for repred special case

    if ( n_stages == 3 )
      labels = labels3;
    else
      labels = labels5;
    
    // by default, requires at least 5 of each 5 epochs to include a trainer
    required_epoch_n = 5;
    if ( param.has( "req-epochs" ) ) required_epoch_n = param.requires_int( "req-epochs" );

    // select /exactly/ N epochs (at random) from each stage/trainer
    // N1, N2, N3, REM, W (or 3 stages if with have '3-stage'
    use_fixed_trainer_req = param.has( "fixed-epochs" );
    if ( use_fixed_trainer_req )
      {
	fixed_trainer_req = param.intvector( "fixed-epochs" );
	if ( fixed_trainer_req.size() != n_stages )
	  Helper::halt( "requiring fixed-epochs=N1,N2,N3,R,W, or NR/R/W epoch counts" );
      }
    
    // swap in fake IDs to SUDS db:  if ids=suds, then suds_0001, suds_0002, suds_0003 etc
    fake_ids = param.has( "ids" ) ? 1 : 0 ;  // rather than bool, use 1, 2, etc as the trainer count
    fake_id_root = fake_ids ? param.value( "ids" ) : "";
    
    
    //
    // channels, w/ sample rates
    //
    // sig=A,B,C
    // lwr=10,10,10
    // upr=20,20,20
    // inc=0.25,0.25,0.25
    // sr=100,100,100
    
    verbose = param.has( "verbose" );
  
    // ultra verbose -- add each TRAINER one at a time and report kappa/MCC
    one_by_one = param.has( "1x1" );

    epoch_lvl_output = param.has( "epoch" );
    
    if ( param.requires( "sig" ) == "*" ) Helper::halt( "requires sig to be set explicitly" );

    siglab = param.strvector( "sig" );
    
    ns = siglab.size();
    
    //
    // channel-specific options can be given
    //
    
    lwr.resize( ns , 0.5 );
    upr.resize( ns , 20 );
    sr.resize( ns , 100 );
    
    if ( param.has( "lwr" ) )
      {
	lwr = param.dblvector( "lwr" );
	if ( lwr.size() != ns ) Helper::halt( "incorrect number of values for lwr" );
      }
    
    if ( param.has( "upr" ) )
      {
	upr = param.dblvector( "upr" );
	if ( upr.size() != ns ) Helper::halt( "incorrect number of values for upr" );
      }
        
    if ( param.has( "sr" ) )
      {
	sr = param.intvector( "sr" );
	if ( sr.size() != ns ) Helper::halt( "incorrect number of values for sr" );
      }
    

    eannot_file = param.has( "eannot" ) ? param.value( "eannot" ) : "" ;

    eannot_ints = param.has( "stage-numbers" );
    
    eannot_prepend = param.has( "prefix" ) ? ( param.value( "prefix" ) + "_" ) : "" ;
    
    mat_dump_file = param.has( "mat" ) ? param.value( "mat" ) : "" ;

  }

  
  //
  // SUDS parameters, needed to be the same across all individuals
  //

  static bool verbose;

  static bool epoch_lvl_output;

  static bool one_by_one;

  static int nc;

  static bool flat_priors;

  static int ns;

  static int n_stages; // 5 or 3 class problem
  
  static int fake_ids;

  static std::string fake_id_root;

  static bool use_fixed_trainer_req;

  static std::vector<int> fixed_trainer_req;
  
  static bool self_classification;

  static double self_classification_prob;

  static double self_classification_kappa;
  
  static std::vector<std::string> siglab;

  static std::vector<double> lwr;

  static std::vector<double> upr;

  static std::vector<int> sr;

  static bool use_repred_weights;

  static bool equal_wgt_in_selected;

  static bool use_mcc;
  
  static bool use_5class_repred;

  static bool use_rem_repred;

  static bool use_kl_weights;
    
  static double wgt_percentile;
  
  static int wgt_exp;

  static bool wgt_mean_normalize;

  static double wgt_mean_th;

  static bool cheat;

  static double denoise_fac;

  static bool standardize_psd;

  static bool standardize_psc;
  
  static bool use_best_guess;
  
  static bool ignore_target_priors;

  static int required_epoch_n;

  static double required_comp_p;
  
  static std::vector<double> outlier_ths;

  // based on trainer mean/SD (averaged), per signal
  static std::vector<double> lwr_h2, upr_h2;
  static std::vector<double> lwr_h3, upr_h3;
  static double hjorth_outlier_th;
  
  static std::string eannot_file;
  static bool        eannot_ints;
  static std::string eannot_prepend;

  static std::string mat_dump_file;

  static bool use_bands;  // no PSC, bands instead 
  
private: 


  // trainer library
  static std::map<std::string,suds_indiv_t*> bank;

  // weight-trainer library
  static std::map<std::string,suds_indiv_t*> wbank;

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
  
  static void make01( Data::Matrix<double> & r ) { 
    const int n = r.dim1();
    const int ns = r.dim2();
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


  static double maxpp( const Data::Vector<double> & r ) { 
    double mx = r[0];    
    for (int j=1;j<suds_t::n_stages;j++) 
      if ( r[j] > mx ) mx = r[j] ;
    return mx;
  }

  static std::string max( const Data::Vector<double> & r , const std::vector<std::string> & labels ) { 
    const int ns = r.size();
    if ( ns != labels.size() )
      Helper::halt( "internal error, max()" );
    int m = 0;
    double mx = r[0];

    for (int j=1; j<ns; j++) 
      if ( r[j] > mx ) { mx = r[j] ; m = j; } 

    return labels[m];        
  }


  static std::vector<std::string> max( const Data::Matrix<double> & r , const std::vector<std::string> & labels ) {

    const int ne = r.dim1();

    std::vector<std::string> p( ne );

    for (int i=0; i<ne; i++) 
      p[i] = max( r.row(i) , labels );

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
	if ( ss == "REM" ) return 0;
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
    if ( ss == "REM" || ss == "R" ) return "R";
    if ( ss == "N1" || ss == "N2" || ss == "N3" || ss == "NR" ) return "NR";
    return "W";
  }

  static std::vector<std::string> NRW( const std::vector<std::string> & ss ) { 
    std::vector<std::string> s( ss.size() );
    for (int i=0;i<ss.size();i++) s[i] = NRW( ss[i] );
    return s;

  }

  // for repred: reduce to just REM versus not , for example
  static std::string Rnot( const std::string & ss ) {     
    if ( ss == "REM" || ss == "R" ) return "R";
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
    if ( s == SUDS_REM ) return "REM";
    if ( s == SUDS_ARTIFACT ) return "BAD";
    if ( s == SUDS_UNKNOWN ) return "?";       
    return "?";
  }
  
  static suds_stage_t type( const std::string & s )
  {
    if ( s == "W" ) return SUDS_WAKE;
    if ( s == "N1" ) return SUDS_N1;
    if ( s == "N2" ) return SUDS_N2;
    if ( s == "N3" ) return SUDS_N3;
    if ( s == "NR" ) return SUDS_NR;
    if ( s == "REM" ) return SUDS_REM;
    if ( s == "BAD" ) return SUDS_ARTIFACT;
    if ( s == "?" ) return SUDS_UNKNOWN;
    return SUDS_UNKNOWN;
  }

  
  static std::map<std::string,std::map<std::string,int> > tabulate( const std::vector<std::string> & a , 
								    const std::vector<std::string> & b , 
								    const bool print = false );

  static void trainer_1x1_evals( const suds_indiv_t & , 
				 const Data::Vector<double> & wgt , 
				 const std::vector<std::string> & );
    
};


#endif
