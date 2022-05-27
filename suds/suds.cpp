
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


#include "suds.h"

#include <vector>
#include <map>
#include <set>
#include <iomanip>

#include "helper/helper.h"
#include "helper/logger.h"
#include "db/db.h"

#include "dirent.h"

#include "stats/eigen_ops.h"
#include "stats/lda.h"
#include "stats/statistics.h"
#include "miscmath/miscmath.h"
#include "miscmath/crandom.h"

#include "edf/edf.h"
#include "edf/slice.h"

#include "dsp/resample.h"
#include "fftw/fftwrap.h"
#include "dsp/mtm/mtm.h"
#include "dsp/tv.h"

extern logger_t logger;

extern writer_t writer;


//
// library version
//

std::string suds_t::suds_lib_version = "SUDS1";

//
// Signals
//

int suds_t::ns;

//
// Run modes
//

// 0, 1, 2, 3  = SUDS, SOAP, RESOAP, PLACE
int suds_t::soap_mode = 0; 
bool suds_t::cache_target = false;
suds_indiv_t suds_t::cached;

//
// LDA/QDA
//

bool suds_t::qda = false;

//
// Model specification
//

suds_model_t suds_t::model;
std::map<std::string,suds_feature_t> suds_model_t::lab2ftr;
std::map<suds_feature_t,std::string> suds_model_t::ftr2lab;

//
// training library bank
// optionally, a second weight-trainer library wbank
//

std::map<std::string,suds_indiv_t*> suds_t::bank;
std::map<std::string,suds_indiv_t*> suds_t::wbank;

// number of SVD components to extract
int suds_t::nc = 0;

// numebr of features
int suds_t::nf = 0;

// fixed spectral slope parameters
std::vector<double> suds_t::slope_range{ 30.0 , 45.0 } ;
double suds_t::slope_th  = 3;
double suds_t::slope_epoch_th = 5;

// SOAP-update threshold
// (of final prediction, of indiv trainer predictions)
double suds_t::soap_global_update_th = -1;
double suds_t::soap_update_th = -1;

// +5 : if stage has << 5 high confidence epochs, then leave as is
// -5 : if stage has << 5 high confidence epochs, set all to missing (i.e. drop that stage)
int suds_t::soap_update_min_epochs = 5;
int suds_t::soap_global_update_min_epochs = 5;


//
// PSD parameters (by default Welch)
//

bool suds_t::use_seg_median;
double suds_t::lwr;
double suds_t::upr;
double suds_t::spectral_resolution;

// ignore for now
bool suds_t::use_mtm; 
double suds_t::mt_tw = 15;
int suds_t::mt_nt = 2 * suds_t::mt_tw - 1 ; 

int suds_t::trim_wake_epochs = -1; 

//
// Hjorth 95% confidence intervals (per channel) from trainers
//

Eigen::ArrayXd suds_t::hjorth1_lwr95;
Eigen::ArrayXd suds_t::hjorth1_upr95;

Eigen::ArrayXd suds_t::hjorth2_lwr95;
Eigen::ArrayXd suds_t::hjorth2_upr95;

Eigen::ArrayXd suds_t::hjorth3_lwr95;
Eigen::ArrayXd suds_t::hjorth3_upr95;

double suds_t::hjorth_outlier_th = 5;


//
// Trainer weighting schemes
//

bool suds_t::use_kl_weights;
bool suds_t::use_soap_weights;
bool suds_t::use_maxpp_weights;
bool suds_t::use_repred_weights;
bool suds_t::use_median_repred_weights;
bool suds_t::use_mcc;
bool suds_t::use_5class_repred;
bool suds_t::use_rem_repred;
double suds_t::wgt_percentile;
int suds_t::wgt_exp;
bool suds_t::equal_wgt_in_selected;
bool suds_t::use_best_guess = true;
// misc
bool suds_t::wgt_mean_normalize;
double suds_t::wgt_mean_th;
bool suds_t::wgt_flip;
bool suds_t::cheat;

std::string suds_t::single_trainer = "";
std::string suds_t::single_wtrainer = "";
bool suds_t::ignore_target_priors = false;

//
// Normalization options
//
double suds_t::denoise_fac;
bool suds_t::standardize_X = true;
bool suds_t::standardize_U = false;
bool suds_t::robust_standardization = false;
double suds_t::winsor1 = 0;
double suds_t::winsor2 = 0;
std::vector<double> suds_t::outlier_ths;

//
// Trainer selection
//

int suds_t::required_epoch_n = 5;
int suds_t::max_epoch_n = -1;
int suds_t::equalize_stages = 0;
double suds_t::required_comp_p = 0.05;
double suds_t::betwithin_ratio = -1 ; 
bool suds_t::self_classification = false;
double suds_t::self_classification_prob = 99;
double suds_t::self_classification_kappa = 0;
bool suds_t::use_fixed_trainer_req;
std::vector<int> suds_t::fixed_trainer_req;

//
// Stage labels
//

// 5 -> N1, N2, N3, R, W
// 3 -> NR, R, W (3-stage)
int suds_t::n_stages = 5; 
std::vector<std::string> suds_t::labels;
std::vector<std::string> suds_t::labels3;
std::vector<std::string> suds_t::labels5;
std::vector<std::string> suds_t::labelsR;

// first classify NR/R/W then N1/N2/N3 if NR is most likely
//  i.e. to avoid R < NR being selected as R just reflecting fact that
//       there are more NR classes
bool suds_t::pick3then5;


//
// Trainer IDs
//

int suds_t::fake_ids = 1 ;
std::string suds_t::fake_id_root;

//
// Output options
//

bool suds_t::verbose = false;
bool suds_t::epoch_lvl_output = false;
bool suds_t::one_by_one = false;
std::string suds_t::eannot_file = "";
bool suds_t::eannot_ints = false;
std::string suds_t::eannot_prepend = "";
bool suds_t::mem_annot = true;
std::string suds_t::mat_dump_file = "";


//
// Elapsed sleep prior model
//

bool suds_t::es_model;
std::string suds_t::es_filename;
Eigen::MatrixXd suds_t::ES_probs;
std::vector<double> suds_t::ES_mins;

bool suds_t::flat_priors;
std::vector<double> suds_t::fixed_priors;
bool suds_t::use_bands;


void suds_t::set_options( param_t & param )
{
  
  // LDA vs QDA (default)?
  if ( param.has( "lda" ) )
    qda = ! Helper::yesno( param.value( "lda" ) ) ; 
  else if ( param.has( "qda" ) )
    qda = Helper::yesno( param.value( "qda" ) ) ;
  
  if ( qda ) logger << "  using QDA for primary predictions\n";
  else logger << "  using LDA for primary predictions\n";

  // use SOAP to fix up predictions? of indiv trainers and/or globally (on final prediction set)?
  if ( param.has( "soap1" ) )
    {
      std::vector<double> p = param.dblvector( "soap1" );
      if ( p.size() != 2 ) Helper::halt( "requires soap1=th,n" );
      soap_update_th = p[0];
      soap_update_min_epochs = round( p[1] );
    }

  if ( param.has( "soap" ) )
    {
      std::vector<double> p = param.dblvector( "soap" );
      if ( p.size() != 2 ) Helper::halt( "requires soap=th,n" );
      soap_global_update_th = p[0];
      soap_global_update_min_epochs = round( p[1] );
    }

  // spectral resolution for Welch
  spectral_resolution = param.has( "segment-sec" ) ? 1 / param.requires_dbl( "segment-sec" ) : 0.25;
  
  // flat priors?
  flat_priors = param.has( "flat-priors" );
  
  // trim leading/trailing wake? -1 means no, otherwise keep N epochs before/after first/final sleep 
  if ( param.has( "trim" ) ) trim_wake_epochs = param.requires_int( "trim" );

  // apply final elapsed stage model; if passed '.', then ignore
  es_model = param.has( "es-model" );
  es_filename = es_model ? param.value( "es-model" ) : "." ;
  if ( es_filename == "." ) es_model = false;
  
  // fixed priors?
  fixed_priors.clear();
  if ( param.has( "fixed-priors" ) )
    fixed_priors = param.dblvector( "fixed-priors" ) ;
  
  // store target (for subsequent SOAP updates?)
  cache_target = param.has( "save" );
  
  // require p<T for each component, in oneway ANOVA a/ stage ( default = 1 );
  required_comp_p = param.has( "pc" ) ? param.requires_dbl( "pc" ) : 0.01;
  if ( param.has( "all-c" ) ) required_comp_p = 99;
  
  // within/total variance ratio
  betwithin_ratio = param.has( "within" ) ? param.requires_dbl( "within" ) : -1 ;
  
  // // smoothing factor (multiple of SD)
  // denoise_fac = param.has( "lambda" ) ? param.requires_dbl( "lambda" ) : 2 ; 
  
  // epoch-level outlier removal for trainers
  if ( param.has( "th" ) ) outlier_ths = param.dblvector( "th" );
  
  // require posterior > threshold rather than most likely call
  self_classification_prob = param.has( "self-prob" ) ? param.requires_dbl( "self-prob" ) : 99;
  
  // require this individual level kappa for self-prediction to keep a trainer
  self_classification_kappa = param.has( "self-kappa" ) ? param.requires_dbl( "self-kappa" ) : 0;
  
  // self-classification threshold? only use trainer epochs where self-classification == observed score
  self_classification =  param.has( "self-kappa" ) || param.has( "self-prob" )  ; 

  // must be within X SD units of trainer distribution for target epoch to be included 
  hjorth_outlier_th = param.has( "th-hjorth" ) ? param.requires_dbl( "th-hjorth" ) : 5 ;
  
  standardize_X  = param.has( "norm-X" ) ? Helper::yesno( param.value( "norm-X" ) ) : true ;
  standardize_U  = param.has( "norm-U" ) ? Helper::yesno( param.value( "norm-U" ) ) : false ;
  
  // use a robust scaler
  robust_standardization = param.has( "robust" );
  // and winsorization (winsor quantiles, q and 1-q)
  if ( robust_standardization )
    {
      std::vector<double> ww = param.dblvector( "robust" );
      if ( ww.size() > 2 ) Helper::halt( "robust requires 1 or 2 numeric args" ); 
      winsor1 = ww.size() == 0 ? 0 : ww[0];
      winsor2 = ww.size() == 2 ? ww[1] : 0 ;
      if ( winsor1 < 0 || winsor1 > 1 ) Helper::halt( "robust parameter(s) should be between 0 and (e.g.) 0.1" );
      if ( winsor1 > 0.5 ) winsor1 = 1 - winsor1;
      
      if ( winsor2 < 0 || winsor2 > 1 ) Helper::halt( "robust parameter(s) should be between 0 and (e.g.) 0.1" );
      if ( winsor2 > 0.5 ) winsor2 = 1 - winsor1;
      
    }
  
  // how to combine predictions across trainers?  by default, use weighted PP (rather than set max to 1 vs 0)
  use_best_guess = param.has( "best-guess" ) ? Helper::yesno( param.value( "best-guess" ) ) : false;
  
  // if target staging present, ignore ; e.g. if it is all 'UNKNOWN'
  ignore_target_priors = param.has( "ignore-prior" );
  
  //
  // Weights 
  //
  
  // weights: take top N % (if 0 use all) based on the weighting
  wgt_percentile = param.has( "wgt-pct" ) ? param.requires_dbl( "wgt-pct" ) : 0 ;
  if ( wgt_percentile < 0 || wgt_percentile > 100 ) Helper::halt( "wgt-pct should be between 0 and 100" );
  
  // among selected, (pct) weight equally
  equal_wgt_in_selected = param.has( "wgt-equal" );
  
  // instead of normalizing to 0..1 range, X / mean( X ) 
  wgt_mean_normalize = param.has( "wgt-mean" );    

  // threshold on the means (only take above average scorers, default = 1) 
  wgt_mean_th = ( param.has( "wgt-mean" ) && param.value( "wgt-mean" ) != "T" ) ? param.requires_dbl( "wgt-mean" ) : 1 ;

  if ( wgt_mean_normalize && ( wgt_percentile > 0 || equal_wgt_in_selected  ) ) 
    Helper::halt( "cannot specify wgt-pct and/or wgt-equal and wgt-mean together" );
  
  // exponential on weight
  wgt_exp = param.has( "wgt-exp" ) ? param.requires_int( "wgt-exp" ) : 0 ;
  if ( wgt_exp < 0 ) Helper::halt( "wgt-exp must be a positive integer" );
  
  if ( param.has( "wgt-exp" ) && ( wgt_mean_normalize || ( wgt_percentile > 0 || equal_wgt_in_selected  ) ) )
    Helper::halt( "cannot specify wgt-exp along with pct, or wgt-mean" );
  
  // wgt1: (do not) use backskip re-weighting
  use_repred_weights = param.has( "wgt-repred" ) ? Helper::yesno( param.value( "wgt-repred" ) ) : true ;
  
  // use median (not mean) repred weights for multiple w-trainers
  use_median_repred_weights = param.has( "wgt-repred-median" );
  
  // use MCC instead of kappa for weighting
  use_mcc = param.has( "wgt-mcc" );
    
  // repred uses 5 classes, not NRW
  use_5class_repred = param.has( "wgt-5" );
    
  // repred uses only REM vs non-REM
  use_rem_repred = param.has( "wgt-rem" ) ;

  // wgt2: use kl_weights
  use_kl_weights = param.has( "wgt-kl" ) ? Helper::yesno( param.value( "wgt-kl" ) ) : false ;

  // SOAP-weights? 
  use_soap_weights = param.has( "wgt-soap" ) ? Helper::yesno( param.value( "wgt-soap" ) ) : false ; 
  if ( use_soap_weights ) use_repred_weights = false;
  
  // MAXPP-weights? 
  use_maxpp_weights = param.has( "wgt-maxpp" ) ? Helper::yesno( param.value( "wgt-maxpp" ) ) : false ; 
  if ( use_maxpp_weights ) use_repred_weights = false;
  
  // turn off 3-then-5 classification
  pick3then5 = param.has( "pick-3-5" ) ? Helper::yesno( param.value( "pick-3-5" ) ) : true ;
    
  // allow cheating (trainer is target)
  cheat = param.has( "cheat" );
    
  // only pick 1 trainer of all db (i.e. for testing purposes w/ mat /verbose output    
  single_trainer = param.has( "single-trainer" ) ? param.value( "single-trainer" ) : "" ;
  single_wtrainer = param.has( "single-wtrainer" ) ? param.value( "single-wtrainer" ) : "" ;

  // debug code: flip weights
  wgt_flip = param.has( "wgt-flip" );
    
  // NR/R/W or N1/N2/N3/R/W classification?
  n_stages = param.has( "3-stage" ) ? 3 : 5;

  labels5 = { "N1" , "N2" , "N3" , "R" , "W" };
  labels3 = { "NR" , "R" , "W" };
  labelsR = { "R" , "NOT" }; // just for repred special case
  
  if ( n_stages == 3 )
    labels = labels3;
  else
    labels = labels5;
    
  // by default, requires at least 5 of each 5 epochs to include a trainer
  required_epoch_n = 10;
  if ( param.has( "req-epochs" ) ) required_epoch_n = param.requires_int( "req-epochs" );
  
  max_epoch_n = -1;
  if ( param.has( "max-epochs" ) ) max_epoch_n = param.requires_int( "max-epochs" );

  // enfore the same # of each epoch 
  // (i.e. drop to the min) 
  equalize_stages = param.has( "equalize-stages" ) ? param.requires_int( "equalize-stages" ) : 0 ; 
  
     
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
  if ( ! param.has( "ids" ) )
    fake_ids = 0; // was set to 1 above (but we can't reset here, as is incremented when used)       
  else
    fake_id_root = fake_ids ? param.value( "ids" ) : "";
  
  verbose = param.has( "verbose" );
  
  // use MTM instead of Welch  : sets 'tw' , and # tapers = 2*tw-1 , e.g. mtm=15
  use_mtm = param.has( "mtm" );
  mt_tw = use_mtm ? param.requires_int( "mtm" ) : 0 ;
  mt_nt = 2 * mt_tw - 1 ; 
    
  // ultra verbose -- add each TRAINER one at a time and report kappa/MCC
  one_by_one = param.has( "1x1" );
  
  epoch_lvl_output = param.has( "epoch" );
  
    
  //
  // Welch param
  //
  
  use_seg_median = ! param.has( "seg-mean" );
  
  // default 'lwr' and 'upr' for SPEC (if not specified
  // otherwise in the model)
  
  lwr = param.has( "lwr" ) ? param.requires_dbl( "lwr" ) : 0.5;
  upr = param.has( "upr" ) ? param.requires_dbl( "upr" ) : 45;
  
  
  //
  // Annotation outputs
  //
  
  eannot_file = param.has( "eannot" ) ? param.value( "eannot" ) : "" ;

  // unless 'eannot' specified, always add predicted stages to in-memory annots
  mem_annot = ! param.has( "eannot" ) ;
  
  eannot_ints = param.has( "stage-numbers" );
  
  eannot_prepend = param.has( "prefix" ) ? ( param.value( "prefix" ) + "_" ) : "s" ;
  
  // sW, sN1, sN2, sN3, sR, s?
  
  // this only works in single-trainer mode
  // shows both predictions:    target | trainer
  //                            trainer | target
  mat_dump_file = param.has( "mat" ) ? param.value( "mat" ) : "" ;
  
  if ( mat_dump_file != "" && single_trainer == "" ) 
    Helper::halt( "can only specifiy verbose 'mat' output in single-trainer or soap mode" ); 
  
}



//
// Primary entry-point for MAKE-SUDS
//

void suds_indiv_t::add_trainer( edf_t & edf , param_t & param )
{

  // build a trainer; returns number of 'valid'/usable stages  
  int n_unique_stages = proc( edf , param , true );
 
  // only include recordings that have all five/three stages included  
  if ( n_unique_stages != suds_t::n_stages ) 
    {
      logger << "  only found " << n_unique_stages
	     << " of " << suds_t::n_stages
	     << " stages, so not adding as a trainer\n";
      return;
    }  

  // fit models
  qda_t qda( y , U );
  qda_model = qda.fit( suds_t::flat_priors );
  
  lda_t lda( y , U );
  lda_model = lda.fit( suds_t::flat_priors );
    
  // save to disk (as text format)
  write( edf , param ); 
  
}


//
// fit LDA, i.e. after reloading 
//

void suds_indiv_t::fit_qlda()
{

  if ( suds_t::qda )
    {
      qda_t qda( y , U );     
      qda_model = qda.fit( suds_t::flat_priors );
    }
  else
    {
      lda_t lda( y , U );     
      lda_model = lda.fit( suds_t::flat_priors );
    }
  
}


//
// make predictions given a trainer's data & model
//

posteriors_t suds_indiv_t::predict( const suds_indiv_t & trainer , const bool use_qda , double * cancor_u , double * cancor_vw )
{
  
  //
  // Project target (this) into trainer space:   U_targ = X_targ * V_trainer * D_trainer^{-1} 
  // subsetting to # of columns
  //
  
  // std::cout << " this    : " <<         id << " nc = " <<         nc << " " <<         W.size() << "\n";
  // std::cout << " trainer : " << trainer.id << " nc = " << trainer.nc << " " << trainer.W.size() << "\n";

  Eigen::MatrixXd trainer_DW = Eigen::MatrixXd::Zero( trainer.nc , trainer.nc );  
  
  for (int i=0;i< trainer.nc; i++)
    trainer_DW(i,i) = 1.0 / trainer.W[i];
  
  // target's projected U matrix, given trainer V and W 

  Eigen::MatrixXd U_projected = X * trainer.V * trainer_DW;
  
  //
  // Canonical correlation of U or V between target and trainer ? 
  //

  if ( cancor_u != NULL || cancor_vw != NULL )
    {

      Eigen::IOFormat fmt1( Eigen::StreamPrecision, Eigen::DontAlignCols );

      Eigen::MatrixXd DW = Eigen::MatrixXd::Zero( nc , nc );

      for (int i=0;i< nc; i++)
	DW(i,i) = 1.0 / W[i];
      
      if ( cancor_u != NULL )
	{
	  Eigen::VectorXd CCAU = eigen_ops::canonical_correlation( U , U_projected );
	  *cancor_u = CCAU.mean();
	  //std::cout << "CCAU " << CCAU.sum() << "\t" << CCAU.transpose() << "\n";
	}

      if ( cancor_vw != NULL )
	{
	  Eigen::VectorXd CCA = eigen_ops::canonical_correlation( V * DW , trainer.V * trainer_DW );
	  *cancor_vw = CCA.mean();
	  //std::cout << "CCA " << CCA.sum() << "\t" << CCA.transpose() << "\n";
	}
    }
  
  
  //
  // Normalize PSC?
  //

  if ( suds_t::standardize_U ) 
    {
      if ( suds_t::robust_standardization )
	eigen_ops::robust_scale( U_projected , true , true , suds_t::winsor2 );
      else
	eigen_ops::scale( U_projected , true , true );
    }


  //
  // smooth U (projected) [ nope - this is now only done at feature level ] 
  //
  
  // if ( suds_t::denoise_fac > 0 ) 
  //   for (int j=0; j< trainer.nc; j++)
  //     {
  // 	double sd = suds_t::standardize_U ? 1 : eigen_ops::sdev( U_projected.col(j) );
  // 	double lambda = suds_t::denoise_fac * sd;
  // 	dsptools::TV1D_denoise( U_projected.col(j) , lambda );
  //     }  

  //
  // Verbose output?
  //
  
  if ( suds_t::mat_dump_file != "" ) 
    {
      
      // Target U
      std::string filename = Helper::expand( suds_t::mat_dump_file ) + ".target.U";
      logger << "  writing target's U matrix to " << filename << "\n";      
      std::ofstream OUT4( filename.c_str() , std::ios::out );      
      OUT4 << U << "\n";
      OUT4.close();

      // Projected U
      filename = Helper::expand( suds_t::mat_dump_file ) + ".projected.U";
      logger << "  writing target's projected U matrix to " << filename << "\n";      
      std::ofstream OUT1( filename.c_str() , std::ios::out );      
      OUT1 << U_projected << "\n";
      OUT1.close();
      
      // Target V
      filename = Helper::expand( suds_t::mat_dump_file ) + ".target.V";
      logger << "  writing target's V matrix to " << filename << "\n";      
      std::ofstream OUT2( filename.c_str() , std::ios::out );      
      OUT2 << V << "\n";
      OUT2.close();

      // Trainer V
      filename = Helper::expand( suds_t::mat_dump_file ) + ".trainer.V";
      logger << "  writing trainer's V matrix to " << filename << "\n";      
      std::ofstream OUT3( filename.c_str() , std::ios::out );      
      OUT3 << trainer.V << "\n";
      OUT3.close();

      // Trainer U
      filename = Helper::expand( suds_t::mat_dump_file ) + ".trainer.U";
      logger << "  writing trainer's U matrix to " << filename << "\n";      
      std::ofstream OUT5( filename.c_str() , std::ios::out );      
      OUT5 << trainer.U << "\n";
      OUT5.close();

      // Trainer sleep stages
      filename = Helper::expand( suds_t::mat_dump_file ) + ".trainer.S";
      logger << "  writing trainer's sleep stages " << filename << "\n";
      std::ofstream OUT6( filename.c_str() , std::ios::out );
      for (int i=0;i<trainer.obs_stage.size();i++)
	OUT6 << suds_t::str( trainer.obs_stage[i] ) << "\n"; 
      OUT6.close();
      
    }

  
  //
  // predict using trainer model
  //
  
  posteriors_t pp;
  
  if ( use_qda )
    pp = posteriors_t( qda_t::predict( trainer.qda_model , U_projected ) );
  else
    pp = posteriors_t( lda_t::predict( trainer.lda_model , U_projected ) );
  
  return pp;
}




//
// Primary scoring routine
//

void suds_t::score( edf_t & edf , param_t & param ) {


  //
  // by this point, bank will be populated with N+ trainers
  //
  

  //
  // create a target 
  //

  suds_indiv_t target( edf.id ) ;
  
  int n_obs = target.proc( edf , param );

  if ( n_obs == 0 ) return;
  

  //
  // Do we have prior staging available for this target?
  //

  bool prior_staging = target.obs_stage.size() != 0 ;

  //
  // for weight training, on use 'self' unless an explicit wdb
  // was specified, in which case use /all/ people in that
  // wdb (even if that is identical to the db)
  //
  
  bool only_self_retrain = use_repred_weights && ! param.has( "wdb" );

  //
  // Does trainer bank contain target?  (if so, plan to skip it)
  //

  bool bank_contains_target = bank.find( target.id ) != bank.end();

  // if 'cheating' (i.e. allow target to be trainer), then say we don't care 
  // about this (i.e. and vectors will be properly +1 sized to accommodate 
  // all trainers)

  if ( suds_t::cheat ) bank_contains_target = false;

  const int bank_size = bank_contains_target ? bank.size() - 1 : bank.size() ; 
    
  //
  // save weights for each trainer, based on re-predicting
  //

  // either kappa or MCC (if 'mcc' option)
  Eigen::ArrayXd wgt_mean = Eigen::ArrayXd::Zero( bank_size ) ;
  Eigen::ArrayXd wgt_median = Eigen::ArrayXd::Zero( bank_size ) ;

  // kappa; not used in any calcs
  Eigen::ArrayXd wgt_max = Eigen::ArrayXd::Zero( bank_size ) ;
  Eigen::ArrayXd wgt_n50 = Eigen::ArrayXd::Zero( bank_size ) ;

  // soap weights
  Eigen::ArrayXd wgt_soap = Eigen::ArrayXd::Zero( bank_size ) ;

  // maxpp weights
  Eigen::ArrayXd wgt_maxpp = Eigen::ArrayXd::Zero( bank_size ) ;

  //
  // Store Kappa3 for each trainer (valid w/ prior staging only)
  //
  
  Eigen::ArrayXd k3_prior = Eigen::ArrayXd::Zero( bank_size ) ;
  

  //
  // Stats on trainers
  //

  std::map<std::string,double> nr_trainer; // number of unique imputed stages 
  std::map<std::string,std::map<std::string,double> > stg_cnt_trainer;  // cnt of each SS for target given this trainer

  //
  // Stats on weight trainers
  //

  std::map<std::string,double> wtrainer_mean_k3;
  std::map<std::string,int> wtrainer_count_k3;
  
  bool w0 = use_soap_weights;
  bool w1 = wbank.size() > 0 && use_repred_weights;
  bool w2 = use_kl_weights;
  bool w3 = use_maxpp_weights;
  
  if ( w0 && w1 ) Helper::halt( "cannot use both SOAP-weights and repred-weights\n" );
  if      ( w1 && w2 ) logger << "  using mean of repred-weights & KL-weights\n";
  else if ( w0 && w2 ) logger << "  using mean of SOAP-weights & KL-weights\n";
  else if ( w1 ) logger << "  using repred-weights only\n";
  else if ( w2 ) logger << "  using KL-weights only\n";
  else if ( w0 ) logger << "  using SOAP-weights only\n";
  else if ( w3 ) logger << "  using MAXPP-weights only\n";
  else logger << "  not applying any weights\n";


  //
  // dump trainer predictions
  //

  
  std::map<trkap_t,std::vector<suds_stage_t> > alltpreds;
  
  const bool dump_trainer_preds = param.has( "dump-preds" );
  
  //
  // iterate over trainers
  //
  
  std::map<std::string,suds_indiv_t*>::const_iterator tt = bank.begin();

  int cntr = 0;

  while ( tt != bank.end() )
    {
      
      if ( (cntr+1) % 50 == 0 ) logger << "   ... " << (cntr+1) << "/" << bank.size() << " trainers\n";

      //
      // Extract this one trainer
      //

      const suds_indiv_t * trainer = tt->second;

      
      //
      // Skip self? [ nb. does not increate cntr ] 
      //
    
      if ( trainer->id == target.id && ! suds_t::cheat ) { ++tt; continue; } 
      
      //
      // Primary prediction call here.  i.e. predict target given
      // trainer, after projecting target X into the trainer-defined
      // space ( i.e. this generates target.U_projected based on
      // trainer, and then uses it to predict target classes, given
      // the trainer model )
      //

      double cancor_u = 0 , cancor_vw = 0; 
      
      posteriors_t prediction = target.predict( *trainer , suds_t::qda , &cancor_u , &cancor_vw );

      
      //
      // Update these predictions by RESOAP-ing on unambiguous epochs?
      // (probably not a good idea!)
      //

      if ( suds_t::soap_update_th > 0 )
	{
	  int changed = target.resoap_update_pp( &prediction.cl , 
						 &prediction.pp ,
						 suds_t::labels , 
						 false  ); // F --> not in global mode      
	  std::cerr << "  TRAINER ...  changed " << changed << " epochs\n";      
	}
      
      
      //
      // Save predictions
      //
           
      target.add( trainer->id , prediction , &cancor_u , &cancor_vw );
      

      // we can likely remove/change this next step: prediction.cl is
      // stored above for this trainer, but use this slot as a temp so
      // we can make the kappa comparison below; also it is used to
      // get the final epoch number below when making the final set;
      // we should go through and tidy this usage up, as this probably
      // not necessary to store twice (but no biggie)

      target.prd_stage = suds_t::type( prediction.cl );   

      
      //
      // Reweighting (using individuals specified in wbank, if any) 
      //
      // Consider that target's predicted stages (from this one particular trainer)
      // are in fact the real/observed stages for this target.    Now, the 'target'
      // stages and target model is used to predict other people (i.e. called 'weight trainers', 
      // and they are effectively targets in this context)
      //
      // Requires at least 2 predicted stages (of sufficient N) have been predicted by the trainer
      // before doing this step 
 
      //  trainer --> target                         : using trainer model (to define U)
      //              target ---> weight trainer1    : using target model (to define U)
      //              target ---> weight trainer2
      //              target ---> weight trainer3


      //
      // Now consider how well this predicts all the weight-trainers
      // i.e. where we also have true stage information
      //
      
      double max_kappa = 0;
      double mean_kappa = 0;
      std::vector<double> track_median_kappa;
      int n_kappa50 = 0;
      int n_kappa_all = 0;
      
      std::map<std::string,int> counts;
      for (int i=0;i<prediction.cl.size();i++) 
	counts[ prediction.cl[ i ] ]++;

      int nr= 0 ; 
      std::map<std::string,int>::const_iterator cc = counts.begin();
      while ( cc != counts.end() )
	{
	  if ( cc->second >= suds_t::required_epoch_n ) ++nr;	  
	  // save output for this trainer: stage epoch counts
	  stg_cnt_trainer[ trainer->id ][ cc->first ] += cc->second;
	  
	  ++cc;
	}

      // save for output
      nr_trainer[ trainer->id ] = nr;
      
      

      //
      // If prior staging is available, report on kappa for this single trainer
      //

      double k3;

      if ( prior_staging )
	{
	  // obs_stage for valid epochs only
	  double kappa3 =  MiscMath::kappa( NRW( str( target.prd_stage ) ) , 
					    NRW( str( target.obs_stage_valid ) ) , 
					    suds_t::str( SUDS_UNKNOWN ) );	  
	  
	  k3_prior[ cntr ] = kappa3;
	  
	  // tmp
	  k3 = kappa3;
	  
	}


      //
      // Dump trainer outputs?
      //
      
      if ( dump_trainer_preds )
	alltpreds[ trkap_t( trainer->id , ( prior_staging ? k3 : 0 ) ) ] = target.prd_stage ; 
      

      
      //
      // Single-trainer verbose matrix dump mode: rename output root
      // so we see trainer --> trainer
      //  then     target --> trainer  (always back to same)  w/ .repred tag
      //

      if ( mat_dump_file != "" ) 
	mat_dump_file += ".repred";
    

      //
      // Loop of re-prediction targets
      //
      
      bool okay_to_fit_model = nr > 1;

      if ( okay_to_fit_model )
	{
	  
	  //
	  // Generate model for prediction based on 'dummy' target (imputed) stages
	  // but U is based on the target's own SVD (i.e. not projected into trainer space);  
	  // Thus we use target.U, which is the original for the target, based on their own data
	  //
	  
	  // set target model for use w/ all different weight-trainers
	  
	  if ( false  && suds_t::qda )  // always use LDA  (rather than QDA) for re-prediction, as some cells may be small
	    {
	      qda_t qda( prediction.cl , target.U ) ;
	      target.qda_model = qda.fit( suds_t::flat_priors );	      	      
	    }
	  else
	    {
	      lda_t lda( prediction.cl , target.U ) ;
	      target.lda_model = lda.fit( suds_t::flat_priors );
	    }
	  
	  
	  //
	  // Consider one or more weight-trainers for this trainer
	  //   Generally: P_C|B|A
	  //   Or, only considering self:   P_A|B|A
	  //
	  //   where A = trainer, B = target, C is weight-trainer (may have C == A as above)

	  std::map<std::string,suds_indiv_t*>::iterator ww = wbank.begin();
	  while ( use_repred_weights && ww != wbank.end() )
	    {
	      
	      suds_indiv_t * weight_trainer = ww->second;
	      
	      // only use self-training
	      if ( only_self_retrain )
		{
		  if ( trainer->id != weight_trainer->id ) { ++ww; continue; } 
		}
	      
	      // do not use target as a weight-trainer (unless we are 'cheating' ;-) 
	      if ( weight_trainer->id == target.id && ! suds_t::cheat ) { ++ww; continue; } 
	      
	      // always use LDA
	      const bool use_qda = false;
	      posteriors_t reprediction( weight_trainer->predict( target , use_qda ) );
	      
	      weight_trainer->prd_stage = suds_t::type( reprediction.cl );
	      
	      // obs_stage for predicted/valid epochs only
	      
	      double kappa = 0 ; 
	      if ( use_5class_repred ) 
		kappa = MiscMath::kappa( reprediction.cl , 
					 str( weight_trainer->obs_stage ) , 
					 suds_t::str( SUDS_UNKNOWN )  ) ;
	      else if ( use_rem_repred ) 
		kappa = MiscMath::kappa( Rnot( reprediction.cl ) , 
					 Rnot( str( weight_trainer->obs_stage ) ) , 
					 suds_t::str( SUDS_UNKNOWN )  );
	      else
		kappa = MiscMath::kappa( NRW( reprediction.cl ) , 
					 NRW( str( weight_trainer->obs_stage ) ) , 
					 suds_t::str( SUDS_UNKNOWN )  );
	      
	      // swap in MCC instead of kappa?
	      if ( use_mcc )
		{		  
		  double macro_f1 = 0 , macro_precision = 0 , macro_recall = 0 , acc = 0;
		  double wgt_f1 = 0 , wgt_precision = 0 , wgt_recall = 0 , mcc = 0;
		  std::vector<double> precision, recall, f1;
		  
		  if ( use_5class_repred )
		    acc = MiscMath::accuracy( str( weight_trainer->obs_stage ) , 
					      reprediction.cl , 
					      suds_t::str( SUDS_UNKNOWN ) , 
					      &suds_t::labels5 ,
					      &precision, &recall, &f1,
					      &macro_precision, &macro_recall, &macro_f1 ,
					      &wgt_precision, &wgt_recall, &wgt_f1 , &mcc);
		  else if ( use_rem_repred ) // just accuracy on REM
		    acc = MiscMath::accuracy( Rnot( str( weight_trainer->obs_stage ) ) , 
					      Rnot( reprediction.cl ) , 
					      suds_t::str( SUDS_UNKNOWN ) , 
					      &suds_t::labelsR ,
					      &precision, &recall, &f1,
					      &macro_precision, &macro_recall, &macro_f1 ,
					      &wgt_precision, &wgt_recall, &wgt_f1 , &mcc);
		  else
		    acc = MiscMath::accuracy( NRW( str( weight_trainer->obs_stage ) ) , 
					      NRW( reprediction.cl ) , 
					      suds_t::str( SUDS_UNKNOWN ) , 
					      &suds_t::labels3 ,
					      &precision, &recall, &f1,
					      &macro_precision, &macro_recall, &macro_f1 ,
					      &wgt_precision, &wgt_recall, &wgt_f1 , &mcc);
		  
		  // swap in MCC
		  kappa = mcc;
		}
	      
	      
	      ++n_kappa_all;
	      if ( kappa > 0.5 ) n_kappa50++;
	      if ( kappa > max_kappa ) max_kappa = kappa;
	      mean_kappa +=  kappa  ;
	      track_median_kappa.push_back( kappa );

	      //
	      // Verbose outputs?
	      //

	      if ( suds_t::verbose ) 
		{
		  wtrainer_mean_k3[ weight_trainer->id ] += kappa;
		  wtrainer_count_k3[ weight_trainer->id ]++;
		}
	      
	      //
	      // For single trainer verbose output mode only:
	      //
	      
	      if ( suds_t::single_wtrainer != "" && suds_t::mat_dump_file != "" )
		{
		  // re-predicted wtrainer : PP, predicted class
		  
		  std::string filename = Helper::expand( suds_t::mat_dump_file ) + ".wtrainer.pp";
		  logger << "  writing wtrainer's PP | target matrix to " << filename << "\n";
		  std::ofstream OUT1( filename.c_str() , std::ios::out );

		  // header
		  std::vector<std::string> labels = target.qda_model.labels;
		  if ( labels.size() == 0 ) labels = target.lda_model.labels;		  
		  if ( labels.size() != reprediction.pp.cols() ) 
		    Helper::halt( "internal error" );
		  
		  for (int i=0; i<reprediction.pp.cols(); i++)
		    OUT1 << labels[i] << " ";
		  
		  OUT1 << "\n";
		  OUT1 << reprediction.pp << "\n";
		  OUT1.close();
		  
		  filename = Helper::expand( suds_t::mat_dump_file ) + ".wtrainer.pred";
		  logger << "  writing wtrainer's predicted stages | target matrix to " << filename << "\n";
		  if ( weight_trainer->epochs.size() != reprediction.cl.size() ) 
		    Helper::halt( "internal error" );

		  std::ofstream OUT2( filename.c_str() , std::ios::out );
		  for (int i=0; i<reprediction.cl.size(); i++) OUT2 << weight_trainer->epochs[i] << "\t" << reprediction.cl[i] << "\n";
		  OUT2.close();

		}

	      //
	      // Next weight trainer
	      //
	      
	      ++ww;
	    }
	  
	}
      
     
      //
      // Trainer weights
      //

      if ( use_repred_weights && wbank.size() > 0 && okay_to_fit_model ) 
	{
	  //	  std::cout << " cntr, etc " << cntr << " " << wgt_mean.size() << " " <<  ( mean_kappa ) / (double)n_kappa_all << "\n";
	  wgt_max[ cntr ] = max_kappa;
	  wgt_mean[ cntr ] = ( mean_kappa ) / (double)n_kappa_all ;
	  wgt_median[ cntr ] = track_median_kappa.size() == 1 ? mean_kappa : MiscMath::median( track_median_kappa );
	  wgt_n50[ cntr ] = n_kappa50;
	}

      //
      // Weights are baed on median (over epochs) of the max PP (over stages)
      //  
      if ( use_maxpp_weights ) 
	{	  
	  wgt_maxpp[ cntr ] = suds_t::median_maxpp( prediction.pp ) ;	  
	}
      
      
      //
      // SOAP-based trainer weights
      //
      // i.e. given a prediction for target B from trainer A , P_B|A
      //   just do SOAP procedure, as if these were the real values for the target
      //   (and using the target's own U, not the trainer-projected value)

      if ( use_soap_weights )
	{

	  // from above, we already have the model fit:

	  //  i.e. these lines of code above:
	  // lda_t lda( prediction.cl , target.U ) ;
	  // target.model = lda.fit( suds_t::flat_priors );

	  // following the self-evaluation (SOAP) procedure, we get kappa
	  // as follows:

	  // always use LDA for SOAP

	  posteriors_t prediction1;
	  if ( false && suds_t::qda )
	    prediction1 = posteriors_t( qda_t::predict( target.qda_model , target.U  ) );
	  else
	    prediction1 = posteriors_t( lda_t::predict( target.lda_model , target.U  ) );
	  
	  double kappa1 = 0 ;
	  
	  if ( use_5class_repred )
	    kappa1 = MiscMath::kappa( prediction1.cl , 
				      prediction.cl , 
				      suds_t::str( SUDS_UNKNOWN ) );
	  else if ( use_rem_repred )
	    kappa1 = MiscMath::kappa( Rnot( prediction1.cl ) , 
				      Rnot( prediction.cl ) , 
				      suds_t::str( SUDS_UNKNOWN ) );
	  else 
	    kappa1 = MiscMath::kappa( NRW( prediction1.cl ) , 
				      NRW( prediction.cl ) , 
				      suds_t::str( SUDS_UNKNOWN ) );

	  wgt_soap[ cntr ] = kappa1;
	  
	}

      
      //
      // Next trainer
      //

      ++cntr;
      ++tt;

    }


  
  //
  // Derive weights for each trainer based on KL divergence from trainer stage distribition to the mean
  // over all trainers
  //

  // normalized from 0..1 

  Eigen::ArrayXd wgt_kl;
  if ( use_kl_weights )
    wgt_kl = eigen_ops::unit_scale( target.wgt_kl() );
  
  

  //
  // Output all weights, and generate 'final' weights
  //
  
  Eigen::ArrayXd wgt = Eigen::ArrayXd::Zero( bank_size );
  std::vector<std::string> used_trainers;
  std::map<std::string,double> twgts;
  
  tt = bank.begin();
  cntr = 0;

  while ( tt != bank.end() )
    {
      
      const suds_indiv_t * trainer = tt->second;
      
      // skip if target is in the trainer bank [ i.e. do not inc. cntr ]
      if ( trainer->id == target.id && ! suds_t::cheat ) { ++tt; continue; } 
      
      writer.level( trainer->id , "TRAINER" );

      writer.value( "NS" , nr_trainer[ trainer->id ] );

      int sum = 0;
      for (int j=0;j<labels.size();j++)
	sum += stg_cnt_trainer[ trainer->id ][ labels[j] ] ;


      for (int j=0;j<labels.size();j++)
	writer.value( "N_" + labels[j] , stg_cnt_trainer[ trainer->id ][labels[j] ] / (double) sum ) ;
      
      if ( prior_staging )
	writer.value( "K3" , k3_prior[ cntr ] );
    
      // canconical corrs
      writer.value( "CCA_U" , target.cancor_u[ trainer->id ] );
      writer.value( "CCA_VW" , target.cancor_vw[ trainer->id ] );
      
      // only output final WGT (below)
      if ( 0 )
	{
	  if ( use_kl_weights )
	    writer.value( "WGT_KL"   , wgt_kl[ cntr ] );
	  
	  if ( use_soap_weights )
	    writer.value( "WGT_SOAP"  , wgt_soap[ cntr ] );
	  
	  if ( use_repred_weights && wbank.size() > 0 ) 
	    {
	      writer.value( "WGT_N50"  , wgt_n50[ cntr ] );
	      writer.value( "WGT_MAX"  , wgt_max[ cntr ] );
	      writer.value( "WGT_MED"  , wgt_median[ cntr ] );
	      writer.value( "WGT_MEAN" , wgt_mean[ cntr ] ); // normalized	  
	    }
	}
      
      //
      // define 'final' weight: if weight trainers exist, 
      // using WGT_MEAN, otherwise WGT_KL
      //
      
      bool w0 = use_soap_weights;
      bool w1 = wbank.size() > 0 && use_repred_weights;
      bool w2 = use_kl_weights;
      bool w3 = use_maxpp_weights;
      
      if ( w1 && w2 )
	wgt[ cntr ] = ( ( use_median_repred_weights ? wgt_median[ cntr ] : wgt_mean[ cntr ] ) +  wgt_kl[ cntr ] ) / 2.0 ; 
      else if ( w0 && w2 ) 
	wgt[ cntr ] = ( wgt_soap[ cntr ] +  wgt_kl[ cntr ] ) / 2.0 ; 
      else if ( w1 )
	wgt[ cntr ] = use_median_repred_weights ? wgt_median[ cntr ] : wgt_mean[ cntr ] ;
      else if ( w2 )
	wgt[ cntr ] = wgt_kl[ cntr ] ;
      else if ( w0 )
	wgt[ cntr ] = wgt_soap[ cntr ];
      else if ( w3 )
	wgt[ cntr ] = wgt_maxpp[ cntr ];
      else
	wgt[ cntr ] = 1 ; 

      // TMP KLUDGE
      //wgt [ cntr ] = target.cancor_u[ trainer->id ] ;
      // TMP KLUDGE END
	
      
      if ( dump_trainer_preds )
	twgts[ trainer->id ] = wgt[ cntr ] ;
      
      used_trainers.push_back( trainer->id );
			     
      ++tt;
      ++cntr;
    }
  writer.unlevel( "TRAINER" );


  //
  // Verbose output: Dump individual trainer predictions 
  //

  if ( dump_trainer_preds )
    {
      logger << "  writing epoch-level individual-trainer predictions to " << param.value( "dump-preds" ) << "\n";
      target.dump_trainer_epoch_matrix( edf , alltpreds , twgts, param.value( "dump-preds" ) ) ; 
    }
  
  //
  // Verbose output: mean weight trainer values
  //

  if ( suds_t::verbose && use_repred_weights && wbank.size() > 0 )
    {
      std::map<std::string,suds_indiv_t*>::iterator ww = wbank.begin();
      while ( ww != wbank.end() )
	{	  
	  suds_indiv_t * weight_trainer = ww->second;
	  
	  if ( weight_trainer->id != target.id ) 
	    {
	      writer.level( weight_trainer->id , "WTRAINER" );
	      double m = wtrainer_mean_k3[ weight_trainer->id ] / (double)wtrainer_count_k3[ weight_trainer->id ];
	      writer.value( "K3" , m );
	    }
	  ++ww;
	}
      writer.unlevel( "WTRAINER" );
            
    }


  //
  // Normalize wgt / truncate at percentile?
  //

  bool has_wgt = ( wbank.size() > 0 && use_repred_weights ) || use_kl_weights || use_soap_weights || use_maxpp_weights;
  
  if ( has_wgt && suds_t::wgt_mean_normalize ) 
    {
      logger << "  normalizing weights by the trainer mean\n";
      double mean_wgt = wgt.mean();
      for ( int i=0; i<wgt.size();i++)
	{	  
	  if ( wgt[i] < 0 ) wgt[i] = 0;
	  else wgt[i] /= mean_wgt;       
	  // default threshold = 1 (i.e. only take values above the mean)
	  if ( wgt[i] < suds_t::wgt_mean_th ) wgt[i] = 0;
	}      
    }
  

  //
  // Unit scale exponential 
  //

  if ( has_wgt && wgt_exp > 1 ) 
    {
      // get MAX --> set to 1.0
      double max = 0;
      for (int i=0;i<wgt.size();i++)
	{
	  if ( wgt[i] < 0 ) wgt[i] = 0; 
	  else if ( wgt[i] > max ) max = wgt[i];
	}
      
      for (int i=0;i<wgt.size();i++)
	{
	  wgt[i] /= max;
	  wgt[i] = pow( wgt[i] , wgt_exp );
	}
      
    }

  //
  // Testing only: flip weights?
  //

  if ( suds_t::wgt_flip )
    {
      logger << "  debug code: flipping weights\n";
      for (int i=0;i<wgt.size();i++)
	wgt[i] = 1 - wgt[i] ;
    }

  //
  // Percentile based scaling (subsetting) 
  //

  if ( has_wgt && suds_t::wgt_percentile > 0 ) 
    {
      
      // get value X = top N% and set to 0/1 if below/above X
      std::vector<double> cc = eigen_ops::copy_array( wgt ) ;
      double threshold = MiscMath::percentile( cc , 1.0 - suds_t::wgt_percentile / 100.0 ) ;
      
      // binarize wgt (if only 1 or 2 trainers, then assign equal weight) 
      if ( wgt.size() < 3 || suds_t::equal_wgt_in_selected ) 
	{
	  for (int i=0;i<wgt.size();i++)
	    wgt[i] = wgt[i] >= threshold ? 1 : 0 ; 
	}
      else
	{
	  for (int i=0;i<wgt.size();i++)
	    wgt[i] = wgt[i] >= threshold ? wgt[i] : 0 ;

	  // unit-scale between /threshold/ and max
	  wgt = eigen_ops::unit_scale( wgt , threshold , 1.0 );
	}
      
    }

  //
  // Output final trainer weights
  //

  for (int t=0; t<used_trainers.size(); t++)
    {
      writer.level( used_trainers[t] , "TRAINER" );
      writer.value( "WGT" , wgt[t] );
    }
  writer.unlevel( "TRAINER" );
  
  
  //
  // Construct for reporting epoch-level stats below (final, and optionally per-trainer)
  //
 
  std::map<int,int> e2e;
  for (int i=0; i<target.epochs.size(); i++) e2e[target.epochs[i]] = i ;  
  const int ne_all = edf.timeline.num_epochs();


  //
  // Construct (weighted) posterior probabilities
  //    
  
  const int ne = target.prd_stage.size();

  // target.prd_stage.clear();
  // target.prd_stage.resize( SUDS_UNKNOWN );

  Eigen::MatrixXd pp = Eigen::MatrixXd::Zero( ne , suds_t::n_stages );

  int ntrainers = 0;
  double tot_wgt = 0;
  int tot_unwgt = 0;

  std::map<std::string,Eigen::MatrixXd >::iterator ii = target.target_posteriors.begin();
  while ( ii != target.target_posteriors.end() )
    {
      
      // get posteriors from this trainer 

      Eigen::MatrixXd & m = ii->second;

      // force 0/1 encoding? i.e. 100% weight placed on most likely
      
      if ( suds_t::use_best_guess ) 
	suds_t::make01( m );
      
      double w = wgt[ ntrainers ];

      tot_wgt += w;

      if ( w > 0 ) 
	tot_unwgt++;

      if ( pp.rows() != m.rows() || pp.cols() != m.cols() )
	Helper::halt( "internal error in compiling posteriors across trainers" );

      // accumulate final posterior set
      
      for (int i=0;i<ne;i++)
	for (int j=0;j<suds_t::n_stages;j++)
	  pp(i,j) += w * m(i,j);

      // verbose output
      if ( suds_t::verbose )
	{
	  const suds_indiv_t * trainer = bank.find( ii->first )->second ;	  
	  writer.level( trainer->id , "TRAINER" );
	  
	  for (int i=0;i<ne_all;i++)
	    {
	      int e = -1;
	      if ( e2e.find( i ) != e2e.end() ) e = e2e[i];
	      if ( e != -1 ) 
		{
		  writer.epoch( edf.timeline.display_epoch( i ) );
		  std::string predss1 = max_inrow( m.row(e) , suds_t::labels );
		  writer.value( "PRED" , predss1 );		  
		  
		  double pp_nr = 0;
		  bool has_nr = false;
		  for (int j=0;j<labels.size();j++)
		    {
		      if ( labels[j] == "NR" ) has_nr = true;
		      if ( labels[j] == "N1" || labels[j] == "N2" || labels[j] == "N3" ) pp_nr += m(e,j);
		      writer.value( "PP_" + labels[j] , m(e,j) );
		    }
		  
		  // automatically aggregate N1+N2+N3 under the 5-class model (or whatever NREM stages are present)
		  if ( ! has_nr )
		    writer.value( "PP_NR" , pp_nr );
		  
		}
	    }
	  writer.unepoch();
	}

      // next trainer
      ++ntrainers;
      ++ii;
    }
  
  if ( suds_t::verbose )
    writer.unlevel( "TRAINER" );

  if ( ntrainers == 0 ) 
    Helper::halt( "no valid trainers, quitting" );
  
  if ( has_wgt && suds_t::wgt_percentile > 0 ) 
    logger << "  constructed posteriors using top "
	   << suds_t::wgt_percentile << " percentile, "
	   << (int)tot_unwgt << " (of " << ntrainers << ") trainers (weighted N = " << tot_wgt << ")\n";
  else if ( has_wgt ) 
    logger << "  constructed posteriors using " << ntrainers << " trainers (weighted N = " << tot_wgt << ")\n";
  else
    logger << "  constructed posteriors using " << ntrainers << " trainers\n";


  //
  // Normalize (weighted) posteriors to sum to 1.0, and get MAP
  //

  double mean_maxpp = 0;

  for (int i=0;i<ne;i++)
    {
      // normalize
      for (int j=0;j<suds_t::n_stages;j++) // 5 or 3 stages
        pp(i,j) /= (double)tot_wgt;

      // track level of confidence for MAP
      mean_maxpp += suds_t::maxpp( pp.row(i) );
      
    }
  mean_maxpp /= (double)ne;


  //
  // Revised estimates based on ES model?
  //

  if ( suds_t::es_model )
    {
      std::vector<std::string> current_prediction;
      
      for (int i=0;i<ne_all;i++)
	{
	  int e = -1;
	  if ( e2e.find( i ) != e2e.end() ) e = e2e[i];
	  if ( e != -1 ) 
	    current_prediction.push_back( max_inrow( pp.row(e) , suds_t::labels ) );
	}
  
      // does nothing is ES model is already attached
      suds_t::read_elapsed_stages( es_filename );
      logger << "  applying ES model to revised final predictions\n";
      // update
      pp = suds_t::apply_es_model( pp , current_prediction );
    }
  
  
  //
  // Report epoch-level stats
  //
 
  // std::map<int,int> e2e;
  // for (int i=0; i<target.epochs.size(); i++) e2e[target.epochs[i]] = i ;  
  // const int ne_all = edf.timeline.num_epochs();

  std::vector<std::string> final_prediction;

  for (int i=0;i<ne_all;i++)
    {
      int e = -1;
      if ( e2e.find( i ) != e2e.end() ) e = e2e[i];
      if ( e != -1 ) 
	{
	  // most likely value
	  std::string predss = max_inrow( pp.row(e) , suds_t::labels );
	  //writer.value( "PRED" , predss );
	  final_prediction.push_back( predss );
	}
    }


  //
  // Update final predictions based on SOAP?
  //

  //  target.summarize_kappa( final_prediction , true );

  if ( suds_t::soap_global_update_th > 0 )
    {
      //std::vector<std::string> ss = suds_t::str( target.obs_stage_valid );
      int changed = target.resoap_update_pp( &final_prediction ,
					     &pp ,
					     suds_t::labels , 
					     true  ); // T --> in global mode      
      logger << "  changed " << changed << " epochs\n";      
    }

  
  //
  // All done w/ SUDS ... now output summaries
  //

  
  //
  // Output epoch-level information
  //
  
  target.summarize_epochs( pp , suds_t::labels , ne_all , edf );
  

  //
  // Summarize staging
  //

  const double epoch_sec = edf.timeline.epoch_length();

  const int bad_epochs = target.summarize_stage_durations( pp , suds_t::labels , ne_all , epoch_sec );
  
  writer.value( "BAD_N" , bad_epochs );
  writer.value( "BAD_P" , bad_epochs/(double) ne_all );

  //
  // Confusion matrics and kappa w/ observed staging
  //

  if ( prior_staging )
    {
      
      //
      // report kappa w/ observed 
      //

      target.summarize_kappa( final_prediction , true );
            
      //
      // also, given correlations between weights and trainer kappas
      //
      
      writer.value( "R_WGT" ,
		    Statistics::correlation( eigen_ops::copy_array( wgt ) ,
					     eigen_ops::copy_array( k3_prior) ) ); 
      
    }
  

  
  //
  // Final SOAP evaluation of /predicted/ stages
  //
  
  // following the self-evaluation (SOAP) procedure, we get kappa
  // as follows:
  
  // build model based on predicted stages
  
  // assume putative 'y' and 'U' will have been constructed, and 'nve' set
  // i.e. this will be called after proc(), or from near the end of proc()

  std::set<std::string> nstages;
  for (int i=0; i<final_prediction.size(); i++) 
    nstages.insert( final_prediction[i] );

  if ( nstages.size() > 1 ) 
    {
      bool valid = false;

      // always use LDA w/ SOAP
      if ( false && suds_t::qda )
	{
	  qda_t self_qda( final_prediction , target.U );
	  qda_model_t self_model = self_qda.fit( suds_t::flat_priors );
	  if ( self_model.valid )
	    {
	      // get predictions: SOAP model (fitting to self)
	      qda_posteriors_t soap_final_prediction = qda_t::predict( self_model , target.U );	      
	      double kappa5 = MiscMath::kappa( soap_final_prediction.cl , final_prediction , suds_t::str( SUDS_UNKNOWN ) );
	      double kappa3 = MiscMath::kappa( NRW( soap_final_prediction.cl ) , NRW( final_prediction ) , suds_t::str( SUDS_UNKNOWN ) );	      
	      writer.value( "SOAP" , kappa5 );
	      writer.value( "SOAP3" , kappa3 );
	    }
	}
      else
	{
	  lda_t self_lda( final_prediction , target.U );
	  lda_model_t self_model = self_lda.fit( suds_t::flat_priors );
	  if ( self_model.valid )
	    {
	      // get predictions: SOAP model (fitting to self)
	      posteriors_t soap_final_prediction = lda_t::predict( self_model , target.U );	      
	      double kappa5 = MiscMath::kappa( soap_final_prediction.cl , final_prediction , suds_t::str( SUDS_UNKNOWN ) );
	      double kappa3 = MiscMath::kappa( NRW( soap_final_prediction.cl ) , NRW( final_prediction ) , suds_t::str( SUDS_UNKNOWN ) );	      
	      writer.value( "SOAP" , kappa5 );
	      writer.value( "SOAP3" , kappa3 );
	    }
	}
      
    }
  
  
  //
  // Misc other output
  //

  writer.value( "MAXPP" , mean_maxpp );



  //
  // Verbose 1-by-1 trainer additions
  //

  if ( suds_t::one_by_one )
    {
      if ( ! prior_staging ) Helper::halt( "need prior staging data for 1x1" );

      // create vector of obs for predicted epochs only (i.e. good ones)
      std::vector<std::string> obs;
      for (int i=0;i<ne_all; i++)
        {
          int e = -1;
          if ( e2e.find( i ) != e2e.end() ) e = e2e[i];
          if ( e == -1 ) continue;          
	  obs.push_back( str( target.obs_stage[i] ) );
	}
            
      suds_t::trainer_1x1_evals( target , wgt , obs );
   }


  

  //
  // Verbose output?
  //

  if ( suds_t::mat_dump_file != "" ) 
    {

      // more versbose INFO including staging 
      // nb. use 'orig' name, w/out .repred tag, so remove last 7 characters

      if ( mat_dump_file.substr( mat_dump_file.size() - 7 ) == ".repred" )
	mat_dump_file = mat_dump_file.substr( 0 , mat_dump_file.size() - 7 );
      
      std::string filename = Helper::expand( mat_dump_file ) ;
      std::ofstream OUT1( filename.c_str() , std::ios::out );
      
      logger << "  writing target epoch-wise matrix to " << filename << "\n";
      OUT1 << "ID\tE";

      for (int i=0;i<target.X.cols();i++)
	OUT1 << "\t" << "X" << (i+1);

      for (int i=0;i<target.U.cols();i++)
	OUT1 << "\t" << "U" << (i+1);

      if ( suds_t::n_stages == 5 )
	{
	  OUT1 << "\tPP_N1"
	       << "\tPP_N2"
	       << "\tPP_N3"
	       << "\tPP_R"
	       << "\tPP_W";      
	}
      else
	{
	  OUT1 << "\tPP_NR"
	       << "\tPP_R"
	       << "\tPP_W";      
	}
      
      OUT1 << "\tPRD";

      if ( prior_staging ) OUT1 << "\tOBS";
      OUT1 << "\n";
      
      // each row/epoch
      for (int i=0;i<ne_all; i++)
	{

	  int e = -1;
	  if ( e2e.find( i ) != e2e.end() ) e = e2e[i];	  
	  
	  if ( e == -1 ) continue;

	  // only display good lines
	  OUT1 << target.id << "\t"
	       << edf.timeline.display_epoch( i ) ;
	  
	  for (int j=0;j<target.X.cols();j++)
	    OUT1 << "\t" << target.X(e,j); 
	  
	  for (int j=0;j<target.U.cols();j++)
	    OUT1 << "\t" << target.U(e,j); 
	  
	  for (int j=0;j<pp.cols();j++)
	    OUT1 << "\t" << pp(e,j); 
	  
	  OUT1 << "\t" << final_prediction[e] ;

	  if ( prior_staging ) 
	    OUT1 << "\t" << str( target.obs_stage[i] ) ;

	  OUT1 << "\n";

	}

      OUT1.close();
    }


  //
  // Add stage predictions to in-memory annotation class
  //
  
  if ( suds_t::mem_annot )
    {

      // epochs[] contains the codes of epochs actually present in the model/valid
      std::map<int,int> e2e;
      for (int i=0; i<target.epochs.size(); i++) e2e[target.epochs[i]] = i ;
      const int ne_all = edf.timeline.num_epochs();

      std::set<std::string> opreds;
      
      for (int i=0; i < ne_all; i++)
	{
	  int e = -1;
	  if ( e2e.find( i ) != e2e.end() ) e = e2e[i];
	  
	  // epoch interval
	  interval_t interval = edf.timeline.epoch( i );

	  // no prediction?
	  if ( e == -1 ) 
	    {
	      std::string predss = suds_t::str( SUDS_UNKNOWN );
	      annot_t * a = edf.timeline.annotations.add( predss );
	      instance_t * instance = a->add( "." , interval , "." );
	    }
	  else
	    {
	      // most likely value
	      std::string predss = suds_t::eannot_prepend + suds_t::max_inrow( pp.row(e) , labels );
	      annot_t * a = edf.timeline.annotations.add( predss );
	      instance_t * instance = a->add( "." , interval , "." );
	    }
	}
    }

  
  //
  // Write .eannot file?
  //
  
  if ( suds_t::eannot_file != "" )
    {
      // expecting this will have individual wild-cards, but these will have 
      // been expanded already;  this is for home folder encoding ~ 
      std::string filename = Helper::expand( suds_t::eannot_file );

      logger << "\n  writing .eannot stage annotations "
	     << ( suds_t::eannot_ints ? "(as integeres) " : "" )
	     << " to " << filename << "\n";

      std::ofstream OUT1( filename.c_str() , std::ios::out );


      // make sure we output all epochs
      for (int i=0;i<ne_all;i++)
	{
	  
	  int e = -1;
	  if ( e2e.find( i ) != e2e.end() ) e = e2e[i];
	  if ( e != -1 )
	    {
	      if ( suds_t::eannot_ints )
		OUT1 << suds_t::num( final_prediction[e] ) << "\n";
	      else
		OUT1 << final_prediction[e] << "\n";	      
	    }
	  else // could not score
	    {
	      if ( suds_t::eannot_ints ) OUT1 << suds_t::num( "?" ) << "\n";
	      else OUT1 << "?" << "\n";
	    }
	}
            
      OUT1.close();      
    }
  
}




void suds_indiv_t::add( const std::string & trainer_id , const posteriors_t & prediction ,
			double * cu , double * cvw )
{
  
  target_posteriors[ trainer_id ] = prediction.pp ;
  
  target_predictions[ trainer_id ] = suds_t::type( prediction.cl );

  if ( cu != NULL ) cancor_u[ trainer_id ] = *cu;

  if ( cvw != NULL ) cancor_vw[ trainer_id ] = *cvw;
  
}



Eigen::ArrayXd suds_indiv_t::wgt_kl() const { 


  // returned weights
  const int nt = target_predictions.size();  

  Eigen::ArrayXd W = Eigen::ArrayXd::Zero( nt );

  if ( nt == 0 ) return W;

  Eigen::MatrixXd Q = Eigen::MatrixXd::Zero( nt , suds_t::n_stages ) ;  

  int r = 0;
  std::map<std::string,std::vector<suds_stage_t> >::const_iterator ii = target_predictions.begin();
  while ( ii != target_predictions.end() ) 
    {

      const double ne = ii->second.size();

      if ( suds_t::n_stages == 5 ) 
	for (int e=0; e<ne; e++) 
	  {	    
	    if      ( ii->second[e] == SUDS_N1 ) Q(r,0)++;
	    else if ( ii->second[e] == SUDS_N2 ) Q(r,1)++;
	    else if ( ii->second[e] == SUDS_N3 ) Q(r,2)++;
	    else if ( ii->second[e] == SUDS_REM ) Q(r,3)++;
	    else if ( ii->second[e] == SUDS_WAKE ) Q(r,4)++;
	  }
      else
	for (int e=0; e<ne; e++) 
	  {	    
	    if      ( ii->second[e] == SUDS_NR ) Q(r,0)++;
	    else if ( ii->second[e] == SUDS_REM ) Q(r,1)++;
	    else if ( ii->second[e] == SUDS_WAKE ) Q(r,2)++;
	  }
	

      // normalize
      for (int s=0;s<suds_t::n_stages;s++) Q(r,s) /= ne;
      
      // next trainer
      ++ii;
      ++r;
    }

  
  // Means

  Eigen::ArrayXd P = Q.colwise().mean();
  
  const double KL_EPS = 1e-6;

  // divergence for each trainer from the mean

  r = 0;
  ii = target_predictions.begin();
  while ( ii != target_predictions.end() ) 
    {      
      // negative KL
      double ss = 0;
      for ( int s = 0 ; s < suds_t::n_stages ; s++ )
	if ( Q(r,s) > KL_EPS ) ss += P[s] * log( P[s] / Q(r,s) );  

      W[r] = -ss;
      
      ++r;
      ++ii;
    }
  
  return W;
}





void suds_t::trainer_1x1_evals( const suds_indiv_t & target , 
				const Eigen::ArrayXd & wgt, 
				const std::vector<std::string> & obs_stages )
{
 
  //
  // Get ordering of trainers
  //
  
  struct trainer_ord_t { 
    trainer_ord_t( double w , const std::string & id ) : w(w) , id(id) { } 
    double w;
    std::string id;
    bool operator< ( const trainer_ord_t & rhs ) const 
    {
      // nb. revserved order, to get highest weighted individuals first
      if ( w < rhs.w ) return false;
      if ( w > rhs.w ) return true;
      return id < rhs.id;
    }
  };
  

  std::set<trainer_ord_t> otrainers;

  int ntrainers = 0;

  std::map<std::string,Eigen::MatrixXd >::const_iterator ii = target.target_posteriors.begin();
  while ( ii != target.target_posteriors.end() )
    {
      otrainers.insert( trainer_ord_t( wgt[ ntrainers ] , ii->first ) );
      ++ntrainers;
      ++ii;
    }
  
  
  //
  // Define starting PP matrix 
  //

  const int ne = target.prd_stage.size();

  Eigen::MatrixXd pp = Eigen::MatrixXd::Zero( ne , suds_t::n_stages );

  double cum_wgt = 0; 

  int nt = 0;


  //
  // Construct (weighted) posterior probabilities
  //    
  
  
  std::set<trainer_ord_t>::const_iterator oo = otrainers.begin();
  while ( oo != otrainers.end() )
    {
      // trainer ranking
      ++nt;

      // get posteriors from this trainer
      Eigen::MatrixXd m = target.target_posteriors.find( oo->id )->second;
      
      // force 0/1 encoding? i.e. 100% weight placed on most likely                                                                                                                                  
      if ( suds_t::use_best_guess ) suds_t::make01( m );
      
      // update PP (weighted)
      if ( oo->w > 0 ) 
	{
	  for (int i=0;i<ne;i++)
	    for (int j=0;j<suds_t::n_stages;j++)
	      pp(i,j) += oo->w * m(i,j);
	}

      // don't worry about normalization, as we are just taking the max per row
      
      //
      // Get predicted (most likely) class per epoch
      //

      std::vector<std::string> current_prediction;
      
      for (int i=0;i<ne;i++)
	current_prediction.push_back( max_inrow( pp.row(i) , suds_t::labels ) );	
      
      //
      // Eval
      //

      if ( current_prediction.size() != obs_stages.size() )
	Helper::halt( "internal error w/ 1x1" );
      
      double kappa = MiscMath::kappa( current_prediction , obs_stages , suds_t::str( SUDS_UNKNOWN )  );
      double kappa3 = MiscMath::kappa( suds_t::NRW( current_prediction ) , suds_t::NRW( obs_stages ) , suds_t::str( SUDS_UNKNOWN ) );
      
      //
      // track cumulative weighted N
      //

      cum_wgt += oo->w;

      //
      // Output 
      //
      
      writer.level( nt , "NTRAINER" );
      writer.value( "TRAINER" , oo->id );
      writer.value( "WGT" , oo->w );
      writer.value( "CUM_WGT" , cum_wgt );
      writer.value( "K" , kappa );
      writer.value( "K3" , kappa3 );

      // next best trainer
      
      ++oo;
    }

  writer.unlevel( "NTRAINER" );

}


std::vector<double> suds_indiv_t::get_priors( const std::vector<double> & p ) const
{
  // take N1 N2 N3 R W priors and rescale to whatever categories exist for this person
  std::vector<double> dummy;
  return dummy;
}



Eigen::MatrixXd suds_t::add_time_track( const int nr , const int tt )
{

  if ( nr <= 0 || tt <= 0 ) Helper::halt( "internal error in add_time_track()" );

  Eigen::MatrixXd T = Eigen::MatrixXd::Zero( nr , tt );

  for (int r=0; r<nr; r++) 
    for ( int c=0; c<tt; c++)
      T(r,c) = pow( ( r / (double)nr ) - 0.5 , c+1 );   

  return T;

}


