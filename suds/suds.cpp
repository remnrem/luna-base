
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


// notes:
// https://link.springer.com/article/10.1007/s10618-019-00638-y


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
#include "dsp/tv.h"

extern logger_t logger;

extern writer_t writer;


// 0, 1, 2 = SUDS, SOAP, RESOAP
int suds_t::soap_mode = 0; 
bool suds_t::cache_target = false;
suds_indiv_t suds_t::cached;

bool suds_t::copy_db_mode = false;
bool suds_t::verbose = false;
bool suds_t::epoch_lvl_output = false;
bool suds_t::one_by_one = false;
std::map<std::string,suds_indiv_t*> suds_t::bank;
std::map<std::string,suds_indiv_t*> suds_t::wbank;
int suds_t::nc;
int suds_t::ns;
bool suds_t::flat_priors;
bool suds_t::use_bands;
std::vector<std::string> suds_t::siglab;
std::vector<double> suds_t::lwr;
std::vector<double> suds_t::upr;
std::vector<int> suds_t::sr;

bool suds_t::use_kl_weights;
bool suds_t::use_soap_weights;
bool suds_t::use_repred_weights;
bool suds_t::use_mcc;
bool suds_t::use_5class_repred;
bool suds_t::use_rem_repred;
double suds_t::wgt_percentile;
int suds_t::wgt_exp;
bool suds_t::equal_wgt_in_selected;
bool suds_t::wgt_mean_normalize;
double suds_t::wgt_mean_th;

bool suds_t::cheat;
std::string suds_t::single_trainer = "";
double suds_t::denoise_fac;
bool suds_t::standardize_psd = true;
bool suds_t::standardize_psc = false;
bool suds_t::robust_standardization = false;
double suds_t::winsor1 = 0;
double suds_t::winsor2 = 0;
bool suds_t::use_best_guess = true;
bool suds_t::ignore_target_priors = false;
std::vector<double> suds_t::outlier_ths;
int suds_t::required_epoch_n = 5;
double suds_t::required_comp_p = 0.05;
bool suds_t::self_classification = false;
double suds_t::self_classification_prob = 99;
double suds_t::self_classification_kappa = 0;

// 5 -> N1, N2, N3, R, W
// 3 -> NR, R, W (3-stage)
int suds_t::n_stages = 5; 
std::vector<std::string> suds_t::labels;
std::vector<std::string> suds_t::labels3;
std::vector<std::string> suds_t::labels5;
std::vector<std::string> suds_t::labelsR;

bool suds_t::use_fixed_trainer_req;
std::vector<int> suds_t::fixed_trainer_req;
int suds_t::fake_ids;
std::string suds_t::fake_id_root;

std::string suds_t::eannot_file = "";
bool suds_t::eannot_ints = false;
std::string suds_t::eannot_prepend = "";
std::string suds_t::mat_dump_file = "";

std::vector<double> suds_t::lwr_h2, suds_t::upr_h2;
std::vector<double> suds_t::lwr_h3, suds_t::upr_h3;
double suds_t::hjorth_outlier_th = 5;

// 1) Epoch-level PSD for 1+ channels (cs_EEG, etc) 
//    These can be given arbitrary labels for targets, but then
//    targets must match these

// 2) PSC to extract N components from the epoch-level

// 3) Outlier removal

// 4) If manual staging data present, fit LDA
//    -->  save in library

// 5) If manual staging data not present, predict using each of N trainers
//    (after projecting epoch-level PSD into the trainer space) 
//   --> evaluate how well this trainer worked by fitting LDA to predicted stages, as seeing how well that
//       manages to impute other trainers (with known staging)



//
// Self evaluation of staging/signals using SUDS ('SOAP')
//

void suds_indiv_t::evaluate( edf_t & edf , param_t & param )
{

  // track ID (needed if caching for RESOAP)
  id = edf.id;

  // this impacts whether epochs w/ missing values are dropped or not  
  suds_t::soap_mode = 1;

  // ensure we do not call self_classify() from proc
  suds_t::self_classification = false;

  // verbose output
  bool epoch_level_output = param.has( "epoch" );

  //
  // assume that we have manual staging ('true') 
  //

  int n_unique_stages = proc( edf , param , true );
  
  //
  // Cache for RESOAP?
  //

  if ( suds_t::cache_target ) 
    {

      logger << "\n  caching " << id << " for a subsequent RESOAP\n";

      // copy
      suds_t::cached = *this;
      
    }


  //
  // Perhaps no observed stages?
  //

  if ( n_unique_stages < 2 )
    {
      logger << "  *** fewer than 2 non-missing stages for this individual, cannot complete SOAP\n";
      return;
    }


  //
  // fit LDA, and extract posteriors (pp) 
  //

  Eigen::MatrixXd pp;
  
  int dummy = self_classify( NULL , &pp );
  
  if ( dummy == 0 ) 
    {
      logger << "  *** not enough data/variability to fit LDA\n";
      return;
    }

  //
  // output stage probabilities 
  //

  const double epoch_sec = edf.timeline.epoch_length();

  const int ne_all = edf.timeline.num_epochs();

  std::vector<std::string> final_pred = suds_t::max( pp , model.labels );

  summarize_kappa( final_pred , true );

  summarize_stage_durations( pp , model.labels , ne_all , epoch_sec );
  
  if ( epoch_level_output )
    summarize_epochs( pp , model.labels , ne_all , edf );


  //
  // Output annotations (of discordant epochs)
  //

  if ( param.has( "annot" ) )
    {
      const std::string annot_folder = param.has("annot-dir") ? param.value( "annot-dir" ) : "./";      
      write_annots( annot_folder , param.value( "annot" ) , pp , model.labels , ne_all , edf );
    }

}


void suds_indiv_t::resoap_alter1( edf_t & edf , int epoch , suds_stage_t stage )
{
  
  
  // actual number of epochs
  
  const int ne_all = edf.timeline.num_epochs();
  
  // for SOAP, number of included epochs based on good signal data 
  // (i.e. as all epochs are included for SOAP targets, irrespective for 
  // ovserved stage being known or not)
  
  const int ne_included = obs_stage_valid.size();
  
  // nb. 'epoch' is 1-based , given by the user
  if ( epoch < 1 || epoch > ne_all ) 
    Helper::halt( "bad epoch value, outside range" );
  
  // some epochs may be skipped, e.g. due to signal outliers
  
  // valid epochs  : std::vector<int> epochs;   
  // all stages    : std::vector<suds_stage_t> obs_stage; 
  // valid stages  : std::vector<suds_stage_t> obs_stage_valid; 
  // same, but str : std::vector<std::string> y;   (send to lda_t() 
  
  // we need to update only y[i]; so the 'original' is kept in obs_stage(_valid)
    
  //
  // Update 'y' and check we have en
  //
  
  // epochs[] contains the codes of epochs actually present in the model/valid                                                                                                                      
  //  (in 0, ... encoding)
  
  bool updated = false;
  
  for (int i=0; i < epochs.size(); i++) 
    {
      
      // internal epoch number = epochs[i]
      // display epoch = edf.timeline.display_epoch( i )
      // nb. in SOAP context, with no restructuring of the EDF, the display epoch
      // will typically be +1 the internal epoch, i.e. not expecting discontinous
      // codes;  the user is expected to 
      
      // for y and obs_stage_valid
      int e0 = i;
      
      // for obs_stage
      int e1 = epochs[i];
      
      // for user disokay epoch (1-based)
      int e2 = edf.timeline.display_epoch( e1 );

      // update this single 'observed' stage
      if ( epoch == e2 ) 
	{
	  logger << "  changing epoch " << epoch << " from " << y[e0] << " to " << suds_t::str( stage ) << "\n";
	  y[e0] = suds_t::str( stage );
	  // obs_stage_valid[ e0 ] = stage;
	  // obs_stage[ e1 ] = stage;
	  updated = true;
	}
      
      // track what we have     
    }
  
  if ( ! updated ) 
    logger << "  no updates made: did not find epoch " << epoch << " (with valid signal data)\n";
  
  
}


void suds_indiv_t::resoap_pickN( edf_t & edf , int pick )
{
  // for evaluation of SOAP only: 

  // pick N each of 'labels', using the original observed stages
  if ( obs_stage_valid.size() != y.size() )
    Helper::halt( "cannot use RESOAP pick without original staging" );
  
  // first scrub
  for (int i=0; i < suds_t::cached.y.size(); i++)
    suds_t::cached.y[i] = suds_t::str( SUDS_UNKNOWN );
  
  const int nss = suds_t::labels.size();

  std::map<std::string,int> scounts;

  // N or more  versus exactly N
  bool exact = pick < 0 ;
  if ( exact ) pick = -pick;

  const int n = y.size();

  // Yates-Fisher shuffle to get a random ordering
  std::vector<int> a( n );
  CRandom::random_draw( a );
  
  std::set<std::string> done;
  for (int i=0; i<n; i++)
    {
      
      int p = a[i]; // random draw

      std::string ss = suds_t::str( obs_stage_valid[p] );
      if ( ss == "?" ) continue;
	       
      if ( exact )
	{
	  // only add if fewer than needed?
	  int c = scounts[ ss ];	  
	  if ( c < pick )
	    {
	      y[p] = ss;
	      ++scounts[ ss];	      
	    }
	}
      else
	{	    
	  y[p] = ss;
	  ++scounts[ ss ];
	}
      
      // done for this stage?
      if ( scounts[ y[p] ] == pick )
	done.insert( y[p] );

      // all done?
      if ( done.size() == nss ) break;	
    }

}


void suds_indiv_t::resoap( edf_t & edf , bool epoch_level_output )
{

  logger << "  re-SOAPing...\n";

  //
  // this impacts format of epoch-level output
  //

  suds_t::soap_mode = 2;
  
  //
  // Count "observed" stages
  //
  
  const int n = y.size();
  
  std::map<std::string,int> ycounts;
  for (int i=0; i<n; i++) ++ycounts[ y[i] ];
  

  //
  // requires at leasrt two stages w/ at least 3 observations, and has to 
  // be greater than the number of PSCs
  //
  
  const int required_n = 3;
  
  logger << "  epoch counts:";
  int s = 0;
  int t = 0;
  int tt = 0;
  std::map<std::string,int>::const_iterator yy = ycounts.begin();
  while ( yy != ycounts.end() )
    {
      logger << " " << yy->first << ":" << yy->second;     
      tt += yy->second;
      if ( yy->first != "?" && yy->second >= required_n ) 
	{
	  ++s;
	  t += yy->second;
	}
      ++yy;
    }
  logger << "\n";

  writer.value( "S" , s );
  writer.value( "OBS_N" , t ); // need at least 3 of each for 't'
  writer.value( "OBS_P" , t/(double)tt );

  
  bool okay = s >= 2 ;
  
  // for p predictors, require at least p+2 observations
  if ( ! ( t > nc+1 ) ) okay = false;
  
  if ( ! okay )
    {
      logger << "  not enough non-missing stages for LDA with " << nc << " predictors\n";
      writer.value( "FIT" , 0 );
      return;
    }
  
  
  //
  // Re-fit the LDA
  //
  
  Eigen::MatrixXd pp;
  
  int dummy = self_classify( NULL , &pp );

  if ( dummy == 0 )
    {
      logger << "  LDA model could not converge with the current stage proposal\n";
      writer.value( "FIT" , 0 );      
      return;
    }


  //
  // Model okay
  //

  writer.value( "FIT" , 1 );


  //
  // output stage probabilities 
  //
  
  const double epoch_sec = edf.timeline.epoch_length();
  
  std::vector<std::string> final_pred = suds_t::max( pp , model.labels );
  
  summarize_kappa( final_pred , true );
  
  // actual number of epochs
  const int ne_all = edf.timeline.num_epochs();

  summarize_stage_durations( pp , model.labels , ne_all , epoch_sec );

  if ( epoch_level_output )
    summarize_epochs( pp , model.labels , ne_all , edf );
  
}


void suds_indiv_t::add_trainer( edf_t & edf , param_t & param )
{

  // build a trainer; returns number of 'valid'/usable stages  
  int n_unique_stages = proc( edf , param , true );
 
  // only include recordings that have all five/three stages included  
  if ( n_unique_stages != suds_t::n_stages ) 
    {
      logger << "  only found " << n_unique_stages << " of " << suds_t::n_stages << " stages, so not adding as a trainer\n";
      return;
    }  
  
  // save to disk: text or binary format?
  if ( param.has( "text" ) )
    write( edf , param ); 
  else
    binary_write( edf , param ); 

}


int suds_indiv_t::self_classify( std::vector<bool> * included , Eigen::MatrixXd * pp )
{

  if ( ! trainer )
    Helper::halt( "can only self-classify trainers (those w/ observed staging" );
  
  // assume putative 'y' and 'U' will have been constructed, and 'nve' set
  // i.e. this will be called after proc(), or from near the end of proc() 

  //
  // fit the LDA to self
  //
  
  fit_lda();
  
  if ( ! model.valid )
    return 0;
  
  //
  // get predictions
  //

  lda_posteriors_t prediction = lda_t::predict( model , U );

  
  // save posteriors?
  if ( pp != NULL ) *pp = prediction.pp ;

  
  //
  // In SOAP mode, all done (we only needed the PP)
  //
  
  if ( suds_t::soap_mode || included == NULL )
    return 1;  // SOAP only cares about a non-zero return value
  
  //
  // Get kappa 
  //

  double kappa = MiscMath::kappa( prediction.cl , y , suds_t::str( SUDS_UNKNOWN )  );

  included->resize( nve , false );
  
  //
  // Optionally, ask whether trainer passes self-classification kappa threshold.  If not
  // make all epochs 'bad', i.e. so that this trainer will not be used
  //
  
  if ( suds_t::self_classification_kappa <= 1 )
    {
      if ( kappa < suds_t::self_classification_kappa )
	{
	  logger << "  trainer does not meet SOAP kappa " << kappa << " < " << suds_t::self_classification_kappa << "\n";
	  return 0;  // all 'included' false at this point
	}      
    }

  //
  // Determine 'bad' epochs
  //

  int okay = 0;

  // hard calls

  if ( suds_t::self_classification_prob > 1 )
    {
      for (int i=0;i<nve;i++)
	{
	  (*included)[i] = prediction.cl[i] == y[i] ; 
	  if ( (*included)[i] ) ++okay;
	}
    }
  else
    {
      logger << "  using threshold of PP > " << suds_t::self_classification_prob << "\n";

      // map labels to slots in PP matrix (this might be non-standard, e.g. no REM) for a
      // given trainer, and so we cannot assume canonical slot positions

      std::map<std::string,int> label2slot;
      for (int j=0;j<model.labels.size();j++)
	label2slot[ model.labels[j] ] = j ;
      
      // check PP for the observated stage
      for (int i=0;i<nve;i++)
	{
	  std::map<std::string,int>::const_iterator ii = label2slot.find( y[i] );
	  if ( ii == label2slot.end() )
	    Helper::halt( "internal error in suds_indiv_t::self_classify() , unrecognized label" );
	  if ( prediction.pp( i , ii->second ) >= suds_t::self_classification_prob )
	    {
	      (*included)[i] = true;
	      ++okay;
	    }
	}
    }

  return okay;
}


int suds_indiv_t::proc( edf_t & edf , param_t & param , bool is_trainer )
{

  //
  // Is this individual a trainer (with known stages) or no?
  //

  trainer = is_trainer;

  //
  // Initial/total number of components to extract
  // from PSC (although we may only retain nc2 <= nc
  // for the LDA)
  //
  
  nc = suds_t::nc;

  //
  // Number of signals
  //

  const int ns = suds_t::ns;

  
  //
  // Signals from this EDF
  //

  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) );

  if ( signals.size() != ns ) 
    {
      logger << "  ** warning: could not find specified signals\n";
      return 0;
    }
  
  //
  // Resample as needed
  //
  
  for (int s=0;s<ns;s++)
    {
      
      if ( edf.header.is_annotation_channel( signals(s) ) )
	Helper::halt( "cannot specificy annotation channel: " + signals.label(s) );
      
      if ( edf.header.sampling_freq( signals(s) ) != suds_t::sr[s] ) 
	dsptools::resample_channel( edf, signals(s) , suds_t::sr[s] );

    }


  //
  // Epoch 
  //

  const int ne = edf.timeline.first_epoch();

  // nb. below:
  //
  //   ne     total number of epochs
  //
  //   nge    num of epochs with 'valid' staging (i.e. no UNKNOWN, etc)
  //
  //   nve    of nge, number that are a) not statistical outliers for 1+ PSC [ stored in output ]
  //                and optionally, b) correctly self-classified 
								      

  //
  // PSD
  //

  double fft_segment_size = param.has( "segment-sec" ) 
    ? param.requires_dbl( "segment-sec" ) : 4 ;
  
  double fft_segment_overlap = param.has( "segment-overlap" ) 
    ? param.requires_dbl( "segment-overlap" ) : 2 ;
  
  if ( edf.timeline.epoch_length() <= ( fft_segment_size + fft_segment_overlap ) )
    {
      fft_segment_overlap = 0;
      fft_segment_size = edf.timeline.epoch_length();
    }

  window_function_t window_function = WINDOW_TUKEY50;	   
  if      ( param.has( "no-window" ) ) window_function = WINDOW_NONE;
  else if ( param.has( "hann" ) ) window_function = WINDOW_HANN;
  else if ( param.has( "hamming" ) ) window_function = WINDOW_HAMMING;
  else if ( param.has( "tukey50" ) ) window_function = WINDOW_TUKEY50;
  
  
  //
  // Get stage information (for trainers only)
  //

  std::vector<bool> retained( ne , true );

  bool has_prior_staging = false;

  // for SUDS trainers, always load observed stages
  // for SUDS target, load unless told not to

  // for SOAP mode, always load unless told not to ( here, 'target' is 'trainer' and so trainer is T)
  // but in this case, set obs_stage, just make it all missing (i.e. as trainers need to have this
  // populated

  if ( suds_t::soap_mode && suds_t::ignore_target_priors )
    {
      has_prior_staging = false; // nb. this is set back to 'true'  after the following 
      // section, i.e. to do read from annotations, but we need to say we have stages
      // (albeit all unknown ones) for downstream code) 
      obs_stage.resize( ne , SUDS_UNKNOWN );
    }
  else if ( trainer )
    {
      edf.timeline.annotations.make_sleep_stage();
      
      if ( ! edf.timeline.hypnogram.construct( &edf.timeline , param , false ) )
	{
	  if ( suds_t::soap_mode ) return 0; // okay to skip for SOAP
	  // but flag as major prob if a trainer
	  Helper::halt( "problem extracting stage information for trainer" );
	}

      // total number of epochs does not match?
      if ( ne != edf.timeline.hypnogram.stages.size() )
	Helper::halt( "problem extracting stage information for trainer" );
      
      has_prior_staging = true;
      
    }
  else if ( ! suds_t::ignore_target_priors )
    {

      // for targets, manual/prior staging may exist, in which case we'll want to track it for comparison
      // unless we've been explicitly told to ignore it (ignore-prior --> suds_t::ignore_target_priors )

      edf.timeline.annotations.make_sleep_stage();

      has_prior_staging = edf.timeline.hypnogram.construct( &edf.timeline , param , false ) ;
      
      if ( has_prior_staging )
	{
	  // total number of epochs does not match?
	  if ( ne != edf.timeline.hypnogram.stages.size() )
	    Helper::halt( "problem extracting stage information for trainer" );
	}
      
    }

  
  // nb. overkill to have two staging classifications, but keep for now,
  // i.e. we may want to extend one type only

  // number of good (retained) epochs

  int nge = 0 ;

  if ( has_prior_staging )
    {
      
      obs_stage.resize( ne , SUDS_UNKNOWN );
      
      for (int ss=0; ss < ne ; ss++ )
	{
	  if ( edf.timeline.hypnogram.stages[ ss ] == UNSCORED
	       || edf.timeline.hypnogram.stages[ ss ] == LIGHTS_ON
	       || edf.timeline.hypnogram.stages[ ss ] == MOVEMENT
	       || edf.timeline.hypnogram.stages[ ss ] == UNKNOWN ) obs_stage[ss] = SUDS_UNKNOWN;
	  
	  else if ( edf.timeline.hypnogram.stages[ ss ] == WAKE ) obs_stage[ss] = SUDS_WAKE;
	  else if ( edf.timeline.hypnogram.stages[ ss ] == NREM1 ) obs_stage[ss] = suds_t::n_stages == 3 ? SUDS_NR : SUDS_N1;
	  else if ( edf.timeline.hypnogram.stages[ ss ] == NREM2 ) obs_stage[ss] = suds_t::n_stages == 3 ? SUDS_NR : SUDS_N2;
	  else if ( edf.timeline.hypnogram.stages[ ss ] == NREM3
		    || edf.timeline.hypnogram.stages[ ss ] == NREM4 )  obs_stage[ss] = suds_t::n_stages == 3 ? SUDS_NR : SUDS_N3;
	  
	  else if ( edf.timeline.hypnogram.stages[ ss ] == REM ) obs_stage[ss] = SUDS_REM;

	  
	  // expand retained class to exclude unknown staging info
	  // note: even for targets, this means we will not try to guess epochs
	  // that have an UNKNOWN value for the 

	  // note: however, if in 'self-classification' mode
	  // (i.e. running SOAP), then we allow these, i.e. to enable
	  // 'auto-completion' staging (i.e. if some small proportion
	  // is staged, and we want to complete the rest)
	  
	  // lda_t will ignore '?' in fitting the model; i.e. so we can process
	  // just the 'full' dataset (as will be used in the prediction part) 
	  // and that way we are sure that we are doing the SOAP part right

	  if ( suds_t::soap_mode )
	    {
	      ++nge;
	    }
	  else
	    {	      
	      if ( obs_stage[ss] == SUDS_UNKNOWN ) 
		retained[ss] = false; 
	      else ++nge;
	    }

	}
    }
  else 
    {

      //
      // for target individuals without staging, include all epochs
      //
      
      nge = ne;
      
    }

  
  //
  // See note above
  //
  
  if ( suds_t::soap_mode && suds_t::ignore_target_priors )
    {
      has_prior_staging = true;
    }
  
  //
  // for QC, estimate Hjorth parameters (only 2nd and 3rd used) over
  // epochs (for each signal) 
  //
  
  h2.resize( nge , ns );
  h3.resize( nge , ns );


  //
  // For band power analysis (liklely EEG only, as a test condition),
  // track frequecy of column
  //

  std::vector<double> frq;
  Eigen::MatrixXd R; // --> B, raw power (PSD is 10log10(R))
  
  //
  // iterate over (retained) epochs
  //
    
  int en = 0 , en_good = 0;

  edf.timeline.first_epoch();
  
  epochs.clear();
  
  while ( 1 ) 
    {
      
      // select epoch
      int epoch = edf.timeline.next_epoch();      	  
      
      if ( epoch == -1 ) break;
      
      if ( en == ne ) Helper::halt( "internal error: over-counted epochs" );

      // retained? if not, skip
      if ( ! retained[ en ] ) 
	{
	  ++en;
	  continue;
	}

      // col counter for PSD aggregation matrix
      int col = 0;
      std::vector<double> firstrow;
      std::vector<double> firstrow2; // if use_bands, for R matrix (raw power)
      
      // iterate over signals      
      for (int s = 0 ; s < ns; s++ )
	{
	  // get data
	  interval_t interval = edf.timeline.epoch( epoch );	  
	  slice_t slice( edf , signals(s) , interval );
	  std::vector<double> * d = slice.nonconst_pdata();

	  // mean centre epoch
	  MiscMath::centre( d );
	  
	  // pwelch() to obtain full PSD
	  const double overlap_sec = fft_segment_overlap;
	  const double segment_sec  = fft_segment_size;
	  const int total_points = d->size();
	  const int segment_points = segment_sec * suds_t::sr[s];
	  const int noverlap_points  = overlap_sec * suds_t::sr[s];
	  
	  // implied number of segments
	  int noverlap_segments = floor( ( total_points - noverlap_points) 
					 / (double)( segment_points - noverlap_points ) );
	  
	  PWELCH pwelch( *d , 
			 suds_t::sr[s] , 
			 segment_sec , 
			 noverlap_segments , 
			 window_function );
	  	  	  
	  // using bin_t, 1 means no binning
	  bin_t bin( suds_t::lwr[s] , suds_t::upr[s] , 1 ); 
	  
	  bin.bin( pwelch.freq , pwelch.psd );

	  bool has_zeros = false;

	  for ( int i = 0 ; i < bin.bfa.size() ; i++ )
	    {
	      
	      if ( bin.bfa[i] >= suds_t::lwr[s] && bin.bfb[i] <= suds_t::upr[s] )
		{
		  
		  // fudge: for now, if find 0 power value, set to a small value;
		  // and ensure that mobility = 0 so that this epoch will be removed
		  // this may skew the first SVD / outlier removal, but should not 
		  // be too bad.... really should remove these epochs first.
		  
		  if ( bin.bspec[i] <= 0 ) 
		    {
		      has_zeros = true;
		      bin.bspec[i] = 1e-4; // -40dB
		    }
		  
		  if ( en_good == 0 ) firstrow.push_back(  10*log10( bin.bspec[i] ) );
		  else PSD( en_good , col ) = 10*log10( bin.bspec[i] ) ; 		  
		  
		  if ( suds_t::use_bands )
		    {
		      if ( en_good == 0 ) firstrow2.push_back( bin.bspec[i] );
		      else R( en_good , col ) = bin.bspec[i] ;

		      // only track on first epoch 
		      if ( en_good == 0 )
			frq.push_back( bin.bfa[i] );
		      
		    }
		  ++col;		  
		}	      
	    }
	  

	  // Hjorth parameters

	  double activity = 0 , mobility = 0 , complexity = 0;
	  MiscMath::hjorth( d , &activity , &mobility , &complexity );
	  h2(en_good,s) = has_zeros ? 0 : mobility ; // ensure epoch removed if any 0 in the PSD
	  h3(en_good,s) = complexity ;
	  
	} // next signal

    
      // store/shape output if first go around
      nbins = col;
      if ( en_good == 0 )
	{
	  PSD.resize( nge , col );
	  for (int i=0;i<col;i++) PSD(0,i) = firstrow[i];

	  if ( suds_t::use_bands )
	    {
	      R.resize( nge , col );
	      for (int i=0;i<col;i++) R(0,i) = firstrow2[i];
	    }
	  
	}

      
      // increase epoch-number
      ++en;
      ++en_good;
      epochs.push_back( epoch );
    
    } // next epoch
  
  // all done: check

  if ( en_good != nge ) Helper::halt( "internal error: under-counted epochs" );


  //
  // Collapse PSD to bands instead:  
  //

  Eigen::MatrixXd B;
  
  if ( suds_t::use_bands )
    {

      // SLOW DELTA THETA ALPHA SIGMA BETA GAMMA
      // (up to) 7 bands

      double lowf = frq[0] , uprf = frq[0];

      for (int f=0;f<frq.size();f++)
	{
	  if ( frq[f] < lowf ) lowf = frq[f];
	  if ( frq[f] > uprf ) uprf = frq[f];
	}

      // ASSUME... DELTA .. BETA, but GAMMA/SLOW may be absent 
      int nbands = 5;
      bool has_slow = lowf < globals::freq_band[ SLOW  ].second ;
      bool has_gamma = uprf >=  globals::freq_band[ GAMMA  ].first ;;
      if ( has_slow ) ++nbands;
      if ( has_gamma ) ++nbands;      

      B = Eigen::MatrixXd::Zero( nge , nbands );

      const int ncol = frq.size();

      if ( ncol != PSD.cols() ) Helper::halt( "problem" );
      
      for (int j=0; j<ncol ; j++ )
	{
	  int b = -1;
	  if      ( frq[j] >= globals::freq_band[ SLOW  ].first && frq[j] < globals::freq_band[ SLOW  ].second ) b = has_slow ? 0 : 0;
	  else if ( frq[j] >= globals::freq_band[ DELTA ].first && frq[j] < globals::freq_band[ DELTA ].second ) b = has_slow ? 1 : 0;
	  else if ( frq[j] >= globals::freq_band[ THETA ].first && frq[j] < globals::freq_band[ THETA ].second ) b = has_slow ? 2 : 1;
	  else if ( frq[j] >= globals::freq_band[ ALPHA ].first && frq[j] < globals::freq_band[ ALPHA ].second ) b = has_slow ? 3 : 2;
	  else if ( frq[j] >= globals::freq_band[ SIGMA ].first && frq[j] < globals::freq_band[ SIGMA ].second ) b = has_slow ? 4 : 3;
	  else if ( frq[j] >= globals::freq_band[ BETA  ].first && frq[j] < globals::freq_band[ BETA  ].second ) b = has_slow ? 5 : 4;
	  else if ( frq[j] >= globals::freq_band[ GAMMA ].first && frq[j] < globals::freq_band[ GAMMA ].second ) b = has_slow ? 6 : 5;
	  
	  if ( b != -1 )
	    for (int i=0; i<nge; i++) B(i,b) += R(i,j);

	}
    }

  
  //
  // Rescale PSD?
  //

  if ( suds_t::standardize_psd )
    {
      if ( suds_t::robust_standardization )
	{
	  logger << "  robust standardizing PSD";
	  if ( suds_t::winsor1 > 0 ) logger << ", winsorizing at " << suds_t::winsor1;
	  logger << "\n";
	  eigen_ops::robust_scale( PSD , suds_t::winsor1 );
	}
      else
	{
	  logger << "  standardizing PSD\n";
	  eigen_ops::scale( PSD , true );
	}
    }
  

  //
  // Get PSC initially (we look for outliers and then remove epochs, and redo the SVD)
  //

  // mean-centre columns if not already done via standardization
  
  if ( ! suds_t::standardize_psd )
    eigen_ops::scale( U , false );

  //
  // SVD
  //

  // W.resize( nbins ); 
  // V.resize( nbins , nbins );
  
  Eigen::BDCSVD<Eigen::MatrixXd> svd( PSD , Eigen::ComputeThinU | Eigen::ComputeThinV );
  U = svd.matrixU();
  V = svd.matrixV();
  W = svd.singularValues();
  
  //
  // Outliers/smoothing
  //
  
  std::vector<bool> valid( nge , true );

  // track reasons for exclusion
  int nout_flat = 0;
  int nout_hjorth = 0;
  int nout_stat = 0;
  
  //
  // Exclusions based on H==0 parameters
  //
  
  for ( int s=0;s<ns;s++)
    {
      for (int i=0;i<nge;i++) 
	{
	  if ( h2( i, s ) < 1e-8 ) { valid[i] = false; nout_flat++; }
	  else if ( h3( i , s ) < 1e-8 ) { valid[i] = false; nout_flat++; } 
	}
    }

  //
  // For targets only, threshold epochs based on per-signal Hjorth
  // from trainers
  // 

  if ( ! trainer )
    {
      logger << "  removing epochs +/-" << suds_t::hjorth_outlier_th << " SD units from the H2 & H3 trainer means\n";

      for ( int s=0;s<ns;s++)
	{
	  for (int i=0;i<nge;i++)
	    {
	      if ( h2( i, s ) <= suds_t::lwr_h2[s] || h2(i,s) >= suds_t::upr_h2[s] ||
		   h3( i, s ) <= suds_t::lwr_h3[s] || h3(i,s) >= suds_t::upr_h3[s] ) { valid[i] = false; nout_hjorth++; } 
	    }
	}
    }


  //
  // Component-based epoch-outlier removal (after removing flat lines)
  //

  for (int o=0;o< suds_t::outlier_ths.size();o++)
    {
      logger << "  removing epochs +/-" << suds_t::outlier_ths[o] << " from PSC means\n";
      for ( int j=0;j<nc;j++)
	{
	  std::vector<double> x;
	  for (int i=0;i<nge;i++) if ( valid[i] ) x.push_back( U(i,j) );
	  if ( x.size() < 2 ) Helper::halt( "no epochs left" );
	  double mean = MiscMath::mean( x );
	  double sd = MiscMath::sdev( x , mean );
	  double lwr = mean - suds_t::outlier_ths[o] * sd;
	  double upr = mean + suds_t::outlier_ths[o] * sd;
	  int c = 0;
	  for (int i=0;i<nge;i++)
	    {
	      if ( valid[i] )
		{
		  if ( x[c] < lwr || x[c] > upr ) { valid[i] = false; nout_stat++; } 
		  ++c;
		}
	    }
	}
    }

  int included = 0;
  for (int i=0;i<nge;i++)
    if ( valid[i] ) ++included;

  logger << "  of " << ne << " total epochs, valid staging for " << nge
         << ", and of those " << included << " passed outlier removal\n";

  logger << "  outliers counts (flat, Hjorth, components = " << nout_flat << ", " << nout_hjorth << ", " << nout_stat << ")\n";


  //
  // Remove bad epochs and repeat (SVD and smoothing)
  //

  // nve = number of valid epochs ( ne > nge > nve ) 
  
  nve = included;

  Eigen::MatrixXd PSD2 = PSD;
  PSD.resize( nve , nbins );
  std::vector<int> epochs2 = epochs;
  epochs.clear();
  
  int r = 0;
  for (int i=0;i<PSD2.rows() ; i++)
    {      
      if ( valid[i] )
	{
	  for (int j=0;j<nbins;j++)
	    PSD(r,j) = PSD2(i,j);
	  
	  epochs.push_back( epochs2[i] );
	  
	  ++r;
	}
    }
  
  // only retain nve obs labels from obs_stage[ne] originals

  if ( has_prior_staging )
    {
      obs_stage_valid.clear();
      
      r = 0;
      for (int i=0;i<ne;i++)
	{
	  if ( retained[i] )
	    {
	      if ( valid[r] )
		obs_stage_valid.push_back( obs_stage[ i ] );
	      ++r;
	    }
	}
    }

  //
  // optional, band-power per epoch tracking
  //

  if ( suds_t::use_bands )
    {
      // and make as log
      int nbands = B.cols();
      Eigen::MatrixXd B2 = B;      
      B.resize( nve , nbands );
      
      int r = 0;
      for (int i=0;i<B2.rows() ; i++)
	{      
	  if ( valid[i] )
	    {
	      for (int j=0;j< nbands; j++)
		B(r,j) = 10*log10( B2(i,j) );
	      ++r;
	    }
	}
    }

  
  //
  // splice out bad epochs for Hjorth parameters
  //
  
  Eigen::MatrixXd hh2 = h2;
  Eigen::MatrixXd hh3 = h3;
  h2.resize( nve , ns ); 
  h3.resize( nve , ns );

  for (int s=0;s<ns;s++)
    {
      int r = 0;      
      for (int i=0; i < valid.size(); i++ )
	if ( valid[i] )
	  {
	    h2(r,s) = hh2(i,s);
	    h3(r,s) = hh3(i,s);
	    ++r;
	  }
    }



  //
  // Rescale PSD?
  //

  if ( suds_t::standardize_psd )
    {

      if ( suds_t::robust_standardization )
	{
	  logger << "  robust re-standardizing PSD after removing bad epochs\n";
	  // nb. not repeating winsorization of PSD here
	  eigen_ops::robust_scale( PSD , 0 );
	  if ( suds_t::use_bands ) eigen_ops::robust_scale( B , 0 );
	}
      else
	{      
	  logger << "  re-standardizing PSD after removing bad epochs\n";	  
	  eigen_ops::scale( PSD , true );
	  if ( suds_t::use_bands ) eigen_ops::scale( B , true );      
	}
    }
  else // just ensure we mean-center in any case
    {
      // mean-center columns (PSD)  
      eigen_ops::scale( PSD , false );
      
      if ( suds_t::use_bands ) eigen_ops::scale( B , false );
      
    }

  

  //
  // Get PSC (post outlier removal)
  //

  // W.resize( nbins ); 
  // V.resize( nbins , nbins );

  Eigen::BDCSVD<Eigen::MatrixXd> svd2( PSD , Eigen::ComputeThinU | Eigen::ComputeThinV );
  U = svd2.matrixU();
  V = svd2.matrixV();
  W = svd2.singularValues();
  

  //
  // Standardize PSC 
  //

  if ( suds_t::standardize_psc )
    {
      if ( suds_t::robust_standardization )
	{
	  logger << "  robust standardizing PSC\n";
	  eigen_ops::robust_scale( U , 0 ); // no repeated winsorization here
	}
      else
	{
	  logger << "  standardizing PSC\n";
	  eigen_ops::scale( U , true );
	}  
    }
  
  //
  // Smooth PSCs?
  //
  
  if ( suds_t::denoise_fac > 0 ) 
    {      
      logger << "  smoothing PSCs lambda=" << suds_t::denoise_fac << " * SD\n";
      
      for (int j=0;j<nc;j++)
	{
	  //std::vector<double> * col = U.col_nonconst_pointer(j)->data_nonconst_pointer();
	  double sd = suds_t::standardize_psc ? 1 : eigen_ops::sdev( U.col(j) );
	  double lambda = suds_t::denoise_fac * sd;	  
	  dsptools::TV1D_denoise( U.col(j) , lambda );
	}
      
      if ( suds_t::use_bands )
	{
	  int nbands = B.cols();
	  for (int j=0; j<nbands; j++)
	    {
	      //std::vector<double> * col = B.col_nonconst_pointer(j)->data_nonconst_pointer();
	      double sd = eigen_ops::sdev( B.col(j) );
	      double lambda = suds_t::denoise_fac * sd;
	      dsptools::TV1D_denoise( B.col(j) , lambda );
	    }
	  
	}

    }
  
  
  //
  // For trainers, optionally only retain PSCs (or bands) that are significantly
  // associated with observed stage in this individual
  //

  if ( trainer && suds_t::required_comp_p < 1 && ! ( suds_t::soap_mode && suds_t::ignore_target_priors ) ) 
    {

      // pull out currently retained epochs
      std::vector<std::string> ss_str;
      int c = 0;
      for ( int i = 0 ; i < ne ; i++ )
        {
          if ( retained[i] )
            {
              if ( valid[c] )
                ss_str.push_back( suds_t::str( obs_stage[i] ) );
              ++c;
            }
        }

      std::set<int> incl_comp;
      for (int j=0;j<nc;j++)
	{	  
	  // standardize column
	  Eigen::VectorXd c = U.col(j);
	  eigen_ops::scale( c , true );
	  
	  double pv = Statistics::anova( ss_str  , eigen_ops::copy_vector( c ) );
	  if ( pv >= 0 && pv <  suds_t::required_comp_p  ) incl_comp.insert( j );
	  writer.level( "PSC_" + Helper::int2str( j+1 ) , "VAR");
	  writer.value( "PV", pv );
	  writer.value( "INC" , pv >= 0 && pv <  suds_t::required_comp_p ); 
	}
      
      //
      // Optionally, compare to band-power association w/ stage for the same set of epochs
      //

      if ( suds_t::use_bands )
	{	  
	 
	  int nbands = B.cols();

	  // fix, for now just assume always has SLOW band... 
	  std::vector<std::string> bands7 = { "SLOW" , "DELTA" , "THETA" , "ALPHA" , "SIGMA" , "BETA" , "GAMMA" } ;
	  std::vector<std::string> bands6 = { "SLOW" , "DELTA" , "THETA" , "ALPHA" , "SIGMA" , "BETA" } ;
	  std::vector<std::string> bands = nbands == 7 ? bands7 : bands6;
	  
	  for (int j=0;j< nbands ;j++)
	    {	      
	      
	      Eigen::VectorXd c = B.col(j);
	      
	      if ( 1 ) // no variance check here now (might need to add back in)
		{
		  eigen_ops::scale( c , true );
		  double pv = Statistics::anova( ss_str  , eigen_ops::copy_vector( c ) );
		  writer.level( bands[j] , "VAR" );
		  writer.value( "PV" , pv  );
		  writer.value( "INC" , pv >= 0 && pv <  suds_t::required_comp_p ); 
		}
	      else
		{		  
                  writer.level( bands[j] , "VAR" );
		  writer.value( "INC" , 0 );
		  
		}
	      
	    }
	}

      writer.unlevel( "VAR" );

      
      //
      // no usable components --> no usable epochs... quit out (this trainer will be ignored)
      //
      
      if ( incl_comp.size() == 0 )
	{
	  logger << "  0 components associated with stage at p<" << suds_t::required_comp_p << ", bailing\n";
	  return 0;
	}
      
      
      //
      // and prune U and V down here
      //
	  
      const int nc2 = incl_comp.size();
      std::vector<bool> incl( nc );
      for (int j=0;j<nc;j++) incl[j] = incl_comp.find( j ) != incl_comp.end();
      
      Eigen::MatrixXd U2 = U;
      U.resize( nve , nc2 );
      for (int i=0;i<nve;i++)
	{
	  int cc = 0;
	  for (int j=0;j<nc;j++)
	    if ( incl[j] ) U(i,cc++) = U2(i,j);
	}
      
      Eigen::ArrayXd W2 = W;
      W.resize( nc2 );
      int cc = 0;
      for (int j=0;j<nc;j++)
	if ( incl[j] ) W[ cc++ ] = W2[ j ];
      
      Eigen::MatrixXd VV = V;
      V.resize( VV.rows() , nc2 );
      for (int i=0;i<VV.rows();i++)
	{
	  int cc = 0;
	  for (int j=0;j<nc;j++)
	    if ( incl[j] ) V(i,cc++) = VV(i,j);
	}
      
          
      //
      // set new 'nc' for this individual (which may be less than suds_t::nc)
      //
      
      logger << "  retaining " << incl_comp.size() << " of " << nc 
	     << " PSCs, based on ANOVA p<" << suds_t::required_comp_p << "\n";

      nc = incl_comp.size();
      
    }

  //
  // Make variables for LDA: shrink down to 'nc' (if not already done by the above
  // component selection step)
  //
  
  if ( U.cols() != nc )
    {
      Eigen::MatrixXd U2 = U;
      U.resize( nve , nc );
      for (int i=0;i<nve;i++)
	for (int j=0;j<nc;j++)
	  U(i,j) = U2(i,j);
      
      W.conservativeResize( nc );
      
      Eigen::MatrixXd VV = V;
      V.resize( VV.rows() , nc );
      for (int i=0;i<VV.rows();i++)
	for (int j=0;j<nc;j++)
	  V(i,j) = VV(i,j);
    }


 
  //
  // Re-Standardize PSC 
  //
  
  if ( suds_t::standardize_psc )
    {

      if ( suds_t::robust_standardization )
	{
	  logger << "  robust re-standardizing PSC";
	  if ( suds_t::winsor2 > 0 ) logger << ", winsorizing at " << suds_t::winsor2;
	  logger << "\n";
	  eigen_ops::robust_scale( U , suds_t::winsor2 );
	}
      else
	{
	  logger << "  re-standardizing PSC\n";
	  eigen_ops::scale( U , true );
	}  
    }
    

  //
  // make class labels ( trainer only )
  //
  
  if ( trainer )
    {
      
      y.clear();
      
      int c = 0;
      for ( int i = 0 ; i < ne ; i++ )
	{
	  if ( retained[i] )
	    {
	      if ( valid[c] )
		y.push_back( suds_t::str( obs_stage[i] ) );
	      ++c;
	    }
	}

      counts.clear();
      for (int i=0;i<y.size();i++) counts[y[i]]++;
      std::map<std::string,int>::const_iterator cc = counts.begin();
      logger << "  epoch counts:";
      while ( cc != counts.end() )
	{
	  logger << " " << cc->first << ":" << cc->second ;
	  ++cc;
	}
      logger << "\n";
    }


  //
  // fit model based based only on band power
  //

  if ( suds_t::use_bands )
    {
      
      // fit based on PSC
      lda_t lda1( y , U );      
      lda_model_t m1 = lda1.fit( suds_t::flat_priors );
      lda_posteriors_t prediction1 = lda_t::predict( m1 , U );
      double kappa1 = MiscMath::kappa( prediction1.cl , y , suds_t::str( SUDS_UNKNOWN ) );
      
      // fit based on band power
      lda_t lda2( y , B );
      lda_model_t m2 = lda2.fit( suds_t::flat_priors );
      lda_posteriors_t prediction2 = lda_t::predict( m2 , B );
      double kappa2 = MiscMath::kappa( prediction2.cl , y , suds_t::str( SUDS_UNKNOWN ) );

      writer.value( "K_PSC" , kappa1 );
      writer.value( "K_BAND" , kappa2);
        
      if ( suds_t::epoch_lvl_output )
	{

	  std::map<int,int> e2e;
	  for (int i=0; i< epochs.size(); i++) e2e[ epochs[i] ] = i ;  
	  const int ne_all = edf.timeline.num_epochs();
	  //if ( prediction1.cl.size() != ne_all ) Helper::halt( "need to adjust in bands lda_t() " );
	  // std::cout << "ne all = " << ne_all << "\n";
	  // std::cout << prediction1.cl.size()  << "\n";
	  // std::cout << epochs.size() << "\n";
	  
	  for (int i=0; i< ne_all; i++)
	    {
	      int e = -1;
	      if ( e2e.find( i ) != e2e.end() ) e = e2e[i];
	      if ( e == -1 ) continue;

	      // show display epoch 'i'
	      writer.epoch( edf.timeline.display_epoch( i ) ) ;

	      //std::cout << "check " << i << " " << edf.timeline.display_epoch( i ) << " " << e << " " << prediction1.cl.size() << "\n";
	      
	      writer.value( "PRED_PSC" , prediction1.cl[e]  );
	      writer.value( "PRED_BAND" , prediction2.cl[e]  );

	      //     	  std::cout << "PSC\n";
	      // std::cout << U.print() << "\n";
	      
	      // std::cout << "PSC POSTERIORS\n";
	      // std::cout << prediction1.pp.print() << "\n";	  
	      
	      // std::cout << "BANDS\n";
	      // std::cout << B.print() << "\n";
	      
	      // std::cout << "BANDS POSTERIORS\n";
	      // std::cout << prediction2.pp.print() << "\n";	  
	      
	      
	    }
	  writer.unepoch();
	  
	}
           
    }


  //
  // Verbose output of PSC  (todo)
  //

  
  //
  // Attempt self-classification, to remove epochs that aren't well self-classified (i.e.
  // do not fit the model) and possible to reject a trainer, if their kappa is sufficiently
  // poor?
  //
  
  if ( trainer && suds_t::self_classification )
    {

      std::vector<bool> okay ;

      int nve2 = self_classify( &okay );
      
      if ( nve2 == 0 )
	{
	  logger << "  trainer not valid based on self-classification thresholds\n";
	  return 0;
	}
      

      //
      // Subset epochs:
      //
      
      //   U  PSD  epochs  y  h2  h3
      
      Eigen::MatrixXd UU = U;      
      U.resize( nve2 , nc );      

      Eigen::MatrixXd PSD2 = PSD;
      PSD.resize( nve2 , nbins );

      std::vector<int> epochs2 = epochs;
      epochs.clear();

      std::vector<suds_stage_t> obs_stage_valid2 = obs_stage_valid;
      if ( has_prior_staging )
	obs_stage_valid.clear();

      Eigen::MatrixXd hh2 = h2;
      h2.resize( nve , ns );
      
      Eigen::MatrixXd hh3 = h3;
      h3.resize( nve , ns );
            
      int r = 0;
      for (int i=0;i< nve; i++)
	{      
	  if ( okay[i] )
	    {

	      // U
	      for (int j=0;j<nc;j++)
		U(r,j) = UU(i,j);

	      // X (PSD)
	      for (int j=0;j<nbins;j++)
		PSD(r,j) = PSD2(i,j);

	      // Epoch tracking
	      epochs.push_back( epochs2[i] );

	      // nb. take from already-pruned set, obs_stage_valid[] ) 
	      if ( has_prior_staging )
		{
		  // std::cout << "hmm " << i << " " 
		  // 	    << obs_stage_valid.size() << " " << obs_stage_valid2.size() << " " << nve << "\n";
		  obs_stage_valid.push_back( obs_stage_valid2[i] );
		}

	      // Hjorth (per signal)
	      for (int s=0;s<ns;s++)
		{
		  h2(r,s) = hh2(i,s);
		  h3(r,s) = hh3(i,s);		  
		}

	      // next good epoch
	      ++r;
	    }
	}
    

      //
      // Redo labels
      //

      std::vector<std::string> yy = y;
      y.clear();
      for ( int i = 0 ; i < nve ; i++ )
	if ( okay[i] ) y.push_back( yy[i] );

      // just in case...?
      if ( y.size() != obs_stage_valid.size() ) Helper::halt( "internal error in proc()" );

      //
      // update nve
      //

      nve = nve2;
      
      //
      // recount stages
      //

      counts.clear();
      for (int i=0;i<y.size();i++) counts[y[i]]++;
      std::map<std::string,int>::const_iterator cc = counts.begin();
      logger << "  updated epoch counts:";
      while ( cc != counts.end() )
	{
	  logger << " " << cc->first << ":" << cc->second ;
	  ++cc;
	}
      logger << "\n";

      logger << "  final count of valid epochs is " << nve << "\n";
    }

  
  //
  // Summarize mean/SD for per-signal Hjorth parameters
  //

  mean_h2 = h2.colwise().mean();
  mean_h3 = h3.colwise().mean();

  sd_h2 = ((h2.array().rowwise() - mean_h2 ).square().colwise().sum()/(h2.rows()-1)).sqrt();
  sd_h3 = ((h3.array().rowwise() - mean_h3 ).square().colwise().sum()/(h3.rows()-1)).sqrt();

  
  // for trainers, returns number of observed stages w/ at least suds_t::required_epoch_n 
  // -- i.e. should be suds_t::n_stages

  int nr = 0;
  std::map<std::string,int>::const_iterator cc = counts.begin();
  while ( cc != counts.end() )
    {
      if ( cc->first != "?" && cc->second >= suds_t::required_epoch_n ) ++nr;      
      ++cc;
    }

  return trainer ? nr : nve ;
  
}



void suds_indiv_t::write( edf_t & edf , param_t & param ) const
{

  const std::string folder = Helper::expand( param.requires( "db" ) );

  const int ns = suds_t::ns;
    
  //
  // Store as an epoch-level EDF
  //

  // Save: 
  //  ID
  //  U [ nve , nc ]  D [ nc ]  V [ nc , nc ] 
  //  stages [ nve ]  epochs[ nve ]

  // create output folder if it does not exist
  std::string syscmd = globals::mkdir_command + " " + folder ;
  int retval = system( syscmd.c_str() );


  // for saving trainers: use EDF ID, or a fake ID?  (e.g. 'ids=suds')
  std::string suds_id = suds_t::fake_ids ? suds_t::fake_id_root + "_" + Helper::int2str( suds_t::fake_ids++ ) : edf.id;
  
  std::string filename = folder + globals::folder_delimiter + suds_id;

  logger << "  writing trainer data to " << filename << "\n";
  
  std::ofstream OUT1( filename.c_str() , std::ios::out );

  // file version code
  OUT1 << "SUDS\t1\n";
  
  OUT1 << "ID\t" << suds_id << "\n"
       << "N_VALID_EPOCHS\t" << nve << "\n"
       << "N_X\t" << nbins << "\n"
       << "N_SIGS\t" << ns << "\n"
       << "N_COMP\t" << nc << "\n";

  // channels , SR [ for comparability w/ other data ] 

  for (int s=0;s<ns;s++)
    {
      OUT1 << "\nCH\t" << suds_t::siglab[s] << "\n"
	   << "SR\t" << suds_t::sr[s] << "\n"
	   << "LWR\t" << suds_t::lwr[s] << "\n"
	   << "UPR\t" << suds_t::upr[s] << "\n"	   
	   << "H2_MN\t" << mean_h2[s] << "\n"
	   << "H2_SD\t" << sd_h2[s] << "\n"
	   << "H3_MN\t" << mean_h3[s] << "\n"
	   << "H3_SD\t" << sd_h3[s] << "\n";
    }


  // stages
  OUT1 << "\nN_STAGES\t" << counts.size() << "\n";
  
  std::map<std::string,int>::const_iterator ss = counts.begin();
  while ( ss != counts.end() )
    {
      OUT1 << ss->first << "\t" << ss->second << "\n";
      ++ss;
    }
  
  // W
  OUT1 << "\nW[" << nc << "]";
  for (int j=0;j<nc;j++)
    OUT1 << " " << W[j];
  OUT1 << "\n";

  // V
  OUT1 << "\nV[" << nbins << "," << nc << "]";
  for (int i=0;i<nbins;i++)
    for (int j=0;j<nc;j++)
      OUT1 << " " << V(i,j);
  OUT1 << "\n";
  
  // stages
  OUT1 << "\nEPOCH_STAGE";
  for (int i=0;i<nve;i++)
    OUT1 << " " << epochs[i] << " " << y[i] ;
  OUT1 << "\n\n";

  // U (to reestimate LDA model upon loading, i.e
  //  to use lda.predict() on target   ; only needs to be nc rather than nbins
  OUT1 << "U[" << nve << "," << nc << "]";
  for (int i=0;i<nve;i++)
    for (int j=0;j<nc;j++)
      OUT1 << " " << U(i,j);
  OUT1 << "\n\n";
  
  // X, RAW DATA (e.g. mean-centered PSD, but possibly other things)
  // i.e. if this trainer is being used as a 'weight trainer',
  // i.e. will project this individuals raw data into the target space
  OUT1 << "X[" << nve << "," << nbins << "]";
  for (int i=0;i<nve;i++)
    for (int j=0;j<nbins;j++)
      OUT1 << " " << PSD(i,j);
  OUT1 << "\n\n";
  
  OUT1.close();

  //
  // All done
  //

}


void suds_indiv_t::write( const std::string & filename ) const
{

  const int ns = suds_t::ns;
  
  logger << "  writing trainer data to " << filename << "\n";
  
  std::ofstream OUT1( filename.c_str() , std::ios::out );

  // file version code
  OUT1 << "SUDS\t1\n";
  
  OUT1 << "ID\t" << id << "\n"
       << "N_VALID_EPOCHS\t" << nve << "\n"
       << "N_X\t" << nbins << "\n"
       << "N_SIGS\t" << ns << "\n"
       << "N_COMP\t" << nc << "\n";

  // channels , SR [ for comparability w/ other data ] 

  for (int s=0;s<ns;s++)
    {
      OUT1 << "\nCH\t" << suds_t::siglab[s] << "\n"
	   << "SR\t" << suds_t::sr[s] << "\n"
	   << "LWR\t" << suds_t::lwr[s] << "\n"
	   << "UPR\t" << suds_t::upr[s] << "\n"	   
	   << "H2_MN\t" << mean_h2[s] << "\n"
	   << "H2_SD\t" << sd_h2[s] << "\n"
	   << "H3_MN\t" << mean_h3[s] << "\n"
	   << "H3_SD\t" << sd_h3[s] << "\n";
    }


  // stages
  OUT1 << "\nN_STAGES\t" << counts.size() << "\n";
  
  std::map<std::string,int>::const_iterator ss = counts.begin();
  while ( ss != counts.end() )
    {
      OUT1 << ss->first << "\t" << ss->second << "\n";
      ++ss;
    }
  
  // W
  OUT1 << "\nW[" << nc << "]";
  for (int j=0;j<nc;j++)
    OUT1 << " " << W[j];
  OUT1 << "\n";

  // V
  OUT1 << "\nV[" << nbins << "," << nc << "]";
  for (int i=0;i<nbins;i++)
    for (int j=0;j<nc;j++)
      OUT1 << " " << V(i,j);
  OUT1 << "\n";
  
  // stages
  OUT1 << "\nEPOCH_STAGE";
  for (int i=0;i<nve;i++)
    OUT1 << " " << epochs[i] << " " << y[i] ;
  OUT1 << "\n\n";

  // U (to reestimate LDA model upon loading, i.e
  //  to use lda.predict() on target   ; only needs to be nc rather than nbins
  OUT1 << "U[" << nve << "," << nc << "]";
  for (int i=0;i<nve;i++)
    for (int j=0;j<nc;j++)
      OUT1 << " " << U(i,j);
  OUT1 << "\n\n";
  
  // X, RAW DATA (e.g. mean-centered PSD, but possibly other things)
  // i.e. if this trainer is being used as a 'weight trainer',
  // i.e. will project this individuals raw data into the target space
  OUT1 << "X[" << nve << "," << nbins << "]";
  for (int i=0;i<nve;i++)
    for (int j=0;j<nbins;j++)
      OUT1 << " " << PSD(i,j);
  OUT1 << "\n\n";
  
  OUT1.close();

  //
  // All done
  //

}


void suds_indiv_t::bwrite( std::ofstream & O , const std::string & s ) 
{
  uint8_t l = s.size();
  O.write( (char*)( &l ), sizeof(uint8_t) );
  O.write( s.c_str(), l );
}

void suds_indiv_t::bwrite( std::ofstream & O , int i ) 
{
  O.write( (char*)( &i ), sizeof(int) );
}

void suds_indiv_t::bwrite( std::ofstream & O , double d ) 
{
  O.write( (char*)( &d ), sizeof(double) );
}

std::string suds_indiv_t::bread_str( std::ifstream & I )
{
  uint8_t len;
  I.read( (char*)( &len ), sizeof(uint8_t) );
  std::vector<char> b( len );
  I.read( &b[0] , len );
  std::string s( b.begin() , b.end() );
  return s;
}

int suds_indiv_t::bread_int( std::ifstream & I )
{
  int i;
  I.read( (char*)( &i ), sizeof(int) );
  return i;
}

double suds_indiv_t::bread_dbl( std::ifstream & I )
{
  double d;
  I.read( (char*)( &d ), sizeof(double) );
  return d;
}


void suds_indiv_t::binary_write( edf_t & edf , param_t & param ) const
{
  
  // same as write(), except the file is binary 
  
  const std::string folder = Helper::expand( param.requires( "db" ) );
  const int ns = suds_t::ns;
  
  // create output folder if it does not exist
  std::string syscmd = globals::mkdir_command + " " + folder ;
  int retval = system( syscmd.c_str() );

  // for saving trainers: EDF ID versus a dummy ID (e.g. 'ids=suds')
  std::string suds_id = suds_t::fake_ids ? suds_t::fake_id_root + "_" + Helper::int2str( suds_t::fake_ids++ ) : edf.id;
  
  std::string filename = folder + globals::folder_delimiter + suds_id ;
  
  logger << "  writing binary-format trainer data to " << filename << "\n";
  
  std::ofstream OUT1( filename.c_str() , std::ios::binary | std::ios::out );

  //
  // file version code
  //

  bwrite( OUT1 , "SUDS1" );
  bwrite( OUT1 , suds_id );
  bwrite( OUT1 , nve );
  bwrite( OUT1 , nbins );
  bwrite( OUT1 , ns );
  bwrite( OUT1 , nc );

  // channels , SR [ for comparability w/ other data ] 

  for (int s=0;s<ns;s++)
     {
       bwrite( OUT1 , suds_t::siglab[s] );
       bwrite( OUT1 , suds_t::sr[s] );
       bwrite( OUT1 , suds_t::lwr[s] );
       bwrite( OUT1 , suds_t::upr[s] );
       bwrite( OUT1 , mean_h2[s] );
       bwrite( OUT1 , sd_h2[s]  );
       bwrite( OUT1 , mean_h3[s] );
       bwrite( OUT1 , sd_h3[s]  );
     }

  // stages (N)
  bwrite( OUT1 , (int)counts.size() );
  
  std::map<std::string,int>::const_iterator ss = counts.begin();
  while ( ss != counts.end() )
    {
      bwrite( OUT1 , ss->first );
      bwrite( OUT1 , ss->second );
      ++ss;
    }
  
  // W [nc]
  for (int j=0;j<nc;j++)
    bwrite( OUT1 , W[j] );
  
   // V [nbins x nc ]
  for (int i=0;i<nbins;i++)
    for (int j=0;j<nc;j++)
      bwrite( OUT1 , V(i,j) );
  
  
  // stages (nve)
  for (int i=0;i<nve;i++)
    {
      bwrite( OUT1 , epochs[i] );
      bwrite( OUT1 , y[i] );
    }
      
  // U (to reestimate LDA model upon loading, i.e
  //  to use lda.predict() on target   ; only needs to be nc rather than nbins
  // U [ nve x nc ] 
  for (int i=0;i<nve;i++)
    for (int j=0;j<nc;j++)
      bwrite( OUT1 , U(i,j) );
  
  // X, RAW DATA (e.g. mean-centered PSD, but possibly other things)
  // i.e. if this trainer is being used as a 'weight trainer',
  // i.e. will project this individuals raw data into the target space
  //  X [ nve x nbins ]
  for (int i=0;i<nve;i++)
    for (int j=0;j<nbins;j++)
      bwrite( OUT1 , PSD(i,j) );
  
  OUT1.close();
  
  //
  // All done
  //
  
}


void suds_indiv_t::binary_write( const std::string & filename ) const
{
  
  // same as binary_write(), except a different interface (for
  // use w/ --copy-suds command;   in future we should make this
  // a single function...)

  const int ns = suds_t::ns;
  
  logger << "  writing binary-format trainer data to " << filename << "\n";
  
  std::ofstream OUT1( filename.c_str() , std::ios::binary | std::ios::out );

  //
  // file version code
  //

  bwrite( OUT1 , "SUDS1" );
  bwrite( OUT1 , id );
  bwrite( OUT1 , nve );
  bwrite( OUT1 , nbins );
  bwrite( OUT1 , ns );
  bwrite( OUT1 , nc );

  // channels , SR [ for comparability w/ other data ] 

  for (int s=0;s<ns;s++)
     {
       bwrite( OUT1 , suds_t::siglab[s] );
       bwrite( OUT1 , suds_t::sr[s] );
       bwrite( OUT1 , suds_t::lwr[s] );
       bwrite( OUT1 , suds_t::upr[s] );
       bwrite( OUT1 , mean_h2[s] );
       bwrite( OUT1 , sd_h2[s]  );
       bwrite( OUT1 , mean_h3[s] );
       bwrite( OUT1 , sd_h3[s]  );
     }

  // stages (N)
  bwrite( OUT1 , (int)counts.size() );
  
  std::map<std::string,int>::const_iterator ss = counts.begin();
  while ( ss != counts.end() )
    {
      bwrite( OUT1 , ss->first );
      bwrite( OUT1 , ss->second );
      ++ss;
    }
  
  // W [nc]
  for (int j=0;j<nc;j++)
    bwrite( OUT1 , W[j] );
  
   // V [nbins x nc ]
  for (int i=0;i<nbins;i++)
    for (int j=0;j<nc;j++)
      bwrite( OUT1 , V(i,j) );
  
  
  // stages (nve)
  for (int i=0;i<nve;i++)
    {
      bwrite( OUT1 , epochs[i] );
      bwrite( OUT1 , y[i] );
    }
      
  // U (to reestimate LDA model upon loading, i.e
  //  to use lda.predict() on target   ; only needs to be nc rather than nbins
  // U [ nve x nc ] 
  for (int i=0;i<nve;i++)
    for (int j=0;j<nc;j++)
      bwrite( OUT1 , U(i,j) );
  
  // X, RAW DATA (e.g. mean-centered PSD, but possibly other things)
  // i.e. if this trainer is being used as a 'weight trainer',
  // i.e. will project this individuals raw data into the target space
  //  X [ nve x nbins ]
  for (int i=0;i<nve;i++)
    for (int j=0;j<nbins;j++)
      bwrite( OUT1 , PSD(i,j) );
  
  OUT1.close();
  
  //
  // All done
  //
  
}


void suds_indiv_t::binary_reload( const std::string & filename , bool load_rawx )
{

  std::ifstream IN1( filename.c_str() , std::ios::binary | std::ios::in );

  std::string dummy;
  std::string suds = bread_str( IN1 );

  if ( suds != "SUDS1" )
    Helper::halt( "bad file format for " + filename );
  
  id = bread_str( IN1 );
  nve = bread_int( IN1 );
  nbins = bread_int( IN1 );
  int this_ns = bread_int( IN1 );
  int this_nc = bread_int( IN1 );

  if ( this_nc == 0 )
    Helper::halt( "0 PSCs for " + filename );

  // reading or just copying?
  if ( suds_t::copy_db_mode )
    {
      // these will not get used, so can set to whatever;
      // they will be set for this person being read, and will
      // be available for the next operation, which is to write
      // the data back for the same person, in a different format (binary/text)
      suds_t::ns = this_ns;
      suds_t::siglab.resize( this_ns );
      suds_t::sr.resize( this_ns );
      suds_t::lwr.resize( this_ns );
      suds_t::upr.resize( this_ns );
    }
  else if (this_ns != suds_t::ns )
    Helper::halt( "different trainer ns=" + Helper::int2str( this_ns )
		  + " in " + filename
		  + ", expecting " + Helper::int2str( suds_t::ns ) ) ; 
  
  // set 'nc' for this individual
  nc = this_nc;

  // ns should be fixed across all individuals (unless in copy mode)
  const int ns = suds_t::ns;
  
  mean_h2.resize(ns); mean_h3.resize(ns);
  sd_h2.resize(ns); sd_h3.resize(ns);
  
  for (int s=0;s<ns;s++)
    {
      std::string this_siglab = bread_str( IN1 );
      int this_sr = bread_int( IN1 );      

      double this_lwr =  bread_dbl( IN1 );
      double this_upr =  bread_dbl( IN1 );
      
      mean_h2[s] = bread_dbl( IN1 );
      sd_h2[s] = bread_dbl( IN1 );
      
      mean_h3[s] = bread_dbl( IN1 );
      sd_h3[s] = bread_dbl( IN1 );

      if( suds_t::copy_db_mode )
	{
	  suds_t::siglab[s] = this_siglab;
	  suds_t::sr[s] = this_sr;
	  suds_t::lwr[s] = this_lwr;
	  suds_t::upr[s] = this_upr;
	}
      else
	{      
	  if ( this_siglab != suds_t::siglab[s] ) Helper::halt( "different signals: " + this_siglab
								+ ", but expecting " + suds_t::siglab[s] );
	  if ( this_sr != suds_t::sr[s] ) Helper::halt( "different SR: " + Helper::int2str( this_sr )
							+ ", but expecting " + Helper::int2str( suds_t::sr[s] ) );
	  if ( this_lwr != suds_t::lwr[s] ) Helper::halt( "different lower-freq: " + Helper::dbl2str( this_lwr )
							  + ", but expecting " + Helper::dbl2str( suds_t::lwr[s] ) );
	  if ( this_upr != suds_t::upr[s] ) Helper::halt( "different upper-freq: " + Helper::dbl2str( this_upr )
							  + ", but expecting " + Helper::dbl2str( suds_t::upr[s] )) ;      
	}
    }

  // stages
  int nstages = bread_int( IN1 );

  for (int i=0;i<nstages;i++)
    {
      std::string sname = bread_str( IN1 );
      int scnt = bread_int( IN1 );
      counts[ sname ] = scnt;
    }
  
  
  // check that these equal suds_t values?
  
  // W [ only nc ]
  W.resize( nc );
  for (int j=0;j<nc;j++)
    W[j] = bread_dbl( IN1 );

  // V [ only nc cols ] 
  V.resize( nbins, nc );
  for (int i=0;i<nbins;i++)
    for (int j=0;j<nc;j++)
      V(i,j) = bread_dbl( IN1 );

  // stages
  y.resize( nve );
  epochs.resize( nve );
  for (int i=0;i<nve;i++)
    {
      epochs[i] = bread_int( IN1 );
      y[i] = bread_str( IN1 );
    }

  obs_stage = suds_t::type( y );
  
  // U (to reestimate LDA model upon loading, i.e
  //  to use lda.predict() on target 
  U.resize( nve , nc );
  for (int i=0;i<nve;i++)
    for (int j=0;j<nc;j++)
      U(i,j) = bread_dbl( IN1 );
  
  // X values (e.g. mean-centered PSD)
  // i.e. if this trainer is being used as a 'weight trainer',
  // i.e. will project this individuals raw data into the target space
  
  if ( load_rawx )
    {      
      // PSD      
      PSD.resize( nve , nbins );
      for (int i=0;i<nve;i++)
	for (int j=0;j<nbins;j++)
	  PSD(i,j) = bread_dbl( IN1 );
    }
  
  IN1.close();
  
  //
  // All done
  //

}


void suds_indiv_t::reload( const std::string & filename , bool load_rawx )
{
  
  std::ifstream IN1( filename.c_str() , std::ios::in );

  std::string dummy;
  std::string suds;
  int version;
  IN1 >> suds >> version;
  
  if ( suds != "SUDS" )
    Helper::halt( "bad file format for " + filename );
  
  int this_ns, this_nc;

  IN1 >> dummy >> id 
      >> dummy >> nve 
      >> dummy >> nbins 
      >> dummy >> this_ns 
      >> dummy >> this_nc ;

  if ( this_nc == 0 )
    Helper::halt( "0 PSCs for " + filename );
  
  
  if ( suds_t::copy_db_mode )
    {
      suds_t::ns = this_ns;
      suds_t::siglab.resize( this_ns );
      suds_t::sr.resize( this_ns );
      suds_t::lwr.resize( this_ns );
      suds_t::upr.resize( this_ns );
    }
  else if (this_ns != suds_t::ns )
    Helper::halt( "different trainer ns=" + Helper::int2str( this_ns )
		  + " in " + filename
		  + ", expecting " + Helper::int2str( suds_t::ns ) ) ; 

  // set 'nc' for this individual
  nc = this_nc;

  // ns should be fixed across all individuals
  const int ns = suds_t::ns;
  
  mean_h2.resize( ns ); mean_h3.resize( ns );
  sd_h2.resize( ns ); sd_h3.resize( ns );
  
  for (int s=0;s<ns;s++)
    {
      std::string this_siglab;
      double this_lwr, this_upr;
      int this_sr;
      IN1 >> dummy >> this_siglab
	  >> dummy >> this_sr 
	  >> dummy >> this_lwr
	  >> dummy >> this_upr	  
	  >> dummy >> mean_h2[s] >> dummy >> sd_h2[s]
	  >> dummy >> mean_h3[s] >> dummy >> sd_h3[s];


      if( suds_t::copy_db_mode )
        {
	  suds_t::siglab[s] = this_siglab;
	  suds_t::sr[s] = this_sr;
	  suds_t::lwr[s] = this_lwr;
	  suds_t::upr[s] = this_upr;
        }
      else
        {
	  if ( this_siglab != suds_t::siglab[s] ) Helper::halt( "different signals: " + this_siglab
								+ ", but expecting " + suds_t::siglab[s] );
	  if ( this_sr != suds_t::sr[s] ) Helper::halt( "different SR: " + Helper::int2str( this_sr )
							+ ", but expecting " + Helper::int2str( suds_t::sr[s] ) );
	  if ( this_lwr != suds_t::lwr[s] ) Helper::halt( "different lower-freq: " + Helper::dbl2str( this_lwr )
							  + ", but expecting " + Helper::dbl2str( suds_t::lwr[s] ) );
	  if ( this_upr != suds_t::upr[s] ) Helper::halt( "different upper-freq: " + Helper::dbl2str( this_upr )
							  + ", but expecting " + Helper::dbl2str( suds_t::upr[s] )) ;      
	}
    }

  // stages
  int nstages;
  IN1 >> dummy >> nstages;
  for (int i=0;i<nstages;i++)
    {
      std::string sname;
      int scnt;
      IN1 >> sname >> scnt;
      counts[ sname ] = scnt;
    }
  
  
  // check that these equal suds_t values?
  
  // W [ only nc ]
  IN1 >> dummy;
  W.resize( nc );
  for (int j=0;j<nc;j++)
    IN1 >> W[j] ;

  // V [ only nc cols ] 
  IN1 >> dummy;
  V.resize( nbins, nc );
  for (int i=0;i<nbins;i++)
    for (int j=0;j<nc;j++)
      IN1 >> V(i,j);

  // stages
  IN1 >> dummy;
  y.resize( nve );
  epochs.resize( nve );
  for (int i=0;i<nve;i++)
    IN1 >> epochs[i] >> y[i] ;
  obs_stage = suds_t::type( y );
  
  // U (to reestimate LDA model upon loading, i.e
  //  to use lda.predict() on target 
  IN1 >> dummy;
  U.resize( nve , nc );
  for (int i=0;i<nve;i++)
    for (int j=0;j<nc;j++)
      IN1 >> U(i,j) ;	

  // X values (e.g. mean-centered PSD)
  // i.e. if this trainer is being used as a 'weight trainer',
  // i.e. will project this individuals raw data into the target space
  
  if ( load_rawx )
    {      
      // PSD
      IN1 >> dummy;
      PSD.resize( nve , nbins );
      for (int i=0;i<nve;i++)
	for (int j=0;j<nbins;j++)
	  IN1 >> PSD(i,j) ;
    }

  IN1.close();
  
  //
  // All done
  //

}


void suds_t::attach_db( const std::string & folder0 , bool binary , bool read_psd )
{

  const std::string folder = Helper::expand( folder0 );
  
  std::map<std::string,suds_indiv_t*> * b = read_psd ? &wbank : &bank ;
    
  // already done?
  if ( b->size() > 0 ) return;

  if ( bank.size() == 0 && wbank.size() == 0 ) 
    logger << "  attaching training data from " << folder << " ...\n";

  // find all files in this folder
  std::vector<std::string> trainer_ids;

  DIR *dir;
  struct dirent *ent;
  if ( (dir = opendir ( folder.c_str() ) )  != NULL )
    {
      /* print all the files and directories within directory */
      while ( (ent = readdir (dir)) != NULL) {

#ifdef _DIRENT_HAVE_D_TYPE
        if ( ent->d_type != DT_REG ) continue;
#endif

        std::string fname = ent->d_name;
        if (  fname == "." || fname == ".." ) continue;

	// in 'single-trainer' test mode, only attach one trainer
	// from all 
	
	if ( single_trainer != "" && single_trainer != fname ) continue;

        // otherwise, we have not yet come across this person, so load
        // up...

	trainer_ids.push_back( fname );
	
      }
      closedir (dir);
    }
  else
    {
      Helper::halt( "could not open directory " + folder );      
    }
  
  //
  // for primary trainers only (! read_psd ) track H2 and H3 distributions
  //

  Eigen::MatrixXd h2m( trainer_ids.size() , suds_t::ns );
  Eigen::MatrixXd h2sd( trainer_ids.size() , suds_t::ns );
  Eigen::MatrixXd h3m( trainer_ids.size() , suds_t::ns );
  Eigen::MatrixXd h3sd( trainer_ids.size() , suds_t::ns );

  
  //
  // load each
  //
  
  for ( int i=0; i<trainer_ids.size() ; i++)
    {

      // is this person the other bank?  we always attach wdb first
      // (i.e. to make sure PSD data are included, if these are
      // needed).  So, check against wbank in that case:

      if ( i % 50 == 0 ) logger << "\n ";
      if ( i % 10 == 0 ) logger << " ";
      logger << ".";

      bool already_loaded =  ( ! read_psd ) && wbank.find( trainer_ids[i] ) != wbank.end() ; 

      suds_indiv_t * trainer = NULL;
      
      if ( already_loaded )
	{
	  // copy pointer over, to already-loaded (wdb) version
	  trainer = wbank[  trainer_ids[i]  ];
	}
      else
	{
	  
	  // otherwise, we need to read in and process as necessary

	  trainer = new suds_indiv_t;
	  
	  if ( binary ) 
	    trainer->binary_reload( folder + globals::folder_delimiter + trainer_ids[i] , read_psd );
	  else
	    trainer->reload( folder + globals::folder_delimiter + trainer_ids[i] , read_psd ); 

	  trainer->fit_lda();
	  
	}

      // store in the relevant bank:
	  
      (*b)[ trainer_ids[i] ] = trainer;

      // for primary trainers only, calculate Hjorth parameters
      
      if ( ! read_psd ) 
	{
	  for (int s=0;s<suds_t::ns;s++)
	    {
	      h2m(i,s) = trainer->mean_h2[s] ;
	      h2sd(i,s) = trainer->sd_h2[s] ;
	      h3m(i,s) = trainer->mean_h3[s] ;
	      h3sd(i,s) = trainer->sd_h3[s] ;
	    }
	}
          
    }
     
  
  logger << "\n"
	 << "  attached " << b->size() << " trainers ("
	 << ( binary ? "binary" : "text" ) << " format, "
	 << ( read_psd ? "with spectra" : "w/out spectra" )
	 << ") from " << folder << "\n";


  //
  // from primary trainers only, track signal-wise Hjorth limits
  //
  

  if ( ! read_psd )
    {
      suds_t::lwr_h2.resize( suds_t::ns );
      suds_t::upr_h2.resize( suds_t::ns );
      suds_t::lwr_h3.resize( suds_t::ns );
      suds_t::upr_h3.resize( suds_t::ns );
      
      Eigen::ArrayXd mean_h2m = h2m.colwise().mean();
      Eigen::ArrayXd mean_h2sd = h2sd.colwise().mean();
      Eigen::ArrayXd mean_h3m = h3m.colwise().mean();
      Eigen::ArrayXd mean_h3sd = h3sd.colwise().mean();

      for (int s=0;s<suds_t::ns;s++)
	{
	  suds_t::lwr_h2[s] = mean_h2m[s] - suds_t::hjorth_outlier_th * mean_h2sd[s];
	  suds_t::upr_h2[s] = mean_h2m[s] + suds_t::hjorth_outlier_th * mean_h2sd[s];

	  suds_t::lwr_h3[s] = mean_h3m[s] - suds_t::hjorth_outlier_th * mean_h3sd[s];
	  suds_t::upr_h3[s] = mean_h3m[s] + suds_t::hjorth_outlier_th * mean_h3sd[s];

	  if ( suds_t::lwr_h2[s] < 0 ) suds_t::lwr_h2[s] = 0;
	  if ( suds_t::lwr_h3[s] < 0 ) suds_t::lwr_h3[s] = 0;
	  
	  logger << "  thresholding " << suds_t::siglab[s] 
		 << " on H2: " << suds_t::lwr_h2[s] << " - " << suds_t::upr_h2[s] 
		 << " and H3: " << suds_t::lwr_h3[s] << " - " << suds_t::upr_h3[s] << "\n";
	}
    }

}



void suds_t::copy_db( const std::string & folder1 , 
		      const std::string & folder2 , 
		      bool from_text )
{
  
  // will set SUDS variables to whatever is just read, i.e. no checking with 
  // the attached EDF, as there is no attached EDF
  
  suds_t::copy_db_mode = true;
  
  // text folder
  std::string tfolder = from_text ? folder1 : folder2;
  
  // binary folder
  std::string bfolder = from_text ? folder2 : folder1;

  logger << "  copying from " << ( from_text ? "[text]" : "[binary]" ) << " " << folder1 << "\n"
	 << "            to " << ( !from_text ? "[text]" : "[binary]" ) << " " << folder2 << "\n";

  //
  // copy from one folder to another, either text->binary or binary->text
  //

  // find all files in this folder (tfolder)
  std::vector<std::string> trainer_ids;

  DIR *dir;
  struct dirent *ent;
  if ( (dir = opendir ( folder1.c_str() ) )  != NULL )
    {
      /* print all the files and directories within directory */
      while ( (ent = readdir (dir)) != NULL) {

#ifdef _DIRENT_HAVE_D_TYPE
        if ( ent->d_type != DT_REG ) continue;
#endif

        std::string fname = ent->d_name;
        if (  fname == "." || fname == ".." ) continue;

        // otherwise, we have not yet come across this person, so load
        // up...

	trainer_ids.push_back( fname );
	
      }
      closedir (dir);
    }
  else
    {
      Helper::halt( "could not open directory " + folder1 );      
    }
  

  //
  // Ensure target folder exists
  //

  std::string syscmd = globals::mkdir_command + " " + folder2 ;
  int retval = system( syscmd.c_str() );
  
  //
  // Read/write each trainer
  //

  for ( int i=0; i< trainer_ids.size() ; i++)
    {
      
      suds_indiv_t * trainer = new suds_indiv_t;
      
      //
      // read text --> write binary 
      //
      
      if ( from_text )
	{
	  trainer->reload( tfolder + globals::folder_delimiter + trainer_ids[i] , true ); 

	  trainer->binary_write( bfolder + globals::folder_delimiter + trainer_ids[i] );
	}
      else
	{
	  //
	  // read  binary --> write text
	  //
	  
	  trainer->binary_reload( bfolder + globals::folder_delimiter + trainer_ids[i] , true ); 

	  trainer->write( tfolder + globals::folder_delimiter + trainer_ids[i] );

	}
	  
      
    }

  logger << "  copied " << trainer_ids.size() << " trainers from " 
	 << ( from_text ? "text to binary" : "binary to text" )
	 << " format\n"
	 << "   from folder: " << folder1 << "\n"
	 << "     to folder: " << folder2 << "\n";

}




// fit LDA, i.e. after reloading U

void suds_indiv_t::fit_lda()
{

  lda_t lda( y , U );

  model = lda.fit( suds_t::flat_priors );
    
}


//
// make predictions given a different individual's signal data
//

lda_posteriors_t suds_indiv_t::predict( const suds_indiv_t & trainer )
{

  //
  // Project target (this) into trainer space:   U_targ = X_targ * V_trainer * D_trainer^{-1} 
  // subsetting to # of columns
  //

  // std::cout << " targ = " << id << " triner = " << trainer.id << "\n";
  // std::cout << " tr nc = " << trainer.nc << " " << trainer.W.size() << "\n";
  // std::cout << "trainer W \n" << trainer.W << "\n";

  Eigen::MatrixXd trainer_DW = Eigen::MatrixXd::Zero( trainer.nc , trainer.nc );  

  for (int i=0;i< trainer.nc; i++)
    trainer_DW(i,i) = 1.0 / trainer.W[i];
  
  U_projected = PSD * trainer.V * trainer_DW;

  //
  // Normalize PSC?
  //

  if ( suds_t::standardize_psc ) 
    {
      if ( suds_t::robust_standardization )
	{
	  eigen_ops::robust_scale( U_projected , suds_t::winsor2 );
	}
      else
	{
	  eigen_ops::scale( U_projected , true );
	}
    }


  //
  // smooth U (projected)
  //
  
  if ( suds_t::denoise_fac > 0 ) 
    for (int j=0; j< trainer.nc; j++)
      {
	double sd = suds_t::standardize_psc ? 1 : eigen_ops::sdev( U_projected.col(j) );
	double lambda = suds_t::denoise_fac * sd;	
	dsptools::TV1D_denoise( U_projected.col(j) , lambda );
      }  

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
    
    }


  //
  // predict using trainer model
  //

  lda_posteriors_t pp = lda_t::predict( trainer.model , U_projected ) ;

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
  // about this (i.e. amnd vectors will be properly +1 sized to accomodate 
  // all trainers

  if ( suds_t::cheat ) bank_contains_target = false;

  const int bank_size = bank_contains_target ? bank.size() - 1 : bank.size() ; 
    
  //
  // save weights for each trainer, based on re-predicting
  //

  // either kappa or MCC (if 'mcc' option)
  Eigen::ArrayXd wgt_mean = Eigen::ArrayXd::Zero( bank_size ) ;

  // kappa; not used in any calcs
  Eigen::ArrayXd wgt_max = Eigen::ArrayXd::Zero( bank_size ) ;
  Eigen::ArrayXd wgt_n50 = Eigen::ArrayXd::Zero( bank_size ) ;

  // soap weights
  Eigen::ArrayXd wgt_soap = Eigen::ArrayXd::Zero( bank_size ) ;

  

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

  if ( w1 && w2 ) logger << "  using mean of repred-weights & KL-weights\n";
  else if ( w1 ) logger << "  using repred-weights only\n";
  else if ( w2 ) logger << "  using KL-weights only\n";
  else if ( w0 ) logger << "  using SOAP-weights only\n";
  else logger << "  not applying any weights\n";
  
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
      
      //      std::cout << "considering " << trainer->id << "\n";

      //
      // Predict target given trainer, after project target PSD into trainer-defined space 
      // ( i.e. this generates target.U_projected based on trainer, and then uses it to 
      //        predict target class given the trainer model )
      //
      
      lda_posteriors_t prediction = target.predict( *trainer );

      
      //
      // Save predictions
      //
           
      target.add( trainer->id , prediction );

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
      int n_kappa50 = 0;
      int n_kappa_all = 0;
      
      std::map<std::string,int> counts;
      for (int i=0;i<prediction.cl.size();i++) counts[ prediction.cl[ i ] ]++;
      int nr= 0 ; 
      std::map<std::string,int>::const_iterator cc = counts.begin();
      while ( cc != counts.end() )
	{
	  if ( cc->second >= suds_t::required_epoch_n ) 
	    { 
	      ++nr;
	      //std::cout << "  inc " << cc->first << "\n";
	    }
	  
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
	  double kappa3 =  MiscMath::kappa( NRW( str( target.prd_stage ) ) , NRW( str( target.obs_stage_valid ) ) , suds_t::str( SUDS_UNKNOWN ) );	  
	  
	  k3_prior[ cntr ] = kappa3;
	  
	  // tmp
	  k3 = kappa3;
	  
	}
      

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
	  // but U basd on the target's own SVD (i.e. not projected into trainer space);  
	  // This we use target.U, which is the original for the target, based on their own data
	  // (we ignore the U_projected which is based on the trainer model)
	  //

	  lda_t lda( prediction.cl , target.U ) ;
      
	  // set target model for use w/ all different weight-trainers

	  target.model = lda.fit( suds_t::flat_priors );



	  //
	  // Consider one or more weight-trainers for this trainer
	  //   Generally: P_C|B|A
	  //   Or, only considering seflf:   P_A|B|A
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
	      
	      lda_posteriors_t reprediction = weight_trainer->predict( target );
	      
	      weight_trainer->prd_stage = suds_t::type( reprediction.cl );
	      
	      // obs_stage for predicted/valid epochs only

	      double kappa = 0 ; 
	      if ( use_5class_repred ) 
		kappa = MiscMath::kappa( reprediction.cl , str( weight_trainer->obs_stage ) , suds_t::str( SUDS_UNKNOWN )  ) ;
	      else if ( use_rem_repred ) 
		kappa = MiscMath::kappa( Rnot( reprediction.cl ) , Rnot( str( weight_trainer->obs_stage ) ) , suds_t::str( SUDS_UNKNOWN )  );
	      else
		kappa = MiscMath::kappa( NRW( reprediction.cl ) , NRW( str( weight_trainer->obs_stage ) ) , suds_t::str( SUDS_UNKNOWN )  );
	      
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
	      
	      if ( suds_t::verbose ) 
		{
		  wtrainer_mean_k3[ weight_trainer->id ] += kappa;
		  wtrainer_count_k3[ weight_trainer->id ]++;
		}
	      
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
	  wgt_n50[ cntr ] = n_kappa50;
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

	  lda_posteriors_t prediction1 = lda_t::predict( target.model , target.U );


	  double kappa1 = 0 ;

	  if ( use_5class_repred )
	    kappa1 = MiscMath::kappa( prediction1.cl , prediction.cl , suds_t::str( SUDS_UNKNOWN ) );
	  else if ( use_rem_repred )
	    kappa1 = MiscMath::kappa( Rnot( prediction1.cl ) , Rnot( prediction.cl ) , suds_t::str( SUDS_UNKNOWN ) );
	  else 
	    kappa1 = MiscMath::kappa( NRW( prediction1.cl ) , NRW( prediction.cl ) , suds_t::str( SUDS_UNKNOWN ) );

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

      if ( w1 && w2 )
	wgt[ cntr ] = ( wgt_mean[ cntr ] +  wgt_kl[ cntr ] ) / 2.0 ; 
      else if ( w1 )
	wgt[ cntr ] = wgt_mean[ cntr ] ;
      else if ( w2 )
	wgt[ cntr ] = wgt_kl[ cntr ] ;
      else if ( w0 )
	wgt[ cntr ] = wgt_soap[ cntr ];
      else
	wgt[ cntr ] = 1 ; 

      used_trainers.push_back( trainer->id );
			     
      ++tt;
      ++cntr;
    }
  writer.unlevel( "TRAINER" );


  
  //
  // Verbose output: mean weight trainer values
  //

  if ( suds_t::verbose && use_repred_weights && wbank.size() > 0 )
    {
      std::map<std::string,suds_indiv_t*>::iterator ww = wbank.begin();
      while ( ww != wbank.end() )
	{	  
	  suds_indiv_t * weight_trainer = ww->second;
	  writer.level( weight_trainer->id , "WTRAINER" );
	  double m = wtrainer_mean_k3[ weight_trainer->id ] / (double)wtrainer_count_k3[ weight_trainer->id ];
	  writer.value( "K3" , m );
	  ++ww;
	}
      writer.unlevel( "WTRAINER" );
            
    }


  //
  // Normalize wgt / truncate at percentile?
  //

  bool has_wgt = ( wbank.size() > 0 && use_repred_weights ) || use_kl_weights || use_soap_weights ;

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
  // Percentile based sclaing (subsetting) 
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

  target.summarize_epochs( pp , suds_t::labels , ne_all , edf );
  

  //
  // Summarize staging
  //

  const double epoch_sec = edf.timeline.epoch_length();

  target.summarize_stage_durations( pp , suds_t::labels , ne_all , epoch_sec );


  
  //
  // Confiusion matrics and kappa w/ observed staging
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
      
      writer.value( "R_K3_WGT" ,
		    Statistics::correlation( eigen_ops::copy_array( wgt ) ,
					     eigen_ops::copy_array( k3_prior) ) ); 
      
      // TODO ; not critical code, so lazy conversion to use Stats correl function, but should add eigen_ops::correlation()

      if ( 0 )
	{
	  
	  if ( k3_prior.size() > 2 )
	    {
	      if ( use_kl_weights )
		writer.value( "R_K3_KL" , Statistics::correlation( eigen_ops::copy_array( wgt_kl ) , eigen_ops::copy_array( k3_prior) ) ); 
	      
	      if ( use_soap_weights )
		writer.value( "R_K3_SOAP" , Statistics::correlation( eigen_ops::copy_array( wgt_soap ) , eigen_ops::copy_array( k3_prior ) ) ); 
	      
	      if ( use_repred_weights && wbank.size() > 0 ) 
		{
		  //	      writer.value( "R_K3_MAX" , Statistics::correlation( eigen_ops::copy_array( wgt_max ) , eigen_ops::copy_array( k3_prior ) ) ); 
		  writer.value( "R_K3_WGT" , Statistics::correlation( eigen_ops::copy_array(wgt_mean) , eigen_ops::copy_array(k3_prior) ) ); 
		  //writer.value( "R_K3_N50" , Statistics::correlation( eigen_ops::copy_array(wgt_n50)  , eigen_ops::copy_array(k3_prior) ) ); 
	      
		  if ( use_kl_weights )
		    {
		      writer.value( "R_MEAN_KL" , Statistics::correlation( eigen_ops::copy_array( wgt_mean ) , eigen_ops::copy_array( wgt_kl ) ) );
		      writer.value( "R_K3_CMB" , Statistics::correlation( eigen_ops::copy_array(wgt) , eigen_ops::copy_array(wgt_kl) ) );		  
		    }
		  
		}
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

      for (int i=0;i<target.PSD.cols();i++)
	OUT1 << "\t" << "X" << (i+1);

      for (int i=0;i<target.U.cols();i++)
	OUT1 << "\t" << "PSC" << (i+1);

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
	  
	  for (int j=0;j<target.PSD.cols();j++)
	    OUT1 << "\t" << target.PSD(e,j); 
	  
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




void suds_indiv_t::add( const std::string & trainer_id , const lda_posteriors_t & prediction )
{
  
  target_posteriors[ trainer_id ] = prediction.pp ;

  target_predictions[ trainer_id ] = suds_t::type( prediction.cl );
  
}



std::map<std::string,std::map<std::string,int> > suds_t::tabulate( const std::vector<std::string> & a , 
								   const std::vector<std::string> & b , 
								   const bool print  )
{
  std::map<std::string,std::map<std::string,int> > res;
  
  const int n = a.size();
  
  if ( n != b.size() ) 
    Helper::halt( "internal error: unequal vectors in tabulate()" );

  // includes unknown stages SUDS_UNKNOWN, '?' in table
  // (but these should be removed from kappa and other stats)

  std::set<std::string> uniq;
  for (int i=0;i<n;i++)
    {
      //std::cout << "CHECK\t" << a[i] << "\t" << b[i] << "\n";
      res[ a[i] ][ b[i] ]++;
      uniq.insert( a[i] );
      uniq.insert( b[i] );
    }

  std::map<std::string,double> rows, cols;
  double tot = 0;
  std::set<std::string>::const_iterator uu = uniq.begin();
  while ( uu != uniq.end() )
    {
      std::set<std::string>::const_iterator jj = uniq.begin();
      while ( jj != uniq.end() )
	{
	  if ( res.find( *uu ) == res.end() )
	    res[ *uu ][ *jj ] = 0;
	  else
	    {
	      std::map<std::string,int> & rjj = res.find(*uu)->second;
	      if ( rjj.find( *jj ) == rjj.end() )
		res[ *uu ][ *jj ] = 0;
	    }	      

	  // col/row marginals
	  rows[ *uu ] += res[ *uu ][ *jj ] ;
	  cols[ *jj ] += res[ *uu ][ *jj ] ;
	  tot += res[ *uu ][ *jj ];
	  
	  ++jj;
	}
      ++uu;
    }
  
  
  if ( print )
    {

      logger << "\t  Obs:\n\t";
      std::set<std::string>::const_iterator uu = uniq.begin();
      while ( uu != uniq.end() )
	{
	  logger << "\t" << *uu;
	  ++uu;
	}
      logger << "\tTot\n";	
      
      logger << "  Pred:";
      uu = uniq.begin();
      while ( uu != uniq.end() )
	{
	  logger << "\t" << *uu;
	  std::set<std::string>::const_iterator jj = uniq.begin();
	  while ( jj != uniq.end() )
	    {
	      logger << "\t" << res[ *uu ][ *jj ];
	      ++jj;
	    }
	  // row sums
	  logger << "\t" << rows[ *uu ]/tot;
	  logger << "\n";
	  ++uu;
	}
      // col sums
      logger << "\tTot:";
      std::set<std::string>::const_iterator jj = uniq.begin();
      while ( jj != uniq.end() )
	{
	  logger << "\t" << cols[ *jj ]/tot;
	  ++jj;
	}
      logger << "\t1.00\n";
    }

  
  return res;
}



Eigen::ArrayXd suds_indiv_t::wgt_kl() const { 

  // returned weights
  const int nt = target_predictions.size();  

  Eigen::ArrayXd W( nt );

  if ( nt == 0 ) return W;

  Eigen::MatrixXd Q( nt , suds_t::n_stages ) ;  

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

  // divergence for each trainer from the mean
  r = 0;
  ii = target_predictions.begin();
  while ( ii != target_predictions.end() ) 
    {      
      // negative KL
      double ss = 0;
      for ( int s = 0 ; s < suds_t::n_stages ; s++ )
	if ( Q(r,s) > 0 ) ss += P[s] * log( P[s] / Q(r,s) );  
      W[r] = -ss;
      
      ++r;
      ++ii;
    }
    
  return W;
}


void suds_indiv_t::write_annots( const std::string & annot_folder , const std::string & aname , 
				 const Eigen::MatrixXd & pp , const std::vector<std::string> & labels , 
				 int ne_all , edf_t & edf )
{

  bool prior_staging = obs_stage.size() != 0 ;
  if ( ! prior_staging ) return;

  const std::string delim = annot_folder[ annot_folder.size() - 1 ] != '/' ? "/" : "";
  
  if ( annot_folder != "" && annot_folder != "./" ) 
    {
      std::string syscmd = globals::mkdir_command + " " + annot_folder ;
      int retval = system( syscmd.c_str() );
    }
  
  // annot label
  annot_t * a_disc3 = edf.timeline.annotations.add( aname + "_disc3" );
  a_disc3->description = "SOAP NR/R/W discordance";

  annot_t * a_disc5 = NULL;
  if ( suds_t::n_stages == 5 )
    {
      a_disc5 = edf.timeline.annotations.add( aname + "_disc5" );
      a_disc5->description = "SOAP N1/N2/N3/R/W discordance";
    }

  annot_t * a_unscr = edf.timeline.annotations.add( aname + "_unscr" );
  a_unscr->description = "SOAP unscored epoch";

  const std::string a_filename3 = annot_folder + delim + aname + "_disc3.annot";
  const std::string a_filename5 = annot_folder + delim + aname + "_disc5.annot";
  const std::string a_filenameU = annot_folder + delim + aname + "_unscr.annot";

  logger << "  writing NR/R/W discordant epochs to " << a_filename3 << "\n";
  if ( suds_t::n_stages == 5 )
    logger << "  writing N1/N2/N3/R/W discordant epochs to " << a_filename5 << "\n";  
  logger << "  writing unscored epochs to " << a_filenameU << "\n";
   
  // epochs[] contains the codes of epochs actually present in the model/valid
  std::map<int,int> e2e;
  for (int i=0; i < epochs.size(); i++) 
    e2e[ epochs[i] ] = i ;  
  
  for ( int i=0;i<ne_all;i++)
    {
      int e = -1;
      if ( e2e.find( i ) != e2e.end() ) e = e2e[i];
            
      // track interval
      interval_t interval = edf.timeline.epoch( i );

      // value found in scoring?
      if ( e != -1 ) 
	{
	  std::string predss = suds_t::max_inrow( pp.row(e) , labels );

	  if ( suds_t::n_stages == 5 )
	    {
	      if ( predss !=  suds_t::str( obs_stage[i] ) )
		a_disc5->add( suds_t::str( obs_stage[i] ) + "->" + predss , interval , "." );
	      
	      if ( suds_t::NRW( predss ) != suds_t::NRW( suds_t::str( obs_stage[i] ) ) )
		a_disc3->add( suds_t::NRW( suds_t::str( obs_stage[i] ) ) + "->" + suds_t::NRW( predss ) , interval , "." );
	    }
	  else
	    {
	      if ( predss !=  suds_t::str( obs_stage[i] ) )
		a_disc3->add( suds_t::str( obs_stage[i] ) + "->" + predss , interval , "." );
	    }
	}
      else
       	{
	  a_unscr->add( "." , interval , "." );
	}
      
    }

  a_disc3->save( a_filename3 );

  if ( a_disc5 ) 
    a_disc5->save( a_filename5 );

  a_unscr->save( a_filenameU );

}

void suds_indiv_t::summarize_epochs( const Eigen::MatrixXd & pp , // posterior probabilities
				     const std::vector<std::string> & labels, // column labels
				     int ne_all , edf_t & edf ) // total number of epochs in EDF
{
  // output epoch-level results: most likely stage, PP, observed stage, flag for discordance, missing/unknown

  bool prior_staging = obs_stage.size() != 0 ;
  
  // epochs[] contains the codes of epochs actually present in the model/valid
  std::map<int,int> e2e;
  for (int i=0; i< epochs.size(); i++) e2e[ epochs[i] ] = i ;  

  for (int i=0;i<ne_all;i++)
    {

      int e = -1;

      if ( e2e.find( i ) != e2e.end() ) e = e2e[i];
      
      writer.epoch( edf.timeline.display_epoch( i ) );

      if ( e != -1 ) 
	{
	  
	  writer.value( "INC" , 1 );

	  double pp_nr = 0;
	  bool has_nr = false;
	  for (int j=0;j<labels.size();j++)
	    {
	      if ( labels[j] == "NR" ) has_nr = true;
	      if ( labels[j] == "N1" || labels[j] == "N2" || labels[j] == "N3" ) pp_nr += pp(e,j); 
	      writer.value( "PP_" + labels[j] , pp(e,j) );
	    }

	  // automatically aggregate N1+N2+N3 under the 5-class model (or whatever NREM stages are present)
	  if ( ! has_nr )
	    writer.value( "PP_NR" , pp_nr );
	
	  // most likely value
	  std::string predss = suds_t::max_inrow( pp.row(e) , labels );
	  writer.value( "PRED" , predss );

	  if ( prior_staging )
	    {
	      // discordance if prior/obs staging available
	      bool disc = obs_stage[i] != SUDS_UNKNOWN && predss !=  suds_t::str( obs_stage[i] ) ;
	      writer.value( "DISC" , disc );

	      // collapse 5->3 ?
	      if ( suds_t::n_stages == 5 )
		{
		  bool disc3 =  obs_stage[i] != SUDS_UNKNOWN && suds_t::NRW( predss ) != suds_t::NRW( suds_t::str( obs_stage[i] ) ) ;
		  writer.value( "DISC3" , disc3 );
		}
	      
	      writer.value( "PRIOR" ,  suds_t::str( obs_stage[i] ) );
	      
	      if ( suds_t::soap_mode == 2 ) 
		writer.value( "PROPOSAL" ,  y[e] );
	      
	    }
	}
      else
	{
	  writer.value( "INC" , 0 );
	  
	  // lookup from all stages
	  if ( prior_staging )
	    writer.value( "PRIOR" ,  suds_t::str( obs_stage[i] ) );	  
	}
      
    }

  writer.unepoch();
  
}

void suds_indiv_t::summarize_stage_durations( const Eigen::MatrixXd & pp , const std::vector<std::string> & labels, int ne_all , double epoch_sec )
{
  
  bool prior_staging = obs_stage.size() != 0 ;
 
  std::map<std::string,double> prd_dur;   // sum of PP
  std::map<std::string,double> prd2_dur;  // based on most likely
  std::map<std::string,double> obs_dur;   // obserevd  (if present)... but based on same epochs as used in the staging (i.e. removing some outliers) 

  std::map<int,int> e2e;
  for (int i=0; i<epochs.size(); i++) e2e[ epochs[i] ] = i ;  
  

  //
  // Get labels -> slots
  //
  
  int n1_slot, n2_slot, n3_slot, nr_slot, rem_slot, wake_slot;
  n1_slot = n2_slot = n3_slot = nr_slot = rem_slot = wake_slot = -1;

  for (int i=0;i< labels.size(); i++)
    {
      if      ( labels[i] == "N1" ) n1_slot = i;
      else if ( labels[i] == "N2" ) n2_slot = i;
      else if ( labels[i] == "N3" ) n3_slot = i;
      else if ( labels[i] == "NR" ) nr_slot = i;
      else if ( labels[i] == "R" ) rem_slot = i;
      else if ( labels[i] == "REM" ) rem_slot = i;
      else if ( labels[i] == "W" ) wake_slot = i;      
    }

  // for output, below
  std::string rem_out = suds_t::n_stages == 5 ? "REM" : "R";
  
  double unknown = 0;
  
  //
  // Aggregate over epochs
  //

  
  for (int i=0;i<ne_all;i++)
    {
      
      int e = -1;
      
      if ( e2e.find( i ) != e2e.end() ) e = e2e[i];
      
      if ( e != -1 ) 
	{
	
	  // most likely value
	  std::string predss = suds_t::max_inrow( pp.row(e) , labels );

	  // track stage duration (based on probabilistic calls)
	  // nb. we do not assume all five/three stages are present here

	  if ( n1_slot != -1 ) prd_dur[ "N1" ]  += pp(e,n1_slot) * epoch_sec ;
	  if ( n2_slot != -1 ) prd_dur[ "N2" ]  += pp(e,n2_slot) * epoch_sec ;
	  if ( n3_slot != -1 ) prd_dur[ "N3" ]  += pp(e,n3_slot) * epoch_sec ;
	  if ( nr_slot != -1 ) prd_dur[ "NR" ]  += pp(e,nr_slot) * epoch_sec ;
	  if ( rem_slot != -1 ) prd_dur[ rem_out ]  += pp(e,rem_slot) * epoch_sec ;
	  if ( wake_slot != -1 ) prd_dur[ "W" ]  += pp(e,wake_slot) * epoch_sec ;
	  
	  // duration based on MAP estimate
	  prd2_dur[ predss ] += epoch_sec;

	  // comparable OBS duration
	  if ( prior_staging )	    
	    obs_dur[ suds_t::str( obs_stage[i] ) ] += epoch_sec; 

	}
      else
	{
	  // track extent of 'bad' epochs
	  unknown += epoch_sec;
	}

    }
  
  //
  // Report stage durations (in minutes)
  //

  for (int s=0; s < suds_t::labels.size(); s++) 
    {
      writer.level( suds_t::labels[ s ] , globals::stage_strat );      
      writer.value( "DUR_PRD" , prd_dur[ suds_t::labels[ s ] ] / 60.0 );

      // alternate estimates, based on most likely predicted epoch
      if ( suds_t::verbose )
	writer.value( "DUR_PRD2" , prd_dur[ suds_t::labels[ s ] ] / 60.0 );

    }  
  
  // unknown/missed epochs
  writer.level( suds_t::str( SUDS_UNKNOWN ) , globals::stage_strat );
  writer.value( "DUR_OBS" , unknown / 60.0 );

  // and done
  writer.unlevel( globals::stage_strat );


  //
  // estimates of observed stage duration (based on comparable set of epochs)
  //

  if ( prior_staging )
    {
      std::map<std::string,double>::const_iterator ss = obs_dur.begin();
      while ( ss != obs_dur.end() )
	{
	  writer.level( ss->first , globals::stage_strat );
	  writer.value( "DUR_OBS" , ss->second / 60.0 );
	  ++ss;
	}
      writer.unlevel( globals::stage_strat );
    }
  
}



void suds_indiv_t::summarize_kappa( const std::vector<std::string> & prd , const bool to_console )
{
  
  if ( to_console )
    logger << std::fixed << std::setprecision(2);
  
  // original reporting (5 or 3 level)
  double kappa = MiscMath::kappa( prd , suds_t::str( obs_stage_valid ) , suds_t::str( SUDS_UNKNOWN ) );

  // accuracy, precision/recall, kappa:   nb. ordering: 'truth' first, then 'predicted' 
  double macro_f1 = 0 , macro_precision = 0 , macro_recall = 0;
  double wgt_f1 = 0 , wgt_precision = 0 , wgt_recall = 0 , mcc = 0;
  std::vector<double> precision, recall, f1;

  double acc = MiscMath::accuracy( suds_t::str( obs_stage_valid ) , prd ,
				   suds_t::str( SUDS_UNKNOWN ) , 
				   &suds_t::labels ,
				   &precision, &recall, &f1,
				   &macro_precision, &macro_recall, &macro_f1 ,
				   &wgt_precision, &wgt_recall, &wgt_f1 , &mcc);
  
  writer.value( "K" , kappa );
  writer.value( "ACC" , acc );

  writer.value( "F1" , macro_f1 );
  writer.value( "MCC" , mcc );
  writer.value( "PREC" , macro_precision );
  writer.value( "RECALL" , macro_recall );
  
  writer.value( "F1_WGT" , wgt_f1 );
  writer.value( "PREC_WGT" , wgt_precision );
  writer.value( "RECALL_WGT" , wgt_recall );

  for ( int l=0;l<suds_t::labels.size();l++)
    {
      writer.level( suds_t::labels[l] , globals::stage_strat );
      writer.value( "F1" , f1[l] );
      writer.value( "PREC" , precision[l] );
      writer.value( "RECALL" , recall[l] );
    }
  writer.unlevel( globals::stage_strat );
  
  if ( to_console ) 
    {
      logger << "\n  Confusion matrix: " << suds_t::n_stages
	     << "-level classification: kappa = " << kappa << ", acc = " << acc << ", MCC = " << mcc << "\n";
      suds_t::tabulate(  prd , suds_t::str( obs_stage_valid ) , true );
    }      
  
  // collapse 5->3?
  if ( suds_t::n_stages == 5 )
    {
      
      double kappa3 = MiscMath::kappa( suds_t::NRW( prd ) , suds_t::NRW( suds_t::str( obs_stage_valid ) ) , suds_t::str( SUDS_UNKNOWN ) );

      // nb. 'truth' / pred
      double macro_f1 = 0 , macro_precision = 0 , macro_recall = 0;
      double wgt_f1 = 0 , wgt_precision = 0 , wgt_recall = 0 , mcc = 0;
      std::vector<double> precision, recall, f1;
      std::vector<std::string> lab3 = { "NR" , "R" , "W" };
      
      double acc3 = MiscMath::accuracy( suds_t::NRW( suds_t::str( obs_stage_valid ) ) , suds_t::NRW( prd ) ,
					suds_t::str( SUDS_UNKNOWN ) , 
					&lab3 ,
					&precision, &recall, &f1,
					&macro_precision, &macro_recall, &macro_f1 ,
					&wgt_precision, &wgt_recall, &wgt_f1 , &mcc );
      
      writer.value( "K3" , kappa3 );
      writer.value( "ACC3" , acc3 );
      
      writer.value( "F13" , macro_f1 );
      writer.value( "MCC3" , mcc );
      writer.value( "PREC3", macro_precision );
      writer.value( "RECALL3" , macro_recall );

      if ( to_console )
	{
	  logger << "\n  Confusion matrix: 3-level classification: kappa = " << kappa3 << ", acc = " << acc3 << ", MCC = " << mcc << "\n";
	  suds_t::tabulate(  suds_t::NRW( prd ) , suds_t::NRW( suds_t:: str( obs_stage_valid ) ) , true );
	}
    }
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
