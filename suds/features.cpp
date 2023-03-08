
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
#include "pdc/pdc.h"

extern logger_t logger;

extern writer_t writer;



int suds_indiv_t::proc( edf_t & edf , param_t & param , bool is_trainer )
{

  //
  // helper struct to connect the proc modules
  //

  suds_helper_t helper (edf , param );


  //
  // Is this individual a trainer (i.e. with known stages) or no?
  //

  trainer = is_trainer;
  
  //
  // Initial/total number of components to extract
  // from PSC (although we may only retain nc2 <= nc
  // for the LDA)
  //
  
  nc = suds_t::nc;

  int rc = 0;


  //
  // Check required signals
  //

  rc = proc_check_channels( &helper );
  if ( rc == 0 ) return 0;
  
  //
  // For trainers, get the observed stages
  //

  rc = proc_extract_observed_stages( &helper ) ;
  if ( rc == 0 ) return 0;

  //
  // Build feature matrix given a model
  //
  
  rc = proc_build_feature_matrix( &helper );
  if ( rc == 0 ) return 0;  
  
  //
  // epoch-level QC (also performs an initial SVD) ( nge --> nve ) 
  //

  rc = proc_initial_svd_and_qc( &helper );
  if ( rc == 0 ) return 0;

  //
  // populate 'y'
  //
  
  rc = proc_class_labels( &helper );
  if ( rc == 0 ) return 0;


  //
  // re-do SVD on dataset w/ bad epochs removed
  //

  rc = proc_main_svd( &helper );
  if ( rc == 0 ) return 0;


  //
  // For SUDS trainers, drop epochs that are not well-classified (i.e. outliers in the current model)
  //
  
  rc = proc_prune_rows( &helper );
  if ( rc == 0 ) return 0;

  
  //
  // re-do (third time) main SVD on final dataset
  //

  rc = proc_main_svd( &helper );
  if ( rc == 0 ) return 0;

    
  //
  // For SUDS trainers, drop any components that do not track well w/ stage
  //
  
  rc = proc_prune_cols( &helper );
  if ( rc == 0 ) return 0;
 

  //
  // some final metrics
  //
  
  return proc_coda( &helper );
  
}


int suds_indiv_t::proc_check_channels( suds_helper_t * helper )
{


  //
  // Signals (and optionally, resampling)
  //
  
  helper->ns = suds_t::model.chs.size();

  // check that all model channels are also present
  // in the EDF
  std::vector<std::string> slabs;
  std::vector<int> slots;

  std::map<std::string,suds_channel_t>::const_iterator ss =  suds_t::model.chs.begin(); 
  while ( ss != suds_t::model.chs.end() )
    {
      int slot = helper->edf.header.signal( ss->first );
      if ( slot == -1 ) Helper::halt( "could not find " + ss->first );
      
      if ( helper->edf.header.is_annotation_channel( slot ) )
	Helper::halt( "cannot specificy annotation channel: " + ss->first );
      
      // need to resample?
      if ( helper->edf.header.sampling_freq( slot ) != ss->second.sr )
        dsptools::resample_channel( helper->edf, slot , ss->second.sr );
      
      // build signal_list_t
      helper->signals.add( slot , ss->first );
      
      ++ss;
    }
  
  return 1;
}



int suds_indiv_t::proc_extract_observed_stages( suds_helper_t * helper )
{


  //
  // Epoch data:
  //
  //   ne     total number of epochs
  //
  //   nge    num of epochs with 'valid' staging (i.e. no UNKNOWN, etc)
  //
  //   nve    of nge, number that are
  //               a) not statistical outliers for any components
  //               b) optionally, correctly self-classified 

  helper->ne = helper->edf.timeline.first_epoch();



  //
  // Get stage information (for trainers only)
  //
  
  helper->retained.resize( helper->ne , true );
  
  helper->has_prior_staging = false;

  // for SUDS trainers, always load observed stages
  // for SUDS target, load unless told not to

  // for SOAP mode, always load unless told not to ( here, 'target' is 'trainer' and so trainer is T)
  // but in this case, set obs_stage, just make it all missing (i.e. as trainers need to have this
  // populated

  if ( suds_t::soap_mode && suds_t::ignore_target_priors )
    {
      helper->has_prior_staging = false;
      // nb. this is set back to 'true'  after the following 
      // section, i.e. to do read from annotations, but we need to say we have stages
      // (albeit all unknown ones) for downstream code) 
      obs_stage.resize( helper->ne , SUDS_UNKNOWN );
    }
  else if ( trainer )
    {
      helper->edf.timeline.annotations.make_sleep_stage( helper->edf.timeline );
      
      if ( ! helper->edf.timeline.hypnogram.construct( &(helper->edf.timeline) , helper->param , false ) )
	{
	  if ( suds_t::soap_mode ) return 0; // okay to skip for SOAP
	  // but flag as major prob if a trainer
	  Helper::halt( "problem extracting stage information for trainer" );
	}
      
      // total number of epochs does not match?
      if ( helper->ne != helper->edf.timeline.hypnogram.stages.size() )
	Helper::halt( "problem extracting stage information for trainer" );
      
      helper->has_prior_staging = true;
      
    }
  else if ( ! suds_t::ignore_target_priors )
    {

      // for targets, manual/prior staging may exist, in which case we'll want to track it for comparison
      // unless we've been explicitly told to ignore it (ignore-prior --> suds_t::ignore_target_priors )

      helper->edf.timeline.annotations.make_sleep_stage( helper->edf.timeline );

      helper->has_prior_staging = helper->edf.timeline.hypnogram.construct( &(helper->edf.timeline) , helper->param , false ) ;
      
      if ( helper->has_prior_staging )
	{
	  // total number of epochs does not match?
	  if ( helper->ne != helper->edf.timeline.hypnogram.stages.size() )
	    Helper::halt( "problem extracting stage information for trainer" );
	}
      
    }

  // nb. overkill to have two staging classifications, but keep for now,
  // i.e. we may want to extend one type only

  //
  // number of good (retained) epochs
  //
  
  helper->nge = 0 ;

  if ( helper->has_prior_staging )
    {
      
      obs_stage.resize( helper->ne , SUDS_UNKNOWN );
      
      for (int ss=0; ss < helper->ne ; ss++ )
	{
	  if ( helper->edf.timeline.hypnogram.stages[ ss ] == UNSCORED
	       || helper->edf.timeline.hypnogram.stages[ ss ] == MOVEMENT
	       || helper->edf.timeline.hypnogram.stages[ ss ] == UNKNOWN )
	    obs_stage[ss] = SUDS_UNKNOWN;
	  
	  else if ( helper->edf.timeline.hypnogram.stages[ ss ] == LIGHTS_ON )
	    obs_stage[ss] = SUDS_LIGHTS;
	  
	  else if ( helper->edf.timeline.hypnogram.stages[ ss ] == WAKE )
	    obs_stage[ss] = SUDS_WAKE;

	  else if ( helper->edf.timeline.hypnogram.stages[ ss ] == NREM1 )
	    obs_stage[ss] = suds_t::n_stages == 3 ? SUDS_NR : SUDS_N1;

	  else if ( helper->edf.timeline.hypnogram.stages[ ss ] == NREM2 )
	    obs_stage[ss] = suds_t::n_stages == 3 ? SUDS_NR : SUDS_N2;

	  else if ( helper->edf.timeline.hypnogram.stages[ ss ] == NREM3
		    || helper->edf.timeline.hypnogram.stages[ ss ] == NREM4 )
	    obs_stage[ss] = suds_t::n_stages == 3 ? SUDS_NR : SUDS_N3;
	  
	  else if ( helper->edf.timeline.hypnogram.stages[ ss ] == REM )
	    obs_stage[ss] = SUDS_REM;

	  
	  // expand retained class to exclude unknown staging info
	  // note: even for targets, this means we will not try to guess epochs
	  // that have an UNKNOWN value for the 

	  // note: however, if in 'self-classification' mode
	  // (i.e. running SOAP), then we allow these, i.e. to enable
	  // 'auto-completion' staging (i.e. if some small proportion
	  // is staged, and we want to complete the rest); we do not do this for L
	  // epochs however, they should be completely ignored
	       
	  // lda_t will ignore '?' in fitting the model; i.e. so we can process
	  // just the 'full' dataset (as will be used in the prediction part) 
	  // and that way we are sure that we are doing the SOAP part right

	  if ( suds_t::soap_mode )
	    {
	      if ( obs_stage[ss] == SUDS_LIGHTS )
                helper->retained[ss] = false;
	      else
		++helper->nge;
	    }
	  else
	    {	      
	      if ( obs_stage[ss] == SUDS_UNKNOWN ||  obs_stage[ss] == SUDS_LIGHTS ) 
		helper->retained[ss] = false; 
	      else ++helper->nge;
	    }

	  // next epoch	 
	}


      //
      // trim leading/trailing wake epochs? (only for SUDS?)
      //
      
      //if ( suds_t::soap_mode == 0 && suds_t::trim_wake_epochs >= 0 ) 
      if ( suds_t::trim_wake_epochs >= 0 ) 
	{
	  int first_sleep = -1;
	  for (int ss=0; ss < helper->ne; ss++)
	    {
	      if ( obs_stage[ss] == SUDS_N1 ||
		   obs_stage[ss] == SUDS_N2 ||
		   obs_stage[ss] == SUDS_N3 ||
		   obs_stage[ss] == SUDS_NR ||
		   obs_stage[ss] == SUDS_REM ) 
		{
		  first_sleep = ss;
		  break;
		}
	    }
	  
	  int last_sleep = helper->ne - 1;
	  for (int ss = helper->ne - 1; ss >=0 ; ss--)
	    {
	      if ( obs_stage[ss] == SUDS_N1 ||
		   obs_stage[ss] == SUDS_N2 ||
		   obs_stage[ss] == SUDS_N3 ||
		   obs_stage[ss] == SUDS_NR ||
		   obs_stage[ss] == SUDS_REM ) 
		{
		  last_sleep = ss;
		  break;
		}
	    }
	  
	  // trim front
	  if ( first_sleep > 0 ) 
	    {
	      //         *
	      // 0 1 2 3 4
	      // if allow 2
	      // X X Y Y S

	      first_sleep -= suds_t::trim_wake_epochs + 1 ; 	      
	      int t = 0;
	      // note, inclusive counting up to X
	      for (int ss=0; ss<= first_sleep; ss++)
		{
		  obs_stage[ss] = SUDS_UNKNOWN; 
		  helper->retained[ss] = false;
		  --helper->nge;
		  ++t;
		}
	      if ( t ) logger << "  trimmed " << t << " leading wake epochs\n";
	      helper->trimmed += t;
	    }
	  
	  // trim end
	  if ( last_sleep < helper->ne - 1 ) 
	    {
	      // * *               
	      // 4 5 6 7 8 9
	      //         X X
	      
	      last_sleep += suds_t::trim_wake_epochs + 1 ;
	      int t=0;
	      for (int ss= helper->ne - 1 ; ss >= last_sleep; ss--)
		{
		  obs_stage[ss] = SUDS_UNKNOWN;
                  helper->retained[ss] = false;
                  --helper->nge;
		  ++t;
		}
	      if ( t ) logger << "  trimmed " << t << " trailing wake epochs\n";
	      helper->trimmed += t;
	    }
	  
	} // end of wake trimming option


    }
  else 
    {

      //
      // for target individuals without staging, include all epochs
      //
      
      helper->nge = helper->ne;
      
    }

  
  //
  // See note above
  //
  
  if ( suds_t::soap_mode && suds_t::ignore_target_priors )
    {
      helper->has_prior_staging = true;
    }
  
  return 1;

}




int suds_indiv_t::proc_build_feature_matrix( suds_helper_t * helper )
{
  
  
  //
  // PSD (Welch) parameters 
  //

  double fft_segment_size = helper->param.has( "segment-sec" ) 
    ? helper->param.requires_dbl( "segment-sec" ) : 4 ;
  
  double fft_segment_overlap = helper->param.has( "segment-overlap" ) 
    ? helper->param.requires_dbl( "segment-overlap" ) : 2 ;
  
  if ( helper->edf.timeline.epoch_length() < fft_segment_size )
    Helper::halt( "Welch segment size (segment-sec) cannot be greater than epoch length" );
  
  window_function_t window_function = WINDOW_TUKEY50;	   
  if      ( helper->param.has( "no-window" ) ) window_function = WINDOW_NONE;
  else if ( helper->param.has( "hann" ) ) window_function = WINDOW_HANN;
  else if ( helper->param.has( "hamming" ) ) window_function = WINDOW_HAMMING;
  else if ( helper->param.has( "tukey50" ) ) window_function = WINDOW_TUKEY50;


  logger << "  applying Welch with " << fft_segment_size << "s segments ("
	 << fft_segment_overlap << "s overlap), using "
	 << ( suds_t::use_seg_median ? "median" : "mean" )
	 << " over segments\n";  
  


  //
  // Size feature matrix X
  //

  nf = suds_t::nf;

  X.resize( helper->nge , nf );

  logger << "  expecting " << nf << " features (for " << helper->nge << " epochs) and " << helper->ns << " channels\n";

  
  //
  // for QC, estimate Hjorth parameters over
  // epochs, for each signal
  //
  
  h1 = Eigen::MatrixXd::Zero( helper->nge , helper->ns );
  h2 = Eigen::MatrixXd::Zero( helper->nge , helper->ns );
  h3 = Eigen::MatrixXd::Zero( helper->nge , helper->ns );
  
  //
  // Track bad epochs 
  //
  
  std::set<int> bad_epochs; // of 0..nge-1 encoding 
  

  //
  // iterate over (retained) epochs
  //
  
  int en = 0 , en_good = 0;
  
  helper->edf.timeline.first_epoch();
  
  epochs.clear();
  
  while ( 1 ) 
    {
      
      //
      // select epoch
      //


      int epoch = helper->edf.timeline.next_epoch();      	  
      
      if ( epoch == -1 ) break;
      
      if ( en == helper->ne ) Helper::halt( "internal error: over-counted epochs" );
      
      
      //
      // retained? if not, skip
      //
      
      if ( ! helper->retained[ en ] ) 
   	{
   	  ++en;
   	  continue;
   	}
      
      //
      // Process this epoch, signal-by-signal, then feature-spec by feature-spec.
      //
      
      interval_t interval = helper->edf.timeline.epoch( epoch );
      
      //
      // is this a bad epoch?
      //
      
      bool bad_epoch = false;

      //
      // Iterate over signals
      //

      for (int s = 0 ; s < helper->ns; s++ )
	{

	  //
	  // Skip if a bad epoch has been flagged already for a previous channel
	  //

	  if ( bad_epoch ) continue;

	  //
	  // Get data
	  //

	  helper->siglab = helper->signals.label(s) ;
	  
	  slice_t slice( helper->edf , helper->signals(s) , interval );
	  
	  const int sr = helper->edf.header.sampling_freq( helper->signals(s) ); 
	  
	  //
	  // get data & mean-center
	  //
	  
	  std::vector<double> * d = slice.nonconst_pdata();
	  
	  
	  //
	  // mean centre epoch
	  //
	  
	  double mean = MiscMath::centre( d );
	  
	  
	  //
	  // extract these channel-specific features 
	  //     
	  
	  const bool do_mean = suds_t::model.has( suds_feature_t::SUDS_MEAN , helper->siglab ) ;
	  
	  const bool do_spectral =
	    suds_t::model.has( suds_feature_t::SUDS_LOGPSD , helper->siglab ) ||
	    suds_t::model.has( suds_feature_t::SUDS_RELPSD , helper->siglab ) ||
	    suds_t::model.has( suds_feature_t::SUDS_SLOPE , helper->siglab ) ||
	    suds_t::model.has( suds_feature_t::SUDS_CVPSD , helper->siglab );
	  
	  const bool do_skew = suds_t::model.has( suds_feature_t::SUDS_SKEW , helper->siglab );
	  
	  const bool do_kurt = suds_t::model.has( suds_feature_t::SUDS_KURTOSIS , helper->siglab );
	  
	  const bool do_hjorth = suds_t::model.has( suds_feature_t::SUDS_HJORTH , helper->siglab );
	  
	  const bool do_pe = suds_t::model.has( suds_feature_t::SUDS_PE , helper->siglab );
	  
	  const bool do_pfd = suds_t::model.has( suds_feature_t::SUDS_FD , helper->siglab );
	  	  	  

	  //
	  // PSD (Welch)
	  //
	  
	  if ( do_spectral )
	    {
	      
	      //
	      // Get spectrum via Welch
	      //
	      
	      const double overlap_sec = fft_segment_overlap;
	      const double segment_sec  = fft_segment_size;
	      const int total_points = d->size();
	      const int segment_points = segment_sec * sr;
	      const int noverlap_points  = overlap_sec * sr;
	      
	      // implied number of segments
	      int noverlap_segments = floor( ( total_points - noverlap_points) 
					     / (double)( segment_points - noverlap_points ) );
	      
	      // also calculate SD over segments for this channel?
	      const bool get_segment_sd = suds_t::model.has( suds_feature_t::SUDS_CVPSD , helper->siglab ) ;
	      
	      PWELCH pwelch( *d , 
			     sr , 
			     segment_sec , 
			     noverlap_segments , 
			     window_function ,
			     suds_t::use_seg_median ,
			     get_segment_sd );
	      
	      // using bin_t, 1 means no binning
	      bin_t bin( suds_t::lwr , suds_t::upr , 1 ); 
	      bin.bin( pwelch.freq , pwelch.psd );	      
	      
	      //
	      // check for zero power values in the 0.5 to 45 Hz range, and flag if so
	      //  -- we will not include this epoch
	      //
	      
	      for ( int i = 0 ; i < bin.bfa.size() ; i++ )
		{
		  if ( bin.bfb[i] > suds_t::upr ) break;
		  if ( bin.bspec[i] <= 0 && bin.bfa[i] >= suds_t::lwr ) 
		    {
		      bad_epoch  = true;		       
		      bin.bspec[i] = 1e-4 ; // set to -40dB as a fudge		   
		    }
		}

	      // track that this is bad / to be removed below?
	      if ( bad_epoch ) bad_epochs.insert( en_good );
	      
	      //
	      // log-PSD?
	      //
	      
	      if ( suds_t::model.has( suds_feature_t::SUDS_LOGPSD , helper->siglab ) && ! bad_epoch )
		{
		  
		  std::vector<int> cols = suds_t::model.cols( suds_feature_t::SUDS_LOGPSD , helper->siglab ) ;
		  const int ncols = cols.size();
		  
		  // this *should* map exactly onto the number of bins between the lwr and upr bounds
		  
		  suds_spec_t spec = suds_t::model.fcmap[ suds_feature_t::SUDS_LOGPSD ][ helper->siglab ];
		  
		  // these have been checked and will be present/valid 
		  const double lwr = spec.arg[ "lwr" ];
		  const double upr = spec.arg[ "upr" ];
		  
		  int b = 0;
		  
		  for ( int i = 0 ; i < bin.bfa.size() ; i++ )
		    {
		      if (  bin.bfa[i] >= lwr && bin.bfa[i] <= upr )
			{
			  //std::cout << " b " << b << " " << bin.bfa[i]  << " " << ncols <<  " " << lwr << " " << upr << "\n";
			  if ( b == ncols ) Helper::halt( "internal error... bad sizes for SPEC" );
			  
			  // save log-scaled power
			  X( en_good , cols[b] ) = 10*log10( bin.bspec[i] ) ; 
			  
			  // next feature column
			  ++b;			   
			}
		    }
		}
	      
	      
	      //
	      // rel-PSD?
	      //
	      
	      if ( suds_t::model.has( suds_feature_t::SUDS_RELPSD , helper->siglab ) && ! bad_epoch )
		{
		  std::vector<int> cols = suds_t::model.cols( suds_feature_t::SUDS_RELPSD , helper->siglab ) ;
		  const int ncols = cols.size();
		  
		  suds_spec_t spec = suds_t::model.fcmap[ suds_feature_t::SUDS_RELPSD ][ helper->siglab ];
		  const double lwr = spec.arg[ "lwr" ];
		  const double upr = spec.arg[ "upr" ];
		  
		  const double zlwr = spec.arg[ "z-lwr" ];
		  const double zupr = spec.arg[ "z-upr" ];
		  
		  // get normalization factor
		  double norm = 0;
		  for ( int i = 0 ; i < bin.bfa.size() ; i++ )
		    {
		      if ( bin.bfa[i] > zupr ) break;		       
		      if ( bin.bfa[i] >= zlwr ) norm += bin.bspec[i] ;
		    }
		  // sanity check
		  if ( norm == 0 )
		    {
		      bad_epoch = true;
		      norm = 1e-4;
		    }
		  
		  int b = 0;				   
		  for ( int i = 0 ; i < bin.bfa.size() ; i++ )
		    {
		      if (  bin.bfa[i] >= lwr && bin.bfa[i] <= upr )
			{
			  if ( b == ncols ) Helper::halt( "internal error... bad sizes for VSPEC" );
			  X( en_good , cols[b] ) = log( bin.bspec[i] / norm ) ; 
			  ++b;			   
			}
		    }
		}
	      
	      
	      //
	      // cv-PSD?
	      //
	      
	      if ( suds_t::model.has( suds_feature_t::SUDS_CVPSD , helper->siglab ) && ! bad_epoch )
		{
		  
		  std::vector<int> cols = suds_t::model.cols( suds_feature_t::SUDS_CVPSD , helper->siglab ) ;
		  const int ncols = cols.size();
		  
		  suds_spec_t spec = suds_t::model.fcmap[ suds_feature_t::SUDS_CVPSD ][ helper->siglab ];
		  const double lwr = spec.arg[ "lwr" ];
		  const double upr = spec.arg[ "upr" ];
		  
		  int b = 0;
		  
		  for ( int i = 0 ; i < pwelch.freq.size() ; i++ )
		    {
		      if (  pwelch.freq[i] >= lwr && pwelch.freq[i] <= upr )
			{
			  if ( b == ncols ) Helper::halt( "internal error... bad sizes for VSPEC" );
			  
			  // save CV of PSD
			  X( en_good , cols[b] ) = pwelch.psdsd[i];
			  
			  // next feature column
			  ++b;			   
			}
		    }
		  
		}
	      
	      //
	      // Spectral slope?
	      //
	      
	      if ( suds_t::model.has( suds_feature_t::SUDS_SLOPE , helper->siglab ) && ! bad_epoch )
		{
		   
		  double bslope = 0, bn = 0;
		  
		  bool okay = spectral_slope_helper( pwelch.psd , 
						     pwelch.freq , 
						     suds_t::slope_range ,
						     suds_t::slope_th , 
						     false ,  // do not output value
						     &bslope , &bn ); 
		  if ( ! okay ) bad_epoch = true;
		  
		  // will be exactly size == 1 
		  std::vector<int> cols = suds_t::model.cols( suds_feature_t::SUDS_SLOPE , helper->siglab ) ;
		  
		  // save slope
		  X( en_good , cols[0] ) = bslope;
		   
		}
	      
	    }
	  

	   //
	   // Time domain features
	   //

	   if ( do_mean && ! bad_epoch )
	     {
	       std::vector<int> cols = suds_t::model.cols( suds_feature_t::SUDS_MEAN , helper->siglab ) ;
	       X( en_good , cols[0] ) = mean; // calculated above when mean-centering
	     }
	   
	   if ( do_skew && ! bad_epoch )
	     {
	       std::vector<int> cols = suds_t::model.cols( suds_feature_t::SUDS_SKEW , helper->siglab ) ;
	       X( en_good , cols[0] ) = MiscMath::skewness( *d , 0 , MiscMath::sdev( *d , 0 ) );
	     }

	   if ( do_kurt && ! bad_epoch )
	     {
	       std::vector<int> cols = suds_t::model.cols( suds_feature_t::SUDS_KURTOSIS , helper->siglab ) ;
	       X( en_good , cols[0] ) = MiscMath::kurtosis0( *d ); // assumes mean-centered
	     }
	   
	   // fractal dimension
	   if ( do_pfd && ! bad_epoch )
	     {
	       std::vector<int> cols = suds_t::model.cols( suds_feature_t::SUDS_FD , helper->siglab ) ;
	       X( en_good , cols[0] ) = MiscMath::petrosian_FD( *d );
	     }
	   
	   // permutation entropy
	   if ( do_pe && ! bad_epoch )
	     {
	       std::vector<int> cols = suds_t::model.cols( suds_feature_t::SUDS_PE , helper->siglab ) ;
	       
	       int sum1 = 1;
	       std::vector<double> pd3 = pdc_t::calc_pd( *d , 3 , 1 , &sum1 );
	       std::vector<double> pd4 = pdc_t::calc_pd( *d , 4 , 1 , &sum1 );
	       std::vector<double> pd5 = pdc_t::calc_pd( *d , 5 , 1 , &sum1 );
	       std::vector<double> pd6 = pdc_t::calc_pd( *d , 6 , 1 , &sum1 );
	       std::vector<double> pd7 = pdc_t::calc_pd( *d , 7 , 1 , &sum1 );
	       
	       X( en_good , cols[0] ) = pdc_t::permutation_entropy( pd3 );
	       X( en_good , cols[1] ) = pdc_t::permutation_entropy( pd4 );
	       X( en_good , cols[2] ) = pdc_t::permutation_entropy( pd5 );
	       X( en_good , cols[3] ) = pdc_t::permutation_entropy( pd6 );
	       X( en_good , cols[4] ) = pdc_t::permutation_entropy( pd7 );
	       
	     }

	   //
	   // Hjorth parameters: these are always calculated for (trainer) QC, but 
	   // they may also be added as explicit features
	   //

	   if ( ! bad_epoch )
	     {
	       double activity = 0 , mobility = 0 , complexity = 0;
	       MiscMath::hjorth( d , &activity , &mobility , &complexity );

	       h1( en_good , s ) = activity ;
	       h2( en_good , s ) = mobility ;
	       h3( en_good , s ) = complexity ; 
	       
	       if ( do_hjorth ) // only take H2 and H3
		 {
		   std::vector<int> cols = suds_t::model.cols( suds_feature_t::SUDS_HJORTH , helper->siglab ) ;
		   X( en_good , cols[0] ) = mobility;
		   X( en_good , cols[1] ) = complexity;
		 }
	     }

	   //
	   // Next signal
	   //
	}


      //
      // track that this was a bad epoch (for at least one signal/metric)
      //
      
      if ( bad_epoch ) 
	bad_epochs.insert( en_good ) ;
      
      //
      // increase epoch-number
      //
      
      ++en;
      ++en_good;
      epochs.push_back( epoch );
      
    } // next epoch
  
  
   // --------------------------------------------------------------------------------
   //
   // Finished adding in signal-based features, epoch-by-epoch
   //
   // --------------------------------------------------------------------------------

   
   // --------------------------------------------------------------------------------
   //
   // Winsorize/rescale metrics before smoothing; this ignores any invariant (0s) cols
   // i.e. time-tracks, smoothed cols not yet filled
   //
   // --------------------------------------------------------------------------------
      
   if ( suds_t::standardize_X )
     {
       if ( suds_t::robust_standardization )
	 {
	   logger << "  robust standardizing X";
	   if ( suds_t::winsor1 > 0 ) logger << ", winsorizing at " << suds_t::winsor1;
	   logger << "\n";

	   // okay to have blank cols at this point ( final true arg ) 
	   eigen_ops::robust_scale( X , true , true , suds_t::winsor1 , true );
	   
	   // if ( ! eigen_ops::robust_scale( X , true , true , suds_t::winsor1 , true ) )
	   //   {
	   //     logger << "  one or more features with no variability, quitting\n";
	   //     return 0;
	   //   }
	 }
       else
	 {
	   logger << "  standardizing X\n";

	   // allow invariant cols ( final true arg ) at this poitn
	   eigen_ops::scale( X , true , true , true ); 

	   // if ( ! eigen_ops::scale( X , true , true , true ) ) 
	   //   {
	   //     logger <<"  one or more features with no variability, quitting\n";
           //     return 0;
	   //   }
	   
	 }
     }


   // --------------------------------------------------------------------------------
   //
   // Final features: time-tracks
   //
   // --------------------------------------------------------------------------------

   if ( suds_t::model.has( suds_feature_t::SUDS_TIME , "." ) )
     {
       suds_spec_t spec = suds_t::model.fcmap[ suds_feature_t::SUDS_TIME ][ "." ];
       const int order = spec.arg[ "order" ];
       if ( order != 0 ) 
	 {
	   logger << "  adding " << order << " time-tracks\n";
	   Eigen::MatrixXd tt = suds_t::add_time_track( X.rows() , order );
	   
	   std::vector<int> cols = suds_t::model.cols( suds_feature_t::SUDS_TIME , "." ) ;
	   if ( cols.size() != order ) Helper::halt( "internal error in column spec. for TIME" );
	   for (int i=0; i<order; i++)
	     X.col( cols[i] ) = tt.col(i); 
	 }
     }

   
   // --------------------------------------------------------------------------------
   //
   // Temporal smoothing/denoising over epochs
   //
   // --------------------------------------------------------------------------------
   
   const bool do_smooth = suds_t::model.has( suds_feature_t::SUDS_SMOOTH , "." );
   const bool do_smooth2 = suds_t::model.has( suds_feature_t::SUDS_SMOOTH2 , "." );
   const bool do_denoise = suds_t::model.has( suds_feature_t::SUDS_DENOISE , "." );
   const bool do_denoise2 = suds_t::model.has( suds_feature_t::SUDS_DENOISE2 , "." );

   if ( do_smooth || do_smooth2 || do_denoise || do_denoise2 )
     {
       int chk = (int)do_smooth + (int)do_smooth2 + (int)do_denoise + (int)do_denoise2 ;
       if ( chk != 1 ) Helper::halt( "can only apply one of SMOOTH, SMOOTH2, DENOISE and DENOISE2" );

       suds_feature_t ftr = do_smooth ? suds_feature_t::SUDS_SMOOTH
	 : ( do_smooth2 ? suds_feature_t::SUDS_SMOOTH2
	     : ( do_denoise ? suds_feature_t::SUDS_DENOISE
		 : suds_feature_t::SUDS_DENOISE2  ) ) ;
       
       suds_spec_t spec = suds_t::model.fcmap[ ftr ][ "." ];
       
       // args
       const double lambda = do_denoise || do_denoise2 ? spec.arg[ "lambda" ] : 0 ;
       const int hwin = do_smooth || do_smooth2 ? spec.arg[ "half-window" ] : 0 ;
       const int fwin = 1 + 2 * hwin;

       // in duplicate mode::
       // any SMOOTH or DENOISE command should have copied all columns *prior* to that
       // argument;  therefore, length of cols should equal first element
       // orig:    0 1 2 3 4
       // smooth:            5 6 7 8 9 
       // // i.e. length 5, 

       // in replace mode:
       // orig:    0 1 2 3 4
       // smooth:  0 1 2 3 4
       
       std::vector<int> cols = suds_t::model.cols( ftr , "." ) ;
       if ( cols.size() == 0 ) Helper::halt( "invalid DENOISE/SMOOTH" );

       if ( ( do_denoise2 || do_smooth2 ) && cols.size() != cols[0] ) Helper::halt( "internal error in DENOISE/SMOOTH col spec." );
	   
       // original column is cols[q] - cols[0], or simply 'q' (as we doublted all cols up to this point)
       const int q = cols.size();

       logger << "  applying " << ( do_denoise || do_denoise2 ? "TV-denoising" : "smoothing" )
	      << " to " << q << " features, "
	      << ( do_denoise2 || do_smooth2 ? "duplicating" : "over-writing" )
	      << " original features\n";
      
       for (int i=0; i<q; i++)
	 {
	   if ( do_denoise || do_denoise2 )
	     {
	       const double sd = eigen_ops::sdev( X.col( i ) );
	       X.col( cols[i] ) = X.col( i );	   
	       dsptools::TV1D_denoise( X.col( cols[i] ) , lambda * sd );
	     }
	   else
	     {	       
	       // median-filter
	       // for (int i=0; i<q; i++)
	       //		 X.col( cols[i] ) = eigen_ops::median_filter( X.col( i ) , fwin );
	       
	       // moving average	       
	       X.col( cols[i] ) = eigen_ops::moving_average( X.col( i ) , fwin );
	     }
	 }
     }
   

   // --------------------------------------------------------------------------------
   //
   // Mean-centre all columns if not already done via standardization
   //
   // --------------------------------------------------------------------------------
   
   if ( ! suds_t::standardize_X )
     eigen_ops::scale( X , true , false );
   

   // --------------------------------------------------------------------------------
   //
   // All done creating the X feature matrix 
   //
   // --------------------------------------------------------------------------------
   
   logger << "  final feature matrix X: " << nf << " features over " << X.rows() << " epochs\n";



   // --------------------------------------------------------------------------------
   //
   // Apply feature weights
   //
   // --------------------------------------------------------------------------------

   if ( 0 ) 
     for (int c=0; c<suds_t::nf; c++)
       X.col(c) *= suds_t::model.W[c];
   
   
   // --------------------------------------------------------------------------------
   //
   // Dump?
   //
   // --------------------------------------------------------------------------------

   if ( 0 ) 
     {
       std::vector<std::string> ll = suds_t::model.labels();
       for (int ii=0;ii<ll.size();ii++)
	 std::cout << ( ii != 0 ? "\t" : "" ) << ll[ii] ;
       std::cout << "\n";
       
       std::cout << X << "\n";
     }

   return 1;
}



int suds_indiv_t::proc_initial_svd_and_qc( suds_helper_t * helper )
{
  
  // --------------------------------------------------------------------------------
  //
  // Get PSC initially (we look for outliers and then remove epochs, and redo the SVD)
  //
  // --------------------------------------------------------------------------------
  
  
  Eigen::BDCSVD<Eigen::MatrixXd> svd( X , Eigen::ComputeThinU | Eigen::ComputeThinV );
  U = svd.matrixU();
  U.conservativeResize( Eigen::NoChange , suds_t::nc );

  V = svd.matrixV();
  V.conservativeResize( Eigen::NoChange , suds_t::nc );

  W = svd.singularValues();
  W.conservativeResize( suds_t::nc );

  // std::cout << " NC = " << suds_t::nc << "\n";
  // std::cout << " sizes = " << U.rows() <<  " " << U.cols() << " ... " << V.rows() << " " << V.cols() << "\n";
  
  
  // --------------------------------------------------------------------------------
  //
  // Outliers in PSC space? 
  //
  // --------------------------------------------------------------------------------
   
  helper->valid.resize( helper->nge , true );
  
  // track reasons for exclusion
  std::set<int> nout_flat;
  std::set<int> nout_hjorth;
  std::set<int> nout_stat;
  std::set<int> nout_tot;

   
   //
   // Exclusions based on H==0 parameters
   //
   
   for ( int s=0;s<helper->ns;s++)
     {
       for (int i=0;i<helper->nge;i++) 
	 {
	   if      ( h1( i , s ) < 1e-8 ) { helper->valid[i] = false; nout_flat.insert(i); }
	   else if ( h2( i , s ) < 1e-8 ) { helper->valid[i] = false; nout_flat.insert(i) ; }
	   else if ( h3( i , s ) < 1e-8 ) { helper->valid[i] = false; nout_flat.insert(i) ; }
	 }
     }
   

   // --------------------------------------------------------------------------------
   //
   // For targets only, threshold epochs based on per-signal Hjorth
   // from trainers
   // 
   // --------------------------------------------------------------------------------
   
   if ( ! trainer )
     {
       logger << "  removing epochs +/-" << suds_t::hjorth_outlier_th << " SD units from Hjorth parameter trainer means\n";
 
       for ( int s=0;s<helper->ns;s++)
	 {
   	  for (int i=0;i<helper->nge;i++)
   	    {
   	      if ( h1( i, s ) <= suds_t::hjorth1_lwr95[s] || h1(i,s) >= suds_t::hjorth1_upr95[s] ||
		   h2( i, s ) <= suds_t::hjorth2_lwr95[s] || h2(i,s) >= suds_t::hjorth2_upr95[s] ||
		   h3( i, s ) <= suds_t::hjorth3_lwr95[s] || h3(i,s) >= suds_t::hjorth3_upr95[s] )
		{
		  helper->valid[i] = false;
		  nout_hjorth.insert(i);
		}
   	    }
   	}
     }


   // --------------------------------------------------------------------------------
   //
   // Component-based epoch-outlier removal (after removing flat lines)
   //
   // --------------------------------------------------------------------------------
   
   for (int o=0;o< suds_t::outlier_ths.size();o++)
     {

       logger << "  ";
       if ( o != 0 ) logger << "(repeatedly) ";
       logger << "removing epochs +/-" << suds_t::outlier_ths[o] << " from U means\n";
       
       for ( int j=0;j<nc;j++)
	 {
	   std::vector<double> x;
	   for (int i=0;i<helper->nge;i++) if ( helper->valid[i] ) x.push_back( U(i,j) );
	   if ( x.size() < 2 ) Helper::halt( "no epochs left" );
	   double mean = MiscMath::mean( x );
	   double sd = MiscMath::sdev( x , mean );
	   double lwr = mean - suds_t::outlier_ths[o] * sd;
	   double upr = mean + suds_t::outlier_ths[o] * sd;
	   int c = 0;
	   for (int i=0;i<helper->nge;i++)
	     {
	       if ( helper->valid[i] )
		 {
		   if ( x[c] < lwr || x[c] > upr ) { helper->valid[i] = false; nout_stat.insert(i); } 
		   ++c;
		 }
	     }
	 }
     }
      
   
   // --------------------------------------------------------------------------------
   //
   // Summarize dropped epochs and remove 
   //
   // --------------------------------------------------------------------------------
   
   int included = 0;

   for (int i=0;i<helper->nge;i++)
     if ( helper->valid[i] ) ++included;
   
   logger << "  of " << helper->ne << " total epochs, valid staging for " << helper->nge
          << ", and of those " << included << " passed outlier removal\n";
      
   std::set<int>::const_iterator oo = nout_flat.begin();
   while ( oo != nout_flat.end() ) { nout_tot.insert( *oo ); ++oo; } 
   oo = nout_hjorth.begin();
   while ( oo != nout_hjorth.end() ) { nout_tot.insert( *oo ); ++oo; } 
   oo = nout_stat.begin();
   while ( oo != nout_stat.end() ) { nout_tot.insert( *oo ); ++oo; } 
   
   logger << "  outlier counts: flat, Hjorth, components, trimmed -> total : "
   	 << nout_flat.size() << ", "
	  << nout_hjorth.size() << ", "
	  << nout_stat.size() << ", "
	  << helper->trimmed << " -> "
	  << nout_tot.size() + helper->trimmed << "\n";
   

   // --------------------------------------------------------------------------------
   //
   // Check we have enough data left
   // 
   // --------------------------------------------------------------------------------
   
   if ( included <= 20 ) 
     {
       logger << "  fewer than 20 epochs left after pruning... quitting\n";
       return 0;
     }


   // --------------------------------------------------------------------------------
   //
   // Remove bad epochs and repeat (SVD and smoothing)
   //
   // --------------------------------------------------------------------------------
   
   // nve = number of valid epochs ( ne > nge > nve ) 
  
   nve = included;

   Eigen::MatrixXd X2 = X;
   X.resize( nve , nf );
   std::vector<int> epochs2 = epochs;
   epochs.clear();
   
   int r = 0;
   for (int i=0;i<X2.rows() ; i++)
     {      
       if ( helper->valid[i] )
	 {	   
	   for (int j=0;j<nf;j++)
	      X(r,j) = X2(i,j);
	   
	   epochs.push_back( epochs2[i] );
	   
	   ++r;
	 }
     }
   

   // only retain nve obs labels from obs_stage[ne] originals

   if ( helper->has_prior_staging )
     {
       obs_stage_valid.clear();
    
       r = 0;
       for (int i=0;i<helper->ne;i++)
   	{
   	  if ( helper->retained[i] )
   	    {
   	      if ( helper->valid[r] )
   		obs_stage_valid.push_back( obs_stage[ i ] );
   	      ++r;
   	    }
   	}
     }

   
   //
   // splice out bad epochs for Hjorth parameters
   //

   Eigen::MatrixXd hh1 = h1;
   Eigen::MatrixXd hh2 = h2;
   Eigen::MatrixXd hh3 = h3;
   h1.resize( nve , helper->ns );
   h2.resize( nve , helper->ns ); 
   h3.resize( nve , helper->ns );

   for (int s=0;s<helper->ns;s++)
     {
       int r = 0;      
       for (int i=0; i < helper->valid.size(); i++ )
	 if ( helper->valid[i] )
	   {
	     h1(r,s) = hh1(i,s);
	     h2(r,s) = hh2(i,s);
	     h3(r,s) = hh3(i,s);
	     ++r;
	   }
     }
   
   
   
   // --------------------------------------------------------------------------------
   //
   // Rescale PSD?
   //
   // --------------------------------------------------------------------------------

   if ( suds_t::standardize_X )
     {
       
       if ( suds_t::robust_standardization )
	 {
	   logger << "  robust re-standardizing X after removing bad epochs\n";
	   // nb. not repeating winsorization of X here
	   
    	  if ( ! eigen_ops::robust_scale( X , true , true , 0 ) ) 
    	    {
    	      logger << "  one or more features with no variability, quitting\n";
    	      return 0;
    	    }

    	}
       else
	 {      
	   logger << "  re-standardizing X after removing bad epochs\n";	  
	   
	   if ( ! eigen_ops::scale( X , true , true ) )
	     {
	       logger << "  one or more features with no variability, quitting\n";
	       return 0;
	     }
	   
	 }
     }
   else // just mean-center X
     {
       eigen_ops::scale( X , true , false );
     }



   // --------------------------------------------------------------------------------
   //
   // Re-weight
   //
   // --------------------------------------------------------------------------------

   if ( 0 ) 
     for (int c=0; c<suds_t::nf; c++)
       X.col(c) *= suds_t::model.W[c];
   

   return 1;
}



int suds_indiv_t::proc_main_svd( suds_helper_t * helper )
{
  
  // --------------------------------------------------------------------------------
  //
  // Get PSC (post outlier removal)
  //
  // --------------------------------------------------------------------------------
  
  // W.resize( nbins ); 
  // V.resize( nbins , nbins );
  
  Eigen::BDCSVD<Eigen::MatrixXd> svd2( X , Eigen::ComputeThinU | Eigen::ComputeThinV );

  U = svd2.matrixU();
  V = svd2.matrixV();
  W = svd2.singularValues();

  U.conservativeResize( Eigen::NoChange , suds_t::nc );
  V.conservativeResize( Eigen::NoChange , suds_t::nc );
  W.conservativeResize( suds_t::nc );
  
  //  std::cout << " (MAIN) sizes = " << U.rows() <<  " " << U.cols() << " ... " << V.rows() << " " << V.cols() << "\n";

  //
  // Standardize U
  //
  
  if ( suds_t::standardize_U )
    {
      if ( suds_t::robust_standardization )
	{
	  logger << "  robust standardizing U\n";
	  
	  // no repeated winsorization here
	  if ( ! eigen_ops::robust_scale( U , true , true , 0 ) )
	    {
	      logger <<"  one or more features with no variability, quitting\n";
	      return 0;
	    }
	}
      else
	{
	  logger << "  standardizing U\n";
	  if ( ! eigen_ops::scale( U , true , true ) )
	    {
	      logger <<"  one or more features with no variability, quitting\n";
	      return 0;
	    }
	}  
    }
  
  return 1;
}



 
int suds_indiv_t::proc_prune_cols( suds_helper_t * helper )
{
  
  // --------------------------------------------------------------------------------
  //
  // For trainers, optionally only retain PSCs (or bands) that are significantly
  // associated with observed stage in this individual AND/OR do not have any stages
  // with greater within-stage variance than between stage variance
  //
  // --------------------------------------------------------------------------------
  
  if ( trainer && ( suds_t::required_comp_p < 1 || suds_t::betwithin_ratio > 0 ) && ! ( suds_t::soap_mode && suds_t::ignore_target_priors ) ) 
    {

      const bool do_anova = suds_t::required_comp_p < 1 ;
      const bool do_bw    = suds_t::betwithin_ratio > 0 ;
      
      //
      // pull out currently retained epochs
      //

      std::set<int> incl_comp;
      
      for (int j=0;j<nc;j++)
	{	  
	  // standardize column
	  Eigen::VectorXd c = U.col(j);
	  eigen_ops::scale( c , true , true );
	  
	  bool okay = true;
	  
	  writer.level( "PSC_" + Helper::int2str( j+1 ) , "VAR");
	  
	  if ( do_anova )
	    {
	      double pv = Statistics::anova( y  , eigen_ops::copy_vector( c ) );
	      writer.value( "PV", pv );	       
	      if ( pv < 0 || pv > suds_t::required_comp_p ) okay = false;
	    }
	  
	  // may have signif stage/group differences, but check that no one stage has a big variance difference also 
	  if ( do_bw )
	    {
	      // nb. c is standardized
	      double wb = eigen_ops::between_within_group_variance( y , c );
	      writer.value( "WMAX", wb );
	      if ( wb > suds_t::betwithin_ratio ) okay = false;	       
	    }
	  
	  if ( okay ) incl_comp.insert( j );
	  writer.value( "INC" , okay );
	  
	}
      
      writer.unlevel( "VAR" );

      // no usable components? --> no usable epochs...
      // quit out (this trainer will be ignored)
      
       if ( incl_comp.size() == 0 )
	 {
	   logger << "  0 p<" << suds_t::required_comp_p << " stage-associated components, bailing\n";
	   return 0; // i.e. no good epochs
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


   
   // --------------------------------------------------------------------------------
   //
   // Make variables for LDA: shrink down to 'nc' (if not already done by the above
   // component selection step)
   //
   // --------------------------------------------------------------------------------
   
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

   return 1;
}


int suds_indiv_t::proc_class_labels( suds_helper_t * helper )
{
  
  // --------------------------------------------------------------------------------
  //
  // make class labels ( trainer only )
  //
  // --------------------------------------------------------------------------------
  
  if ( trainer )
     {
       
       y.clear();
       
       int c = 0;
       for ( int i = 0 ; i < helper->ne ; i++ )
	 {
	   if ( helper->retained[i] )
	     {
	       if ( helper->valid[c] )
		 y.push_back( suds_t::str( obs_stage[i] ) );
	       ++c;
	     }
	 }

       counts.clear();
       for (int i=0;i<y.size();i++) counts[y[i]]++;

       if ( ! suds_t::ignore_target_priors )
	 {
	   std::map<std::string,int>::const_iterator cc = counts.begin();
	   logger << "  epoch counts:";
	   while ( cc != counts.end() )
	     {
	       logger << " " << cc->first << ":" << cc->second ;
	       ++cc;
	     }
	   logger << "\n";
	 }
     }
  
  return 1;
}






int suds_indiv_t::proc_prune_rows( suds_helper_t * helper ) 
{

  // this only applies for trainers

  if ( ! trainer ) return 1;
  
  // retained == size ne,  nge +pos elements
  // valid    == size nge, nve +pos elements
  // okay     == size nve, nve2 +pos elements

  int n0 = 0 , n1 = 0;
  for (int i=0;i<helper->valid.size(); i++) 
    if ( helper->valid[i] ) ++n1 ; else ++n0;

  std::vector<bool> okay( nve , true );
  
  // --------------------------------------------------------------------------------
  //
  // Self-classification (i.e. SOAP) to remove epochs that aren't well self-classified
  //  - possibly reject a trainer, if their SOAP kappa is poor
  //
  // --------------------------------------------------------------------------------
    
   if ( trainer && suds_t::self_classification )
     {
          
       // this returns the number of 'good' epochs 
       
       // std::cout << " nve " << nve << "\n";
       // std::cout << " y " << y.size() << "\n";
       // std::cout << " U " << U.rows() << " " << U.cols() << "\n";
       // std::cout << " epochs " << epochs.size() << "\n";

       int nve2 = self_classify( &okay );
       
       // std::cout << " done , nve2 = " << nve2 << "\n";
       // std::cout << " okay size = " << okay.size() << "\n";
       // std::cout << " retained size = " << helper->retained.size() << "\n";
       // std::cout << " valid size = " << helper->valid.size() << "\n";
       // std::cout << " nve2 = " << nve << " " << nve2 << "\n";
       
       if ( nve2 == 0 )
	 {
	   logger << "  trainer not valid based on self-classification thresholds\n";
	   return 0;
	 }
       
       //
       // track # of ambiguous epochs flagged here
       //
       
       helper->ambig = nve - nve2; 
       
     }
   
   
   // --------------------------------------------------------------------------------
   //
   // Impose a max number of epochs per stage? 
   //
   // --------------------------------------------------------------------------------
   
   if ( helper->has_prior_staging && suds_t::max_epoch_n != -1 )
     {
       
       // reset this
       helper->trimmed = 0;

       std::map<suds_stage_t,std::vector<int> > cnts;
              
       int cc = 0 , cc2 = 0;
       for (int i=0;i<helper->ne;i++)
	 {
	   if ( helper->retained[i] )
	     {
	       // track counts in valid index space
	       if ( helper->valid[cc] )
		 {
		   // and check was not flagged above
		   // store okay[] idx
		   if ( okay[ cc2 ] )
		     cnts[ obs_stage[ i ] ].push_back( cc2 );
		   ++cc2;
		 }
	       ++cc;
	     }
	 }
       
       std::map<suds_stage_t,std::vector<int> >::const_iterator qq = cnts.begin();
       while ( qq != cnts.end() )
	 {
	   if ( qq->second.size() > suds_t::max_epoch_n )
	     {
	       logger << "  reducing " << suds_t::str( qq->first ) 
		      << " from " << qq->second.size() 
		      << " to " << suds_t::max_epoch_n << " epochs\n";
	       int tot = qq->second.size();
	       int rem = tot - suds_t::max_epoch_n;
	       while ( rem ) 
		 {
		   int pick = CRandom::rand( tot );
		   if ( okay[ qq->second[ pick ] ] )
		     {
		       okay[ qq->second[ pick ] ] = false;
		       helper->trimmed++;
		       --rem;
		     }
		 }	      
	     }
	   ++qq;
	 }       
     }
   
   
   // --------------------------------------------------------------------------------
   //
   // Drop epochs that were either ambig and/or extraneous
   //
   // --------------------------------------------------------------------------------       

   
   int nve2 = 0;

   for (int i=0;i<okay.size();i++)
     if ( okay[i] ) ++nve2;
   
   // std::cout << " nve2 = " << nve2 << "\n";
   // std::cout << "  --> ambig, trim " << helper->ambig << " " << helper->trimmed << "\n";

   //   U X epochs  y  h1 h2  h3
   
   Eigen::MatrixXd UU = U;      
   U.resize( nve2 , nc );      
   
   Eigen::MatrixXd XX = X;      
   X.resize( nve2 , X.cols() );      
   
   std::vector<int> epochs2 = epochs;
   epochs.clear();
   
   std::vector<suds_stage_t> obs_stage_valid2 = obs_stage_valid;
   if ( helper->has_prior_staging )
     obs_stage_valid.clear();
   
   Eigen::MatrixXd hh1 = h1;
   h1.resize( nve2 , helper->ns );
   
   Eigen::MatrixXd hh2 = h2;
   h2.resize( nve2 , helper->ns );
   
   Eigen::MatrixXd hh3 = h3;
   h3.resize( nve2 , helper->ns );

   std::vector<std::string> yy = y;
   y.clear();
   
   int r = 0;
   for (int i=0;i < nve; i++)
     {      
       if ( okay[i] )
	 {
	   // U
	   for (int j=0;j<nc;j++)
	     U(r,j) = UU(i,j);
	   
	   // X original features
	   for (int j=0;j<nf;j++)
	     X(r,j) = XX(i,j);
	   
	   // Epoch tracking
	   epochs.push_back( epochs2[i] );
	   
	   // nb. take from already-pruned set, obs_stage_valid[] ) 
	   if ( helper->has_prior_staging )
	     obs_stage_valid.push_back( obs_stage_valid2[i] );	   

	   // labels (as text)
	   y.push_back( yy[i] );

	   // Hjorth (per signal)
	   for (int s=0;s<helper->ns;s++)
	     {
	       h1(r,s) = hh1(i,s);
	       h2(r,s) = hh2(i,s);
	       h3(r,s) = hh3(i,s);		  
	     }
	   
	   // next good epoch
	   ++r;
	 }
     }
        
   //
   // update nve
   //
   
   nve = nve2;
   
     
   
   
   // --------------------------------------------------------------------------------
   //
   // Recount stages
   //
   // --------------------------------------------------------------------------------
   
   if (  trainer && suds_t::self_classification ) 
     logger << "  removed " << helper->ambig << " epochs (posterior < " << suds_t::self_classification_prob << ")\n";
   
   if (  helper->has_prior_staging && suds_t::max_epoch_n != -1 ) 
     logger << "  removed " << helper->trimmed << " epochs to satisfy max-epoch requirements\n";
   
   report_epoch_counts( "final" );

   if ( ! suds_t::ignore_target_priors ) 
     logger << "  final count of valid epochs is " << nve << "\n";
      
   return nve <= 10 ? 0 : 1;
}

 
void suds_indiv_t::report_epoch_counts( const std::string & l )
{
  counts.clear();
  for (int i=0;i<y.size();i++) counts[y[i]]++;
  std::map<std::string,int>::const_iterator cc = counts.begin();

  if ( ! suds_t::ignore_target_priors )
    {
      if ( l == "" ) 
	logger << "  epoch counts:";
      else
	logger << "  " << l << " epoch counts:";
      while ( cc != counts.end() )
	{
	  logger << " " << cc->first << ":" << cc->second ;
	  ++cc;
	}
      logger << "\n";  
    }
}

int suds_indiv_t::proc_coda( suds_helper_t * helper )
{
   
  // --------------------------------------------------------------------------------
  //
  // Summarize mean/SD for per-signal Hjorth parameters
  //
  // --------------------------------------------------------------------------------
   
   mean_h1 = h1.colwise().mean();
   mean_h2 = h2.colwise().mean();
   mean_h3 = h3.colwise().mean();

   sd_h1 = ((h1.array().rowwise() - mean_h1 ).square().colwise().sum()/(h1.rows()-1)).sqrt();
   sd_h2 = ((h2.array().rowwise() - mean_h2 ).square().colwise().sum()/(h2.rows()-1)).sqrt();
   sd_h3 = ((h3.array().rowwise() - mean_h3 ).square().colwise().sum()/(h3.rows()-1)).sqrt();


   // --------------------------------------------------------------------------------
   //
   // for trainers, returns number of observed stages w/ at least suds_t::required_epoch_n 
   // -- i.e. should be suds_t::n_stages
   //
   // --------------------------------------------------------------------------------

   int nr = 0;

   std::map<std::string,int>::const_iterator cc = counts.begin();
   while ( cc != counts.end() )
     {       
       if ( cc->first != "?" && cc->second >= suds_t::required_epoch_n ) ++nr;      
       ++cc;
     }

   return trainer ? nr : nve ;
   
}


void suds_indiv_t::dump_stage_associations( const std::string & filename )
{
  //logger << "  dumping stage <-> feature/component associations to " << filename << "\n";

  std::ofstream O1( Helper::expand( filename ).c_str() , std::ios::out );
  
  // 5-level ANOVA
  //  each versus all other 2-level tests
  //  ?? stage N versus stage Y tests (25/2 pairs)
  
  // grousp  ; 'y'
  //  xxxxxx
  const int n = y.size();
  std::vector<double> is_n1( n );
  std::vector<double> is_n2( n );
  std::vector<double> is_n3( n );
  std::vector<double> is_r( n );
  std::vector<double> is_w( n );

  for (int i=0; i<n; i++)
    {
      is_n1[i] = y[i] == "N1" ;
      is_n2[i] = y[i] == "N2" ;
      is_n3[i] = y[i] == "N3" ;
      is_r[i] = y[i] == "R" ;
      is_w[i] = y[i] == "W" ;
    }
  
  //
  // for each feature
  //
  
  const std::vector<std::string> vars = suds_t::model.labels();

  const int nf = vars.size();

  if ( nf != X.cols() )
    Helper::halt( "internal error in suds_indiv_t::dump_stage_associations()" );

  //
  // header
  //
  
  O1 << "VAR\tU\tP\tF\tPS\tFS"
     << "\tN1\tN2\tN3\tR\tW";

  O1 << "\tN1_N2" 
     << "\tN1_N3"
     << "\tN1_R"
     << "\tN1_W"
     << "\tN2_N3"
     << "\tN2_R"
     << "\tN2_W"
     << "\tN3_R"
     << "\tN3_W"
     << "\tR_W";
    
  // correlation with var(row) and U components
  
  for (int u=0; u<U.cols(); u++)
    O1 << "\tU" << u+1 ;
  O1 << "\n";
  
  //
  // features
  //

  const int nn = nf + U.cols();
    
  for (int c=0; c<nn; c++)
    {

      const int is_U = c >= nf;

      // pick row
      std::vector<double> xx = eigen_ops::copy_vector( is_U ? U.col( c - nf ) : X.col(c) );
      const int n = xx.size();
      
      double f = 0 , b = 0 , w = 0;

      double pv = Statistics::anova( y  , xx , &f , &b , &w );

      O1 << ( is_U ? "U" + Helper::int2str( c - nf + 1 ) : vars[c] )
	 << "\t" << is_U ;
      
      if ( pv < 0 ) O1  << "\tNA\tNA";
      else O1 << "\t" << pv
	      << "\t" << f;

      //
      // sleep-stage (4-way test)
      //

      std::vector<double> xs;
      std::vector<std::string> ss;
      for (int k=0; k<n; k++)
	{
	  if ( y[k] != "W" ) 
	    {
	      xs.push_back( xx[k] );
	      ss.push_back( y[k] );
	    }
	}

      f = 0 , b = 0 , w = 0;
      
      pv = Statistics::anova( ss  , xs , &f , &b , &w );
      
      if ( pv < 0 ) O1  << "\tNA\tNA";
      else O1 << "\t" << pv
              << "\t" << f;
      
      
      //
      // stage-specific
      //

      O1 << "\t" << ( Statistics::correlation( xx , is_n1 ) )
	 << "\t" << ( Statistics::correlation( xx , is_n2 ) )
	 << "\t" << ( Statistics::correlation( xx , is_n3 ) )
	 << "\t" << ( Statistics::correlation( xx , is_r ) ) 
	 << "\t" << ( Statistics::correlation( xx , is_w ) );

      std::vector<std::vector<double> *> p( 5 );
      p[0] = &is_n1;
      p[1] = &is_n2;
      p[2] = &is_n3;
      p[3] = &is_r;
      p[4] = &is_w;
	
      // pairwise
      for (int i=0; i<4; i++)
	for (int j=i+1;j<5;j++)
	  {
	    std::vector<double> xf, xs;
	    int g1 = 0 , g2 = 0;
	    for (int k=0; k<n; k++)
	      {
		if ( (*p[i])[k] > 0.5 || (*p[j])[k] > 0.5 )
		  {
		    xf.push_back( xx[k] );
		    xs.push_back( (*p[i])[k] );
		    g1 += (*p[i])[k];
		    g2 += (*p[j])[k];		    
		  }
	      }

	    if ( g1 > 10 && g2 > 10 )
	      O1 << "\t" << ( Statistics::correlation( xf , xs ) ) ;
	    else
	      O1 << "\t.";

	  }

      // Corr w/ U
      for (int u=0; u<U.cols(); u++)
	O1  << "\t" <<  ( Statistics::correlation( xx , eigen_ops::copy_vector( U.col(u) ) ) );
      
      // all done
      O1 << "\n";
    }  

  O1.close();
}

