
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
  // Signals (and optionally, resampling)
  //
  
  const int ns = suds_t::model.chs.size();

  // check that all model channels are also present
  // in the EDF
  std::vector<std::string> slabs;
  std::vector<int> slots;

  signal_list_t signals;
  int si = 0;

  std::map<std::string,suds_channel_t>::const_iterator ss =  suds_t::model.chs.begin(); 
  while ( ss != suds_t::model.chs.end() )
    {
      int slot = edf.header.signal( ss->first );
      if ( slot == -1 ) Helper::halt( "could not find " + ss->first );

      if ( edf.header.is_annotation_channel( slot ) )
	Helper::halt( "cannot specificy annotation channel: " + ss->first );
      
      // need to resample?
      if ( edf.header.sampling_freq( slot ) != ss->second.sr )
        dsptools::resample_channel( edf, slot , ss->second.sr );

      // build signal_list_t
      signals.add( si , ss->first );

      ++si;
      ++ss;
    }


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

  const int ne = edf.timeline.first_epoch();
  
  //
  // PSD (Welch) parameters 
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


  logger << "  applying Welch with " << fft_segment_size << "s segments ("
	 << fft_segment_overlap << "s overlap), using "
	 << ( suds_t::use_seg_median ? "median" : "mean" )
	 << " over segments\n";  
  
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
      has_prior_staging = false;
      // nb. this is set back to 'true'  after the following 
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

  //
  // number of good (retained) epochs
  //
  
  int nge = 0 ;

  if ( has_prior_staging )
    {
      
      obs_stage.resize( ne , SUDS_UNKNOWN );
      
      for (int ss=0; ss < ne ; ss++ )
	{
	  if ( edf.timeline.hypnogram.stages[ ss ] == UNSCORED
	       || edf.timeline.hypnogram.stages[ ss ] == LIGHTS_ON
	       || edf.timeline.hypnogram.stages[ ss ] == MOVEMENT
	       || edf.timeline.hypnogram.stages[ ss ] == UNKNOWN )
	    obs_stage[ss] = SUDS_UNKNOWN;
	  
	  else if ( edf.timeline.hypnogram.stages[ ss ] == WAKE )
	    obs_stage[ss] = SUDS_WAKE;

	  else if ( edf.timeline.hypnogram.stages[ ss ] == NREM1 )
	    obs_stage[ss] = suds_t::n_stages == 3 ? SUDS_NR : SUDS_N1;

	  else if ( edf.timeline.hypnogram.stages[ ss ] == NREM2 )
	    obs_stage[ss] = suds_t::n_stages == 3 ? SUDS_NR : SUDS_N2;

	  else if ( edf.timeline.hypnogram.stages[ ss ] == NREM3
		    || edf.timeline.hypnogram.stages[ ss ] == NREM4 )
	    obs_stage[ss] = suds_t::n_stages == 3 ? SUDS_NR : SUDS_N3;
	  
	  else if ( edf.timeline.hypnogram.stages[ ss ] == REM )
	    obs_stage[ss] = SUDS_REM;

	  
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
  // Size feature matrix X
  //

  nf = suds_t::nf;

  X.resize( nge , nf );

  logger << "  expecting " << nf << " features (for " << nge << " epochs) and " << ns << " channels\n";

  
  //
  // for QC, estimate Hjorth parameters over
  // epochs, for each signal
  //
  
  h1 = Eigen::MatrixXd::Zero( nge , ns );
  h2 = Eigen::MatrixXd::Zero( nge , ns );
  h3 = Eigen::MatrixXd::Zero( nge , ns );
  
  //
  // Track bad epochs 
  //
  
  std::set<int> bad_epochs; // of 0..nge-1 encoding 
  

  //
  // iterate over (retained) epochs
  //
  
  int en = 0 , en_good = 0;
  
  edf.timeline.first_epoch();
  
  epochs.clear();
  
  while ( 1 ) 
    {
      
      //
      // select epoch
      //


      int epoch = edf.timeline.next_epoch();      	  
      
      if ( epoch == -1 ) break;
      
      if ( en == ne ) Helper::halt( "internal error: over-counted epochs" );
      
      
      //
      // retained? if not, skip
      //
      
      if ( ! retained[ en ] ) 
   	{
   	  ++en;
   	  continue;
   	}
      
      //
      // Process this epoch, signal-by-signal, then feature-spec by feature-spec.
      //
      
      interval_t interval = edf.timeline.epoch( epoch );
      
      //
      // is this a bad epoch?
      //
      
      bool bad_epoch = false;

      //
      // Iterate over signals
      //

      for (int s = 0 ; s < ns; s++ )
	{

	  //
	  // Skip if a bad epoch has been flagged already for a previous channel
	  //

	  if ( bad_epoch ) continue;

	  //
	  // Get data
	  //

	  const std::string siglab = signals.label(s) ;
	  
	  slice_t slice( edf , signals(s) , interval );
	  
	  const int sr = edf.header.sampling_freq( signals(s) ); 

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
	  
	  const bool do_mean = suds_t::model.has( suds_feature_t::SUDS_MEAN , siglab ) ;
	  
	  const bool do_spectral =
	    suds_t::model.has( suds_feature_t::SUDS_LOGPSD , siglab ) ||
	    suds_t::model.has( suds_feature_t::SUDS_RELPSD , siglab ) ||
	    suds_t::model.has( suds_feature_t::SUDS_SLOPE , siglab ) ||
	    suds_t::model.has( suds_feature_t::SUDS_CVPSD , siglab );
	  
	  const bool do_skew = suds_t::model.has( suds_feature_t::SUDS_SKEW , siglab );
	  
	  const bool do_kurt = suds_t::model.has( suds_feature_t::SUDS_KURTOSIS , siglab );
	  
	  const bool do_hjorth = suds_t::model.has( suds_feature_t::SUDS_HJORTH , siglab );
	  
	  const bool do_pe = suds_t::model.has( suds_feature_t::SUDS_PE , siglab );
	  
	  const bool do_pfd = suds_t::model.has( suds_feature_t::SUDS_FD , siglab );
	  	  	  

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
	      const bool get_segment_sd = suds_t::model.has( suds_feature_t::SUDS_CVPSD , siglab ) ;
	      
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
	      
	      if ( suds_t::model.has( suds_feature_t::SUDS_LOGPSD , siglab ) && ! bad_epoch )
		{
		  
		  std::vector<int> cols = suds_t::model.cols( suds_feature_t::SUDS_LOGPSD , siglab ) ;
		  const int ncols = cols.size();
		  
		  // this *should* map exactly onto the number of bins between the lwr and upr bounds
		  
		  suds_spec_t spec = suds_t::model.fcmap[ suds_feature_t::SUDS_LOGPSD ][ siglab ];
		  
		  // these have been checked and will be present/valid 
		  const double lwr = spec.arg[ "lwr" ];
		  const double upr = spec.arg[ "upr" ];
		  
		  int b = 0;
		  
		  for ( int i = 0 ; i < bin.bfa.size() ; i++ )
		    {
		      if (  bin.bfa[i] >= lwr && bin.bfa[i] <= upr )
			{
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
	      
	      if ( suds_t::model.has( suds_feature_t::SUDS_RELPSD , siglab ) && ! bad_epoch )
		{
		  std::vector<int> cols = suds_t::model.cols( suds_feature_t::SUDS_RELPSD , siglab ) ;
		  const int ncols = cols.size();
		  
		  suds_spec_t spec = suds_t::model.fcmap[ suds_feature_t::SUDS_RELPSD ][ siglab ];
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
	      
	      if ( suds_t::model.has( suds_feature_t::SUDS_CVPSD , siglab ) && ! bad_epoch )
		{
		  
		  std::vector<int> cols = suds_t::model.cols( suds_feature_t::SUDS_CVPSD , siglab ) ;
		  const int ncols = cols.size();
		  
		  suds_spec_t spec = suds_t::model.fcmap[ suds_feature_t::SUDS_CVPSD ][ siglab ];
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
	      
	      if ( suds_t::model.has( suds_feature_t::SUDS_SLOPE , siglab ) && ! bad_epoch )
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
		  std::vector<int> cols = suds_t::model.cols( suds_feature_t::SUDS_SLOPE , siglab ) ;
		  
		  // save slope
		  X( en_good , cols[0] ) = bslope;
		   
		}
	      
	    }
	  

	   //
	   // Time domain features
	   //

	   if ( do_mean && ! bad_epoch )
	     {
	       std::vector<int> cols = suds_t::model.cols( suds_feature_t::SUDS_MEAN , siglab ) ;
	       X( en_good , cols[0] ) = mean; // calculated above when mean-centering
	     }
	   
	   if ( do_skew && ! bad_epoch )
	     {
	       std::vector<int> cols = suds_t::model.cols( suds_feature_t::SUDS_SKEW , siglab ) ;
	       X( en_good , cols[0] ) = MiscMath::skewness( *d , 0 , MiscMath::sdev( *d , 0 ) );
	     }

	   if ( do_kurt && ! bad_epoch )
	     {
	       std::vector<int> cols = suds_t::model.cols( suds_feature_t::SUDS_KURTOSIS , siglab ) ;
	       X( en_good , cols[0] ) = MiscMath::kurtosis0( *d ); // assumes mean-centered
	     }
	   
	   // fractal dimension
	   if ( do_pfd && ! bad_epoch )
	     {
	       std::vector<int> cols = suds_t::model.cols( suds_feature_t::SUDS_FD , siglab ) ;
	       X( en_good , cols[0] ) = MiscMath::petrosian_FD( *d );
	     }
	   
	   // permutation entropy
	   if ( do_pe && ! bad_epoch )
	     {
	       std::vector<int> cols = suds_t::model.cols( suds_feature_t::SUDS_PE , siglab ) ;
	       
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
		   std::vector<int> cols = suds_t::model.cols( suds_feature_t::SUDS_HJORTH , siglab ) ;
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


   // --------------------------------------------------------------------------------
   //
   // Get PSC initially (we look for outliers and then remove epochs, and redo the SVD)
   //
   // --------------------------------------------------------------------------------
   
   
   Eigen::BDCSVD<Eigen::MatrixXd> svd( X , Eigen::ComputeThinU | Eigen::ComputeThinV );
   U = svd.matrixU();
   V = svd.matrixV();
   W = svd.singularValues();

   
   // --------------------------------------------------------------------------------
   //
   // Outliers in PSC space? 
   //
   // --------------------------------------------------------------------------------
   
   std::vector<bool> valid( nge , true );
   
   // track reasons for exclusion
   std::set<int> nout_flat;
   std::set<int> nout_hjorth;
   std::set<int> nout_stat;
   std::set<int> nout_tot;

   
   //
   // Exclusions based on H==0 parameters
   //
   
   for ( int s=0;s<ns;s++)
     {
       for (int i=0;i<nge;i++) 
	 {
	   if      ( h1( i , s ) < 1e-8 ) { valid[i] = false; nout_flat.insert(i); }
	   else if ( h2( i , s ) < 1e-8 ) { valid[i] = false; nout_flat.insert(i) ; }
	   else if ( h3( i , s ) < 1e-8 ) { valid[i] = false; nout_flat.insert(i) ; }
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

       for ( int s=0;s<ns;s++)
	 {
   	  for (int i=0;i<nge;i++)
   	    {
   	      if ( h1( i, s ) <= suds_t::hjorth1_lwr95[s] || h1(i,s) >= suds_t::hjorth1_upr95[s] ||
		   h2( i, s ) <= suds_t::hjorth2_lwr95[s] || h2(i,s) >= suds_t::hjorth2_upr95[s] ||
		   h3( i, s ) <= suds_t::hjorth3_lwr95[s] || h3(i,s) >= suds_t::hjorth3_upr95[s] )
		{
		  valid[i] = false;
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
		   if ( x[c] < lwr || x[c] > upr ) { valid[i] = false; nout_stat.insert(i); } 
		   ++c;
		 }
	     }
	 }
     }
   
   
   
   // --------------------------------------------------------------------------------
   //
   // Impose a max number of epochs per stage? 
   //
   // --------------------------------------------------------------------------------
   
   if ( has_prior_staging && suds_t::max_epoch_n != -1 )
     {
   
       std::map<suds_stage_t,std::vector<int> > counts;
   
       int cc = 0;
       for (int i=0;i<ne;i++)
	 {
	   if ( retained[i] )
	     {
	       // track counts in valid index space
	       if ( valid[cc] )
		 counts[ obs_stage[ i ] ].push_back( cc );
	       ++cc;
	     }
	 }
       
       std::map<suds_stage_t,std::vector<int> >::const_iterator qq = counts.begin();
       while ( qq != counts.end() )
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
		   if ( valid[ qq->second[ pick ] ] )
		     {
		       valid[ qq->second[ pick ] ] = false;
		       --rem;
		     }
		 }	      
	     }
	   ++qq;
	 }
       
     }
   
   
   // --------------------------------------------------------------------------------
   //
   // Summarize dropped epochs and remove 
   //
   // --------------------------------------------------------------------------------
   
   int included = 0;

   for (int i=0;i<nge;i++)
     if ( valid[i] ) ++included;
   
   logger << "  of " << ne << " total epochs, valid staging for " << nge
          << ", and of those " << included << " passed outlier removal\n";
   
   
   std::set<int>::const_iterator oo = nout_flat.begin();
   while ( oo != nout_flat.end() ) { nout_tot.insert( *oo ); ++oo; } 
   oo = nout_hjorth.begin();
   while ( oo != nout_hjorth.end() ) { nout_tot.insert( *oo ); ++oo; } 
   oo = nout_stat.begin();
   while ( oo != nout_stat.end() ) { nout_tot.insert( *oo ); ++oo; } 
   
   logger << "  outliers counts (flat, Hjorth, components, total = "
   	 << nout_flat.size() << ", "
	  << nout_hjorth.size() << ", "
	  << nout_stat.size() << ", "
	  << nout_tot.size() << ")\n";
   

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
       if ( valid[i] )
	 {	   
	   for (int j=0;j<nf;j++)
	      X(r,j) = X2(i,j);
	   
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
   // splice out bad epochs for Hjorth parameters
   //

   Eigen::MatrixXd hh1 = h1;
   Eigen::MatrixXd hh2 = h2;
   Eigen::MatrixXd hh3 = h3;
   h1.resize( nve , ns );
   h2.resize( nve , ns ); 
   h3.resize( nve , ns );

   for (int s=0;s<ns;s++)
     {
       int r = 0;      
       for (int i=0; i < valid.size(); i++ )
	 if ( valid[i] )
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
      
   for (int c=0; c<suds_t::nf; c++)
     X.col(c) *= suds_t::model.W[c];

   
   
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
	   eigen_ops::scale( c , true , true );

	   bool okay = true;

	   writer.level( "PSC_" + Helper::int2str( j+1 ) , "VAR");
	   
	   if ( do_anova )
	     {
	       double pv = Statistics::anova( ss_str  , eigen_ops::copy_vector( c ) );
	       writer.value( "PV", pv );	       
	       if ( pv < 0 || pv > suds_t::required_comp_p ) okay = false;
	     }

	   // may have signif stage/group differences, but check that no one stage has a big variance difference also 
	   if ( do_bw )
	     {
	       // nb. c is standardized
	       double wb = eigen_ops::between_within_group_variance( ss_str , c );
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
   

   // --------------------------------------------------------------------------------
   //
   // Re-Standardize U 
   //
   // --------------------------------------------------------------------------------
   
   if ( suds_t::standardize_U )
     {
       
       if ( suds_t::robust_standardization )
	 {
	   logger << "  robust re-standardizing U";
	   if ( suds_t::winsor2 > 0 ) logger << ", winsorizing at " << suds_t::winsor2;
	   logger << "\n";
	   eigen_ops::robust_scale( U , true, true, suds_t::winsor2 );
	 }
       else
   	{
   	  logger << "  re-standardizing U\n";
   	  eigen_ops::scale( U , true , true );
   	}  
     }

   
   // --------------------------------------------------------------------------------
   //
   // make class labels ( trainer only )
   //
   // --------------------------------------------------------------------------------
   
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
   

   // --------------------------------------------------------------------------------
   //
   // Self-classification (i.e. SOAP) to remove epochs that aren't well self-classified
   //  - possibly reject a trainer, if their SOAP kappa is poor
   //
   // --------------------------------------------------------------------------------
   
   if ( trainer && suds_t::self_classification )
     {

       std::vector<bool> okay ;

       // this returns the number of 'good' epochs 

       int nve2 = self_classify( &okay );
       
       if ( nve2 == 0 )
	 {
	   logger << "  trainer not valid based on self-classification thresholds\n";
	   return 0;
	 }
       
       //
       // Subset epochs
       //
       
       //   U X epochs  y  h1 h2  h3
    
       Eigen::MatrixXd UU = U;      
       U.resize( nve2 , nc );      
       
       Eigen::MatrixXd XX = X;      
       X.resize( nve2 , X.cols() );      
       
       std::vector<int> epochs2 = epochs;
       epochs.clear();
       
       std::vector<suds_stage_t> obs_stage_valid2 = obs_stage_valid;
       if ( has_prior_staging )
	 obs_stage_valid.clear();
       
       Eigen::MatrixXd hh1 = h1;
       h1.resize( nve , ns );

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

      	      // X original features
   	      for (int j=0;j<nf;j++)
   		X(r,j) = XX(i,j);
	      
   	      // Epoch tracking
   	      epochs.push_back( epochs2[i] );
	      
   	      // nb. take from already-pruned set, obs_stage_valid[] ) 
   	      if ( has_prior_staging )
		obs_stage_valid.push_back( obs_stage_valid2[i] );

	      
	      // Hjorth (per signal)
	      for (int s=0;s<ns;s++)
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
       // Redo labels
       //
       
       std::vector<std::string> yy = y;
       y.clear();
       for ( int i = 0 ; i < nve ; i++ )
	 if ( okay[i] ) y.push_back( yy[i] );
       
       // just in case...?
       if ( y.size() != obs_stage_valid.size() ) 
	 Helper::halt( "internal error in proc()" );
       
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

