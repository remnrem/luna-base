
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

#include "pops/indiv.h"

#include "helper/helper.h"
#include "helper/logger.h"
#include "stats/eigen_ops.h"
#include "db/db.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "dsp/resample.h"

#include "fftw/fftwrap.h"
#include "pdc/pdc.h"

extern logger_t logger;
extern writer_t writer;

// - what people know
//  - different modalities / ways to study sleep
//  - applications
//  - brain,
//  - heart
//  - lung
  
pops_indiv_t::pops_indiv_t( edf_t & edf ,
			    param_t & param )
{

  const bool training_mode = param.has( "train" );
  
  const bool dump_features = param.has( "dump" );

  trainer = training_mode;
  
  // training (1) : make level-1 stats, stages, save (binary features, BFTR)
  
  // training (2) : load prior data (BFTR files) for each individual in training, (+ in validation)
  //    -> then compile all indivs to make a super-set (still track indiv. level columns)
  //    -> make derived (level 2) metrics (some of these are done across all, e.g. SVD; some within (e.g. NORM)
  //    -> then fit LGBM model, and save
  
  // prediction (1) [ all indiv_t level ] 
  //  - bring in new EDFs
  //  - level-1 stats
  //  - derive level-2 stats for that one individual
  //  - load model
  //  - make prediction
  
  // get any staging
  staging( edf , param );


  
  //
  // training mode: derive level-1 stats and then quit
  //
  
  if ( training_mode )
    {
      level1( edf );
      save1( edf.id , param.requires( "data" ) );      
    }

  
  //
  // Predict: level 1 & 2 states, then fit 
  //

  if ( ! training_mode )
    {

      level1( edf );

      level2();
      
      logger << "  final feature matrix: " << X1.rows() << " rows (epochs) and " << X1.cols() << " columns (features)\n";
      
      if ( dump_features ) 
	{
	  std::string dfile = Helper::expand( param.value( "dump" ) );
	  logger << "  dumping feature matrix to " << dfile << "\n";
	  std::ofstream O1( dfile.c_str() , std::ios::out );
	  O1 << "SS";
	  std::vector<std::string> labels = pops_t::specs.select_labels();
	  for (int i=0; i<labels.size(); i++) O1 << "\t" << labels[i];
	  O1 << "\n";
	  for (int i=0; i<X1.rows(); i++)
	    {
	      O1 << pops_t::label( (pops_stage_t)S[i] );
	      for (int j=0; j<X1.cols(); j++)
		O1 << "\t" << X1(i,j);
	      O1 << "\n";
	    }
	  O1.close();
	}


      if ( ! pops_t::lgbm_model_loaded )
	{
	  pops_t::lgbm.load_model( param.requires( "model" ) );
	  if ( param.has( "config" ) ) 
	    pops_t::lgbm.load_config( param.value( "config" ) );	  
	  pops_t::lgbm_model_loaded = true;
	}
      
      predict();
      
      SHAP();

      summarize();
    }
  
}

void pops_indiv_t::staging( edf_t & edf , param_t & param )
{

  // calculate ne and staging, if present  
  ne = edf.timeline.first_epoch();
  
  // get staging
  edf.timeline.annotations.make_sleep_stage();
  
  bool valid_training = edf.timeline.hypnogram.construct( &(edf.timeline) , param , false );
  
  // trainer
  if ( trainer && ! valid_training )
    Helper::halt( "no valid staging for trainer " + edf.id );
  
  // check epochs line up, if staging present
  if ( valid_training && ne != edf.timeline.hypnogram.stages.size() )    
    Helper::halt( "problem extracting stage information for trainer" );

  // store staging information here
  S.resize( ne , POPS_UNKNOWN );
  E.resize( ne );
  
  // for targets w/ no existing staging, all done
  if ( ! valid_training ) return;

  // convert 
      
  for (int ss=0; ss < ne; ss++ )
    {

      // track 0-based epoch numbers
      E[ss] = ss;
      
      if ( edf.timeline.hypnogram.stages[ ss ] == UNSCORED
	   || edf.timeline.hypnogram.stages[ ss ] == LIGHTS_ON
	   || edf.timeline.hypnogram.stages[ ss ] == MOVEMENT
	   || edf.timeline.hypnogram.stages[ ss ] == UNKNOWN )
	S[ss] = POPS_UNKNOWN;
      
      else if ( edf.timeline.hypnogram.stages[ ss ] == WAKE )
	S[ss] = POPS_WAKE;
      
      else if ( edf.timeline.hypnogram.stages[ ss ] == NREM1 )
	S[ss] = pops_opt_t::n_stages == 3 ? POPS_N1 : POPS_N1;
      
      else if ( edf.timeline.hypnogram.stages[ ss ] == NREM2 )
	S[ss] = pops_opt_t::n_stages == 3 ? POPS_N1 : POPS_N2;
      
      else if ( edf.timeline.hypnogram.stages[ ss ] == NREM3
		|| edf.timeline.hypnogram.stages[ ss ] == NREM4 )
	S[ss] = pops_opt_t::n_stages == 3 ? POPS_N1 : POPS_N3;
      
      else if ( edf.timeline.hypnogram.stages[ ss ] == REM )
	S[ss] = POPS_REM;

          
    } // next epoch 

  //
  // copy original staging  (i.e. as S is set to POPS_UNKNOWN for bad signals, trimming etc)
  //
  
  Sorig = S;
  
  //
  // trim leading/trailing wake epochs? 
  //
  
  if ( pops_opt_t::trim_wake_epochs >= 0 ) 
    {
      int first_sleep = -1;
      for (int ss=0; ss < ne; ss++)
	{
	  if ( S[ss] == POPS_N1 ||
	       S[ss] == POPS_N2 ||
	       S[ss] == POPS_N3 ||
	       S[ss] == POPS_REM ) 
	    {
	      first_sleep = ss;
	      break;
	    }
	}
      
      int last_sleep = ne - 1;
      for (int ss = ne - 1; ss >=0 ; ss--)
	{
	  if ( S[ss] == POPS_N1 ||
	       S[ss] == POPS_N2 ||
	       S[ss] == POPS_N3 ||
	       S[ss] == POPS_REM ) 
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
	  
	  first_sleep -= pops_opt_t::trim_wake_epochs + 1 ; 	      
	  int t = 0;
	  // note, inclusive counting up to X
	  for (int ss=0; ss<= first_sleep; ss++)
	    {
	      S[ss] = POPS_UNKNOWN; 
	      ++t;
	    }
	  if ( t ) logger << "  trimmed " << t << " leading wake epochs\n";
	  
	}
	  
      // trim end
      if ( last_sleep < ne - 1 ) 
	{
	  // * *               
	  // 4 5 6 7 8 9
	  //         X X
	  
	  last_sleep += pops_opt_t::trim_wake_epochs + 1 ;
	  int t=0;
	  for (int ss= ne - 1 ; ss >= last_sleep; ss--)
	    {
	      S[ss] = POPS_UNKNOWN;
	      ++t;
	    }
	  if ( t ) logger << "  trimmed " << t << " trailing wake epochs\n";
	}
	  
    } // end of wake trimming option

  
}


void pops_indiv_t::level1( edf_t & edf )
{

  //
  // score level-1 factors --> X1
  //

  X1.resize( ne , pops_t::specs.n1 );

  logger << "  expecting " << pops_t::specs.n1
	 << " level-1 features (for " << ne
	 << " epochs) and " << pops_t::specs.ns << " signals\n";
    
  //
  // PSD (Welch) parameters 
  //

  double fft_segment_size = pops_opt_t::fft_seg_sec ;  // 4
  double fft_segment_overlap = pops_opt_t::fft_inc_sec ; // 2
  
  if ( edf.timeline.epoch_length() <= ( fft_segment_size + fft_segment_overlap ) )
    {
      fft_segment_overlap = 0;
      fft_segment_size = edf.timeline.epoch_length();
    }

  window_function_t window_function = WINDOW_TUKEY50;	   
    
  logger << "  applying Welch with " << fft_segment_size << "s segments ("
	 << fft_segment_overlap << "s overlap), using "
	 << ( pops_opt_t::welch_median ? "median" : "mean" )
	 << " over segments\n";  
  

  //
  // check signals present in EDF
  //

  std::vector<std::string> slabs;
  std::vector<int> slots;
  const bool silent_signal_search = true;

  signal_list_t signals;
  
  std::map<std::string,pops_channel_t>::const_iterator ss =  pops_t::specs.chs.begin(); 
  while ( ss != pops_t::specs.chs.end() )
    {
      
      // primary?
      int slot = edf.header.signal( ss->first , silent_signal_search );

      // match on an alias?
      if ( slot == -1 ) 
	{

	  const std::set<std::string> & aliases = ss->second.aliases;

	  std::set<std::string>::const_iterator aa = aliases.begin();
	  while ( aa != aliases.end() )
	    {
	      slot = edf.header.signal( *aa , silent_signal_search );
	      if ( slot != -1 ) break;
	      ++aa;
	    }
	}

      // still no match?
      if ( slot == -1 ) 
	Helper::halt( "could not find " + ss->first + " (or any specified aliases)" );
      
      if ( edf.header.is_annotation_channel( slot ) )
	Helper::halt( "cannot specificy annotation channel: " + ss->first );
      
      // need to resample?
      if ( edf.header.sampling_freq( slot ) != ss->second.sr )
        dsptools::resample_channel( edf, slot , ss->second.sr );
      
      // need to rescale?
      if ( Helper::toupper( edf.header.phys_dimension[ slot ] ) != Helper::toupper( ss->second.unit ) )
	{
	  logger << "  rescaling " << ss->first 
		 << " from " << edf.header.phys_dimension[ slot ] 
		 << " to " << ss->second.unit << "\n";
	  edf.rescale( slot , ss->second.unit );
	}

      // build signal_list_t
      signals.add( slot , ss->first );
      
      ++ss;
    }


  //
  // iterate over epochs
  //
  
  int en = 0 ;
  
  edf.timeline.first_epoch();
  
  while ( 1 ) 
    {
      
      int epoch = edf.timeline.next_epoch();      	  
      
      if ( epoch == -1 ) break;
      
      if ( en == ne ) Helper::halt( "internal error: over-counted epochs" );
      
      //
      // skip?
      //
      
      if ( S[ en ] == POPS_UNKNOWN )
   	{
   	  ++en;
   	  continue;
   	}
      
      //
      // Process epoch: signal-by-signal, then feature-spec by feature-spec.
      //
      
      interval_t interval = edf.timeline.epoch( epoch );
      
      bool bad_epoch = false;

      //
      // Iterate over signals
      //

      const int ns = pops_t::specs.ns;
      
      for (int s = 0 ; s < ns; s++ )
	{

	  //
	  // Skip if flagged for a prior channel?
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
	  
	  double mean = MiscMath::centre( d );
	  
	  
	  //
	  // extract these channel-specific features 
	  //     
	  
	  const bool do_mean = pops_t::specs.has( pops_feature_t::POPS_MEAN , siglab ) ;
	  
	  const bool do_spectral =
	    pops_t::specs.has( pops_feature_t::POPS_LOGPSD , siglab ) ||
	    pops_t::specs.has( pops_feature_t::POPS_RELPSD , siglab ) ||	    
	    pops_t::specs.has( pops_feature_t::POPS_CVPSD , siglab ) || 
	    pops_t::specs.has( pops_feature_t::POPS_BANDS , siglab ) ||
	    pops_t::specs.has( pops_feature_t::POPS_RBANDS , siglab ) ||	    
	    pops_t::specs.has( pops_feature_t::POPS_VBANDS , siglab ) || 	  
	    pops_t::specs.has( pops_feature_t::POPS_SLOPE , siglab );

	  const bool do_skew = pops_t::specs.has( pops_feature_t::POPS_SKEW , siglab );
	  
	  const bool do_kurt = pops_t::specs.has( pops_feature_t::POPS_KURTOSIS , siglab );
	  
	  const bool do_hjorth = pops_t::specs.has( pops_feature_t::POPS_HJORTH , siglab );
	  
	  const bool do_pe = pops_t::specs.has( pops_feature_t::POPS_PE , siglab );
	  
	  const bool do_pfd = pops_t::specs.has( pops_feature_t::POPS_FD , siglab );
	  	  	  

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
	      const bool get_segment_sd = pops_t::specs.has( pops_feature_t::POPS_CVPSD , siglab ) ;
	      
	      PWELCH pwelch( *d , 
			     sr , 
			     segment_sec , 
			     noverlap_segments , 
			     window_function ,
			     pops_opt_t::welch_median , 
			     get_segment_sd );
	      
	      // using bin_t, 1 means no binning
	      bin_t bin( pops_opt_t::lwr , pops_opt_t::upr , 1 ); 
	      bin.bin( pwelch.freq , pwelch.psd );	      
	      
	      //
	      // check for zero power values in the 0.5 to 45 Hz range, and flag if so
	      //  -- we will not include this epoch
	      //
	      
	      for ( int i = 0 ; i < bin.bfa.size() ; i++ )
		{
		  if ( bin.bfb[i] > pops_opt_t::upr ) break;
		  if ( bin.bspec[i] <= 0 && bin.bfa[i] >= pops_opt_t::lwr ) 
		    {
		      bad_epoch  = true;		       
		      bin.bspec[i] = 1e-4 ; // set to -40dB as a fudge		   
		    }
		}

	      
	      //
	      // track that this is bad / to be removed below?
	      //

	      if ( bad_epoch ) S[ en ] = POPS_UNKNOWN;
	      
	      //
	      // log-PSD?
	      //
	      
	      if ( pops_t::specs.has( pops_feature_t::POPS_LOGPSD , siglab ) && ! bad_epoch )
		{
		  
		  std::vector<int> cols = pops_t::specs.cols( pops_feature_t::POPS_LOGPSD , siglab ) ;
		  const int ncols = cols.size();
		  
		  // this *should* map exactly onto the number of bins between the lwr and upr bounds
		  pops_spec_t spec = pops_t::specs.fcmap[ pops_feature_t::POPS_LOGPSD ][ siglab ];

		  
		  // these have been checked and will be present/valid 
		  const double lwr = spec.narg( "lwr" );
		  const double upr = spec.narg( "upr" );

		  int b = 0;
		  
		  for ( int i = 0 ; i < bin.bfa.size() ; i++ )
		    {
		      if (  bin.bfa[i] >= lwr && bin.bfa[i] <= upr )
			{
			  if ( b == ncols ) Helper::halt( "internal error... bad sizes for SPEC" );

			  // save log-scaled power
			  X1( en , cols[b] ) = 10*log10( bin.bspec[i] ) ; 
			  
			  // next feature column
			  ++b;			   
			}
		    }
		}
	      
	      
	      //
	      // rel-PSD?
	      //
	      
	      if ( pops_t::specs.has( pops_feature_t::POPS_RELPSD , siglab ) && ! bad_epoch )
		{
		  std::vector<int> cols = pops_t::specs.cols( pops_feature_t::POPS_RELPSD , siglab ) ;
		  const int ncols = cols.size();
		  
		  pops_spec_t spec = pops_t::specs.fcmap[ pops_feature_t::POPS_RELPSD ][ siglab ];
		  const double lwr = spec.narg( "lwr" );
		  const double upr = spec.narg( "upr" );
		  
		  const double zlwr = spec.narg( "z-lwr" );
		  const double zupr = spec.narg( "z-upr" );
		  
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
			  X1( en , cols[b] ) = log( bin.bspec[i] / norm ) ; 
			  ++b;			   
			}
		    }
		}
	      
	      
	      //
	      // cv-PSD?
	      //
	      
	      if ( pops_t::specs.has( pops_feature_t::POPS_CVPSD , siglab ) && ! bad_epoch )
		{
		  
		  std::vector<int> cols = pops_t::specs.cols( pops_feature_t::POPS_CVPSD , siglab ) ;
		  const int ncols = cols.size();
		  
		  pops_spec_t spec = pops_t::specs.fcmap[ pops_feature_t::POPS_CVPSD ][ siglab ];
		  const double lwr = spec.narg( "lwr" );
		  const double upr = spec.narg( "upr" );
		  
		  int b = 0;
		  
		  for ( int i = 0 ; i < pwelch.freq.size() ; i++ )
		    {
		      if (  pwelch.freq[i] >= lwr && pwelch.freq[i] <= upr )
			{
			  if ( b == ncols ) Helper::halt( "internal error... bad sizes for VSPEC" );
			  
			  // save CV of PSD
			  X1( en , cols[b] ) = pwelch.psdsd[i];
			  
			  // next feature column
			  ++b;			   
			}
		    }
		  
		}
	      

	      //
	      // Band power? (abs, rel or CV)
	      //
	      
	      const bool do_bands = pops_t::specs.has( pops_feature_t::POPS_BANDS , siglab )  
		|| pops_t::specs.has( pops_feature_t::POPS_RBANDS , siglab )  
		|| pops_t::specs.has( pops_feature_t::POPS_VBANDS , siglab );

	      if ( do_bands && ! bad_epoch )
		{
		  
		  // fixed 6 bands: 
		  if (  pops_t::specs.has( pops_feature_t::POPS_BANDS , siglab ) 
			|| pops_t::specs.has( pops_feature_t::POPS_RBANDS , siglab ) )
		    {

		      double p_slow = pwelch.psdsum( SLOW );
		      double p_delta = pwelch.psdsum( DELTA );
		      double p_theta = pwelch.psdsum( THETA );
		      double p_alpha = pwelch.psdsum( ALPHA );
		      double p_sigma = pwelch.psdsum( SIGMA );
		      double p_beta = pwelch.psdsum( BETA );

		      if ( pops_t::specs.has( pops_feature_t::POPS_BANDS , siglab ) )
			{
			  std::vector<int> cols = pops_t::specs.cols( pops_feature_t::POPS_BANDS , siglab ) ;
			  const int ncols = cols.size();
			  if ( ncols != 6 ) Helper::halt( "internal error in bands" );
			  // abs power PSD bands
			  X1( en , cols[0] ) = log( p_slow );
			  X1( en , cols[1] ) = log( p_delta );
			  X1( en , cols[2] ) = log( p_theta );
			  X1( en , cols[3] ) = log( p_alpha );
			  X1( en , cols[4] ) = log( p_sigma );
			  X1( en , cols[5] ) = log( p_beta );
			}
		      
		      if ( pops_t::specs.has( pops_feature_t::POPS_RBANDS , siglab ) )
			{
			  std::vector<int> cols = pops_t::specs.cols( pops_feature_t::POPS_RBANDS , siglab ) ;
                          const int ncols = cols.size();
                          if ( ncols != 6 ) Helper::halt( "internal error in rbands" );
			  const double p_total = p_slow + p_delta + p_theta + p_alpha + p_sigma + p_beta;
			  // rel power PSD bands
                          X1( en , cols[0] ) = p_slow / p_total ;
                          X1( en , cols[1] ) = p_delta / p_total ;
                          X1( en , cols[2] ) = p_theta / p_total ;
                          X1( en , cols[3] ) = p_alpha / p_total ;
                          X1( en , cols[4] ) = p_sigma / p_total ;
                          X1( en , cols[5] ) = p_beta / p_total ;
                        }
		    }


		  // VBANDS 

		  if ( pops_t::specs.has( pops_feature_t::POPS_VBANDS , siglab ) )
		    {		      
		      std::vector<int> cols = pops_t::specs.cols( pops_feature_t::POPS_VBANDS , siglab ) ;
		      const int ncols = cols.size();
		      if ( ncols != 6 ) Helper::halt( "internal error in vbands" );
		      // save CV of PSD bands
		      X1( en , cols[0] ) = pwelch.psdsdsum( SLOW );
		      X1( en , cols[1] ) = pwelch.psdsdsum( DELTA );
		      X1( en , cols[2] ) = pwelch.psdsdsum( THETA );
		      X1( en , cols[3] ) = pwelch.psdsdsum( ALPHA );
		      X1( en , cols[4] ) = pwelch.psdsdsum( SIGMA );
		      X1( en , cols[5] ) = pwelch.psdsdsum( BETA );
		    }


		}

	      

	      //
	      // Spectral slope?
	      //
	      
	      if ( pops_t::specs.has( pops_feature_t::POPS_SLOPE , siglab ) && ! bad_epoch )
		{
		   
		  double bslope = 0, bn = 0;
		  
		  bool okay = spectral_slope_helper( pwelch.psd , 
						     pwelch.freq , 
						     pops_opt_t::slope_range ,
						     pops_opt_t::slope_th , 
						     false ,  // do not output value
						     &bslope , &bn ); 
		  if ( ! okay ) bad_epoch = true;
		  
		  // will be exactly size == 1 
		  std::vector<int> cols = pops_t::specs.cols( pops_feature_t::POPS_SLOPE , siglab ) ;
		  
		  // save slope
		  X1( en , cols[0] ) = bslope;
		   
		}
	      
	    }
	  

	   //
	   // Time domain features
	   //

	   if ( do_mean && ! bad_epoch )
	     {
	       std::vector<int> cols = pops_t::specs.cols( pops_feature_t::POPS_MEAN , siglab ) ;
	       X1( en , cols[0] ) = mean; // calculated above when mean-centering
	     }
	   
	   if ( do_skew && ! bad_epoch )
	     {
	       std::vector<int> cols = pops_t::specs.cols( pops_feature_t::POPS_SKEW , siglab ) ;
	       X1( en , cols[0] ) = MiscMath::skewness( *d , 0 , MiscMath::sdev( *d , 0 ) );
	     }

	   if ( do_kurt && ! bad_epoch )
	     {
	       std::vector<int> cols = pops_t::specs.cols( pops_feature_t::POPS_KURTOSIS , siglab ) ;
	       X1( en , cols[0] ) = MiscMath::kurtosis0( *d ); // assumes mean-centered
	     }
	   
	   // fractal dimension
	   if ( do_pfd && ! bad_epoch )
	     {
	       std::vector<int> cols = pops_t::specs.cols( pops_feature_t::POPS_FD , siglab ) ;
	       X1( en , cols[0] ) = MiscMath::petrosian_FD( *d );
	     }
	   
	   // permutation entropy
	   if ( do_pe && ! bad_epoch )
	     {
	       std::vector<int> cols = pops_t::specs.cols( pops_feature_t::POPS_PE , siglab ) ;
	       
	       int sum1 = 1;
	       pops_spec_t spec = pops_t::specs.fcmap[ pops_feature_t::POPS_PE ][ siglab ];
	       int n1 = spec.narg( "from" );
	       int n2 = spec.narg( "to" );
	       if ( cols.size() != n2 - n1 + 1 ) 
		 Helper::halt("internal error in PE cols" );

	       int k = 0;
	       for (int j=n1; j<=n2; j++)
		 {
		   std::vector<double> pd = pdc_t::calc_pd( *d , j , 1 , &sum1 );		   
		   X1( en , cols[k++] ) = pdc_t::permutation_entropy( pd );
		 }	       
	     }

	   //
	   // Hjorth parameters: these are always calculated for (trainer) QC, but 
	   // they may also be added as explicit features
	   //
	   
	   if ( do_hjorth && ! bad_epoch )
	     {
	       double activity = 0 , mobility = 0 , complexity = 0;
	       MiscMath::hjorth( d , &activity , &mobility , &complexity );

	       // use all 3 parameters (log-scaling H1)
	       std::vector<int> cols = pops_t::specs.cols( pops_feature_t::POPS_HJORTH , siglab ) ;
	       X1( en , cols[0] ) = activity > 0 ? log( activity ) : log( 0.0001 ) ;
	       X1( en , cols[1] ) = mobility;
	       X1( en , cols[2] ) = complexity;	     
	     }

	   //
	   // Next signal
	   //
	}


      //
      // track that this was a bad epoch (for at least one signal/metric)
      //
      
      if ( bad_epoch ) S[en] = POPS_UNKNOWN;

      //
      // next epoch
      //
      
      ++en;
      
    } // next epoch


  //
  // Epoch-level outlier removal (at lvl1 stage) 
  //

  // get any lvl1 blocks that have been flagged by an OUTLIER command
  std::map<std::string,pops_spec_t> fm = pops_t::specs.fcmap[ pops_feature_t::POPS_EPOCH_OUTLIER ];
  
  std::map<std::string,pops_spec_t>::const_iterator ii = fm.begin();
  while ( ii != fm.end() )
    {
      // block name
      const std::string & blk = ii->first;
      // get param
      const pops_spec_t & spec = ii->second;
      
      // get lvl1 columns only (
      std::vector<int> c = pops_t::specs.block_cols( blk , pops_t::specs.n1 );
      // std::cout << " --> " << blk << " " 
      // 		<< c.size() << "\n";

      double th = spec.narg( "th" );

      // copy staging
      std::vector<int> S2 = S;
      
      // this amends S2
      for (int j=0;j<c.size();j++)
	{
	  pops_t::outliers( X1.col(c[j]) , th , S , &S2 );
	}

      // after doing one round for this block, update S
      S = S2;
      
      // do proc
      ++ii;
    }


  //
  // Prune out bad rows
  //

  // need to change S, E and X1
  std::vector<int> reslot;
  for (int i=0; i<S.size(); i++)
    {
      if ( S[i] != POPS_UNKNOWN )
	reslot.push_back( i );
    }

  int good = reslot.size();

  // 0 1 2 3 4 5 6 7 8
  // 0 1 . 2 . . 3 4 5

  // 0 1 2 3 4 5 6 7 8
  // 0 1 . 2 . . 3 4 5

  // 0 1 3 6 7 8
  
  for (int i=0; i<good; i++)
    {
      if ( reslot[i] != i )
	{
	  S[i] = S[ reslot[i] ];
	  E[i] = E[ reslot[i] ];
	  X1.row(i) = X1.row( reslot[i] );
	}
    }

  logger << "  pruning rows from " << ne << " to " << good << " epochs\n";

  // final update
  S.resize( good );
  E.resize( good );
  X1.conservativeResize( good , Eigen::NoChange );
  ne = good;
  
  //
  // all done
  //
}


void pops_indiv_t::level2()
{

  // co-opt pops_t::level2() to do this (i.e. same
  // code as used for trainers.   the only difference is 
  // that the SVD W/V will be read from the file, and a 
  // project done
  
  // need to set up duplicates in pops_t and hen copy back
  // bit of a kludge, but this is better than using
  // a duplicated copy of core level 2 features (i.e. if
  // we add stuff

  // expand X1 to include space for level-2 features                                                                                                                                          
  X1.conservativeResize( Eigen::NoChange , pops_t::specs.na );

  pops_t pops;
  pops.from_single_target( *this );
  pops.level2( false ); // false --> not training sample
  pops.copy_back( this );


}


void pops_indiv_t::predict()
{
  P = pops_t::lgbm.predict( X1 );
}

void pops_indiv_t::SHAP()
{
  
  Eigen::MatrixXd SHAP = pops_t::lgbm.SHAP_values( X1 );

  const int n_classes = pops_opt_t::n_stages;
  const int n_features = pops_t::specs.nf;
  
  //
  // always report get means 
  //

  Eigen::MatrixXd C = Eigen::MatrixXd::Zero( n_features , n_classes );
  
  std::vector<std::string> labels = pops_t::specs.select_labels();
  if ( labels.size() != n_features )
    Helper::halt( "internal error in getting labels" );

  int p = 0;
  for (int c = 0 ; c < n_classes ; c++)
    {
      if ( pops_opt_t::n_stages == 5 ) 
	writer.level( pops_t::labels5[ c ] , globals::stage_strat );
      else
	writer.level( pops_t::labels3[ c ] , globals::stage_strat );
      
      for (int r = 0; r < n_features ; r++)
	{
	  writer.level( labels[r] , globals::feature_strat );

	  C(r,c) = SHAP.col(p++).cwiseAbs().mean();
	  writer.value( "SHAP" , C(r,c) );
	}
      writer.unlevel( globals::feature_strat );
    }
  writer.unlevel( globals::stage_strat );


  //
  // Verbose mode: epoch level SHAP (i.e. do for single indiv)
  //
  
  if ( pops_opt_t::epoch_level_SHAP )
    {
      logger << "  reporting epoch-level SHAP values...\n";
      
      int p = 0;
      for (int c = 0 ; c < n_classes ; c++)
	{	  
	  if ( pops_opt_t::n_stages == 5 ) 
	    writer.level( pops_t::labels5[ c ] , globals::stage_strat );
	  else
	    writer.level( pops_t::labels3[ c ] , globals::stage_strat );
	  
	  for (int r = 0; r < n_features ; r++)
	    {
	      writer.level( labels[r] , globals::feature_strat );
	      
	      // epoch	     
	      for (int e = 0; e < ne ; e++) 
		{
		  writer.epoch( e+1 );
		  writer.value( "SHAP" , SHAP(e,p) );
		}
	      writer.unepoch();
	      
	      // next column for SHAP
	      p++;
	      
	    } // next feature
	  writer.unlevel( globals::feature_strat );
	}
      
      writer.unlevel( globals::stage_strat );
    }
  
}

void pops_indiv_t::summarize()
{
  
  std::map<int,double> dur_obs, dur_obs_orig, dur_predf, dur_pred1;
  std::vector<int> preds;
  
  int slp_lat_obs = -1 , slp_lat_prd = -1;
  int rem_lat_obs = -1 , rem_lat_prd = -1;
  
  //
  // epoch-level output (posteriors & predictions)
  //  
  
  double avg_pmax = 0;

  for (int e=0; e<ne; e++)
    {
      
      writer.epoch( E[e] + 1 );
      
      // format always: W R N1 N2 N3
      writer.value( "PP_W" , P(e,0) );   // 0 W 
      writer.value( "PP_R" , P(e,1) );   // 1 R
      writer.value( "PP_N1" , P(e,2) );  // 2 NR
      writer.value( "PP_N2" , P(e,3) );
      writer.value( "PP_N3" , P(e,4) );
      
      // prior
      writer.value( "PRIOR" , pops_t::label( (pops_stage_t)S[e] ) );
      
      // original - note, uses other index back to the orignal epoch-count
      //writer.value( "ORIG" , pops_t::label( (pops_stage_t)Sorig[ E[e] ] ) );
      
      // predicted (original)
      int predx;
      double pmax = P.row(e).maxCoeff(&predx);
      writer.value( "CONF" , pmax );
      avg_pmax += pmax;
      preds.push_back( predx );
      writer.value( "PRED" , pops_t::labels5[ predx ] ) ; 
  
      // slp/rem latency
      if ( slp_lat_obs == -1 && S[e] != POPS_WAKE && S[e] != POPS_UNKNOWN ) 
	slp_lat_obs = E[e];
      
      if ( slp_lat_prd == -1 && predx != POPS_WAKE )
	slp_lat_prd = E[e];

      if ( rem_lat_obs == -1 && S[e] == POPS_REM )
	rem_lat_obs = E[e] - slp_lat_obs ; 
      
      if ( rem_lat_prd == -1 && predx == POPS_REM )
	rem_lat_prd = E[e] - slp_lat_prd ; 

	
      // durations
      dur_obs[ S[e] ]++;
      dur_pred1[ predx ]++;

      for (int ss=0; ss< pops_opt_t::n_stages; ss++)
	dur_predf[ ss ] += P(e,ss);
    }

  writer.unepoch();

  //
  // Summaries
  //  

  // durations in minutes, so get scaling factor
  const double fac = pops_opt_t::epoch_len / 60.0 ; 
    
  // 5-class stats
  pops_stats_t stats( S, preds , 5 );

  // 3-class stats
  pops_stats_t stats3( pops_t::NRW( S ) , pops_t::NRW( preds ) , 3 );

  
  // outputs

  writer.value( "K" , stats.kappa );
  writer.value( "K3" , stats3.kappa );

  writer.value( "ACC" , stats.acc );
  writer.value( "ACC3" , stats3.acc );

  writer.value( "CONF" , avg_pmax / (double)ne );
    
  writer.value( "MCC" , stats.mcc );
  writer.value( "MCC3" , stats3.mcc );

  writer.value( "F1" , stats.macro_f1 );
  writer.value( "PREC" , stats.macro_precision );
  writer.value( "RECALL" , stats.macro_recall );
  
  writer.value( "F1_WGT" , stats.avg_weighted_f1 );
  writer.value( "PREC_WGT" , stats.avg_weighted_precision );
  writer.value( "RECALL_WGT" , stats.avg_weighted_recall );
  
  writer.value( "F13" , stats3.macro_f1 );
  writer.value( "PREC3" , stats3.macro_precision );
  writer.value( "RECALL3" , stats3.macro_recall );

  //
  // stage specific precision/recall
  //

  for ( int l=0;l<pops_opt_t::n_stages; l++)
    {
      writer.level( pops_t::labels5[l] , globals::stage_strat );
      writer.value( "F1" , stats.f1[l] );
      writer.value( "PREC" , stats.precision[l] );
      writer.value( "RECALL" , stats.recall[l] );
    }
  writer.unlevel( globals::stage_strat );


  //
  // sleep and REM latencies
  //

  if ( slp_lat_obs >= 0 ) 
    writer.value( "SLP_LAT_OBS" , slp_lat_obs * fac );
  if ( slp_lat_prd >= 0 ) 
    writer.value( "SLP_LAT_PRD" , slp_lat_prd * fac );
  if ( rem_lat_obs >= 0 )
    writer.value( "REM_LAT_OBS" , rem_lat_obs * fac );
  if ( rem_lat_prd >= 0 )
    writer.value( "REM_LAT_PRD" , rem_lat_prd * fac );


  //
  // Stage level durations
  //

  
  // unknown : dropped epochs going from 
  for (int e=0; e<Sorig.size(); e++)
    dur_obs_orig[ Sorig[e] ]++;
  
  for (int ss=0; ss < pops_opt_t::n_stages ; ss++ )
    {
      writer.level( pops_t::label( (pops_stage_t)ss ) , "SS" ); 
      writer.value( "OBS" ,  fac * dur_obs[ss] );
      writer.value( "ORIG" , fac * dur_obs_orig[ss] );
      writer.value( "PRF" ,  fac * dur_predf[ss] );
      writer.value( "PR1" ,  fac * dur_pred1[ss] );
    }
 
  int masked = Sorig.size() - S.size();
  writer.level( pops_t::label( POPS_UNKNOWN ) , "SS" );
  writer.value( "OBS" ,  fac * masked );
  writer.value( "ORIG" , fac * dur_obs_orig[ POPS_UNKNOWN ] );
  writer.value( "PRF" ,  fac * masked );
  writer.value( "PR1" ,  fac * masked );
  writer.unlevel( "SS" );

  //
  // Confusion matrix, to console
  //

  logger << "  \n  Final kappa = " << stats.kappa << "; 3-class kappa = " << stats3.kappa << "\n";
  logger << "  Confusion matrix: \n";
  std::map<int,std::map<int,int> > table = pops_t::tabulate( S , preds , true );
  logger << "\n";

}





#endif