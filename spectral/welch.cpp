
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

#include "spectral/welch.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "eval.h"
#include "annot/annot.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "db/db.h"
#include "fftw/fftwrap.h"
#include "dsp/mse.h"
#include "miscmath/dynam.h"
#include "fftw/bandaid.h"

extern writer_t writer;
extern logger_t logger;


annot_t * spectral_power( edf_t & edf , 
			  const std::string & signal_label , 
			  param_t & param )
{

  
  // Report full spectrum as well as band power
  const bool show_spectrum = param.has( "spectrum" ) || param.has("epoch-spectrum" );

  // do not report bands
  const bool bands = param.has( "band" ) ? param.yesno( "band" ) : true ; 
  
  // Report dB scale ?
  const bool dB = param.has( "dB" );
  
  // Min required SR to report
  const double min_sr = param.has( "min-sr" ) ? param.requires_dbl( "min-sr" ) : 50 ;

  // Mean center data first? 
  const bool mean_centre_epoch = param.has( "center" ) || param.has( "centre" )
    || param.has( "mean-center" ) || param.has( "mean-centre" );
    
  // detrend signal first?
  const bool remove_linear_trend = param.has( "detrend" );

  if ( mean_centre_epoch && remove_linear_trend )
    Helper::halt( "cannot specify both mean-center and detrend" );

  
  // Spectrum bin width (0 means no binning, default)

  //double bin_width = param.has( "bin" ) ? param.requires_dbl( "bin" ) : 0;
  const double bin_fac = param.has( "fac" ) ? param.requires_int( "fac" ) : 1;

  // Band power per-epoch
  const bool show_epoch = param.has( "epoch" ) || param.has("epoch-spectrum" );
  
  // report variance of PSD across epochs (SD)
  const bool aggregate_psd_sd       = param.has( "sd" );
  const double aggregate_psd_th     = param.has( "th" ) ? param.requires_dbl( "th" ) : 0 ; 
  const bool aggregate_psd_med      = param.has( "median" );

  // get kurtosis of band power (raw - kurt3, or excess kurtosis - kurt )  
  const bool calc_kurt  = param.has( "kurt" ) || param.has( "kurtosis" )
    || param.has( "kurt3" ) || param.has( "kurtosis3" ) ;

  const double kurt_adj = param.has( "kurt3" ) || param.has( "kurtosis3" ) ? +3 : 0 ; 

  // output ratios of band power values
  const bool calc_ratio = param.has( "ratio" );
  if ( calc_ratio && param.empty( "ratio" ) )
    Helper::halt( "cannot have empty ratio arg" );
  // ratio=ALPHA/BETA,THETA/DELTA,...
  const std::string ratios = calc_ratio ? param.value( "ratio" ) : "" ;
  const int ratio_plus1 = param.has( "ratio1" );
  
  // Characterize dynamics: of all epoch-level stats created, get the H1, H2, H3, linear and exponential trend
  // Hjorth stats just consider all points concatenated
  // trend lines are based on observed E numbers 

  const bool calc_dynamics = param.has( "dynamics" );
  
  // Verbose output: full spectrum per epoch

  const bool show_epoch_spectrum = param.has( "epoch-spectrum" );


  //
  // add new signals
  //  --> prefix_CH_N ... where N = 1,2,3, that correspond to Fs in range                                              
  //  --> note, PSD analysis is always done epochwise, which means that this
  //            should be run w/ small (e.g. 4 sec) epoch size and inc 2,
  //    i.e. then the time-series will be "per-epoch" which means 0.5 Hz SR
  //         for the new channel, etc
  //

  const bool new_sigs = param.has( "add" ) || param.has( "add-spectrum" );
  
  const bool new_sigs_relpow = param.has( "add-relpower" );

  if ( new_sigs_relpow && dB ) 
    Helper::halt( "cannot combine add-relpower and dB" );
  
  const bool new_spec_sigs = param.has( "add-spectrum" );
  
  const std::string new_sig_prefix = ( new_sigs && param.empty("add" ) ) ? "" : param.value( "add" ) ;

  std::set<std::string> new_sigs_skip_bands;
  if ( new_sigs && param.has( "skip-bands" ) )
    new_sigs_skip_bands = param.strset( "skip-bands" );
  
  // some extra requirements if adding a new signal
  if ( new_sigs )
    {

      if ( ! edf.header.continuous )
	Helper::halt( "currently, can only specify 'add' with continuous recordings" );
      
      if ( edf.timeline.generic_epochs() )
	Helper::halt( "can not have generic epochs with 'add'" );
      
      if ( edf.timeline.epoch_any_offset() )
	Helper::halt( "cannot use 'add' with any EPOCH offset (e.g. from align)" );

      // to save the signal, the EDF record size must be consistent w/ the new SR of the signal
      // i.e. to have same number of samples per EDF record, cannot have SR = 0.5 Hz with 1-second records, etc

      // perhaps for now, enfore that EDF record size is 1 second, and always output at one-second intervals
      if ( edf.header.record_duration_tp != globals::tp_1sec )
	Helper::halt( "currently, must have 1-second EDF records (use RECORD-SIZE)" );
      if ( edf.timeline.epoch_increment_tp() != globals::tp_1sec )
	Helper::halt( "currently, must have 1-second epoch increment (use EPOCH inc=1 len=4)" );
    }

  
  //
  // for PSD, if adding a new signal (based on EPOCH-level outputs)
  // we cannot have any gaps - i.e. to make life easier with uniform
  // SR of the new signals ;  we can relax this later, but for now
  // the work flow (for infraslow osc. analysis) would be to run
  // round 1 on continuous signal, then epoch (e.g. by stage)
  // then run second FFT; there, the EPOCH'ing will take care of
  // the lengths of segments, if use e.g. EPOCH annot=N2 ... i.e.
  // the epoch size will be extended as appropriate; note: we might want
  // to add an option for weighed average by window /epoch length in that case
  // as those will not all be equal
  //

  
  
  // peak diagnostics
  const bool peak_diagnostics = param.has( "peaks" )
    || param.has( "epoch-peaks" ) || param.has( "peaks-epoch" )
    || param.has( "peaks-verbose") || param.has( "peaks-frq" ) ;
  const int peak_median_filter_n = param.has( "peaks-window" ) ? param.requires_int( "peaks-window" ) : 11 ; 
  const bool verbose_peaks = param.has( "peaks-verbose" );
  const bool peak_per_epoch = param.has( "epoch-peaks" ) || param.has( "peaks-epoch" );
  std::vector<double> peak_range(2); peak_range[0] = 0 ; peak_range[1] = 99999;
  if ( param.has( "peaks-frq" ) ) peak_range = param.dblvector( "peaks-frq" );
  if ( peak_range.size() != 2 || peak_range[0] >= peak_range[1] ) Helper::halt( "bad peaks-frq=lwr,upr" );

  //
  // spectral slope
  //

  const bool spectral_slope = param.has( "slope" );
  const std::vector<double> slope_range = param.dblvector( "slope" );
  const bool spectral_slope_show_epoch = param.has( "epoch-slope" ) || param.has( "slope-epoch" );
 
  if ( spectral_slope )
    {
      if ( slope_range.size() != 2 ||  
	   slope_range[0] >= slope_range[1] ||
	   slope_range[0] <= 0 ||
	   slope_range[1] <= 0 )
	Helper::halt( "expecting slope=lwr,upr" );
    }

  // outlier threshold to remove individual PSD points when calculating a single slope
  const double slope_outlier = param.has( "slope-th" ) ? param.requires_dbl( "slope-th" ) : 3 ;

  // threshold to reove epochs when summarizing slopes over all epochs 
  const double slope_th2     = param.has( "slope-th2" ) ? param.requires_dbl( "slope-th2" ) : 3 ;
  
  // truncate spectra
  double min_power = param.has( "min" ) ? param.requires_dbl( "min" ) : 0.5 ;
  double max_power = param.has( "max" ) ? param.requires_dbl( "max" ) : 25 ;

  // check that slope=X,Y or peaks=X,Y does not necessitate and expanded range
  if ( slope_range.size() == 2 ) 
    {
      if ( min_power > slope_range[0] ) min_power = slope_range[0];
      if ( max_power < slope_range[1] ) max_power = slope_range[1];
    }

  if ( param.has( "peaks-frq" ) )
    {
      if ( min_power > peak_range[0] ) min_power = peak_range[0];
      if ( max_power < peak_range[1] ) max_power = peak_range[1];
    }

  // Calculate MSE
  const bool calc_mse = param.has( "mse" ); 

  //
  // cache PSD for other analyses (e.g. PSC, ASYMM)
  //
  
  const bool cache_data = param.has( "cache" );
  
  const std::string cache_name = cache_data ? param.requires( "cache" ) : "" ;
  
  const bool cache_epochs = param.has( "cache-epochs" );
  
  const bool cache_bands = param.has( "cache-bands" ) ? param.yesno( "cache-bands" ) : true ;
  
  const bool cache_spectrum = param.has( "cache-spectra" ) ? param.yesno( "cache-spectra" ) : false ;
  
  if ( ( cache_epochs || cache_epochs || cache_spectrum ) && ! cache_data ) 
    Helper::halt( "must specify cache=name with cache-epochs, cache-bands or cache-spectra" );
  
  cache_t<double> * cache = NULL ; 

  if ( cache_data )
    cache = edf.timeline.cache.find_num( cache_name );

  // i.e. to only make cache values
  const bool suppress_output = param.has( "silent" );
 
	  
  //
  // Alter PWELCH sliding window parameters
  //

  double fft_segment_size = param.has( "segment-sec" ) 
    ? param.requires_dbl( "segment-sec" ) : 4 ;

  // allow both versions for MTM and backwards compatibility
  double fft_segment_overlap = 2;
  if ( param.has( "segment-inc" ) )
    fft_segment_overlap = param.requires_dbl( "segment-inc" );
  else if ( param.has( "segment-overlap" ) )
    fft_segment_overlap = param.requires_dbl( "segment-overlap" ); 


  //
  // If adding a signal, require that we set segment size/inc to equal epoch size/inc
  //

  if ( new_sigs )
    {

      fft_segment_size = edf.timeline.epoch_length();
      fft_segment_overlap = 0;
      logger << "  with 'add', using epoch duration to set segment-sec=" << fft_segment_size << " and forcing segement-inc=0\n";
    }  
  
  
  //
  // Option to average adjacent points in the power spectra (default = T)
  //
  
  const bool average_adj = param.has( "average-adj" ) ;

  
  //
  // Window function
  //
  
  window_function_t window_function = WINDOW_TUKEY50;	   
  if      ( param.has( "no-window" ) ) window_function = WINDOW_NONE;
  else if ( param.has( "hann" ) ) window_function = WINDOW_HANN;
  else if ( param.has( "hamming" ) ) window_function = WINDOW_HAMMING;
  else if ( param.has( "tukey50" ) ) window_function = WINDOW_TUKEY50;

  //
  // Median vs mean to get eppch PSD (i.e. averaging over segments in Welch)
  //

  const bool use_seg_median = param.has( "segment-median" );

  //
  // Return intra-segment CVs and associated stats
  //

  const bool calc_seg_sd = param.has( "segment-sd" );
  
  //
  // Use nextpow2 for NFFT
  //

  const bool use_nextpow2 = param.has( "pow2" );

  //
  // change power band definitions on-the-fly
  //

  bandaid_t bandaid;

  bandaid.define_bands( param );
  
  
  //
  // Attach signals
  //
  
  signal_list_t signals = edf.header.signal_list( signal_label );  
  
  const int ns = signals.size();

  
  //
  // Obtain sampling freqs (Hz)
  //
  
  std::vector<double> Fs = edf.header.sampling_freq( signals );
  
  //
  // Set first epoch
  //
  
  edf.timeline.first_epoch();
    

  //
  // Check segment lengths (i.e. epoch sizes are fixed)
  //
  
  if ( ! edf.timeline.generic_epochs() )
    {
      if ( edf.timeline.epoch_length() <= ( fft_segment_size + fft_segment_overlap ) )
	{
	  fft_segment_overlap = 0;
	  fft_segment_size = edf.timeline.epoch_length();
	}
    }

  
  //
  // Initiate output
  //
  

  bool epoch_level_output = show_epoch || show_epoch_spectrum || peak_per_epoch || spectral_slope_show_epoch ;

  
  
  //
  // Get each signal
  //
  
  logger << "  calculating PSD from " << min_power << " to " << max_power << " for " << ns << " signals\n"; 

  for (int s = 0 ; s < ns; s++ )
    {
      
      //
      // only consider data tracks
      //
      
      if ( edf.header.is_annotation_channel( signals(s) ) ) 
	continue;

      //
      // bad SR ( require at least 50 Hz)
      //

      if ( Fs[s] < min_sr ) continue;

      //
      // reset bandaid
      //

      bandaid.init();
      
      //
      // Stratify output by channel
      //

      writer.level( signals.label(s) , globals::signal_strat );
      
      
      //
      // get high, low and total power.
      //
      
      int total_epochs = 0;  
      
      // store frequencies from epoch-level analysis
      std::vector<double> freqs;

      std::vector<double> epochs;

      // track F results      
      std::map<int,std::vector<double> > track_freq;
      std::map<int,std::vector<double> > track_freq_logged;
      
      std::map<int,std::vector<double> > track_segcv;
    
      // store spectral slope per epoch for this channel?
      std::vector<double> slopes;
      std::vector<double> slopes_intercept;
      std::vector<double> slopes_rsq;
      
      
      //
      // Set first epoch
      //
      
      edf.timeline.first_epoch();


      //
      // for each each epoch 
      //
      
      while ( 1 ) 
	{
	  
	  int epoch = edf.timeline.next_epoch();      
	  
	  if ( epoch == -1 ) break;
	  
	  interval_t interval = edf.timeline.epoch( epoch );
	  
         
	  //
	  // Need to check segment length?
	  //

	  if ( edf.timeline.generic_epochs() )
	    {
	      // here, all epochs need to have the same segment length, so will skip here
	      // if the epoch is too short
	      //if ( edf.timeline.epoch_length() < ( fft_segment_size + fft_segment_overlap ) )
	      if ( edf.timeline.epoch_length() < fft_segment_size )
		{
		  logger << "  *** skipping epoch " << interval.as_string() << ", too short given segment-sec\n";
		  continue;
		}
	    }
	  

	  //
	  // okay to process
	  //

	  ++total_epochs;
	  
	  //
	  // stratify output by epoch?
	  //

	  if ( epoch_level_output )
	    writer.epoch( edf.timeline.display_epoch( epoch ) );
	  
	  //
	  // Get data
	  //

	  slice_t slice( edf , signals(s) , interval );
	   
	   std::vector<double> * d = slice.nonconst_pdata();

	   //
	   // mean centre epoch?
	   //

	   if ( mean_centre_epoch ) 
	     MiscMath::centre( d );
	   else if ( remove_linear_trend )
	     MiscMath::detrend( d );
	   
	   //
	   // pwelch() to obtain full PSD
	   //
	   
	   const double overlap_sec = fft_segment_overlap;
	   const double segment_sec  = fft_segment_size;
	   
	   const int total_points = d->size();
	   const int segment_points = segment_sec * Fs[s];
	   const int noverlap_points  = overlap_sec * Fs[s];
	   
	   // implied number of segments
	   int noverlap_segments = floor( ( total_points - noverlap_points) 
					  / (double)( segment_points - noverlap_points ) );
	   	   
 	   // logger << "total_points = " << total_points << "\n";
 	   // logger << "nooverlap_segments = " << noverlap_segments << "\n";
 	   // logger << "noverlap_points = " << noverlap_points << "\n";
 	   // logger << "segment_points = " << segment_points << "\n\n";
	   	   
	   PWELCH pwelch( *d , 
			  Fs[s] , 
			  segment_sec , 
			  noverlap_segments , 
			  window_function , 
			  use_seg_median,
			  calc_seg_sd,
			  average_adj ,
			  use_nextpow2 );
	   
	   bandaid.track_bands_per_epoch( pwelch.psdsum( SLOW ),
					  pwelch.psdsum( DELTA ),
					  pwelch.psdsum( THETA ),
					  pwelch.psdsum( ALPHA ),
					  pwelch.psdsum( SIGMA ),
					  pwelch.psdsum( LOW_SIGMA ),
					  pwelch.psdsum( HIGH_SIGMA ),
					  pwelch.psdsum( BETA ),
					  pwelch.psdsum( GAMMA ),
					  pwelch.psdsum( DENOM ) );
	   
	   // double this_slowwave   = pwelch.psdsum( SLOW )  ;      /// globals::band_width( SLOW );
	   // double this_delta      = pwelch.psdsum( DELTA ) ;      /// globals::band_width( DELTA );
	   //  double this_theta      = pwelch.psdsum( THETA ) ;      /// globals::band_width( THETA );
	   //  double this_alpha      = pwelch.psdsum( ALPHA ) ;      /// globals::band_width( ALPHA );
	   //  double this_sigma      = pwelch.psdsum( SIGMA ) ;      /// globals::band_width( SIGMA );
	   //  double this_low_sigma  = pwelch.psdsum( LOW_SIGMA ) ;  /// globals::band_width( LOW_SIGMA );
	   //  double this_high_sigma = pwelch.psdsum( HIGH_SIGMA ) ; /// globals::band_width( HIGH_SIGMA );
	   //  double this_beta       = pwelch.psdsum( BETA )  ;      /// globals::band_width( BETA );
	   //  double this_gamma      = pwelch.psdsum( GAMMA ) ;      /// globals::band_width( GAMMA );]
	   //  double this_total      = pwelch.psdsum( DENOM ) ;      /// globals::band_width( DENOM );

	   
	   
	   //
	   // track epoch numbers (for dynam_t)
	   //
	   
	   epochs.push_back( epoch );

	   //
	   // Epoch-level output
	   //

	   
	   if ( show_epoch || ( cache_epochs && cache_bands ) )
	     {
	       if ( bands && bandaid.total > 0 )
		 {
		   writer.level( globals::band( SLOW ) , globals::band_strat );
		   if ( show_epoch && ! suppress_output ) {
		     writer.value( "PSD" , dB ? 10*log10( bandaid.slow ) : bandaid.slow  );
		     writer.value( "RELPSD" , bandaid.slow / bandaid.total );
		   }
		   if ( cache_data && cache_epochs && cache_bands )
		     cache->add( ckey_t( "PSD" , writer.faclvl() ) , dB ? 10*log10( bandaid.slow ) : bandaid.slow );
		   
		   writer.level( globals::band( DELTA ) , globals::band_strat );
		   if ( show_epoch && ! suppress_output ) {
		     writer.value( "PSD" , dB ? 10*log10( bandaid.delta ) : bandaid.delta );
		     writer.value( "RELPSD" , bandaid.delta / bandaid.total );
		   }
		   if ( cache_data && cache_epochs && cache_bands )
		     cache->add( ckey_t( "PSD" , writer.faclvl() ) , dB ? 10*log10( bandaid.delta ) : bandaid.delta );
		   
		   writer.level( globals::band( THETA ) , globals::band_strat );
		   if ( show_epoch && ! suppress_output ) {
		     writer.value( "PSD" , dB ? 10*log10( bandaid.theta ) : bandaid.theta  );
		     writer.value( "RELPSD" , bandaid.theta / bandaid.total );
		   }
		   if ( cache_data && cache_epochs && cache_bands )
		     cache->add( ckey_t( "PSD" , writer.faclvl() ) , dB ? 10*log10( bandaid.theta ) : bandaid.theta );
		   
		   writer.level( globals::band( ALPHA ) , globals::band_strat );
		   if ( show_epoch && ! suppress_output ) {
		     writer.value( "PSD" , dB ? 10*log10( bandaid.alpha ) : bandaid.alpha );
		     writer.value( "RELPSD" , bandaid.alpha / bandaid.total );
		   }
		   if ( cache_data && cache_epochs && cache_bands )
		     cache->add( ckey_t( "PSD" , writer.faclvl() ) , dB ? 10*log10( bandaid.alpha ) : bandaid.alpha );
		   
		   writer.level( globals::band( SIGMA ) , globals::band_strat );
		   if ( show_epoch && ! suppress_output ) {
		     writer.value( "PSD" , dB ? 10*log10( bandaid.sigma ) : bandaid.sigma );
		     writer.value( "RELPSD" , bandaid.sigma / bandaid.total );
		   }
		   if ( cache_data && cache_epochs && cache_bands )
		     cache->add( ckey_t( "PSD" , writer.faclvl() ) , dB ? 10*log10( bandaid.sigma ) : bandaid.sigma );
		   
		   writer.level( globals::band( LOW_SIGMA ) , globals::band_strat );
		   if ( show_epoch && ! suppress_output ) {
		     writer.value( "PSD" , dB ? 10*log10( bandaid.low_sigma ) : bandaid.low_sigma );
		     writer.value( "RELPSD" , bandaid.low_sigma / bandaid.total );
		   }
		   if ( cache_data && cache_epochs && cache_bands )
		     cache->add( ckey_t( "PSD" , writer.faclvl() ) , dB ? 10*log10( bandaid.low_sigma ) : bandaid.low_sigma );
		   
		   writer.level( globals::band( HIGH_SIGMA ) , globals::band_strat );
		   if ( show_epoch && ! suppress_output ) {
		     writer.value( "PSD" , dB ? 10*log10( bandaid.high_sigma ) : bandaid.high_sigma );
		     writer.value( "RELPSD" , bandaid.high_sigma / bandaid.total );
		   }
		   if ( cache_data && cache_epochs && cache_bands )
		     cache->add( ckey_t( "PSD" , writer.faclvl() ) , dB ? 10*log10( bandaid.high_sigma ) : bandaid.high_sigma );
		   

		   writer.level( globals::band( BETA ) , globals::band_strat );
		   if ( show_epoch && ! suppress_output ) {
		     writer.value( "PSD" , dB ? 10*log10( bandaid.beta ) : bandaid.beta  );
		     writer.value( "RELPSD" , bandaid.beta / bandaid.total );
		   }
		   if ( cache_data && cache_epochs && cache_bands )
		     cache->add( ckey_t( "PSD" , writer.faclvl() ) , dB ? 10*log10( bandaid.beta ) : bandaid.beta );
		   
		   writer.level( globals::band( GAMMA ) , globals::band_strat );
		   if ( show_epoch && ! suppress_output ) {		     
		     writer.value( "PSD" , dB ? 10*log10( bandaid.gamma ) : bandaid.gamma );
		     writer.value( "RELPSD" , bandaid.gamma / bandaid.total );
		   }
		   if ( cache_data && cache_epochs && cache_bands )
		     cache->add( ckey_t( "PSD" , writer.faclvl() ) , dB ? 10*log10( bandaid.gamma ) : bandaid.gamma );
		   
		   writer.level( globals::band( DENOM ) , globals::band_strat );
		   if ( show_epoch && ! suppress_output ) {
		     writer.value( "PSD" , dB ? 10*log10( bandaid.total ) : bandaid.total );				   
		   }
		   if ( cache_data && cache_epochs && cache_bands )
		     cache->add( ckey_t( "PSD" , writer.faclvl() ) , dB ? 10*log10( bandaid.total ) : bandaid.total );
		   
		   writer.unlevel( globals::band_strat );
		 }
	       else if ( bands && cache_data && cache_epochs && cache_bands && ! dB )
		 {
	       // need to enter 0 in this case for cache
	       //  nb. only doing this in non-dB mode (i.e. for ASYMM)
	       writer.level( globals::band( SLOW ) , globals::band_strat );
	       cache->add( ckey_t( "PSD" , writer.faclvl() ) , 0 );
	       
	       writer.level( globals::band( DELTA ) , globals::band_strat );
	       cache->add( ckey_t( "PSD" , writer.faclvl() ) , 0 );
	       
	       writer.level( globals::band( THETA ) , globals::band_strat );
	       cache->add( ckey_t( "PSD" , writer.faclvl() ) , 0 );
	       
	       writer.level( globals::band( ALPHA ) , globals::band_strat );
	       cache->add( ckey_t( "PSD" , writer.faclvl() ) , 0 );
	       
	       writer.level( globals::band( SIGMA ) , globals::band_strat );
	       cache->add( ckey_t( "PSD" , writer.faclvl() ) , 0 );
	       
	       writer.level( globals::band( LOW_SIGMA ) , globals::band_strat );
	       cache->add( ckey_t( "PSD" , writer.faclvl() ) , 0 );

	       writer.level( globals::band( HIGH_SIGMA ) , globals::band_strat );
	       cache->add( ckey_t( "PSD" , writer.faclvl() ) , 0 );

	       writer.level( globals::band( BETA ) , globals::band_strat );
	       cache->add( ckey_t( "PSD" , writer.faclvl() ) , 0 );
	       
	       writer.level( globals::band( GAMMA ) , globals::band_strat );
	       cache->add( ckey_t( "PSD" , writer.faclvl() ) , 0 );
	       
	       writer.level( globals::band( DENOM ) , globals::band_strat );
	       cache->add( ckey_t( "PSD" , writer.faclvl() ) , 0 );
	       
	       writer.unlevel( globals::band_strat );
	       
	     }
	}
      
    
	   //
	   // track over entire spectrum (track on first encounter)
	   //
	   
	   if( freqs.size() == 0 ) 
	     {
	       freqs = pwelch.freq;
	     }
	   
	   
	   if ( freqs.size() == pwelch.psd.size() )
	     {
	       
	       //
	       // accumulate for entire night means; store as dB and raw
	       //
	       
	       if ( show_spectrum || spectral_slope || new_sigs )
		 for (int f=0;f<pwelch.psd.size();f++)
		   {
		     track_freq[ f ].push_back( pwelch.psd[f] );
		     // only track if power > 0
		     if ( pwelch.psd[f] > 0 ) 
		       track_freq_logged[ f ].push_back( 10*log10( pwelch.psd[f] ) );
		   }
	       

	       //
	       // Segment-level stats?
	       //

	       if ( calc_seg_sd ) 
		 for (int f=0;f<pwelch.psd.size();f++)
		   {		   
		     if ( pwelch.psd[f] > 0 ) 
		       track_segcv[ f ].push_back( pwelch.psdsd[f] );
		   }
	       

	       //
	       // epoch-level output?
	       //
	       
	       
	       if ( show_epoch_spectrum || ( cache_epochs && cache_spectrum ) )
		 {		 
		   
		   std::vector<double> f0;
		   
		   // using bin_t 	      
		   bin_t bin( min_power , max_power , bin_fac );
		   bin.bin( freqs , pwelch.psd );
		   
		   bin_t binsd( min_power , max_power , bin_fac );
		   if ( calc_seg_sd )
		     binsd.bin( freqs, pwelch.psdsd );
		   
		   
		   for ( int i = 0 ; i < bin.bfa.size() ; i++ )
		     {		     
		       f0.push_back( ( bin.bfa[i] + bin.bfb[i] ) / 2.0 );		       
		       writer.level( f0[ f0.size()-1 ] , globals::freq_strat );
		       
		       //writer.level( bin.bfa[i] , globals::freq_strat );
		       if ( show_epoch_spectrum && ! suppress_output ) 
			 if ( bin.bspec[i] > 0 || ! dB ) 
			   writer.value( "PSD" , dB? 10*log10( bin.bspec[i] ) : bin.bspec[i] );
		       
		       // cache epoch-level PSD?
		       if ( cache_epochs && cache_spectrum )
			 cache->add( ckey_t( "PSD" , writer.faclvl() ) , dB? 10*log10( bin.bspec[i] ) : bin.bspec[i] );
		       
		       if ( show_epoch_spectrum &&  ! suppress_output ) 
			 if ( bin.nominal[i] != "" )
			   writer.value( "INT" , bin.nominal[i] );

		       // CV?
		       if ( show_epoch_spectrum && ! suppress_output )
			 if ( calc_seg_sd )
			   writer.value( "CV" , binsd.bspec[i] );

		     }
		   writer.unlevel( globals::freq_strat );
		 }

	       
	       //
	       // epoch-level peakedness
	       //
	       
	       if ( peak_per_epoch )
		 {
		   peakedness( pwelch.psd , pwelch.freq , peak_median_filter_n , peak_range , false );
		 }

	       //
	       // epoch-level spectral slope?
	       //

	       if ( spectral_slope ) 
		 {		   
		   
		   double es1 = 0 , intercept = 0 , rsq = 0 ;
		   
		   bool okay = spectral_slope_helper( pwelch.psd ,
						      pwelch.freq ,
						      slope_range ,
						      slope_outlier ,
						      spectral_slope_show_epoch , 
						      &es1 , NULL, &intercept, &rsq );
		   
		   if ( okay )
		     {
		       slopes.push_back( es1 );
		       slopes_intercept.push_back( intercept );
		       slopes_rsq.push_back( rsq );
		     }
		 }
	       
	     }
	   else
	     logger << " *** warning:: skipped a segment: different NFFT/internal problem ... \n";
	   
      
	   //
	   // end of epoch-level strata
	   //

	   if ( epoch_level_output )
	     writer.unepoch();
	   
	   //
	   // next epoch
	   //

	}
      
      
      //
      // Output
      //
      
      const int n = freqs.size();      

      bool okay = total_epochs > 0 ;
      
      if ( ! suppress_output )       
	writer.value( "NE" , total_epochs );
      
      //
      // report full spectrum, or calculate statistics based on the full 
      // spectrum
      //

      if ( okay && ( show_spectrum || peak_diagnostics || spectral_slope ) ) 
	{	  
	      	  
	  //
	  // get mean power across epochs
	  //
	  
	  if ( track_freq.size() != freqs.size() ) 
	    {
	      std::cerr << "track_freq = " << track_freq.size() << " vs freqs = " << freqs.size() << "\n";
	      Helper::halt( "internal error psd_t" );
	    }
	  
	  std::vector<double> means, medians, sds;
	  
	  std::vector<double> cv_means, cv_medians, cv_sds; 

	  int ne_valid = dB ? track_freq_logged[0].size() : track_freq[0].size();
	  int ne_min = ne_valid;
	
	  for (int f=0;f<n;f++) 
	    {
	      
	      //	      std::cout << " F = " << f << "\n";

	      // wanting to get stats of dB or raw?
	      const std::vector<double> & yy = dB ? track_freq_logged[f] : track_freq[f] ;
	      
	      // any outlier removal of epochs?  	      
	      std::vector<double> xx = aggregate_psd_th > 0 && ne_valid > 2 ? 
		MiscMath::outliers( &yy , aggregate_psd_th ) : yy ; 
	      
	      // track min size
	      if ( xx.size() < ne_min ) 
		ne_min = xx.size();
	      
	      const double epoch_mean = MiscMath::mean( xx ) ;
	      const double epoch_sd   = MiscMath::sdev( xx ) ; 
	      	      
	      means.push_back( epoch_mean );
	      
	      if ( aggregate_psd_sd && xx.size() > 2 )
		sds.push_back( epoch_sd ) ;
	      
	      if ( aggregate_psd_med && xx.size() > 2 )
		medians.push_back(  MiscMath::median( xx ) );

	      // segment CV tracking?
	      
	      if ( calc_seg_sd )
		{
		  const std::vector<double> & yy = track_segcv[f];
		  
		  // any outlier removal of epochs?  	      
		  std::vector<double> xx = aggregate_psd_th > 0 && ne_valid > 2 ? 
		    MiscMath::outliers( &yy , aggregate_psd_th ) : yy ; 

		  //		  std::cout << " xx s = " << xx.size() << "\n";

		  cv_means.push_back( MiscMath::mean( xx ) );
		  
		  if ( xx.size() > 2 )
		    cv_sds.push_back( MiscMath::sdev( xx ) );
		  
		  if ( xx.size() > 2 )
		    cv_medians.push_back(  MiscMath::median( xx ) );	      

		} 

	    }
	
	  bin_t bin( min_power , max_power , bin_fac );	  
	  bin.bin( freqs , means );
	  //	  std::cout << " frq means " << freqs.size() << " " << means.size() <<"\n";

	  bin_t bin_med( min_power , max_power , bin_fac );	  
	  if ( aggregate_psd_med && ne_min > 2 )
	    bin_med.bin( freqs , medians );

	  bin_t bin_sds( min_power , max_power , bin_fac );	  
	  if ( aggregate_psd_sd && ne_min > 2 ) 
	    bin_sds.bin( freqs , sds );

	  // segment CV
	  bin_t cv_bin( min_power , max_power , bin_fac );
          bin_t cv_bin_med( min_power , max_power , bin_fac );
          bin_t cv_bin_sds( min_power , max_power , bin_fac );
	  
	  if ( calc_seg_sd )
	    {

	      // std::cout <<"freqs " << freqs.size() << " " << " " << cv_means.size() << " " << cv_medians.size() 
	      //  		<< " " << cv_sds.size() << "\n";
	      
	      cv_bin.bin( freqs , cv_means );
	      
	      if ( ne_min > 2 )
		cv_bin_med.bin( freqs , cv_medians );
	      
	      if ( ne_min > 2 )
		cv_bin_sds.bin( freqs , cv_sds );

	    }
	  
	  //
	  // Get total power
	  //

	  double tot_pow_denom = 0; 
	  for ( int i = 0 ; i < bin.bfa.size() ; i++ )
	    tot_pow_denom += dB ? pow( 10 , bin.bspec[i]/10.0) : bin.bspec[i] ;
	  
	  
	  //
	  // Output raw and relative power spectra
	  //

	  std::vector<double> f0;

	  for ( int i = 0 ; i < bin.bfa.size() ; i++ ) 
	    {
	  
	      f0.push_back( ( bin.bfa[i] + bin.bfb[i] ) / 2.0 );
	      
	      double x = bin.bspec[i] ;
	      
	      if ( show_spectrum ) 
		{
		  writer.level( ( bin.bfa[i] + bin.bfb[i] ) / 2.0 , globals::freq_strat );
		  //writer.level( bin.bfa[i] , globals::freq_strat );
		  
		  // these stats will aleady be dB-scaled if requested
		  if ( ! suppress_output )
		    {
		      writer.value( "PSD" , x );
		      writer.value( "RELPSD" , x / tot_pow_denom );
		      
		      if ( aggregate_psd_med && ne_min > 2 )
			writer.value( "PSD_MD" , bin_med.bspec[i]  );
		      
		      if ( aggregate_psd_sd && ne_min > 2 )
			{
			  writer.value( "PSD_SD" , bin_sds.bspec[i]  );

			  // also give CV assuming log-normal distribution
			  if ( dB ) 
			    {
			      // need SD on natural log scale
			      // bin_sds.bspec[i] is on 10 * log10( X ) scale
			      const double lnsd = log( 10.0 ) * bin_sds.bspec[i] / 10.0 ;
			      const double cv = sqrt( exp( lnsd ) - 1.0 ) ; 
			      writer.value( "PSD_CV" , cv );
			    }
			}
		      

		      if ( calc_seg_sd )
			{
			  writer.value( "SEGCV_MN" , cv_bin.bspec[i]  );
			  writer.value( "SEGCV_MD" , cv_bin_med.bspec[i]  );
			  writer.value( "SEGCV_SD" , cv_bin_sds.bspec[i]  );
			}

		      if ( bin.nominal[i] != "" )
			writer.value( "INT" , bin.nominal[i] );
		    }
		}
	      
	    }
	  
	  if ( show_spectrum )
	    writer.unlevel( globals::freq_strat );
	  
	  
	  std::vector<double> raw_mean_psd = bin.bspec;
	  
	  // kludge... urgh
	  if ( dB ) 
	    {	      
	      for (int i=0; i<raw_mean_psd.size(); i++) 
		raw_mean_psd[i] = pow( 10 , raw_mean_psd[i]/10.0) ;	      
	    }
	  
	  //
	  // Report metrics on the PSD: expecting a raw PSD
	  //
	  
	  if ( peak_diagnostics )
	    peakedness( raw_mean_psd , f0 , peak_median_filter_n , peak_range , verbose_peaks ); 

	  
	  //
	  // spectral slope? (expecting a raw PSD) 
	  //
	
	  if ( spectral_slope ) 
	    {		   
	      spectral_slope_helper( raw_mean_psd , 
				     f0 , 
				     slope_range , 
				     slope_outlier );
	    }
	  
	  
	}

      
      //
      // output spectral slope based on distribution of epoch-level slopes?
      //

      if ( spectral_slope && ! suppress_output ) 
	{
	  if ( slopes.size() > 2 )
	    {

	      // slope
	      std::vector<double> s2 = MiscMath::outliers( &slopes , slope_th2 );
	      if ( s2.size() != 0 ) 
		{
		  double s_mean = MiscMath::mean( s2 );
		  double s_med  = MiscMath::median( s2 );
		  double s_sd   = MiscMath::sdev( s2 , s_mean );		  
		  writer.value( "SPEC_SLOPE_MN" , s_mean );
		  writer.value( "SPEC_SLOPE_MD" , s_med );
		  writer.value( "SPEC_SLOPE_SD" , s_sd );		  
		}
	      
	      // intercept
	      std::vector<double> i2 = MiscMath::outliers( &slopes_intercept , slope_th2 );
	      if ( i2.size() != 0 ) 
		{
		  double i_mean = MiscMath::mean( i2 );
		  double i_med  = MiscMath::median( i2 );
		  double i_sd   = MiscMath::sdev( i2 , i_mean );
		  writer.value( "SPEC_INTERCEPT_MN" , i_mean );
		  writer.value( "SPEC_INTERCEPT_MD" , i_med );
		  writer.value( "SPEC_INTERCEPT_SD" , i_sd );
		}

	      // R-sq
	      std::vector<double> rsq2 = MiscMath::outliers( &slopes_rsq , slope_th2 );
	      if ( rsq2.size() != 0 )
		{
		  double rsq_mean = MiscMath::mean( rsq2 );
		  double rsq_med  = MiscMath::median( rsq2 );	      
		  writer.value( "SPEC_RSQ_MN" , rsq_mean );
		  writer.value( "SPEC_RSQ_MD" , rsq_med );
		}

	    }
	}

      //
      // mean total power
      //

      double mean_total_power = MiscMath::mean( bandaid.track_band[ DENOM ] );
      

      //
      // by band 
      //      

      if ( bands )
	{
	  std::vector<frequency_band_t>::const_iterator bi = bandaid.bands.begin();
	  while ( bi != bandaid.bands.end() )
	    {	   
	      
	      if ( okay ) 
		{
		  double p = MiscMath::mean( bandaid.track_band[ *bi ] );
		  writer.level( globals::band( *bi ) , globals::band_strat );
		  
		  if ( ! suppress_output ) {		
		    writer.value( "PSD" , dB ? 10*log10(p) : p  );
		    writer.value( "RELPSD" , p / mean_total_power );
		  }
		  
		}
	      
	      ++bi;
	    }
	  
	  writer.unlevel( globals::band_strat );
	}
          
      
      //
      // Dynamics?
      //
      

      if ( calc_dynamics )
	{
	  
	  // do we have any _CYCLE epoch-annotations ?
	  
	  bool has_cycles = edf.timeline.epoch_annotation( "_NREMC_1" ) 
	    || edf.timeline.epoch_annotation( "_NREMC_2" )
	    || edf.timeline.epoch_annotation( "_NREMC_3" )
	    || edf.timeline.epoch_annotation( "_NREMC_4" )
	    || edf.timeline.epoch_annotation( "_NREMC_5" )
	    || edf.timeline.epoch_annotation( "_NREMC_6" )
	    || edf.timeline.epoch_annotation( "_NREMC_7" )
	    || edf.timeline.epoch_annotation( "_NREMC_8" )
	    || edf.timeline.epoch_annotation( "_NREMC_9" )
	    || edf.timeline.epoch_annotation( "_NREMC_10" );
	  
	  std::vector<std::string> cycle;
	  
	  if ( has_cycles )
	    {
	      
	      for (int e=0;e<epochs.size(); e++)
		{
		  
		  std::string c = "."; // null
		  
		  // nb. uses current epoch encoding
		  // take up to 10 cycles
		  if      ( edf.timeline.epoch_annotation( "_NREMC_1" , epochs[e] ) ) c = "C1";
		  else if ( edf.timeline.epoch_annotation( "_NREMC_2" , epochs[e] ) ) c = "C2";
		  else if ( edf.timeline.epoch_annotation( "_NREMC_3" , epochs[e] ) ) c = "C3";
		  else if ( edf.timeline.epoch_annotation( "_NREMC_4" , epochs[e] ) ) c = "C4";
		  else if ( edf.timeline.epoch_annotation( "_NREMC_5" , epochs[e] ) ) c = "C5";
		  else if ( edf.timeline.epoch_annotation( "_NREMC_6" , epochs[e] ) ) c = "C6";
		  else if ( edf.timeline.epoch_annotation( "_NREMC_7" , epochs[e] ) ) c = "C7";
		  else if ( edf.timeline.epoch_annotation( "_NREMC_8" , epochs[e] ) ) c = "C8";
		  else if ( edf.timeline.epoch_annotation( "_NREMC_9" , epochs[e] ) ) c = "C9";
		  else if ( edf.timeline.epoch_annotation( "_NREMC_10" , epochs[e] ) ) c = "C10";
		  
		  cycle.push_back( c );
		  
		}
	    }
	  

	  //
	  // band power 
	  //

	  if ( bands )
	    {
	      std::map<frequency_band_t,std::vector<double> >::const_iterator ii = bandaid.track_band.begin();
	      
	      while ( ii != bandaid.track_band.end() )
		{	      
		  writer.level( globals::band( ii->first ) , globals::band_strat );
		  
		  if ( has_cycles )
		    dynam_report_with_log( ii->second , epochs , &cycle );
		  else
		    dynam_report_with_log( ii->second , epochs );
		  
		  ++ii;
		}
	      
	      writer.unlevel( globals::band_strat ); 
	    }
	  

	  //
	  // full spectra?
	  //
	  
	  if ( show_spectrum )
	    {

	      
	      std::map<int,std::vector<double> >::const_iterator ii = track_freq.begin();
	      
	      while ( ii != track_freq.end() )
		{
		  
		  if ( freqs[ ii->first ] > max_power )  { ++ii; continue; } 
		  
		  writer.level( freqs[ ii->first ] , globals::freq_strat );		  
		  
		  if ( has_cycles )
		    dynam_report_with_log( ii->second , epochs , &cycle );
		  else
		    dynam_report_with_log( ii->second , epochs );
		  
		  ++ii;
		  
		}
	      
	      writer.unlevel( globals::freq_strat );
	      
	    }
	  
	}
    
      
      //
      // Output a new signal?
      //

      if ( new_sigs )
	{
	 
	  
	  // note cannot have gaps in recording if adding a new channel
	  // (can be EDF+ but all epochs must be included)
	  // with 4/2 overlap, we'll need to pad w/ a zero at the end,
	  // so force the vector size to be the expected given the same rate
	  // i.e. increment step

	  const double psr = 1.0 / edf.timeline.epoch_inc() ; // e.g. 0.5 for 2s steps
	  
	  // expected
	  const int expected = edf.header.nr * edf.header.record_duration * psr;
	  
	  // we have this many epochs
	  const int obs = edf.timeline.num_epochs();

	  const int diff = expected - obs;
	  
	  // if diff is even, pad at beginning and end
	  // if odd, pad (N-1)/2+1 at start, (N-1)/2 at end
	  int pad1 = 0 , pad2 = 0;
	  if ( diff > 0 )
	    {
	      if ( diff % 2 == 0 )
		pad1 = pad2 = diff / 2;
	      else
		{
		  pad1 = ( diff-1 ) / 2 + 1 ;
		  pad2 = pad1 - 1 ;
		}	      
	    }
	  
	  
	  // but given SR of 1/(epoch-inc), we will typically have fewer (i.e. at the end of the recording,
	  // we might miss 1 sample with epoch len = 4s, epoch step = 2s
	  // just zero-pad these final points
	  
	  //
	  // band power
	  //

	  if ( bands )
	    {
	      
	      // are we storing relative power? 
	      std::vector<double> denom( obs , 0 );  // total denom per interval

	      if ( new_sigs_relpow )
		{
		  
		  std::vector<frequency_band_t>::const_iterator bi = bandaid.bands.begin();
		  while ( bi != bandaid.bands.end() )
		    {	   		  
		      
		      // skip?
		      if ( new_sigs_skip_bands.find( globals::band( *bi ) ) != new_sigs_skip_bands.end() )
			{
			  ++bi;
			  continue;
			}
		      
		      std::vector<double> vec = bandaid.track_band[ *bi ];
		           	      
		      // accumulate: note, this ignores if bands are overlapping or incomplete
		      // no zero-padding here either - that is done last (jsut before insert into edf)
		      for (int i=0; i<vec.size(); i++)
			denom[ i ] += vec[i] ; 
		      
		      ++bi;
		    }

		}
	      
	      //
	      // construct and add channels
	      //

	      int bidx=1;
	      
	      std::vector<frequency_band_t>::const_iterator bi = bandaid.bands.begin();
	      while ( bi != bandaid.bands.end() )
		{	   
		  
		  std::vector<double> vec = bandaid.track_band[ *bi ];
		  
		  // skip?
		  if ( new_sigs_skip_bands.find( globals::band( *bi ) ) != new_sigs_skip_bands.end() )
		    {
		      ++bi;
		      continue;
		    }
		  
		  if ( s == 0 )
		    {
		      freq_range_t fr = globals::freq_band[ *bi ];
		      logger << "   - B" << bidx << " --> " << globals::band( *bi )
			     << " ( " << fr.first << " - " << fr.second << " Hz )"
			     << "\n"; 
		    }
		  
		  
		  // log scale?
		  if ( dB )
		    for (int i=0; i<vec.size(); i++)
		      vec[i] = vec[i] > 0 ? 10 * log10( vec[i] ) : -999;

		
		  // use relpow?
		  if ( new_sigs_relpow )
		    {
		      for (int i=0; i<vec.size(); i++)
			if ( denom[ i ] > 0 ) 
			  vec[i] /= denom[ i ];		      
		    }

		  
		  // zero-pad as needed
		  vec.insert( vec.begin(), pad1 , 0 );
		  for (int j=0;j<pad2;j++) vec.push_back( 0 );
		  
		  // label
		  const std::string slab = new_sig_prefix + signals.label(s) + "_B" + Helper::int2str( bidx );
		  
		  // add signal (will always be 1 Hz)
		  edf.add_signal( slab , 1 , vec );
		  
		  ++bidx;
		  ++bi;
		}
	      
	      if ( bidx > 1 )
		logger << "  for " << signals.label(s) << " added " << bidx-1 << " band signals"
		       << ", padding " << obs << " samples to " << expected << " (adding " << pad1 << " leading and " << pad2 << " trailing 0s)\n";
												
	      
	    }

		  
	  //
	  // spectrum-based channels?
	  //

	  if ( new_spec_sigs )
	    {
	      
	      // add channels
	      std::map<int,std::vector<double> >::const_iterator ii = track_freq.begin();
	      
	      if ( track_freq.size() > 0 ) 
		logger << "  and adding " << track_freq.size() << " signals (" << freqs[ ii->first ] << "Hz - " << max_power << "Hz)"
		       << ", padding " << obs << " samples to " << expected << " (adding " << pad1 << " leading and " << pad2 << " trailing 0s)\n";
	      
	      while ( ii != track_freq.end() )
		{
		  
		  if ( freqs[ ii->first ] > max_power )  { ++ii; continue; }
		  
		  if ( s == 0 )
		    logger << "   - F" << ii->first << " --> " << freqs[ ii->first ] << "Hz\n";

		  std::vector<double> vec = ii->second;
		  
		  if ( dB ) 
		    for (int i=0; i<vec.size(); i++)
		      vec[i] = vec[i] > 0 ? 10 * log10( vec[i] ) : -999;

		  // zero-pad as needed
		  vec.insert( vec.begin(), pad1 , 0 );
		  for (int j=0;j<pad2;j++) vec.push_back( 0 );
		  
		  // label
		  const std::string slab = new_sig_prefix + signals.label(s) + "_F" + Helper::int2str( ii->first ) ;
		  
		  // add signal (will always be 1 Hz
		  edf.add_signal( slab , 1 , vec );
		  
		  ++ii;
		  
		}
	    }
	  
	}
      

      //
      // Multi-scale entropy
      //
      

      if ( calc_mse )
	{
	  
	  int mse_lwr_scale = 1;
	  int mse_upr_scale = 10;
	  int mse_inc_scale = 2;
	  int mse_m = 2;
	  double mse_r = 0.15;

	  mse_t mse( mse_lwr_scale , mse_upr_scale , mse_inc_scale ,
		     mse_m , mse_r );

	  
	  std::map<frequency_band_t,std::vector<double> >::const_iterator ii = bandaid.track_band.begin();

	  while ( ii != bandaid.track_band.end() )
	    {
	      
	      writer.level( globals::band( ii->first ) , globals::band_strat );
	      
	      std::map<int,double> mses = mse.calc( ii->second );
	      
	      std::map<int,double>::const_iterator jj = mses.begin();
	      while ( jj != mses.end() )
		{
		  writer.level( jj->first , "SCALE" );
		  writer.value( "MSE" , jj->second );
		  ++jj;
		}
	      writer.unlevel( "SCALE" );
	      ++ii;
	    }
	  writer.unlevel( globals::band_strat );
	}
     
      //
      // Band-power ratios 
      //
      
      if ( bands && calc_ratio )
	{
	  // ratio=ALPHA/BETA,THETA/DELTA,...

	  std::vector<std::string> r = Helper::parse( Helper::toupper( ratios ) , ',' );
	  
	  bool done_any = false; 
	  
	  std::vector<std::string>::const_iterator rr = r.begin();
	  while ( rr != r.end() )
	    {
	      // A/B
	      std::vector<std::string> tok = Helper::parse( *rr , '/' );
	      if ( tok.size() != 2 ) Helper::halt( "bad format for PSD ratio: " + *rr );
	      	      
	      frequency_band_t b1 = globals::band( tok[0] );
	      frequency_band_t b2 = globals::band( tok[1] );

	      // calculate both mean of epoch-ratios,
	      // as well as ratio of means of epoch-power
	      // optionally (ratio1) add +1 to denom, so it is always defined

	      if ( b1 != UNKNOWN_BAND && b2 != UNKNOWN_BAND )
		{
		  const std::vector<double> & p1 = bandaid.track_band[ b1 ];
		  const std::vector<double> & p2 = bandaid.track_band[ b2 ];

		  if ( p1.size() != p2.size() ) Helper::halt( "internal error" );

		  std::vector<double> rat;
		  double pw1 = 0 , pw2 = 0;
		  for (int i=0;i<p1.size(); i++)
		    {
		      rat.push_back( p1[i] / ( ratio_plus1 + p2[i] ) );
		      pw1 += p1[i];
		      pw2 += p2[i];
		    }

		  if ( rat.size() > 0 )
		    {
		      const double rmean = MiscMath::mean( rat );
		      const double rmedian = MiscMath::median( rat );
		      writer.level( tok[0] , "B1" );
		      writer.level( tok[1] , "B2" );
		      writer.value( "RATIO" , rmean );
		      writer.value( "RATIO_MN" , pw1 / ( ratio_plus1 + pw2 ) );
		      writer.value( "RATIO_MD" , rmedian );		  
		      done_any = true;
		    }
		}
	      ++rr;
	    }
	  
	  if ( done_any )
	    {
	      writer.unlevel( "B2" );
	      writer.unlevel( "B1" );
	    }
	  
	}

      
      //
      // Band-power kurtosis  (redundant now...) 
      //

      if ( bands && calc_kurt )
	{	  
	  std::map<frequency_band_t,std::vector<double> >::const_iterator ii = bandaid.track_band.begin();	  
	  while ( ii != bandaid.track_band.end() )
	    {	      
	      writer.level( globals::band( ii->first ) , globals::band_strat );	      
	      // get kurtosis of dB-scaled values
	      //  track_band is always raw values
	      double k = MiscMath::kurtosis( MiscMath::dB( ii->second ) ) + kurt_adj ;
	      writer.value( "KURT" , k );
	      ++ii;
	    }
	  writer.unlevel( globals::band_strat );
	}


      
      //
      // next signal
      //
    } 


  writer.unlevel( globals::signal_strat );	   
      
  
  // ignore return annot_t * 
  return NULL;

}



