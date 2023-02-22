
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

#include "spectral.h"
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


extern writer_t writer;
extern logger_t logger;



annot_t * spectral_power( edf_t & edf , 
			  const std::string & signal_label , 
			  param_t & param )
{

  
  // Report full spectrum as well as band power
  const bool show_spectrum = param.has( "spectrum" ) || param.has("epoch-spectrum" );

  // Report dB scale ?
  const bool dB = param.has( "dB" );
  
  // Mean center data first? 
  const bool mean_centre_epoch = param.has( "center" ) || param.has( "centre" );
  
  // Spectrum bin width (0 means no binning, default)

  //double bin_width = param.has( "bin" ) ? param.requires_dbl( "bin" ) : 0;
  const double bin_fac = param.has( "fac" ) ? param.requires_int( "fac" ) : 1;

  // Band power per-epoch
  const bool show_epoch = param.has( "epoch" ) || param.has("epoch-spectrum" );
  
  // report variance of PSD across epochs (SD)
  const bool aggregate_psd_sd       = param.has( "sd" );
  const double aggregate_psd_th     = param.has( "th" ) ? param.requires_dbl( "th" ) : 0 ; 
  const bool aggregate_psd_med      = param.has( "median" );

  // Characterize dynamics: of all epoch-level stats created, get the H1, H2, H3, linear and exponential trend
  // Hjorth stats just consider all points concatenated
  // trend lines are based on observed E numbers 

  const bool calc_dynamics = param.has( "dynamics" );
  
  // Verbose output: full spectrum per epoch

  const bool show_epoch_spectrum = param.has( "epoch-spectrum" );

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

  double fft_segment_overlap = param.has( "segment-overlap" ) 
     ? param.requires_dbl( "segment-overlap" ) : 2 ;
       
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
  // Return SD of segments (actually CV)
  //

  const bool calc_seg_sd = param.has( "segment-sd" );
  
  //
  // Use nextpow2 for NFFT
  //

  const bool use_nextpow2 = param.has( "pow2" );


  //
  // User defined 'TOTAL' ? (allow change o the fly)
  //

  globals::freq_band[ DENOM ] = globals::freq_band[ TOTAL ];

  if ( param.has( "total" ) )
    {
      std::vector<std::string> f = Helper::parse( param.value( "total" ) , ",-" );
      if ( f.size() != 2 ) Helper::halt( "expecting band=lower,upper" );
      double f0, f1;
      if ( ! Helper::str2dbl( f[0] , &f0 ) ) Helper::halt( "expecting numeric for total power range" );
      if ( ! Helper::str2dbl( f[1] , &f1 ) ) Helper::halt( "expecting numeric for total power range" );
      if ( f0 >= f1 ) Helper::halt( "expecting band=lower,upper" );
      if ( f0 < 0 || f1 < 0 ) Helper::halt( "negative frequencies specified" );

      // update
      logger << "  setting total power (denominator for RELPSD) to " << f0 << " to " << f1 << " Hz\n"; 
      globals::freq_band[ DENOM ] = freq_range_t( f0 , f1 ) ;
    }
  
  //
  // Define standard band summaries
  //
  
  std::vector<frequency_band_t> bands;
  bands.push_back( SLOW );
  bands.push_back( DELTA );
  bands.push_back( THETA );
  bands.push_back( ALPHA );
  bands.push_back( SIGMA );
  if ( 0 )
    {
      bands.push_back( LOW_SIGMA );
      bands.push_back( HIGH_SIGMA );
    }
  bands.push_back( BETA );
  bands.push_back( GAMMA );
  bands.push_back( DENOM );

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
  // Check segment lengths
  //

  if ( edf.timeline.epoch_length() <= ( fft_segment_size + fft_segment_overlap ) )
    {
      fft_segment_overlap = 0;
      fft_segment_size = edf.timeline.epoch_length();
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

      if ( Fs[s] < 50 ) continue;

      
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

      // band and F results
      std::map<frequency_band_t,std::vector<double> > track_band;
      std::map<int,std::vector<double> > track_freq;
      std::map<int,std::vector<double> > track_freq_logged;

      
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
	  

	  ++total_epochs;

	  
	  interval_t interval = edf.timeline.epoch( epoch );
	  
	  // stratify output by epoch?
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
	   
	   
// 	   logger << "total_points = " << total_points << "\n";
// 	   logger << "nooverlap_segments = " << noverlap_segments << "\n";
// 	   logger << "noverlap_points = " << noverlap_points << "\n";
// 	   logger << "segment_points = " << segment_points << "\n";
	   
// 	   std::cout << "about to fly...\n";
// 	   std::cout << "Fs = " << Fs[s] << "\n";

	   PWELCH pwelch( *d , 
			  Fs[s] , 
			  segment_sec , 
			  noverlap_segments , 
			  window_function , 
			  use_seg_median,
			  calc_seg_sd,
			  average_adj ,
			  use_nextpow2 );
	   

	   double this_slowwave   = pwelch.psdsum( SLOW )  ;      /// globals::band_width( SLOW );
	   double this_delta      = pwelch.psdsum( DELTA ) ;      /// globals::band_width( DELTA );
	   double this_theta      = pwelch.psdsum( THETA ) ;      /// globals::band_width( THETA );
	   double this_alpha      = pwelch.psdsum( ALPHA ) ;      /// globals::band_width( ALPHA );
	   double this_sigma      = pwelch.psdsum( SIGMA ) ;      /// globals::band_width( SIGMA );
	   double this_low_sigma  = pwelch.psdsum( LOW_SIGMA ) ;  /// globals::band_width( LOW_SIGMA );
	   double this_high_sigma = pwelch.psdsum( HIGH_SIGMA ) ; /// globals::band_width( HIGH_SIGMA );
	   double this_beta       = pwelch.psdsum( BETA )  ;      /// globals::band_width( BETA );
	   double this_gamma      = pwelch.psdsum( GAMMA ) ;      /// globals::band_width( GAMMA );]
	   double this_total      = pwelch.psdsum( DENOM ) ;      /// globals::band_width( DENOM );
	   
	   //
	   // track epoch-level band-power statistics
	   //
	   
	   track_band[ SLOW  ].push_back( this_slowwave );
	   track_band[ DELTA ].push_back( this_delta );
	   track_band[ THETA ].push_back( this_theta );
	   track_band[ ALPHA ].push_back( this_alpha );
	   track_band[ SIGMA ].push_back( this_sigma );
	   track_band[ BETA  ].push_back( this_beta );
	   track_band[ GAMMA ].push_back( this_gamma );
	   track_band[ DENOM ].push_back( this_total );

	   if ( 0 )
	     {
	       track_band[ LOW_SIGMA ].push_back( this_low_sigma );
	       track_band[ HIGH_SIGMA ].push_back( this_high_sigma );
	     }
	   
	   //
	   // track epoch numbers (for dynam_t)
	   //
	   
	   epochs.push_back( epoch );

	   //
	   // Epoch-level output
	   //

	   if ( show_epoch || ( cache_epochs && cache_bands ) )
	     {
	       
	       double this_total =  this_slowwave
		 + this_delta
		 + this_theta
		 + this_alpha
		 + this_sigma  
		 + this_beta
		 + this_gamma;
	       
	       if ( this_total > 0 )
		 {
		   writer.level( globals::band( SLOW ) , globals::band_strat );
		   if ( show_epoch && ! suppress_output ) {
		     writer.value( "PSD" , dB ? 10*log10( this_slowwave ) : this_slowwave  );
		     writer.value( "RELPSD" , this_slowwave / this_total );
		   }
		   if ( cache_data && cache_epochs && cache_bands )
		     cache->add( ckey_t( "PSD" , writer.faclvl() ) , dB ? 10*log10( this_slowwave ) : this_slowwave );
		   
		   writer.level( globals::band( DELTA ) , globals::band_strat );
		   if ( show_epoch && ! suppress_output ) {
		     writer.value( "PSD" , dB ? 10*log10( this_delta ) : this_delta );
		     writer.value( "RELPSD" , this_delta / this_total );
		   }
		   if ( cache_data && cache_epochs && cache_bands )
		     cache->add( ckey_t( "PSD" , writer.faclvl() ) , dB ? 10*log10( this_delta ) : this_delta );
		   
		   writer.level( globals::band( THETA ) , globals::band_strat );
		   if ( show_epoch && ! suppress_output ) {
		     writer.value( "PSD" , dB ? 10*log10( this_theta ) : this_theta  );
		     writer.value( "RELPSD" , this_theta / this_total );
		   }
		   if ( cache_data && cache_epochs && cache_bands )
		     cache->add( ckey_t( "PSD" , writer.faclvl() ) , dB ? 10*log10( this_theta ) : this_theta );
		   
		   writer.level( globals::band( ALPHA ) , globals::band_strat );
		   if ( show_epoch && ! suppress_output ) {
		     writer.value( "PSD" , dB ? 10*log10( this_alpha ) : this_alpha );
		     writer.value( "RELPSD" , this_alpha / this_total );
		   }
		   if ( cache_data && cache_epochs && cache_bands )
		     cache->add( ckey_t( "PSD" , writer.faclvl() ) , dB ? 10*log10( this_alpha ) : this_alpha );
		   
		   writer.level( globals::band( SIGMA ) , globals::band_strat );
		   if ( show_epoch && ! suppress_output ) {
		     writer.value( "PSD" , dB ? 10*log10( this_sigma ) : this_sigma );
		     writer.value( "RELPSD" , this_sigma / this_total );
		   }
		   if ( cache_data && cache_epochs && cache_bands )
		     cache->add( ckey_t( "PSD" , writer.faclvl() ) , dB ? 10*log10( this_sigma ) : this_sigma );

		   writer.level( globals::band( BETA ) , globals::band_strat );
		   if ( show_epoch && ! suppress_output ) {
		     writer.value( "PSD" , dB ? 10*log10( this_beta ) : this_beta  );
		     writer.value( "RELPSD" , this_beta / this_total );
		   }
		   if ( cache_data && cache_epochs && cache_bands )
		     cache->add( ckey_t( "PSD" , writer.faclvl() ) , dB ? 10*log10( this_beta ) : this_beta );

		   writer.level( globals::band( GAMMA ) , globals::band_strat );
		   if ( show_epoch && ! suppress_output ) {		     
		     writer.value( "PSD" , dB ? 10*log10( this_gamma ) : this_gamma );
		     writer.value( "RELPSD" , this_gamma / this_total );
		   }
		   if ( cache_data && cache_epochs && cache_bands )
		     cache->add( ckey_t( "PSD" , writer.faclvl() ) , dB ? 10*log10( this_gamma ) : this_gamma );
		   
		   writer.level( globals::band( DENOM ) , globals::band_strat );
		   if ( show_epoch && ! suppress_output ) {
		     writer.value( "PSD" , dB ? 10*log10( this_total ) : this_total );				   
		   }
		   if ( cache_data && cache_epochs && cache_bands )
		     cache->add( ckey_t( "PSD" , writer.faclvl() ) , dB ? 10*log10( this_total ) : this_total );
		   
		   writer.unlevel( globals::band_strat );
		   
		 }
	       else if ( cache_data && cache_epochs && cache_bands && ! dB )
		 {
		   // need to enter 0 in this case for cache
		   //  nb. only doing this in non-dB mode (i.e. for ASYMM(
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
	     

	       
	   // std::cout << "freqs.size() = " << freqs[s].size() << "\n";
	   // std::cout << "pwelch.size() = " << pwelch.psd.size() << "\n";
	   
	   if ( freqs.size() == pwelch.psd.size() )
	     {
	       
	       //
	       // accumulate for entire night means; store as dB and raw
	       //
	       
	       if ( show_spectrum || spectral_slope )
		 for (int f=0;f<pwelch.psd.size();f++)
		   {
		     track_freq[ f ].push_back( pwelch.psd[f] );
		     // only track if power > 0
		     if ( pwelch.psd[f] > 0 ) 
		       track_freq_logged[ f ].push_back( 10*log10( pwelch.psd[f] ) );
		   }

	       //
	       // epoch-level output?
	       //
	       
	       if ( show_epoch_spectrum || ( cache_epochs && cache_spectrum ) )
		 {		 
		   
		   // using bin_t 	      
		   bin_t bin( min_power , max_power , bin_fac );
		   bin.bin( freqs , pwelch.psd );
		   
		   bin_t binsd( min_power , max_power , bin_fac );
		   if ( calc_seg_sd )
		     binsd.bin( freqs, pwelch.psdsd );
		   
		   std::vector<double> f0;
		   
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
	  
	  int ne_valid = dB ? track_freq_logged[0].size() : track_freq[0].size();
	  int ne_min = ne_valid;

	  for (int f=0;f<n;f++) 
	    {
	      // wanting to get stats of dB or raw?
	      const std::vector<double> & yy = dB ? track_freq_logged[f] : track_freq[f] ;
	      
	      // any outlier removal of epochs?  	      
	      std::vector<double> xx = aggregate_psd_th > 0 && ne_valid > 2 ? 
		MiscMath::outliers( &yy , aggregate_psd_th ) : yy ; 
	      
	      // track min size
	      if ( xx.size() < ne_min ) 
		ne_min = xx.size();
	      
	      means.push_back( MiscMath::mean( xx ) );
	      
	      if ( aggregate_psd_sd && xx.size() > 2 )
		sds.push_back( MiscMath::sdev( xx ) );
	      
	      if ( aggregate_psd_med && xx.size() > 2 )
		medians.push_back(  MiscMath::median( xx ) );
	    }
	  
	  bin_t bin( min_power , max_power , bin_fac );	  
	  bin.bin( freqs , means );

	  bin_t bin_med( min_power , max_power , bin_fac );	  
	  if ( aggregate_psd_med && ne_min > 2 )
	    bin_med.bin( freqs , medians );

	  bin_t bin_sds( min_power , max_power , bin_fac );	  
	  if ( aggregate_psd_sd && ne_min > 2 ) 
	    bin_sds.bin( freqs , sds );

	  
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
		      
		      if ( aggregate_psd_med && ne_min > 2 )
			writer.value( "PSD_MD" , bin_med.bspec[i]  );
		      
		      if ( aggregate_psd_sd && ne_min > 2 )
			writer.value( "PSD_SD" , bin_sds.bspec[i]  );
		      
		      if ( bin.nominal[i] != "" )
			writer.value( "INT" , bin.nominal[i] );
		    }
		}

	      //
	      // Cache summary spectra?
	      //

	      if ( cache_data && cache_spectrum )
		cache->add( ckey_t( "PSD" , writer.faclvl() ) , x );
	      
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

      double mean_total_power = MiscMath::mean( track_band[ DENOM ] );
      

      //
      // by band 
      //      
      
      std::vector<frequency_band_t>::const_iterator bi = bands.begin();
      while ( bi != bands.end() )
	{	   

	  if ( okay ) 
	    {
	      double p = MiscMath::mean( track_band[ *bi ] );

	      writer.level( globals::band( *bi ) , globals::band_strat );

	      if ( ! suppress_output ) {		
		writer.value( "PSD" , dB ? 10*log10(p) : p  );
		writer.value( "RELPSD" , p / mean_total_power );
	      }
	      
	      // Cache summary bands?
	      if ( cache_data && cache_bands )
		cache->add( ckey_t( "PSD" , writer.faclvl() ) , dB ? 10*log10(p) : p );
	      
	    }
	  
 	  ++bi;
	}
      
      writer.unlevel( globals::band_strat );


      
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

	  std::map<frequency_band_t,std::vector<double> >::const_iterator ii = track_band.begin();
	  
	  while ( ii != track_band.end() )
	    {	      
	      writer.level( globals::band( ii->first ) , globals::band_strat );
	      
	      if ( has_cycles )
		dynam_report_with_log( ii->second , epochs , &cycle );
	      else
		dynam_report_with_log( ii->second , epochs );
	      
	      ++ii;
	    }
	  
	  writer.unlevel( globals::band_strat ); 


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

	  
	  std::map<frequency_band_t,std::vector<double> >::const_iterator ii = track_band.begin();

	  while ( ii != track_band.end() )
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
      // next signal
      //
    } 


  writer.unlevel( globals::signal_strat );	   
      
  
  // ignore return annot_t * 
  return NULL;

}



