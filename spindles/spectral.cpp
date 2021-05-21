
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
  
  // Report dull spectrum as well as band power
  bool show_spectrum = param.has( "spectrum" ) || param.has("epoch-spectrum" );

  // Report dB scale ?
  bool dB = param.has( "dB" );
  
  // Mean center data first? 
  bool mean_centre_epoch = param.has( "center" ) || param.has( "centre" );
  
  // Spectrum bin width (0 means no binning, default)

  //double bin_width = param.has( "bin" ) ? param.requires_dbl( "bin" ) : 0;
  double bin_fac = param.has( "fac" ) ? param.requires_int( "fac" ) : 1;

  // Band power per-epoch

  bool show_epoch = param.has( "epoch" ) || param.has("epoch-spectrum" );

  // Characterize dynamics: of all epoch-level stats created, get the H1, H2, H3, linear and exponential trend
  // Hjorth stats just consider all points concatenated
  // trend lines are based on observed E numbers 

  bool calc_dynamics = param.has( "dynamics" );
  
  // Verbose output: full spectrum per epoch

  bool show_epoch_spectrum = param.has( "epoch-spectrum" );

  // peak diagnostics
  bool peak_diagnostics = param.has( "peaks" );
  int peak_median_filter_n = param.has( "peaks-window" ) ? param.requires_int( "peaks-window" ) : 11 ; 
  
  // truncate spectra
  double min_power = param.has( "min" ) ? param.requires_dbl( "min" ) : 0.5 ;
  double max_power = param.has( "max" ) ? param.requires_dbl( "max" ) : 25 ;

  // Calculate MSE
  bool calc_mse = param.has( "mse" ); 

  // cache PSD for other analyses (e.g. PSC)
  const bool cache_data = param.has( "cache-metrics" );
  const std::string cache_name = cache_data ? param.value( "cache-metrics" ) : "" ;


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
  
  bool average_adj = param.has( "average-adj" ) ;

  
  //
  // Window function
  //
  
  window_function_t window_function = WINDOW_TUKEY50;	   
  if      ( param.has( "no-window" ) ) window_function = WINDOW_NONE;
  else if ( param.has( "hann" ) ) window_function = WINDOW_HANN;
  else if ( param.has( "hamming" ) ) window_function = WINDOW_HAMMING;
  else if ( param.has( "tukey50" ) ) window_function = WINDOW_TUKEY50;

  //
  // Use nextpow2 for NFFT
  //

  bool use_nextpow2 = param.has( "pow2" );
  
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
  bands.push_back( TOTAL );

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
  

  bool epoch_level_output = show_epoch || show_epoch_spectrum ;
  
  
  
  //
  // Get each signal
  //
  

  logger << "  calculating PSD for " << ns << " signals\n"; 

  for (int s = 0 ; s < ns; s++ )
    {
      

      //
      // only consider data tracks
      //
      
      if ( edf.header.is_annotation_channel( signals(s) ) ) 
	continue;
      

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
			  average_adj ,
			  use_nextpow2 );
	   
	   //	   std::cout << "done\n";


	   double this_slowwave   = pwelch.psdsum( SLOW )  ;      /// globals::band_width( SLOW );
	   double this_delta      = pwelch.psdsum( DELTA ) ;      /// globals::band_width( DELTA );
	   double this_theta      = pwelch.psdsum( THETA ) ;      /// globals::band_width( THETA );
	   double this_alpha      = pwelch.psdsum( ALPHA ) ;      /// globals::band_width( ALPHA );
	   double this_sigma      = pwelch.psdsum( SIGMA ) ;      /// globals::band_width( SIGMA );
	   double this_low_sigma  = pwelch.psdsum( LOW_SIGMA ) ;  /// globals::band_width( LOW_SIGMA );
	   double this_high_sigma = pwelch.psdsum( HIGH_SIGMA ) ; /// globals::band_width( HIGH_SIGMA );
	   double this_beta       = pwelch.psdsum( BETA )  ;      /// globals::band_width( BETA );
	   double this_gamma      = pwelch.psdsum( GAMMA ) ;      /// globals::band_width( GAMMA );]
	   double this_total      = pwelch.psdsum( TOTAL ) ;      /// globals::band_width( TOTAL );
	   
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
	   track_band[ TOTAL ].push_back( this_total );

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

	   if ( show_epoch )
	     {
	       
	       double this_total =  this_slowwave
		 + this_delta
		 + this_theta
		 + this_alpha
		 + this_sigma  
		 + this_beta
		 + this_gamma;
	       
	       writer.level( globals::band( SLOW ) , globals::band_strat );
	       writer.value( "PSD" , dB ? 10*log10( this_slowwave ) : this_slowwave  );
	       writer.value( "RELPSD" , this_slowwave / this_total );
	       
	       writer.level( globals::band( DELTA ) , globals::band_strat );
	       writer.value( "PSD" , dB ? 10*log10( this_delta ) : this_delta );
	       writer.value( "RELPSD" , this_delta / this_total );

	       writer.level( globals::band( THETA ) , globals::band_strat );
	       writer.value( "PSD" , dB ? 10*log10( this_theta ) : this_theta  );
	       writer.value( "RELPSD" , this_theta / this_total );

	       writer.level( globals::band( ALPHA ) , globals::band_strat );
	       writer.value( "PSD" , dB ? 10*log10( this_alpha ) : this_alpha );
	       writer.value( "RELPSD" , this_alpha / this_total );

	       writer.level( globals::band( SIGMA ) , globals::band_strat );
	       writer.value( "PSD" , dB ? 10*log10( this_sigma ) : this_sigma );
	       writer.value( "RELPSD" , this_sigma / this_total );

	       if ( 0 )
		 {
		   writer.level( globals::band( LOW_SIGMA ) , globals::band_strat );
		   writer.value( "PSD" , dB ? 10*log10( this_low_sigma ) : this_low_sigma  );
		   writer.value( "RELPSD" , this_low_sigma / this_total );
		   
		   writer.level( globals::band( HIGH_SIGMA ) , globals::band_strat );
		   writer.value( "PSD" , dB ? 10*log10( this_high_sigma ) : this_high_sigma );
		   writer.value( "RELPSD" , this_high_sigma / this_total );
		 }
	       
	       writer.level( globals::band( BETA ) , globals::band_strat );
	       writer.value( "PSD" , dB ? 10*log10( this_beta ) : this_beta  );
	       writer.value( "RELPSD" , this_beta / this_total );

	       writer.level( globals::band( GAMMA ) , globals::band_strat );
	       writer.value( "PSD" , dB ? 10*log10( this_gamma ) : this_gamma );
	       writer.value( "RELPSD" , this_gamma / this_total );

	       writer.level( globals::band( TOTAL ) , globals::band_strat );
	       writer.value( "PSD" , dB ? 10*log10( this_total ) : this_total );
	       
	       writer.unlevel( globals::band_strat );
	       
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
	       
	       // accumulate for entire night means
	       if ( show_spectrum )
		 for (int f=0;f<pwelch.psd.size();f++)
		   track_freq[ f ].push_back( pwelch.psd[f] );
	       
	       // epoch-level output?
	       
	       if ( show_epoch_spectrum )
		 {		 
		   
		   // using bin_t 	      
		   bin_t bin( min_power , max_power , bin_fac );
		   
		   bin.bin( freqs , pwelch.psd );
		   
		   for ( int i = 0 ; i < bin.bfa.size() ; i++ )
		     {		     
		       writer.level( ( bin.bfa[i] + bin.bfb[i] ) / 2.0 , globals::freq_strat );
		       //writer.level( bin.bfa[i] , globals::freq_strat );
		       writer.value( "PSD" , dB? 10*log10( bin.bspec[i] ) : bin.bspec[i] );
		       if ( bin.nominal[i] != "" )
			 writer.value( "INT" , bin.nominal[i] );
		     }
		   writer.unlevel( globals::freq_strat );
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
      
      writer.var( "NE" , "Number of epochs" );
      
      writer.value( "NE" , total_epochs );

      
      if ( show_spectrum )
	{	  
	  
	  if ( okay )
	    {	  
	      
	      //
	      // cache spectrum? (e.g. for PSC)
	      //
	      
	      cache_t<double> * cache = NULL ; 
	      if ( cache_data )
		cache = edf.timeline.cache.find_num( cache_name );

	      //
	      // get mean power across epochs
	      //
	      
	      if ( track_freq.size() != freqs.size() ) 
		Helper::halt( "internal error psd_t" );
	      
	      std::vector<double> means;
	      for (int f=0;f<n;f++) 
		means.push_back( MiscMath::mean( track_freq[f] ) );

	      bin_t bin( min_power , max_power , bin_fac );

	      bin.bin( freqs , means );

	      std::vector<double> f0;
	      
	      for ( int i = 0 ; i < bin.bfa.size() ; i++ ) 
		{

		  writer.level( ( bin.bfa[i] + bin.bfb[i] ) / 2.0 , globals::freq_strat );
		  //writer.level( bin.bfa[i] , globals::freq_strat );
		  
		  f0.push_back( ( bin.bfa[i] + bin.bfb[i] ) / 2.0 );

		  double x = dB ? 10*log10( bin.bspec[i] ) : bin.bspec[i] ;
		  writer.value( "PSD" , x );

		  if ( cache )
		    cache->add( ckey_t( "PSD" , writer.faclvl() ) , x );
		  
		  if ( bin.nominal[i] != "" )
		    writer.value( "INT" , bin.nominal[i] );
		}
	      writer.unlevel( globals::freq_strat );


	      //
	      // Report metrics on the PSD
	      //

	      if ( peak_diagnostics )
		{
		  double m1, m2;
		  
		  // detrended / smoothed / difference
		  std::vector<double> shape1, shape2, shape3;
		  
		  std::vector<double> logged = bin.bspec;
		  for (int i=0; i<logged.size(); i++)
		    logged[i] = 10*log10( logged[i] );
		  
		  psd_shape_metrics( f0 ,
				     logged , 
				     peak_median_filter_n , 
				     &m1, &m2 ,
				     &shape1, &shape2, &shape3);
		  
		  writer.value( "PK" , m1 );
		  writer.value( "SPK" , m2 );
		  
		  for (int i=0; i<f0.size(); i++)
		    {		      
		      writer.level( f0[i] , globals::freq_strat );
		      writer.value( "DT" , shape1[i] );
		      writer.value( "SM" , shape2[i] );
		      writer.value( "DF" , shape3[i] );		      
		    }
		  writer.unlevel( globals::freq_strat );		  

		}
	      
	    }
	  
	}



      //
      // mean total power
      //

      double mean_total_power = MiscMath::mean( track_band[ TOTAL ] );
      

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
	      writer.value( "PSD" , dB ? 10*log10(p) : p  );
	      writer.value( "RELPSD" , p / mean_total_power );
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
	  
	  bool has_cycles = edf.timeline.epoch_annotation( "_NREMC_1" );
	  
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
	      
	      writer.var( "MSE" , "Multiscale entropy" );
	      
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


