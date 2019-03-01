
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
  
  // Spectrum bin width (0.5 Hz by default)

  double bin_width = param.has( "bin" ) ? param.requires_dbl( "bin" ) : 0.5;

  // Band power per-epoch

  bool show_epoch = param.has( "epoch" ) || param.has("epoch-spectrum" );

  // Verbose output: full spectrum per epoch

  bool show_epoch_spectrum = param.has( "epoch-spectrum" );

  // truncate spectra
  double max_power = param.has( "max" ) ? param.requires_dbl( "max" ) : 20 ;

  // Calculate MSE

  bool calc_mse = param.has( "mse" ); 


  //
  // Alter PWELCH sliding window parameters
  //

  double fft_segment_size = param.has( "segment-sec" ) 
    ? param.requires_dbl( "segment-sec" ) : 4 ;

  double fft_segment_overlap = param.has( "segment-overlap" ) 
    ? param.requires_dbl( "segment-overlap" ) : 2 ;


  //
  // Option to average adjacent points in the power spectra (default = Y)
  //
  
  bool average_adj = param.has( "no-averaging" ) ? false : true;
  
  //
  // Window function
  //
  
  window_function_t window_function = WINDOW_TUKEY50;	   
  if      ( param.has( "no-window" ) ) window_function = WINDOW_NONE;
  else if ( param.has( "hann" ) ) window_function = WINDOW_HANN;
  else if ( param.has( "hamming" ) ) window_function = WINDOW_HAMMING;
  else if ( param.has( "tukey50" ) ) window_function = WINDOW_TUKEY50;

  
  //
  // Spectrum summaries
  //
  
  bool show_user_ranges = param.has( "ranges" ); // lwr, upr, inc

  bool show_user_ranges_per_epoch = param.has( "epoch-ranges" ); // boolean

  if ( show_user_ranges_per_epoch && ! show_user_ranges ) 
    Helper::halt( "need to specify ranges=lwr,upr,inc with epoch-ranges" );

  std::map<freq_range_t,double> user_ranges;
  if ( show_user_ranges )
    {
      std::vector<std::string> tok = param.strvector( "ranges" );
      if ( tok.size() != 3 ) Helper::halt( "expecting 3 values for ranges=lwr,upr,inc" );
      double lwr = 0 , upr = 0 , inc = 0;
      if ( ! Helper::str2dbl( tok[0] , &lwr ) ) Helper::halt( "bad ranges" );
      if ( ! Helper::str2dbl( tok[1] , &upr ) ) Helper::halt( "bad ranges" );
      if ( ! Helper::str2dbl( tok[2] , &inc ) ) Helper::halt( "bad ranges" );
      if ( lwr >= upr ) Helper::halt( "bad ranges" );      
      //      std::cout << lwr << " " << upr << " " << inc << "\n";
      
      for (double f = lwr; f<= upr; f += inc ) 
	user_ranges[ freq_range_t( f , f+inc ) ] = 0;      

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
  bands.push_back( LOW_SIGMA );
  bands.push_back( HIGH_SIGMA );
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
  // get high, low and total power.
  //

  int total_epochs = 0;  

  std::map<int,std::vector<double> > freqs;

  std::map<int,std::map<frequency_band_t,double> > bandmean;
  std::map<int,std::map<freq_range_t,double> > rangemean;  
  std::map<int,std::map<int,double> > freqmean; 
  
  std::map<int,std::map<frequency_band_t,std::vector<double> > > track_mse;
  

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
  

  bool epoch_level_output = show_epoch || show_user_ranges_per_epoch || show_epoch_spectrum ;
  
    
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
       // Get each signal
       //
       
       for (int s = 0 ; s < ns; s++ )
	 {
	   
	   //
	   // only consider data tracks
	   //

	   if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
	   
	   //
	   // Stratify output by channel
	   //
	   
	   if ( epoch_level_output )
	     writer.level( signals.label(s) , globals::signal_strat );
	   
	   //
	   // Get data
	   //

	   slice_t slice( edf , signals(s) , interval );
	   
	   std::vector<double> * d = slice.nonconst_pdata();
	   
	   // mean centre just this epoch
	   if ( mean_centre_epoch ) 
	     MiscMath::centre( d );
	   
	   //
	   // pwelch() to obtain full PSD
	   //
	   
	   // Fixed parameters:: use 4-sec segments with 2-second
	   // overlaps 
	  	   
	   const double overlap_sec = fft_segment_overlap;
	   const double segment_sec  = fft_segment_size;
	   
	   const int total_points = d->size();
	   const int segment_points = segment_sec * Fs[s];
	   const int noverlap_points  = overlap_sec * Fs[s];
	   
	   // implied number of segments
	   int noverlap_segments = floor( ( total_points - noverlap_points) 
					  / (double)( segment_points - noverlap_points ) );
	   
	   
// 	    logger << "total_points = " << total_points << "\n";
// 	    logger << "nooverlap_segments = " << noverlap_segments << "\n";
// 	    logger << "noverlap_points = " << noverlap_points << "\n";
// 	    logger << "segment_points = " << segment_points << "\n";
	   
	   PWELCH pwelch( *d , 
			  Fs[s] , 
			  segment_sec , 
			  noverlap_segments , 
			  window_function , 
			  average_adj );

	   double this_slowwave   = pwelch.psdsum( SLOW )  ;/// globals::band_width( SLOW );
	   double this_delta      = pwelch.psdsum( DELTA ) ;/// globals::band_width( DELTA );
	   double this_theta      = pwelch.psdsum( THETA ) ;/// globals::band_width( THETA );
	   double this_alpha      = pwelch.psdsum( ALPHA ) ;/// globals::band_width( ALPHA );
	   double this_sigma      = pwelch.psdsum( SIGMA ) ;/// globals::band_width( SIGMA );
	   double this_low_sigma  = pwelch.psdsum( LOW_SIGMA ) ;/// globals::band_width( LOW_SIGMA );
	   double this_high_sigma = pwelch.psdsum( HIGH_SIGMA ) ;/// globals::band_width( HIGH_SIGMA );
	   double this_beta       = pwelch.psdsum( BETA )  ;/// globals::band_width( BETA );
	   double this_gamma      = pwelch.psdsum( GAMMA ) ;/// globals::band_width( GAMMA );]
	   double this_total      = pwelch.psdsum( TOTAL ) ;/// globals::band_width( TOTAL );
	   
	   //
	   // Track average over epochs
	   //
	    
	   bandmean[ s ][ SLOW ]       +=  this_slowwave;
	   bandmean[ s ][ DELTA ]      +=  this_delta;
	   bandmean[ s ][ THETA ]      +=  this_theta;
	   bandmean[ s ][ ALPHA ]      +=  this_alpha;
	   bandmean[ s ][ LOW_SIGMA ]  += this_low_sigma;
	   bandmean[ s ][ HIGH_SIGMA ] += this_high_sigma;
	   bandmean[ s ][ SIGMA ]      +=  this_sigma;	   
	   bandmean[ s ][ BETA ]       +=  this_beta;
	   bandmean[ s ][ GAMMA ]      +=  this_gamma;
	   bandmean[ s ][ TOTAL ]      +=  this_total;

	   //
	   // MSE for bands
	   //

	   track_mse[ s ][ SLOW  ].push_back( this_slowwave );
	   track_mse[ s ][ DELTA ].push_back( this_delta );
	   track_mse[ s ][ THETA ].push_back( this_theta );
	   track_mse[ s ][ ALPHA ].push_back( this_alpha );
	   track_mse[ s ][ SIGMA ].push_back( this_sigma );
	   track_mse[ s ][ BETA  ].push_back( this_beta );
	   track_mse[ s ][ GAMMA ].push_back( this_gamma );

	 
	   //
	   // Epoch-level output
	   //

	   //
	   // User-defined ranges
	   //

	   if ( show_user_ranges )
	     {
	       
	       rangemean[ s ] = user_ranges;
	       pwelch.psdsum( &rangemean[s] );
	       
	       if ( show_user_ranges_per_epoch )
		 {
		   
		   const std::map<freq_range_t,double> & x = rangemean[ s ];
		   std::map<freq_range_t,double>::const_iterator ii = x.begin();
		   while ( ii != x.end() )
		     {
		       std::stringstream ss;
		       ss << ii->first.first << ".." << ii->first.second ;
		       
		       writer.level( ss.str() , "FRQRANGE" );

		       // track numeric range mid-point
		       writer.var( "FRQMID" , "Frequency range mid-point" );
		       writer.value( "FRQMID" , ( ii->first.second  + ii->first.first ) / 2.0  );
		       
		       // power value
		       writer.var( "PSD" , "Spectral band power" );
		       writer.value( "PSD" , ii->second );
						
		       ++ii;
		     }

		   writer.unlevel( "FRQRANGE" );

		   // display
		 }
	       
	     }
	   

	   //
	   // detailed, per-EPOCH outout
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
	       writer.value( "PSD" , this_slowwave );
	       writer.value( "RELPSD" , this_slowwave / this_total );
	       
	       writer.level( globals::band( DELTA ) , globals::band_strat );
	       writer.value( "PSD" , this_delta );
	       writer.value( "RELPSD" , this_delta / this_total );

	       writer.level( globals::band( THETA ) , globals::band_strat );
	       writer.value( "PSD" , this_theta );
	       writer.value( "RELPSD" , this_theta / this_total );

	       writer.level( globals::band( ALPHA ) , globals::band_strat );
	       writer.value( "PSD" , this_alpha );
	       writer.value( "RELPSD" , this_alpha / this_total );

	       writer.level( globals::band( SIGMA ) , globals::band_strat );
	       writer.value( "PSD" , this_sigma );
	       writer.value( "RELPSD" , this_sigma / this_total );

	       writer.level( globals::band( LOW_SIGMA ) , globals::band_strat );
	       writer.value( "PSD" , this_low_sigma );
	       writer.value( "RELPSD" , this_low_sigma / this_total );

	       writer.level( globals::band( HIGH_SIGMA ) , globals::band_strat );
	       writer.value( "PSD" , this_high_sigma );
	       writer.value( "RELPSD" , this_high_sigma / this_total );

	       writer.level( globals::band( BETA ) , globals::band_strat );
	       writer.value( "PSD" , this_beta );
	       writer.value( "RELPSD" , this_beta / this_total );

	       writer.level( globals::band( GAMMA ) , globals::band_strat );
	       writer.value( "PSD" , this_gamma );
	       writer.value( "RELPSD" , this_gamma / this_total );

	       writer.level( globals::band( TOTAL ) , globals::band_strat );
	       writer.value( "PSD" , this_total );
	       
	       writer.unlevel( globals::band_strat );
	       
	     }
	 
	   
       //
       // track over entire spectrum
       //
	   
       if( freqs.find( s ) == freqs.end() )
	 {
	   freqs[s] = pwelch.freq;
	 }
       
       // std::cout << "freqs.size() = " << freqs[s].size() << "\n";
       // std::cout << "pwelch.size() = " << pwelch.psd.size() << "\n";

       if ( freqs[s].size() == pwelch.psd.size() )
	 {

	   // accumulate for entire night means
	   for (int f=0;f<pwelch.psd.size();f++)
	     freqmean[ s ][ f ] += pwelch.psd[f];

	   // epoch-level output?
	   
	   if ( show_epoch_spectrum )
	     {		 

	      // using bin_t 	      
	      bin_t bin( bin_width , max_power , Fs[s] );
	      
	      bin.bin( freqs[s] , pwelch.psd );
	      
	      for ( int i = 0 ; i < bin.bfa.size() ; i++ ) 
		{
		  //writer.level( Helper::dbl2str( bin.bfa[i] ) + "-" + Helper::dbl2str( bin.bfb[i] ) ,  globals::freq_strat  );
		  writer.level( ( bin.bfa[i] + bin.bfb[i] ) / 2.0 , globals::freq_strat );
		  writer.value( "PSD" , bin.bspec[i] );
		}
	       writer.unlevel( globals::freq_strat );
	     }
	   
	 }
       else
	 logger << " *** warning:: skipped a segment: different NFFT/internal problem ... \n";

       
       if ( epoch_level_output )
	 writer.unlevel( globals::signal_strat );
       
       
	 } // next signal
       
       if ( epoch_level_output )
	 writer.unepoch();
       
     } // next EPOCH
    
   
  //
  // Output summary power curve
  //
  
  
  for (int s=0; s < ns ; s++ )
    {

      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
      
      const int n = freqs[s].size();      
      
      writer.level( signals.label(s) , globals::signal_strat );
      
      writer.var( "NE" , "Number of epochs" );
      writer.value( "NE" , total_epochs );
      
      if ( show_spectrum )
	{
	  
	  if ( total_epochs > 0 ) 
	    {	  

	      // get mean power across epochs
	      if ( freqmean[s].size() != freqs[s].size() ) Helper::halt( "internal error psd_t" );
	      std::vector<double> means;
	      for (int f=0;f<n;f++) means.push_back( freqmean[ s ][ f ] / (double)total_epochs );
	      
	      // using bin_t 	      
	      bin_t bin( bin_width , max_power , Fs[s] );

	      bin.bin( freqs[s] , means );

	      for ( int i = 0 ; i < bin.bfa.size() ; i++ ) 
		{
		  //writer.level( Helper::dbl2str( bin.bfa[i] ) + "-" + Helper::dbl2str( bin.bfb[i] ) ,  globals::freq_strat  );
		  writer.level( ( bin.bfa[i] + bin.bfb[i] ) / 2.0 , globals::freq_strat );
		  writer.value( "PSD" , bin.bspec[i] );
		}
	       writer.unlevel( globals::freq_strat );

// 	      for (int f=0;f<n;f++)
// 		{	      		  
// 		  double x = freqmean[ s ][ f ] / (double)total_epochs;	      
// 		  writer.level( freqs[s][f] , globals::freq_strat );		  
// 		  writer.value( "PSD" , x ); 		  
// 		}
// 	      writer.unlevel( globals::freq_strat );

	      
	    }
	  
	}

      
      bool okay = total_epochs > 0 ;
      
      // get total power
      double total_power = bandmean[ s ][ TOTAL ] / (double)total_epochs;

//       std::vector<frequency_band_t>::const_iterator bi = bands.begin();
//       while ( bi != bands.end() )
// 	{
// 	  if ( okay ) 
// 	    {
// 	      double p = 
// 	      total_power += p;
// 	    }
// 	  ++bi;
// 	}
      
      
      // by band 
      std::vector<frequency_band_t>::const_iterator bi = bands.begin();
      while ( bi != bands.end() )
	{	   

	  if ( okay ) 
	    {
	      double p = bandmean[ s ][ *bi ] / (double)total_epochs;

	      writer.level( globals::band( *bi ) , globals::band_strat );
	      writer.value( "PSD" , p );
	      writer.value( "RELPSD" , p / total_power );
	    }
	  
 	  ++bi;
	}
      
      writer.unlevel( globals::band_strat );


      // by user-defined ranges
      if ( show_user_ranges )
	{
	  const std::map<freq_range_t,double> & x = rangemean[ s ];
	  std::map<freq_range_t,double>::const_iterator ii = x.begin();
	  while ( ii != x.end() )
	    {
	      std::stringstream ss;
	      ss << ii->first.first << ".." << ii->first.second;

	      writer.level( ss.str() , "FRQRANGE" );
	      writer.value( "FRQMID" , ( ii->first.second  + ii->first.first ) / 2.0 ); 
	      writer.value( "PSD" , ii->second / (double)total_epochs );	      
	      ++ii;
	    }
	  writer.unlevel( "FRQRANGE" );
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

	  
	  const std::map<frequency_band_t,std::vector<double> > & x = track_mse[ s ];
	  std::map<frequency_band_t,std::vector<double> >::const_iterator ii = x.begin();
	  while ( ii != x.end() )
	    {
	      
	      writer.level( globals::band( ii->first ) , globals::band_strat );
	      
// 	      for (int ll=0;ll<10;ll++) std::cout << " " << (ii->second)[ll];
// 	      std::cout << "\n";

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
      
      writer.unlevel( globals::signal_strat );
      
    } // summary for next signal
  
  // redundant, ignore this for now...
  annot_t * a = edf.timeline.annotations.add( "Sigma power" );  
  return a;

}


