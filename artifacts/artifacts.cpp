
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

#include "artifacts.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "annot/annot.h"

#include "helper/helper.h"
#include "helper/logger.h"
#include "eval.h"
#include "db/db.h"

#include "dsp/resample.h"
#include "dsp/mse.h"
#include "dsp/lzw.h"

#include "miscmath/miscmath.h"
#include "fftw/fftwrap.h"

#include <iostream>
#include <fstream>

extern writer_t writer;

extern logger_t logger;

annot_t * brunner_artifact_detection( edf_t & edf , 
				      const std::string & signal_label , 
				      const std::string & filename )
{
  

  Helper::halt("brunner artifact detection: not fully implemented" );
  
  bool write_file = filename != "";
  

  //
  // attach signal
  //
  
  signal_list_t signals = edf.header.signal_list( signal_label );  
  const int ns = signals.size();

  
  // filter channel
  
  //
  // Brunner et al. (1996) power spectra:
  //

  // FFT on 4-sec intervals, with Hamming window
  // 0.25Hz bins collapsed to
  //  0.5Hz bins between 0.25 and 20.0Hz
  //  1Hz bins between 20.25 and 32.0Hz
  // -->results in 52 bins per 4-sec epoch
  
  std::vector<double> blwr, bupr;
  double f = 0.0; 
  
  for (int i = 0 ; i < 52 ; i++ )
    {
      const double w = f < 20 ? 0.5 : 1.0;
      blwr.push_back( f );
      bupr.push_back( f + w );
      f += w;
    }
  
  
  //
  // EPOCHS: 4-second, non-overlapping epochs
  //

  int ne = edf.timeline.set_epoch( 4 , 4 ); 
  
  
  //
  // Store results for highest band
  //
  
  std::vector<double> y;

  //
  // Iterate over epcohs
  //
    
  while ( 1 ) 
    {
      
      int epoch = edf.timeline.next_epoch();      
      
      if ( epoch == -1 ) break;
      
      interval_t interval = edf.timeline.epoch( epoch );
      
      //
      // Get data 
      //
      
      // TODO; iterator over timepoints
      const int s = 0;

      slice_t slice( edf , signals(s) , interval );

      std::vector<double> * d = slice.nonconst_pdata();
      
      //
      // Apply a Hamming window
      //
      
      MiscMath::hamming_window( d );
      
      
      //
      // FFT
      //

      const int Fs = 125; // FIXED sampling rate for now
      const int n = d->size();  // number of input points

      
      FFT fft(n,Fs);      
      
      fft.apply( *d );
      
      int retn = fft.cutoff;
      

      //
      // get power spectra
      //

      std::vector<double> bands = fft.power_bands(blwr,bupr);
      
      //
      // Get average (?sum) of highest bands : 26.25 -- 32.0
      //

      y.push_back( bands[51] + bands[50] + bands[49] + bands[48] + bands[47] + bands[46] );

    }



  //
  // Find bad EPOCHS
  //

  const double th = 4;

  const int ny = y.size();
  std::vector<bool> reject;

  for (int i=0; i<ny;i++)
    {
      // 3 min window, 4-sec epochs --> 45 windows
      // if target window is center, then take 22 each side
      int lwr = i - 22;
      int upr = i + 22;
      if ( lwr < 0 ) lwr = 0;
      if ( upr > ny-1) upr = ny-1;
      std::vector<double> tmp;
      for (int j=lwr;j<=upr;j++) tmp.push_back(y[j]);
      double median = median_destroy( &tmp[0] , tmp.size() );
      reject.push_back(  y[i] > th * median );      
    }
  
  
  //
  // Optionally open .annot file for writing
  //

  if ( write_file )
    {

      std::ofstream F( filename.c_str() , std::ios::out );
      
      F << "NAME\tBrunner\n"
	<< "DESC\tArtifact detection\n"
	<< "TYPE\tBINARY\n"
	<< "EPOCH\t4\t4\n"
	<< "COLS\treject\n";
      
      for (int e=0;e<reject.size();e++)
	F << "E" << "\t" 
	  << e+1 << "\t"
	  << reject[e] << "\n";
      
      F.close();
    }
  

  //
  // Create an epoch-based 'annot_t' object to return
  //
  

  annot_t * a = edf.timeline.annotations.add( "Brunner" );

//   a->set_description( "Brunner et al (1996) automatic artifact detection" );
//   a->add_col( "Brunner" , ATYPE_BINARY );
  
//   for (int epoch=0;epoch<reject.size();epoch++)
//     {
//       bool_event_t event( "brunner" , reject[ epoch ] );
//       a->add( edf.timeline.epoch( epoch ) , &event );
//     }
  
  return a;
}



annot_t * buckelmuller_artifact_detection( edf_t & edf , 
					   param_t & param , 
					   const std::string & signal_label , 
					   const double delta_threshold , 
					   const double beta_threshold , 
					   const double delta_lwr , 
					   const double delta_upr ,
					   const double beta_lwr , 
					   const double beta_upr ,
					   const std::string & filename )
{


  //
  // parameters
  //
  
  bool set_mask = ! param.has( "no-mask" );
  bool verbose = param.has("verbose") || param.has( "epoch" );
  
  //
  // attach signal
  //
  
  signal_list_t signals = edf.header.signal_list( signal_label );  

  const int ns = signals.size();

  //
  // Sampling frequencies
  //

  std::vector<double> Fs = edf.header.sampling_freq( signals );


  //
  // Point to first epoch (assume 30 seconds, but could be different)
  //
  
  int ne = edf.timeline.first_epoch();

 
  //
  // Store per-epoch power
  //

  std::vector<std::vector<double> > delta(ns);
  std::vector<std::vector<double> > beta(ns);

    
  //
  // for each each epoch 
  //

  int cnt = 0;
  std::vector<int> track_epochs;
  
  while ( 1 ) 
    {
      
      int epoch = edf.timeline.next_epoch();      
       
      if ( epoch == -1 ) break;
      
      track_epochs.push_back( epoch );
      
      //
      // Get data 
      //

      interval_t interval = edf.timeline.epoch( epoch );

      for ( int s=0; s<ns; s++ )
	{

	  //
	  // only consider data tracks
	  //
	  
	  if ( edf.header.is_annotation_channel( signals(s) ) ) continue;	  
	  
	  slice_t slice( edf , signals(s) , interval );
	  
	  std::vector<double> * d = slice.nonconst_pdata();

	  const std::vector<uint64_t> * tp = slice.ptimepoints();
	  
	  //
	  // Mean-centre 30-second window
	  //

	  MiscMath::centre( d );
            
	  //
	  // Apply PWELCH to this epoch
	  //

	  // aim to get 10 windows of 4 seconds in 30sec epoch
	  
	  int noverlap_segments = 10;         
	  int segment_size_sec = 4;
	  
	  PWELCH pwelch( *d , Fs[s] , segment_size_sec , noverlap_segments );
      
	  // track power bands     
	  
	  delta[s].push_back( pwelch.psdsum( DELTA ) );
	  beta[s].push_back( pwelch.psdsum( beta_lwr  , beta_upr ) );
	  

	} // next signal

    } // next epoch



  
  //
  // Report for each signal
  //
  
  std::vector<std::vector<double> > delta_average(ns);
  std::vector<std::vector<double> > beta_average(ns);
  
  for (int s=0;s<ns;s++)
    {

      //
      // only consider data tracks
      //
      
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;

      //
      // Output stratifier
      //

      writer.level( signals.label(s) , globals::signal_strat );


      //
      // Make running averages
      //
      
      delta_average[s] = MiscMath::moving_average( delta[s] , 15 );
      beta_average[s]  = MiscMath::moving_average( beta[s]  , 15 );
      
      int total = 0, altered = 0;

      for (int e=0;e<ne;e++)
	{
	  
	  double dfac = delta[s][e] / delta_average[s][e];
	  double bfac = beta[s][e] / beta_average[s][e];
	  
	  bool dmask = dfac > delta_threshold ;
	  bool bmask = bfac > beta_threshold ;
	  bool mask  = dmask || bmask;
	  
	  if ( verbose ) 
	    {
	      
	      writer.epoch( edf.timeline.display_epoch( track_epochs[e] ) );

	      writer.var( "DELTA" , "Delta power" );
	      writer.var( "DELTA_AVG" , "Local average delta power" );
	      writer.var( "DELTA_FAC" , "Relative delta power factor" );
	      writer.var( "BETA" , "Beta power" );
	      writer.var( "BETA_AVG" , "Local average beta power" );
	      writer.var( "BETA_FAC" , "Relative beta power factor" );
	      writer.var( "DELTA_MASK" , "Masked based on delta power" );
	      writer.var( "BETA_MASK" , "Masked based on beta power" );
	      writer.var( "MASK" , "Masked" );		

	      writer.value( "DELTA" , delta[s][e] );
	      writer.value( "DELTA_AVG" , delta_average[s][e] );
	      writer.value( "DELTA_FAC" , dfac );

	      writer.value( "BETA" , beta[s][e] );
	      writer.value( "BETA_AVG" , beta_average[s][e] );
	      writer.value( "BETA_FAC" , bfac );
	      
	      writer.value( "DELTA_MASK" , dmask );
	      writer.value( "BETA_MASK" , bmask );
	      writer.value( "MASK" , mask );

	    }
	  
	  //
	  // Mask this epoch?
	  //
	  
	  if ( set_mask && mask )
	    {
	      if ( !edf.timeline.masked(e) ) ++altered;
	      edf.timeline.set_epoch_mask( e );
	      ++total;
	    }


	  writer.unepoch();
	  
	}

      if ( set_mask )
	logger << " masked " << total << " of " << ne << " epochs, altering " << altered << "\n";
      
      // # masked (i.e. actually), # masked, # total
      // altered , 

      writer.var( "FLAGGED_EPOCHS" , "Number of epochs failing Buckelmueller" );
      writer.var( "ALTERED_EPOCHS" , "Number of epochs actually masked" ); 
      writer.var( "TOTAL_EPOCHS" , "Number of epochs tested" );
      
      writer.value( "FLAGGED_EPOCHS" , total );
      writer.value( "ALTERED_EPOCHS" , altered );
      writer.value( "TOTAL_EPOCHS" , ne );

      writer.unlevel( globals::signal_strat );
      
    } // next signal    
  
  

  //
  // For now, do not return any annot_t
  //
  
  return NULL;


  //
  // In future, we may routinely expand all functions to return annotations
  //
  
  annot_t * a = edf.timeline.annotations.add( "Buckelmuller" );
  
  if ( a == NULL ) std::cout << "is null;\n";

  a->description = "Buckelmuller et al (2006) automatic artifact detection" ;

  for (int s=0;s<ns;s++)
    {
      //
      // only consider data tracks
      //
      
      if ( edf.header.is_annotation_channel( signals(s) ) ) 
	continue;
      
      for (int e=0;e<ne;e++)
	{
	  double dfac = delta[s][e] / delta_average[s][e];
	  double bfac = beta[s][e] / beta_average[s][e];
	  
	  bool reject =  dfac > delta_threshold || bfac > beta_threshold;
	  
	  if ( reject ) 
	    a->add( "buckelmuller:" + signals.label(s)  , edf.timeline.epoch(e) );
	}
      
    }
  
  return a;
  
}



//
// SIGSTATS
//

void  rms_per_epoch( edf_t & edf , param_t & param )
{

  
  // Hjorth parameters: H1, H2, H3
  // Optional: RMS, % clipped signals
  // Turning rate 

  
  // exclude EPOCH is more than 5% of points are clipped
  const double clip_threshold = 0.05;
  
  std::string signal_label = param.requires( "sig" );  

  bool verbose = param.has( "verbose" ) || param.has( "epoch") ;

  bool calc_rms = param.has( "rms" );

  bool calc_clipped = param.has( "clipped" );
  

  //
  // Calculate channel-level statistics?
  //

  bool cstats = param.has( "cstats" ) ;

  double ch_th = cstats ? param.requires_dbl( "cstats" ) : 0 ;

  bool cstats_all = ! param.has( "cstats-unmasked-only" );
  
  
  //
  // Optionally calculate turning rate
  //
  
  bool turning_rate = param.has( "tr" )  || param.has( "tr-epoch" ) || param.has( "tr-d" ) || param.has( "tr-smooth" );

  double tr_epoch_sec = 1.0;
  int    tr_d = 4;
  int    tr_epoch_smooth = 30; // +1 is added afterwards
  
  if ( turning_rate ) 
    {
      if ( param.has( "tr-epoch" ) ) tr_epoch_sec = param.requires_dbl( "tr-epoch" );
      if ( param.has( "tr-d" ) ) tr_d = param.requires_int( "tr-d" );
      if ( param.has( "tr-smooth" ) ) tr_epoch_smooth = param.requires_int( "tr-smooth" );
      logger << " calculating turning rate: d="<< tr_d << " for " << tr_epoch_sec << "sec epochs, smoothed over " << tr_epoch_smooth << " epochs\n";
    }


  //
  // allow for iterative outlier detection, i.e. with multiple
  // comma-delimited thresholds
  //
  
  bool apply_mask = false;

  std::vector<double> th;

  if ( param.has( "threshold" ) )
    {
      apply_mask = true;
      th = param.dblvector( "threshold" );
    }
  else if ( param.has( "th" ) )
    {
      apply_mask = true;
      th = param.dblvector( "th" );
    }
  
  int th_nlevels = th.size();
  
  //
  // channel/epoch (chep) masks; this will not set any full epoch-masks, but it will 
  // populate timeline.chep masks (i.e. that can be subsequently turned into lists of bad channels/epochs
  // and also used to interpolate signals
  //

  bool chep_mask = param.has( "chep" );

  // cannot apply both chep and mask options
  if ( chep_mask )
    {
      // check a masking threshold was specified above

      if ( ! apply_mask )
	Helper::halt( "no threshold ('th') specified for chep" );
      
      // cannot apply epoch mask and chep mask together,
      // so now set this to F anyway
      
      apply_mask = false;
    }
  

  //
  // Calculate per-EPOCH, and also signal-wide, the signal RMS 
  // Also calculate the proportion of flat/clipped signals (calculated per-EPOCH)
  //
  
  //
  // Attach signals
  //
  
  signal_list_t signals = edf.header.signal_list( signal_label );  

  const int ns = signals.size();

  //
  // Store per-epoch statistics
  //

  std::vector<int>      n( ns , 0 );

  std::vector<double> rms( ns , 0 );
  std::vector<double> clipped( ns , 0 );
  std::vector<double> mean_activity( ns , 0 );
  std::vector<double> mean_mobility( ns , 0 );
  std::vector<double> mean_complexity( ns , 0 );
  std::vector<double> mean_turning_rate( ns , 0 ); // nb. this is based on turning-rate epoch size


  //
  // Track original data, if calculating outliers, e_*, i.e. EPOCH level
  //
  
  std::vector<std::vector<double> > e_rms;
  std::vector<std::vector<double> > e_clp;
  std::vector<std::vector<double> > e_act;
  std::vector<std::vector<double> > e_mob;
  std::vector<std::vector<double> > e_cmp;
  std::vector<std::vector<double> > e_epoch;
  std::vector<std::vector<double> > e_tr;

  if ( apply_mask || cstats || chep_mask )
    {
      e_rms.resize( ns );
      e_clp.resize( ns );
      e_act.resize( ns );
      e_mob.resize( ns );
      e_cmp.resize( ns );
      e_epoch.resize( ns );
    }

  if ( turning_rate )
    e_tr.resize( ns );
  
  //
  // Point to first epoch (assume 30 seconds, but could be different)
  //
  
  int ne = edf.timeline.first_epoch();

  if ( ne == 0 ) return;
  
  
  //
  // For each signal
  //
  
  for (int s=0;s<ns;s++)
    {

      //
      // only consider data tracks
      //
      
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
      

      //
      // output stratifier
      //

      writer.level( signals.label(s) , globals::signal_strat );

      //
      // reset to first epoch
      //
      
      edf.timeline.first_epoch();

      //
      // Get sampling rate
      //
      
      int sr = edf.header.sampling_freq( s );
      
      //
      // for each each epoch 
      //
      
      while ( 1 ) 
	{
	  
	  //
	  // Get next epoch
	  //
	  
	  int epoch = edf.timeline.next_epoch();
	  
	  if ( epoch == -1 ) break;
	  
	  interval_t interval = edf.timeline.epoch( epoch );
	  	  
	  slice_t slice( edf , signals(s) , interval );
	  
	  std::vector<double> * d = slice.nonconst_pdata();
	  
	  //
	  // Mean-centre 30-second window, calculate RMS
	  //

	  MiscMath::centre( d );
	  
	  double x = calc_rms ? MiscMath::rms( *d ) : 0 ;
	  	  
	  double c = calc_clipped ? MiscMath::clipped( *d ) : 0 ; 
	  

	  //
	  // Hjorth parameters
	  //

	  double activity = 0 , mobility = 0 , complexity = 0;
	  
	  MiscMath::hjorth( d , &activity , &mobility , &complexity );


	  //
	  // Turning rate
	  //
	  
	  double turning_rate_mean = 0;
	  if ( turning_rate )
	    {

	      std::vector<double> subepoch_tr;
	      
	      turning_rate_mean = MiscMath::turning_rate( d , sr, tr_epoch_sec , tr_d , &subepoch_tr );
	      
	      for (int i=0;i<subepoch_tr.size();i++)
		e_tr[s].push_back( subepoch_tr[i] );
	    }

	  //
	  // Verbose output
	  //
	  
	  if ( verbose )
	    {

	      writer.epoch( edf.timeline.display_epoch( epoch ) );

	     
	      //
	      // Report calculated values
	      //
  
	      writer.value( "H1" , activity , "Epoch Hjorth parameter 1: activity (variance)" );
	      writer.value( "H2" , mobility , "Epoch Hjorth parameter 2: mobility" );
	      writer.value( "H3" , complexity , "Epoch Hjorth parameter 3: complexity" );

	      if ( calc_rms )
		writer.value( "RMS" , x , "Epoch root mean square (RMS)" );
	      
	      if ( calc_rms )
		writer.value( "CLIP" , c , "Proportion of epoch with clipped signal" );

	      if ( turning_rate ) 
		writer.value( "TR" , turning_rate_mean , "Turning rate mean per epoch" );
	    }


	  //
	  // Tot up for individual-level means
	  //

	  if ( calc_rms )
	    rms[s] += x;

	  if ( calc_clipped )
	    clipped[s] += c;
	  
	  mean_activity[s] += activity;
	  mean_mobility[s] += mobility;
	  mean_complexity[s] += complexity;

	  n[s]   += 1;
	  
	  
	  //
	  // Track for thresholding, or channel-stats ?
	  //

	  if ( apply_mask || cstats || chep_mask ) 
	    {
	      if ( calc_rms ) 
		e_rms[s].push_back( x );

	      if ( calc_clipped)
		e_clp[s].push_back( c );

	      e_act[s].push_back( activity );
	      e_mob[s].push_back( mobility );
	      e_cmp[s].push_back( complexity );

	      e_epoch[s].push_back( epoch );
	    }
	
	  
	  //
	  // Next epoch
	  //

	} 

      if ( verbose ) 
	writer.unepoch();
      
      //
      // Next signal
      //
      
    } 

  if ( verbose ) // did we have any output
    writer.unlevel( globals::signal_strat );

    
  
  //
  // Find outliers and mask epochs?
  //
  
  if ( apply_mask || chep_mask ) 
    {
      
      for (int s=0;s<ns;s++)
	{
	  
	  
	  //
	  // only consider data tracks
	  //
	  
	  if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
      

	  int cnt_rms = 0 , cnt_clp = 0 , cnt_act = 0 , cnt_mob = 0 , cnt_cmp = 0;
	  
	  if ( n[s] < 2 ) continue;
	  
	  writer.level( signals.label(s) , globals::signal_strat );

	  // for each iteration of outlier pruning
	  // track which epochs we are dropping here
	  std::vector<bool> dropped( n[s] , false ); 
	  
	  // total number of epochs dropped
	  int total = 0;
	  int altered = 0;

	  for (int o=0;o<th_nlevels;o++)
	    {

	      const int ns = n[s];

	      std::vector<double> act_rms;
	      std::vector<double> act_act;
	      std::vector<double> act_mob;
	      std::vector<double> act_cmp;

	      // calculate current mean for RMS and Hjorth parameters
	      // (clipping based on a single fixed threshold, not statistically)

	      for (int j=0;j<ns;j++)
		{
		  if ( ! dropped[j] ) 
		    {
		      if ( calc_rms )
			act_rms.push_back( e_rms[s][j] );

		      act_act.push_back( e_act[s][j] );
		      act_mob.push_back( e_mob[s][j] );
		      act_cmp.push_back( e_cmp[s][j] );
		    }
		}
	      
	      double mean_rms = calc_rms ? MiscMath::mean( act_rms ) : 0 ;
	      double mean_act = MiscMath::mean( act_act );
	      double mean_mob = MiscMath::mean( act_mob );
	      double mean_cmp = MiscMath::mean( act_cmp );
	      
	      double sd_rms = calc_rms ? MiscMath::sdev( act_rms , mean_rms ) : 1;
	      double sd_act = MiscMath::sdev( act_act , mean_act ); 
	      double sd_mob = MiscMath::sdev( act_mob , mean_mob ); 
	      double sd_cmp = MiscMath::sdev( act_cmp , mean_cmp ); 
	  
	      double lwr_rms = calc_rms ? mean_rms - th[o] * sd_rms : 0 ;
	      double lwr_act = mean_act - th[o] * sd_act;
	      double lwr_mob = mean_mob - th[o] * sd_mob;
	      double lwr_cmp = mean_cmp - th[o] * sd_cmp;
	  
	      double upr_rms = calc_rms ? mean_rms + th[o] * sd_rms : 0 ; 
	      double upr_act = mean_act + th[o] * sd_act;
	      double upr_mob = mean_mob + th[o] * sd_mob;
	      double upr_cmp = mean_cmp + th[o] * sd_cmp;
	      
	      logger << " RMS/Hjorth filtering " << edf.header.label[ signals(s) ] << ", threshold +/-" << th[o] << " SDs";
	      
	      const int ne = e_epoch[s].size();
	      
	      int total_this_iteration = 0;
	      
	      for (int ei=0; ei<ne; ei++)
		{
		  
		  // track which epochs were actually used
		  const int e = e_epoch[s][ei];
		  
		  // skip if this epoch is already masked
		  if ( dropped[ei] ) continue;
		  
		  bool set_mask = false;
		  
		  if ( calc_rms && ( e_rms[s][ei] < lwr_rms || e_rms[s][ei] > upr_rms ) ) 
		    {
		      set_mask = true;
		      cnt_rms++;
		    }
		  
		  // For clipping, use a fixed threshold
		  if ( calc_clipped && e_clp[s][ei] > clip_threshold ) 
		    {
		      set_mask = true;
		      cnt_clp++;
		    }
		  		  
		  if ( e_act[s][ei] < lwr_act || e_act[s][ei] > upr_act ) 
		    {
		      set_mask = true;
		      cnt_act++;
		    }
		  
		  if ( e_mob[s][ei] < lwr_mob || e_mob[s][ei] > upr_mob ) 
		    {
		      set_mask = true;
		      cnt_mob++;
		    }
		  
		  if ( e_cmp[s][ei] < lwr_cmp || e_cmp[s][ei] > upr_cmp )
		    {
		      set_mask = true;
		      cnt_cmp++;
		    }
		  
		  

		  //
		  // full mask
		  //
		  
		  if ( set_mask ) 
		    {
		      
		      if ( chep_mask ) // channel/epoch mask rather than standard epoch mask?
			{
			  edf.timeline.set_chep_mask( e , signals(s) );
			  // do not set 'dropped' here... 
			  // i.e. we'll consider all epochs for subsequent
			  //      signals
			  ++total_this_iteration;
			  ++total;
			}
		      else
			{			  
			  if ( !edf.timeline.masked(e) ) ++altered;
			  edf.timeline.set_epoch_mask( e );		  
			  dropped[ei] = true;
			  ++total_this_iteration;
			  ++total;
			}
		    }
		  
		} // next epoch	  
	      
	      logger << ": removed " << total_this_iteration 
		     << " of " << act_rms.size() 
		     << " epochs of iteration " << o+1 << "\n";

	    } // next outlier iteration
	  
	  
	  //
	  // report final epoch-level masking
	  //
	  
	  if ( verbose )
	    {	      
	      const int ne = e_epoch[s].size();
	      for (int ei = 0 ; ei < ne ; ei++ )		
		{
		  writer.epoch( edf.timeline.display_epoch( e_epoch[s][ei] ) );
		  writer.value( "MASK" , dropped[ei] ? 1 : 0 , "Masked epoch? (1=Y)" ); 
		}
	      writer.unepoch();
	    }
	  
	  logger << " Overall, masked " << total << " of " << ne << " epochs:"
		 << " ACT:" << cnt_act 
		 << " MOB:" << cnt_mob 
		 << " CMP:" << cnt_cmp ;
	  if ( calc_rms ) logger << " RMS:" << cnt_rms ;
	  if ( calc_clipped ) logger << " CLP:" << cnt_clp;
	  logger << "\n";

	  
	  if ( calc_rms ) writer.value( "CNT_RMS" , cnt_rms , "Epochs failing RMS filter" );
	  if ( calc_clipped) writer.value( "CNT_CLP" , cnt_clp , "Epochs failing CLIP filter" );
	  writer.value( "CNT_ACT" , cnt_act , "Epochs failing H1 filter" );
	  writer.value( "CNT_MOB" , cnt_mob , "Epochs failing H2 filter" );
	  writer.value( "CNT_CMP" , cnt_cmp , "Epochs failing H3 filter" );

	  writer.value( "FLAGGED_EPOCHS" , total,    "Number of epochs failing SIGSTATS" );
	  writer.value( "ALTERED_EPOCHS" , altered , "Number of epochs actually masked, i.e. not already masked" );
	  writer.value( "TOTAL_EPOCHS"   , ne,       "Number of epochs tested" );
	
	
	} // next signal      
      
      writer.unlevel( globals::signal_strat );

    }
  

  //
  // Channel level stats (i.e. between channel comparisons, within each epoch)
  //  per epoch
  //   then, averaged all epochs
  //         only for unmasked epochs
  // this only makes sense if we have multiple channels
  //
  
  if ( cstats && ns > 2 )
    {

      logger << "  calculating between-channel statistics, ";
      if ( cstats_all ) logger << "based on all epochs\n";
      else logger << "based only on unmasked epochs\n";
      logger << "  threshold (for P_H1, P_H2, P_H3) is " << ch_th << " SD units\n";
      
      // mean over epochs (Z score of this channel versus all others)
      std::vector<double> m_ch_h1(ns), m_ch_h2(ns), m_ch_h3(ns);

      // number/proportion of epochs where channel has |Z| > ch_th;
      std::vector<double> t_ch_h1(ns), t_ch_h2(ns), t_ch_h3(ns);

      int ne_actual = 0;
      
      for (int ei=0; ei<ne; ei++)
	{

	  // only consider unmasked epochs here?
	  if ( ! cstats_all )
	    if ( edf.timeline.masked( ei ) ) continue;
	  
	  std::vector<double> tmp_h1(ns), tmp_h2(ns), tmp_h3(ns);
	  for (int s=0;s<ns;s++)
	    {
	      tmp_h1[s] = e_act[s][ei];
	      tmp_h2[s] = e_mob[s][ei];
	      tmp_h3[s] = e_cmp[s][ei];	      
	    }

	  // normalize
	  tmp_h1 = MiscMath::Z( tmp_h1 );
	  tmp_h2 = MiscMath::Z( tmp_h2 );
	  tmp_h3 = MiscMath::Z( tmp_h3 );

	  // accumulate
	  for (int s=0;s<ns;s++)
	    {
	      tmp_h1[s] = fabs( tmp_h1[s] );
	      tmp_h2[s] = fabs( tmp_h2[s] );
	      tmp_h3[s] = fabs( tmp_h3[s] );
	      
	      if ( tmp_h1[s] > ch_th ) ++t_ch_h1[s];
	      if ( tmp_h2[s] > ch_th ) ++t_ch_h2[s];
	      if ( tmp_h3[s] > ch_th ) ++t_ch_h3[s];

	      m_ch_h1[s] += tmp_h1[s];
	      m_ch_h2[s] += tmp_h2[s];
	      m_ch_h3[s] += tmp_h3[s];
	      
	    }

	  // track epochs included in analysis
	  
	  ++ne_actual;

	  // next epoch
	}
      
      // normalize by number of epochs

      for (int s=0;s<ns;s++)
	{
	  
	  if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
	  
	  m_ch_h1[s] /= (double)ne_actual;
	  m_ch_h2[s] /= (double)ne_actual;
	  m_ch_h3[s] /= (double)ne_actual;
	  
	  t_ch_h1[s] /= (double)ne_actual;
	  t_ch_h2[s] /= (double)ne_actual;
	  t_ch_h3[s] /= (double)ne_actual;

	  writer.level( signals.label(s) , globals::signal_strat );

	  // all epochs
	  writer.value( "Z_H1" , m_ch_h1[s] );
	  writer.value( "Z_H2" , m_ch_h2[s] );
	  writer.value( "Z_H3" , m_ch_h3[s] );
	  
	  writer.value( "P_H1" , t_ch_h1[s] );
	  writer.value( "P_H2" , t_ch_h2[s] );
	  writer.value( "P_H3" , t_ch_h3[s] );

	  // only unmasked epochs	  
	  
	}
      
      writer.unlevel( globals::signal_strat );

    }
      
  

  //
  // Turning rate sub-epoch level reporting (including smoothing over sub-epochs)
  //

  if ( turning_rate )
    {

      for (int s=0;s<ns;s++)
	{
	  int sr = edf.header.sampling_freq( s );
	  
	  // how many units (in # of sub-epoch units);  +1 means includes self
	  int winsize = 1 + tr_epoch_smooth / tr_epoch_sec ; 

	  logger << "sz = " << e_tr[s].size() << " " << winsize << "\n";
	  e_tr[s] = MiscMath::moving_average( e_tr[s] , winsize );

	  // output
	  writer.level( signals.label(s) , globals::signal_strat );
	  for (int i=0;i<e_tr[s].size();i++)
	    {
	      writer.level( i+1 , "SUBEPOCH" );
	      writer.value( "TR" , e_tr[s][i] );
	    }
	  writer.unlevel( "SUBEPOCH" );
	}
      writer.unlevel( globals::signal_strat );
    }


  
  //
  // Individual level summary
  //
  
  for (int s=0;s<ns;s++)
    {

      //
      // only consider data tracks
      //
      
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
      
      writer.level( signals.label(s) , globals::signal_strat );
      
      writer.value( "H1"   , mean_activity[s] / (double)n[s] , "Mean H1 statistic" );
      writer.value( "H2"   , mean_mobility[s] / (double)n[s] , "Mean H2 statistic" );
      writer.value( "H3"   , mean_complexity[s] / (double)n[s] , "Mean H3 statistic" );

      if ( calc_clipped )
	writer.value( "CLIP" , clipped[s] / (double)n[s] , "Mean CLIP statistic" );
      if ( calc_rms )
	writer.value( "RMS"   , rms[s] / (double)n[s] , "Mean RMS statistic" );      
    }

  writer.unlevel( globals::signal_strat );

}




void  mse_per_epoch( edf_t & edf , param_t & param )
{
  

  //
  // Calculate MSE per epoch and average
  //

  //
  // MSE parameters
  //

  int    m  = param.has("m") ? param.requires_int("m") : 2;  

  double r  = param.has("r") ? param.requires_dbl("r") : 0.15;
  
  std::vector<int> scale;
  if ( param.has("s") ) 
    {
      scale = param.intvector("s");
      if ( scale.size() != 3 ) Helper::halt( "mse s=lwr,upr,inc" );
    }
  else
    {
      scale.push_back(1);
      scale.push_back(10);
      scale.push_back(2);
    }

  //
  // Output
  //
    
  bool verbose = param.has( "verbose" ) ;

  
  //
  // Attach signal(s)
  //

  std::string signal_label = param.requires( "sig" );  
  
  signal_list_t signals = edf.header.signal_list( signal_label );  

  const int ns = signals.size();

  
  
  //
  // For each signal  
  //
  
  for (int s=0;s<ns;s++)
    {

      //
      // only consider data tracks
      //
      
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;


      logger << " estimating MSE for " << signals.label(s) << "\n";
      
      //
      // output stratifier
      //
      
      writer.level( signals.label(s) , globals::signal_strat );	  
      
      //
      // to track overall mean over epochs: scale -> vector of per-epoch MSEs
      //
      
      std::map<int,std::vector<double> > all_mses;


      //
      // Point to first epoch (assume 30 seconds, but could be different)
      //
      
      int ne = edf.timeline.first_epoch();
      
      if ( ne == 0 ) return;
      

      //
      // for each each epoch 
      //

      while ( 1 ) 
	{
	  
	  //
	  // Get next epoch
	  //
	  
	  int epoch = edf.timeline.next_epoch();
	  
	  if ( epoch == -1 ) break;
	  
	  interval_t interval = edf.timeline.epoch( epoch );
	  	  
	  //
	  // get data
	  //

	  slice_t slice( edf , signals(s) , interval );
	  
	  std::vector<double> * d = slice.nonconst_pdata();

	  
	  //
	  // Mean-centre 30-second window, calculate RMS
	  //

	  mse_t mse( scale[0] , scale[1] , scale[2] , m , r );
	  
	  std::map<int,double> mses = mse.calc( *d );
	  
	  //
	  // track
	  //
	  

	  if ( verbose )
	    writer.epoch( edf.timeline.display_epoch( epoch ) );
	  
	  std::map<int,double>::const_iterator ii = mses.begin();
	  while ( ii != mses.end() )
	    {
	      const int & scale = ii->first ; 
	      
	      all_mses[  scale ].push_back( ii->second ); 
	      
	      // verbose output?

	      if ( verbose )
		{
		  writer.level( scale , "SCALE" );
		  writer.value( "MSE" , ii->second );      	  
		}
	      
	      // next scale
	      ++ii;
	    }
	  
	  if ( verbose )
	    writer.unlevel( "SCALE" );
	      
	} // next epoch
  
      if ( verbose ) 
	writer.unepoch();
      
      
      // overall means
      
      std::map<int,std::vector<double> >::const_iterator ii = all_mses.begin();
      while ( ii != all_mses.end() )
	{
	  double mse = 0;
	  const int & scale = ii->first;
	  const std::vector<double> & x = ii->second;
	  for (int i=0;i<x.size();i++) mse += x[i];
	  mse /= (double)x.size();
	  
	  writer.level( scale , "SCALE" );
	  writer.value( "MSE" , mse );
	  
	  ++ii;
	}
      writer.unlevel( "SCALE" );
      
    } // next signal
  
  writer.unlevel( globals::signal_strat );
  
}




void lzw_per_epoch( edf_t & edf , param_t & param )
{
  
  //
  // Calculate LZW per epoch and average
  //

  //
  // LZW parameters
  //

  int nbins   = param.has( "nbins" ) ? param.requires_int( "nbins" ) : 20 ;
  int nsmooth = param.has( "nsmooth" ) ? param.requires_int( "nsmooth" ) : 1 ;

  bool epoched = param.has( "epoch" );

     
  //
  // Attach signal(s)
  //

  std::string signal_label = param.requires( "sig" );  
  
  signal_list_t signals = edf.header.signal_list( signal_label );  

  const int ns = signals.size();


  //
  // Point to first epoch (assume 30 seconds, but could be different)
  //
  
  int ne = edf.timeline.first_epoch();

  if ( ne == 0 ) return;
  
  
  //
  // For each signal  
  //
  
  for (int s=0;s<ns;s++)
    {
      
      //
      // only consider data tracks
      //
      
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;


      //
      // output stratifier
      //
      
      writer.level( signals.label(s) , globals::signal_strat );	  
      

      //
      // per-epoch, or whole-signal calculation?
      //
      
      if ( ! epoched ) 
	{

	  // get all data
	  interval_t interval = edf.timeline.wholetrace();
	  
	  slice_t slice( edf , s , interval );
	  const std::vector<double> * d = slice.pdata();
	  const int sz = d->size();

	  // designed for per-epoch data, but just use first 
	  // slot for entire signal
	  std::vector<std::vector<double> > track_lzw;
	  track_lzw.push_back( *d );
	  	  
	  // coarse-grain signal
	  coarse_t c( track_lzw , nbins , nsmooth );
	  
	  // compress	  
	  lzw_t lzw( c );
      
	  // index	  
	  double index = lzw.size(0) / (double)track_lzw[0].size();

	  // output
	  writer.value( "LZW" , index );
	  	  
	}
     

      //
      // Epoch level analyses
      //
     
      if ( epoched )
	{
	  
	  
	  int ne = edf.timeline.first_epoch();
	  
	  if ( ne == 0 ) return;

	  // track all signals here

	  std::vector<std::vector<double> > track_lzw;
	  std::vector<int> track_e;

	  //
	  // for each each epoch 
	  //
	  
	  while ( 1 ) 
	    {
	      
	      // Get next epoch
	  
	      int epoch = edf.timeline.next_epoch();
	  
	      if ( epoch == -1 ) break;
	  
	      interval_t interval = edf.timeline.epoch( epoch );
	      
	      // get data
	      
	      slice_t slice( edf , signals(s) , interval );
	      
	      std::vector<double> * d = slice.nonconst_pdata();
	      
	      // lzw_t class is designed for per-epoch data to be taken 
	      // all in one structure	      
	      track_lzw.push_back( *d );
	      track_e.push_back( epoch );
	      
	    } // next epoch	      
	      	      
	  // coarse-grain signal
	  coarse_t c( track_lzw , nbins , nsmooth );
	  
	  // compress	  
	  lzw_t lzw( c );
      
	  // index	  
	  for (int e=0; e<track_e.size(); e++)
	    {

	      double index = lzw.size(e) / (double)track_lzw[e].size();
	      
	      // output
	      writer.epoch( edf.timeline.display_epoch( track_e[e] ) );
	      writer.value( "LZW" , index );
	      writer.unepoch();
			     
	    }
	  
	}
      
      writer.unlevel( globals::signal_strat );
    } // next signal
  
}



void    spike_signal( edf_t & edf , int s1 , int s2 , double wgt , const std::string & ns )
{

  if ( s1 == s2 ) return;
  if ( edf.header.is_annotation_channel(s1) ) Helper::halt( "annotation channel specified for SPIKE" ) ;
  if ( edf.header.is_annotation_channel(s2) ) Helper::halt( "annotation channel specified for SPIKE" ) ;

  bool append_new_channel = ns != "";
  
  interval_t interval = edf.timeline.wholetrace();
  
  // Right now, need a similar sampling rate
  
  const int Fs1 = edf.header.sampling_freq( s1 );
  const int Fs2 = edf.header.sampling_freq( s2 );
  
  const std::string label1 = edf.header.label[s1];
  const std::string label2 = edf.header.label[s2];

  if ( Fs1 != Fs2 ) 
  {
    logger << "Note: resampling " << label2 << " to " << Fs1 << " to match " << label1 << "\n";
    dsptools::resample_channel( edf, s2 , Fs1 );
  }

  slice_t slice1( edf , s1 , interval );
  const std::vector<double> * d1 = slice1.pdata();
  const int sz1 = d1->size();
 
  slice_t slice2( edf , s2 , interval );
  const std::vector<double> * d2 = slice2.pdata();
  const int sz2 = d2->size();
  
  if ( sz1 != sz2 ) Helper::halt( "problem in SPIKE, unequal channel lengths" );
      
  // apply SPIKE
  std::vector<double> spiked( sz1 , 0);
  for (int i=0;i<sz1;i++)
    spiked[i] = (*d1)[i] + wgt * (*d2)[i];
  
  
  // either UPDATE
  if ( append_new_channel )
    {
      const std::string label = edf.header.label[s1] + "-spike-" + edf.header.label[s1] + "-wgt-" + Helper::dbl2str( wgt );      
      edf.add_signal( label , Fs1 , spiked ); 
    }
  else // ... else UPDATE
    edf.update_signal( s1 , &spiked );

}
