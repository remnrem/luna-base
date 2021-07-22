
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


#include "spindles.h"
#include "mspindles.h"
#include "plot-spindles.h"
#include "propag.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "eval.h"

#include "annot/annot.h"
#include "intervals/intervals.h"

#include "miscmath/dynam.h"
#include "cwt/cwt.h"
#include "fftw/fftwrap.h"
#include "miscmath/miscmath.h"
#include "artifacts/artifacts.h"
#include "dsp/fir.h"
#include "dsp/hilbert.h"
#include "dsp/slow-waves.h"
#include "defs/defs.h"
#include "db/db.h"

#include "helper/logger.h"
#include "helper/helper.h"

// output
extern writer_t writer;
extern logger_t logger;


annot_t * spindle_wavelet( edf_t & edf , param_t & param )
{

  //
  // Optionally collate all spindles called here, either across all channels:
  //

  mspindles_t mspindles( &edf );

  // ... or within channels only

  std::map<std::string,mspindles_t> ch2mspindles; 
  

  //
  // Wavelet parameters
  //
  
  // center frequencies for wavelets
  std::vector<double> frq;  
  
  if ( param.has( "fc" ) )
    {
      frq = param.dblvector( "fc" );
    }
  else if ( param.has("fc-lower") )
    {
      
      const double   fc_lower = param.requires_dbl("fc-lower") ;

      const double   fc_upper = param.requires_dbl("fc-upper") ;

      const double   fc_step = param.requires_dbl("fc-step"); 
	  
      for (double fc = fc_lower ; fc <= fc_upper ; fc += fc_step ) frq.push_back( fc );
    }
  else
    frq.push_back( 13.5 );


  //
  // Alternate specication, w/ FWHM?
  //
  
  const bool alt_spec = param.has( "fwhm" ) ;

  std::vector<double> fwhm;
  
  if ( alt_spec ) fwhm = param.dblvector( "fwhm" );
  if ( fwhm.size() == 1 && frq.size()>1 ) fwhm.resize( frq.size() , fwhm[0] );

  //
  // number of cycles
  //

  const int num_cycles = param.has("cycles" ) ? param.requires_int( "cycles" ) : 7 ;

  if ( param.has( "cycles" ) && alt_spec ) Helper::halt( "use either fwhm or cycles" );
  
  //
  // Detection parameters
  //
  
  // empirical threshold determination

  bool estimate_empirical_thresholds = param.has( "empirical" ) || param.has( "set-empirical" );
  
  bool use_empirical_thresholds = param.has( "set-empirical" );

  bool verbose_empirical = param.has("verbose-empirical");

  // use local peak-finding threshold method ( still uses th/min0 and th2/min, and max
  const bool use_zpks = param.has( "zpks" );
  const double zpks_window_sec = use_zpks ? param.requires_dbl( "zpks" ) : 0;
  const double zpks_influence = param.has( "influence" ) ? param.requires_dbl( "influence" ) : 0.01;

  // default multiplicative threshold for spindle core, = 4.5 
  double   multiplicative_threshold = param.has( "th" ) ? param.requires_dbl( "th" ) : 4.5 ;
  
  // default multiplicative threshold for spindle core + flank
  double   boundary_threshold       = param.has( "th2" ) ? param.requires_dbl( "th2" ) : 2.0   ;

  // (optional) upper bound for core amplitude threshold, e.g. 4.5 < x < 10  if th-max=10
  double   maximal_threshold        = param.has( "th-max" ) ? param.requires_dbl( "th-max" ) : -9 ;
  
  // minimum spindle core duration (relates to 'th')
  const double   min0_dur_sec              = param.has( "min0" ) ? param.requires_dbl( "min0" ) : 0.3; 
  
  // minimum spindle duration (core + flanking) (relates to 'th2' boundary threshold)
  const double   min_dur_sec              = param.has( "min" ) ? param.requires_dbl( "min" ) : 0.5; 
  
  // maximum spindle duration (core + flanking)
  const double   max_dur_sec              = param.has( "max" ) ? param.requires_dbl( "max" ) : 3.0;
  
  // default 0.1 seconds smoothing of CWT coefficients prior to thresholding
  const double   moving_window_sec        = param.has( "win" ) ? param.requires_dbl( "win" ) : 0.1;
  
  // merge spindles that are within 0.5 sec, by default
  const uint64_t spindle_merge_tp         = ( param.has( "merge" ) ? 
					      param.requires_dbl( "merge" ) : 0.5 ) * globals::tp_1sec; 

  // use CWT median instead of mean when determining thresholds
  const bool     use_median               = param.has( "median" ) ;
  
  // instead of study-wide threshold, adopt a local epoch-based norms
  const int      epoch_norm_sec           = param.has( "local" ) ? param.requires_dbl( "local" ) : 0 ;


  // convert duration thresholds in 'tp' units
  const uint64_t min0_dur_tp = min0_dur_sec * globals::tp_1sec;
  const uint64_t min_dur_tp  = min_dur_sec * globals::tp_1sec;
  const uint64_t max_dur_tp  = max_dur_sec * globals::tp_1sec;
  

  //
  // Analysis/output parameters
  //

  // epoch-level output
  const bool     show_epoch_level         = param.has( "epoch" );

  // spindle-level output
  const bool     show_spindle_level       = param.has( "per-spindle" );
  
  
  // verbose display of all CWT coefficients
  const bool     show_cwt_coeff           = param.has( "show-coef" );

  // detect slow waves and estimate ITPC etc for spindle start/peak/stop and slow waves
  const bool     sw_coupling              = param.has( "so" );

  // show SPINDLES in sample-pints
  const bool     show_sample_points       = param.has( "sp" );

  // intra-spindle frequency changes via HT
  const bool     ht_chirp                 = param.has( "if" );
  
  const double   ht_chirp_frq             = param.has( "if-frq" ) ? param.requires_dbl( "if-frq" ) : 2.0;
  const double   ht_chirp_frq2            = param.has( "if-frq2" ) ? param.requires_dbl( "if-frq2" ) : 0 ;
  const double   ht_chirp_frq_emp         = param.has( "if-frq-emp" ) ? param.requires_dbl( "if-frq-emp" ) : 0 ;

  const int      ht_bins                  = 5; // divide spindle interval into 'n' equal size bins
  const bool     ht_verbose               = param.has( "verbose-if" );
  
  // time-locked signal means
  const bool     tlocking                 = param.has( "tlock" );

  const bool     verbose_time_phase_locking = param.has( "verbose-coupling" );

  // generate a feature file of spindles
  const bool     save_annots                = param.has( "annot" );
  
  // show verbose ENRICH output
  const bool     enrich_output            = param.has( "enrich" );


  //
  // Caches:
  //   - spindle peaks (cache-peaks)
  //   - wavelet power (cache)
  //   - metrics (DENS, etc)
  //
  
  const bool cache_data                   = param.has( "cache" );
  const std::string cache_name = cache_data ? param.value( "cache" ) : "" ;
  
  const bool cache_peaks                  = param.has( "cache-peaks" );
  const std::string cache_peaks_name = cache_peaks ? param.value( "cache-peaks" ) : "";

  cache_t<double> * cache_metrics = param.has( "cache-metrics" ) ? edf.timeline.cache.find_num( param.value( "cache-metrics" ) ) : NULL ;

  
  //
  // Spindle propagation
  //

  const bool do_prop = param.has( "prop" );
  
  sp_props_t props;
  
  //
  // Intersection of multiple wavelets/spindles/channels   ( by default, do not merge across channels)
  //
  
  
  const bool do_channel_merge = param.has( "collate-within-channel" ) ;

  const bool do_merge = do_channel_merge || param.has( "collate" ); 

  
  if ( do_merge )
    {
      // merge two spindles if do_mergeion / union > t
      mspindles.interval_th = param.has( "th-interval" ) ? param.requires_dbl( "th-interval" ) : 0.0;
      
      // for spindles of same frequency, different channels: merge two spindles if intersection / union > t
      mspindles.cross_ch_interval_th = param.has( "th-interval-cross-channel" ) ? param.requires_dbl( "th-interval-cross-channel" ) : 0.0;
      
      // for spindles of different frequencies, same channels: merge two spindles if intersection / union > t
      mspindles.within_ch_interval_th = param.has( "th-interval-within-channel" ) ? param.requires_dbl( "th-interval-within-channel" ) : 0.0;
      
      // merge if within this frequency range (Hz)
      mspindles.frq_th = param.has( "th-frq" ) ? param.requires_dbl( "th-frq" ) : 2; 
      
      // add window around each spindle interval
      mspindles.window = param.has( "window" ) ? param.requires_dbl( "window" ) : 0 ;

      // add MSP_START_HMS and MSP_STOP_HMS to the output
      mspindles.hms = param.has( "hms" );

      // verbose output: for each MSPINDLE list all contributing spindles
      mspindles.per_spindle_verbosity = param.has( "list-all-spindles" );
    }

 

  //
  // Output
  //

  bool hms = param.has("hms");

  clocktime_t starttime( edf.header.starttime );
  if ( ! starttime.valid ) 
    {
      logger << " ** could not find valid start-time in EDF header **\n";
      hms = false;
    }
  



  //
  // Signals
  // 
  
  std::string signal_label = param.requires( "sig" );   
  
  // list of signals
  signal_list_t signals = edf.header.signal_list( signal_label );  
  
  // number of signals
  const int ns = signals.size();

  // nothing to do...
  if ( ns == 0 ) return NULL;

  // sampling rate
  std::vector<double> Fs = edf.header.sampling_freq( signals );
  
 
  //
  // Set up annotation
  //
  
  std::string sp_label = "spindles";

  if ( save_annots ) 
    {
      if ( param.value( "annot" ) != "" )
	sp_label = param.value( "annot" );       
    }

  
  //
  // Add new channels?
  //

  bool add_channels = param.has( "add-channels" );


  //
  // Per-spindle characterization? 
  //

  bool characterize = ! param.has( "no-spindle-stats" );   


  //
  // Draw spindles
  //

  bool spindle_pdf = param.has( "pdf" );
  


  //
  // Save annotation file ( x-ID.annot)
  //

  std::string annotfile = "";

  if ( param.has( "out" ) )
    {
      annotfile = param.value("out");
      logger << "  writing annotation files [" << annotfile << "]\n";
    }

  
  //
  // Keeping track of 'all' overlaps?
  //
  
  std::set<feature_t> all_spindles;
  
  
  //
  // For each signal, over the whole signal
  //
  
  interval_t interval = edf.timeline.wholetrace(); 
  
  for (int s = 0 ; s < ns ; s++ ) 
    {
      
      //
      // Only consider raw signal channels
      //
      
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
      
      //
      // Output
      //
      
      writer.level( signals.label(s) , globals::signal_strat );
      

      //
      // Pull all data
      //
      
      slice_t slice( edf , signals(s) , interval );
      
      const std::vector<double> * d = slice.pdata();
      
      const std::vector<uint64_t> * tp = slice.ptimepoints();
      
      const int np0 = d->size();

      uint64_t dt = 1.0/Fs[s] * globals::tp_1sec; // time in tp-units
      
      double dt_minutes = (double)dt / ( 60 * globals::tp_1sec ) ;
      
      double t_minutes = d->size() * dt_minutes; // total trace time in minutes

      //
      // Run CWT 
      //

      CWT cwt;
      
      cwt.set_sampling_rate( Fs[s] );

      
      for (int fi=0;fi<frq.size();fi++)
	{
	  if ( alt_spec )
	    cwt.alt_add_wavelet( frq[fi] , fwhm[fi] , 10 );  // f( Fc , FWHM , 10 seconds window (fixed number of cycles ) 
	  else
	    cwt.add_wavelet( frq[fi] , num_cycles );  // f( Fc , number of cycles ) 
	}
      
      cwt.load( d );

      cwt.run();


      //
      // Run baseline FFT on the entire signal 
      //
      
      std::map<freq_range_t,double> baseline_fft;

      do_fft( d , Fs[s] , &baseline_fft );


      //
      // This is only populated if we are considering multiple frequencies
      //
      
      std::map<double,std::vector<spindle_t> > f2int;      


      //
      // Set up for optional slow-wave coupling
      //

      hilbert_t * p_hilbert = NULL ;
      
      slow_waves_t * p_sw = NULL ;
    
      
      if ( sw_coupling )
	{
	  
	  // freq
	  double flwr = param.has( "f-lwr" ) ? param.requires_dbl( "f-lwr" ) : 0.2;
	  double fupr = param.has( "f-upr" ) ? param.requires_dbl( "f-upr" ) : 4.5;
	  
	  // time
 	  double tlwr = param.has( "t-lwr" ) ? param.requires_dbl( "t-lwr" ) : 0;
	  double tupr = param.has( "t-upr" ) ? param.requires_dbl( "t-upr" ) : 0;
	  
	  // time of negative peak only
	  double t_neg_lwr = param.has( "t-neg-lwr" ) ? param.requires_dbl( "t-neg-lwr" ) : 0;
	  double t_neg_upr = param.has( "t-neg-upr" ) ? param.requires_dbl( "t-neg-upr" ) : 0;

	  // absolute magnitude
	  double uV_neg = param.has( "uV-neg" ) ? param.requires_dbl( "uV-neg" ) : 0 ;
	  double uV_p2p = param.has( "uV-p2p" ) ? param.requires_dbl( "uV-p2p" ) : 0 ;

	  // relative magnitude 
	  double mag  = param.has( "mag" ) ? param.requires_dbl( "mag" ) : 0 ;
	  bool   use_mean = param.has( "so-mean" ); 
	  bool   ignore_neg_peak = param.has( "ignore-neg-peak" );

	  // for full wave detection, count based on positive-to-negative zero-crossings 
	  // (i.e. negative wave first), or the other way

	  bool use_alternate_neg2pos_zero_crossing = param.has( "neg2pos" ) ;
	  
	  slow_wave_type so_type = SO_FULL;
	  if      ( param.has( "half-wave" ) ) so_type = SO_HALF;
	  else if ( param.has( "negative-half-wave" ) ) so_type = SO_NEGATIVE_HALF;
	  else if ( param.has( "positive-half-wave" ) ) so_type = SO_POSITIVE_HALF;
	  
	  double ripple = 0.01;
	  double tw = 0.5;
	  // filter-Hilbert raw signal for SWs
	  p_hilbert = new hilbert_t( *d , Fs[s] , flwr , fupr , ripple , tw );
	  
	  std::vector<double> ph_peak;
	  
	  // are spindles in slow-waves?
	  std::vector<bool> sw_peak;
	  
	  // find slow-waves	      
	  p_sw = new slow_waves_t( *d , *tp , Fs[s] , mag, ignore_neg_peak , use_mean , uV_neg , uV_p2p , flwr, fupr, 
				   tlwr, tupr, t_neg_lwr , t_neg_upr ,  
				   use_alternate_neg2pos_zero_crossing , so_type );

	  // and phase
	  p_sw->phase_slow_waves();
	  
	  // and display (& potentially cache)
	  p_sw->display_slow_waves( param.has( "verbose" ) , &edf , cache_metrics );
	  

	  if ( verbose_time_phase_locking ) 
	    {
	      
	      //
	      // Time-locked average of ( raw | BPF filtered ) EEG 
	      //
	      
	      
	      // unfiltered EEG (i.e prior to SO BPF)
	      const std::vector<double> * dsig = d;
	      
	      // an alternative would be to plot the BPF version
	      //const std::vector<double> * dsig = p_sw->p_filtered() ;
	      
	      // 36 bins = 10-degree bins; 18 = 20-deg bins
	      int nbins = 36/2;
	      std::vector<double> pl_eeg     = p_sw->phase_locked_averaging( dsig , nbins );
	      
	      if ( pl_eeg.size() > 0 ) // no SO detected anyway
		{
		  writer.var( "SOPL_EEG" , "Slow wave phase-locked average EEG" );
		  
		  double inc = 360 / (double)nbins;
		  double ph = inc/2.0; // use mid-point of range
		  
		  for (int j=0;j<nbins;j++)
		    {
		      writer.level( ph  , "PHASE" );
		      writer.value( "SOPL_EEG" , pl_eeg[j] );
		      ph += inc;
		    }
		  writer.unlevel( "PHASE" );
		}
	      
	      
	      //
	      // Time-locked signal averaging, to negative peak of SO
	      //
	      
	      std::vector<double> tl_eeg = p_sw->time_locked_averaging( dsig , Fs[s] , 1 , 1 );
	      
	      writer.var( "SOTL_EEG" , "Slow wave time-locked average EEG" );
	      
	      int sz = tl_eeg.size();
	      
	      if ( sz > 0 ) // 0 if no SW were detected in the first place
		{
		  int sz2 = - (sz-1)/2;
		  
		  for (int j=0;j<sz;j++)
		    {
		      writer.level( sz2 , "SP" );
		      writer.value( "SOTL_EEG" , tl_eeg[j] );		  
		      ++sz2;
		    }
		  writer.unlevel( "SP" );
		}
	      
	    }
	  
	}

      
      //
      // Now consider results for each Fc separately
      //
          
      for (int fi=0; fi<frq.size(); fi++)
	{
	  
	  logger << "\n  detecting spindles around F_C " << frq[fi] << "Hz for " << signals.label(s) << "\n";

	  if ( alt_spec ) 
	    logger << "  wavelet with FWHM(T) " << fwhm[fi] << "\n";       
	  else
	    logger << "  wavelet with " << num_cycles << " cycles\n";       

	  logger << "  smoothing window = " << moving_window_sec << "s\n";

	  //
	  // Output stratifier: F_C
	  //

	  writer.level( frq[fi] , globals::freq_strat );
	  
	  //
	  // Get results for this F_C
	  //

	  const std::vector<double> & results = cwt.results(fi);
      


	  //
	  // Get a moving average of the result, 0.1 windows, and get mean
	  //
	  
	  int window_points = moving_window_sec * Fs[s];
	  
	  if ( window_points % 2 == 0 ) ++window_points;
	  
	  const std::vector<double> averaged 
	    = MiscMath::moving_average( results , window_points );

	  const double mean = use_median ? MiscMath::median( averaged ) : MiscMath::mean( averaged );

	  
		  
	  //
	  // Find empirical threshold
	  //
	  
	  if ( estimate_empirical_thresholds )
	    {
	      std::vector<double> adj_vals(  averaged.size() );
	      for (int i=0;i<adj_vals.size();i++) adj_vals[i] = averaged[i] / mean; // might be median...
	      std::map<double,double> tvals;
	      double thf = 0;
	      double empirical_threshold = MiscMath::threshold( adj_vals , 0.25 , 20 , 0.25 , &thf, &tvals );

	      if ( verbose_empirical ) 
		{
		  for (int i=0; i<adj_vals.size();i++)
		    std::cout << "AV\t" << adj_vals[i] << "\n";
		}

	      logger << "  estimated empirical thresholds as " << empirical_threshold << "\n";

	      if ( use_empirical_thresholds )
		{
		  logger << "  setting thresholds to empirical value, " << empirical_threshold << "\n";
		  multiplicative_threshold = empirical_threshold;
		  boundary_threshold = multiplicative_threshold * 0.5;
		  maximal_threshold = multiplicative_threshold * 10;
		}

	      writer.var( "EMPTH" , "Empirical threshold" );
	      writer.value( "EMPTH" , empirical_threshold );
	      
	      writer.var( "EMPF" , "Empirical threshold frequency" );
	      writer.value( "EMPF" , thf );

	      // redundant code, but fine to leave as is for now...
	      double t_mean = MiscMath::mean( averaged );
	      double t_median = MiscMath::median( averaged );
	      // std::vector<double> mean_adj_vals(  averaged.size() );
	      // std::vector<double> median_adj_vals(  averaged.size() );
	      // for (int i=0;i<adj_vals.size();i++) mean_adj_vals[i] = averaged[i] / t_mean; 
	      // for (int i=0;i<adj_vals.size();i++) median_adj_vals[i] = averaged[i] / t_median; 
	      
	      if ( t_median > 0 ) 
		writer.value( "MEAN_OVER_MEDIAN" , t_mean / t_median );
	      
	      // writer.value( "ADJMEAN_MEAN" , MiscMath::mean( mean_adj_vals ) );
	      // writer.value( "ADJMEAN_MEDIAN" , MiscMath::median( mean_adj_vals ) );
	      
	      // writer.value( "ADJMEDIAN_MEAN" , MiscMath::mean( median_adj_vals ) );
	      // writer.value( "ADJMEDIAN_MEDIAN" , MiscMath::median( median_adj_vals ) );


	      std::map<double,double>::const_iterator tt = tvals.begin();
	      while ( tt != tvals.end() ) 
		{
		  writer.level( tt->first , "TH" );
		  writer.value("SIGMAB" , tt->second );
		  ++tt;
		}
	      writer.unlevel( "TH" );
	      
	      
	      // if just estimating thresholds, we should skip the actual spindle detection part
	      // go on to next frequency/channel

	      if ( ! use_empirical_thresholds ) 
		continue;
	      
	    }
	
	  // report thresholds

	  logger << "  detection thresholds (core, flank, max)  = " << multiplicative_threshold << ", "
		    << boundary_threshold;
	  if ( maximal_threshold > 0 ) 	 
	    logger << ", " << maximal_threshold;
	  logger << "x" << "\n";
	  
	  logger << "  core duration threshold (core, min, max) = " << min0_dur_sec << ", " << min_dur_sec << ", " << max_dur_sec << "s" << "\n";	  

	  // Set up threshold values as a matrix;  typically, these will use the same mean, 
	  // and so every value will be identical, but allow for the case where we have
	  // e.g. a local sliding window average instead of the whole-night baseline
	  const int sz = averaged.size();
	  
	  // required core to be a spindle	  
	  std::vector<double> threshold( sz , multiplicative_threshold * mean );

	  // optional upper threshold (i.e. to enable looks a sub-threshold spindles)
	  // i.e. reject spindle if max CWT power is above X*mean
	  // -9 means no maximal threshold
	  
	  std::vector<double> threshold_max( sz , maximal_threshold > 0 ? mean * maximal_threshold : -9 );
      
          // start/stop boundaries
          std::vector<double> threshold2( sz, boundary_threshold * mean );


	  //
	  // Adjust thresholds based on local averages?
	  //
	  
	  if ( epoch_norm_sec ) 
	    {

	      // number of points to average
	      // epoch_norm is in # of seconds
	      // note: this specifies the window size (i.e. +/- 50% of the total range)
	      // so set local=60 to have at least one epoch (30s) each side 
	      
	      int window_points = epoch_norm_sec * Fs[s];
	      
	      if ( window_points % 2 == 0 ) ++window_points;

	      // q? okay to use mean here.
	      const std::vector<double> reaveraged =
		( 0 && use_median ) ?
		MiscMath::median_filter( averaged , window_points ) 
	        : MiscMath::moving_average( averaged , window_points );

	      
	      for (int p = 0 ; p < sz ; p++ )
		{
		  threshold[p] = multiplicative_threshold * reaveraged[p] ;
		  threshold2[p] = boundary_threshold * reaveraged[p] ;
		}
	      
	      if ( maximal_threshold > 0 )
		for (int p = 0 ; p < sz ; p++ )
		  threshold_max[p] = maximal_threshold * reaveraged[p];
	    }
	  
	

	  //
	  // Verbose signal display with thresholds
	  //

	  if ( show_cwt_coeff )
	    {
	      
	      writer.var( "RAWCWT" , "Raw CWT coefficient" );
	      writer.var( "CWT" , "CWT coefficient" );
	      writer.var( "CWT_TH" , "CWT primary threshold" );
	      writer.var( "CWT_TH2" , "CWT secondary threshold" );
	      writer.var( "CWT_THMAX" , "CWT maximum threshold" );
	      
	      int np = cwt.points();
	      int nf = cwt.freqs();      
	      if ( np != np0 ) Helper::halt( "internal problem in cwt()" );
	      
	      for (int ti=0;ti<np;ti++)
		{		  
		  writer.interval( interval_t( (*tp)[ti] , (*tp)[ti] ) );
		  writer.value( "RAWCWT" , cwt.raw_result(fi,ti)  );
		  writer.value( "CWT" , cwt.result(fi,ti) );
		  writer.value( "CWT_TH" , threshold[ti] );
		  writer.value( "CWT_TH2" , threshold2[ti] );
		  writer.value( "CWT_THMAX" , threshold_max[ti] );
		}
	      writer.uninterval();
	    }
	  
 		  	  
	  //
	  // Find above threshold regions
	  //
	  
	  std::vector<interval_t> spindles1;
	  std::vector<int> spindles1_start; // sample-points
	  std::vector<int> spindles1_stop; // sample-points

	  if ( use_zpks )
	    {
	      
	      logger << "  robust detection of local peaks\n"; 

	      // all durations for smoothedZ() in sample-points:

	      // prior window to consider ( in sample-points) 
	      const int lag_sp = Fs[s] * zpks_window_sec ; 
	      const uint64_t min0_dur_sp = min0_dur_sec * Fs[s];
	      const uint64_t min_dur_sp  = min_dur_sec * Fs[s];
	      const uint64_t max_dur_sp  = max_dur_sec * Fs[s];
	      	      
	      const bool ignore_negatives = true; 
	      
	      std::vector<interval_t> spindles0;
	      
	      std::vector<int> pk = MiscMath::smoothedZ( averaged , lag_sp , multiplicative_threshold , zpks_influence , min0_dur_sp , 
							 maximal_threshold > 0 ? maximal_threshold : 0 , 
							 boundary_threshold , min_dur_sp , ignore_negatives , &spindles0 );
	      
	      // check total duration <= max_dur_sp
	      
	      for (int i=0; i<spindles0.size(); i++)
		{
		  const interval_t & s = spindles0[i];
		  if ( s.stop - s.start > max_dur_sp ) continue;
		  
		  // tp: use +1 end encoding
		  spindles1.push_back( interval_t( (*tp)[s.start] , (*tp)[s.stop] + 1LLU  ) );

		  // sp: already uses +1 end encoding
		  spindles1_start.push_back( s.start );
		  spindles1_stop.push_back( s.stop );
		}

	    }
	  else
	    {
	      
	      logger << "  basic " << ( use_median ? "median" : "mean" ) << "-based multiplicative threshold rule\n";

	      int start = 0;
	      int stop = 0;
	      int scnt = 0;
	      
	      uint64_t dt = 1.0/(double)Fs[s] * globals::tp_1sec; // time in tp-units
	  
	      if ( averaged.size() != tp->size() ) Helper::halt( "internal error in cwt()\n" );

	      for (int i=0; i<averaged.size(); i++)
		{
		  
		  if ( averaged[i] > threshold[i] ) 
		    {		  
		      
		      if ( scnt == 0 ) start = i;
		      stop = i+1; 
		      ++scnt;
		    }
		  else
		    {
		      
		      if ( scnt > 0 ) 
			{
			  
			  uint64_t start_tp = (*tp)[ start ];
			  uint64_t stop_tp  = (*tp)[ stop ]; // is 1-past end already
			  uint64_t dur_tp   = stop_tp - start_tp ;
			  
			  // does peak area meet duration requirements?
			  
			  if ( dur_tp > min0_dur_tp && dur_tp < max_dur_tp ) 
			    {
			      
			      // core is identified as a spindle, but now extend 
			      // to define boundaries using a lower threshold
			      
			      // prior
			      int j = start;
			      while ( 1 )
				{
				  --j;
				  if ( j <= 0 ) break;
				  if ( averaged[j] < threshold2[j] ) break;
				  start = j;
				}
			      
			      // after
			      j = stop;
			      while ( 1 )
				{
				  ++j;
				  if ( j >= averaged.size() ) break;
				  if ( averaged[j] < threshold2[j] ) break;
				  stop = j+1; // one past end encoding
				}
			      
			      // re-adjusted start/stop times
			      start_tp = (*tp)[ start ];
			      stop_tp  = (*tp)[ stop ];
			      
			      //
			      // Some final checks on whether we should call a spindle here
			      //
			      
			      bool okay = true;
			      
			      //
			      // check that expanded spindle meets the broader definition
			      //
			      
			      uint64_t dur2_tp = stop_tp - start_tp;
			      
			      if ( dur2_tp < min_dur_tp ) okay = false;
			      if ( dur2_tp > max_dur_tp ) okay = false;			  
			      
			      //
			      // check for any max. threshold condition
			      //
			      
			      if ( maximal_threshold > 0 )			  
				{
				  for ( int j = start ; j <= stop ; j++ )
				    {
				      if ( averaged[j] > threshold_max[j] ) { okay = false; break; }
				    }			      
				}
			      
			  
			      //
			      // save this spindle ?
			      //
			      
			      if ( okay )
				{
				  spindles1.push_back( interval_t( start_tp , stop_tp ) );
				  spindles1_start.push_back( start );
				  spindles1_stop.push_back( stop );
				}
			      
			    }
			  
			  scnt = 0;
			}
		  
		    }
		}
	  
	    }


	  
	  //
	  // Merge rule
	  //
	  
	  // in original implementation: if end of two adjacent spindles
	  // within 1-sec, then discard the second; here, instead, we will
	  // merge;  but still apply the 3sec rule
	  
	  std::vector<spindle_t> spindles;
	  
	  int nspindles_premerge = spindles1.size();

	  if ( spindles1.size() > 0 ) 
	    {
	      
	      bool extending = false;
	      
	      uint64_t previous_start = spindles1[0].start;
	      uint64_t previous_stop  = spindles1[0].stop;	  
	      
	      int previous_start_sp = spindles1_start[0];
	      int previous_stop_sp = spindles1_stop[0];

	      
	      for (int i=1;i<spindles1.size();i++)
		{
		  
		  uint64_t this_start = spindles1[i].start;
		  uint64_t this_stop  = spindles1[i].stop;
		  
		  int this_start_sp   = spindles1_start[i];
		  int this_stop_sp   = spindles1_stop[i];
		  
		  // merge? if both ends are within one second (or 
		  // if the first spindle ends 
		  
		  // overlap?
		  if ( this_start < previous_stop ) { extending = true; }

		  // too near?
		  else if ( this_start - previous_stop < spindle_merge_tp ) { extending = true; }
		  
		  // this next spindle is sufficiently far away, so add the previous one
		  else 
		    {				  
		      // does it still meet max duration criterion (i.e. if extended)?
		      if ( previous_stop - previous_start + 1 < max_dur_tp )
			{
			  
			  spindles.push_back( spindle_t( previous_start , previous_stop , 
							 previous_start_sp , previous_stop_sp ) );
			  
			}
		      extending = false;
		    }
		  
		  // update what was last seen
		  
		  if ( ! extending ) 
		    {
		      previous_start = this_start;
		      previous_start_sp = this_start_sp;
		    }

		  previous_stop  = this_stop;
		  previous_stop_sp  = this_stop_sp;
		  
		}
	      
	      
	      // 
	      // and finally add the last spindle (this may either be
	      // an extended one, or not, should not matter
	      //
	      
	      if ( previous_stop - previous_start + 1 < max_dur_tp ) 
		{
		  
		  spindles.push_back( spindle_t( previous_start , previous_stop , 
						previous_start_sp , previous_stop_sp ) );
		  
// 		  spindles.push_back( interval_t( previous_start , previous_stop ) );
// 		  spindles_start.push_back( previous_start_sp );
// 		  spindles_stop.push_back( previous_stop_sp );
		}

	    } 


	  int nspindles_postmerge = spindles.size();

	  logger << "  merged nearby intervals: from " 
		 << spindles1.size() << " to " 
		 << spindles.size() << " unique events\n";
	  

	  
	  //
	  // Track whether each SP is in a spindle or not
	  // Currently, only used when 'if' option
	  //
	  
	  std::vector<bool> in_spindle;
	  
	  if ( ht_chirp ) 
	    {
	      in_spindle.resize( averaged.size() , false );
	      for (int i=0;i<spindles.size();i++) 
		for (int j=spindles[i].start_sp; j <= spindles[i].stop_sp; j++) 
		  in_spindle[j] = true;
	    }
	  
	  


	  //
	  // Characterisation (and display) of each spindle
	  //

	  bool bandpass_filtered_status = false;
	  
	  double total_duration = 0;
	  
	  bool some_data = ! interval.empty();
	  

	  //
	  // Threshold-corrected CWT amplitude (1=at threshold)
	  //
	  
	  std::vector<double> averaged_corr = averaged;
	  for ( int i=0;i<averaged_corr.size();i++) averaged_corr[i] /= threshold[i];
	  


	  //
	  // Track some CH/F level output (i.e. so can all be sent together, given new -t demands...)
	  //
	  
	  
 	  std::map<std::string,double> means;



	  

	  //
	  // Calculate additional spindle parameters, and final spindle-level QC (Q scores, PASS/ENRICH)
	  //	  
	  
	  std::map<double,double> locked;

	  if ( characterize && some_data ) 
	    {	  
	      
	      const double window_f = 4;  // +/- 2 Hz around each peak
	      
	      characterize_spindles( edf , param , signals(s) , 
				     bandpass_filtered_status , 				     
				     frq[fi] , window_f , 
				     "wavelet-" + Helper::dbl2str(frq[fi]) , 
				     &averaged_corr ,    // pass as input threshold-normed CWT
				     d ,                 // original EEG signal 
				     &spindles ,         // this will be annotated/reduced   
				     ( hms ? &starttime : NULL) , 
				     &baseline_fft, 
				     tlocking ? &locked : NULL      // mean signal around spindle troughs
				     );
	      
	      
	    }


	  //
	  // Get mean spindle parameters for this channel/frequency 
	  //
	  	  
	  if ( characterize ) 
	    spindle_stats( spindles , means );
	  
	  //
	  // If we're doing CHIRP analyses, get the mean observed spindle frequency (which we will
	  // use of the basis for the filter-Hilbert characterisation of IF (note, this itself is still
	  // based on some BPF around the target frequency, but not too narrow
	  //
	  
	  double observed_frq = means["FRQ"];
	  
	  logger << "  observed spindle frequency is " << observed_frq << "\n";


	  //
	  // Save spindle locations and CWT coefficients for subsequent propagation analyses?
	  //
	  
	  if ( do_prop ) 
	    {
	      // only adds time-points once -- but will check that size matches on 
	      // subsequent goes around
	      props.add_tp( *tp );
	      
	      // add the actuall dta
	      props.add( frq[fi] , signals.label(s) , spindles , averaged_corr );
	    }
	  
	  

	  
	  //
	  // Optionally, transform of spindle frequencies (+/- 2 Hz ) to get IF
	  //
	  
	  hilbert_t * p_chirp_hilbert = NULL ;
	  std::vector<double> * p_chirp_if = NULL;
	  std::vector<int> * p_chirp_bin = NULL;

	  if ( ht_chirp ) 
	    {
	      // filter-Hilbert raw signal for spindle frequencies

	      // get IF from this;   given noise and issues w/ this approach, 
	      // for now just apply a simple window around the target frequency
	      // and ignore estimates outside of that range...
	      
	      // if BPF has transition frequencies  F-H and F+H
	      // then filter IF on F-2H and F+2H...
		
	      // default BPF window is +/- 2Hz 
	      // 'or' if if-frq2 is set, then use window if-frq .. if-frq2
	      // 'or' if if-frq-emp=X then use observed spindle freq +/- X

	      // set a broad transition widwth 
	      double ripple = 0.02;
	      double tw = 2; 
	      
	      if ( ht_chirp_frq_emp > 0 ) 
		p_chirp_hilbert = new hilbert_t( *d , Fs[s] , observed_frq - ht_chirp_frq_emp  , observed_frq + ht_chirp_frq_emp , ripple , tw );
	      else if ( ht_chirp_frq2 > 0 ) 
		p_chirp_hilbert = new hilbert_t( *d , Fs[s] , ht_chirp_frq  , ht_chirp_frq2 , ripple , tw );
	      else
		p_chirp_hilbert = new hilbert_t( *d , Fs[s] , frq[fi] - ht_chirp_frq  , frq[fi] + ht_chirp_frq , ripple , tw );
	      
	      p_chirp_if = new std::vector<double>;
	      *p_chirp_if = p_chirp_hilbert->instantaneous_frequency( Fs[s] );
	      
	      p_chirp_bin = new std::vector<int>;
	      p_chirp_bin->resize( d->size() , -1 );
	      
	      std::vector<double> isf( ht_bins , 0 );
	      std::vector<int> isfn( ht_bins , 0 );
	      
	      const int nspindles = spindles.size();
	      
	      double ht_lwr =  frq[fi] - 2 * ht_chirp_frq;
	      double ht_upr =  frq[fi] + 2 * ht_chirp_frq;
	      
	      // or use empirical? (add 1 Hz window around here)
	      if ( ht_chirp_frq_emp > 0 ) 
		{
		  ht_lwr = observed_frq - ht_chirp_frq_emp - 1;
                  ht_upr = observed_frq + ht_chirp_frq_emp + 1;
		}
	      // or swap in a fixed band?
	      else if ( ht_chirp_frq2 > 0 ) 
		{
		  ht_lwr = ht_chirp_frq;
		  ht_upr = ht_chirp_frq2;
		}
	      
	      for (int i=0;i<nspindles;i++)
		{
		  int b0 = spindles[i].start_sp;
		  int b1 = spindles[i].stop_sp;
		  double if_spindle = 0;

		  int if_n = 0;
		  int slower_points = 0;

		  std::vector<double> xdata;
		  		  
		  // b0, b1 are start/stop the sample points
		  double denom = ( b1-b0+1 ) / (double)ht_bins;
		  for (int j = b0 ; j <= b1 ; j++ ) 
		    {
		      
		      // is this okay to use? 
		      if ( (*p_chirp_if)[j] >= ht_lwr && (*p_chirp_if)[j] <= ht_upr )
			{
			  // round( prop / (1/B) ) for B bins
			  // so last point (prop==1.0) will be 'over' so we need to 
			  // add to last bin
			  const int idx = j==b1 ? ht_bins-1 : (j-b0) /(double)denom;
		
			  if_spindle += (*p_chirp_if)[j];
			  if_n++;
			  
			  isf[ idx ] += (*p_chirp_if)[j];
			  isfn[ idx ]++;		      
			  (*p_chirp_bin)[j] = idx; // store for use in sw/chirp HT analysis below
			  
			  xdata.push_back( (*p_chirp_if)[j] );
			  
			}
		      else
			{
			  // mark as if not 'in a spindle' so it will not be used for the recording-wide stats below
			  in_spindle[j] = false;
			  // count positive but slower number of sample points
			  if ( (*p_chirp_if)[j] > 0 ) ++slower_points;
			}
		    }
		  
		  //
		  // Mean IF for this spindle
		  //
		  
		  spindles[i].if_spindle = if_spindle;
		  
		}
	      
	      // output by bin
	      for (int j=0;j<ht_bins;j++)
		{
		  writer.level( j+1  , "RELLOC" );
		  writer.value( "IF" , isf[j] / (double)isfn[j] );		  
		}
	      writer.unlevel( "RELLOC" );

	    }

	
	  //
	  // Cache spindle infor
	  //

	  if ( cache_data )
	    {
	      cache_t<double> * cache_num = edf.timeline.cache.find_num( cache_name );
	      cache_num->add( ckey_t( "spindle-wavelet-power" , writer.faclvl() ) , averaged_corr );
	      
	      // cache_t<uint64_t> * cache_tp = edf.timeline.cache_tp( cache_name );
	      // cache_tp->add( ckey_t( "spindle-peaks" , writer.faclvl() ) , averaged_corr );	      

	    }


	  //
	  // Cache spindle peaks
	  //
	  
	  if ( cache_peaks )
            {
              cache_t<int> * cache = edf.timeline.cache.find_int( cache_peaks_name );

	      std::vector<int> peaks;
	      for (int i=0; i<spindles.size(); i++)
		{
		  int p = spindles[i].start_sp + spindles[i].peak_sp;
		  //		  std::cout << "p = " << spindles[i].start_sp  << " " << spindles[i].peak_sp << "\n";
		  peaks.push_back(p);
		}
              cache->add( ckey_t( "points" , writer.faclvl() ) , peaks );
	    }

	  
	  //
	  // Optional slow-wave coupling?
	  //
	  
	  if ( sw_coupling )
	    {

	      // slow waves have already been detected for this channel (in 'so')

	      std::vector<double> ph_peak;
	      
	      // are spindles in slow-waves?

 	      std::vector<bool> sw_peak;

	      if ( verbose_time_phase_locking )
		{

		  //
		  // Phase-locked average of spindle power w.r.t. SO phase
		  //	      
		  
		  // 36 bins = 10-degree bins; 18 = 20-deg bins
		  int nbins = 36/2;
		  
		  std::vector<double> pl_spindle = p_sw->phase_locked_averaging( &averaged_corr , nbins );
		  
		  if ( pl_spindle.size() > 0 ) 
		    {
		      
		      writer.var( "SOPL_CWT" , "Slow wave phase-locked average spindle power" );
		      
		      double inc = 360 / (double)nbins;
		      double ph = inc/2.0; // use mid-point of range
		      
		      for (int j=0;j<nbins;j++)
			{
			  writer.level( ph  , "PHASE" );
			  writer.value( "SOPL_CWT" , pl_spindle[j] );
			  ph += inc;
			}
		      writer.unlevel( "PHASE" );
		    }
		  
	      
		  //
		  // Time-locked averaging
		  //
		  
		  // +1/-1 1 second
		  std::vector<double> tl_spindle = p_sw->time_locked_averaging( &averaged_corr , Fs[s] , 1 , 1 );
		  
		  writer.var( "SOTL_CWT" , "Slow wave time-locked average spindle power" );
		  
		  int sz = tl_spindle.size();
		  
		  if ( sz > 0 )
		    {
		      int sz2 = - (sz-1)/2;
		      
		      for (int j=0;j<sz;j++)
			{
			  writer.level( sz2 , "SP" );
			  writer.value( "SOTL_CWT" , tl_spindle[j] );
			  ++sz2;
			}
		      writer.unlevel( "SP" );
		    }
		  
		}
	      
	      //
	      //  Seed on spindles, consider SO phase
	      //
	      
	      const int nspindles = spindles.size();
	      
	      std::vector<int> sw_spindles_start; std::map<int,int> swmap_start;
	      std::vector<int> sw_spindles_stop; std::map<int,int> swmap_stop;
	      std::vector<int> sw_spindles_peak; std::map<int,int> swmap_peak;
	      std::vector<int> all_spindles_peak;
	      std::vector<double> nearest_sw;
	      std::vector<int> nearest_sw_number;
	      std::vector<uint64_t> spindle_peak;


	      int sw_spin_count = 0;
	      for (int i=0;i<nspindles;i++)
		{
		  int b0 = spindles[i].start_sp;
		  int b1 = spindles[i].stop_sp;
		  double mx = averaged[b0];
		  int mxi = b0;
		  for (int j = b0+1 ; j <= b1 ; j++ )
		    if ( averaged[j] > mx ) { mx = averaged[j] ; mxi = j; } 
		  
		  bool any = false;
		  // is this feature in a slow-wave?
		  if ( p_sw->in_slow_wave(b0) ) 
		    { 
		      any=true; 
		      sw_spindles_start.push_back( b0 ); 
		      swmap_start[i]=sw_spindles_start.size()-1; 
		    }
		  
		  if ( p_sw->in_slow_wave(b1) ) 
		    { 
		      any=true; 
		      sw_spindles_stop.push_back( b1 ); 
		      swmap_stop[i]=sw_spindles_stop.size()-1; 
		    }
		  
		  // is spindle peak in SO ? 
		  if ( p_sw->in_slow_wave(mxi) ) 
		    { any=true;
		      sw_spindles_peak.push_back( mxi ); 
		      swmap_peak[i]=sw_spindles_peak.size()-1; 
		    }
		  
		  if ( any ) ++sw_spin_count;
		  
		  // record all peaks 
		  all_spindles_peak.push_back( mxi ); 
		  
		  // second distance to nearest SW (secs)
		  int sw_num = 0;
		  nearest_sw.push_back( p_sw->nearest(mxi , &sw_num ) ); 
		  nearest_sw_number.push_back( sw_num + 1 ); // make 1-based for visual output
		  spindle_peak.push_back( (*tp)[ mxi ] );
		}


	      //
	      // Proportion of spindles in a SO
	      //
	      
	      std::string analysis_label = "wavelet-" + Helper::dbl2str(frq[fi]) ;
	      
	     
	      //
	      // SW-phase for each spindle in a SO
	      //
	      
	      if ( all_spindles_peak.size() > 0 ) 
		{
		  
		  //
		  // restrict SO/spindle coupling calculations to
		  // spindles that occur within a SO? Or use all
		  // spindles?
		  //

		  std::vector<bool> so_mask;		  

		  // default is to use mask
		  bool use_mask = ! param.has( "all-spindles" );
		  
		  if ( use_mask ) 
		    so_mask = p_sw->sp_in_sw_vec();

		  //
		  // report coupling overlap by SO phase
		  //
		  
		  bool stratify_by_so_phase_bin = param.has ( "stratify-by-phase" );
		  
		  //
		  // Within-epoch permutation (default)
		  //

		  bool eperm = ! param.has( "perm-whole-trace" );
		  
		  double epoch_sec = 0 ;
		  int sr = Fs[s];
		
		  if ( eperm ) 
		    {
		      if ( ! edf.timeline.epoched() ) edf.timeline.ensure_epoched();
		      epoch_sec = edf.timeline.epoch_length();
		    }

		  // use permutation for ITPC values?

		  int nreps = 0;

		  if ( param.has( "nreps" ) ) nreps = param.requires_int( "nreps" );

		  if ( nreps != 0 && nreps < 10 ) Helper::halt( "nreps must be 10+" );
		  

		  //
		  // Perform spindle/SO coupling analysis
		  //

		  itpc_t itpc  = p_hilbert->phase_events( all_spindles_peak , 
							  use_mask ? &so_mask : NULL , nreps ,  // optional, mask events/spindles
							  sr , 
							  epoch_sec ,  // optionally, within-epoch shuffle
							  stratify_by_so_phase_bin // optional, overlap by SO-phase bin
							  );
		  
		  sw_peak = itpc.event_included;
		  
		  ph_peak = itpc.phase;
		  
		  //
		  // Gather output (but don't send to writer until later, i.e. need to 
		  // group all CH/F strata output together for -t mode
		  //

		  // ITPC magnitude of coupling
		  
		  means[ "COUPL_MAG" ] =  itpc.itpc.obs ;
		  //writer.value( "COUPL_MAG"      , itpc.itpc.obs );

		  if ( use_mask )
		    means[ "COUPL_OVERLAP" ] = itpc.ninc.obs ;
		  
		  if ( nreps ) 
		    {
		      means[  "COUPL_MAG_EMP"  ] = itpc.itpc.p ;
		      means[  "COUPL_MAG_NULL" ] = itpc.itpc.mean ;
		      means[  "COUPL_MAG_Z"    ] = ( itpc.itpc.obs - itpc.itpc.mean ) / itpc.itpc.sd ;
		      
		      // proportion of spdinles that overlap a SO 
		      // unless itpc-so was set, this will be meaningless
		      // so only report is a 'mask' was set		  
		      
		      if ( use_mask ) 
			{
			  			  
			  means[ "COUPL_OVERLAP_EMP" ] = itpc.ninc.p ;

			  means[ "COUPL_OVERLAP_NULL" ] = itpc.ninc.mean ;
			  
			  if ( itpc.ninc.sd > 0 )  		    			
			    means[ "COUPL_OVERLAP_Z" ] = ( itpc.ninc.obs - itpc.ninc.mean ) / itpc.ninc.sd ;
			}
		      
		    }

		  //
		  // mean angle; no empirical test results; -9 means no events observed, so set to missing
		  //
		  
		  if ( itpc.angle.obs > -9 ) 
		    means[ "COUPL_ANGLE" ] = itpc.angle.obs ;
		  
		  //
		  // asymptotic significance of coupling test; under
		  // the null, give mean rate of 'significant'
		  // (P<0.05) coupling
		  //
		  
		  means[ "COUPL_PV" ] = itpc.pv.obs ;
		  
		  if ( nreps ) 
		    means[ "COUPL_SIGPV_NULL" ] = itpc.sig.mean ; 

		  
		  //
		  // phase-bin stratified overlap/counts
		  //
		  
		  if ( nreps && 
		       stratify_by_so_phase_bin )
		    {
		      
		      const int nbins = 18;
		      for (int b = 0 ; b < nbins ; b++ ) 
			{
			  writer.level( b * 20 + 10 , "PHASE" );
			  writer.value( "COUPL_OVERLAP"      , itpc.phasebin[b].obs );
			  writer.value( "COUPL_OVERLAP_EMP"    , itpc.phasebin[b].p );
			  //writer.value( "COUPL_OVERLAP_NULL" , itpc.phasebin[b].mean );
			  if ( itpc.phasebin[b].sd > 0 ) 
			    {
			      double Z = ( itpc.phasebin[b].obs - itpc.phasebin[b].mean ) / itpc.phasebin[b].sd ; 
			      writer.value( "COUPL_OVERLAP_Z" , Z );
			    }
			}
		      writer.unlevel( "PHASE" );
		    }

		}
	        
	  
	      //
	      // Individual PEAKS, to be output later, in characterize() 
	      //
	      
	      for (int i=0;i<nspindles;i++)
		{
		  // peak_sec
		  spindles[i].peak_sec = spindle_peak[i] * globals::tp_duration ;
		  
		  if ( nearest_sw_number[i] != 0 ) // is now 1-basewd 
		    {
		      spindles[i].so_nearest = nearest_sw[i] ;
		      spindles[i].so_nearest_num = nearest_sw_number[i] ;		      
		    }
		  
		  if ( sw_peak[ i ] )
		    spindles[i].so_phase_peak =  MiscMath::as_angle_0_pos2neg( ph_peak[ i ] ) ;		  		  		      
		  else 
		    spindles[i].so_phase_peak =  -9;
		  
		}
	      
	    	     	      	      
	      //
	      // Optional, consideration of spindle chirp as a function of SO phase
	      //
	      
	      if ( ht_chirp ) 
		{
		  

		  // look at IF in spindle range, i.e. but not specifically detecting individual spindles
		  // i.e. parallel to phase and time locked SW analyses above

		  int nbins = 36/2;
	      
		  
		  //
		  // Look at spindle IF as a function of SW phase 		  
		  //

		  std::vector<double> pl_chirp = p_sw->phase_locked_averaging( p_chirp_if , nbins , &in_spindle );
	      	      	      	      
		  double inc = 360 / (double)nbins;
		  double ph = inc/2.0; // use mid-point of range
		  
		  for (int j=0;j<nbins;j++)
		    {
		      writer.level( ph  , "PHASE" );
		      writer.value( "IF" , pl_chirp[j] );
		      ph += inc;
		    }
		  writer.unlevel( "PHASE" );
		  
		  // time-locked SO spindle IF -- code not used -- phase-locked analysis above should 
		  // be sufficient

		  if ( 0 ) 
		    {
		      // +1/-1 1 second
		      std::vector<double> tl_chirp = p_sw->time_locked_averaging( p_chirp_if , Fs[s] , 1 , 1 );
		      
		      //		      writer.var( "IF" , "Slow wave time-locked average spindle frequency/chirp" );
		      
		      int sz = tl_chirp.size();
		      if ( sz > 0 ) 
			{
			  int sz2 = - (sz-1)/2;
			  
			  for (int j=0;j<sz;j++)
			    {
			      writer.level( sz2 , "SP" );
			      writer.value( "IF" , tl_chirp[j] );
			      ++sz2;
			    }
			  writer.unlevel( "SP" );
			}
		    }

		  
		  //
		  // Spindle IF as a function of slow wave phase and position in the spindle
		  //
		  
		  for (int h=0;h<ht_bins;h++)
		    {

		      writer.level( h+1  , "RELLOC" );
		      
		      // black out 
		      std::vector<bool> in_spindle_and_bin = in_spindle;
		      for (int i=0; i < in_spindle_and_bin.size(); i++) 
			if ( (*p_chirp_bin)[i] != h ) in_spindle_and_bin[i] = false;
		      
		      int ccc = 0;
		      for (int i=0;i<in_spindle_and_bin.size(); i++) if ( in_spindle_and_bin[i] ) ++ccc;
		      //logger << "h etc = " << h << " " << in_spindle_and_bin.size() << " " << ccc << "\n";
		      
		      std::vector<double> pl_chirp = p_sw->phase_locked_averaging( p_chirp_if , nbins , &in_spindle_and_bin );
		      
		      double inc = 360 / (double)nbins;
		      double ph = inc/2.0; // use mid-point of range
		      
		      for (int j=0;j<nbins;j++)
			{
			  writer.level( ph  , "PHASE" );
			  writer.value( "IF" , pl_chirp[j] );
			  ph += inc;
			}
		      writer.unlevel( "PHASE" );
		    }
		  writer.unlevel( "RELLOC" );
		  
		}
	    
	    
	    } // end of SW-coupling code
	  
	  

	  //
	  // Per-spindle level output
	  //

	  if ( show_spindle_level )
	    per_spindle_output( &spindles , param , ( hms ? &starttime : NULL) , &baseline_fft );



	  //
	  // plot spindles?  
	  //
	  
	  if ( spindle_pdf && some_data )
	    {
	      std::string analysis_label = "wavelet-" + Helper::dbl2str(frq[fi]) ;
	      std::string fname = param.value( "pdf" ) + "-" + signals.label(s) + "-" + analysis_label + ".pdf";
	      logger << "  writing PDF of spindle traces to " << fname << "\n";
	      std::map<uint64_t,double> avgmap;
	      for (int j=0;j<averaged.size();j++) avgmap[ (*tp)[j] ] = averaged[j] ;
	      draw_spindles( edf , param , fname , signals(s) , spindles, &avgmap );
	    }

	  
	  //
	  // output: time-locked signal averaging 
	  //
	  
	  if ( characterize & tlocking && spindles.size() > 0 )
	    {
	      
	      double tlock_min = locked.begin()->first;
	      double tlock_max = (--locked.end())->first;
	      
	      writer.var( "TLOCK" , "Average EEG amplitude time-locked to spindle peak" );
	      
	      std::map<double,double>::const_iterator ll = locked.begin();
	      while ( ll != locked.end() )
		{
		  writer.level( ll->first , "MSEC" );
		  writer.value( "TLOCK" , ll->second );
		  ++ll;
		}
	      writer.unlevel( "MSEC" );

	    }
	

	  //
	  // Estiamte of spindle density to console
	  //

	  bool empty = spindles.size() == 0; 
	  
	  if ( ! empty )
	    logger << "  estimated spindle density is " << spindles.size() / t_minutes << "\n";
  


	  //
	  // Save for an 'intersection' command?
	  //

	  if ( do_merge &&  spindles.size() > 0 )
	    {
	    
	      // channel specific info? 

	      if ( do_channel_merge )
		{

		  // on first instance, assign same parameters and pointer to edf_t
		  if ( ch2mspindles.find( signals.label(s) ) == ch2mspindles.end() )
		    ch2mspindles[ signals.label(s) ] = mspindles; 
		  
		  ch2mspindles[ signals.label(s) ].add( spindles , 
							Fs[s] , averaged.size() , 
							frq[fi] , 		
							signals(s) , 
							signals.label(s) + ":" + Helper::dbl2str( frq[fi] ) ); 
		  
		}
	      else
		{
		  
		  mspindles.add( spindles , 
				 Fs[s] , averaged.size() , 
				 frq[fi] , 		
				 signals(s) , 
				 signals.label(s) + ":" + Helper::dbl2str( frq[fi] ) );
		  
		}

	    }


	 	  

	  //
	  // Per-EPOCH summary and test of over-dispersion
	  //
	  

	  edf.timeline.first_epoch();
	      
	  std::vector<int> epoch_counts;
	  
	  while ( 1 ) 
	    {
	      
	      int epoch = edf.timeline.next_epoch();      
	      
	      if ( epoch == -1 ) break;
              
	      interval_t interval = edf.timeline.epoch( epoch );
	      
	      
	      const int nsp = spindles.size();
	      
	      int sp_epoch = 0;
	      
	      for (int i=0 ; i<nsp; i++)
		{
		  
		  // is this spindle included?
		  
		  if ( ! spindles[i].include ) continue;
		  
		  // dummy interval for just starting point of spindle 
		  // i.e. for purpose of assigning to EPOCH (so as not
		  // to double-count spindles that overlap epoch boundaries
		  
		  interval_t spstart( spindles[i].tp.start , spindles[i].tp.start );
		  
		  if ( interval.overlaps( spstart ) )
		    ++sp_epoch;
		  else if ( spstart.is_after( interval ) )
		    break; // spindles are in order, so can skip
		}
	      
	      // record
	      epoch_counts.push_back( sp_epoch );
	      

	      //
	      // per-epoch output
	      //

	      if ( show_epoch_level )
		{	      
		  writer.epoch( edf.timeline.display_epoch( epoch ) );		  
		  // per-epoch spindle count
		  writer.value( "N" , sp_epoch );
		}
	      
	    }


	  // close out any epoch-level output
	  if ( show_epoch_level )
	    writer.unepoch();
	  
	      
	  // 
	  // Test for over-dispersion of spindle counts
	  //
	  
	  if ( ! empty )
	    {
	      double pval = 0;
	      double stat = MiscMath::overdispersion( epoch_counts , &pval );
	      
	      writer.var( "DSPERSION" ,"Spindle epoch-dispersion index" );
	      writer.var( "DISPERSION_P" , "Spindle epoch-dispersion index p-value" );
	      writer.var( "NE" , "Number of epochs for spindle detection" );
	      
	      means[ "DISPERSION" ] = stat ;
	      means[ "DISPERSION_P" ] = pval ;
	      means[ "NE" ] = (int)epoch_counts.size() ;
	    }
	  
	      
	  

	  //
	  // Main output
	  //

	  
	  //
	  // Output over all epochs
	  //
	  
	  if ( characterize && ! empty )
	    {
	      
	      writer.var( "N01"     , "Number of spindles prior to merging" );
	      writer.var( "N02"     , "Number of spindles prior to QC" );
	      writer.var( "N"       , "Final number of spindles" );
	      
	      writer.var( "MINS"  , "Number of minutes for spindle detection" );	      
	      writer.var( "DENS"  , "Spindle density (per minute)" );
	      writer.var( "AMP"   , "Mean spindle amplitude" );
	      writer.var( "DUR"   , "Mean spindle duration" );
	      writer.var( "FWHM"  , "Mean spindle FWHM" );
	      writer.var( "NOSC"  , "Mean spindle number of oscillations" );
	      writer.var( "FRQ"   , "Mean spindle frequency (zero-crossing method)" );
	      writer.var( "FFT"   , "Mean spindle frequency (FFT)" );	      
	      writer.var( "SYMM"  , "Mean spindle symmetry index" );
	      writer.var( "SYMM2"  , "Mean spindle folded symmetry index" );
	      writer.var( "CHIRP" , "Mean spindle chirp index" );
	      
	      
	      writer.value( "N01" , nspindles_premerge );  // original
	      writer.value( "N02" , nspindles_postmerge ); // post merging
	      writer.value( "N" ,  (int)spindles.size()  ) ;    // post merging and QC	    
	      writer.value( "MINS" , t_minutes );
	      writer.value( "DENS" , spindles.size() / t_minutes );
	      
	      writer.value( "ISA_S" , means[ "ISA_PER_SPINDLE" ] );
	      writer.value( "ISA_M" , means[ "ISA_TOTAL" ] / t_minutes );
	      writer.value( "ISA_T" , means[ "ISA_TOTAL" ] );
	      writer.value( "Q"     , means[ "Q" ] );

	      writer.value( "AMP" , means["AMP"] );
	      writer.value( "DUR" , means["DUR"] );
	      writer.value( "FWHM" , means["FWHM"] );
	      writer.value( "NOSC" , means["NOSC"] );
	      writer.value( "FRQ" , means["FRQ"] );
	      writer.value( "FFT" , means["FFT"] );
	      writer.value( "SYMM" , means["SYMM"] );
	      writer.value( "SYMM2" , means["SYMM2"] );
	      writer.value( "CHIRP" , means["CHIRP"] );
	      writer.value( "CHIRPF" , means["CHIRPF"] );

	      if ( globals::devel )
		{
		  writer.value( "F_POS" , means["FPOS"] );
		  writer.value( "F_NEG" , means["FNEG"] );
		  writer.value( "F_ALL" , means["FALL"] );
		  writer.value( "F_DIF" , means["FPOS"] - means["FNEG"] ); 

		  writer.value( "B_POS" , means["BPOS"] );
		  writer.value( "B_NEG" , means["BNEG"] );
		  writer.value( "B_ALL" , means["BALL"] );
		  writer.value( "B_DIF" , means["BPOS"] - means["BNEG"] ); 		  

		  writer.value( "V_POS" , means["VPOS"] );
		  writer.value( "V_NEG" , means["VNEG"] );
		  writer.value( "V_ALL" , means["VALL"] );
		  writer.value( "V_DIF" , means["VPOS"] - means["VNEG"]  );
		}
	      
	      // cache main metrics also?
	      if ( cache_metrics )
		{
		  std::map<std::string,std::string> faclvl = writer.faclvl() ;
		  cache_metrics->add( ckey_t( "DENS" ,  faclvl ) , spindles.size() / t_minutes );
		  cache_metrics->add( ckey_t( "AMP" ,   faclvl ) , means["AMP"]  );
		  cache_metrics->add( ckey_t( "DUR" ,   faclvl ) , means["DUR"]  );
		  cache_metrics->add( ckey_t( "ISA_S" , faclvl ) , means["ISA_PER_SPINDLE"] );
		  cache_metrics->add( ckey_t( "CHIRP" , faclvl ) , means["CHIRP"]  );	      
		}

	      writer.value( "DISPERSION" , means[ "DISPERSION" ] );
	      writer.value( "DISPERSION_P" , means[ "DISPERSION_P" ] );
	      writer.value( "NE" , means[ "NE" ] );

	      write_if_exists( "COUPL_MAG" , means );
	      write_if_exists( "COUPL_MAG_EMP" , means );
	      write_if_exists( "COUPL_MAG_Z" , means );
	      write_if_exists( "COUPL_MAG_NULL" , means );

	      write_if_exists( "COUPL_OVERLAP" , means );
	      write_if_exists( "COUPL_OVERLAP_EMP" , means );
	      write_if_exists( "COUPL_OVERLAP_NULL" , means );
	      write_if_exists( "COUPL_OVERLAP_Z" , means );

	      write_if_exists( "COUPL_ANGLE" , means );
	      write_if_exists( "COUPL_PV" , means );

	      write_if_exists( "COUPL_SIGPV_NULL" , means );

	    }




	  
	  //
	  // Adding new signals?
	  //
	  
	  if ( add_channels )
	    {
	      
	      const int n1 = spindles.size();
	      
	      // make 0/1 for spindle call 
	      std::vector<double> is_spindle( averaged.size() , 0 );
	      for (int i=0;i<n1;i++)
		{
		  int start = spindles[i].start_sp;
		  int stop  = spindles[i].stop_sp;		  
		  for (int j=start;j<=stop;j++) is_spindle[j] = 1.0;
		}
	      
	      // normalize CWT values
	      std::vector<double> zt = MiscMath::Z( results );
	      
	      edf.add_signal( "CWT-raw-" + Helper::dbl2str(frq[fi]) + "-" + signals.label(s) , Fs[s] , zt );
	      
	      // show averaged values, but only above threshold
 	      std::vector<double> copy = averaged;
 	      for (int i=0;i<copy.size();i++) 
		{
		  copy[i] /= threshold[i] ;
		  if ( copy[i] < 1 ) copy[i] = 0; 
		}
 	      edf.add_signal( "CWT-avg-" + Helper::dbl2str(frq[fi]) + "-" + signals.label(s) , Fs[s] , copy );	      
	      
	      edf.add_signal( "spindle-" + Helper::dbl2str(frq[fi]) + "-" + signals.label(s) , Fs[s] , is_spindle );
	      
	    }

	
	  //
	  // Record as an .annot file?
	  //

	  if ( save_annots )
	    {

	      // annot label
	      const std::string analysis_label = Helper::dbl2str(frq[fi]) ;

	      const std::string aname = sp_label;// + "-" + analysis_label;
		
	      annot_t * a = edf.timeline.annotations.add( aname );
	      a->description = "Spindle intervals";
	      
	      logger << "  creating annotation class: " << aname 
		     << ", instance: " << analysis_label 
		     << ", channel: " << signals.label(s) << "\n";

	      // use F_C as instance label
	      for (int i=0;i<spindles.size();i++)
		instance_t * instance = a->add( analysis_label , spindles[i].tp , signals.label(s) );
	      
	    }
	  

	  //
	  // Clean-up at the spindle/F level
	  //
	  

	  if ( ht_chirp ) 
	    {
	      delete p_chirp_hilbert;
	      delete p_chirp_if;
	      delete p_chirp_bin;
	      
	      p_chirp_hilbert = NULL;
	      p_chirp_if = NULL;
	      p_chirp_bin = NULL;
	      
	    }

	  
	  //
	  // Next wavelet frequency
	  //
	  
	  writer.unlevel( globals::freq_strat );
	  
	}
            
      
      //
      // Clean up
      //

      if ( sw_coupling )
	{

	  delete p_sw;
	  delete p_hilbert;
	  
	  p_sw = NULL;
	  p_hilbert = NULL;
	}

      
      //
      // Next signal
      //
      
      writer.unlevel( globals::signal_strat );
      
    }
  
  

  //
  // Spindle propagation analysis?
  //

  if ( do_prop ) 
    {
      
      const double w = 1.0; // 1 sec window

      // all channels
      std::set<std::string> c;
      for (int s = 0 ; s < ns ; s++ )
	{
	  if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
	  c.insert( signals.label(s) );
	}
      
      // freqs one-at-a-time
      for ( int fi=0; fi<frq.size(); fi++ )
	{
	  
	  writer.level( frq[fi] , globals::freq_strat );

	  std::set<double> f;
	  f.insert( frq[fi] );	  
	  
	  std::vector<double> avgs;

	  // take each channel as seed
	  for (int s = 0 ; s < ns ; s++ )
	    {
	      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
	      
	      // do analysis
	      avgs.push_back( props.analyse( f , c , signals.label(s) , w ) );
	      
	      // report average
	    }
	  
	  // scale seed values
	  if ( avgs.size() )
	    {
	      // report both original means, and rescaled versions
	      // i.e. in raw seconds, plus relative to earliest/latest
	      std::vector<double> orig = avgs;

	      double seed_min = avgs[0];
	      double seed_max = avgs[0];
	      for ( int s=0; s<avgs.size(); s++)
		{
		  if ( avgs[s] < seed_min ) seed_min = avgs[s];
		  if ( avgs[s] > seed_max ) seed_max = avgs[s];
		}

	      // scale
	      if ( seed_max - seed_min > 0 )
		for ( int s=0; s<avgs.size(); s++)
		  avgs[s] = ( avgs[s] - seed_min ) / ( seed_max - seed_min );
	      
	      // report
	      int ss = 0;
	      for (int s = 0 ; s < ns ; s++ )
		{
		  if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
		  // track which channel we are seeding on
		  writer.level(  signals.label(s) , "SEED" );
		  writer.value( "T" , orig[ss] ); // nb. ss not s
		  writer.value( "R" , avgs[ss] ); // nb. ss not s
		  ++ss;
		}
	      writer.unlevel( "SEED" );
	    } 
	} // next freq bin      
     
      writer.unlevel( globals::freq_strat );
    }


  //
  // Collation of spindles across any frequencies/channels 
  //

  if ( do_merge )
    {
      
      if ( do_channel_merge )  // stratify by channel
	{
	  
	  std::map<std::string,mspindles_t>::iterator mm = ch2mspindles.begin();
	  while ( mm != ch2mspindles.end() )
	    {

	      // output stratified by channel
	      writer.level( mm->first , globals::signal_strat );

	      mspindles_t & ms = mm->second;	      
	      ms.collate();
	      ms.output( signals );
	      
	      writer.unlevel( globals::signal_strat );

	      // no plots yet...
	      ++mm;
	    }	  
	  
	}

      else // merge across all channels 
	{
	  
	  // collate all
	  mspindles.collate();
	  
	  // some output
	  mspindles.output( signals );
	  
	  // plot merged spindles?
	  mspindles.plot( "mspindles.pdf" );
	}

    }


  //
  // If we added new channels, then we need to save a new EDF
  //

  if ( add_channels ) 
    proc_write( edf , param );
    
  
  
  return NULL;
  
}




void characterize_spindles( edf_t & edf , 
			    param_t & param , 
			    const int s0 ,
			    bool bandpass_filtered , 
			    const double target_f , 
			    const double window_f ,
			    const std::string & analysis_label ,
			    const std::vector<double> * averaged ,
			    const std::vector<double> * original_signal , 
			    std::vector<spindle_t>    * spindles ,
			    clocktime_t               * starttime , 
			    std::map<freq_range_t,double> * baseline ,
			    std::map<double,double> * locked 
			    )   // input, is 
{
 
  //
  // Copy key output modes
  //

  const bool enrich_output = param.has( "enrich" );

  //
  // Create a copy of this signal, if it does not already exist
  //
  
  const std::string signal_label = edf.header.label[s0];
  const std::string new_label = signal_label + "_BP" + "_" + Helper::dbl2str( target_f ) + "_" + Helper::dbl2str( window_f );

  if ( ! edf.header.has_signal( new_label ) )
    {

      // copy the existing signal

      edf.copy_signal( edf.header.label[s0] , new_label );
      
      const int s = edf.header.signal( new_label );

      // and do we need to band-pass filter this new signal?
	
      if ( ! bandpass_filtered ) 
	{
	  
	  logger << "  filtering at " << target_f - window_f * 0.5  << " to " << target_f + window_f * 0.5 << "\n";
	  
	  // default above is 4 Hz, i.e. +/- 2 Hz    ~ 9-13Hz for slow spindles,    13-17Hz for fast spindles

	  //  ripple = 0.02 tw=4 bandpass=9,13 
	  //  ripple = 0.02 tw=4 bandpass=13,17 
	  
	  dsptools::apply_fir( edf , s , fir_t::BAND_PASS ,
			       1 , // 1 = Kaiser window
			       0.02 , 4 , // ripple , transition width
			       target_f - window_f * 0.5 ,
			       target_f + window_f * 0.5 ,
			       0, // order (if not Kaiser win)
			       fir_t::RECTANGULAR , // large window, no window
			       true  // use FFT convolution
			       );
	  
	}
      
    }


  const int s = edf.header.signal( new_label );

   //
   // Output
   //
   
   const int n = spindles->size();

   //
   // Spindle-level QC filters (set default at 0, i.e spindle-activity must be more likely)
   //

   bool qc_q = true;
   double qc_qmin = 0 , qc_qmax = -1;
   if ( param.has( "q" ) ) { qc_q = true; qc_qmin = param.requires_dbl( "q" ); } 
   if ( param.has( "q-max" ) ) { qc_q = true; qc_qmax = param.requires_dbl( "q-max" ); }
   
   //
   // track if we QC any spindles out at this step
   //

   bool removed_some = false;


   //
   // Iterate over each spindle
   //
   
   for (int i=0;i<n;i++)
     {

       spindle_t * spindle = &(*spindles)[ i ];
       
       //
       // duration
       //
       
       spindle->dur = (spindle->tp.stop - spindle->tp.start + 1 ) / (double)globals::tp_1sec;
      
       
       //
       // pull out band-pass filtered data for actual spindle
       //
       
       slice_t slice( edf , s , spindle->tp );
       
       std::vector<double> d = *slice.pdata();

       std::vector<uint64_t> tp = *slice.ptimepoints();
       
       const int Fs = edf.header.sampling_freq( s );
       
       const double period_sec = 1.0/(double)Fs;

       const uint64_t period = period_sec * globals::tp_1sec;

       const int npoints = d.size();
       
      
       //
       // ISA (scale by SR)
       //
       
       spindle->isa = 0;
       
       if ( averaged != NULL )
	 {
	   
	  int start = spindle->start_sp;
	  int stop  = spindle->stop_sp;
	  
	  for (int s=start;s<=stop;s++)
	    spindle->isa += (*averaged)[s] ;
	 }
       
       spindle->isa /= (double)Fs;

       
       //
       // Sanity check
       //
       
       if ( npoints < 2 ) 
	 {
	   // next spindle
	   continue;	   
	 }

       
       
       //
       // FWHM estimate of duration
       //
       
       spindle->fwhm = 0;
       
       if ( averaged != NULL )
	 {
	  
	  int start = spindle->start_sp;
	  int stop  = spindle->stop_sp;
	  
	  // get max (CWT are all positive) 
	  double mx = 0;
	  int mxi = start;
	  
	  for (int s=start;s<=stop;s++)
	    {	      
	      if ( (*averaged)[s] > mx ) 
		{
		  mx = (*averaged)[s] ;
		  mxi = s;
		}
	    }
	  
	  // Move out until we hit 50% drop
	  int lwr = mxi , upr = mxi ; 
	  double half = mx / 2.0;
	  while ( 1 ) // lower
	    {
	      if ( lwr < 0 ) break;
	      if ( (*averaged)[lwr] <= half ) break;
	      --lwr;
	    }
	  
	  while ( 1 ) // upper
	    {
	      if ( upr == averaged->size() ) break;
	      if ( (*averaged)[upr] <= half ) break;
	      ++upr;
	    }
	  
	  // *assumes* a contiguous segment... could be problematic 
	  spindle->fwhm = (1.0/(double)Fs) * ( upr - lwr + 1 );
	  
	}
      

      //
      // Get max/avarage of the statistic
      //
      
      spindle->max_stat = 0;
      spindle->mean_stat = 0;
      
      if ( averaged != NULL )
	{
	  
	  int start = spindle->start_sp;
	  int stop  = spindle->stop_sp;

	  // get max (CWT are all positive) 
	  
	  double sum = 0;
	  
	  for (int s=start;s<=stop;s++)
	    {
	      double x = (*averaged)[s] ; 
	      if ( x > spindle->max_stat ) spindle->max_stat = x ;
	      sum += x ;	      
	    }
	  spindle->mean_stat = sum / (double)(stop - start + 1);
	}
      
      

      //
      // Find largest peak-to-peak amplitude
      //
       
      // rule of thumb: call something a peak if the surrounding +/i 2
      // points are all smaller, or all larger; for tied values, extend to the next 
      // variable ones, and track only the last one as the peak
      
      std::vector<int> peak;
      
      for (int p=2;p<npoints-2;p++)
	{
	  
	  // tied w/ the next point? 
	  if ( d[p] == d[p+1] ) continue;

	  // peaks  

	  int gt = 0 , lt = 0;
	  
	  // forwards
	  if      ( d[p] < d[p+1] ) ++lt;
	  else if ( d[p] > d[p+1] ) ++gt;
	  
	  if      ( d[p] < d[p+2] ) ++lt;
	  else if ( d[p] > d[p+2] ) ++gt;

	  // backwards, which might require we skip ties previously skipped
	  int bck = 0;
	  bool skip = false;
	  while ( 1 ) 
	    {
	      if ( p - bck < 1 ) { skip = true; break; } 
	      if ( d[p] == d[p-bck] ) ++bck;
	      else break;
	    }

	  if ( skip ) continue;
	  
	  if      ( d[p] < d[p-bck] ) ++lt;
	  else if ( d[p] > d[p-bck] ) ++gt;
	  
	  if      ( d[p] < d[p-bck-1] ) ++lt;
	  else if ( d[p] > d[p-bck-1] ) ++gt;
	  
	  if ( gt == 4 ) peak.push_back(p); 
	  else if ( lt == 4 ) peak.push_back(p);
	  
	}
      
      // something strange?  bail
      if ( peak.size() < 2 ) 
	{
	  logger << " *** warning: spindle w/ only a single peak... should not happen... bailing" << "\n";
	  continue;
	}

      //
      // Zero-crossings, in seconds, with linear interpolation between points
      //

      std::vector<double> zc;  // duration of i to i+1 zero-crossing (i.e. half-waves)
      std::vector<bool> zcp;   // T : pos-neg; F : neg-pos

      for (int p=0;p<npoints-1;p++)
        {
	  const bool pos2neg = d[p] >= 0 && d[p+1] < 0 ;
	  const bool neg2pos = d[p] <= 0 && d[p+1] > 0 ;
	  if ( ! ( pos2neg || neg2pos ) ) continue;

	  const double s1 = p * period_sec ;
	  const double frac = fabs( d[p] ) / ( fabs( d[p] ) + fabs( d[p+1] ) );

	  zc.push_back( s1 + frac * period_sec );
	  zcp.push_back( pos2neg );
	}
      
      //
      // pos/neg halfwave freqs given ZCs
      //
      
      spindle->posf = 0;
      spindle->negf = 0;
      int posc = 0 , negc = 0;

      // duration of halfwaves
      std::vector<double> wpos, wneg, wall; 
      std::vector<double> tpos, tneg, tall; 

      for (int z=0; z<zc.size()-1; z++)
	{
	  // transform to frequency 1/SW
	  double w = 1.0 / ( 2 * ( zc[z+1] - zc[z] ) );
	  double t = ( zc[z] + zc[z+1] ) / 2.0;
	  
	  // neg-halfwave
	  if ( zcp[z] )
	    {
	      spindle->negf += w ;
	      ++negc;
	      wneg.push_back(w);
	      tneg.push_back(t);
	    }
	  else
	    {
	      spindle->posf += w ;
	      ++posc;
	      wpos.push_back( w );
	      tpos.push_back( t );
	    }
	  
	  wall.push_back( w );
	  tall.push_back( t );
	}
      
      spindle->posf /= (double)posc;
      spindle->negf /= (double)negc;


      //
      // Slope of frequency implied by ZC / stratified by POS and NEG halfwaves
      //

      // Y = frequencxy implied by halfwave
      // X = time midpoint of bounding ZC (elapsed seconds from spindle start)
      dynam_t zall( wall , tall ); 
      dynam_t zpos( wpos , tpos ); 
      dynam_t zneg( wneg , tneg ); 

      zall.linear_trend( &spindle->allb , NULL );
      zpos.linear_trend( &spindle->posb , NULL );
      zneg.linear_trend( &spindle->negb , NULL );
      
      zall.mean_variance( &spindle->allf , &spindle->allv );
      zpos.mean_variance( &spindle->posf , &spindle->posv );
      zneg.mean_variance( &spindle->negf , &spindle->negv );

      
      //
      // Simple spindle 'chirp' metrics
      //  - contrast of first vs second half of spindle
      //  - based on peak-to-peak durations (both pos + neg, so F = 1/2T)
      //
      
      double int1 = 0 , int2 = 0;   // mean duration (in sample-points) peaks within each half
      int   cint1 = 0 , cint2 = 0;  // number of peak-to-peak intervals in each half (pos-neg and neg-pos)
      
      for (int pi = 0 ; pi < peak.size() ; pi++)
	{
	  // simple first/second half chirp
	  double pos = peak[pi] / (double)(npoints-1);
	  if ( pos < 0.5 ) 
	    { 
	      if ( pi > 0 ) { int1 += peak[pi] - peak[pi-1]; cint1++; } 
	    }
	  else if ( pos > 0.5 ) 
	    {
	      if ( pi < peak.size() -1 ) { int2 += peak[pi+1] - peak[pi]; cint2++; } 
	    }
	}
      
      
      // assume we will always have at least 2 peaks in each half
      // i.e. this was a detected spindle, but just in case...
      // just in case give a invalid code

      spindle->chirp      = -99999;
      spindle->chirp_fdif = -99999;
      
      bool valid_chirp = cint1 > 1 && cint2 > 1 ; 

      if ( valid_chirp )
	{
	  // go mean from peak-to-peak duration in sample points, to implied frequency, Hz
	  double f1 = 1.0 / ( 2 * ( period_sec * int1/(double)cint1 ) );
	  double f2 = 1.0 / ( 2 * ( period_sec * int2/(double)cint2 ) );

	  // +ve means getting faster: absolute diffference (Hz)
	  spindle->chirp_fdif = f2 - f1 ; 
	  
	  // old CHIRP definition: log scaled ratio
	  spindle->chirp = log( ( int1/(double)cint1 )  / (int2/(double)cint2  )  );
	  
	  std::cout << "   chirp = " << f1 << " " << f2 << " " << spindle->chirp_fdif << " " << spindle->chirp << "\n";
	}


      //
      // Max peak-to-peak, i.e. amplitude
      //
      
      double max_p2p = 0;
      double max_p2p_idx = 0;

      //
      // Lowest trough (i.e. index of location of spindle 'peak')
      //
      

      double lowest     = peak[0];
      int    lowest_idx = 0;

      for (int k=1;k<peak.size();k++)
	{	  
	  
	  const double t = fabs( d[ peak[k] ] - d[ peak[k-1] ] );
	  
	  if ( t > max_p2p ) 
	    {
	      max_p2p = t;
	      max_p2p_idx = ( (peak[k]+peak[k-1] )/(double)2.0) / (double)npoints; // mean, standardarized	      
	    }
	  
	  if ( d[ peak[k] ] < lowest )
	    {
	      lowest     = d[ peak[k] ];
	      lowest_idx = peak[k] ;  
	    }	  

	}
      
      // spindle 'peak' defined as lowest trough
      spindle->peak_sp = lowest_idx; 
      
      // spindle symmetry (based on mid-point of largest peak-to-trough)
      spindle->symm = max_p2p_idx;
      
      // folded symmetry index (i.e. 0 = mid-way; 1 = 0 )
      spindle->symm2 = 2.0 * fabs( spindle->symm - 0.5 ) ;

      // spindle amp
      spindle->amp = max_p2p;



      
      //
      // FFT for modal spindle frequency of spindle
      // (performed on bandpass filtered data)
      //

      real_FFT fft( npoints , MiscMath::nextpow2( npoints ) , Fs , WINDOW_HANN );     
      fft.apply( d );
      int cutoff = fft.cutoff;
      
      double max = 0;
      spindle->fft = 0;
           
      // skip DC component
      for (int j=1;j<fft.cutoff;j++)
	{
	  if ( fft.X[j] > max ) 
       	    {
       	      max = fft.X[j];
       	      spindle->fft = fft.frq[j];
       	    }
	}
     
      spindle->nosc = peak.size() / (double)2.0;
      spindle->frq  = spindle->nosc / (double)spindle->dur;
      
      //
      // FFT on original data, compared to baseline
      //

     
      if ( baseline )
	{  
	  
	  slice_t slice0( edf , s0 , spindle->tp );

	  // copy over for ranges
	  std::map<freq_range_t,double> spindle_fft = *baseline;	  
	  
	  // fixed at:
	  // 0.5..4
	  // 4..8

	  // 10..13.5  <slow spindles> 
	  // 13.5..16  <fast spindles>

	  // 20..30
	  
	  do_fft( slice0.nonconst_pdata() , Fs , &spindle_fft );
	  
	  // calculate enrichment (log10-scale), so set min to v. low...
	  double q_spindle = -999 , q_baseline = -999;

	  std::map<freq_range_t,double>::const_iterator ff = spindle_fft.begin();
	  while ( ff != spindle_fft.end() )
	    {

	      const double & baseline_band_power = (*baseline)[ ff->first ] ;
	      const double & spindle_band_power = ff->second;
	      const freq_range_t & band = ff->first;
	      
	      // relative enrichment (to baseline)  [ log scale ]
	      double re = spindle_band_power - baseline_band_power;
	      
	      // relative enrichment (to baseline)
	      //double re = spindle_band_power / baseline_band_power;
	      
	      // store
	      spindle->enrich[ ff->first ] = re;
	      
	      // calculate overall q score
	      // take 'spindle' as the two middle categories

	      // quality score: 10..16 is spindle range
	      // 
	      if ( band.first <= 16 && band.second >= 10 )
		{
		  // i.e. get largest of slow and fast bands
		  if ( re > q_spindle ) q_spindle = re;
		}
	      else
		{
		  // i.e. get largest of non-spindle bands
		  if ( re > q_baseline) q_baseline = re;
		}
	      ++ff;
	    }
	  
	  // relative relative enrichment [ log scale ]
	  spindle->qual = q_spindle - q_baseline ;

	  // relative relative enrichment [ log scale ]
	  //spindle->qual = q_spindle / q_baseline ;

	  // i.e. max( B_S / B_S_0 ) / max( B_NS / B_NS_0 ) 
	  // or logs
	  
	  // QUAL filter? 
	  
	  if ( qc_q ) 
	    {
	      if ( spindle->qual < qc_qmin ) spindle->include = false;
	      if ( qc_qmax > 0 && spindle->qual > qc_qmax ) spindle->include = false;
	    }

	}

      
      //
      // [ REMOVE ] Optional, time-locked analysis?  [ for QC+ spindles only ] 
      //
      
      if ( locked && spindle->include )
	{
	  
	  // use original signal, plus a window (+/-) 1.5 seconds around center)
	  const double window_sec = 2.0;
	  const int    nbins      = window_sec * Fs ;
	  const double window_left = - ( window_sec / 2.0 );

	  // start point (left of window) for peak minus half window
	  int orig_sp = spindle->start_sp + lowest_idx - ( window_sec / 2.0 ) * Fs ;
	  
	  //	  logger << "os = " <<orig_sp << "\n";
	  
	  uint64_t centre = tp[ lowest_idx ];

	  interval_t i0;
	  i0.set_window( centre , window_sec * globals::tp_1sec );
	  
	  slice_t slice0( edf , s0 , i0 );
	  
	  const std::vector<double> * d0 = slice0.pdata();
	  const std::vector<uint64_t> * tp0 = slice0.ptimepoints(); 
	  
	  const uint64_t step_tp = globals::tp_1sec  * ( window_sec / (double)nbins ) ;
	  const double   step_sec = ( window_sec / (double)nbins );
	  
	  if ( orig_sp >= 0 ) 
	    {
	      for (int l=0;l<d0->size();l++)
		{
		  const uint64_t pos = (*tp0)[l] - (*tp0)[0];
		  const int      bin = pos / step_tp ; 
		  const double   fbin = window_left + bin * step_sec ; 
		  
		  // weight by CWT for spindle...
		  
		  if ( 0 ) 
		    std::cout << "TL\t" 
			      << edf.id << "\t"  
			      << target_f << "\t"
			      << orig_sp << "\t"
			      << i << "\t"
			      << l << "\t" 
			      << (*averaged)[orig_sp] << "\t"
			      << (*d0)[l] << "\t"
			      << (*averaged)[orig_sp] * (*d0)[l] << "\n";
		  
		  ++orig_sp;
		  
		  (*locked)[ fbin ] += (*d0)[l]; 
		}
	      
	    }
	}
      
      
      
      //
      // Next spindle
      //

      if ( ! spindle->include ) removed_some = true;
      
      
     }



   //
   // Prune spindle list?
   //

   if ( removed_some ) 
     {
       std::vector<spindle_t> copy_spindles = *spindles;
       spindles->clear();
       for (int i=0;i<copy_spindles.size();i++)
	 if ( copy_spindles[i].include ) spindles->push_back( copy_spindles[i] );
       logger << "  QC'ed spindle list from " << copy_spindles.size() << " to " << spindles->size() << "\n";
     }
      
   
   //
   // Denominator for mean of spindle-locked average signal
   //
   
   if ( locked ) 
     {
       std::map<double,double>::iterator ll = locked->begin();
       while ( ll != locked->end() )
	 {
	   ll->second /= (double)spindles->size();
	   ++ll;
	 }
     }


   //
   // Remove tmp channel we created
   //
   
   if (  edf.header.has_signal( new_label ) )
     {
       int s = edf.header.signal( new_label );
       edf.drop_signal( s );	  
     }

}


void per_spindle_output( std::vector<spindle_t>    * spindles ,
			 param_t & param , 
			 clocktime_t               * starttime , 
			 std::map<freq_range_t,double> * baseline )
{
 
  const bool enrich_output = param.has( "enrich" );

  const int n = spindles->size();
  
   //
   // Per-spindle output
   //
  
   for (int i = 0 ; i < spindles->size(); i++ ) 
     {
       
       spindle_t * spindle = &(*spindles)[i];
       
       writer.level( i+1 , "SPINDLE" );  // 1-based spindle count
       
       writer.value( "START"  , spindle->tp.start * globals::tp_duration );
       writer.value( "STOP"   , spindle->tp.stop * globals::tp_duration );
       
       writer.value( "START_SP"  , spindle->start_sp );
       writer.value( "STOP_SP"   , spindle->stop_sp  );
       
       if ( starttime != NULL )
	 {

	   double tp1_sec =  spindle->tp.start / (double)globals::tp_1sec;
	   clocktime_t present1 = *starttime;
	   present1.advance_seconds( tp1_sec );
	   // add down to 1/100th of a second
	   double tp1_extra = tp1_sec - (long)tp1_sec;
	   
	   double tp2_sec =  spindle->tp.stop / (double)globals::tp_1sec;
	   clocktime_t present2 = *starttime;
	   present2.advance_seconds( tp2_sec );
	   double tp2_extra = tp2_sec - (long)tp2_sec;
	   
	   writer.value( "START_HMS"  , present1.as_string() +  Helper::dbl2str_fixed( tp1_extra , globals::time_format_dp ).substr(1) );
	   writer.value( "STOP_HMS"   , present2.as_string() +  Helper::dbl2str_fixed( tp2_extra , globals::time_format_dp ).substr(1) );
	   
	 }
       
       writer.value( "AMP"    , spindle->amp      );
       writer.value( "DUR"    , spindle->dur      );
       writer.value( "FWHM"   , spindle->fwhm     );
       writer.value( "NOSC"   , spindle->nosc     );
       writer.value( "FRQ"    , spindle->frq      );
       writer.value( "FFT"    , spindle->fft      );
       writer.value( "SYMM"   , spindle->symm     );
       writer.value( "SYMM2"  , spindle->symm2    );
       writer.value( "ISA"    , spindle->isa      );


	 
       if ( globals::devel )
	 {
	   writer.value( "F_POS"    , spindle->posf      );
	   writer.value( "F_NEG"    , spindle->negf      );
	   writer.value( "F_ALL"    , spindle->allf      );
	   writer.value( "F_DIF"    , spindle->posf - spindle->negf );

	   writer.value( "B_POS"    , spindle->posb      );
	   writer.value( "B_NEG"    , spindle->negb      );
	   writer.value( "B_ALL"    , spindle->allb      );
	   writer.value( "B_DIF"    , spindle->posb - spindle->negb );

	   writer.value( "V_POS"    , spindle->posv      );
	   writer.value( "V_NEG"    , spindle->negv      );
	   writer.value( "V_ALL"    , spindle->allv      );
	   writer.value( "V_DIF"    , spindle->posv - spindle->negv );
	   
	 }
       
       if ( spindle->chirp > -99998 ) 
	 writer.value( "CHIRP"  , spindle->chirp );

       if ( spindle->chirp_fdif > -99998 ) 
	 writer.value( "CHIRPF"  , spindle->chirp_fdif );
	      
       writer.value( "MAXSTAT" , spindle->max_stat );
       writer.value( "MEANSTAT" , spindle->mean_stat );
       
       if ( param.has( "so" ) )
	 {
	   writer.value( "PEAK" , spindle->peak_sec );

	   if ( spindle->so_nearest_num != 0 ) 
	     {
	       writer.value( "SO_NEAREST" , spindle->so_nearest );
	       writer.value( "SO_NEAREST_NUM" , spindle->so_nearest_num );
	     }
	   
	   if ( spindle->so_phase_peak >= 0 ) 
	     writer.value( "SO_PHASE_PEAK" , spindle->so_phase_peak );
	 }

       if ( param.has( "if" ) )
	 writer.value( "IF" , spindle->if_spindle );
       
      
       //
       // Enrichment relative to the baseline
       //
       
       if ( baseline ) 
	 {
	   
	   writer.value( "Q"      , spindle->qual     );
	   writer.value( "PASS"   , spindle->include  );
	   
	   if ( enrich_output )
	     {
	       std::map<freq_range_t,double>::const_iterator bb = spindle->enrich.begin();
	       while ( bb != spindle->enrich.end() )
		 {
		   writer.level( globals::print( bb->first ) , globals::band_strat );
		   writer.value( "ENRICH" , bb->second );
		   ++bb;
		 }
	       writer.unlevel( globals::band_strat );
	     }
	 }

     }

   // end of per-spindle output
   writer.unlevel( "SPINDLE" );
    
}






void do_fft( const std::vector<double> * d , const int Fs , std::map<freq_range_t,double> * freqs )
{

  // Fixed parameters:: use 4-sec segments with 2-second
  // overlaps and Hanning window
  
  double overlap_sec  = 2;
  double segment_sec  = 4;
  double length_sec = d->size() / (double)Fs;

  // check length
  if ( length_sec <= ( segment_sec + overlap_sec ) )
    {
      overlap_sec = 0;
      segment_sec = length_sec;
    }
  
  const int total_points = d->size();
  const int segment_points = segment_sec * Fs;
  const int noverlap_points  = overlap_sec * Fs;
	   
  int noverlap_segments = floor( ( total_points - noverlap_points) 
				 / (double)( segment_points - noverlap_points ) );
  
  //  std::cout << "total_points " << total_points << " " << segment_points << " " << noverlap_points << " " << noverlap_segments << "\n";

  PWELCH pwelch( *d , 
		 Fs , 
		 segment_sec , 
		 noverlap_segments , 
		 WINDOW_HANN );
  
  freqs->clear();

  // 1 .. 25 
  // for (double f=0.5 ; f <= 25.5 ; f++ )
  //   (*freqs)[ freq_range_t( f , f+1 ) ] = 0;
  
  (*freqs)[ freq_range_t( 0.5   , 4    ) ] = 0 ;
  (*freqs)[ freq_range_t( 4     , 8    ) ] = 0 ;
  (*freqs)[ freq_range_t( 10    , 13.5 ) ] = 0 ;
  (*freqs)[ freq_range_t( 13.5  , 16 )   ] = 0 ;
  (*freqs)[ freq_range_t( 20    ,  30 ) ] = 0 ;
  
  // populate
  pwelch.psdmean( freqs );
  
  // log-scale
  std::map<freq_range_t,double>::iterator ff = freqs->begin();
  while ( ff != freqs->end() )
    {
      ff->second = log10( ff->second ); 
      //ff->second =  ff->second ; 
      ++ff;
    }



}



void spindle_stats( const std::vector<spindle_t> & spindles , std::map<std::string,double> & results ) 
{

  double dur = 0 , fwhm = 0 , amp = 0 , nosc = 0 , frq = 0 , fft = 0 , symm = 0 , 
    symm2 = 0, chirp = 0 , chirp_fdif = 0 , isa = 0 , qual = 0 ;

  double negf = 0 , posf = 0 , allf = 0;
  double negb = 0 , posb = 0 , allb = 0;
  double negv = 0 , posv = 0 , allv = 0;

  int denom = 0;

  std::map<freq_range_t,double> enrich; // versus baseline
  
  std::vector<spindle_t>::const_iterator ii = spindles.begin();
  while ( ii != spindles.end() )
    {

      // at this point, all should be included, but just in case
      // keep this in
      if ( ! ii->include ) { ++ii; continue; }

      ++denom; // a QC+ spindle
      dur += ii->dur;
      amp += ii->amp;
      fwhm += ii->fwhm;
      nosc += ii->nosc;
      frq += ii->frq;
      fft += ii->fft;
      symm += ii->symm;
      symm2 += ii->symm2;

      chirp += ii->chirp;
      chirp_fdif += ii->chirp_fdif;

      negf += ii->negf;
      posf += ii->posf;
      allf += ii->allf;
      
      negb += ii->negb;
      posb += ii->posb;
      allb += ii->allb;
      
      negv += ii->negv;
      posv += ii->posv;
      allv += ii->allv;

      isa += ii->isa;
      qual += ii->qual;

      // relative enrichment compared to baseline
      std::map<freq_range_t,double>::const_iterator ss = ii->enrich.begin();
      while ( ss != ii->enrich.end() )
	{
	  enrich[ ss->first ] += ss->second;
	  ++ss;
	}

      ++ii;
    }
  

  results[ "AMP" ]      = amp /(double)denom;
  results[ "TOTDUR" ]   = dur;
  results[ "DUR" ]      = dur / (double)denom;
  results[ "FWHM" ]     = fwhm / (double)denom;
  results[ "NOSC" ]     = nosc / (double)denom;
  results[ "FRQ" ]      = frq / (double)denom;
  results[ "FFT" ]      = fft / (double)denom;
  results[ "SYMM" ]     = symm / (double)denom;
  results[ "SYMM2" ]    = symm2 / (double)denom;

  results[ "CHIRP" ]    = chirp / (double)denom;
  results[ "CHIRPF" ]    = chirp_fdif / (double)denom;

  results[ "FNEG" ]    = negf / (double)denom;
  results[ "FPOS" ]    = posf / (double)denom;
  results[ "FALL" ]    = allf / (double)denom;

  results[ "BNEG" ]    = negb / (double)denom;
  results[ "BPOS" ]    = posb / (double)denom;
  results[ "BALL" ]    = allb / (double)denom;

  results[ "VNEG" ]    = negv / (double)denom;
  results[ "VPOS" ]    = posv / (double)denom;
  results[ "VALL" ]    = allv / (double)denom;

  results[ "Q" ]        = qual / (double)denom;
  
  results[ "ISA_PER_SPINDLE" ] = isa / (double)denom;
  results[ "ISA_TOTAL" ] = isa;
  
  // relative enrichment compared to baseline
  std::map<freq_range_t,double>::iterator ee = enrich.begin();
  while ( ee != enrich.end() )
    {
      results[ "E" + globals::print( ee->first ) ] = ee->second / (double)denom;
      ++ee;
    }
  

}



annot_t * spindle_bandpass( edf_t & edf , param_t & param )
{


  //
  // Attach signals
  // 

  std::string signal_label = param.requires( "sig" );   

  signal_list_t signals = edf.header.signal_list( signal_label );  
  
  const int ns = signals.size();
  
  //
  // Obtain sampling freqs (Hz)
  //
  
  std::vector<double> Fs = edf.header.sampling_freq( signals );

  
  //
  // Annotations to save
  //
  
  annot_t * a = edf.timeline.annotations.add( "spindles-v2" );
  a->description = "Martin et al. spindles" ;
  

  //
  // For each signal
  //
 
  for (int s = 0 ; s < ns ; s++ ) 
    {

      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;

      //
      // Based on: Martin et al. "Topography of age-related changes in
      // sleep spindles", Neurobio Aging 34(2), 2013, pp 468-476
      //
      // Method 'A4' from Warby et al.
      //
      
      // [# Band-pass filter EEG, calculate RMS in sliding windows and apply a
      // constant threshold. Detect a spindle if the RMS exceeds a constant
      // threshold for 0.3-3 s.]
      
      // Parameters

      const double p_resolution    = 0.25;
      const double p_percentile    = 95;
      const double p_window_length = Fs[s] * p_resolution; 
      
      // 1. Bandpass filter signal from C3-M2 in the 11-15 Hz band
      // 2. Calculate the RMS of the bandpass filtered signal with a time
      //    resolution of 25 ms using a time window of 25 ms [# no overlap]
      // 3. threshold <- 95th percentile of RMS signal [# only S2+S3+S4]
      // 4. if ( RMS > threshold  &&  0.3s <= duration above threshold <= 3s )
      //    then [Detect spindle]

      
      //
      // Filter entire signal
      //
      
      // ripple = 0.005 , transition width (Hz) = 0.5 Hz 
      dsptools::apply_fir( edf , signals(s) , fir_t::BAND_PASS ,
			   1 , // Kaiser window
			   0.02 , 0.5 , // ripple , TW
			   10 , 16			   
			   );

      //
      // Get windows of 0.25seconds, no overlap (i.e. advance by 0.25)
      //

      int ne = edf.timeline.set_epoch( p_resolution , p_resolution );

      //
      // Aggregate RMS per window 
      //

      std::vector<double> rms;

            
      //
      // Get data
      //
      
      while ( 1 ) 
	{
	  
	  int epoch = edf.timeline.next_epoch();      
	  
	  if ( epoch == -1 ) break;
	  
	  interval_t interval = edf.timeline.epoch( epoch );
	  
	  //
	  // Get data
	  //

	  slice_t slice( edf , signals(s) , interval );
	  
	  const std::vector<double> * signal = slice.pdata();
	  
	  //
	  // Calculate RMS for window
	  //

	  double t = MiscMath::rms( *signal );
	  
	  rms.push_back(t);
	  
	} // next 0.25 window

      
      //
      // Get threshold (95th percentile)
      //
      
      const int n_bins = rms.size();
      
      const int t95 = n_bins * ( p_percentile / 100.0 );
      
      const double threshold = MiscMath::kth_smallest_preserve( rms , t95 );
      
      const uint64_t bin_ms = p_resolution * globals::tp_1sec;
      
      std::vector<spindle_t> spindles;
      
      uint64_t start = 0 , stop = 0;

      // count of current spindle 'length' (in 0.25s windows
      int scnt = 0;
      
      // for 0.3 to 3s duration, means at least 2 bins,
      // but not more than 12
      
      for (int i=0;i<n_bins;i++)
	{
	  
	  if ( rms[i] >= threshold )
	    {
	      if ( scnt == 0 ) // start of putative spindle
		{
		  scnt = 1;
		  start = i * bin_ms;
		}
	      else // continue a window
		{
		  ++scnt;
		  stop = (i+1) * bin_ms - 1;
		}
	    }
	  else
	    {
	      
	      if ( scnt ) // close a window?
		{
		  if ( scnt >= 2 && scnt <= 12 ) 
		    {
		      spindles.push_back( spindle_t( start , stop , 0 , 0 ) );
		    }
		}
	      scnt = 0; // reset scnt in any case
	    }
	}
    
      //
      // Characterisation of each spindle
      //
      
      bool bandpass_filtered = true;
      
      characterize_spindles( edf , param , signals(s) , bandpass_filtered ,  13 , 4 , // i.e. 11 to 15 
			     "bandpass" ,
			     NULL , NULL , &spindles , NULL , NULL );
      
      std::map<std::string,double> means;

      spindle_stats( spindles , means );

      //
      // Save in the annotation class
      //

      std::vector<spindle_t>::const_iterator ii = spindles.begin();
      while ( ii != spindles.end() )
	{	  
	  const interval_t & spindle = ii->tp;
	  a->add(  signals.label(s)  , spindle , signals.label(s) );      	  	  
	  ++ii;
	}


      //
      // Per-spindle level output 
      //
      
      if ( false ) 
	{
	  
	  std::vector<spindle_t>::const_iterator ii = spindles.begin();
	  int cnt = 0;
	  
	  while ( ii != spindles.end() )
	    {	  
	      
	      const spindle_t & spindle = *ii;	  
	      
	      writer.level( ++cnt , "SPINDLE" );
	      writer.var( "SINGLE_SP_START" , "Single spindle start time-point" );
	      writer.var( "SINGLE_SP_STOP" , "Single spindle stop time-point" );
	      writer.var( "SINGLE_SP_DUR" , "Single spindle stop time-point" );

	      writer.value( "SINGLE_SP_START" , spindle.tp.start * globals::tp_duration );
	      writer.value( "SINGLE_SP_STOP" , spindle.tp.stop * globals::tp_duration );
	      writer.value( "SINGLE_SP_DUR" , (spindle.tp.stop-spindle.tp.start+1)/(double)globals::tp_1sec );

	      ++ii;
	    }
	  writer.unlevel( "SPINDLE" );
	}
      
      const double t_minutes = ( n_bins * p_resolution ) / 60.0 ; 

      bool empty = spindles.size() == 0; 

      if ( empty ) 
	std::cout << "INDIV" << "\t"
		  << edf.id << "\t" 
		  << "[" << globals::current_tag << "]\t"
		  << signals.label(s) << "\t"
		  << 0 << "\t"
		  << t_minutes << "\t"
		  << 0 << "\t"
		  << 0 << "\t" 
		  << "NA\t"
		  << "NA\t"
		  << "NA\t"
		  << "NA\t"
		  << "NA\t"
		  << "NA\n";     
      else
	std::cout << "INDIV" << "\t"
		  << edf.id << "\t" 
		  << "[" << globals::current_tag << "]\t"
		  << signals.label(s) << "\t"
		  << spindles.size() << "\t"
		  << t_minutes << "\t"
		  << spindles.size() / t_minutes << "\t"
		  << means["TOTDUR"] << "\t" 
		  << means["AMP"] << "\t"
		  << means["DUR"] << "\t"
		  << means["NOSC"] << "\t"
		  << means["FRQ"] << "\t"
		  << means["FFT"] << "\t"
// 		  << means["TREND"] << "\t"
// 		  << means["ABSTREND"] << "\t"
		  << means["SYMM"] << "\t"
		  << means["SYMM2"] << "\n"; 
      

    } // Next signal

  
  //  a->save( "spindle.annot" );
  
  return a;
  
}

// helper function
void write_if_exists( const std::string & s , const std::map<std::string,double> & means ) 
{ 
  std::map<std::string,double>::const_iterator ss = means.find( s );
  if ( ss != means.end() ) writer.value( s , ss->second );
}


