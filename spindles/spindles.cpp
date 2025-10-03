
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


// tmp hack:
bool  prune;
double prune_rate;

bool  prune_fixed;
int  prune_n;


#include "spindles.h"
#include "mspindles.h"
#include "plot-spindles.h"
#include "propag.h"
#include "dynamics/qdynam.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "param.h"

#include "annot/annot.h"
#include "intervals/intervals.h"
#include "miscmath/crandom.h"
#include "miscmath/dynam.h"

#include "cwt/cwt.h"
#include "fftw/fftwrap.h"
#include "miscmath/miscmath.h"
#include "artifacts/artifacts.h"
#include "dsp/fir.h"
#include "dsp/hilbert.h"
#include "dsp/spline.h"
#include "spindles/slow-waves.h"
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
  // Signals
  // 
  
  std::string signal_label = param.requires( "sig" );   
  
  // list of signals (data only)

  const bool no_annotations = true;

  signal_list_t signals = edf.header.signal_list( signal_label , no_annotations );  
  
  // number of signals
  const int ns = signals.size();
  
  // nothing to do...
  if ( ns == 0 ) return NULL;
  
  // sampling rate
  std::vector<double> Fs = edf.header.sampling_freq( signals );
  
  

  //
  // Optionally collate all spindles called here, either across all channels:
  //

  mspindles_t mspindles( &edf );

  // ... or within channels only

  std::map<std::string,mspindles_t> ch2mspindles; 

  //
  // Optionally, read in pre-computed spindles from an annotation (i.e. rather than use
  // the CWT);  we still want F set, as this is used in the characterisation stage;  we also
  // want to stratify by channel;
  //

  const bool using_precomputed = param.has( "precomputed" ) ;

  // channel -> precomputed spindle list

  std::map<std::string,std::vector<spindle_t> > precomputed;

  if ( using_precomputed )
    {
      
      // check that all channels have the same SR if using 'precomputed'       
      
      const std::string from_annot = param.value( "precomputed" );
      
      // get sample-points for each channel

      for (int s=0; s<ns; s++)
	{
	  const bool only_extract_this_channel = true;

	  std::vector<interval_t> sample_points , time_points ; 

	  int orig_n;
	  
	  // get sample points, and times
	  int n = edf.timeline.annot2sp(edf, from_annot, only_extract_this_channel,
					&sample_points ,
					&time_points ,
					&orig_n , 
					signals.label(s) );

	  //
	  // build the spindle list
	  //
	  
	  std::vector<spindle_t> sp;

	  // nb. have to remove 1 SP from end 
	  for (int i=0; i<n; i++)
	    sp.push_back( spindle_t( time_points[i].start , time_points[i].stop ,
				     sample_points[i].start , sample_points[i].stop - 1 ) );
	  
	  //
	  // add to the list
	  //
	  
	  precomputed[ signals.label(s) ] = sp; 
	  
	  logger << "  mapped " << sp.size() << " of " 
		 << orig_n << " spindles for "
		 << signals.label(s) << "\n";
	  
	}

    }
  
  
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

  // winsorize cwt-amplitude
  const bool winsorize              = param.has( "winsor" );
  const double winsorize_th         = winsorize ? param.requires_dbl( "winsor" ) : 0;

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
  // Epoch-level dynamics?
  //
  
  const bool calc_dynamics = param.has( "dynam" );
  
  //
  // Analysis/output parameters
  //

  // epoch-level output
  const bool     show_epoch_level         = calc_dynamics || param.has( "epoch" );

  // spindle-level output
  const bool     show_spindle_level       = param.has( "per-spindle" );
  
  
  // verbose display of all CWT coefficients
  const bool     show_cwt_coeff           = param.has( "show-coef" );


  //
  // Misc
  //

  prune = param.has( "prune" );
  prune_rate  = prune ? param.requires_dbl( "prune" ) : 0 ;

  prune_fixed = param.has( "fixed" ) ;
  prune_n     = prune_fixed ? param.requires_int( "fixed" ) : 0 ;

  if ( prune && prune_fixed ) Helper::halt( "cannot specify both prune and fixed" );
  
  
  //
  // Coupling
  //

  // detect slow waves and estimate ITPC etc for spindle start/anchor/stop and slow waves
  const bool     sw_coupling              = param.has( "so" );
  
  // by default, SW coupling is for same channel as the original; but we can allow this to vary
  //  (if null "", below use same channel)
  const std::string sw_channel = ( ! param.empty( "so" ) ) ? param.value( "so" ) : "" ; 

  const bool  external_sw_channel = sw_channel != "";
  
  // by default, SW phase coupling anchor is with "peak" (point of max CWT amplitude)
  //  always calculate this, and use this in the primary spindle level output 

  // however, we can also evaluate SO-coupling w/ different anchors -- this
  // is sent to a new strata: F CH ANCHOR 
  // alternatively, consider position metrics (-1 start , 0 middle , +1 stop)
  // 0 spindle peak, -1 start, +1 end :: choice of coupling metric
  
  // alternatively, add second offset 
  // here we take a vector of anchors and offsets
  
  // anchor based on -1 .. +1 for spindle start/stop
  std::vector<double> sw_coupling_anchor = param.dblvector( "anchor" );
  const bool sw_coupling_has_anchors = sw_coupling_anchor.size();
  
  // scale -1 to +1 becomes 0 to 1
  for (int i=0; i<sw_coupling_anchor.size(); i++)
    {
      sw_coupling_anchor[i] = ( 1.0 + sw_coupling_anchor[i] ) / 2.0 ;
      
      if ( sw_coupling_anchor[i] < 0 || sw_coupling_anchor[i] > 1 ) 
	Helper::halt( "expecting anchor values between -1 and 1" );
    }
  
  // optional offset - if specified, must must 'anchor' length
  // if only offset if specified, then this is applied to max peak
  std::vector<double> sw_coupling_offset  = param.dblvector( "offset" );
  const bool sw_coupling_has_offsets = sw_coupling_offset.size();
  
  if ( sw_coupling_has_anchors && sw_coupling_has_offsets ) 
    {
      if ( sw_coupling_anchor.size() != sw_coupling_offset.size() )
	Helper::halt( "anchor and offset must be of similar length, if both specified" ); 
    }
  
  // if not specified, set to zero to match anchor length
  if ( sw_coupling_has_anchors && ! sw_coupling_has_offsets ) 
    sw_coupling_offset.resize( sw_coupling_anchor.size() , 0 );

  // or, add a single offset, which will match w/ the peak-ampltiude (default) anchor
  if ( ! ( sw_coupling_has_anchors || sw_coupling_has_offsets ) )
    sw_coupling_offset.push_back( 0 );
  
  
  //
  // Non-SO coupling: we may have a channel we want to consider phase
  // and look at spindle coupling, but not detect canonical SO -
  // e.g. a non-EEG channel.  here, we specify by the phase=<channel>
  // option instead
  //

  const bool phase_coupling = param.has( "phase" ) ;
  if ( phase_coupling && sw_coupling )
    Helper::halt( "cannot specify both 'so' and 'phase' options together" );
  if ( phase_coupling && param.empty( "phase" ) ) Helper::halt( "expecting phase=channel" );
  const std::string phase_channel = phase_coupling ? param.value( "phase" ) : "";
  double phase_lwr = 1, phase_upr = 4;
  if ( phase_coupling ) 
    {
      std::vector<double> pp = param.dblvector( "phase-band" );
      if ( pp.size() != 2 ) Helper::halt( "expecting phase-band=lwr,upr" );
      phase_lwr = pp[0];
      phase_upr = pp[1];
      if ( phase_lwr >= phase_upr ) Helper::halt( "expecting phase-band=lwr,upr" );      
    }
  const double fir_tw = param.has( "tw" ) ? param.requires_dbl( "tw" ) : 1 ; 
  const double fir_ripple = param.has( "ripple" ) ? param.requires_dbl( "ripple" ) : 0.02 ; 


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
 
  

  //
  // Caches:
  //   - spindle peaks (cache-peaks); spindle mid-points (cache-mids)
  //   - wavelet power (cache)
  //   - metrics (DENS, etc)
  //
  
  const bool cache_data                   = param.has( "cache" );
  const std::string cache_name = cache_data ? param.value( "cache" ) : "" ;
  
  const bool cache_peaks                  = param.has( "cache-peaks" );
  const std::string cache_peaks_name = cache_peaks ? param.value( "cache-peaks" ) : "";
  const bool cache_amp_peak  = param.has( "amp-peaks" ) ;

  const bool cache_mids                  = param.has( "cache-mids" );
  const std::string cache_mids_name = cache_mids ? param.value( "cache-mids" ) : "";
  
  const bool cache_peaks_sec             = param.has( "cache-peaks-sec" );
  const std::string cache_peaks_sec_name = cache_peaks_sec ? param.value( "cache-peaks-sec" ) : "";


  //
  // Spindle-level QC filters 
  //

  spindle_qc_t q( param );

  //
  // Spindle propagation
  //

  const bool verbose_prop = param.has( "prop-verbose" ) || param.has( "verbose-prop" ) ;

  const bool do_prop = verbose_prop || param.has( "prop" );
  
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
  // Set up annotation
  //
  
  std::string sp_label = "SP";
  
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
  // Get SO/Phase channel if explicitly specified?  Only need to do once
  //
  
  std::vector<double> dxslow;
  std::vector<uint64_t> tpxslow;
  int srslow = 0;

  if ( sw_channel != "" || phase_channel != "" )
    {
      std::string ext_ch =  phase_channel == "" ? sw_channel : phase_channel ;
      const int slow_slot = edf.header.signal( ext_ch ) ;
      if ( slow_slot == -1 ) Helper::halt( "could not find " + ext_ch );
      
      if ( edf.header.is_annotation_channel( slow_slot ) )
	Helper::halt( ext_ch + " is an annotation channel" );

      srslow = edf.header.sampling_freq( slow_slot );
      
      // get data
      slice_t slice( edf , slow_slot , edf.timeline.wholetrace() );
      dxslow = *slice.nonconst_pdata();
      tpxslow = *slice.nonconst_ptimepoints();      
      logger << "  attached SO/phase channel " << ext_ch << " SR=" << srslow << "Hz\n";
    }
  
  
   
  //
  // For each signal, over the whole signal
  //
  
  interval_t interval = edf.timeline.wholetrace(); 
  
  for (int s = 0 ; s < ns ; s++ ) 
    {
      
      //
      // Only consider raw signal channels, with at least 50 Hz SR
      //
      
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;

      if ( Fs[s] < 50 )
	{
	  logger << "  skipping channel " <<  signals.label(s) << " with sample rate of " << Fs[s] << " < 50 Hz\n";
	  continue;
	}
      
      //
      // If specified, SO/phase channel should have the same sample rate (?necessary?)
      //

      if ( srslow != 0 && srslow != Fs[s] )
	{
	  Helper::halt( "sample rate for sig and so/phase channels must be the same; use RESAMPLE" );
	}
      
      //
      // Output
      //
      
      writer.level( signals.label(s) , globals::signal_strat );

      //
      // track epoch-level dynamics?
      //
      
      qdynam_t qd;
      if ( calc_dynamics )
	qd.init( edf , param );

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
      // Set 
      //

      std::vector<double> * dslow = srslow == 0 ? (std::vector<double> *)d : &dxslow ; 

      std::vector<uint64_t> * tpslow = srslow == 0 ? (std::vector<uint64_t> *)tp : &tpxslow ;

      
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
      // Baseline CWT on signal
      //

      q.init( d , Fs[s] );

      //
      // This is only populated if we are considering multiple frequencies
      //
      
      std::map<double,std::vector<spindle_t> > f2int;      


      //
      // Set up for optional (slow-wave) coupling
      //
      
      hilbert_t * p_hilbert = NULL ;
      
      slow_waves_t * p_sw = NULL ;
     
      
      if ( sw_coupling || phase_coupling )
	{
	  
	  // get SW param from the options
	  slow_wave_param_t sw_par( param );
	  
	  // do not explicitly detect SOs if in phase mode
	  if ( phase_coupling ) sw_par.skip = true;
	  
	  // set channel label
	  sw_par.ch = external_sw_channel ? sw_channel : signals.label(s);
	  
	  // filter-Hilbert raw signal for SWs
	  p_hilbert = sw_coupling ? 
	    new hilbert_t( *dslow , Fs[s] , sw_par.f_lwr , sw_par.f_upr , sw_par.fir_ripple , sw_par.fir_tw ) :
	    new hilbert_t( *dslow , Fs[s] , phase_lwr , phase_upr , fir_ripple , fir_tw );

	  // find slow-waves (note: inefficient as filters slow channel twice; but if p_hilbert may
	  // be pointing at a different channel, can't always use the same variant... but should fix
	  p_sw = new slow_waves_t( *dslow , *tpslow , Fs[s] , sw_par ) ;
	  
	  // and phase
	  p_sw->phase_slow_waves();

	  // and display 
	  if ( sw_coupling ) 
	    p_sw->display_slow_waves( param.has( "verbose" ) , &edf );
	  
	  //
	  // display SO characteristics 
	  //

	  if ( sw_coupling && verbose_time_phase_locking ) 
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
	    logger << "  wavelet with FWHM(T) " << fwhm[fi];       
	  else
	    logger << "  wavelet with " << num_cycles << " cycles";       

	  logger << ", smoothing window " << moving_window_sec << "s\n";

	  //
	  // Output stratifier: F_C
	  //

	  writer.level( frq[fi] , globals::freq_strat );
	  
	  //
	  // Get results for this F_C
	  //

	  std::vector<double> results = cwt.results(fi);


	  //
	  // Add to Q-filter (before any winsorizing or smoothing, etc)
	  //

	  q.init_spindle( results );
	  
	  
	  //
	  // Optionally winsorize
	  //

	  double emp_winsorize = 0;
	  std::vector<bool> winsorized_sp;
	  if ( winsorize )
	    {
	      logger << "  winsorizing at " << winsorize_th << "\n";
	      winsorized_sp.resize( results.size() , false );
	      MiscMath::winsorize( &results , winsorize_th ,
				   1 ,  // 1 = one-sided, not that it should matter
				   NULL , // ignore lower
				   &emp_winsorize , // save upper winsor-th
				   &winsorized_sp // store which samples winsorized
				   );  
	    }
	  
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
	      
	      if ( t_median > 0 ) 
		writer.value( "MEAN_OVER_MEDIAN" , t_mean / t_median );
	      
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

	  logger << "  detection thresholds (core, flank, max)  = "
		 << multiplicative_threshold << ", "
		 << boundary_threshold;
	  if ( maximal_threshold > 0 ) 	 
	    logger << ", " << maximal_threshold;
	  logger << "x" << "\n";
	  
	  logger << "  core duration threshold (core, min, max) = " << min0_dur_sec << ", "
		 << min_dur_sec << ", " << max_dur_sec << "s" << "\n";	  

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
	  // Find above threshold regions (detect given CWT,or use a precomputed list)
	  //
	  
	  std::vector<interval_t> spindles1;
	  std::vector<int> spindles1_start; // sample-points
	  std::vector<int> spindles1_stop; // sample-points

	  if ( using_precomputed )
	    {
	      std::vector<spindle_t> prespins = precomputed[ signals.label(s) ];
	      const int n = prespins.size();

	      for (int i=0;i<n; i++)
		{
		  spindles1.push_back( prespins[i].tp );
		  spindles1_start.push_back( prespins[i].start_sp );
		  spindles1_stop.push_back( prespins[i].stop_sp );
		}
	      logger << "  attached " << n <<  " precomputed spindles...\n";
	      
	    }	  
	  else if ( use_zpks )
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
		  
		  // tp: use +1 end encoding (i.e. as s.stop is +1 end in SP space already)
		  spindles1.push_back( interval_t( (*tp)[s.start] , (*tp)[s.stop] ) );
		  
		  // sp: smoothedZ returns 1 past end encoding; adjust here for sp
		  spindles1_start.push_back( s.start );
		  spindles1_stop.push_back( s.stop - 1 );
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
		      stop = i; // end encoding for SP
		      ++scnt;
		    }
		  else
		    {
		      
		      if ( scnt > 0 ) 
			{
			  
			  uint64_t start_tp = (*tp)[ start ];
			  uint64_t stop_tp  = (*tp)[ stop + 1 ] ; 
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
				  //stop = j+1; // end-encoding	  
				  stop = j; // end-encoding
				}
			      
			      // re-adjusted start/stop times
			      start_tp = (*tp)[ start ];
			      stop_tp  = (*tp)[ stop + 1 ]; // 1-TP past end encoding
			      
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
		  
		  // merge? if both ends are within X seconds (default = 0.5 ) 
		  // but if using a precomputed spindle list, take exactly as is
		  // i.e. could have completely overlapping events, but that is fine

		  // overlap?
		  if ( (!using_precomputed) && this_start < previous_stop ) { extending = true; }
		  
		  // too near?
		  else if ( (!using_precomputed) && this_start - previous_stop < spindle_merge_tp ) { extending = true; }
		  
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
		  
		}

	    } 


	  int nspindles_postmerge = spindles.size();
	  
	  logger << "  merged nearby intervals: from " 
		 << spindles1.size() << " to " 
		 << spindles.size() << " unique events\n";
	  


	  

	  // ------------------------------------------------------------
	  // 
	  // Spindle detection over -- now do spindle analysis 
	  //  - only requires 'spindles[]  and averaged[] 
	  //
	  // ------------------------------------------------------------
	  
	  //
	  // Track whether each SP is in a spindle or not
	  // Currently, only used when 'if' option
	  //
	  
	  std::vector<bool> in_spindle;
	  std::vector<bool> in_pos_hw; // (neg2pos ZC - pos2neg ZC)
	  std::vector<bool> in_neg_hw; // (pos2neg ZC - neg2pos ZC)
	  std::vector<bool> in_pos_slope; // ( beak peak - pos peak ) 
	  std::vector<bool> in_neg_slope; // ( beak peak - pos peak ) 

	  if ( ht_chirp ) 
	    {
	      in_spindle.resize( averaged.size() , false );
	      for (int i=0;i<spindles.size();i++) 
		for (int j=spindles[i].start_sp; j <= spindles[i].stop_sp; j++) 
		  in_spindle[j] = true;

	      // these get editted later, i.e. to zero out additional terms 
	      // based on exact popsition in the spindle 
	      in_pos_hw = in_spindle;
	      in_neg_hw = in_spindle;
	      in_pos_slope = in_spindle;
	      in_neg_slope = in_spindle;	     
	    }
	  
	  std::vector<bool> * p_in_pos_hw = ht_chirp ? &in_pos_hw : NULL ; 
	  std::vector<bool> * p_in_neg_hw = ht_chirp ? &in_neg_hw : NULL ; 
	  std::vector<bool> * p_in_pos_slope = ht_chirp ? &in_pos_slope : NULL ; 
	  std::vector<bool> * p_in_neg_slope = ht_chirp ? &in_neg_slope : NULL ; 
	 

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
	  for ( int i=0;i<averaged_corr.size();i++)
	    averaged_corr[i] /= threshold[i];
	  

	  
	  
	  
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
	      
	      characterize_spindles( edf , param , signals(s) , 
				     bandpass_filtered_status , 				     
				     frq[fi] ,
				     "wavelet-" + Helper::dbl2str(frq[fi]) , 
				     &averaged_corr ,    // pass as input threshold-normed CWT
				     d ,                 // original EEG signal 
				     &spindles ,         // this will be annotated/reduced   
				     ( hms ? &starttime : NULL) , 
				     &q, 
				     tlocking ? &locked : NULL ,     // mean signal around spindle troughs
				     p_in_pos_hw , p_in_neg_hw , p_in_pos_slope , p_in_neg_slope , 
				     winsorize ? &winsorized_sp : NULL 
				     );	      
	      
	    }
	  

	  
	  //
	  // Get mean spindle parameters for this channel/frequency 
	  //

	  if ( characterize ) 
	    spindle_stats( spindles , means );

		  
	  //
	  // Verbose signal display with thresholds
	  //

	  if ( show_cwt_coeff )
	    {
	      
	      int np = cwt.points();
	      int nf = cwt.freqs();      
	      if ( np != np0 ) Helper::halt( "internal problem in cwt()" );
	      
	      // add spindle track for initial/putative spindles (pre-merge, QC)
	      std::vector<bool> sputative(np, false);
	      if ( spindles1.size() != 0 )
		{
		  
		  interval_t curr = spindles1[0];
		  int si = 0;
		  int ti = 0;
		  int ns = spindles1.size();
		  while ( ti < np )
		    {
		      // skip to next event? 
		      if ( (*tp)[ti] >= curr.stop )
			{
			  ++si;
			  if ( si == ns ) break;
			  curr = spindles1[si];
			}

		      if ( (*tp)[ti] >= curr.start && (*tp)[ti] < curr.stop )
			sputative[ti] = true;

		      ++ti;		      
		    }
		}
	      
	      std::vector<bool> sfinal(np, false);
	      if ( spindles.size() != 0 )
		{
		  
		  interval_t curr = spindles[0].tp;
		  int si = 0;
		  int ti = 0;
		  int ns = spindles.size();
		  while ( ti < np )
		    {
		      // skip to next event? 
		      if ( (*tp)[ti] >= curr.stop )
			{
			  ++si;
			  if ( si == ns ) break;
			  curr = spindles[si].tp;
			}

		      if ( (*tp)[ti] >= curr.start && (*tp)[ti] < curr.stop )
			sfinal[ti] = true;

		      ++ti;		      
		    }
		}
	      
	      
	      for (int ti=0;ti<np;ti++)
		{		  
		  writer.interval( interval_t( (*tp)[ti] , (*tp)[ti] ) );
		  writer.value( "SEC" , (*tp)[ti] * globals::tp_duration );
		  writer.value( "RAWCWT" , cwt.raw_result(fi,ti)  );
		  writer.value( "CWT" , cwt.result(fi,ti)  );
		  writer.value( "AVG" , averaged[ti] );
		  writer.value( "AVG_CORR" , averaged_corr[ti] );
		  writer.value( "CWT_TH" , threshold[ti] );
		  writer.value( "CWT_TH2" , threshold2[ti] );
		  writer.value( "CWT_THMAX" , threshold_max[ti] );
		  writer.value( "PUTATIVE" , (int)sputative[ti] );
		  writer.value( "SPINDLE" , (int)sfinal[ti] );
		}
	      writer.uninterval();
	    }
	  


	  
	  //
	  // If we're doing CHIRP analyses, get the mean observed spindle frequency (which we will
	  // use of the basis for the filter-Hilbert characterisation of IF (note, this itself is still
	  // based on some BPF around the target frequency, but not too narrow
	  //
	  
	  double observed_frq = means["FRQ"];
	  
	  logger << "  estimated spindle frequency is " << observed_frq << "\n";


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
	      // 'or' if if-frq-emp=X then use observed spindle freq +/- X  <<- best option

	      // set a broad transition width 
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
		  isf[j] /= (double)isfn[j] ;
		  writer.value( "IF" , isf[j] );
		}
	      writer.unlevel( "RELLOC" );
	      
	      // variance of RELLOC freq
	      double vbin = MiscMath::sdev( isf );
	      writer.value( "FVAR" , vbin );

	      double fmin = 0;
	      double fmax = 0;	      
	      MiscMath::minmax( isf , &fmin , &fmax );
	      writer.value( "FRNG" , fmax - fmin );
	      
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
	  // Cache spindle peaks: in sample points 
	  //  'peaks' or literal mid-points
	  //
	  
	  if ( cache_peaks )
            {
              cache_t<int> * cache = edf.timeline.cache.find_int( cache_peaks_name );

	      std::vector<int> peaks;
	      for (int i=0; i<spindles.size(); i++)
		{
		  int p = spindles[i].start_sp + spindles[i].max_trough_sp;
		  peaks.push_back(p);
		}
              cache->add( ckey_t( "points" , writer.faclvl() ) , peaks );
	    }
	  
	  if ( cache_mids )
	    {
              cache_t<int> * cache = edf.timeline.cache.find_int( cache_mids_name );

	      std::vector<int> peaks;
	      for (int i=0; i<spindles.size(); i++)
		{
		  int p = spindles[i].start_sp + ( spindles[i].stop_sp - spindles[i].start_sp ) / 2;
		  peaks.push_back(p);
		}
              cache->add( ckey_t( "points" , writer.faclvl() ) , peaks );
	    }


	  //
	  // Cache spindle peaks: in seconds
	  //
	  
	  if ( cache_peaks_sec )
            {
	      logger << "  caching peaks (seconds) in " << cache_peaks_sec_name << "\n";
              cache_t<double> * cache = edf.timeline.cache.find_num( cache_peaks_sec_name );
	      
	      std::vector<double> peaks;
	      for (int i=0; i<spindles.size(); i++)
		{
		  spindle_t & spindle = spindles[i];
		  if ( cache_amp_peak )
		    {
		      double peak_amp_sec = spindle.peak_amp_sp / (double)( spindle.stop_sp - spindle.start_sp ) ;
		      peak_amp_sec = spindle.tp.start + peak_amp_sec * ( spindle.tp.stop - spindle.tp.start ) ;
		      peaks.push_back( peak_amp_sec * globals::tp_duration );
		    }
		  else
		    {
		      peaks.push_back( spindle.tp_mid * globals::tp_duration );
		    }
		      
		}
              cache->add( ckey_t( "seconds" , writer.faclvl() ) , peaks );
	    }

	  //
	  // Align SO events/phase w/ these spindles (sw_coupling)
	  // Alternatively, align the phase bins of a specified channel w/ spindles (phase_coupling)
	  //  i.e. in the latter case, we do not try to detect SO discretely
	  //
	  
	  if ( sw_coupling || phase_coupling )
	    {
	      
	      //
	      // slow waves have already been detected for this channel (in 'so')
	      //

	      
	      //
	      // SO-stratifed analysis (i.e. not COUPL analyses per se)
	      //
	      
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
	      // Primary spindle/SO coupling & overlap analysis
	      //


	      // always iterate over offset vector
	      // either:
	      //   anchor = PEAK AMPLITUDE, run for each 
	      
	      for ( int a=0; a<sw_coupling_offset.size(); a++)
		{
		  const bool use_peak_amplitde = ! sw_coupling_has_anchors;
		  
		  std::string anchor_label = use_peak_amplitde ? 
		    "A" : ( "T" + Helper::dbl2str( sw_coupling_anchor[a] ) );		  
		  if ( sw_coupling_offset[a] < -0.0000001 ) 
		    anchor_label += "minus" + Helper::dbl2str( sw_coupling_offset[a] ).substr(1) + "s";
		  else if ( sw_coupling_offset[a] > 0.0000001 ) 
		    anchor_label += "plus" + Helper::dbl2str( sw_coupling_offset[a] ) + "s";

		  // only use more complex output form is > 1 anchor
		  
		  if ( sw_coupling_has_anchors ) 
		    {
		      logger << "  running coupling analysis for anchor tag " << anchor_label << "\n";
		      
		      writer.level( anchor_label , globals::anchor_strat );		      
		    }

		  
		  // estimate SO phase at spindle 'anchor' point
		  // and a boolean for whether it falls within a SO
		  
		  // typically, anchor is "peak" or point of max wavelet
		  // alternatively, if 'couple' arg is set, it can be 
		  // positional (start/middle/stop)
		  
		  // not all spindles will have a valid anchor (i.e. if it can be outside
		  // the spindle body - due to discontinuities), and so we have to track
		  // which subset of spindles are included here [ in anchor_idx ] 
	      
		  std::map<int,int> anchor_idx;
	      
		  // SO phase for this spindle 
		  std::vector<double> ph_anchor;
		  
		  // was this spindle anchor in a SO?
		  std::vector<bool> sw_anchor;
		  		  
		  // all spindles
		  const int nspindles = spindles.size();
	      
		  std::vector<int> sw_spindles_anchor; 
		  std::vector<int> valid_spindles_anchor;
		  
		  std::vector<double> nearest_sw;
		  std::vector<int> nearest_sw_number;
		  std::vector<uint64_t> spindle_anchor;
		  
		  for (int i=0;i<nspindles;i++)
		    {
		      // spindle start/stop in sample points
		      int b0 = spindles[i].start_sp;
		      int b1 = spindles[i].stop_sp;
		      
		      // define spindle anchor 'ba' in one of a few ways, below		  
		      int ba = b0;
		      
		      // positional: anchor defined as range 0..1 from start/stop
		      if ( ! use_peak_amplitde )
			{		      
			  // n.b. b1 is *inclusive* 
			  int blen = b1 - b0 + 1 ; 
			  ba = b0 + sw_coupling_anchor[a] * blen;
			}
		      else // maximal CWT amplitude
			{
			  double mx = averaged[b0];
			  int mxi = b0;
			  for (int j = b0+1 ; j <= b1 ; j++ )
			    if ( averaged[j] > mx ) { mx = averaged[j] ; ba = j; } 		      
			}		  
		      
		      //
		      // any offset to the anchor?
		      //
		      
		      const int offset_sp = sw_coupling_offset[a] * Fs[s]; 
		      
		      //
		      // test for offset inducing a discontinuity (this includes if
		      // we go off the start/end of the map 
		      //  we know that 'ba' will start being somewhere within the spindle
		      // body, and that every spindle will be continuous
		      // so just test ba vs ba+offset
		      
		      bool valid_anchor = true;
		      if ( offset_sp < 0 ) 
			valid_anchor = ! timeline_t::discontinuity( *tp , Fs[s] , ba + offset_sp , ba );
		      else if ( offset_sp > 0 ) 
			valid_anchor = ! timeline_t::discontinuity( *tp , Fs[s] , ba , ba + offset_sp );
		      
		      //
		      // If this is a valid anchor, include in the SO-phase/overlap analysis
		      //
		      
		      if ( valid_anchor ) 
			{
			  
			  const int valid_idx = anchor_idx.size();
			  
			  // track that this spindle is in the SO phase/overlap analysis
			  anchor_idx[ i ] = valid_idx;
			  
			  // calculate anchor (sp)
			  ba += offset_sp ;
			  
			  // spindle anchor in SO ? 
			  if ( p_sw->in_slow_wave(ba) ) 
			    sw_spindles_anchor.push_back( ba ); 
			  
			  // record all (valid) anchors 
			  valid_spindles_anchor.push_back( ba ); 
			  
			  // time (seconds) to nearest SO from this anchor
			  int sw_num = 0;
			  nearest_sw.push_back( p_sw->nearest(ba , &sw_num ) ); 
			  nearest_sw_number.push_back( sw_num + 1 ); // make 1-based for visual output
			  spindle_anchor.push_back( (*tp)[ ba ] );
			}
		      
		      //
		      // next spindle
		      //
		    }
	      
		  //
		  // SO-phase for each spindle in a SO
		  //
		  
		  if ( valid_spindles_anchor.size() > 0 ) 
		    {
		      
		      //
		      // restrict SO/spindle coupling calculations to
		      // spindles that occur within a SO? Or use all
		      // spindles?
		      //
		      
		      std::vector<bool> so_mask;		  
		      
		      // default is to use mask in SO mode; in phase mode, do not
		      bool use_mask = phase_coupling ? false : ! param.has( "all-spindles" );
		      
		      // if non-SO coupling (i.e. other band phases specified, not detected SO) , this is not needed
		      if ( phase_coupling ) use_mask = false;
		      
		      if ( use_mask ) 
			so_mask = p_sw->sp_in_sw_vec();
		      

		      //
		      // report coupling overlap by SO/external channel phase
		      //
		      
		      bool stratify_by_so_phase_bin = param.has ( "stratify-by-phase" ) || phase_coupling ; 
		      
		      
		      //
		      // Within-epoch permutation (default)
		      //
		      
		      bool eperm = ! param.has( "perm-whole-trace" );
		      
		      double epoch_sec = 0 ;
		      
		      int sr = Fs[s];
		      
		      if ( eperm ) 
			{
			  if ( edf.timeline.generic_epochs() )
			    {
			      if ( param.has( "fixed-epoch-dur" ) )
				epoch_sec = param.requires_dbl( "fixed-epoch-dur" );
			      else
				Helper::halt( "cannot run within-epoch permutation with generic (differently-sized) epochs: add 'perm-whole-trace' " );
			    }
			  else 
			    {
			      if ( ! edf.timeline.epoched() ) edf.timeline.ensure_epoched();
			      epoch_sec = edf.timeline.epoch_length();
			    }
			}
		      
		      //
		      // Permutation for ITPC values?
		      //
		      
		      int nreps = 0;
		      
		      if ( param.has( "nreps" ) ) nreps = param.requires_int( "nreps" );
		      
		      if ( nreps != 0 && nreps < 10 ) Helper::halt( "nreps must be 10+" );
		      
		      
		      //
		      // Perform spindle/SO coupling analysis
		      //
		      
		      itpc_t itpc = p_hilbert->phase_events( valid_spindles_anchor , 
							     use_mask ? &so_mask : NULL , 
							     nreps , 
							     sr , 
							     epoch_sec , // opt: within-epoch shuffle
							     stratify_by_so_phase_bin // opt: overlap by SO-phase bin
							     );
		      
		      //
		      // Repeat for all spindles? - i.e. always run twice ratehr than use 'all-spindles' option
		      //
		      
		      itpc_t itpc_all;
		      if ( use_mask ) 
			itpc_all = p_hilbert->phase_events( valid_spindles_anchor ,
							    NULL ,  // i.e. same as above, but no mask
							    nreps ,
							    sr ,
							    epoch_sec , // opt: within-epoch shuffle                               
							    false      //  opt: overlap by SO-phase bin               
							    );
		      
		      
		      // n.b. these only include all valid spindles
		      // i.e. not necessarily all spindles
		      //  anchor_idx[] maps this set to the originals
		      
		      sw_anchor = itpc.event_included;
		      ph_anchor = itpc.phase;
		      
		      //
		      // if >1 We can output now, as this is stratified by ANCHOR + F + CH
		      //  i.e. so does not need to be interleaved with the F + CH outputs
		      //  in the case of -t mode
		      //
		      
		      // ITPC magnitude of coupling
		      
		      // i.e. by default overlap with mid-point
		      writer.value( "COUPL_ANCHOR" , (int)anchor_idx.size() );
		      
		      writer.value( "COUPL_MAG" , itpc.itpc.obs );
		      
		      

		      //
		      // gross overlap 
		      //
		      
		      if ( use_mask )
			{
			  writer.value( "COUPL_OVERLAP" , itpc.ninc.obs );
			  writer.value( "CDENS" , itpc.ninc.obs / t_minutes ) ;
			  writer.value( "UDENS" , ( spindles.size() - itpc.ninc.obs ) / t_minutes ) ;
			  
			  // phase coupling stats based on all spindles (not just those overlapping a SO)
			  writer.value( "COUPL_ALL_MAG" , itpc_all.itpc.obs );
			  
			}
		      
		      
		      //
		      // phase coupling
		      //

		      if ( nreps ) 
			{
			  
			  writer.value( "COUPL_MAG_EMP" , itpc.itpc.p );
			  writer.value( "COUPL_MAG_NULL" , itpc.itpc.mean );
			  writer.value( "COUPL_MAG_Z" , ( itpc.itpc.obs - itpc.itpc.mean ) / itpc.itpc.sd );
			  
			  // proportion of spindles that overlap a SO 
			  // unless itpc-so was set, this will be meaningless
			  // so only report is a 'mask' was set		  
			  if ( use_mask ) 
			    {
			      writer.value( "COUPL_OVERLAP_EMP" , itpc.ninc.p );
			      writer.value( "COUPL_OVERLAP_NULL" , itpc.ninc.mean );
			      if ( itpc.ninc.sd > 0 )  		    			
				writer.value( "COUPL_OVERLAP_Z" , ( itpc.ninc.obs - itpc.ninc.mean ) / itpc.ninc.sd );

			      // phase coupling based on all spindles (not just overlapping a SO)
			      writer.value( "COUPL_ALL_MAG_EMP" , itpc_all.itpc.p );
			      writer.value( "COUPL_ALL_MAG_NULL" , itpc_all.itpc.mean );
			      writer.value( "COUPL_ALL_MAG_Z" , ( itpc_all.itpc.obs - itpc_all.itpc.mean ) / itpc_all.itpc.sd );
			      
			      
			    }
			}
		      
		      //
		      // mean angle; no empirical test results; -9 means no events observed, so set to missing
		      //
		      
		      if ( itpc.angle.obs > -9 ) 
			writer.value( "COUPL_ANGLE" , itpc.angle.obs );
		      
		      if ( use_mask && itpc_all.angle.obs > -9 )
			writer.value( "COUPL_ALL_ANGLE" , itpc_all.angle.obs );
		      

		      //
		      // asymptotic significance of coupling test; under
		      // the null, give mean rate of 'significant'
		      // (P<0.05) coupling
		      //
		      
		      writer.value( "COUPL_PV" , itpc.pv.obs ) ;
		      if ( use_mask ) 
			writer.value( "COUPL_ALL_PV" , itpc_all.pv.obs ) ;
		      
		      if ( nreps ) 
			{
			  writer.value( "COUPL_SIGPV_NULL" , itpc.sig.mean ) ; 
			  if ( use_mask ) 
			    writer.value( "COUPL_ALL_SIGPV_NULL" , itpc_all.sig.mean ) ; 
			}
		      

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
			       writer.value( "COUPL_OVERLAP_EXP"   , itpc.phasebin[b].mean );
			       writer.value( "COUPL_OVERLAP_EMP"    , itpc.phasebin[b].p );

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
		   // Individual spindle anchors (to be output later) in characterize() 
		   // n.b. have to map from valid set to originals [ anchor_idx ]
		   //
		   //  ** only set these for the first run considered
		  //
		  
		  if ( a == 0 )
		    {
		      for (int i=0;i<nspindles;i++)
			{
			  
			  const bool valid = anchor_idx.find( i ) != anchor_idx.end();
			  
			  if ( valid )
			    {
			      
			      const int aidx = anchor_idx[ i ];
			      
			      // track spindle anchor position (for output only)
			      spindles[i].anchor_sec = spindle_anchor[aidx] * globals::tp_duration ;
			      
			      if ( nearest_sw_number[ aidx ] != 0 ) // is now 1-basewd 
				{
				  spindles[i].so_nearest = nearest_sw[ aidx ] ;
				  spindles[i].so_nearest_num = nearest_sw_number[ aidx ] ;		      
				}
			      else
				spindles[i].so_nearest_num = 0;
			      
			      if ( sw_anchor[ aidx ] )
				spindles[i].so_phase_anchor = MiscMath::as_angle_0_pos2neg( ph_anchor[ aidx ] ) ;
			      else 
				spindles[i].so_phase_anchor = -9;
			      
			    }
			  else // for spindles w/ no valid anchors...
			    {
			      // set anchor to NA
			      spindles[i].anchor_sec = -9;
			      spindles[i].so_phase_anchor = -9;
			      spindles[i].so_nearest_num = 0; // as is 1-based for lookup 
			    }
			}
		    }
		}
	      

	      //
	      // only using complex anchor-statifier if >1 anchor specified
	      //

	      if ( sw_coupling_has_anchors ) 
		writer.unlevel( globals::anchor_strat );
	      
	      //
	      // End of COUPL analyses
	      //

	      
	      
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

		  double inc = 360 / (double)nbins;
		  double ph = inc/2.0; // use mid-point of range
		  std::vector<double> so_angle( nbins );


		  std::vector<int> pl_chirp_n ; 
		  std::vector<double> pl_chirp = p_sw->phase_locked_averaging( p_chirp_if , nbins , 
									       &in_spindle , &pl_chirp_n );
		  
		  // returned empty is no SW
		  if ( pl_chirp.size() == nbins )
		    {
		      for (int j=0;j<nbins;j++)
			{
			  writer.level( ph  , "PHASE" );
			  writer.value( "IF" , pl_chirp[j] );
			  writer.value( "N" , pl_chirp_n[j] );
			  so_angle[j] = clocs_t::deg2rad( ph );
			  ph += inc;
			}
		      writer.unlevel( "PHASE" );
		  		      
		      // Circular-linear correlation between SO phase and spindle frequency		  
		      double r_cl = Statistics::circular_linear_correlation( so_angle , pl_chirp );
		      writer.value( "R_PHASE_IF" , r_cl );
		    }
		  

		  //
		  // Repeat SO/IF analysis, but stratify by positive/negative HWs / slopes
		  //
		  
		  if ( 0 ) 
		    {

		      bool ran_typed = false;
		      
		      // POS HW:
		      
		      pl_chirp = p_sw->phase_locked_averaging( p_chirp_if , nbins , &in_pos_hw );		  
		      
		      if ( pl_chirp.size() == nbins )		    
			{
			  ran_typed = true;
			  writer.level( "POS_HW" , "TYPE" );
			  ph = inc/2.0; // use mid-point of range 
			  for (int j=0;j<nbins;j++)
			    {
			      writer.level( ph  , "PHASE" );
			      writer.value( "IF" , pl_chirp[j] );		      
			      ph += inc;
			    }
			  writer.unlevel( "PHASE" );
			  double r_cl = Statistics::circular_linear_correlation( so_angle , pl_chirp );
			  writer.value( "R_PHASE_IF" , r_cl );
			}
		      
		      // NEG HW:	  
		      
		      pl_chirp = p_sw->phase_locked_averaging( p_chirp_if , nbins , &in_neg_hw );		  
		      
		      if ( pl_chirp.size() == nbins)
			{
			  ran_typed = true;
			  writer.level( "NEG_HW" , "TYPE" );
			  ph = inc/2.0; // use mid-point of range
			  for (int j=0;j<nbins;j++)
			    {
			      writer.level( ph  , "PHASE" );
			      writer.value( "IF" , pl_chirp[j] );		      
			      ph += inc;
			    }
			  writer.unlevel( "PHASE" );
			  double r_cl = Statistics::circular_linear_correlation( so_angle , pl_chirp );
			  writer.value( "R_PHASE_IF" , r_cl );
			}
		      
		      // POS SLOPE:
		      
		      pl_chirp = p_sw->phase_locked_averaging( p_chirp_if , nbins , &in_pos_slope );		 
		      if ( pl_chirp.size() == nbins)
			{
			  ran_typed = true;
			  writer.level( "POS_SLOPE" , "TYPE" );
			  ph = inc/2.0; // use mid-point of range
			  for (int j=0;j<nbins;j++)
			    {
			      writer.level( ph  , "PHASE" );
			      writer.value( "IF" , pl_chirp[j] );		      
			      ph += inc;
			    }
			  writer.unlevel( "PHASE" );
			  double r_cl = Statistics::circular_linear_correlation( so_angle , pl_chirp );
			  writer.value( "R_PHASE_IF" , r_cl );
			}
		      
		      // NEG SLOPE:
		      
		      pl_chirp = p_sw->phase_locked_averaging( p_chirp_if , nbins , &in_neg_slope );		  
		      if ( pl_chirp.size() == nbins)
			{
			  ran_typed = true;
			  writer.level( "NEG_SLOPE" , "TYPE" );
			  ph = inc/2.0; // use mid-point of range
			  for (int j=0;j<nbins;j++)
			    {
			      writer.level( ph  , "PHASE" );
			      writer.value( "IF" , pl_chirp[j] );		      
			      ph += inc;
			    }
			  writer.unlevel( "PHASE" );
			  double r_cl = Statistics::circular_linear_correlation( so_angle , pl_chirp );
			  writer.value( "R_PHASE_IF" , r_cl );
			  
			}		  
		      
		      if ( ran_typed ) 
			writer.unlevel( "TYPE" );
		      
		    }
		  
		  
		  //
		  // time-locked SO spindle IF
		  //  in general, phase-locked analysis above should be sufficient		  
		  //  but include this for now
		  //

		  if ( globals::devel ) 
		    {
		      // +/- from negative peak 1 second
		      std::vector<double> tl_chirp = p_sw->time_locked_averaging( p_chirp_if , Fs[s] , 1 , 1 );
		      
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
		  // Spindle IF as a function of slow wave phase AND position in the spindle
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
		      
		      std::vector<int> pl_chirp_n ; 
		      std::vector<double> pl_chirp = p_sw->phase_locked_averaging( p_chirp_if , nbins , 
										   &in_spindle_and_bin , 
										   &pl_chirp_n );
		      
		      
		      if ( pl_chirp.size() == nbins ) 
			{
			  double inc = 360 / (double)nbins;
			  double ph = inc/2.0; // use mid-point of range
			  
			  std::vector<double> so_angle;
			  std::vector<double> pl_chirp2;
			  
			  for (int j=0;j<nbins;j++)
			    {
			      writer.level( ph  , "PHASE" );
			      writer.value( "IF" , pl_chirp[j] );
			      writer.value( "N" , pl_chirp_n[j] );
			      
			      if ( pl_chirp_n[j] > 0 )
				{
				  so_angle.push_back( clocs_t::deg2rad( ph ) );
				  pl_chirp2.push_back( pl_chirp[j] );
				}
			      
			      ph += inc;
			    }
			  writer.unlevel( "PHASE" );
			  
			  // cicular-linear correlation of PHASE and F (i.e. conditional on RELLOC bin)
			  if ( pl_chirp2.size() >= 5 )
			    {
			      double r_cl = Statistics::circular_linear_correlation( so_angle , pl_chirp );
			      if ( r_cl > -8 ) 
				writer.value( "R_PHASE_IF" , r_cl );
			    }
			  
			}
		    }		  
		  
		  // next spindle bin
		  writer.unlevel( "RELLOC" );
		  
		}
	    
	    
	    } // end of all SW-based options
	  
	  
	  
	  //
	  // Per-spindle level output
	  //

	  if ( show_spindle_level )
	    per_spindle_output( &spindles ,
				param ,
				( hms ? &starttime : NULL) , 
				q );
				


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
	  
	  if ( characterize && tlocking && spindles.size() > 0 )
	    {
	      
	      double tlock_min = locked.begin()->first;
	      double tlock_max = (--locked.end())->first;		  
		  
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

	  const bool has_coupling = param.has( "so" );
	  
	  while ( 1 ) 
	    {
	      
	      int epoch = edf.timeline.next_epoch();      
	      
	      if ( epoch == -1 ) break;
              
	      interval_t interval = edf.timeline.epoch( epoch );
	      
	      
	      const int nsp = spindles.size();
	      
	      int sp_epoch = 0;

	      // 'average' spindle_t for epoch-level stats
	      spindle_t avg(0,0,0,0);
	      avg.amp = avg.dur = avg.frq = avg.nosc = avg.isa = 0;
	      avg.symm = avg.symm2 = avg.chirp = avg.frq_h1 = avg.frq_h2 = 0;
	      avg.so_nearest = 0;
	      avg.norm_amp_max = avg.norm_amp_mean = 0;
	      	      
	      for (int i=0 ; i<nsp; i++)
		{
		  
		  // is this spindle included?
		  
		  if ( ! spindles[i].include ) continue;
		  
		  // dummy interval for just starting point of spindle 
		  // i.e. for purpose of assigning to EPOCH (so as not
		  // to double-count spindles that overlap epoch boundaries
		  
		  interval_t spstart( spindles[i].tp.start , spindles[i].tp.start );
		  
		  if ( interval.overlaps( spstart ) )
		    {
		      ++sp_epoch;

		      avg.amp += spindles[i].amp;
		      avg.dur += spindles[i].dur;
		      avg.frq += spindles[i].frq;
		      avg.nosc += spindles[i].nosc;
		      avg.isa += spindles[i].isa;
		      avg.symm += spindles[i].symm;
		      avg.symm2 += spindles[i].symm2;
		      avg.chirp += spindles[i].chirp;
		      avg.frq_h1 += spindles[i].frq_h1;
		      avg.frq_h2 += spindles[i].frq_h2;
		      avg.norm_amp_max += spindles[i].norm_amp_max;
		      avg.norm_amp_mean += spindles[i].norm_amp_mean;

		      if ( has_coupling )
			avg.so_nearest += spindles[i].so_nearest <= 1e-6; // i.e. PROP SO OVERLAP
					      
		    }
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

		  if ( sp_epoch )
		    {
		      
		      avg.amp /= (double)sp_epoch ;
		      avg.dur /= (double)sp_epoch ;
		      avg.frq /= (double)sp_epoch ;
		      avg.nosc /= (double)sp_epoch ;
		      avg.isa /= (double)sp_epoch ;
		      avg.symm /= (double)sp_epoch ;
		      avg.symm2 /= (double)sp_epoch;
		      avg.chirp /= (double)sp_epoch ;
		      avg.frq_h1 /= (double)sp_epoch ;
		      avg.frq_h2 /= (double)sp_epoch ;
		      avg.so_nearest /= (double)sp_epoch ;
		      avg.norm_amp_max /= (double)sp_epoch ;
		      avg.norm_amp_mean /= (double)sp_epoch ;

		      writer.value( "AMP" , avg.amp  );
		      writer.value( "ACT_MX" , avg.norm_amp_max  );
		      writer.value( "ACT_MN" , avg.norm_amp_mean  );
		      writer.value( "DUR" , avg.dur  );
		      writer.value( "FRQ" , avg.frq  );
		      writer.value( "NOSC" , avg.nosc );
		      writer.value( "ISA" , avg.isa );
		      writer.value( "SYMM" , avg.symm  );
		      writer.value( "SYMM2" , avg.symm2 );
		      writer.value( "CHIRP" , avg.chirp  );
		      writer.value( "FRQ1" , avg.frq_h1  );
		      writer.value( "FRQ2" , avg.frq_h2  );		      
		      if ( has_coupling	)
			writer.value( "COUPL_OVERLAP" , avg.so_nearest );

		      //
		      // track metrics dynamics (if observed 1+ spindle)
		      //

		      if ( calc_dynamics )
			{
			  const int e = edf.timeline.display_epoch( epoch ) - 1;
			  
			  qd.add( writer.faclvl_notime() , "AMP" , e  , avg.amp );
			  qd.add( writer.faclvl_notime() , "ACT_MX" , e  , avg.norm_amp_max );
			  qd.add( writer.faclvl_notime() , "ACT_MN" , e  , avg.norm_amp_mean );
			  qd.add( writer.faclvl_notime() , "DUR" , e  , avg.dur );
			  qd.add( writer.faclvl_notime() , "FRQ" , e  , avg.frq );
			  qd.add( writer.faclvl_notime() , "NOSC" , e  , avg.nosc );
			  qd.add( writer.faclvl_notime() , "ISA" , e  , avg.isa );
			  qd.add( writer.faclvl_notime() , "SYMM" , e  , avg.symm );
			  qd.add( writer.faclvl_notime() , "SYMM2" , e  , avg.symm2 );
			  qd.add( writer.faclvl_notime() , "CHIRP" , e  , avg.chirp );
			  qd.add( writer.faclvl_notime() , "FRQ1" , e  , avg.frq_h1 );
			  qd.add( writer.faclvl_notime() , "FRQ2" , e  , avg.frq_h2 );
			  if ( has_coupling ) 
			    qd.add( writer.faclvl_notime() , "COUPL_OVERLAP" , e  , avg.so_nearest );
			  
			}
		      
		    }
		  

		  // N can be tracked if 0 
		  if ( calc_dynamics )
		    {
		      const int e = edf.timeline.display_epoch( epoch ) - 1;
		      qd.add( writer.faclvl_notime() , "N" , e  , sp_epoch );		      
		    }
		  
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
	      
	      writer.var( "DISPERSION" ,"Spindle epoch-dispersion index" );
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
	      writer.var( "AMP"   , "Mean spindle amplitude" );
	      writer.var( "DUR"   , "Mean spindle duration" );
	      writer.var( "FWHM"  , "Mean spindle FWHM" );
	      writer.var( "NOSC"  , "Mean spindle number of oscillations" );
	      writer.var( "FRQ"   , "Mean spindle frequency (zero-crossing method)" );
	      writer.var( "FFT"   , "Mean spindle frequency (FFT)" );	      
	      writer.var( "SYMM"  , "Mean spindle symmetry index" );
	      writer.var( "SYMM2"  , "Mean spindle folded symmetry index" );
	      writer.var( "CHIRP" , "Mean spindle chirp index" );

	      writer.value( "CWT_TH" , mean ); // mean or median
	      if ( winsorize ) 
		{
		  writer.value( "WIN_TH" , emp_winsorize ); // 
		  writer.value( "WIN1" , means[ "W_ANY" ] );
		  writer.value( "WIN" , means[ "W" ] );
		}
	      
	      writer.value( "N01" , nspindles_premerge );  // original
	      writer.value( "N02" , nspindles_postmerge ); // post merging
	      writer.value( "N" ,  (int)spindles.size()  ) ;    // post merging and QC	    

	      writer.value( "P01" , (int)spindles.size() / (double)nspindles_premerge ); // relative final vs orig (impact of merge)
	      writer.value( "P02" , (int)spindles.size() / (double)nspindles_postmerge ); // reltive final vs post-merge (impact of Q)
	      
	      writer.value( "MINS" , t_minutes );
	      writer.value( "DENS" , spindles.size() / t_minutes );
	      
	      writer.value( "ISA_S" , means[ "ISA_PER_SPINDLE" ] );
	      writer.value( "ISA_M" , means[ "ISA_TOTAL" ] / t_minutes );
	      writer.value( "ISA_T" , means[ "ISA_TOTAL" ] );
	      writer.value( "Q"     , means[ "Q" ] );

	      
	      writer.value( "AMP" , means["AMP"] );
	      writer.value( "ACT_MX" , means["ACT_MX"] );
	      writer.value( "ACT_MN" , means["ACT_MN"] );
	      writer.value( "DUR" , means["DUR"] );
	      writer.value( "FWHM" , means["FWHM"] );
	      writer.value( "NOSC" , means["NOSC"] );

	      writer.value( "SYMM" , means["SYMM"] );
	      writer.value( "SYMM2" , means["SYMM2"] );
	      writer.value( "SYMM_AMP" , means["SYMM_AMP"] );
	      writer.value( "SYMM_TROUGH" , means["SYMM_TROUGH"] );
	      writer.value( "SEC_P2P" , means["SEC_P2P"] );
	      writer.value( "SEC_AMP" , means["SEC_AMP"] );
	      writer.value( "SEC_TROUGH" , means["SEC_TROUGH"] );
	      
	      // spindle frequency
	      writer.value( "FRQ" , means["FRQ"] );
	      writer.value( "FFT" , means["FFT"] );	      
	      writer.value( "CHIRP" , means["CHIRP"] );	      
	      writer.value( "FVAR2" , means["FVAR2"] ); // spindle variance
	      writer.value( "FRNG2" , means["FRNG2"] ); // spindle freq range
	      writer.value( "FRQ1" , means["FRQ1"] );
	      writer.value( "FRQ2" , means["FRQ2"] );
	      
	      
	      if ( globals::devel )
		{
		  
		  //
		  // positive spindle components
		  //

		  writer.level( "POS" , "SUBSET" );
		  writer.value( "F_HW" , means["FPOS"] );
		  writer.value( "F_S" , means["FSPOS"] );
		  writer.value( "B_HW" , means["BPOS"] );
		  writer.value( "B_S" , means["BSPOS"] );
		  writer.value( "V" , means["VPOS"] );
		  const double pos_amp = means[ "POSISA_PER_SP" ] ;
		  writer.value( "AMP" , pos_amp );

		  
		  //
		  // negative spindle components
		  //
		  
		  writer.level( "NEG" , "SUBSET" );
		  writer.value( "F_HW" , means["FNEG"] );
		  writer.value( "F_S" , means["FSNEG"] );
		  writer.value( "B_HW" , means["BNEG"] );
		  writer.value( "B_S" , means["BSNEG"] );
		  writer.value( "V" , means["VNEG"] );		  
		  const double neg_amp = means[ "NEGISA_PER_SP" ] ;
		  writer.value( "AMP" , neg_amp );

		  
		  //
		  // Avereged (i.e. all, based on HWs)
		  //

		  writer.level( "ALL" , "SUBSET" );
		  writer.value( "F_HW" , means["FALL"] );
		  writer.value( "B_HW" , means["BALL"] );
		  writer.value( "V" , means["VALL"] );
		  
		  //
		  // Specific difference contrasts (POS vs NEG)
		  //

		  writer.level( "DIF" , "SUBSET" );
		  writer.value( "F_HW" , means["FPOS"] - means["FNEG"] ); 
		  writer.value( "F_S" , means["FSPOS"] - means["FSNEG"] ); 		  
		  writer.value( "B_HW" , means["BPOS"] - means["BNEG"] ); 		  
		  writer.value( "B_S" , means["BSPOS"] - means["BSNEG"] ); 		  
		  writer.value( "V" , means["VPOS"] - means["VNEG"]  );		  
		  writer.value( "AMP" , pos_amp - neg_amp );
		  
		  writer.unlevel( "SUBSET" );
		}
	      
	      writer.value( "DISPERSION" , means[ "DISPERSION" ] );
	      writer.value( "DISPERSION_P" , means[ "DISPERSION_P" ] );
	      writer.value( "NE" , means[ "NE" ] );
	      
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
	      
	      annot_t * a = edf.annotations->add( aname );
	      a->description = "Spindle intervals";
	      
	      logger << "  creating annotation class: " << aname 
		     << ", instance: " << analysis_label 
		     << ", channel: " << signals.label(s) << "\n";
	      

	      // use F_C as instance label
	      for (int i=0;i<spindles.size();i++)
		{
		  
		  const spindle_t & spindle = spindles[i];
		  
		  instance_t * instance = a->add( analysis_label , spindle.tp , signals.label(s) );

		  // any meta-data?
		  
		  instance->set( "amp" , spindle.amp );
		  instance->set( "isa" , spindle.isa );
		  instance->set( "frq" , spindle.frq );
		  instance->set( "dur" , spindle.dur );

		  if ( param.has( "so" ) )
		    instance->set( "soc" , (int)( fabs( spindle.so_nearest < 1e-8 ) ) );
		  
		  // index tp
		  //instance->set( "mid", "tp:" + Helper::int2str( spindle.tp_mid ) );
		  instance->set( "rp_mid", (double)( spindle.tp_mid - spindle.tp.start )
				 / double( spindle.tp.stop - spindle.tp.start ) );
		  
		}
	      
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
      // report dynamics for this signal
      //
      
      if ( calc_dynamics )
	qd.proc_all();

      
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

      props.regional( 1 );
      
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
	      avgs.push_back( props.analyse( f , c , signals.label(s) , w , verbose_prop ) );
	      
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
			    const std::string & analysis_label ,
			    const std::vector<double> * averaged ,
			    const std::vector<double> * original_signal , 
			    std::vector<spindle_t>    * spindles ,
			    clocktime_t               * starttime , 
			    spindle_qc_t * q , 
			    std::map<double,double>   * locked ,
			    std::vector<bool> * p_in_pos_hw , 
			    std::vector<bool> * p_in_neg_hw , 
			    std::vector<bool> * p_in_pos_slope , 
			    std::vector<bool> * p_in_neg_slope ,
			    std::vector<bool> * p_winsp 
			    )   // input, is 
{
 

  //
  // Filtering parameters
  //

  
  // default, bandpass width 4 Hz (i.e. +/- 2 Hz around the peak frequency
  const double window_f = param.has( "char-bp" ) ? param.requires_dbl( "char-bp" ) : 4 ;

  // ripple and transition width
  const double window_ripple = param.has( "char-ripple" ) ? param.requires_dbl( "char-ripple" ) : 0.02 ;
  const double window_tw = param.has( "char-tw" ) ? param.requires_dbl( "char-tw" ) : 4 ;
  const int hann = param.has( "hann" ) ? param.requires_int( "hann" ) : 0 ; 

  std::vector<double> window_ripple_vec( 1 , window_ripple );
  std::vector<double> window_tw_vec( 1 , window_tw );

  
  //
  // Create a copy of this signal, if it does not already exist
  //
  
  const std::string signal_label = edf.header.label[s0];
  const std::string new_label = signal_label + "_BP" + "_" 
    + Helper::dbl2str( target_f ) + "_" + Helper::dbl2str( window_f );

  if ( ! edf.header.has_signal( new_label ) )
    {

      // copy the existing signal

      edf.copy_signal( edf.header.label[s0] , new_label );
      
      const int s = edf.header.signal( new_label );

      // and do we need to band-pass filter this new signal?
      
      if ( ! bandpass_filtered ) 
	{

	  //
	  // Kaiser window 
	  //

	  if ( hann == 0 ) 
	    {
	      
	      logger << "  filtering at " 
		     << target_f - window_f * 0.5  << " to " 
		     << target_f + window_f * 0.5 << ", ripple=" 
		     << window_ripple << " tw=" << window_tw << "\n";
	      
	      // default above is 4 Hz, i.e. +/- 2 Hz
	      //  ~ 9-13Hz for slow spindles,    13-17Hz for fast spindles
	      
	      // defaults: 
	      //  char-ripple = 0.02 char-tw=4 char-bp=4  (i.e. bandpass=9,13 )
	      //  char-ripple = 0.02 char-tw=4 char-bp=4  (i.e. bandpass=13,17 )
	      
	      dsptools::apply_fir( edf , s , fir_t::BAND_PASS ,
				   1 , // 1 = Kaiser window
				   window_ripple_vec , window_tw_vec , 
				   target_f - window_f * 0.5 ,
				   target_f + window_f * 0.5 ,
				   0, // order (if not Kaiser win)
				   fir_t::RECTANGULAR , // large window, no window
				   true  // use FFT convolution
				   );
	    }
	  
	  //
	  // Hann window 
	  //

	  else
	    {
	      
	      logger << "  filtering at "
		     << target_f - window_f * 0.5  << " to "
		     << target_f + window_f * 0.5 << ", Hann window, order = " << hann << "\n";

	      std::vector<double> zero;
	      
	      dsptools::apply_fir( edf , s , fir_t::BAND_PASS ,
				   2 , // fixed order
				   zero , zero , // Kaiser stuff 
				   target_f - window_f * 0.5 ,
				   target_f + window_f * 0.5 ,
				   hann, // --> even order: odd / symmetric FIR
				   fir_t::HANN , // Hann window
				   true  // use FFT convolution
				   );
	    }
	  
	}
      
    }
  
  const int s = edf.header.signal( new_label );
  
 
  //
  // Iterate over each spindle
  //

  const int n = spindles->size();
  
  // track if we QC any spindles out at this step
  bool removed_some = false;

   for (int i=0;i<n;i++)
     {

       spindle_t * spindle = &(*spindles)[ i ];
       
       
       //
       // pull out band-pass filtered data for actual spindle
       //
       
       slice_t slice( edf , s , spindle->tp );

       // nb. we are going to mean-centre these data too
       // for this particular spindle... avoids some weird cases
       // where we have huge offsets in original signal amplitudes around
       // this spindle, meaning that even after BP filter, we can
       // have spindle intervals with no zero-crossings
       
       std::vector<double> d = MiscMath::centre( *slice.pdata() );

       std::vector<uint64_t> tp = *slice.ptimepoints();
       
       const int Fs = edf.header.sampling_freq( s );
       
       const double period_sec = 1.0/(double)Fs;

       const uint64_t period = period_sec * globals::tp_1sec;

       const int npoints = d.size();


       // Sanity check --> next spindle
       if ( npoints < 2 ) 
	 continue;	   
	 

       //
       // duration
       //
       
       spindle->dur = (spindle->tp.stop - spindle->tp.start ) / (double)globals::tp_1sec;
      

       //
       // ISA (scale by SR) and normalized ACT and peak ampl position (sp-offset)
       //
       
       spindle->isa = 0;
       spindle->norm_amp_max = 0;
       spindle->peak_amp_sp = 0;
       spindle->peak_amp_rel = -1;
       
       if ( averaged != NULL )
	 {
	   
	   int start = spindle->start_sp;
	   int stop  = spindle->stop_sp;
	   
	   for (int s=start;s<=stop;s++)
	     {
	       spindle->isa += (*averaged)[s] ;

	       // track max (ACT_MX)
	       if ( (*averaged)[s] > spindle->norm_amp_max )
		 {
		   spindle->norm_amp_max = (*averaged)[s];
		   spindle->peak_amp_sp = s - start ; // 0-based sp offset                    
       		 }
	     }

	   // ACT_MN (i.e. simple mean, no duration info) 
	   spindle->norm_amp_mean = spindle->isa / (double)(stop-start+1);
	   
	   // ISA_S (i.e. includes duration info as norm by Fs only
	   spindle->isa /= (double)Fs;

	   // relative position of peak-amp
	   spindle->peak_amp_rel = spindle->peak_amp_sp / (double)( stop - start );
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
       // spindle winsorization status
       //
       
       int start = spindle->start_sp;
       int stop  = spindle->stop_sp;

       // marker for no winsorization done = -1
       spindle->winsor = -1; 

       if ( p_winsp != NULL )
	 {
	   spindle->winsor = 0;	   
	   for (int s=start;s<=stop;s++)
	     if ( (*p_winsp)[s] ) spindle->winsor++; 
	 }
     

       
      //
      // Find largest peak-to-peak amplitude
      //
       
      // rule of thumb: call something a peak if the surrounding +/i 2
      // points are all smaller, or all larger; for tied values, extend to the next 
      // variable ones, and track only the last one as the peak
      
      std::vector<int> peak;
      
      
      // track signal max/min and 5% of max/min
      // of the filtered signal

      double fltmax = 0 , fltmin = 0; 
      

      // 0 nothing; 2 pos peak,  -2 negative peak
      //            1 neg2pos ZC -1 pos2neg ZC
      
      // below, we define halfwaves both :: ZC to ZC     (halfwaves)
      //                                    peak-to-peak (slopes)

      std::vector<int> landmark( npoints , 0 );

      for (int p=2;p<npoints-2;p++)
	{
	  
	  // tied w/ the next point? 
	  if ( d[p] == d[p+1] ) continue;

	  // track max/min
	  if ( d[p] < fltmin ) fltmin = d[p];
	  else if ( d[p] > fltmax ) fltmax = d[p];
	  
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
	  
	  if ( gt == 4 ) 
	    {
	      peak.push_back(p); 
	      landmark[p] = 2;
	    }	  
	  else if ( lt == 4 ) 
	    {
	      peak.push_back(p);
	      landmark[p] = -2;
	    }

	}
      
      // something strange?  bail
      if ( peak.size() < 2 ) 
	{
	  logger << " *** warning: spindle w/ only a single peak... should not happen... bailing" << "\n";
	  continue;
	}


      const bool ht_verb = false;

      if ( ht_verb ) 
	{
	  std::cout << "\n\nSPINDLE " << i+1 << " npoints = "
		    << npoints << " ::: " <<  spindle->start_sp << " " <<  spindle->stop_sp << "\n";
	  for (int p=0;p<npoints;p++)
	    std::cout << p << "\t" <<  p * period_sec << "\t" << d[p] << "\t" << (*original_signal)[ spindle->start_sp + p] << "\n";	  
	}


     
      
      //
      // Zero-crossings, in seconds, with linear interpolation between points
      //

      std::vector<double> zc;  // duration of i to i+1 zero-crossing (i.e. half-waves)
      std::vector<bool> zcp;   // T : pos-neg; F : neg-pos
      std::vector<int> zcs;    // track which sample point (is before the ZC) 

      for (int p=0;p<npoints-1;p++)
        {

	  const bool pos2neg = d[p] >= 0 && d[p+1] < 0 ;
	  const bool neg2pos = d[p] <= 0 && d[p+1] > 0 ;	  
	  if ( ! ( pos2neg || neg2pos ) ) continue;
	  
	  // set landmarks
	  if ( pos2neg ) landmark[p] = -1;
	  else landmark[p] = +1;
	  
	  // linear interpolation between two ZC points
	  const double s1 = p * period_sec ;
	  const double frac = fabs( d[p] ) / ( fabs( d[p] ) + fabs( d[p+1] ) );
	  
	  zc.push_back( s1 + frac * period_sec );
	  zcp.push_back( pos2neg );
	  zcs.push_back( p );
	}


      // no zero-crossings?
      if ( zc.size() == 0 )
        {
          logger << " *** warning: spindle w/out any zero-crossings... bailing on this spindle\n";
	  

	  // slice_t slice( edf , s0 , spindle->tp );
	  // std::vector<double> d2 = *slice.pdata();
	  // if ( d2.size() != npoints ) Helper::halt( "hmm" );
	  // for (int jj=0; jj<npoints; jj++)
	  //   std::cout << " " << d[jj] << "\t" << d2[jj] << "\n";
	  // std::cout << "\n\n";

	  continue;
        }


      if ( ht_verb )
	{
	  std::cout << "\nINTERPOLATED ZCs (dur. from i to i+1 ZC):\n";
	  for (int zz=0; zz<zc.size(); zz++)
	    std::cout << zz << "\t" << zc[zz] << "\t" << ( zcp[zz] ? "P2N" : "N2P" ) << " nearest prior SP " << zcs[zz] << "\n"; 	      
	}
      

      //
      // Positive peaks, with cubic interpolation at arbitrary temporal resolution
      //
      
      // time_points
      std::vector<double> t0( npoints );
      for (int p=0; p<npoints; p++) t0[p] = p * period_sec;
      
      tk::spline spline;
      spline.set_points( t0 , d );

      const double tdur = (npoints-1.0) * period_sec  ;
      const double tinc = tdur * 0.001;

      // spline interpolation at this resolution
      std::vector<double> interp;
      for (double tt = 0 ; tt <= tdur ; tt += tinc )
	interp.push_back( spline(tt) );


      if ( ht_verb && 0 ) 
	{
	  std::cout << "\n1000-point interpolated sequence (for peak-to-trough finding)\n";
	  for (int p=0;p<interp.size(); p++) 
	    std::cout << "INTERP\t" << p << "\t" << p * tinc << "\t" << interp[p] << "\n";
	}

      // find peaks/troughs

      const int nint = interp.size();
      std::vector<double> tpk;
      std::vector<bool> ppk; // peak / trough 
      std::vector<double> mpk;  // peak magnitude 
      double int_fltmin = 0, int_fltmax = 0; // track min/max 

      for (int p=2; p < nint-2;p++)
	{
	  
	  // tied w/ the next point? 
	  if ( interp[p] == interp[p+1] ) continue;

	  // peaks  	  
	  int gt = 0 , lt = 0;
	  
	  // forwards
	  if      ( interp[p] < interp[p+1] ) ++lt;
	  else if ( interp[p] > interp[p+1] ) ++gt;
	  
	  if      ( interp[p] < interp[p+2] ) ++lt;
	  else if ( interp[p] > interp[p+2] ) ++gt;
	  
	  // backwards, which might require we skip ties previously skipped
	  int bck = 0;
	  bool skip = false;
	  while ( 1 ) 
	    {
	      if ( p - bck < 1 ) { skip = true; break; } 
	      if ( interp[p] == interp[p-bck] ) ++bck;
	      else break;
	    }

	  if ( skip ) continue;
	  
	  if      ( interp[p] < interp[p-bck] ) ++lt;
	  else if ( interp[p] > interp[p-bck] ) ++gt;
	  
	  if      ( interp[p] < interp[p-bck-1] ) ++lt;
	  else if ( interp[p] > interp[p-bck-1] ) ++gt;
	  
	  if ( gt == 4 ) 
	    {	      
	      tpk.push_back( p * tinc ); 
	      ppk.push_back( true );
	      mpk.push_back( interp[p]  );
	      if ( interp[p] > int_fltmax ) int_fltmax = interp[p] ;
	    }
	  else if ( lt == 4 ) 
	    {
	      tpk.push_back( p * tinc ); 
	      ppk.push_back( false );	      
	      mpk.push_back( interp[p]  );
	      if ( interp[p] < int_fltmin ) int_fltmin = interp[p] ;
	    }
	}
      
      
      //
      // Get intervals from peak-trough, vice-versa
      //

      // want to get two vectors or intervals: for pos-slopes and neg-slopes
      std::vector<double> wpos_slope, wneg_slope;  // duration of slope interval 
      std::vector<double> tpos_slope, tneg_slope;  // time mid-points of slope interval
      
      // for (int p=0; p<tpk.size() ; p++)
      // 	std::cout << "peak = " << p << "/" << tpk.size() << " " << ppk[p] << "\t" << tpk[p] << "\n";

      // by default. needs to be 5% of max peak-to-peak
      const double req_p2p = ( int_fltmax - int_fltmin ) * 0.05;

      for (int p=0; p<tpk.size() - 1 ; p++)
	{

	  // skip if consecutive peaks of the same polarity
	  //  i.e. typically local minima, occassionally happens
	  if ( ppk[p] == ppk[p+1] ) continue;
	  
	  // also, skip if either range is not > % of max peak-to-peak
	  if ( fabs( mpk[p] ) + fabs( mpk[p+1] ) < req_p2p ) 
	    {
	      //std::cerr << "*****************************************************SKIPPY\n";
	      continue;
	    }

	  bool pos_slope = ! ppk[p] ; // from NEG to POS

	  const double w1 = tpk[p];
          const double w2 = tpk[p+1];

          const double w = 1.0 / ( 2 * ( w2 - w1 ) );
          const double t = ( w1 + w2 ) / 2.0 ;

	  // tracsk diff and midpoint
	  if ( pos_slope )
	    {
	      wpos_slope.push_back( w );
	      tpos_slope.push_back( t );
	    }
	  else
	    {
	      wneg_slope.push_back( w );
	      tneg_slope.push_back( t );
	    }
	}


      
      //
      // populate boolean vectors for use in IF analyses? (if 'if' is set) 
      //
      
      if ( p_in_pos_hw != NULL ) 
	{
	  
	  // determine HW status just by polarity of original signal
	  // use landmarks to get slopes
	  
	  // get first landmark to figure out slope
	  bool pos_slope = false;
	  
	  for (int p=0; p<npoints; p++)
	    {
	      if ( landmark[p] == 2 || landmark[p] == -2 ) 
		{
		  pos_slope = landmark[p] == 2;
		  break;
		}
	    }

	  
          int sp = spindle->start_sp;

	  for (int p=0; p<npoints; p++)
	    {

	      // HW status
	      
	      const bool pos_hw = d[p] >= 0;

	      // mask according to state
	      if ( pos_hw ) 
		(*p_in_neg_hw)[sp] = false;
	      else
		(*p_in_pos_hw)[sp] = false;

	      // determine slope
	      if ( landmark[p] == 2 ) pos_slope = false;
	      else if ( landmark[p] == -2 ) pos_slope = true;

	      // mask according to state
	      if ( pos_slope ) (*p_in_neg_slope)[sp] = false;
	      else (*p_in_pos_slope)[sp] = false;
	      	      
	      // std::cout << " det " << i << "\t" << d[p] << " " << landmark[p] << " " <<  (*p_in_pos_hw)[sp] <<  (*p_in_neg_hw)[sp] << " " <<  (*p_in_pos_slope)[sp] <<  (*p_in_neg_slope)[sp] << "\n";

	      // move along spindle sample point
	      ++sp;

	    }
	  //std::cout << "\n npoints = " << npoints << "\n";

	}


      //
      // ISA stratified by pos/neg HWs (currently, yoked to 'if' option)
      //
      
      spindle->posisa = spindle->negisa = 0;
      spindle->possp = spindle->negsp = 0;

      if ( p_in_pos_hw != NULL &&  averaged != NULL )
	{
	  
	  int start = spindle->start_sp;
	  int stop  = spindle->stop_sp;
	  
	  for (int s=start;s<=stop;s++)
	    {
	      if ( (*p_in_pos_hw)[s] )
		{
		  spindle->posisa += (*averaged)[s] ;
		  ++(spindle->possp);
		  // std::cout << " adding pos " << (*averaged)[s]  << "  " << spindle->posisa << " " << spindle->negisa << " / " 
		  // 	    << spindle->possp << " " << spindle->negsp << "\n";
		}
	      else
		{
		  spindle->negisa += (*averaged)[s] ;
		  ++(spindle->negsp);

		  // std::cout << " adding neg " << (*averaged)[s]  << "  " << spindle->posisa << " " << spindle->negisa << " / " 
		  // 	    << spindle->possp << " " << spindle->negsp << "\n";

		}
	    }
	  
	  spindle->posisa /= (double)Fs;
	  spindle->negisa /= (double)Fs;

	  spindle->posisa /= (double)spindle->possp;
	  spindle->negisa /= (double)spindle->negsp;

	}

		  
      //      std::cout << " spindle->posisa = " << spindle->posisa << " " << spindle->negisa << "\n";

      //
      // pos/neg halfwaves freqs; we check to make sure that there
      // is a sufficient peak between ZCs here: e.g. X% of the max
      //
      
            
      const double fltmin_th = 0.05 * fltmin;
      const double fltmax_th = 0.05 * fltmax;
      
      if ( ht_verb ) 
	{	  
	  std::cout << "\nEXTRACTING HWs for F_* and B_* calculations\n";
	  std::cout << "fltmin_th (5%) = " << fltmin_th << "\n"; 
	  std::cout << "fltmax_th (5%) = " << fltmax_th << "\n"; 
	}

      // frequency of halfwaves
      std::vector<double> wpos, wneg, wall; 
      std::vector<double> tpos, tneg, tall; 
      
      for (int z=0; z<zc.size()-1; z++)
	{

	  //	  std::cout << "zcs = " << zcs.size() << " " << zcp.size() << " " << " and " << z << "\n";
	  
	  // neg HW (i.e. initiated by post->neg ZC)
	  const bool neg_halfwave = zcp[z];

	  // is the spanned peak sufficient?
	  bool okay = false;
	  
	  for (int p=zcs[z]; p <= zcs[z+1]; p++)
	    {	      
	      if ( neg_halfwave ) 
		{
		  if ( d[p] < fltmin_th ) { okay = true; break; }
		}
	      else
		{
		  if ( d[p] > fltmax_th ) { okay = true; break; }
		}
	    }
	  
	  
	  // skip?
	  if ( ! okay ) 
	    {
	      if ( ht_verb ) 
		std::cout << "  *** skipping " << ( neg_halfwave ? "neg" : "pos" ) << " HW from (1-based) " << z+1 << " to " << z+2 << " ZCs\n";
	      continue;
	    }
	  
	  // transform to frequency 1/SW
	  double w = 1.0 / ( 2 * ( zc[z+1] - zc[z] ) );
	  double t = ( zc[z] + zc[z+1] ) / 2.0;
	  
	  if ( neg_halfwave )
	    {
	      if ( ht_verb ) std::cout << " +++ adding ZC-defined (NEG HW) = " << zc[z] << " to " << zc[z+1] << " = " << w << " " << t << "\n"; 
	      wneg.push_back(w);
	      tneg.push_back(t);
	    }
	  else
	    {	      
	      if ( ht_verb ) std::cout << " +++ adding ZC-defined (POS HW) = " << zc[z] << " to " << zc[z+1] << " = " << w << " " << t << "\n"; 
	      wpos.push_back( w );
	      tpos.push_back( t );
	    }
	  
	  wall.push_back( w );
	  tall.push_back( t );
	}


      //
      // F min/max
      //

      double fmin = wall[0];
      double fmax = wall[0];
      for (int z=1; z<wall.size(); z++)
	{
	  if ( wall[z] < fmin ) fmin = wall[z];
	  else if ( wall[z] > fmax ) fmax = wall[z];
	}
      
      spindle->frq_range = fmax - fmin;


      //
      // Durations/freqs of slopes -- no need to check things here
      //
          
      std::vector<double> wpos2, wneg2;
      std::vector<double> tpos2, tneg2; 
      
      for (int p=0; p<peak.size()-1; p++)
	{
	  // neg-peak to pos-peak :
	  const bool pos_slope = landmark[ peak[p] ] == -2 ; 
	  
	  // slope duration, in Hz
	  const double w1 = peak[p] * period_sec;
	  const double w2 = peak[p+1] * period_sec;

	  const double w = 1.0 / ( 2 * ( w2 - w1 ) );
	  const double t = ( w1 + w2 ) / 2.0 ; 
	  
	  // neg-halfwave
	  if ( pos_slope )
	    {
	      wpos2.push_back(w);  
	      tpos2.push_back(t);
	    }
	  else
	    {
	      wneg2.push_back(w);	      
	      tneg2.push_back(t);
	    }

	}


      //
      // Slope of frequency (both HWs and slopes) stratified by pos/neg polarity
      //

      // Y = implied frequency
      // X = time midpoint of halfwave/slope

      dynam_t zall( wall , tall ); 
      dynam_t zpos( wpos , tpos ); 
      dynam_t zneg( wneg , tneg ); 

      // if ( i == 10 ) 
      // 	{
      // 	  for (int ww=0; ww<wpos.size(); ww++)
      // 	    std::cout << " pos = " << tpos[ww] << " = " << wpos[ww] << "\n";
	  
      // 	  for (int ww=0; ww<wneg.size(); ww++)
      // 	    std::cout << " neg = " << tneg[ww] << " = " << wneg[ww] << "\n";
      // 	}


      zall.linear_trend( &spindle->allb , NULL );
      zpos.linear_trend( &spindle->posb , NULL );
      zneg.linear_trend( &spindle->negb , NULL );
      
      zall.mean_variance( &spindle->allf , &spindle->allv );
      zpos.mean_variance( &spindle->posf , &spindle->posv );
      zneg.mean_variance( &spindle->negf , &spindle->negv );

      if ( ht_verb )
	{
	  std::cout << "\nFINAL num INTERVALS = POS, NEG, ALL " << zpos.size() << "\t" << zneg.size() << "\t" << zall.size() << "\n";
	  std::cout << "\nFINAL F_POS = " << spindle->posf << "   and F_NEG = " << spindle->negf << " F_ALL = " << spindle->allf << "\n";
	  std::cout << "\nFINAL B_POS = " << spindle->posb << "   and B_NEG = " << spindle->negb << " B_ALL = " << spindle->allb << "\n";
	}

      // if ( i == 10 ) { 
      // 	std::cout << " MEAN = F_POS, F_NEG = " <<  spindle->posf << " " <<  spindle->negf << "\n";
      // }

      // slopes
      
      // old - based on samples
      // dynam_t zpos2( wpos2 , tpos2 ); 
      //  dynam_t zneg2( wneg2 , tneg2 ); 
       
      // new - use interpolated variants
      dynam_t zpos2( wpos_slope , tpos_slope ); 
      dynam_t zneg2( wneg_slope , tneg_slope ); 

      // ignore variance for the slopes (for now)
      zpos2.linear_trend( &spindle->pos2b , NULL );
      zneg2.linear_trend( &spindle->neg2b , NULL );
      
      zpos2.mean_variance( &spindle->pos2f , NULL );
      zneg2.mean_variance( &spindle->neg2f , NULL );
          

      // VERB
      // std::cout << "\n\n .... \n\n";
      // if ( spindle->posf > 1000 ) 
      // 	{

      // 	  std::cout << "d--> " << d.size() << " " << original_signal->size() << "\n";
	  
      // 	  for (int l=0; l<npoints; l++)
      // 	    {
      // 	      std::cout << l << "\t" << landmark[l] << "\t" << d[l] << " " << (*original_signal)[ spindle->start_sp + l] << "\n";
      // 	    }

      // 	  std::cout << " ZC : neg\n";
      // 	  for (int l=0;l<zneg.size();l++) std::cout << wneg[l] << " " << tneg[l] << "\n";

      // 	  std::cout << " ZC : pos\n";
      // 	  for (int l=0;l<zpos.size();l++) std::cout << wpos[l] << " " << tpos[l] << "\n";
	  
      // 	  std::cout << " DONE>...\n";
      // 	}


      
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
      
      spindle->chirp = -99999;
      spindle->frq_h1 = spindle->frq_h2 = 0;

      bool valid_chirp = cint1 > 1 && cint2 > 1 ; 

      if ( valid_chirp )
	{
	  // go mean from peak-to-peak duration in sample points, to implied frequency, Hz
	  spindle->frq_h1 = 1.0 / ( 2 * ( period_sec * int1/(double)cint1 ) );
	  spindle->frq_h2 = 1.0 / ( 2 * ( period_sec * int2/(double)cint2 ) );

	  // +ve means getting faster: absolute diffference (Hz)
	  spindle->chirp = spindle->frq_h2 - spindle->frq_h1 ; 
	  
	  // old CHIRP definition: log scaled ratio
	  //spindle->chirp = log( ( int1/(double)cint1 )  / (int2/(double)cint2  )  );
	  
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
      // 0-based sp offset, similar to spindle->peak_amp_sp
      spindle->max_trough_sp = lowest_idx; 
      
      // track tp of spindle max neg-trough 
      spindle->max_trough_rel = spindle->max_trough_sp / (double)( spindle->stop_sp - spindle->start_sp ) ;      
      spindle->tp_mid = spindle->tp.start + spindle->max_trough_rel * ( spindle->tp.stop - spindle->tp.start ) ;  
      
      // spindle symmetry (based on mid-point of largest peak-to-trough)
      spindle->symm = max_p2p_idx;
      
      // folded symmetry index (i.e. 0 = mid-way; 1 = 0 )
      spindle->symm2 = 2.0 * fabs( spindle->symm - 0.5 ) ;

      // spindle amp (max peak-to-peak)
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
      // Q stats? 
      //
      
      if ( q->qc_q ) 
	{
	  
	  // get Q-score
	  spindle->qual = q->qual( spindle );
	  
	  // failed?
	  if ( spindle->qual < q->qc_qmin ) spindle->include = false;
	  if ( q->qc_qmax > 0 && spindle->qual > q->qc_qmax ) spindle->include = false;
	  
	  if ( q->verbose_all || ( q->verbose_failed && ! spindle->include ) ) 
	    q->output( spindle , i  );
	}
     
          
      
      //
      // [ REMOVE - i.e. do generically ] Optional, time-locked analysis?  [ for QC+ spindles only ] 
      //
      
      if ( locked && spindle->include )
	{
	  
	  // use original signal, plus a window (+/-) 1.5 seconds around center)
	  const double window_sec = 2.0;
	  const int    nbins      = window_sec * Fs ;
	  const double window_left = - ( window_sec / 2.0 );

	  // start point (left of window) for peak minus half window
	  int orig_sp = spindle->start_sp + lowest_idx - ( window_sec / 2.0 ) * Fs ;	  
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
   // Prune spindle list based on Q ?
   //

   if ( removed_some ) 
     {
       std::vector<spindle_t> copy_spindles = *spindles;
       spindles->clear();
       for (int i=0;i<copy_spindles.size();i++)
	 if ( copy_spindles[i].include ) spindles->push_back( copy_spindles[i] );
       logger << "  Q-filtered spindle list from " << copy_spindles.size() << " to " << spindles->size() << "\n";
     }
   else
     logger << "  no spindles removed based on Q-filtering\n";



   //
   // Random pruning (debug) 
   //
   
   if ( prune ) 
     {
       const int n = spindles->size();
       std::vector<spindle_t> spindles0 = *spindles;
       spindles->clear();
       for (int i=0; i<n; i++)
	 {
	   if ( CRandom::rand() <= prune_rate )
	     spindles->push_back( spindles0[i] );
	 }
       
       logger << "  randomly selected " << spindles->size() << " spindles\n";
       
     }
   
   
   if ( prune_fixed )	
     {
       // select a fixed N - if we can't then set to 0
       //  i.e. so we have either exactly N spindles, or none

       const int n = spindles->size();
       
       if ( n < prune_n )
	 {
	   spindles->clear();
	   logger << "  did not achieve " << prune_n << " spindles - removing all\n";
	 }
       else if ( n > prune_n )
	 {
	   
	   std::vector<spindle_t> spindles0 = *spindles;
	   spindles->clear();
	   
	   // take N at random
	   std::vector<int> a( n );
	   CRandom::random_draw( a );
	   for (int i=0; i<prune_n; i++)
	     spindles->push_back( spindles0[ a[i] ] );
	   
	   logger << "  randomly selected fixed value of " << prune_n << " spindles\n";
	   
	 }
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


void per_spindle_output( std::vector<spindle_t> * spindles ,
			 param_t & param , 
			 clocktime_t * starttime ,
			 spindle_qc_t & q )
{
 
  const int n = spindles->size();

  const bool has_coupling = param.has( "so" );

  const bool has_if = param.has( "if" );
  
  //
  // Per-spindle output
  //
  
  for (int i = 0 ; i < spindles->size(); i++ ) 
    {
      
      spindle_t * spindle = &(*spindles)[i];
      
      writer.level( i+1 , "SPINDLE" );  // 1-based spindle count
 
      
      writer.value( "START"   , spindle->tp.start * globals::tp_duration );
      double peak_amp_sec = spindle->peak_amp_sp / (double)( spindle->stop_sp - spindle->start_sp ) ;
      peak_amp_sec = spindle->tp.start + peak_amp_sec * ( spindle->tp.stop - spindle->tp.start ) ;  
      writer.value( "PEAK"    ,  peak_amp_sec * globals::tp_duration );
      writer.value( "TROUGH"  ,  spindle->tp_mid * globals::tp_duration );
      writer.value( "STOP"    , spindle->tp.stop * globals::tp_duration );
      
      writer.value( "START_SP"   , spindle->start_sp );
      writer.value( "PEAK_SP"    ,  spindle->start_sp + spindle->peak_amp_sp );
      writer.value( "TROUGH_SP"  ,  spindle->start_sp + spindle->max_trough_sp );
      writer.value( "STOP_SP"    , spindle->stop_sp  );
      
      if ( spindle->winsor != -1 )
	writer.value( "WINSOR"    , spindle->winsor / (double)( spindle->stop_sp - spindle->start_sp + 1 )  );
       
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
       writer.value( "ACT_MX" , spindle->norm_amp_max      );
       writer.value( "ACT_MN" , spindle->norm_amp_mean      );
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
	   
	   // based on slopes
	   writer.value( "FS_POS"    , spindle->pos2f      );
	   writer.value( "FS_NEG"    , spindle->neg2f      );
	   writer.value( "FS_DIF"    , spindle->pos2f - spindle->neg2f );

	   writer.value( "BS_POS"    , spindle->pos2b      );
	   writer.value( "BS_NEG"    , spindle->neg2b      );	   
	   writer.value( "BS_DIF"    , spindle->pos2b - spindle->neg2b );

	   writer.value( "AMP_POS" , spindle->posisa      );
	   writer.value( "AMP_NEG" , spindle->negisa      );
	   
	   writer.value( "SP_POS" , spindle->possp  );
	   writer.value( "SP_NEG" , spindle->negsp  );
	   
	 }
       
       
       if ( spindle->chirp > -99998 ) 
	 {	   
	   writer.value( "CHIRP"  , spindle->chirp );
	   writer.value( "FRQ1"  , spindle->frq_h1 );
	   writer.value( "FRQ2"  , spindle->frq_h2 );	   
	 }
       
       if ( has_coupling )
	 {
	   // if no valid anchor, make NA in output
	   if ( spindle->anchor_sec >= 0 )
	     writer.value( "ANCHOR" , spindle->anchor_sec );
	   
	   if ( spindle->so_nearest_num != 0 ) 
	     {
	       writer.value( "SO_NEAREST" , spindle->so_nearest );
	       writer.value( "SO_NEAREST_NUM" , spindle->so_nearest_num );
	     }
	   
	   if ( spindle->so_phase_anchor >= 0 ) 
	     writer.value( "SO_PHASE_ANCHOR" , spindle->so_phase_anchor );
	 }

       if ( has_if )
	 writer.value( "IF" , spindle->if_spindle );
       
       // Q-score
       if ( q.qc_q )
	 {	   
	   writer.value( "Q"      , spindle->qual     );
	   writer.value( "PASS"   , spindle->include  );
	 }
       
     }

   // end of per-spindle output
   writer.unlevel( "SPINDLE" );
    
}



void spindle_qc_t::proc( param_t & param )
{

  qc_q = true;

  if ( param.has( "noq" ) )
    qc_q = false;
  
  // thresholds

  if ( qc_q )
    {
      if ( param.has( "q" ) )
	qc_qmin = param.requires_dbl( "q" );
      else
	qc_qmin = 0;
      
      if ( qc_qmin <= -9 )
	qc_q = false;
    }

  // nothing to do?
  if ( ! qc_q )
    {
      logger << "  no Q-score filtering will be applied\n";
      return;      
    }
       
  if ( param.has( "q-max" ) )
    qc_qmax = param.requires_dbl( "q-max" );

  // outputs
  verbose_failed = param.has( "q-verbose" );
  verbose_all = param.has( "q-verbose-all" );
  
  // defaults : 6 Hz  (5 cycles)
  //          : 35 Hz ( 7 cycles) 

  ofc.clear();
  ofc.push_back( 6 );
  ofc.push_back( 35 );
  ocycles.clear();
  ocycles.push_back( 5 );
  ocycles.push_back( 7 );
  
  // target fcs
  if ( param.has( "q-frq" ) )
    {
      ofc.clear();
      ocycles.clear();
      std::vector<double> pp = param.dblvector( "q-frq" );
      if ( pp.size() % 2 != 0 )
	Helper::halt( "expecting q-frq=fc,cycles(,fc,cycles,...}" );
      for (int i=0; i<pp.size(); i++)
	{
	  if ( pp[i] < 0 ) Helper::halt( "expecting q-frq F_cvalues above 0 Hz" );
	  if ( pp[i+1] < 3 ) Helper::halt( "expecting q-frq cycles values 3+ " );
	  ofc.push_back( pp[i] );
	  ocycles.push_back( pp[i+1] );
	}
    }

  // store
  nob = ofc.size();

  logger << "  Q-filtering q<" << qc_qmin ;
  for (int i=0; i<nob; i++)
    logger << ", " << ofc[i] << "Hz (" << ocycles[i] << " cyc)";
  logger << "\n";
  
}

void spindle_qc_t::init( const std::vector<double> * _ts , int _Fs )
{

  if ( ! qc_q ) return;       
  
  // keep a copy of the original time-series
  // (this is *assumed* to still be in scope whenever q.output() is called) 

  ts = _ts;
  Fs = _Fs;
  
  ocwt.resize( nob );
  
  CWT cwt;
  cwt.set_sampling_rate( Fs );
  for (int fi=0; fi<nob; fi++)
    cwt.add_wavelet( ofc[fi] , ocycles[fi] );  // f( Fc , number of cycles )
  cwt.load( ts );
  cwt.run();
  
  for (int fi=0; fi<nob; fi++)
    {
      ocwt[fi] = cwt.results(fi);
      const double om = MiscMath::mean( ocwt[fi] );
      for (int i=0; i<ocwt[fi].size(); i++)
	ocwt[fi][i] /= om;
    }
}


void spindle_qc_t::init_spindle( const std::vector<double> & s )
{
  if ( ! qc_q ) return;

  scwt = s;
  double sm = MiscMath::mean( scwt );
  for (int i=0; i<scwt.size(); i++)
    scwt[i] /= sm;
}


double spindle_qc_t::qual( spindle_t * spindle )
{
  
  double s = 0;
  std::vector<double> o( nob , 0 );
  const int np = spindle->stop_sp - spindle->start_sp + 1;
  
  for (int i=spindle->start_sp; i<= spindle->stop_sp; i++)
    {
      s += scwt[i];
      for (int j=0; j<nob; j++)
	o[j] += ocwt[j][i];
    }

  // means
  s /= (double)np;
  for (int j=0; j<nob; j++)
    o[j] /= (double)np;
  
  // get max
  double omx = o[0];
  for (int j=1; j<nob; j++)
    if ( o[j] > omx ) omx = o[j];
  if ( omx <= 0 ) omx = 1e-12;
  
  return log2( s ) - log2( omx );

}

void spindle_qc_t::output( spindle_t * spindle , const int n )
{
  // verbose outputs
  const int np = spindle->stop_sp - spindle->start_sp + 1;

  int s1 = spindle->start_sp - 2 * Fs; 
  int s2 = spindle->stop_sp + 2 * Fs;
  if ( s1 < 0 ) s1 = 0;
  if ( s2 >= ts->size() ) s2 = ts->size() - 1 ;
  
  for (int i=s1; i<=s2; i++)
    {
      const bool in_spindle = i >= spindle->start_sp && i <= spindle->stop_sp;
      
      std::cout << "QS\t"
		<< n << "\t"	
		<< in_spindle << "\t"
		<< (*ts)[i] << "\t"
		<< scwt[i] ;
      for (int j=0; j<nob; j++)
	std::cout << "\t" << ocwt[j][i];
      std::cout << "\n";
    }
}


void spindle_stats( const std::vector<spindle_t> & spindles , std::map<std::string,double> & results ) 
{

  double dur = 0 , fwhm = 0 , amp = 0 , nosc = 0 , isa = 0 , qual = 0 ;

  double norm_amp_max = 0 , norm_amp_mean = 0;

  double winsor1 = 0 , winsor2 = 0;
  
  // relative measures of spindle "peaks" (0..1)
  double symm = 0 , symm2 = 0, symm_amp = 0 , symm_max_trough = 0;
  
  // in seconds from start
  double symm_sec = 0 , symm_amp_sec = 0 , symm_max_trough_sec = 0;

  // freq/chirp 
  double frq = 0 , fft = 0 , chirp = 0, frq1 = 0 , frq2 = 0 , frq_range = 0;

  // HWs
  double negf = 0 , posf = 0 , allf = 0;
  double negb = 0 , posb = 0 , allb = 0;
  double negv = 0 , posv = 0 , allv = 0;
  double posisa = 0 , negisa = 0 , possp = 0 , negsp = 0;
  
  // slopes
  double neg2f = 0 , pos2f = 0 ;
  double neg2b = 0 , pos2b = 0 ;

  int denom = 0;

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

      // relative symmetry measures
      symm += ii->symm;   // based on max peak-to-peak
      symm2 += ii->symm2; // as above
      symm_amp += ii->peak_amp_rel; 
      symm_max_trough += ii->max_trough_rel;
      
      // in seconds
      symm_sec += ii->symm * ii->dur ;
      symm_amp_sec += ii->peak_amp_rel * ii->dur;
      symm_max_trough_sec += ii->max_trough_rel * ii->dur;

      // winsorization
      winsor1 += ii->winsor > 0 ; // any clipping for this spindle?
      winsor2 += ii->winsor / (double)(ii->stop_sp - ii->start_sp + 1 );
      
      // freq
      frq += ii->frq;
      fft += ii->fft;
      chirp += ii->chirp;      
      frq1 += ii->frq_h1;
      frq2 += ii->frq_h2;
      frq_range += ii->frq_range;

      negf += ii->negf;
      posf += ii->posf;
      allf += ii->allf;
      
      negb += ii->negb;
      posb += ii->posb;
      allb += ii->allb;
      
      negv += ii->negv;
      posv += ii->posv;
      allv += ii->allv;
      
      // based on slopes
      neg2f += ii->neg2f;
      pos2f += ii->pos2f;
            
      neg2b += ii->neg2b;
      pos2b += ii->pos2b;

      isa += ii->isa;
      posisa += ii->posisa;
      negisa += ii->negisa;

      norm_amp_max += ii->norm_amp_max;
      norm_amp_mean += ii->norm_amp_mean;
      
      possp += ii->possp;
      negsp += ii->negsp;
      
      qual += ii->qual;
      
      ++ii;
    }
  

  results[ "AMP" ]      = amp /(double)denom;

  results[ "ACT_MX" ]      = norm_amp_max /(double)denom;
  results[ "ACT_MN" ]      = norm_amp_mean /(double)denom;

  results[ "TOTDUR" ]   = dur;
  results[ "DUR" ]      = dur / (double)denom;
  results[ "FWHM" ]     = fwhm / (double)denom;
  results[ "NOSC" ]     = nosc / (double)denom;
  results[ "FRQ" ]      = frq / (double)denom;
  results[ "FFT" ]      = fft / (double)denom;

  results[ "SYMM" ]     = symm / (double)denom;
  results[ "SYMM2" ]    = symm2 / (double)denom;
  results[ "SYMM_AMP" ] = symm_amp / (double)denom;
  results[ "SYMM_TROUGH" ]    = symm_max_trough / (double)denom;

  results[ "SEC_P2P" ]    = symm_sec / (double)denom;
  results[ "SEC_AMP" ]    = symm_amp_sec / (double)denom;
  results[ "SEC_TROUGH" ] = symm_max_trough_sec / (double)denom;

  results[ "CHIRP" ]   = chirp / (double)denom;
  results[ "FRQ1" ]    = frq1 / (double)denom;
  results[ "FRQ2" ]    = frq2 / (double)denom;
  results[ "FVAR2" ]    = allv / (double)denom;
  results[ "FRNG2" ]    = frq_range / (double)denom;

  results[ "FNEG" ]    = negf / (double)denom;
  results[ "FPOS" ]    = posf / (double)denom;
  results[ "FALL" ]    = allf / (double)denom;

  results[ "BNEG" ]    = negb / (double)denom;
  results[ "BPOS" ]    = posb / (double)denom;
  results[ "BALL" ]    = allb / (double)denom;

  results[ "VNEG" ]    = negv / (double)denom;
  results[ "VPOS" ]    = posv / (double)denom;
  results[ "VALL" ]    = allv / (double)denom;

  // slopes
  results[ "FSNEG" ]    = neg2f / (double)denom;
  results[ "FSPOS" ]    = pos2f / (double)denom;

  results[ "BSNEG" ]    = neg2b / (double)denom;
  results[ "BSPOS" ]    = pos2b / (double)denom;

  results[ "Q" ]        = qual / (double)denom;
  results[ "W_ANY" ]    = winsor1 / (double)denom;
  results[ "W" ]        = winsor2 / (double)denom;
  
  results[ "ISA_PER_SPINDLE" ] = isa / (double)denom;
  results[ "ISA_TOTAL" ] = isa;
  
  results[ "POSISA_PER_SP" ] = posisa / (double)denom;
  results[ "POSSP" ] = possp / (double)denom;

  results[ "NEGISA_PER_SP" ] = negisa / (double)denom;
  results[ "NEGSP" ] = negsp / (double)denom;
    

}

