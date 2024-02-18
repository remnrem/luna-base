
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

#include "mtm.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "eval.h"
#include "fftw/fftwrap.h"
#include "fftw/bandaid.h"

#include "db/db.h"
#include "helper/helper.h"
#include "helper/logger.h"

extern writer_t writer;
extern logger_t logger; 

#include "stats/Eigen/Dense"
#include "stats/eigen_ops.h"





void mtm::wrapper( edf_t & edf , param_t & param )
{
  
  std::string signal_label = param.requires( "sig" );

  signal_list_t signals = edf.header.signal_list( signal_label );    

  std::vector<double> Fs = edf.header.sampling_freq( signals );
  
  const int ns = signals.size();

  // nb. for efficiency, MTM uses its own segmentation of signals as well as generic
  // epochs;  epochs are optional - either way, segments are specified.  In a typical
  // case by segment and epoch might be 30 seconds, i.e. one segment per epoch
  //  can specify either segment-level output and/or epoch-level output;  in the
  //  above case, they would be the same
  
  // it also gives an average spectrum at the end, which is the average/median of all epoch-level
  // results, if the data are epoched
  
  // nb. using whole trace assumes a continuous EDF; using epochs
  //     allows for gaps (i.e. as an epoch is always continuos)

  //
  // Segment size (within epoch) definitions
  //
  
  const double segment_size_sec = param.has( "segment-sec" ) ? param.requires_dbl( "segment-sec" ) : 30;
  const double segment_step_sec = param.has( "segment-inc" ) ? param.requires_dbl( "segment-inc" ) : segment_size_sec ;
  
  
  //
  // set up band values
  //

  bandaid_t bandaid;

  bandaid.define_bands( param );

  // report bands?
  
  const bool bands = param.has( "band" ) ? param.yesno( "band" ) : true ;

  
  //
  // spectral kurtosis values
  //

  // 'excess' kurtosis, i.e. = 0 for N(0,1)
  const bool spec_kurt  = param.has( "speckurt" ) || param.has( "speckurt3" ) || param.has( "alternate-speckurt" );

  // log-scale spec-kurt (for fbin only right nowe)
  const bool sklog = param.has( "speckurt-fbin-log" );
  
  // 'alternative' definition of spectral kurtosis (i.e. kurt of both time and freq bins, within band)
  const bool kurt_altdef = param.has( "alternate-speckurt" );

  // standardard k, i.e. = 3 for N(0,1)
  const bool spec_kurt3 = param.has( "speckurt3" );

  //
  // ratios of band-powers
  //

  // output ratios of band power values
  const bool calc_ratio = param.has( "ratio" );
  if ( calc_ratio && param.empty( "ratio" ) )
    Helper::halt( "cannot have empty ratio arg" );
  // ratio=ALPHA/BETA,THETA/DELTA,...
  const std::string ratios = calc_ratio ? param.value( "ratio" ) : "" ;
  const int ratio_plus1 = param.has( "ratio1" );


  // running w/ out epoch
  // whole recording -----------------------------|   overall stats
  // |seg1|seg2|seg3|seg4|seg5|seg6|seg7|seg8|seg9|   SEG level stats

  // whole recording -----------------------------|   overall stats 
  // | epoch1 ------| epoch2-------| epoch3 ------|   epoch-level stats
  // |seg1|seg2|seg3|seg1|seg2|seg3|seg1|seg2|seg3|   SEG level stats

  //
  // epoch-wise analysis? (as per standard epoch mechanism)
  //

  const bool epoch_level_output = param.has( "epoch-output" );
  
  const bool epoch_level_output_spectra = param.has( "epoch-spectra" );

  // we might want to compute the whole-night stats epoch-wise, but not care about epoch
  // level output, so allow 'epoch' ( by itself, does not give epoch level output) 
  
  const bool epochwise = param.has( "epoch" ) || epoch_level_output ; 

  //
  // if epochs specified, check that segments will fit (at least one per epoch)
  //

  if ( epochwise )
    {

      edf.timeline.ensure_epoched();
      edf.timeline.first_epoch();

      // if using generic epochs, then we need to check below on an epoch-by-epoch basis
      // otherwise, we can check here

      if ( ! edf.timeline.generic_epochs() )
	{
	  if ( edf.timeline.epoch_length() < segment_size_sec )
	    {
	      logger << "  segment size = " << segment_size_sec << " seconds\n"
		     << "  epoch size   = " << edf.timeline.epoch_length() << " seconds\n";
	      Helper::halt( "segment-sec is larger than epoch duration...\nno segments would be placed, please change parameters" );
	    } 
	}
    }
  

  //  
  // report segment-level? i.e. epochs must contain 1 or more segments
  //
  
  const bool segment_level_output = param.has( "segment-output" ) || param.has( "segment-spectra" );
  
  const bool segment_level_output_spectra = param.has( "segment-spectra" );


  //
  // other misc functions
  //
  
  const bool display_tapers = param.has( "dump-tapers" );
  
  const bool mean_center = param.has( "mean-center" ) || param.has( "mean-centre" );

  const bool remove_linear_trend = param.has( "detrend" );
  
  if ( mean_center && remove_linear_trend )
    Helper::halt( "cannot specify both mean-center and detrend" );
  
  //
  // create new signals?
  // prefix_CH_N ... where N = 1,2,3, that correspond to Fs in range
  //
  
  const bool new_sigs = param.has( "add" ) ;
  
  std::string new_sig_prefix = new_sigs ? param.value( "add" ) : "" ; 

  if ( new_sigs && epochwise )
    Helper::halt( "cannot specify 'add' (to make a new signal) when using epoch-level analysis" );
  

  //
  // MTM parameters (tw or nw)
  //
  
  double npi = 3;
  if ( param.has( "nw" ) ) npi = param.requires_dbl( "nw" );
  else if ( param.has( "tw" ) ) npi = param.requires_dbl( "tw" );
  
  const int nwin = param.has( "t" ) ? param.requires_int( "t" ) : 2*floor(npi)-1 ;


  //
  // Required minimum SR to attempt MTM
  //

  const int min_sr = param.has( "sr" ) ? param.requires_int( "sr" ) : 0 ; 

  //
  // Start/stop times? Only allow this in wholetrace (i.e. not epoch) mode
  //
  
  const bool   restrict_start = param.has( "start" );
  const bool   restrict_stop = param.has( "stop" );
  const double restrict_start_sec = param.has( "start" ) ? param.requires_dbl( "start" ) : 0 ;
  const double restrict_stop_sec = param.has( "stop" ) ? param.requires_dbl( "stop" ) : 0 ;

  if ( ( restrict_start || restrict_stop ) && epochwise )
    Helper::halt( "can only specify start/stop times (to select a subset of segments) when not in epoch-mode" );
  
  //
  // Reporting full spectrum? (default 0.5 to 25 Hz)
  //
  
  double min_f = param.has( "min" ) ? param.requires_dbl( "min" ) : 0.5; 
  double max_f = param.has( "max" ) ? param.requires_dbl( "max" ) : 25;  


  //
  // Spectral slope
  //
  
  const bool spectral_slope = param.has( "slope" );
  const std::vector<double> slope_range = param.dblvector( "slope" );

  if ( spectral_slope && epochwise )
    Helper::halt( "cannot currently do slope and epochwise analysis for MTM\nuse segment-slopes" );
  
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

    
  // output
  
  const bool dB = param.has( "dB" );


  //
  // Channel checks
  //
  
  std::set<int> srs;

  int ns_used = 0;
  
  for (int s = 0 ; s < ns; s++ )
    {

      if ( edf.header.is_annotation_channel( signals(s) ) )
	continue;
      
      if ( min_sr && Fs[s] < min_sr )
	continue;

      ++ns_used;
      
      srs.insert( Fs[s] );
    }

  if ( ns_used == 0 ) return;
  
  if ( spec_kurt && srs.size() != 1 )
    Helper::halt( "all SRs must be similar if using speckurt option" );
  
  //
  // Precompute tapers (for each Fs) 
  //

  logger << "  precomputing " << nwin << " tapers for " << srs.size() << " distinct sample rates\n";
  std::map<int,mtm_t> sr2tapers;
  std::set<int>::const_iterator ss = srs.begin();
  while ( ss != srs.end() )
    {      
      mtm_t mtm( npi , nwin );
      const int segment_size = *ss * segment_size_sec;
      mtm.store_tapers( segment_size );      
      sr2tapers[ *ss ] = mtm;
      ++ss;
    }

  //
  // Epoch-trackers
  //  - in epoch mode, will need to track epoch-level stats to get final versions
  //
  
  // channel->band->epoch->value
  std::vector<std::vector<std::vector<double> > > etrack_bandpower;
  if ( bands ) 
    etrack_bandpower.resize( ns_used );

  // epoch length (for generic case)
  std::vector<double> etrack_length;
  
  std::vector<std::vector<std::vector<double> > > etrack_relbandpower;
  etrack_relbandpower.resize( ns_used );

  // channel->ratio->epoch->value
  std::vector<std::vector<std::vector<double> > > etrack_bandratios;
  if ( bands )
    etrack_bandratios.resize( ns_used );

  // channel->freq (if diff SRs?)  
  std::vector<std::vector<double> > etrack_freqs;
  etrack_freqs.resize( ns_used );
  
  // channel->freq->epoch->value  
  std::vector<std::vector<std::vector<double> > > etrack_power;
  etrack_power.resize( ns_used );
  
  // channel->epoch->value
  // std::vector<std::vector<double> > etrack_slope;
  // etrack_slope.resize( ns_used );

  // band->epoch->value (already averaged over channels)
  //   avg over channels
  std::vector<std::vector<double> > etrack_chavg_speckurt;  
  std::vector<std::vector<double> > etrack_chavg_specskew;
  std::vector<std::vector<double> > etrack_chavg_speccv;
  // channel->band->epoch->value
  std::vector<std::vector<std::vector<double> > > etrack_speckurt;  
  std::vector<std::vector<std::vector<double> > > etrack_specskew;
  std::vector<std::vector<std::vector<double> > > etrack_speccv;
  if ( bands )
    {
      etrack_speckurt.resize( ns_used );
      etrack_specskew.resize( ns_used );
      etrack_speccv.resize( ns_used );
    }
  
  // channel->freqbin->epoch->value
  std::vector<std::vector<std::vector<double> > > etrack_fspeckurt;
  std::vector<std::vector<std::vector<double> > > etrack_fspecskew;
  std::vector<std::vector<std::vector<double> > > etrack_fspeccv;
  if ( bands )
    {
      etrack_fspeckurt.resize( ns_used );
      etrack_fspecskew.resize( ns_used );
      etrack_fspeccv.resize( ns_used );
    }
  
  if ( epochwise )
    logger << "  epochwise analysis, iterating over " << edf.timeline.num_epochs() << " epochs\n";
  else
    logger << "  single analysis of the entire available signal\n";
  
  //
  // Either iterate over epochs, or whole trace
  //

  bool show_initial_log_output = true;
  
  while ( 1 ) 
    {
      
      //
      // which interval to analyse?
      //  -- either one epoch or the whole interval 
      //
      
      int epoch = epochwise ? edf.timeline.next_epoch() : 0 ;
      
      if ( epoch == -1 ) break;
            
      interval_t interval = epochwise ? edf.timeline.epoch( epoch ) : edf.timeline.wholetrace(); 
      
      //
      // Need to check segment length?
      //

      if ( epochwise && edf.timeline.generic_epochs() )
	{
	  // here, all epochs need to have the same segment length, so will skip here
	  // if the epoch is too short
	  if ( edf.timeline.epoch_length() < segment_size_sec )
	    {
	      logger << "  *** skipping epoch " << interval.as_string()
		     << ", too short given segment-sec = " << segment_size_sec << "\n";
	      continue;
	    }
	}
      

      //
      // Track epoch length (function works for standard & generic epochs)
      //
      
      etrack_length.push_back( edf.timeline.epoch_length() );
      
      //
      // Stratify output by epoch, if epoch-level analysis/output
      //
      
      if ( epoch_level_output )
	writer.epoch( edf.timeline.display_epoch( epoch ) );

      //
      // Used if tracking spectral kurtosis (average over channels)
      //
      
      spectral_kurtosis_t skurt( spec_kurt3 );
      
      
      //
      // For this interval, now process all signals
      //

      int ns1 = -1; // actually used channel index
      
      for (int s = 0 ; s < ns; s++ )
	{
	  
	  //
	  // only consider data tracks
	  //
	  
	  if ( edf.header.is_annotation_channel( signals(s) ) )
	    continue;
	  
	  //
	  // Min. required SR?
	  //
	  
	  if ( min_sr && Fs[s] < min_sr )
	    continue;

	  //
	  // Get used-channel index for use below
	  //

	  ++ns1;
	  
	  //
	  // Clear up
	  //
	
	  bandaid.init();

	  
	  //
	  // Stratify output by channel
	  //
	  
	  writer.level( signals.label(s) , globals::signal_strat );
      
	  
	  //
	  // Get data
	  //
      
	  slice_t slice( edf , signals(s) , interval );
	  
	  const std::vector<double> * d = slice.pdata();
	  
	  //
	  // Step size in sample-points
	  //
	  
	  const int segment_size = Fs[s] * segment_size_sec;
	  const int segment_step = Fs[s] * segment_step_sec;
	  const uint64_t delta_tp = globals::tp_1sec / Fs[s] ;
	  
	  //
	  // Get time points (and flags for segments that span discontinuities)
	  //  - also, indicate whether all should be computed (at segment/epoch level)
	  //  - this is for moonlight MTM interactive viewer mainly
	  
	  const std::vector<uint64_t> * tp = slice.ptimepoints();
	  const int np = tp->size();
      
	  std::vector<double> start, stop;
	  std::vector<int> start_sp, stop_sp; // original signal encoding
      
	  // when 'add' option, count number of segments spanning this sample point
	  std::vector<int> addn( np , 0 ); 
	  
	  // get 
	  int nf = -1;
	  Eigen::MatrixXd addX;
	  
	  std::vector<bool> disc, restrict;
	  
	  int p = 0;
	  int nn = 0;
	  int actual = 0;
	  
	  // actual segment size (may differ from requested due to sample rates)
	  // here, check against this);    actually, seems like it will be okay
	  // to use start/stop gap checking w/ epoch-level analysis.  but for now
	  // can leave as is
	  
	  const double segment_sec2 = segment_size / (double)Fs[s];
	  
	  while ( 1 ) {
	    if ( p + segment_size > np ) break;
	    
	    double start_sec = (*tp)[p] * globals::tp_duration;
	    double stop_sec = ( (*tp)[ p + segment_size - 1 ] + delta_tp ) * globals::tp_duration; // '1past'
	    double implied_sec = stop_sec - start_sec;
	
	    start_sp.push_back( p );
	    stop_sp.push_back( p + segment_size - 1 ) ; // 'last point in seg'
	    start.push_back( start_sec );
	    stop.push_back( stop_sec );
	    disc.push_back( fabs( implied_sec - segment_sec2 ) > 0.0001 );
	    
	    bool okay = true;
	    if ( restrict_start && start_sec < restrict_start_sec ) okay = false;
	    if ( restrict_stop && stop_sec > restrict_stop_sec ) okay = false;
	    restrict.push_back( ! okay );
	    if ( okay ) ++actual;
	    ++nn;
	    
	    // std::cout << "seg " << nn << "\t" << p << "\t" << start_sec << "\t" << stop_sec << "\t"
	    // 	  << " sz " << implied_sec << " " << segment_sec2 << " " 
	    // 	  << ( fabs( implied_sec - segment_sec2 ) > 0.001 )
	    //   	  << "\trestrict=" << ! okay << "\n";
	
	    // next segment
	    p += segment_step;
	    
	  }
	  
	  // nothing to do... note; this should not happen in epoch mode, as we've already checked the size
	  // above???  either way, fine to leave as is for now, i.e. is a bad set of parameters if any
	  // epoch is shorter than the segment size
	  
	  if ( actual == 0 )
	    {
	      logger << "  *** no segments to process, leaving MTM...\n";
	      return;
	    }
	  
	  //
	  // call MTM
	  //
      
	  mtm_t mtm( npi , nwin );
	  
	  mtm.dB = dB;
	  mtm.opt_remove_mean = mean_center;
	  mtm.opt_remove_trend = remove_linear_trend;
	  mtm.bandaid = &bandaid;
		  
	  // possibly restrict to a subset of segments?
	  if ( restrict_start || restrict_stop )
	    mtm.restrict = restrict;
	  
	  // actual MTM 
	  mtm_t * precomputed = &(sr2tapers[ Fs[s] ]);
	  mtm.apply( d , Fs[s] , segment_size , segment_step , show_initial_log_output , precomputed );
	  
	  // cannot show channel progress in epoch mode (epoch then channel processed)
	  if ( ! epochwise )
	    {
	      if ( show_initial_log_output ) logger << "  processed channel(s):";
	      logger << " " << signals.label(s) ;
	    }
	  
	  //
	  // count freq bins on first pass if not already done? could add only reduced set perhaps?
	  //   - this is not allowed in epoch-mode
	  
	  if ( new_sigs && nf == -1 ) 
	    {	  
	      nf = 0;
	      for ( int i = 0 ; i < mtm.f.size() ; i++ )
		if ( mtm.f[i] >= min_f && mtm.f[i] <= max_f )
		  ++nf;
	      addX = Eigen::MatrixXd::Zero( np , nf );	  
	    }

	  //
	  // track spec-kurt?
	  //

	  if ( spec_kurt )
	    {
	      // first channel only, add freqs)
	      if ( ns1 == 0 ) skurt.setf( mtm.f );

	      // add data
	      skurt.add( ns1 , mtm.raw_espec );
	    }
	      
	  //
	  // Output: tapers? (only once) 
	  //
	  
	  if ( display_tapers && show_initial_log_output )
	    {	  
	      for( int i=0; i< mtm.tapers.rows(); i++)
		{
		  writer.level( i+1 , "SP" );
		  for(int j=0; j< mtm.tapers.cols(); j++)
		    {
		      writer.level( j+1 , "TAPER" );
		      writer.value( "W" , mtm.tapers(i,j) );
		    }
		  writer.unlevel( "TAPER" );
		}
	      writer.unlevel( "SP" );
	      
	      for(int j=0; j< mtm.lam.size(); j++) 
		{
		  writer.level( j+1 , "TAPER" );
		  writer.value( "LAMBDA" , mtm.lam[j] );
		}
	      writer.unlevel( "TAPER" );
	      
	    }
	  
	  //
	  // track freqs for epoch-level analysis
	  //

	  if ( epochwise && etrack_freqs[ns1].size() == 0 ) 
	    {
	      for ( int i = 0 ; i < mtm.f.size() ; i++ )
		if ( mtm.f[i] >= min_f && mtm.f[i] <= max_f )
		  etrack_freqs[ns1].push_back( mtm.f[i] );
	      // size power values slots
 	      etrack_power[ns1].resize( etrack_freqs[ns1].size() );
	      
	    }
	  
	  //
	  // Track epoch spectra?
	  //
	  
	  if ( epochwise )
	    {
	      int fidx = 0;
	      for ( int i = 0 ; i < mtm.f.size() ; i++ )
                if ( mtm.f[i] >= min_f && mtm.f[i] <= max_f )
                  {
		    // channel -> freq -> epoch/value
		    etrack_power[ ns1 ][ fidx ].push_back( mtm.spec[i] );
		    ++fidx;
		  }
	    }


	  
	  //
	  // Output: averaged spectrum (will be either per epoch, or overall)
	  //
	  
	  if ( ( ! epochwise ) || ( epochwise && epoch_level_output_spectra ) )
	    {
	      for ( int i = 0 ; i < mtm.f.size() ; i++ ) 
		{
		  if ( mtm.f[i] >= min_f && mtm.f[i] <= max_f ) 
		    {
		      writer.level( mtm.f[i] , globals::freq_strat  );
		      writer.value( "MTM" , mtm.spec[i] );
		    }
		}
	      writer.unlevel( globals::freq_strat );
	    }
	  
	  
	  //
	  // bands - output?
	  //

	  if ( bands )
	    {
	      if ( ( ! epochwise ) || ( epochwise && epoch_level_output ) )
		{
		  bandaid.calc_bandpower( mtm.f , mtm.raw_spec );
		  
		  const double mean_total_power = bandaid.fetch( DENOM );
		  
		  std::vector<frequency_band_t>::const_iterator bb = bandaid.bands.begin();
		  while ( bb != bandaid.bands.end() )
		    {
		      writer.level( globals::band( *bb ) , globals::band_strat );
		      double p = bandaid.fetch( *bb );
		      writer.value( "MTM" , dB ? 10*log10(p) : p  );
		      writer.value( "REL" , p / mean_total_power );
		      
		      ++bb;
		    }
		  writer.unlevel( globals::band_strat );
		}
	      
	      //
	      // band - track?
	      //
	      
	      if ( epochwise )
		{
		  // make space as needed?
		  if ( etrack_bandpower[ns1].size() == 0 )
		    {
		      etrack_bandpower[ns1].resize( bandaid.size() );
		      etrack_relbandpower[ns1].resize( bandaid.size() );
		    }
		  
		  // add this epoch
		  bandaid.calc_bandpower( mtm.f , mtm.raw_spec );
		  
		  const double mean_total_power = bandaid.fetch( DENOM );
		  
		  int bi=0;
		  std::vector<frequency_band_t>::const_iterator bb = bandaid.bands.begin();
		  while ( bb != bandaid.bands.end() )
		    {
		      double p = bandaid.fetch( *bb ) ;
		      // dB? scaling (all bandaid_t storage is raw)
		      etrack_bandpower[ns1][bi].push_back( dB ? 10*log10(p) : p  );
		      etrack_relbandpower[ns1][bi].push_back( p / mean_total_power );
		      ++bi;
		      ++bb;
		    }
		}
	      
	      //
	      // Band-power ratios ( gets the ratio per segment, then takes mean/median for this interval (epoch/whole)
	      //
	      
	      if ( calc_ratio )
		{
		  // ratio=ALPHA/BETA,THETA/DELTA,...
		  std::vector<std::string> r = Helper::parse( Helper::toupper( ratios ) , ',' );
		  
		  // are we tracking? if so, make space
		  if ( epochwise && etrack_bandratios[ns1].size() == 0 )
		    etrack_bandratios[ns1].resize( r.size() );
		  
		  bool done_any = false; 
		  
		  int rn = 0;
		  
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
		      
		      if ( b1 == UNKNOWN_BAND || b2 == UNKNOWN_BAND )
			Helper::halt( "unknown band values in ratios" );
		      
		      // i.e. all segments for this epoch/interval
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
			  
			  // track?
			  if ( epochwise )
			    {
			      etrack_bandratios[ns1][rn].push_back( rmean );
			    }
			  
			  // output now?
			  if ( ( ! epochwise ) || epoch_level_output )
			    {
			      const double rmedian = MiscMath::median( rat );
			      
			      writer.level( tok[0] , "B1" );
			      writer.level( tok[1] , "B2" );
			      writer.value( "RATIO" , rmean );
			      writer.value( "RATIO_MN" , pw1 / ( ratio_plus1 + pw2 ) );
			      writer.value( "RATIO_MD" , rmedian );		  
			      done_any = true;
			    }
			}
		      
		      ++rn;
		      ++rr;
		    }
		  
		  if ( done_any )
		    {
		      writer.unlevel( "B2" );
		      writer.unlevel( "B1" );
		    }
		  
		}
	    }

	  
	  //
	  // compute spectral slope?
	  //
	
	  // if dB mode, then based on the mean of db-scaled power (where
	  // we need to convert back to raw); otherwise, just use raw
	  
	  if ( spectral_slope ) 
	    {		   
	      
	      std::vector<double> xx = dB ? mtm.spec : mtm.raw_spec;
	      
	      if ( dB ) 
		for (int i=0; i<xx.size(); i++) 
		  xx[i] = pow( 10 , xx[i]/10.0) ;	      
	      
	      spectral_slope_helper( xx , 
				     mtm.f , 
				     slope_range , 
				     slope_outlier );
	    }
	  
      
      
	  //
	  // store (and optionally ouput) epoch-wise slope?
	  //
	  
	  std::vector<double> slopes;
	  
	  if ( segment_level_output || spectral_slope || new_sigs )
	    {
	      const int nsegs = mtm.espec.size();
	      
	      if ( nsegs != start.size() )
		Helper::halt( "internal error in MTM timing:"
			      + Helper::int2str( nsegs )
			      + " vs "
			      + Helper::int2str( (int)start.size() ) );
	      
	      if ( segment_level_output || new_sigs )
		{
		  for ( int j = 0 ; j < nsegs ; j++)
		    {
		      
		      if ( ! restrict[j] ) 
			{
			  
			  //
			  // add main output
			  //
			  
			  if ( segment_level_output )
			    {
			      writer.level( j+1 , "SEG" );	  
			      writer.value( "START" , start[j] );
			      writer.value( "STOP" , stop[j] );
			      writer.value( "DISC" , (int)disc[j] );
			      
			      // spectrogram-level output?
			      if ( segment_level_output_spectra )
				{
				  for ( int i = 0 ; i < mtm.f.size() ; i++ ) 
				    {
				      if ( mtm.f[i] >= min_f && mtm.f[i] <= max_f ) 
					{
					  writer.level( mtm.f[i] , globals::freq_strat  );
					  writer.value( "MTM" , mtm.espec[j][i] );
					}
				    }
				  writer.unlevel( globals::freq_strat );
				}
			      
			      //
			      // segment-level band-level output?
			      //

			      if ( bands )
				{
				  bandaid.calc_bandpower( mtm.f , mtm.raw_espec[j] );
				  
				  const double mean_total_power = bandaid.fetch( DENOM );
				  
				  std::vector<frequency_band_t>::const_iterator bb = bandaid.bands.begin();
				  while ( bb != bandaid.bands.end() )
				    {
				      writer.level( globals::band( *bb ) , globals::band_strat );
				      double p = bandaid.fetch( *bb );
				      writer.value( "MTM" , dB ? 10*log10(p) : p  );
				      if ( ! ( *bb == TOTAL || *bb == DENOM ) )
					writer.value( "REL" , p / mean_total_power );			  
				      ++bb;
				    }
				  writer.unlevel( globals::band_strat );
				}
			    }

			  
			  //
			  // make new signals?		      
			  //
			  
			  if ( new_sigs ) 
			    {
			      int s1 = start_sp[j];
			      int s2 = stop_sp[j];
			      //std::cout << " s1 s2 = " << s1 << " .. " << s2 << "\n";
			      for (int p=s1; p<=s2; p++)
				{
				  addn[p]++;
				  
				  int fidx = 0;
				  for ( int i = 0 ; i < mtm.f.size() ; i++ )
				    if ( mtm.f[i] >= min_f && mtm.f[i] <= max_f )
				      {
					addX(p,fidx++) += mtm.espec[j][i];
					//std::cout << " adding seg << " << j << " sample " << p << " freq " << fidx-1 <<"=" << mtm.f[i] << " = " << mtm.espec[j][i] << "\n";
				      }
				}
			      
			    }
			}
		      
		    }
		}
	      
	      
	      //
	      // segment level spectral slope? (based on raw power)
	      //
	  
	      if ( spectral_slope )
		{		  
		  
		  for ( int j = 0 ; j < nsegs ; j++)
		    {

		      if ( segment_level_output )
			writer.level( j+1 , "SEG" );
		      
		      if ( ! restrict[j] )
			{
			  double es1 = 0 ;
			  
			  bool okay = spectral_slope_helper( mtm.raw_espec[j] , 
							     mtm.f ,						     
							     slope_range ,
							     slope_outlier ,
							     segment_level_output , 
							     &es1 );
			  
			  if ( okay ) slopes.push_back( es1 );
			}
		    }
		}
	      
	    }
	  
	  if ( segment_level_output ) 
	    writer.unlevel( "SEG" );
      
   
	  
	  //
	  // spectral slope based on distribution of epoch-level slopes?
	  //
	  
	  if ( spectral_slope ) 
	    {
	      if ( slopes.size() > 2 )
		{
		  std::vector<double> s2 = MiscMath::outliers( &slopes , slope_th2 );
		  double s_mean = MiscMath::mean( s2 );
		  double s_med  = MiscMath::median( s2 );
		  double s_sd   = MiscMath::sdev( s2 , s_mean );
		  writer.value( "SPEC_SLOPE_MN" , s_mean );
		  writer.value( "SPEC_SLOPE_MD" , s_med );
		  writer.value( "SPEC_SLOPE_SD" , s_sd );
		}
	    }
	  
	  
	  //
	  // (per-channel) Spectral kurtosis (calc'ed per epoch/interval but also here averaging channels
	  //   -- only for 'primary' definiton (not alternate-speckurt) 
	  //
	  
	  if ( spec_kurt && ! kurt_altdef )
	    {

	      //
	      // output now?
	      //

	      if ( ( ! epochwise ) || epoch_level_output )
		{

		  if ( bands )
		    {
		      std::set<frequency_band_t>::const_iterator bb = skurt.bands.begin();
		      while ( bb != skurt.bands.end() )
			{
			  writer.level( globals::band( *bb ) , globals::band_strat );
			  
			  // current (i.e. latest added) channel = ns1
			  // kurtosis2() is the 'primary' def
			  double spsk , spcv ; 
			  double spku = skurt.kurtosis2( ns1 , *bb , &spcv, &spsk ) ;
			  
			  if ( spku > -900 ) { 
			    writer.value( "SPECCV" , spcv );
			    writer.value( "SPECSKEW" , spsk );
			    writer.value( "SPECKURT" , spku );
			  }
			  
			  ++bb;
			}
		      writer.unlevel( globals::band_strat );
		    }
		}
	      
	      //
	      // track?
	      //
	      
              if ( epochwise )
		{
		  if ( bands )
		    {
		      // bandwise speckurt
		      if ( etrack_speckurt[ns1].size() == 0 ) 
			{
			  etrack_speckurt[ns1].resize( skurt.bands.size() );
			  etrack_specskew[ns1].resize( skurt.bands.size() );
			  etrack_speccv[ns1].resize( skurt.bands.size() );
			}
		      
		      int bn=0;
		      std::set<frequency_band_t>::const_iterator bb = skurt.bands.begin();
		      while ( bb != skurt.bands.end() )
			{
			  
			  double spsk ,spcv ;
			  double spku = kurt_altdef ? skurt.kurtosis( *bb , &spcv, &spsk ) : skurt.kurtosis2( ns1, *bb , &spcv, &spsk ) ;
			  
			  etrack_speckurt[ns1][bn].push_back( spku );
			  etrack_speccv[ns1][bn].push_back( spcv );
			  etrack_specskew[ns1][bn].push_back( spsk );
			  
			  ++bb;
			  ++bn;
			}
		    }
		}
	      

	      //
	      // freqbin speckurt
	      //

	      int nf = etrack_freqs[ns1].size();
	      
	      if ( etrack_fspeckurt[ns1].size() == 0 ) 
		{
		  etrack_fspeckurt[ns1].resize( nf );
		  etrack_fspecskew[ns1].resize( nf );
		  etrack_fspeccv[ns1].resize( nf );
		}

	      // calc spec kurt
	      std::vector<double> fspsk ,fspcv ;
	      std::vector<double> fspku = skurt.kurtosis2_fbin( sklog, ns1, &fspcv, &fspsk ) ;
	      
              int fidx = 0;
              for ( int i = 0 ; i < mtm.f.size() ; i++ )
                if ( mtm.f[i] >= min_f && mtm.f[i] <= max_f )
                  {
                    // channel -> freq -> epoch/value 
                    etrack_fspeckurt[ ns1 ][ fidx ].push_back( fspku[i] );
		    etrack_fspecskew[ ns1 ][ fidx ].push_back( fspsk[i] );
		    etrack_fspeccv[ ns1 ][ fidx ].push_back( fspcv[i] );
		    ++fidx;
                  }
            
	    }
	  
	  

	  //
	  // add new signals?
	  //

	  if ( new_sigs )
	    {
	      
	      int fidx = 0;
	      for ( int i = 0 ; i < mtm.f.size() ; i++ )
		if ( mtm.f[i] >= min_f && mtm.f[i] <= max_f )
		  {
		    
		    const std::string new_sig_label = new_sig_prefix + "_" + signals.label(s) + "_" + Helper::int2str( fidx + 1 );
		    
		    std::vector<double> dat = eigen_ops::copy_array( addX.col( fidx ) );
		    
		    // normalize
		    for (int p=0; p<np; p++) dat[p] /= (double)addn[p];
		    
		    if ( dat.size() != np ) Helper::halt( "internal error in MTM 'add'" );
		    
		    logger << "  adding new signal " << new_sig_label << " ( MTM @ " << mtm.f[i] << " Hz )\n";
		    
		    edf.add_signal( new_sig_label , Fs[s] , dat );
		    
		    ++fidx;
		  }
	    }


	  //
	  // at tis point, we're done w/ console outputs
	  //

	  show_initial_log_output = false;

	  
	} // next signal

      writer.unlevel( globals::signal_strat );


      //
      // Spectral kurtosis (calc'ed per epoch/interval but also here averaging channels
      //
      
      if ( spec_kurt )
	{
	  // average values over channels
	  skurt.average_channels();

	  //
	  // output now?
	  //

	  if ( ( ! epochwise ) || epoch_level_output )
	    {
	      if ( bands )
		{
		  std::set<frequency_band_t>::const_iterator bb = skurt.bands.begin();
		  while ( bb != skurt.bands.end() )
		    {
		      writer.level( globals::band( *bb ) , globals::band_strat );
		      
		      double spsk , spcv ; 
		      double spku = kurt_altdef ? skurt.kurtosis( *bb , &spcv, &spsk ) : skurt.kurtosis2( *bb , &spcv, &spsk ) ;
		      
		      if ( spku > -900 )
			{
			  writer.value( "SPECCV" , spcv );
			  writer.value( "SPECSKEW" , spsk );
			  writer.value( "SPECKURT" , spku );
			}
		      
		      ++bb;
		    }
		  writer.unlevel( globals::band_strat );
		}
	    }

	  //
	  // track?
	  //
	  
	  // first, size up
	  if ( bands )
	    {
	      if ( etrack_chavg_speckurt.size() == 0 )
		{
		  etrack_chavg_speckurt.resize( skurt.bands.size() );
		  etrack_chavg_specskew.resize( skurt.bands.size() );
		  etrack_chavg_speccv.resize( skurt.bands.size() );	      
		}
	      
	      int bn=0;
	      std::set<frequency_band_t>::const_iterator bb = skurt.bands.begin();
	      while ( bb != skurt.bands.end() )
		{
		  
		  double spsk ,spcv ;
		  double spku = kurt_altdef ? skurt.kurtosis( *bb , &spcv, &spsk ) : skurt.kurtosis2( *bb , &spcv, &spsk ) ;
		  
		  etrack_chavg_speckurt[bn].push_back( spku );
		  etrack_chavg_speccv[bn].push_back( spcv );
		  etrack_chavg_specskew[bn].push_back( spsk );
		  
		  ++bb;
		  ++bn;
		}
	    }
	}

      
      //
      // All done / next epoch?
      //
      
      if ( ! epochwise )
	break;
      
    } // next epoch
  
  if ( epoch_level_output )
    writer.unepoch();  
  
  logger << "\n";


  //
  // summarized whole-night outputs (based on epoch-level results)
  //

  if ( epochwise )
    {

      //
      // Output mean spectra
      //

      int ns1 = -1;
      for ( int s=0; s<ns; s++)
	{
	  if ( edf.header.is_annotation_channel( signals(s) ) )
	    continue;
	  if ( min_sr && Fs[s] < min_sr )
	    continue;

	  // output for this channel
	  ++ns1;

          writer.level( signals.label(s) , globals::signal_strat );

	  //
	  // basic MTM spectra
	  //
	  
	  int fn = etrack_freqs[ns1].size();
	  for (int fi=0; fi<fn; fi++)
	    {
	      writer.level( etrack_freqs[ns1][fi]  , globals::freq_strat  );

	      
	      const double pmean = MiscMath::mean( etrack_power[ns1][fi] );

	      const double pmed  = MiscMath::median( etrack_power[ns1][fi] );
	      const double psd   = MiscMath::sdev( etrack_power[ns1][fi] );
	      
	      writer.value( "MTM" , pmean );

	      // allow for variable length epochs, so report weighted mean
	      // -- TOOD - should add this for all outputs...  thus why for
	      //           now we add WMTM as a reminder that this is incomplete...
	      
	      if (edf.timeline.generic_epochs() )
		{
		  const double wpmean = MiscMath::weighted_mean( etrack_power[ns1][fi] , etrack_length );
		  writer.value( "WMTM" , wpmean );
		}
	      
	      if ( etrack_power[ns1][fi].size() > 2 ) 
		{
		  writer.value( "MTM_MD" , pmed );
		  writer.value( "MTM_SD" , psd );
		}
	    }
	  writer.unlevel( globals::freq_strat );

	  //
	  // band power
	  //

	  if ( bands )
	    {
	      
	      int bn = bandaid.size();
	      
	      for (int bi=0; bi<bn; bi++)
		{
		  writer.level( globals::band( bandaid.bands[bi] ) , globals::band_strat  );
		  const double pmean = MiscMath::mean( etrack_bandpower[ns1][bi] );
		  
		  const double rpmean = MiscMath::mean( etrack_relbandpower[ns1][bi] );
		  
		  writer.value( "MTM" , pmean );
		  writer.value( "REL" , rpmean );
		  
		  if ( etrack_bandpower[ns1][bi].size() > 2 ) 
		    {
		      const double pmed  = MiscMath::median( etrack_bandpower[ns1][bi] );
		      const double psd   = MiscMath::sdev( etrack_bandpower[ns1][bi] );
		      writer.value( "MTM_MD" , pmed );
		      writer.value( "MTM_SD" , psd );
		    }
		  
		  if ( etrack_relbandpower[ns1][bi].size() > 2 )
		    {
		      const double rpmed  = MiscMath::median( etrack_relbandpower[ns1][bi] );
		      const double rpsd   = MiscMath::sdev( etrack_relbandpower[ns1][bi] );
		      writer.value( "REL_MD" , rpmed );
		      writer.value( "REL_SD" , rpsd );
		    }
		}
	      writer.unlevel( globals::band_strat );
	    }

	  
	  //
	  // spectral kurtsosis (channel-specific version)
	  //
	
	  if ( spec_kurt && ! kurt_altdef )
	    {

	      // initial dummy, to get bands (defined in constructor)
	      spectral_kurtosis_t skurt;

	      if ( bands )
		{
		  int bn=0;
		  std::set<frequency_band_t>::const_iterator bb = skurt.bands.begin();
		  while ( bb != skurt.bands.end() )
		    {
		      writer.level( globals::band( *bb ) , globals::band_strat  );
		      
		      if ( MiscMath::mean( etrack_speckurt[ns1][bn] )  > -900 )
			{
			  writer.value( "SPECKURT" , MiscMath::mean( etrack_speckurt[ns1][bn] ) );
			  writer.value( "SPECKURT_MD" , MiscMath::median( etrack_speckurt[ns1][bn] ) );
			  
			  writer.value( "SPECSKEW" ,  MiscMath::mean( etrack_specskew[ns1][bn] ) );
			  writer.value( "SPECSKEW_MD" , MiscMath::median( etrack_specskew[ns1][bn] ) );
			  
			  writer.value( "SPECCV" ,    MiscMath::mean( etrack_speccv[ns1][bn] ) );
			  writer.value( "SPECCV_MD" , MiscMath::median( etrack_speccv[ns1][bn] ) );
			}
		      
		      ++bb;
		      ++bn;
		    }
		  writer.unlevel( globals::band_strat );
		}

	      // freq-bin
	      int fn = etrack_freqs[ns1].size();
	      for (int fi=0; fi<fn; fi++)
		{
		  writer.level( etrack_freqs[ns1][fi]  , globals::freq_strat  );

		  const double sk_mn = MiscMath::mean( etrack_fspeckurt[ns1][fi] );
		  const double ss_mn = MiscMath::mean( etrack_fspecskew[ns1][fi] );
		  const double sc_mn = MiscMath::mean( etrack_fspeccv[ns1][fi] );
		  
		  writer.value( "SPECKURT" , MiscMath::mean( etrack_fspeckurt[ns1][fi] ) );
		  writer.value( "SPECSKEW" , MiscMath::mean( etrack_fspecskew[ns1][fi] ) );
		  writer.value( "SPECCV" , MiscMath::mean( etrack_fspeccv[ns1][fi] ) );
		  
		  if ( etrack_power[ns1][fi].size() > 2 )
		    {
		      writer.value( "SPECKURT_MD" , MiscMath::median( etrack_fspeckurt[ns1][fi] ) );
		      writer.value( "SPECSKEW_MD" , MiscMath::median( etrack_fspecskew[ns1][fi] ) );
		      writer.value( "SPECCV_MD" , MiscMath::median( etrack_fspeccv[ns1][fi] ) );
		    }
		}
	      writer.unlevel( globals::freq_strat );
	      
	    }
	  
	  
	  //
	  // band ratios
	  //
	  
	  if ( bands && calc_ratio )
	    {
	      // ratio=ALPHA/BETA,THETA/DELTA,...
	      std::vector<std::string> r = Helper::parse( Helper::toupper( ratios ) , ',' );
	      
	      bool done_any = false; 
	      int rn = 0;
	      
	      std::vector<std::string>::const_iterator rr = r.begin();
	      while ( rr != r.end() )
		{
		  // A/B
		  std::vector<std::string> tok = Helper::parse( *rr , '/' );
		  frequency_band_t b1 = globals::band( tok[0] );
		  frequency_band_t b2 = globals::band( tok[1] );

		  const double rmean = MiscMath::mean( etrack_bandratios[ns1][rn] );
		  
		  writer.level( tok[0] , "B1" );
		  writer.level( tok[1] , "B2" );
		  writer.value( "RATIO" , rmean );
		  

		  if ( etrack_bandratios[ns1][rn].size() > 2 ) 
		    {
		      const double rmed  = MiscMath::median( etrack_bandratios[ns1][rn] );
		      const double rsd   = MiscMath::sdev( etrack_bandratios[ns1][rn] );
		      writer.value( "RATIO_MD" , rmed );
		      writer.value( "RATIO_SD" , rsd );
		    }
		  
		  done_any = true;
		  
		  ++rn;
		  ++rr;
		}
	      
	      if ( done_any )
		{
		  writer.unlevel( "B2" );
		  writer.unlevel( "B1" );
		}
	    }

	  
	  
	  
	} // next channel
      
      writer.unlevel( globals::signal_strat );


      //
      // spectral kurtsosis (avg over channels variant alternate-speckurt)
      //
      
      if ( spec_kurt && kurt_altdef && bands )
	{

	  // initial dummy, to get bands (defined in constructor)
	  spectral_kurtosis_t skurt;
	  
	  int bn=0;
	  std::set<frequency_band_t>::const_iterator bb = skurt.bands.begin();
	  while ( bb != skurt.bands.end() )
	    {
	      writer.level( globals::band( *bb ) , globals::band_strat  );

	      if ( MiscMath::mean( etrack_chavg_speckurt[bn] ) > -900 )
		{
		  writer.value( "SPECKURT" , MiscMath::mean( etrack_chavg_speckurt[bn] ) );
		  writer.value( "SPECKURT_MD" , MiscMath::median( etrack_chavg_speckurt[bn] ) );
		  
		  writer.value( "SPECSKEW" ,  MiscMath::mean( etrack_chavg_specskew[bn] ) );
		  writer.value( "SPECSKEW_MD" , MiscMath::median( etrack_chavg_specskew[bn] ) );
		  
		  writer.value( "SPECCV" ,    MiscMath::mean( etrack_chavg_speccv[bn] ) );
		  writer.value( "SPECCV_MD" , MiscMath::median( etrack_chavg_speccv[bn] ) );
		}
	      
	      ++bb;
	      ++bn;
	    }
	  writer.unlevel( globals::band_strat );
	}
      
      
    } // end of summary outputs 

  // all done
}



