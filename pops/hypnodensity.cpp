
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

#include "pops/hypnodensity.h"
#include "pops/pops.h"
#include "pops/indiv.h"
#include "pops/options.h"
#include "pops/spec.h"

#include "edf/edf.h"
#include "annot/annot.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "param.h"
#include "defs/defs.h"

#include <map>
#include <numeric>
#include <cmath>

extern logger_t logger;


void pops_hypnodensity( edf_t & edf , param_t & param )
{

  //
  // Require continuous EDF
  //

  // ORIGINAL — too strict: rejected any EDF+D including single-segment post-EDGER recordings:
  //   if ( ! edf.header.continuous )
  //     Helper::halt( "POPS hypnodensity mode requires a continuous EDF (EDF+D not supported)" );
  //
  // CHANGED: only reject recordings with genuine gaps between segments.  A single-segment
  // EDF+D (e.g. produced by EDGER trimming leading/trailing wake) is effectively continuous;
  // epoch_offset_tp is now applied within each segment in calc_epochs() (see timeline/epochs.cpp).
  // To revert: uncomment the two original lines above and remove the two lines below.
  if ( edf.is_actually_discontinuous() )
    Helper::halt( "POPS hypnodensity mode requires a single-segment EDF (no within-recording gaps)" );


  //
  // Parse parameters
  //

  const int N = param.value( "hypnodensity" ) != "" ? param.requires_int( "hypnodensity" ) : 6;

  if ( N < 1 || 30 % N != 0 )
    Helper::halt( "hypnodensity=N requires N to divide 30 evenly (valid: 1,2,3,5,6,10,15,30)" );

  if ( pops_opt_t::n_stages == 3 )
    Helper::halt( "hypnodensity mode requires a 5-class model (not 3-class)" );

  const double stride_sec = 30.0 / N;

  //
  // Check EDF record duration is compatible with the requested stride before doing anything else
  //
  // add_signal() requires n_samples_per_record = record_duration / stride_sec to be a positive integer.
  // Common issue: 1-second EDF records can't hold sub-Hz signals (e.g. 1s / 5s = 0.2).
  // Use RECORD-SIZE to adjust (e.g. RECORD-SIZE dur=30).
  //

  const double n_spr_d = edf.header.record_duration / stride_sec;
  const int n_spr = (int)std::lround( n_spr_d );
  if ( std::fabs( n_spr_d - n_spr ) > 1e-6 || n_spr < 1 )
    Helper::halt( "EDF record duration (" + Helper::dbl2str( edf.header.record_duration )
		  + "s) is not divisible by hypnodensity stride (" + Helper::dbl2str( stride_sec )
		  + "s); use RECORD-SIZE (e.g. dur=30) before running POPS hypnodensity" );

  // sample index offset: use (N-1)/2 (integer division) rather than N/2
  // This gives floor((N-1)/2) leading edge samples and ceil((N-1)/2) trailing edge samples —
  // the most balanced achievable split, with any 1-sample asymmetry falling at the trailing end.
  // Examples: N=6 -> offset=2 (+10s), N=3 -> offset=1 (+10s, symmetric), N=5 -> offset=2 (+12s, symmetric)
  const int mid_offset = ( N - 1 ) / 2;

  const std::string prefix = param.has( "prefix" ) ? param.value( "prefix" ) : "PP" ;

  // add-nrem123 (individual N1/N2/N3 channels): default true
  const bool add_nrem123 = param.has( "add-nrem123" ) ? param.yesno( "add-nrem123" ) : true ;

  // add-nrem (summed NR = N1+N2+N3 channel): default false
  const bool add_nrem = param.has( "add-nrem" ) ? param.yesno( "add-nrem" ) : false ;

  // equiv channel combination parameters (mirrors existing POPS logic)
  const int equivn = pops_opt_t::equivs.size();
  const double comb_conf = param.has( "conf" ) ? param.requires_dbl( "conf" ) : 0 ;
  const int comb_method = param.has( "mean" ) ? 3 : param.has( "geo" ) ? 2 : 1 ;

  // feature ranges
  const double range_th   = param.has( "ranges-th" )   ? param.requires_dbl( "ranges-th" )   : 4 ;
  const double range_prop = param.has( "ranges-prop" )  ? param.requires_dbl( "ranges-prop" )  : 0.33 ;

  if ( ! add_nrem123 && ! add_nrem )
    Helper::halt( "hypnodensity: at least one of add-nrem123 (default T) or add-nrem must be true" );

  // model iterations
  int num_iter = 0;
  if      ( param.has( "iterations" ) ) num_iter = param.requires_int( "iterations" );
  else if ( param.has( "iter" ) )       num_iter = param.requires_int( "iter" );


  // Fs code for add_signal: negative means n_samples_per_record directly
  const double Fs_code = -(double)n_spr;

  const int total_samples = edf.header.nr * n_spr;

  const uint64_t stride_tp = (uint64_t)( stride_sec * globals::tp_1sec );


  //
  // Load model if not already cached
  //

  std::string model_file = ".";
  if ( param.has( "model" ) )
    model_file = param.value( "model" );
  else if ( pops_opt_t::pops_root != "" )
    model_file = pops_opt_t::pops_root + ".mod";
  if ( model_file != "." )
    model_file = pops_t::update_filepath( model_file );

  if ( model_file == "." )
    Helper::halt( "POPS hypnodensity requires a model file, via model or lib args" );

  if ( pops_t::lgbm_model_loaded != model_file )
    {
      pops_t::lgbm.load_model( model_file );
      pops_t::lgbm_model_loaded = model_file;
    }


  //
  // Load feature ranges if present
  //

  std::string ranges_file = ".";
  if ( param.has( "ranges" ) )
    ranges_file = param.value( "ranges" );
  else if ( pops_opt_t::pops_root != "" && pops_opt_t::if_root_apply_ranges )
    {
      ranges_file = pops_t::update_filepath( pops_opt_t::pops_root + ".ranges" );
      if ( ! Helper::fileExists( ranges_file ) ) ranges_file = ".";
    }

  if ( ranges_file != "." )
    pops_t::read_ranges( ranges_file );


  //
  // Validate equiv channel originals (same check as existing POPS)
  //

  if ( equivn )
    {
      const std::map<std::string,std::string> & em = pops_opt_t::equivs[0];
      for ( auto ee = em.begin(); ee != em.end(); ++ee )
	{
	  if ( pops_t::specs.chs.find( ee->first ) == pops_t::specs.chs.end() )
	    Helper::halt( "could not find root equivalence channel in feature spec: " + ee->first );
	}
    }


  //
  // Result store: output sample index -> 5-stage posterior vector
  //

  std::map<int,Eigen::VectorXd> sample_to_post;


  //
  // Reusable indiv object (minimal init; fields populated per stride below)
  //

  pops_indiv_t indiv( &edf );


  //
  // Per-stride loop
  //

  logger << "  hypnodensity mode: N=" << N
	 << " strides, stride=" << stride_sec << "s"
	 << ", output sample rate=1/" << (int)stride_sec << "Hz"
	 << ", " << n_spr << " sample(s) per EDF record\n";

  if ( equivn )
    {
      logger << "  equivalence mode: " << equivn << " channel set(s); method = ";
      if ( comb_method == 1 ) logger << "most confident per epoch";
      else if ( comb_method == 2 ) logger << "geometric mean";
      else logger << "confidence-weighted mean";
      if ( comb_method != 1 ) logger << ", minimum conf score = " << comb_conf;
      logger << "\n";
    }

  for ( int k = 0; k < N; k++ )
    {
      const uint64_t offset_tp = (uint64_t)( k * stride_sec * globals::tp_1sec );

      edf.timeline.set_epoch( 30.0 , 30.0 , offset_tp );

      logger << "\n  stride " << k << " (offset=" << k*(int)stride_sec
	     << "s): " << edf.timeline.num_epochs() << " epochs\n";

      if ( edf.timeline.num_epochs() == 0 )
	continue;


      //
      // Equivalence channel loop (runs once when equivn == 0)
      //

      std::vector<pops_sol_t> sols;
      int eq1 = 0;

      while ( true )
	{

	  //
	  // Set up equivalence channel swap for this iteration
	  //

	  if ( equivn )
	    {
	      if ( eq1 == equivn ) break;

	      pops_opt_t::equiv_swapins = pops_opt_t::equivs[ eq1 ];
	      ++eq1;

	      pops_opt_t::equiv_label = "";
	      for ( auto ee = pops_opt_t::equiv_swapins.begin();
		    ee != pops_opt_t::equiv_swapins.end(); ++ee )
		{
		  if ( ee != pops_opt_t::equiv_swapins.begin() )
		    pops_opt_t::equiv_label += ";";
		  pops_opt_t::equiv_label += ee->first + "->" + ee->second;
		}
	      logger << "    using equivalence set " << eq1 << "/" << equivn
		     << ": " << pops_opt_t::equiv_label << "\n";

	      // for actual swaps (not self-mapping), check all channels present
	      if ( eq1 != 1 )
		{
		  bool all_present = true;
		  for ( auto & ee : pops_opt_t::equiv_swapins )
		    if ( ! edf.header.has_signal( ee.second ) )
		      { all_present = false; break; }
		  if ( ! all_present )
		    {
		      logger << "  ** skipping equivalence set (missing channel)\n";
		      continue;
		    }
		}
	    }


	  //
	  // Level 1: extract per-epoch features
	  //

	  indiv.level1( edf );

	  // Populate E, S, Sorig (not done by level1 -- normally done by staging())
	  // E: 0-based epoch indices; S/Sorig: all WAKE (ignored in output, needed by combine())
	  indiv.E.resize( indiv.ne );
	  std::iota( indiv.E.begin() , indiv.E.end() , 0 );
	  indiv.S.assign( indiv.ne , POPS_WAKE );
	  indiv.Sorig.assign( indiv.ne , POPS_WAKE );
	  indiv.has_staging = false;


	  //
	  // Level 2: temporal smoothing, normalisation, SVD projection
	  //

	  indiv.level2( ! pops_opt_t::verbose );


	  //
	  // Optionally clamp out-of-range features
	  //

	  if ( ranges_file != "." )
	    indiv.apply_ranges( range_th , range_prop );


	  //
	  // Predict: fills indiv.P (ne x 5) and indiv.PS
	  //

	  indiv.predict( num_iter );


	  //
	  // Save solution for equiv combining
	  //

	  if ( equivn )
	    {
	      pops_sol_t sol;
	      sol.E = indiv.E;
	      sol.S = indiv.S;
	      sol.P = indiv.P;
	      sols.push_back( sol );
	    }


	  if ( equivn == 0 ) break;

	} // end equiv loop


      //
      // Combine equivalence channel solutions (if any)
      //

      if ( equivn && sols.size() > 0 )
	{
	  logger << "  combining " << sols.size()
		 << " equivalence-channel solutions for stride " << k << "\n";
	  indiv.combine( sols , comb_method , comb_conf );
	}

      if ( equivn )
	logger << "\n";


      //
      // Map each epoch's posterior to its output sample index.
      //
      // sidx = window_start_sample + mid_offset
      //      = (k + epoch_m * N) + (N-1)/2
      //
      // This places the posterior at floor((N-1)/2) stride-steps into the window,
      // giving floor((N-1)/2) leading edge samples and ceil((N-1)/2) trailing edge
      // samples — the most balanced split achievable, with any 1-sample asymmetry
      // falling at the trailing end.  All arithmetic is exact integer (no rounding).
      //

      for ( int e = 0; e < (int)indiv.E.size(); e++ )
	{
	  interval_t interval = edf.timeline.epoch( indiv.E[e] );
	  // Compute sample index using pure integer stride arithmetic.
	  // interval.start is always an exact multiple of stride_tp, so the
	  // division is exact (no rounding).  mid_offset is in stride steps.
	  const int window_start_sample = (int)( interval.start / stride_tp );
	  const int sidx = window_start_sample + mid_offset;

	  if ( sidx >= 0 && sidx < total_samples )
	    sample_to_post[ sidx ] = indiv.P.row( e );
	}

    } // end stride loop

  // restore default 30s non-overlapping epochs
  edf.timeline.set_epoch( 30.0 , 30.0 );


  //
  // Build output signal vectors via zero-order hold
  //
  // Leading edge:  samples before the first valid sample — filled with the globally
  //                first valid posterior (temporally nearest).
  // Trailing edge: samples after the last valid sample — filled by carrying forward
  //                the last-seen valid posterior.
  // No-data fallback: uniform 0.2 across all 5 stages (only reached if recording
  //                   produced zero epochs for every stride).
  //

  // Uniform fallback (no valid data at all)
  Eigen::VectorXd uniform = Eigen::VectorXd::Constant( 5 , 0.2 );

  // Globally first valid posterior (used for leading-edge ZOH)
  Eigen::VectorXd global_first_post = uniform;
  const bool have_any = ! sample_to_post.empty();
  if ( have_any )
    global_first_post = sample_to_post.begin()->second;

  const int first_valid_sample = have_any ? sample_to_post.begin()->first  : -1;
  int       last_valid_sample  = have_any ? sample_to_post.rbegin()->first : -1;

  std::vector<double> sigW ( total_samples , 0.0 );
  std::vector<double> sigR ( total_samples , 0.0 );
  std::vector<double> sigN1( total_samples , 0.0 );
  std::vector<double> sigN2( total_samples , 0.0 );
  std::vector<double> sigN3( total_samples , 0.0 );
  std::vector<double> sigNR( total_samples , 0.0 );

  Eigen::VectorXd last_seen_post = global_first_post;

  for ( int m = 0; m < total_samples; m++ )
    {
      Eigen::VectorXd post;

      if ( sample_to_post.count( m ) )
	{
	  post = sample_to_post[ m ];
	  last_seen_post = post;
	}
      else if ( ! have_any )
	{
	  post = uniform;
	}
      else if ( m < first_valid_sample )
	{
	  // leading edge: use globally first valid posterior
	  post = global_first_post;
	}
      else
	{
	  // trailing edge: carry forward last seen valid posterior
	  post = last_seen_post;
	}

      sigW [ m ] = post[ POPS_WAKE ];
      sigR [ m ] = post[ POPS_REM  ];
      sigN1[ m ] = post[ POPS_N1   ];
      sigN2[ m ] = post[ POPS_N2   ];
      sigN3[ m ] = post[ POPS_N3   ];
      sigNR[ m ] = post[ POPS_N1   ] + post[ POPS_N2 ] + post[ POPS_N3 ];
    }


  //
  // Add signals to EDF
  //

  if ( add_nrem123 )
    {
      edf.add_signal( prefix + "_N1" , Fs_code , sigN1 );
      edf.add_signal( prefix + "_N2" , Fs_code , sigN2 );
      edf.add_signal( prefix + "_N3" , Fs_code , sigN3 );
      edf.add_signal( prefix + "_R"  , Fs_code , sigR  );
      edf.add_signal( prefix + "_W"  , Fs_code , sigW  );
    }

  if ( add_nrem )
    edf.add_signal( prefix + "_NR" , Fs_code , sigNR );

  logger << "  added hypnodensity channels with prefix '" << prefix << "'\n";


  //
  // Edge annotation: mark ZOH-filled leading and trailing regions
  //

  if ( first_valid_sample > 0 || last_valid_sample < total_samples - 1 )
    {
      annot_t * aedge = edf.annotations->add( prefix + "/edge" );
      aedge->description = "POPS hypnodensity ZOH-filled edge region";

      if ( first_valid_sample > 0 )
	{
	  // [0 , first_valid_sample * stride_tp)
	  const uint64_t tp_stop = (uint64_t)first_valid_sample * stride_tp;
	  aedge->add( "." , interval_t( 0 , tp_stop ) , "." );
	}

      if ( last_valid_sample >= 0 && last_valid_sample < total_samples - 1 )
	{
	  // [(last_valid_sample+1) * stride_tp , total_duration_tp)
	  const uint64_t tp_start = (uint64_t)( last_valid_sample + 1 ) * stride_tp;
	  aedge->add( "." , interval_t( tp_start , edf.timeline.total_duration_tp ) , "." );
	}

      logger << "  edge annotation '" << prefix << "/edge' added"
	     << " (leading: " << first_valid_sample * (int)stride_sec << "s"
	     << ", trailing from sample " << last_valid_sample + 1 << ")\n";
    }

}

#endif
