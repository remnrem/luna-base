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

//
// COMBINE-EMG : build a single continuous EMG channel from 2+ candidate
// channels (raw or already bipolar-referenced) by scoring short windows
// for artifact/quality (flatline, clipping, excess RMS, line-noise,
// low-frequency drift, burstiness) and switching to whichever candidate
// is cleanest, with hysteresis/minimum-dwell/smoothing to avoid flapping
// and a local-normalization + crossfade + clamp stitching step so the
// output doesn't itself look like a string of artifacts.
//
// Selection never leaves a window unassigned: flat/clip/excess-RMS/
// line-noise are soft (per-candidate, per-window) penalties on a
// per-candidate z-scored composite score, not hard vetoes -- so even if
// every candidate looks bad in a given window, the least-bad one is still
// picked and stitched in (continuously, no gaps in the output).
//

#include "dsp/combine_emg.h"

#include "param.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "db/db.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "annot/annot.h"
#include "defs/defs.h"
#include "miscmath/miscmath.h"
#include "fftw/fftwrap.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

extern writer_t writer;
extern logger_t logger;

namespace {

  // one contiguous EDF segment, mapped onto the concatenated whole-trace
  // sample index space shared by all (same-SR) candidate channels
  struct seg_t {
    interval_t iv;
    int start_idx;
    int n_samples;
  };

  static std::vector<seg_t> build_segments( edf_t & edf , const int sig_n )
  {
    std::vector<seg_t> segs;
    std::set<interval_t> ivs = edf.timeline.segments();
    int offset = 0;
    for ( std::set<interval_t>::const_iterator ss = ivs.begin() ; ss != ivs.end() ; ++ss )
      {
        slice_t slice( edf , sig_n , *ss );
        const std::vector<double> * pdata = slice.pdata();
        const int n = pdata == NULL ? 0 : (int)pdata->size();
        if ( n <= 0 ) continue;
        seg_t s;
        s.iv = *ss;
        s.start_idx = offset;
        s.n_samples = n;
        segs.push_back( s );
        offset += n;
      }
    return segs;
  }

  static std::vector<double> read_whole( edf_t & edf , const int sig_n )
  {
    slice_t slice( edf , sig_n , edf.timeline.wholetrace() );
    const std::vector<double> * pdata = slice.pdata();
    return pdata == NULL ? std::vector<double>() : *pdata;
  }

  // ensure a signal is expressed in microvolts (so amplitude thresholds are
  // meaningful); mirrors dsptools::qc_t::ensure_units()
  static void ensure_uV( edf_t & edf , int s )
  {
    const std::string pdim = Helper::trim( edf.header.phys_dimension[s] );
    if ( Helper::imatch( pdim , "uV" ) ) return;
    const bool is_mV = Helper::imatch( pdim , "mV" );
    const bool is_V  = Helper::imatch( pdim , "V" );
    if ( is_mV || is_V )
      edf.rescale( s , "uV" );
    else
      logger << "  ** warning: " << edf.header.label[s]
             << " physical dimension is '" << ( pdim.empty() ? "(empty)" : pdim )
             << "' -- not a recognised voltage unit (uV/mV/V);"
             << " amplitude thresholds may not apply\n";
  }

  // non-enforcing stage lookup: only looks for existing stage annotations,
  // never invokes/derives automated staging
  static bool any_stage_annotations( edf_t & edf )
  {
    static const char * labels[] = { "N1" , "N2" , "N3" , "N4" , "R" , "REM" , "W" };
    for ( int i = 0 ; i < 7 ; ++i )
      if ( edf.annotations->find( labels[i] ) != NULL ) return true;
    return false;
  }

  static std::vector< std::pair<uint64_t,uint64_t> > stage_intervals( edf_t & edf , const std::vector<std::string> & labels )
  {
    std::vector< std::pair<uint64_t,uint64_t> > ivs;
    for ( int i = 0 ; i < (int)labels.size() ; ++i )
      {
        annot_t * a = edf.annotations->find( labels[i] );
        if ( a == NULL ) continue;
        annot_map_t evts = a->extract( edf.timeline.wholetrace() );
        for ( annot_map_t::const_iterator aa = evts.begin() ; aa != evts.end() ; ++aa )
          ivs.push_back( std::make_pair( aa->first.interval.start , aa->first.interval.stop ) );
      }
    std::sort( ivs.begin() , ivs.end() );
    return ivs;
  }

  static bool tp_in_intervals( const std::vector< std::pair<uint64_t,uint64_t> > & ivs , uint64_t tp )
  {
    if ( ivs.empty() ) return false;
    std::vector< std::pair<uint64_t,uint64_t> >::const_iterator it =
      std::upper_bound( ivs.begin() , ivs.end() ,
                         std::make_pair( tp , std::numeric_limits<uint64_t>::max() ) );
    if ( it != ivs.begin() )
      {
        --it;
        if ( tp >= it->first && tp < it->second ) return true;
      }
    return false;
  }

  static double median_abs_dev( const std::vector<double> & v , const double med )
  {
    std::vector<double> dev( v.size() );
    for ( size_t i = 0 ; i < v.size() ; ++i ) dev[i] = std::fabs( v[i] - med );
    return MiscMath::median( dev );
  }

  // per-window raw quality metrics for a single candidate channel, computed
  // once (phase A) and reused once each candidate's own log-broadband-power
  // distribution is known (phase B), so the composite score compares
  // candidates on a scale-consistent (per-candidate z-scored) basis rather
  // than raw, channel-amplitude-dependent power
  struct win_metrics_t {
    bool flat, clip, hi_rms, hi_ln;
    bool has_spectral;
    double log_broad, line_ratio, drift_ratio, kurt_bonus;
    bool is_rem;
  };

  // one contiguous segment's decision windows + per-candidate raw metrics
  struct seg_work_t {
    seg_t seg;
    std::vector< std::pair<int,int> > windows;
    std::vector< std::vector<win_metrics_t> > M; // [candidate][window]
  };

} // anonymous namespace


combine_emg_t::combine_emg_t( edf_t & edf , param_t & param )
{

  //
  // -----------------------------------------------------------------------
  // Parameters
  // -----------------------------------------------------------------------
  //

  const std::string sigstr = param.requires( "sig" );
  signal_list_t sig0 = edf.header.signal_list( sigstr );
  if ( sig0.size() < 2 )
    Helper::halt( "requires 2+ signals in sig=" );

  const bool auto_pairs = param.yesno( "pairs" , false , true );

  const std::string new_label  = param.has( "new"   ) ? param.value( "new"   ) : "cEMG";
  const std::string annot_name = param.has( "annot" ) ? param.value( "annot" ) : "emg_src";

  const double win_s          = param.has( "win" )             ? param.requires_dbl( "win" )             : 20.0;
  const double flat_sd_th     = param.has( "flat-sd-th" )      ? param.requires_dbl( "flat-sd-th" )      : 0.5;
  const double flat_deriv_prop= param.has( "flat-deriv-prop" ) ? param.requires_dbl( "flat-deriv-prop" ) : 0.8;
  const double flat_deriv_eps = param.has( "flat-deriv-eps" )  ? param.requires_dbl( "flat-deriv-eps" )  : 1e-4;
  const double clip_prop_th   = param.has( "clip-prop" )       ? param.requires_dbl( "clip-prop" )       : 0.01;
  const double rms_th         = param.has( "rms-th" )          ? param.requires_dbl( "rms-th" )          : 300.0;
  const double line_bw        = param.has( "line-bw" )         ? param.requires_dbl( "line-bw" )         : 2.0;
  const double line_lo        = param.has( "line-lo" )         ? param.requires_dbl( "line-lo" )         : 10.0;
  const double line_hi        = param.has( "line-hi" )         ? param.requires_dbl( "line-hi" )         : 100.0;
  const double line_th        = param.has( "line-th" )         ? param.requires_dbl( "line-th" )         : 0.3;
  const double drift_hi       = param.has( "drift-hi" )        ? param.requires_dbl( "drift-hi" )        : 0.5;

  // penalties applied to the (z-scored) composite score when a criterion
  // is flagged -- soft, not absolute exclusions, so a window is never left
  // with no acceptable candidate; graded by how severe/irrecoverable the
  // problem typically is
  const double flat_penalty   = param.has( "flat-penalty" )    ? param.requires_dbl( "flat-penalty" )    : 10.0;
  const double clip_penalty   = param.has( "clip-penalty" )    ? param.requires_dbl( "clip-penalty" )    : 10.0;
  const double rms_penalty    = param.has( "rms-penalty" )     ? param.requires_dbl( "rms-penalty" )     : 2.0;
  const double ln_penalty     = param.has( "ln-penalty" )      ? param.requires_dbl( "ln-penalty" )      : 3.0;

  const double switch_margin  = param.has( "switch-margin" )   ? param.requires_dbl( "switch-margin" )   : 2.5;
  const double min_dwell_s    = param.has( "min-dwell" )       ? param.requires_dbl( "min-dwell" )       : 180.0;
  const int    smooth_win     = param.has( "smooth-win" )      ? param.requires_int( "smooth-win" )      : 5;
  const double xfade_s        = param.has( "xfade" )           ? param.requires_dbl( "xfade" )           : 1.0;
  const bool   do_norm        = param.yesno( "norm" , true , true );
  const double norm_win_s     = param.has( "norm-win" )        ? param.requires_dbl( "norm-win" )        : 60.0;
  const double clamp_z        = param.has( "clamp" )           ? param.requires_dbl( "clamp" )           : 8.0;
  const bool   use_staging_opt= param.yesno( "use-staging" , true , true );
  const double stage_weight   = param.has( "stage-weight" )    ? param.requires_dbl( "stage-weight" )    : 0.5;

  //
  // -----------------------------------------------------------------------
  // Candidate resolution: either sig= is already the candidate pool, or
  // (pairs=T) sig= lists raw electrodes and we auto-derive all C(n,2)
  // bipolar pairs via the same edf_t::reference() used by REFERENCE
  // -----------------------------------------------------------------------
  //

  std::vector<int> cand_n;
  std::vector<std::string> cand_lab;

  if ( auto_pairs )
    {
      const int n0 = sig0.size();
      for ( int i = 0 ; i < n0 ; ++i )
        for ( int j = i+1 ; j < n0 ; ++j )
          {
            signal_list_t s1; s1.signals.push_back( sig0(i) ); s1.signal_labels.push_back( sig0.label(i) );
            signal_list_t r1; r1.signals.push_back( sig0(j) ); r1.signal_labels.push_back( sig0.label(j) );
            const std::string new_lab = sig0.label(i) + "_" + sig0.label(j);
            edf.reference( s1 , r1 , true , new_lab , 0 , false , false );
            const int slot = edf.header.signal( new_lab );
            if ( slot == -1 )
              Helper::halt( "internal error creating pairwise reference " + new_lab );
            cand_n.push_back( slot );
            cand_lab.push_back( new_lab );
          }
      logger << "    pairs=T : derived " << cand_n.size() << " bipolar candidates from "
             << n0 << " raw electrodes\n";
    }
  else
    {
      for ( int i = 0 ; i < sig0.size() ; ++i )
        {
          cand_n.push_back( sig0(i) );
          cand_lab.push_back( sig0.label(i) );
        }
    }

  const int nc = (int)cand_n.size();

  logger << "    candidate channels : ";
  for ( int i = 0 ; i < nc ; ++i ) logger << ( i ? "," : "" ) << cand_lab[i];
  logger << "\n";

  const double fs = edf.header.sampling_freq( cand_n[0] );
  for ( int i = 1 ; i < nc ; ++i )
    if ( std::fabs( edf.header.sampling_freq( cand_n[i] ) - fs ) > 1e-4 )
      Helper::halt( "requires all candidate channels to have the same sample rate" );
  if ( fs <= 0 ) Helper::halt( "invalid sample rate" );

  for ( int i = 0 ; i < nc ; ++i ) ensure_uV( edf , cand_n[i] );

  //
  // -----------------------------------------------------------------------
  // Read whole-recording data for each candidate (once)
  // -----------------------------------------------------------------------
  //

  std::vector< std::vector<double> > cand_data( nc );
  for ( int i = 0 ; i < nc ; ++i ) cand_data[i] = read_whole( edf , cand_n[i] );

  const int N = (int)cand_data[0].size();
  if ( N == 0 ) { logger << "  ** no data read, skipping... \n"; return; }
  for ( int i = 1 ; i < nc ; ++i )
    if ( (int)cand_data[i].size() != N )
      Helper::halt( "internal error: candidate channel length mismatch" );

  std::vector<seg_t> segs = build_segments( edf , cand_n[0] );
  if ( segs.empty() ) { logger << "  ** no segments found, skipping...\n"; return; }

  //
  // -----------------------------------------------------------------------
  // Global representative scale, used to rescale the final z-stitched
  // output back to approximate physical (uV) units
  // -----------------------------------------------------------------------
  //

  double global_scale = 1.0;
  {
    std::vector<double> sigmas( nc );
    for ( int i = 0 ; i < nc ; ++i )
      {
        const double med = MiscMath::median( cand_data[i] );
        sigmas[i] = 1.4826 * median_abs_dev( cand_data[i] , med );
      }
    const double gm = MiscMath::median( sigmas );
    if ( gm > 0 ) global_scale = gm;
  }

  //
  // -----------------------------------------------------------------------
  // Optional, non-enforcing staging (REM figure-of-merit only)
  // -----------------------------------------------------------------------
  //

  const bool has_staging = any_stage_annotations( edf );
  const bool use_staging = use_staging_opt && has_staging;
  std::vector< std::pair<uint64_t,uint64_t> > rem_ivs;
  if ( use_staging )
    {
      std::vector<std::string> rem_labels; rem_labels.push_back( "R" ); rem_labels.push_back( "REM" );
      rem_ivs = stage_intervals( edf , rem_labels );
    }
  logger << "    staging (REM figure-of-merit) : "
         << ( use_staging ? "found, in use" : ( has_staging ? "found, disabled (use-staging=F)" : "not found" ) )
         << "\n";

  //
  // -----------------------------------------------------------------------
  // Window geometry
  // -----------------------------------------------------------------------
  //

  const int win_samples       = std::max( 1 , (int)std::lround( win_s * fs ) );
  const int norm_win_samples  = std::max( win_samples , (int)std::lround( norm_win_s * fs ) );
  const int xfade_samples     = std::max( 0 , (int)std::lround( xfade_s * fs ) );
  const int min_dwell_windows = std::max( 1 , (int)std::lround( min_dwell_s / win_s ) );

  // global accumulators (across all segments) for reporting
  std::vector<int>    g_ne( nc , 0 ) , g_cnt_used( nc , 0 );
  std::vector<int>    g_cnt_flat( nc , 0 ) , g_cnt_clip( nc , 0 ) , g_cnt_hi_rms( nc , 0 ) , g_cnt_hi_ln( nc , 0 );
  std::vector<double> g_sum_score( nc , 0.0 );
  int g_n_switch = 0;

  //
  // -----------------------------------------------------------------------
  // Phase A: build decision windows and raw (uncalibrated) per-candidate
  // metrics for every segment; accumulate each candidate's own log-
  // broadband-power distribution for later z-scoring
  // -----------------------------------------------------------------------
  //

  std::vector<seg_work_t> work( segs.size() );
  std::vector< std::vector<double> > cand_log_broad( nc );

  for ( size_t si = 0 ; si < segs.size() ; ++si )
    {
      const seg_t & seg = segs[si];
      seg_work_t & sw = work[si];
      sw.seg = seg;

      // partition this segment into non-overlapping decision windows;
      // merge an undersized final remainder into the previous window
      std::vector< std::pair<int,int> > & windows = sw.windows;
      {
        int a = seg.start_idx;
        const int seg_end = seg.start_idx + seg.n_samples;
        while ( a < seg_end )
          {
            const int b = std::min( a + win_samples , seg_end );
            windows.push_back( std::make_pair( a , b ) );
            a = b;
          }
        if ( windows.size() >= 2 )
          {
            std::pair<int,int> & last = windows.back();
            if ( ( last.second - last.first ) < win_samples / 2 )
              {
                std::pair<int,int> & prev = windows[ windows.size() - 2 ];
                prev.second = last.second;
                windows.pop_back();
              }
          }
      }

      const int nw = (int)windows.size();
      if ( nw == 0 ) continue;

      sw.M.assign( nc , std::vector<win_metrics_t>( nw ) );

      for ( int w = 0 ; w < nw ; ++w )
        {
          const int a = windows[w].first , b = windows[w].second;
          const int n = b - a;

          bool is_rem = false;
          if ( use_staging )
            {
              const int center_local = ( (a+b)/2 ) - seg.start_idx;
              const uint64_t tmid = seg.iv.start + (uint64_t)std::llround( (double)center_local / fs * globals::tp_1sec );
              is_rem = tp_in_intervals( rem_ivs , tmid );
            }

          for ( int c = 0 ; c < nc ; ++c )
            {
              std::vector<double> d( cand_data[c].begin() + a , cand_data[c].begin() + b );

              const double sd     = MiscMath::sdev( d );
              const double deriv  = MiscMath::flat( d , flat_deriv_eps );
              const bool f_flat   = ( sd < flat_sd_th ) || ( deriv >= flat_deriv_prop );
              const bool f_clip   = MiscMath::clipped( d ) >= clip_prop_th;
              const double rms    = MiscMath::rms( d );
              const bool f_rms    = rms > rms_th;

              double line_ratio = 0.0 , drift_ratio = 0.0 , log_broad = 0.0;
              bool f_ln = false , has_spectral = false;
              {
                const double nyquist = fs / 2.0;
                const double broad_hi = std::min( line_hi , nyquist - 0.5 );
                const double seg_sec = 4.0 , step_sec = 2.0;
                const int seg_pts = (int)( seg_sec * fs );
                const int ovl_pts = (int)( step_sec * fs );
                if ( broad_hi > line_lo && seg_pts <= n && seg_pts > ovl_pts )
                  {
                    const int noverlap_segments = (int)std::floor( ( n - ovl_pts ) / (double)( seg_pts - ovl_pts ) );
                    if ( noverlap_segments >= 1 )
                      {
                        PWELCH pwelch( d , fs , seg_sec , noverlap_segments , WINDOW_TUKEY50 );
                        const double p_broad = pwelch.psdsum( line_lo , broad_hi );
                        if ( p_broad > 0 )
                          {
                            log_broad = std::log( p_broad );
                            has_spectral = true;
                            const double f50_lo = 50.0 - line_bw , f50_hi = 50.0 + line_bw;
                            const double f60_lo = 60.0 - line_bw , f60_hi = 60.0 + line_bw;
                            const double p50 = f50_hi < nyquist ? pwelch.psdsum( f50_lo , f50_hi ) : 0.0;
                            const double p60 = f60_hi < nyquist ? pwelch.psdsum( f60_lo , f60_hi ) : 0.0;
                            const double p_line = std::max( p50 , p60 );
                            line_ratio = p_line / p_broad;
                            f_ln = ( p50 > 0 || p60 > 0 ) && line_ratio > line_th;
                            const double drift_hi_capped = std::min( drift_hi , broad_hi );
                            if ( drift_hi_capped > 0 )
                              drift_ratio = pwelch.psdsum( 0.0 , drift_hi_capped ) / p_broad;
                          }
                      }
                  }
                // both ratios compare power across different bands (not a
                // fraction of a shared total), so they are not naturally
                // bounded to [0,1] -- e.g. drift_ratio compares low-frequency
                // power against a broadband reference that excludes low
                // frequencies, and can be very large for signals with
                // substantial low-frequency content; clip so a single
                // extreme ratio can't dominate the composite score
                line_ratio  = std::min( line_ratio  , 5.0 );
                drift_ratio = std::min( drift_ratio , 5.0 );
              }

              const double kurt = MiscMath::kurtosis( d );
              const double kurt_bonus = std::max( 0.0 , kurt - 3.0 );

              win_metrics_t & mm = sw.M[c][w];
              mm.flat = f_flat; mm.clip = f_clip; mm.hi_rms = f_rms; mm.hi_ln = f_ln;
              mm.has_spectral = has_spectral;
              mm.log_broad = log_broad; mm.line_ratio = line_ratio; mm.drift_ratio = drift_ratio;
              mm.kurt_bonus = kurt_bonus; mm.is_rem = is_rem;

              ++g_ne[c];
              if ( f_flat ) ++g_cnt_flat[c];
              if ( f_clip ) ++g_cnt_clip[c];
              if ( f_rms  ) ++g_cnt_hi_rms[c];
              if ( f_ln   ) ++g_cnt_hi_ln[c];
              if ( has_spectral ) cand_log_broad[c].push_back( log_broad );
            }
        }
    }

  // per-candidate log-broadband-power calibration (median/MAD), so the
  // composite score compares candidates on a common, scale-free basis
  std::vector<double> log_broad_med( nc , 0.0 ) , log_broad_sigma( nc , 1.0 );
  for ( int c = 0 ; c < nc ; ++c )
    {
      if ( cand_log_broad[c].empty() ) continue;
      const double med = MiscMath::median( cand_log_broad[c] );
      double sigma = 1.4826 * median_abs_dev( cand_log_broad[c] , med );
      if ( sigma <= 0 ) sigma = 1.0;
      log_broad_med[c] = med; log_broad_sigma[c] = sigma;
    }

  //
  // -----------------------------------------------------------------------
  // Phase B: score (using the calibrated z-scale), select, smooth, and
  // stitch the output waveform, segment by segment
  // -----------------------------------------------------------------------
  //

  annot_t * a_src = edf.annotations->add( annot_name );

  std::vector<double> out( N , 0.0 );

  for ( size_t si = 0 ; si < segs.size() ; ++si )
    {
      const seg_work_t & sw = work[si];
      const seg_t & seg = sw.seg;
      const std::vector< std::pair<int,int> > & windows = sw.windows;
      const int nw = (int)windows.size();
      if ( nw == 0 ) continue;

      // composite, calibrated score per candidate/window
      std::vector< std::vector<double> > score( nc , std::vector<double>( nw ) );
      for ( int c = 0 ; c < nc ; ++c )
        for ( int w = 0 ; w < nw ; ++w )
          {
            const win_metrics_t & mm = sw.M[c][w];
            // clip to +/-5: an unbounded near-silent window (p_broad -> 0) can
            // otherwise produce an arbitrarily large negative log/z value that
            // would dominate both MEAN_SCORE and the switch decision
            double z_broad = mm.has_spectral ? ( mm.log_broad - log_broad_med[c] ) / log_broad_sigma[c] : 0.0;
            z_broad = std::max( -5.0 , std::min( 5.0 , z_broad ) );
            double s = z_broad - mm.line_ratio - mm.drift_ratio + 0.25 * mm.kurt_bonus;
            if ( mm.is_rem ) s += stage_weight * mm.kurt_bonus;
            s -= flat_penalty * (double)mm.flat;
            s -= clip_penalty * (double)mm.clip;
            s -= rms_penalty  * (double)mm.hi_rms;
            s -= ln_penalty   * (double)mm.hi_ln;
            score[c][w] = s;
            g_sum_score[c] += s;
          }

      //
      // switch selection: hysteresis + minimum dwell (no hard exclusion --
      // the best-scoring candidate is always chosen, so there is always an
      // active candidate and never a gap in the output)
      //

      std::vector<int> active( nw , -1 );
      int cur = -1;
      int windows_since_switch = min_dwell_windows + 1;

      for ( int w = 0 ; w < nw ; ++w )
        {
          int best = 0; double best_score = score[0][w];
          for ( int c = 1 ; c < nc ; ++c )
            if ( score[c][w] > best_score ) { best = c; best_score = score[c][w]; }

          if ( cur == -1 )
            {
              cur = best;
              windows_since_switch = 0;
            }
          else if ( cur != best )
            {
              const bool dwell_ok  = windows_since_switch >= min_dwell_windows;
              const bool margin_ok = ( best_score - score[cur][w] ) > switch_margin;

              if ( dwell_ok && margin_ok )
                {
                  cur = best;
                  windows_since_switch = 0;
                }
              else
                ++windows_since_switch;
            }
          else
            ++windows_since_switch;

          active[w] = cur;
        }

      // smoothing: short majority-vote filter to remove residual single-window flicker
      if ( smooth_win > 1 )
        {
          std::vector<int> smoothed = active;
          const int half = smooth_win / 2;
          for ( int w = 0 ; w < nw ; ++w )
            {
              std::map<int,int> counts;
              for ( int k = std::max(0,w-half) ; k <= std::min(nw-1,w+half) ; ++k )
                counts[ active[k] ]++;
              int best_lab = active[w] , best_cnt = -1;
              for ( std::map<int,int>::const_iterator it = counts.begin() ; it != counts.end() ; ++it )
                if ( it->second > best_cnt ) { best_cnt = it->second; best_lab = it->first; }
              smoothed[w] = best_lab;
            }
          active = smoothed;
        }

      // tallies
      for ( int w = 0 ; w < nw ; ++w )
        {
          ++g_cnt_used[ active[w] ];
          if ( w > 0 && active[w] != active[w-1] ) ++g_n_switch;
        }

      //
      // waveform construction
      //

      std::vector<double> win_med( nw , 0.0 ) , win_sigma( nw , 1.0 );

      for ( int w = 0 ; w < nw ; ++w )
        {
          const int a = windows[w].first , b = windows[w].second;
          const int c = active[w];

          if ( do_norm )
            {
              const int center = (a+b)/2;
              int na = std::max( seg.start_idx , center - norm_win_samples/2 );
              int nb = std::min( seg.start_idx + seg.n_samples , center + norm_win_samples/2 );
              if ( nb <= na ) { na = a; nb = b; }

              std::vector<double> ref( cand_data[c].begin() + na , cand_data[c].begin() + nb );
              const double med = MiscMath::median( ref );
              double sigma = 1.4826 * median_abs_dev( ref , med );
              if ( sigma <= 0 ) sigma = 1.0;

              win_med[w] = med; win_sigma[w] = sigma;

              for ( int k = a ; k < b ; ++k )
                out[k] = ( cand_data[c][k] - med ) / sigma;
            }
          else
            {
              for ( int k = a ; k < b ; ++k )
                out[k] = cand_data[c][k];
            }
        }

      // crossfade at genuine switch points (post-smoothing)
      if ( xfade_samples > 0 )
        {
          for ( int w = 1 ; w < nw ; ++w )
            {
              if ( active[w] == active[w-1] ) continue;

              const int boundary = windows[w].first;
              const int half = xfade_samples / 2;
              const int lo = std::max( seg.start_idx , boundary - half );
              const int hi = std::min( seg.start_idx + seg.n_samples , boundary + half );
              if ( hi <= lo ) continue;

              const int c_out = active[w-1] , c_in = active[w];

              for ( int k = lo ; k < hi ; ++k )
                {
                  const double frac = (double)( k - lo ) / (double)( hi - lo );
                  const double wgt_in = 0.5 - 0.5 * std::cos( M_PI * frac );

                  double z_out , z_in;
                  if ( do_norm )
                    {
                      z_out = ( cand_data[c_out][k] - win_med[w-1] ) / win_sigma[w-1];
                      z_in  = ( cand_data[c_in][k]  - win_med[w]   ) / win_sigma[w];
                    }
                  else
                    {
                      z_out = cand_data[c_out][k];
                      z_in  = cand_data[c_in][k];
                    }
                  out[k] = ( 1.0 - wgt_in ) * z_out + wgt_in * z_in;
                }
            }
        }

      // clamp outliers (only meaningful in normalized/z-unit space)
      if ( do_norm && clamp_z > 0 )
        for ( int k = seg.start_idx ; k < seg.start_idx + seg.n_samples ; ++k )
          out[k] = std::max( -clamp_z , std::min( clamp_z , out[k] ) );

      //
      // provenance annotation: contiguous runs of the active candidate
      //

      {
        int run_start = 0;
        for ( int w = 1 ; w <= nw ; ++w )
          {
            if ( w == nw || active[w] != active[run_start] )
              {
                const int local_a = windows[run_start].first - seg.start_idx;
                const int local_b = windows[w-1].second - seg.start_idx;
                const uint64_t t0 = seg.iv.start + (uint64_t)std::llround( (double)local_a / fs * globals::tp_1sec );
                const uint64_t t1 = seg.iv.start + (uint64_t)std::llround( (double)local_b / fs * globals::tp_1sec );
                a_src->add( cand_lab[ active[run_start] ] , interval_t( t0 , t1 ) , "." );
                run_start = w;
              }
          }
      }

      // next segment
    }

  // rescale normalized output back to approximate physical (uV) units
  if ( do_norm )
    for ( int k = 0 ; k < N ; ++k )
      out[k] *= global_scale;

  //
  // -----------------------------------------------------------------------
  // Write the new combined channel
  // -----------------------------------------------------------------------
  //

  if ( ! edf.init_signal( new_label , fs ) )
    Helper::halt( "could not create new channel '" + new_label + "' (label already exists?)" );
  const int new_slot = edf.header.signal( new_label );
  edf.update_signal( new_slot , &out );

  //
  // -----------------------------------------------------------------------
  // Output
  // -----------------------------------------------------------------------
  //

  const double total_hours = ( N / fs ) / 3600.0;

  for ( int c = 0 ; c < nc ; ++c )
    {
      writer.level( cand_lab[c] , "CH" );
      writer.value( "NE"       , g_ne[c] );
      writer.value( "PCT_USED" , g_ne[c] > 0 ? 100.0 * g_cnt_used[c] / (double)g_ne[c] : 0.0 );
      writer.value( "PCT_FLAT" , g_ne[c] > 0 ? 100.0 * g_cnt_flat[c] / (double)g_ne[c] : 0.0 );
      writer.value( "PCT_CLIP" , g_ne[c] > 0 ? 100.0 * g_cnt_clip[c] / (double)g_ne[c] : 0.0 );
      writer.value( "PCT_HI_RMS" , g_ne[c] > 0 ? 100.0 * g_cnt_hi_rms[c] / (double)g_ne[c] : 0.0 );
      writer.value( "PCT_HI_LN"  , g_ne[c] > 0 ? 100.0 * g_cnt_hi_ln[c] / (double)g_ne[c] : 0.0 );
      if ( g_ne[c] > 0 ) writer.value( "MEAN_SCORE" , g_sum_score[c] / (double)g_ne[c] );
    }
  writer.unlevel( "CH" );

  writer.value( "N_CAND"     , nc );
  writer.value( "N_SWITCH"   , g_n_switch );
  writer.value( "SWITCH_RATE", total_hours > 0 ? g_n_switch / total_hours : 0.0 );
  writer.value( "NEW"        , new_label );
  writer.value( "NEW_SR"     , fs );

  logger << "  complete: " << g_n_switch << " switches, new channel '" << new_label << "'\n\n";

}
