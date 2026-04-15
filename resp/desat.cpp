
//    --------------------------------------------------------------------
//
//    This file is part of Luna.
//
//    LUNA is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    LUNA is distributed in the hope that it will be useful,
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
// DESAT : Oxygen desaturation detector
//
// Designed to be robust to noisy SpO2 signals:
//   - heavy quantization (e.g. integer steps)
//   - ~1 Hz SR (variable)
//   - dropout, spikes, drift
//   - out-of-range values
//
// Algorithm overview:
//   1. Flag out-of-range and spike samples as artifact
//   2. For each contiguous segment, scan forward with a rolling
//      120-second median baseline (excluding artifact and desat samples)
//   3. A desat starts when signal drops >= drop_threshold below baseline
//   4. A desat ends when signal recovers to >= baseline - drop_threshold * recovery_frac
//   5. Events lasting >= min_dur seconds are recorded and annotated
//

#include "resp/desat.h"

#include "param.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "db/db.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "annot/annot.h"
#include "defs/defs.h"
#include "dsp/iir.h"

#include <algorithm>
#include <cmath>
#include <deque>
#include <memory>
#include <numeric>
#include <sstream>
#include <utility>

extern writer_t writer;
extern logger_t logger;


// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

static double median_of( std::vector<double> & v )
{
  // compute median in-place (sorts v)
  if ( v.empty() ) return 0.0;
  const int n = v.size();
  std::sort( v.begin(), v.end() );
  if ( n % 2 == 1 ) return v[ n / 2 ];
  return 0.5 * ( v[ n/2 - 1 ] + v[ n/2 ] );
}


// ---------------------------------------------------------------------------
// Main constructor — processes all segments for one signal
// ---------------------------------------------------------------------------

desat_t::desat_t( edf_t & edf , param_t & param )
{

  // Dispatch to Matlab-style peak-valley mode if requested
  if ( param.has( "mode" ) && param.value( "mode" ) == "matlab" )
    {
      desat_matlab_t dm( edf , param );
      return;
    }

  //
  // -----------------------------------------------------------------------
  // Parameters
  // -----------------------------------------------------------------------
  //

  // Required: signal label
  if ( ! param.has( "sig" ) )
    Helper::halt( "DESAT requires sig=<label>" );

  const std::string sig_label = param.value( "sig" );

  // Artifact: values outside [low_art, 100] are hard artifacts
  const double low_art     = param.has( "low"      ) ? param.requires_dbl( "low"      ) : 50.0;

  // Spike detection: flag sample if |x[i] - x[i-1]| > spike_th (SpO2 units, 0-100 scale)
  const double spike_th    = param.has( "spike"    ) ? param.requires_dbl( "spike"    ) : 10.0;

  // Desaturation threshold: signal must drop >= drop_th below baseline
  const double drop_th     = param.has( "drop"     ) ? param.requires_dbl( "drop"     ) : 3.0;

  // Minimum duration (seconds) for an event to be recorded
  const double min_dur     = param.has( "dur"      ) ? param.requires_dbl( "dur"      ) : 10.0;

  // Baseline window (seconds of prior valid signal used for rolling median)
  const double baseline_sec = param.has( "baseline" ) ? param.requires_dbl( "baseline" ) : 120.0;

  // Minimum number of valid samples required in baseline window before detecting desats
  const int min_bsln       = param.has( "min-bsln" ) ? param.requires_int( "min-bsln" ) : 30;

  // Recovery: desat ends when signal >= baseline_at_onset - drop_th * recovery_frac
  const double recovery    = param.has( "recovery" ) ? param.requires_dbl( "recovery" ) : 0.5;

  // Annotation names
  const std::string desat_annot_name = param.has( "annot"     ) ? param.value( "annot"     ) : "desat";
  const std::string art_annot_name   = param.has( "art-annot" ) ? param.value( "art-annot" ) : "spo2-artifact";

  // SpO2 thresholds for time-below metrics
  // (allow user to specify custom set; default: 90, 88, 85, 80)
  std::vector<double> pct_thresholds = { 90.0, 88.0, 85.0, 80.0 };
  if ( param.has( "pct-th" ) )
    pct_thresholds = param.dblvector( "pct-th" );

  //
  // -----------------------------------------------------------------------
  // Log parameters clearly
  // -----------------------------------------------------------------------
  //

  logger << "\n  DESAT : oxygen desaturation detector\n";
  logger << "    signal             : " << sig_label << "\n";
  logger << "    artifact threshold : SpO2 < " << low_art << "% or > 100%\n";
  logger << "    spike threshold    : |delta| > " << spike_th << " SpO2 units\n";
  logger << "    drop threshold     : " << drop_th << "% below baseline\n";
  logger << "    minimum duration   : " << min_dur << " seconds\n";
  logger << "    baseline window    : " << baseline_sec << " seconds\n";
  logger << "    min baseline N     : " << min_bsln << " valid samples\n";
  logger << "    recovery fraction  : " << recovery << " (recover to baseline - " << drop_th * recovery << "%)\n";
  logger << "    desat annotation   : " << desat_annot_name << "\n";
  logger << "    artifact annotation: " << art_annot_name << "\n";

  //
  // -----------------------------------------------------------------------
  // Get signal
  // -----------------------------------------------------------------------
  //

  int sig_n = edf.header.signal( sig_label );
  if ( sig_n == -1 )
    {
      logger << "  ** could not find signal '" << sig_label << "', skipping DESAT\n";
      return;
    }

  const int fs = edf.header.sampling_freq( sig_n );
  if ( fs <= 0 )
    {
      logger << "  ** signal has zero sample rate, skipping DESAT\n";
      return;
    }

  logger << "    sample rate        : " << fs << " Hz\n";

  //
  // -----------------------------------------------------------------------
  // Get contiguous segments
  // -----------------------------------------------------------------------
  //

  std::set<interval_t> segs = edf.timeline.segments();
  if ( segs.empty() )
    {
      logger << "  ** no segments found, skipping DESAT\n";
      return;
    }

  logger << "    segments           : " << segs.size() << "\n\n";

  //
  // -----------------------------------------------------------------------
  // Annotation objects (created once, populated per-segment)
  // -----------------------------------------------------------------------
  //

  annot_t * desat_annot = edf.annotations->add( desat_annot_name );
  annot_t * art_annot   = edf.annotations->add( art_annot_name );

  //
  // -----------------------------------------------------------------------
  // Per-segment processing — accumulate global stats
  // -----------------------------------------------------------------------
  //

  int    total_desats    = 0;
  double total_t_valid   = 0.0;
  double total_t_art     = 0.0;
  double total_t_desat   = 0.0;
  double sum_spo2        = 0.0;
  long   n_spo2          = 0;
  double sum_nadir       = 0.0;
  double sum_drop        = 0.0;
  double sum_dur         = 0.0;
  double min_nadir       = 101.0;
  std::vector<double> pct_below_counts( pct_thresholds.size(), 0.0 );
  long   valid_count     = 0;
  int    event_counter   = 0;

  // iterate segments
  int seg_idx = 0;
  std::set<interval_t>::const_iterator si = segs.begin();
  while ( si != segs.end() )
    {
      const interval_t & seg = *si;

      //
      // Pull signal for this segment
      //

      slice_t slice( edf, sig_n, seg );
      const std::vector<double> * pdata = slice.pdata();
      if ( pdata == NULL || pdata->empty() ) { ++si; ++seg_idx; continue; }

      const int N = (int)pdata->size();
      const std::vector<double> & x = *pdata;

      // baseline_n: number of samples in the rolling baseline window
      const int baseline_n = (int)( baseline_sec * fs );
      // min_dur_samps: minimum samples for a qualifying event
      const int min_dur_samps = (int)( min_dur * fs );

      //
      // Step 1: artifact flags
      //
      //  hard_artifact : out-of-range (excluded from everything)
      //  spike         : sudden jump (excluded from baseline, treated as non-signal)
      //

      std::vector<bool> hard_art( N, false );
      std::vector<bool> spike_art( N, false );

      // out-of-range
      for ( int i = 0; i < N; i++ )
        if ( x[i] < low_art || x[i] > 100.0 )
          hard_art[i] = true;

      // spikes: compare consecutive valid (non-hard-artifact) samples;
      // if the jump exceeds spike_th, mark both the "before" and "after" if they
      // form an isolated excursion; otherwise mark the higher-change sample
      for ( int i = 1; i < N; i++ )
        {
          if ( hard_art[i] || hard_art[i-1] ) continue;
          if ( std::fabs( x[i] - x[i-1] ) > spike_th )
            {
              // mark the sample further from its outer neighbour as the spike
              // simple heuristic: if next sample also looks OK, mark current as spike
              if ( i + 1 < N && !hard_art[i+1] &&
                   std::fabs( x[i] - x[i+1] ) > spike_th )
                spike_art[i] = true;   // isolated spike (x[i] differs from both neighbours)
              else if ( i >= 2 && !hard_art[i-2] &&
                        std::fabs( x[i-1] - x[i-2] ) > spike_th )
                spike_art[i-1] = true; // x[i-1] differs from both neighbours
              else
                spike_art[i] = true;   // just mark the later one
            }
        }

      // combined bad-sample flag
      std::vector<bool> bad( N, false );
      for ( int i = 0; i < N; i++ )
        bad[i] = hard_art[i] || spike_art[i];

      //
      // Step 2: emit artifact annotations for contiguous bad runs
      //
      {
        bool in_art = false;
        int  art_start = 0;
        for ( int i = 0; i <= N; i++ )
          {
            bool b = (i < N) ? bad[i] : false;
            if ( !in_art && b ) { in_art = true; art_start = i; }
            else if ( in_art && !b )
              {
                // convert sample indices to absolute time points
                uint64_t tp_start = seg.start
                  + (uint64_t)round( (double)art_start / fs * globals::tp_1sec );
                uint64_t tp_stop  = seg.start
                  + (uint64_t)round( (double)i         / fs * globals::tp_1sec );
                art_annot->add( ".", interval_t( tp_start, tp_stop ), "." );
                in_art = false;
              }
          }
      }

      //
      // Step 3: desat detection (forward scan with rolling median baseline)
      //
      //  window_vals  : circular buffer of recent valid, non-desat sample values
      //                 (aligned to window of last baseline_n samples, but only
      //                  valid samples actually stored)
      //
      //  We store (index, value) pairs so we can evict samples that have aged out.
      //

      std::deque< std::pair<int,double> > bsln_buf;  // (sample_index, value) of valid baseline samples
      std::vector<bool> in_desat_flag( N, false );   // mark samples as being in (tentative) desat

      bool   in_desat       = false;
      int    desat_start_i  = -1;
      double desat_baseline = 0.0;
      double desat_nadir    = 100.0;

      std::vector<desat_event_t> seg_events;

      for ( int i = 0; i < N; i++ )
        {
          // -- evict samples that have aged out of the baseline window
          while ( !bsln_buf.empty() && ( i - bsln_buf.front().first ) > baseline_n )
            bsln_buf.pop_front();

          // -- compute baseline from buffer (copy values, sort for median)
          bool   bsln_valid = false;
          double baseline   = 0.0;
          if ( (int)bsln_buf.size() >= min_bsln )
            {
              std::vector<double> tmp;
              tmp.reserve( bsln_buf.size() );
              for ( const auto & p : bsln_buf )
                tmp.push_back( p.second );
              baseline   = median_of( tmp );
              bsln_valid = true;
            }

          // -- current sample handling
          if ( bad[i] )
            {
              // bad sample: don't update desat state, don't push to baseline
              // but if we're in a desat, the desat "pauses" — we neither
              // end it nor update the nadir; bad samples are not annotated as desat
              // (they already get the artifact annotation)
            }
          else
            {
              // valid sample

              if ( !in_desat )
                {
                  // check for desat onset
                  if ( bsln_valid && x[i] < baseline - drop_th )
                    {
                      in_desat       = true;
                      desat_start_i  = i;
                      desat_baseline = baseline;
                      desat_nadir    = x[i];
                      in_desat_flag[i] = true;
                    }
                  else
                    {
                      // normal: push to baseline buffer
                      bsln_buf.push_back( { i, x[i] } );
                    }
                }
              else
                {
                  // currently in desat
                  in_desat_flag[i] = true;
                  if ( x[i] < desat_nadir ) desat_nadir = x[i];

                  // recovery criterion: signal >= baseline_at_onset - drop_th * recovery_frac
                  const double recovery_line = desat_baseline - drop_th * recovery;
                  if ( x[i] >= recovery_line )
                    {
                      // desat ended at this sample
                      int desat_stop_i = i;
                      int dur_samps    = desat_stop_i - desat_start_i;

                      if ( dur_samps >= min_dur_samps )
                        {
                          desat_event_t ev;
                          ev.start_sec = (double)seg.start / globals::tp_1sec
                                         + (double)desat_start_i / fs;
                          ev.stop_sec  = (double)seg.start / globals::tp_1sec
                                         + (double)desat_stop_i  / fs;
                          ev.duration  = (double)dur_samps / fs;
                          ev.nadir     = desat_nadir;
                          ev.baseline  = desat_baseline;
                          ev.drop      = desat_baseline - desat_nadir;
                          ev.seg       = seg_idx;
                          seg_events.push_back( ev );

                          // create annotation
                          uint64_t tp_start = seg.start
                            + (uint64_t)round( (double)desat_start_i / fs * globals::tp_1sec );
                          uint64_t tp_stop  = seg.start
                            + (uint64_t)round( (double)desat_stop_i  / fs * globals::tp_1sec );

                          instance_t * inst = desat_annot->add(
                            ".", interval_t( tp_start, tp_stop ), "." );
                          inst->set( "nadir",    desat_nadir );
                          inst->set( "baseline", desat_baseline );
                          inst->set( "drop",     ev.drop );
                          inst->set( "dur",      ev.duration );
                        }

                      in_desat   = false;
                      desat_nadir = 100.0;
                      // push the recovery sample to baseline
                      bsln_buf.push_back( { i, x[i] } );
                    }
                  // if still in desat: do NOT push to baseline
                }
            }
        } // end per-sample loop

      // handle desat that extends to end of segment
      if ( in_desat )
        {
          int desat_stop_i = N;
          int dur_samps    = desat_stop_i - desat_start_i;
          if ( dur_samps >= min_dur_samps )
            {
              desat_event_t ev;
              ev.start_sec = (double)seg.start / globals::tp_1sec
                             + (double)desat_start_i / fs;
              ev.stop_sec  = (double)seg.start / globals::tp_1sec
                             + (double)desat_stop_i  / fs;
              ev.duration  = (double)dur_samps / fs;
              ev.nadir     = desat_nadir;
              ev.baseline  = desat_baseline;
              ev.drop      = desat_baseline - desat_nadir;
              ev.seg       = seg_idx;
              seg_events.push_back( ev );

              uint64_t tp_start = seg.start
                + (uint64_t)round( (double)desat_start_i / fs * globals::tp_1sec );
              uint64_t tp_stop  = seg.stop;  // use segment end

              instance_t * inst = desat_annot->add(
                ".", interval_t( tp_start, tp_stop ), "." );
              inst->set( "nadir",    desat_nadir );
              inst->set( "baseline", desat_baseline );
              inst->set( "drop",     ev.drop );
              inst->set( "dur",      ev.duration );
            }
        }

      //
      // Step 4: per-segment statistics
      //

      double seg_t_valid  = 0.0;
      double seg_t_art    = 0.0;
      double seg_t_desat  = 0.0;
      double seg_sum_spo2 = 0.0;
      long   seg_n_spo2   = 0;
      std::vector<double> seg_pct_below( pct_thresholds.size(), 0.0 );
      long seg_valid_n = 0;

      for ( int i = 0; i < N; i++ )
        {
          double dt = 1.0 / fs;
          if ( bad[i] )
            {
              seg_t_art += dt;
            }
          else
            {
              seg_t_valid += dt;
              seg_sum_spo2 += x[i];
              seg_n_spo2++;

              for ( int ti = 0; ti < (int)pct_thresholds.size(); ti++ )
                if ( x[i] < pct_thresholds[ti] )
                  seg_pct_below[ti] += dt;

              seg_valid_n++;

              if ( in_desat_flag[i] )
                seg_t_desat += dt;
            }
        }

      // accumulate into totals
      total_desats  += (int)seg_events.size();
      total_t_valid += seg_t_valid;
      total_t_art   += seg_t_art;
      total_t_desat += seg_t_desat;
      sum_spo2      += seg_sum_spo2;
      n_spo2        += seg_n_spo2;
      valid_count   += seg_valid_n;

      for ( int ti = 0; ti < (int)pct_thresholds.size(); ti++ )
        pct_below_counts[ti] += seg_pct_below[ti];

      for ( const auto & ev : seg_events )
        {
          sum_nadir += ev.nadir;
          sum_drop  += ev.drop;
          sum_dur   += ev.duration;
          if ( ev.nadir < min_nadir ) min_nadir = ev.nadir;
        }

      //
      // Step 5: per-event writer output (strat by event counter)
      //

      for ( const auto & ev : seg_events )
        {
          ++event_counter;
          writer.level( event_counter , "DESAT" );
          writer.value( "START"    , ev.start_sec );
          writer.value( "STOP"     , ev.stop_sec );
          writer.value( "DUR"      , ev.duration );
          writer.value( "NADIR"    , ev.nadir );
          writer.value( "BASELINE" , ev.baseline );
          writer.value( "DROP"     , ev.drop );
          writer.value( "SEG"      , ev.seg );
          writer.unlevel( "DESAT" );
        }

      ++seg_idx;
      ++si;
    } // end segment loop


  //
  // -----------------------------------------------------------------------
  // Individual-level output
  // -----------------------------------------------------------------------
  //

  const double mean_spo2  = ( n_spo2 > 0 ) ? ( sum_spo2 / n_spo2 ) : 0.0;
  const double mean_nadir = ( total_desats > 0 ) ? ( sum_nadir / total_desats ) : 0.0;
  const double mean_drop  = ( total_desats > 0 ) ? ( sum_drop  / total_desats ) : 0.0;
  const double mean_dur   = ( total_desats > 0 ) ? ( sum_dur   / total_desats ) : 0.0;

  writer.value( "N"        , total_desats );
  writer.value( "T_VALID"  , total_t_valid );
  writer.value( "T_ART"    , total_t_art );
  writer.value( "T_DESAT"  , total_t_desat );
  writer.value( "MEAN_SPO2", mean_spo2 );
  writer.value( "NADIR_MEAN", mean_nadir );
  writer.value( "NADIR_MIN" , ( total_desats > 0 ) ? min_nadir : 0.0 );
  writer.value( "DROP_MEAN" , mean_drop );
  writer.value( "DUR_MEAN"  , mean_dur );
  writer.value( "DUR_TOTAL" , sum_dur );

  // time-below-threshold as % of valid time
  if ( total_t_valid > 0 )
    for ( int ti = 0; ti < (int)pct_thresholds.size(); ti++ )
      {
        std::ostringstream oss;
        oss << "PCT_LT" << (int)pct_thresholds[ti];
        writer.value( oss.str(), 100.0 * pct_below_counts[ti] / total_t_valid );
      }
  else
    for ( int ti = 0; ti < (int)pct_thresholds.size(); ti++ )
      {
        std::ostringstream oss;
        oss << "PCT_LT" << (int)pct_thresholds[ti];
        writer.value( oss.str(), 0.0 );
      }

  logger << "  DESAT complete: " << total_desats << " events detected"
         << "  (valid time " << total_t_valid << "s"
         << ", artifact " << total_t_art << "s)\n\n";
}


// ---------------------------------------------------------------------------
//
// desat_matlab_t : Matlab-style peak-valley oxygen desaturation detector
//
// Implements the CalcODI / SpO2ArtifactReject approach:
//   1. Artifact: out-of-range + ±neighbor-second expansion
//   2. Nearest-neighbour gap-fill (for peak detection only)
//   3. peakdet (Billauer, delta=0.5) on inverted signal
//   4. Iterative peak-valley skeleton pruning (MagAv_thres)
//   5. Local pre/post peak search within ±baseline-sec
//   6. ODI2/3/4 output; optional sleep-stage gating via N1/N2/N3/R annotations
//
// ---------------------------------------------------------------------------


// Internal: peakdet result
struct desat_pk_t {
  std::vector<int> maxX;  // sample indices of local maxima (in original, non-flipped signal)
  std::vector<int> minX;  // sample indices of local minima
};


// Billauer peakdet on a std::vector<double>.
// flip=true inverts the signal so that signal valleys become detected maxima
// (and the roles of maxX/minX are swapped back by the caller).
static desat_pk_t desat_peakdet( const std::vector<double> & v ,
                                  double delta ,
                                  bool flip )
{
  desat_pk_t r;
  const int n = (int)v.size();
  if ( n == 0 || delta <= 0 ) return r;

  const double sgn = flip ? -1.0 : 1.0;

  double mn =  1e300, mx = -1e300;
  int    mnpos = 0,   mxpos = 0;
  bool   lookformax = true;

  for ( int i = 0; i < n; i++ )
    {
      double th = sgn * v[i];
      if ( th > mx ) { mx = th; mxpos = i; }
      if ( th < mn ) { mn = th; mnpos = i; }

      if ( lookformax )
        {
          if ( th < mx - delta )
            {
              r.maxX.push_back( mxpos );
              mn = th; mnpos = i;
              lookformax = false;
            }
        }
      else
        {
          if ( th > mn + delta )
            {
              r.minX.push_back( mnpos );
              mx = th; mxpos = i;
              lookformax = true;
            }
        }
    }
  return r;
}


// Check whether timepoint tp falls inside any interval in a sorted list of
// [start, stop) pairs (sorted by start).
static bool tp_in_sleep( const std::vector< std::pair<uint64_t,uint64_t> > & ivs ,
                          uint64_t tp )
{
  if ( ivs.empty() ) return false;
  // upper_bound on start
  auto it = std::upper_bound( ivs.begin(), ivs.end(),
                               std::make_pair( tp , (uint64_t)UINT64_MAX ) );
  if ( it != ivs.begin() )
    {
      --it;
      if ( tp >= it->first && tp < it->second ) return true;
    }
  return false;
}


// ---------------------------------------------------------------------------
// Main constructor
// ---------------------------------------------------------------------------

desat_matlab_t::desat_matlab_t( edf_t & edf , param_t & param )
{

  //
  // -----------------------------------------------------------------------
  // Parameters
  // -----------------------------------------------------------------------
  //

  if ( ! param.has( "sig" ) )
    Helper::halt( "DESAT mode=matlab requires sig=<label>" );

  const std::string sig_label = param.value( "sig" );

  // Hard artifact floor — Matlab uses 40; values > 101 are artifact, 100-101 clamped to 100
  const double low_art       = param.has( "low"       ) ? param.requires_dbl( "low"       ) : 40.0;

  // Neighbourhood expansion around artifact samples (seconds) — Matlab: SpO2neighborhoodD_dt=2
  const double neighbor_sec  = param.has( "neighbor"  ) ? param.requires_dbl( "neighbor"  ) : 2.0;

  // Butterworth low-pass filter cutoff (Hz) — Matlab default: 1 Hz, order 2; set lp=0 to disable
  const double lp_cutoff     = param.has( "lp"        ) ? param.requires_dbl( "lp"        ) : 1.0;
  const int    lp_order      = param.has( "lp-order"  ) ? param.requires_int( "lp-order"  ) : 2;

  // Round signal to nearest integer after LP filter — Matlab always does this; use no-round to skip
  const bool   do_round      = ! param.has( "no-round" );

  // Peak-valley pruning threshold (MagAv_thres in Matlab, units: % SpO2)
  const double mag_thres     = param.has( "mag-thres" ) ? param.requires_dbl( "mag-thres" ) : 1.5;

  // Local pre/post baseline search window (seconds each side)
  const double baseline_sec  = param.has( "baseline"  ) ? param.requires_dbl( "baseline"  ) : 120.0;

  // Annotation names
  const std::string desat_annot_name = param.has( "annot"     ) ? param.value( "annot"     ) : "desat";
  const std::string art_annot_name   = param.has( "art-annot" ) ? param.value( "art-annot" ) : "spo2-artifact";

  logger << "\n  DESAT (matlab mode) : oxygen desaturation detector\n";
  logger << "    signal             : " << sig_label << "\n";
  logger << "    artifact threshold : SpO2 < " << low_art << "% or > 101% (100-101 clamped)\n";
  logger << "    neighbor expansion : ±" << neighbor_sec << " s\n";
  if ( lp_cutoff > 0 )
    logger << "    LP filter          : order=" << lp_order << "  cutoff=" << lp_cutoff << " Hz"
           << ( do_round ? "  +round" : "" ) << "\n";
  else
    logger << "    LP filter          : disabled\n";
  logger << "    mag-thres (pruning): " << mag_thres << "% SpO2\n";
  logger << "    baseline window    : ±" << baseline_sec << " s\n";
  logger << "    desat annotation   : " << desat_annot_name << "\n";
  logger << "    artifact annotation: " << art_annot_name << "\n";

  //
  // -----------------------------------------------------------------------
  // Get signal
  // -----------------------------------------------------------------------
  //

  int sig_n = edf.header.signal( sig_label );
  if ( sig_n == -1 )
    {
      logger << "  ** could not find signal '" << sig_label << "', skipping DESAT\n";
      return;
    }

  const int fs = edf.header.sampling_freq( sig_n );
  if ( fs <= 0 )
    {
      logger << "  ** signal has zero sample rate, skipping DESAT\n";
      return;
    }

  logger << "    sample rate        : " << fs << " Hz\n";

  const uint64_t fs_tp = globals::tp_1sec / (uint64_t)fs;  // timepoints per sample

  //
  // -----------------------------------------------------------------------
  // Optional LP filter (applied once to full signal via segments)
  // We apply per-segment to avoid filter edge effects crossing gaps.
  // -----------------------------------------------------------------------
  //

  std::unique_ptr<iir_t> lp_filt;
  if ( lp_cutoff > 0.0 )
    {
      lp_filt.reset( new iir_t );
      lp_filt->init( BUTTERWORTH_LOWPASS , lp_order ,
                     (double)fs , lp_cutoff , 0.0 );
    }

  //
  // -----------------------------------------------------------------------
  // Build sleep-stage intervals (N1, N2, N3, R) for ODI gating
  // If no staging annotations are found, ODI denominator = all valid time.
  // -----------------------------------------------------------------------
  //

  std::vector< std::pair<uint64_t,uint64_t> > sleep_ivs;  // sorted by start
  bool has_staging = false;

  {
    const std::vector<std::string> stage_labels = { "N1","N2","N3","R" };
    for ( const auto & lbl : stage_labels )
      {
        annot_t * a = edf.annotations->find( lbl );
        if ( a == NULL ) continue;
        has_staging = true;
        annot_map_t evts = a->extract( edf.timeline.wholetrace() );
        for ( auto & kv : evts )
          sleep_ivs.push_back( { kv.first.interval.start ,
                                  kv.first.interval.stop  } );
      }
    std::sort( sleep_ivs.begin(), sleep_ivs.end() );
  }

  if ( has_staging )
    logger << "    staging            : found (ODI gated to sleep)\n";
  else
    logger << "    staging            : not found (ODI over all valid time)\n";

  //
  // -----------------------------------------------------------------------
  // Contiguous segments
  // -----------------------------------------------------------------------
  //

  std::set<interval_t> segs = edf.timeline.segments();
  if ( segs.empty() )
    {
      logger << "  ** no segments found, skipping DESAT\n";
      return;
    }

  logger << "    segments           : " << segs.size() << "\n\n";

  //
  // -----------------------------------------------------------------------
  // Annotation objects
  // -----------------------------------------------------------------------
  //

  annot_t * desat_annot = edf.annotations->add( desat_annot_name );
  annot_t * art_annot   = edf.annotations->add( art_annot_name );

  //
  // -----------------------------------------------------------------------
  // Global accumulators
  // -----------------------------------------------------------------------
  //

  double total_t_valid   = 0.0;
  double total_t_art     = 0.0;
  double total_t_sleep_valid = 0.0;  // valid + sleep (ODI denominator)
  double sum_spo2        = 0.0;
  long   n_spo2          = 0;

  // All qualifying events (MagDown >= 2%) across segments
  struct ml_event_t {
    double start_sec;
    double stop_sec;
    double nadir_sec;  // seconds from recording start of nadir sample
    double nadir_val;
    double pre_val;    // SpO2 at pre-peak
    double post_val;   // SpO2 at post-peak
    double mag_down;
    double mag_up;
    int    seg;
    bool   in_sleep;
  };

  std::vector<ml_event_t> all_events;

  int seg_idx = 0;
  for ( const auto & seg : segs )
    {
      //
      // Pull signal for this segment
      //

      slice_t slice( edf, sig_n, seg );
      const std::vector<double> * pdata = slice.pdata();
      if ( pdata == NULL || pdata->empty() ) { ++seg_idx; continue; }

      // work on a mutable copy
      std::vector<double> x( *pdata );
      const int N = (int)x.size();

      //
      // Step 1a: Optional LP filter
      // Guard: skip if cutoff >= Nyquist (fs/2); this happens when signal is
      // already at ~1 Hz and the default 1 Hz cutoff would be invalid.
      //

      if ( lp_filt && lp_cutoff < (double)fs / 2.0 )
        {
          // re-init filter state for each segment to avoid inter-segment leakage
          lp_filt->init( BUTTERWORTH_LOWPASS , lp_order ,
                         (double)fs , lp_cutoff , 0.0 );
          x = lp_filt->apply( x );
        }
      else if ( lp_filt )
        {
          logger << "  ** DESAT (matlab): LP cutoff " << lp_cutoff
                 << " Hz >= Nyquist (" << (double)fs/2.0
                 << " Hz) for this signal, skipping filter\n";
        }

      //
      // Step 1b: Optional rounding to integer SpO2 values
      //

      if ( do_round )
        for ( int i = 0; i < N; i++ )
          x[i] = std::round( x[i] );

      //
      // Step 2: Artifact detection
      //

      const int neighbor_samps = (int)( neighbor_sec * fs );

      // Hard artifact: out-of-range
      // Matches Matlab: clamp 100-101 to 100, flag < low_art or > 101 as artifact
      std::vector<bool> bad( N, false );
      for ( int i = 0; i < N; i++ )
        {
          if ( x[i] > 100.0 && x[i] <= 101.0 ) x[i] = 100.0;  // clamp
          if ( x[i] < low_art || x[i] > 101.0 ) bad[i] = true;
        }

      // Neighbourhood expansion: dilate bad zones by ±neighbor_samps
      if ( neighbor_samps > 0 )
        {
          std::vector<bool> bad2 = bad;
          for ( int i = 0; i < N; i++ )
            {
              if ( !bad[i] ) continue;
              const int lo = std::max( 0,   i - neighbor_samps );
              const int hi = std::min( N-1, i + neighbor_samps );
              for ( int j = lo; j <= hi; j++ ) bad2[j] = true;
            }
          bad = bad2;
        }

      //
      // Step 3: Emit artifact annotations for contiguous bad runs
      //

      {
        bool in_art = false;
        int  art_start = 0;
        for ( int i = 0; i <= N; i++ )
          {
            bool b = ( i < N ) ? bad[i] : false;
            if ( !in_art && b ) { in_art = true; art_start = i; }
            else if ( in_art && !b )
              {
                uint64_t tp_s = seg.start + (uint64_t)art_start * fs_tp;
                uint64_t tp_e = seg.start + (uint64_t)i         * fs_tp;
                art_annot->add( ".", interval_t( tp_s, tp_e ), "." );
                in_art = false;
              }
          }
      }

      //
      // Step 4: Per-segment statistics
      // Done here (before peakdet) so early-exit continues don't skip it.
      //

      {
        const double dt = 1.0 / fs;
        for ( int i = 0; i < N; i++ )
          {
            const uint64_t tp_i = seg.start + (uint64_t)i * fs_tp;
            if ( bad[i] )
              {
                total_t_art += dt;
              }
            else
              {
                total_t_valid += dt;
                sum_spo2 += x[i];
                ++n_spo2;
                bool is_sl = has_staging
                              ? tp_in_sleep( sleep_ivs , tp_i )
                              : true;
                if ( is_sl ) total_t_sleep_valid += dt;
              }
          }
      }

      //
      // Step 5: Nearest-neighbour gap-fill for peak detection
      // Bad samples are replaced with the nearest valid sample value.
      // This mirrors Matlab's fillmissing(SaO2,'nearest').
      // The bad[] mask is retained for stats.
      //

      std::vector<double> xf = x;  // filled version

      // Forward pass: fill from previous valid
      {
        double last_valid = 0.0;
        bool   have_valid = false;
        for ( int i = 0; i < N; i++ )
          {
            if ( !bad[i] ) { last_valid = x[i]; have_valid = true; }
            else if ( have_valid ) xf[i] = last_valid;
          }
      }
      // Backward pass: fill from next valid (handles leading bad runs)
      {
        double next_valid = 0.0;
        bool   have_valid = false;
        for ( int i = N-1; i >= 0; i-- )
          {
            if ( !bad[i] ) { next_valid = x[i]; have_valid = true; }
            else if ( have_valid && xf[i] == x[i] )  // still unfilled
              xf[i] = next_valid;
          }
      }

      //
      // Step 5: peakdet on gap-filled signal (inverted, delta=0.5)
      //
      // With flip=true, peakdet finds signal valleys (local minima) as maxX
      // and signal peaks (local maxima) as minX.  We swap the names back:
      //   SaO2MaxIdx <- p.minX   (local maxima of SpO2)
      //   SaO2MinIdx <- p.maxX   (local minima of SpO2 = desaturation nadirs)
      //

      desat_pk_t p = desat_peakdet( xf , 0.5 , true );

      if ( p.maxX.empty() || p.minX.empty() ) { ++seg_idx; continue; }

      // Rename to match Matlab variable names (swapped because of flip)
      std::deque<int> SaO2MaxIdx( p.minX.begin(), p.minX.end() );
      std::deque<int> SaO2MinIdx( p.maxX.begin(), p.maxX.end() );

      // Ensure sequence starts and ends with a maximum:
      //   [MAX] MIN [MAX] MIN ... [MAX]
      while ( !SaO2MaxIdx.empty() && SaO2MaxIdx.front() < 0 )
        SaO2MaxIdx.pop_front();

      while ( !SaO2MinIdx.empty() && !SaO2MaxIdx.empty() &&
              SaO2MinIdx.back() >= SaO2MaxIdx.back() )
        SaO2MinIdx.pop_back();

      while ( !SaO2MinIdx.empty() && !SaO2MaxIdx.empty() &&
              SaO2MinIdx.front() <= SaO2MaxIdx.front() )
        SaO2MinIdx.pop_front();

      if ( SaO2MaxIdx.size() < 2 || SaO2MinIdx.empty() ) { ++seg_idx; continue; }

      //
      // Step 6: Iterative peak-valley skeleton pruning
      // Merge pairs whose min(MagDown, MagUp) < mag_thres
      //

      {
        auto calc_mags = [&]( std::vector<double> & MD ,
                               std::vector<double> & MU )
          {
            int nm = (int)SaO2MinIdx.size();
            MD.assign( nm, 0.0 );
            MU.assign( nm, 0.0 );
            for ( int i = 0; i < nm; i++ )
              {
                MD[i] = xf[ SaO2MaxIdx[i]   ] - xf[ SaO2MinIdx[i] ];
                MU[i] = xf[ SaO2MaxIdx[i+1] ] - xf[ SaO2MinIdx[i] ];
              }
          };

        std::vector<double> MagDown, MagUp;
        calc_mags( MagDown, MagUp );

        while ( !SaO2MaxIdx.empty() && !SaO2MinIdx.empty() )
          {
            // Find minimum of MagDown
            int md_i = (int)( std::min_element( MagDown.begin(), MagDown.end() )
                               - MagDown.begin() );
            double md_v = MagDown[ md_i ];

            // Find minimum of MagUp
            int mu_i = (int)( std::min_element( MagUp.begin(), MagUp.end() )
                               - MagUp.begin() );
            double mu_v = MagUp[ mu_i ];

            double minMag = std::min( md_v, mu_v );
            if ( minMag > mag_thres ) break;

            if ( md_v <= mu_v )  // collapse around the smallest drop
              {
                SaO2MaxIdx.erase( SaO2MaxIdx.begin() + md_i );
                SaO2MinIdx.erase( SaO2MinIdx.begin() + md_i );
              }
            else                 // collapse around the smallest recovery
              {
                SaO2MaxIdx.erase( SaO2MaxIdx.begin() + mu_i + 1 );
                SaO2MinIdx.erase( SaO2MinIdx.begin() + mu_i );
              }

            if ( SaO2MaxIdx.size() < 2 || SaO2MinIdx.empty() ) break;
            calc_mags( MagDown, MagUp );
          }
      }

      if ( SaO2MaxIdx.size() < 2 || SaO2MinIdx.empty() ) { ++seg_idx; continue; }

      //
      // Step 7: Refine pre/post baseline to local maximum within ±baseline_sec
      //

      const int bsln_samps = (int)( baseline_sec * fs );
      const int n_ev = (int)SaO2MinIdx.size();

      // Initialise pre/post indices from adjacent max positions
      std::vector<int> SpO2prei ( n_ev );
      std::vector<int> SpO2posti( n_ev );
      for ( int i = 0; i < n_ev; i++ )
        {
          SpO2prei [i] = SaO2MaxIdx[i];
          SpO2posti[i] = SaO2MaxIdx[i+1];
        }

      for ( int i = 0; i < n_ev; i++ )
        {
          int nadir = SaO2MinIdx[i];

          // Left search: find last index of maximum within [nadir-ileft, nadir]
          int ileft = std::min( nadir - SaO2MaxIdx[i] , bsln_samps );
          int left_from = nadir - ileft;
          if ( left_from < 0 ) left_from = 0;
          double best_val = xf[ left_from ];
          int    best_idx = left_from;
          for ( int j = left_from; j <= nadir; j++ )
            if ( xf[j] >= best_val ) { best_val = xf[j]; best_idx = j; }
          SpO2prei[i] = best_idx;

          // Right search: find first index of maximum within [nadir, nadir+iright]
          int iright = std::min( SaO2MaxIdx[i+1] - nadir , bsln_samps );
          int right_to = std::min( nadir + iright , N-1 );
          best_val = xf[ nadir ];
          best_idx = nadir;
          for ( int j = nadir; j <= right_to; j++ )
            if ( xf[j] > best_val ) { best_val = xf[j]; best_idx = j; }
          SpO2posti[i] = best_idx;
        }

      //
      // Step 8: Collect events (MagDown >= 2%) and accumulate statistics
      //

      for ( int i = 0; i < n_ev; i++ )
        {
          double mag_down = xf[ SpO2prei[i]  ] - xf[ SaO2MinIdx[i] ];
          double mag_up   = xf[ SpO2posti[i] ] - xf[ SaO2MinIdx[i] ];

          if ( mag_down < 2.0 ) continue;  // Matlab minimum threshold

          // Absolute timepoints
          uint64_t tp_pre   = seg.start + (uint64_t)SpO2prei [i] * fs_tp;
          uint64_t tp_nadir = seg.start + (uint64_t)SaO2MinIdx[i] * fs_tp;
          uint64_t tp_post  = seg.start + (uint64_t)SpO2posti[i] * fs_tp;

          ml_event_t ev;
          ev.start_sec  = (double)tp_pre   / globals::tp_1sec;
          ev.stop_sec   = (double)tp_post  / globals::tp_1sec;
          ev.nadir_sec  = (double)tp_nadir / globals::tp_1sec;
          ev.nadir_val  = xf[ SaO2MinIdx[i] ];
          ev.pre_val    = xf[ SpO2prei[i]  ];
          ev.post_val   = xf[ SpO2posti[i] ];
          ev.mag_down   = mag_down;
          ev.mag_up     = mag_up;
          ev.seg        = seg_idx;
          ev.in_sleep   = has_staging
                           ? tp_in_sleep( sleep_ivs , tp_nadir )
                           : true;  // no staging → treat all as eligible

          all_events.push_back( ev );

          // Annotation: span from pre-peak to post-peak
          instance_t * inst = desat_annot->add(
            ".", interval_t( tp_pre, tp_post ), "." );
          inst->set( "nadir",    ev.nadir_val );
          inst->set( "pre",      ev.pre_val   );
          inst->set( "post",     ev.post_val  );
          inst->set( "mag_down", mag_down );
          inst->set( "mag_up",   mag_up   );
        }

      ++seg_idx;
    }  // end segment loop


  //
  // -----------------------------------------------------------------------
  // Per-event output (stratified by DESAT_M level)
  // -----------------------------------------------------------------------
  //

  int event_counter = 0;
  for ( const auto & ev : all_events )
    {
      ++event_counter;
      writer.level( event_counter , "DESAT_M" );
      writer.value( "START",    ev.start_sec );
      writer.value( "STOP",     ev.stop_sec  );
      writer.value( "DUR",      ev.stop_sec - ev.start_sec );
      writer.value( "NADIR_T",  ev.nadir_sec );
      writer.value( "NADIR",    ev.nadir_val  );
      writer.value( "PRE",      ev.pre_val    );
      writer.value( "POST",     ev.post_val   );
      writer.value( "MAG_DOWN", ev.mag_down   );
      writer.value( "MAG_UP",   ev.mag_up     );
      writer.value( "SEG",      ev.seg        );
      if ( has_staging )
        writer.value( "SLEEP", (int)ev.in_sleep );
      writer.unlevel( "DESAT_M" );
    }


  //
  // -----------------------------------------------------------------------
  // Individual-level summary: ODI2, ODI3, ODI4
  // Denominator: hours of valid signal (sleep only, if staging present)
  // -----------------------------------------------------------------------
  //

  const double denom_hr = total_t_sleep_valid / 3600.0;

  // Count events at each threshold (only those flagged in_sleep for ODI)
  int n2 = 0, n3 = 0, n4 = 0;
  int n2s = 0, n3s = 0, n4s = 0;  // sleep-gated counts
  for ( const auto & ev : all_events )
    {
      if ( ev.mag_down >= 2.0 ) { ++n2; if ( ev.in_sleep ) ++n2s; }
      if ( ev.mag_down >= 3.0 ) { ++n3; if ( ev.in_sleep ) ++n3s; }
      if ( ev.mag_down >= 4.0 ) { ++n4; if ( ev.in_sleep ) ++n4s; }
    }

  writer.value( "N2",  n2 );
  writer.value( "N3",  n3 );
  writer.value( "N4",  n4 );

  if ( denom_hr > 0 )
    {
      writer.value( "ODI2", n2s / denom_hr );
      writer.value( "ODI3", n3s / denom_hr );
      writer.value( "ODI4", n4s / denom_hr );
    }
  else
    {
      writer.value( "ODI2", 0.0 );
      writer.value( "ODI3", 0.0 );
      writer.value( "ODI4", 0.0 );
    }

  writer.value( "T_VALID",       total_t_valid );
  writer.value( "T_ART",         total_t_art   );
  writer.value( "T_SLEEP_VALID", total_t_sleep_valid );
  writer.value( "MEAN_SPO2",     n_spo2 > 0 ? sum_spo2 / n_spo2 : 0.0 );

  logger << "  DESAT (matlab) complete:"
         << "  N2=" << n2 << "  N3=" << n3 << "  N4=" << n4
         << "  (valid " << total_t_valid << "s"
         << ", sleep-valid " << total_t_sleep_valid << "s"
         << ", artifact " << total_t_art << "s)\n\n";
}
