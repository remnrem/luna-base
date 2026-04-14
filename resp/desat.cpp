
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

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>

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
