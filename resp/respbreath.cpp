
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
// RESPBREATH : Respiratory breath segmentation
//
// Algorithm summary
// -----------------
//   For each contiguous EDF segment and each input channel:
//     1. Bandpass-filter the raw signal (IIR Butterworth).
//     2. Apply light box-car smoothing to stabilise extrema.
//     3. Infer or apply polarity so inspiration is upward.
//     4. Estimate rolling local noise/amplitude (rolling MAD).
//     5. Detect alternating local extrema (troughs and peaks).
//     6. Enforce physiologic timing constraints and minimum prominence.
//     7. Construct breath events (trough → peak → trough).
//     8. Score each breath's confidence.
//     9. Mark artifact intervals where timing cannot be estimated.
//  Multi-channel:
//     10. Fuse channel results: choose the locally best-quality channel
//         and use secondary channels to rescue timing and boost confidence.
//  Output:
//     11. Emit per-breath and per-segment writer output.
//     12. Emit BREATH, INSP, EXP, RESP_ART annotations.
//

#include "resp/respbreath.h"

#include "param.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "db/db.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "annot/annot.h"
#include "defs/defs.h"
#include "dsp/iir.h"
#include "dsp/resample.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>
#include <deque>
#include <limits>

extern writer_t writer;
extern logger_t logger;


// ===========================================================================
// Internal helper types
// ===========================================================================

// Candidate extremum (before alternation enforcement)
struct extremum_t {
  int    idx;       // sample index within segment
  double val;       // filtered+smoothed value
  bool   is_max;    // true=peak, false=trough
  double prominence;
};

// Per-channel results for one contiguous segment
struct ch_seg_result_t {
  std::vector<breath_event_t>   breaths;
  std::vector<breath_artifact_t> artifacts;
  double quality_score;   // [0,1] used for fusion channel selection
};


// ===========================================================================
// Internal helpers (static)
// ===========================================================================

// Simple box-car moving average (causal, centred over ±half_win samples)
static std::vector<double> smooth_boxcar( const std::vector<double> & x , int half_win )
{
  const int n = (int)x.size();
  if ( n == 0 ) return {};

  // prefix sum for O(n) box-car smoothing
  std::vector<double> psum( n + 1, 0.0 );
  for ( int i = 0; i < n; i++ ) psum[i+1] = psum[i] + x[i];

  std::vector<double> out( n );
  for ( int i = 0; i < n; i++ )
    {
      int lo = std::max( 0, i - half_win );
      int hi = std::min( n - 1, i + half_win );
      out[i] = ( psum[hi+1] - psum[lo] ) / (double)( hi - lo + 1 );
    }
  return out;
}

// Median of a small sorted copy
static double vec_median( std::vector<double> v )
{
  if ( v.empty() ) return 0.0;
  std::sort( v.begin(), v.end() );
  int n = (int)v.size();
  if ( n % 2 == 1 ) return v[n/2];
  return 0.5 * ( v[n/2-1] + v[n/2] );
}

// Rolling MAD with a 1-second step so it adapts quickly at amplitude
// transitions (e.g. hypopneas), while remaining fast: O(N/step * W*logW).
// For 8h at 50 Hz with step=50 and win=250: ~28k blocks, well under 1s.
static std::vector<double> rolling_mad( const std::vector<double> & x , int win_samps , int step_samps )
{
  const int n = (int)x.size();
  std::vector<double> out( n, 0.0 );
  if ( n == 0 ) return out;

  const int half = win_samps / 2;

  for ( int i = 0; i < n; i += step_samps )
    {
      int lo = std::max( 0, i - half );
      int hi = std::min( n - 1, i + half );

      std::vector<double> w;
      w.reserve( hi - lo + 1 );
      for ( int j = lo; j <= hi; j++ ) w.push_back( x[j] );

      double med = vec_median( w );

      std::vector<double> ad;
      ad.reserve( w.size() );
      for ( double v : w ) ad.push_back( std::fabs( v - med ) );

      double m = vec_median( ad );

      int fill_end = std::min( n, i + step_samps );
      for ( int j = i; j < fill_end; j++ ) out[j] = m;
    }
  return out;
}

// Infer signal polarity: returns +1 if inspiration is already upward,
// -1 if signal should be flipped.
// Heuristic: compute mean of the top-quartile excursions (peaks) and
// bottom-quartile excursions (troughs) of the smoothed signal; if
// |mean_top| > |mean_bottom| the signal is already positive-up.
// This works well for nasal pressure and effort belts; falls back to +1
// if the signal is flat or uncertain.
static double infer_polarity( const std::vector<double> & x )
{
  if ( x.empty() ) return 1.0;

  // subtract mean
  double mean = 0.0;
  for ( double v : x ) mean += v;
  mean /= (double)x.size();

  double sum_pos = 0.0, sum_neg = 0.0;
  int n_pos = 0, n_neg = 0;
  for ( double v : x )
    {
      double d = v - mean;
      if ( d > 0 ) { sum_pos += d; ++n_pos; }
      else         { sum_neg -= d; ++n_neg; }
    }
  double avg_pos = ( n_pos > 0 ) ? ( sum_pos / n_pos ) : 0.0;
  double avg_neg = ( n_neg > 0 ) ? ( sum_neg / n_neg ) : 0.0;

  // If positive and negative excursions are similar, uncertain → keep as-is
  if ( avg_pos < 1e-12 && avg_neg < 1e-12 ) return 1.0;

  // nasal pressure: inspiration creates a large positive excursion
  // effort belts: depends on channel; we can't always tell from statistics alone
  // Conservative: if neither side strongly dominates, keep +1
  return ( avg_pos >= avg_neg ) ? 1.0 : -1.0;
}

// Detect local extrema in the signal.
// A sample x[i] is a local maximum if it is the maximum within a
// neighbourhood of half_win samples on each side.
// Only report the best (argmax/argmin) within each local plateau.
static std::vector<extremum_t> find_extrema( const std::vector<double> & x ,
                                               int half_win )
{
  const int n = (int)x.size();
  std::vector<extremum_t> out;
  if ( n < 3 ) return out;

  // For each sample decide if it is a strict local max or min
  // using a neighbourhood window
  for ( int i = half_win; i < n - half_win; i++ )
    {
      double v = x[i];
      bool is_max = true, is_min = true;
      for ( int k = i - half_win; k <= i + half_win; k++ )
        {
          if ( k == i ) continue;
          if ( x[k] >= v ) is_max = false;
          if ( x[k] <= v ) is_min = false;
        }
      if ( is_max || is_min )
        {
          extremum_t e;
          e.idx        = i;
          e.val        = v;
          e.is_max     = is_max;
          e.prominence = 0.0;  // computed properly in breath construction loop
          out.push_back( e );
        }
    }
  return out;
}

// Enforce strict trough-peak alternation in a list of extrema.
// When two consecutive extrema are the same type, keep the one with
// higher prominence (for maxima: larger value; for minima: smaller value).
static void enforce_alternation( std::vector<extremum_t> & ex )
{
  if ( ex.empty() ) return;

  bool changed = true;
  while ( changed )
    {
      changed = false;
      for ( int i = 1; i < (int)ex.size(); i++ )
        {
          if ( ex[i-1].is_max == ex[i].is_max )
            {
              // same type: keep the more prominent one
              int drop = ( ex[i-1].prominence >= ex[i].prominence ) ? i : i - 1;
              ex.erase( ex.begin() + drop );
              changed = true;
              break;
            }
        }
    }
}

// Score a breath's confidence on [0,1].
// Factors:
//   - I:E ratio (deviation from physiologic range)
//   - cycle duration (must be within bounds)
//   - inspiratory amplitude vs recent history
static double score_confidence( const breath_event_t & b ,
                                  double min_cycle , double max_cycle ,
                                  double min_half_cycle , double max_half_cycle ,
                                  double recent_amp_median )
{
  double score = 1.0;

  // Penalise timing outliers (soft, not hard — hard rejection is done earlier)
  if ( b.t_tot < min_cycle || b.t_tot > max_cycle )
    score *= 0.2;
  else if ( b.t_tot < min_cycle * 1.5 || b.t_tot > max_cycle * 0.75 )
    score *= 0.7;

  if ( b.t_insp < min_half_cycle || b.t_insp > max_half_cycle )
    score *= 0.4;
  if ( b.t_exp  < min_half_cycle || b.t_exp  > max_half_cycle )
    score *= 0.4;

  // Penalise implausible I:E (expect 0.3 to 0.8 roughly)
  if ( b.t_tot > 0 )
    {
      double ie_frac = b.t_insp / b.t_tot;
      if ( ie_frac < 0.15 || ie_frac > 0.85 ) score *= 0.5;
    }

  // Penalise very small amplitude relative to recent good breaths
  if ( recent_amp_median > 1e-12 )
    {
      double amp_frac = b.amp_insp / recent_amp_median;
      if ( amp_frac < 0.1 ) score *= 0.1;
      else if ( amp_frac < 0.3 ) score *= 0.5;
    }

  return std::max( 0.0, std::min( 1.0, score ) );
}

// Merge overlapping / nearby artifact intervals into a single list,
// and drop intervals shorter than min_art_sec.
static std::vector<breath_artifact_t> merge_artifacts(
    const std::vector<breath_artifact_t> & raw ,
    double merge_gap_sec ,
    double min_art_sec )
{
  if ( raw.empty() ) return raw;

  // sort by start
  std::vector<breath_artifact_t> sorted = raw;
  std::sort( sorted.begin(), sorted.end(),
             []( const breath_artifact_t & a, const breath_artifact_t & b )
             { return a.start_sec < b.start_sec; } );

  std::vector<breath_artifact_t> merged;
  breath_artifact_t cur = sorted[0];

  for ( int i = 1; i < (int)sorted.size(); i++ )
    {
      if ( sorted[i].start_sec <= cur.end_sec + merge_gap_sec )
        {
          // extend
          if ( sorted[i].end_sec > cur.end_sec )
            {
              cur.end_sec  = sorted[i].end_sec;
              cur.tp_end   = sorted[i].tp_end;
            }
        }
      else
        {
          if ( cur.end_sec - cur.start_sec >= min_art_sec )
            merged.push_back( cur );
          cur = sorted[i];
        }
    }
  if ( cur.end_sec - cur.start_sec >= min_art_sec )
    merged.push_back( cur );

  return merged;
}


// ===========================================================================
// Per-channel processing for one contiguous segment
// ===========================================================================

static ch_seg_result_t process_channel_segment(
    const std::vector<double> & raw ,
    int fs ,
    uint64_t seg_start_tp ,   // absolute timepoint of segment start
    double   seg_start_sec ,  // seconds
    int seg_idx ,
    int ch_idx ,
    // parameters
    double hp_hz ,
    double lp_hz ,
    double smooth_sec ,
    double min_half_cycle ,
    double max_half_cycle ,
    double min_cycle ,
    double max_cycle ,
    double prom_z ,
    double amp_rel ,
    double recover_sec ,
    double merge_art_sec ,
    double min_art_sec ,
    bool   use_hilbert ,
    double flip_polarity , // +1 or -1
    bool   auto_flip ,
    double conf_thresh ,
    double flatline_thresh_z ,  // MAD-relative threshold for flatline
    double min_seg_sec          // minimum segment length to attempt analysis
    )
{
  ch_seg_result_t result;
  result.quality_score = 0.0;

  const int N = (int)raw.size();
  if ( N < (int)( min_seg_sec * fs ) )
    {
      // Too short: mark entire segment as artifact
      breath_artifact_t art;
      art.start_sec = seg_start_sec;
      art.end_sec   = seg_start_sec + (double)N / fs;
      art.tp_start  = seg_start_tp;
      art.tp_end    = seg_start_tp + (uint64_t)round( (double)N / fs * globals::tp_1sec );
      art.source_ch = ch_idx;
      art.reason    = "segment_too_short";
      result.artifacts.push_back( art );
      return result;
    }

  
  // -----------------------------------------------------------------------
  // Step 1: Bandpass filter (IIR Butterworth)
  // -----------------------------------------------------------------------

  std::vector<double> filt;
  try {
    iir_t iir;
    iir.init( BUTTERWORTH_BANDPASS , 2 , (double)fs , hp_hz , lp_hz );
    filt = iir.apply( raw );
  } catch (...) {
    // filter init failed (e.g. frequency too low for sample rate)
    // mark entire segment as artifact
    breath_artifact_t art;
    art.start_sec = seg_start_sec;
    art.end_sec   = seg_start_sec + (double)N / fs;
    art.tp_start  = seg_start_tp;
    art.tp_end    = seg_start_tp + (uint64_t)round( (double)N / fs * globals::tp_1sec );
    art.source_ch = ch_idx;
    art.reason    = "filter_failed";
    result.artifacts.push_back( art );
    return result;
  }

  // -----------------------------------------------------------------------
  // Step 2: Light smoothing
  // -----------------------------------------------------------------------

  const int smooth_samps = std::max( 1, (int)round( smooth_sec * fs / 2.0 ) );
  std::vector<double> sig = smooth_boxcar( filt, smooth_samps );

  // -----------------------------------------------------------------------
  // Step 3: Polarity
  // -----------------------------------------------------------------------

  double polarity = flip_polarity;
  if ( auto_flip )
    polarity = infer_polarity( sig );

  if ( polarity < 0 )
    for ( double & v : sig ) v = -v;

  // -----------------------------------------------------------------------
  // Step 4: Rolling MAD for local noise / amplitude estimate
  // -----------------------------------------------------------------------

  // Use ~10 second windows for noise estimation
  const int mad_win  = std::max( 3,  (int)round( 5.0 * fs ) );  // 5s window
  const int mad_step = std::max( 1,  (int)round( 1.0 * fs ) );  // 1s step
  std::vector<double> mad = rolling_mad( sig, mad_win, mad_step );

  // Global signal MAD (for flatline detection)
  double global_mad = vec_median( mad );

  // -----------------------------------------------------------------------
  // Step 5: Flatline / all-artifact detection
  // -----------------------------------------------------------------------

  // If the entire segment is nearly flat, mark it as artifact
  if ( global_mad < 1e-12 )
    {
      breath_artifact_t art;
      art.start_sec = seg_start_sec;
      art.end_sec   = seg_start_sec + (double)N / fs;
      art.tp_start  = seg_start_tp;
      art.tp_end    = seg_start_tp + (uint64_t)round( (double)N / fs * globals::tp_1sec );
      art.source_ch = ch_idx;
      art.reason    = "flatline";
      result.artifacts.push_back( art );
      return result;
    }

  // -----------------------------------------------------------------------
  // Step 6: Optional Hilbert envelope (for future use / amplitude aid)
  // -----------------------------------------------------------------------
  // TODO: if use_hilbert is true, compute hilbert_t envelope of the bandpass
  // signal and use it as an alternative amplitude/noise estimate.
  // For now, we rely on rolling MAD.

  if ( use_hilbert ) Helper::halt( "use-hilbert not implemented" );
  
  // -----------------------------------------------------------------------
  // Step 7: Find candidate extrema
  // -----------------------------------------------------------------------

  // neighbourhood window for extrema detection: half a minimum half-cycle
  const int ext_half_win = std::max( 2, (int)round( min_half_cycle * fs * 0.4 ) );

  std::vector<extremum_t> extr = find_extrema( sig, ext_half_win );

  if ( extr.size() < 2 )
    {
      breath_artifact_t art;
      art.start_sec = seg_start_sec;
      art.end_sec   = seg_start_sec + (double)N / fs;
      art.tp_start  = seg_start_tp;
      art.tp_end    = seg_start_tp + (uint64_t)round( (double)N / fs * globals::tp_1sec );
      art.source_ch = ch_idx;
      art.reason    = "no_extrema";
      result.artifacts.push_back( art );
      return result;
    }

  // -----------------------------------------------------------------------
  // Step 8+9: Enforce strict trough-peak alternation
  // (prominence is checked per-triplet in Step 10 using actual peak-trough
  // amplitude, which is DC-invariant and physically meaningful)
  // -----------------------------------------------------------------------

  enforce_alternation( extr );

  // Make sure we start with a trough (if inspiration is upward, a breath
  // starts at a trough)
  while ( !extr.empty() && extr[0].is_max )
    extr.erase( extr.begin() );

  // -----------------------------------------------------------------------
  // Step 10: Construct breath candidates
  //   Breath i = trough[i] → peak[i] → trough[i+1]
  //   Requires alternating: trough, peak, trough, peak, ...
  // -----------------------------------------------------------------------

  const int min_half_samps = (int)round( min_half_cycle * fs );
  const int max_half_samps = (int)round( max_half_cycle * fs );
  const int min_cycle_samps = (int)round( min_cycle * fs );
  const int max_cycle_samps = (int)round( max_cycle * fs );

  // Rolling median amplitude from recent good breaths (for amp_rel threshold)
  std::deque<double> recent_amps;
  const int recent_n = 5;  // look back at last 5 good breaths
  double recent_amp_median = 0.0;

  std::vector<breath_event_t> breaths;
  std::vector<breath_artifact_t> raw_arts;

  // Iterate over trough-peak-trough triplets
  for ( int i = 0; i + 2 < (int)extr.size(); i += 2 )
    {
      // extr[i]   = trough T_i   (is_max == false)
      // extr[i+1] = peak P_i     (is_max == true)
      // extr[i+2] = trough T_{i+1} (is_max == false)

      if ( extr[i].is_max || !extr[i+1].is_max || extr[i+2].is_max )
        {
          // Alternation violated — skip one step
          // This should not happen after enforce_alternation, but guard anyway
          continue;
        }

      const extremum_t & T0 = extr[i];
      const extremum_t & P  = extr[i+1];
      const extremum_t & T1 = extr[i+2];

      int dt_insp  = P.idx  - T0.idx;
      int dt_exp   = T1.idx - P.idx;
      int dt_cycle = T1.idx - T0.idx;

      // Hard timing constraints
      if ( dt_insp  < min_half_samps || dt_insp  > max_half_samps ) continue;
      if ( dt_exp   < min_half_samps || dt_exp   > max_half_samps ) continue;
      if ( dt_cycle < min_cycle_samps || dt_cycle > max_cycle_samps ) continue;

      // Amplitude constraint relative to recent good breaths
      double amp_insp = P.val - T0.val;
      double amp_exp  = P.val - T1.val;
      double amp_sym  = P.val - 0.5 * ( T0.val + T1.val );

      if ( amp_insp <= 0.0 ) continue;  // peak must be above trough

      // Prominence: peak-to-trough amplitude normalised by local MAD at the peak.
      // Using amp_sym (peak - mean adjacent troughs) makes this DC-invariant.
      double noise_at_peak = ( mad[ P.idx ] > 1e-15 ) ? mad[ P.idx ] : 1e-15;
      double prominence    = amp_sym / noise_at_peak;
      if ( prominence < prom_z ) continue;

      if ( recent_amp_median > 1e-12 )
        if ( amp_insp < amp_rel * recent_amp_median ) continue;

      // Build breath event
      breath_event_t b;
      b.start_sec  = seg_start_sec + (double)T0.idx / fs;
      b.peak_sec   = seg_start_sec + (double)P.idx  / fs;
      b.end_sec    = seg_start_sec + (double)T1.idx / fs;
      b.t_insp     = (double)dt_insp  / fs;
      b.t_exp      = (double)dt_exp   / fs;
      b.t_tot      = (double)dt_cycle / fs;
      b.amp_insp   = amp_insp;
      b.amp_exp    = amp_exp;
      b.amp_sym    = amp_sym;
      b.source_ch  = ch_idx;
      b.seg_idx    = seg_idx;
      b.fused      = false;
      b.n_supporting_ch = 1;

      b.tp_start   = seg_start_tp + (uint64_t)round( (double)T0.idx / fs * globals::tp_1sec );
      b.tp_peak    = seg_start_tp + (uint64_t)round( (double)P.idx  / fs * globals::tp_1sec );
      b.tp_end     = seg_start_tp + (uint64_t)round( (double)T1.idx / fs * globals::tp_1sec );

      // Score confidence
      b.confidence = score_confidence( b,
                                        min_cycle, max_cycle,
                                        min_half_cycle, max_half_cycle,
                                        recent_amp_median );
      b.low_conf   = ( b.confidence < conf_thresh );

      // Update rolling amplitude history from good breaths
      if ( !b.low_conf )
        {
          recent_amps.push_back( amp_insp );
          if ( (int)recent_amps.size() > recent_n ) recent_amps.pop_front();
          std::vector<double> tmp( recent_amps.begin(), recent_amps.end() );
          recent_amp_median = vec_median( tmp );
        }

      breaths.push_back( b );
    }

  // -----------------------------------------------------------------------
  // Step 11: Detect artifact intervals
  //   Regions not covered by accepted breaths are candidate artifacts.
  //   Also mark regions where there were no extrema to work with.
  // -----------------------------------------------------------------------

  // Build a per-sample "covered" flag
  std::vector<bool> covered( N, false );
  for ( const breath_event_t & b : breaths )
    {
      int i0 = (int)round( (b.start_sec - seg_start_sec) * fs );
      int i1 = (int)round( (b.end_sec   - seg_start_sec) * fs );
      i0 = std::max( 0, std::min( N-1, i0 ) );
      i1 = std::max( 0, std::min( N-1, i1 ) );
      for ( int k = i0; k <= i1; k++ ) covered[k] = true;
    }

  // Any run of uncovered samples >= min_art_sec → artifact
  {
    bool in_gap = false;
    int  gap_start = 0;

    // We also want to catch the region before first breath and after last
    for ( int i = 0; i <= N; i++ )
      {
        bool unc = ( i < N ) ? !covered[i] : false;
        if ( !in_gap && unc ) { in_gap = true; gap_start = i; }
        else if ( in_gap && !unc )
          {
            double dur = (double)( i - gap_start ) / fs;
            if ( dur >= min_art_sec )
              {
                breath_artifact_t art;
                art.start_sec = seg_start_sec + (double)gap_start / fs;
                art.end_sec   = seg_start_sec + (double)i / fs;
                art.tp_start  = seg_start_tp + (uint64_t)round( (double)gap_start / fs * globals::tp_1sec );
                art.tp_end    = seg_start_tp + (uint64_t)round( (double)i         / fs * globals::tp_1sec );
                art.source_ch = ch_idx;
                art.reason    = "uncovered";
                raw_arts.push_back( art );
              }
            in_gap = false;
          }
      }
  }

  // Flag low-confidence breaths that are adjacent to artifact regions
  for ( breath_event_t & b : breaths )
    {
      for ( const breath_artifact_t & a : raw_arts )
        {
          if ( std::fabs( b.start_sec - a.end_sec   ) < recover_sec ||
               std::fabs( b.end_sec   - a.start_sec ) < recover_sec )
            b.near_artifact = true;
        }
    }

  // -----------------------------------------------------------------------
  // Step 12: Compute channel quality score for fusion
  // -----------------------------------------------------------------------

  double total_sec   = (double)N / fs;
  double covered_sec = 0.0;
  double high_conf_n = 0.0;
  for ( const breath_event_t & b : breaths )
    {
      covered_sec += b.t_tot;
      if ( !b.low_conf ) high_conf_n += 1.0;
    }

  double frac_covered   = ( total_sec > 0 ) ? ( covered_sec / total_sec ) : 0.0;
  double frac_high_conf = ( breaths.empty() ) ? 0.0 : ( high_conf_n / breaths.size() );

  result.quality_score = 0.5 * frac_covered + 0.5 * frac_high_conf;

  // -----------------------------------------------------------------------
  // Merge artifact intervals
  // -----------------------------------------------------------------------

  result.artifacts = merge_artifacts( raw_arts, merge_art_sec, min_art_sec );
  result.breaths   = breaths;

  return result;
}


// ===========================================================================
// Multi-channel fusion
// ===========================================================================

static std::vector<breath_event_t> fuse_channels(
    const std::vector<ch_seg_result_t> & ch_results ,
    const std::vector<std::string>      & ch_labels ,
    int primary_ch ,
    double fuse_window_sec ,
    double conf_thresh )
{
  const int nch = (int)ch_results.size();
  if ( nch == 0 ) return {};
  if ( nch == 1 ) return ch_results[0].breaths;

  // Choose locally best channel (highest quality score)
  int best_ch = primary_ch;
  double best_q = ch_results[ primary_ch ].quality_score;
  for ( int c = 0; c < nch; c++ )
    if ( ch_results[c].quality_score > best_q )
      { best_q = ch_results[c].quality_score; best_ch = c; }

  // Start from the best channel's breath list
  std::vector<breath_event_t> fused = ch_results[ best_ch ].breaths;

  // For each fused breath, check how many other channels have a matching breath
  for ( breath_event_t & fb : fused )
    {
      int n_support = 1;  // best_ch always supports itself

      for ( int c = 0; c < nch; c++ )
        {
          if ( c == best_ch ) continue;
          for ( const breath_event_t & ob : ch_results[c].breaths )
            {
              // Match by peak time within fuse_window
              if ( std::fabs( ob.peak_sec - fb.peak_sec ) <= fuse_window_sec )
                {
                  ++n_support;
                  // Boost confidence if secondary agrees
                  fb.confidence = std::min( 1.0, fb.confidence + 0.15 );
                  break;
                }
            }
        }

      fb.n_supporting_ch = n_support;
      fb.fused = ( n_support > 1 );
      if ( n_support > 1 ) fb.low_conf = ( fb.confidence < conf_thresh );
    }

  // Rescue: if a breath in a secondary channel has no match in the primary
  // breath list but high confidence, add it (this handles temporary primary
  // channel dropouts)
  for ( int c = 0; c < nch; c++ )
    {
      if ( c == best_ch ) continue;
      for ( const breath_event_t & ob : ch_results[c].breaths )
        {
          if ( ob.low_conf ) continue;

          bool matched = false;
          for ( const breath_event_t & fb : fused )
            if ( std::fabs( ob.peak_sec - fb.peak_sec ) <= fuse_window_sec )
              { matched = true; break; }

          if ( !matched )
            {
              // Rescued breath from secondary channel
              breath_event_t rb = ob;
              rb.source_ch = c;
              rb.fused     = true;
              rb.n_supporting_ch = 1;
              fused.push_back( rb );
            }
        }
    }

  // Re-sort by start time
  std::sort( fused.begin(), fused.end(),
             []( const breath_event_t & a , const breath_event_t & b )
             { return a.start_sec < b.start_sec; } );

  return fused;
}


// ===========================================================================
// Summary statistics
// ===========================================================================

struct breath_summary_t {
  int    n_breath;
  double breath_rate;     // breaths/min
  double mean_ttot;
  double median_ttot;
  double sd_ttot;
  double cv_ttot;
  double median_tinsp;
  double median_texp;
  double ie_ratio;
  double median_amp;
  double lowconf_pct;
  double artifact_pct;   // fraction of segment time covered by artifact
  int    n_fused;
  double fused_pct;
};

static breath_summary_t compute_summary( const std::vector<breath_event_t> & breaths ,
                                          const std::vector<breath_artifact_t> & arts ,
                                          double total_sec )
{
  breath_summary_t s;
  s.n_breath     = (int)breaths.size();
  s.breath_rate  = ( total_sec > 0 ) ? ( s.n_breath / total_sec * 60.0 ) : 0.0;
  s.mean_ttot    = 0.0;
  s.median_ttot  = 0.0;
  s.sd_ttot      = 0.0;
  s.cv_ttot      = 0.0;
  s.median_tinsp = 0.0;
  s.median_texp  = 0.0;
  s.ie_ratio     = 0.0;
  s.median_amp   = 0.0;
  s.lowconf_pct  = 0.0;
  s.artifact_pct = 0.0;
  s.n_fused      = 0;
  s.fused_pct    = 0.0;

  if ( breaths.empty() ) return s;

  std::vector<double> ttot_v, tinsp_v, texp_v, amp_v;
  int n_low = 0;
  int n_fused = 0;

  for ( const breath_event_t & b : breaths )
    {
      ttot_v.push_back( b.t_tot );
      tinsp_v.push_back( b.t_insp );
      texp_v.push_back( b.t_exp );
      amp_v.push_back( b.amp_sym );
      if ( b.low_conf ) ++n_low;
      if ( b.fused    ) ++n_fused;
    }

  // mean / sd of ttot
  double sum = 0.0;
  for ( double v : ttot_v ) sum += v;
  s.mean_ttot = sum / ttot_v.size();

  double sumsq = 0.0;
  for ( double v : ttot_v ) sumsq += ( v - s.mean_ttot ) * ( v - s.mean_ttot );
  s.sd_ttot = ( ttot_v.size() > 1 ) ? std::sqrt( sumsq / ( ttot_v.size() - 1 ) ) : 0.0;
  s.cv_ttot = ( s.mean_ttot > 0 ) ? ( s.sd_ttot / s.mean_ttot ) : 0.0;

  s.median_ttot  = vec_median( ttot_v );
  s.median_tinsp = vec_median( tinsp_v );
  s.median_texp  = vec_median( texp_v );
  s.median_amp   = vec_median( amp_v );
  s.ie_ratio     = ( s.median_texp > 0 ) ? ( s.median_tinsp / s.median_texp ) : 0.0;

  s.lowconf_pct  = 100.0 * n_low  / breaths.size();
  s.n_fused      = n_fused;
  s.fused_pct    = 100.0 * n_fused / breaths.size();

  // artifact fraction of total time
  double art_sec = 0.0;
  for ( const breath_artifact_t & a : arts ) art_sec += ( a.end_sec - a.start_sec );
  s.artifact_pct = ( total_sec > 0 ) ? ( 100.0 * art_sec / total_sec ) : 0.0;

  return s;
}


// ===========================================================================
// Main constructor
// ===========================================================================

respbreath_t::respbreath_t( edf_t & edf , param_t & param )
{
  //
  // -------------------------------------------------------------------------
  // Parameters
  // -------------------------------------------------------------------------
  //

  std::string signal_label = param.requires( "sig" );

  const bool no_annotations = true;
    
  signal_list_t signals = edf.header.signal_list( signal_label , no_annotations );

  const int ns = signals.size();
  
  if ( ns == 0 )
    {
      logger << "  no matching channels found, bailing\n";
      return;
    }


  // enforce <= 50 Hz signals

  for (int s=0; s<ns; s++)
    {
      if ( edf.header.sampling_freq( signals(s) ) > 50 )
	dsptools::resample_channel( edf, signals(s) , 50 );
    }
  
  
  // Primary channel selection: first listed, unless override given
  std::string primary_label = param.has( "primary" ) ? param.value( "primary" ) : "first";

  // Polarity
  const std::string flip_str = param.has( "flip" ) ? param.value( "flip" ) : "auto";
  bool   auto_flip    = ( flip_str == "auto" );
  double flip_polarity = 1.0;
  if ( flip_str == "yes" || flip_str == "1" ) { auto_flip = false; flip_polarity = -1.0; }
  if ( flip_str == "no"  || flip_str == "0" ) { auto_flip = false; flip_polarity =  1.0; }

  // Filter parameters
  const double hp_hz      = param.has( "hp"     ) ? param.requires_dbl( "hp"     ) : 0.03;
  const double lp_hz      = param.has( "lp"     ) ? param.requires_dbl( "lp"     ) : 1.5;

  // Smoothing
  const double smooth_sec = param.has( "smooth" ) ? param.requires_dbl( "smooth" ) : 0.20;

  // Timing constraints
  const double min_half_cycle = param.has( "min-half-cycle" ) ? param.requires_dbl( "min-half-cycle" ) : 0.40;
  const double max_half_cycle = param.has( "max-half-cycle" ) ? param.requires_dbl( "max-half-cycle" ) : 6.00;
  const double min_cycle      = param.has( "min-cycle"      ) ? param.requires_dbl( "min-cycle"      ) : 0.80;
  const double max_cycle      = param.has( "max-cycle"      ) ? param.requires_dbl( "max-cycle"      ) : 12.0;

  // Quality thresholds
  const double prom_z         = param.has( "prom-z"   ) ? param.requires_dbl( "prom-z"   ) : 1.0;
  const double amp_rel        = param.has( "amp-rel"  ) ? param.requires_dbl( "amp-rel"  ) : 0.15;
  const double conf_thresh    = param.has( "conf"     ) ? param.requires_dbl( "conf"     ) : 0.40;

  // Artifact handling
  const double recover_sec   = param.has( "recover"   ) ? param.requires_dbl( "recover"   ) : 2.0;
  const double merge_art_sec = param.has( "merge-art" ) ? param.requires_dbl( "merge-art" ) : 1.0;
  const double min_art_sec   = param.has( "min-art"   ) ? param.requires_dbl( "min-art"   ) : 2.0;
  const double min_seg_sec   = param.has( "min-seg"   ) ? param.requires_dbl( "min-seg"   ) : 5.0;
  const double flatline_z    = 0.01;  // internal only

  // Multi-channel options
  const bool   do_fuse     = param.has( "fuse" ) ? ( param.value( "fuse" ) == "yes" || param.value( "fuse" ) == "1" ) : true;
  const double fuse_window = param.has( "fuse-window" ) ? param.requires_dbl( "fuse-window" ) : 0.75;

  // Hilbert
  const bool use_hilbert = param.has( "use-hilbert" )
    ? ( param.value( "use-hilbert" ) == "yes" || param.value( "use-hilbert" ) == "1" )
    : false;

  // Annotation labels
  const std::string annot_breath = param.has( "annot-breath" ) ? param.value( "annot-breath" ) : "BREATH";
  const std::string annot_insp   = param.has( "annot-insp"   ) ? param.value( "annot-insp"   ) : "INSP";
  const std::string annot_exp    = param.has( "annot-exp"    ) ? param.value( "annot-exp"    ) : "EXP";
  const std::string annot_art    = param.has( "annot-art"    ) ? param.value( "annot-art"    ) : "RESP_ART";

  const bool emit_insp_exp = param.has( "annot-insp" ) || param.has( "annot-exp" );
  const bool verbose       = param.has( "verbose" ) && ( param.value( "verbose" ) == "yes" );

  //
  // -------------------------------------------------------------------------
  // Log parameters
  // -------------------------------------------------------------------------
  //
  
  logger << "\n  respiratory breath segmentation\n";
  logger << "    channels           : ";
  for ( int ci = 0; ci < ns; ci++ )
    logger << signals.label(ci) << ( ci + 1 < ns ? "," : "" );
  logger << "\n";
  logger << "    bandpass           : " << hp_hz << " – " << lp_hz << " Hz\n";
  logger << "    cycle constraints  : " << min_cycle << " – " << max_cycle << " s\n";
  logger << "    half-cycle         : " << min_half_cycle << " – " << max_half_cycle << " s\n";
  logger << "    prom_z / amp_rel   : " << prom_z << " / " << amp_rel << "\n";
  logger << "    polarity           : " << flip_str << "\n";
  logger << "    fusion             : " << ( do_fuse ? "yes" : "no" ) << "\n";
  logger << "    breath annotation  : " << annot_breath << "\n";
  logger << "    artifact annotation: " << annot_art << "\n\n";

  //
  // -------------------------------------------------------------------------
  // Resolve signal indices
  // -------------------------------------------------------------------------
  //

  std::vector<int> sig_ns;
  std::vector<std::string> valid_labels;

  for ( int s=0; s<ns; s++ )
    {

      const std::string lbl = signals.label(s);
      
      int sn = edf.header.signal( lbl );

      if ( sn == -1 )
        {
          logger << "  ** could not find signal '" << lbl << "', skipping\n";
          continue;
        }

      int fs_check = edf.header.sampling_freq( sn );

      if ( fs_check <= 0 )
        {
          logger << "  ** signal '" << lbl << "' has zero sample rate, skipping\n";
          continue;
        }

      sig_ns.push_back( sn );

      valid_labels.push_back( lbl );
    }

  if ( sig_ns.empty() )
    {
      logger << "  ** no valid signals found, aborting\n";
      return;
    }

  // Determine primary channel index
  int primary_ch = 0;
  if ( primary_label != "first" && primary_label != "auto" )
    {
      for ( int ci = 0; ci < (int)valid_labels.size(); ci++ )
        if ( valid_labels[ci] == primary_label )
          { primary_ch = ci; break; }
    }

  //
  // -------------------------------------------------------------------------
  // Create annotation objects
  // -------------------------------------------------------------------------
  //

  annot_t * ann_breath = edf.annotations->add( annot_breath );
  annot_t * ann_art    = edf.annotations->add( annot_art );
  annot_t * ann_insp   = emit_insp_exp ? edf.annotations->add( annot_insp ) : nullptr;
  annot_t * ann_exp    = emit_insp_exp ? edf.annotations->add( annot_exp  ) : nullptr;

  //
  // -------------------------------------------------------------------------
  // Get contiguous segments
  // -------------------------------------------------------------------------
  //

  std::set<interval_t> segs = edf.timeline.segments();
  if ( segs.empty() )
    {
      logger << "  ** RESPBREATH: no segments found\n";
      return;
    }

  logger << "    segments           : " << segs.size() << "\n\n";

  //
  // -------------------------------------------------------------------------
  // Accumulate global stats
  // -------------------------------------------------------------------------
  //

  std::vector<breath_event_t>   all_breaths;
  std::vector<breath_artifact_t> all_artifacts;
  double total_rec_sec = 0.0;

  int seg_idx = 0;
  int breath_counter = 0;

  //
  // -------------------------------------------------------------------------
  // Per-segment loop
  // -------------------------------------------------------------------------
  //

  for ( const interval_t & seg : segs )
    {
      const double seg_start_sec = (double)seg.start / globals::tp_1sec;
      const double seg_end_sec   = (double)seg.stop  / globals::tp_1sec;
      const double seg_dur_sec   = seg_end_sec - seg_start_sec;
      total_rec_sec += seg_dur_sec;

      // Process each channel for this segment
      std::vector<ch_seg_result_t> ch_results;
      ch_results.reserve( sig_ns.size() );

      for ( int ci = 0; ci < (int)sig_ns.size(); ci++ )
        {
          int sig_n = sig_ns[ci];
          int fs    = edf.header.sampling_freq( sig_n );

          slice_t slice( edf, sig_n, seg );
          const std::vector<double> * pdata = slice.pdata();
          if ( pdata == nullptr || pdata->empty() )
            {
              ch_seg_result_t empty;
              empty.quality_score = 0.0;
              // mark whole seg as artifact for this channel
              breath_artifact_t art;
              art.start_sec = seg_start_sec;
              art.end_sec   = seg_end_sec;
              art.tp_start  = seg.start;
              art.tp_end    = seg.stop;
              art.source_ch = ci;
              art.reason    = "no_data";
              empty.artifacts.push_back( art );
              ch_results.push_back( empty );
              continue;
            }

          ch_seg_result_t res = process_channel_segment(
              *pdata , fs ,
              seg.start , seg_start_sec ,
              seg_idx , ci ,
              hp_hz , lp_hz ,
              smooth_sec ,
              min_half_cycle , max_half_cycle ,
              min_cycle , max_cycle ,
              prom_z , amp_rel ,
              recover_sec , merge_art_sec , min_art_sec ,
              use_hilbert ,
              flip_polarity , auto_flip ,
              conf_thresh ,
              flatline_z ,
              min_seg_sec );

          ch_results.push_back( res );
        }

      // -----------------------------------------------------------------------
      // Multi-channel fusion (or just use single channel result)
      // -----------------------------------------------------------------------

      std::vector<breath_event_t>    seg_breaths;
      std::vector<breath_artifact_t> seg_arts;

      if ( do_fuse && ch_results.size() > 1 )
        {
          seg_breaths = fuse_channels( ch_results, valid_labels,
                                        primary_ch, fuse_window, conf_thresh );
          // Use artifacts from best-quality channel for annotation
          double best_q = -1.0;
          int best_c = 0;
          for ( int c = 0; c < (int)ch_results.size(); c++ )
            if ( ch_results[c].quality_score > best_q )
              { best_q = ch_results[c].quality_score; best_c = c; }
          seg_arts = ch_results[ best_c ].artifacts;
        }
      else
        {
          // Single channel (or fusion disabled): use primary channel
          seg_breaths = ch_results[ primary_ch ].breaths;
          seg_arts    = ch_results[ primary_ch ].artifacts;
        }

      // -----------------------------------------------------------------------
      // (no per-segment writer output; individual-level summary covers all segments combined)

      // -----------------------------------------------------------------------
      // Per-breath writer output and annotation emission
      // -----------------------------------------------------------------------

      for ( breath_event_t & b : seg_breaths )
        {
          b.breath_idx = breath_counter++;

          // Writer output strat by CH x BREATH
          const std::string ch_label = ( b.source_ch >= 0 && b.source_ch < (int)valid_labels.size() )
            ? valid_labels[ b.source_ch ] : "?";
          writer.level( ch_label      , "CH"     );
          writer.level( b.breath_idx  , "N" );
          writer.value( "START"    , b.start_sec  );
          writer.value( "PEAK"     , b.peak_sec   );
          writer.value( "END"      , b.end_sec    );
          writer.value( "TINSP"    , b.t_insp     );
          writer.value( "TEXP"     , b.t_exp      );
          writer.value( "TTOT"     , b.t_tot      );
          writer.value( "AMP_INSP" , b.amp_insp   );
          writer.value( "AMP_EXP"  , b.amp_exp    );
          writer.value( "AMP_SYM"  , b.amp_sym    );
          writer.value( "CONF"     , b.confidence );
          writer.value( "LOW_CONF" , (int)b.low_conf );
          writer.value( "FUSED"    , (int)b.fused    );
          writer.value( "N_SUPPORT", b.n_supporting_ch );
          writer.unlevel( "N"  );
          writer.unlevel( "CH" );

          // BREATH annotation (start → end)
          instance_t * inst = ann_breath->add(
              ".", interval_t( b.tp_start, b.tp_end ), "." );
          inst->set( "peak"      , b.peak_sec   );
          inst->set( "tinsp"     , b.t_insp     );
          inst->set( "texp"      , b.t_exp      );
          inst->set( "ttot"      , b.t_tot      );
          inst->set( "conf"  , b.confidence );
          inst->set( "fused" , (int)b.fused );

          // Optional INSP / EXP sub-interval annotations
          if ( emit_insp_exp && ann_insp && ann_exp )
            {
              ann_insp->add( ".", interval_t( b.tp_start, b.tp_peak ), "." );
              ann_exp->add(  ".", interval_t( b.tp_peak,  b.tp_end  ), "." );
            }
        }

      // -----------------------------------------------------------------------
      // Artifact annotation emission
      // -----------------------------------------------------------------------

      for ( const breath_artifact_t & a : seg_arts )
        {
          instance_t * inst = ann_art->add(
              ".", interval_t( a.tp_start, a.tp_end ), "." );
          inst->set( "reason", a.reason );
          if ( a.source_ch >= 0 && a.source_ch < (int)valid_labels.size() )
            inst->set( "channel", valid_labels[ a.source_ch ] );
        }

      // Accumulate for global summary
      for ( const breath_event_t   & b : seg_breaths ) all_breaths.push_back( b );
      for ( const breath_artifact_t & a : seg_arts    ) all_artifacts.push_back( a );

      ++seg_idx;
    } // end segment loop


  //
  // -------------------------------------------------------------------------
  // Global (individual-level) summary output
  // -------------------------------------------------------------------------
  //

  const breath_summary_t gsumm = compute_summary( all_breaths, all_artifacts, total_rec_sec );

  writer.value( "N_BREATH"    , gsumm.n_breath    );
  writer.value( "BREATH_RATE" , gsumm.breath_rate );
  writer.value( "MEAN_TTOT"   , gsumm.mean_ttot   );
  writer.value( "MEDIAN_TTOT" , gsumm.median_ttot );
  writer.value( "SD_TTOT"     , gsumm.sd_ttot     );
  writer.value( "CV_TTOT"     , gsumm.cv_ttot     );
  writer.value( "MEDIAN_TINSP", gsumm.median_tinsp );
  writer.value( "MEDIAN_TEXP" , gsumm.median_texp  );
  writer.value( "IE_RATIO"    , gsumm.ie_ratio     );
  writer.value( "MEDIAN_AMP"  , gsumm.median_amp   );
  writer.value( "LOWCONF_PCT" , gsumm.lowconf_pct  );
  writer.value( "ARTIFACT_PCT", gsumm.artifact_pct );
  writer.value( "FUSED_PCT"   , gsumm.fused_pct    );

  // Per-channel support breakdown
  if ( valid_labels.size() > 1 )
    {
      std::vector<int> ch_primary_count( valid_labels.size(), 0 );
      for ( const breath_event_t & b : all_breaths )
        if ( b.source_ch >= 0 && b.source_ch < (int)valid_labels.size() )
          ch_primary_count[ b.source_ch ]++;

      for ( int ci = 0; ci < (int)valid_labels.size(); ci++ )
        {
          writer.level( valid_labels[ci] , "CH" );
          double pct = ( all_breaths.empty() ) ? 0.0
                       : 100.0 * ch_primary_count[ci] / all_breaths.size();
          writer.value( "PRIMARY_USED_PCT", pct );
          writer.unlevel( "CH" );
        }
    }

  logger << "  RESPBREATH complete : " << all_breaths.size() << " breaths detected"
         << " across " << seg_idx << " segment(s)"
         << " (" << all_artifacts.size() << " artifact interval(s))\n\n";
}
