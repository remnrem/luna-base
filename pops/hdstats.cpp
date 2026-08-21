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
//  HDSTATS: analysis of hypnodensity signals (posterior-probability EDF channels)
//
//  Organizes output around three core domains:
//    MIXEDNESS            - how diffuse is the posterior at each sample
//    INSTABILITY          - how much does the posterior move over time
//    TRANSITION STRUCTURE - what happens near likely state boundaries
//
//    i.e. stable periods  vs.  unstable (mixed/intermediate/transitional) periods
//

#include "hdstats.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "annot/annot.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "defs/defs.h"
#include "db/db.h"
#include "param.h"

#include <cmath>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <map>
#include <set>
#include <limits>

extern logger_t logger;
extern writer_t writer;

namespace {

std::string make_invalid_channel_label( const std::string & prefix )
{
  return prefix + "_INVALID";
}

std::string infer_invalid_channel_label( const std::vector<std::string> & ch )
{
  if ( ch.empty() ) return "PP_INVALID";

  std::string prefix = "";
  for ( int i = 0; i < (int)ch.size(); i++ )
    {
      const std::string & label = ch[i];
      const std::string::size_type pos = label.rfind( '_' );
      if ( pos == std::string::npos || pos == 0 ) return "PP_INVALID";
      const std::string cur = label.substr( 0 , pos );
      if ( i == 0 )
        prefix = cur;
      else if ( cur != prefix )
        return "PP_INVALID";
    }

  return make_invalid_channel_label( prefix );
}

std::string infer_invalid_channel_label( const std::string & signal_label )
{
  const std::string::size_type pos = signal_label.rfind( '_' );
  if ( pos == std::string::npos || pos == 0 ) return "";
  return make_invalid_channel_label( signal_label.substr( 0 , pos ) );
}

}


// ============================================================
// Internal types
// ============================================================

struct hd_params_t
{
  std::vector<std::string> ch;      // channel names [W, N1, N2, N3, R]
  bool do_3state;
  bool do_hd_metrics;
  double window_sec;                // transition window half-width (seconds)
  double lag_sec;                   // lag for TV_lag (seconds)
  double stable_min_sec;            // min seconds of contiguous flank stability for STABLE events
  double conf_th;                   // confidence threshold for frac_below metric
  double hd_smooth_sec;             // smoothing window for soft sleep/REM trajectories
  double hd_onset_win_min;          // forward integrated-area window for onset detection
  double hd_sleep_mass_th;          // mean sleep probability threshold for sleep onset
  double hd_rem_mass_th;            // mean REM probability threshold for REM onset
  double hd_offset_win_min;         // backward integrated-area window for final-sleep detection
  double hd_offset_mass_th;         // mean sleep probability threshold for final-sleep detection
  std::string annot_name;           // annotation class for stratification (empty = none)
  bool verbose;                     // emit per-sample HDSIG table
  int min_events_profile;           // minimum events needed to emit aligned profile
  int min_samples_region;           // minimum samples needed to emit REGION summaries
  std::string emit_annot_label;     // if non-empty, emit transition annotations under this class
  double stable_ps_th;              // per-sample PS must be < this to qualify as stable
  double stable_conf_th;            // per-sample C must be > this to qualify as stable
  double min_shift;                 // minimum L1 posterior shift for TR event acceptance (0=disabled)
  double shift_win_sec;             // pre/post window (seconds) for shift computation
  bool emit_annot_typed_classes;    // if true, encode directional transition type in point-event class name
  double bout_min_sec;              // minimum stable-run span (seconds) to qualify as a bout
  double bout_gap_sec;               // max within-run gap (seconds) tolerated before splitting a bout
  double bout_smooth_sec;           // centered moving-average half-window (seconds) applied to PS/CONF before the bout stability test (0=off)
  std::string ch_valid;             // channel name for PP_INVALID invalidity flag (optional)
  std::string emit_sig_label;       // if non-empty, emit PS/H/CONF per-sample signals under this channel-name prefix
  bool do_emit_hypnogram;           // if true, emit hard-called (argmax) stage run annotations
  std::string hypno_prefix;         // optional prefix for those annotation classes, e.g. "hd_"
};


// ============================================================
// Layer 1: data loading and validation
// ============================================================

struct hd_data_t
{
  int N;
  double Fs;
  std::vector<std::vector<double>> p;   // p[5][N], stages: W N1 N2 N3 R
  std::vector<std::vector<double>> p3;  // p3[3][N], stages: W NREM R
  std::vector<uint64_t> tp;             // timepoints
  std::vector<bool> valid;              // per-sample validity (false = NA/ZOH-filled)
  int n_invalid_lead  = 0;             // consecutive invalid samples at start
  int n_invalid_trail = 0;             // consecutive invalid samples at end
  int n_invalid_mid   = 0;             // invalid samples in the interior

  bool load( edf_t & edf, const hd_params_t & par )
  {
    const int K = 5;

    // Resolve channels and verify all exist
    std::vector<int> sidx( K );
    for ( int k = 0; k < K; k++ )
      {
	signal_list_t sl = edf.header.signal_list( par.ch[k] );
	if ( sl.size() == 0 )
	  Helper::halt( "channel not found: " + par.ch[k] );
	if ( edf.header.is_annotation_channel( sl(0) ) )
	  Helper::halt( par.ch[k] + " is an annotation channel, not a signal" );
	sidx[k] = sl(0);
      }

    // Verify same sample rate for all channels
    Fs = edf.header.sampling_freq( sidx[0] );
    for ( int k = 1; k < K; k++ )
      {
	double fk = edf.header.sampling_freq( sidx[k] );
	if ( std::fabs( fk - Fs ) > 1e-6 )
	  Helper::halt( "all hypnodensity channels must share the same sample rate; "
			+ par.ch[k] + " has " + Helper::dbl2str(fk) + " Hz vs "
			+ Helper::dbl2str(Fs) + " Hz for " + par.ch[0] );
      }

    // Read the whole trace for all channels
    interval_t whole = edf.timeline.wholetrace();

    p.resize( K );
    for ( int k = 0; k < K; k++ )
      {
	slice_t sl( edf, sidx[k], whole );
	const std::vector<double> * d = sl.pdata();
	p[k] = *d;
	if ( k == 0 )
	  {
	    N = (int)p[0].size();
	    const std::vector<uint64_t> * tps = sl.ptimepoints();
	    tp = *tps;
	  }
	else if ( (int)p[k].size() != N )
	  Helper::halt( "channel length mismatch for " + par.ch[k] );
      }

    if ( N == 0 )
      Helper::halt( "no samples found" );

    // Load optional PP_INVALID channel to seed the validity mask.
    // Any sample with PP_INVALID >= 0.5 (i.e. flagged invalid) is masked out.
    // If channel is absent all samples start as valid.
    valid.assign( N, true );
    {
      signal_list_t sl_v = edf.header.signal_list( par.ch_valid );
      if ( sl_v.size() > 0 && ! edf.header.is_annotation_channel( sl_v(0) ) )
        {
          slice_t sv( edf, sl_v(0), whole );
          const std::vector<double> * dv = sv.pdata();
          if ( (int)dv->size() == N )
            {
              for ( int i = 0; i < N; i++ )
                valid[i] = (*dv)[i] < 0.5;   // 0=genuine, 1=ZOH-filled/invalid
              logger << "  loaded invalidity channel '" << par.ch_valid << "'\n";
            }
          else
            logger << "  invalidity channel '" << par.ch_valid << "' length mismatch — ignoring\n";
        }
    }

    // Validate and optionally renormalize row sums.
    // Samples already marked invalid (PP_INVALID=1) are zeroed out and skipped.
    // Samples marked valid but with implausible row sums are warned and marked invalid.
    int n_renorm = 0, n_bad_rowsum = 0;
    for ( int i = 0; i < N; i++ )
      {
        if ( ! valid[i] )
          {
            for ( int k = 0; k < K; k++ ) p[k][i] = 0.0;
            continue;
          }

        bool sample_nonfinite = false;
	double rowsum = 0.0;
	for ( int k = 0; k < K; k++ )
          {
            if ( ! std::isfinite( p[k][i] ) ) sample_nonfinite = true;
            rowsum += p[k][i];
          }

	if ( sample_nonfinite || ! std::isfinite( rowsum ) || rowsum < 0.5 || rowsum > 1.5 )
	  {
	    logger << "  warning: invalid posterior row at sample " << i
		   << " (rowsum=" << rowsum << ") — marking invalid\n";
	    valid[i] = false;
	    ++n_bad_rowsum;
	    for ( int k = 0; k < K; k++ ) p[k][i] = 0.0;
	    continue;
	  }

	if ( std::fabs( rowsum - 1.0 ) > 1e-4 )
	  {
	    ++n_renorm;
	    for ( int k = 0; k < K; k++ ) p[k][i] /= rowsum;
	  }

	// Clamp to [0,1]
	for ( int k = 0; k < K; k++ )
	  {
	    if ( p[k][i] < 0.0 ) p[k][i] = 0.0;
	    if ( p[k][i] > 1.0 ) p[k][i] = 1.0;
	  }
      }

    if ( n_renorm > 0 )
      logger << "  renormalized " << n_renorm << " of " << N << " samples\n";
    if ( n_bad_rowsum > 0 )
      logger << "  " << n_bad_rowsum << " samples had implausible row sums and were marked invalid\n";

    // Classify invalid samples as leading / trailing / middle
    {
      int lo = 0;
      while ( lo < N && ! valid[lo] ) ++lo;
      n_invalid_lead = lo;

      int hi = N - 1;
      while ( hi >= lo && ! valid[hi] ) --hi;
      n_invalid_trail = ( N - 1 ) - hi;

      for ( int i = lo; i <= hi; i++ )
        if ( ! valid[i] ) ++n_invalid_mid;
    }

    // Build 3-state collapsed posteriors: W, NREM=N1+N2+N3, R
    if ( par.do_3state )
      {
	p3.resize( 3, std::vector<double>( N, 0.0 ) );
	for ( int i = 0; i < N; i++ )
	  {
	    p3[0][i] = p[0][i];               // W
	    p3[1][i] = p[1][i] + p[2][i] + p[3][i];  // NREM
	    p3[2][i] = p[4][i];               // R
	  }
      }

    return true;
  }
};


// ============================================================
// Layer 2: per-sample derived signals
// ============================================================

struct hd_derived_t
{
  // 5-state signals
  std::vector<double> H, C, Mg;          // entropy, confidence, margin
  std::vector<double> TV, TV_lag;        // motion metrics
  std::vector<int>    argmax5;

  // 3-state signals (only populated if do_3state)
  std::vector<double> H3, C3, Mg3;
  std::vector<double> TV3, TV3_lag;
  std::vector<int>    argmax3;

  void compute( const hd_data_t & dat, const hd_params_t & par )
  {
    const int N = dat.N;
    const int lag = std::max( 1, (int)std::round( par.lag_sec * dat.Fs ) );

    H.assign( N, 0.0 );
    C.assign( N, 0.0 );
    Mg.assign( N, 0.0 );
    TV.assign( N, 0.0 );
    TV_lag.assign( N, 0.0 );
    argmax5.assign( N, 0 );

    for ( int i = 0; i < N; i++ )
      {
	// Entropy H(t) = -sum p_k log(p_k)
	double h = 0.0;
	for ( int k = 0; k < 5; k++ )
	  {
	    double pk = dat.p[k][i];
	    if ( pk > 0.0 ) h -= pk * std::log( pk );
	  }
	H[i] = h;

	// Find max and second-max for confidence and margin
	double mx1 = -1.0, mx2 = -1.0;
	int amax = 0;
	for ( int k = 0; k < 5; k++ )
	  {
	    if ( dat.p[k][i] > mx1 ) { mx2 = mx1; mx1 = dat.p[k][i]; amax = k; }
	    else if ( dat.p[k][i] > mx2 ) { mx2 = dat.p[k][i]; }
	  }
	C[i]      = mx1;
	Mg[i]     = ( mx2 >= 0.0 ) ? mx1 - mx2 : mx1;
	argmax5[i] = amax;

	// Motion metrics (step from previous sample)
	if ( i > 0 )
	  {
	    double tv = 0.0;
	    for ( int k = 0; k < 5; k++ )
	      tv += std::fabs( dat.p[k][i] - dat.p[k][i-1] );
	    TV[i] = 0.5 * tv;
	  }

	// Longer-lag motion
	if ( i >= lag )
	  {
	    double tv = 0.0;
	    for ( int k = 0; k < 5; k++ )
	      tv += std::fabs( dat.p[k][i] - dat.p[k][i-lag] );
	    TV_lag[i] = 0.5 * tv;
	  }
      }

    if ( ! par.do_3state ) return;

    // 3-state derived signals
    const int lag3 = lag;  // same lag in samples
    H3.assign( N, 0.0 );
    C3.assign( N, 0.0 );
    Mg3.assign( N, 0.0 );
    TV3.assign( N, 0.0 );
    TV3_lag.assign( N, 0.0 );
    argmax3.assign( N, 0 );

    for ( int i = 0; i < N; i++ )
      {
	double h = 0.0;
	for ( int k = 0; k < 3; k++ )
	  {
	    double pk = dat.p3[k][i];
	    if ( pk > 0.0 ) h -= pk * std::log( pk );
	  }
	H3[i] = h;

	double mx1 = -1.0, mx2 = -1.0;
	int amax = 0;
	for ( int k = 0; k < 3; k++ )
	  {
	    if ( dat.p3[k][i] > mx1 ) { mx2 = mx1; mx1 = dat.p3[k][i]; amax = k; }
	    else if ( dat.p3[k][i] > mx2 ) { mx2 = dat.p3[k][i]; }
	  }
	C3[i]      = mx1;
	Mg3[i]     = ( mx2 >= 0.0 ) ? mx1 - mx2 : mx1;
	argmax3[i]  = amax;

	if ( i > 0 )
	  {
	    double tv = 0.0;
	    for ( int k = 0; k < 3; k++ )
	      tv += std::fabs( dat.p3[k][i] - dat.p3[k][i-1] );
	    TV3[i] = 0.5 * tv;
	  }

	if ( i >= lag3 )
	  {
	    double tv = 0.0;
	    for ( int k = 0; k < 3; k++ )
	      tv += std::fabs( dat.p3[k][i] - dat.p3[k][i-lag3] );
	    TV3_lag[i] = 0.5 * tv;
	  }
      }
  }
};


// ============================================================
// Layer 3: transition detection
// ============================================================

struct hd_event_t
{
  int    idx;         // boundary between samples idx-1 and idx
  int    from_st;     // dominant left state
  int    to_st;       // dominant right state
  double shift;       // 0.5 * L1(left,right)
  double peak_tv;     // local TV summary around boundary
  bool   is_dir;      // definite X->Y transition
  bool   is_stable;   // stable-X to stable-Y transition
  hd_event_t() : idx(0), from_st(-1), to_st(-1), shift(0.0), peak_tv(0.0), is_dir(false), is_stable(false) {}
};

static std::vector<double> make_linear_weights( const int n )
{
  std::vector<double> w( n , 1.0 );
  if ( n <= 1 ) return w;
  for ( int i = 0 ; i < n ; i++ )
    w[i] = ( i + 1 ) / (double)n;
  return w;
}

static std::vector<double> weighted_support(
    const std::vector<std::vector<double>> & P,
    int lo,
    int hi,
    bool reverse_weights )
{
  const int K = P.size();
  std::vector<double> out( K , 0.0 );
  if ( hi < lo ) return out;
  const int n = hi - lo + 1;
  const std::vector<double> w = make_linear_weights( n );
  double wsum = 0.0;
  for ( int j = 0 ; j < n ; j++ )
    {
      const int si = lo + j;
      const double ww = reverse_weights ? w[n - 1 - j] : w[j];
      wsum += ww;
      for ( int k = 0 ; k < K ; k++ )
        out[k] += ww * P[k][si];
    }
  if ( wsum > 0.0 )
    for ( int k = 0 ; k < K ; k++ ) out[k] /= wsum;
  return out;
}

static int dominant_state( const std::vector<double> & p )
{
  int amax = -1;
  double mx = -1.0;
  for ( int k = 0 ; k < (int)p.size() ; k++ )
    if ( p[k] > mx ) { mx = p[k]; amax = k; }
  return amax;
}

static double dominant_margin( const std::vector<double> & p , const int amax )
{
  if ( amax < 0 || amax >= (int)p.size() ) return 0.0;
  double mx2 = -1.0;
  for ( int k = 0 ; k < (int)p.size() ; k++ )
    if ( k != amax && p[k] > mx2 ) mx2 = p[k];
  return p[amax] - std::max( 0.0 , mx2 );
}

// Centered moving average, used only to smooth PS/CONF for the bout-formation
// stability test -- so a single noisy sample doesn't reset the whole
// stable-run clock the way it does for the (unsmoothed) IS_STABLE/TRANS path.
static std::vector<double> smooth_series( const std::vector<double> & x, int half_win )
{
  const int N = (int)x.size();
  if ( half_win <= 0 || N == 0 ) return x;
  std::vector<double> prefix( N + 1, 0.0 );
  for ( int i = 0; i < N; i++ ) prefix[i+1] = prefix[i] + x[i];
  std::vector<double> out( N );
  for ( int i = 0; i < N; i++ )
    {
      const int lo = std::max( 0,     i - half_win );
      const int hi = std::min( N - 1, i + half_win );
      out[i] = ( prefix[hi+1] - prefix[lo] ) / (double)( hi - lo + 1 );
    }
  return out;
}

// STABLE mask: contiguous runs where TV < tv_th AND C > conf_th for >= min_run_samples
static void build_stable_mask(
    const std::vector<double> & tv,
    const std::vector<double> & conf,
    int N,
    int min_run_samples,
    double tv_th,
    double conf_th,
    std::vector<bool> & is_stable )
{
  is_stable.assign( N, false );
  std::vector<bool> q( N, false );
  for ( int i = 0; i < N; i++ )
    q[i] = tv[i] < tv_th && conf[i] > conf_th;
  int i = 0;
  while ( i < N )
    {
      if ( !q[i] ) { ++i; continue; }
      int lo = i;
      while ( i < N && q[i] ) ++i;
      if ( i - lo >= min_run_samples )
        for ( int j = lo; j < i; j++ ) is_stable[j] = true;
    }
}

// A bout: a long, settled, single-dominant-stage run (see bout-min/bout-gap).
struct hd_bout_t { int lo, hi, dom; };

// Find bouts within context: contiguous is_stable runs, merged across gaps of
// <= bout_gap_sec, kept only if the merged span is >= bout_min_sec.
static std::vector<hd_bout_t> compute_bouts(
    const std::vector<bool> & context_mask,
    const std::vector<bool> & is_stable,
    const std::vector<int> & argmax,
    const double bout_min_sec,
    const double bout_gap_sec,
    const double Fs )
{
  const int bout_min_smp = std::max( 1, (int)std::round( bout_min_sec * Fs ) );
  const int bout_gap_smp = std::max( 0, (int)std::round( bout_gap_sec * Fs ) );

  // Step 1: find raw stable runs within context
  std::vector<hd_bout_t> raw;
  const int Ns = (int)context_mask.size();
  int si = 0;
  while ( si < Ns )
    {
      if ( !context_mask[si] || !is_stable[si] ) { ++si; continue; }
      int lo = si;
      while ( si < Ns && context_mask[si] && is_stable[si] ) ++si;
      std::map<int,int> scnt;
      for ( int j = lo; j < si; j++ ) ++scnt[ argmax[j] ];
      int dom = -1, best = -1;
      for ( const auto & kv : scnt )
        if ( kv.second > best ) { best = kv.second; dom = kv.first; }
      hd_bout_t r; r.lo = lo; r.hi = si - 1; r.dom = dom;
      raw.push_back( r );
    }

  // Step 2: merge adjacent same-stage runs separated by <= bout_gap_smp
  std::vector<hd_bout_t> merged;
  for ( const auto & r : raw )
    {
      if ( !merged.empty() &&
           r.dom == merged.back().dom &&
           r.lo - merged.back().hi - 1 <= bout_gap_smp )
        merged.back().hi = r.hi;   // extend span
      else
        merged.push_back( r );
    }

  // Step 3: require span >= bout_min_sec to qualify as a bout
  std::vector<hd_bout_t> bouts;
  for ( const auto & r : merged )
    if ( r.hi - r.lo + 1 >= bout_min_smp )
      bouts.push_back( r );

  return bouts;
}

// Bout-to-bout transitions: adjacent bouts with different dominant stages.
// The gap between them is always non-empty (see build_stable_mask: any state
// change forces at least one non-stable sample immediately after it, so
// adjacent bouts never directly touch), but can otherwise be anywhere from a
// single sample to many minutes long -- any intervening run too short to
// itself qualify as a bout is simply absorbed into the gap. Anchor at the
// highest-shift TRANS event already detected within the gap (the same
// windowed, pre/post-averaged shift measure TRANS itself uses to place its
// own events) rather than raw per-sample TV: TV is a single-step (lag-1)
// posterior derivative, noisy enough that its peak can land on an unrelated
// local blip somewhere in a long, messy gap instead of the actual crossing.
// Falls back to the peak-TV sample only for the rare gap with no TRANS event
// in it at all. Reuses hd_event_t so the existing profile-shape machinery
// needs no duplication.
static std::vector<hd_event_t> compute_bout_transitions( const std::vector<hd_bout_t> & bouts,
                                                           const std::vector<hd_event_t> & events,
                                                           const std::vector<double> & tv )
{
  std::vector<hd_event_t> out;
  for ( int ri = 1; ri < (int)bouts.size(); ri++ )
    if ( bouts[ri].dom != bouts[ri-1].dom )
      {
        const int gap_lo = bouts[ri-1].hi + 1;
        const int gap_hi = bouts[ri].lo   - 1;

        int peak = -1;
        double best_shift = -1.0;
        for ( const auto & te : events )
          if ( te.idx >= gap_lo && te.idx <= gap_hi && te.shift > best_shift )
            {
              best_shift = te.shift;
              peak = te.idx;
            }

        if ( peak < 0 )
          {
            peak = gap_lo;
            for ( int i = gap_lo + 1; i <= gap_hi; i++ )
              if ( tv[i] > tv[peak] ) peak = i;
          }

        hd_event_t e;
        e.idx      = peak;
        e.from_st  = bouts[ri-1].dom;
        e.to_st    = bouts[ri].dom;
        e.is_dir   = true;
        e.is_stable = true;
        out.push_back( e );
      }
  return out;
}

// Count of TRANS events whose index falls within a bout-to-bout gap with a
// matching stage pair — i.e. detected crossings corroborated by genuine
// settled-block-to-settled-block movement, as opposed to crossings inside
// stretches that never consolidate on one or both sides. Speaks to fragmentation:
// a low match rate means most detected switching isn't anchored to real bouts.
static int count_trans_bout_matches(
    const std::vector<hd_event_t> & trans_events,
    const std::vector<hd_bout_t> & bouts )
{
  int n_matched = 0;
  for ( const auto & e : trans_events )
    for ( int ri = 1; ri < (int)bouts.size(); ri++ )
      {
        if ( bouts[ri].dom == bouts[ri-1].dom ) continue;
        if ( e.from_st != bouts[ri-1].dom || e.to_st != bouts[ri].dom ) continue;
        if ( e.idx >= bouts[ri-1].hi && e.idx <= bouts[ri].lo )
          {
            ++n_matched;
            break;
          }
      }
  return n_matched;
}

struct hd_trans_t
{
  std::vector<hd_event_t> events5;
  std::vector<hd_event_t> bout_events5;
  std::vector<bool>       is_trans5;
  std::vector<bool>       is_bout_trans5;
  std::vector<bool>       is_stable5;
  std::vector<bool>       is_bout_stable5;   // smoothed PS/CONF stability test, bout-formation only

  std::vector<hd_event_t> events3;
  std::vector<hd_event_t> bout_events3;
  std::vector<bool>       is_trans3;
  std::vector<bool>       is_bout_trans3;
  std::vector<bool>       is_stable3;
  std::vector<bool>       is_bout_stable3;   // smoothed PS/CONF stability test, bout-formation only

  static std::vector<hd_event_t> detect_boundary_events(
      const std::vector<std::vector<double>> & P,
      const std::vector<double> & tv,
      const std::vector<bool> & is_stable,
      const std::vector<int> & argmax,
      const hd_params_t & par,
      const double Fs )
  {
    const int N = argmax.size();
    const int stable_smp = std::max( 1, (int)std::round( par.stable_min_sec * Fs ) );
    const int shift_smp  = std::max( 1, (int)std::round( par.shift_win_sec  * Fs ) );
    const int win_smp    = std::max( 1, (int)std::round( par.window_sec     * Fs ) );
    std::vector<hd_event_t> candidates;

    const double dir_margin_th = std::max( 0.0 , 0.5 * par.min_shift );

    for ( int i = 1 ; i < N ; i++ )
      {
        const int pre_lo  = std::max( 0, i - shift_smp );
        const int pre_hi  = i - 1;
        const int post_lo = i;
        const int post_hi = std::min( N - 1, i + shift_smp - 1 );
        if ( pre_hi < pre_lo || post_hi < post_lo ) continue;

        const std::vector<double> left  = weighted_support( P, pre_lo, pre_hi, true );
        const std::vector<double> right = weighted_support( P, post_lo, post_hi, false );
        const int left_st  = dominant_state( left );
        const int right_st = dominant_state( right );

        double l1 = 0.0;
        for ( int k = 0 ; k < (int)P.size() ; k++ )
          l1 += std::fabs( right[k] - left[k] );
        const double shift = 0.5 * l1;

        if ( shift < par.min_shift ) continue;

        hd_event_t e;
        e.idx = i;
        e.shift = shift;
        e.peak_tv = std::max( tv[ std::max( 0, i - 1 ) ] , tv[ std::min( N - 1, i ) ] );

        const double left_margin  = dominant_margin( left, left_st );
        const double right_margin = dominant_margin( right, right_st );

        if ( left_st >= 0 && right_st >= 0 &&
             left_st != right_st &&
             left_margin >= dir_margin_th &&
             right_margin >= dir_margin_th )
          {
            e.from_st = left_st;
            e.to_st = right_st;
            e.is_dir = true;
          }

        // Only directional X->Y boundaries count as transition events.
        if ( ! e.is_dir ) continue;

        if ( e.is_dir &&
             i - stable_smp >= 0 &&
             i + stable_smp - 1 < N )
          {
            bool left_ok = true, right_ok = true;
            for ( int j = i - stable_smp ; j < i ; j++ )
              if ( ! is_stable[j] || argmax[j] != e.from_st ) { left_ok = false; break; }
            for ( int j = i ; j < i + stable_smp ; j++ )
              if ( ! is_stable[j] || argmax[j] != e.to_st ) { right_ok = false; break; }
            e.is_stable = left_ok && right_ok;
          }

        candidates.push_back( e );
      }

    if ( candidates.empty() ) return candidates;

    // Require each accepted event to be a salient local maximum in left/right shift
    // among nearby candidates with the same directional label.
    std::vector<hd_event_t> peak_candidates;
    for ( int ci = 0 ; ci < (int)candidates.size() ; ci++ )
      {
        const hd_event_t & cur = candidates[ci];
        bool is_peak = true;
        for ( int cj = 0 ; cj < (int)candidates.size() ; cj++ )
          {
            if ( ci == cj ) continue;
            const hd_event_t & oth = candidates[cj];
            if ( std::abs( oth.idx - cur.idx ) > shift_smp ) continue;
            if ( oth.from_st != cur.from_st || oth.to_st != cur.to_st ) continue;
            if ( oth.shift > cur.shift )
              {
                is_peak = false;
                break;
              }
          }
        if ( is_peak ) peak_candidates.push_back( cur );
      }

    if ( peak_candidates.empty() ) return peak_candidates;

    const int min_gap = std::max( shift_smp , win_smp );
    std::vector<hd_event_t> events;
    hd_event_t best = peak_candidates[0];

    for ( int ci = 1 ; ci < (int)peak_candidates.size() ; ci++ )
      {
        const hd_event_t & cur = peak_candidates[ci];
        if ( cur.idx - best.idx < min_gap &&
             cur.from_st == best.from_st &&
             cur.to_st == best.to_st )
          {
            if ( cur.shift > best.shift ||
                 ( cur.shift == best.shift && cur.peak_tv > best.peak_tv ) )
              best = cur;
          }
        else
          {
            events.push_back( best );
            best = cur;
          }
      }

    events.push_back( best );
    return events;
  }

  static void build_trans_mask(
      const std::vector<hd_event_t> & events,
      int N,
      int window_samples,
      std::vector<bool> & is_trans )
  {
    is_trans.assign( N, false );
    for ( const auto & e : events )
      {
        int lo = std::max( 0,     e.idx - window_samples );
        int hi = std::min( N - 1, e.idx + window_samples );
        for ( int i = lo; i <= hi; i++ ) is_trans[i] = true;
      }
  }

  void detect( const hd_data_t & dat,
	       const hd_derived_t & der,
	       const hd_params_t & par )
  {
    const int N              = dat.N;
    const int win_smp        = std::max( 1, (int)std::round( par.window_sec       * dat.Fs ) );
    const int stable_smp     = std::max( 1, (int)std::round( par.stable_min_sec   * dat.Fs ) );
    const int bout_smooth_smp = std::max( 0, (int)std::round( par.bout_smooth_sec * dat.Fs ) );

    // 5-state detection
    {
      build_stable_mask( der.TV, der.C, N, stable_smp,
                         par.stable_ps_th, par.stable_conf_th, is_stable5 );
      events5 = detect_boundary_events( dat.p, der.TV, is_stable5, der.argmax5, par, dat.Fs );
      build_trans_mask( events5, N, win_smp, is_trans5 );

      // TRANS takes priority: clear stable for any peri-event sample
      for ( int i = 0; i < N; i++ )
        if ( is_trans5[i] ) is_stable5[i] = false;

      // Bouts are formed on a separately smoothed PS/CONF stability test --
      // scoped to bout-formation only, so a single noisy sample doesn't reset
      // the whole stable-run clock the way it still does for IS_STABLE/TRANS.
      const std::vector<double> ps_smooth5   = smooth_series( der.TV, bout_smooth_smp );
      const std::vector<double> conf_smooth5 = smooth_series( der.C,  bout_smooth_smp );
      build_stable_mask( ps_smooth5, conf_smooth5, N, stable_smp,
                         par.stable_ps_th, par.stable_conf_th, is_bout_stable5 );
      for ( int i = 0; i < N; i++ )
        if ( is_trans5[i] ) is_bout_stable5[i] = false;

      // BOUT_TRANS: per-sample halo (±window_sec), same construction as
      // IS_TRANS but anchored on bout-to-bout transitions instead of events.
      const std::vector<bool> all_mask( N, true );
      const std::vector<hd_bout_t> bouts5 = compute_bouts( all_mask, is_bout_stable5, der.argmax5,
                                                             par.bout_min_sec, par.bout_gap_sec, dat.Fs );
      bout_events5 = compute_bout_transitions( bouts5, events5, der.TV );
      build_trans_mask( bout_events5, N, win_smp, is_bout_trans5 );
    }

    if ( ! par.do_3state ) return;

    // 3-state detection (independent from 5-state)
    {
      build_stable_mask( der.TV3, der.C3, N, stable_smp,
                         par.stable_ps_th, par.stable_conf_th, is_stable3 );
      events3 = detect_boundary_events( dat.p3, der.TV3, is_stable3, der.argmax3, par, dat.Fs );
      build_trans_mask( events3, N, win_smp, is_trans3 );

      // TRANS takes priority: clear stable for any peri-event sample
      for ( int i = 0; i < N; i++ )
        if ( is_trans3[i] ) is_stable3[i] = false;

      const std::vector<double> ps_smooth3   = smooth_series( der.TV3, bout_smooth_smp );
      const std::vector<double> conf_smooth3 = smooth_series( der.C3,  bout_smooth_smp );
      build_stable_mask( ps_smooth3, conf_smooth3, N, stable_smp,
                         par.stable_ps_th, par.stable_conf_th, is_bout_stable3 );
      for ( int i = 0; i < N; i++ )
        if ( is_trans3[i] ) is_bout_stable3[i] = false;

      const std::vector<bool> all_mask3( N, true );
      const std::vector<hd_bout_t> bouts3 = compute_bouts( all_mask3, is_bout_stable3, der.argmax3,
                                                             par.bout_min_sec, par.bout_gap_sec, dat.Fs );
      bout_events3 = compute_bout_transitions( bouts3, events3, der.TV3 );
      build_trans_mask( bout_events3, N, win_smp, is_bout_trans3 );
    }
  }
};


// ============================================================
// Layer 4: summary engine
// ============================================================

// Percentile from a vector (value in [0,1])
static double percentile( std::vector<double> v, double p )
{
  if ( v.empty() ) return std::numeric_limits<double>::quiet_NaN();
  std::sort( v.begin(), v.end() );
  double idx = p * ( (double)v.size() - 1.0 );
  int lo = (int)idx;
  int hi = lo + 1;
  if ( hi >= (int)v.size() ) return v.back();
  double frac = idx - lo;
  return v[lo] * (1.0 - frac) + v[hi] * frac;
}

// Pearson correlation between two same-length vectors
static double pearson( const std::vector<double> & x,
		       const std::vector<double> & y )
{
  if ( x.size() != y.size() || x.empty() ) return std::numeric_limits<double>::quiet_NaN();
  int n = (int)x.size();
  double mx = 0.0, my = 0.0;
  for ( int i = 0; i < n; i++ ) { mx += x[i]; my += y[i]; }
  mx /= n; my /= n;
  double num = 0.0, dx2 = 0.0, dy2 = 0.0;
  for ( int i = 0; i < n; i++ )
    {
      double dx = x[i] - mx, dy = y[i] - my;
      num += dx * dy;
      dx2 += dx * dx;
      dy2 += dy * dy;
    }
  double denom = std::sqrt( dx2 * dy2 );
  return ( denom > 0.0 ) ? num / denom : std::numeric_limits<double>::quiet_NaN();
}

static std::vector<double> moving_average_centered( const std::vector<double> & x , int width )
{
  const int N = x.size();
  std::vector<double> out( N , 0.0 );
  if ( N == 0 ) return out;

  width = std::max( 1 , width );
  const int half = width / 2;

  std::vector<double> prefix( N + 1 , 0.0 );
  for ( int i = 0 ; i < N ; i++ ) prefix[i+1] = prefix[i] + x[i];

  for ( int i = 0 ; i < N ; i++ )
    {
      const int lo = std::max( 0 , i - half );
      const int hi = std::min( N - 1 , i + half );
      out[i] = ( prefix[hi+1] - prefix[lo] ) / (double)( hi - lo + 1 );
    }

  return out;
}

static std::vector<double> window_mean_forward( const std::vector<double> & x , int width )
{
  const int N = x.size();
  std::vector<double> out( N , 0.0 );
  if ( N == 0 ) return out;

  width = std::max( 1 , width );
  std::vector<double> prefix( N + 1 , 0.0 );
  for ( int i = 0 ; i < N ; i++ ) prefix[i+1] = prefix[i] + x[i];

  for ( int i = 0 ; i < N ; i++ )
    {
      const int hi = std::min( N , i + width );
      out[i] = ( prefix[hi] - prefix[i] ) / (double)( hi - i );
    }

  return out;
}

static std::vector<double> window_mean_backward( const std::vector<double> & x , int width )
{
  const int N = x.size();
  std::vector<double> out( N , 0.0 );
  if ( N == 0 ) return out;

  width = std::max( 1 , width );
  std::vector<double> prefix( N + 1 , 0.0 );
  for ( int i = 0 ; i < N ; i++ ) prefix[i+1] = prefix[i] + x[i];

  for ( int i = 0 ; i < N ; i++ )
    {
      const int lo = std::max( 0 , i - width + 1 );
      out[i] = ( prefix[i+1] - prefix[lo] ) / (double)( i - lo + 1 );
    }

  return out;
}

struct hd_hypno_metrics_t
{
  bool valid_onset = false;
  bool valid_offset = false;
  bool valid_window = false;
  bool valid_rem = false;

  int onset_idx = -1;
  int offset_idx = -1;
  int rem_idx = -1;

  double total_min = 0.0;
  double dt_min = 0.0;
  double hd_tst = 0.0;
  double hd_spt = 0.0;
  double hd_waso = 0.0;
  double hd_sme = 0.0;
  double hd_sol = 0.0;
  double hd_rem_lat = 0.0;
  double hd_rem_lat2 = 0.0;

  std::map<std::string,double> mins;

  void compute( const hd_data_t & dat , const hd_params_t & par )
  {
    mins.clear();
    const int N = dat.N;
    if ( N == 0 ) return;

    dt_min = 1.0 / dat.Fs / 60.0;
    total_min = N * dt_min;

    std::vector<double> pW( N , 0.0 ) , pR( N , 0.0 ) , pNR( N , 0.0 ) , pS( N , 0.0 );
    for ( int i = 0 ; i < N ; i++ )
      {
        pW[i] = dat.p[0][i];
        pR[i] = dat.p[4][i];
        pNR[i] = dat.p[1][i] + dat.p[2][i] + dat.p[3][i];
        pS[i] = pNR[i] + pR[i];

        mins["W"]  += dat.p[0][i] * dt_min;
        mins["N1"] += dat.p[1][i] * dt_min;
        mins["N2"] += dat.p[2][i] * dt_min;
        mins["N3"] += dat.p[3][i] * dt_min;
        mins["R"]  += dat.p[4][i] * dt_min;
        mins["NR"] += pNR[i] * dt_min;
        mins["S"]  += pS[i] * dt_min;
      }

    hd_tst = mins["S"];

    const int smooth_smp = std::max( 1 , (int)std::round( par.hd_smooth_sec * dat.Fs ) );
    const int onset_smp  = std::max( 1 , (int)std::round( par.hd_onset_win_min * 60.0 * dat.Fs ) );
    const int offset_smp = std::max( 1 , (int)std::round( par.hd_offset_win_min * 60.0 * dat.Fs ) );

    const std::vector<double> pS_sm = moving_average_centered( pS , smooth_smp );
    const std::vector<double> pR_sm = moving_average_centered( pR , smooth_smp );
    const std::vector<double> a_sleep = window_mean_forward( pS_sm , onset_smp );
    const std::vector<double> a_rem   = window_mean_forward( pR_sm , onset_smp );
    const std::vector<double> b_sleep = window_mean_backward( pS_sm , offset_smp );

    for ( int i = 0 ; i < N ; i++ )
      if ( a_sleep[i] >= par.hd_sleep_mass_th )
        {
          onset_idx = i;
          valid_onset = true;
          hd_sol = i * dt_min;
          break;
        }

    if ( ! valid_onset ) return;

    for ( int i = N - 1 ; i >= onset_idx ; i-- )
      if ( b_sleep[i] >= par.hd_offset_mass_th )
        {
          offset_idx = i;
          valid_offset = true;
          break;
        }

    if ( valid_offset && offset_idx > onset_idx )
      {
        valid_window = true;
        hd_spt = ( offset_idx - onset_idx ) * dt_min;

        double sleep_mass_window = 0.0;
        for ( int i = onset_idx ; i <= offset_idx ; i++ )
          {
            hd_waso += pW[i] * dt_min;
            sleep_mass_window += pS[i] * dt_min;
          }
        hd_sme = hd_spt > 0 ? 100.0 * sleep_mass_window / hd_spt : 0.0;
      }

    for ( int i = onset_idx ; i < N ; i++ )
      if ( a_rem[i] >= par.hd_rem_mass_th )
        {
          rem_idx = i;
          valid_rem = true;
          hd_rem_lat = ( i - onset_idx ) * dt_min;

          double sleep_mass_to_rem = 0.0;
          for ( int j = onset_idx ; j <= i ; j++ ) sleep_mass_to_rem += pS[j] * dt_min;
          hd_rem_lat2 = sleep_mass_to_rem;
          break;
        }
  }
};

struct hd_region_stats_t
{
  int n = 0;
  bool valid = false;
  // Domain A: mixedness
  double mean_H = 0, p90_H = 0;
  double mean_C = 0, frac_C_below = 0, mean_Mg = 0;
  // Domain B: instability
  double mean_TV = 0, p90_TV = 0, mean_TV_lag = 0;
  double corr_H_TV = 0;
  // Mean posterior probability for each state across included samples.
  std::vector<double> mean_p;
};

static hd_region_stats_t compute_region(
    const hd_data_t & dat,
    const hd_derived_t & der,
    const std::vector<bool> & mask,   // which samples to include
    const hd_params_t & par,
    bool threestate )
{
  hd_region_stats_t rs;

  std::vector<double> vH, vC, vMg, vTV, vTV_lag;
  std::vector<double> sum_p( threestate ? 3 : 5 , 0.0 );
  int n_p = 0;

  const std::vector<double> & H_ref       = threestate ? der.H3      : der.H;
  const std::vector<double> & C_ref       = threestate ? der.C3      : der.C;
  const std::vector<double> & Mg_ref      = threestate ? der.Mg3     : der.Mg;
  const std::vector<double> & TV_ref      = threestate ? der.TV3     : der.TV;
  const std::vector<double> & TVlag_ref   = threestate ? der.TV3_lag : der.TV_lag;
  const std::vector<std::vector<double>> & P_ref = threestate ? dat.p3 : dat.p;
  const int K = threestate ? 3 : 5;

  int N = (int)mask.size();
  for ( int i = 0; i < N; i++ )
    {
      if ( !mask[i] ) continue;
      vH.push_back( H_ref[i] );
      vC.push_back( C_ref[i] );
      vMg.push_back( Mg_ref[i] );
      vTV.push_back( TV_ref[i] );
      vTV_lag.push_back( TVlag_ref[i] );
      for ( int k = 0; k < K; k++ ) sum_p[k] += P_ref[k][i];
      ++n_p;
    }

  rs.n = (int)vH.size();
  if ( rs.n == 0 ) return rs;
  rs.valid = true;

  // Means
  auto vmean = [](const std::vector<double> & v) {
    if (v.empty()) return 0.0;
    return std::accumulate( v.begin(), v.end(), 0.0 ) / v.size();
  };
  rs.mean_H      = vmean( vH );
  rs.p90_H       = percentile( vH, 0.9 );
  rs.mean_C      = vmean( vC );
  rs.mean_Mg     = vmean( vMg );
  rs.mean_TV     = vmean( vTV );
  rs.p90_TV      = percentile( vTV, 0.9 );
  rs.mean_TV_lag = vmean( vTV_lag );
  rs.mean_p.assign( K, std::numeric_limits<double>::quiet_NaN() );
  if ( n_p > 0 )
    for ( int k = 0; k < K; k++ )
      rs.mean_p[k] = sum_p[k] / n_p;

  // Fraction with confidence below threshold
  int n_below = 0;
  for ( double v : vC ) if ( v < par.conf_th ) ++n_below;
  rs.frac_C_below = (double)n_below / rs.n;

  // Entropy-TV correlation
  rs.corr_H_TV = pearson( vH, vTV );

  return rs;
}


// State-set mixing: which combinations to report, and how the density is computed.
//
// Label order follows state index order (W, N1, N2, N3, R / W, NR, R) so labels
// are always emitted the same way regardless of iteration order.
static std::vector<std::vector<int>> hd_mix_combinations( const bool threestate )
{
  if ( threestate )
    return { {0,1}, {0,2}, {1,2}, {0,1,2} };
  return { {0,1}, {0,2}, {0,3}, {0,4}, {1,2}, {1,3}, {1,4}, {2,3}, {2,4}, {3,4},
           {1,2,3}, {0,1,4} };
}

// HD_MINS (under MIX): expected minutes jointly ambiguous across a set of states.
// For each sample, its contribution is the minimum posterior across the set, i.e.
// a fractional/soft occupancy weight, integrated exactly as HD_MINS is elsewhere
// in this file (no threshold — every sample contributes, however small).
static double compute_mix_mins(
    const std::vector<std::vector<double>> & P,
    const std::vector<bool> & mask,
    const std::vector<int> & states,
    const double Fs )
{
  const double dt_min = 1.0 / Fs / 60.0;
  double sum_min = 0.0;
  for ( int i = 0; i < (int)mask.size(); i++ )
    {
      if ( !mask[i] ) continue;
      double mn = 1.0;
      for ( int k : states )
        if ( P[k][i] < mn ) mn = P[k][i];
      sum_min += mn * dt_min;
    }
  return sum_min;
}


// Transition shape statistics and aligned profiles
struct hd_trans_stats_t
{
  int n_events = 0;
  double density  = 0;         // events per hour
  double mean_width_sec = 0;
  double mean_peak_H    = 0;
  double mean_min_C     = 0;
  double mean_TV_area   = 0;
  bool valid = false;
};

struct hd_profile_t
{
  int n_events = 0;
  std::vector<double> offsets;   // seconds relative to event center; for discrete
                                 // hard/motion events, the center is the midpoint
                                 // between samples i-1 and i, not the later sample i
  std::vector<double> H, C, Mg, TV;
  std::vector<std::vector<double>> P;   // mean posterior profiles, one vector per state
  bool valid = false;
};

// Compute transition shape stats and aligned profiles for events restricted to a region mask
static void compute_trans_shape(
    const hd_derived_t & der,
    const std::vector<hd_event_t> & all_events,
    const std::vector<bool> & region_mask,
    const std::vector<bool> & is_stable,
    const hd_data_t & dat,
    const hd_params_t & par,
    double Fs,
    bool threestate,
    bool require_dir,
    bool require_stable,
    hd_trans_stats_t & tshape,
    hd_profile_t & profile )
{
  const int N = (int)region_mask.size();
  const int win_smp = std::max( 1, (int)std::round( par.window_sec * Fs ) );

  // Filter events to those whose center sample falls within the region
  std::vector<const hd_event_t *> local_events;
  for ( const auto & e : all_events )
    if ( e.idx >= 0 && e.idx < N &&
         region_mask[e.idx] &&
         ( !require_dir || e.is_dir ) &&
         ( !require_stable || e.is_stable ) )
      local_events.push_back( &e );

  tshape.n_events = (int)local_events.size();
  tshape.valid    = false;

  // Density: count region samples for duration estimate
  int n_region = 0;
  for ( bool b : region_mask ) if (b) ++n_region;
  double dur_hr = (double)n_region / Fs / 3600.0;
  tshape.density = ( dur_hr > 0 ) ? tshape.n_events / dur_hr : 0.0;

  if ( tshape.n_events == 0 ) return;
  tshape.valid = true;

  const std::vector<double> & H_ref  = threestate ? der.H3  : der.H;
  const std::vector<double> & C_ref  = threestate ? der.C3  : der.C;
  const std::vector<double> & TV_ref = threestate ? der.TV3 : der.TV;
  const std::vector<double> & Mg_ref = threestate ? der.Mg3 : der.Mg;
  const std::vector<std::vector<double>> & P_ref = threestate ? dat.p3 : dat.p;
  const int K = threestate ? 3 : 5;

  auto centered_tv = [&]( int si ) {
    if ( N <= 1 ) return TV_ref[si];
    if ( si <= 0 ) return TV_ref[1];
    if ( si >= N - 1 ) return TV_ref[N - 1];
    return 0.5 * ( TV_ref[si] + TV_ref[si + 1] );
  };

  // Per-event scalars
  double sum_width = 0, sum_peak_H = 0, sum_min_C = 0, sum_TV_area = 0;

  // Profile accumulation: 2*win_smp + 1 offset positions
  int profile_len = 2 * win_smp + 1;
  std::vector<double> sum_H( profile_len, 0.0 ), sum_C( profile_len, 0.0 );
  std::vector<double> sum_Mg( profile_len, 0.0 ), sum_TV( profile_len, 0.0 );
  std::vector<std::vector<double>> sum_P( K , std::vector<double>( profile_len , 0.0 ) );
  std::vector<int>    cnt( profile_len, 0 );

  // Baseline H: stable-core samples only; fall back to full region if none available
  std::vector<double> all_H_region, stable_H_region;
  for ( int i = 0; i < N; i++ )
    if ( region_mask[i] )
      {
        all_H_region.push_back( H_ref[i] );
        if ( is_stable[i] ) stable_H_region.push_back( H_ref[i] );
      }
  const std::vector<double> & baseline_src = stable_H_region.empty() ? all_H_region : stable_H_region;
  double baseline_H = percentile( baseline_src, 0.5 );
  double H_width_th = baseline_H + 0.5 * ( percentile( all_H_region, 0.9 ) - baseline_H );

  int n_events_profile = 0;

  for ( const auto * ep : local_events )
    {
      int ci = ep->idx;
      int lo = std::max( 0, ci - win_smp );
      int hi = std::min( N - 1, ci + win_smp );

      // Peak entropy and min confidence in window
      double peak_H = 0.0, min_C = 1.0, tv_area = 0.0;
      for ( int i = lo; i <= hi; i++ )
	{
	  if ( H_ref[i] > peak_H ) peak_H = H_ref[i];
	  if ( C_ref[i] < min_C  ) min_C  = C_ref[i];
	  // TV area by trapezoid sum (step size = 1 sample)
	  if ( i > lo ) tv_area += 0.5 * ( TV_ref[i-1] + TV_ref[i] );
	}

      // Transition width: duration where H > H_width_th
      int width_cnt = 0;
      for ( int i = lo; i <= hi; i++ )
	if ( H_ref[i] > H_width_th ) ++width_cnt;
      double width_sec = width_cnt / Fs;

      sum_peak_H  += peak_H;
      sum_min_C   += min_C;
      sum_TV_area += tv_area / Fs;  // convert to seconds
      sum_width   += width_sec;

      // Accumulate profile only for events with a full window
      if ( lo == ci - win_smp && hi == ci + win_smp )
	{
	  ++n_events_profile;
	  for ( int j = 0; j < profile_len; j++ )
	    {
	      int si = lo + j;
	      sum_H[j]  += H_ref[si];
	      sum_C[j]  += C_ref[si];
	      sum_Mg[j] += Mg_ref[si];
	      sum_TV[j] += centered_tv( si );
	      for ( int k = 0 ; k < K ; k++ )
		sum_P[k][j] += P_ref[k][si];
	      cnt[j]++;
	    }
	}
    }

  int ne = tshape.n_events;
  tshape.mean_peak_H    = sum_peak_H  / ne;
  tshape.mean_min_C     = sum_min_C   / ne;
  tshape.mean_TV_area   = sum_TV_area / ne;
  tshape.mean_width_sec = sum_width   / ne;

  // Profile
  profile.n_events = n_events_profile;
  profile.valid    = n_events_profile >= par.min_events_profile;

  if ( profile.valid )
    {
      profile.offsets.resize( profile_len );
      profile.H.resize( profile_len );
      profile.C.resize( profile_len );
      profile.Mg.resize( profile_len );
      profile.TV.resize( profile_len );
      profile.P.resize( K );
      for ( int k = 0 ; k < K ; k++ ) profile.P[k].resize( profile_len );
      for ( int j = 0; j < profile_len; j++ )
	{
	  // Events are indexed at sample i, but both hard argmax changes and TV[i]
	  // motion steps represent the boundary between samples i-1 and i.  Shift
	  // the aligned profile by +0.5 sample so OFFSET=0 is the transition midpoint.
	  profile.offsets[j] = ( (double)( j - win_smp ) + 0.5 ) / Fs;
	  profile.H[j]  = cnt[j] > 0 ? sum_H[j]  / cnt[j] : std::numeric_limits<double>::quiet_NaN();
	  profile.C[j]  = cnt[j] > 0 ? sum_C[j]  / cnt[j] : std::numeric_limits<double>::quiet_NaN();
	  profile.Mg[j] = cnt[j] > 0 ? sum_Mg[j] / cnt[j] : std::numeric_limits<double>::quiet_NaN();
	  profile.TV[j] = cnt[j] > 0 ? sum_TV[j] / cnt[j] : std::numeric_limits<double>::quiet_NaN();
	  for ( int k = 0 ; k < K ; k++ )
	    profile.P[k][j] = cnt[j] > 0 ? sum_P[k][j] / cnt[j] : std::numeric_limits<double>::quiet_NaN();
	}
    }
}


// ============================================================
// Output helpers
// ============================================================

// Core mixedness/instability summary: valid at any granularity (top-level, ANNOT,
// STATE, SS, or within a STABLE/UNSTABLE region) since it's built from means and
// percentiles, which stay well-defined under range restriction.
static void write_region( const hd_region_stats_t & rs, bool threestate )
{
  if ( !rs.valid ) return;
  writer.value( "N",          rs.n );
  writer.value( "H",          rs.mean_H );
  writer.value( "H_P90",      rs.p90_H );
  writer.value( "CONF",       rs.mean_C );
  writer.value( "FRAC_LOW_CONF", rs.frac_C_below );
  writer.value( "MG",         rs.mean_Mg );
  writer.value( "PS",         rs.mean_TV );
  writer.value( "PS_P90",     rs.p90_TV );
  writer.value( "PS_LAG",     rs.mean_TV_lag );
  if ( threestate )
    {
	if ( rs.mean_p.size() >= 3 )
	{
	  writer.value( "P_W",   rs.mean_p[0] );
	  writer.value( "P_NR",  rs.mean_p[1] );
	  writer.value( "P_R",   rs.mean_p[2] );
	}
    }
  else
    {
      if ( rs.mean_p.size() >= 5 )
	{
	  writer.value( "P_W",   rs.mean_p[0] );
	  writer.value( "P_N1",  rs.mean_p[1] );
	  writer.value( "P_N2",  rs.mean_p[2] );
	  writer.value( "P_N3",  rs.mean_p[3] );
	  writer.value( "P_R",   rs.mean_p[4] );
	}
    }
}

// CORR_H_PS: entropy/PS correlation. Only valid at bare (unstratified-by-REGION)
// granularity — within a STABLE or UNSTABLE bucket, both entropy and PS are
// artificially range-restricted by the same thresholds that define the bucket,
// which destabilizes a correlation coefficient (unlike the means/percentiles
// above, which remain well-defined under range restriction).
static void write_region_corr( const hd_region_stats_t & rs )
{
  if ( !rs.valid ) return;
  writer.value( "CORR_H_PS", rs.corr_H_TV );
}

static void write_trans_stats( const hd_trans_stats_t & ts )
{
  writer.value( "TRANS_N",         ts.n_events );
  writer.value( "TRANS_DENS",      ts.density );
  if ( ts.valid )
    {
      writer.value( "TRANS_DUR",      ts.mean_width_sec );
      writer.value( "TRANS_PEAK_H",   ts.mean_peak_H );
      writer.value( "TRANS_MIN_CONF", ts.mean_min_C );
      writer.value( "TRANS_TOT_PS",   ts.mean_TV_area );
    }
}

static std::string hd_state_label( const bool threestate , const int st )
{
  if ( threestate )
    {
      if ( st == 0 ) return "W";
      if ( st == 1 ) return "NR";
      if ( st == 2 ) return "R";
    }
  else
    {
      if ( st == 0 ) return "W";
      if ( st == 1 ) return "N1";
      if ( st == 2 ) return "N2";
      if ( st == 3 ) return "N3";
      if ( st == 4 ) return "R";
    }
  return "?";
}

static std::string hd_run_state_label( const std::vector<int> & argmax ,
				       const int lo ,
				       const int hi ,
				       const bool threestate )
{
  if ( lo > hi || lo < 0 || hi >= (int)argmax.size() ) return ".";

  std::map<int,int> counts;
  for ( int i = lo ; i <= hi ; i++ ) ++counts[ argmax[i] ];
  if ( counts.empty() ) return ".";

  int best_st = -1;
  int best_n = -1;
  for ( std::map<int,int>::const_iterator cc = counts.begin() ; cc != counts.end() ; ++cc )
    if ( cc->second > best_n )
      {
	best_st = cc->first;
	best_n = cc->second;
      }

  if ( best_st < 0 ) return ".";
  return hd_state_label( threestate , best_st );
}

static bool hd_make_pp_nrem( edf_t & edf , param_t & param )
{
  const std::string out_label = param.has( "NR" ) ? param.value( "NR" ) : "PP_NR";
  const std::vector<std::string> in_labels = {
    param.has( "N1" ) ? param.value( "N1" ) : "PP_N1",
    param.has( "N2" ) ? param.value( "N2" ) : "PP_N2",
    param.has( "N3" ) ? param.value( "N3" ) : "PP_N3"
  };

  if ( edf.header.has_signal( out_label ) )
    {
      logger << "  HDSTATS emit-nr-pp: " << out_label << " already exists, skipping\n";
      return false;
    }

  interval_t whole = edf.timeline.wholetrace();
  std::vector<double> nr;
  std::vector<double> invalid;
  double ref_sr = 0.0;
  std::string ref_label = ".";
  int n_found = 0;
  bool have_invalid = false;
  const std::string out_invalid_label =
    param.has( "VALID" ) ? param.value( "VALID" ) : out_label + "_INVALID";

  for ( int i = 0 ; i < (int)in_labels.size() ; i++ )
    {
      signal_list_t sl = edf.header.signal_list( in_labels[i] );
      if ( sl.size() == 0 ) continue;
      if ( edf.header.is_annotation_channel( sl(0) ) )
	Helper::halt( "HDSTATS emit-nr-pp: " + in_labels[i] + " is an annotation channel, not a signal" );

      const int slot = sl(0);
      const double sr = edf.header.sampling_freq( slot );
      slice_t s( edf , slot , whole );
      const std::vector<double> * d = s.pdata();

      if ( n_found == 0 )
	{
	  ref_sr = sr;
	  ref_label = in_labels[i];
	  nr = *d;
          invalid.assign( nr.size() , 0.0 );
	}
      else
	{
	  if ( std::fabs( sr - ref_sr ) > 1e-6 )
	    Helper::halt( "HDSTATS emit-nr-pp: sample-rate mismatch for "
			  + in_labels[i] + " (" + Helper::dbl2str(sr) + " Hz) vs "
			  + ref_label + " (" + Helper::dbl2str(ref_sr) + " Hz)" );
	  if ( d->size() != nr.size() )
	    Helper::halt( "HDSTATS emit-nr-pp: length mismatch for " + in_labels[i] );
	  for ( int j = 0 ; j < (int)nr.size() ; j++ ) nr[j] += (*d)[j];
	}

      const std::string invalid_label = infer_invalid_channel_label( in_labels[i] );
      if ( invalid_label != "" )
        {
          signal_list_t sl_invalid = edf.header.signal_list( invalid_label );
          if ( sl_invalid.size() > 0 && ! edf.header.is_annotation_channel( sl_invalid(0) ) )
            {
              slice_t s_invalid( edf , sl_invalid(0) , whole );
              const std::vector<double> * d_invalid = s_invalid.pdata();
              if ( d_invalid->size() == invalid.size() )
                {
                  have_invalid = true;
                  for ( int j = 0 ; j < (int)invalid.size() ; j++ )
                    if ( (*d_invalid)[j] >= 0.5 ) invalid[j] = 1.0;
                }
            }
        }

      ++n_found;
    }

  if ( n_found == 0 )
    Helper::halt( "HDSTATS emit-nr-pp: could not find any source NREM posterior channels; "
		  "require at least one of "
		  + in_labels[0] + ", " + in_labels[1] + ", " + in_labels[2] );

  for ( int j = 0 ; j < (int)nr.size() ; j++ )
    {
      bool bad = have_invalid && invalid[j] >= 0.5;
      if ( ! bad && ( ! std::isfinite( nr[j] ) || nr[j] < 0.0 ) )
        {
          bad = true;
          invalid[j] = 1.0;
          have_invalid = true;
        }

      if ( bad ) nr[j] = 0.0;
    }

  edf.add_signal( out_label , ref_sr , nr );
  if ( have_invalid && ! edf.header.has_signal( out_invalid_label ) )
    edf.add_signal( out_invalid_label , ref_sr , invalid );
  logger << "  HDSTATS emit-nr-pp: added " << out_label
	 << " from " << n_found << " available NREM posterior channel(s)\n";
  return true;
}

static void hd_emit_signals( edf_t & edf,
			      const hd_data_t & dat,
			      const hd_derived_t & der,
			      const hd_params_t & par )
{
  auto add_sig = [&]( const std::string & label, const std::vector<double> & x )
    {
      if ( edf.header.has_signal( label ) )
	{
	  logger << "  HDSTATS emit-sigs: " << label << " already exists, skipping\n";
	  return;
	}
      edf.add_signal( label, dat.Fs, x );
    };

  add_sig( par.emit_sig_label + "_PS",   der.TV );
  add_sig( par.emit_sig_label + "_H",    der.H );
  add_sig( par.emit_sig_label + "_CONF", der.C );

  if ( par.do_3state )
    {
      add_sig( par.emit_sig_label + "_PS3",   der.TV3 );
      add_sig( par.emit_sig_label + "_H3",    der.H3 );
      add_sig( par.emit_sig_label + "_CONF3", der.C3 );
    }

  logger << "  HDSTATS emit-sigs: added PS/H/CONF signal(s) under prefix '" << par.emit_sig_label << "'\n";
}

static void write_hd_hypno_metrics( const hd_hypno_metrics_t & hm )
{
  if ( hm.valid_onset ) writer.value( "HD_SOL" , hm.hd_sol );
  if ( hm.valid_window )
    {
      writer.value( "HD_SPT" , hm.hd_spt );
      writer.value( "HD_WASO" , hm.hd_waso );
      writer.value( "HD_SME" , hm.hd_sme );
    }
  if ( hm.valid_rem )
    {
      writer.value( "HD_REM_LAT" , hm.hd_rem_lat );
      writer.value( "HD_REM_LAT2" , hm.hd_rem_lat2 );
    }

  static const std::vector<std::string> stages = { "W", "N1", "N2", "N3", "R", "NR", "S" };
  for ( const auto & st : stages )
    {
      std::map<std::string,double>::const_iterator ii = hm.mins.find( st );
      if ( ii == hm.mins.end() ) continue;
      writer.level( st , "SS" );
      writer.value( "HD_MINS" , ii->second );
      // W and S: % of total recording time; sleep sub-stages: % of TST
      const bool use_total = ( st == "W" || st == "S" );
      const double denom = use_total ? hm.total_min : hm.hd_tst;
      writer.value( "HD_PCT" , denom > 0 ? 100.0 * ii->second / denom : 0.0 );
      writer.value( "HD_DENS" , hm.hd_spt > 0 ? ii->second / hm.hd_spt : 0.0 );
      writer.unlevel( "SS" );
    }
}

static std::vector<bool> make_stage_mask(
    const std::vector<bool> & region_mask,
    const std::vector<int> & argmax,
    const int stage )
{
  std::vector<bool> mask( region_mask.size(), false );
  for ( int i = 0; i < (int)region_mask.size(); i++ )
    mask[i] = region_mask[i] && argmax[i] == stage;
  return mask;
}

// Write HDSTATS rows for one context (global or per-stratum).
// Emits a bare (unstratified) summary plus REGION strata: STABLE, UNSTABLE.
// Emits STATE stratum only when do_3state is true.
// Emits SS stratum for primary summary metrics only, based on argmax stage.
static void write_hdstats(
    const hd_derived_t & der,
    const hd_trans_t & tr,
    const std::vector<bool> & region_mask5,   // which samples are in this context (5-state)
    const std::vector<bool> & region_mask3,   // same for 3-state (may be same as region_mask5)
    const hd_params_t & par,
    const hd_data_t & dat,
    const hd_hypno_metrics_t * hm )
{
  if ( hm != NULL && par.do_hd_metrics )
    write_hd_hypno_metrics( *hm );

  auto make_all   = [&]( const std::vector<bool> & rm ) { return rm; };
  auto make_stable = [&]( const std::vector<bool> & rm,
			   const std::vector<bool> & is_stable ) {
    std::vector<bool> m( rm.size() );
    for ( int i = 0; i < (int)rm.size(); i++ ) m[i] = rm[i] && is_stable[i];
    return m;
  };
  // Helper: emit one state (5 or 3)
  auto emit_state = [&]( bool threestate ) {

    const std::vector<bool> & rmask       = threestate ? region_mask3 : region_mask5;
    const std::vector<bool> & is_stable   = threestate ? tr.is_stable3 : tr.is_stable5;
    const std::vector<bool> & is_bout_stable = threestate ? tr.is_bout_stable3 : tr.is_bout_stable5;
    const std::vector<int>  & argmax      = threestate ? der.argmax3    : der.argmax5;
    const int K = threestate ? 3 : 5;

    auto emit_context = [&]( const std::vector<bool> & context_mask,
			     const bool emit_top_level_transition_stats,
			     const bool emit_profiles_and_pairs )
      {
	// Helper: emit a transition profile vector under the current writer context.
	// Per-state posteriors are only meaningful once events are restricted to a
	// single pair (BOUT_TRANS,SEC) — pooled across mixed pair types they average
	// together unrelated state identities, so the generic SEC table omits them.
	auto emit_profile_rows = [&]( const hd_profile_t & pf, const bool emit_posteriors ) {
	  for ( int j = 0; j < (int)pf.offsets.size(); j++ )
	    {
	      writer.level( pf.offsets[j], globals::sec_strat );
	      writer.value( "H",  pf.H[j] );
	      writer.value( "CONF", pf.C[j] );
	      writer.value( "MG", pf.Mg[j] );
	      writer.value( "PS", pf.TV[j] );
	      if ( emit_posteriors )
		{
		  if ( threestate )
		    {
		      writer.value( "P_W",  pf.P[0][j] );
		      writer.value( "P_NR", pf.P[1][j] );
		      writer.value( "P_R",  pf.P[2][j] );
		    }
		  else
		    {
		      writer.value( "P_W",  pf.P[0][j] );
		      writer.value( "P_N1", pf.P[1][j] );
		      writer.value( "P_N2", pf.P[2][j] );
		      writer.value( "P_N3", pf.P[3][j] );
		      writer.value( "P_R",  pf.P[4][j] );
		    }
		}
	      writer.unlevel( globals::sec_strat );
	    }
	};

	// HDSTATS table: bare (unstratified-by-STABLE) summary, plus a binary
	// STABLE/UNSTABLE breakdown. UNSTABLE merges what were previously the
	// separate TRANS (near a crossing) and unlabeled remainder buckets —
	// the finer distinction is still available via the event-based
	// TRANS/SEC tables below, which are unaffected by this.
	{
	  auto make_unstable = [&]( const std::vector<bool> & rm ) {
	    std::vector<bool> m( rm.size() );
	    for ( int i = 0; i < (int)rm.size(); i++ )
	      m[i] = rm[i] && !is_stable[i];
	    return m;
	  };

	  hd_region_stats_t rs_all      = compute_region( dat, der, make_all(context_mask),               par, threestate );
	  hd_region_stats_t rs_stable   = compute_region( dat, der, make_stable(context_mask, is_stable),  par, threestate );
	  hd_region_stats_t rs_unstable = compute_region( dat, der, make_unstable(context_mask),           par, threestate );

	  const bool ok_all      = rs_all.valid      && rs_all.n      >= par.min_samples_region;
	  const bool ok_stable   = rs_stable.valid   && rs_stable.n   >= par.min_samples_region;
	  const bool ok_unstable = rs_unstable.valid && rs_unstable.n >= par.min_samples_region;

	  // Bare summary: no REGION token, so this lands wherever we're already
	  // nested (top-level/ANNOT/STATE for the main pass, bare SS=<stage> for
	  // the stage-conditioned pass). CORR_H_PS and PCT_STABLE only make sense
	  // here, not within a STABLE/UNSTABLE bucket (see write_region_corr).
	  if ( ok_all )
	    {
	      write_region( rs_all, threestate );
	      write_region_corr( rs_all );
	      writer.value( "PCT_STABLE", 100.0 * rs_stable.n / rs_all.n );
	    }

	  if ( ok_stable )
	    {
	      writer.level( 1, "STABLE" );
	      write_region( rs_stable, threestate );
	      writer.unlevel( "STABLE" );
	    }

	  if ( ok_unstable )
	    {
	      writer.level( 0, "STABLE" );
	      write_region( rs_unstable, threestate );
	      writer.unlevel( "STABLE" );
	    }

	  // Stable-vs-unstable contrasts (no REGION stratum; top-level context only)
	  if ( emit_top_level_transition_stats && ok_stable && ok_unstable )
	    {
	      writer.value( "H_RATIO_UN_ST",   rs_stable.mean_H  > 0
			    ? rs_unstable.mean_H  / rs_stable.mean_H  : std::numeric_limits<double>::quiet_NaN() );
	      writer.value( "CONF_DIFF_UN_ST", rs_unstable.mean_C - rs_stable.mean_C );
	      writer.value( "PS_RATIO_UN_ST",  rs_stable.mean_TV > 0
			    ? rs_unstable.mean_TV / rs_stable.mean_TV : std::numeric_limits<double>::quiet_NaN() );
	    }
	}

	// HDTRANS table (transition shape)
	{
	  const std::vector<hd_event_t> & evts = threestate ? tr.events3 : tr.events5;
	  hd_trans_stats_t tshape;
	  hd_profile_t     profile;
	  compute_trans_shape( der, evts, context_mask, is_stable, dat, par, dat.Fs, threestate,
                               false, false, tshape, profile );

	  // Bouts: long, settled, single-dominant-stage runs (see bout-min/bout-gap).
	  // Reused below both to count bout-to-bout transitions (N_BOUT_TRANS) and,
	  // via compute_bout_transitions(), to anchor every OFFSET profile — OFFSET
	  // no longer has its own independent flank-qualification notion.
	  const std::vector<hd_bout_t> bouts = compute_bouts( context_mask, is_bout_stable, argmax,
	                                                       par.bout_min_sec, par.bout_gap_sec, dat.Fs );
	  const std::vector<hd_event_t> bout_events = compute_bout_transitions( bouts, evts, threestate ? der.TV3 : der.TV );

	  // TRANS events restricted to this context (matches how TRANS_N itself, via
	  // compute_trans_shape's region_mask filtering, is scoped).
	  std::vector<hd_event_t> evts_in_context;
	  for ( const auto & e : evts )
	    if ( e.idx >= 0 && e.idx < (int)context_mask.size() && context_mask[e.idx] )
	      evts_in_context.push_back( e );

	  if ( emit_top_level_transition_stats )
	    {
	      write_trans_stats( tshape );

	      const int n_matched = count_trans_bout_matches( evts_in_context, bouts );
	      writer.value( "TRANS_N_MATCHED",   n_matched );
	      writer.value( "TRANS_PCT_MATCHED", tshape.n_events > 0 ? 100.0 * n_matched / tshape.n_events : 0.0 );

	      int n_in_region = 0;
	      for ( bool b : context_mask ) if (b) ++n_in_region;
	      double dur_hr = (double)n_in_region / dat.Fs / 3600.0;
	      writer.value( "BOUT_TRANS_N",    (int)bout_events.size() );
	      writer.value( "BOUT_TRANS_DENS", dur_hr > 0 ? (int)bout_events.size() / dur_hr : 0.0 );

	      // MIX table: state-set mixing, in the same MINS/DENS units as the HD_*
	      // hypnodensity metrics (top-level context only; not under SS, matching
	      // how HD_MINS/HD_DENS themselves are only available where hm != NULL).
	      {
		const std::vector<std::vector<double>> & P_ref = threestate ? dat.p3 : dat.p;
		for ( const auto & combo : hd_mix_combinations( threestate ) )
		  {
		    std::string label = hd_state_label( threestate, combo[0] );
		    for ( int ci = 1; ci < (int)combo.size(); ci++ )
		      label += "_" + hd_state_label( threestate, combo[ci] );
		    const double mix_mins = compute_mix_mins( P_ref, context_mask, combo, dat.Fs );
		    writer.level( label, "MIX" );
		    writer.value( "HD_MINS", mix_mins );
		    if ( hm != NULL )
		      writer.value( "HD_DENS", hm->hd_spt > 0 ? mix_mins / hm->hd_spt : 0.0 );
		    writer.unlevel( "MIX" );
		  }
	      }
	    }

	  if ( ! emit_profiles_and_pairs ) return;

          // Generic OFFSET: aligned profile pooled across all bout-to-bout transitions.
          hd_trans_stats_t tshape_bout;
          hd_profile_t     profile_bout;
          compute_trans_shape( der, bout_events, context_mask, is_stable, dat, par, dat.Fs, threestate,
                               false, false, tshape_bout, profile_bout );
	  if ( profile_bout.valid )
	    emit_profile_rows( profile_bout, false );

	  // TRANS: pair-specific summary statistics still use directional events (unchanged).
	  std::set<std::pair<int,int>> trans_pairs;
	  for ( const auto & e : evts )
	    if ( e.is_dir && e.from_st >= 0 && e.to_st >= 0 && e.from_st != e.to_st )
	      trans_pairs.insert( std::make_pair( e.from_st , e.to_st ) );

	  for ( const auto & tp : trans_pairs )
	    {
	      std::vector<hd_event_t> evts_pair;
	      for ( const auto & e : evts )
		if ( e.from_st == tp.first && e.to_st == tp.second )
		  evts_pair.push_back( e );

	      if ( evts_pair.empty() ) continue;

	      hd_trans_stats_t tshape_pair;
              hd_profile_t     profile_pair_dir;
	      compute_trans_shape( der, evts_pair, context_mask, is_stable, dat, par, dat.Fs, threestate,
                                   true, false, tshape_pair, profile_pair_dir );

	      std::vector<hd_event_t> evts_pair_in_context;
	      for ( const auto & e : evts_pair )
		if ( e.idx >= 0 && e.idx < (int)context_mask.size() && context_mask[e.idx] )
		  evts_pair_in_context.push_back( e );
	      const int n_matched_pair = count_trans_bout_matches( evts_pair_in_context, bouts );

	      const std::string trans_label = hd_state_label( threestate , tp.first ) + "to"
		+ hd_state_label( threestate , tp.second );

	      writer.level( trans_label , "TRANS" );
	      write_trans_stats( tshape_pair );
	      writer.value( "TRANS_N_MATCHED",   n_matched_pair );
	      writer.value( "TRANS_PCT_MATCHED", tshape_pair.n_events > 0 ? 100.0 * n_matched_pair / tshape_pair.n_events : 0.0 );
	      writer.unlevel( "TRANS" );
	    }

	  // BOUT: pair-specific aligned profiles, one per observed bout-to-bout stage pair.
	  std::set<std::pair<int,int>> bout_pairs;
	  for ( const auto & e : bout_events )
	    bout_pairs.insert( std::make_pair( e.from_st , e.to_st ) );

	  for ( const auto & tp : bout_pairs )
	    {
	      std::vector<hd_event_t> bout_events_pair;
	      for ( const auto & e : bout_events )
		if ( e.from_st == tp.first && e.to_st == tp.second )
		  bout_events_pair.push_back( e );

	      if ( bout_events_pair.empty() ) continue;

	      hd_trans_stats_t tshape_bout_pair;
              hd_profile_t     profile_bout_pair;
	      compute_trans_shape( der, bout_events_pair, context_mask, is_stable, dat, par, dat.Fs, threestate,
                                   false, false, tshape_bout_pair, profile_bout_pair );

	      const std::string bout_label = hd_state_label( threestate , tp.first ) + "to"
		+ hd_state_label( threestate , tp.second );

	      writer.level( bout_label , "BOUT_TRANS" );
	      if ( profile_bout_pair.valid )
		emit_profile_rows( profile_bout_pair, true );
	      writer.unlevel( "BOUT_TRANS" );
	    }
	}
      };

    // Overall context output remains unchanged and includes detailed transition profiles.
    emit_context( rmask, true, true );

    // Stage-conditioned output is REGION-only; do not emit separate top-level
    // transition or ratio summaries under SS.
    for ( int st = 0; st < K; st++ )
      {
	const std::vector<bool> stage_mask = make_stage_mask( rmask, argmax, st );
	bool any = false;
	for ( int i = 0; i < (int)stage_mask.size(); i++ )
	  if ( stage_mask[i] ) { any = true; break; }
	if ( ! any ) continue;
	writer.level( hd_state_label( threestate, st ), "SS" );
	emit_context( stage_mask, false, false );
	writer.unlevel( "SS" );
      }
  };

  // NSS is always emitted, even when only 5-state is computed, so table shape
  // never depends on the 3state parameter.
  writer.level( "5", "NSS" );
  emit_state( false );
  writer.unlevel( "NSS" );

  if ( par.do_3state )
    {
      writer.level( "3", "NSS" );
      emit_state( true );
      writer.unlevel( "NSS" );
    }
}

static void set_hdstats_offset_outputs( bool enabled )
{
  // Generic SEC: no per-state posteriors (pooled across mixed pair types).
  const std::vector<std::string> sec_tables = {
    "NSS,SEC",
    "ANNOT,NSS,SEC"
  };
  const std::vector<std::string> sec_vars = {
    "H", "CONF", "MG", "PS"
  };

  // Pair-specific BOUT_TRANS,SEC: posteriors are meaningful once restricted to one pair.
  const std::vector<std::string> bout_sec_tables = {
    "NSS,BOUT_TRANS,SEC",
    "ANNOT,NSS,BOUT_TRANS,SEC"
  };
  const std::vector<std::string> bout_sec_vars = {
    "H", "CONF", "MG", "PS", "P_W", "P_NR", "P_N1", "P_N2", "P_N3", "P_R"
  };

  for ( int i = 0; i < (int)sec_tables.size(); i++ )
    for ( int j = 0; j < (int)sec_vars.size(); j++ )
      globals::cmddefs().register_var( "HDSTATS", sec_tables[i], sec_vars[j], enabled );

  for ( int i = 0; i < (int)bout_sec_tables.size(); i++ )
    for ( int j = 0; j < (int)bout_sec_vars.size(); j++ )
      globals::cmddefs().register_var( "HDSTATS", bout_sec_tables[i], bout_sec_vars[j], enabled );
}

struct hdstats_offset_outputs_scope_t
{
  hdstats_offset_outputs_scope_t()  { set_hdstats_offset_outputs( true ); }
  ~hdstats_offset_outputs_scope_t() { set_hdstats_offset_outputs( false ); }
};


// ============================================================
// Layer 5: top-level command
// ============================================================

void proc_hdstats( edf_t & edf, param_t & param )
{
  // Detection/bout logic throughout this file assumes samples 0..N-1 are
  // contiguous and evenly spaced in time; a genuine gap would silently treat
  // non-adjacent samples as neighbors (see hd_data_t::load()).
  if ( edf.is_actually_discontinuous() )
    Helper::halt( "HDSTATS requires a single-segment EDF (no within-recording gaps)" );

  hdstats_offset_outputs_scope_t offset_outputs_scope;

  const bool emit_nr_pp =
    ( param.has( "emit-nr-pp" ) && param.yesno( "emit-nr-pp" ) );
  const bool emit_stats = param.has( "emit-stats" ) ? param.yesno( "emit-stats" ) : true;

  if ( emit_nr_pp )
    hd_make_pp_nrem( edf , param );

  if ( ! emit_stats ) return;

  // Parse parameters
  hd_params_t par;

  // Channel names for W, N1, N2, N3, R
  par.ch.resize( 5 );
  par.ch[0] = param.has( "W" )  ? param.value( "W" )  : "PP_W";
  par.ch[1] = param.has( "N1" ) ? param.value( "N1" ) : "PP_N1";
  par.ch[2] = param.has( "N2" ) ? param.value( "N2" ) : "PP_N2";
  par.ch[3] = param.has( "N3" ) ? param.value( "N3" ) : "PP_N3";
  par.ch[4] = param.has( "R" )  ? param.value( "R" )  : "PP_R";

  par.do_hd_metrics    = param.has( "hd-metrics" ) ? param.yesno( "hd-metrics" ) : true;
  par.do_3state       = param.yesno( "3state", true );
  par.window_sec      = param.has( "window"    ) ? param.requires_dbl( "window"     ) : 60.0;
  par.lag_sec         = param.has( "lag"       ) ? param.requires_dbl( "lag"        ) : 30.0;
  par.stable_min_sec  = param.has( "stable-min") ? param.requires_dbl( "stable-min" ) : 30.0;
  par.conf_th         = param.has( "conf-th"   ) ? param.requires_dbl( "conf-th"    ) : 0.8;
  par.hd_smooth_sec   = param.has( "hd-smooth" ) ? param.requires_dbl( "hd-smooth" ) : 30.0;
  par.hd_onset_win_min = param.has( "hd-onset-win" ) ? param.requires_dbl( "hd-onset-win" ) : 10.0;
  par.hd_sleep_mass_th = param.has( "hd-sleep-mass-th" ) ? param.requires_dbl( "hd-sleep-mass-th" ) : 0.60;
  par.hd_rem_mass_th   = param.has( "hd-rem-mass-th" ) ? param.requires_dbl( "hd-rem-mass-th" ) : 0.30;
  par.hd_offset_win_min = param.has( "hd-offset-win" ) ? param.requires_dbl( "hd-offset-win" ) : 10.0;
  par.hd_offset_mass_th = param.has( "hd-offset-mass-th" ) ? param.requires_dbl( "hd-offset-mass-th" ) : 0.60;
  par.annot_name      = param.has( "annot"     ) ? param.value( "annot" )              : "";
  par.verbose         = param.yesno( "verbose" );
  par.min_events_profile = param.has( "min-events" ) ? param.requires_int( "min-events" ) : 3;
  par.min_samples_region = param.has( "min-samples" ) ? param.requires_int( "min-samples" ) : 20;
  par.emit_annot_typed_classes = param.has( "emit-annots-by-type" ) ? param.yesno( "emit-annots-by-type" ) : false;
  par.emit_annot_label   =
    param.has( "emit-annots" ) ? ( ( ! param.empty("emit-annots" ) ) ? param.value( "emit-annots" ) : "hd" )
    : ( par.emit_annot_typed_classes ? "hd" : "" );
  par.stable_ps_th    = param.has( "stable-ps"   ) ? param.requires_dbl( "stable-ps"   ) : 0.05;
  par.stable_conf_th  = param.has( "stable-conf" ) ? param.requires_dbl( "stable-conf" ) : 0.70;
  par.min_shift       = param.has( "min-shift"   ) ? param.requires_dbl( "min-shift"   ) : 0.10;
  par.shift_win_sec   = param.has( "shift-win"   ) ? param.requires_dbl( "shift-win"   ) : 30.0;
  par.bout_min_sec    = param.has( "bout-min"    ) ? param.requires_dbl( "bout-min"    ) : 120.0;
  par.bout_gap_sec    = param.has( "bout-gap"    ) ? param.requires_dbl( "bout-gap"    ) : 30.0;
  par.bout_smooth_sec = param.has( "bout-smooth" ) ? param.requires_dbl( "bout-smooth" ) : 15.0;
  par.ch_valid        = param.has( "VALID" ) ? param.value( "VALID" ) : infer_invalid_channel_label( par.ch );
  par.emit_sig_label  =
    param.has( "emit-sigs" ) ? ( ( ! param.empty( "emit-sigs" ) ) ? param.value( "emit-sigs" ) : "hd" ) : "";
  par.do_emit_hypnogram = param.has( "emit-hypnogram" );
  par.hypno_prefix    =
    ( par.do_emit_hypnogram && ! param.empty( "emit-hypnogram" ) )
    ? param.value( "emit-hypnogram" ) + "_" : "";

  logger << "  loading hypnodensity channels ["
	 << par.ch[0] << ", " << par.ch[1] << ", " << par.ch[2] << ", "
	 << par.ch[3] << ", " << par.ch[4] << "]\n";

  // Layer 1: load data
  hd_data_t dat;
  if ( ! dat.load( edf, par ) ) return;

  logger << "  N=" << dat.N << " samples at Fs=" << dat.Fs << " Hz ("
	 << dat.N / dat.Fs / 3600.0 << " hours)\n";

  // Layer 2: derive per-sample signals
  hd_derived_t der;
  der.compute( dat, par );

  if ( ! par.emit_sig_label.empty() )
    hd_emit_signals( edf, dat, der, par );

  // Layer 3: detect transitions
  hd_trans_t tr;
  tr.detect( dat, der, par );

  logger << "  detected " << tr.events5.size() << " transition events (5-state)\n";
  if ( par.do_3state )
    logger << "  detected " << tr.events3.size() << " transition events (3-state)\n";

  // Optionally emit transition events as annotations
  if ( ! par.emit_annot_label.empty() && ! dat.tp.empty() )
    {
      const uint64_t sample_tp = (uint64_t)std::round( (double)globals::tp_1sec / dat.Fs );

      // Helper: emit contiguous runs of a boolean mask as whole-interval annotations
      auto emit_runs = [&]( annot_t * a ,
			    const std::vector<bool> & mask ,
			    const std::string & inst ,
			    const std::vector<int> * argmax = NULL ,
			    const bool threestate = false )
        {
          int i = 0;
          while ( i < dat.N )
            {
              if ( !mask[i] ) { ++i; continue; }
              int lo = i;
              while ( i < dat.N && mask[i] ) ++i;
              int hi = i - 1;
              const std::string run_inst =
		argmax == NULL ? inst : hd_run_state_label( *argmax , lo , hi , threestate );
              a->add( run_inst , interval_t( dat.tp[lo] , dat.tp[hi] + sample_tp ) , "." );
            }
        };

      // Emit TRANS events as point events; use instance label "from->to" or "." if states unknown
      const uint64_t half_sample_tp = sample_tp / 2;
      for ( const auto & e : tr.events5 )
	{
	  const std::string inst =
	    ( e.is_dir && e.from_st >= 0 && e.to_st >= 0 )
	    ? hd_state_label( false, e.from_st ) + "to" + hd_state_label( false, e.to_st )
	    : ".";
	  const std::string cls =
	    ( par.emit_annot_typed_classes && inst != "." )
	    ? par.emit_annot_label + "_trans_event_" + inst
	    : par.emit_annot_label + "_trans_event";
	  uint64_t tp0 = dat.tp[ e.idx ] >= half_sample_tp ? dat.tp[ e.idx ] - half_sample_tp : dat.tp[ e.idx ];
	  edf.annotations->add( cls )->add( inst , interval_t( tp0 , tp0 ) , "." );
	}

      // Emit BOUT_TRANS transitions as point events, anchored at the bout-gap
      // midpoint; from/to states are always known (bouts are always dominant-state runs)
      for ( const auto & e : tr.bout_events5 )
	{
	  const std::string inst = hd_state_label( false, e.from_st ) + "to" + hd_state_label( false, e.to_st );
	  const std::string cls =
	    par.emit_annot_typed_classes
	    ? par.emit_annot_label + "_bout_trans_event_" + inst
	    : par.emit_annot_label + "_bout_trans_event";
	  uint64_t tp0 = dat.tp[ e.idx ] >= half_sample_tp ? dat.tp[ e.idx ] - half_sample_tp : dat.tp[ e.idx ];
	  edf.annotations->add( cls )->add( inst , interval_t( tp0 , tp0 ) , "." );
	}

      // Emit stable, trans, and bout_trans whole intervals (5-state) -- the same
      // three independent flags as the TIME table's STABLE/TRANS/BOUT_TRANS
      annot_t * atr_stable = edf.annotations->add( par.emit_annot_label + "_stable" );
      emit_runs( atr_stable , tr.is_stable5 , "stable" , &der.argmax5 , false );

      annot_t * atr_trans = edf.annotations->add( par.emit_annot_label + "_trans" );
      emit_runs( atr_trans , tr.is_trans5 , "trans" );

      annot_t * atr_bout_trans = edf.annotations->add( par.emit_annot_label + "_bout_trans" );
      emit_runs( atr_bout_trans , tr.is_bout_trans5 , "bout_trans" );

      // If 3-state mode, also emit 3-state transitions, stable, trans, and bout_trans intervals
      if ( par.do_3state )
	{
	  for ( const auto & e : tr.events3 )
	    {
	      const std::string inst =
		( e.is_dir && e.from_st >= 0 && e.to_st >= 0 )
		? hd_state_label( true, e.from_st ) + "to" + hd_state_label( true, e.to_st )
		: ".";
	      const std::string cls =
		( par.emit_annot_typed_classes && inst != "." )
		? par.emit_annot_label + "_3_trans_event_" + inst
		: par.emit_annot_label + "_3_trans_event";
	      uint64_t tp0 = dat.tp[ e.idx ] >= half_sample_tp ? dat.tp[ e.idx ] - half_sample_tp : dat.tp[ e.idx ];
	      edf.annotations->add( cls )->add( inst , interval_t( tp0 , tp0 ) , "." );
	    }

	  for ( const auto & e : tr.bout_events3 )
	    {
	      const std::string inst = hd_state_label( true, e.from_st ) + "to" + hd_state_label( true, e.to_st );
	      const std::string cls =
		par.emit_annot_typed_classes
		? par.emit_annot_label + "_3_bout_trans_event_" + inst
		: par.emit_annot_label + "_3_bout_trans_event";
	      uint64_t tp0 = dat.tp[ e.idx ] >= half_sample_tp ? dat.tp[ e.idx ] - half_sample_tp : dat.tp[ e.idx ];
	      edf.annotations->add( cls )->add( inst , interval_t( tp0 , tp0 ) , "." );
	    }

	  annot_t * atr3_stable = edf.annotations->add( par.emit_annot_label + "_3_stable" );
	  emit_runs( atr3_stable , tr.is_stable3 , "stable" , &der.argmax3 , true );

	  annot_t * atr3_trans = edf.annotations->add( par.emit_annot_label + "_3_trans" );
	  emit_runs( atr3_trans , tr.is_trans3 , "trans" );

	  annot_t * atr3_bout_trans = edf.annotations->add( par.emit_annot_label + "_3_bout_trans" );
	  emit_runs( atr3_bout_trans , tr.is_bout_trans3 , "bout_trans" );
	}

      logger << "  emitting transition annotations under '" << par.emit_annot_label << "'\n";
    }

  // Optionally emit hard-called (argmax) stage runs as annotations, one class
  // per stage (e.g. N1, N2, N3, R, W, or prefix_N1 etc if a prefix was given).
  // Pooled to conventional 30s epochs (mean posterior over each epoch's
  // sub-samples, then argmax) rather than called at the native sample
  // resolution -- otherwise, whenever the source posteriors are sampled
  // faster than 30s (e.g. 5s blocks), sub-epoch flicker produces a mass of
  // very short runs that don't line up with any single 30s epoch. This
  // assumes samples 0..N-1 are contiguous and start at t=0 (true for a
  // continuous EDF); it does not account for timeline gaps -- see the
  // detect()/compute_bouts() machinery elsewhere in this file, which shares
  // the same sample-index-contiguity assumption throughout.
  if ( par.do_emit_hypnogram && ! dat.tp.empty() )
    {
      const uint64_t sample_tp = (uint64_t)std::round( (double)globals::tp_1sec / dat.Fs );
      const int epoch_smp = std::max( 1, (int)std::round( 30.0 * dat.Fs ) );
      const int n_epochs  = ( dat.N + epoch_smp - 1 ) / epoch_smp;

      std::vector<int>  epoch_argmax( n_epochs, -1 );
      std::vector<bool> epoch_valid( n_epochs, false );
      std::vector<int>  epoch_lo( n_epochs ), epoch_hi( n_epochs );

      for ( int e = 0; e < n_epochs; e++ )
	{
	  const int lo = e * epoch_smp;
	  const int hi = std::min( dat.N, lo + epoch_smp ) - 1;
	  epoch_lo[e] = lo;
	  epoch_hi[e] = hi;

	  double sum[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
	  int n_valid = 0;
	  for ( int i = lo; i <= hi; i++ )
	    if ( dat.valid[i] )
	      {
		for ( int k = 0; k < 5; k++ ) sum[k] += dat.p[k][i];
		++n_valid;
	      }
	  if ( n_valid == 0 ) continue;

	  epoch_valid[e] = true;
	  int best = 0;
	  for ( int k = 1; k < 5; k++ )
	    if ( sum[k] > sum[best] ) best = k;
	  epoch_argmax[e] = best;
	}

      int e = 0;
      while ( e < n_epochs )
	{
	  if ( ! epoch_valid[e] ) { ++e; continue; }
	  const int lo_e = e;
	  const int st = epoch_argmax[e];
	  while ( e < n_epochs && epoch_valid[e] && epoch_argmax[e] == st ) ++e;
	  const int hi_e = e - 1;
	  const std::string cls = par.hypno_prefix + hd_state_label( false, st );
	  edf.annotations->add( cls )->add(
	    "." , interval_t( dat.tp[ epoch_lo[lo_e] ] , dat.tp[ epoch_hi[hi_e] ] + sample_tp ) , "." );
	}
      logger << "  emitting hard-called stage annotations (pooled to 30s epochs) under prefix '"
	     << par.hypno_prefix << "'\n";
    }

  hd_hypno_metrics_t hm;
  if ( par.do_hd_metrics )
    hm.compute( dat, par );

  // All-samples mask: valid samples only
  std::vector<bool> all_mask( dat.N, false );
  for ( int i = 0; i < dat.N; i++ ) all_mask[i] = dat.valid[i];

  // Emit validity summary
  {
    const int n_valid = (int)std::count( dat.valid.begin(), dat.valid.end(), true );
    writer.value( "N_VALID",   n_valid );
    writer.value( "PCT_VALID", dat.N > 0 ? 100.0 * n_valid / dat.N : 0.0 );
    if ( dat.n_invalid_lead  > 0 ) writer.value( "N_INVALID_LEAD",  dat.n_invalid_lead );
    if ( dat.n_invalid_trail > 0 ) writer.value( "N_INVALID_TRAIL", dat.n_invalid_trail );
    if ( dat.n_invalid_mid   > 0 ) writer.value( "N_INVALID_MID",   dat.n_invalid_mid );
  }

  // Global summary
  write_hdstats( der, tr, all_mask, all_mask, par, dat, par.do_hd_metrics ? &hm : NULL );

  // Optional HDSIG verbose output
  if ( par.verbose )
    {
      const double sec0 = dat.tp.empty() ? 0.0 : (double)dat.tp[0] / (double)globals::tp_1sec;
      for ( int i = 0; i < dat.N; i++ )
	{
	  double t_sec = dat.tp.empty()
	    ? i / dat.Fs
	    : (double)dat.tp[i] / (double)globals::tp_1sec - sec0;

	  writer.level( t_sec, "TIME" );
	  writer.value( "H",        der.H[i] );
	  writer.value( "CONF",     der.C[i] );
	  writer.value( "MG",       der.Mg[i] );
	  writer.value( "PS",       der.TV[i] );
	  writer.value( "PS_LAG",   der.TV_lag[i] );
	  writer.value( "ARGMAX",   der.argmax5[i] );
	  writer.value( "IS_TRANS",      (int)tr.is_trans5[i] );
	  writer.value( "IS_BOUT_TRANS", (int)tr.is_bout_trans5[i] );
	  writer.value( "IS_STABLE",     (int)tr.is_stable5[i] );
	  if ( par.do_3state )
	    {
	      writer.value( "H3",          der.H3[i] );
	      writer.value( "CONF3",       der.C3[i] );
	      writer.value( "PS3",         der.TV3[i] );
	      writer.value( "PS3_LAG",     der.TV3_lag[i] );
	      writer.value( "ARGMAX3",     der.argmax3[i] );
	      writer.value( "IS_TRANS3",      (int)tr.is_trans3[i] );
	      writer.value( "IS_BOUT_TRANS3", (int)tr.is_bout_trans3[i] );
	      writer.value( "IS_STABLE3",     (int)tr.is_stable3[i] );
	    }
	  writer.unlevel( "TIME" );
	}
    }

  // Stratification by annotation
  if ( par.annot_name.empty() ) return;

  annot_t * annot = (*edf.annotations)( par.annot_name );
  if ( annot == NULL )
    {
      logger << "  annotation '" << par.annot_name << "' not found; skipping stratification\n";
      return;
    }
  if ( annot->empty() )
    Helper::halt( " annotation '" + par.annot_name + "' is empty" );

  // Collect all unique level IDs in this annotation
  std::set<std::string> levels;
  {
    annot_map_t::const_iterator ii = annot->interval_events.begin();
    while ( ii != annot->interval_events.end() )
      {
	levels.insert( ii->first.id );
	++ii;
      }
  }

  // For each level, build a sample mask and compute summary
  for ( const std::string & lv : levels )
    {
      std::vector<bool> mask5( dat.N, false );

      annot_map_t::const_iterator ii = annot->interval_events.begin();
      while ( ii != annot->interval_events.end() )
	{
	  if ( ii->first.id == lv )
	    {
	      uint64_t a_start = ii->first.interval.start;
	      uint64_t a_stop  = ii->first.interval.stop;

	      // Map annotation interval to sample indices using binary search on tp[]
	      if ( !dat.tp.empty() )
		{
		  // First sample >= a_start
		  int lo = (int)( std::lower_bound( dat.tp.begin(), dat.tp.end(), a_start )
				  - dat.tp.begin() );
		  // Last sample < a_stop
		  int hi = (int)( std::lower_bound( dat.tp.begin(), dat.tp.end(), a_stop )
				  - dat.tp.begin() );
		  for ( int i = lo; i < hi && i < dat.N; i++ )
		    mask5[i] = true;
		}
	    }
	  ++ii;
	}

      // 3-state uses the same time mask (same sample positions)
      writer.level( lv, globals::annot_strat );
      write_hdstats( der, tr, mask5, mask5, par, dat, NULL );
      writer.unlevel( globals::annot_strat );
    }
}
