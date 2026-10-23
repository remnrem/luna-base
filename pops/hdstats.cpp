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
#include "hdstats-pc.h"
#include "hdstats-model.h"

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
#include <array>
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
  bool qc;                          // emit operational and fit diagnostics
  int min_events_profile;           // minimum events needed to emit aligned profile
  int min_samples_region;           // minimum samples needed to emit REGION summaries
  std::string emit_annot_label;     // if non-empty, emit transition annotations under this class
  double stable_ps_th;              // per-sample PS must be < this to qualify as stable
  double stable_conf_th;            // per-sample C must be > this to qualify as stable
  double min_shift;                 // minimum L1 posterior shift for TR event acceptance (0=disabled)
  double shift_win_sec;             // pre/post window (seconds) for shift computation
  bool emit_annot_typed_classes;    // if true, encode directional transition type in point-event class name
  bool emit_all_trans;              // if true, pad TRANS/BOUT_TRANS pair summaries with unobserved directions
  double bout_min_sec;              // minimum stable-run span (seconds) to qualify as a bout
  double bout_gap_sec;               // max within-run gap (seconds) tolerated before splitting a bout
  double bout_smooth_sec;           // centered moving-average half-window (seconds) applied to PS/CONF before the bout stability test (0=off)
  std::string ch_valid;             // channel name for PP_INVALID invalidity flag (optional)
  std::string emit_sig_label;       // if non-empty, emit PS/H/CONF per-sample signals under this channel-name prefix
  bool do_emit_hypnogram;           // if true, emit hard-called (argmax) stage run annotations
  std::string hypno_prefix;         // optional prefix for those annotation classes, e.g. "hd_"
  std::string emit_spt_label;       // root for soft-SPT annotations; empty disables
  std::string emit_cycle_label;     // root for complete soft-ultradian cycle annotations; empty disables
  bool do_longrange;                // global long-range hypnodensity dynamics
  std::vector<double> lr_lags_sec;  // physical lags (seconds)
  std::vector<double> lr_windows_sec; // non-overlapping physical windows (seconds)
  double lr_min_duration_sec;       // minimum contiguous valid-core duration
  std::string lr_policy;             // gap-aware (default) or strict
  bool lr_fail_on_empty;             // retain an optional fatal empty-panel policy
  std::string analysis;              // spt (default), fixed, or recording
  double fixed_duration_sec;         // fixed analysis duration, anchored at sleep onset
  double analysis_flank_pre_sec;     // physical time retained before hypnodensity onset
  double analysis_flank_post_sec;    // physical time retained after hypnodensity offset
  std::string pc_out;                // optional level-1 HDSTATS feature directory
  std::string pc_proj;               // optional fixed reference PCA model
  double pc_min_loading;             // minimum absolute oriented loading for linked FTR output
  bool do_ultradian;                 // soft ultradian depth trajectory / phase
  double ultra_fwhm_sec;             // Gaussian FWHM in physical seconds
  double ultra_prominence;           // hysteretic peak/trough prominence
  double ultra_min_sleep_support;    // Gaussian-weighted local sleep mass [0,1]
  std::array<double,5> ultra_weights; // W,N1,N2,N3,R depth weights
  bool emit_ultradian;               // write U_PHASE/U_AMP EDF signals
  bool ultra_diagnostics;            // also write U_DEPTH/U_AMP_RAW EDF signals
  bool do_ultra_condition;           // joint time + harmonic phase descriptions
  bool ultra_condition_weighted;     // optional U_AMP-weighted sensitivity fit
  double ultra_max_condition;        // max standardized UVAR design condition number
  bool do_ultra_regions;             // optional four-region U_PHASE summaries
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
// Soft ultradian organization
//
// This is deliberately a turning-point coordinate, rather than a band-pass
// phase estimate.  U_DEPTH is the sleep-conditional Gaussian average of the
// configured posterior depth contrast.  A small state-machine applies a
// prominence (hysteresis) threshold to U_DEPTH and only commits an extremum
// once its reversal has been observed.  Thus neither recording edges nor an
// unsupported (usually wake) gap can manufacture a half-cycle.
// ============================================================

struct hd_ultradian_extremum_t
{
  int idx;
  bool trough;
};

struct hd_ultradian_t
{
  std::vector<double> depth;       // NaN where sleep support is inadequate
  std::vector<double> phase;       // NaN where no qualified half-cycle exists
  std::vector<double> amp;         // sleep-weighted excursion amplitude
  std::vector<double> raw_amp;     // unattenuated half-cycle amplitude
  std::vector<bool> supported;
  std::vector<hd_ultradian_extremum_t> extrema;
  int n_half_cycles = 0;
  int n_complete_cycles = 0;

  void compute( const hd_data_t & dat , const hd_params_t & par )
  {
    const int n = dat.N;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    depth.assign( n, nan );
    phase.assign( n, nan );
    amp.assign( n, 0.0 );
    raw_amp.assign( n, 0.0 );
    supported.assign( n, false );
    extrema.clear();
    n_half_cycles = n_complete_cycles = 0;
    if ( n == 0 ) return;

    // sigma = FWHM / (2 sqrt(2 log(2))); truncating at 3 sigma loses only
    // 0.27% of the nominal kernel mass while keeping the native-rate method
    // inexpensive for the 5 s/30 s posterior series for which HDSTATS is
    // intended.
    const double sigma_sec = par.ultra_fwhm_sec / 2.3548200450309493;
    const int radius = std::max( 1, (int)std::ceil( 3.0 * sigma_sec * dat.Fs ) );
    std::vector<double> kernel( radius + 1 );
    for ( int d = 0; d <= radius; ++d )
      {
        const double x = d / dat.Fs / sigma_sec;
        kernel[d] = std::exp( -0.5 * x * x );
      }

    for ( int i = 0; i < n; ++i )
      {
        double numerator = 0.0, sleep_mass = 0.0, available_mass = 0.0;
        const int lo = std::max( 0, i - radius );
        const int hi = std::min( n - 1, i + radius );
        for ( int j = lo; j <= hi; ++j )
          {
            if ( ! dat.valid[j] ) continue;
            const double k = kernel[ std::abs( j - i ) ];
            const double sleep = 1.0 - dat.p[0][j];
            const double d = par.ultra_weights[0] * dat.p[0][j]
              + par.ultra_weights[1] * dat.p[1][j]
              + par.ultra_weights[2] * dat.p[2][j]
              + par.ultra_weights[3] * dat.p[3][j]
              + par.ultra_weights[4] * dat.p[4][j];
            numerator += k * d;
            sleep_mass += k * sleep;
            available_mass += k;
          }

        // Normalising by available observations avoids boundary padding
        // artifacts.  The support ratio intentionally remains separate from
        // the depth itself: no sleep evidence is not intermediate depth.
        const double support = available_mass > 0.0 ? sleep_mass / available_mass : 0.0;
        if ( dat.valid[i] && support >= par.ultra_min_sleep_support && sleep_mass > 0.0 )
          {
            depth[i] = numerator / sleep_mass;
            supported[i] = std::isfinite( depth[i] );
          }
      }

    // Scan each contiguous supported interval separately.  This hysteretic
    // extrema finder is a prominence procedure: a trough/peak is accepted
    // only after a reversal of at least ultra_prominence.  It gives a single,
    // chronological, alternating sequence and deliberately does not add
    // unconfirmed extrema at either supported-run boundary.
    int run_lo = 0;
    while ( run_lo < n )
      {
        while ( run_lo < n && ! supported[run_lo] ) ++run_lo;
        if ( run_lo == n ) break;
        int run_hi = run_lo;
        while ( run_hi + 1 < n && supported[run_hi+1] ) ++run_hi;

        int low_i = run_lo, high_i = run_lo;
        double low = depth[run_lo], high = depth[run_lo];
        int direction = 0; // +1 rising / tracking a peak; -1 falling / trough
        std::vector<hd_ultradian_extremum_t> local;

        for ( int i = run_lo + 1; i <= run_hi; ++i )
          {
            const double x = depth[i];
            if ( direction == 0 )
              {
                if ( x < low ) { low = x; low_i = i; }
                if ( x > high ) { high = x; high_i = i; }
                if ( high - low >= par.ultra_prominence )
                  {
                    if ( low_i < high_i )
                      {
                        local.push_back( { low_i, true } );
                        direction = 1; // the rise has confirmed the trough
                        high = x; high_i = i;
                      }
                    else
                      {
                        local.push_back( { high_i, false } );
                        direction = -1; // the fall has confirmed the peak
                        low = x; low_i = i;
                      }
                  }
              }
            else if ( direction > 0 )
              {
                if ( x > high ) { high = x; high_i = i; }
                if ( high - x >= par.ultra_prominence )
                  {
                    local.push_back( { high_i, false } );
                    direction = -1;
                    low = x; low_i = i;
                  }
              }
            else
              {
                if ( x < low ) { low = x; low_i = i; }
                if ( x - low >= par.ultra_prominence )
                  {
                    local.push_back( { low_i, true } );
                    direction = 1;
                    high = x; high_i = i;
                  }
              }
          }

        // The state machine can use the first observation as a provisional
        // baseline while it determines direction.  It is not retained as a
        // detected extremum.  Edge portions are handled below as explicitly
        // partial half-cycles, rather than being silently discarded.
        if ( ! local.empty() && local.front().idx == run_lo ) local.erase( local.begin() );

        auto wrap_phase = []( double p ) {
          while ( p < 0.0 ) p += 2.0 * M_PI;
          while ( p >= 2.0 * M_PI ) p -= 2.0 * M_PI;
          return p;
        };
        auto assign_phase = [&]( const int lo, const int hi,
                                 const double phase_lo, const double phase_hi,
                                 const double excursion ) {
          if ( hi <= lo ) return;
          for ( int i = lo; i <= hi; ++i )
            {
              if ( ! supported[i] ) continue;
              const double f = ( i - lo ) / (double)( hi - lo );
              phase[i] = wrap_phase( phase_lo + f * ( phase_hi - phase_lo ) );
              raw_amp[i] = excursion;
              amp[i] = excursion * ( 1.0 - dat.p[0][i] );
            }
        };

        // A detected interior extremum can support the observed portion of a
        // half-cycle up to a supported-run edge.  The edge is *not* inserted
        // into extrema: its phase span is scaled by the observed depth change
        // relative to the neighbouring confirmed excursion (or prominence if
        // this is the sole extremum).  This keeps wake gaps unbridged while
        // retaining first/last partial trajectories instead of requiring an
        // unobservable reversal outside the recording.
        if ( ! local.empty() )
          {
            const hd_ultradian_extremum_t & first = local.front();
            if ( first.idx > run_lo )
              {
                const double observed = std::fabs( depth[first.idx] - depth[run_lo] );
                const double reference = local.size() > 1
                  ? std::fabs( depth[local[1].idx] - depth[first.idx] )
                  : std::max( par.ultra_prominence, observed );
                const double span = M_PI * std::min( 1.0, reference > 0.0 ? observed / reference : 0.0 );
                const double at_extremum = first.trough ? 0.0 : M_PI;
                assign_phase( run_lo, first.idx, at_extremum - span, at_extremum, observed );
              }
          }

        // Consecutive, confirmed extrema are the ordinary full half-cycles.
        for ( int k = 0; k + 1 < (int)local.size(); ++k )
          {
            const hd_ultradian_extremum_t & a = local[k];
            const hd_ultradian_extremum_t & b = local[k+1];
            if ( a.trough == b.trough || b.idx <= a.idx ) continue;
            const double excursion = std::fabs( depth[b.idx] - depth[a.idx] );
            assign_phase( a.idx, b.idx, a.trough ? 0.0 : M_PI,
                          a.trough ? M_PI : 2.0 * M_PI, excursion );
            ++n_half_cycles;
          }

        if ( ! local.empty() )
          {
            const hd_ultradian_extremum_t & last = local.back();
            if ( last.idx < run_hi )
              {
                const double observed = std::fabs( depth[run_hi] - depth[last.idx] );
                const double reference = local.size() > 1
                  ? std::fabs( depth[last.idx] - depth[local[local.size()-2].idx] )
                  : std::max( par.ultra_prominence, observed );
                const double span = M_PI * std::min( 1.0, reference > 0.0 ? observed / reference : 0.0 );
                const double at_extremum = last.trough ? 0.0 : M_PI;
                assign_phase( last.idx, run_hi, at_extremum,
                              last.trough ? at_extremum + span : at_extremum + span,
                              observed );
              }
          }

        // If a supported run contains a substantial monotonic excursion but
        // no interior reversal (e.g. the recording starts/ends mid-cycle),
        // retain that observed partial component.  No phase is assigned
        // outside the two observed depth endpoints.
        if ( local.empty() )
          {
            int lo_extreme = run_lo, hi_extreme = run_lo;
            for ( int i = run_lo + 1; i <= run_hi; ++i )
              {
                if ( depth[i] < depth[lo_extreme] ) lo_extreme = i;
                if ( depth[i] > depth[hi_extreme] ) hi_extreme = i;
              }
            const double excursion = depth[hi_extreme] - depth[lo_extreme];
            if ( excursion >= par.ultra_prominence && lo_extreme != hi_extreme )
              {
                if ( lo_extreme < hi_extreme ) assign_phase( lo_extreme, hi_extreme, 0.0, M_PI, excursion );
                else assign_phase( hi_extreme, lo_extreme, M_PI, 2.0 * M_PI, excursion );
                ++n_half_cycles;
              }
          }

        for ( int k = 0; k + 2 < (int)local.size(); ++k )
          if ( local[k].trough && ! local[k+1].trough && local[k+2].trough )
            ++n_complete_cycles;

        extrema.insert( extrema.end(), local.begin(), local.end() );
        run_lo = run_hi + 1;
      }
  }
};

// Apply the inverse of the class weights used to train a five-state model.
// The public option is deliberately in the same order as CALIBRATION:
// W,R,N1,N2,N3.  HDSTATS stores its channels as W,N1,N2,N3,R, hence the map.
static int hd_apply_weight_correction( hd_data_t * dat,
                                       const hd_params_t & par,
                                       const std::array<double,5> & weights )
{
  const int hdstats_to_weight[] = { 0, 2, 3, 4, 1 };
  int n_corrected = 0;

  for ( int i = 0; i < dat->N; ++i )
    {
      if ( ! dat->valid[i] ) continue;

      double denom = 0.0;
      for ( int k = 0; k < 5; ++k )
        denom += dat->p[k][i] / weights[ hdstats_to_weight[k] ];

      // load() has already rejected invalid posterior rows.  Retain this
      // guard so that malformed input can never be written back as a valid
      // posterior when update-pp=T.
      if ( ! std::isfinite( denom ) || denom <= 0.0 )
        {
          dat->valid[i] = false;
          for ( int k = 0; k < 5; ++k ) dat->p[k][i] = 0.0;
          continue;
        }

      for ( int k = 0; k < 5; ++k )
        dat->p[k][i] = ( dat->p[k][i] / weights[ hdstats_to_weight[k] ] ) / denom;
      ++n_corrected;
    }

  // load() creates these before the correction, so rebuild them afterwards.
  if ( par.do_3state )
    for ( int i = 0; i < dat->N; ++i )
      {
        dat->p3[0][i] = dat->p[0][i];
        dat->p3[1][i] = dat->p[1][i] + dat->p[2][i] + dat->p[3][i];
        dat->p3[2][i] = dat->p[4][i];
      }

  return n_corrected;
}

// Update only valid rows.  This intentionally preserves source values for
// PP_INVALID-marked, -9, or otherwise invalid rows, as well as the validity
// channel itself.  Use the common posterior [-9,1] physical range on write:
// an input channel can have an artificially narrow range (e.g. a nearly
// constant probability), which is not safe to retain after correction.
static void hd_update_posterior_channels( edf_t & edf,
                                          const hd_data_t & dat,
                                          const hd_params_t & par )
{
  const interval_t whole = edf.timeline.wholetrace();

  for ( int k = 0; k < 5; ++k )
    {
      signal_list_t sl = edf.header.signal_list( par.ch[k] );
      if ( sl.size() == 0 )
        Helper::halt( "HDSTATS update-pp: channel disappeared: " + par.ch[k] );

      const int slot = sl(0);
      slice_t s( edf, slot, whole );
      std::vector<double> d = *s.pdata();
      if ( (int)d.size() != dat.N )
        Helper::halt( "HDSTATS update-pp: channel length changed: " + par.ch[k] );

      for ( int i = 0; i < dat.N; ++i )
        if ( dat.valid[i] ) d[i] = dat.p[k][i];

      int16_t dmin = -32768, dmax = 32767;
      const double pmin = -9.0, pmax = 1.0;
      edf.update_signal( slot, &d, &dmin, &dmax, &pmin, &pmax );
    }
}


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

static std::vector<double> smooth_series( const std::vector<double> & x,
                                          int half_win,
                                          const std::vector<bool> & allowed )
{
  const int N = (int)x.size();
  if ( (int)allowed.size() != N ) return smooth_series( x, half_win );
  if ( half_win <= 0 || N == 0 ) return x;
  std::vector<double> out( N, 0.0 );
  for ( int i = 0; i < N; ++i )
    {
      if ( ! allowed[i] ) continue;
      double sum = 0.0;
      int n = 0;
      for ( int j = std::max( 0, i - half_win ); j <= std::min( N - 1, i + half_win ); ++j )
        if ( allowed[j] ) { sum += x[j]; ++n; }
      if ( n > 0 ) out[i] = sum / n;
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
    std::vector<bool> & is_stable,
    const std::vector<bool> * allowed = NULL )
{
  is_stable.assign( N, false );
  std::vector<bool> q( N, false );
  for ( int i = 0; i < N; i++ )
    q[i] = ( allowed == NULL || (*allowed)[i] ) && tv[i] < tv_th && conf[i] > conf_th;
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
      const double Fs,
      const std::vector<bool> * allowed = NULL )
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
        bool support_ok = true;
        if ( allowed != NULL )
          for ( int j = pre_lo; j <= post_hi; ++j )
            if ( ! (*allowed)[j] ) { support_ok = false; break; }
        if ( ! support_ok ) continue;

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
	       const hd_params_t & par,
               const std::vector<bool> * allowed = NULL )
  {
    const int N              = dat.N;
    const int win_smp        = std::max( 1, (int)std::round( par.window_sec       * dat.Fs ) );
    const int stable_smp     = std::max( 1, (int)std::round( par.stable_min_sec   * dat.Fs ) );
    const int bout_smooth_smp = std::max( 0, (int)std::round( par.bout_smooth_sec * dat.Fs ) );

    // 5-state detection
    {
      build_stable_mask( der.TV, der.C, N, stable_smp,
                         par.stable_ps_th, par.stable_conf_th, is_stable5, allowed );
      events5 = detect_boundary_events( dat.p, der.TV, is_stable5, der.argmax5, par, dat.Fs, allowed );
      build_trans_mask( events5, N, win_smp, is_trans5 );
      if ( allowed != NULL )
        for ( int i = 0; i < N; ++i ) if ( ! (*allowed)[i] ) is_trans5[i] = false;

      // TRANS takes priority: clear stable for any peri-event sample
      for ( int i = 0; i < N; i++ )
        if ( is_trans5[i] ) is_stable5[i] = false;

      // Bouts are formed on a separately smoothed PS/CONF stability test --
      // scoped to bout-formation only, so a single noisy sample doesn't reset
      // the whole stable-run clock the way it still does for IS_STABLE/TRANS.
      const std::vector<double> ps_smooth5   = allowed == NULL ? smooth_series( der.TV, bout_smooth_smp ) : smooth_series( der.TV, bout_smooth_smp, *allowed );
      const std::vector<double> conf_smooth5 = allowed == NULL ? smooth_series( der.C, bout_smooth_smp ) : smooth_series( der.C, bout_smooth_smp, *allowed );
      build_stable_mask( ps_smooth5, conf_smooth5, N, stable_smp,
                         par.stable_ps_th, par.stable_conf_th, is_bout_stable5, allowed );
      for ( int i = 0; i < N; i++ )
        if ( is_trans5[i] ) is_bout_stable5[i] = false;

      // BOUT_TRANS: per-sample halo (±window_sec), same construction as
      // IS_TRANS but anchored on bout-to-bout transitions instead of events.
      const std::vector<bool> all_mask = allowed == NULL ? std::vector<bool>( N, true ) : *allowed;
      const std::vector<hd_bout_t> bouts5 = compute_bouts( all_mask, is_bout_stable5, der.argmax5,
                                                             par.bout_min_sec, par.bout_gap_sec, dat.Fs );
      bout_events5 = compute_bout_transitions( bouts5, events5, der.TV );
      build_trans_mask( bout_events5, N, win_smp, is_bout_trans5 );
    }

    if ( ! par.do_3state ) return;

    // 3-state detection (independent from 5-state)
    {
      build_stable_mask( der.TV3, der.C3, N, stable_smp,
                         par.stable_ps_th, par.stable_conf_th, is_stable3, allowed );
      events3 = detect_boundary_events( dat.p3, der.TV3, is_stable3, der.argmax3, par, dat.Fs, allowed );
      build_trans_mask( events3, N, win_smp, is_trans3 );
      if ( allowed != NULL )
        for ( int i = 0; i < N; ++i ) if ( ! (*allowed)[i] ) is_trans3[i] = false;

      // TRANS takes priority: clear stable for any peri-event sample
      for ( int i = 0; i < N; i++ )
        if ( is_trans3[i] ) is_stable3[i] = false;

      const std::vector<double> ps_smooth3   = allowed == NULL ? smooth_series( der.TV3, bout_smooth_smp ) : smooth_series( der.TV3, bout_smooth_smp, *allowed );
      const std::vector<double> conf_smooth3 = allowed == NULL ? smooth_series( der.C3, bout_smooth_smp ) : smooth_series( der.C3, bout_smooth_smp, *allowed );
      build_stable_mask( ps_smooth3, conf_smooth3, N, stable_smp,
                         par.stable_ps_th, par.stable_conf_th, is_bout_stable3, allowed );
      for ( int i = 0; i < N; i++ )
        if ( is_trans3[i] ) is_bout_stable3[i] = false;

      const std::vector<bool> all_mask3 = allowed == NULL ? std::vector<bool>( N, true ) : *allowed;
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

// The analysis interval is a mask on the original clock, never a compacted
// signal.  SPT mode uses the existing sustained sleep onset/final-sleep logic
// and optionally retains physical-time flanks around that interval. Fixed mode
// is an exact-duration window whose first sample is the sustained sleep onset.
struct hd_analysis_t
{
  std::string mode = "recording";
  bool fallback_to_recording = false;
  bool fixed_ok = true;
  double fixed_available_sec = std::numeric_limits<double>::quiet_NaN();
  int lo = 0;
  int hi = -1;
  std::vector<bool> physical_mask;
  std::vector<bool> valid_mask;
  int n_invalid = 0;
};

static hd_analysis_t hd_make_analysis( const hd_data_t & dat,
                                        const hd_hypno_metrics_t & hm,
                                        const hd_params_t & par )
{
  hd_analysis_t out;
  out.hi = dat.N - 1;
  if ( par.analysis == "fixed" )
    {
      out.mode = "fixed";
      if ( ! hm.valid_onset )
        out.fixed_ok = false;
      else
        {
          const int duration = std::max( 1, (int)std::round( par.fixed_duration_sec * dat.Fs ) );
          out.lo = hm.onset_idx;
          out.fixed_available_sec = ( dat.N - out.lo ) / dat.Fs;
          out.hi = out.lo + duration - 1;
          if ( out.hi >= dat.N )
            {
              out.fixed_ok = false;
              // Keep the interval valid for any optional QC calculation; it
              // represents what is available, never a shortened analysis.
              out.hi = dat.N - 1;
            }
        }
    }
  else if ( par.analysis == "spt" )
    {
      if ( hm.valid_window )
        {
          const int pre = std::max( 0, (int)std::round( par.analysis_flank_pre_sec * dat.Fs ) );
          const int post = std::max( 0, (int)std::round( par.analysis_flank_post_sec * dat.Fs ) );
          out.mode = "spt";
          out.lo = std::max( 0, hm.onset_idx - pre );
          out.hi = std::min( dat.N - 1, hm.offset_idx + post );
        }
      else
        {
          out.fallback_to_recording = true;
          logger << "  warning: HDSTATS analysis=spt could not define a sustained sleep interval; using recording\n";
        }
    }
  out.physical_mask.assign( dat.N, false );
  out.valid_mask.assign( dat.N, false );
  for ( int i = out.lo; i <= out.hi; ++i )
    {
      out.physical_mask[i] = true;
      out.valid_mask[i] = dat.valid[i];
      if ( ! dat.valid[i] ) ++out.n_invalid;
    }
  return out;
}

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

static std::vector<bool> make_stage_mask( const std::vector<bool> &,
                                          const std::vector<int> &,
                                          const int );
static std::string hd_state_label( const bool, const int );

// Level-1 HDSTATS-PC features include the hard-coded structured panel inputs,
// including zero-filled directional densities. Time/SEC profiles are excluded.
// Keys encode their factors and use the selected analysis interval.
static void hd_pc_add( hdstats_pc::feature_map_t * out,
                       const std::string & key,
                       const double value )
{
  (*out)[key] = value;
}

// The ultradian PC block uses confidence and PS at both state resolutions.
// The ordinary HDSTATS tables retain a broader set of UQ/UVAR diagnostics.
// Every key is
// present even when ultradian phase cannot be estimated, so a PC model has a
// stable schema and treats such a record as unavailable rather than silently
// ascribing it zero phase modulation.
static const std::vector<std::string> & hd_pc_ultra_metrics()
{
  static const std::vector<std::string> metrics = { "CONF", "CONF3", "PS", "PS3" };
  return metrics;
}

static bool hd_pc_ultra_metric( const std::string & metric )
{
  const std::vector<std::string> & metrics = hd_pc_ultra_metrics();
  return std::find( metrics.begin(), metrics.end(), metric ) != metrics.end();
}

static void hd_pc_add_ultradian_placeholders( hdstats_pc::feature_map_t * out )
{
  const double na = std::numeric_limits<double>::quiet_NaN();
  for ( const std::string & metric : hd_pc_ultra_metrics() )
    {
      for ( const std::string & uq : std::vector<std::string>{ "D1", "D2", "L1", "L2" } )
        {
          hd_pc_add( out, hdstats_pc::key( "U_MEAN", "UQ/" + uq + ";UMETRIC/" + metric ), na );
        }
      for ( const std::string & statistic : std::vector<std::string>{
              "U_SLOPE", "U_MOD", "U_MOD_COS", "U_MOD_SIN", "U_R2", "U_PARTIAL_R2" } )
        hd_pc_add( out, hdstats_pc::key( statistic, "UVAR/" + metric ), na );
    }
}

static hdstats_pc::feature_map_t hd_pc_features(
    const hd_data_t & dat,
    const hd_derived_t & der,
    const hd_trans_t & tr,
    const hd_hypno_metrics_t & hm,
    const hd_analysis_t & analysis,
    const hd_params_t & par )
{
  hdstats_pc::feature_map_t out;
  const std::vector<bool> & mask = analysis.valid_mask;
  // Match compute_trans_shape(): count only events centered on valid analysis
  // samples and use those same samples for the density denominator.
  const double hours = std::count( mask.begin(), mask.end(), true ) / dat.Fs / 3600.0;
  auto event_density = [&]( const std::vector<hd_event_t> & events ) {
    int n = 0;
    for ( const auto & e : events )
      if ( e.idx >= 0 && e.idx < (int)mask.size() && mask[e.idx] ) ++n;
    return hours > 0.0 ? n / hours : 0.0;
  };

  const hd_region_stats_t rs5 = compute_region( dat, der, mask, par, false );
  // Structured panels use exactly the selected analysis interval and the
  // same statistics as the normal HDSTATS tables. Their schema lives in C++;
  // docs/hdstats/hdstats.pc is only a readable inventory.
  const double na = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> mins( 7, 0.0 ); // W N1 N2 N3 R NR S
  for ( int i = 0; i < dat.N; ++i )
    if ( mask[i] )
      {
        mins[0] += dat.p[0][i] / dat.Fs / 60.0;
        mins[1] += dat.p[1][i] / dat.Fs / 60.0;
        mins[2] += dat.p[2][i] / dat.Fs / 60.0;
        mins[3] += dat.p[3][i] / dat.Fs / 60.0;
        mins[4] += dat.p[4][i] / dat.Fs / 60.0;
      }
  mins[5] = mins[1] + mins[2] + mins[3]; mins[6] = mins[5] + mins[4];
  const std::vector<std::string> labels = { "W", "N1", "N2", "N3", "R", "NR", "S" };
  auto add = [&]( const std::string & var, const std::string & stratum, const double value ) {
    hd_pc_add( &out, hdstats_pc::key( var, stratum ), value );
  };
  const double analysis_min = ( analysis.hi - analysis.lo + 1 ) / dat.Fs / 60.0;
  const double sleep_min = mins[6];
  const bool selected_window = analysis.mode != "recording";
  for ( const bool three : { false, true } )
    {
      if ( three && !par.do_3state ) continue;
      const std::string s = three ? "NSS/3" : "NSS/5";
      const hd_region_stats_t & all = three ? compute_region( dat, der, mask, par, true ) : rs5;
      const std::vector<bool> & stable_mask = three ? tr.is_stable3 : tr.is_stable5;
      const std::vector<bool> & bout_stable = three ? tr.is_bout_stable3 : tr.is_bout_stable5;
      const std::vector<int> & argmax = three ? der.argmax3 : der.argmax5;
      const std::vector<hd_event_t> & events = three ? tr.events3 : tr.events5;
      const std::vector<hd_event_t> & bout_events = three ? tr.bout_events3 : tr.bout_events5;
      const bool all_ok = all.valid && all.n >= par.min_samples_region;
      auto masked = [&]( const bool want_stable ) {
        std::vector<bool> result( mask.size(), false );
        for ( int i = 0; i < (int)mask.size(); ++i ) result[i] = mask[i] && ( stable_mask[i] == want_stable );
        return result;
      };
      const hd_region_stats_t st = compute_region( dat, der, masked(true), par, three );
      const hd_region_stats_t un = compute_region( dat, der, masked(false), par, three );
      const bool st_ok = st.valid && st.n >= par.min_samples_region;
      const bool un_ok = un.valid && un.n >= par.min_samples_region;
      add( "H", s, all_ok ? all.mean_H : na );
      add( "H_P90", s, all_ok ? all.p90_H : na );
      add( "CONF", s, all_ok ? all.mean_C : na );
      add( "FRAC_LOW_CONF", s, all_ok ? all.frac_C_below : na );
      add( "MG", s, all_ok ? all.mean_Mg : na );
      add( "PS", s, all_ok ? all.mean_TV : na );
      add( "PS_P90", s, all_ok ? all.p90_TV : na );
      add( "PS_LAG", s, all_ok ? all.mean_TV_lag : na );
      add( "CORR_H_PS", s, all_ok ? all.corr_H_TV : na );
      add( "PCT_STABLE", s, all_ok ? 100.0 * st.n / all.n : na );
      add( "H_RATIO_UN_ST", s, st_ok && un_ok && st.mean_H > 0 ? un.mean_H / st.mean_H : na );
      add( "CONF_DIFF_UN_ST", s, st_ok && un_ok ? un.mean_C - st.mean_C : na );
      add( "PS_RATIO_UN_ST", s, st_ok && un_ok && st.mean_TV > 0 ? un.mean_TV / st.mean_TV : na );
      for ( const auto & pair : std::vector<std::pair<std::string,const hd_region_stats_t *>>{
            {"1",&st}, {"0",&un} } )
        {
          const bool ok = pair.first == "1" ? st_ok : un_ok;
          const std::string t = s + ";STABLE/" + pair.first;
          add( "H", t, ok ? pair.second->mean_H : na );
          add( "CONF", t, ok ? pair.second->mean_C : na );
          add( "PS", t, ok ? pair.second->mean_TV : na );
        }

      const std::vector<std::vector<double>> & P = three ? dat.p3 : dat.p;
      const std::vector<std::vector<int>> mixes = hd_mix_combinations( three );
      for ( const auto & combo : mixes )
        {
          std::string label = hd_state_label( three, combo[0] );
          for ( int i = 1; i < (int)combo.size(); ++i ) label += "_" + hd_state_label( three, combo[i] );
          const double mix_mins = compute_mix_mins( P, mask, combo, dat.Fs );
          add( "HD_MINS", s + ";MIX/" + label, mix_mins );
          const double mix_denom = selected_window ? analysis_min : hm.hd_spt;
          add( "HD_DENS", s + ";MIX/" + label, mix_denom > 0 ? mix_mins / mix_denom : 0.0 );
        }

      const auto bouts = compute_bouts( mask, bout_stable, argmax, par.bout_min_sec, par.bout_gap_sec, dat.Fs );
      const auto scoped_bouts = compute_bout_transitions( bouts, events, three ? der.TV3 : der.TV );
      auto emit_family = [&]( const std::vector<hd_event_t> & family, const std::string & prefix ) {
        hd_trans_stats_t stats;
        hd_profile_t profile;
        compute_trans_shape( der, family, mask, stable_mask, dat, par, dat.Fs, three,
                             false, false, stats, profile );
        add( prefix + "_N", s, stats.n_events );
        add( prefix + "_DENS", s, stats.density );
        add( prefix + "_DUR", s, stats.valid ? stats.mean_width_sec : na );
        add( prefix + "_PEAK_H", s, stats.valid ? stats.mean_peak_H : na );
        add( prefix + "_MIN_CONF", s, stats.valid ? stats.mean_min_C : na );
        add( prefix + "_TOT_PS", s, stats.valid ? stats.mean_TV_area : na );
        return stats.n_events;
      };
      const int n_trans = emit_family( events, "TRANS" );
      emit_family( scoped_bouts, "BOUT_TRANS" );
      std::vector<hd_event_t> scoped_events;
      for ( const auto & e : events ) if ( e.idx >= 0 && e.idx < (int)mask.size() && mask[e.idx] )
        scoped_events.push_back( e );
      const int n_matched = count_trans_bout_matches( scoped_events, bouts );
      add( "TRANS_N_MATCHED", s, n_matched );
      add( "TRANS_PCT_MATCHED", s, n_trans > 0 ? 100.0 * n_matched / n_trans : 0.0 );
      const int K = three ? 3 : 5;
      for ( int from = 0; from < K; ++from ) for ( int to = 0; to < K; ++to ) if ( from != to )
        {
          const std::string pair = hd_state_label( three, from ) + "to" + hd_state_label( three, to );
          auto pair_density = [&]( const std::vector<hd_event_t> & family, const bool directional ) {
            int count = 0;
            for ( const auto & e : family )
              if ( ( !directional || e.is_dir ) && e.from_st == from && e.to_st == to &&
                   e.idx >= 0 && e.idx < (int)mask.size() && mask[e.idx] ) ++count;
            return hours > 0 ? count / hours : 0.0;
          };
          add( "TRANS_DENS", s + ";TRANS/" + pair, pair_density( events, true ) );
          add( "BOUT_TRANS_DENS", s + ";BOUT_TRANS/" + pair, pair_density( scoped_bouts, false ) );
        }
    }

  if ( hm.valid_onset ) add( "HD_SOL", ".", hm.hd_sol );
  add( "HD_SPT", ".", selected_window ? analysis_min : ( hm.valid_window ? hm.hd_spt : na ) );
  add( "HD_WASO", ".", selected_window ? mins[0] : ( hm.valid_window ? hm.hd_waso : na ) );
  add( "HD_SME", ".", selected_window ? ( analysis_min > 0 ? 100.0 * sleep_min / analysis_min : 0.0 )
                          : ( hm.valid_window ? hm.hd_sme : na ) );
  const bool rem_in_analysis = hm.valid_rem &&
    ( ! selected_window || ( hm.rem_idx >= analysis.lo && hm.rem_idx <= analysis.hi ) );
  add( "HD_REM_LAT", ".", rem_in_analysis ? hm.hd_rem_lat : na );
  add( "HD_REM_LAT2", ".", rem_in_analysis ? hm.hd_rem_lat2 : na );
  for ( int k = 0; k < 5; ++k )
    {
      const std::string st = labels[k], stratum = "SS/" + st;
      const double value = selected_window ? mins[k] : hm.mins.at( st );
      const double denom = selected_window ? ( k == 0 ? analysis_min : sleep_min )
                               : ( k == 0 ? hm.total_min : hm.hd_tst );
      add( "HD_MINS", stratum, value );
      add( "HD_PCT", stratum, denom > 0 ? 100.0 * value / denom : 0.0 );
      add( "HD_DENS", stratum, ( selected_window ? analysis_min : hm.hd_spt ) > 0
           ? value / ( selected_window ? analysis_min : hm.hd_spt ) : 0.0 );
    }
  for ( const auto & stage : std::vector<std::pair<bool,int>>{
          { true, 0 }, { false, 2 }, { true, 2 } } )
    if ( stage.first ? par.do_3state : true )
      {
        const bool three = stage.first;
        const int k = stage.second;
        const std::string label = hd_state_label( three, k );
        const std::vector<bool> stage_mask = make_stage_mask(
          mask, three ? der.argmax3 : der.argmax5, k );
        const hd_region_stats_t region = compute_region( dat, der, stage_mask, par, three );
        const bool ok = region.valid && region.n >= par.min_samples_region;
        const std::string s = "SS/" + label + ";NSS/" + ( three ? "3" : "5" );
        add( "P_" + label, s, ok ? region.mean_p[k] : na );
        add( "PS", s, ok ? region.mean_TV : na );
      }
  hd_pc_add_ultradian_placeholders( &out );
  return out;
}

static void write_trans_stats( const hd_trans_stats_t & ts,
                               const std::string & prefix,
                               const bool emit_na_for_empty = false )
{
  writer.value( prefix + "_N",    ts.n_events );
  writer.value( prefix + "_DENS", ts.density );
  if ( ts.valid )
    {
      writer.value( prefix + "_DUR",      ts.mean_width_sec );
      writer.value( prefix + "_PEAK_H",   ts.mean_peak_H );
      writer.value( prefix + "_MIN_CONF", ts.mean_min_C );
      writer.value( prefix + "_TOT_PS",   ts.mean_TV_area );
    }
  else if ( emit_na_for_empty )
    {
      const double na = std::numeric_limits<double>::quiet_NaN();
      writer.value( prefix + "_DUR",      na );
      writer.value( prefix + "_PEAK_H",   na );
      writer.value( prefix + "_MIN_CONF", na );
      writer.value( prefix + "_TOT_PS",   na );
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

// Missing phase must not be stored as zero: zero is the trough convention.
// EDF itself has no NA representation, so U_PHASE uses -1 radians and the
// companion U_PHASE_VALID channel makes the convention explicit.  U_DEPTH's
// diagnostic sentinel is chosen below the observed depth range.
static void hd_emit_ultradian_signals( edf_t & edf,
                                       const hd_data_t & dat,
                                       const hd_ultradian_t & ultra,
                                       const hd_params_t & par )
{
  auto add_sig = [&]( const std::string & label, const std::vector<double> & x,
                      const double pmin, const double pmax )
    {
      if ( edf.header.has_signal( label ) )
        {
          logger << "  HDSTATS emit-ultradian: " << label << " already exists, skipping\n";
          return;
        }
      edf.add_signal( label, dat.Fs, x, pmin, pmax );
    };

  std::vector<double> p( dat.N, -1.0 ), pvalid( dat.N, 0.0 );
  std::vector<double> ua( dat.N, 0.0 ), ura( dat.N, 0.0 );
  double depth_min = std::numeric_limits<double>::infinity();
  double depth_max = -std::numeric_limits<double>::infinity();
  for ( int i = 0; i < dat.N; ++i )
    {
      if ( std::isfinite( ultra.phase[i] ) )
        { p[i] = ultra.phase[i]; pvalid[i] = 1.0; }
      ua[i] = ultra.amp[i];
      ura[i] = ultra.raw_amp[i];
      if ( std::isfinite( ultra.depth[i] ) )
        {
          depth_min = std::min( depth_min, ultra.depth[i] );
          depth_max = std::max( depth_max, ultra.depth[i] );
        }
    }

  add_sig( "U_PHASE", p, -1.0, 2.0 * M_PI );
  add_sig( "U_PHASE_VALID", pvalid, 0.0, 1.0 );
  double amp_max = 1.0;
  for ( int i = 0; i < dat.N; ++i ) amp_max = std::max( amp_max, std::max( ua[i], ura[i] ) );
  add_sig( "U_AMP", ua, 0.0, amp_max );

  if ( par.ultra_diagnostics )
    {
      const double sentinel = std::isfinite( depth_min )
        ? depth_min - std::max( 1.0, depth_max - depth_min + 1.0 ) : -9.0;
      std::vector<double> ud( dat.N, sentinel );
      for ( int i = 0; i < dat.N; ++i )
        if ( std::isfinite( ultra.depth[i] ) ) ud[i] = ultra.depth[i];
      add_sig( "U_DEPTH", ud, sentinel, std::isfinite( depth_max ) ? depth_max : 1.0 );
      add_sig( "U_AMP_RAW", ura, 0.0, amp_max );
    }

  logger << "  HDSTATS emit-ultradian: added U_PHASE, U_PHASE_VALID and U_AMP"
         << ( par.ultra_diagnostics ? ", U_DEPTH and U_AMP_RAW" : "" ) << "\n";
}

// These are intentionally labelled as HDSTATS soft-ultradian cycles, rather
// than canonical scored NREM cycles. A cycle requires three consecutive,
// confirmed extrema (trough -> peak -> trough) within one contiguous,
// selected, valid supported run. Thus it neither invents partial cycles at
// record/window edges nor bridges invalid/wake-supported gaps.
static int hd_emit_ultradian_cycles( edf_t & edf,
                                     const hd_data_t & dat,
                                     const hd_ultradian_t & ultra,
                                     const hd_analysis_t & analysis,
                                     const hd_params_t & par )
{
  if ( par.emit_cycle_label.empty() || dat.tp.empty() ) return 0;

  const uint64_t sample_tp = (uint64_t)std::round( (double)globals::tp_1sec / dat.Fs );
  const std::string class_name = par.emit_cycle_label + "_cycle";
  int n_emitted = 0;
  for ( int k = 0; k + 2 < (int)ultra.extrema.size(); ++k )
    {
      const hd_ultradian_extremum_t & a = ultra.extrema[k];
      const hd_ultradian_extremum_t & b = ultra.extrema[k+1];
      const hd_ultradian_extremum_t & c = ultra.extrema[k+2];
      if ( ! a.trough || b.trough || ! c.trough || a.idx >= b.idx || b.idx >= c.idx ) continue;
      if ( a.idx < analysis.lo || c.idx > analysis.hi ) continue;

      bool contiguous_selected_support = true;
      for ( int i = a.idx; i <= c.idx; ++i )
        if ( ! analysis.valid_mask[i] || ! ultra.supported[i] )
          { contiguous_selected_support = false; break; }
      if ( ! contiguous_selected_support ) continue;

      ++n_emitted;
      edf.annotations->add( class_name )->add(
        "C" + std::to_string( n_emitted ),
        interval_t( dat.tp[a.idx], dat.tp[c.idx] + sample_tp ), "." );
    }

  return n_emitted;
}

struct hd_ultra_fit_t
{
  bool ok = false;
  std::string status = "NO_PHASE";
  std::vector<double> beta;
  double sse = std::numeric_limits<double>::quiet_NaN();
  double design_condition = std::numeric_limits<double>::quiet_NaN();
};

// Condition number of the weighted design after centering and scaling every
// non-intercept column. This makes the QC invariant to units (hours versus
// radians) and tests the real question: can time, cos(phase), and sin(phase)
// be distinguished over the observed phase support?
static double hd_ultra_standardized_condition( const std::vector<std::vector<double>> & x,
                                                const std::vector<double> & w )
{
  const int n = (int)x.size();
  const int p = n == 0 ? 0 : (int)x[0].size();
  if ( p == 0 || (int)w.size() != n ) return std::numeric_limits<double>::infinity();
  double sum_w = 0.0;
  std::vector<double> mean( p, 0.0 ), ss( p, 0.0 );
  for ( int i = 0; i < n; ++i )
    {
      if ( !std::isfinite(w[i]) || w[i] <= 0.0 || (int)x[i].size() != p )
        return std::numeric_limits<double>::infinity();
      sum_w += w[i];
      for ( int j = 1; j < p; ++j )
        {
          if ( !std::isfinite(x[i][j]) ) return std::numeric_limits<double>::infinity();
          mean[j] += w[i] * x[i][j];
        }
    }
  if ( !std::isfinite(sum_w) || sum_w <= 0.0 ) return std::numeric_limits<double>::infinity();
  for ( int j = 1; j < p; ++j ) mean[j] /= sum_w;
  for ( int i = 0; i < n; ++i )
    for ( int j = 1; j < p; ++j )
      { const double d = x[i][j] - mean[j]; ss[j] += w[i] * d * d; }
  for ( int j = 1; j < p; ++j )
    if ( !std::isfinite(ss[j]) || ss[j] <= 1e-14 * sum_w )
      return std::numeric_limits<double>::infinity();

  std::vector<std::vector<double>> g( p, std::vector<double>( p, 0.0 ) );
  for ( int i = 0; i < n; ++i )
    {
      std::vector<double> z( p, 0.0 );
      z[0] = 1.0 / std::sqrt(sum_w);
      for ( int j = 1; j < p; ++j ) z[j] = ( x[i][j] - mean[j] ) / std::sqrt(ss[j]);
      for ( int j = 0; j < p; ++j )
        for ( int k = 0; k <= j; ++k ) g[j][k] += w[i] * z[j] * z[k];
    }
  for ( int j = 0; j < p; ++j ) for ( int k = j + 1; k < p; ++k ) g[j][k] = g[k][j];

  // Jacobi eigensolver is ample for this fixed 4-by-4 Gram matrix.
  for ( int iter = 0; iter < 100; ++iter )
    {
      int a = 0, b = 1;
      double largest = 0.0;
      for ( int j = 0; j < p; ++j ) for ( int k = j + 1; k < p; ++k )
        if ( std::fabs(g[j][k]) > largest ) { largest = std::fabs(g[j][k]); a = j; b = k; }
      if ( largest < 1e-12 ) break;
      const double theta = 0.5 * std::atan2( 2.0 * g[a][b], g[b][b] - g[a][a] );
      const double c = std::cos(theta), s = std::sin(theta);
      const double gaa = g[a][a], gbb = g[b][b], gab = g[a][b];
      g[a][a] = c*c*gaa - 2.0*s*c*gab + s*s*gbb;
      g[b][b] = s*s*gaa + 2.0*s*c*gab + c*c*gbb;
      g[a][b] = g[b][a] = 0.0;
      for ( int j = 0; j < p; ++j )
        if ( j != a && j != b )
          {
            const double gja = g[j][a], gjb = g[j][b];
            g[j][a] = g[a][j] = c*gja - s*gjb;
            g[j][b] = g[b][j] = s*gja + c*gjb;
          }
    }
  double lo = g[0][0], hi = g[0][0];
  for ( int j = 1; j < p; ++j ) { lo = std::min(lo,g[j][j]); hi = std::max(hi,g[j][j]); }
  return lo > 0.0 && std::isfinite(lo) && std::isfinite(hi)
    ? std::sqrt(hi / lo) : std::numeric_limits<double>::infinity();
}

// Weighted least squares for a deliberately tiny design matrix.  Numerical
// rank alone is insufficient: a partial monotonic phase ramp can be full-rank
// yet almost indistinguishable from elapsed time.
static hd_ultra_fit_t hd_ultra_fit_linear( const std::vector<std::vector<double>> & x,
                                           const std::vector<double> & y,
                                           const std::vector<double> & w,
                                           const double max_condition )
{
  hd_ultra_fit_t out;
  const int n = (int)y.size();
  const int p = x.empty() ? 0 : (int)x[0].size();
  if ( n < 2 * p || p == 0 ) { out.status = "INSUFFICIENT_N"; return out; }
  out.design_condition = hd_ultra_standardized_condition( x, w );
  if ( !std::isfinite(out.design_condition) || out.design_condition > max_condition )
    { out.status = "ILL_CONDITIONED"; return out; }
  std::vector<std::vector<double>> a( p, std::vector<double>( p + 1, 0.0 ) );
  for ( int i = 0; i < n; ++i )
    {
      if ( ! std::isfinite(y[i]) || !std::isfinite(w[i]) || w[i] <= 0.0 )
        { out.status = "INVALID_INPUT"; return out; }
      for ( int j = 0; j < p; ++j )
        {
          if ( !std::isfinite(x[i][j]) ) { out.status = "INVALID_INPUT"; return out; }
          a[j][p] += w[i] * x[i][j] * y[i];
          for ( int k = 0; k < p; ++k ) a[j][k] += w[i] * x[i][j] * x[i][k];
        }
    }

  double max_pivot = 0.0;
  for ( int j = 0; j < p; ++j )
    for ( int k = 0; k < p; ++k ) max_pivot = std::max( max_pivot, std::fabs(a[j][k]) );
  if ( max_pivot <= 0.0 || !std::isfinite(max_pivot) ) { out.status = "SINGULAR"; return out; }

  for ( int c = 0; c < p; ++c )
    {
      int pivot = c;
      for ( int r = c + 1; r < p; ++r )
        if ( std::fabs(a[r][c]) > std::fabs(a[pivot][c]) ) pivot = r;
      const double pv = std::fabs(a[pivot][c]);
      // This threshold is intentionally conservative.  It is a design-QC
      // rule, not an inferential p-value or a required number of cycles.
      if ( pv / max_pivot < 1e-8 || !std::isfinite(pv) ) { out.status = "ILL_CONDITIONED"; return out; }
      if ( pivot != c ) std::swap( a[pivot], a[c] );
      const double d = a[c][c];
      for ( int j = c; j <= p; ++j ) a[c][j] /= d;
      for ( int r = 0; r < p; ++r )
        if ( r != c )
          {
            const double f = a[r][c];
            for ( int j = c; j <= p; ++j ) a[r][j] -= f * a[c][j];
          }
    }

  out.beta.assign( p, 0.0 );
  for ( int j = 0; j < p; ++j ) out.beta[j] = a[j][p];
  out.sse = 0.0;
  for ( int i = 0; i < n; ++i )
    {
      double fitted = 0.0;
      for ( int j = 0; j < p; ++j ) fitted += x[i][j] * out.beta[j];
      const double e = y[i] - fitted;
      out.sse += w[i] * e * e;
    }
  out.ok = std::isfinite(out.sse);
  out.status = out.ok ? "OK" : "NUMERIC_FAILURE";
  return out;
}

// Per-metric joint Y ~ time + cos(U_PHASE) + sin(U_PHASE) descriptions.  These
// are global/analysis-interval outputs, analogous to the existing long-range
// trend summaries, rather than thousands of heavily overlapping strata.
static void hd_write_ultra_conditioned_metrics( const hd_data_t & dat,
                                                 const hd_derived_t & der,
                                                 const hd_ultradian_t & ultra,
                                                 const hd_analysis_t & analysis,
                                                 const hd_params_t & par,
                                                 hdstats_pc::feature_map_t * pc_features = NULL,
                                                 const bool emit_output = true )
{
  struct metric_t { const char * label; const std::vector<double> * x; int gap; };
  std::vector<metric_t> metrics = {
    { "H", &der.H, 0 }, { "CONF", &der.C, 0 }, { "MG", &der.Mg, 0 },
    { "PS", &der.TV, 1 }, { "PS_LAG", &der.TV_lag,
      std::max(1, (int)std::round(par.lag_sec * dat.Fs)) }
  };
  if ( par.do_3state )
    {
      metrics.push_back( { "H3", &der.H3, 0 } );
      metrics.push_back( { "CONF3", &der.C3, 0 } );
      metrics.push_back( { "MG3", &der.Mg3, 0 } );
      metrics.push_back( { "PS3", &der.TV3, 1 } );
      metrics.push_back( { "PS3_LAG", &der.TV3_lag,
                           std::max(1, (int)std::round(par.lag_sec * dat.Fs)) } );
    }

  const double nan = std::numeric_limits<double>::quiet_NaN();
  for ( const auto & metric : metrics )
    {
      std::vector<int> idx;
      std::vector<double> raw_w;
      for ( int i = analysis.lo; i <= analysis.hi; ++i )
        {
          if ( !analysis.valid_mask[i] || !std::isfinite(ultra.phase[i])
               || !std::isfinite((*(metric.x))[i]) ) continue;
          if ( metric.gap > 0 && ( i < metric.gap || !dat.valid[i-metric.gap] ) ) continue;
          const double wi = par.ultra_condition_weighted ? ultra.amp[i] : 1.0;
          if ( !std::isfinite(wi) || wi <= 0.0 ) continue;
          idx.push_back(i); raw_w.push_back(wi);
        }

      std::string status = "NO_PHASE";
      double slope = nan, modulation = nan, peak = nan, r2 = nan, partial_r2 = nan, design_condition = nan;
      double phase_cos = nan, phase_sin = nan;
      if ( !idx.empty() )
        {
          double mean_t = 0.0, sum_w = 0.0;
          for ( int q = 0; q < (int)idx.size(); ++q )
            { mean_t += raw_w[q] * idx[q] / dat.Fs / 3600.0; sum_w += raw_w[q]; }
          mean_t = sum_w > 0.0 ? mean_t / sum_w : 0.0;
          // Normalise positive amplitude weights to unit mean: numerical
          // scaling changes neither coefficient estimates nor R2, and makes
          // the optional weighting definition explicit.
          const double mean_w = sum_w / idx.size();
          std::vector<double> y, w;
          std::vector<std::vector<double>> x_time, x_full;
          double ybar = 0.0;
          for ( int q = 0; q < (int)idx.size(); ++q )
            {
              const int i = idx[q];
              const double wi = raw_w[q] / mean_w;
              const double t = i / dat.Fs / 3600.0 - mean_t;
              y.push_back( (*(metric.x))[i] ); w.push_back(wi);
              ybar += wi * y.back();
              x_time.push_back( { 1.0, t } );
              x_full.push_back( { 1.0, t, std::cos(ultra.phase[i]), std::sin(ultra.phase[i]) } );
            }
          double wsum = std::accumulate(w.begin(), w.end(), 0.0);
          ybar = wsum > 0.0 ? ybar / wsum : 0.0;
          double sst = 0.0;
          for ( int q = 0; q < (int)y.size(); ++q ) sst += w[q] * (y[q]-ybar) * (y[q]-ybar);
          hd_ultra_fit_t time_fit = hd_ultra_fit_linear( x_time, y, w, par.ultra_max_condition );
          hd_ultra_fit_t full_fit = hd_ultra_fit_linear( x_full, y, w, par.ultra_max_condition );
          design_condition = full_fit.design_condition;
          status = full_fit.status;
          if ( full_fit.ok && time_fit.ok && sst > 1e-12 && time_fit.sse > 1e-12 )
            {
              slope = full_fit.beta[1];
              phase_cos = full_fit.beta[2];
              phase_sin = full_fit.beta[3];
              modulation = std::sqrt( full_fit.beta[2]*full_fit.beta[2]
                                      + full_fit.beta[3]*full_fit.beta[3] );
              peak = std::atan2( full_fit.beta[3], full_fit.beta[2] );
              if ( peak < 0.0 ) peak += 2.0 * M_PI;
              r2 = 1.0 - full_fit.sse / sst;
              partial_r2 = ( time_fit.sse - full_fit.sse ) / time_fit.sse;
              if ( partial_r2 < 0.0 && partial_r2 > -1e-10 ) partial_r2 = 0.0;
              status = "OK";
            }
          else if ( full_fit.ok && time_fit.ok && sst <= 1e-12 ) status = "ZERO_VARIANCE";
        }

      if ( emit_output )
        {
          writer.level( metric.label, "UVAR" );
          writer.value( "U_SLOPE", slope );
          writer.value( "U_MOD", modulation );
          writer.value( "U_PEAK", peak );
          writer.value( "U_R2", r2 );
          writer.value( "U_PARTIAL_R2", partial_r2 );
          if ( par.qc )
            {
              writer.level( 1, "QC" );
              writer.value( "U_DESIGN_COND", design_condition );
              writer.value( "U_STATUS", status );
              writer.unlevel( "QC" );
            }
          writer.unlevel( "UVAR" );
        }

      // A raw circular peak cannot be a sensible scalar PC feature: phases
      // just below 2*pi and just above zero are adjacent physiologically but
      // far apart numerically.  The two signed harmonic coefficients retain
      // both U_MOD and U_PEAK without that discontinuity.  QC failures remain
      // missing rather than being recoded as no modulation.
      if ( pc_features != NULL && hd_pc_ultra_metric( metric.label ) && status == "OK" )
        {
          const std::string s = "UVAR/" + std::string( metric.label );
          hd_pc_add( pc_features, hdstats_pc::key( "U_SLOPE", s ), slope );
          hd_pc_add( pc_features, hdstats_pc::key( "U_MOD", s ), modulation );
          hd_pc_add( pc_features, hdstats_pc::key( "U_MOD_COS", s ), phase_cos );
          hd_pc_add( pc_features, hdstats_pc::key( "U_MOD_SIN", s ), phase_sin );
          hd_pc_add( pc_features, hdstats_pc::key( "U_R2", s ), r2 );
          hd_pc_add( pc_features, hdstats_pc::key( "U_PARTIAL_R2", s ), partial_r2 );
        }
    }
}

// Complementary descriptive view of the continuous harmonic fit.  These are
// phase regions, not scored sleep-cycle states: each observation remains on
// the same soft U_PHASE coordinate used by UVAR.
static void hd_write_ultra_phase_regions( const hd_data_t & dat,
                                           const hd_derived_t & der,
                                           const hd_ultradian_t & ultra,
                                           const hd_analysis_t & analysis,
                                           const hd_params_t & par,
                                           hdstats_pc::feature_map_t * pc_features = NULL,
                                           const bool emit_output = true )
{
  struct metric_t { const char * label; const std::vector<double> * x; int gap; };
  std::vector<metric_t> metrics = {
    { "H", &der.H, 0 }, { "CONF", &der.C, 0 }, { "MG", &der.Mg, 0 },
    { "PS", &der.TV, 1 }, { "PS_LAG", &der.TV_lag,
      std::max(1, (int)std::round(par.lag_sec * dat.Fs)) }
  };
  if ( par.do_3state )
    {
      metrics.push_back( { "H3", &der.H3, 0 } );
      metrics.push_back( { "CONF3", &der.C3, 0 } );
      metrics.push_back( { "MG3", &der.Mg3, 0 } );
      metrics.push_back( { "PS3", &der.TV3, 1 } );
      metrics.push_back( { "PS3_LAG", &der.TV3_lag,
                           std::max(1, (int)std::round(par.lag_sec * dat.Fs)) } );
    }
  static const char * region_label[] = { "D1", "D2", "L1", "L2" };

  for ( const auto & metric : metrics )
    {
      double sum[4] = { 0.0, 0.0, 0.0, 0.0 };
      int n[4] = { 0, 0, 0, 0 };
      for ( int i = analysis.lo; i <= analysis.hi; ++i )
        {
          if ( !analysis.valid_mask[i] || !std::isfinite(ultra.phase[i])
               || !std::isfinite((*(metric.x))[i]) ) continue;
          if ( metric.gap > 0 && ( i < metric.gap || !dat.valid[i-metric.gap] ) ) continue;
          int q = (int)std::floor( ultra.phase[i] / ( M_PI / 2.0 ) );
          if ( q < 0 ) continue;
          if ( q > 3 ) q = 3; // defensively accepts round-off just below 2 pi
          sum[q] += (*(metric.x))[i];
          ++n[q];
        }
      for ( int q = 0; q < 4; ++q )
        if ( n[q] >= par.min_samples_region )
          {
            const double mean = sum[q] / n[q];
            if ( emit_output )
              {
                writer.level( region_label[q], "UQ" );
                writer.level( metric.label, "UMETRIC" );
                writer.value( "U_MEAN", mean );
                writer.unlevel( "UMETRIC" );
                writer.unlevel( "UQ" );
              }
            if ( pc_features != NULL && hd_pc_ultra_metric( metric.label ) )
              {
                hd_pc_add( pc_features, hdstats_pc::key( "U_MEAN", std::string("UQ/") + region_label[q] + ";UMETRIC/" + metric.label ), mean );
              }
          }
    }
}

static void write_hd_hypno_metrics( const hd_hypno_metrics_t & hm,
                                    const hd_data_t & dat,
                                    const hd_analysis_t & analysis )
{
  if ( hm.valid_onset ) writer.value( "HD_SOL" , hm.hd_sol );
  if ( analysis.mode != "recording" )
    {
      const double dt_min = 1.0 / dat.Fs / 60.0;
      const double analysis_min = ( analysis.hi - analysis.lo + 1 ) * dt_min;
      double wake_min = 0.0, sleep_min = 0.0;
      std::map<std::string,double> mins;
      for ( int i = analysis.lo; i <= analysis.hi; ++i )
        if ( dat.valid[i] )
          {
            const double w = dat.p[0][i];
            const double nr = dat.p[1][i] + dat.p[2][i] + dat.p[3][i];
            const double r = dat.p[4][i];
            wake_min += w * dt_min;
            sleep_min += ( nr + r ) * dt_min;
            mins["W"] += w * dt_min; mins["N1"] += dat.p[1][i] * dt_min;
            mins["N2"] += dat.p[2][i] * dt_min; mins["N3"] += dat.p[3][i] * dt_min;
            mins["R"] += r * dt_min; mins["NR"] += nr * dt_min; mins["S"] += ( nr + r ) * dt_min;
          }
      writer.value( "HD_SPT", analysis_min );
      writer.value( "HD_WASO", wake_min );
      writer.value( "HD_SME", analysis_min > 0.0 ? 100.0 * sleep_min / analysis_min : 0.0 );
      if ( hm.valid_rem && hm.rem_idx >= analysis.lo && hm.rem_idx <= analysis.hi )
        {
          writer.value( "HD_REM_LAT" , hm.hd_rem_lat );
          writer.value( "HD_REM_LAT2" , hm.hd_rem_lat2 );
        }
      static const std::vector<std::string> stages = { "W", "N1", "N2", "N3", "R", "NR", "S" };
      for ( const auto & st : stages )
        {
          writer.level( st, "SS" );
          const double value = mins[st];
          writer.value( "HD_MINS", value );
          const bool use_total = st == "W" || st == "S";
          const double denom = use_total ? analysis_min : sleep_min;
          writer.value( "HD_PCT", denom > 0.0 ? 100.0 * value / denom : 0.0 );
          writer.value( "HD_DENS", analysis_min > 0.0 ? value / analysis_min : 0.0 );
          writer.unlevel( "SS" );
        }
      return;
    }
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
    const hd_hypno_metrics_t * hm,
    const hd_analysis_t * analysis )
{
  if ( hm != NULL && par.do_hd_metrics )
    write_hd_hypno_metrics( *hm, dat, *analysis );

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
	      writer.value( "N", pf.n_events );
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
	  hd_trans_stats_t tshape_bout;
	  hd_profile_t     profile_bout;
	  if ( emit_top_level_transition_stats || emit_profiles_and_pairs )
	    compute_trans_shape( der, bout_events, context_mask, is_stable, dat, par, dat.Fs, threestate,
	                         false, false, tshape_bout, profile_bout );

	  // TRANS events restricted to this context (matches how TRANS_N itself, via
	  // compute_trans_shape's region_mask filtering, is scoped).
	  std::vector<hd_event_t> evts_in_context;
	  for ( const auto & e : evts )
	    if ( e.idx >= 0 && e.idx < (int)context_mask.size() && context_mask[e.idx] )
	      evts_in_context.push_back( e );

	  if ( emit_top_level_transition_stats )
	    {
	      write_trans_stats( tshape, "TRANS" );

	      const int n_matched = count_trans_bout_matches( evts_in_context, bouts );
	      writer.value( "TRANS_N_MATCHED",   n_matched );
	      writer.value( "TRANS_PCT_MATCHED", tshape.n_events > 0 ? 100.0 * n_matched / tshape.n_events : 0.0 );

	      write_trans_stats( tshape_bout, "BOUT_TRANS" );

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
		{
		  const double denom = analysis->mode != "recording"
		    ? ( analysis->hi - analysis->lo + 1 ) / dat.Fs / 60.0
		    : hm->hd_spt;
		  writer.value( "HD_DENS", denom > 0 ? mix_mins / denom : 0.0 );
		}
		    writer.unlevel( "MIX" );
		  }
	      }
	    }

	  if ( ! emit_profiles_and_pairs ) return;

	  // Generic aligned profiles, pooled across every event family.  The ALL
	  // level keeps the family explicit while retaining one table per family.
	  writer.level( "ALL", "TRANS" );
	  if ( profile.valid )
	    emit_profile_rows( profile, false );
	  writer.unlevel( "TRANS" );

	  writer.level( "ALL", "BOUT_TRANS" );
	  if ( profile_bout.valid )
	    emit_profile_rows( profile_bout, false );
	  writer.unlevel( "BOUT_TRANS" );

	  // TRANS: pair-specific summary statistics still use directional events (unchanged).
	  std::set<std::pair<int,int>> trans_pairs;
	  for ( const auto & e : evts )
	    if ( e.is_dir && e.from_st >= 0 && e.to_st >= 0 && e.from_st != e.to_st )
	      trans_pairs.insert( std::make_pair( e.from_st , e.to_st ) );
	  if ( par.emit_all_trans )
	    for ( int from = 0; from < K; from++ )
	      for ( int to = 0; to < K; to++ )
		if ( from != to ) trans_pairs.insert( std::make_pair( from, to ) );

	  for ( const auto & tp : trans_pairs )
	    {
	      std::vector<hd_event_t> evts_pair;
	      for ( const auto & e : evts )
		if ( e.from_st == tp.first && e.to_st == tp.second )
		  evts_pair.push_back( e );

	      if ( evts_pair.empty() && !par.emit_all_trans ) continue;

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
	      const bool empty_pair = tshape_pair.n_events == 0;
	      write_trans_stats( tshape_pair, "TRANS", par.emit_all_trans && empty_pair );
	      writer.value( "TRANS_N_MATCHED",   n_matched_pair );
	      writer.value( "TRANS_PCT_MATCHED", tshape_pair.n_events > 0
		    ? 100.0 * n_matched_pair / tshape_pair.n_events
		    : ( par.emit_all_trans ? std::numeric_limits<double>::quiet_NaN() : 0.0 ) );
	      if ( profile_pair_dir.valid )
		emit_profile_rows( profile_pair_dir, true );
	      writer.unlevel( "TRANS" );
	    }

	  // BOUT: pair-specific aligned profiles, one per observed bout-to-bout stage pair.
	  std::set<std::pair<int,int>> bout_pairs;
	  for ( const auto & e : bout_events )
	    bout_pairs.insert( std::make_pair( e.from_st , e.to_st ) );
	  if ( par.emit_all_trans )
	    for ( int from = 0; from < K; from++ )
	      for ( int to = 0; to < K; to++ )
		if ( from != to ) bout_pairs.insert( std::make_pair( from, to ) );

	  for ( const auto & tp : bout_pairs )
	    {
	      std::vector<hd_event_t> bout_events_pair;
	      for ( const auto & e : bout_events )
		if ( e.from_st == tp.first && e.to_st == tp.second )
		  bout_events_pair.push_back( e );

	      if ( bout_events_pair.empty() && !par.emit_all_trans ) continue;

	      hd_trans_stats_t tshape_bout_pair;
              hd_profile_t     profile_bout_pair;
	      compute_trans_shape( der, bout_events_pair, context_mask, is_stable, dat, par, dat.Fs, threestate,
                                   false, false, tshape_bout_pair, profile_bout_pair );

	      const std::string bout_label = hd_state_label( threestate , tp.first ) + "to"
		+ hd_state_label( threestate , tp.second );

	      writer.level( bout_label , "BOUT_TRANS" );
	      write_trans_stats( tshape_bout_pair, "BOUT_TRANS",
				 par.emit_all_trans && tshape_bout_pair.n_events == 0 );
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

// Long-range summaries are deliberately global: unlike the main HDSTATS
// summaries, they describe the temporal ordering of one continuous trace and
// therefore cannot be meaningfully restricted to arbitrary annotation masks.
static double hd_sample_sd( const std::vector<double> & x , const int lo , const int hi )
{
  const int n = hi - lo + 1;
  if ( n < 2 ) return std::numeric_limits<double>::quiet_NaN();
  double sum = 0.0;
  for ( int i = lo; i <= hi; i++ ) sum += x[i];
  const double mean = sum / n;
  double ss = 0.0;
  for ( int i = lo; i <= hi; i++ )
    {
      const double d = x[i] - mean;
      ss += d * d;
    }
  return std::sqrt( ss / ( n - 1 ) );
}

static double hd_sample_sd_values( const std::vector<double> & x );

// Prefix moments make a complete window's mean and sample SD O(1).  Long-
// range BETWEEN_SD/WITHIN_SD evaluates every possible non-overlapping tiling
// origin; rescanning every sample of every window would otherwise make that
// exact boundary-invariant calculation O(N * window-width).
struct hd_prefix_moments_t
{
  std::vector<long double> sum;
  std::vector<long double> sumsq;

  explicit hd_prefix_moments_t( const std::vector<double> & x )
  {
    const int n = (int)x.size();
    sum.assign( n + 1, 0.0L );
    sumsq.assign( n + 1, 0.0L );
    for ( int i = 0; i < n; ++i )
      {
        const long double v = x[i];
        sum[i+1] = sum[i] + v;
        sumsq[i+1] = sumsq[i] + v * v;
      }
  }

  void window( const int lo, const int hi, double * mean, double * sd ) const
  {
    const int n = hi - lo + 1;
    const long double s = sum[hi+1] - sum[lo];
    const long double ss = sumsq[hi+1] - sumsq[lo];
    const long double m = s / n;
    long double centered = ss - s * m;
    // Floating cancellation may turn an exactly constant window into a
    // microscopic negative value.  The reference two-pass implementation
    // returns zero in that case.
    if ( centered < 0.0L && centered > -1e-18L * ( std::fabs(ss) + 1.0L ) ) centered = 0.0L;
    if ( centered < 0.0L ) centered = 0.0L;
    *mean = (double)m;
    *sd = n < 2 ? std::numeric_limits<double>::quiet_NaN()
      : std::sqrt( (double)( centered / ( n - 1 ) ) );
  }
};

struct hd_lr_tiling_stats_t
{
  double sum_nwin = 0.0;
  double sum_within_sd = 0.0;
  double sum_between_sd = 0.0;
  int n_offsets_used = 0;
};

// Preserve the old statistic exactly at the definition level: each offset
// retains only wholly valid blocks, requires four blocks, then contributes
// one (unweighted) mean WITHIN_SD and BETWEEN_SD.  Only the repeated scans of
// an already-known complete block are replaced with prefix lookups.
static hd_lr_tiling_stats_t hd_lr_tiling_stats( const hd_data_t & dat,
                                                const std::vector<double> & x,
                                                const int core_lo,
                                                const int core_hi,
                                                const int width,
                                                const bool ps_trace,
                                                const bool gap_aware )
{
  hd_lr_tiling_stats_t out;
  const int n = dat.N;
  std::vector<int> invalid_prefix( n + 1, 0 ), ps_break_prefix( n + 1, 0 );
  for ( int i = 0; i < n; ++i )
    {
      invalid_prefix[i+1] = invalid_prefix[i] + ( dat.valid[i] ? 0 : 1 );
      const bool ps_break = i > 0 && i != core_lo && !dat.valid[i-1];
      ps_break_prefix[i+1] = ps_break_prefix[i] + ( ps_break ? 1 : 0 );
    }
  const hd_prefix_moments_t moments( x );

  for ( int offset = 0; offset < width; ++offset )
    {
      std::vector<double> means;
      double sum_within = 0.0;
      for ( int lo = core_lo + offset; lo + width - 1 <= core_hi; lo += width )
        {
          const int hi = lo + width - 1;
          if ( invalid_prefix[hi+1] != invalid_prefix[lo] ) continue;
          if ( gap_aware && ps_trace && ps_break_prefix[hi+1] != ps_break_prefix[lo] ) continue;
          double mean, sd;
          moments.window( lo, hi, &mean, &sd );
          means.push_back( mean );
          sum_within += sd;
        }
      const int nwin = (int)means.size();
      if ( nwin < 4 ) continue;
      out.sum_nwin += nwin;
      out.sum_within_sd += sum_within / nwin;
      out.sum_between_sd += hd_sample_sd_values( means );
      ++out.n_offsets_used;
    }
  return out;
}

// Long-range missing-data bookkeeping is intentionally separate from the
// ordinary HDSTATS validity summary.  In particular, runs retain their source
// sample indexes: they are never concatenated into an artificial time axis.
struct hd_lr_quality_t
{
  std::vector<std::pair<int,int>> runs;
  int first_valid = -1;
  int last_valid = -1;
  int n_valid = 0;
  int n_invalid = 0;
  int n_internal_gaps = 0;
  int longest_internal_gap = 0;
};

// Timing summaries are only meaningful when the posterior assigns a
// minimally appreciable amount of mass to the quantity being timed.  This is
// deliberately expressed in posterior-mass seconds, so it has the same
// interpretation for raw stage probabilities and conditional series.
static const double hd_lr_timing_min_mass_sec = 120.0;
static const double hd_lr_timing_min_denom = 1e-8;

struct hd_lr_timing_t
{
  double com = std::numeric_limits<double>::quiet_NaN();
  double com_m = std::numeric_limits<double>::quiet_NaN();
  double t50 = std::numeric_limits<double>::quiet_NaN();
  double t50_m = std::numeric_limits<double>::quiet_NaN();
  double mass_m = 0.0;
  int n_valid = 0;
  bool estimable = false;
  std::string na_reason = "NO_ELIGIBLE_OBSERVATIONS";
};

static hd_lr_quality_t hd_lr_quality( const hd_data_t & dat,
                                      const std::vector<bool> * allowed = NULL )
{
  hd_lr_quality_t q;
  int i = 0;
  while ( i < dat.N )
    {
      const bool usable = dat.valid[i] && ( allowed == NULL || (*allowed)[i] );
      if ( ! usable )
        {
          // Outside an analysis interval is intentionally not an invalid
          // observation for its long-range quality accounting.
          if ( allowed == NULL || (*allowed)[i] ) ++q.n_invalid;
          ++i;
          continue;
        }
      const int lo = i;
      while ( i < dat.N && dat.valid[i] && ( allowed == NULL || (*allowed)[i] ) )
        { ++q.n_valid; ++i; }
      q.runs.push_back( std::make_pair( lo, i - 1 ) );
    }
  if ( ! q.runs.empty() )
    {
      q.first_valid = q.runs.front().first;
      q.last_valid = q.runs.back().second;
      for ( int r = 1; r < (int)q.runs.size(); ++r )
        {
          const int gap = q.runs[r].first - q.runs[r-1].second - 1;
          ++q.n_internal_gaps;
          q.longest_internal_gap = std::max( q.longest_internal_gap, gap );
        }
    }
  return q;
}

// Centre of mass and T50 retain source-clock positions across excluded gaps:
// a retained sample after a gap is not silently moved earlier in the night.
// x may contain NaNs to mark conditionally ineligible samples.
static hd_lr_timing_t hd_lr_timing( const std::vector<double> & x,
                                    const hd_data_t & dat,
                                    const hd_lr_quality_t & quality )
{
  hd_lr_timing_t out;
  if ( (int)x.size() != dat.N || quality.first_valid < 0 || quality.last_valid < quality.first_valid )
    return out;

  double mass = 0.0;
  double moment_min = 0.0;
  for ( const auto & run : quality.runs )
    for ( int i = run.first; i <= run.second; ++i )
      {
        const double v = x[i];
        if ( !std::isfinite( v ) || v < 0.0 ) continue;
        ++out.n_valid;
        mass += v;
        // A sample represents an interval, not an instantaneous observation:
        // use its midpoint for the centre of mass.
        moment_min += v * ( i - quality.first_valid + 0.5 ) / dat.Fs / 60.0;
      }

  out.mass_m = mass / dat.Fs / 60.0;
  if ( out.n_valid == 0 ) return out;
  if ( mass / dat.Fs < hd_lr_timing_min_mass_sec )
    {
      out.na_reason = "INSUFFICIENT_POSTERIOR_MASS";
      return out;
    }

  out.com_m = moment_min / mass;
  const double duration_min = ( quality.last_valid - quality.first_valid + 1 ) / dat.Fs / 60.0;
  if ( duration_min <= 0.0 )
    {
      out.na_reason = "INSUFFICIENT_TIME_VARIATION";
      return out;
    }
  out.com = out.com_m / duration_min;

  const double half_mass = 0.5 * mass;
  double cumulative = 0.0;
  for ( const auto & run : quality.runs )
    for ( int i = run.first; i <= run.second; ++i )
      {
        const double v = x[i];
        if ( !std::isfinite( v ) || v < 0.0 ) continue;
        if ( cumulative + v >= half_mass )
          {
            // Treat posterior mass as uniform inside the sample interval,
            // which gives a continuous, sampling-rate-independent T50.
            const double frac = v > 0.0 ? ( half_mass - cumulative ) / v : 0.0;
            out.t50_m = ( i - quality.first_valid + frac ) / dat.Fs / 60.0;
            out.t50 = out.t50_m / duration_min;
            out.estimable = true;
            out.na_reason = "OK";
            return out;
          }
        cumulative += v;
      }

  // This should be unreachable after a finite positive mass calculation, but
  // avoids reporting a partially populated result if a future input type
  // changes the numerical assumptions above.
  out.na_reason = "NO_CUMULATIVE_HALF_MASS";
  return out;
}

static void hd_write_lr_timing( const hd_lr_timing_t & timing,
                                const bool emit_output,
                                const bool emit_qc )
{
  if ( emit_output )
    {
      writer.value( "COM", timing.com );
      writer.value( "COM_M", timing.com_m );
      writer.value( "T50", timing.t50 );
      writer.value( "T50_M", timing.t50_m );
    }
  if ( emit_qc )
    {
      writer.level( 1, "QC" );
      writer.value( "TIMING_MASS_M", timing.mass_m );
      writer.value( "TIMING_NVALID", timing.n_valid );
      writer.value( "TIMING_ESTIMABLE", timing.estimable ? 1 : 0 );
      writer.value( "TIMING_NA_REASON", timing.na_reason );
      writer.unlevel( "QC" );
    }
}

static double hd_sample_sd_values( const std::vector<double> & x )
{
  return x.size() < 2 ? std::numeric_limits<double>::quiet_NaN()
    : hd_sample_sd( x, 0, (int)x.size() - 1 );
}

static void hd_write_lr_status_values( const bool estimable, const std::string & reason )
{
  writer.value( "ESTIMABLE", estimable ? 1 : 0 );
  if ( ! estimable ) writer.value( "NA_REASON", reason );
}

static void hd_write_lr_status( const bool estimable, const std::string & reason,
                                const bool emit_qc )
{
  if ( !emit_qc ) return;
  writer.level( 1, "QC" );
  hd_write_lr_status_values( estimable, reason );
  writer.unlevel( "QC" );
}

// A companion, metric-granular row removes any ambiguity when a shared
// result row contains (for example) an estimable PS but a zero-variance SIM.
// Its inherited strata make it directly joinable to the numerical row.
static void hd_write_lr_metric_quality( const std::string & metric,
                                        const std::string & count_name,
                                        const int count,
                                        const bool estimable,
                                        const std::string & reason,
                                        const bool emit_qc )
{
  if ( !emit_qc ) return;
  writer.level( metric, "METRIC" );
  writer.level( 1, "QC" );
  writer.value( count_name, count );
  hd_write_lr_status_values( estimable, reason );
  writer.unlevel( "QC" );
  writer.unlevel( "METRIC" );
}

static double hd_pair_corr( const std::vector<double> & x , const int lo , const int hi , const int lag )
{
  const int n = hi - lo + 1 - lag;
  if ( n < 2 ) return std::numeric_limits<double>::quiet_NaN();

  double sx = 0.0, sy = 0.0;
  for ( int i = lo; i <= hi - lag; i++ )
    {
      sx += x[i];
      sy += x[i+lag];
    }
  const double mx = sx / n;
  const double my = sy / n;
  double sxx = 0.0, syy = 0.0, sxy = 0.0;
  for ( int i = lo; i <= hi - lag; i++ )
    {
      const double dx = x[i] - mx;
      const double dy = x[i+lag] - my;
      sxx += dx * dx;
      syy += dy * dy;
      sxy += dx * dy;
    }
  return sxx > 0.0 && syy > 0.0 ? sxy / std::sqrt( sxx * syy )
    : std::numeric_limits<double>::quiet_NaN();
}

static void write_hdstats_longrange(
    const hd_data_t & dat,
    const hd_derived_t & der,
    const hd_params_t & par,
    const hd_lr_quality_t & quality,
    hdstats_pc::feature_map_t * pc_features = NULL,
    const bool emit_output = true )
{
  const int core_lo = quality.first_valid;
  const int core_hi = quality.last_valid;
  const int n_core = quality.n_valid;
  const double NaN = std::numeric_limits<double>::quiet_NaN();

  auto emit_representation = [&]( const bool threestate ) {
    const std::vector<std::vector<double>> & P = threestate ? dat.p3 : dat.p;
    const std::vector<double> & H  = threestate ? der.H3  : der.H;
    const std::vector<double> & PS = threestate ? der.TV3 : der.TV;
    const int K = threestate ? 3 : 5;

    writer.level( threestate ? "3" : "5", "NSS" );

    std::vector<hd_lr_timing_t> stage_timing( K );
    for ( int k = 0; k < K; ++k )
      stage_timing[k] = hd_lr_timing( P[k], dat, quality );

    // Lagged stage persistence, full-vector displacement, and a centered
    // multivariate similarity.  Centering makes SIM a temporal similarity,
    // rather than merely a restatement of overall stage occupancy.
    for ( double requested_lag_sec : par.lr_lags_sec )
      {
        const int lag = std::max( 1, (int)std::round( requested_lag_sec * dat.Fs ) );
        const double lag_sec = lag / dat.Fs;

        std::vector<double> mean( K, 0.0 );
        for ( int k = 0; k < K; k++ )
          for ( const auto & run : quality.runs )
            for ( int i = run.first; i <= run.second; ++i ) mean[k] += P[k][i];
        if ( n_core > 0 )
          for ( int k = 0; k < K; k++ ) mean[k] /= n_core;

        double mean_l1 = 0.0;
        double dot = 0.0, ss0 = 0.0, ss1 = 0.0;
        std::vector<std::vector<double>> acf0( K ), acf1( K );
        int npair = 0;
        // Strict mode has always required a core of at least 3 * lag samples;
        // its shortest accepted trace therefore has 2 * lag retained pairs.
        // Preserve that documented minimum in gap-aware mode after excluding
        // every pair that would cross a valid-run boundary.
        const int min_lag_pairs = 2 * lag;
        for ( const auto & run : quality.runs )
          for ( int i = run.first; i + lag <= run.second; ++i )
            {
              ++npair;
              double l1 = 0.0;
              for ( int k = 0; k < K; k++ )
                {
                  l1 += std::fabs( P[k][i+lag] - P[k][i] );
                  const double x0 = P[k][i] - mean[k];
                  const double x1 = P[k][i+lag] - mean[k];
                  dot += x0 * x1;
                  ss0 += x0 * x0;
                  ss1 += x1 * x1;
                  acf0[k].push_back( P[k][i] );
                  acf1[k].push_back( P[k][i+lag] );
                }
              mean_l1 += 0.5 * l1;
            }

        writer.level( lag_sec, "LRLAG" );
        if ( par.qc )
          {
            writer.level( 1, "QC" );
            writer.value( "NPAIR", npair );
            writer.unlevel( "QC" );
          }
        const bool ps_ok = npair >= min_lag_pairs;
        const bool sim_ok = ps_ok && ss0 > 0.0 && ss1 > 0.0;
        if ( emit_output )
          {
            writer.value( "PS", ps_ok ? mean_l1 / npair : NaN );
            writer.value( "SIM", sim_ok ? dot / std::sqrt( ss0 * ss1 ) : NaN );
          }
        if ( pc_features != NULL )
          {
            const std::string s = "NSS/" + std::string( threestate ? "3" : "5" )
              + ";LRLAG/" + Helper::dbl2str( lag_sec );
            hd_pc_add( pc_features, hdstats_pc::key( "PS", s ), ps_ok ? mean_l1 / npair : NaN );
            hd_pc_add( pc_features, hdstats_pc::key( "SIM", s ), sim_ok ? dot / std::sqrt( ss0 * ss1 ) : NaN );
          }
        hd_write_lr_status( ps_ok && sim_ok, npair < min_lag_pairs ? "INSUFFICIENT_PAIRS" :
                            "ZERO_VARIANCE", par.qc );
        hd_write_lr_metric_quality( "PS", "NPAIR", npair, ps_ok,
                                    "INSUFFICIENT_PAIRS", par.qc );
        hd_write_lr_metric_quality( "SIM", "NPAIR", npair, sim_ok,
                                    npair < min_lag_pairs ? "INSUFFICIENT_PAIRS" : "ZERO_VARIANCE", par.qc );
        for ( int k = 0; k < K; k++ )
          {
            writer.level( hd_state_label( threestate, k ), "SS" );
            const int n = (int)acf0[k].size();
            double sx = 0.0, sy = 0.0, sxx = 0.0, syy = 0.0, sxy = 0.0;
            for ( int j = 0; j < n; ++j ) { sx += acf0[k][j]; sy += acf1[k][j]; }
            const double mx = n ? sx / n : 0.0, my = n ? sy / n : 0.0;
            for ( int j = 0; j < n; ++j )
              {
                const double dx = acf0[k][j] - mx, dy = acf1[k][j] - my;
                sxx += dx * dx; syy += dy * dy; sxy += dx * dy;
              }
            const bool acf_ok = n >= min_lag_pairs && sxx > 0.0 && syy > 0.0;
            if ( par.qc )
              {
                writer.level( 1, "QC" );
                writer.value( "NPAIR", n );
                writer.unlevel( "QC" );
              }
            if ( emit_output ) writer.value( "ACF", acf_ok ? sxy / std::sqrt( sxx * syy ) : NaN );
            if ( pc_features != NULL )
              hd_pc_add( pc_features, hdstats_pc::key( "ACF",
                           "SS/" + hd_state_label(threestate,k) + ";NSS/"
                           + ( threestate ? "3" : "5" ) + ";LRLAG/" + Helper::dbl2str(lag_sec) ),
                         acf_ok ? sxy / std::sqrt( sxx * syy ) : NaN );
            hd_write_lr_status( acf_ok, n < min_lag_pairs ? "INSUFFICIENT_PAIRS" : "ZERO_VARIANCE", par.qc );
            hd_write_lr_metric_quality( "ACF", "NPAIR", n, acf_ok,
                                        n < min_lag_pairs ? "INSUFFICIENT_PAIRS" : "ZERO_VARIANCE", par.qc );
            writer.unlevel( "SS" );
          }
        writer.unlevel( "LRLAG" );
      }

    // Slopes are per full-recording normalized time unit (from the start to
    // the end of the valid core).  They are intentionally not annualized or
    // epoch-count dependent.
    std::vector<const std::vector<double> *> trend_x;
    std::vector<std::string> trend_label;
    for ( int k = 0; k < K; k++ )
      {
        trend_x.push_back( &P[k] );
        trend_label.push_back( hd_state_label( threestate, k ) );
      }
    trend_x.push_back( &H );  trend_label.push_back( "H" );
    trend_x.push_back( &PS ); trend_label.push_back( "PS" );

    for ( int q = 0; q < (int)trend_x.size(); q++ )
      {
        const std::vector<double> & x = *trend_x[q];
        std::vector<int> idx;
        for ( const auto & run : quality.runs )
          for ( int i = run.first; i <= run.second; ++i )
            // PS at the first sample of a later run would be a step from an
            // invalid sample in the ordinary per-sample trace, so it is not a
            // valid long-range observation.
            if ( par.lr_policy == "strict" || q < K + 1 || i == core_lo || dat.valid[i-1] ) idx.push_back( i );
        double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
        for ( int i : idx )
          {
            const double u = core_hi > core_lo ? (double)( i - core_lo ) / ( core_hi - core_lo ) : 0.0;
            sx += u;
            sy += x[i];
            sxx += u * u;
            sxy += u * x[i];
          }
        const int n = (int)idx.size();
        const double den = n * sxx - sx * sx;
        const bool slope_ok = n >= 2 && den != 0.0;
        writer.level( trend_label[q], "VAR" );
        const double slope = slope_ok ? ( n * sxy - sx * sy ) / den : NaN;
        if ( emit_output ) writer.value( "SLOPE", slope );
        if ( q < K ) hd_write_lr_timing( stage_timing[q], emit_output, par.qc );
        if ( pc_features != NULL )
          {
            const std::string s = "NSS/" + std::string( threestate ? "3" : "5" )
              + ";VAR/" + trend_label[q];
            hd_pc_add( pc_features, hdstats_pc::key( "SLOPE", s ), slope );
            if ( q < K && !( threestate && trend_label[q] == "NR" ) )
              {
                hd_pc_add( pc_features, hdstats_pc::key( "T50", s ), stage_timing[q].t50 );
                hd_pc_add( pc_features, hdstats_pc::key( "COM", s ), stage_timing[q].com );
              }
          }
        hd_write_lr_status( slope_ok, "INSUFFICIENT_TIME_VARIATION", par.qc );
        hd_write_lr_metric_quality( "SLOPE", "NVALID", n, slope_ok,
                                    "INSUFFICIENT_TIME_VARIATION", par.qc );
        writer.unlevel( "VAR" );
      }

    // N3 within NREM and REM within sleep separate the timing of those stages
    // from changing wake/NREM opportunity.  They are inherently five-state
    // quantities and are therefore not emitted for the collapsed 3-state
    // representation.  Samples with essentially zero conditioning mass are
    // ineligible rather than assigned a numerically unstable ratio.
    if ( ! threestate )
      {
        std::vector<double> n3_nrem( dat.N, NaN );
        std::vector<double> r_sleep( dat.N, NaN );
        for ( const auto & run : quality.runs )
          for ( int i = run.first; i <= run.second; ++i )
            {
              const double nrem = P[1][i] + P[2][i] + P[3][i];
              const double sleep = 1.0 - P[0][i];
              if ( nrem > hd_lr_timing_min_denom ) n3_nrem[i] = P[3][i] / nrem;
              if ( sleep > hd_lr_timing_min_denom ) r_sleep[i] = P[4][i] / sleep;
            }

        const auto emit_conditional_timing = [&]( const std::string & label,
                                                   const std::vector<double> & x ) {
          writer.level( label, "VAR" );
          const hd_lr_timing_t timing = hd_lr_timing( x, dat, quality );
          hd_write_lr_timing( timing, emit_output, par.qc );
          if ( pc_features != NULL )
            {
              hd_pc_add( pc_features, hdstats_pc::key( "T50", "NSS/5;VAR/" + label ), timing.t50 );
              hd_pc_add( pc_features, hdstats_pc::key( "COM", "NSS/5;VAR/" + label ), timing.com );
            }
          writer.unlevel( "VAR" );
        };
        emit_conditional_timing( "N3_NREM", n3_nrem );
        emit_conditional_timing( "R_SLEEP", r_sleep );
      }

    // WITHIN_SD/BETWEEN_SD tile the valid core into non-overlapping blocks of
    // a fixed width. A single tiling is anchored arbitrarily at the first
    // valid sample, so its block boundaries - and therefore both statistics -
    // depend on exactly where that happens to fall, for no physiological
    // reason. To remove that dependence, every one of the `width` possible
    // tiling phases is evaluated (each still built from the same handful of
    // genuinely non-overlapping blocks as a single tiling would use) and the
    // reported WITHIN_SD/BETWEEN_SD are the mean across phases: an unbiased
    // reduction of phase noise, not a different quantity (unlike naively
    // sliding a single window by one sample and pooling all the resulting,
    // heavily autocorrelated means into one SD, which *would* measure
    // something else - see the derivation notes alongside this command).
    for ( double requested_window_sec : par.lr_windows_sec )
      {
        const int width = std::max( 1, (int)std::round( requested_window_sec * dat.Fs ) );
        const double window_sec = width / dat.Fs;
        for ( int q = 0; q < (int)trend_x.size(); q++ )
          {
            const std::vector<double> & x = *trend_x[q];

            const hd_lr_tiling_stats_t tiling = hd_lr_tiling_stats(
              dat, x, core_lo, core_hi, width, q == K + 1,
              par.lr_policy == "gap-aware" );
            const double sum_nwin = tiling.sum_nwin;
            const double sum_within_sd = tiling.sum_within_sd;
            const double sum_between_sd = tiling.sum_between_sd;
            const int n_offsets_used = tiling.n_offsets_used;

            const bool within_ok  = n_offsets_used > 0 && width >= 2;
            const bool between_ok = n_offsets_used > 0;
            // NWIN is now the mean block count per phase (same scale/meaning
            // as the old single-tiling count); N_OFFSETS records how many of
            // the `width` candidate phases actually contributed.
            const int nwin_mean = n_offsets_used > 0 ? (int)std::round( sum_nwin / n_offsets_used ) : 0;
            writer.level( window_sec, "LRWIN" );
            writer.level( trend_label[q], "VAR" );
            if ( par.qc )
              {
                writer.level( 1, "QC" );
                writer.value( "NWIN", nwin_mean );
                writer.value( "N_OFFSETS", n_offsets_used );
                writer.unlevel( "QC" );
              }
            const double within_sd  = within_ok  ? sum_within_sd  / n_offsets_used : NaN;
            const double between_sd = between_ok ? sum_between_sd / n_offsets_used : NaN;
            if ( emit_output )
              {
                writer.value( "WITHIN_SD", within_sd );
                writer.value( "BETWEEN_SD", between_sd );
              }
            if ( pc_features != NULL )
              {
                const std::string s = "NSS/" + std::string( threestate ? "3" : "5" )
                  + ";VAR/" + trend_label[q] + ";LRWIN/" + Helper::dbl2str(window_sec);
                hd_pc_add( pc_features, hdstats_pc::key( "WITHIN_SD", s ), within_sd );
                hd_pc_add( pc_features, hdstats_pc::key( "BETWEEN_SD", s ), between_sd );
              }
            hd_write_lr_status( within_ok && between_ok,
                                n_offsets_used == 0 ? "NO_ELIGIBLE_WINDOWS" : "INSUFFICIENT_WINDOWS", par.qc );
            hd_write_lr_metric_quality( "WITHIN_SD", "NWIN", nwin_mean, within_ok,
                                        n_offsets_used == 0 ? "NO_ELIGIBLE_WINDOWS" : "INSUFFICIENT_WINDOW_SAMPLES", par.qc );
            hd_write_lr_metric_quality( "BETWEEN_SD", "NWIN", nwin_mean, between_ok,
                                        n_offsets_used == 0 ? "NO_ELIGIBLE_WINDOWS" : "INSUFFICIENT_WINDOWS", par.qc );
            writer.unlevel( "VAR" );
            writer.unlevel( "LRWIN" );
          }
      }

    writer.unlevel( "NSS" );
  };

  emit_representation( false );
  if ( par.do_3state ) emit_representation( true );
}

static void set_hdstats_offset_outputs( bool enabled )
{
  // Pooled (ALL) and pair-specific profiles share a table. Per-state
  // posteriors are emitted only for named directional-pair levels.
  const std::vector<std::string> trans_sec_tables = {
    "NSS,TRANS,SEC",
    "ANNOT,NSS,TRANS,SEC"
  };
  const std::vector<std::string> bout_sec_tables = {
    "NSS,BOUT_TRANS,SEC",
    "ANNOT,NSS,BOUT_TRANS,SEC"
  };
  const std::vector<std::string> sec_vars = {
    "N", "H", "CONF", "MG", "PS", "P_W", "P_NR", "P_N1", "P_N2", "P_N3", "P_R"
  };

  for ( int i = 0; i < (int)trans_sec_tables.size(); i++ )
    for ( int j = 0; j < (int)sec_vars.size(); j++ )
      globals::cmddefs().register_var( "HDSTATS", trans_sec_tables[i], sec_vars[j], enabled );
  for ( int i = 0; i < (int)bout_sec_tables.size(); i++ )
    for ( int j = 0; j < (int)sec_vars.size(); j++ )
      globals::cmddefs().register_var( "HDSTATS", bout_sec_tables[i], sec_vars[j], enabled );
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
  const bool emit_ultradian_only = param.has( "emit-ultradian" ) && param.yesno( "emit-ultradian" );

  if ( emit_nr_pp )
    hd_make_pp_nrem( edf , param );

  // weights= can also be used with emit-stats=F to make an in-memory
  // posterior correction before a later command (typically WRITE).
  if ( ! emit_stats && ! param.has( "weights" ) && ! emit_ultradian_only ) return;

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
  par.qc              = param.has( "qc" ) ? param.yesno( "qc" ) : false;
  par.min_events_profile = param.has( "min-events" ) ? param.requires_int( "min-events" ) : 3;
  par.min_samples_region = param.has( "min-samples" ) ? param.requires_int( "min-samples" ) : 20;
  par.emit_annot_typed_classes = param.has( "emit-annots-by-type" ) ? param.yesno( "emit-annots-by-type" ) : false;
  par.emit_all_trans = param.yesno( "emit-all-trans" );
  par.emit_annot_label   =
    param.has( "emit-annots" ) ? ( ( ! param.empty("emit-annots" ) ) ? param.value( "emit-annots" ) : "hd" )
    : ( par.emit_annot_typed_classes ? "hd" : "" );
  // Reuse an explicitly selected transition root; otherwise HDSTATS' default
  // generated-annotation root is "hd". emit-spt=F/N/0 disables this trio.
  par.emit_spt_label = par.emit_annot_label.empty() ? "hd" : par.emit_annot_label;
  if ( param.has( "emit-spt" ) )
    {
      const std::string setting = param.value( "emit-spt", true );
      if ( setting == "F" || setting == "N" || setting == "0" )
        par.emit_spt_label = "";
      else if ( setting == "T" || setting == "Y" || setting == "1" || param.empty( "emit-spt" ) )
        { /* retain the inherited/default root */ }
      else
        par.emit_spt_label = param.value( "emit-spt" );
    }
  // Unlike emit-spt, cycle annotation emission is opt-in. T selects the
  // conventional "hd" root; any other non-boolean value is a caller-supplied
  // root (e.g. emit-cycles=foo creates foo_cycle).
  par.emit_cycle_label = "";
  if ( param.has( "emit-cycles" ) )
    {
      const std::string setting = param.value( "emit-cycles", true );
      if ( setting == "T" || setting == "Y" || setting == "1" || param.empty( "emit-cycles" ) )
        par.emit_cycle_label = "hd";
      else if ( setting != "F" && setting != "N" && setting != "0" )
        par.emit_cycle_label = param.value( "emit-cycles" );
    }
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
  par.do_longrange = param.has( "longrange" ) ? param.yesno( "longrange" ) : true;
  par.lr_lags_sec = param.has( "lr-lags" )
    ? param.dblvector( "lr-lags" )
    : std::vector<double>{ 600.0, 1800.0, 3600.0, 5400.0 };
  par.lr_windows_sec = param.has( "lr-windows" )
    ? param.dblvector( "lr-windows" )
    : std::vector<double>{ 1800.0, 3600.0 };
  par.lr_min_duration_sec = param.has( "lr-min-duration" )
    ? param.requires_dbl( "lr-min-duration" ) : 4.5 * 3600.0;
  par.lr_policy = ( ! param.has( "gap-aware" ) || param.yesno( "gap-aware" ) ) ? "gap-aware" : "strict";
  par.lr_fail_on_empty = param.has( "fail-on-empty" ) ? param.yesno( "fail-on-empty" ) : false;
  par.fixed_duration_sec = 0.0;
  if ( param.has( "fixed" ) )
    {
      const double fixed_hours = param.requires_dbl( "fixed" );
      if ( ! std::isfinite( fixed_hours ) || fixed_hours <= 0.0 )
        Helper::halt( "HDSTATS fixed= must be a finite positive duration in hours" );
      par.analysis = "fixed";
      par.fixed_duration_sec = fixed_hours * 3600.0;
    }
  else
    par.analysis = ( ! param.has( "spt" ) || param.yesno( "spt" ) ) ? "spt" : "recording";
  const double analysis_flank_sec = param.has( "spt-flank" ) ? param.requires_dbl( "spt-flank" ) : 300.0;
  par.analysis_flank_pre_sec = param.has( "spt-flank-pre" ) ? param.requires_dbl( "spt-flank-pre" ) : analysis_flank_sec;
  par.analysis_flank_post_sec = param.has( "spt-flank-post" ) ? param.requires_dbl( "spt-flank-post" ) : analysis_flank_sec;
  par.pc_out = param.has( "pc-out" ) ? param.value( "pc-out" ) : "";
  par.pc_proj = param.has( "pc-proj" ) ? Helper::expand( param.value( "pc-proj" ) ) : "";
  const bool emit_model_desc = param.yesno( "model-desc" );
  if ( emit_model_desc && par.pc_proj.empty() )
    Helper::halt( "HDSTATS model-desc requires pc-proj=<HDSTATS-PC V3 model>" );
  par.pc_min_loading = param.has( "pc-min-loading" ) ? param.requires_dbl( "pc-min-loading" ) : 0.10;
  if ( ! std::isfinite(par.pc_min_loading) || par.pc_min_loading < 0.0 )
    Helper::halt( "pc-min-loading must be a finite non-negative value" );
  const bool emit_main_metrics = par.pc_proj.empty() ||
    ( param.has( "pc-all" ) && param.yesno( "pc-all" ) );
  par.do_ultradian = param.has( "ultradian" ) ? param.yesno( "ultradian" ) : true;
  par.ultra_fwhm_sec = param.has( "ultra-fwhm" ) ? param.requires_dbl( "ultra-fwhm" ) : 600.0;
  par.ultra_prominence = param.has( "ultra-prominence" ) ? param.requires_dbl( "ultra-prominence" ) : 1.0;
  par.ultra_min_sleep_support = param.has( "ultra-min-sleep-support" )
    ? param.requires_dbl( "ultra-min-sleep-support" ) : 0.25;
  par.ultra_weights = {{ 0.0, -1.0, 1.0, 2.0, -2.0 }};
  if ( param.has( "ultra-weights" ) )
    {
      const std::vector<double> w = param.dblvector( "ultra-weights" );
      if ( w.size() != 5 )
        Helper::halt( "HDSTATS ultra-weights= requires five finite values in W,N1,N2,N3,R order" );
      for ( int k = 0; k < 5; ++k )
        {
          if ( ! std::isfinite( w[k] ) )
            Helper::halt( "HDSTATS ultra-weights= requires five finite values in W,N1,N2,N3,R order" );
          par.ultra_weights[k] = w[k];
        }
    }
  par.emit_ultradian = param.has( "emit-ultradian" ) && param.yesno( "emit-ultradian" );
  par.ultra_diagnostics = param.has( "ultra-diagnostics" ) && param.yesno( "ultra-diagnostics" );
  par.do_ultra_condition = param.has( "ultra-condition" ) ? param.yesno( "ultra-condition" ) : true;
  par.ultra_condition_weighted = param.has( "ultra-condition-weighted" )
    && param.yesno( "ultra-condition-weighted" );
  par.ultra_max_condition = param.has( "ultra-condition-max" )
    ? param.requires_dbl( "ultra-condition-max" ) : 12.0;
  par.do_ultra_regions = param.has( "ultra-regions" ) ? param.yesno( "ultra-regions" ) : true;

  const bool do_weight_correction = param.has( "weights" );
  const bool update_pp = param.yesno( "update-pp", false );
  std::array<double,5> weights;
  if ( do_weight_correction )
    {
      const std::vector<double> w = param.dblvector( "weights" );
      if ( w.size() != 5 )
        Helper::halt( "HDSTATS weights= requires five finite positive values in W,R,N1,N2,N3 order" );
      for ( int k = 0; k < 5; ++k )
        {
          if ( ! std::isfinite( w[k] ) || w[k] <= 0.0 )
            Helper::halt( "HDSTATS weights= requires five finite positive values in W,R,N1,N2,N3 order" );
          weights[k] = w[k];
        }
    }
  else if ( update_pp )
    Helper::halt( "HDSTATS update-pp=T requires weights=" );

  if ( par.analysis_flank_pre_sec < 0.0 || par.analysis_flank_post_sec < 0.0 )
    Helper::halt( "analysis flanks must be non-negative seconds" );

  if ( par.do_ultradian )
    {
      if ( ! std::isfinite( par.ultra_fwhm_sec ) || par.ultra_fwhm_sec <= 0.0 )
        Helper::halt( "ultra-fwhm must be positive seconds" );
      if ( ! std::isfinite( par.ultra_prominence ) || par.ultra_prominence <= 0.0 )
        Helper::halt( "ultra-prominence must be positive" );
      if ( ! std::isfinite( par.ultra_min_sleep_support )
           || par.ultra_min_sleep_support <= 0.0 || par.ultra_min_sleep_support > 1.0 )
        Helper::halt( "ultra-min-sleep-support must be in (0,1]" );
      if ( !std::isfinite( par.ultra_max_condition ) || par.ultra_max_condition <= 1.0 )
        Helper::halt( "ultra-condition-max must be greater than 1" );
    }

  if ( par.do_longrange )
    {
      if ( ! par.annot_name.empty() )
        Helper::halt( "HDSTATS longrange requires the global trace; use longrange=F with annot=" );
      if ( par.lr_min_duration_sec <= 0.0 )
        Helper::halt( "lr-min-duration must be positive" );
      for ( double x : par.lr_lags_sec )
        if ( x <= 0.0 ) Helper::halt( "all lr-lags must be positive seconds" );
      for ( double x : par.lr_windows_sec )
        if ( x <= 0.0 ) Helper::halt( "all lr-windows must be positive seconds" );
    }

  auto hd_join_secs = []( const std::vector<double> & v ) {
    std::string s;
    for ( int i = 0; i < (int)v.size(); ++i )
      s += ( i ? "," : "" ) + Helper::dbl2str( v[i] );
    return s;
  };

  logger << "  HDSTATS: analysis=" << par.analysis
         << ", 3state=" << ( par.do_3state ? "T" : "F" )
         << ", window=" << par.window_sec << "s"
         << ", lag=" << par.lag_sec << "s"
         << ", stable-min=" << par.stable_min_sec << "s\n";

  logger << "  ULTRADIAN: enabled=" << ( par.do_ultradian ? "T" : "F" );
  if ( par.do_ultradian )
    logger << ", fwhm=" << par.ultra_fwhm_sec << "s"
           << ", prominence=" << par.ultra_prominence
           << ", sleep-support>=" << par.ultra_min_sleep_support
           << "\n             edf=" << ( par.emit_ultradian ? "T" : "F" )
           << ", diagnostics=" << ( par.ultra_diagnostics ? "T" : "F" )
           << ", uvar=" << ( par.do_ultra_condition ? "T" : "F" )
           << ( par.ultra_condition_weighted ? "/weighted" : "/unweighted" )
           << " (max-cond=" << par.ultra_max_condition << ")"
           << ", uq=" << ( par.do_ultra_regions ? "T" : "F" );
  logger << "\n";

  logger << "  LONG-RANGE: enabled=" << ( par.do_longrange ? "T" : "F" );
  if ( par.do_longrange )
    logger << ", lags=" << hd_join_secs( par.lr_lags_sec ) << "s"
           << ", windows=" << hd_join_secs( par.lr_windows_sec ) << "s"
           << "\n              policy=" << par.lr_policy
           << ", min-duration=" << par.lr_min_duration_sec / 3600.0 << "h"
           << ", block-SD=all tiling offsets";
  logger << "\n";

  logger << "  POSTERIORS: W=" << par.ch[0] << ", N1=" << par.ch[1]
         << ", N2=" << par.ch[2] << ", N3=" << par.ch[3]
         << ", R=" << par.ch[4] << "\n";

  // Layer 1: load data
  hd_data_t dat;
  if ( ! dat.load( edf, par ) ) return;

  if ( do_weight_correction )
    {
      const int n_corrected = hd_apply_weight_correction( &dat, par, weights );
      const int n_invalid = dat.N - n_corrected;
      logger << "  HDSTATS weights: applied inverse class-weight correction "
             << "(W,R,N1,N2,N3="
             << weights[0] << "," << weights[1] << "," << weights[2] << ","
             << weights[3] << "," << weights[4] << ") to "
             << n_corrected << " valid samples; left " << n_invalid
             << " invalid sample(s) unchanged\n";
      if ( update_pp )
        {
          hd_update_posterior_channels( edf, dat, par );
          logger << "  HDSTATS update-pp: updated posterior channels in place ["
                 << par.ch[0] << ", " << par.ch[1] << ", " << par.ch[2] << ", "
                 << par.ch[3] << ", " << par.ch[4]
                 << "]; invalid samples and " << par.ch_valid << " were preserved\n";
        }
    }

  // A compact signal-only invocation is useful for EDF QC and avoids
  // generating the normal HDSTATS tables.  It still uses the same validated
  // posterior rows (and optional correction) as a regular analysis.
  if ( ! emit_stats )
    {
      if ( par.do_ultradian && par.emit_ultradian )
        {
          hd_ultradian_t ultra;
          ultra.compute( dat, par );
          hd_emit_ultradian_signals( edf, dat, ultra, par );
          logger << "  ultradian: " << ultra.n_half_cycles << " qualifying half-cycle(s), "
                 << ultra.n_complete_cycles << " strict complete cycle(s)\n";
        }
      return;
    }

  logger << "  N=" << dat.N << " samples at Fs=" << dat.Fs << " Hz ("
	 << dat.N / dat.Fs / 3600.0 << " hours)\n";

  // Layer 2: derive per-sample signals
  hd_derived_t der;
  der.compute( dat, par );

  hd_ultradian_t ultra;
  if ( par.do_ultradian )
    {
      ultra.compute( dat, par );
      if ( par.emit_ultradian ) hd_emit_ultradian_signals( edf, dat, ultra, par );
      logger << "  ultradian: " << ultra.n_half_cycles << " qualifying half-cycle(s), "
             << ultra.n_complete_cycles << " strict complete cycle(s)\n";
    }

  if ( ! par.emit_sig_label.empty() )
    hd_emit_signals( edf, dat, der, par );

  // The sleep-period selector relies on the existing robust hypnodensity
  // onset/final-sleep logic even when the optional HD timing table is off.
  hd_hypno_metrics_t hm;
  hm.compute( dat, par );
  const hd_analysis_t analysis = hd_make_analysis( dat, hm, par );
  const bool make_pc_features = ! par.pc_out.empty() || ! par.pc_proj.empty();
  if ( par.analysis == "fixed" )
    {
      writer.value( "HD_FIXED_AVAILABLE_HOURS", analysis.fixed_available_sec / 3600.0 );
      writer.value( "HD_FIXED_OK", analysis.fixed_ok ? 1 : 0 );
      if ( ! analysis.fixed_ok )
        {
          if ( hm.valid_onset )
            logger << "  warning: HDSTATS fixed="
                   << Helper::dbl2str( par.fixed_duration_sec / 3600.0 )
                   << " hours skips this recording: only "
                   << Helper::dbl2str( analysis.fixed_available_sec / 3600.0 )
                   << " hours are available after sustained sleep onset\n";
          else
            logger << "  warning: HDSTATS fixed= skips this recording: could not define sustained sleep onset\n";
          if ( par.qc )
            {
              writer.level( 1, "QC" );
              const int n_valid = (int)std::count( dat.valid.begin(), dat.valid.end(), true );
              writer.value( "N_VALID", n_valid );
              writer.value( "PCT_VALID", dat.N > 0 ? 100.0 * n_valid / dat.N : 0.0 );
              if ( dat.n_invalid_lead  > 0 ) writer.value( "N_INVALID_LEAD",  dat.n_invalid_lead );
              if ( dat.n_invalid_trail > 0 ) writer.value( "N_INVALID_TRAIL", dat.n_invalid_trail );
              if ( dat.n_invalid_mid   > 0 ) writer.value( "N_INVALID_MID",   dat.n_invalid_mid );
              writer.value( "HD_ANALYSIS", "fixed" );
              writer.value( "HD_ANALYSIS_START_SEC", hm.valid_onset ? analysis.lo / dat.Fs : std::numeric_limits<double>::quiet_NaN() );
              writer.value( "HD_ANALYSIS_STOP_SEC", hm.valid_onset ? ( analysis.hi + 1 ) / dat.Fs : std::numeric_limits<double>::quiet_NaN() );
              writer.value( "HD_ANALYSIS_SEC", hm.valid_onset ? ( analysis.hi - analysis.lo + 1 ) / dat.Fs : 0.0 );
              writer.value( "HD_ANALYSIS_INVALID_SEC", analysis.n_invalid / dat.Fs );
              writer.value( "HD_FIXED_SKIPPED", 1 );
              writer.unlevel( "QC" );
            }
          return;
        }
    }
  if ( analysis.mode != "recording" )
    logger << "  HDSTATS analysis interval: " << ( analysis.hi - analysis.lo + 1 ) / dat.Fs
           << " seconds (samples " << analysis.lo << "-" << analysis.hi << ")\n";

  if ( ! par.emit_spt_label.empty() && ! dat.tp.empty() )
    {
      const uint64_t sample_tp = (uint64_t)std::round( (double)globals::tp_1sec / dat.Fs );
      if ( hm.valid_window )
        {
          const uint64_t onset_tp = dat.tp[hm.onset_idx];
          const uint64_t final_wake_tp = dat.tp[hm.offset_idx] + sample_tp;
          edf.annotations->add( par.emit_spt_label + "_sleep_onset" )->add(
            ".", interval_t( onset_tp, onset_tp ), "." );
          edf.annotations->add( par.emit_spt_label + "_final_wake" )->add(
            ".", interval_t( final_wake_tp, final_wake_tp ), "." );
          edf.annotations->add( par.emit_spt_label + "_spt" )->add(
            ".", interval_t( onset_tp, final_wake_tp ), "." );
          logger << "  emitted soft SPT annotations under '" << par.emit_spt_label << "'\n";
        }
      if ( analysis.mode == "fixed" )
        {
          edf.annotations->add( par.emit_spt_label + "_fixed" )->add(
            ".", interval_t( dat.tp[analysis.lo], dat.tp[analysis.hi] + sample_tp ), "." );
          logger << "  emitted fixed analysis annotation under '" << par.emit_spt_label << "'\n";
        }
    }

  if ( par.do_ultradian && ! par.emit_cycle_label.empty() )
    {
      const int n_cycles = hd_emit_ultradian_cycles( edf, dat, ultra, analysis, par );
      logger << "  emitted " << n_cycles << " complete soft-ultradian cycle annotation(s) under '"
             << par.emit_cycle_label << "'\n";
    }

  // Strict mode retains the established contiguous-core contract exactly.
  // Gap-aware mode instead records every valid run and evaluates each metric
  // against its own usable subset without joining source-time across a gap.
  hd_lr_quality_t lr_quality;
  bool longrange_available = true;
  if ( par.do_longrange )
    {
      lr_quality = hd_lr_quality( dat, &analysis.physical_mask );
      int lr_core_lo = par.lr_policy == "strict" ? analysis.lo : lr_quality.first_valid;
      int lr_core_hi = par.lr_policy == "strict" ? analysis.hi : lr_quality.last_valid;

      if ( par.lr_policy == "gap-aware" )
        {
          if ( lr_quality.n_internal_gaps > 0 )
            logger << "  warning: HDSTATS longrange excluded "
                   << lr_quality.n_invalid / dat.Fs << " invalid seconds (longest internal gap "
                   << lr_quality.longest_internal_gap / dat.Fs << " seconds; "
                   << lr_quality.runs.size() << " valid runs)\n";
          if ( lr_quality.n_valid == 0 && par.lr_fail_on_empty )
            Helper::halt( "HDSTATS longrange has no valid observations" );
          if ( lr_quality.n_valid == 0 && make_pc_features )
            longrange_available = false;
        }
      else
        {
          while ( lr_core_lo <= analysis.hi && ! dat.valid[lr_core_lo] ) ++lr_core_lo;
          while ( lr_core_hi >= lr_core_lo && ! dat.valid[lr_core_hi] ) --lr_core_hi;
          std::string reason;
          if ( lr_core_lo > lr_core_hi ) reason = "HDSTATS longrange requires valid hypnodensity samples";
          if ( reason.empty() )
            for ( int i = lr_core_lo; i <= lr_core_hi; i++ )
              if ( ! dat.valid[i] )
                { reason = "HDSTATS longrange requires one contiguous valid core (no internal invalid samples)"; break; }
          const int lr_n = reason.empty() ? lr_core_hi - lr_core_lo + 1 : 0;
          if ( reason.empty() && lr_n / dat.Fs < par.lr_min_duration_sec )
            reason = "HDSTATS longrange requires at least "
              + Helper::dbl2str( par.lr_min_duration_sec / 3600.0 )
              + " hours of contiguous valid hypnodensity data";
          if ( reason.empty() ) for ( double x : par.lr_lags_sec )
            {
              const int lag = std::max( 1, (int)std::round( x * dat.Fs ) );
              if ( lr_n < 3 * lag )
                { reason = "HDSTATS longrange core is too short for lr-lags=" + Helper::dbl2str(x) + " seconds"; break; }
            }
          if ( reason.empty() ) for ( double x : par.lr_windows_sec )
            {
              const int width = std::max( 1, (int)std::round( x * dat.Fs ) );
              if ( lr_n / width < 4 )
                { reason = "HDSTATS longrange core is too short for lr-windows=" + Helper::dbl2str(x) + " seconds"; break; }
            }
          if ( !reason.empty() )
            {
              if ( !make_pc_features || par.lr_fail_on_empty ) Helper::halt(reason);
              longrange_available = false;
              logger << "  warning: " << reason << "; long-range PC inputs unavailable\n";
            }
          else logger << "  longrange valid core: " << lr_n / dat.Fs / 3600.0 << " hours\n";
        }
    }

  // Layer 3: detect transitions.  In an onset-selected mode no stability, support window,
  // event or bout is allowed to borrow a sample outside the analysis interval.
  hd_trans_t tr;
  tr.detect( dat, der, par, analysis.mode != "recording" ? &analysis.valid_mask : NULL );

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

  // All substantive summaries are limited to the selected physical interval.
  // HD_SOL remains referenced to recording start by design.
  const std::vector<bool> & all_mask = analysis.valid_mask;
  hdstats_pc::feature_map_t pc_features;
  if ( make_pc_features ) pc_features = hd_pc_features( dat, der, tr, hm, analysis, par );

  // Operational validity/analysis diagnostics are opt-in QC, separate from
  // the ordinary night-level phenotype summaries below.
  double twt_min = 0.0;
  for ( int i = analysis.lo; i <= analysis.hi; ++i )
    if ( analysis.valid_mask[i] ) twt_min += dat.p[0][i] / dat.Fs / 60.0;
  if ( par.qc )
  {
    writer.level( 1, "QC" );
    const int n_valid = (int)std::count( dat.valid.begin(), dat.valid.end(), true );
    writer.value( "N_VALID",   n_valid );
    writer.value( "PCT_VALID", dat.N > 0 ? 100.0 * n_valid / dat.N : 0.0 );
    if ( dat.n_invalid_lead  > 0 ) writer.value( "N_INVALID_LEAD",  dat.n_invalid_lead );
    if ( dat.n_invalid_trail > 0 ) writer.value( "N_INVALID_TRAIL", dat.n_invalid_trail );
    if ( dat.n_invalid_mid   > 0 ) writer.value( "N_INVALID_MID",   dat.n_invalid_mid );
    writer.value( "HD_ANALYSIS", analysis.mode );
    writer.value( "HD_ANALYSIS_START_SEC", analysis.lo / dat.Fs );
    writer.value( "HD_ANALYSIS_STOP_SEC", ( analysis.hi + 1 ) / dat.Fs );
    writer.value( "HD_ANALYSIS_SEC", ( analysis.hi - analysis.lo + 1 ) / dat.Fs );
    writer.value( "HD_EXCLUDED_PRE_SEC", analysis.lo / dat.Fs );
    writer.value( "HD_EXCLUDED_POST_SEC", ( dat.N - analysis.hi - 1 ) / dat.Fs );
    writer.value( "HD_ANALYSIS_INVALID_SEC", analysis.n_invalid / dat.Fs );
    writer.value( "HD_ANALYSIS_FALLBACK", analysis.fallback_to_recording ? 1 : 0 );
    writer.unlevel( "QC" );
  }
  if ( emit_main_metrics ) writer.value( "HD_TWT", twt_min );

  // Ultradian coverage is reported independently of excursion amplitude:
  // unsupported depth (insufficient local sleep evidence) is distinct from a
  // supported but low-amplitude half-cycle.
  if ( par.do_ultradian )
    {
      int n_supported = 0, n_phase = 0;
      int n_phase_runs = 0;
      bool previous_phase = false;
      double sleep_mass = 0.0, amp_mass = 0.0;
      double covered_sleep_mass = 0.0, raw_amp_mass = 0.0;
      double phase_advance = 0.0;
      for ( int i = analysis.lo; i <= analysis.hi; ++i )
        if ( analysis.valid_mask[i] )
          {
            const double sleep = 1.0 - dat.p[0][i];
            sleep_mass += sleep;
            amp_mass += ultra.amp[i]; // unsupported / no-half-cycle is zero by definition
            if ( ultra.supported[i] ) ++n_supported;
            const bool defined = std::isfinite( ultra.phase[i] );
            if ( defined )
              {
                ++n_phase;
                covered_sleep_mass += sleep;
                raw_amp_mass += ultra.raw_amp[i] * sleep;
                if ( !previous_phase ) ++n_phase_runs;
                else
                  {
                    // U_PHASE is wrapped at troughs.  Only accumulate within
                    // an uninterrupted defined run; never bridge a gap.
                    const double previous = ultra.phase[i-1];
                    double dphi = ultra.phase[i] - previous;
                    if ( dphi < -M_PI ) dphi += 2.0 * M_PI;
                    if ( dphi > 0.0 && std::isfinite(dphi) ) phase_advance += dphi;
                  }
              }
            previous_phase = defined;
          }
        else previous_phase = false;
      const int n_analysis_valid = (int)std::count( analysis.valid_mask.begin(), analysis.valid_mask.end(), true );
      const double nan = std::numeric_limits<double>::quiet_NaN();
      // Keep the ordinary whole-night panel intentionally compact.  These
      // three are the primary ultradian phenotypes; the supporting counts and
      // conservative strict-turning-point diagnostic are verbose QC only.
      if ( emit_main_metrics ) writer.value( "U_MAG", sleep_mass > 0.0 ? amp_mass / sleep_mass : nan );
      if ( emit_main_metrics ) writer.value( "U_COVERAGE", sleep_mass > 0.0 ? covered_sleep_mass / sleep_mass : nan );
      // A valid night with no phase has zero supported progression.  Missing
      // input, by contrast, cannot assess an effective cycle count.
      if ( emit_main_metrics ) writer.value( "U_NCYC", n_analysis_valid > 0 ? phase_advance / (2.0 * M_PI) : nan );
      if ( make_pc_features )
        {
          hd_pc_add( &pc_features, "U_MAG", sleep_mass > 0.0 ? amp_mass / sleep_mass : nan );
          hd_pc_add( &pc_features, "U_COVERAGE", sleep_mass > 0.0 ? covered_sleep_mass / sleep_mass : nan );
          hd_pc_add( &pc_features, "U_NCYC", n_analysis_valid > 0 ? phase_advance / (2.0 * M_PI) : nan );
        }
      if ( par.qc && par.verbose )
        {
          writer.level( 1, "QC" );
          writer.value( "U_N_HALF_CYCLES", ultra.n_half_cycles );
          // This intentionally excludes edge-partial and gap-split components;
          // U_NCYC above is the primary fractional phase-progression measure.
          writer.value( "U_STRICT_NCYC", ultra.n_complete_cycles );
          writer.value( "U_N_SUPPORTED", n_supported );
          writer.value( "U_N_PHASE", n_phase );
          writer.value( "U_PCT_SUPPORTED", n_analysis_valid > 0 ? 100.0 * n_supported / n_analysis_valid : 0.0 );
          writer.value( "U_PCT_PHASE", n_analysis_valid > 0 ? 100.0 * n_phase / n_analysis_valid : 0.0 );
          writer.value( "U_VALID_SEC", n_analysis_valid / dat.Fs );
          writer.value( "U_SLEEP_SEC", sleep_mass / dat.Fs );
          writer.value( "U_PHASE_SEC", n_phase / dat.Fs );
          writer.value( "U_N_PHASE_RUNS", n_phase_runs );
          writer.value( "U_MEAN_RAW", covered_sleep_mass > 0.0 ? raw_amp_mass / covered_sleep_mass : nan );
          writer.unlevel( "QC" );
        }
    }

  // HDSTATS-PC has a fixed ultradian schema, so calculate its selected UVAR
  // terms whenever a PC line/projection is requested.  The ordinary UVAR
  // table remains controlled by ultra-condition=.
  if ( par.do_ultradian && ( par.do_ultra_condition || make_pc_features ) )
    {
      if ( par.qc && par.do_ultra_condition )
        {
          writer.level( 1, "QC" );
          writer.value( "U_CONDITION_WEIGHTED", par.ultra_condition_weighted ? 1 : 0 );
          writer.unlevel( "QC" );
        }
      hd_write_ultra_conditioned_metrics( dat, der, ultra, analysis, par,
                                          make_pc_features ? &pc_features : NULL,
                                          par.do_ultra_condition && emit_main_metrics );
    }

  // As above, the PC UQ features are calculated even when the descriptive
  // output table is suppressed with ultra-regions=F.
  if ( par.do_ultradian && ( par.do_ultra_regions || make_pc_features ) )
    hd_write_ultra_phase_regions( dat, der, ultra, analysis, par,
                                  make_pc_features ? &pc_features : NULL,
                                  par.do_ultra_regions && emit_main_metrics );

  // Global summary
  if ( emit_main_metrics )
    write_hdstats( der, tr, all_mask, all_mask, par, dat, par.do_hd_metrics ? &hm : NULL, &analysis );

  if ( par.do_longrange )
    {
      if ( par.qc )
        {
          writer.level( 1, "QC" );
          writer.value( "LR_POLICY", par.lr_policy );
          writer.value( "LR_N_VALID_RUNS", (int)lr_quality.runs.size() );
          writer.value( "LR_VALID_SEC", lr_quality.n_valid / dat.Fs );
          writer.value( "LR_INVALID_SEC", lr_quality.n_invalid / dat.Fs );
          writer.value( "LR_LONGEST_INTERNAL_GAP_SEC", lr_quality.longest_internal_gap / dat.Fs );
          writer.value( "LR_N_INTERNAL_GAPS", lr_quality.n_internal_gaps );
          if ( par.lr_policy == "strict" && longrange_available )
            writer.value( "LR_CORE_SEC", ( lr_quality.last_valid - lr_quality.first_valid + 1 ) / dat.Fs );
          writer.unlevel( "QC" );
        }
      if ( longrange_available )
        write_hdstats_longrange( dat, der, par, lr_quality,
                                 make_pc_features ? &pc_features : NULL,
                                 emit_main_metrics );
    }

  // Level-1 HDSTATS-PC exchange: one compact line per individual.  This is
  // independent of STOUT/destrat and can be concatenated directly for the
  // level-2 reference fit. A projection is evaluated from the same feature
  // map, so it always shares the current HDSTATS parameterization.
  if ( ! par.pc_out.empty() || ! par.pc_proj.empty() )
    {
      if ( ! par.pc_out.empty() )
        {
          hdstats_pc::write_features( par.pc_out, edf.id, pc_features );
          logger << "  wrote HDSTATS-PC features for " << edf.id << "\n";
        }
      if ( ! par.pc_proj.empty() )
        {
          // This uses the portable V3 reader directly, rather than the
          // projection cache, so API clients and HDSTATS interpret optional
          // model annotations with exactly the same parser.
          if ( emit_model_desc )
            {
              hdstats_model::model_t model;
              std::string error;
              if ( ! hdstats_model::read( par.pc_proj, &model, &error ) ) Helper::halt( error );
              writer.level( "MODEL_DESC", "MODEL" );
              writer.value( "N_PANEL", (int)model.panels.size() );
              writer.value( "N_DESC", (int)model.descriptions.size() );
              writer.value( "N_EXT_ASSOC", (int)model.external_associations.size() );
              writer.value( "N_PC_CORR", (int)model.pc_correlations.size() );
              writer.value( "N_FTR_CORR", (int)model.feature_correlations.size() );
              for ( const auto & description : model.descriptions )
                {
                  writer.level( description.entity_type, "ENTITY" );
                  writer.level( description.entity_id, "MODEL_ID" );
                  writer.value( "DESC", description.text );
                  writer.unlevel( "MODEL_ID" );
                  writer.unlevel( "ENTITY" );
                }
              for ( const auto & association : model.external_associations )
                {
                  writer.level( association.entity_type, "ENTITY" );
                  writer.level( association.entity_id, "MODEL_ID" );
                  writer.level( association.factor, "FACTOR" );
                  writer.value( "SCORE", association.score );
                  writer.unlevel( "FACTOR" );
                  writer.unlevel( "MODEL_ID" );
                  writer.unlevel( "ENTITY" );
                }
              for ( const auto & correlation : model.pc_correlations )
                {
                  writer.level( correlation.first, "PC1" );
                  writer.level( correlation.second, "PC2" );
                  writer.value( "CORR", correlation.value );
                  writer.unlevel( "PC2" );
                  writer.unlevel( "PC1" );
                }
              for ( const auto & correlation : model.feature_correlations )
                {
                  writer.level( correlation.first, "FTR1" );
                  writer.level( correlation.second, "FTR2" );
                  writer.value( "CORR", correlation.value );
                  writer.unlevel( "FTR2" );
                  writer.unlevel( "FTR1" );
                }
              writer.unlevel( "MODEL" );
            }
          // Domain-coded PC IDs deliberately avoid PANEL/PC structural missingness:
          // every declared component is emitted, with SCORE=NA when its source
          // panel cannot be projected for this individual.
          for ( const auto & component : hdstats_pc::project_components(
                  par.pc_proj, pc_features, par.pc_min_loading ) )
            {
              writer.level( component.id, "PC" );
              writer.value( "SCORE", component.score );
              for ( const auto & contribution : component.contributions )
                {
                  writer.level( contribution.feature, "FTR" );
                  writer.value( "VAL", contribution.value );
                  writer.value( "RAW", contribution.raw_value );
                  writer.unlevel( "FTR" );
                }
              writer.unlevel( "PC" );
            }
        }
    }

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
	  if ( par.do_ultradian )
	    {
	      if ( std::isfinite( ultra.depth[i] ) ) writer.value( "U_DEPTH", ultra.depth[i] );
	      if ( std::isfinite( ultra.phase[i] ) ) writer.value( "U_PHASE", ultra.phase[i] );
	      writer.value( "U_PHASE_VALID", std::isfinite( ultra.phase[i] ) ? 1 : 0 );
	      writer.value( "U_AMP", ultra.amp[i] );
	      writer.value( "U_AMP_RAW", ultra.raw_amp[i] );
	    }
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
		    mask5[i] = analysis.valid_mask[i];
		}
	    }
	  ++ii;
	}

      // 3-state uses the same time mask (same sample positions)
      writer.level( lv, globals::annot_strat );
      if ( emit_main_metrics )
        write_hdstats( der, tr, mask5, mask5, par, dat, NULL, &analysis );
      writer.unlevel( globals::annot_strat );
    }
}
