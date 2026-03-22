
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

#include "actig.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "param.h"
#include "annot/annot.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "db/db.h"
#include "defs/defs.h"
#include "miscmath/miscmath.h"

#include <cmath>
#include <numeric>
#include <algorithm>

extern writer_t writer;
extern logger_t logger;


// -----------------------------------------------------------------------
//
// Helper: find sliding window of W bins with min or max mean
//  (wraps around for circadian profiles)
//
// -----------------------------------------------------------------------

// Median of a (copied) vector
static double vec_median( std::vector<double> v )
{
  if ( v.empty() ) return 0.0;
  const size_t n = v.size();
  std::sort( v.begin(), v.end() );
  return ( n % 2 == 0 ) ? ( v[n/2 - 1] + v[n/2] ) / 2.0 : v[n/2];
}

static double vec_mean( const std::vector<double> & v )
{
  if ( v.empty() ) return 0.0;
  double s = 0;
  for ( double x : v ) s += x;
  return s / v.size();
}

static double vec_sd( const std::vector<double> & v )
{
  if ( v.size() < 2 ) return 0.0;
  const double mn = vec_mean( v );
  double s = 0;
  for ( double x : v ) { double d = x - mn; s += d * d; }
  return std::sqrt( s / ( v.size() - 1 ) );
}

static double vec_mad( const std::vector<double> & v , const double med )
{
  if ( v.empty() ) return 0.0;
  std::vector<double> abs_dev( v.size() );
  for ( size_t i = 0 ; i < v.size() ; i++ )
    abs_dev[i] = std::fabs( v[i] - med );
  return vec_median( abs_dev );
}

struct bout_counts_t {
  std::vector<int> sleep_epochs;
  std::vector<int> wake_epochs;
};

static bout_counts_t extract_bouts( const std::vector<bool> & is_gap ,
                                    const std::vector<bool> & is_wake ,
                                    const std::vector<bool> & in_period ,
                                    int first , int last )
{
  bout_counts_t b;
  if ( first == -1 ) return b;
  int i = first;
  while ( i <= last )
    {
      if ( is_gap[i] || !in_period[i] ) { ++i; continue; }
      const bool wake = is_wake[i];
      int j = i + 1;
      while ( j <= last && !is_gap[j] && in_period[j] && is_wake[j] == wake ) ++j;
      if ( wake ) b.wake_epochs.push_back( j - i );
      else        b.sleep_epochs.push_back( j - i );
      i = j;
    }
  return b;
}

// percentile of bout durations in minutes; returns 0 if empty
static double pct_min( const std::vector<int> & runs , double q , double epoch_sec )
{
  if ( runs.empty() ) return 0.0;
  std::vector<double> mins;
  mins.reserve( runs.size() );
  for ( int r : runs ) mins.push_back( r * epoch_sec / 60.0 );
  return MiscMath::percentile( mins , q );
}

// Shannon conditional entropy H(X_t | X_{t-1}) for a 2-state (sleep/wake)
// sequence within in_period, skipping gaps and out-of-period epochs.
// Returns entropy in bits; 0 if fewer than 2 valid epochs.
static double transition_entropy( const std::vector<bool> & is_gap ,
                                  const std::vector<bool> & is_wake ,
                                  const std::vector<bool> & in_period ,
                                  int first , int last )
{
  if ( first == -1 ) return 0.0;
  int n[2][2] = { {0,0}, {0,0} };
  int prev = -1;
  for ( int i = first ; i <= last ; i++ )
    {
      if ( is_gap[i] || !in_period[i] ) continue;
      const int s = is_wake[i] ? 1 : 0;
      if ( prev >= 0 ) n[prev][s]++;
      prev = s;
    }
  const int total = n[0][0] + n[0][1] + n[1][0] + n[1][1];
  if ( total == 0 ) return 0.0;
  double h = 0.0;
  for ( int i = 0 ; i < 2 ; i++ )
    {
      const int row = n[i][0] + n[i][1];
      if ( row == 0 ) continue;
      const double pi = row / (double)total;
      for ( int j = 0 ; j < 2 ; j++ )
        if ( n[i][j] > 0 )
          {
            const double pij = n[i][j] / (double)row;
            h -= pi * pij * std::log2( pij );
          }
    }
  return h;
}

// Shannon entropy of the empirical bout-length distribution (bits).
// Measures variability of bout lengths; 0 if empty or all same length.
static double rl_entropy( const std::vector<int> & bouts )
{
  if ( bouts.size() < 2 ) return 0.0;
  std::map<int,int> cnt;
  for ( int r : bouts ) cnt[r]++;
  const double N = (double)bouts.size();
  double h = 0.0;
  for ( auto & kv : cnt )
    {
      const double p = kv.second / N;
      h -= p * std::log2( p );
    }
  return h;
}

// Transition probabilities P(sleep→wake) and P(wake→sleep) within a period.
// Uses the ML estimator with last-epoch correction and a Bayesian lambda
// smoothing factor (paper: Springer BMC Med Res Methodol 2024; default λ=1,
// corresponding to a uniform / Laplace-smoothing prior).
//
// Formula (s=1): TP_SW = (n_SW + λ) / (T_sleep − I(last=sleep) + λ)
//                TP_WS = (n_WS + λ) / (T_wake  − I(last=wake)  + λ)
//
// Gap epochs reset the transition chain (no transition counted across a gap)
// but all valid in-period epochs contribute to epoch counts.
struct tp_pair_t { double tp_sw = 0.0 , tp_ws = 0.0; };

static tp_pair_t transition_probs( const std::vector<bool> & is_gap ,
                                   const std::vector<bool> & is_wake ,
                                   const std::vector<bool> & in_period ,
                                   int first , int last ,
                                   double lambda )
{
  tp_pair_t tp;
  if ( first == -1 ) return tp;
  int n_sw = 0 , n_ws = 0;
  int t_sleep = 0 , t_wake = 0;
  int last_state = -1;   // state of last valid (non-gap, in-period) epoch
  int prev = -1;         // previous state in current unbroken chain; reset on gap
  for ( int i = first ; i <= last ; i++ )
    {
      if ( is_gap[i] || !in_period[i] ) { prev = -1; continue; }
      const int s = is_wake[i] ? 1 : 0;
      if ( s == 0 ) t_sleep++; else t_wake++;
      if ( prev == 0 && s == 1 ) n_sw++;
      if ( prev == 1 && s == 0 ) n_ws++;
      prev = s;
      last_state = s;
    }
  // MLE last-epoch correction: last in-period epoch is observed but has no
  // successor within the period, so subtract it from the at-risk denominator.
  const double denom_S = t_sleep - ( last_state == 0 ? 1 : 0 ) + lambda;
  const double denom_W = t_wake  - ( last_state == 1 ? 1 : 0 ) + lambda;
  if ( denom_S > 0 ) tp.tp_sw = ( n_sw + lambda ) / denom_S;
  if ( denom_W > 0 ) tp.tp_ws = ( n_ws + lambda ) / denom_W;
  return tp;
}

static double first_day_boundary_sec( const clocktime_t & startdatetime ,
                                      const int anchor_hr ,
                                      const double total_sec )
{
  double first_boundary_sec = total_sec;

  if ( startdatetime.valid )
    {
      const double start_tod_sec = startdatetime.h * 3600.0
        + startdatetime.m * 60.0
        + startdatetime.s;
      const double anchor_sec = anchor_hr * 3600.0;
      first_boundary_sec = anchor_sec - start_tod_sec;
      if ( first_boundary_sec <= 0 )
        first_boundary_sec += 86400.0;
    }

  return first_boundary_sec;
}

struct actig_day_qc_t
{
  double act_med = 0.0;
  double act_iqr = 0.0;
  double act_sd = 0.0;
  double act_mad = 0.0;
  double act_min = 0.0;
  double act_max = 0.0;
  double act_p05 = 0.0;
  double act_p95 = 0.0;
  double flat_frac = 0.0;
  double lowvar_cv = 0.0;
  double nearfloor_frac = 0.0;
  double nonzero_frac = 0.0;
  int active_epoch_n = 0;
  double longest_quiet_run_min = 0.0;
  double longest_lowvar_run_min = 0.0;
  double longest_sleep_bout_min = 0.0;
  double longest_wake_bout_min = 0.0;
  int wake_run_n = 0;
  int sleep_run_n = 0;
  int tech_flag_n = 0;
  bool flag_flat = false;
  bool flag_lowvar = false;
  bool flag_nearfloor = false;
  bool flag_allsleep = false;
  bool flag_allwake = false;
  bool warn_highsleep = false;
  bool warn_longsleep = false;
  bool warn_lowwakeruns = false;
  bool warn_implausible = false;
  bool day_warn = false;
  bool day_excluded_tech = false;
  bool day_excluded_extreme = false;
  bool day_excluded = false;
  bool day_ok = false;
};


// Fine-step sliding L5/M10 window over a 24 h profile.
//
// The profile is treated as piecewise-constant with bin width bin_min.
// The window is placed at positions t = 0, step_hr, 2*step_hr, ...
// and its mean is computed as the exact weighted integral over the
// piecewise-constant profile, handling circular wrap-around.
//
// step_min == 60  reproduces the traditional integer-hour-step behaviour.
//
// bin_min : minutes per profile bin (1 for the minute-resolution NP profile).
//   The profile must cover exactly 24 hours (profile.size() * bin_min == 1440).
//   window_hours and result_onset_hr are always in fractional hours.
//   step_min controls search resolution; meaningful down to bin_min.
static void find_window_fine( const std::vector<double> & profile ,
			      double bin_min ,
			      int    window_hours ,
			      bool   find_min ,
			      double * result_mean ,
			      double * result_onset_hr ,   // fractional hours [0,24)
			      int    step_min = 1 )
{
  const int    p      = (int)profile.size();
  const double bhr    = bin_min / 60.0;          // bin width in hours
  const double W      = (double)window_hours;    // window width in hours

  if ( p <= 0 || bhr <= 0 || W <= 0 || W >= 24.0 )
    {
      *result_mean     = 0;
      *result_onset_hr = 0;
      return;
    }

  const double step_hr = step_min / 60.0;
  const int    n_steps = (int)std::round( 24.0 / step_hr );

  double best       = find_min ? 1e300 : -1e300;
  double best_onset = 0.0;

  for ( int si = 0 ; si < n_steps ; si++ )
    {
      const double t   = si * step_hr;  // window onset (fractional hours)
      const double end = t + W;

      double sum = 0.0;

      for ( int h = 0 ; h < p ; h++ )
	{
	  const double h0 = h * bhr , h1 = h0 + bhr;
	  double overlap;

	  if ( end <= 24.0 )
	    {
	      overlap = std::max( 0.0 , std::min( h1 , end ) - std::max( h0 , t ) );
	    }
	  else
	    {
	      // window wraps past midnight: [t,24) ∪ [0, end-24)
	      double ow1 = std::max( 0.0 , std::min( h1 , 24.0    ) - std::max( h0 , t   ) );
	      double ow2 = std::max( 0.0 , std::min( h1 , end-24.0 ) - std::max( h0 , 0.0 ) );
	      overlap = ow1 + ow2;
	    }

	  sum += profile[h] * overlap;
	}

      double m = sum / W;

      if ( ( find_min && m < best ) || ( !find_min && m > best ) )
	{
	  best       = m;
	  best_onset = t;
	}
    }

  *result_mean     = best;
  *result_onset_hr = best_onset;
}


// -----------------------------------------------------------------------
//
// ACTIG : nonparametric circadian metrics and wake/sleep scoring
//
//  Gap handling: gaps are identified via timestamp discontinuities in
//  ptimepoints() after wholetrace() extraction.  Any epoch bin whose
//  sample count falls below gap-min-pct % of the expected count is
//  flagged as a gap and excluded from all calculations.  Gaps must be
//  communicated to ACTIG via the standard MASK/RE mechanism prior to
//  calling this command; no annotation-based gap parameter is needed.
//
// -----------------------------------------------------------------------

void actig::actig( edf_t & edf , param_t & param )
{

  //
  // Signal (optional in prescored mode)
  //

  const bool do_score        = param.has( "score" );
  const bool do_score_period = param.has( "score-period" );
  const bool do_any_score    = do_score || do_score_period;
  const bool do_prescored    = param.has( "prescored" );
  const bool do_verbose      = param.has( "verbose" );
  const bool do_investigate  = param.has( "investigate" );

  // Epoch-level prescored labels (S/W by default)
  std::string presc_sleep_label = "S";
  std::string presc_wake_label  = "W";
  if ( do_prescored && ! param.empty( "prescored" ) )
    {
      std::vector<std::string> tok = Helper::parse( param.value( "prescored" ) , "," );
      if ( tok.size() != 2 || tok[0] == "" || tok[1] == "" )
	Helper::halt( "ACTIG: prescored must be empty or exactly two comma-delimited labels: sleep-epoch,wake-epoch" );
      presc_sleep_label = tok[0];
      presc_wake_label  = tok[1];
    }

  if ( do_score && do_score_period )
    Helper::halt( "ACTIG: cannot specify both score and score-period" );
  if ( do_any_score && do_prescored )
    Helper::halt( "ACTIG: cannot specify both score/score-period and prescored" );

  // need_signal: signal required only when actively scoring
  // have_signal: signal available for NP metrics and gap detection
  //   (required for scoring, or optionally supplied in prescored/annotation mode)
  const bool need_signal = do_any_score;
  // sig=* is injected by Luna as a default; only treat sig as explicitly provided
  // when it names a specific channel (not the wildcard)
  const bool sig_explicit = param.has( "sig" ) && ! param.empty( "sig" )
                            && param.value( "sig" ) != "*";
  const bool have_signal = do_any_score || sig_explicit;
  const std::string siglab = need_signal ? param.requires( "sig" )
                           : ( sig_explicit ? param.value( "sig" ) : "" );

  int slot = -1;
  double Fs = 0;

  if ( have_signal )
    {
      signal_list_t signals = edf.header.signal_list( siglab );

      if ( signals.size() != 1 )
	Helper::halt( "ACTIG requires exactly one signal" );

      if ( edf.header.is_annotation_channel( signals(0) ) )
	Helper::halt( "ACTIG: cannot use an annotation channel" );

      slot = signals(0);
      Fs = edf.header.sampling_freq( slot );

      if ( Fs <= 0 )
	Helper::halt( "ACTIG: invalid sampling rate" );
    }


  //
  // Options
  //

  // epoch length for activity binning (default 60s = 1 minute)
  const double epoch_sec = param.has( "epoch" ) ? param.requires_dbl( "epoch" ) : 60.0;
  if ( epoch_sec <= 0 )
    Helper::halt( "ACTIG: epoch must be > 0 seconds" );

  // bin size for NP metrics (default 60 minutes)
  const double np_bin_min = param.has( "bin" ) ? param.requires_dbl( "bin" ) : 60.0;
  if ( np_bin_min <= 0 )
    Helper::halt( "ACTIG: bin must be > 0 minutes" );

  // Restrict NP metrics to complete anchored 24 h days, dropping
  // leading and trailing partial days when a valid start datetime exists.
  const bool np_full_days = ! param.has( "np-full-days" ) || param.yesno( "np-full-days" );

  // L5/M10 window sizes (hours)
  const int l_hours = param.has( "l" ) ? param.requires_int( "l" ) : 5;
  const int m_hours = param.has( "m" ) ? param.requires_int( "m" ) : 10;

  // NP sliding-window step: minutes per L5/M10 window position (default 1 min)
  // Use np-traditional to revert to original 60-min step (integer-hour onsets)
  const bool np_traditional = param.has( "np-traditional" );
  const int  np_step_min    = np_traditional ? 60
    : ( param.has( "np-step" ) ? param.requires_int( "np-step" ) : 1 );
  if ( np_step_min < 1 || np_step_min > 60 || 60 % np_step_min != 0 )
    Helper::halt( "ACTIG: np-step must be 1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 30, or 60" );

  // use sum rather than mean for epoch binning (e.g., for raw counts)
  const bool use_sum = param.has( "sum" );

  // scoring method
  const std::string method = param.has( "method" ) ? param.value( "method" ) : "luna";

  // threshold for threshold method
  const double thresh = param.has( "thresh" ) ? param.requires_dbl( "thresh" ) : -1;

  // labels for wake/sleep/gap output annotations
  // score-period defaults to WP/SP; score defaults to W/S
  const std::string wake_label  = param.has( "wake" )    ? param.value( "wake" )
                                : ( do_score_period ? "WP" : "W" );
  const std::string sleep_label = param.has( "sleep" )   ? param.value( "sleep" )
                                : ( do_score_period ? "SP" : "S" );

  // Sleep/wake period window annotations (SP/WP by default).
  // These define the sleep/wake windows for fragmentation metrics.
  // Looked up regardless of scoring mode; falls back to empirical if absent.
  const std::string period_sleep_label = param.has( "sleep-period" ) ? param.value( "sleep-period" ) : "SP";
  const std::string period_wake_label  = param.has( "wake-period" )  ? param.value( "wake-period" )  : "WP";

  // Bayesian smoothing for transition probability estimates (λ in ML formula).
  // λ=1 is a uniform/Laplace prior; avoids zero-count issues in short recordings.
  const double tp_lambda = param.has( "tp-lambda" ) ? param.requires_dbl( "tp-lambda" ) : 1.0;

  // Day boundary anchor hour for NP and per-day summaries.
  const int day_anchor_hr = param.has( "day-anchor" ) ? param.requires_int( "day-anchor" ) : 12;
  if ( day_anchor_hr < 0 || day_anchor_hr > 23 )
    Helper::halt( "ACTIG: day-anchor must be 0..23" );

  // minimum % of valid epochs per NP bin for the bin to be included
  const double gap_min_pct = param.has( "gap-min-pct" ) ? param.requires_dbl( "gap-min-pct" ) : 50.0;
  if ( gap_min_pct < 0 || gap_min_pct > 100 )
    Helper::halt( "ACTIG: gap-min-pct must be 0..100" );

  // minimum valid minutes per day for inclusion in day averages (default 16 h)
  const double day_min_valid_min = param.has( "day-min-valid" ) ? param.requires_dbl( "day-min-valid" ) : 960.0;
  if ( day_min_valid_min < 0 )
    Helper::halt( "ACTIG: day-min-valid must be >= 0" );

  // Day-level QC
  const bool qc_day               = ! param.has( "qc-day" )              || param.yesno( "qc-day" );
  const bool qc_exclude_flat      = ! param.has( "qc-exclude-flat" )     || param.yesno( "qc-exclude-flat" );
  const bool qc_exclude_lowvar    = ! param.has( "qc-exclude-lowvar" )   || param.yesno( "qc-exclude-lowvar" );
  const bool qc_exclude_nearfloor = ! param.has( "qc-exclude-nearfloor" )|| param.yesno( "qc-exclude-nearfloor" );
  const bool qc_warn_implausible  = ! param.has( "qc-warn-implausible" ) || param.yesno( "qc-warn-implausible" );
  const bool qc_warn_longsleep    = ! param.has( "qc-warn-longsleep" )   || param.yesno( "qc-warn-longsleep" );
  const bool qc_warn_highsleep    = ! param.has( "qc-warn-highsleep" )   || param.yesno( "qc-warn-highsleep" );
  const std::string qc_exclude_label  = param.has( "qc-exclude-out" )  ? param.value( "qc-exclude-out" )  : "QC_EXCLUDED";
  const std::string qc_extreme_label  = param.has( "qc-extreme-out" )  ? param.value( "qc-extreme-out" )  : "QC_EXTREME";
  const bool qc_exclude_allsleep      = ! param.has( "qc-exclude-allsleep" ) || param.yesno( "qc-exclude-allsleep" );
  const bool qc_exclude_allwake       = ! param.has( "qc-exclude-allwake" )  || param.yesno( "qc-exclude-allwake" );
  const double qc_allsleep_th         = param.has( "qc-allsleep-th" ) ? param.requires_dbl( "qc-allsleep-th" ) : 0.98;
  const double qc_allwake_th          = param.has( "qc-allwake-th" )  ? param.requires_dbl( "qc-allwake-th" )  : 0.98;

  const double qc_flat_frac_th      = param.has( "qc-flat-frac-th" )      ? param.requires_dbl( "qc-flat-frac-th" )      : 0.80;
  const double qc_flat_delta_th     = param.has( "qc-flat-delta-th" )     ? param.requires_dbl( "qc-flat-delta-th" )     : 0.0;
  const double qc_lowvar_frac_th    = param.has( "qc-lowvar-frac-th" )    ? param.requires_dbl( "qc-lowvar-frac-th" )    : 0.80;
  const double qc_lowvar_cv_th      = param.has( "qc-lowvar-cv-th" )      ? param.requires_dbl( "qc-lowvar-cv-th" )      : 0.05;
  const double qc_nearfloor_frac_th = param.has( "qc-nearfloor-frac-th" ) ? param.requires_dbl( "qc-nearfloor-frac-th" ) : 0.85;
  const double qc_nearfloor_q_th    = param.has( "qc-nearfloor-q-th" )    ? param.requires_dbl( "qc-nearfloor-q-th" )    : 0.05;
  const int    qc_min_active_epochs = param.has( "qc-min-active-epochs" ) ? param.requires_int( "qc-min-active-epochs" ) : 24;

  const double qc_warn_sleep_pct    = param.has( "qc-warn-sleep-pct" )    ? param.requires_dbl( "qc-warn-sleep-pct" )    : 0.85;
  const double qc_warn_longsleep_h  = param.has( "qc-warn-longsleep-h" )  ? param.requires_dbl( "qc-warn-longsleep-h" )  : 16.0;
  const double qc_warn_max_sleep_h  = param.has( "qc-warn-max-sleep-h" )  ? param.requires_dbl( "qc-warn-max-sleep-h" )  : 20.0;
  const int    qc_warn_low_wakeruns = param.has( "qc-warn-low-wakeruns" ) ? param.requires_int( "qc-warn-low-wakeruns" ) : 1;

  auto check_frac01 = [&]( const double x , const std::string & name )
  {
    if ( x < 0.0 || x > 1.0 )
      Helper::halt( "ACTIG: " + name + " must be in 0..1" );
  };
  check_frac01( qc_flat_frac_th      , "qc-flat-frac-th" );
  check_frac01( qc_lowvar_frac_th    , "qc-lowvar-frac-th" );
  check_frac01( qc_nearfloor_frac_th , "qc-nearfloor-frac-th" );
  check_frac01( qc_nearfloor_q_th    , "qc-nearfloor-q-th" );
  check_frac01( qc_warn_sleep_pct    , "qc-warn-sleep-pct" );
  if ( qc_min_active_epochs < 0 ) Helper::halt( "ACTIG: qc-min-active-epochs must be >= 0" );
  if ( qc_warn_low_wakeruns < 0 ) Helper::halt( "ACTIG: qc-warn-low-wakeruns must be >= 0" );
  if ( qc_warn_longsleep_h < 0.0 ) Helper::halt( "ACTIG: qc-warn-longsleep-h must be >= 0" );
  if ( qc_warn_max_sleep_h < qc_warn_longsleep_h )
    Helper::halt( "ACTIG: qc-warn-max-sleep-h must be >= qc-warn-longsleep-h" );

  // Threshold method smoothing: triangular window half-width in minutes (0 = off)
  const double thresh_smooth_min  = param.has( "thresh-smooth" ) ? param.requires_dbl( "thresh-smooth" ) : 0.0;
  if ( thresh_smooth_min < 0 )
    Helper::halt( "ACTIG: thresh-smooth must be >= 0" );

  // Luna method parameters
  // score-period uses coarser defaults to produce contiguous period blocks
  const double luna_smooth_min    = param.has( "smooth" )      ? param.requires_dbl( "smooth" )
                                  : ( do_score_period ? 30.0  :  1.0  );
  const double luna_burst_min     = param.has( "burst" )       ? param.requires_dbl( "burst" )
                                  : ( do_score_period ? 10.0  :  3.0  );
  const double luna_burst_z       = param.has( "burst-z" )     ? param.requires_dbl( "burst-z" )     :  0.5;
  const double luna_quiet_z       = param.has( "quiet-z" )     ? param.requires_dbl( "quiet-z" )     : -0.5;
  const double luna_active_frac   = param.has( "active-frac" ) ? param.requires_dbl( "active-frac" ) :  0.20;
  const double luna_min_sleep_min = param.has( "min-sleep" )   ? param.requires_dbl( "min-sleep" )
                                  : ( do_score_period ? 60.0  :  1.0  );
  const double luna_max_gap_min   = param.has( "max-gap" )     ? param.requires_dbl( "max-gap" )
                                  : ( do_score_period ? 10.0  :  1.0  );
  const double luna_min_wake_min  = param.has( "min-wake" )    ? param.requires_dbl( "min-wake" )
                                  : ( do_score_period ? 20.0  :  1.0  );

  // add diagnostic channels (_Z, _L, _D, _CS, _SW) to the EDF (luna method only)
  const bool do_channels = param.has( "channels" );

  // Local sleep debt
  const bool        do_debt          = param.has( "debt" );
  const std::string debt_target_str  = param.has( "debt-target" )   ? param.value( "debt-target" )        : "";
  const int         debt_R           = param.has( "debt-recent" )    ? param.requires_int( "debt-recent" )  : 2;
  const int         debt_B           = param.has( "debt-base" )      ? param.requires_int( "debt-base" )    : 7;
  const int         debt_min_base    = param.has( "debt-min-base" )  ? param.requires_int( "debt-min-base" ): 3;
  const int         debt_min_z       = param.has( "debt-min-z" )     ? param.requires_int( "debt-min-z" )   : 5;
  const double      debt_w           = param.has( "debt-w" )         ? param.requires_dbl( "debt-w" )       : 0.5;

  // All-nights normalisation (Z-score target night vs all other valid nights)
  const bool        do_norm          = param.has( "norm" );
  const std::string norm_target_str  = param.has( "norm-target" )    ? param.value( "norm-target" )
                                     : ( param.has( "debt-target" )  ? param.value( "debt-target" ) : "" );
  const int         norm_min         = param.has( "norm-min" )        ? param.requires_int( "norm-min" )    : 5;
  const double      norm_w           = param.has( "norm-w" )          ? param.requires_dbl( "norm-w" )      : 0.5;


  if ( do_verbose )
    {
      logger << "  settings: sig=" << ( need_signal ? siglab : "<prescored>" );
      if ( need_signal ) logger << ", sr=" << Fs << " Hz";
      logger << ", epoch=" << epoch_sec << "s"
	     << ", bin=" << np_bin_min << "m\n";
      logger << "  settings: L=" << l_hours << "h"
	     << ", M=" << m_hours << "h"
	     << ", score=" << ( do_score ? "yes" : ( do_score_period ? "period" : "no" ) )
	     << ", method=" << method
	     << ", sum=" << ( use_sum ? "yes" : "no" ) << "\n";
      logger << "  day-min-valid=" << day_min_valid_min << " min\n";
      logger << "  day QC: " << ( qc_day ? "enabled" : "disabled" )
	     << ", tech flags require >=2 of [flat, lowvar, nearfloor]"
	     << ", warnings do not change inclusion\n";
    }


  //
  // Extract full signal + timestamps
  //
  //  Gaps (from prior MASK/RE or EDF+D discontinuities) appear as
  //  timestamp jumps in ptimepoints().  We detect them during binning.
  //

  interval_t whole = edf.timeline.wholetrace( true );

  std::vector<double>   sig_data;
  std::vector<uint64_t> sig_tp;
  const std::vector<double>   * d  = &sig_data;
  const std::vector<uint64_t> * tp = &sig_tp;

  if ( have_signal )
    {
      slice_t slice( edf , slot , whole );

      if ( slice.pdata()->empty() )
	{
	  logger << "  no data extracted\n";
	  return;
	}

      if ( slice.pdata()->size() != slice.ptimepoints()->size() )
	Helper::halt( "ACTIG: internal error: data/timepoint size mismatch" );

      sig_data = *slice.pdata();
      sig_tp   = *slice.ptimepoints();
    }


  //
  // Time-anchored epoch binning
  //
  //  Instead of binning by sequential index (which collapses gaps),
  //  we assign each sample to a bin based on its absolute timestamp
  //  offset from the recording start.  Bins with fewer than
  //  gap-min-pct % of expected samples are flagged as gaps.
  //

  const uint64_t recording_start_tp = whole.start;
  const uint64_t recording_total_tp = whole.stop - whole.start;
  const double   total_sec          = recording_total_tp / (double)globals::tp_1sec;
  const clocktime_t startdatetime( edf.header.startdate , edf.header.starttime );
  const double start_tod_sec = startdatetime.valid
    ? ( startdatetime.h * 3600.0 + startdatetime.m * 60.0 + startdatetime.s )
    : 0.0;

  const uint64_t epoch_tp = (uint64_t)( epoch_sec * globals::tp_1sec );
  if ( epoch_tp == 0 )
    Helper::halt( "ACTIG: epoch_sec too small for time-point resolution" );

  const int n_epochs = (int)( recording_total_tp / epoch_tp );

  if ( n_epochs == 0 )
    {
      logger << "  no epochs after binning\n";
      return;
    }

  // minimum samples required per epoch for it to be considered valid
  const int expected_spe = have_signal ? (int)( Fs * epoch_sec + 1e-9 ) : 1;
  const int min_spe      = have_signal ? std::max( 1 , (int)( gap_min_pct / 100.0 * expected_spe ) ) : 1;

  std::vector<double> epochs( n_epochs , 0.0 );
  std::vector<int>    epoch_count( n_epochs , 0 );
  std::vector<bool>   is_gap( n_epochs , have_signal );   // start all as gap; valid ones cleared below

  int n_valid_epochs = 0 , n_gap_epochs = 0;
  if ( have_signal )
    {
      for ( int i = 0 ; i < (int)d->size() ; i++ )
	{
	  if ( (*tp)[i] < recording_start_tp ) continue;  // safety guard
	  const uint64_t offset_tp = (*tp)[i] - recording_start_tp;
	  const int ebin = (int)( offset_tp / epoch_tp );
	  if ( ebin >= 0 && ebin < n_epochs )
	    {
	      epochs[ ebin ] += (*d)[i];
	      epoch_count[ ebin ]++;
	    }
	}

      for ( int e = 0 ; e < n_epochs ; e++ )
	{
	  if ( epoch_count[e] >= min_spe )
	    {
	      epochs[e] = use_sum ? epochs[e] : epochs[e] / epoch_count[e];
	      is_gap[e] = false;
	      n_valid_epochs++;
	    }
	  else
	    {
	      epochs[e] = 0.0;   // zero-fill so gap epochs never contribute
	      n_gap_epochs++;
	    }
	}
    }

  else
    {
      // prescored with no signal: no signal-based gap detection;
      // EDF+D holes will be flagged below via the record structure.
      n_valid_epochs = n_epochs;
      n_gap_epochs = 0;
      if ( do_verbose )
	logger << "  prescored mode: no signal provided; NP metrics will not be computed\n";
    }

  // EDF+D holes: epoch bins that fall between retained records after MASK+RE.
  // Detected from the EDF record structure (rec2tp) — works regardless of
  // whether a signal was provided.  Keep is_gap[e]=true for scoring (never
  // bridge across holes); EDF+D holes are already marked by EXCLUDED annotations.
  std::vector<bool> is_edf_hole( n_epochs , false );
  if ( !edf.header.continuous )
    {
      // Start with all bins as holes; clear any bin covered by a retained record.
      std::fill( is_edf_hole.begin() , is_edf_hole.end() , true );
      const uint64_t rec_dur_tp = edf.header.record_duration_tp;
      for ( auto & kv : edf.timeline.rec2tp )
	{
	  const uint64_t rec_start = kv.second;
	  if ( rec_start < recording_start_tp ) continue;
	  const uint64_t rec_stop  = rec_start + rec_dur_tp;
	  const int e0 = (int)( ( rec_start - recording_start_tp ) / epoch_tp );
	  const int e1 = (int)( ( rec_stop  - recording_start_tp - 1 ) / epoch_tp );
	  for ( int e = e0 ; e <= e1 && e < n_epochs ; e++ )
	    is_edf_hole[e] = false;
	}
      // In no-signal mode, is_gap[] starts all-false; mark holes as gap so
      // they are excluded from TST/TWT and the prescored annotation lookup.
      if ( !have_signal )
	for ( int e = 0 ; e < n_epochs ; e++ )
	  if ( is_edf_hole[e] ) { is_gap[e] = true; n_gap_epochs++; n_valid_epochs--; }
    }

  logger << "  recording: " << total_sec << " sec ("
	 << total_sec / 3600.0 << " hrs), "
	 << n_epochs << " epoch bins ("
	 << n_valid_epochs << " valid, "
	 << n_gap_epochs   << " gap)\n";


  // ---------------------------------------------------------------
  //
  // INVESTIGATE mode: signal pre-screen / parameter suggestion report
  //
  // ---------------------------------------------------------------

  if ( do_investigate )
    {
      if ( ! have_signal )
	Helper::halt( "ACTIG investigate requires sig=<channel>" );

      // Collect valid epoch values
      std::vector<double> valid_vals;
      valid_vals.reserve( n_valid_epochs );
      for ( int e = 0; e < n_epochs; e++ )
	if ( ! is_gap[e] ) valid_vals.push_back( epochs[e] );

      const int nv = (int) valid_vals.size();
      if ( nv == 0 ) { logger << "  INVESTIGATE: no valid epochs found\n"; return; }

      // Sort for percentile queries
      std::vector<double> sv = valid_vals;
      std::sort( sv.begin(), sv.end() );
      auto pct = [&]( double p ) -> double {
	const int idx = std::max( 0, std::min( (int)( p * nv ), nv - 1 ) );
	return sv[ idx ];
      };

      // Zero fraction
      const int n_zero = (int)std::count_if( valid_vals.begin(), valid_vals.end(),
					     []( double x ){ return x == 0.0; } );
      const double zero_frac = n_zero / (double) nv;

      // Overall mean / SD / CV
      const double g_mean = vec_mean( valid_vals );
      const double g_sd   = vec_sd( valid_vals );
      const double g_cv   = g_mean > 1e-9 ? g_sd / g_mean : 0.0;

      // log(1+x) space — mirrors luna scoring exactly
      std::vector<double> ly( nv );
      for ( int i = 0; i < nv; i++ ) ly[i] = std::log( 1.0 + valid_vals[i] );

      double ly_med = vec_median( ly );
      if ( ly_med < 1e-9 )
	{
	  // active-sub-distribution correction (same as luna method)
	  std::vector<double> act_ly;
	  for ( double v : ly ) if ( v > 1e-9 ) act_ly.push_back( v );
	  if ( ! act_ly.empty() ) ly_med = vec_median( act_ly );
	}
      const double ly_mad  = vec_mad( ly, ly_med );
      // z-score for a zero epoch: log(1+0) = 0
      const double z_zero  = ( ly_mad > 1e-9 ) ? ( 0.0 - ly_med ) / ly_mad
	                                         : std::numeric_limits<double>::quiet_NaN();
      const double active_frac_obs = ( nv - n_zero ) / (double) nv;

      // Zero-run and nonzero-run analysis
      std::vector<int> zero_runs, nonzero_runs;
      {
	int cur_z = 0, cur_nz = 0;
	for ( int e = 0; e < n_epochs; e++ )
	  {
	    if ( is_gap[e] )
	      {
		if ( cur_z  > 0 ) { zero_runs.push_back( cur_z );      cur_z  = 0; }
		if ( cur_nz > 0 ) { nonzero_runs.push_back( cur_nz );  cur_nz = 0; }
		continue;
	      }
	    if ( epochs[e] == 0.0 )
	      { if ( cur_nz > 0 ) { nonzero_runs.push_back( cur_nz ); cur_nz = 0; } cur_z++; }
	    else
	      { if ( cur_z  > 0 ) { zero_runs.push_back( cur_z );    cur_z  = 0; } cur_nz++; }
	  }
	if ( cur_z  > 0 ) zero_runs.push_back( cur_z );
	if ( cur_nz > 0 ) nonzero_runs.push_back( cur_nz );
      }

      auto run_med_min = [&]( const std::vector<int> & r ) -> double {
	if ( r.empty() ) return 0.0;
	std::vector<double> rd( r.begin(), r.end() );
	return vec_median( rd ) * epoch_sec / 60.0;
      };
      auto run_max_min = [&]( const std::vector<int> & r ) -> double {
	if ( r.empty() ) return 0.0;
	return *std::max_element( r.begin(), r.end() ) * epoch_sec / 60.0;
      };

      // Flat-frac: fraction of consecutive same-valued valid epoch pairs
      int flat_pair_n = 0, flat_pair_match_n = 0;
      {
	double prev = std::numeric_limits<double>::quiet_NaN();
	for ( int e = 0; e < n_epochs; e++ )
	  {
	    if ( is_gap[e] ) { prev = std::numeric_limits<double>::quiet_NaN(); continue; }
	    if ( std::isfinite( prev ) )
	      {
		flat_pair_n++;
		if ( std::fabs( epochs[e] - prev ) <= qc_flat_delta_th ) flat_pair_match_n++;
	      }
	    prev = epochs[e];
	  }
      }
      const double global_flat_frac = flat_pair_n > 0 ? flat_pair_match_n / (double) flat_pair_n : 0.0;

      // Nearfloor: fraction at/below P(qc_nearfloor_q_th) + eps
      const double nf_q_val = pct( qc_nearfloor_q_th );
      const double nf_eps   = std::max( std::fabs( pct( 0.95 ) - nf_q_val ) * 1e-6 , 1e-9 );
      const double nf_th    = nf_q_val + nf_eps;
      int n_nearfloor = 0;
      for ( double x : valid_vals ) if ( x <= nf_th ) ++n_nearfloor;
      const double global_nearfloor_frac = n_nearfloor / (double) nv;

      // ---- Derive report items ----

      const bool warn_flat      = global_flat_frac      >= qc_flat_frac_th;
      const bool warn_nearfloor = global_nearfloor_frac >= qc_nearfloor_frac_th;
      const bool warn_lowvar    = g_cv                  <  qc_lowvar_cv_th;
      const int  n_qc_triggers  = (int)warn_flat + (int)warn_nearfloor + (int)warn_lowvar;
      const bool qc_exclusion_likely = n_qc_triggers >= 2;

      const double max_z_min = run_max_min( zero_runs );

      // z-score check: will zero epochs qualify as sleep?
      const bool zeros_qualify = std::isfinite( z_zero ) && ( z_zero < luna_quiet_z );

      // Non-zero run stats (for score-period min-wake suggestion)
      const double nz_med_min  = run_med_min( nonzero_runs );
      const double nz_max_min  = run_max_min( nonzero_runs );

      // Zero-run P75 (for score-period min-sleep suggestion)
      double zero_run_p75_min = 0.0;
      if ( ! zero_runs.empty() )
	{
	  std::vector<double> zrd( zero_runs.begin(), zero_runs.end() );
	  std::sort( zrd.begin(), zrd.end() );
	  const int idx = std::min( (int)( 0.75 * zrd.size() ), (int)zrd.size() - 1 );
	  zero_run_p75_min = zrd[idx] * epoch_sec / 60.0;
	}

      // score-period: suggest min-sleep = 2 x P75 of zero runs, floored at 20 min
      // idea: the minimum sleep period should be longer than most short zero bursts
      const double sp_min_sleep_suggest = std::max( 20.0, 2.0 * zero_run_p75_min );
      // round to nearest 5 min
      const double sp_min_sleep = std::round( sp_min_sleep_suggest / 5.0 ) * 5.0;

      // score-period: suggest min-wake = median nonzero run, floored at 10 min
      const double sp_min_wake_suggest = std::max( 10.0, nz_med_min );
      const double sp_min_wake = std::round( sp_min_wake_suggest / 5.0 ) * 5.0;

      // Common optional flags (quiet-z, active-frac, QC)
      std::string common_flags;
      if ( ! zeros_qualify && std::isfinite( z_zero ) )
	common_flags += " quiet-z=" + Helper::dbl2str( std::floor( z_zero ) - 1.0 , 1 );
      if ( active_frac_obs < luna_active_frac && active_frac_obs > 0 )
	common_flags += " active-frac=" + Helper::dbl2str( active_frac_obs * 0.5 , 2 );
      if ( qc_exclusion_likely )
	{
	  if ( warn_flat )      common_flags += " qc-exclude-flat=F";
	  if ( warn_nearfloor ) common_flags += " qc-exclude-nearfloor=F";
	  if ( warn_lowvar )    common_flags += " qc-exclude-lowvar=F";
	}

      const std::string base = "  ACTIG sig=" + siglab
	+ ( Fs < 1.0 && expected_spe > 1 ? " sum" : "" );
      const std::string cmd_score    = base + " score" + common_flags;
      const std::string cmd_period   = base + " score-period" + common_flags
	+ " min-sleep=" + Helper::dbl2str( sp_min_sleep , 0 )
	+ " min-wake="  + Helper::dbl2str( sp_min_wake  , 0 );

      // ---- Print report ----
      logger << "\n";
      logger << "  ============================================================\n";
      logger << "  ACTIG INVESTIGATE: " << siglab << "\n";
      logger << "  ============================================================\n";
      logger << "\n";
      logger << "  Fs=" << Fs << " Hz, epoch=" << epoch_sec << " s ("
	     << (int)( Fs * epoch_sec + 0.5 ) << " sample/epoch), "
	     << nv << " valid epochs, " << n_gap_epochs << " gap\n";
      logger << "\n";

      // Signal distribution (compact)
      logger << "  Distribution  zeros=" << Helper::dbl2str( zero_frac * 100.0 , 1 ) << "%"
	     << "  P25=" << Helper::dbl2str( pct(0.25) , 2 )
	     << "  P50=" << Helper::dbl2str( pct(0.50) , 2 )
	     << "  P75=" << Helper::dbl2str( pct(0.75) , 2 )
	     << "  P95=" << Helper::dbl2str( pct(0.95) , 2 )
	     << "  max=" << Helper::dbl2str( pct(1.00) , 2 ) << "\n";
      logger << "  Mean=" << Helper::dbl2str( g_mean , 2 )
	     << "  SD=" << Helper::dbl2str( g_sd , 2 )
	     << "  CV=" << Helper::dbl2str( g_cv , 3 ) << "\n";
      logger << "\n";

      // Zero-run summary
      if ( ! zero_runs.empty() )
	logger << "  Zero runs: " << zero_runs.size()
	       << " runs, median=" << Helper::dbl2str( run_med_min( zero_runs ) , 1 ) << " min"
	       << ", max=" << Helper::dbl2str( max_z_min , 1 ) << " min\n";
      else
	logger << "  Zero runs: none\n";
      logger << "\n";

      // Issues and actions
      int issue_n = 0;
      logger << "  --- Issues ---\n";

      // 1. Zero-epoch scoring
      if ( std::isfinite( z_zero ) )
	{
	  if ( zeros_qualify )
	    {
	      logger << "  [OK] Zero-valued epochs score as sleep candidates\n"
		     << "       (z=" << Helper::dbl2str( z_zero , 2 )
		     << " < quiet-z=" << luna_quiet_z << "; no change needed)\n";
	    }
	  else
	    {
	      ++issue_n;
	      logger << "  [" << issue_n << "] Zero epochs score as WAKE, not sleep\n"
		     << "      z(zero)=" << Helper::dbl2str( z_zero , 2 )
		     << " >= quiet-z=" << luna_quiet_z << "\n"
		     << "      Fix: add  quiet-z=" << Helper::dbl2str( std::floor( z_zero ) - 1.0 , 1 )
		     << "  or use  method=threshold thresh=1\n";
	    }
	}

      // 2. Non-wear / long zero runs
      if ( max_z_min >= 240.0 )
	{
	  ++issue_n;
	  logger << "  [" << issue_n << "] Very long zero runs (max " << Helper::dbl2str( max_z_min , 0 )
		 << " min = " << Helper::dbl2str( max_z_min / 60.0 , 1 ) << " hrs)\n"
		 << "      Zero epochs score as sleep — but extended zeros are often non-wear,\n"
		 << "      not genuine sleep. There is no reliable way to distinguish the two\n"
		 << "      from activity counts alone. Options:\n"
		 << "        a) Accept zeros as sleep (device records during sleep with no movement)\n"
		 << "        b) Annotate confirmed non-wear periods before scoring and use MASK\n"
		 << "        c) Use method=threshold thresh=1 for a simple zero=sleep rule,\n"
		 << "           then inspect scored output manually\n";
	}
      else if ( max_z_min >= 60.0 )
	{
	  ++issue_n;
	  logger << "  [" << issue_n << "] Long zero runs present (max " << Helper::dbl2str( max_z_min , 0 )
		 << " min)\n"
		 << "      May be non-wear or genuine sleep; inspect scored output\n";
	}

      // 3. Low sampling rate
      if ( Fs < 1.0 && expected_spe > 1 )
	{
	  ++issue_n;
	  logger << "  [" << issue_n << "] Sub-Hz sample rate (" << Fs << " Hz, "
		 << expected_spe << " samples/epoch)\n"
		 << "      Add  sum  so epoch values are totals, not per-sample means\n";
	}

      // 4. active-frac mismatch
      if ( active_frac_obs < luna_active_frac && active_frac_obs > 0 )
	{
	  ++issue_n;
	  logger << "  [" << issue_n << "] Active (nonzero) fraction " << Helper::dbl2str( active_frac_obs * 100.0 , 1 )
		 << "% < active-frac=" << luna_active_frac << "\n"
		 << "      The burst-density gate uses active-frac to cap the fraction of\n"
		 << "      high-activity epochs in a local window before marking as wake.\n"
		 << "      With so few nonzero epochs this threshold may block sleep scoring.\n"
		 << "      Suggest: active-frac=" << Helper::dbl2str( active_frac_obs * 0.5 , 2 ) << "\n";
	}

      // 5. QC exclusions
      if ( qc_exclusion_likely )
	{
	  ++issue_n;
	  logger << "  [" << issue_n << "] Day QC flags likely to exclude whole days\n"
		 << "      Global estimates:";
	  if ( warn_flat )
	    logger << "  flat=" << Helper::dbl2str( global_flat_frac * 100.0 , 1 ) << "%";
	  if ( warn_nearfloor )
	    logger << "  nearfloor=" << Helper::dbl2str( global_nearfloor_frac * 100.0 , 1 ) << "%";
	  if ( warn_lowvar )
	    logger << "  CV=" << Helper::dbl2str( g_cv , 3 );
	  logger << "\n"
		 << "      Excluded days appear as all-WAKE in output, not as gaps.\n"
		 << "      Add these params to disable the offending flags:\n";
	  if ( warn_flat )      logger << "        qc-exclude-flat=F\n";
	  if ( warn_nearfloor ) logger << "        qc-exclude-nearfloor=F\n";
	  if ( warn_lowvar )    logger << "        qc-exclude-lowvar=F\n";
	}
      else if ( n_qc_triggers == 1 )
	{
	  logger << "  [OK] Day QC: signal looks mostly OK (1 marginal check, but 2 needed to drop a day)\n";
	}
      else
	{
	  logger << "  [OK] Day QC: signal looks fine, days should not be dropped by quality checks\n";
	}

      if ( issue_n == 0 )
	logger << "  No issues detected\n";

      logger << "\n";

      // Suggested commands
      logger << "  --- Suggested commands ---\n";
      logger << "\n";
      logger << "  Epoch-level W/S scoring:\n";
      logger << cmd_score << "\n";
      logger << "\n";
      logger << "  Sleep/wake period windows (SP/WP):\n";
      logger << cmd_period << "\n";
      if ( ! zero_runs.empty() )
	logger << "  (min-sleep=" << Helper::dbl2str( sp_min_sleep , 0 )
	       << " from 2 x P75 zero-run=" << Helper::dbl2str( zero_run_p75_min , 1 ) << " min"
	       << ";  min-wake=" << Helper::dbl2str( sp_min_wake , 0 )
	       << " from median nonzero-run=" << Helper::dbl2str( nz_med_min , 1 ) << " min)\n";
      if ( max_z_min >= 240.0 )
	{
	  logger << "\n";
	  logger << "  Simple zero=sleep threshold alternative:\n";
	  logger << base << " score method=threshold thresh=1\n";
	}
      logger << "\n";
      logger << "  ============================================================\n";
      logger << "\n";
      return;
    }


  // ---------------------------------------------------------------
  //
  // Nonparametric circadian metrics
  //
  // ---------------------------------------------------------------

  if ( have_signal && ! do_score_period )
    {

  //
  // Aggregate 1-min epochs into NP bins (default 60 min)
  // A NP bin is valid if >= gap-min-pct of its constituent epochs are valid.
  //

  const int epochs_per_bin = (int)( np_bin_min * 60.0 / epoch_sec );

  if ( epochs_per_bin <= 0 )
    Helper::halt( "ACTIG: NP bin size must be >= epoch length" );

  const int bins_per_day = (int)( 1440.0 / np_bin_min );
  if ( bins_per_day <= 0 )
    Helper::halt( "ACTIG: bin must be <= 1440 minutes" );

  const double np_first_boundary_sec = first_day_boundary_sec( startdatetime , day_anchor_hr , total_sec );
  int np_epoch_start = 0;
  int np_epoch_stop  = n_epochs;
  int np_day_count   = 0;

  if ( np_full_days && startdatetime.valid )
    {
      if ( total_sec > np_first_boundary_sec )
        np_day_count = (int)( ( total_sec - np_first_boundary_sec ) / 86400.0 );

      const double np_window_start_sec = np_first_boundary_sec;
      const double np_window_stop_sec  = np_first_boundary_sec + np_day_count * 86400.0;

      np_epoch_start = (int)std::ceil( np_window_start_sec / epoch_sec - 1e-9 );
      np_epoch_stop  = (int)std::floor( np_window_stop_sec / epoch_sec + 1e-9 );

      if ( np_epoch_start < 0 ) np_epoch_start = 0;
      if ( np_epoch_stop > n_epochs ) np_epoch_stop = n_epochs;
      if ( np_epoch_stop < np_epoch_start ) np_epoch_stop = np_epoch_start;
    }

  const int np_epochs_available = np_epoch_stop - np_epoch_start;
  const int n_bins = np_epochs_available / epochs_per_bin;

  if ( n_bins <= 0 )
    {
      logger << "  no NP bins (recording too short for current bin settings)\n";
      return;
    }

  if ( n_bins < 24 )
    logger << "  warning: fewer than 24 NP bins (" << n_bins
	   << "), metrics may be unreliable\n";

  const int min_valid_epe = std::max( 1 , (int)( gap_min_pct / 100.0 * epochs_per_bin ) );

  std::vector<double> bins( n_bins , 0.0 );
  std::vector<bool>   is_gap_bin( n_bins , true );

  for ( int b = 0 ; b < n_bins ; b++ )
    {
      double sum   = 0.0;
      int valid_n  = 0;
      const int e0 = np_epoch_start + b * epochs_per_bin;
      const int e1 = e0 + epochs_per_bin;
      for ( int e = e0 ; e < e1 ; e++ )
	if ( !is_gap[e] ) { sum += epochs[e]; valid_n++; }

      if ( valid_n >= min_valid_epe )
	{
	  bins[b] = sum / valid_n;
	  is_gap_bin[b] = false;
	}
    }

  int n_valid_bins = 0;
  for ( int b = 0 ; b < n_bins ; b++ )
    if ( !is_gap_bin[b] ) n_valid_bins++;


  //
  // Grand mean and total variance (valid bins only)
  //

  const int n_full_days = n_bins / bins_per_day;
  const int np_dropped_bins = n_bins % bins_per_day;

  int np_n_bins = n_bins;
  if ( np_full_days )
    {
      np_n_bins = n_bins - np_dropped_bins;
      if ( np_n_bins == 0 )
        logger << "  NP full-days restriction left no complete anchored 24 h days; NP metrics will not be computed\n";
      else if ( do_verbose )
        {
          if ( startdatetime.valid )
            {
              const int dropped_leading_epochs = np_epoch_start;
              const int dropped_trailing_epochs = n_epochs - np_epoch_stop
                + np_dropped_bins * epochs_per_bin;
              if ( dropped_leading_epochs > 0 || dropped_trailing_epochs > 0 )
                logger << "  NP full-days restriction kept " << n_full_days
                       << " anchored 24 h day(s) using day-anchor=" << day_anchor_hr
                       << ", dropping " << dropped_leading_epochs << " leading epoch(s)"
                       << " and " << dropped_trailing_epochs << " trailing epoch(s)\n";
            }
          else if ( np_n_bins != n_bins )
            logger << "  NP full-days restriction dropped trailing partial day: "
                   << np_dropped_bins << " NP bins ("
                   << ( np_dropped_bins * epochs_per_bin ) << " epochs)\n";
        }
    }
  const int np_epochs_used = np_n_bins * epochs_per_bin;
  const int np_epoch_used_stop = np_epoch_start + np_epochs_used;

  int np_n_valid_epochs = 0;
  for ( int e = np_epoch_start ; e < np_epoch_used_stop ; e++ )
    if ( !is_gap[e] ) ++np_n_valid_epochs;

  double grand_mean = 0;
  int    gm_n = 0;
  for ( int i = 0 ; i < np_n_bins ; i++ )
    if ( !is_gap_bin[i] ) { grand_mean += bins[i]; gm_n++; }
  if ( gm_n > 0 ) grand_mean /= gm_n;

  double total_var = 0;
  for ( int i = 0 ; i < np_n_bins ; i++ )
    if ( !is_gap_bin[i] )
      { double diff = bins[i] - grand_mean; total_var += diff * diff; }


  //
  // IS (Interdaily Stability)
  //
  //  Requires at least 2 days of valid bin coverage.
  //  Slot means computed from valid bins only.
  //

  double IS = 0;

  int np_n_valid_bins = 0;
  for ( int b = 0 ; b < np_n_bins ; b++ )
    if ( !is_gap_bin[b] ) np_n_valid_bins++;

  if ( np_n_valid_bins >= 2 * bins_per_day && total_var > 0 )
    {
      std::vector<double> slot_mean( bins_per_day , 0 );
      std::vector<int>    slot_n( bins_per_day , 0 );

      for ( int i = 0 ; i < np_n_bins ; i++ )
	if ( !is_gap_bin[i] )
	  {
	    int s = i % bins_per_day;
	    slot_mean[s] += bins[i];
	    slot_n[s]++;
	  }

      double between_var = 0;
      for ( int h = 0 ; h < bins_per_day ; h++ )
	{
	  if ( slot_n[h] > 0 ) slot_mean[h] /= slot_n[h];
	  double diff = slot_mean[h] - grand_mean;
	  between_var += slot_n[h] * diff * diff;
	}

      // OLD IS = ( (double)n_valid_bins / bins_per_day ) * between_var / total_var;
      IS = between_var / total_var;
    }


  //
  // IV (Intradaily Variability)
  //
  //  Cross-gap bin pairs are excluded from the numerator sum.
  //  Denominator uses valid bin count and valid-only total variance.
  //

  double IV = 0;

  if ( np_n_valid_bins >= 2 && total_var > 0 )
    {
      double diff_var = 0;
      int    diff_n   = 0;

      for ( int i = 1 ; i < np_n_bins ; i++ )
	if ( !is_gap_bin[i] && !is_gap_bin[i-1] )
	  {
	    double diff = bins[i] - bins[i-1];
	    diff_var += diff * diff;
	    diff_n++;
	  }

      if ( diff_n > 0 )
	IV = ( (double)np_n_valid_bins * diff_var ) / ( (double)diff_n * total_var );
    }


  //
  // L5 / M10 and RA
  //
  //  Use a 1440-point minute-of-day profile across the retained NP window,
  //  matching the nparACT approach of computing L5/M10 from minute averages
  //  rather than from the coarser hourly NP bins.
  //  Slots with no valid data remain 0.
  //

  double L5_val = 0 , M10_val = 0 , RA = 0;
  double L5_onset_hr = 0.0 , M10_onset_hr = 0.0;

  if ( np_n_bins > 0 )
    {
      const int minute_slots_per_day = 1440;
      std::vector<double> profile( minute_slots_per_day , 0.0 );
      std::vector<double> profile_w( minute_slots_per_day , 0.0 );

      for ( int e = np_epoch_start ; e < np_epoch_used_stop ; e++ )
	if ( !is_gap[e] )
	  {
	    const double start_min = ( ( e - np_epoch_start ) * epoch_sec ) / 60.0;
	    const double end_min   = ( ( e - np_epoch_start + 1 ) * epoch_sec ) / 60.0;

	    int m0 = (int)std::floor( start_min );
	    int m1 = (int)std::ceil( end_min );

	    for ( int m = m0 ; m < m1 ; m++ )
	      {
		const double slot0 = m;
		const double slot1 = m + 1.0;
		const double overlap = std::max( 0.0 , std::min( end_min , slot1 ) - std::max( start_min , slot0 ) );
		if ( overlap <= 0.0 ) continue;

		int md = m % minute_slots_per_day;
		if ( md < 0 ) md += minute_slots_per_day;

		profile[md]   += epochs[e] * overlap;
		profile_w[md] += overlap;
	      }
	  }
      int valid_profile_slots = 0;
      for ( int m = 0 ; m < minute_slots_per_day ; m++ )
	if ( profile_w[m] > 0.0 )
	  {
	    profile[m] /= profile_w[m];
	    valid_profile_slots++;
	  }

      if ( valid_profile_slots > 0 )
	{
	  find_window_fine( profile , 1.0 , l_hours , true  , &L5_val  , &L5_onset_hr  , np_step_min );
	  find_window_fine( profile , 1.0 , m_hours , false , &M10_val , &M10_onset_hr , np_step_min );

	  double denom = M10_val + L5_val;
	  RA = denom > 0 ? ( M10_val - L5_val ) / denom : 0;
	}
    }


  //
  // Clock-time labels and numeric onset values for L5/M10
  //
  // L5_onset_hr / M10_onset_hr are fractional hours from the profile
  // origin (= retained NP window start time of day).  Add that origin to
  // convert to a 24 h clock position, then extract HH:MM:SS and
  // decimal minutes-from-midnight for downstream analysis.
  //

  const double np_profile_origin_sec = startdatetime.valid
    ? std::fmod( startdatetime.h * 3600.0
                 + startdatetime.m * 60.0
                 + startdatetime.s
                 + np_epoch_start * epoch_sec , 86400.0 )
    : 0.0;
  const double np_profile_origin_hr = np_profile_origin_sec / 3600.0;

  // helper: fractional profile-hours → HH:MM:SS string and minutes-from-midnight
  auto onset_to_hms_and_min = [&]( double onset_hr ,
				   std::string * hms_str ,
				   double      * min_from_midnight )
  {
    double clock_hr = std::fmod( np_profile_origin_hr + onset_hr , 24.0 );
    if ( clock_hr < 0.0 ) clock_hr += 24.0;

    int hh = (int)clock_hr;
    int mm = (int)std::round( ( clock_hr - hh ) * 60.0 );
    if ( mm == 60 ) { hh = ( hh + 1 ) % 24; mm = 0; }

    *hms_str          = ( hh < 10 ? "0" : "" ) + Helper::int2str( hh ) + ":"
                      + ( mm < 10 ? "0" : "" ) + Helper::int2str( mm ) + ":00";
    *min_from_midnight = hh * 60.0 + mm;
  };

  std::string L5_hms , M10_hms;
  double L5_min_midnight = 0.0 , M10_min_midnight = 0.0;

  onset_to_hms_and_min( L5_onset_hr  , &L5_hms  , &L5_min_midnight  );
  onset_to_hms_and_min( M10_onset_hr , &M10_hms , &M10_min_midnight );


  //
  // Output NP metrics
  //

  if ( have_signal )
    {
      writer.value( "NP_IS" , IS );
      writer.value( "NP_IV" , IV );

      writer.value( "NP_L5"  , L5_val );
      writer.value( "NP_M10" , M10_val );
      writer.value( "NP_RA" , RA );

      writer.value( "NP_L5_ONSET"     , L5_hms );
      writer.value( "NP_M10_ONSET"    , M10_hms );
      writer.value( "NP_L5_ONSET_MIN" , L5_min_midnight );
      writer.value( "NP_M10_ONSET_MIN", M10_min_midnight );

      if ( do_verbose )
	{
	  logger << "  NP metrics: IS=" << IS
		 << " IV=" << IV
		 << " L" << l_hours << "=" << L5_val  << " (" << L5_hms  << ", " << L5_min_midnight  << " min)"
		 << " M" << m_hours << "=" << M10_val << " (" << M10_hms << ", " << M10_min_midnight << " min)"
		 << " RA=" << RA
		 << "\n";
	  if ( ! np_traditional )
	    logger << "  NP L5/M10 onset resolution: " << np_step_min << " min"
		   << " (use np-traditional for classic 1-hour resolution)\n";
	}
    }
  else if ( do_verbose )
    logger << "  NP metrics: no signal; IS/IV/L5/M10/RA not computed\n";

  writer.value( "NP_NE"    , np_full_days ? np_epochs_used : n_epochs );
  writer.value( "NP_NBINS" , np_n_bins );
  writer.value( "NP_NDAYS" , n_full_days );

  if ( have_signal )
    {
      if ( ! np_full_days )
	np_n_valid_epochs = n_valid_epochs;

      writer.value( "NP_NE_VALID"    , np_n_valid_epochs );
      writer.value( "NP_NBINS_VALID" , np_n_valid_bins );
      if ( do_verbose )
	logger << "  NP bins: " << np_n_bins << " total, "
	       << np_n_valid_bins << " valid\n";
    }

  } // if ( have_signal ) -- end NP section


  // ---------------------------------------------------------------
  //
  // Wake/sleep scoring
  //
  // ---------------------------------------------------------------

  //
  // Score each epoch as wake (true) or sleep (false), either from:
  //  - ACTIG scoring (if score requested), or
  //  - existing annotations (if prescored requested).
  // Gap epochs are never scored; is_wake[e] is irrelevant when is_gap[e].
  //

  std::vector<bool> is_wake( n_epochs , true );
  std::string score_method = method;

  // If no scoring mode and no explicit prescored, check whether epoch
  // annotations already exist in memory (e.g. from a prior ACTIG score run).
  // If so, treat as implicit prescored; otherwise return after NP metrics.
  bool do_prescored_eff = do_prescored;
  if ( ! do_any_score && ! do_prescored )
    {
      if ( do_debt )
	Helper::halt( "ACTIG: debt requires score or prescored to be specified" );
      if ( do_norm )
	Helper::halt( "ACTIG: norm requires score or prescored to be specified" );
      annot_t * impl_w = edf.timeline.annotations->find( presc_wake_label );
      annot_t * impl_s = edf.timeline.annotations->find( presc_sleep_label );
      if ( impl_w == NULL || impl_s == NULL )
	{
	  logger << "  no score mode and no epoch annotations ("
		 << presc_sleep_label << "/" << presc_wake_label
		 << ") found; skipping TST and fragmentation outputs\n";
	  return;
	}
      logger << "  no score mode but epoch annotations ("
	     << presc_sleep_label << "/" << presc_wake_label
	     << ") found; treating as prescored\n";
      do_prescored_eff = true;
    }

  if ( do_prescored_eff )
    {
      annot_t * a_wake_in  = edf.timeline.annotations->find( presc_wake_label );
      annot_t * a_sleep_in = edf.timeline.annotations->find( presc_sleep_label );
      if ( a_wake_in != NULL && a_sleep_in != NULL )
	{
	  int n_w = 0 , n_s = 0 , n_conflict = 0 , n_unset = 0;
	  for ( int i = 0 ; i < n_epochs ; i++ )
	    {
	      if ( is_gap[i] ) continue;
	      double start_sec = i * epoch_sec;
	      double stop_sec  = ( i + 1 ) * epoch_sec;
	      if ( stop_sec > total_sec ) stop_sec = total_sec;
	      interval_t ival( Helper::sec2tp( start_sec ) ,
			       Helper::sec2tp( stop_sec ) );
	      const bool has_w = ! a_wake_in->extract( ival ).empty();
	      const bool has_s = ! a_sleep_in->extract( ival ).empty();
	      if ( has_w && has_s )
		{
		  is_wake[i] = true;
		  ++n_conflict;
		  ++n_w;
		}
	      else if ( has_w )
		{
		  is_wake[i] = true;
		  ++n_w;
		}
	      else if ( has_s )
		{
		  is_wake[i] = false;
		  ++n_s;
		}
	      else
		{
		  is_wake[i] = true;
		  ++n_unset;
		}
	    }
	  score_method = "annots";
	  logger << "  prescored mode: using existing annotations sleep=" << presc_sleep_label
		 << ", wake=" << presc_wake_label
		 << " (wake-marked epochs=" << n_w
		 << ", sleep-marked epochs=" << n_s
		 << ", conflicts=" << n_conflict
		 << ", unmarked defaulted-to-wake=" << n_unset << ")\n";
	}
      else
	{
	  Helper::halt( "ACTIG: prescored mode requires both annotations: sleep="
			+ presc_sleep_label + ", wake=" + presc_wake_label );
	}
    }

  // Cole-Kripke weights are calibrated for 60-second epochs.
  else if ( ( method == "cole" || method == "threshold" ) &&
	    std::fabs( epoch_sec - 60.0 ) > 1e-6 )
    Helper::halt( "ACTIG: cole and threshold methods require epoch=60" );

  if ( do_any_score && ( method == "threshold" || method == "thresh" ) )
    {

      double t = thresh;

      // if no threshold given, use median of valid epochs only
      if ( t < 0 )
	{
	  std::vector<double> valid_vals;
	  valid_vals.reserve( n_valid_epochs );
	  for ( int i = 0 ; i < n_epochs ; i++ )
	    if ( !is_gap[i] ) valid_vals.push_back( epochs[i] );

	  if ( !valid_vals.empty() )
	    {
	      std::sort( valid_vals.begin() , valid_vals.end() );
	      t = valid_vals[ valid_vals.size() / 2 ];
	    }
	  else t = 0;

	  if ( do_verbose )
	    logger << "  threshold method: median of valid epochs = " << t << "\n";
	}

      // Optional triangular-window smoothing before threshold
      const int thresh_smooth_half = thresh_smooth_min > 0
	? std::max( 1, (int)std::round( thresh_smooth_min * 60.0 / epoch_sec ) ) : 0;

      if ( thresh_smooth_half > 0 )
	{
	  std::vector<double> smoothed( n_epochs , 0.0 );
	  for ( int i = 0 ; i < n_epochs ; i++ )
	    {
	      if ( is_gap[i] ) continue;
	      double wsum = 0.0, wtot = 0.0;
	      for ( int j = -thresh_smooth_half ; j <= thresh_smooth_half ; j++ )
		{
		  const int idx = i + j;
		  const double w = thresh_smooth_half + 1 - std::abs( j );
		  const double v = ( idx < 0 || idx >= n_epochs || is_gap[idx] ) ? 0.0 : epochs[idx];
		  wsum += w * v;
		  wtot += w;
		}
	      smoothed[i] = wtot > 0 ? wsum / wtot : 0.0;
	    }
	  for ( int i = 0 ; i < n_epochs ; i++ )
	    if ( !is_gap[i] ) is_wake[i] = ( smoothed[i] >= t );
	}
      else
	{
	  for ( int i = 0 ; i < n_epochs ; i++ )
	    if ( !is_gap[i] ) is_wake[i] = ( epochs[i] >= t );
	}

    }
  else if ( do_any_score && method == "cole" )
    {

      //
      // Cole-Kripke (1992)
      //   Score = 0.001 * sum( w[j] * A[n+offset[j]] )
      //   Gap neighbors are zero-substituted (conservative: biases toward
      //   sleep, avoids spurious wake calls at gap edges).
      //

      const double scale    = 0.001;
      const int n_weights   = 7;
      const int offsets[7]  = { -4, -3, -2, -1, 0, 1, 2 };
      const double weights[7] = { 404, 598, 326, 441, 1408, 508, 350 };

      const double wake_thresh = param.has( "cole-thresh" )
	? param.requires_dbl( "cole-thresh" ) : 1.0;

      for ( int i = 0 ; i < n_epochs ; i++ )
	{
	  if ( is_gap[i] ) continue;   // leave gap epochs unscored

	  double score = 0;
	  for ( int w = 0 ; w < n_weights ; w++ )
	    {
	      int idx = i + offsets[w];
	      if ( idx < 0 ) idx = 0;
	      if ( idx >= n_epochs ) idx = n_epochs - 1;
	      // zero-substitute gap neighbors
	      double val = is_gap[idx] ? 0.0 : epochs[idx];
	      score += weights[w] * val;
	    }
	  score *= scale;
	  is_wake[i] = ( score >= wake_thresh );
	}

    }
  else if ( do_any_score && method == "luna" )
    {

      // ---------------------------------------------------------------
      // Luna generic actigraphy sleep/wake detector
      //
      // Device-agnostic, epoch-size-agnostic algorithm:
      //   1. Partition observed epochs into contiguous blocks.
      //   2. Per-block log transform + robust (median/MAD) normalisation.
      //   3. Centred moving-median smoothing → L[i].
      //   4. Local burst-density estimate → D[i].
      //   5. Candidate sleep: L[i] < quiet_threshold_z AND D[i] < max_active_fraction.
      //   6. Run-length persistence: min sleep duration, gap filling,
      //      wake-break requirement.
      // All steps are confined to contiguous blocks; gaps are never bridged.
      // ---------------------------------------------------------------

      const int n_smooth     = std::max( 1, (int)std::round( luna_smooth_min    * 60.0 / epoch_sec ) );
      const int n_burst_w    = std::max( 1, (int)std::round( luna_burst_min     * 60.0 / epoch_sec ) );
      const int n_min_sleep  = std::max( 1, (int)std::round( luna_min_sleep_min * 60.0 / epoch_sec ) );
      const int n_gap_fill   = std::max( 0, (int)std::round( luna_max_gap_min   * 60.0 / epoch_sec ) );
      const int n_wake_break = std::max( 1, (int)std::round( luna_min_wake_min  * 60.0 / epoch_sec ) );

      if ( do_verbose )
	{
	  logger << "  luna: smooth=" << luna_smooth_min << "m (n=" << n_smooth << ")"
		 << ", burst=" << luna_burst_min << "m (n=" << n_burst_w << ")"
		 << ", burst-z=" << luna_burst_z
		 << ", quiet-z=" << luna_quiet_z
		 << ", active-frac=" << luna_active_frac << "\n";
	  logger << "  luna: min-sleep=" << luna_min_sleep_min << "m (n=" << n_min_sleep << ")"
		 << ", max-gap=" << luna_max_gap_min << "m (n=" << n_gap_fill << ")"
		 << ", min-wake=" << luna_min_wake_min << "m (n=" << n_wake_break << ")\n";
	}

      // Epoch-level diagnostic arrays (allocated only when channels requested)
      std::vector<double> eg_Z, eg_L, eg_D, eg_CS;
      if ( do_channels )
	{
	  eg_Z.assign ( n_epochs, 0.0 );
	  eg_L.assign ( n_epochs, 0.0 );
	  eg_D.assign ( n_epochs, 0.0 );
	  eg_CS.assign( n_epochs, 0.0 );
	}

      // Step 1: Identify contiguous blocks (maximal runs of non-gap epochs)
      std::vector<std::pair<int,int>> blocks;
      {
	int i = 0;
	while ( i < n_epochs )
	  {
	    if ( is_gap[i] ) { i++; continue; }
	    int j = i;
	    while ( j < n_epochs && !is_gap[j] ) j++;
	    blocks.push_back( { i, j - 1 } );
	    i = j;
	  }
      }

      if ( do_verbose )
	logger << "  luna: " << blocks.size() << " contiguous block(s)\n";

      for ( auto & blk : blocks )
	{
	  const int bs  = blk.first;
	  const int be  = blk.second;
	  const int bsz = be - bs + 1;

	  // Step 2: Log transform
	  std::vector<double> y( bsz );
	  for ( int i = 0; i < bsz; i++ )
	    y[i] = std::log( 1.0 + epochs[ bs + i ] );

	  // Step 3: Robust normalisation (median / MAD)
	  //
	  // When ≥50% of epochs have zero activity (non-wear or heavily
	  // inactive periods dominate the block), the global median collapses
	  // to zero.  Because log(1+x) ≥ 0 always, every z-score would then
	  // be ≥ 0 and no epoch could ever pass the quiet-z < 0 criterion.
	  // Fix: if the global median is effectively zero, centre instead on
	  // the median of strictly positive log-values (the "active sub-
	  // distribution").  This automatically equals the P(f_zero × 100)
	  // percentile of all data — the higher the zero fraction f_zero,
	  // the higher the reference percentile — so quiet/sleep epochs
	  // always receive negative z-scores regardless of zero-inflation.
	  double med = vec_median( y );
	  if ( med < 1e-9 )
	    {
	      std::vector<double> act_y;
	      act_y.reserve( bsz );
	      for ( double v : y ) if ( v > 1e-9 ) act_y.push_back( v );
	      if ( ! act_y.empty() ) med = vec_median( act_y );
	    }

	  std::vector<double> abs_dev( bsz );
	  for ( int i = 0; i < bsz; i++ )
	    abs_dev[i] = std::fabs( y[i] - med );

	  double mad = vec_median( abs_dev );

	  if ( mad < 1e-9 )
	    {
	      // Fallback 1: IQR-based scale
	      std::vector<double> sy = y;
	      std::sort( sy.begin(), sy.end() );
	      const double q1  = sy[ bsz / 4 ];
	      const double q3  = sy[ std::min( bsz - 1, 3 * bsz / 4 ) ];
	      const double iqr = q3 - q1;
	      if ( iqr > 1e-9 )
		mad = iqr / 1.3489795;   // IQR/1.349 ≈ σ for a normal distribution
	      else
		{
		  // Fallback 2: standard deviation
		  double mn = 0;
		  for ( double v : y ) mn += v;
		  mn /= bsz;
		  double var = 0;
		  for ( double v : y ) { double d = v - mn; var += d * d; }
		  const double sd = bsz > 1 ? std::sqrt( var / ( bsz - 1 ) ) : 0.0;
		  mad = sd > 1e-9 ? sd : 1e-9;
		}
	    }

	  std::vector<double> z( bsz );
	  for ( int i = 0; i < bsz; i++ )
	    z[i] = ( y[i] - med ) / mad;

	  // Step 4: Centred moving-median smoothing → L
	  std::vector<double> L( bsz );
	  {
	    const int half = n_smooth / 2;
	    for ( int i = 0; i < bsz; i++ )
	      {
		const int lo = std::max( 0, i - half );
		const int hi = std::min( bsz - 1, i + half );
		std::vector<double> win( z.begin() + lo, z.begin() + hi + 1 );
		L[i] = vec_median( win );
	      }
	  }

	  // Step 5: Local burst density → D
	  //   D[i] = fraction of epochs in window where z[j] > burst_threshold_z
	  std::vector<double> D( bsz );
	  {
	    const int half = n_burst_w / 2;
	    for ( int i = 0; i < bsz; i++ )
	      {
		const int lo    = std::max( 0, i - half );
		const int hi    = std::min( bsz - 1, i + half );
		const int n_win = hi - lo + 1;
		int n_bursts = 0;
		for ( int j = lo; j <= hi; j++ )
		  if ( z[j] > luna_burst_z ) n_bursts++;
		D[i] = (double)n_bursts / n_win;
	      }
	  }

	  // Step 6: Candidate sleep
	  std::vector<bool> cand( bsz );
	  for ( int i = 0; i < bsz; i++ )
	    cand[i] = ( L[i] < luna_quiet_z ) && ( D[i] < luna_active_frac );

	  // Diagnostic: summarise distribution and candidate counts
	  {
	    int n_cand_quiet = 0 , n_cand_burst = 0 , n_cand = 0;
	    double y_min = y[0] , y_max = y[0] , z_min = z[0] , z_max = z[0];
	    double D_mean = 0;
	    for ( int i = 0; i < bsz; i++ )
	      {
		y_min = std::min( y_min , y[i] ); y_max = std::max( y_max , y[i] );
		z_min = std::min( z_min , z[i] ); z_max = std::max( z_max , z[i] );
		D_mean += D[i];
		if ( L[i] < luna_quiet_z )  n_cand_quiet++;
		if ( D[i] < luna_active_frac ) n_cand_burst++;
		if ( cand[i] ) n_cand++;
	      }
	    D_mean /= bsz;
	    if ( do_verbose )
	      logger << "  luna diag block [" << bs << "," << be << "] bsz=" << bsz
		     << ": log(act) range=[" << y_min << "," << y_max << "]"
		     << " med=" << vec_median(y) << " z-range=[" << z_min << "," << z_max << "]"
		     << " D_mean=" << D_mean
		     << " n_quiet(L<" << luna_quiet_z << ")=" << n_cand_quiet
		     << " n_burst_ok(D<" << luna_active_frac << ")=" << n_cand_burst
		     << " n_cand=" << n_cand << "\n";
	  }

	  // Step 7: Run-length persistence rules

	  // Pass 1: Fill short non-candidate gaps sandwiched inside sleep
	  std::vector<bool> filled = cand;
	  {
	    int i = 0;
	    while ( i < bsz )
	      {
		if ( !filled[i] )
		  {
		    int j = i;
		    while ( j < bsz && !filled[j] ) j++;
		    const int rlen = j - i;
		    const bool before = ( i > 0   && filled[i - 1] );
		    const bool after  = ( j < bsz && filled[j]     );
		    if ( rlen <= n_gap_fill && before && after )
		      for ( int k = i; k < j; k++ ) filled[k] = true;
		    i = j;
		  }
		else i++;
	      }
	  }

	  // Post pass-1 count
	  {
	    int n1 = 0; for ( int i=0;i<bsz;i++) if(filled[i]) n1++;
	    if ( do_verbose )
	      logger << "  luna diag: after gap-fill (pass1): " << n1 << " sleep epochs\n";
	  }

	  // Pass 2: Remove sleep runs shorter than n_min_sleep
	  {
	    int i = 0;
	    while ( i < bsz )
	      {
		if ( filled[i] )
		  {
		    int j = i;
		    while ( j < bsz && filled[j] ) j++;
		    const int rlen = j - i;
		    if ( rlen < n_min_sleep )
		      for ( int k = i; k < j; k++ ) filled[k] = false;
		    i = j;
		  }
		else i++;
	      }
	  }

	  // Post pass-2 count
	  {
	    int n2 = 0; for ( int i=0;i<bsz;i++) if(filled[i]) n2++;
	    if ( do_verbose )
	      logger << "  luna diag: after min-sleep (pass2): " << n2 << " sleep epochs\n";
	  }

	  // Pass 3: Wake-persistence state machine
	  //   Once in sleep, stay asleep until n_wake_break consecutive
	  //   non-sleep epochs accumulate.  Sleep always terminates at the
	  //   block boundary.  Any pending "short wake" that never reaches
	  //   n_wake_break before a new sleep epoch is retracted (filled).
	  std::vector<bool> final_sleep( bsz, false );
	  {
	    bool in_sleep       = false;
	    int  wake_run_start = -1;

	    for ( int i = 0; i <= bsz; i++ )   // i == bsz is end-of-block sentinel
	      {
		const bool cur    = ( i < bsz ) && filled[i];
		const bool at_end = ( i >= bsz );

		if ( cur )
		  {
		    // Sleep epoch: enter / remain in sleep, cancel pending wake.
		    in_sleep       = true;
		    wake_run_start = -1;
		    final_sleep[i] = true;
		  }
		else if ( in_sleep )
		  {
		    // Non-sleep epoch while in sleep state.
		    if ( wake_run_start < 0 )
		      wake_run_start = i;

		    const int run_end = ( i < bsz ) ? i : bsz - 1;
		    const int rlen    = run_end - wake_run_start + 1;

		    if ( rlen >= n_wake_break || at_end )
		      {
			// Confirmed wake (or block ended): retroactively
			// un-fill any epochs we had temporarily marked sleep.
			for ( int k = wake_run_start; k < i && k < bsz; k++ )
			  final_sleep[k] = false;
			in_sleep       = false;
			wake_run_start = -1;
		      }
		    else if ( i < bsz )
		      {
			// Short non-sleep: still in sleep (tentatively filled).
			final_sleep[i] = true;
		      }
		  }
		// !in_sleep and !cur: final_sleep[i] stays false.
	      }
	  }

	  // Accumulate epoch-level diagnostic values
	  if ( do_channels )
	    for ( int i = 0; i < bsz; i++ )
	      {
		eg_Z [ bs + i ] = z[i];
		eg_L [ bs + i ] = L[i];
		eg_D [ bs + i ] = D[i];
		eg_CS[ bs + i ] = cand[i] ? 1.0 : 0.0;
	      }

	  // Map back to global is_wake[]
	  for ( int i = 0; i < bsz; i++ )
	    is_wake[ bs + i ] = !final_sleep[i];

	}  // end per-block loop

      // Add diagnostic channels to the EDF if requested
      if ( do_channels )
	{
	  // Build epoch-level SW array from is_wake[]: sleep=1, wake=0.
	  // Gap epochs are never observed (EDF+D has no samples there),
	  // so their value is irrelevant; leave at 0.
	  std::vector<double> eg_SW( n_epochs, 0.0 );
	  for ( int e = 0; e < n_epochs; e++ )
	    if ( !is_gap[e] && !is_wake[e] ) eg_SW[e] = 1.0;

	  // Helper: expand epoch-level values to sample-level using
	  // the same timestamp-to-bin mapping used during binning.
	  auto expand = [&]( const std::vector<double> & epoch_vals ) -> std::vector<double>
	  {
	    std::vector<double> sig( d->size(), 0.0 );
	    for ( int i = 0; i < (int)d->size(); i++ )
	      {
		if ( (*tp)[i] < recording_start_tp ) continue;
		const uint64_t offset_tp = (*tp)[i] - recording_start_tp;
		const int ebin = (int)( offset_tp / epoch_tp );
		if ( ebin >= 0 && ebin < n_epochs )
		  sig[i] = epoch_vals[ ebin ];
	      }
	    return sig;
	  };

	  // For sub-1 Hz signals, Fs * record_duration rounds to 0 in add_signal().
	  // Use the negative-Fs backdoor: pass -(n_samples_per_record) directly.
	  const int    n_spr   = edf.header.n_samples[ slot ];
	  const double ch_Fs   = ( Fs >= 1.0 ) ? Fs : -(double)n_spr;

	  edf.add_signal( siglab + "_Z"  , ch_Fs , expand( eg_Z  ) );
	  edf.add_signal( siglab + "_L"  , ch_Fs , expand( eg_L  ) );
	  edf.add_signal( siglab + "_D"  , ch_Fs , expand( eg_D  ) );
	  edf.add_signal( siglab + "_CS" , ch_Fs , expand( eg_CS ) );
	  edf.add_signal( siglab + "_SW" , ch_Fs , expand( eg_SW ) );

	  logger << "  luna: added channels: "
		 << siglab << "_Z (normalised), "
		 << siglab << "_L (smoothed), "
		 << siglab << "_D (burst density), "
		 << siglab << "_CS (candidate sleep), "
		 << siglab << "_SW (sleep=1, wake=0)\n";
	}

    }  // end luna method

  else if ( do_any_score )
    {
      Helper::halt( "ACTIG: unknown scoring method '" + method
		    + "' (use: luna, cole, threshold)" );
    }


  //
  // Create WAKE, SLEEP, and GAP annotations from contiguous runs
  //
  //  Epoch state: 0 = gap, 1 = wake, 2 = sleep
  //

  annot_t * a_wake  = NULL;
  annot_t * a_sleep = NULL;
  if ( do_any_score )
    {
      if ( edf.timeline.annotations->find( wake_label ) != NULL )
	Helper::halt( "ACTIG: score output label '" + wake_label
		      + "' already exists; drop or remap it before scoring" );
      if ( edf.timeline.annotations->find( sleep_label ) != NULL )
	Helper::halt( "ACTIG: score output label '" + sleep_label
		      + "' already exists; drop or remap it before scoring" );

      a_wake  = edf.timeline.annotations->add( wake_label );
      a_sleep = edf.timeline.annotations->add( sleep_label );

      a_wake->description  = "Actigraphy-scored wake";
      a_sleep->description = "Actigraphy-scored sleep";
      a_wake->file = a_sleep->file = "ACTIG";
    }

  auto epoch_state = [&]( int i ) -> int {
    if ( is_gap[i]  ) return 0;
    if ( is_wake[i] ) return 1;
    return 2;
  };

  int wake_inst = 1 , sleep_inst = 1;
  int total_wake = 0 , total_sleep = 0;

  int run_start = 0;
  for ( int i = 1 ; i <= n_epochs ; i++ )
    {
      if ( i == n_epochs || epoch_state(i) != epoch_state(run_start) )
	{
	  double start_sec = run_start * epoch_sec;
	  double stop_sec  = i * epoch_sec;
	  if ( stop_sec > total_sec ) stop_sec = total_sec;

	  interval_t ival( Helper::sec2tp( start_sec ) ,
			   Helper::sec2tp( stop_sec ) );
	  const int run_len = i - run_start;
	  const int st = epoch_state( run_start );

	  if ( st == 1 )
	    {
	      if ( do_any_score ) a_wake->add( Helper::int2str( wake_inst ) , ival , "." );
	      ++wake_inst;
	      total_wake += run_len;
	    }
	  else
	    {
	      if ( do_any_score ) a_sleep->add( Helper::int2str( sleep_inst ) , ival , "." );
	      ++sleep_inst;
	      total_sleep += run_len;
	    }

	  run_start = i;
	}
    }


  //
  // score-period: annotations written; skip all metric computation and output
  //

  if ( do_score_period )
    {
      logger << "  score-period: SP/WP annotations written; all metrics skipped\n";
      return;
    }


  //
  // Overall summary stats
  //  TST + TWT are computed over valid (non-gap) epochs only.
  //  sleep_pct denominator is valid epochs, not total epochs.
  //

  const double total_sleep_min     = total_sleep * epoch_sec / 60.0;
  const double total_wake_min      = total_wake  * epoch_sec / 60.0;
  const double sleep_pct           = n_valid_epochs > 0
    ? 100.0 * total_sleep / (double)n_valid_epochs : 0;

  // Period annotations (S/W): define sleep/wake windows for fragmentation.
  // Looked up regardless of scoring mode; fall back to empirical if absent.
  annot_t * a_sleep_period = edf.timeline.annotations->find( period_sleep_label );
  annot_t * a_wake_period  = edf.timeline.annotations->find( period_wake_label  );
  const bool have_sleep_period_annot = ( a_sleep_period != NULL );
  const bool have_wake_period_annot  = ( a_wake_period  != NULL );

  // Per-epoch flags: does this epoch overlap an S- or W-period annotation?
  std::vector<bool> in_sleep_period( n_epochs , false );
  std::vector<bool> in_wake_period ( n_epochs , false );
  for ( int i = 0 ; i < n_epochs ; i++ )
    {
      const double start_sec = i * epoch_sec;
      double stop_sec = ( i + 1 ) * epoch_sec;
      if ( stop_sec > total_sec ) stop_sec = total_sec;
      const interval_t ival( Helper::sec2tp( start_sec ) ,
			      Helper::sec2tp( stop_sec ) );
      if ( have_sleep_period_annot )
	in_sleep_period[i] = ! a_sleep_period->extract( ival ).empty();
      if ( have_wake_period_annot )
	in_wake_period[i]  = ! a_wake_period->extract(  ival ).empty();
    }
  if ( have_sleep_period_annot )
    logger << "  period annotations: using '" << period_sleep_label
	   << "' for sleep-window definition\n";
  if ( have_wake_period_annot )
    logger << "  period annotations: using '" << period_wake_label
	   << "' for wake-window definition\n";

  writer.value( "SCORE_TST"   , total_sleep_min );
  writer.value( "SCORE_TWT"   , total_wake_min );
  if ( total_sleep + total_wake > 0 )
    writer.value( "SCORE_SLEEP_PCT" , sleep_pct );



  // -----------------------------------------------------------------------
  //
  // Per-day summary
  //
  //  Epochs are clock-time-anchored (bin e corresponds to recording_start
  //  + e * epoch_sec), so day partitioning is straightforward.
  //  Gap epochs are tracked separately; they do not contribute to TST/TWT.
  //  A day is INCLUDED in cross-day averages only if valid minutes >= day-min-valid.
  //
  // -----------------------------------------------------------------------

  double first_boundary_sec = first_day_boundary_sec( startdatetime , day_anchor_hr , total_sec );

  auto accumulate_day_seconds = [&]( std::vector<double> * x ,
				     double start_sec , double stop_sec )
  {
    while ( start_sec < stop_sec )
      {
	int day_idx = 0;
	double day_stop = first_boundary_sec;

	if ( start_sec >= first_boundary_sec )
	  {
	    day_idx = 1 + (int)( ( start_sec - first_boundary_sec ) / 86400.0 );
	    day_stop = first_boundary_sec + day_idx * 86400.0;
	  }

	if ( day_stop > total_sec ) day_stop = total_sec;

	double seg_stop = stop_sec < day_stop ? stop_sec : day_stop;

	if ( (size_t)day_idx >= x->size() )
	  x->resize( day_idx + 1 , 0.0 );

	(*x)[ day_idx ] += seg_stop - start_sec;
	start_sec = seg_stop;
      }
  };

  auto day_index_for_sec = [&]( double sec ) -> int
  {
    if ( sec < first_boundary_sec ) return 0;
    return 1 + (int)( ( sec - first_boundary_sec ) / 86400.0 );
  };

  std::vector<double> sleep_sec_day;
  std::vector<double> wake_sec_day;
  std::vector<double> sleep_sec_sp_day;   // TST within SP window
  std::vector<double> wake_sec_sp_day;    // TWT within SP window
  std::vector<double> sleep_sec_wp_day;   // TST within WP window
  std::vector<double> wake_sec_wp_day;    // TWT within WP window

  std::vector<int> day_epoch_n;        // valid (non-gap) epochs per day
  std::vector<int> day_wake_epoch_n;
  std::vector<int> day_trans_to_w;
  std::vector<int> day_first_idx;      // first valid epoch index in day
  std::vector<int> day_last_idx;       // last valid epoch index in day

  auto ensure_day = [&]( int d )
  {
    if ( d < 0 ) return;
    if ( (size_t)(d+1) > day_epoch_n.size() )
      {
	day_epoch_n.resize( d+1 , 0 );
	day_wake_epoch_n.resize( d+1 , 0 );
	day_trans_to_w.resize( d+1 , 0 );
	day_first_idx.resize( d+1 , -1 );
	day_last_idx.resize( d+1 , -1 );
      }
  };

  // Main per-epoch accumulation loop
  for ( int i = 0 ; i < n_epochs ; i++ )
    {
      const double start_sec = i * epoch_sec;
      const int di = day_index_for_sec( start_sec );
      ensure_day( di );

      if ( is_gap[i] ) continue;

      const double stop_sec = std::min( (i+1) * epoch_sec , total_sec );

      day_epoch_n[di]++;
      if ( is_wake[i] ) day_wake_epoch_n[di]++;
      if ( day_first_idx[di] == -1 ) day_first_idx[di] = i;
      day_last_idx[di] = i;

      if ( is_wake[i] )
	accumulate_day_seconds( &wake_sec_day , start_sec , stop_sec );
      else
	accumulate_day_seconds( &sleep_sec_day , start_sec , stop_sec );

      // SP/WP-stratified TST/TWT
      if ( in_sleep_period[i] )
	{
	  if ( is_wake[i] ) accumulate_day_seconds( &wake_sec_sp_day  , start_sec , stop_sec );
	  else              accumulate_day_seconds( &sleep_sec_sp_day , start_sec , stop_sec );
	}
      else if ( in_wake_period[i] )
	{
	  if ( is_wake[i] ) accumulate_day_seconds( &wake_sec_wp_day  , start_sec , stop_sec );
	  else              accumulate_day_seconds( &sleep_sec_wp_day , start_sec , stop_sec );
	}
    }

  // Transitions (skip cross-gap)
  for ( int i = 1 ; i < n_epochs ; i++ )
    if ( !is_gap[i] && !is_gap[i-1] && !is_wake[i-1] && is_wake[i] )
      {
	const int di = day_index_for_sec( i * epoch_sec );
	ensure_day( di );
	day_trans_to_w[di]++;
      }

  const size_t n_days = std::max( { sleep_sec_day.size() ,
				    wake_sec_day.size() ,
				    day_epoch_n.size() } );

  double sum_sleep_sec = 0.0;
  double sum_wake_sec  = 0.0;
  double sum_sleep_sec_sp = 0.0 , sum_wake_sec_sp = 0.0;
  double sum_sleep_sec_wp = 0.0 , sum_wake_sec_wp = 0.0;
  int    n_included_days = 0;
  int    n_qc_ok_days = 0;
  int    n_qc_excluded_days = 0;
  int    n_qc_warn_days = 0;
  int    n_included_postqc_days = 0;

  // Accumulators for day-averaged fragmentation metrics
  double sum_frag_sfi           = 0.0;
  double sum_frag_mi_act        = 0.0;
  double sum_frag_fi_act        = 0.0;
  double sum_frag_sfi_act       = 0.0;
  double sum_frag_wake_bout_med  = 0.0;
  double sum_frag_wake_bout_p90  = 0.0;
  double sum_frag_wake_bout_max  = 0.0;
  double sum_frag_sleep_bout_med = 0.0;
  double sum_frag_sleep_bout_p10 = 0.0;
  double sum_frag_sleep_bout_max = 0.0;
  double sum_frag_trans_ent_S    = 0.0;
  double sum_frag_rl_ent_se_S    = 0.0;
  double sum_frag_rl_ent_we_S    = 0.0;
  double sum_frag_tp_sw_S        = 0.0;
  double sum_frag_tp_ws_S        = 0.0;
  double sum_frag_trans_ent_W    = 0.0;
  double sum_frag_rl_ent_se_W    = 0.0;
  double sum_frag_rl_ent_we_W    = 0.0;
  double sum_frag_tp_sw_W        = 0.0;
  double sum_frag_tp_ws_W        = 0.0;
  int    n_frag_dayavg           = 0;
  int    n_frag_W_dayavg         = 0;

  std::vector<actig_day_qc_t> day_qc( n_days );
  annot_t * a_qc_excluded = NULL;
  annot_t * a_qc_extreme  = NULL;
  writer.numeric_factor( "DAY" );
  if ( qc_day && qc_exclude_label != "" )
    {
      a_qc_excluded = edf.timeline.annotations->add( qc_exclude_label );
      a_qc_excluded->description = "ACTIG QC-excluded day (signal quality)";
      a_qc_excluded->file = "ACTIG";
    }
  if ( qc_day && qc_extreme_label != "" && ( qc_exclude_allsleep || qc_exclude_allwake ) )
    {
      a_qc_extreme = edf.timeline.annotations->add( qc_extreme_label );
      a_qc_extreme->description = "ACTIG QC-excluded day (extreme sleep/wake ratio)";
      a_qc_extreme->file = "ACTIG";
    }

  auto join_flags = []( const std::vector<std::string> & x ) -> std::string
  {
    std::string s;
    for ( size_t i = 0 ; i < x.size() ; i++ )
      {
	if ( i ) s += ",";
	s += x[i];
      }
    return s;
  };

  auto minutes_from_epochs = [&]( int run_epochs ) -> double
  {
    return run_epochs * epoch_sec / 60.0;
  };

  for ( size_t d = 0 ; d < n_days ; d++ )
    {
      const double sleep_sec    = d < sleep_sec_day.size()  ? sleep_sec_day[d]   : 0.0;
      const double wake_sec     = d < wake_sec_day.size()   ? wake_sec_day[d]    : 0.0;
      const double sleep_sec_sp = d < sleep_sec_sp_day.size() ? sleep_sec_sp_day[d] : 0.0;
      const double wake_sec_sp  = d < wake_sec_sp_day.size()  ? wake_sec_sp_day[d]  : 0.0;
      const double sleep_sec_wp = d < sleep_sec_wp_day.size() ? sleep_sec_wp_day[d] : 0.0;
      const double wake_sec_wp  = d < wake_sec_wp_day.size()  ? wake_sec_wp_day[d]  : 0.0;
      const int day_epochs      = d < day_epoch_n.size()    ? day_epoch_n[d]     : 0;
      const int day_wake_epochs = d < day_wake_epoch_n.size()? day_wake_epoch_n[d]: 0;

      int    day_trans_win  = 0;
      double day_sfi_luna_h = 0.0;

      const double valid_min_day   = day_epochs * epoch_sec / 60.0;
      const double day_total_valid = sleep_sec + wake_sec;
      const double day_sleep_pct   = day_total_valid > 0
	? 100.0 * sleep_sec / day_total_valid : 0.0;

      const bool included = ( valid_min_day >= day_min_valid_min );

      // Per-day fragmentation (gap-aware, same logic as overall)
      int    day_imm_bout_n = 0 , day_imm1_bout_n = 0;
      double day_mi_pct = 0 , day_imm1_pct = 0 , day_sfi_act = 0;
      bool   day_has_sleep_window = false;
      bool   day_has_wake_window  = false;
      int    day_frag_wake_bout_n = 0;
      double day_frag_wake_bout_med = 0 , day_frag_wake_bout_p90 = 0 , day_frag_wake_bout_max = 0;
      double day_frag_sleep_bout_med = 0 , day_frag_sleep_bout_p10 = 0 , day_frag_sleep_bout_max = 0;
      double day_trans_ent_S = 0 , day_rl_ent_se_S = 0 , day_rl_ent_we_S = 0;
      double day_trans_ent_W = 0 , day_rl_ent_se_W = 0 , day_rl_ent_we_W = 0;
      double day_tp_sw_S = 0 , day_tp_ws_S = 0;
      double day_tp_sw_W = 0 , day_tp_ws_W = 0;
      actig_day_qc_t qc;
      bool lowvar_frac_flag = false;

      if ( d < day_first_idx.size() && day_first_idx[d] != -1 )
	{
	  const int b = day_first_idx[d];
	  const int e = day_last_idx[d];
	  std::vector<double> day_epochs_valid;
	  day_epochs_valid.reserve( day_epochs );
	  int flat_pair_n = 0 , flat_pair_match_n = 0;
	  int prev_valid_idx = -1;

	  // QC data collection (valid epochs only)
	  for ( int i = b ; i <= e ; i++ )
	    if ( !is_gap[i] )
	      {
		day_epochs_valid.push_back( epochs[i] );
		if ( prev_valid_idx != -1 )
		  {
		    flat_pair_n++;
		    if ( std::fabs( epochs[i] - epochs[ prev_valid_idx ] ) <= qc_flat_delta_th )
		      flat_pair_match_n++;
		  }
		prev_valid_idx = i;
	      }

	  // Build in_S_day: which epochs within [b,e] are "in the sleep period".
	  // Source is S-period annotations if present; otherwise empirical
	  // (first..last non-gap sleep epoch within the day's range).
	  // Fragmentation code consumes this vector directly — no further branching.
	  std::vector<bool> in_S_day( n_epochs , false );
	  if ( have_sleep_period_annot )
	    {
	      for ( int i = b ; i <= e ; i++ )
		in_S_day[i] = in_sleep_period[i];
	    }
	  else
	    {
	      int fs = -1 , ls = -1;
	      for ( int i = b ; i <= e ; i++ )
		if ( !is_gap[i] && !is_wake[i] )
		  {
		    if ( fs == -1 ) fs = i;
		    ls = i;
		  }
	      if ( fs != -1 )
		for ( int i = fs ; i <= ls ; i++ )
		  in_S_day[i] = true;
	    }

	  // Sleep window bounds (for extract_bouts range and day_has_sleep_window)
	  int first_sleep_day = -1 , last_sleep_day = -1;
	  for ( int i = b ; i <= e ; i++ )
	    if ( in_S_day[i] )
	      {
		if ( first_sleep_day == -1 ) first_sleep_day = i;
		last_sleep_day = i;
	      }

	  if ( first_sleep_day != -1 )
	    {
	      // MI + SFI: only over in-S-period, non-gap epochs.
	      // Gaps and out-of-S-period epochs are skipped; prev_in_S tracks
	      // the last valid epoch so transitions are never counted across them.
	      int wake_in_win = 0 , valid_in_win = 0 , sleep_in_win = 0;
	      int prev_in_S = -1;
	      for ( int i = first_sleep_day ; i <= last_sleep_day ; i++ )
		{
		  if ( is_gap[i] || !in_S_day[i] ) continue;
		  valid_in_win++;
		  if ( is_wake[i] ) wake_in_win++;
		  else              sleep_in_win++;
		  if ( prev_in_S != -1 && !is_wake[prev_in_S] && is_wake[i] )
		    ++day_trans_win;
		  prev_in_S = i;
		}
	      day_mi_pct = valid_in_win > 0
		? 100.0 * wake_in_win / (double)valid_in_win : 0.0;
	      const double day_win_sleep_h = sleep_in_win * epoch_sec / 3600.0;
	      day_sfi_luna_h = day_win_sleep_h > 0 ? day_trans_win / day_win_sleep_h : 0.0;

	      day_has_sleep_window = true;
	      const bout_counts_t day_bouts = extract_bouts( is_gap , is_wake , in_S_day ,
                                                             first_sleep_day , last_sleep_day );
	      day_imm_bout_n  = (int)day_bouts.sleep_epochs.size();
	      day_imm1_bout_n = 0;
	      for ( int r : day_bouts.sleep_epochs ) if ( r <= 1 ) ++day_imm1_bout_n;

	      day_imm1_pct = day_imm_bout_n > 0
		? 100.0 * day_imm1_bout_n / (double)day_imm_bout_n : 0.0;
	      day_sfi_act = day_mi_pct + day_imm1_pct;

	      day_frag_wake_bout_n   = (int)day_bouts.wake_epochs.size();
	      day_frag_wake_bout_med = pct_min( day_bouts.wake_epochs  , 0.50 , epoch_sec );
	      day_frag_wake_bout_p90 = pct_min( day_bouts.wake_epochs  , 0.90 , epoch_sec );
	      day_frag_wake_bout_max = day_bouts.wake_epochs.empty()  ? 0.0
		: *std::max_element( day_bouts.wake_epochs.begin()  , day_bouts.wake_epochs.end() )
		  * epoch_sec / 60.0;
	      day_frag_sleep_bout_med = pct_min( day_bouts.sleep_epochs , 0.50 , epoch_sec );
	      day_frag_sleep_bout_p10 = pct_min( day_bouts.sleep_epochs , 0.10 , epoch_sec );
	      day_frag_sleep_bout_max = day_bouts.sleep_epochs.empty() ? 0.0
		: *std::max_element( day_bouts.sleep_epochs.begin() , day_bouts.sleep_epochs.end() )
		  * epoch_sec / 60.0;

	      // Entropy metrics for S period
	      day_trans_ent_S = transition_entropy( is_gap , is_wake , in_S_day ,
						    first_sleep_day , last_sleep_day );
	      day_rl_ent_se_S = rl_entropy( day_bouts.sleep_epochs );
	      day_rl_ent_we_S = rl_entropy( day_bouts.wake_epochs  );

	      // Transition probabilities for S period (ML + Bayesian λ correction)
	      const tp_pair_t tp_S = transition_probs( is_gap , is_wake , in_S_day ,
						       first_sleep_day , last_sleep_day , tp_lambda );
	      day_tp_sw_S = tp_S.tp_sw;
	      day_tp_ws_S = tp_S.tp_ws;
	    }

	  // Build in_W_day: W-period annotation if present, otherwise complement
	  // of S period (all non-gap, non-S-period epochs within the day).
	  std::vector<bool> in_W_day( n_epochs , false );
	  if ( have_wake_period_annot )
	    {
	      for ( int i = b ; i <= e ; i++ )
		in_W_day[i] = in_wake_period[i];
	    }
	  else
	    {
	      for ( int i = b ; i <= e ; i++ )
		if ( !is_gap[i] && !in_S_day[i] )
		  in_W_day[i] = true;
	    }

	  // W window bounds
	  int first_wake_day = -1 , last_wake_day = -1;
	  for ( int i = b ; i <= e ; i++ )
	    if ( in_W_day[i] )
	      {
		if ( first_wake_day == -1 ) first_wake_day = i;
		last_wake_day = i;
	      }

	  // Entropy metrics for W period
	  if ( first_wake_day != -1 )
	    {
	      const bout_counts_t W_bouts = extract_bouts( is_gap , is_wake , in_W_day ,
							   first_wake_day , last_wake_day );
	      day_trans_ent_W = transition_entropy( is_gap , is_wake , in_W_day ,
						    first_wake_day , last_wake_day );
	      day_rl_ent_se_W = rl_entropy( W_bouts.sleep_epochs );
	      day_rl_ent_we_W = rl_entropy( W_bouts.wake_epochs  );

	      // Transition probabilities for W period
	      const tp_pair_t tp_W = transition_probs( is_gap , is_wake , in_W_day ,
						       first_wake_day , last_wake_day , tp_lambda );
	      day_tp_sw_W = tp_W.tp_sw;
	      day_tp_ws_W = tp_W.tp_ws;
	      day_has_wake_window = true;
	    }

	  if ( !day_epochs_valid.empty() )
	    {
	      qc.act_min = *std::min_element( day_epochs_valid.begin() , day_epochs_valid.end() );
	      qc.act_max = *std::max_element( day_epochs_valid.begin() , day_epochs_valid.end() );
	      qc.act_med = vec_median( day_epochs_valid );
	      qc.act_sd  = vec_sd( day_epochs_valid );
	      qc.act_mad = vec_mad( day_epochs_valid , qc.act_med );
	      const double q25 = MiscMath::percentile( day_epochs_valid , 0.25 );
	      const double q75 = MiscMath::percentile( day_epochs_valid , 0.75 );
	      qc.act_iqr = q75 - q25;
	      qc.act_p05 = MiscMath::percentile( day_epochs_valid , qc_nearfloor_q_th );
	      qc.act_p95 = MiscMath::percentile( day_epochs_valid , 0.95 );

	      const double mn = vec_mean( day_epochs_valid );
	      const double cv_denom = std::max( std::fabs( mn ) , 1e-9 );
	      qc.lowvar_cv = qc.act_sd / cv_denom;
	      qc.flat_frac = flat_pair_n > 0 ? flat_pair_match_n / (double)flat_pair_n : 0.0;

	      const double scale_eps = std::max( std::fabs( qc.act_p95 - qc.act_p05 ) * 1e-6 , 1e-9 );
	      const double nearfloor_th = qc.act_p05 + scale_eps;
	      const double nonzero_th   = scale_eps;
	      const double active_th    = qc.act_med + qc.act_mad;

	      int nearfloor_n = 0;
	      int nonzero_n = 0;
	      int lowvar_state_n = 0;
	      int quiet_run = 0 , quiet_run_best = 0;
	      int lowvar_run = 0 , lowvar_run_best = 0;
	      const int lowvar_window = std::max( 3 , (int)std::round( 60.0 * 60.0 / epoch_sec ) );

	      for ( int i = b ; i <= e ; i++ )
		{
		  if ( is_gap[i] )
		    {
		      quiet_run = 0;
		      lowvar_run = 0;
		      continue;
		    }

		  if ( epochs[i] <= nearfloor_th )
		    {
		      nearfloor_n++;
		      quiet_run++;
		      if ( quiet_run > quiet_run_best ) quiet_run_best = quiet_run;
		    }
		  else
		    quiet_run = 0;

		  if ( epochs[i] > nonzero_th ) nonzero_n++;
		  if ( epochs[i] > active_th ) qc.active_epoch_n++;

		  std::vector<double> win;
		  win.reserve( lowvar_window );
		  for ( int j = i ; j <= e && (int)win.size() < lowvar_window ; j++ )
		    {
		      if ( is_gap[j] ) break;
		      win.push_back( epochs[j] );
		    }

		  bool lowvar_state = false;
		  if ( win.size() >= 3 )
		    {
		      const double win_mn = vec_mean( win );
		      const double win_cv = vec_sd( win ) / std::max( std::fabs( win_mn ) , 1e-9 );
		      lowvar_state = win_cv <= qc_lowvar_cv_th;
		    }
		  if ( lowvar_state )
		    {
		      lowvar_state_n++;
		      lowvar_run++;
		      if ( lowvar_run > lowvar_run_best ) lowvar_run_best = lowvar_run;
		    }
		  else
		    lowvar_run = 0;
		}

	      const double lowvar_frac = lowvar_state_n / (double)day_epochs_valid.size();
	      qc.nearfloor_frac = nearfloor_n / (double)day_epochs_valid.size();
	      qc.nonzero_frac = nonzero_n / (double)day_epochs_valid.size();
	      qc.longest_quiet_run_min = minutes_from_epochs( quiet_run_best );
	      qc.longest_lowvar_run_min = minutes_from_epochs( lowvar_run_best );
	      lowvar_frac_flag = lowvar_frac >= qc_lowvar_frac_th;
	    }

	  int cur_sleep_run = 0 , cur_wake_run = 0;
	  int best_sleep_run = 0 , best_wake_run = 0;
	  int prev_state = -1;
	  for ( int i = b ; i <= e ; i++ )
	    {
	      if ( is_gap[i] )
		{
		  cur_sleep_run = 0;
		  cur_wake_run = 0;
		  prev_state = -1;
		  continue;
		}

	      const int state = is_wake[i] ? 1 : 0;
	      if ( state != prev_state )
		{
		  if ( state == 1 ) qc.wake_run_n++;
		  else qc.sleep_run_n++;
		}

	      if ( state == 1 )
		{
		  cur_wake_run++;
		  cur_sleep_run = 0;
		  if ( cur_wake_run > best_wake_run ) best_wake_run = cur_wake_run;
		}
	      else
		{
		  cur_sleep_run++;
		  cur_wake_run = 0;
		  if ( cur_sleep_run > best_sleep_run ) best_sleep_run = cur_sleep_run;
		}

	      prev_state = state;
	    }

	  qc.longest_sleep_bout_min = minutes_from_epochs( best_sleep_run );
	  qc.longest_wake_bout_min  = minutes_from_epochs( best_wake_run );
	}

      if ( qc_day )
	{
	  // Signal-quality tech flags are only meaningful when an activity
	  // signal is present; skip them in prescored-no-signal mode so that
	  // an all-zero epochs array does not spuriously exclude every day.
	  if ( have_signal )
	    {
	      qc.flag_flat = qc_exclude_flat
		&& qc.flat_frac >= qc_flat_frac_th;
	      qc.flag_lowvar = qc_exclude_lowvar
		&& lowvar_frac_flag
		&& qc.lowvar_cv <= qc_lowvar_cv_th
		&& qc.active_epoch_n < qc_min_active_epochs;
	      qc.flag_nearfloor = qc_exclude_nearfloor
		&& qc.nearfloor_frac >= qc_nearfloor_frac_th
		&& qc.active_epoch_n < qc_min_active_epochs;
	    }

	  qc.tech_flag_n = ( qc.flag_flat ? 1 : 0 )
	    + ( qc.flag_lowvar ? 1 : 0 )
	    + ( qc.flag_nearfloor ? 1 : 0 );

	  qc.warn_highsleep = qc_warn_highsleep
	    && day_total_valid > 0
	    && ( day_sleep_pct / 100.0 ) >= qc_warn_sleep_pct;
	  qc.warn_longsleep = qc_warn_longsleep
	    && qc.longest_sleep_bout_min >= qc_warn_longsleep_h * 60.0;
	  if ( qc_warn_longsleep && qc.longest_sleep_bout_min >= qc_warn_max_sleep_h * 60.0 )
	    qc.warn_longsleep = true;
	  qc.warn_lowwakeruns = qc.wake_run_n <= qc_warn_low_wakeruns;
	  qc.warn_implausible = qc_warn_implausible
	    && ( day_sleep_pct / 100.0 ) >= qc_warn_sleep_pct
	    && qc.longest_sleep_bout_min >= qc_warn_longsleep_h * 60.0
	    && qc.active_epoch_n < qc_min_active_epochs;

	  // Extreme sleep/wake ratio flags (applied to included days only,
	  // requires day_total_valid > 0 so sleep_pct is meaningful)
	  if ( day_total_valid > 0 )
	    {
	      qc.flag_allsleep = qc_exclude_allsleep
		&& ( day_sleep_pct / 100.0 ) >= qc_allsleep_th;
	      qc.flag_allwake  = qc_exclude_allwake
		&& ( day_sleep_pct / 100.0 ) <= ( 1.0 - qc_allwake_th );
	    }
	  qc.day_excluded_extreme = qc.flag_allsleep || qc.flag_allwake;

	  qc.day_warn = qc.warn_highsleep || qc.warn_longsleep
	    || qc.warn_lowwakeruns || qc.warn_implausible;
	  qc.day_excluded_tech = qc.tech_flag_n >= 2;
	}

      qc.day_excluded = ( !included ) || qc.day_excluded_tech || qc.day_excluded_extreme;
      qc.day_warn = included && !qc.day_excluded_tech && !qc.day_excluded_extreme && qc.day_warn;
      qc.day_ok = included && !qc.day_excluded_tech && !qc.day_excluded_extreme && !qc.day_warn;
      day_qc[d] = qc;

      const int day_num = (int)d + 1;
      const std::string day_name = "day" + Helper::int2str( day_num );
      double day_start_sec = 0.0;
      double day_stop_sec = total_sec;
      if ( d == 0 )
	day_stop_sec = std::min( first_boundary_sec , total_sec );
      else
	{
	  day_start_sec = first_boundary_sec + ( d - 1 ) * 86400.0;
	  day_stop_sec = std::min( first_boundary_sec + d * 86400.0 , total_sec );
	}

      if ( qc.day_excluded_tech && a_qc_excluded != NULL && day_stop_sec > day_start_sec )
	a_qc_excluded->add( "." ,
			    interval_t( Helper::sec2tp( day_start_sec ) ,
					Helper::sec2tp( day_stop_sec ) ) ,
			    "." );

      if ( qc.day_excluded_extreme && a_qc_extreme != NULL && day_stop_sec > day_start_sec )
	a_qc_extreme->add( "." ,
			   interval_t( Helper::sec2tp( day_start_sec ) ,
				       Helper::sec2tp( day_stop_sec ) ) ,
			   "." );

      writer.level( day_num , "DAY" );
      writer.value( "SCORE_TST"     , sleep_sec / 60.0 );
      writer.value( "SCORE_TWT"     , wake_sec  / 60.0 );
      writer.value( "SCORE_TST_SP"  , sleep_sec_sp / 60.0 );
      writer.value( "SCORE_TWT_SP"  , wake_sec_sp  / 60.0 );
      writer.value( "SCORE_TST_WP"  , sleep_sec_wp / 60.0 );
      writer.value( "SCORE_TWT_WP"  , wake_sec_wp  / 60.0 );
      writer.value( "VALID_MIN"         , valid_min_day );
      writer.value( "INCLUDED"          , included ? 1 : 0 );
      writer.value( "SCORE_EPOCH_N"     , day_epochs );
      writer.value( "SCORE_WAKE_EPOCH_N", day_wake_epochs );
      writer.value( "QC_DAY_OK"         , qc.day_ok       ? 1 : 0 );
      writer.value( "QC_DAY_EXCLUDED"   , qc.day_excluded ? 1 : 0 );
      writer.value( "QC_DAY_WARN"       , qc.day_warn     ? 1 : 0 );
      if ( day_total_valid > 0 )
	writer.value( "SCORE_SLEEP_PCT" , day_sleep_pct );
      if ( day_has_sleep_window )
	{
	  writer.value( "FRAG_SFI"         , day_sfi_luna_h );
	  writer.value( "FRAG_SFI_N"       , day_trans_win );
	  writer.value( "FRAG_SFI_ACT"     , day_sfi_act );
	  writer.value( "FRAG_MI"          , day_mi_pct );
	  writer.value( "FRAG_FI"          , day_imm1_pct );
	  writer.value( "FRAG_WAKE_BOUT_N"   , day_frag_wake_bout_n );
	  writer.value( "FRAG_WAKE_BOUT_MED" , day_frag_wake_bout_med );
	  writer.value( "FRAG_WAKE_BOUT_P90" , day_frag_wake_bout_p90 );
	  writer.value( "FRAG_WAKE_BOUT_MAX" , day_frag_wake_bout_max );
	  writer.value( "FRAG_SLEEP_BOUT_MED", day_frag_sleep_bout_med );
	  writer.value( "FRAG_SLEEP_BOUT_P10", day_frag_sleep_bout_p10 );
	  writer.value( "FRAG_SLEEP_BOUT_MAX", day_frag_sleep_bout_max );
	  writer.value( "FRAG_TRANS_ENT_S"       , day_trans_ent_S );
	  writer.value( "FRAG_RL_ENT_SE_S"       , day_rl_ent_se_S );
	  writer.value( "FRAG_RL_ENT_WE_S"       , day_rl_ent_we_S );
	  writer.value( "FRAG_TP_SW_S"           , day_tp_sw_S );
	  writer.value( "FRAG_TP_WS_S"           , day_tp_ws_S );
	}
      if ( day_has_wake_window )
	{
	  writer.value( "FRAG_TRANS_ENT_W" , day_trans_ent_W );
	  writer.value( "FRAG_RL_ENT_SE_W" , day_rl_ent_se_W );
	  writer.value( "FRAG_RL_ENT_WE_W" , day_rl_ent_we_W );
	  writer.value( "FRAG_TP_SW_W"     , day_tp_sw_W );
	  writer.value( "FRAG_TP_WS_W"     , day_tp_ws_W );
	}
      if ( day_epochs > 0 )
	{
	  writer.value( "QC_TECH_FLAG_N"          , qc.tech_flag_n );
	  writer.value( "QC_FLAT_FRAC"            , qc.flat_frac );
	  writer.value( "QC_LOWVAR_CV"            , qc.lowvar_cv );
	  writer.value( "QC_NEARFLOOR_FRAC"       , qc.nearfloor_frac );
	  writer.value( "QC_ACTIVE_EPOCH_N"       , qc.active_epoch_n );
	  writer.value( "QC_LONGEST_QUIET_RUN_MIN"  , qc.longest_quiet_run_min );
	  writer.value( "QC_LONGEST_LOWVAR_RUN_MIN" , qc.longest_lowvar_run_min );
	  writer.value( "QC_ACT_MED"              , qc.act_med );
	  writer.value( "QC_ACT_MAD"              , qc.act_mad );
	  writer.value( "QC_ACT_SD"               , qc.act_sd );
	  writer.value( "QC_ACT_P05"              , qc.act_p05 );
	  writer.value( "QC_ACT_P95"              , qc.act_p95 );
	  writer.value( "QC_LONGEST_SLEEP_BOUT_MIN" , qc.longest_sleep_bout_min );
	  writer.value( "QC_LONGEST_WAKE_BOUT_MIN"  , qc.longest_wake_bout_min );
	  writer.value( "QC_FLAG_FLAT"            , qc.flag_flat           ? 1 : 0 );
	  writer.value( "QC_FLAG_LOWVAR"          , qc.flag_lowvar         ? 1 : 0 );
	  writer.value( "QC_FLAG_NEARFLOOR"       , qc.flag_nearfloor      ? 1 : 0 );
	  writer.value( "QC_FLAG_ALLSLEEP"        , qc.flag_allsleep       ? 1 : 0 );
	  writer.value( "QC_FLAG_ALLWAKE"         , qc.flag_allwake        ? 1 : 0 );
	  writer.value( "QC_DAY_EXCLUDED_EXTREME" , qc.day_excluded_extreme? 1 : 0 );
	  writer.value( "QC_WARN_HIGHSLEEP"       , qc.warn_highsleep  ? 1 : 0 );
	  writer.value( "QC_WARN_LONGSLEEP"       , qc.warn_longsleep  ? 1 : 0 );
	  writer.value( "QC_WARN_LOWWAKERUNS"     , qc.warn_lowwakeruns ? 1 : 0 );
	  writer.value( "QC_WARN_IMPLAUSIBLE"     , qc.warn_implausible ? 1 : 0 );
	}
      writer.unlevel( "DAY" );

      if ( included )
	{
	  n_included_days++;
	}
      if ( qc.day_ok )
	{
	  sum_sleep_sec += sleep_sec;
	  sum_wake_sec  += wake_sec;
	  sum_sleep_sec_sp += sleep_sec_sp;  sum_wake_sec_sp += wake_sec_sp;
	  sum_sleep_sec_wp += sleep_sec_wp;  sum_wake_sec_wp += wake_sec_wp;
	  n_included_postqc_days++;
	  n_qc_ok_days++;
	}
      else if ( qc.day_warn )
	{
	  sum_sleep_sec += sleep_sec;
	  sum_wake_sec  += wake_sec;
	  sum_sleep_sec_sp += sleep_sec_sp;  sum_wake_sec_sp += wake_sec_sp;
	  sum_sleep_sec_wp += sleep_sec_wp;  sum_wake_sec_wp += wake_sec_wp;
	  n_included_postqc_days++;
	  n_qc_warn_days++;
	}
      if ( ( qc.day_ok || qc.day_warn ) && day_has_sleep_window )
	{
	  sum_frag_sfi           += day_sfi_luna_h;
	  sum_frag_mi_act        += day_mi_pct;
	  sum_frag_fi_act        += day_imm1_pct;
	  sum_frag_sfi_act       += day_sfi_act;
	  sum_frag_wake_bout_med += day_frag_wake_bout_med;
	  sum_frag_wake_bout_p90 += day_frag_wake_bout_p90;
	  sum_frag_wake_bout_max += day_frag_wake_bout_max;
	  sum_frag_sleep_bout_med+= day_frag_sleep_bout_med;
	  sum_frag_sleep_bout_p10+= day_frag_sleep_bout_p10;
	  sum_frag_sleep_bout_max+= day_frag_sleep_bout_max;
	  sum_frag_trans_ent_S   += day_trans_ent_S;
	  sum_frag_rl_ent_se_S   += day_rl_ent_se_S;
	  sum_frag_rl_ent_we_S   += day_rl_ent_we_S;
	  sum_frag_tp_sw_S       += day_tp_sw_S;
	  sum_frag_tp_ws_S       += day_tp_ws_S;
	  ++n_frag_dayavg;
	}
      if ( ( qc.day_ok || qc.day_warn ) && day_has_wake_window )
	{
	  sum_frag_trans_ent_W   += day_trans_ent_W;
	  sum_frag_rl_ent_se_W   += day_rl_ent_se_W;
	  sum_frag_rl_ent_we_W   += day_rl_ent_we_W;
	  sum_frag_tp_sw_W       += day_tp_sw_W;
	  sum_frag_tp_ws_W       += day_tp_ws_W;
	  ++n_frag_W_dayavg;
	}
      if ( qc.day_excluded ) n_qc_excluded_days++;

      if ( qc_day )
	{
	  std::vector<std::string> tech_names;
	  std::vector<std::string> warn_names;
	  if ( qc.flag_flat ) tech_names.push_back( "flat" );
	  if ( qc.flag_lowvar ) tech_names.push_back( "lowvar" );
	  if ( qc.flag_nearfloor ) tech_names.push_back( "nearfloor" );
	  if ( qc.warn_highsleep ) warn_names.push_back( "highsleep" );
	  if ( qc.warn_longsleep ) warn_names.push_back( "longsleep" );
	  if ( qc.warn_lowwakeruns ) warn_names.push_back( "lowwakeruns" );
	  if ( qc.warn_implausible ) warn_names.push_back( "possible_nonwear" );
	  if ( do_verbose )
	    logger << "  day " << day_name
		   << ": valid=" << valid_min_day << "m"
		   << ", excluded=" << ( qc.day_excluded ? "yes" : "no" )
		   << ", tech-flags=" << qc.tech_flag_n
		   << " [" << join_flags( tech_names ) << "]"
		   << ", warn=" << warn_names.size()
		   << " [" << join_flags( warn_names ) << "]\n";
	}
    }

  if ( n_days > 0 )
    {
      writer.value( "DAY_N"          , (int)n_days );
      writer.value( "DAY_N_INCLUDED" , n_included_days );
      writer.value( "DAY_N_EXCLUDED" , (int)n_days - n_included_days );
      writer.value( "DAY_N_QC_OK" , n_qc_ok_days );
      writer.value( "DAY_N_QC_EXCLUDED" , n_qc_excluded_days );
      writer.value( "DAY_N_QC_WARN" , n_qc_warn_days );
      writer.value( "DAY_N_INCLUDED_POSTQC" , n_included_postqc_days );

      if ( n_included_postqc_days > 0 )
	{
	  const double total_scored_sec = sum_sleep_sec + sum_wake_sec;
	  const double n_pqc = (double)n_included_postqc_days;
	  writer.value( "SCORE_TST_DAYAVG"      , ( sum_sleep_sec    / n_pqc ) / 60.0 );
	  writer.value( "SCORE_TWT_DAYAVG"      , ( sum_wake_sec     / n_pqc ) / 60.0 );
	  writer.value( "SCORE_TST_SP_DAYAVG"   , ( sum_sleep_sec_sp / n_pqc ) / 60.0 );
	  writer.value( "SCORE_TWT_SP_DAYAVG"   , ( sum_wake_sec_sp  / n_pqc ) / 60.0 );
	  writer.value( "SCORE_TST_WP_DAYAVG"   , ( sum_sleep_sec_wp / n_pqc ) / 60.0 );
	  writer.value( "SCORE_TWT_WP_DAYAVG"   , ( sum_wake_sec_wp  / n_pqc ) / 60.0 );
	  writer.value( "SCORE_SLEEP_PCT_DAYAVG" ,
			total_scored_sec > 0
			? 100.0 * sum_sleep_sec / total_scored_sec : 0.0 );
	}
      if ( n_frag_dayavg > 0 )
	{
	  const double n = (double)n_frag_dayavg;
	  writer.value( "FRAG_SFI"            , sum_frag_sfi           / n );
	  writer.value( "FRAG_MI"             , sum_frag_mi_act        / n );
	  writer.value( "FRAG_FI"             , sum_frag_fi_act        / n );
	  writer.value( "FRAG_SFI_ACT"        , sum_frag_sfi_act       / n );
	  writer.value( "FRAG_WAKE_BOUT_MED"  , sum_frag_wake_bout_med / n );
	  writer.value( "FRAG_WAKE_BOUT_P90"  , sum_frag_wake_bout_p90 / n );
	  writer.value( "FRAG_WAKE_BOUT_MAX"  , sum_frag_wake_bout_max / n );
	  writer.value( "FRAG_SLEEP_BOUT_MED" , sum_frag_sleep_bout_med/ n );
	  writer.value( "FRAG_SLEEP_BOUT_P10" , sum_frag_sleep_bout_p10/ n );
	  writer.value( "FRAG_SLEEP_BOUT_MAX" , sum_frag_sleep_bout_max/ n );
	  writer.value( "FRAG_TRANS_ENT_S"    , sum_frag_trans_ent_S   / n );
	  writer.value( "FRAG_RL_ENT_SE_S"   , sum_frag_rl_ent_se_S   / n );
	  writer.value( "FRAG_RL_ENT_WE_S"   , sum_frag_rl_ent_we_S   / n );
	  writer.value( "FRAG_TP_SW_S"       , sum_frag_tp_sw_S       / n );
	  writer.value( "FRAG_TP_WS_S"       , sum_frag_tp_ws_S       / n );
	  writer.value( "FRAG_N"              , n_frag_dayavg );
	}
      if ( n_frag_W_dayavg > 0 )
	{
	  const double nw = (double)n_frag_W_dayavg;
	  writer.value( "FRAG_TRANS_ENT_W"    , sum_frag_trans_ent_W   / nw );
	  writer.value( "FRAG_RL_ENT_SE_W"    , sum_frag_rl_ent_se_W   / nw );
	  writer.value( "FRAG_RL_ENT_WE_W"    , sum_frag_rl_ent_we_W   / nw );
	  writer.value( "FRAG_TP_SW_W"        , sum_frag_tp_sw_W       / nw );
	  writer.value( "FRAG_TP_WS_W"        , sum_frag_tp_ws_W       / nw );
	}
    }

  if ( do_verbose && n_days > 0 )
    {
      logger << "  per-day output: " << n_days << " day(s), "
	     << n_included_days << " included (>= "
	     << day_min_valid_min << " valid min)"
	     << ( startdatetime.valid ? "" : " [no valid start datetime; treating entire record as one day]" )
	     << "\n";
      if ( qc_day )
	logger << "  per-day QC: " << n_days << " day(s), "
	       << n_included_postqc_days << " included, "
	       << n_qc_excluded_days << " excluded, "
	       << n_qc_warn_days << " warning-only\n";
    }

  logger << "  scored " << n_epochs << " epoch bins: "
	 << total_sleep << " sleep (" << total_sleep_min << " min), "
	 << total_wake  << " wake ("  << total_wake_min  << " min), "
	 << n_gap_epochs << " gap (" << n_gap_epochs * epoch_sec / 60.0 << " min), "
	 << sleep_pct << "% sleep (of valid)\n";


  // -----------------------------------------------------------------------
  //
  // Local sleep debt
  //
  //  Computes how much the K nights immediately before a target night
  //  deviated from that individual's own baseline sleep over the prior
  //  B nights.  Both quantity (TST) and continuity (SFI) are summarised.
  //  Relative and z-score metrics are only emitted when sufficient
  //  baseline data are available (debt-min-base / debt-min-z).
  //
  //  Reference: within-person z-score approach used in actigraphy + PSG
  //  studies (e.g. Drake et al. 2013; Kroenke et al. 2019).
  //
  // -----------------------------------------------------------------------

  if ( do_debt )
    {

      if ( debt_target_str.empty() )
	Helper::halt( "ACTIG: debt requires debt-target (day number or YYYY-MM-DD date)" );

      if ( debt_R < 1 )
	Helper::halt( "ACTIG: debt-recent must be >= 1" );
      if ( debt_B < 1 )
	Helper::halt( "ACTIG: debt-base must be >= 1" );
      if ( debt_w < 0.0 || debt_w > 1.0 )
	Helper::halt( "ACTIG: debt-w must be in 0..1" );


      // ------------------------------------------------------------------
      // Resolve target to a 0-based day index
      // ------------------------------------------------------------------

      int target_idx = -1;

      int trial_int = 0;
      const bool debt_target_is_int = debt_target_str.find_first_not_of( "0123456789" ) == std::string::npos
	                            && Helper::str2int( debt_target_str , &trial_int );
      if ( debt_target_is_int )
	{
	  // Integer: user-supplied 1-based day number
	  if ( trial_int < 1 || trial_int > (int)n_days )
	    Helper::halt( "ACTIG: debt-target day " + debt_target_str
			  + " out of range (1.." + Helper::int2str( (int)n_days ) + ")" );
	  target_idx = trial_int - 1;
	}
      else
	{
	  // Date string: requires valid EDF start datetime
	  if ( !startdatetime.valid )
	    Helper::halt( "ACTIG: debt-target as date requires a valid EDF start date/time" );

	  // Detect format: length-10 with 4th char '-' → YYYY-MM-DD (YMD); else DMY
	  const date_format_t dfmt =
	    ( debt_target_str.size() == 10 && debt_target_str[4] == '-' ) ? YMD : DMY;

	  const date_t target_date( debt_target_str , dfmt );
	  const date_t start_date ( edf.header.startdate );  // EDF header: DD.MM.YY → DMY

	  // Offset from recording start to the anchor time (noon by default)
	  // of the target date.  The ACTIG day that starts at that anchor
	  // contains the PSG night.
	  const double target_sec = (double)target_date.diff( start_date ) * 86400.0
	    + (double)day_anchor_hr * 3600.0
	    - start_tod_sec;

	  target_idx = day_index_for_sec( target_sec );

	  if ( target_idx < 0 || target_idx >= (int)n_days )
	    {
	      logger << "  debt: warning: debt-target date " << debt_target_str
		     << " maps to day " << ( target_idx + 1 )
		     << " which is outside the recording (1.."
		     << (int)n_days << "); skipping debt outputs\n";
	      target_idx = -1;
	    }
	}

      writer.value( "DEBT_MAPPED" , target_idx >= 0 ? 1 : 0 );

      if ( target_idx >= 0 )
	{

      logger << "  debt: target=day" << ( target_idx + 1 )
	     << ", recent=" << debt_R
	     << ", base=" << debt_B
	     << ", min-base=" << debt_min_base
	     << ", min-z=" << debt_min_z
	     << ", w(TST)=" << debt_w << "\n";


      // ------------------------------------------------------------------
      // Per-day helpers  (post-QC inclusion = INCLUDED and not QC excluded)
      // ------------------------------------------------------------------

      auto day_is_included = [&]( int d ) -> bool
      {
	if ( d < 0 || d >= (int)n_days ) return false;
	return d < (int)day_qc.size() ? day_qc[d].day_ok || day_qc[d].day_warn : false;
      };

      auto day_tst_min = [&]( int d ) -> double
      {
	return ( d < (int)sleep_sec_day.size() ? sleep_sec_day[d] : 0.0 ) / 60.0;
      };

      // SFI = sleep→wake transitions per hour of sleep; 0 if no sleep
      auto day_sfi_val = [&]( int d ) -> double
      {
	const double sh = ( d < (int)sleep_sec_day.size() ? sleep_sec_day[d] : 0.0 ) / 3600.0;
	const int    tr = ( d < (int)day_trans_to_w.size() ? day_trans_to_w[d] : 0 );
	return sh > 0 ? tr / sh : 0.0;
      };


      // ------------------------------------------------------------------
      // Collect values for the two windows
      //   recent   : [ target - R  ..  target - 1 ]
      //   baseline : [ target - R - B  ..  target - R - 1 ]
      // Only INCLUDED days contribute; days with no sleep are excluded
      // from the SFI window (but not the TST window).
      // ------------------------------------------------------------------

      std::vector<double> recent_tst, recent_sfi;
      std::vector<double> base_tst,   base_sfi;

      for ( int d = target_idx - debt_R ; d <= target_idx - 1 ; d++ )
	if ( day_is_included(d) )
	  {
	    recent_tst.push_back( day_tst_min(d) );
	    const double s = day_sfi_val(d);
	    if ( s > 0 ) recent_sfi.push_back( s );
	  }

      for ( int d = target_idx - debt_R - debt_B ; d <= target_idx - debt_R - 1 ; d++ )
	if ( day_is_included(d) )
	  {
	    base_tst.push_back( day_tst_min(d) );
	    const double s = day_sfi_val(d);
	    if ( s > 0 ) base_sfi.push_back( s );
	  }

      const int n_recent = (int)recent_tst.size();
      const int n_base   = (int)base_tst.size();

      logger << "  debt: " << n_recent << " recent night(s), "
	     << n_base << " baseline night(s) with valid data\n";


      // ------------------------------------------------------------------
      // Output
      // ------------------------------------------------------------------

      writer.value( "DEBT_TARGET_DAY" , target_idx + 1 );
      writer.value( "DEBT_N_RECENT"   , n_recent );
      writer.value( "DEBT_N_BASE"     , n_base );

      // Per-window absolute summaries
      if ( n_recent > 0 )
	{
	  writer.value( "DEBT_TST_RECENT"  , vec_mean( recent_tst ) );
	  writer.value( "DEBT_FRAG_RECENT" ,
			recent_sfi.empty() ? 0.0 : vec_mean( recent_sfi ) );
	}

      if ( n_base > 0 )
	{
	  writer.value( "DEBT_TST_BASE"  , vec_median( base_tst ) );
	  writer.value( "DEBT_FRAG_BASE" ,
			base_sfi.empty() ? 0.0 : vec_median( base_sfi ) );
	}

      // Relative metrics: require n_base >= debt_min_base
      if ( n_recent > 0 && n_base >= debt_min_base )
	{
	  const double r_tst     = vec_mean( recent_tst );
	  const double b_tst_med = vec_median( base_tst );
	  const double delta     = b_tst_med - r_tst;
	  writer.value( "DEBT_TST_DELTA" , delta );
	  writer.value( "DEBT_TST_REL"   , b_tst_med > 0 ? delta / b_tst_med : 0.0 );
	}

      // Z-score metrics + composite: require n_base >= debt_min_z
      if ( n_recent > 0 && n_base >= debt_min_z )
	{
	  const double r_tst    = vec_mean( recent_tst );
	  const double b_tst_mn = vec_mean( base_tst );
	  const double b_tst_sd = vec_sd( base_tst );

	  double tst_z    = 0.0;
	  bool   has_tst_z  = false;

	  if ( b_tst_sd > 0 )
	    {
	      tst_z     = ( r_tst - b_tst_mn ) / b_tst_sd;
	      has_tst_z = true;
	      writer.value( "DEBT_TST_Z" , tst_z );
	    }

	  // Fragmentation Z: additionally require SFI data in both windows
	  double frag_z   = 0.0;
	  bool   has_frag_z = false;

	  if ( !recent_sfi.empty() && (int)base_sfi.size() >= debt_min_z )
	    {
	      const double r_sfi    = vec_mean( recent_sfi );
	      const double b_sfi_mn = vec_mean( base_sfi );
	      const double b_sfi_sd = vec_sd( base_sfi );

	      if ( b_sfi_sd > 0 )
		{
		  frag_z     = ( r_sfi - b_sfi_mn ) / b_sfi_sd;
		  has_frag_z = true;
		  writer.value( "DEBT_FRAG_Z" , frag_z );
		}
	    }

	  // Composite DEBT_INDEX:  w * (−TST_Z) + (1−w) * FRAG_Z
	  // Higher = more deprived.  Falls back gracefully if only one
	  // component is available.
	  if ( has_tst_z || has_frag_z )
	    {
	      double index;
	      if      ( has_tst_z && has_frag_z )
		index = debt_w * ( -tst_z ) + ( 1.0 - debt_w ) * frag_z;
	      else if ( has_tst_z  )
		index = -tst_z;
	      else
		index = frag_z;
	      writer.value( "DEBT_INDEX" , index );
	    }
	}

	}  // end if target_idx >= 0

    }  // end do_debt


  // -----------------------------------------------------------------------
  //
  // All-nights normalisation
  //
  //  Z-scores the target night's TST (and SFI) against all other valid
  //  nights in the recording.  Unlike debt, no prior baseline window is
  //  required, so this works even when actigraphy starts at PSG night 1.
  //
  //  Reference pool: all days that pass QC (day_ok or day_warn) other
  //  than the target night itself.  Outputs are suppressed when the
  //  reference pool has fewer than norm-min nights (default 5).
  //
  // -----------------------------------------------------------------------

  if ( do_norm )
    {

      if ( norm_target_str.empty() )
	Helper::halt( "ACTIG: norm requires norm-target (day number or YYYY-MM-DD date)" );

      // ------------------------------------------------------------------
      // Resolve target to a 0-based day index (same logic as do_debt)
      // ------------------------------------------------------------------

      int target_idx = -1;

      int trial_int = 0;
      const bool norm_target_is_int = norm_target_str.find_first_not_of( "0123456789" ) == std::string::npos
	                            && Helper::str2int( norm_target_str , &trial_int );
      if ( norm_target_is_int )
	{
	  if ( trial_int < 1 || trial_int > (int)n_days )
	    Helper::halt( "ACTIG: norm-target day " + norm_target_str
			  + " out of range (1.." + Helper::int2str( (int)n_days ) + ")" );
	  target_idx = trial_int - 1;
	}
      else
	{
	  if ( !startdatetime.valid )
	    Helper::halt( "ACTIG: norm-target as date requires a valid EDF start date/time" );

	  const date_format_t dfmt =
	    ( norm_target_str.size() == 10 && norm_target_str[4] == '-' ) ? YMD : DMY;

	  const date_t target_date( norm_target_str , dfmt );
	  const date_t start_date ( edf.header.startdate );

	  const double target_sec = (double)target_date.diff( start_date ) * 86400.0
	    + (double)day_anchor_hr * 3600.0
	    - start_tod_sec;

	  target_idx = day_index_for_sec( target_sec );

	  if ( target_idx < 0 || target_idx >= (int)n_days )
	    {
	      logger << "  norm: warning: norm-target date " << norm_target_str
		     << " maps to day " << ( target_idx + 1 )
		     << " which is outside the recording (1.."
		     << (int)n_days << "); skipping norm outputs\n";
	      target_idx = -1;
	    }
	}

      writer.value( "NORM_MAPPED" , target_idx >= 0 ? 1 : 0 );

      if ( target_idx >= 0 )
	{

      // ------------------------------------------------------------------
      // Per-day helpers
      // ------------------------------------------------------------------

      auto day_is_included_n = [&]( int d ) -> bool
      {
	if ( d < 0 || d >= (int)n_days ) return false;
	return d < (int)day_qc.size() ? day_qc[d].day_ok || day_qc[d].day_warn : false;
      };

      auto day_tst_min_n = [&]( int d ) -> double
      {
	return ( d < (int)sleep_sec_day.size() ? sleep_sec_day[d] : 0.0 ) / 60.0;
      };

      auto day_sfi_val_n = [&]( int d ) -> double
      {
	const double sh = ( d < (int)sleep_sec_day.size() ? sleep_sec_day[d] : 0.0 ) / 3600.0;
	const int    tr = ( d < (int)day_trans_to_w.size() ? day_trans_to_w[d] : 0 );
	return sh > 0 ? tr / sh : 0.0;
      };

      // ------------------------------------------------------------------
      // Collect all valid nights except the target
      // ------------------------------------------------------------------

      std::vector<double> ref_tst, ref_sfi;

      for ( int d = 0; d < (int)n_days; d++ )
	{
	  if ( d == target_idx ) continue;
	  if ( !day_is_included_n(d) ) continue;
	  ref_tst.push_back( day_tst_min_n(d) );
	  const double s = day_sfi_val_n(d);
	  if ( s > 0 ) ref_sfi.push_back( s );
	}

      const int n_ref = (int)ref_tst.size();

      logger << "  norm: target=day" << ( target_idx + 1 )
	     << ", " << n_ref << " reference night(s), min=" << norm_min << "\n";

      // ------------------------------------------------------------------
      // Output
      // ------------------------------------------------------------------

      const double tgt_tst = day_tst_min_n( target_idx );
      const double tgt_sfi = day_sfi_val_n( target_idx );

      writer.value( "NORM_TARGET_DAY" , target_idx + 1 );
      writer.value( "NORM_N"          , n_ref );
      writer.value( "NORM_TST"        , tgt_tst );
      writer.value( "NORM_FRAG"       , tgt_sfi );

      if ( n_ref >= norm_min )
	{
	  const double mn_tst = vec_mean( ref_tst );
	  const double sd_tst = vec_sd( ref_tst );

	  writer.value( "NORM_TST_MN" , mn_tst );
	  writer.value( "NORM_TST_SD" , sd_tst );

	  double tst_z    = 0.0;
	  bool   has_tst_z  = false;

	  if ( sd_tst > 0 )
	    {
	      tst_z     = ( tgt_tst - mn_tst ) / sd_tst;
	      has_tst_z = true;
	      writer.value( "NORM_TST_Z" , tst_z );
	    }

	  double frag_z   = 0.0;
	  bool   has_frag_z = false;

	  if ( (int)ref_sfi.size() >= norm_min )
	    {
	      const double mn_sfi = vec_mean( ref_sfi );
	      const double sd_sfi = vec_sd( ref_sfi );

	      writer.value( "NORM_FRAG_MN" , mn_sfi );
	      writer.value( "NORM_FRAG_SD" , sd_sfi );

	      if ( sd_sfi > 0 )
		{
		  frag_z     = ( tgt_sfi - mn_sfi ) / sd_sfi;
		  has_frag_z = true;
		  writer.value( "NORM_FRAG_Z" , frag_z );
		}
	    }

	  // Composite index: w*(-TST_Z) + (1-w)*FRAG_Z  (higher = more deprived)
	  if ( has_tst_z || has_frag_z )
	    {
	      double index;
	      if      ( has_tst_z && has_frag_z )
		index = norm_w * ( -tst_z ) + ( 1.0 - norm_w ) * frag_z;
	      else if ( has_tst_z )
		index = -tst_z;
	      else
		index = frag_z;
	      writer.value( "NORM_INDEX" , index );
	    }
	}

	}  // end if target_idx >= 0

    }  // end do_norm

}  // end actig::actig()


// -----------------------------------------------------------------------
//
// DAYS : create day and hourly clock-time annotations
//
// -----------------------------------------------------------------------

void actig::days( edf_t & edf , param_t & param )
{

  //
  // Requires a valid EDF start date and time
  //

  clocktime_t startdatetime( edf.header.startdate , edf.header.starttime );

  if ( ! startdatetime.valid )
    Helper::halt( "DAYS requires a valid EDF start date and time" );

  //
  // Recording duration in seconds
  //

  const double total_sec = edf.timeline.total_duration_tp / (double)globals::tp_1sec;

  if ( total_sec <= 0 ) return;


  //
  // Options
  //

  const int anchor_hr = param.has( "anchor" ) ? param.requires_int( "anchor" ) : 12;

  if ( anchor_hr < 0 || anchor_hr > 23 )
    Helper::halt( "DAYS anchor must be 0..23" );

  const bool add_days = param.yesno( "days" , true , true );

  const std::string day_prefix = param.has( "prefix" ) ? param.value( "prefix" ) : "day";

  const bool add_hours = param.has( "hours" );

  const std::string hr_prefix = param.has( "hour-prefix" ) ? param.value( "hour-prefix" ) : "";

  const bool add_halves = param.has( "halves" );

  // Arbitrary bin duration
  const bool add_bins = param.has( "bin" );
  double bin_dur_sec = 3600.0;
  std::string bin_prefix = "b";
  int bin_ndigits = 2;

  if ( add_bins )
    {
      std::string binval = param.value( "bin" );
      double mult = 1.0; // default: seconds
      if ( ! binval.empty() )
	{
	  char last = binval[ binval.size() - 1 ];
	  if ( last == 's' || last == 'S' ) { mult = 1.0;    binval = binval.substr( 0, binval.size() - 1 ); }
	  else if ( last == 'm' || last == 'M' ) { mult = 60.0;   binval = binval.substr( 0, binval.size() - 1 ); }
	  else if ( last == 'h' || last == 'H' ) { mult = 3600.0; binval = binval.substr( 0, binval.size() - 1 ); }
	}
      double bval = 0;
      if ( ! Helper::str2dbl( binval , &bval ) || bval <= 0 )
	Helper::halt( "DAYS bin: invalid duration '" + param.value( "bin" ) + "'" );
      bin_dur_sec = bval * mult;
      if ( bin_dur_sec <= 0 || bin_dur_sec > 86400 )
	Helper::halt( "DAYS bin duration must be > 0 and <= 86400 seconds" );

      bin_prefix = param.has( "bin-prefix" ) ? param.value( "bin-prefix" ) : "b";
      const int n_bins_per_day = (int)std::ceil( 86400.0 / bin_dur_sec );
      bin_ndigits = n_bins_per_day > 999 ? 4 : ( n_bins_per_day > 99 ? 3 : 2 );
    }

  // Clock-window annotations
  const bool add_windows = param.has( "window" );
  const bool add_weekend = param.has( "weekend" );
  const std::string weekend_label = param.has( "weekend-label" ) ? param.value( "weekend-label" ) : "WEEKEND";

  const bool verbose = param.has( "verbose" );
  if ( verbose ) writer.numeric_factor( "DAY" );


  //
  // Compute first day boundary
  //

  double start_tod_sec = startdatetime.h * 3600.0
    + startdatetime.m * 60.0
    + startdatetime.s;

  double first_boundary_sec = first_day_boundary_sec( startdatetime , anchor_hr , total_sec );


  //
  // Build day annotations
  //

  int day_num = 1;
  int n_full_days = 0;
  int n_partial_days = 0;

  double day_start_sec = 0;
  double day_stop_sec  = first_boundary_sec;

  if ( day_stop_sec > total_sec )
    day_stop_sec = total_sec;

  annot_t * weekend_annot = NULL;
  int n_weekend_annots = 0;
  if ( add_weekend && weekend_label != "" )
    {
      weekend_annot = edf.timeline.annotations->add( weekend_label );
      weekend_annot->description = weekend_label;
      weekend_annot->file = "DAYS";
    }

  while ( day_start_sec < total_sec )
    {

      const std::string name = day_prefix
	+ ( day_num < 10 ? "0" : "" )
	+ Helper::int2str( day_num );

      uint64_t tp_start = Helper::sec2tp( day_start_sec );
      uint64_t tp_stop  = Helper::sec2tp( day_stop_sec );

      if ( add_days )
	{
	  annot_t * a = edf.timeline.annotations->add( name );
	  a->description = name;
	  a->file = "DAYS";
	  a->add( "." , interval_t( tp_start , tp_stop ) , "." );
	}

      clocktime_t ds = startdatetime;
      ds.advance_seconds( day_start_sec );

      const bool is_weekend = date_t( ds.as_date_string() , DMY ).is_weekend();

      if ( weekend_annot != NULL && is_weekend )
	{
	  weekend_annot->add( "." , interval_t( tp_start , tp_stop ) , "." );
	  ++n_weekend_annots;
	}

      if ( verbose )
	{
	  clocktime_t de = startdatetime;
	  de.advance_seconds( day_stop_sec );

	  writer.level( day_num , "DAY" );
	  writer.value( "START_SEC" , day_start_sec );
	  writer.value( "STOP_SEC" , day_stop_sec );
	  writer.value( "DUR_SEC" , day_stop_sec - day_start_sec );
	  writer.value( "START_HMS" , ds.as_string( ':' ) );
	  writer.value( "STOP_HMS" , de.as_string( ':' ) );
	  writer.value( "START_DATE" , ds.as_date_string( '/' ) );
	  writer.value( "WEEKEND" , is_weekend ? 1 : 0 );
	}


      if ( day_stop_sec - day_start_sec >= 86400.0 - 1.0 )
	++n_full_days;
      else
	++n_partial_days;

      ++day_num;
      day_start_sec = day_stop_sec;
      day_stop_sec  = day_start_sec + 86400.0;

      if ( day_stop_sec > total_sec )
	day_stop_sec = total_sec;
    }

  if ( verbose )
    writer.unlevel( "DAY" );

  if ( add_days )
    {
      logger << "  DAYS: created " << (n_full_days + n_partial_days) << " day annotations"
	     << " (" << n_full_days << " full, " << n_partial_days << " partial)";
      if ( anchor_hr != 12 )
	logger << " (anchor=" << anchor_hr << "h)";
      logger << "\n";
    }
  if ( weekend_annot != NULL )
    logger << "  DAYS: created " << n_weekend_annots << " weekend annotations\n";


  //
  // Hourly annotations: 00h, 01h, ..., 23h
  //

  if ( add_hours )
    {

      int n_hr_annots = 0;

      for (int hr = 0; hr < 24; hr++)
	{
	  std::string hr_str = ( hr < 10 ? "0" : "" ) + Helper::int2str( hr );
	  const std::string name = hr_prefix + hr_str + "h";

	  annot_t * a = edf.timeline.annotations->add( name );
	  a->description = name;
	  a->file = "DAYS";

	  double hr_tod_sec = hr * 3600.0;

	  double offset = hr_tod_sec - start_tod_sec;
	  if ( offset < 0 )
	    offset += 86400.0;

	  double hr_start = offset;
	  int inst = 1;

	  while ( hr_start < total_sec )
	    {
	      double hr_stop = hr_start + 3600.0;
	      if ( hr_stop > total_sec )
		hr_stop = total_sec;

	      if ( hr_stop > hr_start )
		{
		  a->add( Helper::int2str( inst ) ,
			  interval_t( Helper::sec2tp( hr_start ) ,
				      Helper::sec2tp( hr_stop ) ) ,
			  "." );
		  ++n_hr_annots;
		  ++inst;
		}

	      hr_start += 86400.0;
	    }
	}

      logger << "  DAYS: created " << n_hr_annots << " hourly annotations\n";
    }


  //
  // AM/PM halves
  //

  if ( add_halves )
    {

      int n_half_annots = 0;

      for (int half = 0; half < 2; half++)
	{
	  const std::string name = half == 0 ? "AM" : "PM";
	  double half_tod_sec = half * 43200.0;

	  annot_t * a = edf.timeline.annotations->add( name );
	  a->description = name;
	  a->file = "DAYS";

	  double offset = half_tod_sec - start_tod_sec;
	  if ( offset < 0 )
	    offset += 86400.0;

	  double h_start = offset;
	  int inst = 1;

	  while ( h_start < total_sec )
	    {
	      double h_stop = h_start + 43200.0;
	      if ( h_stop > total_sec )
		h_stop = total_sec;

	      if ( h_stop > h_start )
		{
		  a->add( Helper::int2str( inst ) ,
			  interval_t( Helper::sec2tp( h_start ) ,
				      Helper::sec2tp( h_stop ) ) ,
			  "." );
		  ++n_half_annots;
		  ++inst;
		}

	      h_start += 86400.0;
	    }
	}

      logger << "  DAYS: created " << n_half_annots << " AM/PM annotations\n";
    }


  //
  // Arbitrary clock-time bin annotations
  //

  if ( add_bins )
    {
      int n_bin_annots = 0;
      const int n_bins_per_day = (int)std::ceil( 86400.0 / bin_dur_sec );

      // helper: seconds-past-midnight --> HHMMSS string (allows h=24 for end-of-day)
      auto tod_hhmmss = []( double sec ) -> std::string {
	int h = (int)( sec / 3600.0 );
	int m = (int)( ( sec - h * 3600.0 ) / 60.0 );
	int s = (int)( sec - h * 3600.0 - m * 60.0 );
	return Helper::zero_pad( h , 2 ) + Helper::zero_pad( m , 2 ) + Helper::zero_pad( s , 2 );
      };

      for ( int b = 0; b < n_bins_per_day; b++ )
	{
	  // clock-time offset of bin b from midnight
	  double bin_tod_sec = b * bin_dur_sec;
	  double bin_end_tod = std::min( (b + 1) * bin_dur_sec , 86400.0 );

	  const std::string name = bin_prefix
	    + Helper::zero_pad( b , bin_ndigits )
	    + "_" + tod_hhmmss( bin_tod_sec )
	    + "_" + tod_hhmmss( bin_end_tod );

	  annot_t * a = edf.timeline.annotations->add( name );
	  a->description = name;
	  a->file = "DAYS";

	  double offset = bin_tod_sec - start_tod_sec;
	  if ( offset < 0 ) offset += 86400.0;

	  double b_start = offset;
	  int inst = 1;

	  while ( b_start < total_sec )
	    {
	      double b_stop = b_start + bin_dur_sec;
	      if ( b_stop > total_sec ) b_stop = total_sec;

	      if ( b_stop > b_start )
		{
		  a->add( Helper::int2str( inst ) ,
			  interval_t( Helper::sec2tp( b_start ) ,
				      Helper::sec2tp( b_stop ) ) ,
			  "." );
		  ++n_bin_annots;
		  ++inst;
		}

	      b_start += 86400.0;
	    }
	}

      logger << "  DAYS: created " << n_bin_annots << " bin annotations"
	     << " (" << bin_dur_sec << "s bins, " << n_bins_per_day << " per day)\n";
    }


  //
  // Arbitrary clock-window annotations (repeating each day)
  //

  if ( add_windows )
    {
      const std::vector<std::string> wspecs = param.strvector( "window" , "," );
      int auto_label_idx = 1;

      for ( int wi = 0; wi < (int)wspecs.size(); wi++ )
	{
	  const std::string & spec = wspecs[wi];
	  std::string label;
	  std::string timerange;

	  // Format: label:HH:MM-HH:MM  or  HH:MM-HH:MM
	  if ( ! spec.empty() && ! isdigit( (unsigned char)spec[0] ) )
	    {
	      // has a label prefix before the first ':'
	      const size_t colon_pos = spec.find( ':' );
	      if ( colon_pos == std::string::npos )
		Helper::halt( "DAYS window: bad format '" + spec + "'; expected [label:]HH:MM-HH:MM" );
	      label     = spec.substr( 0 , colon_pos );
	      timerange = spec.substr( colon_pos + 1 );
	    }
	  else
	    {
	      label     = "W" + Helper::int2str( auto_label_idx++ );
	      timerange = spec;
	    }

	  // Parse HH:MM-HH:MM  (find '-' after at least pos 3 to skip past first HH:MM)
	  const size_t dash_pos = timerange.find( '-' , 3 );
	  if ( dash_pos == std::string::npos )
	    Helper::halt( "DAYS window: bad time range '" + timerange + "'; expected HH:MM-HH:MM" );

	  const std::string start_str = timerange.substr( 0 , dash_pos );
	  const std::string end_str   = timerange.substr( dash_pos + 1 );

	  int sh = 0, sm = 0, eh = 0, em = 0;
	  if ( sscanf( start_str.c_str() , "%d:%d" , &sh , &sm ) != 2 ||
	       sscanf( end_str.c_str()   , "%d:%d" , &eh , &em ) != 2 )
	    Helper::halt( "DAYS window: could not parse times in '" + spec + "'" );

	  if ( sh < 0 || sh > 23 || sm < 0 || sm > 59 ||
	       eh < 0 || eh > 23 || em < 0 || em > 59 )
	    Helper::halt( "DAYS window: invalid time values in '" + spec + "'" );

	  const double win_start_tod = sh * 3600.0 + sm * 60.0;
	  const double win_end_tod   = eh * 3600.0 + em * 60.0;

	  // duration; handle midnight-crossing windows
	  double win_dur_sec;
	  if ( win_end_tod > win_start_tod )
	    win_dur_sec = win_end_tod - win_start_tod;
	  else
	    win_dur_sec = 86400.0 - win_start_tod + win_end_tod;

	  if ( win_dur_sec <= 0 )
	    Helper::halt( "DAYS window: zero-duration window in '" + spec + "'" );

	  annot_t * a = edf.timeline.annotations->add( label );
	  a->description = label;
	  a->file = "DAYS";

	  double offset = win_start_tod - start_tod_sec;
	  if ( offset < 0 ) offset += 86400.0;

	  double w_start = offset;
	  int inst = 1;
	  int n_this_window = 0;

	  while ( w_start < total_sec )
	    {
	      double w_stop = w_start + win_dur_sec;
	      if ( w_stop > total_sec ) w_stop = total_sec;

	      if ( w_stop > w_start )
		{
		  a->add( Helper::int2str( inst ) ,
			  interval_t( Helper::sec2tp( w_start ) ,
				      Helper::sec2tp( w_stop ) ) ,
			  "." );
		  ++n_this_window;
		  ++inst;
		}

	      w_start += 86400.0;
	    }

	  logger << "  DAYS: created " << n_this_window
		 << " window annotations for '" << label << "' ("
		 << start_str << "-" << end_str << ")\n";
	}
    }

}
