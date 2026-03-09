
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

#include <cmath>
#include <numeric>
#include <algorithm>
#include <memory>

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


// Fine-step sliding L5/M10 window over the 24-bin hourly profile.
//
// The profile is treated as piecewise-constant (each bin = 1 hour).
// The window is placed at positions t = 0, step_hr, 2*step_hr, ...
// and its mean is computed as the exact weighted integral over the
// piecewise-constant profile, handling circular wrap-around.
//
// step_min == 60  reproduces the traditional integer-hour-step behaviour.
//
static void find_window_fine( const std::vector<double> & hourly_profile ,
			      int    window_hours ,
			      bool   find_min ,
			      double * result_mean ,
			      double * result_onset_hr ,   // fractional hours [0,24)
			      int    step_min = 1 )
{
  const int p = 24;

  if ( (int)hourly_profile.size() != p || window_hours <= 0 || window_hours >= p )
    {
      *result_mean     = 0;
      *result_onset_hr = 0;
      return;
    }

  const double step_hr = step_min / 60.0;
  const int    n_steps = (int)std::round( 24.0 / step_hr );
  const double W       = (double)window_hours;

  double best       = find_min ? 1e300 : -1e300;
  double best_onset = 0.0;

  for ( int si = 0 ; si < n_steps ; si++ )
    {
      const double t   = si * step_hr;  // window onset (fractional hours from profile start)
      const double end = t + W;          // window end (may exceed 24)

      double sum = 0.0;

      for ( int h = 0 ; h < p ; h++ )
	{
	  const double h0 = h , h1 = h + 1.0;
	  double overlap;

	  if ( end <= 24.0 )
	    {
	      // no circular wrap
	      overlap = std::max( 0.0 , std::min( h1 , end ) - std::max( h0 , t ) );
	    }
	  else
	    {
	      // window wraps past midnight: [t,24) ∪ [0, end-24)
	      double ow1 = std::max( 0.0 , std::min( h1 , 24.0   ) - std::max( h0 , t      ) );
	      double ow2 = std::max( 0.0 , std::min( h1 , end-24.0) - std::max( h0 , 0.0   ) );
	      overlap = ow1 + ow2;
	    }

	  sum += hourly_profile[h] * overlap;
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

  const bool do_score = param.has( "score" );
  const bool do_prescored = param.has( "prescored" );

  std::string presc_sleep_label = "S";
  std::string presc_wake_label  = "W";
  if ( do_prescored && ! param.empty( "prescored" ) )
    {
      std::vector<std::string> tok = Helper::parse( param.value( "prescored" ) , "," );
      if ( tok.size() != 2 || tok[0] == "" || tok[1] == "" )
	Helper::halt( "ACTIG: prescored must be empty or exactly two comma-delimited labels: sleep,wake" );
      presc_sleep_label = tok[0];
      presc_wake_label  = tok[1];
    }

  if ( do_score && do_prescored )
    Helper::halt( "ACTIG: cannot specify both score and prescored" );

  const bool need_signal = ! do_prescored;
  const std::string siglab = need_signal ? param.requires( "sig" ) : "";

  int slot = -1;
  double Fs = 0;

  if ( need_signal )
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

  // labels for wake/sleep/gap annotations
  const std::string wake_label  = param.has( "wake" )  ? param.value( "wake" )  : "W";
  const std::string sleep_label = param.has( "sleep" ) ? param.value( "sleep" ) : "S";
  const std::string gap_label   = param.has( "gap-out" ) ? param.value( "gap-out" ) : "ACTIG_G";

  // day boundary anchor hour for per-day TST/WASO summaries
  const int day_anchor_hr = param.has( "day-anchor" ) ? param.requires_int( "day-anchor" ) : 12;
  if ( day_anchor_hr < 0 || day_anchor_hr > 23 )
    Helper::halt( "ACTIG: day-anchor must be 0..23" );

  // gap: minimum % of expected samples for an epoch to be considered valid
  const double gap_min_pct = param.has( "gap-min-pct" ) ? param.requires_dbl( "gap-min-pct" ) : 50.0;
  if ( gap_min_pct < 0 || gap_min_pct > 100 )
    Helper::halt( "ACTIG: gap-min-pct must be 0..100" );

  // minimum valid minutes per day for inclusion in day averages (default 16 h)
  const double day_min_valid_min = param.has( "day-min-valid" ) ? param.requires_dbl( "day-min-valid" ) : 960.0;
  if ( day_min_valid_min < 0 )
    Helper::halt( "ACTIG: day-min-valid must be >= 0" );

  // Luna method parameters
  const double luna_smooth_min    = param.has( "smooth" )      ? param.requires_dbl( "smooth" )      : 10.0;
  const double luna_burst_min     = param.has( "burst" )       ? param.requires_dbl( "burst" )       : 10.0;
  const double luna_burst_z       = param.has( "burst-z" )     ? param.requires_dbl( "burst-z" )     :  0.5;
  const double luna_quiet_z       = param.has( "quiet-z" )     ? param.requires_dbl( "quiet-z" )     : -0.5;
  const double luna_active_frac   = param.has( "active-frac" ) ? param.requires_dbl( "active-frac" ) :  0.20;
  const double luna_min_sleep_min = param.has( "min-sleep" )   ? param.requires_dbl( "min-sleep" )   : 15.0;
  const double luna_max_gap_min   = param.has( "max-gap" )     ? param.requires_dbl( "max-gap" )     :  2.0;
  const double luna_min_wake_min  = param.has( "min-wake" )    ? param.requires_dbl( "min-wake" )    :  5.0;

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


  std::string endorsed_opts, not_endorsed_opts;
  auto add_opt = []( std::string * x , const std::string & s )
  {
    if ( ! x->empty() ) *x += ", ";
    *x += s;
  };
  auto track_opt = [&]( bool on , const std::string & name )
  {
    if ( on ) add_opt( &endorsed_opts , name );
    else add_opt( &not_endorsed_opts , name );
  };

  track_opt( param.has( "epoch" ) , "epoch" );
  track_opt( param.has( "bin" ) , "bin" );
  track_opt( param.has( "l" ) , "l" );
  track_opt( param.has( "m" ) , "m" );
  track_opt( param.has( "np-step" )        , "np-step" );
  track_opt( param.has( "np-traditional" ) , "np-traditional" );
  track_opt( param.has( "sum" ) , "sum" );
  track_opt( param.has( "score" ) , "score" );
  track_opt( param.has( "prescored" ) , "prescored" );
  track_opt( param.has( "method" ) , "method" );
  track_opt( param.has( "thresh" ) , "thresh" );
  track_opt( param.has( "cole-thresh" ) , "cole-thresh" );
  track_opt( param.has( "wake" ) , "wake" );
  track_opt( param.has( "sleep" ) , "sleep" );
  track_opt( param.has( "gap-out" ) , "gap-out" );
  track_opt( param.has( "day-anchor" ) , "day-anchor" );
  track_opt( param.has( "gap-min-pct" ) , "gap-min-pct" );
  track_opt( param.has( "day-min-valid" ) , "day-min-valid" );
  track_opt( param.has( "smooth" )       , "smooth" );
  track_opt( param.has( "burst" )        , "burst" );
  track_opt( param.has( "burst-z" )      , "burst-z" );
  track_opt( param.has( "quiet-z" )      , "quiet-z" );
  track_opt( param.has( "active-frac" )  , "active-frac" );
  track_opt( param.has( "min-sleep" )    , "min-sleep" );
  track_opt( param.has( "max-gap" )      , "max-gap" );
  track_opt( param.has( "min-wake" )     , "min-wake" );
  track_opt( param.has( "channels" )     , "channels" );
  track_opt( param.has( "debt" )         , "debt" );
  track_opt( param.has( "debt-target" )  , "debt-target" );
  track_opt( param.has( "debt-recent" )  , "debt-recent" );
  track_opt( param.has( "debt-base" )    , "debt-base" );
  track_opt( param.has( "debt-min-base" ), "debt-min-base" );
  track_opt( param.has( "debt-min-z" )   , "debt-min-z" );
  track_opt( param.has( "debt-w" )       , "debt-w" );

  logger << "  settings: sig=" << ( need_signal ? siglab : "<prescored>" );
  if ( need_signal ) logger << ", sr=" << Fs << " Hz";
  logger << ", epoch=" << epoch_sec << "s"
	 << ", bin=" << np_bin_min << "m\n";
  logger << "  settings: L=" << l_hours << "h"
	 << ", M=" << m_hours << "h"
	 << ", score=" << ( do_score ? "yes" : "no" )
	 << ", method=" << method
	 << ", sum=" << ( use_sum ? "yes" : "no" ) << "\n";
  logger << "  gap settings: gap-min-pct=" << gap_min_pct
	 << "%, day-min-valid=" << day_min_valid_min << " min\n";
  logger << "  options endorsed: "
	 << ( endorsed_opts.empty() ? "none" : endorsed_opts ) << "\n";
  logger << "  options not endorsed: "
	 << ( not_endorsed_opts.empty() ? "none" : not_endorsed_opts ) << "\n";


  //
  // Extract full signal + timestamps
  //
  //  Gaps (from prior MASK/RE or EDF+D discontinuities) appear as
  //  timestamp jumps in ptimepoints().  We detect them during binning.
  //

  interval_t whole = edf.timeline.wholetrace( true );

  std::unique_ptr<slice_t> slice;
  std::vector<double> empty_d;
  std::vector<uint64_t> empty_tp;
  const std::vector<double>   * d  = &empty_d;
  const std::vector<uint64_t> * tp = &empty_tp;
  if ( need_signal )
    {
      slice.reset( new slice_t( edf , slot , whole ) );
      d = slice->pdata();
      tp = slice->ptimepoints();

      if ( d->empty() )
	{
	  logger << "  no data extracted\n";
	  return;
	}

      if ( d->size() != tp->size() )
	Helper::halt( "ACTIG: internal error: data/timepoint size mismatch" );
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
  const int expected_spe = need_signal ? (int)( Fs * epoch_sec + 1e-9 ) : 1;
  const int min_spe      = need_signal ? std::max( 1 , (int)( gap_min_pct / 100.0 * expected_spe ) ) : 1;

  std::vector<double> epochs( n_epochs , 0.0 );
  std::vector<int>    epoch_count( n_epochs , 0 );
  std::vector<bool>   is_gap( n_epochs , need_signal );   // prescored mode: no signal gaps

  int n_valid_epochs = 0 , n_gap_epochs = 0;
  if ( need_signal )
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
      n_valid_epochs = n_epochs;
      n_gap_epochs = 0;
      logger << "  prescored mode: skipping signal-based gap detection and NP signal extraction\n";
    }

  logger << "  recording: " << total_sec << " sec ("
	 << total_sec / 3600.0 << " hrs), "
	 << n_epochs << " epoch bins ("
	 << n_valid_epochs << " valid, "
	 << n_gap_epochs   << " gap)\n";


  // ---------------------------------------------------------------
  //
  // Nonparametric circadian metrics
  //
  // ---------------------------------------------------------------

  //
  // Aggregate 1-min epochs into NP bins (default 60 min)
  // A NP bin is valid if >= gap-min-pct of its constituent epochs are valid.
  //

  const int epochs_per_bin = (int)( np_bin_min * 60.0 / epoch_sec );

  if ( epochs_per_bin <= 0 )
    Helper::halt( "ACTIG: NP bin size must be >= epoch length" );

  const int n_bins = n_epochs / epochs_per_bin;

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
      for ( int e = b * epochs_per_bin ; e < (b+1) * epochs_per_bin ; e++ )
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

  double grand_mean = 0;
  int    gm_n = 0;
  for ( int i = 0 ; i < n_bins ; i++ )
    if ( !is_gap_bin[i] ) { grand_mean += bins[i]; gm_n++; }
  if ( gm_n > 0 ) grand_mean /= gm_n;

  double total_var = 0;
  for ( int i = 0 ; i < n_bins ; i++ )
    if ( !is_gap_bin[i] )
      { double diff = bins[i] - grand_mean; total_var += diff * diff; }


  //
  // IS (Interdaily Stability)
  //
  //  Requires at least 2 days of valid bin coverage.
  //  Slot means computed from valid bins only.
  //

  const int bins_per_day = (int)( 1440.0 / np_bin_min );
  if ( bins_per_day <= 0 )
    Helper::halt( "ACTIG: bin must be <= 1440 minutes" );
  const int n_full_days = n_bins / bins_per_day;

  double IS = 0;

  if ( n_valid_bins >= 2 * bins_per_day && total_var > 0 )
    {
      std::vector<double> slot_mean( bins_per_day , 0 );
      std::vector<int>    slot_n( bins_per_day , 0 );

      for ( int i = 0 ; i < n_bins ; i++ )
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

      IS = ( (double)n_valid_bins / bins_per_day ) * between_var / total_var;
    }


  //
  // IV (Intradaily Variability)
  //
  //  Cross-gap bin pairs are excluded from the numerator sum.
  //  Denominator uses valid bin count and valid-only total variance.
  //

  double IV = 0;

  if ( n_valid_bins >= 2 && total_var > 0 )
    {
      double diff_var = 0;
      int    diff_n   = 0;

      for ( int i = 1 ; i < n_bins ; i++ )
	if ( !is_gap_bin[i] && !is_gap_bin[i-1] )
	  {
	    double diff = bins[i] - bins[i-1];
	    diff_var += diff * diff;
	    diff_n++;
	  }

      if ( diff_n > 0 )
	IV = ( (double)n_valid_bins * diff_var ) / ( (double)diff_n * total_var );
    }


  //
  // L5 / M10 and RA
  //
  //  Profile means computed from valid bins only.
  //  Slots with no valid data remain 0 (noted via NP_NE_GAP).
  //

  double L5_val = 0 , M10_val = 0 , RA = 0;
  double L5_onset_hr = 0.0 , M10_onset_hr = 0.0;

  if ( bins_per_day == 24 && n_valid_bins >= bins_per_day )
    {
      std::vector<double> profile( 24 , 0 );
      std::vector<int>    profile_n( 24 , 0 );

      for ( int i = 0 ; i < n_bins ; i++ )
	if ( !is_gap_bin[i] )
	  {
	    int h = i % 24;
	    profile[h] += bins[i];
	    profile_n[h]++;
	  }
      for ( int h = 0 ; h < 24 ; h++ )
	if ( profile_n[h] > 0 ) profile[h] /= profile_n[h];

      find_window_fine( profile , l_hours , true  , &L5_val  , &L5_onset_hr  , np_step_min );
      find_window_fine( profile , m_hours , false , &M10_val , &M10_onset_hr , np_step_min );

      double denom = M10_val + L5_val;
      RA = denom > 0 ? ( M10_val - L5_val ) / denom : 0;
    }


  //
  // Clock-time labels and numeric onset values for L5/M10
  //
  // L5_onset_hr / M10_onset_hr are fractional hours from the profile
  // origin (= recording start, rounded to the hour).  Add start_hour
  // to convert to a 24 h clock position, then extract HH:MM:SS and
  // decimal minutes-from-midnight for downstream analysis.
  //

  clocktime_t startdatetime( edf.header.startdate , edf.header.starttime );

  double start_tod_sec = startdatetime.valid
    ? ( startdatetime.h * 3600.0 + startdatetime.m * 60.0 + startdatetime.s )
    : 0 ;

  int start_hour = (int)( start_tod_sec / 3600.0 );

  // helper: fractional profile-hours → HH:MM:SS string and minutes-from-midnight
  auto onset_to_hms_and_min = [&]( double onset_hr ,
				   std::string * hms_str ,
				   double      * min_from_midnight )
  {
    double clock_hr = std::fmod( start_hour + onset_hr , 24.0 );
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

  writer.value( "NP_IS" , IS );
  writer.value( "NP_IV" , IV );

  writer.value( "NP_L"  + Helper::int2str( l_hours ) , L5_val );
  writer.value( "NP_M"  + Helper::int2str( m_hours ) , M10_val );
  writer.value( "NP_RA" , RA );

  writer.value( "NP_L" + Helper::int2str( l_hours ) + "_ONSET"     , L5_hms );
  writer.value( "NP_M" + Helper::int2str( m_hours ) + "_ONSET"     , M10_hms );
  writer.value( "NP_L" + Helper::int2str( l_hours ) + "_ONSET_MIN" , L5_min_midnight );
  writer.value( "NP_M" + Helper::int2str( m_hours ) + "_ONSET_MIN" , M10_min_midnight );

  writer.value( "NP_NE"        , n_epochs );
  writer.value( "NP_NE_VALID"  , n_valid_epochs );
  writer.value( "NP_NE_GAP"    , n_gap_epochs );
  writer.value( "NP_NBINS"     , n_bins );
  writer.value( "NP_NBINS_VALID" , n_valid_bins );
  writer.value( "NP_NDAYS"     , n_full_days );

  logger << "  NP metrics: IS=" << IS
	 << " IV=" << IV
	 << " L" << l_hours << "=" << L5_val  << " (" << L5_hms  << ", " << L5_min_midnight  << " min)"
	 << " M" << m_hours << "=" << M10_val << " (" << M10_hms << ", " << M10_min_midnight << " min)"
	 << " RA=" << RA
	 << "\n";
  if ( ! np_traditional )
    logger << "  NP L5/M10 onset resolution: " << np_step_min << " min"
	   << " (use np-traditional for classic 1-hour resolution)\n";
  logger << "  NP bins: " << n_bins << " total, "
	 << n_valid_bins << " valid, "
	 << (n_bins - n_valid_bins) << " gap\n";


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

  if ( do_prescored )
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
  else if ( ! do_score )
    {
      logger << "  score not requested; skipping wake/sleep annotations and TST outputs\n";
      return;
    }

  // Cole-Kripke weights are calibrated for 60-second epochs.
  else if ( ( method == "cole" || method == "threshold" ) &&
	    std::fabs( epoch_sec - 60.0 ) > 1e-6 )
    Helper::halt( "ACTIG: cole and threshold methods require epoch=60" );

  if ( do_score && ( method == "threshold" || method == "thresh" ) )
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

	  logger << "  threshold method: median of valid epochs = " << t << "\n";
	}

      for ( int i = 0 ; i < n_epochs ; i++ )
	if ( !is_gap[i] ) is_wake[i] = ( epochs[i] >= t );

    }
  else if ( do_score && method == "cole" )
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
  else if ( do_score && method == "luna" )
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

      logger << "  luna: smooth=" << luna_smooth_min << "m (n=" << n_smooth << ")"
	     << ", burst=" << luna_burst_min << "m (n=" << n_burst_w << ")"
	     << ", burst-z=" << luna_burst_z
	     << ", quiet-z=" << luna_quiet_z
	     << ", active-frac=" << luna_active_frac << "\n";
      logger << "  luna: min-sleep=" << luna_min_sleep_min << "m (n=" << n_min_sleep << ")"
	     << ", max-gap=" << luna_max_gap_min << "m (n=" << n_gap_fill << ")"
	     << ", min-wake=" << luna_min_wake_min << "m (n=" << n_wake_break << ")\n";

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

  else if ( do_score )
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
  annot_t * a_gap_a = NULL;
  if ( do_score )
    {
      a_wake  = edf.timeline.annotations->add( wake_label );
      a_sleep = edf.timeline.annotations->add( sleep_label );
      a_gap_a = edf.timeline.annotations->add( gap_label );

      a_wake->description  = "Actigraphy-scored wake";
      a_sleep->description = "Actigraphy-scored sleep";
      a_gap_a->description = "Actigraphy gap / non-wear";
      a_wake->file = a_sleep->file = a_gap_a->file = "ACTIG";
    }

  auto epoch_state = [&]( int i ) -> int {
    if ( is_gap[i]  ) return 0;
    if ( is_wake[i] ) return 1;
    return 2;
  };

  int wake_inst = 1 , sleep_inst = 1 , gap_inst = 1;
  int total_wake = 0 , total_sleep = 0 , total_gap_epochs = 0;

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

	  if ( st == 0 )
	    {
	      if ( do_score ) a_gap_a->add( Helper::int2str( gap_inst ) , ival , "." );
	      ++gap_inst;
	      total_gap_epochs += run_len;
	    }
	  else if ( st == 1 )
	    {
	      if ( do_score ) a_wake->add( Helper::int2str( wake_inst ) , ival , "." );
	      ++wake_inst;
	      total_wake += run_len;
	    }
	  else
	    {
	      if ( do_score ) a_sleep->add( Helper::int2str( sleep_inst ) , ival , "." );
	      ++sleep_inst;
	      total_sleep += run_len;
	    }

	  run_start = i;
	}
    }


  //
  // Overall summary stats
  //  TST + WASO are computed over valid (non-gap) epochs only.
  //  sleep_pct denominator is valid epochs, not total epochs.
  //

  const double total_sleep_min     = total_sleep * epoch_sec / 60.0;
  const double total_wake_min      = total_wake  * epoch_sec / 60.0;
  const double total_gap_min       = total_gap_epochs * epoch_sec / 60.0;
  const double sleep_pct           = n_valid_epochs > 0
    ? 100.0 * total_sleep / (double)n_valid_epochs : 0;

  // SFI: sleep->wake transitions; skip transitions that cross a gap
  int trans_to_w = 0;
  for ( int i = 1 ; i < n_epochs ; i++ )
    if ( !is_gap[i] && !is_gap[i-1] && !is_wake[i-1] && is_wake[i] )
      ++trans_to_w;

  const double total_sleep_h = total_sleep * epoch_sec / 3600.0;
  const double sfi_luna_h    = total_sleep_h > 0 ? trans_to_w / total_sleep_h : 0;

  // SFI_ACT: sleep window = [first valid sleep epoch, last valid sleep epoch]
  int first_sleep_epoch = -1 , last_sleep_epoch = -1;
  for ( int i = 0 ; i < n_epochs ; i++ )
    if ( !is_gap[i] && !is_wake[i] )
      {
	if ( first_sleep_epoch == -1 ) first_sleep_epoch = i;
	last_sleep_epoch = i;
      }

  int    imm_bout_n = 0 , imm1_bout_n = 0;
  double mi_pct = 0 , imm1_pct = 0 , sfi_act = 0;

  if ( first_sleep_epoch != -1 )
    {
      // MI: wake % within sleep window, over valid epochs only
      int wake_in_win = 0 , valid_in_win = 0;
      for ( int i = first_sleep_epoch ; i <= last_sleep_epoch ; i++ )
	if ( !is_gap[i] ) { valid_in_win++; if ( is_wake[i] ) wake_in_win++; }
      mi_pct = valid_in_win > 0 ? 100.0 * wake_in_win / (double)valid_in_win : 0;

      // IMM bouts: consecutive sleep; gap breaks a bout
      int i = first_sleep_epoch;
      while ( i <= last_sleep_epoch )
	{
	  if ( is_gap[i] || is_wake[i] ) { ++i; continue; }
	  int j = i + 1;
	  while ( j <= last_sleep_epoch && !is_gap[j] && !is_wake[j] ) ++j;
	  const int run_epochs = j - i;
	  ++imm_bout_n;
	  if ( run_epochs <= 1 ) ++imm1_bout_n;
	  i = j;
	}

      imm1_pct = imm_bout_n > 0 ? 100.0 * imm1_bout_n / (double)imm_bout_n : 0;
      sfi_act  = mi_pct + imm1_pct;
    }

  writer.value( "SCORE_TST_MIN"     , total_sleep_min );
  writer.value( "SCORE_WASO_MIN"    , total_wake_min );
  writer.value( "SCORE_GAP_MIN"     , total_gap_min );
  writer.value( "SCORE_SLEEP_PCT"   , sleep_pct );
  writer.value( "FRAG_SFI"         , sfi_luna_h );
  writer.value( "FRAG_SFI_N"       , trans_to_w );
  writer.value( "FRAG_SFI_ACT"     , sfi_act );
  writer.value( "FRAG_MI_PCT"      , mi_pct );
  writer.value( "FRAG_IMM1_PCT"    , imm1_pct );
  writer.value( "FRAG_IMM_BOUT_N"  , imm_bout_n );
  writer.value( "FRAG_IMM1_BOUT_N" , imm1_bout_n );
  writer.value( "SCORE_METHOD"     , score_method );
  writer.value( "SCORE_WAKE_RUN_N" , wake_inst - 1 );
  writer.value( "SCORE_SLEEP_RUN_N", sleep_inst - 1 );
  writer.value( "SCORE_GAP_RUN_N"  , gap_inst - 1 );


  // -----------------------------------------------------------------------
  //
  // Per-day summary
  //
  //  Epochs are clock-time-anchored (bin e corresponds to recording_start
  //  + e * epoch_sec), so day partitioning is straightforward.
  //  Gap epochs are tracked separately; they do not contribute to TST/WASO.
  //  A day is INCLUDED in cross-day averages only if valid minutes >= day-min-valid.
  //
  // -----------------------------------------------------------------------

  double first_boundary_sec = total_sec;

  if ( startdatetime.valid )
    {
      const double anchor_sec = day_anchor_hr * 3600.0;
      first_boundary_sec = anchor_sec - start_tod_sec;
      if ( first_boundary_sec <= 0 )
	first_boundary_sec += 86400.0;
    }

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

  std::vector<int> day_epoch_n;        // valid (non-gap) epochs per day
  std::vector<int> day_gap_epoch_n;    // gap epochs per day
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
	day_gap_epoch_n.resize( d+1 , 0 );
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

      if ( is_gap[i] )
	{
	  // gap epochs: count separately, do not contribute to wake/sleep
	  day_gap_epoch_n[di]++;
	  continue;
	}

      const double stop_sec = std::min( (i+1) * epoch_sec , total_sec );

      day_epoch_n[di]++;
      if ( is_wake[i] ) day_wake_epoch_n[di]++;
      if ( day_first_idx[di] == -1 ) day_first_idx[di] = i;
      day_last_idx[di] = i;

      if ( is_wake[i] )
	accumulate_day_seconds( &wake_sec_day , start_sec , stop_sec );
      else
	accumulate_day_seconds( &sleep_sec_day , start_sec , stop_sec );
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
  int    n_included_days = 0;

  for ( size_t d = 0 ; d < n_days ; d++ )
    {
      const double sleep_sec    = d < sleep_sec_day.size()  ? sleep_sec_day[d]   : 0.0;
      const double wake_sec     = d < wake_sec_day.size()   ? wake_sec_day[d]    : 0.0;
      const int day_epochs      = d < day_epoch_n.size()    ? day_epoch_n[d]     : 0;
      const int day_gap_epochs  = d < day_gap_epoch_n.size()? day_gap_epoch_n[d] : 0;
      const int day_wake_epochs = d < day_wake_epoch_n.size()? day_wake_epoch_n[d]: 0;
      const int day_trans       = d < day_trans_to_w.size() ? day_trans_to_w[d]  : 0;

      const double day_sleep_h     = sleep_sec / 3600.0;
      const double day_sfi_luna_h  = day_sleep_h > 0 ? day_trans / day_sleep_h : 0.0;

      const double valid_min_day   = day_epochs * epoch_sec / 60.0;
      const double gap_min_day     = day_gap_epochs * epoch_sec / 60.0;
      const double total_day_min   = valid_min_day + gap_min_day;
      const double valid_pct_day   = total_day_min > 0
	? 100.0 * valid_min_day / total_day_min : 0.0;
      const double day_total_valid = sleep_sec + wake_sec;
      const double day_sleep_pct   = day_total_valid > 0
	? 100.0 * sleep_sec / day_total_valid : 0.0;

      const bool included = ( valid_min_day >= day_min_valid_min );

      // Per-day fragmentation (gap-aware, same logic as overall)
      int    day_imm_bout_n = 0 , day_imm1_bout_n = 0;
      double day_mi_pct = 0 , day_imm1_pct = 0 , day_sfi_act = 0;

      if ( d < day_first_idx.size() && day_first_idx[d] != -1 )
	{
	  const int b = day_first_idx[d];
	  const int e = day_last_idx[d];

	  int first_sleep_day = -1 , last_sleep_day = -1;
	  for ( int i = b ; i <= e ; i++ )
	    if ( !is_gap[i] && !is_wake[i] )
	      {
		if ( first_sleep_day == -1 ) first_sleep_day = i;
		last_sleep_day = i;
	      }

	  if ( first_sleep_day != -1 )
	    {
	      int wake_in_win = 0 , valid_in_win = 0;
	      for ( int i = first_sleep_day ; i <= last_sleep_day ; i++ )
		if ( !is_gap[i] ) { valid_in_win++; if ( is_wake[i] ) wake_in_win++; }
	      day_mi_pct = valid_in_win > 0
		? 100.0 * wake_in_win / (double)valid_in_win : 0.0;

	      int i = first_sleep_day;
	      while ( i <= last_sleep_day )
		{
		  if ( is_gap[i] || is_wake[i] ) { ++i; continue; }
		  int j = i + 1;
		  while ( j <= last_sleep_day && !is_gap[j] && !is_wake[j] ) ++j;
		  const int run_epochs = j - i;
		  ++day_imm_bout_n;
		  if ( run_epochs <= 1 ) ++day_imm1_bout_n;
		  i = j;
		}

	      day_imm1_pct = day_imm_bout_n > 0
		? 100.0 * day_imm1_bout_n / (double)day_imm_bout_n : 0.0;
	      day_sfi_act = day_mi_pct + day_imm1_pct;
	    }
	}

      const std::string day_name = "day" + Helper::int2str( (int)d + 1 );
      writer.level( day_name , "DAY" );
      writer.value( "SCORE_TST_MIN"     , sleep_sec / 60.0 );
      writer.value( "SCORE_WASO_MIN"    , wake_sec  / 60.0 );
      writer.value( "SCORE_GAP_MIN"     , gap_min_day );
      writer.value( "VALID_MIN"         , valid_min_day );
      writer.value( "VALID_PCT"         , valid_pct_day );
      writer.value( "INCLUDED"          , included ? 1 : 0 );
      writer.value( "SCORE_SLEEP_PCT"   , day_sleep_pct );
      writer.value( "FRAG_SFI"         , day_sfi_luna_h );
      writer.value( "FRAG_SFI_N"       , day_trans );
      writer.value( "FRAG_SFI_ACT"     , day_sfi_act );
      writer.value( "FRAG_MI_PCT"      , day_mi_pct );
      writer.value( "FRAG_IMM1_PCT"    , day_imm1_pct );
      writer.value( "FRAG_IMM_BOUT_N"  , day_imm_bout_n );
      writer.value( "FRAG_IMM1_BOUT_N" , day_imm1_bout_n );
      writer.value( "SCORE_EPOCH_N"    , day_epochs );
      writer.value( "SCORE_WAKE_EPOCH_N", day_wake_epochs );
      writer.unlevel( "DAY" );

      if ( included )
	{
	  sum_sleep_sec += sleep_sec;
	  sum_wake_sec  += wake_sec;
	  n_included_days++;
	}
    }

  if ( n_days > 0 )
    {
      writer.value( "DAY_N"          , (int)n_days );
      writer.value( "DAY_N_INCLUDED" , n_included_days );
      writer.value( "DAY_N_EXCLUDED" , (int)n_days - n_included_days );

      if ( n_included_days > 0 )
	{
	  const double total_scored_sec = sum_sleep_sec + sum_wake_sec;
	  writer.value( "SCORE_TST_DAYAVG_MIN"   , ( sum_sleep_sec / n_included_days ) / 60.0 );
	  writer.value( "SCORE_WASO_DAYAVG_MIN"  , ( sum_wake_sec  / n_included_days ) / 60.0 );
	  writer.value( "SCORE_SLEEP_PCT_DAYAVG" ,
			total_scored_sec > 0
			? 100.0 * sum_sleep_sec / total_scored_sec : 0.0 );
	}
    }

  if ( n_days > 0 )
    logger << "  per-day output: " << n_days << " day(s), "
	   << n_included_days << " included (>= "
	   << day_min_valid_min << " valid min)"
	   << ( startdatetime.valid ? "" : " [no valid start datetime; 24h from recording start]" )
	   << "\n";

  logger << "  scored " << n_epochs << " epoch bins: "
	 << total_sleep << " sleep (" << total_sleep_min << " min), "
	 << total_wake  << " wake ("  << total_wake_min  << " min), "
	 << total_gap_epochs << " gap (" << total_gap_min << " min), "
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
      if ( Helper::str2int( debt_target_str , &trial_int ) )
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
	    Helper::halt( "ACTIG: debt-target date " + debt_target_str
			  + " maps to day " + Helper::int2str( target_idx + 1 )
			  + " which is outside the recording (1.."
			  + Helper::int2str( (int)n_days ) + ")" );
	}

      logger << "  debt: target=day" << ( target_idx + 1 )
	     << ", recent=" << debt_R
	     << ", base=" << debt_B
	     << ", min-base=" << debt_min_base
	     << ", min-z=" << debt_min_z
	     << ", w(TST)=" << debt_w << "\n";


      // ------------------------------------------------------------------
      // Per-day helpers  (INCLUDED = valid minutes >= day-min-valid)
      // ------------------------------------------------------------------

      auto day_is_included = [&]( int d ) -> bool
      {
	if ( d < 0 || d >= (int)n_days ) return false;
	const double vm = ( d < (int)day_epoch_n.size() ? day_epoch_n[d] : 0 )
	  * epoch_sec / 60.0;
	return vm >= day_min_valid_min;
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

    }  // end do_debt

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

  const std::string day_prefix = param.has( "prefix" ) ? param.value( "prefix" ) : "day";

  const bool add_hours = param.has( "hours" );

  const std::string hr_prefix = param.has( "hour-prefix" ) ? param.value( "hour-prefix" ) : "";

  const bool add_halves = param.has( "halves" );

  const bool verbose = param.has( "verbose" );


  //
  // Compute first day boundary
  //

  double start_tod_sec = startdatetime.h * 3600.0
    + startdatetime.m * 60.0
    + startdatetime.s;

  double anchor_sec = anchor_hr * 3600.0;

  double first_boundary_sec = anchor_sec - start_tod_sec;
  if ( first_boundary_sec <= 0 )
    first_boundary_sec += 86400.0;


  //
  // Build day annotations
  //

  int day_num = 1;

  double day_start_sec = 0;
  double day_stop_sec  = first_boundary_sec;

  if ( day_stop_sec > total_sec )
    day_stop_sec = total_sec;

  while ( day_start_sec < total_sec )
    {

      const std::string name = day_prefix + Helper::int2str( day_num );

      uint64_t tp_start = Helper::sec2tp( day_start_sec );
      uint64_t tp_stop  = Helper::sec2tp( day_stop_sec );

      annot_t * a = edf.timeline.annotations->add( name );
      a->description = name;
      a->file = "DAYS";

      a->add( "." , interval_t( tp_start , tp_stop ) , "." );

      if ( verbose )
	{
	  clocktime_t ds = startdatetime;
	  ds.advance_seconds( day_start_sec );

	  clocktime_t de = startdatetime;
	  de.advance_seconds( day_stop_sec );

	  writer.level( name , "DAY" );
	  writer.value( "START_SEC" , day_start_sec );
	  writer.value( "STOP_SEC" , day_stop_sec );
	  writer.value( "DUR_SEC" , day_stop_sec - day_start_sec );
	  writer.value( "START_HMS" , ds.as_string( ':' ) );
	  writer.value( "STOP_HMS" , de.as_string( ':' ) );
	  writer.value( "START_DATE" , ds.as_date_string( '/' ) );
	}

      logger << "  " << name
	     << " : " << day_start_sec << " - " << day_stop_sec
	     << " sec (" << (day_stop_sec - day_start_sec) / 3600.0 << " hrs)\n";

      ++day_num;
      day_start_sec = day_stop_sec;
      day_stop_sec  = day_start_sec + 86400.0;

      if ( day_stop_sec > total_sec )
	day_stop_sec = total_sec;
    }

  if ( verbose )
    writer.unlevel( "DAY" );

  logger << "  DAYS: created " << (day_num - 1) << " day annotations";
  if ( anchor_hr != 12 )
    logger << " (anchor=" << anchor_hr << "h)";
  logger << "\n";


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

}
