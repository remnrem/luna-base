
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
//  Organizes output around four conceptual domains:
//    A. MIXEDNESS  - how diffuse is the posterior at each sample
//    B. INSTABILITY - how much does the posterior move over time
//    C. TRANSITION STRUCTURE - what happens near likely state boundaries
//    D. CONTEXT - how A/B/C vary across user-defined annotation strata
//
//  The primary scientific distinction preserved throughout is:
//    stable mixed/intermediate states  vs.  transitional instability
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

extern logger_t logger;
extern writer_t writer;


// ============================================================
// Internal types
// ============================================================

enum hd_trans_method_t { HD_HARD, HD_MOTION, HD_BOTH };

struct hd_params_t
{
  std::vector<std::string> ch;      // channel names [W, N1, N2, N3, R]
  bool do_3state;
  hd_trans_method_t method;
  double motion_th;                 // TV threshold for motion-based detection
  double window_sec;                // transition window half-width (seconds)
  double lag_sec;                   // lag for TV_lag (seconds)
  double stable_min_sec;            // min seconds from any event to qualify as stable-core
  double conf_th;                   // confidence threshold for frac_below metric
  std::string annot_name;           // annotation class for stratification (empty = none)
  bool verbose;                     // emit per-sample HDSIG table
  int min_events_profile;           // minimum events needed to emit aligned profile
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
  std::vector<uint64_t> tp;            // timepoints (in Luna units)

  bool load( edf_t & edf, const hd_params_t & par )
  {
    const int K = 5;

    // Resolve channels and verify all exist
    std::vector<int> sidx( K );
    for ( int k = 0; k < K; k++ )
      {
	signal_list_t sl = edf.header.signal_list( par.ch[k] );
	if ( sl.size() == 0 )
	  Helper::halt( "HDSTATS: channel not found: " + par.ch[k] );
	if ( edf.header.is_annotation_channel( sl(0) ) )
	  Helper::halt( "HDSTATS: " + par.ch[k] + " is an annotation channel, not a signal" );
	sidx[k] = sl(0);
      }

    // Verify same sample rate for all channels
    Fs = edf.header.sampling_freq( sidx[0] );
    for ( int k = 1; k < K; k++ )
      {
	double fk = edf.header.sampling_freq( sidx[k] );
	if ( std::fabs( fk - Fs ) > 1e-6 )
	  Helper::halt( "HDSTATS: all hypnodensity channels must share the same sample rate; "
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
	  Helper::halt( "HDSTATS: channel length mismatch for " + par.ch[k] );
      }

    if ( N == 0 )
      Helper::halt( "HDSTATS: no samples found" );

    // Validate and optionally renormalize row sums
    int n_renorm = 0;
    for ( int i = 0; i < N; i++ )
      {
	double rowsum = 0.0;
	for ( int k = 0; k < K; k++ ) rowsum += p[k][i];

	if ( rowsum < 0.5 || rowsum > 1.5 )
	  Helper::halt( "HDSTATS: row sum " + Helper::dbl2str(rowsum)
			+ " at sample " + Helper::int2str(i) + " is implausible" );

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
      logger << "  HDSTATS: renormalized " << n_renorm << " of " << N << " samples\n";

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
  std::vector<double> TV, L2, TV_lag;    // motion metrics
  std::vector<double> mix_wn1;           // min(W, N1)
  std::vector<double> mix_n2n3;          // min(N2, N3)
  std::vector<double> mix_rn1;           // min(R, N1)
  std::vector<int>    argmax5;

  // 3-state signals (only populated if do_3state)
  std::vector<double> H3, C3, Mg3;
  std::vector<double> TV3, TV3_lag;
  std::vector<double> mix3_wnrem;        // min(W, NREM)
  std::vector<double> mix3_nremr;        // min(NREM, R)
  std::vector<int>    argmax3;

  void compute( const hd_data_t & dat, const hd_params_t & par )
  {
    const int N = dat.N;
    const int lag = std::max( 1, (int)std::round( par.lag_sec * dat.Fs ) );

    H.assign( N, 0.0 );
    C.assign( N, 0.0 );
    Mg.assign( N, 0.0 );
    TV.assign( N, 0.0 );
    L2.assign( N, 0.0 );
    TV_lag.assign( N, 0.0 );
    mix_wn1.assign( N, 0.0 );
    mix_n2n3.assign( N, 0.0 );
    mix_rn1.assign( N, 0.0 );
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

	// Pairwise mixing
	mix_wn1[i]  = std::min( dat.p[0][i], dat.p[1][i] );
	mix_n2n3[i] = std::min( dat.p[2][i], dat.p[3][i] );
	mix_rn1[i]  = std::min( dat.p[4][i], dat.p[1][i] );

	// Motion metrics (step from previous sample)
	if ( i > 0 )
	  {
	    double tv = 0.0, l2sq = 0.0;
	    for ( int k = 0; k < 5; k++ )
	      {
		double d = dat.p[k][i] - dat.p[k][i-1];
		tv   += std::fabs( d );
		l2sq += d * d;
	      }
	    TV[i] = 0.5 * tv;
	    L2[i] = std::sqrt( l2sq );
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
    mix3_wnrem.assign( N, 0.0 );
    mix3_nremr.assign( N, 0.0 );
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

	mix3_wnrem[i] = std::min( dat.p3[0][i], dat.p3[1][i] );
	mix3_nremr[i] = std::min( dat.p3[1][i], dat.p3[2][i] );

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
  int    idx;         // sample index of event center
  int    from_st;     // argmax state just before (hard method; -1 if motion-only)
  int    to_st;       // argmax state just after (hard method; -1 if motion-only)
  double peak_tv;     // TV value at event (motion method; NaN if hard-only)
};

// Detect transitions by hard argmax changes
static std::vector<hd_event_t> detect_hard(
    const std::vector<int> & argmax,
    const std::vector<double> & tv,
    int N )
{
  std::vector<hd_event_t> events;
  for ( int i = 1; i < N; i++ )
    {
      if ( argmax[i] != argmax[i-1] )
	{
	  hd_event_t e;
	  e.idx     = i;
	  e.from_st = argmax[i-1];
	  e.to_st   = argmax[i];
	  e.peak_tv = tv[i];
	  events.push_back( e );
	}
    }
  return events;
}

// Detect transitions by TV peaks above threshold, then merge events within min_gap samples
static std::vector<hd_event_t> detect_motion(
    const std::vector<double> & tv,
    const std::vector<int> & argmax,
    int N,
    double threshold,
    int min_gap )
{
  // Find local peaks above threshold
  std::vector<int> peaks;
  for ( int i = 1; i < N - 1; i++ )
    {
      if ( tv[i] >= threshold &&
	   tv[i] >= tv[i-1] &&
	   tv[i] >= tv[i+1] )
	peaks.push_back( i );
    }
  // Handle edge: last sample can be a peak
  if ( N > 1 && tv[N-1] >= threshold && tv[N-1] >= tv[N-2] )
    peaks.push_back( N - 1 );

  // Merge peaks within min_gap: keep the higher-TV peak in each cluster
  std::vector<hd_event_t> events;
  int prev = -1;
  int prev_peak_idx = -1;
  double prev_peak_tv = 0.0;

  for ( int pi = 0; pi < (int)peaks.size(); pi++ )
    {
      int cur = peaks[pi];
      double cur_tv = tv[cur];

      if ( prev_peak_idx >= 0 && ( cur - prev_peak_idx ) < min_gap )
	{
	  // Same cluster: keep higher-TV peak
	  if ( cur_tv > prev_peak_tv )
	    {
	      prev_peak_idx = cur;
	      prev_peak_tv  = cur_tv;
	    }
	}
      else
	{
	  // Emit the previous cluster's representative
	  if ( prev_peak_idx >= 0 )
	    {
	      hd_event_t e;
	      e.idx     = prev_peak_idx;
	      e.from_st = ( prev_peak_idx > 0 ) ? argmax[prev_peak_idx - 1] : -1;
	      e.to_st   = argmax[prev_peak_idx];
	      e.peak_tv = prev_peak_tv;
	      events.push_back( e );
	    }
	  prev_peak_idx = cur;
	  prev_peak_tv  = cur_tv;
	}
    }
  // Emit the last cluster
  if ( prev_peak_idx >= 0 )
    {
      hd_event_t e;
      e.idx     = prev_peak_idx;
      e.from_st = ( prev_peak_idx > 0 ) ? argmax[prev_peak_idx - 1] : -1;
      e.to_st   = argmax[prev_peak_idx];
      e.peak_tv = prev_peak_tv;
      events.push_back( e );
    }

  return events;
}

// Given an event list and window/stable sizes, compute per-sample boolean masks
static void build_masks(
    const std::vector<hd_event_t> & events,
    int N,
    int window_samples,
    int stable_min_samples,
    std::vector<bool> & is_trans,
    std::vector<bool> & is_stable )
{
  is_trans.assign( N, false );
  is_stable.assign( N, true );  // initially all stable; will clear near events

  for ( const auto & e : events )
    {
      int lo_t = std::max( 0,     e.idx - window_samples );
      int hi_t = std::min( N - 1, e.idx + window_samples );
      int lo_s = std::max( 0,     e.idx - stable_min_samples );
      int hi_s = std::min( N - 1, e.idx + stable_min_samples );

      for ( int i = lo_t; i <= hi_t; i++ ) is_trans[i]  = true;
      for ( int i = lo_s; i <= hi_s; i++ ) is_stable[i] = false;
    }
}

struct hd_trans_t
{
  std::vector<hd_event_t> events5;
  std::vector<bool>       is_trans5;
  std::vector<bool>       is_stable5;

  std::vector<hd_event_t> events3;
  std::vector<bool>       is_trans3;
  std::vector<bool>       is_stable3;

  void detect( const hd_data_t & dat,
	       const hd_derived_t & der,
	       const hd_params_t & par )
  {
    const int N = dat.N;
    const int win_smp    = std::max( 1, (int)std::round( par.window_sec     * dat.Fs ) );
    const int stable_smp = std::max( 1, (int)std::round( par.stable_min_sec * dat.Fs ) );
    const int min_gap    = win_smp;  // merge motion events within one window width

    // 5-state detection
    std::vector<hd_event_t> hard5, motion5;
    if ( par.method == HD_HARD || par.method == HD_BOTH )
      hard5 = detect_hard( der.argmax5, der.TV, N );
    if ( par.method == HD_MOTION || par.method == HD_BOTH )
      motion5 = detect_motion( der.TV, der.argmax5, N, par.motion_th, min_gap );

    if ( par.method == HD_HARD )        events5 = hard5;
    else if ( par.method == HD_MOTION ) events5 = motion5;
    else
      {
	// Merge hard + motion events, sort by index, re-merge within min_gap
	events5 = hard5;
	events5.insert( events5.end(), motion5.begin(), motion5.end() );
	std::sort( events5.begin(), events5.end(),
		   [](const hd_event_t & a, const hd_event_t & b){ return a.idx < b.idx; } );
	// Deduplicate / merge nearby
	std::vector<hd_event_t> merged;
	int prev_idx = -999999;
	double prev_tv = 0.0;
	for ( int ei = 0; ei < (int)events5.size(); ei++ )
	  {
	    if ( events5[ei].idx - prev_idx < min_gap )
	      {
		if ( events5[ei].peak_tv > merged.back().peak_tv )
		  merged.back() = events5[ei];
	      }
	    else
	      {
		merged.push_back( events5[ei] );
		prev_idx = events5[ei].idx;
		prev_tv  = events5[ei].peak_tv;
	      }
	  }
	events5 = merged;
      }

    build_masks( events5, N, win_smp, stable_smp, is_trans5, is_stable5 );

    if ( ! par.do_3state ) return;

    // 3-state detection (independent from 5-state)
    std::vector<hd_event_t> hard3, motion3;
    if ( par.method == HD_HARD || par.method == HD_BOTH )
      hard3 = detect_hard( der.argmax3, der.TV3, N );
    if ( par.method == HD_MOTION || par.method == HD_BOTH )
      motion3 = detect_motion( der.TV3, der.argmax3, N, par.motion_th, min_gap );

    if ( par.method == HD_HARD )        events3 = hard3;
    else if ( par.method == HD_MOTION ) events3 = motion3;
    else
      {
	events3 = hard3;
	events3.insert( events3.end(), motion3.begin(), motion3.end() );
	std::sort( events3.begin(), events3.end(),
		   [](const hd_event_t & a, const hd_event_t & b){ return a.idx < b.idx; } );
	std::vector<hd_event_t> merged;
	for ( int ei = 0; ei < (int)events3.size(); ei++ )
	  {
	    if ( !merged.empty() && events3[ei].idx - merged.back().idx < min_gap )
	      {
		if ( events3[ei].peak_tv > merged.back().peak_tv )
		  merged.back() = events3[ei];
	      }
	    else
	      merged.push_back( events3[ei] );
	  }
	events3 = merged;
      }

    build_masks( events3, N, win_smp, stable_smp, is_trans3, is_stable3 );
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

struct hd_region_stats_t
{
  int n = 0;
  bool valid = false;
  // Domain A: mixedness
  double mean_H = 0, sd_H = 0, p90_H = 0;
  double mean_C = 0, frac_C_below = 0, mean_Mg = 0;
  // Domain B: instability
  double mean_TV = 0, sd_TV = 0, p90_TV = 0, mean_TV_lag = 0;
  double corr_H_TV = 0;
  // Domain C: pairwise mixing (a=WN1/WNREM, b=N2N3, c=RN1/NREMR)
  double mean_mix_a = 0, mean_mix_b = 0, mean_mix_c = 0;
};

static hd_region_stats_t compute_region(
    const hd_derived_t & der,
    const std::vector<bool> & mask,   // which samples to include
    const hd_params_t & par,
    bool threestate )
{
  hd_region_stats_t rs;

  std::vector<double> vH, vC, vMg, vTV, vTV_lag, vMix_a, vMix_b, vMix_c;

  const std::vector<double> & H_ref       = threestate ? der.H3      : der.H;
  const std::vector<double> & C_ref       = threestate ? der.C3      : der.C;
  const std::vector<double> & Mg_ref      = threestate ? der.Mg3     : der.Mg;
  const std::vector<double> & TV_ref      = threestate ? der.TV3     : der.TV;
  const std::vector<double> & TVlag_ref   = threestate ? der.TV3_lag : der.TV_lag;
  const std::vector<double> & mix_a_ref   = threestate ? der.mix3_wnrem : der.mix_wn1;
  const std::vector<double> & mix_c_ref   = threestate ? der.mix3_nremr : der.mix_rn1;
  // mix_b is N2/N3 mixing — only meaningful in 5-state
  const std::vector<double> & mix_b_ref   = der.mix_n2n3;

  int N = (int)mask.size();
  for ( int i = 0; i < N; i++ )
    {
      if ( !mask[i] ) continue;
      vH.push_back( H_ref[i] );
      vC.push_back( C_ref[i] );
      vMg.push_back( Mg_ref[i] );
      vTV.push_back( TV_ref[i] );
      vTV_lag.push_back( TVlag_ref[i] );
      vMix_a.push_back( mix_a_ref[i] );
      if ( !threestate ) vMix_b.push_back( mix_b_ref[i] );
      vMix_c.push_back( mix_c_ref[i] );
    }

  rs.n = (int)vH.size();
  if ( rs.n == 0 ) return rs;
  rs.valid = true;

  // Means
  auto vmean = [](const std::vector<double> & v) {
    if (v.empty()) return 0.0;
    return std::accumulate( v.begin(), v.end(), 0.0 ) / v.size();
  };
  auto vsd = [&vmean](const std::vector<double> & v) {
    if (v.size() < 2) return 0.0;
    double m = vmean(v);
    double s = 0.0;
    for (auto x : v) s += (x-m)*(x-m);
    return std::sqrt( s / (v.size()-1) );
  };

  rs.mean_H      = vmean( vH );
  rs.sd_H        = vsd( vH );
  rs.p90_H       = percentile( vH, 0.9 );
  rs.mean_C      = vmean( vC );
  rs.mean_Mg     = vmean( vMg );
  rs.mean_TV     = vmean( vTV );
  rs.sd_TV       = vsd( vTV );
  rs.p90_TV      = percentile( vTV, 0.9 );
  rs.mean_TV_lag = vmean( vTV_lag );
  rs.mean_mix_a  = vmean( vMix_a );
  rs.mean_mix_b  = threestate ? std::numeric_limits<double>::quiet_NaN() : vmean( vMix_b );
  rs.mean_mix_c  = vmean( vMix_c );

  // Fraction with confidence below threshold
  int n_below = 0;
  for ( double v : vC ) if ( v < par.conf_th ) ++n_below;
  rs.frac_C_below = (double)n_below / rs.n;

  // Entropy-TV correlation
  rs.corr_H_TV = pearson( vH, vTV );

  return rs;
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
    const hd_data_t & dat,
    const hd_params_t & par,
    double Fs,
    bool threestate,
    hd_trans_stats_t & tshape,
    hd_profile_t & profile )
{
  const int N = (int)region_mask.size();
  const int win_smp = std::max( 1, (int)std::round( par.window_sec * Fs ) );

  // Filter events to those whose center sample falls within the region
  std::vector<const hd_event_t *> local_events;
  for ( const auto & e : all_events )
    if ( e.idx >= 0 && e.idx < N && region_mask[e.idx] )
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

  // Baseline H for width estimation: mean H in stable-core of this region
  // (We use a simple approximation: mean H in region samples that are not transition-zone)
  // This is just the p50 of H across the region for convenience
  std::vector<double> all_H_region;
  for ( int i = 0; i < N; i++ )
    if ( region_mask[i] ) all_H_region.push_back( H_ref[i] );
  double baseline_H = percentile( all_H_region, 0.5 );

  // Entropy threshold for width measurement: halfway between baseline and max
  // We'll also compute width based on TV threshold
  double H_width_th = baseline_H + 0.5 * ( percentile( all_H_region, 0.9 ) - baseline_H );
  double TV_width_th = par.motion_th * 0.5;

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

static void write_region( const hd_region_stats_t & rs, bool threestate )
{
  if ( !rs.valid ) return;
  writer.value( "N",          rs.n );
  writer.value( "MEAN_H",     rs.mean_H );
  writer.value( "SD_H",       rs.sd_H );
  writer.value( "P90_H",      rs.p90_H );
  writer.value( "MEAN_C",     rs.mean_C );
  writer.value( "FRAC_C_LT",  rs.frac_C_below );
  writer.value( "MEAN_MG",    rs.mean_Mg );
  writer.value( "MEAN_TV",    rs.mean_TV );
  writer.value( "SD_TV",      rs.sd_TV );
  writer.value( "P90_TV",     rs.p90_TV );
  writer.value( "MEAN_TV_LAG", rs.mean_TV_lag );
  writer.value( "CORR_H_TV",  rs.corr_H_TV );
  writer.value( "MEAN_MIX_A", rs.mean_mix_a );
  if ( !threestate )
    writer.value( "MEAN_MIX_B", rs.mean_mix_b );
  writer.value( "MEAN_MIX_C", rs.mean_mix_c );
}

static void write_trans_stats( const hd_trans_stats_t & ts )
{
  writer.value( "N_TRANS",      ts.n_events );
  writer.value( "TRANS_DENS",   ts.density );
  if ( ts.valid )
    {
      writer.value( "MEAN_TRANS_W",    ts.mean_width_sec );
      writer.value( "MEAN_PEAK_H",     ts.mean_peak_H );
      writer.value( "MEAN_MIN_C",      ts.mean_min_C );
      writer.value( "MEAN_TV_AREA",    ts.mean_TV_area );
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

// Write HDSTATS rows for one context (global or per-stratum).
// Emits REGION strata: ALL, STABLE, TRANS.
// Emits STATE stratum only when do_3state is true.
static void write_hdstats(
    const hd_derived_t & der,
    const hd_trans_t & tr,
    const std::vector<bool> & region_mask5,   // which samples are in this context (5-state)
    const std::vector<bool> & region_mask3,   // same for 3-state (may be same as region_mask5)
    const hd_params_t & par,
    const hd_data_t & dat )
{
  auto make_all   = [&]( const std::vector<bool> & rm ) { return rm; };
  auto make_stable = [&]( const std::vector<bool> & rm,
			   const std::vector<bool> & is_stable ) {
    std::vector<bool> m( rm.size() );
    for ( int i = 0; i < (int)rm.size(); i++ ) m[i] = rm[i] && is_stable[i];
    return m;
  };
  auto make_trans = [&]( const std::vector<bool> & rm,
			  const std::vector<bool> & is_trans ) {
    std::vector<bool> m( rm.size() );
    for ( int i = 0; i < (int)rm.size(); i++ ) m[i] = rm[i] && is_trans[i];
    return m;
  };

  // Helper: emit one state (5 or 3)
  auto emit_state = [&]( bool threestate ) {

    const std::vector<bool> & rmask     = threestate ? region_mask3 : region_mask5;
    const std::vector<bool> & is_trans  = threestate ? tr.is_trans3  : tr.is_trans5;
    const std::vector<bool> & is_stable = threestate ? tr.is_stable3 : tr.is_stable5;

    // HDSTATS table (REGION strata)
    {
      hd_region_stats_t rs_all    = compute_region( der, make_all(rmask),                par, threestate );
      hd_region_stats_t rs_stable = compute_region( der, make_stable(rmask, is_stable),  par, threestate );
      hd_region_stats_t rs_trans  = compute_region( der, make_trans(rmask,  is_trans),   par, threestate );

      writer.level( "ALL",   "REGION" );
      write_region( rs_all, threestate );
      writer.unlevel( "REGION" );

      writer.level( "STABLE", "REGION" );
      write_region( rs_stable, threestate );
      writer.unlevel( "REGION" );

      writer.level( "TRANS", "REGION" );
      write_region( rs_trans, threestate );
      writer.unlevel( "REGION" );

      // Stable-vs-transition ratios (no REGION stratum)
      if ( rs_all.valid )
	{
	  if ( rs_stable.valid && rs_trans.valid )
	    {
	      writer.value( "H_RATIO_TR_ST",  rs_stable.mean_H > 0
			    ? rs_trans.mean_H / rs_stable.mean_H : std::numeric_limits<double>::quiet_NaN() );
	      writer.value( "CONF_DIFF_TR_ST", rs_trans.mean_C - rs_stable.mean_C );
	      writer.value( "TV_RATIO_TR_ST",  rs_stable.mean_TV > 0
			    ? rs_trans.mean_TV / rs_stable.mean_TV : std::numeric_limits<double>::quiet_NaN() );
	    }
	}
    }

    // HDTRANS table (transition shape)
    {
      const std::vector<hd_event_t> & evts = threestate ? tr.events3 : tr.events5;
      hd_trans_stats_t tshape;
      hd_profile_t     profile;
      compute_trans_shape( der, evts, rmask, dat, par, dat.Fs, threestate, tshape, profile );
      write_trans_stats( tshape );

      // HDTRANS_PROFILE table
      if ( profile.valid )
	{
	  for ( int j = 0; j < (int)profile.offsets.size(); j++ )
	    {
	      writer.level( profile.offsets[j], "OFFSET" );
	      writer.value( "H",  profile.H[j] );
	      writer.value( "C",  profile.C[j] );
	      writer.value( "MG", profile.Mg[j] );
	      writer.value( "TV", profile.TV[j] );
	      writer.unlevel( "OFFSET" );
	    }
	}

      // Transition-pair specific outputs for definite argmax changes only.
      std::set<std::pair<int,int>> trans_pairs;
      for ( const auto & e : evts )
	if ( e.from_st >= 0 && e.to_st >= 0 && e.from_st != e.to_st )
	  trans_pairs.insert( std::make_pair( e.from_st , e.to_st ) );

      for ( const auto & tp : trans_pairs )
	{
	  std::vector<hd_event_t> evts_pair;
	  for ( const auto & e : evts )
	    if ( e.from_st == tp.first && e.to_st == tp.second )
	      evts_pair.push_back( e );

	  if ( evts_pair.empty() ) continue;

	  hd_trans_stats_t tshape_pair;
	  hd_profile_t     profile_pair;
	  compute_trans_shape( der, evts_pair, rmask, dat, par, dat.Fs, threestate, tshape_pair, profile_pair );

	  const std::string trans_label = hd_state_label( threestate , tp.first ) + "->"
	    + hd_state_label( threestate , tp.second );

	  writer.level( trans_label , "TRANS" );
	  write_trans_stats( tshape_pair );

	  if ( profile_pair.valid )
	    {
	      for ( int j = 0; j < (int)profile_pair.offsets.size(); j++ )
		{
		  writer.level( profile_pair.offsets[j], "OFFSET" );
		  writer.value( "H",  profile_pair.H[j] );
		  writer.value( "C",  profile_pair.C[j] );
		  writer.value( "MG", profile_pair.Mg[j] );
		  writer.value( "TV", profile_pair.TV[j] );

		  if ( threestate )
		    {
		      writer.value( "P_W",  profile_pair.P[0][j] );
		      writer.value( "P_NR", profile_pair.P[1][j] );
		      writer.value( "P_R",  profile_pair.P[2][j] );
		    }
		  else
		    {
		      writer.value( "P_W",  profile_pair.P[0][j] );
		      writer.value( "P_N1", profile_pair.P[1][j] );
		      writer.value( "P_N2", profile_pair.P[2][j] );
		      writer.value( "P_N3", profile_pair.P[3][j] );
		      writer.value( "P_R",  profile_pair.P[4][j] );
		    }
		  writer.unlevel( "OFFSET" );
		}
	    }

	  writer.unlevel( "TRANS" );
	}
    }
  };

  if ( par.do_3state )
    {
      writer.level( "5", "STATE" );
      emit_state( false );
      writer.unlevel( "STATE" );

      writer.level( "3", "STATE" );
      emit_state( true );
      writer.unlevel( "STATE" );
    }
  else
    {
      emit_state( false );
    }
}


// ============================================================
// Layer 5: top-level command
// ============================================================

void proc_hdstats( edf_t & edf, param_t & param )
{
  // Parse parameters
  hd_params_t par;

  // Channel names for W, N1, N2, N3, R
  par.ch.resize( 5 );
  par.ch[0] = param.has( "W" )  ? param.value( "W" )  : "PP_W";
  par.ch[1] = param.has( "N1" ) ? param.value( "N1" ) : "PP_N1";
  par.ch[2] = param.has( "N2" ) ? param.value( "N2" ) : "PP_N2";
  par.ch[3] = param.has( "N3" ) ? param.value( "N3" ) : "PP_N3";
  par.ch[4] = param.has( "R" )  ? param.value( "R" )  : "PP_R";

  par.do_3state       = param.yesno( "3state" );
  par.window_sec      = param.has( "window"    ) ? param.requires_dbl( "window"     ) : 60.0;
  par.lag_sec         = param.has( "lag"       ) ? param.requires_dbl( "lag"        ) : 30.0;
  par.stable_min_sec  = param.has( "stable-min") ? param.requires_dbl( "stable-min" ) : 60.0;
  par.motion_th       = param.has( "motion-th" ) ? param.requires_dbl( "motion-th"  ) : 0.1;
  par.conf_th         = param.has( "conf-th"   ) ? param.requires_dbl( "conf-th"    ) : 0.8;
  par.annot_name      = param.has( "annot"     ) ? param.value( "annot" )              : "";
  par.verbose         = param.yesno( "verbose" );
  par.min_events_profile = param.has( "min-events" ) ? param.requires_int( "min-events" ) : 3;

  std::string method_str = param.has( "transition" ) ? param.value( "transition" ) : "motion";
  if      ( method_str == "hard"   ) par.method = HD_HARD;
  else if ( method_str == "motion" ) par.method = HD_MOTION;
  else if ( method_str == "both"   ) par.method = HD_BOTH;
  else Helper::halt( "HDSTATS: transition must be hard, motion, or both" );

  logger << "  HDSTATS: loading hypnodensity channels ["
	 << par.ch[0] << ", " << par.ch[1] << ", " << par.ch[2] << ", "
	 << par.ch[3] << ", " << par.ch[4] << "]\n";

  // Layer 1: load data
  hd_data_t dat;
  if ( ! dat.load( edf, par ) ) return;

  logger << "  HDSTATS: N=" << dat.N << " samples at Fs=" << dat.Fs << " Hz ("
	 << dat.N / dat.Fs / 3600.0 << " hours)\n";

  // Layer 2: derive per-sample signals
  hd_derived_t der;
  der.compute( dat, par );

  // Layer 3: detect transitions
  hd_trans_t tr;
  tr.detect( dat, der, par );

  logger << "  HDSTATS: detected " << tr.events5.size() << " transition events (5-state)\n";
  if ( par.do_3state )
    logger << "  HDSTATS: detected " << tr.events3.size() << " transition events (3-state)\n";

  // All-samples mask (true everywhere)
  std::vector<bool> all_mask( dat.N, true );

  // Global summary
  write_hdstats( der, tr, all_mask, all_mask, par, dat );

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
	  writer.value( "C",        der.C[i] );
	  writer.value( "MG",       der.Mg[i] );
	  writer.value( "TV",       der.TV[i] );
	  writer.value( "TV_LAG",   der.TV_lag[i] );
	  writer.value( "MIX_A",    der.mix_wn1[i] );
	  writer.value( "MIX_B",    der.mix_n2n3[i] );
	  writer.value( "MIX_C",    der.mix_rn1[i] );
	  writer.value( "ARGMAX",   der.argmax5[i] );
	  writer.value( "IS_TRANS", (int)tr.is_trans5[i] );
	  writer.value( "IS_STABLE",(int)tr.is_stable5[i] );
	  if ( par.do_3state )
	    {
	      writer.value( "H3",        der.H3[i] );
	      writer.value( "C3",        der.C3[i] );
	      writer.value( "TV3",       der.TV3[i] );
	      writer.value( "ARGMAX3",   der.argmax3[i] );
	      writer.value( "IS_TRANS3", (int)tr.is_trans3[i] );
	    }
	  writer.unlevel( "TIME" );
	}
    }

  // Stratification by annotation
  if ( par.annot_name.empty() ) return;

  annot_t * annot = (*edf.annotations)( par.annot_name );
  if ( annot == NULL )
    {
      logger << "  HDSTATS: annotation '" << par.annot_name << "' not found; skipping stratification\n";
      return;
    }

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
      write_hdstats( der, tr, mask5, mask5, par, dat );
      writer.unlevel( globals::annot_strat );
    }
}
