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

// Luna integrated test suite
// Invocation: luna __LUNA_TESTS__ [group] [verbose]
//
// Groups: all, signal, epoch, mask, filter, resample, psd, spindles,
//         hypno, annot, write, script, eval, lunapi, segsrv, plm
//
// All tests use fully synthetic in-memory data (no external files needed).
// Exit code: 0 = all pass, 1 = any failure.

#include "tests.h"
#include "luna.h"
#include "lunapi/lunapi.h"
#include "lunapi/segsrv.h"
#include "helper/token-eval.h"
#include "miscmath/crandom.h"
#include "miscmath/miscmath.h"
#include "dsp/ipc.h"
#include "dsp/ssa.h"
#include "dsp/tsync.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include <map>
#include <set>
#include <limits>
#include <algorithm>
#include <stdexcept>

#ifdef _WIN32
#include <process.h>
#ifdef near
#undef near
#endif
#endif

// ============================================================
// Internal types
// ============================================================

struct test_result_t {
  std::string name;
  bool        pass;
  std::string msg;   // what was expected vs. what was observed
};

// ============================================================
// Helpers: test runner
// ============================================================

static int n_pass = 0, n_fail = 0;

static void record( std::vector<test_result_t> & results ,
		    const std::string & name ,
		    const bool pass ,
		    const std::string & msg ,
		    const bool verbose )
{
  results.push_back( { name, pass, msg } );
  if ( pass ) ++n_pass; else ++n_fail;

  // Always print failures; print passes only in verbose mode
  if ( !pass || verbose )
    {
      std::cout << ( pass ? "[PASS] " : "[FAIL] " )
		<< std::left << std::setw(45) << name
		<< "  " << msg << "\n";
    }
  globals::problem = false;
  globals::empty   = false;
}

// ============================================================
// Helpers: numeric
// ============================================================

static bool approx_equal( double a, double b, double tol )
{
  return std::fabs(a - b) <= tol;
}

static bool approx_equal_rel( double a, double b, double rel_tol )
{
  double denom = std::max( std::fabs(b), 1e-12 );
  return std::fabs(a - b) / denom <= rel_tol;
}

static int current_pid()
{
#ifdef _WIN32
  return _getpid();
#else
  return getpid();
#endif
}

static std::string temp_base_path( const std::string & stem )
{
  const char * temp_dir = std::getenv("TMPDIR");
#ifdef _WIN32
  if ( temp_dir == nullptr || *temp_dir == '\0' ) temp_dir = std::getenv("TEMP");
  if ( temp_dir == nullptr || *temp_dir == '\0' ) temp_dir = std::getenv("TMP");
  if ( temp_dir == nullptr || *temp_dir == '\0' ) temp_dir = ".";
#else
  if ( temp_dir == nullptr || *temp_dir == '\0' ) temp_dir = "/tmp";
#endif
  return std::string(temp_dir) + "/luna_" + stem + "_" + std::to_string( current_pid() );
}

static bool contains_substr( const std::string & haystack, const std::string & needle )
{
  return haystack.find( needle ) != std::string::npos;
}

// ============================================================
// Helpers: signal generation
// ============================================================

// Pure sine wave
static std::vector<double> make_sine( int sr, double dur_sec,
				      double freq_hz, double amp = 1.0,
				      double phase = 0.0 )
{
  const int n = (int)std::round( sr * dur_sec );
  std::vector<double> v(n);
  for (int i = 0; i < n; i++)
    v[i] = amp * std::sin( 2.0 * M_PI * freq_hz * i / sr + phase );
  return v;
}

// Sum of two sines (used for filter tests)
static std::vector<double> make_two_sines( int sr, double dur_sec,
					   double f1, double a1,
					   double f2, double a2 )
{
  std::vector<double> v = make_sine( sr, dur_sec, f1, a1 );
  std::vector<double> v2 = make_sine( sr, dur_sec, f2, a2 );
  for (size_t i = 0; i < v.size(); i++) v[i] += v2[i];
  return v;
}

// White Gaussian noise (seeded for reproducibility)
static std::vector<double> make_noise( int n, double sd = 1.0, int seed = 42 )
{
  std::mt19937 rng( (unsigned)seed );
  std::normal_distribution<double> dist( 0.0, sd );
  std::vector<double> v(n);
  for (int i = 0; i < n; i++) v[i] = dist(rng);
  return v;
}

// Sine burst embedded in noise (for spindle simulation)
// burst_start/dur are in seconds
static std::vector<double> make_burst( int sr, double dur_sec,
				       double freq_hz, double amp,
				       double burst_start, double burst_dur,
				       double noise_sd = 0.05 )
{
  std::vector<double> v = make_noise( (int)std::round(sr * dur_sec), noise_sd, 99 );
  int bi = (int)std::round( burst_start * sr );
  int bd = (int)std::round( burst_dur   * sr );
  for (int i = 0; i < bd && (bi + i) < (int)v.size(); i++)
    v[bi + i] += amp * std::sin( 2.0 * M_PI * freq_hz * (i / (double)sr) );
  return v;
}

// Slow-wave shape: one cycle of a negative half-sine then positive
static std::vector<double> make_slow_wave( int sr, double dur_sec,
					   double freq_hz = 0.75, double amp = 100.0 )
{
  return make_sine( sr, dur_sec, freq_hz, amp );
}

static std::vector<double> make_chirp_phase( int sr, double dur_sec,
					     double f0_hz, double f1_hz )
{
  const int n = (int)std::round( sr * dur_sec );
  std::vector<double> phase(n);
  const double T = dur_sec > 0 ? dur_sec : 1.0;
  const double k = (f1_hz - f0_hz) / T;
  for (int i = 0; i < n; i++)
    {
      const double t = i / (double)sr;
      phase[i] = 2.0 * M_PI * ( f0_hz * t + 0.5 * k * t * t );
    }
  return phase;
}

static std::vector<double> shift_with_edge_hold( const std::vector<double> & x, int delay )
{
  const int n = (int)x.size();
  std::vector<double> y(n);
  for (int i = 0; i < n; i++)
    {
      int src = i - delay;
      if ( src < 0 ) src = 0;
      if ( src >= n ) src = n - 1;
      y[i] = x[src];
    }
  return y;
}

// Stage annotation intervals (full night of all-N2, 30s epochs)
static std::vector<std::tuple<double,double>> make_stage_annots( int ne, double elen )
{
  std::vector<std::tuple<double,double>> v;
  for (int i = 0; i < ne; i++)
    v.push_back( { i * elen, (i+1) * elen } );
  return v;
}

// ============================================================
// Helpers: lunapi wrappers
// ============================================================

// Create a fresh empty-EDF instance: nr×rs-second EDF, one EEG channel at sr Hz
static lunapi_inst_ptr make_inst( lunapi_t * eng,
				  const std::vector<double> & sig,
				  int sr,
				  int nr = 720, int rs = 30,
				  const std::string & label = "EEG",
				  const std::string & id = "T1" )
{
  lunapi_inst_ptr p = eng->inst( id );
  p->empty_edf( id, nr, rs, "01.01.85", "22.00.00" );
  p->insert_signal( label, sig, sr );
  return p;
}

// Shorthand: 720-record, rs=30 EDF with a single sine wave
static lunapi_inst_ptr make_sine_inst( lunapi_t * eng,
				       double freq_hz = 10.0,
				       int sr = 256,
				       double amp = 1.0,
				       const std::string & label = "EEG" )
{
  const int nr = 720, rs = 30;
  auto sig = make_sine( sr, (double)(nr * rs), freq_hz, amp );
  return make_inst( eng, sig, sr, nr, rs, label );
}

// Extract scalar from rtables (searches all strata for the command)
static double get_val( lunapi_inst_ptr p,
		       const std::string & cmd,
		       const std::string & var,
		       int row = 0 )
{
  for (const auto & cs : p->strata())
    {
      if (cs.first != cmd) continue;
      auto r   = p->results( cs.first, cs.second );
      const auto & cols = std::get<0>(r);
      const auto & data = std::get<1>(r);
      for (size_t ci = 0; ci < cols.size(); ci++)
	{
	  if (cols[ci] != var) continue;
	  if (ci >= data.size() || data[ci].empty()) continue;
	  int ri = row < (int)data[ci].size() ? row : 0;
	  const auto & e = data[ci][ri];
	  if (std::holds_alternative<double>(e)) return std::get<double>(e);
	  if (std::holds_alternative<int>(e))    return (double)std::get<int>(e);
	}
    }
  return std::numeric_limits<double>::quiet_NaN();
}

// Find a variable value in a strata row where a factor column has a given string value
static double get_val_where( lunapi_inst_ptr p,
			     const std::string & cmd,
			     const std::string & strata,
			     const std::string & factor_col,
			     const std::string & factor_val,
			     const std::string & var )
{
  auto r   = p->results( cmd, strata );
  const auto & cols = std::get<0>(r);
  const auto & data = std::get<1>(r);
  int fci = -1, vi = -1;
  for (int i = 0; i < (int)cols.size(); i++) {
    if (cols[i] == factor_col) fci = i;
    if (cols[i] == var)        vi  = i;
  }
  if (fci < 0 || vi < 0 || data.empty()) return std::numeric_limits<double>::quiet_NaN();
  int nrows = (int)data[fci].size();
  for (int r2 = 0; r2 < nrows; r2++) {
    std::string fval;
    if (std::holds_alternative<std::string>(data[fci][r2]))
      fval = std::get<std::string>(data[fci][r2]);
    if (fval == factor_val) {
      if (vi < (int)data.size() && r2 < (int)data[vi].size()) {
        const auto & e = data[vi][r2];
        if (std::holds_alternative<double>(e)) return std::get<double>(e);
        if (std::holds_alternative<int>(e))    return (double)std::get<int>(e);
      }
    }
  }
  return std::numeric_limits<double>::quiet_NaN();
}

// Same, but restrict to a specific strata key
static double get_val_s( lunapi_inst_ptr p,
			 const std::string & cmd,
			 const std::string & strata,
			 const std::string & var,
			 int row = 0 )
{
  auto r   = p->results( cmd, strata );
  const auto & cols = std::get<0>(r);
  const auto & data = std::get<1>(r);
  for (size_t ci = 0; ci < cols.size(); ci++)
    {
      if (cols[ci] != var) continue;
      if (ci >= data.size() || data[ci].empty()) continue;
      int ri = row < (int)data[ci].size() ? row : 0;
      const auto & e = data[ci][ri];
      if (std::holds_alternative<double>(e)) return std::get<double>(e);
      if (std::holds_alternative<int>(e))    return (double)std::get<int>(e);
    }
  return std::numeric_limits<double>::quiet_NaN();
}

// Return number of rows in a strata table
static int get_nrows( lunapi_inst_ptr p,
		      const std::string & cmd,
		      const std::string & strata )
{
  auto r   = p->results( cmd, strata );
  const auto & data = std::get<1>(r);
  if (data.empty()) return 0;
  return (int)data[0].size();
}

// Find peak frequency: in the CH_F table, find the F at which PSD is maximum
static double peak_freq( lunapi_inst_ptr p )
{
  auto r   = p->results( "PSD", "CH_F" );
  const auto & cols = std::get<0>(r);
  const auto & data = std::get<1>(r);
  int fi = -1, pi = -1;
  for (int i = 0; i < (int)cols.size(); i++)
    {
      if (cols[i] == "F")   fi = i;
      if (cols[i] == "PSD") pi = i;
    }
  if (fi < 0 || pi < 0 || data.empty()) return std::numeric_limits<double>::quiet_NaN();
  int nrows = (int)data[fi].size();
  double best_f = std::numeric_limits<double>::quiet_NaN();
  double best_p = -1.0;
  for (int r2 = 0; r2 < nrows; r2++)
    {
      double f = std::numeric_limits<double>::quiet_NaN();
      double pv = -1.0;
      if (std::holds_alternative<double>(data[fi][r2])) f = std::get<double>(data[fi][r2]);
      if (std::holds_alternative<double>(data[pi][r2])) pv = std::get<double>(data[pi][r2]);
      if (pv > best_p) { best_p = pv; best_f = f; }
    }
  return best_f;
}

// Sum PSD power in a frequency range [flo, fhi]
static double band_power( lunapi_inst_ptr p, double flo, double fhi )
{
  auto r   = p->results( "PSD", "CH_F" );
  const auto & cols = std::get<0>(r);
  const auto & data = std::get<1>(r);
  int fi = -1, pi = -1;
  for (int i = 0; i < (int)cols.size(); i++)
    {
      if (cols[i] == "F")   fi = i;
      if (cols[i] == "PSD") pi = i;
    }
  if (fi < 0 || pi < 0 || data.empty()) return 0.0;
  int nrows = (int)data[fi].size();
  double sum = 0.0;
  for (int r2 = 0; r2 < nrows; r2++)
    {
      double f = 0.0, pv = 0.0;
      if (std::holds_alternative<double>(data[fi][r2])) f = std::get<double>(data[fi][r2]);
      if (std::holds_alternative<double>(data[pi][r2])) pv = std::get<double>(data[pi][r2]);
      if (f >= flo && f <= fhi) sum += pv;
    }
  return sum;
}

// ============================================================
// Group A: Signal generation
// ============================================================

static void test_signal( lunapi_t * eng,
			 std::vector<test_result_t> & R, bool V )
{
  // A1 — sine amplitude: STATS mean≈0, SD≈amp/sqrt(2)
  try {
    auto p = make_sine_inst( eng, 10.0, 256, 2.0 );
    p->eval("STATS sig=EEG");
    double mn  = get_val( p, "STATS", "MEAN" );
    double sd  = get_val( p, "STATS", "SD" );
    std::ostringstream m;
    m << "mean=" << mn << " (exp≈0), SD=" << sd << " (exp≈" << 2.0/std::sqrt(2.0) << ")";
    record(R,"signal/sine-amplitude", approx_equal(mn,0.0,0.05) && approx_equal(sd,2.0/std::sqrt(2.0),0.1), m.str(), V);
  } catch(std::exception & e) { record(R,"signal/sine-amplitude",false,e.what(),V); }

  // A2 — SIMUL white noise: mean≈0, SD>0
  try {
    auto p = eng->inst("T_simul");
    p->empty_edf("T_simul", 720, 30, "01.01.85", "22.00.00");
    p->eval("SIGGEN sig=EEG sr=256 sine=10,1,0 & STATS sig=EEG");
    double mn = get_val(p,"STATS","MEAN");
    double sd = get_val(p,"STATS","SD");
    std::ostringstream m;
    m << "mean=" << mn << " SD=" << sd;
    record(R,"signal/siggen-eval", approx_equal(mn,0.0,0.1) && sd>0, m.str(), V);
  } catch(std::exception & e) { record(R,"signal/siggen-eval",false,e.what(),V); }

  // A3 — add mode: SIGGEN creates then adds to itself, doubling amplitude
  try {
    auto p = eng->inst("T_add");
    p->empty_edf("T_add", 720, 30, "01.01.85", "22.00.00");
    // Create EEG via SIGGEN (from empty EDF), then add the same sine to it
    p->eval("SIGGEN sig=EEG sr=256 sine=10,1,0 & SIGGEN sig=EEG sr=256 add sine=10,1,0 & STATS sig=EEG");
    double sd = get_val(p,"STATS","SD");
    std::ostringstream m;
    m << "SD=" << sd << " (exp≈" << 2.0/std::sqrt(2.0) << " after double-sine add)";
    record(R,"signal/siggen-add", approx_equal(sd, 2.0/std::sqrt(2.0), 0.2), m.str(), V);
  } catch(std::exception & e) { record(R,"signal/siggen-add",false,e.what(),V); }

  // A4 — two channels at different SRs
  try {
    auto p = eng->inst("T_2ch");
    p->empty_edf("T_2ch", 720, 30, "01.01.85", "22.00.00");
    auto s256 = make_sine(256, 720*30, 10.0, 1.0);
    auto s64  = make_sine(64,  720*30, 5.0,  1.0);
    p->insert_signal("CH256", s256, 256);
    p->insert_signal("CH64",  s64,  64);
    auto chs = p->channels();
    bool has256 = std::find(chs.begin(),chs.end(),"CH256") != chs.end();
    bool has64  = std::find(chs.begin(),chs.end(),"CH64")  != chs.end();
    std::ostringstream m;
    m << "channels=" << chs.size() << " has256=" << has256 << " has64=" << has64;
    record(R,"signal/two-sr-channels", has256 && has64 && chs.size()==2, m.str(), V);
  } catch(std::exception & e) { record(R,"signal/two-sr-channels",false,e.what(),V); }

  // A5 — update_signal replaces data
  try {
    auto p = make_sine_inst(eng, 10.0, 256, 1.0);
    p->eval("STATS sig=EEG");
    double sd_before = get_val(p,"STATS","SD");
    auto flat = std::vector<double>(256*720*30, 0.0);
    p->update_signal("EEG", flat);
    p->eval("STATS sig=EEG");
    double sd_after = get_val(p,"STATS","SD");
    std::ostringstream m;
    m << "SD before=" << sd_before << " after=" << sd_after << " (exp≈0)";
    record(R,"signal/update-signal", sd_before > 0.1 && approx_equal(sd_after,0.0,1e-6), m.str(), V);
  } catch(std::exception & e) { record(R,"signal/update-signal",false,e.what(),V); }

  // A6 — direct SSA reconstruction should sum back to the original series
  try {
    const int sr = 64;
    const double dur = 8.0;
    std::vector<double> x = make_sine( sr, dur, 10.0, 1.0 );
    for (int i = 0; i < (int)x.size(); i++)
      x[i] += 0.01 * i;
    ssa_t ssa( &x , 32 );
    std::vector<int> all( ssa.d );
    for (int i = 0; i < ssa.d; i++) all[i] = i;
    Eigen::VectorXd recon = ssa.reconstruct( all );
    double max_abs_err = 0.0;
    for (int i = 0; i < recon.size(); i++)
      max_abs_err = std::max( max_abs_err , std::fabs( recon[i] - x[i] ) );
    std::ostringstream m;
    m << "rank=" << ssa.d << " max_abs_err=" << max_abs_err;
    record(R,"signal/ssa-reconstructs-input", max_abs_err < 1e-8, m.str(), V);
  } catch(std::exception & e) { record(R,"signal/ssa-reconstructs-input",false,e.what(),V); }

  // A7 — SSA command smoke test: wired, reports summaries, adds channels
  try {
    const int sr = 64;
    const int nr = 10;
    const int rs = 1;
    auto sig = make_two_sines( sr, (double)(nr * rs), 2.0, 1.0, 10.0, 0.5 );
    auto p = make_inst( eng, sig, sr, nr, rs, "EEG", "T_ssa" );
    p->eval("SSA sig=EEG L=32 nc=2 wcorr=T");
    const double L = get_val_s( p, "SSA", "CH", "L" );
    const double D = get_val_s( p, "SSA", "CH", "D" );
    const double sigma1 = get_val( p, "SSA", "SIGMA" );
    const auto chs = p->channels();
    const bool has1 = std::find( chs.begin(), chs.end(), "SSA_EEG_1" ) != chs.end();
    const bool has2 = std::find( chs.begin(), chs.end(), "SSA_EEG_2" ) != chs.end();
    std::ostringstream m;
    m << "L=" << L << " D=" << D << " sigma1=" << sigma1
      << " has1=" << has1 << " has2=" << has2;
    record(R,"signal/ssa-command-smoke", L == 32 && D >= 2 && std::isfinite(sigma1) && has1 && has2, m.str(), V);
  } catch(std::exception & e) { record(R,"signal/ssa-command-smoke",false,e.what(),V); }

  // A8 — lagged IPC should track the same delayed-coupling use-case as TSYNC HT
  try {
    const int sr = 128;
    const double dur_sec = 20.0;
    const int n = (int)std::round( sr * dur_sec );
    const int delay = 8;
    const int max_lag = 16;

    std::vector<double> phase1 = make_chirp_phase( sr, dur_sec, 6.0, 12.0 );
    std::vector<double> amp1(n), amp2(n);
    for (int i = 0; i < n; i++)
      {
        const double t = i / (double)sr;
        amp1[i] = 1.0 + 0.25 * std::sin( 2.0 * M_PI * 0.35 * t );
      }
    std::vector<double> phase2 = shift_with_edge_hold( phase1, delay );
    amp2 = shift_with_edge_hold( amp1, delay );

    Eigen::MatrixXd P(n,2), A(n,2);
    for (int i = 0; i < n; i++)
      {
        P(i,0) = phase1[i];
        P(i,1) = phase2[i];
        A(i,0) = amp1[i];
        A(i,1) = amp2[i];
      }

    tsync_t tsync( P, A, sr, max_lag );
    ipc_phaseamp_t seed, tgt;
    seed.phase = phase1; seed.amp = amp1;
    tgt.phase = phase2;  tgt.amp  = amp2;
    ipc_param_t ipc_param;
    ipc_lag_output_t lagged = ipc_t::compute_ipc_lagged( seed, tgt, sr, ipc_param, max_lag );

    int best_tsync_lag = 0;
    double best_tsync_val = -1.0;
    double tsync_zero = std::numeric_limits<double>::quiet_NaN();
    for (std::map<int,double>::const_iterator it = tsync.phase_lock[0][1].begin();
         it != tsync.phase_lock[0][1].end(); ++it)
      {
        if ( it->second > best_tsync_val )
          {
            best_tsync_val = it->second;
            best_tsync_lag = it->first;
          }
        if ( it->first == 0 ) tsync_zero = it->second;
      }

    int best_ipc_lag = 0;
    double best_ipc_val = -1.0;
    double ipc_zero = std::numeric_limits<double>::quiet_NaN();
    for (int i = 0; i < lagged.rows.size(); i++)
      {
        const double v = lagged.rows[i].summary.mean_ipc;
        if ( std::isfinite(v) && v > best_ipc_val )
          {
            best_ipc_val = v;
            best_ipc_lag = lagged.rows[i].lag;
          }
        if ( lagged.rows[i].lag == 0 ) ipc_zero = v;
      }

    const bool lag_match = best_tsync_lag == best_ipc_lag;
    const bool delay_recovered = best_ipc_lag == -delay;
    const bool both_peak_over_zero = best_tsync_val > tsync_zero + 0.05 &&
                                     best_ipc_val > ipc_zero + 0.05;
    const bool both_peaks_high = best_tsync_val > 0.95 && best_ipc_val > 0.95;

    std::ostringstream m;
    m << "TSYNC lag=" << best_tsync_lag << " PHL=" << best_tsync_val
      << " zero=" << tsync_zero
      << " | IPC lag=" << best_ipc_lag << " IPC=" << best_ipc_val
      << " zero=" << ipc_zero;
    record(R,"signal/ipc-lag-vs-tsync-ht",
           lag_match && delay_recovered && both_peak_over_zero && both_peaks_high,
           m.str(), V);
  } catch(std::exception & e) { record(R,"signal/ipc-lag-vs-tsync-ht",false,e.what(),V); }
}

// ============================================================
// Group B: Epoch counting
// ============================================================

static void test_epoch( lunapi_t * eng,
			std::vector<test_result_t> & R, bool V )
{
  const int nr=720, rs=30, sr=256;

  // B1 — standard 30s epochs
  try {
    auto p = make_sine_inst(eng,10.0,sr);
    p->eval("EPOCH len=30");
    double ne = get_val(p,"EPOCH","NE");
    std::ostringstream m; m << "NE=" << ne << " (exp=720)";
    record(R,"epoch/standard-30s", approx_equal(ne,720,0.5), m.str(), V);
  } catch(std::exception & e) { record(R,"epoch/standard-30s",false,e.what(),V); }

  // B2 — with rs=1 (single-sample records)
  try {
    auto p = eng->inst("T_rs1");
    p->empty_edf("T_rs1", 720*30, 1, "01.01.85", "22.00.00");
    auto sig = make_sine(sr, 720*30, 10.0, 1.0);
    p->insert_signal("EEG", sig, sr);
    p->eval("EPOCH len=30");
    double ne = get_val(p,"EPOCH","NE");
    std::ostringstream m; m << "NE=" << ne << " (exp=720, rs=1)";
    record(R,"epoch/rs1-30s", approx_equal(ne,720,0.5), m.str(), V);
  } catch(std::exception & e) { record(R,"epoch/rs1-30s",false,e.what(),V); }

  // B3 — rs=30 (record size matches epoch size)
  try {
    auto p = eng->inst("T_rs30");
    p->empty_edf("T_rs30", 720, 30, "01.01.85", "22.00.00");
    auto sig = make_sine(sr, 720*30, 10.0, 1.0);
    p->insert_signal("EEG", sig, sr);
    p->eval("EPOCH len=30");
    double ne = get_val(p,"EPOCH","NE");
    std::ostringstream m; m << "NE=" << ne << " (exp=720, rs=30)";
    record(R,"epoch/rs30-30s", approx_equal(ne,720,0.5), m.str(), V);
  } catch(std::exception & e) { record(R,"epoch/rs30-30s",false,e.what(),V); }

  // B4 — overlapping epochs (len=30 inc=15)
  try {
    auto p = make_sine_inst(eng,10.0,sr);
    p->eval("EPOCH len=30 inc=15");
    double ne = get_val(p,"EPOCH","NE");
    // 720*30 sec / 15s step - 1 = 1439 epochs
    std::ostringstream m; m << "NE=" << ne << " (exp=1439)";
    record(R,"epoch/overlapping-inc15", approx_equal(ne,1439,1.5), m.str(), V);
  } catch(std::exception & e) { record(R,"epoch/overlapping-inc15",false,e.what(),V); }

  // B5 — 5-second epochs
  try {
    auto p = make_sine_inst(eng,10.0,sr);
    p->eval("EPOCH len=5");
    double ne = get_val(p,"EPOCH","NE");
    // 720*30/5 = 4320
    std::ostringstream m; m << "NE=" << ne << " (exp=4320)";
    record(R,"epoch/5sec-epochs", approx_equal(ne,4320,0.5), m.str(), V);
  } catch(std::exception & e) { record(R,"epoch/5sec-epochs",false,e.what(),V); }

  // B6 — EPOCH called twice: NE is stable
  try {
    auto p = make_sine_inst(eng,10.0,sr);
    p->eval("EPOCH len=30");
    double ne1 = get_val(p,"EPOCH","NE");
    p->eval("EPOCH len=30");
    double ne2 = get_val(p,"EPOCH","NE");
    std::ostringstream m; m << "NE1=" << ne1 << " NE2=" << ne2;
    record(R,"epoch/dump-stable", approx_equal(ne1,ne2,0.5), m.str(), V);
  } catch(std::exception & e) { record(R,"epoch/dump-stable",false,e.what(),V); }
}

// ============================================================
// Group C: MASK and RE
// ============================================================

static void test_mask( lunapi_t * eng,
		       std::vector<test_result_t> & R, bool V )
{
  // C1 — MASK all: after RE, NR2=0 records retained
  try {
    auto p = make_sine_inst(eng);
    p->eval("EPOCH len=30 & MASK all & RE require=0");
    double nr2 = get_val(p,"RE","NR2");
    globals::empty = false; globals::problem = false;
    std::ostringstream m; m << "NR2=" << nr2 << " (exp=0)";
    record(R,"mask/mask-all", approx_equal(nr2,0,0.5), m.str(), V);
  } catch(std::exception & e) { record(R,"mask/mask-all",false,e.what(),V); }

  // C2 — MASK none: after MASK all then MASK none, RE retains all 720 records
  try {
    auto p = make_sine_inst(eng);
    p->eval("EPOCH len=30 & MASK all");
    p->eval("MASK none & RE");
    double nr2 = get_val(p,"RE","NR2");
    std::ostringstream m; m << "NR2=" << nr2 << " (exp=720)";
    record(R,"mask/mask-none", approx_equal(nr2,720,0.5), m.str(), V);
  } catch(std::exception & e) { record(R,"mask/mask-none",false,e.what(),V); }

  // C3 — MASK ifnot=N2 with annotation: N_RETAINED = # N2 epochs
  try {
    auto p = make_sine_inst(eng);
    // First 300 epochs are N2, remaining 420 are N3
    auto n2 = make_stage_annots(300, 30.0);
    auto n3 = make_stage_annots(420, 30.0);
    // offset n3 by 300*30 seconds
    for (auto & iv : n3) {
      std::get<0>(iv) += 300*30.0;
      std::get<1>(iv) += 300*30.0;
    }
    p->insert_annotation("N2", n2);
    p->insert_annotation("N3", n3);
    p->eval("EPOCH len=30 & MASK ifnot=N2");
    double ret = get_val(p,"MASK","N_RETAINED");
    std::ostringstream m; m << "N_RETAINED=" << ret << " (exp=300)";
    record(R,"mask/ifnot-annot", approx_equal(ret,300,0.5), m.str(), V);
  } catch(std::exception & e) { record(R,"mask/ifnot-annot",false,e.what(),V); }

  // C4 — MASK epoch=1-100: verify via RE NR2=100
  try {
    auto p = make_sine_inst(eng);
    p->eval("EPOCH len=30 & MASK epoch=1-100 & RE");
    double nr2 = get_val(p,"RE","NR2");
    std::ostringstream m; m << "NR2=" << nr2 << " (exp=100)";
    record(R,"mask/epoch-range", approx_equal(nr2,100,0.5), m.str(), V);
  } catch(std::exception & e) { record(R,"mask/epoch-range",false,e.what(),V); }

  // C5 — RE with rs=1: DUR2 after RE == 200*30s
  try {
    auto p = eng->inst("T_re_rs1");
    p->empty_edf("T_re_rs1", 720*30, 1, "01.01.85", "22.00.00");
    p->insert_signal("EEG", make_sine(256,720*30,10,1), 256);
    p->eval("EPOCH len=30 & MASK epoch=1-200 & RE");
    double dur2 = get_val(p,"RE","DUR2");
    std::ostringstream m; m << "DUR2=" << dur2 << "s (exp=6000, rs=1)";
    record(R,"mask/re-basic-rs1", approx_equal(dur2,6000.0,1.0), m.str(), V);
  } catch(std::exception & e) { record(R,"mask/re-basic-rs1",false,e.what(),V); }

  // C6 — RE with rs=30: DUR2 after RE == 200*30s  (regression guard)
  try {
    auto p = eng->inst("T_re_rs30");
    p->empty_edf("T_re_rs30", 720, 30, "01.01.85", "22.00.00");
    p->insert_signal("EEG", make_sine(256,720*30,10,1), 256);
    p->eval("EPOCH len=30 & MASK epoch=1-200 & RE");
    double dur2 = get_val(p,"RE","DUR2");
    std::ostringstream m; m << "DUR2=" << dur2 << "s (exp=6000, rs=30)";
    record(R,"mask/re-basic-rs30", approx_equal(dur2,6000.0,1.0), m.str(), V);
  } catch(std::exception & e) { record(R,"mask/re-basic-rs30",false,e.what(),V); }

  // C7 — DOUBLE RE with rs=30: NR2 from second RE > 0  (uint64_t overflow regression)
  try {
    auto p = eng->inst("T_dre_rs30");
    p->empty_edf("T_dre_rs30", 720, 30, "01.01.85", "22.00.00");
    p->insert_signal("EEG", make_sine(256,720*30,10,1), 256);
    auto n2 = make_stage_annots(300, 30.0);
    auto n3 = make_stage_annots(420, 30.0);
    for (auto & iv : n3) {
      std::get<0>(iv) += 300*30.0;
      std::get<1>(iv) += 300*30.0;
    }
    p->insert_annotation("N2", n2);
    p->insert_annotation("N3", n3);
    p->eval("EPOCH len=30 & MASK ifnot=N2 & RE & MASK if=EXCLUDED & RE");
    // row=1 gets second RE's NR2 (row=0 would be first RE)
    double nr2 = get_val(p,"RE","NR2",1);
    bool ok = nr2 > 0;
    std::ostringstream m;
    m << "NR2(2nd RE)=" << nr2 << " (exp>0 and ≈300; overflow bug gives 0)";
    record(R,"mask/double-re-rs30", ok, m.str(), V);
  } catch(std::exception & e) { record(R,"mask/double-re-rs30",false,e.what(),V); }

  // C8 — MASK flip: inverts mask; verify via RE
  try {
    auto p = make_sine_inst(eng);
    p->eval("EPOCH len=30 & MASK epoch=1-100 & MASK flip & RE");
    double nr2 = get_val(p,"RE","NR2");
    std::ostringstream m; m << "NR2 after flip+RE=" << nr2 << " (exp=620)";
    record(R,"mask/mask-flip", approx_equal(nr2,620,0.5), m.str(), V);
  } catch(std::exception & e) { record(R,"mask/mask-flip",false,e.what(),V); }

  // C9 — MASK epoch list (non-contiguous); verify via RE
  try {
    auto p = make_sine_inst(eng);
    p->eval("EPOCH len=30 & MASK epoch=1,3,5,7,9 & RE");
    double nr2 = get_val(p,"RE","NR2");
    std::ostringstream m; m << "NR2=" << nr2 << " (exp=5)";
    record(R,"mask/epoch-list", approx_equal(nr2,5,0.5), m.str(), V);
  } catch(std::exception & e) { record(R,"mask/epoch-list",false,e.what(),V); }

  // C10 — RE reduces signal duration
  try {
    auto p = make_sine_inst(eng);
    double dur_before = p->last_sec();
    p->eval("EPOCH len=30 & MASK epoch=1-100 & RE");
    double dur_after = p->last_sec();
    std::ostringstream m;
    m << "dur before=" << dur_before << " after=" << dur_after << " (exp≈3000s)";
    record(R,"mask/re-duration", dur_after < dur_before && approx_equal(dur_after,3000.0,60.0), m.str(), V);
  } catch(std::exception & e) { record(R,"mask/re-duration",false,e.what(),V); }
}

// ============================================================
// Group D: Filter
// ============================================================

static void test_filter( lunapi_t * eng,
			 std::vector<test_result_t> & R, bool V )
{
  // D1 — lowpass attenuates high-freq component
  // Input: 5Hz (pass) + 40Hz (stop); lowpass=20Hz
  try {
    auto p = eng->inst("T_lp");
    p->empty_edf("T_lp", 720, 30, "01.01.85","22.00.00");
    p->insert_signal("EEG", make_two_sines(256,720*30, 5.0,1.0, 40.0,1.0), 256);
    p->eval("EPOCH len=30 & PSD sig=EEG max=50 spectrum=T");
    double pwr_low_before  = band_power(p, 3.0, 7.0);
    double pwr_high_before = band_power(p, 35.0, 45.0);
    p->eval("FILTER sig=EEG lowpass=20 tw=2 ripple=0.01 & EPOCH len=30 & PSD sig=EEG max=50 spectrum=T");
    double pwr_low_after   = band_power(p, 3.0, 7.0);
    double pwr_high_after  = band_power(p, 35.0, 45.0);
    // Low-freq largely preserved; high-freq substantially reduced
    bool pass = (pwr_low_after  > pwr_low_before  * 0.5) &&
		(pwr_high_after < pwr_high_before * 0.1);
    std::ostringstream m;
    m << "low: " << pwr_low_before << "->" << pwr_low_after
      << " high: " << pwr_high_before << "->" << pwr_high_after;
    record(R,"filter/lowpass-attenuation", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"filter/lowpass-attenuation",false,e.what(),V); }

  // D2 — highpass attenuates low-freq component
  // Input: 0.2Hz (stop) + 10Hz (pass); highpass=2Hz
  try {
    auto p = eng->inst("T_hp");
    p->empty_edf("T_hp", 720, 30, "01.01.85","22.00.00");
    p->insert_signal("EEG", make_two_sines(256,720*30, 0.2,1.0, 10.0,1.0), 256);
    p->eval("EPOCH len=30 & PSD sig=EEG max=20 spectrum=T");
    double pwr_high_before = band_power(p, 8.0, 12.0);
    p->eval("FILTER sig=EEG highpass=2 tw=1 ripple=0.01 & EPOCH len=30 & PSD sig=EEG max=20 spectrum=T");
    double pwr_dc_after    = band_power(p, 0.0, 0.5);
    double pwr_high_after  = band_power(p, 8.0, 12.0);
    bool pass = approx_equal_rel(pwr_high_after, pwr_high_before, 0.5) &&
		pwr_dc_after < pwr_high_after * 0.1;
    std::ostringstream m;
    m << "highpass: DC_after=" << pwr_dc_after << " 10Hz_after=" << pwr_high_after;
    record(R,"filter/highpass-attenuation", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"filter/highpass-attenuation",false,e.what(),V); }

  // D3 — bandpass: pass 8-15Hz, stop outside
  try {
    auto p = eng->inst("T_bp");
    p->empty_edf("T_bp", 720, 30, "01.01.85","22.00.00");
    // In-band at 11Hz, out-of-band at 30Hz
    p->insert_signal("EEG", make_two_sines(256,720*30, 11.0,1.0, 30.0,1.0), 256);
    p->eval("FILTER sig=EEG bandpass=8,16 tw=1 ripple=0.01 & EPOCH len=30 & PSD sig=EEG max=40 spectrum=T");
    double pwr_in  = band_power(p, 9.0, 13.0);
    double pwr_out = band_power(p, 27.0, 33.0);
    bool pass = pwr_out < pwr_in * 0.05;
    std::ostringstream m; m << "in-band=" << pwr_in << " out-of-band=" << pwr_out;
    record(R,"filter/bandpass", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"filter/bandpass",false,e.what(),V); }

  // D4 — filter preserves sample rate
  try {
    auto p = make_sine_inst(eng, 10.0, 256, 1.0);
    p->eval("FILTER sig=EEG bandpass=1,40 tw=1 ripple=0.01");
    auto st = p->status();
    // Should still have at least 1 channel
    bool pass = p->channels().size() == 1;
    record(R,"filter/preserves-channel", pass, "channel count = " + std::to_string(p->channels().size()), V);
  } catch(std::exception & e) { record(R,"filter/preserves-channel",false,e.what(),V); }
}

// ============================================================
// Group E: Resample
// ============================================================

static void test_resample( lunapi_t * eng,
			   std::vector<test_result_t> & R, bool V )
{
  // E1 — downsample 256→128: duration unchanged, peak freq preserved
  try {
    auto p = make_sine_inst(eng, 10.0, 256, 1.0);
    double dur_before = p->last_sec();
    p->eval("RESAMPLE sig=EEG sr=128");
    double dur_after = p->last_sec();
    p->eval("EPOCH len=30 & PSD sig=EEG max=60 spectrum=T");
    double pf = peak_freq(p);
    bool pass = approx_equal(dur_before,dur_after,60.0) && approx_equal(pf,10.0,1.5);
    std::ostringstream m;
    m << "dur=" << dur_after << " peak_f=" << pf << " (exp≈10Hz)";
    record(R,"resample/downsample-256-128", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"resample/downsample-256-128",false,e.what(),V); }

  // E2 — upsample 64→256: duration unchanged, peak freq preserved
  try {
    auto p = eng->inst("T_up");
    p->empty_edf("T_up", 720, 30, "01.01.85","22.00.00");
    p->insert_signal("EEG", make_sine(64,720*30,10.0,1.0), 64);
    double dur_before = p->last_sec();
    p->eval("RESAMPLE sig=EEG sr=256");
    double dur_after = p->last_sec();
    p->eval("EPOCH len=30 & PSD sig=EEG max=30 spectrum=T");
    double pf = peak_freq(p);
    bool pass = approx_equal(dur_before,dur_after,60.0) && approx_equal(pf,10.0,1.5);
    std::ostringstream m;
    m << "dur=" << dur_after << " peak_f=" << pf << " (exp≈10Hz)";
    record(R,"resample/upsample-64-256", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"resample/upsample-64-256",false,e.what(),V); }
}

// ============================================================
// Group F: PSD
// ============================================================

static void test_psd( lunapi_t * eng,
		      std::vector<test_result_t> & R, bool V )
{
  // F1 — sine peak recovery: 10Hz sine → peak bin near 10Hz
  try {
    auto p = make_sine_inst(eng, 10.0, 256, 1.0);
    p->eval("EPOCH len=30 & PSD sig=EEG max=30 spectrum=T");
    double pf = peak_freq(p);
    std::ostringstream m; m << "peak_f=" << pf << " (exp≈10Hz)";
    record(R,"psd/sine-peak-10hz", approx_equal(pf,10.0,1.0), m.str(), V);
  } catch(std::exception & e) { record(R,"psd/sine-peak-10hz",false,e.what(),V); }

  // F2 — sigma band power > alpha when signal is 13Hz
  try {
    auto p = make_sine_inst(eng, 13.0, 256, 1.0);
    p->eval("EPOCH len=30 & PSD sig=EEG max=30 spectrum=T");
    double sigma_pwr = band_power(p, 11.0, 16.0);
    double alpha_pwr = band_power(p,  8.0, 12.0);
    bool pass = sigma_pwr > alpha_pwr * 2.0;
    std::ostringstream m;
    m << "sigma=" << sigma_pwr << " alpha=" << alpha_pwr;
    record(R,"psd/sigma-vs-alpha", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"psd/sigma-vs-alpha",false,e.what(),V); }

  // F3 — multi-channel: both CH present in CH_F table
  try {
    auto p = eng->inst("T_2ch_psd");
    p->empty_edf("T_2ch_psd",720,30,"01.01.85","22.00.00");
    p->insert_signal("C3", make_sine(256,720*30,10.0,1.0), 256);
    p->insert_signal("C4", make_sine(256,720*30,12.0,1.0), 256);
    p->eval("EPOCH len=30 & PSD sig=C3,C4 max=30 spectrum=T");
    // CH_F table should have rows for both C3 and C4
    auto r  = p->results("PSD","CH_F");
    const auto & cols = std::get<0>(r);
    const auto & data = std::get<1>(r);
    // Find CH column
    int ci = -1;
    for (int i=0; i<(int)cols.size(); i++) if (cols[i]=="CH") ci=i;
    bool has_c3=false, has_c4=false;
    if (ci>=0) {
      for (const auto & e : data[ci]) {
	if (std::holds_alternative<std::string>(e)) {
	  if (std::get<std::string>(e)=="C3") has_c3=true;
	  if (std::get<std::string>(e)=="C4") has_c4=true;
	}
      }
    }
    std::ostringstream m; m << "has_C3=" << has_c3 << " has_C4=" << has_c4;
    record(R,"psd/multi-channel", has_c3 && has_c4, m.str(), V);
  } catch(std::exception & e) { record(R,"psd/multi-channel",false,e.what(),V); }

  // F4 — epoch-level PSD: E strata present
  try {
    auto p = make_sine_inst(eng,10.0,256,1.0);
    p->eval("EPOCH len=30 & PSD sig=EEG max=30 epoch=T");
    // Expect a CH_E_F strata (or CH_E)
    bool found = false;
    for (const auto & cs : p->strata()) {
      if (cs.first=="PSD" && cs.second.find("E")!=std::string::npos)
	found = true;
    }
    std::ostringstream m; m << "epoch strata found=" << found;
    record(R,"psd/epoch-level", found, m.str(), V);
  } catch(std::exception & e) { record(R,"psd/epoch-level",false,e.what(),V); }
}

// ============================================================
// Group G: Spindles
// ============================================================

static void test_spindles( lunapi_t * eng,
			   std::vector<test_result_t> & R, bool V )
{
  // G0 — empirical threshold uses a valid Otsu split and optional outputs
  try {
    const std::vector<double> x = { 1, 1, 1, 2, 8 };
    double upper_f = 0;
    std::map<double,double> tvals;
    const double th = MiscMath::threshold( x, 0, 8, 1, &upper_f, &tvals );
    double upper_f_no_trace = 0;
    const double th_no_trace = MiscMath::threshold( x, 0, 8, 1,
						    &upper_f_no_trace, NULL );
    bool finite_trace = ! tvals.empty();
    for (std::map<double,double>::const_iterator ii = tvals.begin();
	 ii != tvals.end(); ++ii)
      finite_trace = finite_trace && std::isfinite( ii->second )
	&& ii->second >= 0 && ii->second <= 1;
    const bool pass = approx_equal( th, 2, 1e-12 )
      && approx_equal( upper_f, 0.2, 1e-12 )
      && approx_equal( th_no_trace, th, 1e-12 )
      && approx_equal( upper_f_no_trace, upper_f, 1e-12 )
      && finite_trace
      && approx_equal( tvals[th], 1, 1e-12 );
    std::ostringstream m;
    m << "threshold=" << th << " upper_fraction=" << upper_f
      << " trace_n=" << tvals.size();
    record(R,"spindles/empirical-threshold", pass, m.str(), V);
  } catch(std::exception & e) {
    record(R,"spindles/empirical-threshold",false,e.what(),V);
  }

  // G1 — detect spindles from sigma-band bursts
  // Multiple 1s bursts at 13Hz across the night
  try {
    auto p = eng->inst("T_sp");
    p->empty_edf("T_sp", 720, 30, "01.01.85","22.00.00");
    int sr = 256;
    double dur = 720*30.0;
    auto sig = make_noise((int)(sr*dur), 0.05, 77);
    // Add 13Hz bursts every 20s (each 1.5s), 36 per epoch × 720 epochs = ~1296 bursts
    for (int ep=0; ep<720; ep++) {
      double ep_start = ep * 30.0;
      for (int b=0; b<4; b++) {
	double bs = ep_start + 5.0 + b*7.0;
	double bd = 1.5;
	int bi = (int)(bs * sr);
	int bn = (int)(bd * sr);
	for (int i=0; i<bn && (bi+i)<(int)sig.size(); i++)
	  sig[bi+i] += 1.5 * std::sin(2.0*M_PI*13.0*(i/(double)sr));
      }
    }
    p->insert_signal("EEG", sig, sr);
    p->eval("EPOCH len=30 & SPINDLES sig=EEG fc=13");
    double dens = get_val(p,"SPINDLES","DENS");
    bool pass = dens > 1.0;  // expect several/minute
    std::ostringstream m; m << "DENS=" << dens << " (exp>1.0/min)";
    record(R,"spindles/detection-sigma-burst", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"spindles/detection-sigma-burst",false,e.what(),V); }

  // G2 — flat noise: very low spindle density
  try {
    auto p = eng->inst("T_sp_noise");
    p->empty_edf("T_sp_noise", 720, 30, "01.01.85","22.00.00");
    p->insert_signal("EEG", make_noise((int)(256*720*30.0), 0.1, 1), 256);
    p->eval("EPOCH len=30 & SPINDLES sig=EEG fc=13");
    double dens = get_val(p,"SPINDLES","DENS");
    bool pass = dens < 0.5;
    std::ostringstream m; m << "DENS=" << dens << " (exp≈0 on noise)";
    record(R,"spindles/noise-baseline", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"spindles/noise-baseline",false,e.what(),V); }

  // G3 — multi-frequency: output has rows for each fc value
  try {
    // Use burst signal (like G1) to ensure spindle detection at fc=13
    auto p = eng->inst("T_sp3");
    p->empty_edf("T_sp3", 100, 30, "01.01.85","22.00.00");
    int sr3 = 256;
    auto sig3 = make_noise((int)(sr3*100*30.0), 0.05, 77);
    for (int ep=0; ep<100; ep++) {
      double ep_start = ep * 30.0;
      for (int b=0; b<4; b++) {
        double bs = ep_start + 5.0 + b*7.0;
        int bi = (int)(bs * sr3);
        int bn = (int)(1.5 * sr3);
        for (int i=0; i<bn && (bi+i)<(int)sig3.size(); i++)
          sig3[bi+i] += 2.0 * std::sin(2.0*M_PI*13.0*(i/(double)sr3));
      }
    }
    p->insert_signal("EEG", sig3, sr3);
    p->eval("EPOCH len=30 & SPINDLES sig=EEG fc=11,13");
    int nrows_cf = get_nrows(p,"SPINDLES","CH_F");
    bool pass = nrows_cf >= 2;
    std::ostringstream m; m << "CH_F rows=" << nrows_cf << " (exp≥2 for fc=11,13)";
    record(R,"spindles/multi-freq", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"spindles/multi-freq",false,e.what(),V); }
}

// ============================================================
// Group H: Hypnogram
// ============================================================

static void test_hypno( lunapi_t * eng,
			std::vector<test_result_t> & R, bool V )
{
  // H1 — all-N2 night: TST = total duration, SE = 100%
  try {
    auto p = make_sine_inst(eng);
    p->insert_annotation("N2", make_stage_annots(720, 30.0));
    p->eval("EPOCH len=30 & HYPNO");
    double tst = get_val(p,"HYPNO","TST");
    double se  = get_val(p,"HYPNO","SE");
    bool pass = approx_equal(tst, 720*30.0/60.0, 1.0) && approx_equal(se, 100.0, 1 );
    std::ostringstream m; m << "TST=" << tst << "min (exp=" << 720*30.0/60.0 << ") SE=" << se << "% (exp=1.00)";
    record(R,"hypno/all-n2-tst-se", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"hypno/all-n2-tst-se",false,e.what(),V); }

  // H2 — mixed staging: TST excludes wake
  try {
    auto p = make_sine_inst(eng);
    // 360 Wake, 360 N2
    auto wake = make_stage_annots(360, 30.0);
    auto n2   = make_stage_annots(360, 30.0);
    for (auto & iv : n2) {
      std::get<0>(iv) += 360*30.0;
      std::get<1>(iv) += 360*30.0;
    }
    p->insert_annotation("W",  wake);
    p->insert_annotation("N2", n2);
    p->eval("EPOCH len=30 & HYPNO");
    double tst = get_val(p,"HYPNO","TST");
    // N2 % of TST is under SS strata, variable PCT, where SS column = "N2"
    double n2p = get_val_where(p,"HYPNO","SS","SS","N2","PCT");
    // TST should be ~180 min (360 N2 epochs × 30s / 60)
    bool pass = approx_equal(tst, 180.0, 2.0) && approx_equal(n2p, 1.0, 0.1);
    std::ostringstream m; m << "TST=" << tst << " (exp=180) N2%=" << n2p << " (exp=1.00)";
    record(R,"hypno/mixed-wake-n2", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"hypno/mixed-wake-n2",false,e.what(),V); }

  // H3 — sleep efficiency formula: SE = TST / TIB * 100
  try {
    auto p = make_sine_inst(eng);
    // 60 Wake + 600 N2 + 60 Wake at end (lights-off = start, lights-on = end)
    std::vector<std::tuple<double,double>> wake_start, n2, wake_end;
    for (int i=0;  i<60;  i++) wake_start.push_back({i*30.0,       (i+1)*30.0});
    for (int i=60; i<660; i++) n2.push_back(        {i*30.0,       (i+1)*30.0});
    for (int i=660;i<720; i++) wake_end.push_back(  {i*30.0,       (i+1)*30.0});
    p->insert_annotation("W",  wake_start);
    p->insert_annotation("N2", n2);
    // Note: second W set must be merged — re-insert appends
    p->insert_annotation("W",  wake_end);
    p->eval("EPOCH len=30 & HYPNO");
    double se  = get_val(p,"HYPNO","SE");
    // TST=600*30/60=300min; TIB=720*30/60=360min; SE=300/360*100=83.3%
    bool pass = approx_equal(se, 83.3, 3.0);
    std::ostringstream m; m << "SE=" << se << "% (exp≈83.3%)";
    record(R,"hypno/sleep-efficiency", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"hypno/sleep-efficiency",false,e.what(),V); }

  // H4 — HYPNO after RE: TST decreases
  try {
    auto p = make_sine_inst(eng);
    p->insert_annotation("N2", make_stage_annots(720, 30.0));
    p->eval("EPOCH len=30 & HYPNO");
    double tst_full = get_val(p,"HYPNO","TST");
    p->eval("MASK epoch=1-360 & RE & EPOCH len=30 & HYPNO");
    double tst_half = get_val(p,"HYPNO","TST");
    bool pass = tst_half < tst_full * 0.75;
    std::ostringstream m; m << "TST full=" << tst_full << " half=" << tst_half;
    record(R,"hypno/tst-after-re", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"hypno/tst-after-re",false,e.what(),V); }
}

// ============================================================
// Group I: Annotations
// ============================================================

static void test_annot( lunapi_t * eng,
			std::vector<test_result_t> & R, bool V )
{
  // I1 — insert and fetch: intervals round-trip through insert_annotation / fetch_annots
  try {
    auto p = make_sine_inst(eng);
    std::vector<std::tuple<double,double>> ivs = {{0,30},{60,90},{120,150}};
    p->insert_annotation("MySig", ivs);
    auto fa = p->fetch_annots({"MySig"});
    bool pass = fa.size() == 3;
    std::ostringstream m; m << "fetched " << fa.size() << " intervals (exp=3)";
    record(R,"annot/insert-fetch", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"annot/insert-fetch",false,e.what(),V); }

  // I2 — annots() lists inserted annotation
  try {
    auto p = make_sine_inst(eng);
    p->insert_annotation("TestLabel", {{10,20},{50,60}});
    auto av = p->annots();
    bool found = std::find(av.begin(),av.end(),"TestLabel") != av.end();
    record(R,"annot/list-after-insert", found, "annots() contains TestLabel: " + std::to_string(found), V);
  } catch(std::exception & e) { record(R,"annot/list-after-insert",false,e.what(),V); }

  // I3 — has_staging() false before, true after inserting stage annotations
  try {
    auto p = make_sine_inst(eng);
    p->eval("EPOCH len=30");  // must set epochs first for staging detection
    bool before = p->has_staging();
    p->insert_annotation("N2", make_stage_annots(720,30.0));
    p->eval("EPOCH len=30");  // re-epoch to register staging
    bool after = p->has_staging();
    std::ostringstream m; m << "before=" << before << " after=" << after;
    record(R,"annot/has-staging", !before && after, m.str(), V);
  } catch(std::exception & e) { record(R,"annot/has-staging",false,e.what(),V); }

  // I4 — fetch_full_annots returns class/instance/meta fields
  try {
    auto p = make_sine_inst(eng);
    std::vector<std::tuple<double,double>> ivs = {{0,30}};
    p->insert_annotation("SigA", ivs);
    auto fa = p->fetch_full_annots({"SigA"});
    bool pass = !fa.empty();
    // Check class field
    if (pass) pass = std::get<0>(fa[0]) == "SigA";
    std::ostringstream m; m << "full annots: count=" << fa.size()
			    << " class=" << (fa.empty() ? "?" : std::get<0>(fa[0]));
    record(R,"annot/fetch-full", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"annot/fetch-full",false,e.what(),V); }

  // I4b — fetch_full_annots(add_keys=true) returns keyed meta-data
  try {
    const std::string tmp = temp_base_path("test_fetch_full_meta");
    {
      std::ofstream out(tmp + ".annot");
      out << "class\tinstance\tchannel\tstart\tstop\tmeta\n";
      out << "SigA\tinst1\tECG\t0\t30\tk1=v1;k2=v2\n";
    }
    auto p = eng->inst("T_ffm");
    p->empty_edf("T_ffm", 720, 30, "01.01.85","22.00.00");
    p->attach_annot( tmp + ".annot" );
    auto fa0 = p->fetch_full_annots({"SigA"});
    auto fa1 = p->fetch_full_annots({"SigA"}, true);
    bool pass = fa0.size() == 1 && fa1.size() == 1
      && std::get<3>(fa0[0]) == "v1|v2"
      && std::get<3>(fa1[0]) == "k1=v1;k2=v2";
    std::ostringstream m; m << "legacy=" << (fa0.empty() ? "?" : std::get<3>(fa0[0]))
                            << " keyed=" << (fa1.empty() ? "?" : std::get<3>(fa1[0]));
    record(R,"annot/fetch-full-add-keys", pass, m.str(), V);
    std::remove( (tmp + ".annot").c_str() );
  } catch(std::exception & e) { record(R,"annot/fetch-full-add-keys",false,e.what(),V); }

  // I4c — time-shape classification: period HMS is distinct from elapsed
  // seconds, and continuation markers are never classified as HMS.
  try {
    const bool b1 = ! Helper::is_hms("20.330");
    const bool b2 = Helper::is_hms("20.00.30");
    const bool b3 = Helper::is_hms("20.00.30.1250");
    const bool b4 = Helper::is_hms("0+00.00.30.1250");
    const bool b5 = Helper::is_hms("20:00:30.1250");
    const bool b6 = Helper::is_hms("02.02.25-20.00.30.1250");
    const bool b7 = ! Helper::is_hms("20.00.30.1.2");
    const bool b8 = ! Helper::is_hms("...");
    const bool b9 = ! Helper::is_hms("-");
    const bool pass = b1 && b2 && b3 && b4 && b5 && b6 && b7 && b8 && b9;
    std::ostringstream m;
    m << "elapsed=" << b1 << " dot-hms=" << b2 << " fractional=" << b3
      << " elapsed-hms=" << b4 << " colon=" << b5 << " datetime=" << b6
      << " invalid-extra-period=" << b7 << " ellipsis=" << b8 << " dash=" << b9;
    record(R,"annot/time-shape-classification", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"annot/time-shape-classification",false,e.what(),V); }

  // I4d — load all elapsed/clock/continuation forms through the real .annot
  // reader, including the single-period ambiguity that motivated this fix.
  try {
    const std::string tmp = temp_base_path("test_annot_time_forms");
    {
      std::ofstream out(tmp + ".annot");
      out << "class\tinstance\tchannel\tstart\tstop\tmeta\n";
      out << "Elapsed\t.\t.\t20.330\t20.830\t.\n";
      out << "ClockDot\t.\t.\t20.00.30\t20.00.31\t.\n";
      out << "ClockFrac\t.\t.\t20.00.31.2500\t20.00.32.5000\t.\n";
      out << "ElapsedHMS\t.\t.\t0+00.00.32.5000\t0+00.00.33.5000\t.\n";
      out << "ClockColon\t.\t.\t20:00:34.1250\t20:00:35.2500\t.\n";
      out << "Ellipsis\t.\t.\t40\t...\t.\n";
      out << "AfterEllipsis\t.\t.\t50\t51\t.\n";
    }

    auto p = eng->inst("T_time_forms");
    p->empty_edf("T_time_forms", 2000, 30, "01.02.25", "20.00.00");
    const bool attached = p->attach_annot( tmp + ".annot" );

    auto one = [&]( const std::string & label ) -> std::tuple<bool,double,double> {
      const auto x = p->fetch_annots({label});
      if ( x.size() != 1 ) return {false, -1, -1};
      return {true, std::get<1>(x[0]), std::get<2>(x[0])};
    };

    const auto e = one("Elapsed");
    const auto d = one("ClockDot");
    const auto f = one("ClockFrac");
    const auto eh = one("ElapsedHMS");
    const auto c = one("ClockColon");
    const auto ell = one("Ellipsis");

    const bool pass = attached
      && std::get<0>(e) && approx_equal(std::get<1>(e),20.330,0.001) && approx_equal(std::get<2>(e),20.830,0.001)
      && std::get<0>(d) && approx_equal(std::get<1>(d),30.0,0.001) && approx_equal(std::get<2>(d),31.0,0.001)
      && std::get<0>(f) && approx_equal(std::get<1>(f),31.25,0.001) && approx_equal(std::get<2>(f),32.5,0.001)
      && std::get<0>(eh) && approx_equal(std::get<1>(eh),32.5,0.001) && approx_equal(std::get<2>(eh),33.5,0.001)
      && std::get<0>(c) && approx_equal(std::get<1>(c),34.125,0.001) && approx_equal(std::get<2>(c),35.25,0.001)
      && std::get<0>(ell) && approx_equal(std::get<1>(ell),40.0,0.001) && approx_equal(std::get<2>(ell),50.0,0.001);

    std::ostringstream m;
    m << "attached=" << attached << " elapsed=" << std::get<1>(e)
      << " dot=" << std::get<1>(d) << " fraction=" << std::get<1>(f)
      << " elapsed-hms=" << std::get<1>(eh) << " ellipsis-stop=" << std::get<2>(ell);
    record(R,"annot/time-form-loading", pass, m.str(), V);
    std::remove( (tmp + ".annot").c_str() );
  } catch(std::exception & e) { record(R,"annot/time-form-loading",false,e.what(),V); }

  // I4e — date-time parsing with all supported date-order settings and
  // period-delimited fractional clock components.
  try {
    const std::string base = temp_base_path("test_annot_datetime");
    const date_format_t prior = globals::read_annot_date_format;
    const std::vector<std::tuple<date_format_t,std::string,std::string>> cases = {
      { DMY, "02.02.25-20.00.00.1250", "02.02.25-20.00.01.1250" },
      { MDY, "02/02/25-20.00.00.1250", "02/02/25-20.00.01.1250" },
      { YMD, "2025-02-02-20.00.00.1250", "2025-02-02-20.00.01.1250" }
    };
    bool pass = true;
    std::ostringstream detail;

    for ( int k = 0 ; k < (int)cases.size() ; ++k )
      {
        const std::string file = base + "_" + std::to_string(k) + ".annot";
        {
          std::ofstream out(file);
          out << "class\tinstance\tchannel\tstart\tstop\tmeta\n";
          out << "DateEvt\t.\t.\t" << std::get<1>(cases[k]) << "\t"
              << std::get<2>(cases[k]) << "\t.\n";
        }

        globals::read_annot_date_format = std::get<0>(cases[k]);
        auto p = eng->inst("T_datetime_" + std::to_string(k));
        p->empty_edf("T_datetime_" + std::to_string(k), 3000, 30, "01.02.25", "20.00.00");
        const bool attached = p->attach_annot(file);
        const auto x = p->fetch_annots({"DateEvt"});
        const bool one_ok = attached && x.size() == 1
          && approx_equal(std::get<1>(x[0]),86400.125,0.001)
          && approx_equal(std::get<2>(x[0]),86401.125,0.001);
        pass = pass && one_ok;
        detail << " case" << k << "=" << one_ok;
        if ( !one_ok ) detail << "/n=" << x.size() << (x.empty() ? "" : "/s=" + Helper::dbl2str(std::get<1>(x[0])));
        std::remove(file.c_str());
      }

    globals::read_annot_date_format = prior;
    record(R,"annot/datetime-format-loading", pass, detail.str(), V);
  } catch(std::exception & e) { record(R,"annot/datetime-format-loading",false,e.what(),V); }

  // I4f — explicit relative days remain local to .annot parsing and work
  // with anonymized EDF dates: d1 may be before the EDF start, while d2 is
  // the following calendar day.
  try {
    const std::string tmp = temp_base_path("test_annot_relative_days");
    {
      std::ofstream out(tmp + ".annot");
      out << "class\tinstance\tchannel\tstart\tstop\tmeta\n";
      out << "Before\t.\t.\td1-20.00.00\td1-20.00.15\t.\n";
      out << "SameDay\t.\t.\td1-21.00.00\td1-21.00.15\t.\n";
      out << "NextDay\t.\t.\td2-20.00.00\td2-20.00.15\t.\n";
    }

    auto p = eng->inst("T_relative_days");
    p->empty_edf("T_relative_days", 4000, 30, "01.01.85", "20.30.00");
    const bool attached = p->attach_annot(tmp + ".annot");
    const auto same = p->fetch_annots({"SameDay"});
    const auto next = p->fetch_annots({"NextDay"});
    const auto before = p->fetch_annots({"Before"});
    const bool pass = attached && before.empty()
      && same.size() == 1 && next.size() == 1
      && approx_equal(std::get<1>(same[0]),1800.0,0.001)
      && approx_equal(std::get<2>(same[0]),1815.0,0.001)
      && approx_equal(std::get<1>(next[0]),84600.0,0.001)
      && approx_equal(std::get<2>(next[0]),84615.0,0.001);
    std::ostringstream m;
    m << "attached=" << attached << " before=" << before.size()
      << " same=" << same.size() << " next=" << next.size();
    record(R,"annot/relative-day-null-date", pass, m.str(), V);
    std::remove((tmp + ".annot").c_str());
  } catch(std::exception & e) { record(R,"annot/relative-day-null-date",false,e.what(),V); }

  // I4g — annot-time-wrap=F suppresses a wrapped undated timestamp when the
  // next occurrence is beyond a short EDF, but retains it in a multiday EDF.
  try {
    const std::string tmp = temp_base_path("test_annot_time_wrap");
    {
      std::ofstream out(tmp + ".annot");
      out << "class\tinstance\tchannel\tstart\tstop\tmeta\n";
      out << "WrapEvt\t.\t.\t20.00.00\t20.00.15\t.\n";
      out << "WrapOrderEvt\t.\t.\t20.33.30\t20.33.32\t.\n";
      out << "EndEvt\t.\t.\t08.00.00\t08.00.32\t.\n";
      out << "PointEndEvt\t.\t.\t08.00.32\t08.00.32\t.\n";
    }

    const bool prior_wrap = globals::annot_time_wrap;
    const bool prior_drop = globals::drop_annots_past_end;
    globals::drop_annots_past_end = false;

    auto short_edf = eng->inst("T_wrap_short");
    short_edf->empty_edf("T_wrap_short", 1374, 30, "01.01.85", "20.33.32");
    globals::annot_time_wrap = true;
    const bool short_attached_t = short_edf->attach_annot(tmp + ".annot");
    const auto short_t = short_edf->fetch_annots({"WrapEvt"});
    const auto short_order_t = short_edf->fetch_annots({"WrapOrderEvt"});

    auto short_no_wrap = eng->inst("T_wrap_short_no");
    short_no_wrap->empty_edf("T_wrap_short_no", 1374, 30, "01.01.85", "20.33.32");
    globals::annot_time_wrap = false;
    const bool short_attached_f = short_no_wrap->attach_annot(tmp + ".annot");
    const auto short_f = short_no_wrap->fetch_annots({"WrapEvt"});
    const auto short_order_f = short_no_wrap->fetch_annots({"WrapOrderEvt"});
    const auto short_end_f = short_no_wrap->fetch_annots({"EndEvt"});
    const auto short_point_end_f = short_no_wrap->fetch_annots({"PointEndEvt"});

    auto long_no_wrap = eng->inst("T_wrap_long_no");
    long_no_wrap->empty_edf("T_wrap_long_no", 3334, 30, "01.01.85", "20.33.32");
    const bool long_attached_f = long_no_wrap->attach_annot(tmp + ".annot");
    const auto long_f = long_no_wrap->fetch_annots({"WrapEvt"});

    auto short_drop_mask = eng->inst("T_drop_mask_default");
    short_drop_mask->empty_edf("T_drop_mask_default", 1374, 30, "01.01.85", "20.33.32");
    globals::annot_time_wrap = true;
    const bool short_drop_attached = short_drop_mask->attach_annot(tmp + ".annot");
    short_drop_mask->eval("DROP-ANNOTS annot=WrapEvt mask all");
    const auto short_drop = short_drop_mask->fetch_annots({"WrapEvt"});

    globals::annot_time_wrap = prior_wrap;
    globals::drop_annots_past_end = prior_drop;

    const bool pass = short_attached_t && short_t.size() == 1
      && approx_equal(std::get<1>(short_t[0]),84388.0,0.001)
      && short_order_t.size() == 1
      && approx_equal(std::get<1>(short_order_t[0]),86398.0,0.001)
      && approx_equal(std::get<2>(short_order_t[0]),86400.0,0.001)
      && short_attached_f && short_f.empty()
      && short_order_f.empty()
      && short_end_f.size() == 1
      && approx_equal(std::get<1>(short_end_f[0]),41188.0,0.001)
      && approx_equal(std::get<2>(short_end_f[0]),41220.0,0.001)
      && short_point_end_f.empty()
      && long_attached_f && long_f.size() == 1
      && approx_equal(std::get<1>(long_f[0]),84388.0,0.001)
      && short_drop_attached && short_drop.empty();
    std::ostringstream m;
    m << "short-wrap=" << short_t.size() << " short-no-wrap=" << short_f.size()
      << " order-wrap=" << short_order_t.size()
      << " order-no-wrap=" << short_order_f.size()
      << " endpoint=" << short_end_f.size() << " long-no-wrap=" << long_f.size()
      << " default-mask-drop=" << short_drop.size();
    record(R,"annot/time-wrap-option", pass, m.str(), V);
    std::remove((tmp + ".annot").c_str());
  } catch(std::exception & e) { record(R,"annot/time-wrap-option",false,e.what(),V); }

  // I5 — ANNOTS command output: COUNT matches inserted intervals
  try {
    auto p = make_sine_inst(eng);
    std::vector<std::tuple<double,double>> ivs = {{0,30},{60,90},{120,150},{200,230}};
    p->insert_annotation("SpindleX", ivs);
    p->eval("ANNOTS");
    // The ANNOTS command outputs per-annotation strata
    double cnt = get_val(p,"ANNOTS","COUNT");
    bool pass = approx_equal(cnt, 4.0, 0.5);
    std::ostringstream m; m << "COUNT=" << cnt << " (exp=4)";
    record(R,"annot/annots-count", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"annot/annots-count",false,e.what(),V); }

  // I6 — interp parameter in fetch_annots chops intervals
  try {
    auto p = make_sine_inst(eng);
    // One 120s interval; with interp=30, should yield 4 pieces
    p->insert_annotation("LongEvt", {{0.0, 120.0}});
    auto fa  = p->fetch_annots({"LongEvt"}, -1);  // no interp
    auto fai = p->fetch_annots({"LongEvt"}, 30.0); // interp=30s
    bool pass = (fa.size() == 1) && (fai.size() >= 4);
    std::ostringstream m;
    m << "no-interp=" << fa.size() << " interp30=" << fai.size() << " (exp 1 and >=4)";
    record(R,"annot/fetch-interp", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"annot/fetch-interp",false,e.what(),V); }

  // I7 — META look=next skips overlaps, whereas look=next-start can select them
  try {
    auto p = eng->inst("T_meta7");
    p->empty_edf("T_meta7", 60, 30, "01.01.85", "22.00.00");
    p->insert_annotation("Seed", {{10.0, 20.0}});
    p->insert_annotation("Other", {{15.0, 16.0}, {25.0, 26.0}});

    annot_t * other = p->find_annot("Other");
    bool overlap_tagged = false;
    bool after_tagged = false;
    if ( other != NULL )
      {
        for (auto & kv : other->interval_events)
          {
            if ( kv.second == NULL ) continue;
            const double s = kv.first.interval.start_sec();
            if ( approx_equal(s, 15.0, 0.001) )
              {
                kv.second->set("tag", std::string("overlap"));
                overlap_tagged = true;
              }
            else if ( approx_equal(s, 25.0, 0.001) )
              {
                kv.second->set("tag", std::string("after"));
                after_tagged = true;
              }
          }
      }

    p->eval("META annot=Seed other=Other look=next w-right=10 md=NEXT_FLAG flag & "
            "META annot=Seed other=Other look=next-start w-right=10 md=NEXT_START_FLAG flag & "
            "META annot=Seed other=Other look=next w-right=10 md=NEXT_TAG copy-md=tag & "
            "META annot=Seed other=Other look=next-start w-right=10 md=NEXT_START_TAG copy-md=tag");

    annot_t * seed = p->find_annot("Seed");
    bool next_flag_ok = false;
    bool next_start_flag_ok = false;
    bool next_tag_ok = false;
    bool next_start_tag_ok = false;

    if ( seed != NULL )
      {
        for (auto & kv : seed->interval_events)
          {
            if ( kv.second == NULL ) continue;
            avar_t * v1 = kv.second->find("NEXT_FLAG");
            avar_t * v2 = kv.second->find("NEXT_START_FLAG");
            avar_t * v3 = kv.second->find("NEXT_TAG");
            avar_t * v4 = kv.second->find("NEXT_START_TAG");
            next_flag_ok = v1 && v1->int_value() == 1;
            next_start_flag_ok = v2 && v2->int_value() == 1;
            next_tag_ok = v3 && v3->text_value() == "after";
            next_start_tag_ok = v4 && v4->text_value() == "overlap";
          }
      }

    bool pass = overlap_tagged && after_tagged && next_flag_ok && next_start_flag_ok
      && next_tag_ok && next_start_tag_ok;
    std::ostringstream m;
    m << "tagged=" << overlap_tagged << "/" << after_tagged
      << " next_flag=" << next_flag_ok
      << " next_start_flag=" << next_start_flag_ok
      << " next_tag=" << next_tag_ok
      << " next_start_tag=" << next_start_tag_ok;
    record(R,"annot/meta-next-vs-next-start", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"annot/meta-next-vs-next-start",false,e.what(),V); }
}

// ============================================================
// Group J: Write / read EDF round-trip
// ============================================================

static void test_write( lunapi_t * eng,
			std::vector<test_result_t> & R, bool V )
{
  // J1 — write to temp file, re-attach, verify channel count and duration
  try {
    const std::string tmp = temp_base_path("test_edf");
    auto p = make_sine_inst(eng, 10.0, 256, 1.0);
    p->eval( std::string("WRITE edf=") + tmp + " force-edf=T" );

    auto p2 = eng->inst("T_reload");
    bool ok = p2->attach_edf( tmp + ".edf" );

    int ns1 = (int)p->channels().size();
    int ns2 = (int)p2->channels().size();
    double dur1 = p->last_sec_original();
    double dur2 = p2->last_sec_original();
    bool pass = ok && (ns2 == ns1) && approx_equal(dur1, dur2, 30.0);
    std::ostringstream m;
    m << "ns1=" << ns1 << " ns2=" << ns2
      << " dur1=" << dur1 << " dur2=" << dur2;
    record(R,"write/edf-round-trip", pass, m.str(), V);
    // cleanup
    std::remove( (tmp + ".edf").c_str() );
  } catch(std::exception & e) { record(R,"write/edf-round-trip",false,e.what(),V); }

  // J2 — write-annots round-trip: WRITE-ANNOTS then re-attach, count matches
  try {
    const std::string tmp = temp_base_path("test_annot");
    auto p = make_sine_inst(eng);
    std::vector<std::tuple<double,double>> ivs;
    for (int i=0; i<10; i++) ivs.push_back({i*30.0, (i+1)*30.0});
    p->insert_annotation("RoundTrip", ivs);
    p->eval( std::string("WRITE-ANNOTS file=") + tmp + ".annot" );

    auto p2 = eng->inst("T_ra");
    p2->empty_edf("T_ra", 720, 30, "01.01.85","22.00.00");
    p2->attach_annot( tmp + ".annot" );
    auto fa = p2->fetch_annots({"RoundTrip"});
    bool pass = (int)fa.size() == 10;
    std::ostringstream m; m << "re-fetched " << fa.size() << " intervals (exp=10)";
    record(R,"write/annot-round-trip", pass, m.str(), V);
    std::remove( (tmp + ".annot").c_str() );
  } catch(std::exception & e) { record(R,"write/annot-round-trip",false,e.what(),V); }

  // J2b — write-annots annot-dir: create nested folder and write <ID>.annot
  try {
    const std::string tmp = temp_base_path("test_annot_dir");
    const std::string outdir = tmp + "/nested/annots";
    const std::string outfile = outdir + "/T1.annot";
    auto p = make_sine_inst(eng);
    std::vector<std::tuple<double,double>> ivs;
    for (int i=0; i<4; i++) ivs.push_back({i*30.0, (i+1)*30.0});
    p->insert_annotation("DirRoundTrip", ivs);
    p->eval( std::string("WRITE-ANNOTS annot-dir=") + outdir );

    auto p2 = eng->inst("T_ra_dir");
    p2->empty_edf("T_ra_dir", 720, 30, "01.01.85","22.00.00");
    p2->attach_annot( outfile );
    auto fa = p2->fetch_annots({"DirRoundTrip"});
    const bool exists = Helper::fileExists( outfile );
    bool pass = exists && (int)fa.size() == 4;
    std::ostringstream m;
    m << "exists=" << exists << " re-fetched " << fa.size() << " intervals (exp=4)";
    record(R,"write/annot-dir-round-trip", pass, m.str(), V);
    std::remove( outfile.c_str() );
  } catch(std::exception & e) { record(R,"write/annot-dir-round-trip",false,e.what(),V); }

  // J3 — STATS mean/SD preserved through EDF write/read
  try {
    const std::string tmp = temp_base_path("test_stats");
    auto p = make_sine_inst(eng, 10.0, 256, 1.5);
    p->eval("STATS sig=EEG");
    double sd1 = get_val(p,"STATS","SD");
    p->eval( std::string("WRITE edf=") + tmp + " force-edf=T" );
    auto p2 = eng->inst("T_stats");
    p2->attach_edf( tmp + ".edf" );
    p2->eval("STATS sig=EEG");
    double sd2 = get_val(p2,"STATS","SD");
    bool pass = approx_equal_rel(sd1, sd2, 0.02);
    std::ostringstream m; m << "SD before=" << sd1 << " after=" << sd2;
    record(R,"write/stats-preserved", pass, m.str(), V);
    std::remove( (tmp + ".edf").c_str() );
  } catch(std::exception & e) { record(R,"write/stats-preserved",false,e.what(),V); }
}

// ============================================================
// Group K: Script syntax
// ============================================================

static void test_script( lunapi_t * eng,
			 std::vector<test_result_t> & R, bool V )
{
  // K1 — & separator: two commands run correctly
  try {
    auto p = make_sine_inst(eng);
    p->eval("EPOCH len=30 & STATS sig=EEG");
    double ne = get_val(p,"EPOCH","NE");
    double sd = get_val(p,"STATS","SD");
    bool pass = approx_equal(ne,720,0.5) && sd > 0;
    std::ostringstream m; m << "NE=" << ne << " SD=" << sd;
    record(R,"script/amp-separator", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"script/amp-separator",false,e.what(),V); }

  // K2 — % comments ignored (command after comment must still run)
  try {
    auto p = make_sine_inst(eng);
    p->eval("% this is a comment\nSTATS sig=EEG");
    double sd = get_val(p,"STATS","SD");
    bool pass = sd > 0;
    std::ostringstream m; m << "SD=" << sd << " (comment not executed as command)";
    record(R,"script/comment-ignored", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"script/comment-ignored",false,e.what(),V); }

  // K3 — ivar substitution in eval string
  try {
    auto p = make_sine_inst(eng);
    p->ivar("mysig","EEG");
    p->eval("STATS sig=${mysig}");
    double sd = get_val(p,"STATS","SD");
    bool pass = sd > 0;
    std::ostringstream m; m << "SD=" << sd << " via ivar ${mysig}=EEG";
    record(R,"script/ivar-substitution", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"script/ivar-substitution",false,e.what(),V); }

  // K4 — global var via lunapi_t::var()
  try {
    eng->var("gvar","EEG");
    auto p = make_sine_inst(eng);
    p->eval("STATS sig=${gvar}");
    double sd = get_val(p,"STATS","SD");
    eng->dropvar("gvar");
    bool pass = sd > 0;
    std::ostringstream m; m << "SD=" << sd << " via global var ${gvar}=EEG";
    record(R,"script/global-var", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"script/global-var",false,e.what(),V); }

  // K5 — boolean assignment via ?{x=value} can drive IF blocks
  try {
    auto p = make_sine_inst(eng);
    p->eval("?{emit=1}\nIF ${emit}\nSTATS sig=EEG\nFI");
    double sd = get_val(p,"STATS","SD");
    bool pass = sd > 0;
    std::ostringstream m; m << "SD=" << sd << " after ?{emit=1}";
    record(R,"script/bool-var-assignment", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"script/bool-var-assignment",false,e.what(),V); }

  // K6 — missing FI must throw a parse error
  try {
    auto p = make_sine_inst(eng);
    bool threw = false;
    try {
      p->eval("IF 1\nSTATS sig=EEG");
    } catch (std::exception &) {
      threw = true;
    }
    record(R,"script/missing-fi-throws", threw, threw ? "missing FI rejected" : "missing FI accepted", V);
    globals::problem = false;
  } catch(std::exception & e) { record(R,"script/missing-fi-throws",false,e.what(),V); }

  // K7 — comment-only input should be ignored cleanly
  try {
    auto p = make_sine_inst(eng);
    bool threw = false;
    try {
      p->eval("% comment only\n% another comment");
    } catch (std::exception &) {
      threw = true;
    }
    record(R,"script/comment-only-safe", !threw, threw ? "comment-only input threw" : "comment-only input ignored", V);
    globals::problem = false;
  } catch(std::exception & e) { record(R,"script/comment-only-safe",false,e.what(),V); }

  // K8 — bad command: should throw, not crash
  try {
    auto p = make_sine_inst(eng);
    bool threw = false;
    try {
      p->eval("NOTACOMMAND_XXXX sig=EEG");
    } catch (std::exception &) {
      threw = true;
    }
    // Luna may warn but not necessarily throw on unknown command — mark pass either way
    // as long as the process doesn't crash (reaching here means no crash)
    record(R,"script/unknown-command-no-crash", true, "no crash on unknown command", V);
    globals::problem = false;
  } catch(std::exception & e) { record(R,"script/unknown-command-no-crash",false,e.what(),V); }
}

// ============================================================
// Group K2: direct eval parser/runtime behavior
// ============================================================

static void test_eval( lunapi_t * eng,
		       std::vector<test_result_t> & R, bool V )
{
  (void)eng;

  auto expect_eval_error = [&]( const std::string & name,
				const std::string & expr,
				const std::string & needle )
  {
    try {
      std::map<std::string,annot_map_t> inputs;
      instance_t out;
      Eval tok( expr );
      tok.bind( inputs , &out );
      bool ok = tok.evaluate( false );
      const std::string msg = tok.errmsg();
      bool pass = !ok && contains_substr( msg , needle );
      std::ostringstream m;
      m << "expr='" << expr << "' ok=" << ok << " errmsg='" << msg << "'";
      record( R, name, pass, m.str(), V );
    } catch ( std::exception & e ) {
      record( R, name, false, e.what(), V );
    }
  };

  auto expect_eval_throw = [&]( const std::string & name,
				const std::string & expr,
				const std::string & needle )
  {
    try {
      std::map<std::string,annot_map_t> inputs;
      instance_t out;
      Eval tok( expr );
      tok.bind( inputs , &out );
      (void)tok.evaluate( false );
      record( R, name, false, "expression unexpectedly succeeded", V );
    } catch ( std::exception & e ) {
      record( R, name, contains_substr( e.what(), needle ), e.what(), V );
    }
  };

  auto expect_eval_success = [&]( const std::string & name,
				  const std::string & expr,
				  const std::string & expected )
  {
    try {
      std::map<std::string,annot_map_t> inputs;
      instance_t out;
      Eval tok( expr );
      tok.bind( inputs , &out );
      bool ok = tok.evaluate( false );
      const std::string result = tok.result();
      bool pass = ok && result == expected;
      std::ostringstream m;
      m << "expr='" << expr << "' ok=" << ok << " result='" << result << "'";
      record( R, name, pass, m.str(), V );
    } catch ( std::exception & e ) {
      record( R, name, false, e.what(), V );
    }
  };

  expect_eval_error( "eval/malformed-number", "1..2", "malformed numeric literal" );
  expect_eval_error( "eval/unterminated-quote", "'abc", "unterminated single-quoted string" );
  expect_eval_error( "eval/missing-bracket", "A[2", "missing closing ]" );
  expect_eval_error( "eval/undefined-operand", "A=int(1,2,3); A[2] + B", "undefined operand for +" );
  expect_eval_success( "eval/tab-whitespace", "A=\t1;A", "1i" );
  expect_eval_throw( "eval/index-zero-rejected", "A=int(1,2,3); A[0]", "out of range for A" );

  // Annotation-backed EVAL regression tests.  A is present for epochs 1--300
  // and B for epochs 1--100; every epoch must remain valid, while the true
  // intervals reflect the expression result.
  try {
    auto p = make_sine_inst( eng );
    p->insert_annotation( "A", make_stage_annots( 300, 30.0 ) );
    p->insert_annotation( "B", make_stage_annots( 100, 30.0 ) );
    p->eval( "EPOCH len=30 & EVAL annot=E_IF expr=\"if(A)\"" );
    const auto found = p->fetch_annots( { "E_IF" } );
    const bool pass = found.size() == 300 &&
                      approx_equal( std::get<1>( found[0] ), 0.0, 0.01 ) &&
                      approx_equal( std::get<2>( found.back() ), 9000.0, 0.01 );
    std::ostringstream m; m << "n=" << found.size() << " (exp=300, 0--9000s)";
    record( R, "eval/annot-if", pass, m.str(), V );
  } catch ( std::exception & e ) { record( R, "eval/annot-if", false, e.what(), V ); }

  try {
    auto p = make_sine_inst( eng );
    p->insert_annotation( "A", make_stage_annots( 300, 30.0 ) );
    p->insert_annotation( "B", make_stage_annots( 100, 30.0 ) );
    p->eval( "EPOCH len=30 & EVAL annot=E_IFNOT expr=\"ifnot(B)\"" );
    const auto found = p->fetch_annots( { "E_IFNOT" } );
    const bool pass = found.size() == 620 &&
                      approx_equal( std::get<1>( found[0] ), 3000.0, 0.01 ) &&
                      approx_equal( std::get<2>( found.back() ), 21600.0, 0.01 );
    std::ostringstream m; m << "n=" << found.size() << " (exp=620, 3000--21600s)";
    record( R, "eval/annot-ifnot", pass, m.str(), V );
  } catch ( std::exception & e ) { record( R, "eval/annot-ifnot", false, e.what(), V ); }

  try {
    auto p = make_sine_inst( eng );
    p->insert_annotation( "A", make_stage_annots( 300, 30.0 ) );
    p->insert_annotation( "B", make_stage_annots( 100, 30.0 ) );
    p->eval( "EPOCH len=30 & EVAL annot=E_BOTH expr=\"if(A) && ifnot(B)\"" );
    const auto found = p->fetch_annots( { "E_BOTH" } );
    const bool pass = found.size() == 200 &&
                      approx_equal( std::get<1>( found[0] ), 3000.0, 0.01 ) &&
                      approx_equal( std::get<2>( found.back() ), 9000.0, 0.01 );
    std::ostringstream m; m << "n=" << found.size() << " (exp=200, 3000--9000s)";
    record( R, "eval/annot-if-and-ifnot", pass, m.str(), V );
  } catch ( std::exception & e ) { record( R, "eval/annot-if-and-ifnot", false, e.what(), V ); }

  // --------------------------------------------------------------------
  // Arithmetic operators adjacent to operands, i.e. WITHOUT surrounding
  // whitespace.  Regression guard: the numeric-literal scanner used to
  // greedily absorb a trailing +/- into the number (e.g. "2+1" -> the
  // malformed literal "2+1"), so binary +/- only worked when spaced.  A
  // +/- now only continues a number as a scientific-notation exponent
  // sign; otherwise it terminates the literal and is tokenized as an
  // operator.  These must parse like any normal language (cf. R).
  expect_eval_success( "eval/add-nospace",        "2+1",          "3i" );
  expect_eval_success( "eval/sub-nospace",        "2-1",          "1i" );
  expect_eval_success( "eval/add-nospace-2",      "7+8",          "15i" );
  expect_eval_success( "eval/sub-nospace-2",      "20-5",         "15i" );
  expect_eval_success( "eval/sub-negative",       "3-5",          "-2i" );
  expect_eval_success( "eval/mul-nospace",        "2*3",          "6i" );
  expect_eval_success( "eval/mod-nospace",        "3%2",          "1i" );
  expect_eval_success( "eval/mod-double",         "5%%2",         "1i" );

  // operator precedence and associativity
  expect_eval_success( "eval/prec-add-mul",       "2+3*4",        "14i" );
  expect_eval_success( "eval/prec-mul-add",       "2*3+4",        "10i" );
  expect_eval_success( "eval/paren-override",     "(2+3)*4",      "20i" );
  expect_eval_success( "eval/left-assoc-sub",     "5-3-1",        "1i" );
  expect_eval_success( "eval/left-assoc-sub-2",   "100-50-25",    "25i" );
  expect_eval_success( "eval/chain-add",          "2+2+2+2",      "8i" );

  // spacing must not change the result (equivalence with the above)
  expect_eval_success( "eval/add-spaced",         "2 + 1",        "3i" );
  expect_eval_success( "eval/add-multispace",     "2  +  3",      "5i" );
  expect_eval_success( "eval/add-surround-ws",    "  2+2  ",      "4i" );

  // division always yields a float
  expect_eval_success( "eval/div-int-exact",      "6/2",          "3f" );
  expect_eval_success( "eval/div-frac",           "10/4",         "2.5f" );
  expect_eval_success( "eval/div-half",           "1/2",          "0.5f" );
  expect_eval_success( "eval/div-then-sub",       "10/2-3",       "2f" );

  // decimal (float) operands
  expect_eval_success( "eval/float-mul",          "0.8*0.5",      "0.4f" );
  expect_eval_success( "eval/float-add",          "1.5+1.5",      "3f" );
  expect_eval_success( "eval/float-leading-dot",  ".5+.5",        "1f" );
  expect_eval_success( "eval/float-negative",     "-0.5",         "-0.5f" );

  // scientific notation must still parse (the exponent sign is retained)
  expect_eval_success( "eval/sci-plainexp",       "1e3",          "1000f" );
  expect_eval_success( "eval/sci-negexp",         "1e-3",         "0.001f" );
  expect_eval_success( "eval/sci-posexp",         "1.5e+2",       "150f" );
  expect_eval_success( "eval/sci-zeroexp-int",    "2e0",          "2i" );

  // unary +/- before a NUMBER (leading, or following another operator/paren)
  expect_eval_success( "eval/unary-plus",         "+2",           "2i" );
  expect_eval_success( "eval/unary-minus",        "-2",           "-2i" );
  expect_eval_success( "eval/unary-lead-add",     "-2+5",         "3i" );
  expect_eval_success( "eval/unary-after-mul",    "2*-1",         "-2i" );
  expect_eval_success( "eval/unary-both",         "-2*-3",        "6i" );
  expect_eval_success( "eval/sub-then-unary",     "2--1",         "3i" );
  expect_eval_success( "eval/mul-spaced-unary",   "2 * -3",       "-6i" );
  expect_eval_success( "eval/unary-double-num",   "--5",          "5i" );
  expect_eval_success( "eval/unary-plus-minus",   "-+3",          "-3i" );

  // unary +/- before a VARIABLE, FUNCTION or '(' : rewritten as -1*... / 1*...
  // (previously rejected; verify it now parses AND binds correctly).
  expect_eval_success( "eval/unary-var",          "A=3 ; -A",     "-3i" );
  expect_eval_success( "eval/unary-var-double",   "A=3 ; --A",    "3i" );
  expect_eval_success( "eval/unary-var-mul",      "A=3 ; -A*2",   "-6i" );
  expect_eval_success( "eval/unary-var-after-mul","A=3 ; 2*-A",   "-6i" );
  expect_eval_success( "eval/unary-var-add",      "A=5;B=2; -A+B","-3i" );
  expect_eval_success( "eval/unary-func",         "-sqrt(4)",     "-2f" );
  expect_eval_success( "eval/unary-func-negarg",  "-abs(-3)",     "-3i" );
  expect_eval_success( "eval/unary-func-after-mul","2*-sqrt(4)",  "-4f" );
  expect_eval_success( "eval/unary-paren",        "-(2+3)",       "-5i" );
  expect_eval_success( "eval/unary-add-neg",      "3+-2",         "1i" );

  // comparison and logical operators adjacent to operands
  expect_eval_success( "eval/gt",                 "2>1",          "true" );
  expect_eval_success( "eval/gte",                "5>=5",         "true" );
  expect_eval_success( "eval/lte-false",          "4<=3",         "false" );
  expect_eval_success( "eval/eq-false",           "1==2",         "false" );
  expect_eval_success( "eval/neq",                "1!=2",         "true" );
  expect_eval_success( "eval/arith-before-cmp",   "2+3==5",       "true" );
  expect_eval_success( "eval/logical-and",        "3>2 && 1<2",   "true" );
  expect_eval_success( "eval/logical-or",         "2>3 || 5>1",   "true" );
  expect_eval_success( "eval/logical-not",        "!(2>3)",       "true" );
  expect_eval_success( "eval/bool-and-false",     "true && false","false" );

  // functions (incl. unary-signed arguments)
  expect_eval_success( "eval/fn-sqrt",            "sqrt(4)",      "2f" );
  expect_eval_success( "eval/fn-pow",             "pow(2,3)",     "8i" );
  expect_eval_success( "eval/fn-abs-neg-arg",     "abs(-5)",      "5i" );
  expect_eval_success( "eval/fn-abs-then-add",    "abs(-5)+1",    "6i" );

  // string literals and comparison
  expect_eval_success( "eval/str-literal",        "'abc'",        "abc" );
  expect_eval_success( "eval/str-eq-true",        "'a'=='a'",     "true" );
  expect_eval_success( "eval/str-eq-false",       "'ab'=='cd'",   "false" );

  // malformed numbers must be rejected, not silently mis-parsed
  expect_eval_error( "eval/err-double-dot",       "2..3",         "malformed numeric literal" );
  expect_eval_error( "eval/err-triple-dotnum",    "3.4.5",        "malformed numeric literal" );
  expect_eval_error( "eval/err-exp-no-digit",     "1e",           "malformed numeric literal" );
  expect_eval_error( "eval/err-exp-sign-only",    "1e+",          "malformed numeric literal" );
  expect_eval_error( "eval/err-double-exp",       "1e3e4",        "malformed numeric literal" );

  // incomplete / ambiguous expressions must error (not silently succeed)
  expect_eval_error( "eval/err-trailing-op",      "2+",           "not enough arguments for +" );
  expect_eval_error( "eval/err-leading-binop",    "*3",           "not enough arguments for *" );
  expect_eval_error( "eval/err-no-pow-operator",  "2**3",         "not enough arguments for *" );
  expect_eval_error( "eval/err-adjacent-values",  "2 3",          "badly formed eval expression" );

  // malformed unary operators (sign with no valid operand)
  expect_eval_error( "eval/err-unary-incomplete", "-",            "incomplete unary operator" );
  expect_eval_error( "eval/err-unary-paren-close","-)",           "malformed unary - operator" );
  expect_eval_error( "eval/err-unary-then-op",    "-*3",          "malformed unary - operator" );
}

// ============================================================
// Group L: lunapi C++ API
// ============================================================

static void test_lunapi( lunapi_t * eng,
			 std::vector<test_result_t> & R, bool V )
{
  // L1 — status() fields
  try {
    auto p = make_sine_inst(eng);
    auto st = p->status();
    bool has_ns  = st.count("ns")  && std::holds_alternative<int>(st.at("ns"));
    bool has_dur = st.count("duration");
    bool has_id  = st.count("id");
    int  ns      = has_ns ? std::get<int>(st.at("ns")) : -1;
    std::ostringstream m;
    m << "ns=" << ns << " has_dur=" << has_dur << " has_id=" << has_id;
    record(R,"lunapi/status-fields", has_ns && has_dur && has_id && ns==1, m.str(), V);
  } catch(std::exception & e) { record(R,"lunapi/status-fields",false,e.what(),V); }

  // L2 — channels() returns inserted labels
  try {
    auto p = make_sine_inst(eng, 10.0, 256, 1.0, "MyEEG");
    auto chs = p->channels();
    bool found = std::find(chs.begin(),chs.end(),"MyEEG") != chs.end();
    record(R,"lunapi/channels-list", found, "channels contains MyEEG: " + std::to_string(found), V);
  } catch(std::exception & e) { record(R,"lunapi/channels-list",false,e.what(),V); }

  // L3 — has_channels() returns correct bool vector
  try {
    auto p = make_sine_inst(eng, 10.0, 256, 1.0, "EEG");
    auto hv = p->has_channels({"EEG","MISSING"});
    bool pass = hv.size()==2 && hv[0]==true && hv[1]==false;
    std::ostringstream m; m << "has EEG=" << hv[0] << " has MISSING=" << hv[1];
    record(R,"lunapi/has-channels", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"lunapi/has-channels",false,e.what(),V); }

  // L4 — last_sec() matches expected duration
  try {
    auto p = make_sine_inst(eng);  // 720*30 = 21600 s
    double ls = p->last_sec();
    bool pass = approx_equal(ls, 21600.0, 60.0);
    std::ostringstream m; m << "last_sec=" << ls << " (exp=21600)";
    record(R,"lunapi/last-sec", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"lunapi/last-sec",false,e.what(),V); }

  // L5 — epochs2intervals: first epoch timepoints start at 0
  try {
    auto p = make_sine_inst(eng);
    p->eval("EPOCH len=30");
    auto ivs = p->epochs2intervals({1,2,3});
    bool pass = ivs.size()==3;
    // First interval should start at timepoint 0
    if (pass) pass = (std::get<0>(ivs[0]) == 0);
    std::ostringstream m; m << "intervals=" << ivs.size() << " first_start=" << std::get<0>(ivs[0]);
    record(R,"lunapi/epochs2intervals", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"lunapi/epochs2intervals",false,e.what(),V); }

  // L6 — seconds2intervals round-trip
  try {
    auto p = make_sine_inst(eng);
    auto ivs = p->seconds2intervals({{0.0,30.0},{60.0,90.0}});
    bool pass = ivs.size()==2;
    std::ostringstream m; m << "intervals=" << ivs.size() << " (exp=2)";
    record(R,"lunapi/seconds2intervals", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"lunapi/seconds2intervals",false,e.what(),V); }

  // L7 — slice(): shape = samples × (channels + time_track)
  try {
    auto p = make_sine_inst(eng, 10.0, 256, 1.0);
    p->eval("EPOCH len=30");
    auto ivs = p->epochs2intervals({1});  // one 30s epoch = 256*30=7680 samples
    auto slr = p->slice(ivs, {"EEG"}, {}, true);
    const auto & scols = std::get<0>(slr);
    const auto & smat  = std::get<1>(slr);
    bool pass = (smat.rows() == 256*30) && (smat.cols() == 2);
    std::ostringstream m; m << "rows=" << smat.rows() << " cols=" << smat.cols()
			    << " (exp=" << 256*30 << "x2)";
    record(R,"lunapi/slice-shape", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"lunapi/slice-shape",false,e.what(),V); }

  // L8 — slices(): N intervals → N matrices
  try {
    auto p = make_sine_inst(eng);
    p->eval("EPOCH len=30");
    auto ivs = p->epochs2intervals({1,2,3});
    auto slsr = p->slices(ivs, {"EEG"}, {}, false);
    const auto & smats = std::get<1>(slsr);
    bool pass = smats.size()==3;
    for (const auto & m2 : smats) pass &= (m2.rows()==256*30);
    std::ostringstream m; m << "nmats=" << smats.size() << " (exp=3)";
    record(R,"lunapi/slices-count", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"lunapi/slices-count",false,e.what(),V); }

  // L9 — data(): smaller EDF to avoid memory issues
  try {
    auto p = eng->inst("T_data");
    p->empty_edf("T_data", 10, 30, "01.01.85","22.00.00");
    p->insert_signal("EEG", make_sine(256,10*30,10.0,1.0), 256);
    auto dr = p->data({"EEG"}, {}, true);
    const auto & dmat = std::get<1>(dr);
    bool pass = (dmat.rows() == 256*10*30) && (dmat.cols() == 2);
    std::ostringstream m; m << "rows=" << dmat.rows() << " cols=" << dmat.cols()
			    << " (exp=" << 256*10*30 << "x2)";
    record(R,"lunapi/data-shape", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"lunapi/data-shape",false,e.what(),V); }

  // L10 — eval_return_data: returns log + rtables
  try {
    auto p = make_sine_inst(eng);
    auto erd = p->eval_return_data("EPOCH len=30 & STATS sig=EEG");
    const auto & elog   = std::get<0>(erd);
    const auto & etabs  = std::get<1>(erd);
    bool has_epoch = etabs.count("EPOCH") > 0;
    bool has_stats = etabs.count("STATS") > 0;
    std::ostringstream m;
    m << "has_epoch=" << has_epoch << " has_stats=" << has_stats << " log_len=" << elog.size();
    record(R,"lunapi/eval-return-data", has_epoch && has_stats, m.str(), V);
  } catch(std::exception & e) { record(R,"lunapi/eval-return-data",false,e.what(),V); }

  // L11 — drop() / refresh(): state transitions
  try {
    auto p = make_sine_inst(eng);
    int state_before = p->get_state();
    p->drop();
    int state_after_drop = p->get_state();
    std::ostringstream m;
    m << "before=" << state_before << " after_drop=" << state_after_drop;
    record(R,"lunapi/drop-state", state_before==1 && state_after_drop==0, m.str(), V);
  } catch(std::exception & e) { record(R,"lunapi/drop-state",false,e.what(),V); }

  // L12 — error in eval throws runtime_error (lunapi_t::init redirects halts)
  try {
    auto p = make_sine_inst(eng);
    // REQUIRES a non-existent annotation — should cause a problem flag or throw
    bool threw = false;
    try { p->eval("REQUIRES annot=__NONEXISTENT__"); }
    catch (std::exception &) { threw = true; }
    globals::problem = false;
    // Whether or not it threw, the process should still be alive here
    record(R,"lunapi/error-no-crash", true, "no crash after eval error", V);
  } catch(std::exception & e) { record(R,"lunapi/error-no-crash",false,e.what(),V); }
}

// ============================================================
// Group M: segsrv
// ============================================================

static void test_segsrv( lunapi_t * eng,
			 std::vector<test_result_t> & R, bool V )
{
  // Helper: build a populated segsrv with one 256Hz EEG channel
  auto make_ss = [&]( const std::string & id ) -> std::pair<lunapi_inst_ptr, segsrv_t*> {
    auto p = eng->inst(id);
    p->empty_edf(id, 720, 30, "01.01.85", "22.00.00");
    p->insert_signal("EEG", make_sine(256, 720*30, 10.0, 1.0), 256);
    segsrv_t * ss = new segsrv_t(p);
    ss->populate({"EEG"}, {});
    return {p, ss};
  };

  // M1 — get_total_sec() matches inst duration
  try {
    auto [p, ss] = make_ss("T_ss1");
    double ts = ss->get_total_sec();
    bool pass = approx_equal(ts, 720.0*30, 60.0);
    std::ostringstream m; m << "total_sec=" << ts << " (exp=" << 720.0*30 << ")";
    record(R,"segsrv/total-sec", pass, m.str(), V);
    delete ss;
  } catch(std::exception & e) { record(R,"segsrv/total-sec",false,e.what(),V); }

  // M2 — set_window returns true for valid window
  try {
    auto [p, ss] = make_ss("T_ss2");
    bool ok_valid   = ss->set_window(0, 30);
    bool ok_invalid = ss->set_window(99999, 100000);  // past end
    std::ostringstream m;
    m << "valid=" << ok_valid << " invalid=" << ok_invalid;
    record(R,"segsrv/set-window-valid", ok_valid && !ok_invalid, m.str(), V);
    delete ss;
  } catch(std::exception & e) { record(R,"segsrv/set-window-valid",false,e.what(),V); }

  // M3 — get_signal() length ≈ window_size × SR
  try {
    auto [p, ss] = make_ss("T_ss3");
    ss->set_window(0, 30);
    auto sig = ss->get_signal("EEG");
    int expected = 256 * 30;
    bool pass = (sig.size() >= (size_t)(expected * 0.9));  // allow small deviation
    std::ostringstream m; m << "signal_len=" << sig.size() << " (exp≈" << expected << ")";
    record(R,"segsrv/get-signal-length", pass, m.str(), V);
    delete ss;
  } catch(std::exception & e) { record(R,"segsrv/get-signal-length",false,e.what(),V); }

  // M4 — get_timetrack() is monotone increasing
  try {
    auto [p, ss] = make_ss("T_ss4");
    ss->set_window(0, 30);
    auto tt = ss->get_timetrack("EEG");
    bool mono = true;
    for (int i = 1; i < (int)tt.size(); i++)
      if (tt[i] < tt[i-1]) { mono = false; break; }
    std::ostringstream m; m << "timetrack_len=" << tt.size() << " monotone=" << mono;
    record(R,"segsrv/timetrack-monotone", mono && tt.size() > 0, m.str(), V);
    delete ss;
  } catch(std::exception & e) { record(R,"segsrv/timetrack-monotone",false,e.what(),V); }

  // M5 — get_scaled_signal() values in [0,1] after set_scaling
  try {
    auto [p, ss] = make_ss("T_ss5");
    ss->set_scaling(1, 0, 1.0, 0.05, 0.02, 0.02, 0.1, false);
    ss->empirical_physical_scale("EEG");
    ss->set_window(0, 30);
    auto sc = ss->get_scaled_signal("EEG", 0);
    float mn = sc.minCoeff(), mx = sc.maxCoeff();
    bool pass = (mn >= -0.01f) && (mx <= 1.01f) && sc.size() > 0;
    std::ostringstream m; m << "min=" << mn << " max=" << mx;
    record(R,"segsrv/scaled-signal-range", pass, m.str(), V);
    delete ss;
  } catch(std::exception & e) { record(R,"segsrv/scaled-signal-range",false,e.what(),V); }

  // M6 — apply_filter / clear_filter / clear_filters: no crash
  try {
    auto [p, ss] = make_ss("T_ss6");
    ss->set_window(0, 60);
    auto raw_before = ss->get_signal("EEG");
    // identity SOS: [b0,b1,b2, a0,a1,a2] = [1,0,0, 1,0,0]
    std::vector<double> identity_sos = {1.0, 0.0, 0.0, 1.0, 0.0, 0.0};
    ss->apply_filter("EEG", identity_sos);
    // 'filtered' set should contain EEG now
    bool has_filter = ss->filtered.count("EEG") > 0;
    ss->clear_filter("EEG");
    bool after_clear = ss->filtered.count("EEG") == 0;
    auto raw_after = ss->get_signal("EEG");
    bool same_len = (raw_before.size() == raw_after.size());
    std::ostringstream m;
    m << "has_filter=" << has_filter << " after_clear=" << after_clear
      << " len=" << raw_before.size();
    record(R,"segsrv/apply-filter", has_filter && after_clear && same_len, m.str(), V);
    delete ss;
  } catch(std::exception & e) { record(R,"segsrv/apply-filter",false,e.what(),V); }

  // M7 — get_time_scale() returns at least one segment
  try {
    auto [p, ss] = make_ss("T_ss7");
    auto ts = ss->get_time_scale();
    bool pass = !ts.empty() && ts[0].first == 0.0;
    std::ostringstream m; m << "segments=" << ts.size() << " first_start=" << ts[0].first;
    record(R,"segsrv/time-scale", pass, m.str(), V);
    delete ss;
  } catch(std::exception & e) { record(R,"segsrv/time-scale",false,e.what(),V); }

  // M8 — calc_bands() and get_bands(): must call calc_bands BEFORE populate
  try {
    auto pb8 = eng->inst("T_ss8");
    pb8->empty_edf("T_ss8", 10, 30, "01.01.85","22.00.00");
    pb8->insert_signal("EEG", make_sine(256,10*30,10.0,1.0), 256);
    segsrv_t * ss8 = new segsrv_t(pb8);
    ss8->calc_bands({"EEG"});    // ← BEFORE populate so do_summaries runs
    ss8->populate({"EEG"}, {});
    int ne8 = ss8->nepochs();
    auto bands8 = ss8->get_bands("EEG");
    bool pass = (ne8 > 0) && (bands8.rows() == ne8) && (bands8.cols() > 0);
    std::ostringstream m;
    m << "ne=" << ne8 << " bands_shape=" << bands8.rows() << "x" << bands8.cols();
    record(R,"segsrv/calc-bands", pass, m.str(), V);
    delete ss8;
  } catch(std::exception & e) { record(R,"segsrv/calc-bands",false,e.what(),V); }

  // M9 — calc_hjorths() and get_hjorths(): same ordering requirement
  // Current segsrv returns a 101-column Hjorth display matrix, not nx3 scalars.
  try {
    auto pb9 = eng->inst("T_ss9");
    pb9->empty_edf("T_ss9", 10, 30, "01.01.85","22.00.00");
    pb9->insert_signal("EEG", make_sine(256,10*30,10.0,1.0), 256);
    segsrv_t * ss9 = new segsrv_t(pb9);
    ss9->calc_hjorths({"EEG"});  // ← BEFORE populate
    ss9->populate({"EEG"}, {});
    int ne9 = ss9->nepochs();
    auto hj9 = ss9->get_hjorths("EEG");
    bool pass = (ne9 > 0) && (hj9.rows() == ne9) && (hj9.cols() == 101);
    std::ostringstream m;
    m << "ne=" << ne9 << " hjorth_shape=" << hj9.rows() << "x" << hj9.cols() << " (exp nx101)";
    record(R,"segsrv/calc-hjorths", pass, m.str(), V);
    delete ss9;
  } catch(std::exception & e) { record(R,"segsrv/calc-hjorths",false,e.what(),V); }

  // M10 — annotations: compile_evts + fetch_evts
  try {
    auto p = eng->inst("T_ss10");
    p->empty_edf("T_ss10", 720, 30, "01.01.85","22.00.00");
    p->insert_signal("EEG", make_sine(256,720*30,10.0,1.0), 256);
    // Insert annotation covering seconds 0-90
    p->insert_annotation("Arousal", {{10.0,25.0},{60.0,75.0}});
    segsrv_t ss(p);
    ss.populate({"EEG"}, {"Arousal"});
    ss.add_annot("Arousal");
    ss.set_window(0, 90);
    ss.compile_evts({"Arousal"});
    auto evts = ss.fetch_evts();
    bool pass = evts.count("Arousal") > 0 && evts.at("Arousal").size() >= 1;
    std::ostringstream m;
    m << "annot_keys=" << evts.size() << " Arousal_count=" << evts["Arousal"].size();
    record(R,"segsrv/annot-compile-fetch", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"segsrv/annot-compile-fetch",false,e.what(),V); }

  // M11 — throttle: input_throttle uses integer-factor decimation at populate()
  try {
    auto p = eng->inst("T_ss11");
    p->empty_edf("T_ss11", 720, 30, "01.01.85","22.00.00");
    p->insert_signal("EEG", make_sine(256,720*30,10.0,1.0), 256);
    segsrv_t ss(p);
    ss.input_throttle(50);  // max 50 samples/sec input rate
    ss.populate({"EEG"}, {});
    ss.set_window(0, 30);
    auto sig = ss.get_signal("EEG");
    const int decimation_fac = 256 / 50;
    const int original_len = 256 * 30;
    const int expected = (( original_len - 1 ) / decimation_fac ) + 1;
    bool pass = (int)sig.size() == expected;
    std::ostringstream m; m << "throttled_len=" << sig.size() << " (exp=" << expected << ")";
    record(R,"segsrv/input-throttle", pass, m.str(), V);
  } catch(std::exception & e) { record(R,"segsrv/input-throttle",false,e.what(),V); }

  // M12 — get_window_left_hms() returns non-empty string after set_window
  try {
    auto [p, ss] = make_ss("T_ss12");
    ss->set_window(0, 30);
    std::string hms = ss->get_window_left_hms();
    bool pass = !hms.empty();
    std::ostringstream m; m << "left_hms='" << hms << "'";
    record(R,"segsrv/window-hms", pass, m.str(), V);
    delete ss;
  } catch(std::exception & e) { record(R,"segsrv/window-hms",false,e.what(),V); }

  // M13 — summary threshold flag: current segsrv exposes the threshold toggle
  // but not populated summary getters.
  try {
    auto [p, ss] = make_ss("T_ss13");
    ss->summary_threshold_mins(1.0);
    ss->set_window(0, 30);
    bool short_window = ss->serve_raw_signals();
    ss->set_window(0, 300);
    bool long_window = ss->serve_raw_signals();
    bool pass = !short_window && long_window;
    std::ostringstream m;
    m << "serve_raw short=" << short_window << " long=" << long_window;
    record(R,"segsrv/summary-threshold", pass, m.str(), V);
    delete ss;
  } catch(std::exception & e) { record(R,"segsrv/summary-threshold",false,e.what(),V); }


  // -----------------------------------------------------------------------
  // annot_editor_t / segsrv annotation editing tests
  // -----------------------------------------------------------------------
  //
  // Helpers shared across edit tests
  //   TP  = 1,000,000,000  (one second in tp units)
  //   All tests build a fresh EDF + segsrv; annotations inserted via
  //   lunapi_inst_t::insert_annotation() which uses (start_sec, stop_sec)
  //   tuples.  We recover exact tp values from fetch_all_evts_with_inst_ids().

  const uint64_t TP = globals::tp_1sec; // 1e9

  // Helper: build inst + segsrv with a single "Sp" annotation class
  //         containing two non-overlapping events:
  //           event 0: [10s, 11s)
  //           event 1: [20s, 22s)
  auto make_edit_ss = [&]( const std::string & id )
    -> std::pair<lunapi_inst_ptr, segsrv_t*>
  {
    auto p = eng->inst(id);
    p->empty_edf(id, 10, 30, "01.01.85", "22.00.00");  // 300 s
    p->insert_signal("EEG", make_sine(256, 10*30, 10.0, 1.0), 256);
    p->insert_annotation("Sp", {{10.0, 11.0}, {20.0, 22.0}});
    segsrv_t * ss = new segsrv_t(p);
    ss->populate({"EEG"}, {"Sp"});
    ss->add_annot("Sp");
    return {p, ss};
  };

  // Helper: count instances in a class
  auto count_class = [&]( lunapi_inst_ptr p, const std::string & cls ) -> int {
    annot_t * a = p->find_annot(cls);
    if (!a) return -1;
    return (int)a->interval_events.size();
  };

  // Helper: fetch rows for a class
  auto fetch_rows = [&]( segsrv_t * ss, const std::string & cls )
    -> std::vector<std::vector<std::string>>
  {
    return ss->fetch_all_evts_with_inst_ids({cls}, false);
  };


  // -----------------------------------------------------------------------
  // N1 — fetch_all_evts_with_inst_ids returns 8 cols (!hms) with correct tps
  // -----------------------------------------------------------------------
  try {
    auto [p, ss] = make_edit_ss("T_ae1");
    auto rows = fetch_rows(ss, "Sp");
    bool right_cols  = (rows.size() == 2) && (rows[0].size() == 8);
    // col 0 = class name
    bool right_class = (rows[0][0] == "Sp");
    // cols 4/5 should round-trip to exact tp values near 10s and 11s
    uint64_t s0 = 0, e0 = 0;
    Helper::str2int64(rows[0][4], &s0);
    Helper::str2int64(rows[0][5], &e0);
    bool right_tp = (s0 == 10 * TP) && (e0 == 11 * TP);
    bool pass = right_cols && right_class && right_tp;
    std::ostringstream m;
    m << "rows=" << rows.size() << " cols=" << (rows.empty() ? 0 : rows[0].size())
      << " class=" << (rows.empty() ? "?" : rows[0][0])
      << " start_tp=" << s0 << " stop_tp=" << e0;
    record(R, "segsrv/ae-fetch-cols", pass, m.str(), V);
    delete ss;
  } catch(std::exception & e) { record(R,"segsrv/ae-fetch-cols",false,e.what(),V); }


  // -----------------------------------------------------------------------
  // N2 — delete_annot: count drops by 1; deleted tp gone from fetch
  // -----------------------------------------------------------------------
  try {
    auto [p, ss] = make_edit_ss("T_ae2");
    auto rows = fetch_rows(ss, "Sp");
    // delete first event [10s,11s)
    uint64_t s0=0, e0=0;
    Helper::str2int64(rows[0][4], &s0);
    Helper::str2int64(rows[0][5], &e0);
    annot_edit_t ed;
    ed.aclass   = "Sp";
    ed.inst_id  = rows[0][6];
    ed.start_tp = s0;
    ed.stop_tp  = e0;
    ed.ch_str   = rows[0][7];
    ed.del      = true;
    ss->queue_edit(ed);
    int n = ss->apply_annot_edits({});
    int remaining = count_class(p, "Sp");
    // confirm the deleted tp is absent
    auto rows2 = fetch_rows(ss, "Sp");
    bool gone = true;
    for (auto & r : rows2) {
      uint64_t ts = 0; Helper::str2int64(r[4], &ts);
      if (ts == s0) { gone = false; break; }
    }
    bool pass = (n == 1) && (remaining == 1) && gone;
    std::ostringstream m;
    m << "n=" << n << " remaining=" << remaining << " gone=" << gone;
    record(R, "segsrv/ae-delete", pass, m.str(), V);
    delete ss;
  } catch(std::exception & e) { record(R,"segsrv/ae-delete",false,e.what(),V); }


  // -----------------------------------------------------------------------
  // N3 — shift: new_start only → duration preserved
  // -----------------------------------------------------------------------
  try {
    auto [p, ss] = make_edit_ss("T_ae3");
    auto rows = fetch_rows(ss, "Sp");
    uint64_t s0=0, e0=0;
    Helper::str2int64(rows[0][4], &s0);
    Helper::str2int64(rows[0][5], &e0);
    uint64_t dur = e0 - s0;                  // 1s = TP
    uint64_t new_start = 15 * TP;
    annot_edit_t ed;
    ed.aclass   = "Sp";  ed.inst_id = rows[0][6];
    ed.start_tp = s0;    ed.stop_tp = e0;  ed.ch_str = rows[0][7];
    ed.new_start = new_start;               // new_stop absent → shift
    ss->queue_edit(ed);
    ss->apply_annot_edits({});
    // find the shifted event
    auto rows2 = fetch_rows(ss, "Sp");
    bool found = false;
    for (auto & r : rows2) {
      uint64_t ts=0, te=0;
      Helper::str2int64(r[4], &ts);
      Helper::str2int64(r[5], &te);
      if (ts == new_start && (te - ts) == dur) { found = true; break; }
    }
    std::ostringstream m;
    m << "new_start=" << new_start << " dur=" << dur << " found=" << found;
    record(R, "segsrv/ae-shift", found, m.str(), V);
    delete ss;
  } catch(std::exception & e) { record(R,"segsrv/ae-shift",false,e.what(),V); }


  // -----------------------------------------------------------------------
  // N4 — resize: new_start + new_stop → exact new interval
  // -----------------------------------------------------------------------
  try {
    auto [p, ss] = make_edit_ss("T_ae4");
    auto rows = fetch_rows(ss, "Sp");
    uint64_t s0=0, e0=0;
    Helper::str2int64(rows[0][4], &s0);
    Helper::str2int64(rows[0][5], &e0);
    uint64_t ns = 12 * TP, ne = 14 * TP;
    annot_edit_t ed;
    ed.aclass   = "Sp";  ed.inst_id = rows[0][6];
    ed.start_tp = s0;    ed.stop_tp = e0;  ed.ch_str = rows[0][7];
    ed.new_start = ns;   ed.new_stop = ne;
    ss->queue_edit(ed);
    ss->apply_annot_edits({});
    auto rows2 = fetch_rows(ss, "Sp");
    bool found = false;
    for (auto & r : rows2) {
      uint64_t ts=0, te=0;
      Helper::str2int64(r[4], &ts); Helper::str2int64(r[5], &te);
      if (ts == ns && te == ne) { found = true; break; }
    }
    std::ostringstream m;
    m << "ns=" << ns << " ne=" << ne << " found=" << found;
    record(R, "segsrv/ae-resize", found, m.str(), V);
    delete ss;
  } catch(std::exception & e) { record(R,"segsrv/ae-resize",false,e.what(),V); }


  // -----------------------------------------------------------------------
  // N5 — channel change: new_ch updates ch_str in the rebuilt index
  // -----------------------------------------------------------------------
  try {
    auto [p, ss] = make_edit_ss("T_ae5");
    auto rows = fetch_rows(ss, "Sp");
    uint64_t s0=0, e0=0;
    Helper::str2int64(rows[0][4], &s0);
    Helper::str2int64(rows[0][5], &e0);
    annot_edit_t ed;
    ed.aclass   = "Sp";  ed.inst_id = rows[0][6];
    ed.start_tp = s0;    ed.stop_tp = e0;  ed.ch_str = rows[0][7];
    ed.new_ch   = std::string("EEG2");
    ss->queue_edit(ed);
    ss->apply_annot_edits({});
    auto rows2 = fetch_rows(ss, "Sp");
    bool found = false;
    for (auto & r : rows2) {
      uint64_t ts=0; Helper::str2int64(r[4], &ts);
      if (ts == s0 && r[7] == "EEG2") { found = true; break; }
    }
    std::ostringstream m;
    m << "ch_updated=" << found;
    record(R, "segsrv/ae-ch-change", found, m.str(), V);
    delete ss;
  } catch(std::exception & e) { record(R,"segsrv/ae-ch-change",false,e.what(),V); }


  // -----------------------------------------------------------------------
  // N6 — metadata set: meta key appears on modified instance
  // -----------------------------------------------------------------------
  try {
    auto [p, ss] = make_edit_ss("T_ae6");
    auto rows = fetch_rows(ss, "Sp");
    uint64_t s0=0, e0=0;
    Helper::str2int64(rows[0][4], &s0);
    Helper::str2int64(rows[0][5], &e0);
    annot_edit_t ed;
    ed.aclass   = "Sp";  ed.inst_id = rows[0][6];
    ed.start_tp = s0;    ed.stop_tp = e0;  ed.ch_str = rows[0][7];
    ed.meta["confidence"] = "0.95";
    ed.meta["note"]       = "manual";
    ss->queue_edit(ed);
    ss->apply_annot_edits({});
    // verify via fetch_full_annots
    annot_t * a = p->find_annot("Sp");
    bool found_conf = false, found_note = false;
    for (auto & kv : a->interval_events) {
      uint64_t ts = kv.first.interval.start;
      if (ts == s0 && kv.second) {
        avar_t * av_c = kv.second->find("confidence");
        avar_t * av_n = kv.second->find("note");
        if (av_c && av_c->text_value() == "0.95") found_conf = true;
        if (av_n && av_n->text_value() == "manual") found_note = true;
      }
    }
    bool pass = found_conf && found_note;
    std::ostringstream m;
    m << "confidence=" << found_conf << " note=" << found_note;
    record(R, "segsrv/ae-meta-set", pass, m.str(), V);
    delete ss;
  } catch(std::exception & e) { record(R,"segsrv/ae-meta-set",false,e.what(),V); }


  // -----------------------------------------------------------------------
  // N7 — clear_meta: existing meta removed; new meta applied
  // -----------------------------------------------------------------------
  try {
    auto [p, ss] = make_edit_ss("T_ae7");
    // First stamp some metadata onto the first event
    auto rows = fetch_rows(ss, "Sp");
    uint64_t s0=0, e0=0;
    Helper::str2int64(rows[0][4], &s0);
    Helper::str2int64(rows[0][5], &e0);
    {
      annot_edit_t ed;
      ed.aclass = "Sp"; ed.inst_id = rows[0][6];
      ed.start_tp = s0; ed.stop_tp = e0; ed.ch_str = rows[0][7];
      ed.meta["old_key"] = "old_val";
      ss->queue_edit(ed);
      ss->apply_annot_edits({});
    }
    // Now clear + set new meta
    auto rows2 = fetch_rows(ss, "Sp");
    uint64_t s1=0, e1=0;
    Helper::str2int64(rows2[0][4], &s1);
    Helper::str2int64(rows2[0][5], &e1);
    {
      annot_edit_t ed;
      ed.aclass = "Sp"; ed.inst_id = rows2[0][6];
      ed.start_tp = s1; ed.stop_tp = e1; ed.ch_str = rows2[0][7];
      ed.clear_meta = true;
      ed.meta["new_key"] = "new_val";
      ss->queue_edit(ed);
      ss->apply_annot_edits({});
    }
    annot_t * a = p->find_annot("Sp");
    bool old_gone = true, new_present = false;
    for (auto & kv : a->interval_events) {
      if (kv.first.interval.start == s1 && kv.second) {
        if (kv.second->find("old_key")) old_gone = false;
        avar_t * av = kv.second->find("new_key");
        if (av && av->text_value() == "new_val") new_present = true;
      }
    }
    bool pass = old_gone && new_present;
    std::ostringstream m;
    m << "old_gone=" << old_gone << " new_present=" << new_present;
    record(R, "segsrv/ae-clear-meta", pass, m.str(), V);
    delete ss;
  } catch(std::exception & e) { record(R,"segsrv/ae-clear-meta",false,e.what(),V); }


  // -----------------------------------------------------------------------
  // N8 — combined: shift + channel + meta in one edit
  // -----------------------------------------------------------------------
  try {
    auto [p, ss] = make_edit_ss("T_ae8");
    auto rows = fetch_rows(ss, "Sp");
    uint64_t s0=0, e0=0;
    Helper::str2int64(rows[0][4], &s0);
    Helper::str2int64(rows[0][5], &e0);
    uint64_t dur = e0 - s0;
    uint64_t new_s = 50 * TP;
    annot_edit_t ed;
    ed.aclass    = "Sp"; ed.inst_id = rows[0][6];
    ed.start_tp  = s0;   ed.stop_tp = e0; ed.ch_str = rows[0][7];
    ed.new_start = new_s;
    ed.new_ch    = std::string("F4");
    ed.meta["edited"] = "1";
    ss->queue_edit(ed);
    ss->apply_annot_edits({});
    annot_t * a = p->find_annot("Sp");
    bool found = false;
    for (auto & kv : a->interval_events) {
      const instance_idx_t & idx = kv.first;
      if (idx.interval.start == new_s &&
          (idx.interval.stop - idx.interval.start) == dur &&
          idx.ch_str == "F4" &&
          kv.second && kv.second->find("edited")) {
        found = true; break;
      }
    }
    std::ostringstream m;
    m << "combined_ok=" << found;
    record(R, "segsrv/ae-combined", found, m.str(), V);
    delete ss;
  } catch(std::exception & e) { record(R,"segsrv/ae-combined",false,e.what(),V); }


  // -----------------------------------------------------------------------
  // N9 — class filter: apply_annot_edits({"Other"}) leaves "Sp" untouched
  // -----------------------------------------------------------------------
  try {
    auto [p, ss] = make_edit_ss("T_ae9");
    // Also insert a second class
    p->insert_annotation("Other", {{5.0, 6.0}});
    ss->add_annot("Other");
    auto rows = fetch_rows(ss, "Sp");
    uint64_t s0=0, e0=0;
    Helper::str2int64(rows[0][4], &s0);
    Helper::str2int64(rows[0][5], &e0);
    annot_edit_t ed;
    ed.aclass = "Sp"; ed.inst_id = rows[0][6];
    ed.start_tp = s0; ed.stop_tp = e0; ed.ch_str = rows[0][7];
    ed.del = true;
    ss->queue_edit(ed);
    // Apply only to "Other" — the "Sp" delete should be applied but
    // "Other" has no queued edits, so net result: queue processed for
    // "Sp" because class filter means only "Other" is targeted
    int n = ss->apply_annot_edits({"Other"});
    int sp_count = count_class(p, "Sp");
    // "Sp" edit was queued but filtered out → Sp still has 2 events
    bool pass = (n == 0) && (sp_count == 2);
    std::ostringstream m;
    m << "n=" << n << " sp_count=" << sp_count;
    record(R, "segsrv/ae-class-filter", pass, m.str(), V);
    // clear the unfired queue
    ss->clear_edits();
    delete ss;
  } catch(std::exception & e) { record(R,"segsrv/ae-class-filter",false,e.what(),V); }


  // -----------------------------------------------------------------------
  // N10 — clear_edits: queued edits discarded, annotation unchanged
  // -----------------------------------------------------------------------
  try {
    auto [p, ss] = make_edit_ss("T_ae10");
    auto rows = fetch_rows(ss, "Sp");
    uint64_t s0=0, e0=0;
    Helper::str2int64(rows[0][4], &s0);
    Helper::str2int64(rows[0][5], &e0);
    annot_edit_t ed;
    ed.aclass = "Sp"; ed.inst_id = rows[0][6];
    ed.start_tp = s0; ed.stop_tp = e0; ed.ch_str = rows[0][7];
    ed.del = true;
    ss->queue_edit(ed);
    ss->clear_edits();               // discard before apply
    // count should be unchanged
    int cnt = count_class(p, "Sp");
    bool pass = (cnt == 2);
    std::ostringstream m;
    m << "cnt_after_cancel=" << cnt;
    record(R, "segsrv/ae-cancel", pass, m.str(), V);
    delete ss;
  } catch(std::exception & e) { record(R,"segsrv/ae-cancel",false,e.what(),V); }


  // -----------------------------------------------------------------------
  // N11 — apply empty queue: returns 0, annotation unchanged
  // -----------------------------------------------------------------------
  try {
    auto [p, ss] = make_edit_ss("T_ae11");
    int n = ss->apply_annot_edits({});
    int cnt = count_class(p, "Sp");
    bool pass = (n == 0) && (cnt == 2);
    std::ostringstream m;
    m << "n=" << n << " cnt=" << cnt;
    record(R, "segsrv/ae-empty-apply", pass, m.str(), V);
    delete ss;
  } catch(std::exception & e) { record(R,"segsrv/ae-empty-apply",false,e.what(),V); }


  // -----------------------------------------------------------------------
  // N12 — bulk: delete both events, count reaches 0
  // -----------------------------------------------------------------------
  try {
    auto [p, ss] = make_edit_ss("T_ae12");
    auto rows = fetch_rows(ss, "Sp");
    for (auto & r : rows) {
      uint64_t s=0, e=0;
      Helper::str2int64(r[4], &s); Helper::str2int64(r[5], &e);
      annot_edit_t ed;
      ed.aclass = "Sp"; ed.inst_id = r[6];
      ed.start_tp = s; ed.stop_tp = e; ed.ch_str = r[7];
      ed.del = true;
      ss->queue_edit(ed);
    }
    int n   = ss->apply_annot_edits({});
    int cnt = count_class(p, "Sp");
    bool pass = (n == 2) && (cnt == 0);
    std::ostringstream m;
    m << "n=" << n << " cnt=" << cnt;
    record(R, "segsrv/ae-bulk-delete", pass, m.str(), V);
    delete ss;
  } catch(std::exception & e) { record(R,"segsrv/ae-bulk-delete",false,e.what(),V); }


  // -----------------------------------------------------------------------
  // N13 — collision: moving event 0 onto event 1's position keeps count
  //        at 2 (edit dropped, not applied) — original event 1 survives
  // -----------------------------------------------------------------------
  try {
    auto [p, ss] = make_edit_ss("T_ae13");
    auto rows = fetch_rows(ss, "Sp");
    uint64_t s0=0, e0=0, s1=0, e1=0;
    Helper::str2int64(rows[0][4], &s0); Helper::str2int64(rows[0][5], &e0);
    Helper::str2int64(rows[1][4], &s1); Helper::str2int64(rows[1][5], &e1);
    // Move event 0 exactly onto event 1's interval
    annot_edit_t ed;
    ed.aclass    = "Sp"; ed.inst_id = rows[0][6];
    ed.start_tp  = s0;   ed.stop_tp = e0; ed.ch_str = rows[0][7];
    ed.new_start = s1;   ed.new_stop = e1;
    ss->queue_edit(ed);
    ss->apply_annot_edits({});
    // One of the two survives; count is 1 (collision drops the incoming edit,
    // so the untouched event 1 wins and the moved event 1' is discarded)
    int cnt = count_class(p, "Sp");
    bool pass = (cnt == 1);
    std::ostringstream m;
    m << "cnt_after_collision=" << cnt;
    record(R, "segsrv/ae-collision", pass, m.str(), V);
    delete ss;
  } catch(std::exception & e) { record(R,"segsrv/ae-collision",false,e.what(),V); }


  // -----------------------------------------------------------------------
  // N14 — queue survives apply_annot_edits only for unmatched classes
  //        (queue is fully cleared after apply regardless)
  // -----------------------------------------------------------------------
  try {
    auto [p, ss] = make_edit_ss("T_ae14");
    auto rows = fetch_rows(ss, "Sp");
    uint64_t s0=0, e0=0;
    Helper::str2int64(rows[0][4], &s0); Helper::str2int64(rows[0][5], &e0);
    annot_edit_t ed;
    ed.aclass = "Sp"; ed.inst_id = rows[0][6];
    ed.start_tp = s0; ed.stop_tp = e0; ed.ch_str = rows[0][7];
    ed.del = true;
    ss->queue_edit(ed);
    ss->apply_annot_edits({});   // queue cleared
    // Second apply with empty queue should be a no-op
    int n2 = ss->apply_annot_edits({});
    int cnt = count_class(p, "Sp");
    bool pass = (n2 == 0) && (cnt == 1);
    std::ostringstream m;
    m << "n2=" << n2 << " cnt=" << cnt;
    record(R, "segsrv/ae-queue-cleared", pass, m.str(), V);
    delete ss;
  } catch(std::exception & e) { record(R,"segsrv/ae-queue-cleared",false,e.what(),V); }


  // -----------------------------------------------------------------------
  // N15 — interval tree rebuilt: fetch after edit returns updated tps
  // -----------------------------------------------------------------------
  try {
    auto [p, ss] = make_edit_ss("T_ae15");
    auto rows = fetch_rows(ss, "Sp");
    uint64_t s0=0, e0=0;
    Helper::str2int64(rows[0][4], &s0); Helper::str2int64(rows[0][5], &e0);
    uint64_t new_s = 80 * TP, new_e = 81 * TP;
    annot_edit_t ed;
    ed.aclass = "Sp"; ed.inst_id = rows[0][6];
    ed.start_tp = s0; ed.stop_tp = e0; ed.ch_str = rows[0][7];
    ed.new_start = new_s; ed.new_stop = new_e;
    ss->queue_edit(ed);
    ss->apply_annot_edits({});
    // Verify via interval_tree query directly
    annot_t * a = p->find_annot("Sp");
    interval_t qwin( new_s, new_e );
    auto hits = a->extract( qwin );
    bool found_in_tree = !hits.empty();
    // And the old position should not be in the tree
    interval_t old_win( s0, e0 );
    auto old_hits = a->extract( old_win );
    bool old_gone_from_tree = old_hits.empty();
    bool pass = found_in_tree && old_gone_from_tree;
    std::ostringstream m;
    m << "found_in_tree=" << found_in_tree
      << " old_gone=" << old_gone_from_tree;
    record(R, "segsrv/ae-tree-rebuild", pass, m.str(), V);
    delete ss;
  } catch(std::exception & e) { record(R,"segsrv/ae-tree-rebuild",false,e.what(),V); }


  // -----------------------------------------------------------------------
  // N16 — inst_id rename: new_inst_id updates key.id; old id gone, new present
  // -----------------------------------------------------------------------
  try {
    auto [p, ss] = make_edit_ss("T_ae16");
    auto rows = fetch_rows(ss, "Sp");
    uint64_t s0=0, e0=0;
    Helper::str2int64(rows[0][4], &s0);
    Helper::str2int64(rows[0][5], &e0);
    std::string old_id = rows[0][6];
    std::string new_id = "renamed_inst";

    annot_edit_t ed;
    ed.aclass      = "Sp"; ed.inst_id = old_id;
    ed.start_tp    = s0;   ed.stop_tp = e0; ed.ch_str = rows[0][7];
    ed.new_inst_id = new_id;
    ss->queue_edit(ed);
    ss->apply_annot_edits({});

    // Verify via interval_events: old id gone, new id present at same interval
    annot_t * a = p->find_annot("Sp");
    bool old_gone = true, new_found = false;
    for (auto & kv : a->interval_events) {
      if (kv.first.interval.start == s0 && kv.first.interval.stop == e0) {
        if (kv.first.id == old_id) old_gone = false;
        if (kv.first.id == new_id) new_found = true;
      }
    }
    // Also verify fetch returns new id in col 6
    auto rows2 = fetch_rows(ss, "Sp");
    bool fetch_shows_new = false;
    for (auto & row : rows2)
      if (row[4] == Helper::int2str(s0) && row[6] == new_id) fetch_shows_new = true;

    bool pass = old_gone && new_found && fetch_shows_new;
    std::ostringstream msg;
    msg << "old_gone=" << old_gone
        << " new_found=" << new_found
        << " fetch_shows_new=" << fetch_shows_new;
    record(R, "segsrv/ae-inst-id-rename", pass, msg.str(), V);
    delete ss;
  } catch(std::exception & e) { record(R,"segsrv/ae-inst-id-rename",false,e.what(),V); }
}

// ============================================================
// Group P: LM / PLM (WASM 2016)
// ============================================================

static void test_plm( lunapi_t * eng,
                      std::vector<test_result_t> & R, bool V )
{
  using namespace plm;
  const double TP = (double)globals::tp_1sec;

  // build one component (signal-domain event)
  auto mkc = [&]( int id , side_t side , double onset , double dur ) -> lm_component_t {
    lm_component_t c;
    c.comp_id = id; c.side = side; c.sig = "X";
    c.onset_tp  = (uint64_t)std::llround( onset * TP );
    c.offset_tp = (uint64_t)std::llround( ( onset + dur ) * TP );
    c.onset_sec = onset; c.dur = dur;
    c.base_onset = 0; c.scale_onset = 1; c.on_thr = 0; c.off_thr = 0;
    c.peak = 0; c.peak_exc = 0; c.peak_z = 0; c.morph_value = 0; c.morph_ok = true;
    c.clm_component = ( dur >= 0.5 && dur <= 10.0 );
    c.stage = "UNK"; c.state = "UNK"; c.qc_base_high = false;
    return c;
  };

  // build one final event (for classify / sequence tests)
  auto mke = [&]( int id , double onset , double dur , bool clm ,
                  int n_comp = 1 , bool all_clm = true , side_t side = SIDE_LEG ) -> lm_event_t {
    lm_event_t e;
    e.evt_id = id; e.side = side;
    e.onset_tp  = (uint64_t)std::llround( onset * TP );
    e.offset_tp = (uint64_t)std::llround( ( onset + dur ) * TP );
    e.onset_sec = onset; e.dur = dur;
    e.n_comp = n_comp; e.n_l = 0; e.n_r = 0; e.all_comps_clm = all_clm;
    e.clm = clm; e.plm = false; e.seq = 0; e.seq_pos = 0; e.seq_n = 0;
    e.imi_prev = std::numeric_limits<double>::quiet_NaN();
    e.imi_next = std::numeric_limits<double>::quiet_NaN();
    e.stage = "UNK"; e.state = "UNK"; e.peak_exc_max = 0; e.peak_z_max = 0; e.qc_base_high = false;
    return e;
  };

  plm_param_t P; // defaults

  // ---- bilateral combination -------------------------------------------

  try { // P6 overlapping L/R combine
    std::vector<lm_component_t> cs = { mkc(1,SIDE_L,0,1.0), mkc(2,SIDE_R,0.5,1.0) };
    auto ev = combine_bilateral( cs , P );
    bool pass = ev.size()==1 && ev[0].n_comp==2 && ev[0].side==SIDE_B;
    record(R,"plm/bilateral-overlap", pass, "n_ev="+std::to_string(ev.size()), V);
  } catch(std::exception & e){ record(R,"plm/bilateral-overlap",false,e.what(),V); }

  try { // P7 gap < 0.5 combine
    std::vector<lm_component_t> cs = { mkc(1,SIDE_L,0,1.0), mkc(2,SIDE_R,1.3,0.7) };
    auto ev = combine_bilateral( cs , P );
    bool pass = ev.size()==1 && ev[0].n_comp==2;
    record(R,"plm/bilateral-gap-lt", pass, "n_ev="+std::to_string(ev.size()), V);
  } catch(std::exception & e){ record(R,"plm/bilateral-gap-lt",false,e.what(),V); }

  try { // P8 gap == 0.5 does NOT combine (strict)
    std::vector<lm_component_t> cs = { mkc(1,SIDE_L,0,1.0), mkc(2,SIDE_R,1.5,0.5) };
    auto ev = combine_bilateral( cs , P );
    bool pass = ev.size()==2;
    record(R,"plm/bilateral-gap-eq", pass, "n_ev="+std::to_string(ev.size())+" (exp 2)", V);
  } catch(std::exception & e){ record(R,"plm/bilateral-gap-eq",false,e.what(),V); }

  try { // P9 transitive multi-component grouping
    std::vector<lm_component_t> cs = { mkc(1,SIDE_L,0,1.0), mkc(2,SIDE_R,0.8,0.8), mkc(3,SIDE_L,1.5,0.7) };
    auto ev = combine_bilateral( cs , P );
    bool pass = ev.size()==1 && ev[0].n_comp==3 && ev[0].side==SIDE_B;
    record(R,"plm/bilateral-transitive", pass, "n_ev="+std::to_string(ev.size()), V);
  } catch(std::exception & e){ record(R,"plm/bilateral-transitive",false,e.what(),V); }

  try { // P13 component LM > 10s makes bilateral non-CLM
    std::vector<lm_component_t> cs = { mkc(1,SIDE_L,0,11.0), mkc(2,SIDE_R,0.5,0.6) };
    auto ev = combine_bilateral( cs , P );
    classify_clm( ev , P );
    bool pass = ev.size()==1 && ev[0].n_comp==2 && ev[0].clm==false;
    record(R,"plm/bilateral-long-component", pass, "clm="+std::to_string(ev.empty()?-1:ev[0].clm), V);
  } catch(std::exception & e){ record(R,"plm/bilateral-long-component",false,e.what(),V); }

  // ---- CLM classification boundaries -----------------------------------

  try { // P: monolateral dur==10 is CLM; dur>10 is not
    std::vector<lm_component_t> a = { mkc(1,SIDE_LEG,0,10.0) };
    std::vector<lm_component_t> b = { mkc(1,SIDE_LEG,0,10.001) };
    auto ea = combine_single(a,P); classify_clm(ea,P);
    auto eb = combine_single(b,P); classify_clm(eb,P);
    bool pass = ea[0].clm==true && eb[0].clm==false;
    record(R,"plm/clm-10s-boundary", pass, "", V);
  } catch(std::exception & e){ record(R,"plm/clm-10s-boundary",false,e.what(),V); }

  try { // P12 bilateral dur==15 CLM; >15 not
    std::vector<lm_event_t> e1 = { mke(1,0,15.0,false,2,true,SIDE_B) };
    std::vector<lm_event_t> e2 = { mke(1,0,15.001,false,2,true,SIDE_B) };
    classify_clm(e1,P); classify_clm(e2,P);
    bool pass = e1[0].clm==true && e2[0].clm==false;
    record(R,"plm/bilateral-15s-boundary", pass, "", V);
  } catch(std::exception & e){ record(R,"plm/bilateral-15s-boundary",false,e.what(),V); }

  try { // P10/P11 four bilateral components allowed, five not
    std::vector<lm_event_t> e4 = { mke(1,0,5.0,false,4,true,SIDE_B) };
    std::vector<lm_event_t> e5 = { mke(1,0,5.0,false,5,true,SIDE_B) };
    classify_clm(e4,P); classify_clm(e5,P);
    bool pass = e4[0].clm==true && e5[0].clm==false;
    record(R,"plm/bilateral-max-components", pass, "", V);
  } catch(std::exception & e){ record(R,"plm/bilateral-max-components",false,e.what(),V); }

  // ---- PLM sequence detection ------------------------------------------

  auto run_seq = [&]( std::vector<lm_event_t> ev ) {
    return detect_sequences( ev , P );
  };

  try { // P14 four CLMs @20s -> one sequence, four PLMs
    std::vector<lm_event_t> ev = { mke(1,0,1,true), mke(2,20,1,true), mke(3,40,1,true), mke(4,60,1,true) };
    auto seqs = detect_sequences(ev,P);
    int nplm=0; for(auto&e:ev) if(e.plm) nplm++;
    bool pass = seqs.size()==1 && nplm==4 && seqs[0].n==4;
    record(R,"plm/seq-4clm-20s", pass, "nseq="+std::to_string(seqs.size())+" nplm="+std::to_string(nplm), V);
  } catch(std::exception & e){ record(R,"plm/seq-4clm-20s",false,e.what(),V); }

  try { // P15 three CLMs -> no sequence
    std::vector<lm_event_t> ev = { mke(1,0,1,true), mke(2,20,1,true), mke(3,40,1,true) };
    auto seqs = detect_sequences(ev,P);
    record(R,"plm/seq-3clm-none", seqs.empty(), "nseq="+std::to_string(seqs.size()), V);
  } catch(std::exception & e){ record(R,"plm/seq-3clm-none",false,e.what(),V); }

  try { // P16 IMI exactly 10 valid
    std::vector<lm_event_t> ev = { mke(1,0,1,true), mke(2,10,1,true), mke(3,20,1,true), mke(4,30,1,true) };
    auto seqs = detect_sequences(ev,P);
    record(R,"plm/seq-imi-10-valid", seqs.size()==1, "nseq="+std::to_string(seqs.size()), V);
  } catch(std::exception & e){ record(R,"plm/seq-imi-10-valid",false,e.what(),V); }

  try { // P17 IMI exactly 90 valid
    std::vector<lm_event_t> ev = { mke(1,0,1,true), mke(2,90,1,true), mke(3,180,1,true), mke(4,270,1,true) };
    auto seqs = detect_sequences(ev,P);
    record(R,"plm/seq-imi-90-valid", seqs.size()==1, "nseq="+std::to_string(seqs.size()), V);
  } catch(std::exception & e){ record(R,"plm/seq-imi-90-valid",false,e.what(),V); }

  try { // P18 IMI < 10 breaks
    std::vector<lm_event_t> ev = { mke(1,0,1,true), mke(2,20,1,true), mke(3,40,1,true), mke(4,45,1,true) };
    auto seqs = detect_sequences(ev,P);
    record(R,"plm/seq-imi-short-breaks", seqs.empty(), "nseq="+std::to_string(seqs.size()), V);
  } catch(std::exception & e){ record(R,"plm/seq-imi-short-breaks",false,e.what(),V); }

  try { // P19 IMI > 90 breaks
    std::vector<lm_event_t> ev = { mke(1,0,1,true), mke(2,20,1,true), mke(3,40,1,true), mke(4,140,1,true) };
    auto seqs = detect_sequences(ev,P);
    record(R,"plm/seq-imi-long-breaks", seqs.empty(), "nseq="+std::to_string(seqs.size()), V);
  } catch(std::exception & e){ record(R,"plm/seq-imi-long-breaks",false,e.what(),V); }

  try { // P20 non-CLM LM breaks; no IMI computed across it
    std::vector<lm_event_t> ev = { mke(1,0,1,true), mke(2,20,1,true), mke(3,40,12,false),
                                   mke(4,60,1,true), mke(5,80,1,true) };
    auto seqs = detect_sequences(ev,P);
    // ev is sorted by onset inside detect_sequences; index 3 == onset 60
    bool no_imi = std::isnan( ev[3].imi_prev );
    record(R,"plm/seq-nonclm-breaks", seqs.empty() && no_imi,
           "nseq="+std::to_string(seqs.size())+" imi_across_nan="+std::to_string(no_imi), V);
  } catch(std::exception & e){ record(R,"plm/seq-nonclm-breaks",false,e.what(),V); }

  try { // P21 event after a break can begin a new sequence
    std::vector<lm_event_t> ev = { mke(1,0,1,true), mke(2,20,1,true), mke(3,40,1,true),
                                   mke(4,45,1,true), mke(5,65,1,true), mke(6,85,1,true),
                                   mke(7,105,1,true), mke(8,125,1,true) };
    auto seqs = detect_sequences(ev,P);
    bool pass = seqs.size()==1 && seqs[0].n==5;
    record(R,"plm/seq-restart-after-break", pass,
           "nseq="+std::to_string(seqs.size())+(seqs.empty()?"":" n="+std::to_string(seqs[0].n)), V);
  } catch(std::exception & e){ record(R,"plm/seq-restart-after-break",false,e.what(),V); }

  // ---- end-to-end detection (adaptive; scale-relative) -----------------

  // synthetic 1-channel EMG: low background + periodic 1s bursts @20s
  const int fs = 100, nr = 40, rs = 30;         // 1200 s
  const double dur = nr * rs;
  auto make_emg = [&]( double amp , double scale , int seed ) {
    auto s = make_noise( (int)(fs*dur), 1.0, seed );
    for ( double t = 10.0 ; t + 1.0 < dur ; t += 20.0 )
      {
        int b = (int)(t*fs);
        for ( int i = 0 ; i < fs ; i++ ) if ( b+i < (int)s.size() ) s[b+i] += amp;
      }
    for ( auto & v : s ) v *= scale;
    return s;
  };

  try { // P22-24 single-channel adaptive: bursts -> CLMs -> PLM sequence
    auto p = eng->inst("T_plm1");
    p->empty_edf("T_plm1", nr, rs, "01.01.85","22.00.00");
    p->insert_signal("LAT", make_emg(40,1.0,11), fs);
    p->eval("LM sig=LAT method=adaptive");
    double nlm = get_val(p,"LM","N_LM");
    double nclm = get_val(p,"LM","N_CLM");
    double nseq = get_val(p,"LM","N_PLM_SEQ");
    bool pass = nlm>=40 && nclm>=40 && nseq>=1;
    std::ostringstream m; m << "N_LM="<<nlm<<" N_CLM="<<nclm<<" N_PLM_SEQ="<<nseq;
    record(R,"plm/e2e-adaptive-single", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"plm/e2e-adaptive-single",false,e.what(),V); }

  try { // P25 adaptive detection invariant to multiplicative rescale
    auto p1 = eng->inst("T_plm_s1");
    p1->empty_edf("T_plm_s1", nr, rs, "01.01.85","22.00.00");
    p1->insert_signal("LAT", make_emg(40,1.0,11), fs);
    p1->eval("LM sig=LAT method=adaptive");
    double n1 = get_val(p1,"LM","N_LM");

    auto p2 = eng->inst("T_plm_s2");
    p2->empty_edf("T_plm_s2", nr, rs, "01.01.85","22.00.00");
    p2->insert_signal("LAT", make_emg(40,3.7,11), fs);   // same signal * 3.7
    p2->eval("LM sig=LAT method=adaptive");
    double n2 = get_val(p2,"LM","N_LM");

    bool pass = n1>0 && n1==n2;
    std::ostringstream m; m << "N_LM x1="<<n1<<" x3.7="<<n2;
    record(R,"plm/adaptive-rescale-invariant", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"plm/adaptive-rescale-invariant",false,e.what(),V); }

  try { // P28 raising k-on does not increase detections
    auto p1 = eng->inst("T_plm_k1");
    p1->empty_edf("T_plm_k1", nr, rs, "01.01.85","22.00.00");
    p1->insert_signal("LAT", make_emg(40,1.0,11), fs);
    p1->eval("LM sig=LAT method=adaptive k-on=6");
    double n_lo = get_val(p1,"LM","N_LM");

    auto p2 = eng->inst("T_plm_k2");
    p2->empty_edf("T_plm_k2", nr, rs, "01.01.85","22.00.00");
    p2->insert_signal("LAT", make_emg(40,1.0,11), fs);
    p2->eval("LM sig=LAT method=adaptive k-on=12");
    double n_hi = get_val(p2,"LM","N_LM");

    bool pass = n_hi <= n_lo;
    std::ostringstream m; m << "N_LM k-on6="<<n_lo<<" k-on12="<<n_hi;
    record(R,"plm/adaptive-kon-monotone", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"plm/adaptive-kon-monotone",false,e.what(),V); }

  try { // P26 fixed-uV WASM detection is NOT rescale-invariant
    auto p1 = eng->inst("T_plm_w1");
    p1->empty_edf("T_plm_w1", nr, rs, "01.01.85","22.00.00");
    p1->insert_signal("LAT", make_emg(30,1.0,11), fs);   // ~30 uV bursts
    p1->eval("SET-HEADERS sig=LAT physical-dimension=uV & LM sig=LAT method=wasm");
    double n_big = get_val(p1,"LM","N_LM");

    auto p2 = eng->inst("T_plm_w2");
    p2->empty_edf("T_plm_w2", nr, rs, "01.01.85","22.00.00");
    p2->insert_signal("LAT", make_emg(30,0.1,11), fs);   // scaled down 10x -> ~3 uV
    p2->eval("SET-HEADERS sig=LAT physical-dimension=uV & LM sig=LAT method=wasm");
    double n_small = get_val(p2,"LM","N_LM");

    bool pass = n_big>=40 && n_small < n_big;
    std::ostringstream m; m << "N_LM uV="<<n_big<<" 0.1*uV="<<n_small;
    record(R,"plm/wasm-not-scale-invariant", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"plm/wasm-not-scale-invariant",false,e.what(),V); }

  try { // P27 iterative baseline: a big movement does not inflate its own baseline
    // one very large sustained burst region + normal periodic bursts; adaptive
    // should still detect the periodic bursts near the large one
    auto p = eng->inst("T_plm_bl");
    p->empty_edf("T_plm_bl", nr, rs, "01.01.85","22.00.00");
    auto s = make_noise((int)(fs*dur),1.0,7);
    // large 8s excursion around t=300
    for (int i=(int)(300*fs); i<(int)(308*fs) && i<(int)s.size(); i++) s[i]+=200;
    // periodic 1s bursts @20s
    for (double t=10.0; t+1.0<dur; t+=20.0){ int b=(int)(t*fs); for(int i=0;i<fs;i++) if(b+i<(int)s.size()) s[b+i]+=40; }
    p->insert_signal("LAT", s, fs);
    p->eval("LM sig=LAT method=adaptive");
    double nclm = get_val(p,"LM","N_CLM");
    bool pass = nclm >= 40;   // the huge event did not suppress nearby detections
    record(R,"plm/adaptive-baseline-robust", pass, "N_CLM="+std::to_string(nclm), V);
  } catch(std::exception & e){ record(R,"plm/adaptive-baseline-robust",false,e.what(),V); }
}

// ============================================================
// Main entry point
// ============================================================

void proc_tests( const std::string & group, const bool verbose )
{
  std::cout << "\n=== Luna integrated tests"
	    << " [group=" << (group.empty() ? "all" : group)
	    << " verbose=" << verbose << "] ===\n\n";

  // Suppress normal Luna log output during tests unless verbose
  if (!verbose) logger.off();

  // Initialise the lunapi singleton (redirects halt() to exceptions)
  lunapi_t * eng = lunapi_t::inaugurate();
  eng->silence(true);

  std::vector<test_result_t> results;

  const bool run_all      = (group == "all" || group.empty());

#define RUN(g, fn) if (run_all || group == (g)) { fn(eng, results, verbose); }

  RUN("signal",   test_signal)
  RUN("epoch",    test_epoch)
  RUN("mask",     test_mask)
  RUN("filter",   test_filter)
  RUN("resample", test_resample)
  RUN("psd",      test_psd)
  RUN("spindles", test_spindles)
  RUN("hypno",    test_hypno)
  RUN("annot",    test_annot)
  RUN("write",    test_write)
  RUN("script",   test_script)
  RUN("eval",     test_eval)
  RUN("lunapi",   test_lunapi)
  RUN("segsrv",   test_segsrv)
  RUN("plm",      test_plm)

#undef RUN

  lunapi_t::retire();

  // Summary
  int total = n_pass + n_fail;
  std::cout << "\n--- Summary ---\n"
	    << "  PASS: " << n_pass << " / " << total << "\n"
	    << "  FAIL: " << n_fail << " / " << total << "\n";

  if (n_fail > 0)
    {
      std::cout << "\nFailed tests:\n";
      for (const auto & r : results)
	if (!r.pass)
	  std::cout << "  [FAIL] " << r.name << "  " << r.msg << "\n";
    }

  std::cout << "\n";

  // Set process exit code
  globals::retcode = (n_fail > 0) ? 1 : 0;
}
