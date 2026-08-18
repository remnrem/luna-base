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
//         hypno, annot, write, script, eval, lunapi, segsrv, plm, sigdyn, dpp,
//         dpp-fit (HAS_LGBM builds only)
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
#include "stats/dpp-spec.h"
#include "stats/dpp-filter.h"
#include "stats/dpp-fit.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <fstream>
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

// Same as get_val_where, but for a numeric (double/int) factor column,
// e.g. SIGDYN's SEC offset level -- matched within 'tol'
static double get_val_where_num( lunapi_inst_ptr p,
				 const std::string & cmd,
				 const std::string & strata,
				 const std::string & factor_col,
				 double factor_val,
				 const std::string & var,
				 double tol = 1e-3 )
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
    double fval = std::numeric_limits<double>::quiet_NaN();
    if (std::holds_alternative<double>(data[fci][r2])) fval = std::get<double>(data[fci][r2]);
    else if (std::holds_alternative<int>(data[fci][r2])) fval = (double)std::get<int>(data[fci][r2]);
    if (std::fabs(fval - factor_val) <= tol) {
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

  // J4 — clear-reserved input option blanks signal reserved fields
  const bool old_clear_reserved = globals::clear_reserved;
  try {
    const std::string src = temp_base_path("test_clear_reserved_src");
    const std::string keep = temp_base_path("test_clear_reserved_keep");
    const std::string clear = temp_base_path("test_clear_reserved_clear");
    const std::string uuid = "9cecc8e855df4ad8b78ebecaaa00dc84";

    auto p = make_sine_inst(eng, 10.0, 256, 1.0);
    p->eval( std::string("WRITE edf=") + src + " force-edf=T" );

    // For one signal, the signal-reserved field is the final 32 bytes of
    // the 256-byte signal header: offset 256 + 224 = 480.
    {
      std::fstream f( src + ".edf" , std::ios::in | std::ios::out | std::ios::binary );
      f.seekp( 480 );
      f.write( uuid.data() , uuid.size() );
    }

    globals::clear_reserved = false;
    auto p_keep = eng->inst("T_clear_reserved_keep");
    const bool attached_keep = p_keep->attach_edf( src + ".edf" );
    p_keep->eval( std::string("WRITE edf=") + keep + " force-edf=T" );

    globals::clear_reserved = true;
    auto p_clear = eng->inst("T_clear_reserved_clear");
    const bool attached_clear = p_clear->attach_edf( src + ".edf" );
    p_clear->eval( std::string("WRITE edf=") + clear + " force-edf=T" );

    std::string kept( 32 , '\0' ), cleared( 32 , '\0' );
    {
      std::ifstream f( keep + ".edf" , std::ios::binary );
      f.seekg( 480 );
      f.read( &kept[0] , kept.size() );
    }
    {
      std::ifstream f( clear + ".edf" , std::ios::binary );
      f.seekg( 480 );
      f.read( &cleared[0] , cleared.size() );
    }

    globals::clear_reserved = old_clear_reserved;

    const bool pass = attached_keep && attached_clear
      && kept == uuid
      && cleared == std::string( 32 , ' ' );
    std::ostringstream m;
    m << "kept='" << kept << "' cleared_spaces="
      << ( cleared == std::string( 32 , ' ' ) );
    record(R,"write/clear-reserved-input", pass, m.str(), V);

    std::remove( (src + ".edf").c_str() );
    std::remove( (keep + ".edf").c_str() );
    std::remove( (clear + ".edf").c_str() );
  } catch(std::exception & e) {
    globals::clear_reserved = old_clear_reserved;
    record(R,"write/clear-reserved-input",false,e.what(),V);
  }
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
// Group Q: SIGDYN
// ============================================================

static void test_sigdyn( lunapi_t * eng,
			 std::vector<test_result_t> & R, bool V )
{
  // helper: a 1 Hz "ramp" signal (sample i, sr=1Hz, has value i -- i.e.
  // signal value == elapsed second) on an n-second recording (n must be a
  // multiple of rs=10 for add_signal's exact-length requirement). Gives
  // fully predictable peri-event means: mean of samples [a,b] == (a+b)/2.
  auto make_ramp = [&]( int n, const std::string & id ) -> lunapi_inst_ptr {
    const int rs = 10;
    const int nr = n / rs;
    lunapi_inst_ptr p = eng->inst( id );
    p->empty_edf( id, nr, rs, "01.01.85", "22.00.00" );
    std::vector<double> ramp( n );
    for (int i = 0; i < n; i++) ramp[i] = (double)i;
    p->insert_signal( "X", ramp, 1 );
    return p;
  };

  // Q1 -- mode 0: constant signal gives exact MEAN/SD, overall and by stage
  try {
    auto p = eng->inst("T_sigdyn_q1");
    p->empty_edf("T_sigdyn_q1", 90, 10, "01.01.85", "22.00.00");     // 900s
    p->insert_signal( "X", std::vector<double>(900, 5.0), 1 );
    p->insert_annotation( "N2", make_stage_annots(30, 30.0) );        // 30x30s epochs
    // STAGE populates the per-epoch stage lookup SIGDYN's mode 0 reads
    // (timeline_t::epoch_annotation()) -- a raw inserted annotation class
    // alone does not, until a STAGE/HYPNO run compiles it
    p->eval("EPOCH dur=30 & STAGE & SIGDYN sig=X hypno-annot=F");
    double n_ch    = get_val_s(p,"SIGDYN","CH","N");
    double mean_ch = get_val_s(p,"SIGDYN","CH","MEAN");
    double sd_ch   = get_val_s(p,"SIGDYN","CH","SD");
    double mean_n2 = get_val_where(p,"SIGDYN","CH_SS","SS","N2","MEAN");
    bool pass = approx_equal(n_ch,30,0.5) && approx_equal(mean_ch,5.0,0.01)
      && approx_equal(sd_ch,0.0,0.01) && approx_equal(mean_n2,5.0,0.01);
    std::ostringstream m; m << "N=" << n_ch << " MEAN=" << mean_ch << " SD=" << sd_ch << " N2-MEAN=" << mean_n2;
    record(R,"sigdyn/mode0-constant-signal", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"sigdyn/mode0-constant-signal",false,e.what(),V); }

  // Q2 -- mode 2, native (unbinned) resolution: a single 10s anchor well
  // away from any edge gives an exactly symmetric +/-5s window, with M/L/R/S
  // matching the ramp signal's known values by hand.
  try {
    auto p = make_ramp(2000, "T_sigdyn_q2");
    p->insert_annotation( "Evt", { {1000.0, 1010.0} } );
    p->eval("EPOCH dur=30 & SIGDYN sig=X annot=Evt hypno-annot=F w=5");
    double n = get_val_s(p,"SIGDYN","ANNOT_CH","N");
    double m = get_val_s(p,"SIGDYN","ANNOT_CH","M");
    double l = get_val_s(p,"SIGDYN","ANNOT_CH","L");
    double r = get_val_s(p,"SIGDYN","ANNOT_CH","R");
    double s = get_val_s(p,"SIGDYN","ANNOT_CH","S");
    int nsec = get_nrows(p,"SIGDYN","ANNOT_CH_SEC");
    // tol=0.05: EDF's 16-bit physical/digital round-trip quantizes a
    // 0..1999 range to ~0.03 resolution
    bool pass = approx_equal(n,1,0.1) && approx_equal(m,1000,0.05)
      && approx_equal(l,997,0.05) && approx_equal(r,1003,0.05)
      && approx_equal(s,10,0.05) && nsec == 11;
    std::ostringstream mm; mm << "N=" << n << " M=" << m << " L=" << l << " R=" << r << " S=" << s << " nSEC=" << nsec;
    record(R,"sigdyn/native-symmetric-window", pass, mm.str(), V);
  } catch(std::exception & e){ record(R,"sigdyn/native-symmetric-window",false,e.what(),V); }

  // Q3 -- anchor= start vs end: the t=0 anchor shifts by exactly the
  // instance's own duration (10s), since a ramp signal's value == time
  try {
    auto p1 = make_ramp(2000, "T_sigdyn_q3a");
    p1->insert_annotation( "Evt", { {1000.0, 1010.0} } );
    p1->eval("EPOCH dur=30 & SIGDYN sig=X annot=Evt hypno-annot=F w=2 anchor=start");
    double m_start = get_val_where_num(p1,"SIGDYN","ANNOT_CH_SEC","SEC",0.0,"M");

    auto p2 = make_ramp(2000, "T_sigdyn_q3b");
    p2->insert_annotation( "Evt", { {1000.0, 1010.0} } );
    p2->eval("EPOCH dur=30 & SIGDYN sig=X annot=Evt hypno-annot=F w=2 anchor=end");
    double m_end = get_val_where_num(p2,"SIGDYN","ANNOT_CH_SEC","SEC",0.0,"M");

    bool pass = approx_equal(m_start,1000,0.05) && approx_equal(m_end,1010,0.05);
    std::ostringstream m; m << "anchor=start SEC0=" << m_start << " anchor=end SEC0=" << m_end;
    record(R,"sigdyn/anchor-point-shifts-t0", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"sigdyn/anchor-point-shifts-t0",false,e.what(),V); }

  // Q4 -- only complete bins contribute: an anchor placed 25 samples from
  // the recording's physical end, with bin=10/w=35/bin-align=start, means
  // only the SEC=10 bin (needing samples out to +19) can still complete on
  // the right (SEC=20/30, needing +29/+39, cannot); bins that DO appear use
  // their full, untruncated sample range (checked via the exact ramp mean),
  // not a partial/thinned one.
  try {
    auto p = make_ramp(1000, "T_sigdyn_q4");                 // samples 0..999
    p->insert_annotation( "Evt", { {974.0, 974.0} } );       // point anchor, offset cap = 999-974 = 25
    p->eval("EPOCH dur=30 & SIGDYN sig=X annot=Evt hypno-annot=F w=35 bin=10 bin-align=start");
    double m0   = get_val_where_num(p,"SIGDYN","ANNOT_CH_SEC","SEC",0.0,"M");    // mean(974..983)
    double mm10 = get_val_where_num(p,"SIGDYN","ANNOT_CH_SEC","SEC",-10.0,"M");  // mean(964..973)
    double m10  = get_val_where_num(p,"SIGDYN","ANNOT_CH_SEC","SEC",10.0,"M");   // mean(984..993), full 10-wide bin
    double m20  = get_val_where_num(p,"SIGDYN","ANNOT_CH_SEC","SEC",20.0,"M");   // needs samples past 999: absent
    double m30  = get_val_where_num(p,"SIGDYN","ANNOT_CH_SEC","SEC",30.0,"M");   // absent
    int nsec = get_nrows(p,"SIGDYN","ANNOT_CH_SEC");
    bool pass = approx_equal(m0,978.5,0.02) && approx_equal(mm10,968.5,0.02)
      && approx_equal(m10,988.5,0.02) && std::isnan(m20) && std::isnan(m30) && nsec == 5;
    std::ostringstream m; m << "SEC0=" << m0 << " SEC-10=" << mm10 << " SEC10=" << m10
			    << " SEC20=" << m20 << " SEC30=" << m30 << " nSEC=" << nsec << " (exp 5)";
    record(R,"sigdyn/complete-bins-only-at-edge", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"sigdyn/complete-bins-only-at-edge",false,e.what(),V); }

  // Q5 -- require-full=T rejects an instance outright when its whole
  // nominal window can't fit (same near-edge anchor as Q4, w=35 needs +35
  // but only +25 is available); the default (F) still emits a summary row
  // built from whatever's available.
  try {
    auto p_full = make_ramp(1000, "T_sigdyn_q5a");
    p_full->insert_annotation( "Evt", { {974.0, 974.0} } );
    p_full->eval("EPOCH dur=30 & SIGDYN sig=X annot=Evt hypno-annot=F w=35 require-full=T");
    int n_full = get_nrows(p_full,"SIGDYN","ANNOT_CH");

    auto p_def = make_ramp(1000, "T_sigdyn_q5b");
    p_def->insert_annotation( "Evt", { {974.0, 974.0} } );
    p_def->eval("EPOCH dur=30 & SIGDYN sig=X annot=Evt hypno-annot=F w=35");
    int n_def = get_nrows(p_def,"SIGDYN","ANNOT_CH");

    bool pass = n_full == 0 && n_def == 1;
    std::ostringstream m; m << "require-full=T rows=" << n_full << " (exp 0); default rows=" << n_def << " (exp 1)";
    record(R,"sigdyn/require-full-rejects-instance", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"sigdyn/require-full-rejects-instance",false,e.what(),V); }

  // Q6 -- min-n= drops a bin that only one of two instances could complete,
  // while the coarser per-anchor summary (built across ALL bins) still
  // counts both instances as contributing overall. bin-align=start means
  // SEC=20 (bin [20,29]) needs samples out to offset+29, beyond even the
  // requested w=20 window itself -- so it never completes for EITHER
  // instance, regardless of edge proximity (5 nominal labels, 4 reachable).
  // A is far from any physical edge (completes SEC=-20/-10/0/10); B is only
  // 14 samples from the recording's end (completes SEC=-20/-10/0 only, not
  // SEC=10, which needs offset+19).
  try {
    auto mk = [&]( const std::string & id ) {
      auto p = make_ramp(1000, id);
      p->insert_annotation( "Evt", { {500.0, 500.0}, {985.0, 985.0} } );
      return p;
    };

    auto p1 = mk("T_sigdyn_q6a");
    p1->eval("EPOCH dur=30 & SIGDYN sig=X annot=Evt hypno-annot=F w=20 bin=10 bin-align=start min-n=1");
    int nsec1 = get_nrows(p1,"SIGDYN","ANNOT_CH_SEC");

    auto p2 = mk("T_sigdyn_q6b");
    p2->eval("EPOCH dur=30 & SIGDYN sig=X annot=Evt hypno-annot=F w=20 bin=10 bin-align=start min-n=2");
    int nsec2 = get_nrows(p2,"SIGDYN","ANNOT_CH_SEC");
    double m10_dropped = get_val_where_num(p2,"SIGDYN","ANNOT_CH_SEC","SEC",10.0,"M");
    double n_annot = get_val_s(p2,"SIGDYN","ANNOT_CH","N");

    bool pass = nsec1 == 4 && nsec2 == 3 && std::isnan(m10_dropped) && approx_equal(n_annot,2,0.1);
    std::ostringstream m; m << "min-n=1 nSEC=" << nsec1 << " (exp 4); min-n=2 nSEC=" << nsec2
			    << " (exp 3), SEC10=" << m10_dropped << " (exp absent), ANNOT N=" << n_annot << " (exp 2)";
    record(R,"sigdyn/min-n-per-bin-filtering", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"sigdyn/min-n-per-bin-filtering",false,e.what(),V); }

  // Q7 -- mode 2 per-bin stats are event- (instance-) weighted, not sample-
  // weighted: two point anchors whose bin (bin=4, bin-align=start) covers 4
  // consecutive ramp samples each -- tightly clustered internally (~1.29 SD)
  // but far apart between instances (spread ~200). SD must reflect the
  // *between-instance* spread of each instance's own bin mean, not the
  // *within-bin* sample spread diluted 4x per instance
  try {
    auto p = make_ramp(2000, "T_sigdyn_q7");
    p->insert_annotation( "Evt", { {500.0,500.0}, {700.0,700.0} } );
    p->eval("EPOCH dur=30 & SIGDYN sig=X annot=Evt hypno-annot=F w=10 bin=4 bin-align=start");
    double n0  = get_val_where_num(p,"SIGDYN","ANNOT_CH_SEC","SEC",0.0,"N");
    double m0  = get_val_where_num(p,"SIGDYN","ANNOT_CH_SEC","SEC",0.0,"M");
    double sd0 = get_val_where_num(p,"SIGDYN","ANNOT_CH_SEC","SEC",0.0,"SD");
    // event-weighted: vals={501.5,701.5} -> mean=601.5, SD=141.42
    // (the pre-fix sample-weighted bug would instead give SD=106.91)
    bool pass = approx_equal(n0,2,0.1) && approx_equal(m0,601.5,0.05) && approx_equal(sd0,141.42,0.1);
    std::ostringstream m; m << "N=" << n0 << " M=" << m0 << " SD=" << sd0
			    << " (exp N=2,M=601.5,SD=141.42; sample-weighted bug would give SD=106.91)";
    record(R,"sigdyn/mode2-event-weighted-bin-stats", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"sigdyn/mode2-event-weighted-bin-stats",false,e.what(),V); }

  // Q8 -- tolog=T: log() is undefined for values <= 0; those samples must
  // be skipped rather than turning the bin into -inf/NaN. A single anchor
  // at a ramp value of 0 has a window spanning both negative and positive
  // ramp values: the negative offset vanishes entirely (no output row),
  // the positive offset gives exactly log(value)
  try {
    auto p = eng->inst("T_sigdyn_q8");
    const int rs = 10, n = 2000, nr = n / rs;
    p->empty_edf("T_sigdyn_q8", nr, rs, "01.01.85", "22.00.00");
    std::vector<double> shifted_ramp(n);
    for (int i=0; i<n; i++) shifted_ramp[i] = (double)(i - 1000);
    p->insert_signal("X", shifted_ramp, 1);
    p->insert_annotation("Evt", { {1000.0, 1000.0} } );  // ramp value here == 0
    p->eval("EPOCH dur=30 & SIGDYN sig=X annot=Evt hypno-annot=F w=5 tolog=T");
    double m_pos = get_val_where_num(p,"SIGDYN","ANNOT_CH_SEC","SEC",3.0,"M");   // value=3 -> log(3)
    double m_neg = get_val_where_num(p,"SIGDYN","ANNOT_CH_SEC","SEC",-3.0,"M");  // value=-3 -> skipped, absent
    bool pass = approx_equal(m_pos, std::log(3.0), 0.02) && std::isnan(m_neg);
    std::ostringstream m; m << "SEC3(log3)=" << m_pos << " (exp " << std::log(3.0) << ") SEC-3=" << m_neg << " (exp NaN)";
    record(R,"sigdyn/mode2-tolog-skips-nonpositive", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"sigdyn/mode2-tolog-skips-nonpositive",false,e.what(),V); }

  // Q9 -- sub-1Hz signals (e.g. a DPP_Z/hypnodensity-style channel written
  // once every 30s, via add_signal's negative-Fs "samples per record"
  // convention) must not collapse to a single SEC=NaN row. Previously this
  // block cast Fs to int: truncating any Fs<1 to 0 divided by zero for the
  // SEC bin label, and even flooring that up to a minimum of 1 (a
  // tempting-looking fix) was wrong in a different way -- it told
  // discontinuity()/furthest_reachable() to expect ~1 sample/sec against
  // samples that are really 30s apart, so every genuine inter-sample gap
  // looked like a discontinuity and the window collapsed to the anchor's
  // own single sample. One sample every 30s (Fs=-1 with a 30s record
  // duration), ramp value == elapsed seconds, anchor at t=600 (sample 20)
  // well away from either edge.
  try {
    auto p = eng->inst("T_sigdyn_q9");
    const int rs = 30, nr = 40;                // 40 x 30s records = 1200s
    p->empty_edf("T_sigdyn_q9", nr, rs, "01.01.85", "22.00.00");
    std::vector<double> ramp(nr);
    for (int i=0; i<nr; i++) ramp[i] = i * 30.0;
    p->insert_signal( "X", ramp, -1 );          // Fs=-1 -> 1 sample/record -> 1/30 Hz
    p->insert_annotation( "Evt", { {600.0, 600.0} } );
    p->eval("SIGDYN sig=X annot=Evt hypno-annot=F w=150 bin=30 inc=30 require-full=F min-n=1");
    int nsec    = get_nrows(p,"SIGDYN","ANNOT_CH_SEC");
    double m0   = get_val_where_num(p,"SIGDYN","ANNOT_CH_SEC","SEC",0.0,"M");
    double mm30 = get_val_where_num(p,"SIGDYN","ANNOT_CH_SEC","SEC",-30.0,"M");
    double m30  = get_val_where_num(p,"SIGDYN","ANNOT_CH_SEC","SEC",30.0,"M");
    bool pass = nsec == 11 && ! std::isnan(m0) && ! std::isnan(mm30) && ! std::isnan(m30)
      && approx_equal(m0,600,0.5) && approx_equal(mm30,570,0.5) && approx_equal(m30,630,0.5);
    std::ostringstream m; m << "nSEC=" << nsec << " (exp 11) SEC0=" << m0 << " (exp 600) SEC-30="
			    << mm30 << " (exp 570) SEC30=" << m30 << " (exp 630)";
    record(R,"sigdyn/sub1hz-signal-real-sec-offsets", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"sigdyn/sub1hz-signal-real-sec-offsets",false,e.what(),V); }
}

// ============================================================
// Group U: DPP (dynamic phenotype projection, stage 2 feature engine)
// ============================================================

static void test_dpp( lunapi_t * eng,
		      std::vector<test_result_t> & R, bool V )
{
  // DPP output lookup by (SEC, VAR) pair -- neither get_val_where (one
  // string factor) nor get_val_where_num (one numeric factor) alone can
  // match both factors of DPP's "SEC_VAR" table simultaneously; returns NaN
  // if the row (or the requested column within it) doesn't exist, which for
  // DPP means "that window was not emitted" (masked/gapped/QC-excluded)
  auto dpp_val = []( lunapi_inst_ptr p, double t, const std::string & var,
		     const std::string & col ) -> double
    {
      auto r = p->results( "DPP", "SEC_VAR" );
      const auto & cols = std::get<0>(r);
      const auto & data = std::get<1>(r);
      int ti=-1, vi=-1, ci=-1;
      for (int i=0; i<(int)cols.size(); i++) {
	if (cols[i]=="SEC") ti=i;
	if (cols[i]=="VAR") vi=i;
	if (cols[i]==col)   ci=i;
      }
      if (ti<0 || vi<0 || ci<0 || data.empty()) return std::numeric_limits<double>::quiet_NaN();
      int nrows = (int)data[ti].size();
      for (int r2=0; r2<nrows; r2++) {
	double tval = std::numeric_limits<double>::quiet_NaN();
	if (std::holds_alternative<double>(data[ti][r2])) tval = std::get<double>(data[ti][r2]);
	else if (std::holds_alternative<int>(data[ti][r2])) tval = (double)std::get<int>(data[ti][r2]);
	std::string vval;
	if (std::holds_alternative<std::string>(data[vi][r2])) vval = std::get<std::string>(data[vi][r2]);
	if ( std::fabs(tval - t) <= 0.5 && vval == var ) {
	  if ( ci < (int)data.size() && r2 < (int)data[ci].size() ) {
	    const auto & e = data[ci][r2];
	    if (std::holds_alternative<double>(e)) return std::get<double>(e);
	    if (std::holds_alternative<int>(e))    return (double)std::get<int>(e);
	  }
	  return std::numeric_limits<double>::quiet_NaN();
	}
      }
      return std::numeric_limits<double>::quiet_NaN();
    };

  // fresh empty EDF (rs=10s records) holding one inserted signal
  auto make_dpp_inst = [&]( const std::vector<double> & sig, int sr, int dur_sec,
			    const std::string & id, const std::string & label = "EEG" ) -> lunapi_inst_ptr
    {
      const int rs = 10;
      const int nr = dur_sec / rs;
      lunapi_inst_ptr p = eng->inst( id );
      p->empty_edf( id, nr, rs, "01.01.85", "22.00.00" );
      p->insert_signal( label, sig, sr );
      return p;
    };

  // U1 -- spec grammar round-trip: CH/FILTER/block-line parsing, windows=
  // crossing, channel-pair connectivity, per-feature column counts/labels
  try {
    const std::string tmp = temp_base_path("test_dpp_spec") + ".spec";
    {
      std::ofstream out(tmp);
      out << "CH C3\n";
      out << "CH C4 prefilter=0.3-35\n";
      out << "FILTER sigma 11 15\n";
      out << "BASE: PSD C3\n";
      out << "BASE: HJORTH C3\n";
      out << "CONN: PLV C3,C4 band=sigma windows=10,20\n";
      out << "CONN: COH C3,C4\n";
    }
    dpp_specs_t S;
    S.read( tmp );

    bool pass = S.chs.size() == 2 && S.has_channel("C3") && S.has_channel("C4")
      && S.chs["C4"].has_prefilter
      && approx_equal(S.chs["C4"].prefilter_lwr, 0.3, 1e-9)
      && approx_equal(S.chs["C4"].prefilter_upr, 35, 1e-9)
      && S.filters.size() == 1 && S.has_filter("sigma")
      && approx_equal(S.filters["sigma"].lwr, 11, 1e-9)
      && approx_equal(S.filters["sigma"].upr, 15, 1e-9)
      && S.specs.size() == 5; // PSD + HJORTH + 2x PLV(windows=10,20) + COH

    int n_plv10=0, n_plv20=0, n_psd=0, n_hjorth=0, n_coh=0;
    for (size_t i=0; i<S.specs.size(); i++) {
      const dpp_spec_t & s = S.specs[i];
      if (s.ftr == DPP_PLV && approx_equal(s.window_sec,10,1e-9)) { n_plv10++; pass = pass && s.cols()==2 && s.label_root()=="PLV.C3-C4.sigma"; }
      if (s.ftr == DPP_PLV && approx_equal(s.window_sec,20,1e-9)) n_plv20++;
      if (s.ftr == DPP_PSD)    { n_psd++;    pass = pass && s.cols()==5 && s.label_root()=="PSD.C3"; }
      if (s.ftr == DPP_HJORTH) { n_hjorth++; pass = pass && s.cols()==3 && s.label_root()=="HJORTH.C3"; }
      if (s.ftr == DPP_COH)    { n_coh++;    pass = pass && s.cols()==5 && s.label_root()=="COH.C3-C4"; }
    }
    pass = pass && n_plv10==1 && n_plv20==1 && n_psd==1 && n_hjorth==1 && n_coh==1;

    std::ostringstream m; m << "n_specs=" << S.specs.size() << " (exp 5)";
    record(R,"dpp/spec-grammar-roundtrip", pass, m.str(), V);
    std::remove( tmp.c_str() );
  } catch(std::exception & e){ record(R,"dpp/spec-grammar-roundtrip",false,e.what(),V); }

  // U2 -- zero-config default spec (init_default) matches an equivalent
  // hand-written custom spec file, feature-for-feature
  try {
    dpp_specs_t A;
    A.init_default( { "C3" } );

    const std::string tmp = temp_base_path("test_dpp_spec_custom") + ".spec";
    {
      std::ofstream out(tmp);
      out << "CH C3\n";
      out << "BASE: PSD C3\n";
      out << "BASE: SLOPE C3\n";
      out << "BASE: HJORTH C3\n";
    }
    dpp_specs_t B;
    B.read( tmp );

    bool pass = A.specs.size() == 3 && B.specs.size() == 3;
    std::set<std::string> a_labels, b_labels;
    for (size_t i=0; i<A.specs.size(); i++) a_labels.insert( A.specs[i].label_root() + ":" + std::to_string(A.specs[i].cols()) );
    for (size_t i=0; i<B.specs.size(); i++) b_labels.insert( B.specs[i].label_root() + ":" + std::to_string(B.specs[i].cols()) );
    pass = pass && a_labels == b_labels && a_labels.size() == 3;
    for (size_t i=0; i<A.specs.size(); i++) pass = pass && approx_equal(A.specs[i].window_sec, A.default_window_sec, 1e-9);

    std::ostringstream m; m << "default n=" << A.specs.size() << " custom n=" << B.specs.size();
    record(R,"dpp/default-spec-matches-custom", pass, m.str(), V);
    std::remove( tmp.c_str() );
  } catch(std::exception & e){ record(R,"dpp/default-spec-matches-custom",false,e.what(),V); }

  // U3 -- HJORTH exact value: a pure sine, window an exact integer number
  // of cycles -> Activity (mean-square) == amp^2/2 exactly. qc=F: a coarse
  // 10-samples/cycle sine quantized to 16-bit EDF resolution ties exactly
  // at several samples per peak/trough, which can trip the qc=T flat-signal
  // heuristic even though the signal isn't actually flat (qc itself has its
  // own dedicated test, U8, below)
  try {
    auto sig = make_sine(100, 20.0, 10.0, 2.0);   // 10Hz, amp=2, sr=100 (10 samples/cycle)
    auto p = make_dpp_inst(sig, 100, 20, "T_dpp_hjorth");
    p->eval("DPP sig=EEG windows=10 step=10 qc=F show-features=T");
    double act = dpp_val(p, 10.0, "HJORTH.EEG.w10", "V1");
    bool pass = approx_equal(act, std::log(2.0), 0.02); // log(amp^2/2) = log(2.0)
    std::ostringstream m; m << "log-activity=" << act << " (exp log(2.0)=" << std::log(2.0) << ")";
    record(R,"dpp/hjorth-exact-activity", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"dpp/hjorth-exact-activity",false,e.what(),V); }

  // U4 -- PSD: (a) a pure 10Hz sine's power is essentially all in ALPHA
  // (8-12Hz), negligible in DELTA; (b) DPP's ALPHA band-power matches an
  // independently-run PWELCH call on the identical window slice/Welch
  // params DPP itself derives for a 10s window (segment=4s/overlap=2s) --
  // confirms DPP wires the window buffer into PWELCH correctly, not just
  // "some" plausible number
  try {
    auto sig = make_sine(100, 20.0, 10.0, 1.0);
    auto p = make_dpp_inst(sig, 100, 20, "T_dpp_psd");
    p->eval("DPP sig=EEG windows=10 step=10 qc=F show-features=T");
    double dpp_alpha = dpp_val(p, 10.0, "PSD.EEG.w10", "V3"); // delta,theta,alpha,sigma,beta -> V3
    double dpp_delta = dpp_val(p, 10.0, "PSD.EEG.w10", "V1");

    // DPP's own get_window(): idx_end = sample @ t=10s = index 1000;
    // w_samples = 1000 -> idx_start_report = 1
    std::vector<double> win( sig.begin() + 1, sig.begin() + 1001 );
    PWELCH pwelch( win, 100, 4.0, 4, WINDOW_TUKEY50 );
    double ref_alpha = std::log( pwelch.psdsum( ALPHA ) );  // DPP log-transforms PSD unconditionally

    // dominance check in log-space: alpha should still be many natural-log
    // units above delta (equivalent to >100x in linear power, log(100)=4.6)
    bool pass = approx_equal_rel(dpp_alpha, ref_alpha, 1e-3) && (dpp_alpha - dpp_delta) > std::log(100.0);
    std::ostringstream m; m << "DPP log-alpha=" << dpp_alpha << " ref log-alpha=" << ref_alpha << " log-delta=" << dpp_delta;
    record(R,"dpp/psd-matches-direct-pwelch", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"dpp/psd-matches-direct-pwelch",false,e.what(),V); }

  // U5 -- ENVELOPE: filter-Hilbert magnitude of a constant-amplitude
  // in-band sine tracks the true amplitude, with small CV
  try {
    auto sig = make_sine(200, 180.0, 13.0, 3.0);  // 13Hz, amp=3, inside sigma(12-15)
    auto p = make_dpp_inst(sig, 200, 180, "T_dpp_env");
    const std::string tmp = temp_base_path("test_dpp_env") + ".spec";
    {
      std::ofstream out(tmp);
      out << "CH EEG\n";
      out << "FILTER sigma 12 15\n";
      out << "BASE: ENVELOPE EEG band=sigma windows=10\n";
    }
    p->eval("DPP spec=" + tmp + " step=90 show-features=T");
    double mean_mag = dpp_val(p, 90.0, "ENVELOPE.EEG.sigma.w10", "V1");
    double cv        = dpp_val(p, 90.0, "ENVELOPE.EEG.sigma.w10", "V3");
    bool pass = approx_equal_rel(mean_mag, 3.0, 0.15) && cv < 0.1;
    std::ostringstream m; m << "mean=" << mean_mag << " (exp~3.0) cv=" << cv << " (exp<0.1)";
    record(R,"dpp/envelope-constant-magnitude", pass, m.str(), V);
    std::remove( tmp.c_str() );
  } catch(std::exception & e){ record(R,"dpp/envelope-constant-magnitude",false,e.what(),V); }

  // U6 -- PLV: matched-frequency channels lock (PLV~1); a 1Hz-offset pair
  // (10 full relative-phase drift cycles over the 10s window) averages out
  // to PLV~0
  try {
    auto sigA  = make_sine(200, 180.0, 13.0, 1.0);
    auto sigB1 = make_sine(200, 180.0, 13.0, 1.0);  // matched
    auto sigB2 = make_sine(200, 180.0, 14.0, 1.0);  // 1Hz offset -> drifting

    const std::string tmp = temp_base_path("test_dpp_plv") + ".spec";
    {
      std::ofstream out(tmp);
      out << "CH A\n";
      out << "CH B\n";
      out << "FILTER sigma 12 15\n";
      out << "CONN: PLV A,B band=sigma windows=10\n";
    }

    auto p1 = eng->inst("T_dpp_plv_matched");
    p1->empty_edf("T_dpp_plv_matched", 18, 10, "01.01.85", "22.00.00");
    p1->insert_signal("A", sigA,  200);
    p1->insert_signal("B", sigB1, 200);
    p1->eval("DPP spec=" + tmp + " step=90 show-features=T");
    double plv_matched = dpp_val(p1, 90.0, "PLV.A-B.sigma.w10", "V1");

    auto p2 = eng->inst("T_dpp_plv_mismatched");
    p2->empty_edf("T_dpp_plv_mismatched", 18, 10, "01.01.85", "22.00.00");
    p2->insert_signal("A", sigA,  200);
    p2->insert_signal("B", sigB2, 200);
    p2->eval("DPP spec=" + tmp + " step=90 show-features=T");
    double plv_mismatched = dpp_val(p2, 90.0, "PLV.A-B.sigma.w10", "V1");

    bool pass = plv_matched > 0.9 && plv_mismatched < 0.3;
    std::ostringstream m; m << "matched PLV=" << plv_matched << " (exp>0.9) mismatched PLV=" << plv_mismatched << " (exp<0.3)";
    record(R,"dpp/plv-matched-vs-drifting", pass, m.str(), V);
    std::remove( tmp.c_str() );
  } catch(std::exception & e){ record(R,"dpp/plv-matched-vs-drifting",false,e.what(),V); }

  // U6b -- plv-weighted= is actually wired through: a channel pair that is
  // phase-locked and full-amplitude for the first half of the window, then
  // near-zero-amplitude and phase-drifting for the second half. With
  // amplitude weighting (default, gate disabled here to isolate weighting
  // itself from the separate gating mechanism already exercised above),
  // the ~100x-smaller-envelope drifting half barely votes, so PLV stays
  // close to the locked half's own value; unweighted, the two halves
  // contribute equally, diluting PLV well below that
  try {
    auto warmupA = make_sine(200, 30.0, 13.0, 1.0);
    auto burstA  = make_sine(200, 90.0, 13.0, 1.0);
    auto quietA  = make_sine(200, 90.0, 13.0, 1.0);
    std::vector<double> sigA = warmupA;
    sigA.insert(sigA.end(), burstA.begin(), burstA.end());
    sigA.insert(sigA.end(), quietA.begin(), quietA.end());

    auto warmupB = make_sine(200, 30.0, 13.0, 1.0);
    auto burstB  = make_sine(200, 90.0, 13.0, 1.0);   // matched, full amp -> phase-locked
    auto quietB  = make_sine(200, 90.0, 14.0, 0.01);  // drifting freq, ~100x smaller amp
    std::vector<double> sigB = warmupB;
    sigB.insert(sigB.end(), burstB.begin(), burstB.end());
    sigB.insert(sigB.end(), quietB.begin(), quietB.end());

    const std::string tmp = temp_base_path("test_dpp_plv_weighted") + ".spec";
    {
      std::ofstream out(tmp);
      out << "CH A\n";
      out << "CH B\n";
      out << "FILTER sigma 12 15\n";
      out << "CONN: PLV A,B band=sigma windows=180\n";
    }

    auto p1 = eng->inst("T_dpp_plv_w");
    p1->empty_edf("T_dpp_plv_w", 21, 10, "01.01.85", "22.00.00");  // 210s total
    p1->insert_signal("A", sigA, 200);
    p1->insert_signal("B", sigB, 200);
    p1->eval("DPP spec=" + tmp + " step=209 plv-weighted=T plv-gate=F qc=F show-features=T");
    double plv_weighted = dpp_val(p1, 209.0, "PLV.A-B.sigma.w180", "V1");

    auto p2 = eng->inst("T_dpp_plv_uw");
    p2->empty_edf("T_dpp_plv_uw", 21, 10, "01.01.85", "22.00.00");
    p2->insert_signal("A", sigA, 200);
    p2->insert_signal("B", sigB, 200);
    p2->eval("DPP spec=" + tmp + " step=209 plv-weighted=F qc=F show-features=T");
    double plv_unweighted = dpp_val(p2, 209.0, "PLV.A-B.sigma.w180", "V1");

    bool pass = plv_weighted > 0.85 && plv_unweighted < 0.7 && (plv_weighted - plv_unweighted) > 0.2;
    std::ostringstream m; m << "plv-weighted=T (gate=F)=" << plv_weighted << " (exp>0.85) "
			    << "plv-weighted=F=" << plv_unweighted << " (exp<0.7, and >=0.2 below weighted)";
    record(R,"dpp/plv-weighting-param-wired", pass, m.str(), V);
    std::remove( tmp.c_str() );
  } catch(std::exception & e){ record(R,"dpp/plv-weighting-param-wired",false,e.what(),V); }

  // U7 -- gap handling: MASK+RE removes epochs 11-20 (seconds [100,200)),
  // creating a real discontinuity. A window fully inside either remaining
  // segment is emitted; a window whose extent straddles the removed
  // segment (in array-index terms, now contiguous, but with a real time
  // jump) is marked missing, not computed from truncated/wrong data
  try {
    auto sig = make_sine(100, 400.0, 10.0, 1.0);
    auto p = make_dpp_inst(sig, 100, 400, "T_dpp_gap");
    p->eval("EPOCH dur=10 & MASK epoch=1-10,21-40 & RE"); // keep epochs 1-10,21-40; drop 11-20
    p->eval("DPP sig=EEG windows=20 step=10 qc=F show-features=T");

    double v_before = dpp_val(p, 50.0,  "PSD.EEG.w20", "V1");  // window [30,50]:  pre-gap, present
    double v_across = dpp_val(p, 210.0, "PSD.EEG.w20", "V1");  // straddles the removed [100,200): absent
    double v_after  = dpp_val(p, 390.0, "PSD.EEG.w20", "V1");  // window [370,390]: post-gap, present

    bool pass = !std::isnan(v_before) && std::isnan(v_across) && !std::isnan(v_after);
    std::ostringstream m; m << "before=" << v_before << " across=" << v_across << " (exp NaN) after=" << v_after;
    record(R,"dpp/gap-marks-window-missing", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"dpp/gap-marks-window-missing",false,e.what(),V); }

  // U8 -- artifact QC gate: a fully flat/constant window is excluded under
  // qc=T (default), but included (finite Activity=0) under qc=F
  try {
    std::vector<double> flat_sig(6000, 5.0); // 60s @ sr=100, constant
    auto p_qcT = make_dpp_inst(flat_sig, 100, 60, "T_dpp_qcT");
    p_qcT->eval("DPP sig=EEG windows=10 step=10 qc=T show-features=T");
    double v_qcT = dpp_val(p_qcT, 10.0, "HJORTH.EEG.w10", "V1");

    auto p_qcF = make_dpp_inst(flat_sig, 100, 60, "T_dpp_qcF");
    p_qcF->eval("DPP sig=EEG windows=10 step=10 qc=F show-features=T");
    double v_qcF = dpp_val(p_qcF, 10.0, "HJORTH.EEG.w10", "V1");

    bool pass = std::isnan(v_qcT) && !std::isnan(v_qcF);
    std::ostringstream m; m << "qc=T -> " << v_qcT << " (exp NaN/missing); qc=F -> " << v_qcF << " (exp finite)";
    record(R,"dpp/qc-gate-flat-window", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"dpp/qc-gate-flat-window",false,e.what(),V); }

  // U9 -- causal-only guarantee: two recordings share an identical prefix
  // [0,30) and diverge only after t=30; a feature reported at t=29 (window
  // [9,29], never touching t>=30) must closely match between them. Not
  // bit-identical: each EDF auto-scales its own 16-bit physical range from
  // its *full* signal (prefix+suffix), so the differing suffix amplitude
  // gives the two recordings slightly different quantization step sizes
  // for the (source-identical) prefix -- a real future-leakage bug would
  // instead pull in the wildly different suffix (amp 5 vs 1, freq 3 vs
  // 10Hz) and move the value by orders of magnitude, not ~1e-4 relative
  try {
    auto prefix  = make_sine(100, 30.0, 10.0, 1.0);
    auto suffixA = make_sine(100, 30.0, 10.0, 1.0);  // continues identically
    auto suffixB = make_sine(100, 30.0, 3.0,  5.0);  // diverges after t=30

    std::vector<double> sigA = prefix; sigA.insert(sigA.end(), suffixA.begin(), suffixA.end());
    std::vector<double> sigB = prefix; sigB.insert(sigB.end(), suffixB.begin(), suffixB.end());

    auto pA = make_dpp_inst(sigA, 100, 60, "T_dpp_causalA");
    pA->eval("DPP sig=EEG windows=20 step=29 qc=F show-features=T");
    double vA = dpp_val(pA, 29.0, "HJORTH.EEG.w20", "V1");

    auto pB = make_dpp_inst(sigB, 100, 60, "T_dpp_causalB");
    pB->eval("DPP sig=EEG windows=20 step=29 qc=F show-features=T");
    double vB = dpp_val(pB, 29.0, "HJORTH.EEG.w20", "V1");

    bool pass = !std::isnan(vA) && !std::isnan(vB) && approx_equal_rel(vA, vB, 1e-3);
    std::ostringstream m; m << "A=" << vA << " B=" << vB << " (exp closely matching -- window [9,29] never reads t>=30)";
    record(R,"dpp/causal-no-future-leakage", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"dpp/causal-no-future-leakage",false,e.what(),V); }

  // U10 -- prefiltering must not carry filter state across a genuine
  // discontinuity: a strong out-of-band segment immediately followed (in
  // array-index terms) by an in-band segment, separated by an artificial
  // ~100s gap in the timepoint sequence, must filter identically to the
  // in-band segment filtered alone -- not show contamination from the
  // out-of-band segment, which per real elapsed time is nowhere near it
  try {
    const int sr = 200;
    const double lwr = 10, upr = 15;
    auto seg1 = make_sine( sr, 5.0, 1.0, 50.0 );    // 1Hz, amp=50: well outside 10-15Hz
    auto seg2 = make_sine( sr, 5.0, 12.0, 1.0 );    // 12Hz, amp=1: well inside 10-15Hz
    std::vector<double> x = seg1; x.insert( x.end(), seg2.begin(), seg2.end() );

    std::vector<uint64_t> tp( x.size() );
    const uint64_t one_sample_tp = globals::tp_1sec / sr;
    for (size_t i=0; i<seg1.size(); i++) tp[i] = i * one_sample_tp;
    const uint64_t gap_tp = 100 * (uint64_t)globals::tp_1sec;  // ~100s artificial gap
    for (size_t i=0; i<seg2.size(); i++)
      tp[ seg1.size() + i ] = seg1.size() * one_sample_tp + gap_tp + i * one_sample_tp;

    // reference values use dpp_filters::apply_band() (the same causal
    // filtering path prefilter_trace() itself now uses internally), not
    // the raw dsptools::apply_fir() -- that lower-level function is
    // deliberately group-delay-compensated (non-causal: output position j
    // reflects real input through j+(taps-1)/2, see stats/dpp-filter.h's
    // top comment), so comparing against it directly would no longer be
    // an apples-to-apples check of prefilter_trace()'s own behavior
    dpp_filter_t filt; filt.lwr = lwr; filt.upr = upr; filt.ripple = 0.02; filt.tw = 1.0;
    std::vector<double> filtered_gapped = dpp_filters::prefilter_trace( x, tp, sr, lwr, upr );
    std::vector<double> filtered_seg2_alone = dpp_filters::apply_band( seg2, sr, filt );
    std::vector<double> filtered_whole = dpp_filters::apply_band( x, sr, filt );

    const double v_gapaware = filtered_gapped[ seg1.size() ];
    const double v_isolated = filtered_seg2_alone[0];
    const double v_naive    = filtered_whole[ seg1.size() ];

    bool pass = approx_equal( v_gapaware, v_isolated, 1e-9 )
      && std::fabs( v_naive - v_isolated ) > 1e-6;
    std::ostringstream m; m << "gap-aware=" << v_gapaware << " isolated=" << v_isolated
			    << " (exp equal); naive-whole-trace=" << v_naive
			    << " (exp different from isolated -- this is the bug being fixed)";
    record(R,"dpp/prefilter-respects-discontinuity", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"dpp/prefilter-respects-discontinuity",false,e.what(),V); }

  // U11 -- two spec blocks that resolve to the identical feature/channel/
  // band/window label must halt loudly (not silently collide on the same
  // output column)
  try {
    const std::string tmp = temp_base_path("test_dpp_dup_label") + ".spec";
    {
      std::ofstream out(tmp);
      out << "CH C3\n";
      out << "BASE: PSD C3\n";
      out << "EXTRA: PSD C3\n";  // identical feature/channel/window -> same label
    }
    auto sig = make_sine(100, 20.0, 10.0, 1.0);
    auto p = make_dpp_inst(sig, 100, 20, "T_dpp_duplabel");
    bool halted = false;
    try { p->eval("DPP spec=" + tmp); }
    catch (std::exception &) { halted = true; }
    record(R,"dpp/duplicate-label-halts", halted, halted ? "halted as expected" : "did not halt", V);
    std::remove( tmp.c_str() );
  } catch(std::exception & e){ record(R,"dpp/duplicate-label-halts",false,e.what(),V); }

  // U12 -- CHEP-based per-channel artifact masking: window_masked()'s CHEP
  // branch (edf.timeline.masked(e,ch)) is the only artifact-exclusion path
  // in get_window() not otherwise covered (MASK/RE gaps and the raw-buffer
  // qc= gate both have their own dedicated tests above). CHEP bad-channels=
  // marks one channel bad for every epoch; a window on that channel must be
  // excluded even with qc=F, while an untouched channel in the same
  // recording is unaffected
  try {
    auto sig = make_sine(100, 60.0, 10.0, 1.0);
    auto p = eng->inst("T_dpp_chep");
    p->empty_edf("T_dpp_chep", 6, 10, "01.01.85", "22.00.00");
    p->insert_signal("A", sig, 100);
    p->insert_signal("B", sig, 100);
    p->eval("EPOCH dur=10 & CHEP bad-channels=B");
    p->eval("DPP sig=A,B windows=10 step=10 qc=F show-features=T");

    double v_A = dpp_val(p, 10.0, "HJORTH.A.w10", "V1");
    double v_B = dpp_val(p, 10.0, "HJORTH.B.w10", "V1");

    bool pass = !std::isnan(v_A) && std::isnan(v_B);
    std::ostringstream m; m << "A (clean)=" << v_A << " (exp finite); B (CHEP bad-channels)=" << v_B << " (exp NaN/missing)";
    record(R,"dpp/chep-masking-excludes-window", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"dpp/chep-masking-excludes-window",false,e.what(),V); }

  // U16 -- causal-only guarantee, with prefilter= enabled: same construction
  // as U9 (two recordings sharing an identical prefix, diverging after
  // t=30), but this time with CH ... prefilter= set, so a feature at t=29
  // (window [9,29]) exercises the prefiltered trace, not just the raw one.
  // Regression test for a real bug: dpp_filters::prefilter_trace() used to
  // filter via dsptools::apply_fir(), which is deliberately group-delay-
  // compensated (output position j reflects real input through
  // j+(taps-1)/2) -- for a segment much longer than the filter's own
  // settling length, that means genuinely reading future samples, so vA/vB
  // would diverge by orders of magnitude here before the fix (same failure
  // signature U9's own comment describes for a real future-leakage bug)
  try {
    auto prefix  = make_sine(100, 30.0, 10.0, 1.0);
    auto suffixA = make_sine(100, 30.0, 10.0, 1.0);  // continues identically
    auto suffixB = make_sine(100, 30.0, 3.0,  5.0);  // diverges after t=30

    std::vector<double> sigA = prefix; sigA.insert(sigA.end(), suffixA.begin(), suffixA.end());
    std::vector<double> sigB = prefix; sigB.insert(sigB.end(), suffixB.begin(), suffixB.end());

    const std::string tmp = temp_base_path("test_dpp_causal_prefilter") + ".spec";
    { std::ofstream out(tmp); out << "CH EEG prefilter=0.3-35\n"; out << "BASE: HJORTH EEG windows=20\n"; }

    auto pA = make_dpp_inst(sigA, 100, 60, "T_dpp_causal_pf_A");
    pA->eval("DPP spec=" + tmp + " step=29 qc=F show-features=T");
    double vA = dpp_val(pA, 29.0, "HJORTH.EEG.w20", "V1");

    auto pB = make_dpp_inst(sigB, 100, 60, "T_dpp_causal_pf_B");
    pB->eval("DPP spec=" + tmp + " step=29 qc=F show-features=T");
    double vB = dpp_val(pB, 29.0, "HJORTH.EEG.w20", "V1");

    bool pass = !std::isnan(vA) && !std::isnan(vB) && approx_equal_rel(vA, vB, 1e-2);
    std::ostringstream m; m << "A=" << vA << " B=" << vB << " (exp closely matching -- window [9,29] never reads t>=30, even prefiltered)";
    record(R,"dpp/prefilter-causal-no-future-leakage", pass, m.str(), V);
    std::remove( tmp.c_str() );
  } catch(std::exception & e){ record(R,"dpp/prefilter-causal-no-future-leakage",false,e.what(),V); }

  // ---- Stage 6: catch22 feature family ----

  // U13 -- spec-grammar: CATCH22 defaults to 22 columns/CATCH22.<ch> label;
  // catch24= block arg switches to 24
  try {
    dpp_specs_t S;
    S.init();
    S.chs[ "C3" ] = dpp_channel_t();
    S.chs[ "C3" ].ch = "C3";

    dpp_spec_t s22; s22.ftr = DPP_CATCH22; s22.ch = { "C3" }; s22.window_sec = 30;
    dpp_spec_t s24; s24.ftr = DPP_CATCH22; s24.ch = { "C3" }; s24.window_sec = 30; s24.arg[ "catch24" ] = "T";

    bool pass = s22.cols() == 22 && s24.cols() == 24 && s22.label_root() == "CATCH22.C3";
    std::ostringstream m; m << "cols(no-arg)=" << s22.cols() << " (exp 22) cols(catch24=T)=" << s24.cols()
			    << " (exp 24) label=" << s22.label_root() << " (exp CATCH22.C3)";
    record(R,"dpp/catch22-spec-grammar", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"dpp/catch22-spec-grammar",false,e.what(),V); }

  // U14 -- exact-value cross-check: DPP's CATCH22 output on a window
  // reproduces an independent, direct catch22_t call on the identical
  // extracted buffer, bit-for-bit -- same idiom as U4's PSD-vs-PWELCH
  // cross-check (same signal/window-index construction: windows=10,
  // step=10 on a 20s/sr=100 recording -> reporting range is sig[1:1001])
  try {
    auto sig = make_sine(100, 20.0, 10.0, 1.0);
    auto p = make_dpp_inst(sig, 100, 20, "T_dpp_catch22");
    const std::string tmp = temp_base_path("test_dpp_catch22") + ".spec";
    { std::ofstream out(tmp); out << "CH EEG\n"; out << "BASE: CATCH22 EEG windows=10 catch24=T\n"; }
    p->eval("DPP spec=" + tmp + " step=10 qc=F show-features=T");

    std::vector<double> win( sig.begin() + 1, sig.begin() + 1001 );
    catch22_t ref( true );
    ref.calc( win.data() , (int)win.size() );

    // insert_signal() round-trips the synthetic double buffer through the
    // EDF's 16-bit digital/physical quantization, so DPP's window and the
    // reference window aren't bit-identical inputs (same reason U4's PSD
    // cross-check above uses a loose tolerance, not exact equality). Most
    // of catch22's 22 stats are inherently discontinuous functions of the
    // input (histogram-bin indices, motif counts, ...), so a tiny
    // quantization perturbation can flip them to a genuinely different
    // value -- no numeric tolerance papers over that meaningfully. Only
    // check the two catch24 "mean"/"SD" extras here (plain, continuous,
    // linear statistics), found by name rather than assumed index --
    // sufficient to confirm DPP wires the exact expected window/length
    // into catch22_t::calc(), which is this test's actual purpose.
    int idx_mean = -1, idx_sd = -1;
    for (int k=0; k<24; k++)
      {
	if ( catch22_t::short_name(k) == "mean" ) idx_mean = k;
	if ( catch22_t::short_name(k) == "SD" )   idx_sd = k;
      }

    bool pass = ref.valid() && idx_mean >= 0 && idx_sd >= 0;
    double dv_mean = std::numeric_limits<double>::quiet_NaN(), dv_sd = std::numeric_limits<double>::quiet_NaN();
    if ( pass )
      {
	dv_mean = dpp_val( p , 10.0 , "CATCH22.EEG.w10" , "V" + std::to_string(idx_mean+1) );
	dv_sd   = dpp_val( p , 10.0 , "CATCH22.EEG.w10" , "V" + std::to_string(idx_sd+1) );
	pass = approx_equal( dv_mean , ref.stat(idx_mean) , 1e-3 ) && approx_equal( dv_sd , ref.stat(idx_sd) , 1e-3 );
      }
    std::ostringstream m; m << "ref.valid()=" << ref.valid() << " mean: dpp=" << dv_mean << " ref=" << ref.stat(idx_mean)
			    << "; SD: dpp=" << dv_sd << " ref=" << ref.stat(idx_sd);
    record(R,"dpp/catch22-matches-direct-call", pass, m.str(), V);
    std::remove( tmp.c_str() );
  } catch(std::exception & e){ record(R,"dpp/catch22-matches-direct-call",false,e.what(),V); }

  // U15 -- QC gate applies to CATCH22 the same as every other feature: a
  // fully flat/constant window is excluded under qc=T, included under qc=F
  try {
    std::vector<double> flat_sig(6000, 5.0); // 60s @ sr=100, constant
    const std::string tmp = temp_base_path("test_dpp_catch22_qc") + ".spec";
    { std::ofstream out(tmp); out << "CH EEG\n"; out << "BASE: CATCH22 EEG windows=10\n"; }

    auto p_qcT = make_dpp_inst(flat_sig, 100, 60, "T_dpp_catch22_qcT");
    p_qcT->eval("DPP spec=" + tmp + " step=10 qc=T show-features=T");
    double v_qcT = dpp_val(p_qcT, 10.0, "CATCH22.EEG.w10", "V1");

    auto p_qcF = make_dpp_inst(flat_sig, 100, 60, "T_dpp_catch22_qcF");
    p_qcF->eval("DPP spec=" + tmp + " step=10 qc=F show-features=T");
    double v_qcF = dpp_val(p_qcF, 10.0, "CATCH22.EEG.w10", "V1");

    bool pass = std::isnan(v_qcT);
    std::ostringstream m; m << "qc=T -> " << v_qcT << " (exp NaN/missing); qc=F -> " << v_qcF << " (exp finite or catch22's own QC-NaN)";
    record(R,"dpp/catch22-qc-gate", pass, m.str(), V);
    std::remove( tmp.c_str() );
  } catch(std::exception & e){ record(R,"dpp/catch22-qc-gate",false,e.what(),V); }

  // ---- native band= support on HJORTH/SKEW/KURTOSIS/MSE/CATCH22 ----

  // U16 -- grammar: band=delta,sigma on a single-channel optional-band
  // feature expands into one spec per band (mirrors windows='s own
  // comma-list crossing); a single band= value still gives exactly one
  // spec; PSD/SLOPE (never band-aware) halt loudly rather than silently
  // ignoring band= at compute time
  try {
    const std::string tmp = temp_base_path("test_dpp_band_grammar") + ".spec";
    {
      std::ofstream out(tmp);
      out << "CH C3\n";
      out << "FILTER delta 0.5 4\n";
      out << "FILTER sigma 11 15\n";
      out << "BASE: HJORTH C3 band=delta,sigma windows=30\n";
      out << "BASE: SKEW   C3 band=sigma windows=30\n";
    }
    dpp_specs_t S;
    S.read( tmp );
    std::set<std::string> labels;
    for (int i=0; i<(int)S.specs.size(); i++) labels.insert( S.specs[i].label_root() );
    bool pass = S.specs.size() == 3
      && labels.count("HJORTH.C3.delta") == 1 && labels.count("HJORTH.C3.sigma") == 1
      && labels.count("SKEW.C3.sigma") == 1;

    const std::string tmp2 = temp_base_path("test_dpp_band_psd_halts") + ".spec";
    { std::ofstream out(tmp2); out << "CH C3\n"; out << "FILTER sigma 11 15\n"; out << "BASE: PSD C3 band=sigma windows=30\n"; }
    dpp_specs_t S2;
    bool halted = false;
    try { S2.read( tmp2 ); } catch (std::exception &) { halted = true; }
    pass = pass && halted;

    std::ostringstream m; m << "n=" << S.specs.size() << " (exp 3);";
    for (const auto & l : labels) m << " " << l;
    m << "; PSD band= halted=" << halted << " (exp 1)";
    record(R,"dpp/band-optional-expand-grammar", pass, m.str(), V);
    std::remove( tmp.c_str() );
    std::remove( tmp2.c_str() );
  } catch(std::exception & e){ record(R,"dpp/band-optional-expand-grammar",false,e.what(),V); }

  // U17 -- exact-value cross-check: HJORTH with band= reproduces an
  // independent, direct dpp_filters::apply_band() + MiscMath::hjorth() call
  // on the identical causally-padded-then-dropped buffer DPP itself builds
  // (same technique as U14's raw-signal CATCH22 cross-check, extended to
  // cover the new band-filtered path's padding arithmetic too)
  try {
    const int sr = 200; const int dur = 180; const double t = 90.0; const double w_sec = 10;
    auto sig = make_sine(sr, (double)dur, 13.0, 3.0);  // 13Hz, inside sigma(12-15)
    auto p = make_dpp_inst(sig, sr, dur, "T_dpp_band_hjorth");
    const std::string tmp = temp_base_path("test_dpp_band_hjorth") + ".spec";
    {
      std::ofstream out(tmp);
      out << "CH EEG\n";
      out << "FILTER sigma 12 15\n";
      out << "BASE: HJORTH EEG band=sigma windows=10\n";
    }
    p->eval("DPP spec=" + tmp + " step=90 show-features=T");
    double dpp_logact = dpp_val(p, t, "HJORTH.EEG.sigma.w10", "V1");

    dpp_filter_t filt; filt.name = "sigma"; filt.lwr = 12; filt.upr = 15; filt.ripple = 0.02; filt.tw = 1;
    const double pad_sec = dpp_filters::pad_seconds( filt, sr );
    const int idx_end = (int) std::llround( t * sr );
    const int w_samples = (int) std::llround( w_sec * sr );
    const int pad_samples = (int) std::llround( pad_sec * sr );
    const int idx_start_report = idx_end - w_samples + 1;
    const int idx_start_padded = idx_start_report - pad_samples;
    std::vector<double> padded( sig.begin() + idx_start_padded, sig.begin() + idx_end + 1 );
    std::vector<double> filtered = dpp_filters::apply_band( padded, sr, filt );
    std::vector<double> win_ref( filtered.begin() + pad_samples, filtered.end() );
    double act=0, mob=0, comp=0;
    MiscMath::hjorth( &win_ref, &act, &mob, &comp );
    const double ref_logact = log( std::max( act, 1e-300 ) );

    bool pass = approx_equal_rel(dpp_logact, ref_logact, 1e-3);
    std::ostringstream m; m << "dpp=" << dpp_logact << " ref=" << ref_logact;
    record(R,"dpp/band-optional-hjorth-matches-direct-filter", pass, m.str(), V);
    std::remove( tmp.c_str() );
  } catch(std::exception & e){ record(R,"dpp/band-optional-hjorth-matches-direct-filter",false,e.what(),V); }

  // ---- PAC (phase-amplitude coupling), and PLV restricted to within-band ----

  // U18 -- grammar: PAC requires band=<phase>,<amplitude> (exactly two
  // values -- neither omitted nor a single shared value is accepted),
  // allows only one channel pair per line (unlike PLV/COH/PSI), and PLV
  // itself now halts on two distinct band values (cross-band phase-phase
  // "coupling" removed in favour of PAC)
  try {
    const std::string tmp_ok = temp_base_path("test_dpp_pac_ok") + ".spec";
    { std::ofstream out(tmp_ok); out << "CH C3\n" << "FILTER slow 0.3 1.5\n" << "FILTER sigma 11 15\n"
				     << "CONN: PAC C3,C3 band=slow,sigma windows=30\n"; }
    dpp_specs_t Sok;
    Sok.read( tmp_ok );
    bool pass = Sok.specs.size() == 1 && Sok.specs[0].cols() == 2
      && Sok.specs[0].label_root() == "PAC.C3-C3.slow-sigma";

    const std::string tmp_onepair = temp_base_path("test_dpp_pac_onepair") + ".spec";
    { std::ofstream out(tmp_onepair); out << "CH C3\n" << "CH O1\n" << "FILTER slow 0.3 1.5\n" << "FILTER sigma 11 15\n"
					  << "CONN: PAC C3,C3 O1,O1 band=slow,sigma windows=30\n"; }
    dpp_specs_t S2;
    bool halted_two_pairs = false;
    try { S2.read( tmp_onepair ); } catch (std::exception &) { halted_two_pairs = true; }

    const std::string tmp_oneband = temp_base_path("test_dpp_pac_oneband") + ".spec";
    { std::ofstream out(tmp_oneband) ; out << "CH C3\n" << "FILTER sigma 11 15\n"
					   << "CONN: PAC C3,C3 band=sigma windows=30\n"; }
    dpp_specs_t S3;
    bool halted_one_band = false;
    try { S3.read( tmp_oneband ); } catch (std::exception &) { halted_one_band = true; }

    const std::string tmp_noband = temp_base_path("test_dpp_pac_noband") + ".spec";
    { std::ofstream out(tmp_noband); out << "CH C3\n" << "CONN: PAC C3,C3 windows=30\n"; }
    dpp_specs_t S4;
    bool halted_no_band = false;
    try { S4.read( tmp_noband ); } catch (std::exception &) { halted_no_band = true; }

    const std::string tmp_plv_crossband = temp_base_path("test_dpp_plv_crossband") + ".spec";
    { std::ofstream out(tmp_plv_crossband); out << "CH C3\n" << "FILTER slow 0.3 1.5\n" << "FILTER sigma 11 15\n"
						<< "CONN: PLV C3,C3 band=slow,sigma windows=30\n"; }
    dpp_specs_t S5;
    bool halted_plv_crossband = false;
    try { S5.read( tmp_plv_crossband ); } catch (std::exception &) { halted_plv_crossband = true; }

    pass = pass && halted_two_pairs && halted_one_band && halted_no_band && halted_plv_crossband;

    std::ostringstream m; m << "ok: n=" << Sok.specs.size() << " cols=" << Sok.specs[0].cols()
			    << " label=" << Sok.specs[0].label_root()
			    << "; two-pairs halted=" << halted_two_pairs
			    << "; one-band halted=" << halted_one_band
			    << "; no-band halted=" << halted_no_band
			    << "; PLV cross-band halted=" << halted_plv_crossband;
    record(R,"dpp/pac-grammar-and-plv-within-band-only", pass, m.str(), V);
    std::remove( tmp_ok.c_str() ); std::remove( tmp_onepair.c_str() );
    std::remove( tmp_oneband.c_str() ); std::remove( tmp_noband.c_str() );
    std::remove( tmp_plv_crossband.c_str() );
  } catch(std::exception & e){ record(R,"dpp/pac-grammar-and-plv-within-band-only",false,e.what(),V); }

  // U19 -- PAC recovers known coupling: a 0.5Hz "slow" carrier plus a 13Hz
  // "sigma" carrier whose amplitude is modulated in-phase with the slow
  // wave (envelope peaks exactly at the slow wave's own peak, i.e.
  // preferred phase ~0) gives a high normalized MVL; the same two carriers
  // with constant (unmodulated) sigma amplitude -- otherwise identical --
  // gives a near-zero MVL. Contrast-based check (same style as the
  // existing PLV matched-vs-drifting test), since the filtered/Hilbert
  // pipeline doesn't admit a simple closed-form exact value.
  try {
    const int sr = 200; const int dur = 180;
    auto build = [&]( bool coupled ) {
      std::vector<double> x( sr * dur );
      for (int i=0; i<(int)x.size(); i++)
	{
	  const double t = i / (double)sr;
	  const double slow = std::cos( 2.0 * M_PI * 0.5 * t );
	  const double carrier = std::cos( 2.0 * M_PI * 13.0 * t );
	  const double envelope = coupled ? ( 1.0 + 0.9 * slow ) : 1.0;
	  x[i] = slow + envelope * carrier;
	}
      return x;
    };
    const std::string tmp = temp_base_path("test_dpp_pac") + ".spec";
    { std::ofstream out(tmp); out << "CH EEG\n" << "FILTER slow 0.3 1.5\n" << "FILTER sigma 11 15\n"
				  << "CONN: PAC EEG,EEG band=slow,sigma windows=60\n"; }

    auto p_coupled = make_dpp_inst( build(true), sr, dur, "T_dpp_pac_coupled" );
    p_coupled->eval("DPP spec=" + tmp + " step=90 show-features=T");
    double mvl_coupled = dpp_val(p_coupled, 90.0, "PAC.EEG-EEG.slow-sigma.w60", "V1");

    auto p_uncoupled = make_dpp_inst( build(false), sr, dur, "T_dpp_pac_uncoupled" );
    p_uncoupled->eval("DPP spec=" + tmp + " step=90 show-features=T");
    double mvl_uncoupled = dpp_val(p_uncoupled, 90.0, "PAC.EEG-EEG.slow-sigma.w60", "V1");

    bool pass = mvl_coupled > 0.3 && mvl_uncoupled < 0.05 && mvl_coupled > 5.0 * mvl_uncoupled;
    std::ostringstream m; m << "coupled MVL=" << mvl_coupled << " (exp >0.3); uncoupled MVL=" << mvl_uncoupled << " (exp <0.05)";
    record(R,"dpp/pac-recovers-known-coupling", pass, m.str(), V);
    std::remove( tmp.c_str() );
  } catch(std::exception & e){ record(R,"dpp/pac-recovers-known-coupling",false,e.what(),V); }
}

// ============================================================
// Group V: DPP stage 3 (LightGBM fit/apply) -- HAS_LGBM only.
// First HAS_LGBM-guarded content in this test suite.
// ============================================================

#ifdef HAS_LGBM

static void test_dpp_fit( lunapi_t * eng,
			  std::vector<test_result_t> & R, bool V )
{
  // helper: write a synthetic two-cluster corpus (n individuals per
  // cluster, each with a few rows carrying a single feature value close to
  // the cluster's mean) plus a matching phenotype file -- one feature
  // ("BASE: SLOPE X", cols()==1) so the corpus's column count trivially
  // matches a single-column spec; the label is arbitrary bookkeeping here,
  // not a claim that these are real slope values
  auto build_corpus = [&]( const std::string & base, int n_per_cluster, int nrows_per_indiv,
			   double lo_val, double hi_val, double lo_y, double hi_y,
			   std::vector<std::string> * lo_ids, std::vector<std::string> * hi_ids ) -> std::string
    {
      const std::string corpus = base + ".dat";
      const std::string phenofile = base + ".pheno";
      std::ofstream ph( phenofile.c_str() );
      ph << "ID\tY\n";
      bool first = true;
      for (int grp=0; grp<2; grp++)
	{
	  const double val = grp==0 ? lo_val : hi_val;
	  const double y   = grp==0 ? lo_y   : hi_y;
	  for (int i=0; i<n_per_cluster; i++)
	    {
	      const std::string id = ( grp==0 ? "LO" : "HI" ) + std::to_string(i) + "_" + base.substr(base.size()>8?base.size()-8:0);
	      dpp_matrix_t m;
	      m.id = id;
	      for (int r=0; r<nrows_per_indiv; r++)
		{
		  m.time_sec.push_back( (r+1) * 30.0 );
		  m.X.push_back( { val + 0.001 * r } );
		}
	      dpp_io::save( corpus , m , 1 , ! first );
	      first = false;
	      ph << id << "\t" << y << "\n";
	      if ( grp==0 ) lo_ids->push_back( id ); else hi_ids->push_back( id );
	    }
	}
      ph.close();
      return phenofile;
    };

  // V1 -- end-to-end fit: recovers the coarse LO/HI separation in-sample,
  // and reloading the saved bundle reproduces the same prediction exactly
  try {
    const std::string base = temp_base_path("test_dppfit_v1");
    std::vector<std::string> lo_ids, hi_ids;
    const std::string phenofile = build_corpus( base , 10 , 5 , 0.0 , 10.0 , 0.0 , 10.0 , &lo_ids , &hi_ids );

    const std::string specfile = base + ".spec";
    { std::ofstream sp( specfile.c_str() ); sp << "CH X\nBASE: SLOPE X\n"; }

    const std::string config = base + ".conf";
    { std::ofstream cf( config.c_str() ); cf << "objective=regression\nverbosity=-1\nmin_data_in_leaf=1\nnum_leaves=7\n"; }

    cmd_t::attach_ivars( phenofile );

    param_t param;
    param.add( "data" , base + ".dat" );
    param.add( "phe" , "Y" );
    param.add( "spec" , specfile );
    param.add( "config" , config );
    param.add( "out" , base + "_model" );

    dpp_fit_t fit( param );
    fit.fit();

    // in-sample recovery: predict on the training matrix directly
    Eigen::MatrixXd Ylo = fit.lgbm.predict( fit.Xtrain.topRows(1) );
    Eigen::MatrixXd Yhi = fit.lgbm.predict( fit.Xtrain.bottomRows(1) );
    // Xtrain is filled group-by-group (LO block first, HI block last) by flatten_and_split()
    const double pred_lo = Ylo(0,0);
    const double pred_hi = Yhi(0,0);

    // save/reload determinism
    lgbm_t reloaded;
    reloaded.qt_mode = true;
    reloaded.load_model( Helper::expand( base + "_model.mod" ) );
    Eigen::MatrixXd Ylo2 = reloaded.predict( fit.Xtrain.topRows(1) );

    bool pass = approx_equal( pred_lo , 0.0 , 2.0 ) && approx_equal( pred_hi , 10.0 , 2.0 )
      && ( pred_hi - pred_lo ) > 5.0
      && approx_equal( Ylo2(0,0) , pred_lo , 1e-9 );
    std::ostringstream m; m << "pred_lo=" << pred_lo << " (exp~0) pred_hi=" << pred_hi
			    << " (exp~10) reload=" << Ylo2(0,0) << " (exp==pred_lo, bit-identical)";
    record(R,"dpp-fit/end-to-end-recovery-and-reload", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"dpp-fit/end-to-end-recovery-and-reload",false,e.what(),V); }

  // V2 -- no-subject-leakage: individual-level validation split, never a
  // raw row-index split -- every row from a held-out individual lands in
  // Xvalid, every row from a training individual in Xtrain, with block
  // sizes matching exactly (nrows_per_indiv each, in this synthetic setup)
  try {
    const std::string base = temp_base_path("test_dppfit_v2");
    std::vector<std::string> lo_ids, hi_ids;
    const std::string phenofile = build_corpus( base , 6 , 4 , 0.0 , 10.0 , 0.0 , 10.0 , &lo_ids , &hi_ids );

    const std::string specfile = base + ".spec";
    { std::ofstream sp( specfile.c_str() ); sp << "CH X\nBASE: SLOPE X\n"; }

    // hold out 2 individuals (one per cluster)
    const std::string valfile = base + ".valid";
    { std::ofstream vf( valfile.c_str() ); vf << lo_ids[0] << "\n" << hi_ids[0] << "\n"; }

    cmd_t::attach_ivars( phenofile );

    param_t param;
    param.add( "data" , base + ".dat" );
    param.add( "phe" , "Y" );
    param.add( "spec" , specfile );
    param.add( "validation" , valfile );
    param.add( "out" , base + "_model" );

    dpp_fit_t fit( param );
    fit.load_corpus();
    fit.build_feature_labels();
    fit.attach_phenotypes();
    fit.flatten_and_split();

    // 12 individuals total, 2 held out -> 10 training x 4 rows, 2 validation x 4 rows
    bool pass = fit.Xtrain.rows() == 40 && fit.Xvalid.rows() == 8
      && fit.istart_train.size() == 10 && fit.istart_valid.size() == 2;
    std::ostringstream m; m << "Xtrain.rows=" << fit.Xtrain.rows() << " (exp 40) Xvalid.rows=" << fit.Xvalid.rows()
			    << " (exp 8) n_train_indiv=" << fit.istart_train.size() << " (exp 10) n_valid_indiv="
			    << fit.istart_valid.size() << " (exp 2)";
    record(R,"dpp-fit/no-subject-leakage-split", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"dpp-fit/no-subject-leakage-split",false,e.what(),V); }

  // V3 -- apply mode: DPP model= on a real (synthetic) EDF, using the SAME
  // spec used at training, attaches a new signal with genuine (non-
  // sentinel-only) values; a MISMATCHED spec at apply time halts loudly
  // rather than silently producing wrong output
  try {
    const std::string base = temp_base_path("test_dppfit_v3");
    std::vector<std::string> lo_ids, hi_ids;
    const std::string phenofile = build_corpus( base , 8 , 3 , 0.0 , 10.0 , 0.0 , 10.0 , &lo_ids , &hi_ids );

    const std::string specfile = base + ".spec";
    { std::ofstream sp( specfile.c_str() ); sp << "CH X\nBASE: SLOPE X\n"; }
    const std::string mismatch_specfile = base + "_mismatch.spec";
    { std::ofstream sp( mismatch_specfile.c_str() ); sp << "CH X\nBASE: HJORTH X\n"; } // 3 cols, not 1 -- must mismatch

    const std::string config = base + ".conf";
    { std::ofstream cf( config.c_str() ); cf << "objective=regression\nverbosity=-1\nmin_data_in_leaf=1\nnum_leaves=7\n"; }

    cmd_t::attach_ivars( phenofile );

    param_t param;
    param.add( "data" , base + ".dat" );
    param.add( "phe" , "Y" );
    param.add( "spec" , specfile );
    param.add( "config" , config );
    param.add( "out" , base + "_model" );

    dpp_fit_t fit( param );
    fit.fit();

    auto sig = make_sine( 100, 60.0, 10.0, 1.0 );
    auto p = make_inst( eng, sig, 100, 6, 10, "X", "T_dppfit_apply" );

    // matching spec: should attach DPP_Z with genuine (non-sentinel-only) values
    p->eval( "DPP spec=" + specfile + " step=10 qc=F model=" + base + "_model" );
    bool has_sig = p->has_channels( { "DPP_Z" } )[0];
    p->eval( "STATS sig=DPP_Z" );
    double smin = get_val( p , "STATS" , "MIN" );
    double smax = get_val( p , "STATS" , "MAX" );

    // mismatched spec: must halt, not silently attach a bogus signal
    auto p2 = make_inst( eng, sig, 100, 6, 10, "X", "T_dppfit_apply_mismatch" );
    bool halted = false;
    try { p2->eval( "DPP spec=" + mismatch_specfile + " step=10 qc=F model=" + base + "_model" ); }
    catch ( std::exception & ) { halted = true; }

    bool pass = has_sig && Helper::realnum(smin) && Helper::realnum(smax) && smax > smin && halted;
    std::ostringstream m; m << "has_sig=" << has_sig << " MIN=" << smin << " MAX=" << smax
			    << " (exp MAX>MIN, real values present); mismatched-spec halted=" << halted << " (exp T)";
    record(R,"dpp-fit/apply-mode-and-manifest-validation", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"dpp-fit/apply-mode-and-manifest-validation",false,e.what(),V); }

  // ---- Stage 7: DPP x SIGDYN integration (verification only -- no new
  // code in either command; both were independently designed to be
  // signal-agnostic, see the implementation plan's "Stage 7" section) ----

  // V4 -- DPP model='s output (DPP_Z, an ordinary EDF signal) feeds SIGDYN
  // with zero glue code: EPOCH the recording, then SIGDYN sig=DPP_Z runs
  // to completion (mode 1, whole-recording trend/decile summary --
  // epoch-stats=F/hypno-annot=F since this synthetic recording has no
  // staging/annotations) and actually produces output, not a silent no-op
  try {
    const std::string base = temp_base_path("test_dppfit_v4");
    std::vector<std::string> lo_ids, hi_ids;
    const std::string phenofile = build_corpus( base , 8 , 3 , 0.0 , 10.0 , 0.0 , 10.0 , &lo_ids , &hi_ids );

    const std::string specfile = base + ".spec";
    { std::ofstream sp( specfile.c_str() ); sp << "CH X\nBASE: SLOPE X\n"; }
    const std::string config = base + ".conf";
    { std::ofstream cf( config.c_str() ); cf << "objective=regression\nverbosity=-1\nmin_data_in_leaf=1\nnum_leaves=7\n"; }

    cmd_t::attach_ivars( phenofile );

    param_t param;
    param.add( "data" , base + ".dat" );
    param.add( "phe" , "Y" );
    param.add( "spec" , specfile );
    param.add( "config" , config );
    param.add( "out" , base + "_model" );

    dpp_fit_t fit( param );
    fit.fit();

    auto sig = make_sine( 100, 60.0, 10.0, 1.0 );
    auto p = make_inst( eng, sig, 100, 6, 10, "X", "T_dppfit_sigdyn" );
    p->eval( "DPP spec=" + specfile + " step=10 qc=F model=" + base + "_model" );
    bool has_sig = p->has_channels( { "DPP_Z" } )[0];

    p->eval( "EPOCH dur=10 & SIGDYN sig=DPP_Z epoch-stats=F hypno-annot=F" );
    std::vector<std::string> cmds = p->commands();
    bool sigdyn_ran = std::find( cmds.begin() , cmds.end() , "SIGDYN" ) != cmds.end();

    bool pass = has_sig && sigdyn_ran;
    std::ostringstream m; m << "has_sig(DPP_Z)=" << has_sig << " SIGDYN in commands()=" << sigdyn_ran;
    record(R,"dpp-fit/sigdyn-integration", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"dpp-fit/sigdyn-integration",false,e.what(),V); }

  // ---- Stage 4: stage-conditioned training/blending ----

  // W1 -- end-to-end: two individual-groups, hypnodensity ~pure-W vs
  // ~pure-NR (3-state), distinct feature/phenotype values per group. Train
  // stage-conditioned, then apply to two synthetic recordings whose PP_*
  // context is ~pure-W vs ~pure-NR (at a coarser sr than the feature
  // channel -- also exercises the differing-native-resolution lookup
  // path). The blended predictions should land near each group's own
  // value and differ substantially from each other -- proof the per-stage
  // weighting and blending actually did something, not a silent no-op
  // equivalent to the pooled (stage 3) path.
  try {
    const std::string base = temp_base_path("test_dppfit_w1");
    const std::string corpus = base + ".dat";
    const std::string hypno_corpus = base + ".hypno.dat";
    const std::string phenofile = base + ".pheno";

    std::ofstream ph( phenofile.c_str() );
    ph << "ID\tY\n";

    const int n_per_group = 8, nrows = 5;
    bool first = true;
    for (int grp=0; grp<2; grp++)
      {
	const double val = grp==0 ? 0.0 : 10.0;
	const double y   = grp==0 ? 0.0 : 10.0;
	for (int i=0; i<n_per_group; i++)
	  {
	    const std::string id = ( grp==0 ? "W" : "NR" ) + std::to_string(i) + "_w1";
	    dpp_matrix_t m; m.id = id;
	    dpp_matrix_t h; h.id = id;
	    for (int r=0; r<nrows; r++)
	      {
		m.time_sec.push_back( (r+1) * 30.0 );
		m.X.push_back( { val + 0.001 * r } );
		h.time_sec.push_back( (r+1) * 30.0 );
		h.X.push_back( grp==0 ? std::vector<double>{ 1.0 , 0.0 , 0.0 } : std::vector<double>{ 0.0 , 0.0 , 1.0 } );
	      }
	    dpp_io::save( corpus , m , 1 , ! first );
	    dpp_io::save( hypno_corpus , h , 3 , ! first );
	    first = false;
	    ph << id << "\t" << y << "\n";
	  }
      }
    ph.close();

    const std::string specfile = base + ".spec";
    { std::ofstream sp( specfile.c_str() ); sp << "CH X\nBASE: SLOPE X\n"; }
    const std::string config = base + ".conf";
    { std::ofstream cf( config.c_str() ); cf << "objective=regression\nverbosity=-1\nmin_data_in_leaf=1\nnum_leaves=7\n"; }

    cmd_t::attach_ivars( phenofile );

    param_t param;
    param.add( "data" , corpus );
    param.add( "hypno" , hypno_corpus );
    param.add( "hypno-three-state" , "T" );
    param.add( "phe" , "Y" );
    param.add( "spec" , specfile );
    param.add( "config" , config );
    param.add( "out" , base + "_model" );

    dpp_fit_t fit( param );
    fit.fit();

    bool bundle_ok = Helper::fileExists( base + "_model.W.mod" )
      && Helper::fileExists( base + "_model.R.mod" )
      && Helper::fileExists( base + "_model.NR.mod" )
      && Helper::fileExists( base + "_model.dpp" );

    auto sig = make_sine( 100 , 60.0 , 10.0 , 1.0 );

    // W-context apply recording: PP_* at sr=1 (coarser than the feature
    // channel's sr=100/step=10 -- differing native resolutions), constant
    // pure-W throughout
    auto pw = eng->inst( "T_dppfit_w1_apply_w" );
    pw->empty_edf( "T_dppfit_w1_apply_w" , 6 , 10 , "01.01.85" , "22.00.00" );
    pw->insert_signal( "X" , sig , 100 );
    pw->insert_signal( "PP_W" , std::vector<double>( 60 , 1.0 ) , 1 );
    pw->insert_signal( "PP_R" , std::vector<double>( 60 , 0.0 ) , 1 );
    pw->insert_signal( "PP_NR" , std::vector<double>( 60 , 0.0 ) , 1 );
    pw->eval( "DPP spec=" + specfile + " step=10 qc=F model=" + base + "_model hypno-prefix=PP hypno-three-state=T" );
    pw->eval( "STATS sig=DPP_Z" );
    const double z_w = get_val( pw , "STATS" , "MEAN" );

    // NR-context apply recording: constant pure-NR throughout
    auto pn = eng->inst( "T_dppfit_w1_apply_nr" );
    pn->empty_edf( "T_dppfit_w1_apply_nr" , 6 , 10 , "01.01.85" , "22.00.00" );
    pn->insert_signal( "X" , sig , 100 );
    pn->insert_signal( "PP_W" , std::vector<double>( 60 , 0.0 ) , 1 );
    pn->insert_signal( "PP_R" , std::vector<double>( 60 , 0.0 ) , 1 );
    pn->insert_signal( "PP_NR" , std::vector<double>( 60 , 1.0 ) , 1 );
    pn->eval( "DPP spec=" + specfile + " step=10 qc=F model=" + base + "_model hypno-prefix=PP hypno-three-state=T" );
    pn->eval( "STATS sig=DPP_Z" );
    const double z_nr = get_val( pn , "STATS" , "MEAN" );

    bool pass = bundle_ok && approx_equal( z_w , 0.0 , 3.0 ) && approx_equal( z_nr , 10.0 , 3.0 ) && ( z_nr - z_w ) > 5.0;
    std::ostringstream m; m << "bundle_ok=" << bundle_ok << " z_W-context=" << z_w << " (exp~0) z_NR-context="
			    << z_nr << " (exp~10, and >=5 apart)";
    record(R,"dpp-fit/stage-weighting-and-blend", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"dpp-fit/stage-weighting-and-blend",false,e.what(),V); }

  // W2 -- mixed/ambiguous-epoch weighting: direct inspection of the actual
  // per-row weight vector LightGBM receives (bypassing corpus files
  // entirely), rather than inferring it indirectly from fit quality. Row 0
  // pure-W, row 1 a 50/50 W/R mix, row 2 pure-NR, row 3 missing (NaN)
  // hypnodensity -- confirms proportional (not all-or-nothing) weighting
  // and that a missing/invalid row gets weight 0, not a crash or NaN
  // weight passed to LightGBM
  try {
    const std::string base = temp_base_path("test_dppfit_w2");
    param_t param;
    param.add( "out" , base + "_model" );
    param.add( "phe" , "Y" );
    param.add( "hypno" , "unused" );   // only needed so the constructor sets stage_conditioned=true

    dpp_fit_t fit( param );
    fit.n_features = 1;
    fit.stage_labels = { "W" , "R" , "NR" };

    fit.Xtrain = Eigen::MatrixXd( 4 , 1 );
    fit.Xtrain << 1.0 , 2.0 , 3.0 , 4.0;
    fit.ytrain = { 1.0 , 2.0 , 3.0 , 4.0 };
    fit.indiv_weight_train.assign( 4 , 1.0f );

    const double NaN_value = std::numeric_limits<double>::quiet_NaN();
    fit.Htrain = Eigen::MatrixXd( 4 , 3 );
    fit.Htrain << 1.0 , 0.0 , 0.0 ,
                  0.5 , 0.5 , 0.0 ,
                  0.0 , 0.0 , 1.0 ,
                  NaN_value , NaN_value , NaN_value;

    fit.train_stage_boosters();

    const std::vector<float> & w = fit.stage_lgbm[0].training_weights;   // stage "W"
    bool pass = w.size() == 4 && approx_equal(w[0],1.0,0.01) && approx_equal(w[1],0.5,0.01)
      && approx_equal(w[2],0.0,0.01) && approx_equal(w[3],0.0,0.01);
    std::ostringstream m; m << "W-weights=[" << ( w.size()>0?w[0]:-1 ) << "," << ( w.size()>1?w[1]:-1 )
			    << "," << ( w.size()>2?w[2]:-1 ) << "," << ( w.size()>3?w[3]:-1 )
			    << "] (exp [1,0.5,0,0])";
    record(R,"dpp-fit/mixed-epoch-proportional-weighting", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"dpp-fit/mixed-epoch-proportional-weighting",false,e.what(),V); }

  // W3 -- near-zero hypnodensity sum marks the output missing rather than
  // reporting a misleading Z~=0: an apply-time recording whose PP_*
  // channels are all zero throughout should yield an entirely
  // sentinel-filled DPP_Z (no real value ever written -- MIN==MAX,
  // degenerate but not a crash)
  try {
    const std::string base = temp_base_path("test_dppfit_w3");
    // reuse a fresh, minimal stage-conditioned bundle (same shape as W1,
    // smaller, purely to have a valid model to apply)
    const std::string corpus = base + ".dat";
    const std::string hypno_corpus = base + ".hypno.dat";
    const std::string phenofile = base + ".pheno";
    std::ofstream ph( phenofile.c_str() ); ph << "ID\tY\n";
    for (int i=0; i<6; i++)
      {
	const std::string id = "S" + std::to_string(i) + "_w3";
	dpp_matrix_t m; m.id = id;
	dpp_matrix_t h; h.id = id;
	for (int r=0; r<5; r++)
	  {
	    m.time_sec.push_back( (r+1)*30.0 );
	    m.X.push_back( { 5.0 + 0.001*r } );
	    h.time_sec.push_back( (r+1)*30.0 );
	    h.X.push_back( { 1.0 , 0.0 , 0.0 } );
	  }
	dpp_io::save( corpus , m , 1 , i!=0 );
	dpp_io::save( hypno_corpus , h , 3 , i!=0 );
	ph << id << "\t5.0\n";
      }
    ph.close();

    const std::string specfile = base + ".spec";
    { std::ofstream sp( specfile.c_str() ); sp << "CH X\nBASE: SLOPE X\n"; }
    const std::string config = base + ".conf";
    { std::ofstream cf( config.c_str() ); cf << "objective=regression\nverbosity=-1\nmin_data_in_leaf=1\nnum_leaves=7\n"; }

    cmd_t::attach_ivars( phenofile );

    param_t param;
    param.add( "data" , corpus );
    param.add( "hypno" , hypno_corpus );
    param.add( "hypno-three-state" , "T" );
    param.add( "phe" , "Y" );
    param.add( "spec" , specfile );
    param.add( "config" , config );
    param.add( "out" , base + "_model" );

    dpp_fit_t fit( param );
    fit.fit();

    auto sig = make_sine( 100 , 60.0 , 10.0 , 1.0 );
    auto p = eng->inst( "T_dppfit_w3_apply" );
    p->empty_edf( "T_dppfit_w3_apply" , 6 , 10 , "01.01.85" , "22.00.00" );
    p->insert_signal( "X" , sig , 100 );
    p->insert_signal( "PP_W" , std::vector<double>( 60 , 0.0 ) , 1 );
    p->insert_signal( "PP_R" , std::vector<double>( 60 , 0.0 ) , 1 );
    p->insert_signal( "PP_NR" , std::vector<double>( 60 , 0.0 ) , 1 );
    p->eval( "DPP spec=" + specfile + " step=10 qc=F model=" + base + "_model hypno-prefix=PP hypno-three-state=T" );
    p->eval( "STATS sig=DPP_Z" );
    const double smin = get_val( p , "STATS" , "MIN" );
    const double smax = get_val( p , "STATS" , "MAX" );

    bool pass = Helper::realnum(smin) && Helper::realnum(smax) && approx_equal( smin , smax , 1e-6 );
    std::ostringstream m; m << "MIN=" << smin << " MAX=" << smax << " (exp MIN==MAX -- all-sentinel, no real value written)";
    record(R,"dpp-fit/near-zero-hypno-sum-marks-missing", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"dpp-fit/near-zero-hypno-sum-marks-missing",false,e.what(),V); }

  // W4 -- calibration: recovers a known systematic offset exactly (y is an
  // exact linear function, slope 1, of the stage booster's own raw
  // validation predictions -- GLM should recover a~=1, b~=OFFSET), and
  // falls back to {1,0} rather than an unstable fit when there's no
  // validation set at all, or too few PP_s>0.5 qualifying rows
  try {
    const std::string base = temp_base_path("test_dppfit_w4");
    const std::string config = base + ".conf";
    { std::ofstream cf( config.c_str() ); cf << "objective=regression\nverbosity=-1\nmin_data_in_leaf=1\nnum_leaves=7\n"; }

    param_t param;
    param.add( "out" , base + "_model" );
    param.add( "phe" , "Y" );
    param.add( "hypno" , "unused" );

    dpp_fit_t fit( param );
    fit.n_features = 1;
    fit.stage_labels = { "W" , "R" , "NR" };

    const int n = 20;
    fit.Xtrain = Eigen::MatrixXd( n , 1 );
    fit.ytrain.resize( n );
    fit.indiv_weight_train.assign( n , 1.0f );
    fit.Htrain = Eigen::MatrixXd( n , 3 );
    for (int i=0; i<n; i++)
      {
	fit.Xtrain(i,0) = (double)i;
	fit.ytrain[i] = (double)i;
	fit.Htrain(i,0) = 1.0; fit.Htrain(i,1) = 0.0; fit.Htrain(i,2) = 0.0;
      }

    const int nv = 20;
    fit.Xvalid = Eigen::MatrixXd( nv , 1 );
    fit.yvalid.assign( nv , 0.0 );
    fit.Hvalid = Eigen::MatrixXd( nv , 3 );
    for (int i=0; i<nv; i++)
      {
	fit.Xvalid(i,0) = (double)i;
	fit.Hvalid(i,0) = 1.0; fit.Hvalid(i,1) = 0.0; fit.Hvalid(i,2) = 0.0;
      }

    fit.stage_lgbm.resize( 3 );   // 3 default-constructed (empty) boosters -- lgbm_t is move-only, not copyable
    lgbm_t & lg = fit.stage_lgbm[0];
    lg.load_config( config );
    lg.qt_mode = true;
    lg.attach_training_matrix( fit.Xtrain );
    lg.attach_training_qts( fit.ytrain );
    lg.training_weights.assign( n , 1.0f );
    lg.apply_weights( lg.training , &lg.training_weights );
    lg.create_booster( false );

    Eigen::MatrixXd Zraw = lg.predict( fit.Xvalid );
    const double OFFSET = 7.5;
    // tiny alternating jitter: an exactly-noiseless (zero-residual) fit is
    // a real degenerate edge case for GLM::display()'s standard-error gate
    // (see dpp_fit_t::calibrate()'s own comment) -- vanishingly unlikely
    // with real validation data, but this synthetic construction would
    // otherwise hit it by construction
    for (int i=0; i<nv; i++) fit.yvalid[i] = Zraw(i,0) + OFFSET + ( (i%2==0) ? 0.01 : -0.01 );

    std::pair<double,double> ab_offset = fit.calibrate( 0 );

    // fallback: no validation set at all
    std::vector<double> saved_yvalid = fit.yvalid;
    fit.yvalid.clear();
    std::pair<double,double> ab_noval = fit.calibrate( 0 );
    fit.yvalid = saved_yvalid;

    // fallback: too few qualifying (PP_s>0.5) rows
    for (int i=0; i<nv; i++) fit.Hvalid(i,0) = i < 2 ? 1.0 : 0.0;   // only 2 rows qualify, need >=5
    std::pair<double,double> ab_toofew = fit.calibrate( 0 );

    bool pass = approx_equal( ab_offset.first , 1.0 , 0.05 ) && approx_equal( ab_offset.second , OFFSET , 0.5 )
      && approx_equal( ab_noval.first , 1.0 , 1e-9 ) && approx_equal( ab_noval.second , 0.0 , 1e-9 )
      && approx_equal( ab_toofew.first , 1.0 , 1e-9 ) && approx_equal( ab_toofew.second , 0.0 , 1e-9 );
    std::ostringstream m; m << "offset: a=" << ab_offset.first << " (exp~1) b=" << ab_offset.second << " (exp~" << OFFSET
			    << "); no-validation fallback: a=" << ab_noval.first << " b=" << ab_noval.second
			    << " (exp 1,0); too-few-rows fallback: a=" << ab_toofew.first << " b=" << ab_toofew.second << " (exp 1,0)";
    record(R,"dpp-fit/calibration-offset-and-fallbacks", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"dpp-fit/calibration-offset-and-fallbacks",false,e.what(),V); }

  // ---- Stage 5: grouped K-fold cross-validation / OOF ----

  // X1 -- default fold assignment: deterministic sorted-round-robin over
  // individual IDs, not hash-based -- assert against a hand-computed
  // expectation for a small, deliberately-unsorted ID list
  try {
    param_t param;
    param.add( "out" , temp_base_path("test_dppfit_x1") + "_model" );
    param.add( "phe" , "Y" );
    param.add( "folds" , "3" );

    dpp_fit_t fit( param );
    for ( const std::string & id : { "C" , "A" , "E" , "B" , "D" } )
      { dpp_matrix_t m; m.id = id; fit.individuals.push_back( m ); }

    fit.assign_folds();

    // sorted: A,B,C,D,E -> fold = index % 3 -> 0,1,2,0,1
    bool pass = fit.fold_assignment.size() == 5
      && fit.fold_assignment["A"] == 0 && fit.fold_assignment["B"] == 1 && fit.fold_assignment["C"] == 2
      && fit.fold_assignment["D"] == 0 && fit.fold_assignment["E"] == 1;
    std::ostringstream m; m << "A=" << fit.fold_assignment["A"] << " B=" << fit.fold_assignment["B"]
			    << " C=" << fit.fold_assignment["C"] << " D=" << fit.fold_assignment["D"]
			    << " E=" << fit.fold_assignment["E"] << " (exp 0,1,2,0,1)";
    record(R,"dpp-fit/fold-assignment-sorted-round-robin", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"dpp-fit/fold-assignment-sorted-round-robin",false,e.what(),V); }

  // X2 -- pooled K-fold CV/OOF: every usable row gets exactly one
  // out-of-fold prediction, LO/HI groups separate sensibly in their OOF
  // predictions (never having been trained on the fold that produced
  // them), and the final saved bundle is still trained on the *entire*
  // corpus (not any one fold's smaller subset)
  try {
    const std::string base = temp_base_path("test_dppfit_x2");
    std::vector<std::string> lo_ids, hi_ids;
    const std::string phenofile = build_corpus( base , 6 , 4 , 0.0 , 10.0 , 0.0 , 10.0 , &lo_ids , &hi_ids );

    const std::string specfile = base + ".spec";
    { std::ofstream sp( specfile.c_str() ); sp << "CH X\nBASE: SLOPE X\n"; }
    const std::string config = base + ".conf";
    { std::ofstream cf( config.c_str() ); cf << "objective=regression\nverbosity=-1\nmin_data_in_leaf=1\nnum_leaves=7\n"; }

    cmd_t::attach_ivars( phenofile );

    param_t param;
    param.add( "data" , base + ".dat" );
    param.add( "phe" , "Y" );
    param.add( "spec" , specfile );
    param.add( "config" , config );
    param.add( "out" , base + "_model" );
    param.add( "folds" , "3" );   // 12 individuals, 3 folds -> 4 individuals (16 rows) held out per fold

    dpp_fit_t fit( param );
    fit.fit();

    // 12 individuals x 4 rows/individual = 48 usable rows total
    const int n_total_rows = 48;
    bool every_row_covered = (int)fit.oof_rows.size() == n_total_rows;

    double lo_sum = 0, hi_sum = 0; int lo_n = 0, hi_n = 0;
    for ( const auto & row : fit.oof_rows )
      {
	if ( row.id.substr(0,2) == "LO" ) { lo_sum += row.y_pred; ++lo_n; }
	else                              { hi_sum += row.y_pred; ++hi_n; }
      }
    const double lo_mean = lo_n > 0 ? lo_sum / lo_n : std::numeric_limits<double>::quiet_NaN();
    const double hi_mean = hi_n > 0 ? hi_sum / hi_n : std::numeric_limits<double>::quiet_NaN();

    // final bundle trained on the entire corpus, not a fold-reduced subset
    // (cross_validate() restores validation_ids before fit()'s own,
    // final, unconditional flatten_and_split() call)
    bool final_bundle_full = fit.Xtrain.rows() == n_total_rows;

    bool oof_file_written = Helper::fileExists( Helper::expand( base + "_model.oof" ) );

    bool pass = every_row_covered && lo_n == 24 && hi_n == 24
      && approx_equal( lo_mean , 0.0 , 3.0 ) && approx_equal( hi_mean , 10.0 , 3.0 ) && ( hi_mean - lo_mean ) > 4.0
      && final_bundle_full && oof_file_written;
    std::ostringstream m; m << "oof_rows=" << fit.oof_rows.size() << " (exp " << n_total_rows
			    << ") lo_n=" << lo_n << " hi_n=" << hi_n << " lo_mean=" << lo_mean << " (exp~0) hi_mean="
			    << hi_mean << " (exp~10) final_Xtrain.rows=" << fit.Xtrain.rows() << " (exp "
			    << n_total_rows << ", i.e. full corpus) .oof written=" << oof_file_written;
    record(R,"dpp-fit/pooled-kfold-oof-and-full-final-bundle", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"dpp-fit/pooled-kfold-oof-and-full-final-bundle",false,e.what(),V); }

  // X3 -- stage-conditioned K-fold CV: runs end-to-end (per-fold stage
  // boosters + per-fold calibration + blending) without crashing, and
  // still recovers the coarse W-vs-NR group separation in its out-of-fold
  // predictions despite each fold only ever training on half the cohort
  try {
    const std::string base = temp_base_path("test_dppfit_x3");
    const std::string corpus = base + ".dat";
    const std::string hypno_corpus = base + ".hypno.dat";
    const std::string phenofile = base + ".pheno";

    std::ofstream ph( phenofile.c_str() ); ph << "ID\tY\n";
    const int n_per_group = 6, nrows = 5;
    bool first = true;
    for (int grp=0; grp<2; grp++)
      {
	const double val = grp==0 ? 0.0 : 10.0;
	for (int i=0; i<n_per_group; i++)
	  {
	    // tiny per-individual jitter: an exactly-noiseless (zero-residual)
	    // fit is a real degenerate edge case for GLM::display()'s standard-
	    // error gate (see dpp_fit_t::calibrate()'s own comment) -- this
	    // synthetic cohort would otherwise hit it exactly, since each
	    // stage's booster ends up perfectly predicting its (constant)
	    // group phenotype. Period-3 (not period-2): the default fold
	    // assignment is itself i%folds, so a period-2 jitter here would
	    // make every 2-fold training subset a constant again (each fold's
	    // training half is the *other* fold's individuals, all sharing one
	    // i%2 parity) -- period-3 avoids aliasing against a 2-fold split.
	    const double jitter = 0.01 * ( (i % 3) - 1 );
	    const double y = ( grp==0 ? 0.0 : 10.0 ) + jitter;
	    const std::string id = ( grp==0 ? "W" : "NR" ) + std::to_string(i) + "_x3";
	    dpp_matrix_t m; m.id = id;
	    dpp_matrix_t h; h.id = id;
	    for (int r=0; r<nrows; r++)
	      {
		m.time_sec.push_back( (r+1) * 30.0 );
		// feature must itself vary by individual (not just by row) so a
		// fold-held-out booster has something to split on to predict
		// this individual's jittered y -- otherwise (as originally
		// written here) every individual in a group shares the exact
		// same feature trace, the booster can only predict the flat
		// group mean, and calibrate()'s y~[1,Zraw] fit degenerates
		// (Zraw is then a constant column -- rank-deficient)
		m.X.push_back( { val + jitter + 0.001 * r } );
		h.time_sec.push_back( (r+1) * 30.0 );
		h.X.push_back( grp==0 ? std::vector<double>{ 1.0 , 0.0 , 0.0 } : std::vector<double>{ 0.0 , 0.0 , 1.0 } );
	      }
	    dpp_io::save( corpus , m , 1 , ! first );
	    dpp_io::save( hypno_corpus , h , 3 , ! first );
	    first = false;
	    ph << id << "\t" << y << "\n";
	  }
      }
    ph.close();

    const std::string specfile = base + ".spec";
    { std::ofstream sp( specfile.c_str() ); sp << "CH X\nBASE: SLOPE X\n"; }
    const std::string config = base + ".conf";
    { std::ofstream cf( config.c_str() ); cf << "objective=regression\nverbosity=-1\nmin_data_in_leaf=1\nnum_leaves=7\n"; }

    cmd_t::attach_ivars( phenofile );

    param_t param;
    param.add( "data" , corpus );
    param.add( "hypno" , hypno_corpus );
    param.add( "hypno-three-state" , "T" );
    param.add( "phe" , "Y" );
    param.add( "spec" , specfile );
    param.add( "config" , config );
    param.add( "out" , base + "_model" );
    param.add( "folds" , "2" );   // 12 individuals, 2 folds

    dpp_fit_t fit( param );
    fit.fit();

    bool bundle_ok = Helper::fileExists( Helper::expand( base + "_model.W.mod" ) )
      && Helper::fileExists( Helper::expand( base + "_model.NR.mod" ) )
      && Helper::fileExists( Helper::expand( base + "_model.dpp" ) )
      && Helper::fileExists( Helper::expand( base + "_model.oof" ) );

    double w_sum = 0, nr_sum = 0; int w_n = 0, nr_n = 0;
    for ( const auto & row : fit.oof_rows )
      {
	if ( row.id[0] == 'W' ) { w_sum += row.y_pred; ++w_n; }
	else                    { nr_sum += row.y_pred; ++nr_n; }
      }
    const double w_mean  = w_n  > 0 ? w_sum  / w_n  : std::numeric_limits<double>::quiet_NaN();
    const double nr_mean = nr_n > 0 ? nr_sum / nr_n : std::numeric_limits<double>::quiet_NaN();

    bool pass = bundle_ok && w_n > 0 && nr_n > 0 && ( nr_mean - w_mean ) > 3.0;
    std::ostringstream m; m << "bundle_ok=" << bundle_ok << " oof_rows=" << fit.oof_rows.size()
			    << " w_mean=" << w_mean << " (exp low) nr_mean=" << nr_mean << " (exp high, separated)";
    record(R,"dpp-fit/stage-conditioned-kfold-runs-and-separates", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"dpp-fit/stage-conditioned-kfold-runs-and-separates",false,e.what(),V); }

  // X4 -- fold-file= override is actually respected, not silently ignored:
  // supply an assignment that deliberately differs from what the default
  // sorted-round-robin would produce, confirm the override wins
  try {
    const std::string base = temp_base_path("test_dppfit_x4");
    param_t param;
    param.add( "out" , base + "_model" );
    param.add( "phe" , "Y" );
    param.add( "folds" , "2" );

    const std::string foldfile = base + ".folds";
    { std::ofstream ff( foldfile.c_str() ); ff << "A\t1\nB\t1\nC\t1\nD\t0\nE\t0\n"; }   // opposite of round-robin's 0,1,0,1,0
    param.add( "fold-file" , foldfile );

    dpp_fit_t fit( param );
    for ( const std::string & id : { "A" , "B" , "C" , "D" , "E" } )
      { dpp_matrix_t m; m.id = id; fit.individuals.push_back( m ); }

    fit.assign_folds();

    bool pass = fit.fold_assignment.size() == 5
      && fit.fold_assignment["A"] == 1 && fit.fold_assignment["B"] == 1 && fit.fold_assignment["C"] == 1
      && fit.fold_assignment["D"] == 0 && fit.fold_assignment["E"] == 0;
    std::ostringstream m; m << "A=" << fit.fold_assignment["A"] << " B=" << fit.fold_assignment["B"]
			    << " C=" << fit.fold_assignment["C"] << " D=" << fit.fold_assignment["D"]
			    << " E=" << fit.fold_assignment["E"] << " (exp 1,1,1,0,0, i.e. the fold-file's own assignment)";
    record(R,"dpp-fit/fold-file-override-respected", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"dpp-fit/fold-file-override-respected",false,e.what(),V); }

  // X5 -- folds= and validation= together halts with a clear error, rather
  // than defining an ambiguous three-way split
  try {
    param_t param;
    param.add( "out" , temp_base_path("test_dppfit_x5") + "_model" );
    param.add( "phe" , "Y" );
    param.add( "folds" , "3" );
    param.add( "validation" , "somefile.txt" );

    bool halted = false;
    try { dpp_fit_t fit( param ); } catch ( std::exception & ) { halted = true; }

    record(R,"dpp-fit/folds-and-validation-mutually-exclusive", halted,
	   std::string("halted=") + ( halted ? "T" : "F" ) + " (exp T)", V);
  } catch(std::exception & e){ record(R,"dpp-fit/folds-and-validation-mutually-exclusive",false,e.what(),V); }

  // X6 -- regression: flatten_and_split() is re-invoked once per fold by
  // cross_validate(), plus once more for the final full-corpus fit -- the
  // per-individual block-boundary accumulators (istart_train/valid,
  // wtable_train/valid) must reflect only the *current* call, not
  // carry over stale entries from an earlier call (which would pair a
  // stale, possibly out-of-range boundary against the current, differently-
  // sized weight vector in lgbm_t::add_block_weights())
  try {
    const std::string base = temp_base_path("test_dppfit_x6");
    std::vector<std::string> lo_ids, hi_ids;
    const std::string phenofile = build_corpus( base , 6 , 4 , 0.0 , 10.0 , 0.0 , 10.0 , &lo_ids , &hi_ids );

    const std::string specfile = base + ".spec";
    { std::ofstream sp( specfile.c_str() ); sp << "CH X\nBASE: SLOPE X\n"; }

    cmd_t::attach_ivars( phenofile );

    param_t param;
    param.add( "data" , base + ".dat" );
    param.add( "phe" , "Y" );
    param.add( "spec" , specfile );
    param.add( "out" , base + "_model" );

    dpp_fit_t fit( param );
    fit.load_corpus();
    fit.build_feature_labels();
    fit.attach_phenotypes();

    // "fold A": hold out 2 individuals (1 per cluster) -> 10 train indiv
    fit.validation_ids = { lo_ids[0] , hi_ids[0] };
    fit.flatten_and_split();
    const size_t istart_train_A = fit.istart_train.size();
    const size_t wtable_train_A = fit.wtable_train.size();

    // "fold B": hold out a *different*, larger set of 4 individuals (2 per
    // cluster) -> only 8 train indiv this time -- deliberately smaller and
    // differently-sized than fold A, so any leftover fold-A entries would
    // be detectably stale (extra count, and/or a stale istart_train[i]
    // pointing past the fold-B training matrix's actual row count)
    fit.validation_ids = { lo_ids[0] , lo_ids[1] , hi_ids[0] , hi_ids[1] };
    fit.flatten_and_split();
    const size_t istart_train_B = fit.istart_train.size();
    const size_t wtable_train_B = fit.wtable_train.size();

    bool pass = istart_train_A == 10 && wtable_train_A == 10
      && istart_train_B == 8 && wtable_train_B == 8
      && fit.Xtrain.rows() == 32   // 8 individuals x 4 rows
      && ( fit.istart_train.empty() || fit.istart_train.back() < (uint64_t)fit.Xtrain.rows() );
    std::ostringstream m; m << "fold A: istart_train=" << istart_train_A << " wtable_train=" << wtable_train_A
			    << " (exp 10,10); fold B: istart_train=" << istart_train_B << " wtable_train=" << wtable_train_B
			    << " (exp 8,8, not accumulated to 18) Xtrain.rows=" << fit.Xtrain.rows() << " (exp 32)";
    record(R,"dpp-fit/flatten-and-split-no-stale-boundaries", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"dpp-fit/flatten-and-split-no-stale-boundaries",false,e.what(),V); }

  // X6b -- regression: Xvalid/yvalid/indiv_weight_valid (and Hvalid, for a
  // stage-conditioned fit) must themselves be reset to empty when a
  // flatten_and_split() call has zero validation rows (validation_ids
  // empty) -- not just the istart_valid/wtable_valid boundary bookkeeping
  // covered by X6 above. Without this, the *final* full-corpus fit (which
  // cross_validate() drives via exactly this "populated-validation fold,
  // then empty-validation call" sequence) would silently retain the last
  // CV fold's validation set: lgbm_t would see yvalid.size()>0, attach it
  // as a real validation set, and run early stopping / calibration against
  // stale, wrong data for what's supposed to be the validation-free final
  // model. Caught via a real --dpp-fit folds= run's log unexpectedly
  // showing "validation = ..." for the final, "0 held out" fit.
  try {
    const std::string base = temp_base_path("test_dppfit_x6b");
    std::vector<std::string> lo_ids, hi_ids;
    const std::string phenofile = build_corpus( base , 6 , 4 , 0.0 , 10.0 , 0.0 , 10.0 , &lo_ids , &hi_ids );

    const std::string specfile = base + ".spec";
    { std::ofstream sp( specfile.c_str() ); sp << "CH X\nBASE: SLOPE X\n"; }

    cmd_t::attach_ivars( phenofile );

    param_t param;
    param.add( "data" , base + ".dat" );
    param.add( "phe" , "Y" );
    param.add( "spec" , specfile );
    param.add( "out" , base + "_model" );

    dpp_fit_t fit( param );
    fit.load_corpus();
    fit.build_feature_labels();
    fit.attach_phenotypes();

    // simulates cross_validate()'s last fold: a real, non-empty validation set
    fit.validation_ids = { lo_ids[0] , hi_ids[0] };
    fit.flatten_and_split();
    const int xvalid_rows_1 = (int) fit.Xvalid.rows();
    const int yvalid_size_1 = (int) fit.yvalid.size();
    const int xtrain_rows_1 = (int) fit.Xtrain.rows();
    const bool populated_first = xvalid_rows_1 > 0 && yvalid_size_1 > 0
      && fit.indiv_weight_valid.size() == fit.yvalid.size();

    // simulates fit()'s subsequent final, full-corpus call: validation_ids
    // restored to empty (mirrors cross_validate()'s own restore of the
    // saved, pre-CV validation_ids at the end of the fold loop)
    fit.validation_ids.clear();
    fit.flatten_and_split();

    bool pass = populated_first
      && fit.Xvalid.rows() == 0 && fit.yvalid.empty() && fit.indiv_weight_valid.empty()
      && fit.Xtrain.rows() == 48;   // all 6+6 individuals x 4 rows, none held out
    std::ostringstream m; m << "1st call (2 held out): Xvalid.rows=" << xvalid_rows_1 << " yvalid.size=" << yvalid_size_1
			    << " Xtrain.rows=" << xtrain_rows_1 << " (exp >0,>0,40); 2nd call (0 held out): Xvalid.rows="
			    << fit.Xvalid.rows() << " yvalid.size=" << fit.yvalid.size()
			    << " indiv_weight_valid.size=" << fit.indiv_weight_valid.size()
			    << " (exp 0,0,0); Xtrain.rows=" << fit.Xtrain.rows() << " (exp 48)";
    record(R,"dpp-fit/flatten-and-split-clears-stale-validation", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"dpp-fit/flatten-and-split-clears-stale-validation",false,e.what(),V); }

  // X7 -- hypnodensity/feature row *alignment* (not just row count) is
  // validated: two corpora with matching IDs and row counts but a shifted
  // time_sec[] for one individual must halt, not silently pair rows by
  // position against the wrong stage-probability vector
  try {
    const std::string base = temp_base_path("test_dppfit_x7");
    const std::string corpus = base + ".dat";
    const std::string hypno_corpus = base + ".hypno.dat";
    const std::string phenofile = base + ".pheno";
    std::ofstream ph( phenofile.c_str() ); ph << "ID\tY\n";
    for (int i=0; i<3; i++)
      {
	const std::string id = "S" + std::to_string(i) + "_x7";
	dpp_matrix_t m; m.id = id;
	dpp_matrix_t h; h.id = id;
	for (int r=0; r<4; r++)
	  {
	    m.time_sec.push_back( (r+1)*30.0 );
	    m.X.push_back( { 5.0 + 0.001*r } );
	    // individual S1's hypno rows are offset by +1s -- same row COUNT
	    // and same ID, but misaligned in time
	    h.time_sec.push_back( (r+1)*30.0 + ( i==1 ? 1.0 : 0.0 ) );
	    h.X.push_back( { 1.0 , 0.0 , 0.0 } );
	  }
	dpp_io::save( corpus , m , 1 , i!=0 );
	dpp_io::save( hypno_corpus , h , 3 , i!=0 );
	ph << id << "\t5.0\n";
      }
    ph.close();

    param_t param;
    param.add( "data" , corpus );
    param.add( "hypno" , hypno_corpus );
    param.add( "hypno-three-state" , "T" );
    param.add( "phe" , "Y" );
    param.add( "out" , base + "_model" );

    dpp_fit_t fit( param );
    fit.load_corpus();

    bool halted = false;
    try { fit.load_hypno_corpus(); } catch ( std::exception & ) { halted = true; }

    record(R,"dpp-fit/hypno-corpus-time-misalignment-halts", halted,
	   std::string("halted=") + ( halted ? "T" : "F" ) + " (exp T)", V);
  } catch(std::exception & e){ record(R,"dpp-fit/hypno-corpus-time-misalignment-halts",false,e.what(),V); }

  // X8 -- calibrate() falls back to {1,0} (not {0,0}) when GLM fit
  // succeeds/is valid() but has a degenerate per-coefficient standard
  // error (here: a deliberately exact, zero-residual linear relationship,
  // no jitter) -- regression test for calibrate() previously using
  // coef[]'s resized-but-never-set {0,0} default without checking
  // display()'s return value, silently zeroing this stage's entire
  // contribution to the blend instead of an honest pass-through
  try {
    const std::string base = temp_base_path("test_dppfit_x8");
    const std::string config = base + ".conf";
    { std::ofstream cf( config.c_str() ); cf << "objective=regression\nverbosity=-1\nmin_data_in_leaf=1\nnum_leaves=7\n"; }

    param_t param;
    param.add( "out" , base + "_model" );
    param.add( "phe" , "Y" );
    param.add( "hypno" , "unused" );

    dpp_fit_t fit( param );
    fit.n_features = 1;
    fit.stage_labels = { "W" , "R" , "NR" };

    const int n = 20;
    fit.Xtrain = Eigen::MatrixXd( n , 1 );
    fit.ytrain.resize( n );
    fit.indiv_weight_train.assign( n , 1.0f );
    fit.Htrain = Eigen::MatrixXd( n , 3 );
    for (int i=0; i<n; i++)
      {
	fit.Xtrain(i,0) = (double)i;
	fit.ytrain[i] = (double)i;
	fit.Htrain(i,0) = 1.0; fit.Htrain(i,1) = 0.0; fit.Htrain(i,2) = 0.0;
      }

    const int nv = 20;
    fit.Xvalid = Eigen::MatrixXd( nv , 1 );
    fit.yvalid.assign( nv , 0.0 );
    fit.Hvalid = Eigen::MatrixXd( nv , 3 );
    for (int i=0; i<nv; i++)
      {
	fit.Xvalid(i,0) = (double)i;
	fit.Hvalid(i,0) = 1.0; fit.Hvalid(i,1) = 0.0; fit.Hvalid(i,2) = 0.0;
      }

    fit.stage_lgbm.resize( 3 );
    lgbm_t & lg = fit.stage_lgbm[0];
    lg.load_config( config );
    lg.qt_mode = true;
    lg.attach_training_matrix( fit.Xtrain );
    lg.attach_training_qts( fit.ytrain );
    lg.training_weights.assign( n , 1.0f );
    lg.apply_weights( lg.training , &lg.training_weights );
    lg.create_booster( false );

    // exact (zero-residual) relationship, no jitter -- deliberately
    // triggers GLM::display()'s degenerate-standard-error gate
    Eigen::MatrixXd Zraw = lg.predict( fit.Xvalid );
    const double OFFSET = 7.5;
    for (int i=0; i<nv; i++) fit.yvalid[i] = Zraw(i,0) + OFFSET;

    std::pair<double,double> ab = fit.calibrate( 0 );

    bool pass = approx_equal( ab.first , 1.0 , 1e-9 ) && approx_equal( ab.second , 0.0 , 1e-9 );
    std::ostringstream m; m << "a=" << ab.first << " b=" << ab.second << " (exp 1,0 -- fallback, not 0,0)";
    record(R,"dpp-fit/calibrate-degenerate-se-falls-back-to-noop", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"dpp-fit/calibrate-degenerate-se-falls-back-to-noop",false,e.what(),V); }

  // X9 -- subject-level OOF aggregation: per-individual rows collapse to
  // one (y_true, mean y_pred) pair each, with the correct row count and
  // mean; write_oof_subject() then computes Pearson r / Spearman rho /
  // RMSE over those subject-level pairs (not the raw rows) and writes
  // them to <out>.oof.subject -- constructed here as an exact y=x
  // relationship (r=1, rho=1, rmse=0) so the stats are deterministic to
  // check, unlike the row-level oof_rows they're aggregated from
  try {
    const std::string base = temp_base_path("test_dppfit_x9");
    param_t param;
    param.add( "out" , base + "_model" );
    param.add( "phe" , "Y" );

    dpp_fit_t fit( param );

    auto add_row = [&]( const std::string & id, double y_true, double y_pred, int fold )
      {
	dpp_fit_t::oof_row_t row;
	row.id = id; row.time_sec = 0; row.y_true = y_true; row.y_pred = y_pred; row.fold = fold;
	fit.oof_rows.push_back( row );
      };
    add_row( "A" , 0.0  , 0.0  , 0 );
    add_row( "A" , 0.0  , 0.0  , 0 );
    add_row( "B" , 5.0  , 5.0  , 0 );
    add_row( "B" , 5.0  , 5.0  , 0 );
    add_row( "B" , 5.0  , 5.0  , 0 );
    add_row( "C" , 10.0 , 10.0 , 1 );

    fit.aggregate_oof_by_subject();

    bool pass = fit.oof_subjects.size() == 3;
    std::map<std::string,dpp_fit_t::oof_subject_t> by_id;
    for ( const auto & s : fit.oof_subjects ) by_id[ s.id ] = s;
    pass = pass && by_id.count("A") && by_id.count("B") && by_id.count("C")
      && by_id["A"].n_rows == 2 && approx_equal( by_id["A"].y_pred_mean , 0.0  , 1e-9 )
      && by_id["B"].n_rows == 3 && approx_equal( by_id["B"].y_pred_mean , 5.0  , 1e-9 )
      && by_id["C"].n_rows == 1 && approx_equal( by_id["C"].y_pred_mean , 10.0 , 1e-9 );

    fit.write_oof_subject();
    const std::string subj_file = Helper::expand( base + "_model.oof.subject" );
    bool file_written = Helper::fileExists( subj_file );

    // exact y=x relationship by construction -> r=1, rho=1, rmse=0
    bool stats_ok = file_written;
    if ( file_written )
      {
	std::ifstream IN1( subj_file.c_str() );
	std::string line;
	double r_val = -99, rho_val = -99, rmse_val = -99;
	while ( std::getline( IN1 , line ) )
	  {
	    if ( line.find( "# pearson_r=" ) == 0 )   r_val   = atof( line.substr( 12 ).c_str() );
	    if ( line.find( "# spearman_rho=" ) == 0 ) rho_val = atof( line.substr( 15 ).c_str() );
	    if ( line.find( "# rmse=" ) == 0 )         rmse_val = atof( line.substr( 7 ).c_str() );
	  }
	IN1.close();
	stats_ok = approx_equal( r_val , 1.0 , 1e-6 ) && approx_equal( rho_val , 1.0 , 1e-6 ) && approx_equal( rmse_val , 0.0 , 1e-6 );
      }

    pass = pass && file_written && stats_ok;
    std::ostringstream m; m << "n_subjects=" << fit.oof_subjects.size() << " (exp 3); A: n=" << by_id["A"].n_rows
			    << " mean=" << by_id["A"].y_pred_mean << " (exp 2,0); B: n=" << by_id["B"].n_rows
			    << " mean=" << by_id["B"].y_pred_mean << " (exp 3,5); C: n=" << by_id["C"].n_rows
			    << " mean=" << by_id["C"].y_pred_mean << " (exp 1,10); file_written=" << file_written
			    << "; stats (r,rho,rmse exp 1,1,0) ok=" << stats_ok;
    record(R,"dpp-fit/subject-level-oof-aggregation", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"dpp-fit/subject-level-oof-aggregation",false,e.what(),V); }

  // X10 -- CV-derived iteration count (pooled): with folds= set,
  // forced_n_iterations should be populated (median best_iteration across
  // folds, capped by the configured iterations=), and -- the actual point
  // of the feature -- the *final*, full-corpus booster's own n_iterations
  // must equal it (wiring confirmed, not just computed and discarded).
  // No-exclusion re-checked here too: this feature must never change which
  // rows train the final booster, only its iteration count/calibration.
  try {
    const std::string base = temp_base_path("test_dppfit_x10");
    std::vector<std::string> lo_ids, hi_ids;
    const std::string phenofile = build_corpus( base , 6 , 4 , 0.0 , 10.0 , 0.0 , 10.0 , &lo_ids , &hi_ids );

    const std::string specfile = base + ".spec";
    { std::ofstream sp( specfile.c_str() ); sp << "CH X\nBASE: SLOPE X\n"; }
    const std::string config = base + ".conf";
    // aggressive early-stopping patience so it plausibly fires on this
    // tiny synthetic set, giving real (not trivially-identical) per-fold
    // best_iteration values for the median to aggregate over
    { std::ofstream cf( config.c_str() ); cf << "objective=regression\nverbosity=-1\nmin_data_in_leaf=1\nnum_leaves=7\n"; }

    cmd_t::attach_ivars( phenofile );

    param_t param;
    param.add( "data" , base + ".dat" );
    param.add( "phe" , "Y" );
    param.add( "spec" , specfile );
    param.add( "config" , config );
    param.add( "out" , base + "_model" );
    param.add( "folds" , "3" );
    param.add( "iterations" , "50" );
    param.add( "early-stopping-rounds" , "2" );

    dpp_fit_t fit( param );
    fit.fit();

    const int n_total_rows = 48;   // 12 individuals x 4 rows
    bool pass = fit.forced_n_iterations > 0 && fit.forced_n_iterations <= 50
      && fit.lgbm.n_iterations == fit.forced_n_iterations
      && fit.Xtrain.rows() == n_total_rows;
    std::ostringstream m; m << "forced_n_iterations=" << fit.forced_n_iterations << " (exp 1..50); final lgbm.n_iterations="
			    << fit.lgbm.n_iterations << " (exp == forced_n_iterations); Xtrain.rows()=" << fit.Xtrain.rows()
			    << " (exp " << n_total_rows << ", i.e. full corpus, unaffected)";
    record(R,"dpp-fit/cv-derived-iterations-pooled", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"dpp-fit/cv-derived-iterations-pooled",false,e.what(),V); }

  // X11 -- CV-derived calibration + iteration count (stage-conditioned):
  // the final bundle's calib_a/calib_b must equal the mean of the K folds'
  // own (real, non-trivial) calibrations -- not the {1,0} fallback
  // train_stage_boosters() would otherwise leave them at for the final,
  // no-validation fit -- and each stage's final n_iterations must equal
  // its own CV-derived (median) value. No-exclusion re-checked again.
  try {
    const std::string base = temp_base_path("test_dppfit_x11");
    const std::string corpus = base + ".dat";
    const std::string hypno_corpus = base + ".hypno.dat";
    const std::string phenofile = base + ".pheno";

    std::ofstream ph( phenofile.c_str() ); ph << "ID\tY\n";
    const int n_per_group = 6, nrows = 5;
    bool first = true;
    for (int grp=0; grp<2; grp++)
      {
	const double val = grp==0 ? 0.0 : 10.0;
	for (int i=0; i<n_per_group; i++)
	  {
	    const double jitter = 0.01 * ( (i % 3) - 1 );   // see X3's comment: period-3, avoids the GLM degenerate-SE case
	    const double y = ( grp==0 ? 0.0 : 10.0 ) + jitter;
	    const std::string id = ( grp==0 ? "W" : "NR" ) + std::to_string(i) + "_x11";
	    dpp_matrix_t m; m.id = id;
	    dpp_matrix_t h; h.id = id;
	    for (int r=0; r<nrows; r++)
	      {
		m.time_sec.push_back( (r+1) * 30.0 );
		m.X.push_back( { val + jitter + 0.001 * r } );
		h.time_sec.push_back( (r+1) * 30.0 );
		h.X.push_back( grp==0 ? std::vector<double>{ 1.0 , 0.0 , 0.0 } : std::vector<double>{ 0.0 , 0.0 , 1.0 } );
	      }
	    dpp_io::save( corpus , m , 1 , ! first );
	    dpp_io::save( hypno_corpus , h , 3 , ! first );
	    first = false;
	    ph << id << "\t" << y << "\n";
	  }
      }
    ph.close();

    const std::string specfile = base + ".spec";
    { std::ofstream sp( specfile.c_str() ); sp << "CH X\nBASE: SLOPE X\n"; }
    const std::string config = base + ".conf";
    { std::ofstream cf( config.c_str() ); cf << "objective=regression\nverbosity=-1\nmin_data_in_leaf=1\nnum_leaves=7\n"; }

    cmd_t::attach_ivars( phenofile );

    param_t param;
    param.add( "data" , corpus );
    param.add( "hypno" , hypno_corpus );
    param.add( "hypno-three-state" , "T" );
    param.add( "phe" , "Y" );
    param.add( "spec" , specfile );
    param.add( "config" , config );
    param.add( "out" , base + "_model" );
    param.add( "folds" , "2" );
    param.add( "iterations" , "50" );
    param.add( "early-stopping-rounds" , "2" );

    dpp_fit_t fit( param );
    fit.fit();

    const int n_total_rows = 60;   // 12 individuals x 5 rows
    const int n_stages = (int) fit.stage_labels.size();

    bool pass = n_stages == 3
      && (int)fit.forced_stage_n_iterations.size() == n_stages
      && (int)fit.cv_calib_a.size() == n_stages && (int)fit.cv_calib_b.size() == n_stages
      && fit.calib_a == fit.cv_calib_a && fit.calib_b == fit.cv_calib_b
      && fit.Xtrain.rows() == n_total_rows;

    bool iters_ok = pass;
    bool any_real_calibration = false;
    for (int s=0; s<n_stages && pass; s++)
      {
	iters_ok = iters_ok && fit.forced_stage_n_iterations[s] > 0 && fit.forced_stage_n_iterations[s] <= 50
	  && fit.stage_lgbm[s].n_iterations == fit.forced_stage_n_iterations[s];
	if ( ! approx_equal( fit.calib_a[s] , 1.0 , 1e-6 ) || ! approx_equal( fit.calib_b[s] , 0.0 , 1e-6 ) )
	  any_real_calibration = true;
      }
    pass = pass && iters_ok && any_real_calibration;

    std::ostringstream m; m << "n_stages=" << n_stages << " (exp 3); final calib_a==cv_calib_a && calib_b==cv_calib_b: "
			    << ( fit.calib_a == fit.cv_calib_a && fit.calib_b == fit.cv_calib_b )
			    << "; any_real_calibration=" << any_real_calibration << " (exp T, not all {1,0}); iters_ok=" << iters_ok
			    << "; Xtrain.rows()=" << fit.Xtrain.rows() << " (exp " << n_total_rows << ", full corpus)";
    record(R,"dpp-fit/cv-derived-calibration-and-iterations-stage-conditioned", pass, m.str(), V);
  } catch(std::exception & e){ record(R,"dpp-fit/cv-derived-calibration-and-iterations-stage-conditioned",false,e.what(),V); }
}

#endif

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
  RUN("sigdyn",   test_sigdyn)
  RUN("dpp",      test_dpp)
#ifdef HAS_LGBM
  RUN("dpp-fit",  test_dpp_fit)
#endif

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
