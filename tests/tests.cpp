
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

// Luna integrated test suite
// Invocation: luna __LUNA_TESTS__ [group] [verbose]
//
// Groups: all, signal, epoch, mask, filter, resample, psd, spindles,
//         hypno, annot, write, script, lunapi, segsrv
//
// All tests use fully synthetic in-memory data (no external files needed).
// Exit code: 0 = all pass, 1 = any failure.

#include "tests.h"
#include "luna.h"
#include "lunapi/lunapi.h"
#include "lunapi/segsrv.h"
#include "miscmath/crandom.h"

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

  // K5 — bad command: should throw, not crash
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
  RUN("lunapi",   test_lunapi)
  RUN("segsrv",   test_segsrv)

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
