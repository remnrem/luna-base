
// TODO
//  --> |skew| for polarity invariant features / check all other features are polarity invariant
//  --> added phys |median| & IQR (polarity invariant)
//  --> select up to N of each epochs class only

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

#ifdef HAS_LGBM

#include "dsp/ctypes.h"
#include "edf/edf.h"
#include "dsp/resample.h"
#include "edf/slice.h"

#include "dsp/wrappers.h"
#include "dsp/acf.h"

#include "db/db.h"
#include "param.h"

#include <set>
#include <random>

extern logger_t logger;
extern writer_t writer;

// move to MiscMath

#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <limits>
#include <algorithm>
#include <unordered_set>

lgbm_t ctypes_t::lgbm;

std::string ctypes_t::lgbm_model_loaded = "";

std::vector<std::string> ctypes_t::varlist;
std::vector<std::string> ctypes_t::var_root;
std::vector<std::string> ctypes_t::var_hz;
std::vector<std::string> ctypes_t::var_trans;
std::vector<std::string> ctypes_t::clslist;



// ---- helpers ----

std::vector<double> ctypes_t::make_diff(const std::vector<double>& x)
{
  const size_t n = x.size();  

  std::vector<double> d(n); // do not change size
  if ( n == 0 ) return d;
  d[0] = 0; // do not change size

  for (size_t i = 1; i < n; ++i)
    d[i - 1] = x[i] - x[i - 1];
  
  return d;
}


#include <set>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>

static inline bool is_nan_d(double v) { return std::isnan(v); }


double ctypes_t::mean(const std::vector<double>& x)
{
  if (x.empty()) return 0.0;
  return std::accumulate(x.begin(), x.end(), 0.0) / (double)x.size();
}

double ctypes_t::sd_sample(const std::vector<double>& x)
{
  const size_t n = x.size();
  if (n < 2) return 0.0;
  const double mu = mean(x);
  double ss = 0.0;
  for (double v : x) {
    const double d = v - mu;
    ss += d * d;
  }
  return std::sqrt(ss / (double)(n - 1));
}

double ctypes_t::mean_ptr(const double* x, size_t n)
{
  if (n == 0) return 0.0;
  
  double s = 0.0;
  for (size_t i = 0; i < n; ++i)
    s += x[i];
  
  return s / (double)n;
}

double ctypes_t::sd_sample_ptr(const double* x, size_t n)
{
  if (n < 2) return 0.0;
  
  const double mu = mean_ptr(x, n);
  double ss = 0.0;
  
  for (size_t i = 0; i < n; ++i) {
    const double d = x[i] - mu;
    ss += d * d;
  }
  
  return std::sqrt(ss / (double)(n - 1));
}

  
// quantile using nth_element (no interpolation; fast and sufficient for QC)
double ctypes_t::quantile(std::vector<double> xs, double p)
{
  const size_t n = xs.size();
  if (n == 0) return 0.0;
  if (p <= 0.0) return *std::min_element(xs.begin(), xs.end());
  if (p >= 1.0) return *std::max_element(xs.begin(), xs.end());
  const size_t k = (size_t)std::llround(p * (double)(n - 1));
  std::nth_element(xs.begin(), xs.begin() + (long)k, xs.end());
  return xs[k];
}

double ctypes_t::mad(const std::vector<double>& x)
{
  if (x.empty()) return 0.0;
  std::vector<double> xs = x;
  const double med = quantile(xs, 0.5);
  
  std::vector<double> ad;
  ad.reserve(x.size());
  for (double v : x) ad.push_back(std::abs(v - med));
  return quantile(ad, 0.5); // unscaled MAD (constant=1)
}

// sign helper for ZCR; treats 0 as "no sign"
int ctypes_t::sgn(double v, double eps )
{
  if (v > eps) return +1;
  if (v < -eps) return -1;
  return 0;
}

// ---- requested metrics ----

double ctypes_t::line_length(const std::vector<double>& x)
{
  const size_t n = x.size();
  if (n < 2) return 0.0;
  double ll = 0.0;
  for (size_t i = 1; i < n; ++i) ll += std::abs(x[i] - x[i - 1]);
  return ll;
}

// normalized line length: mean(|dx|) / SD(x)  (unitless)
double ctypes_t::line_length_norm_sd(const std::vector<double>& x)
{
  const size_t n = x.size();
  if (n < 2) return 0.0;
  const double s = sd_sample(x);
  if (!(std::isfinite(s) && s > 1e-12)) return 0.0;
  return line_length(x) / ((double)(n - 1) * s);
}

// zero-crossing rate (fraction of sign changes); ignores exact zeros by carrying last nonzero sign
double ctypes_t::zero_cross_rate(const std::vector<double>& x, double eps )
{
  const size_t n = x.size();
  if (n < 2) return 0.0;
  
  int prev = 0;
  // find first non-zero sign
  size_t i0 = 0;
  for (; i0 < n; ++i0) {
    prev = sgn(x[i0], eps);
    if (prev != 0) break;
  }
  if (prev == 0) return 0.0;
  
  size_t crossings = 0;
  size_t steps = 0;
  for (size_t i = i0 + 1; i < n; ++i) {
    const int cur = sgn(x[i], eps);
    if (cur == 0) continue;         // ignore zeros
    ++steps;
    if (cur != prev) ++crossings;
    prev = cur;
  }
  if (steps == 0) return 0.0;
  return (double)crossings / (double)steps;
}

// flatline fraction: fraction of samples with |dx| < eps, where eps is relative to MAD by default
double ctypes_t::flatline_fraction(const std::vector<double>& x)
{


  const size_t n = x.size();


  // Cheap invariants; if these fail you're already in UB land.
  if (n > 0) {
    // Touch first/last to trigger crash earlier if x is invalid
    volatile double a = x.front();
    volatile double b = x.back();
    (void)a; (void)b;
  }


  if (n < 2) return 1.0;

  const double m = mad(x);

  const double eps = (m > 0.0) ? 1e-6 * m : 1e-12;
  
  size_t flat = 0;
  for (size_t i = 1; i < n; ++i)
    if (std::abs(x[i] - x[i - 1]) < eps) ++flat;
  
  return (double)flat / (double)(n - 1);
}

// clip fraction: fraction near empirical rails (quantile-based)
// rails: qlo=0.001, qhi=0.999; band: delta = 0.01*(qhi-qlo)
double ctypes_t::clip_fraction(const std::vector<double>& x)
{
  const size_t n = x.size();
  if (n == 0) return 0.0;
  
  std::vector<double> xs = x;
  const double qlo = quantile(xs, 0.001);
  const double qhi = quantile(xs, 0.999);
  const double range = qhi - qlo;
  if (!(std::isfinite(range) && range > 0.0)) return 0.0;
  
  const double delta = 0.01 * range;
  
  size_t clipped = 0;
  for (double v : x)
    if (v <= qlo + delta || v >= qhi - delta) ++clipped;
  
  return (double)clipped / (double)n;
}

// Quantize to an integer bin: q = round(x / eps)
// Robust to overflow; ignores non-finite values at call sites.
long long ctypes_t::qbin(double v, double eps )
{
  // eps must be > 0; caller enforces
  const double z = v / eps;
  
  // Guard overflow before rounding/casting
  const double hi = (double)std::numeric_limits<long long>::max();
  const double lo = (double)std::numeric_limits<long long>::min();
  
  if (z >= hi) return std::numeric_limits<long long>::max();
  if (z <= lo) return std::numeric_limits<long long>::min();
  
  return (long long)std::llround(z);
}

// 1) unique_frac: fraction of unique quantized values among finite samples
double ctypes_t::unique_frac_q(const std::vector<double>& x, double eps )
{
  if (x.empty() || !(eps > 0.0)) return 0.0;
  
  std::unordered_set<long long> bins;
  bins.reserve(x.size());
  
  int n_valid = 0;
  for (double v : x) {
    if (!std::isfinite(v)) continue;
    ++n_valid;
    bins.insert(qbin(v, eps));
  }
  
  if (n_valid == 0) return 0.0;
  return (double)bins.size() / (double)n_valid;
}

// 2) most_common_value_frac: modal quantized-bin fraction among finite samples
double ctypes_t::most_common_value_frac_q(const std::vector<double>& x, double eps )
{
  if (x.empty() || !(eps > 0.0)) return 0.0;
  
  std::unordered_map<long long, int> counts;
  counts.reserve(x.size());
  
  int n_valid = 0;
  for (double v : x) {
    if (!std::isfinite(v)) continue;
    ++n_valid;
    ++counts[qbin(v, eps)];
  }
  
  if (n_valid == 0) return 0.0;
  
  int max_count = 0;
  for (const auto& kv : counts)
    if (kv.second > max_count) max_count = kv.second;
  
  return (double)max_count / (double)n_valid;
}

// 3) transition_rate: fraction of transitions between successive finite samples,
// where "transition" means quantized bin changes.
double ctypes_t::transition_rate_q(const std::vector<double>& x, double eps )
{
  if (x.size() < 2 || !(eps > 0.0)) return 0.0;
  
  bool have_prev = false;
  long long prev = 0;
  
  int valid_steps = 0;  // number of adjacent finite-sample pairs considered
  int transitions = 0;
  
  for (size_t i = 0; i < x.size(); ++i) {
    const double v = x[i];
    if (!std::isfinite(v)) continue;
    
    const long long cur = qbin(v, eps);
    
    if (!have_prev) {
      prev = cur;
      have_prev = true;
      continue;
    }
    
    ++valid_steps;
    if (cur != prev) ++transitions;
    prev = cur;
  }
  
  if (valid_steps == 0) return 0.0;
  return (double)transitions / (double)valid_steps;
}



// Robust per-second proxy signals:
//  - x1[t] = median of samples in second t  (trend proxy)
//  - a1[t] = median of |samples| in second t (envelope proxy)
//
// Notes:
//  - Works for arbitrary sr (including non-integer like 31.25).
//  - Uses floor(sr) samples per second with an error-accumulator so the
//    long-run rate matches sr (i.e., some seconds get n or n+1 samples).
//  - Does NOT invent samples (no interpolation).
//  - If sr < 2, we treat the input as already "slow": copy x -> x1 and |x| -> a1.
//  - If x contains NaN/Inf, those are ignored within each second. If a second
//    has no finite samples, outputs NaN for that second.


double ctypes_t::median_inplace(std::vector<double>& v) {
  if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
  const size_t n = v.size();
  const size_t mid = n / 2;
  std::nth_element(v.begin(), v.begin() + mid, v.end());
  double m = v[mid];
  if (n % 2 == 0) {
    auto max_it = std::max_element(v.begin(), v.begin() + mid);
    m = 0.5 * (m + *max_it);
  }
  return m;
}

void ctypes_t::make_1s(const std::vector<double>& x,
		      double sr,
		      std::vector<double>* x1,
		      std::vector<double>* a1)
{

  if (!x1 || !a1) return;
  x1->clear();
  a1->clear();
  
  if (!(sr > 0.0) || x.empty()) return;
  
  // If already very slow, don't try to "bin"; just pass through.
  // This keeps comparability and avoids pathological 0/1 samples-per-bin behavior.

  if (sr < 2.0) {
    x1->reserve(x.size());
    a1->reserve(x.size());
    for (double v : x) {
      if (is_finite(v)) {
        x1->push_back(v);
        a1->push_back(std::fabs(v));
      } else {
        x1->push_back(std::numeric_limits<double>::quiet_NaN());
        a1->push_back(std::numeric_limits<double>::quiet_NaN());
      }
    }
    return;
  }

  // Binning policy:
  // Use integer sample counts per "1 second" bin that average to sr over time.
  // Example sr=31.25 => base=31 with +1 sample added in 25% of bins.

  const int base = std::max(1, static_cast<int>(std::floor(sr)));
  const double frac = sr - static_cast<double>(base);
  
  std::vector<double> buf;   buf.reserve(static_cast<size_t>(std::ceil(sr)));
  std::vector<double> abuf;  abuf.reserve(static_cast<size_t>(std::ceil(sr)));
  
  size_t i = 0;
  double err = 0.0;
  
  while (i < x.size()) {
    int n = base;
    err += frac;
    if (err >= 1.0) { n += 1; err -= 1.0; }
    
    // Clamp to remaining samples.
    const size_t remain = x.size() - i;
    if (static_cast<size_t>(n) > remain) n = static_cast<int>(remain);
    
    buf.clear();
    abuf.clear();
    buf.reserve(n);
    abuf.reserve(n);
    
    for (int k = 0; k < n; ++k) {
      const double v = x[i + static_cast<size_t>(k)];
      if (is_finite(v)) {
        buf.push_back(v);
        abuf.push_back(std::fabs(v));
      }
    }

    if (buf.empty()) {
      x1->push_back(std::numeric_limits<double>::quiet_NaN());
      a1->push_back(std::numeric_limits<double>::quiet_NaN());
    } else {
      x1->push_back(median_inplace(buf));
      a1->push_back(median_inplace(abuf));
    }
    
    i += static_cast<size_t>(n);
  }
}



void spectral_stats_t::compute(const std::vector<double>& x, double sr)
{

  const bool slow_sr =  sr < 32 ;
  
  std::vector<double> f, p;

  // expecting either 128 Hz (per-epoch analysis)
  //  or a 1 Hz (whole-night) signal - adjust windows accordinaly
  const double segment_sec = slow_sr ? 128 : 4 ;
  const double overlap_sec = slow_sr ?  96 : 2 ; 
  
  bool okay = dsptools::welch(x, sr, f, p, segment_sec , overlap_sec ); // single PSD call

  if (f.size() != p.size() || f.size() < 3 || ! okay )
    Helper::halt( "welch() returned invalid arrays");
  
  // Keep non-negative freqs only, enforce ascending (assumed)
  // Also guard against negative/NaN power.
  double Ptot = 0.0;
  for (size_t i = 0; i < p.size(); ++i)
    if (f[i] >= 0.0 && std::isfinite(p[i]) && p[i] >= 0.0)
      Ptot += p[i];
  
  if (!(Ptot > 0.0))
    {
      // All zeros / invalid power. Treat as degenerate.
      centroid = edge90 = bandwidth = std::numeric_limits<double>::quiet_NaN();
      flatness = std::numeric_limits<double>::quiet_NaN();
      entropy  = std::numeric_limits<double>::quiet_NaN();
      high = low = std::numeric_limits<double>::quiet_NaN();
      return;
    }
  
  // --- centroid (power-weighted mean frequency) ---
  double num = 0.0;
  for (size_t i = 0; i < p.size(); ++i)
    if (f[i] >= 0.0 && std::isfinite(f[i]) && std::isfinite(p[i]) && p[i] >= 0.0)
      num += f[i] * p[i];
  centroid = num / Ptot;
  
  // --- edge90 (frequency below which 90% of total power lies) ---
  double cum = 0.0;
  const double target = 0.90 * Ptot;
  edge90 = f.back();
  for (size_t i = 0; i < p.size(); ++i)
    {
      if (!(f[i] >= 0.0) || !std::isfinite(p[i]) || p[i] < 0.0) continue;
      cum += p[i];
      if (cum >= target)
	{
	  edge90 = f[i];
	  break;
	}
    }
  
  
  // --- bandwidth (2nd central moment around centroid) ---
  double var = 0.0;
  for (size_t i = 0; i < p.size(); ++i)
    {
      if (!(f[i] >= 0.0) || !std::isfinite(f[i]) || !std::isfinite(p[i]) || p[i] < 0.0) continue;
      const double df = f[i] - centroid;
      var += (df * df) * p[i];
    }
  bandwidth = std::sqrt(var / Ptot);
  
  // --- flatness (geometric mean / arithmetic mean of power) ---
  // Use a floor to avoid log(0). Flatness near 1 => noise-like; near 0 => peaky/tonal.
  const double eps = 1e-30;
    double log_sum = 0.0;
    size_t k = 0;
    for (size_t i = 0; i < p.size(); ++i)
      {
	if (!(f[i] >= 0.0) || !std::isfinite(p[i]) || p[i] < 0.0) continue;
	log_sum += std::log(p[i] + eps);
	++k;
      }
    if (k == 0)
      flatness = 0.0;
    else
      {
	const double geo = std::exp(log_sum / static_cast<double>(k));
      const double ari = Ptot / static_cast<double>(k);
      flatness = (ari > 0.0) ? (geo / ari) : 0.0;
      }

    // --- spectral entropy (normalized Shannon entropy of power distribution) ---
    // Compute on normalized p_i = P(f_i)/Ptot over non-negative freqs.
    double H = 0.0;
    size_t m = 0;
    for (size_t i = 0; i < p.size(); ++i)
    {
      if (!(f[i] >= 0.0) || !std::isfinite(p[i]) || p[i] < 0.0) continue;
      const double pi = p[i] / Ptot;
      if (pi > 0.0) H += -pi * std::log(pi);
      ++m;
    }
    // Normalize to [0,1] by dividing by log(m) (max entropy for uniform distribution)
    if (m > 1)
      entropy = H / std::log(static_cast<double>(m));
    else
      entropy = 0.0;

    // --- highlow: high-band power / low-band power ---
    // Defaults chosen for PSG-like work; adjust as desired.
    // Low:  0.1–2 Hz (resp/drift-ish)
    // High: 20–45 Hz (EMG-ish at Fs>=100)
    // Mid(ref): 2 - 15 Hz 
    // If Nyquist is below 45, this automatically clamps via f grid.
    const double lo1 = 0.1, lo2 = 2.0;
    const double mi1 = 2.0 , mi2 = 15.0;
    const double hi1 = 20.0, hi2 = 45.0;

    double Plo = 0.0, Pmi = 0.0, Phi = 0.0;
    for (size_t i = 0; i < p.size(); ++i)
    {
      if (!(f[i] >= 0.0) || !std::isfinite(f[i]) || !std::isfinite(p[i]) || p[i] < 0.0) continue;
      const double fi = f[i];
      if (fi >= lo1 && fi < lo2) Plo += p[i];
      if (fi >= mi1 && fi < mi2) Pmi += p[i];      
      if (fi >= hi1 && fi < hi2) Phi += p[i];
    }
    // Avoid division-by-zero; return 0 if no low-band power.
    high = (Pmi > 0.0) ? (Phi / Pmi) : 0.0;
    low = (Pmi > 0.0) ? (Plo / Pmi) : 0.0;


    // clamp if slow-SR
    if ( slow_sr )
      {
	high = low = std::numeric_limits<double>::quiet_NaN();
	flatness = edge90 = std::numeric_limits<double>::quiet_NaN();
      }
    
}



ctypes_t::ctypes_t( edf_t & edf , param_t & param)
{

  proc( edf , param );
    
}


//
// compute features, and then output and/or predict 
//

void ctypes_t::proc( edf_t & edf , param_t & param )
{

  //
  // options
  //

  // make predictions (default, w/ libray ct1
  
  const bool make_predictions = param.yesno( "predict" , true, true ); 
  
  const std::string model_path = param.has( "path" ) ? param.value( "path" ) : ".";
  
  const std::string model_lib  = param.has( "lib" ) ? param.value( "lib" ) : "ct1";

  const bool ignore_staging    = param.yesno( "ignore-staging" , false , true );

  // 10 mins (20 epochs of each of four stagea stage )
  // else total 40 mins
  const int sel_num_epochs     = param.has( "num-epoch" ) ? param.requires_int( "num-epoch" ) : ( ignore_staging ? 80 : 20 ) ;

  
  if ( make_predictions )
    attach_model( model_path, model_lib );
	
  
  //
  // signals  < 32 Hz --> resample to 1 Hz
  //         >= 32 Hz --> resample to 128 Hz
  //
  
  const double Fs_thresh = param.has( "fs-min" ) ? param.requires_dbl( "fs-min" ) : 32.0;

  const double Fs_norm   = param.has( "fs" )     ? param.requires_dbl( "fs" )      : 128.0;
  const double Fs_low    = param.has( "fs-low" ) ? param.requires_dbl( "fs-low" )      : 1.0;

  //
  // epoch durations
  //

  const double e1dur = 300.0;
  const double e1inc = 60.0;

  const double e128dur = 30.0;
  const double e128inc = 30.0;


  // to calculate epoch stats:  mean e95 must be > 0.5 Hz OR at least 20% of epochs must have e95 > 0.5
  edge95_mean  = param.has( "edge95-mean" )  ? param.requires_dbl( "edge95-mean" )  : 0.5 ;
  edge95_th    = param.has( "edge95-th" )    ? param.requires_dbl( "edge95-th" )    : 0.5 ;
  edge95_prop  = param.has( "edge95-prop" )  ? param.requires_dbl( "edge95-prop" )  : 0.2 ; 

  // to calculate whole-night stats: mean e95 must be < 2 Hz 
  edge95_mean2 = param.has( "edge95-mean2" ) ? param.requires_dbl( "edge95-mean2" ) : 2 ;
  edge95_th2   = param.has( "edge95-th2" )   ? param.requires_dbl( "edge95-th2" )   : 2 ;
  edge95_prop2 = param.has( "edge95-prop2" ) ? param.requires_dbl( "edge95-prop2" ) : 0.2 ; 
    
  //
  // options
  //

  const bool epoch_output = param.yesno( "epoch" , false , true );

  // if not in predict mode, then we want CH-level output by default
  const bool output = param.yesno( "output" , make_predictions ? false : true , true );
  
  //
  // attach LGBM model as needed
  //

  if ( make_predictions )
    {
      logger << "  attaching LGBM model " << model_lib << " from " << model_path << "\n";
    }

  
  //
  // add temp dummy 128 Hz signal to determine sample->epoch mappings (30 seconds, no overlap)
  //
  
  edf.init_signal( "__d128__" , 128 );
  int dslot = edf.header.signal( "__d128__" );
  
  std::vector<std::pair<int,int> > epoch2samples128;
  std::vector<int> epochs128;
  
  int s1 = 0;
  int s2 = 128 * 30;
  
  edf.timeline.set_epoch( e128dur , e128inc );  
  const int ne = edf.timeline.first_epoch();
  while ( 1 )
    {
      int epoch = edf.timeline.next_epoch();
      if ( epoch == -1 ) break;
      interval_t interval = edf.timeline.epoch( epoch );
      // get smps too from dummy 128-Hz signal
      slice_t slice( edf , dslot , interval , 1 , false, true );
      std::vector<int> * sp = slice.nonconst_psmps();
      if ( ! sp->empty() ) {
	//epoch2samples128.push_back( std::make_pair( sp->front() , sp->back() ) );
	epoch2samples128.push_back( std::make_pair( s1, s2 ) );
	s1 += 128 * 30;
	s2 += 128 * 30;
	epochs128.push_back( edf.timeline.display_epoch( epoch ) ); // track /display/ epoch codes
      }
    }

  // get wholetrace length for a 128 Hz too

  slice_t slice( edf , dslot , edf.timeline.wholetrace() );
  const std::vector<double> * d = slice.pdata();
  const size_t n128 = d->size();
  
  edf.drop_signal( dslot );

  //
  // add temp dummy 1 Hz signal to determine sample->epoch mappings (300 seconds, 60 seconds overlap)
  //
  
  edf.init_signal( "__d1__" , 1 );
  dslot = edf.header.signal( "__d1__" );
  
  std::vector<std::pair<int,int> > epoch2samples1;
  std::vector<int> epochs1;

  // reset epochs
  edf.timeline.set_epoch( e1dur , e1inc );

  const int ne1 = edf.timeline.first_epoch();
  
  s1 = 0;
  s2 = 300;
  
  while ( 1 )
    {
      int epoch = edf.timeline.next_epoch();
      if ( epoch == -1 ) break;
      interval_t interval = edf.timeline.epoch( epoch );
      // get smps too from dummy 1-Hz signal
      slice_t slice( edf , dslot , interval , 1 , false, true );
      std::vector<int> * sp = slice.nonconst_psmps();
      if ( ! sp->empty() ) {
	//std::cout << " originals(1Hz) = " << sp->front()  << " " << sp->back() << "\n";
	//std::cout << " originals(1Hz) = " << s1  << " " << s2 << "\n";
	//epoch2samples1.push_back( std::make_pair( sp->front() , sp->back() ) );
	epoch2samples1.push_back( std::make_pair( s1, s2 ) ) ;
	s1 += 60;
	s2 += 60;
	epochs1.push_back( edf.timeline.display_epoch( epoch ) ); // track /display/ epoch codes
      }
    }

  // get wholetrace length for a 1 Hz too

  slice_t slice1( edf , dslot , edf.timeline.wholetrace() );
  const std::vector<double> * d1 = slice1.pdata();
  const size_t n1 = d1->size();
  //  std::cout << " n1 = " << n1 << "\n";
  edf.drop_signal( dslot );

  
  //
  // get signals
  //

  std::string signal_label = param.requires( "sig" );
    
  const bool NO_ANNOTS = true;
    
  signal_list_t signals = edf.header.signal_list( signal_label , NO_ANNOTS );
  
  const int ns = signals.size();
  
  if ( ns == 0 ) return; 

  std::vector<double> Fs_orig = edf.header.sampling_freq( signals );

    
  //
  // get staging
  //

  s_n1.clear(); s_n2.clear(); s_n3.clear(); s_rem.clear();
  
  // reset to 30-second epochs
  edf.timeline.set_epoch( e128dur , e128inc );

  edf.timeline.first_epoch();

  std::map<int,int> stage;
  std::set<int> selected_epochs;
  
  if ( ! ignore_staging )
    {
      
      edf.annotations->make_sleep_stage( edf.timeline );
      
      // valid?
      bool has_staging = edf.timeline.hypnogram.construct( &(edf.timeline) , param , false );
      
      // valid, but empty?
      if ( has_staging && edf.timeline.hypnogram.empty() )
	has_staging = false;
        
      // exclude gaps from stages epoch count
      if ( has_staging )
	{
	  std::vector<int> ungapped;
	  for (int i=0; i<edf.timeline.hypnogram.stages.size(); i++)
	    if ( edf.timeline.hypnogram.stages[i] != GAP )
	      ungapped.push_back( edf.timeline.hypnogram.stages[i] );
	  
	  // total number of epochs does not match?
	  if ( ne != ungapped.size() )
	    Helper::halt( "problem extracting stage information: " 
			  + Helper::int2str( ne ) + " epochs but found stage info for " 
			  + Helper::int2str( (int)ungapped.size()) );
	  
	  int idx = 0;
	  edf.timeline.first_epoch();
	  while ( 1 )
	    {
	      int epoch = edf.timeline.next_epoch();
	      if ( epoch == -1 ) break;
	      stage[ epoch ] = ungapped[idx++];
	    }
	  // reset
	  edf.timeline.first_epoch();
	}
    }

  
  double num_n1 = std::numeric_limits<double>::quiet_NaN();
  double num_n2 = std::numeric_limits<double>::quiet_NaN();
  double num_n3 = std::numeric_limits<double>::quiet_NaN();
  double num_rem = std::numeric_limits<double>::quiet_NaN();

  if ( stage.size() != 0 )
    {            
      
      num_rem = 0; 
      num_n1 = 0;
      num_n2 = 0;
      num_n3 = 0;
      
      int num_epoch = 0;
      
      int idx = 0;
      
      for (const auto& [e, stg] : stage )
	{
	  if ( stg == NREM1 || stg == NREM2 || stg == NREM3 || stg == NREM4 || stg == REM )
	    {	      
	      const int e1 = edf.timeline.display_epoch( e ); // to matches epochs[] in calc()
	      
	      ++num_epoch;
	      if ( stg == REM )
		s_rem.insert( e1 );
	      else if ( stg == NREM3 || stg == NREM4 )
		s_n3.insert( e1 );
	      else if ( stg == NREM1 )
		s_n1.insert( e1 ) ;
	      else if ( stg == NREM2 )
		s_n2.insert( e1 );
	    }
	  
	  ++idx;
	}
      
      num_n1 = s_n1.size();
      num_n2 = s_n2.size();
      num_n3 = s_n3.size();
      num_rem  = s_rem.size();
      
      int num_s_all = num_n1 + num_n2 + num_n3 + num_rem;

      // get up to K of each stage
      s_rem = select( s_rem , sel_num_epochs );
      s_n1 = select( s_n1 , sel_num_epochs );
      s_n2 = select( s_n2 , sel_num_epochs );
      s_n3 = select( s_n3, sel_num_epochs );

      num_n1 = s_n1.size();
      num_n2 = s_n2.size();
      num_n3 = s_n3.size();
      num_rem = s_rem.size();

      int num_s_sel = num_n1 + num_n2 + num_n3 + num_rem;
      
      selected_epochs.clear();
      for (const auto& e : s_n1 ) selected_epochs.insert( e );
      for (const auto& e : s_n2 ) selected_epochs.insert( e );
      for (const auto& e : s_n3 ) selected_epochs.insert( e );
      for (const auto& e : s_rem ) selected_epochs.insert( e );
            
      // summary 
      logger << "  detected " << num_epoch << " epochs, w/ " << num_s_sel << " of " << num_s_all << " sleep epochs selected\n"
	     << "  (N1 = " << num_n1 << ", N2 = " << num_n2 << ", N3 = " << num_n3 << ", REM = " << num_rem << " selected)\n";
      
    }
  
    
    
  //
  // iterate over channels
  //
  
  logger << "  processing:";
    
  for (int s=0; s<ns; s++)
    {

      writer.level( signals.label(s) , globals::signal_strat );

      logger << " " << signals.label(s);
      if ( s % 10 == 9 ) logger << "\n      ";

      // build primary feature vector
      
      ctypes_ftrs_t ftr;
      
      ftr.num_rem = num_rem;
      ftr.num_n1 = num_n1;
      ftr.num_n2 = num_n2;
      ftr.num_n3 = num_n3;
      ftr.num_any = std::numeric_limits<double>::quiet_NaN();
           
      //
      // original signal
      //

      interval_t interval = edf.timeline.wholetrace();
      slice_t slice( edf , signals(s) , interval );	  
      std::vector<double> * d = slice.nonconst_pdata();

      //
      // from original signal - physical scales
      //
      // take 'as is', w/ exception that if V will always express as uV

      const std::string pdim = edf.header.phys_dimension[ signals(s) ];
      
      bool is_mV = Helper::imatch( pdim , "mV" );
      bool is_V  = Helper::imatch( pdim , "V" );  

      double fac = 1.0;
      if ( is_mV ) fac = 1000;
      else if ( is_V ) fac = 1000000;
      
      double q1 = std::numeric_limits<double>::quiet_NaN();
      double q3 = std::numeric_limits<double>::quiet_NaN();      
      {
	std::vector<double> w = *d;
	ftr.phys_median = fabs( fac * median_inplace( w ) ); // |median| for polarity invariance
      }
      {
	std::vector<double> w = *d;
	q1 = fac * percentile_inplace(w, 0.25);
      }
      {
	std::vector<double> w = *d;
	q3 = fac * percentile_inplace(w, 0.75);
      }      
      ftr.phys_iqr = fac * ( q3 - q1 ) ;
      

      if ( output )
	{
	  writer.value( "phys_median" , ftr.phys_median );
	  writer.value( "phys_iqr" , ftr.phys_iqr );
	  writer.value( "num_any" , ftr.num_any );
	  writer.value( "num_n1" , ftr.num_n1 );
	  writer.value( "num_n2" , ftr.num_n2 );
	  writer.value( "num_n3" , ftr.num_n3 );	  
	  writer.value( "num_rem" , ftr.num_rem ); 
	}
      
      //
      // 1 Hz signal: get whole trace
      //

      // calculate features
      
      calc_1Hz_stats( *d , Fs_orig[s] , epochs1, epoch2samples1, n1, &ftr , epoch_output );
            
      // output features
      if ( output )
	{
	  writer.level( 1 , "F" );
	  for (const auto& [key, ptr] : feature_map )
	    {
	      writer.level( "RAW", "TRANS" );
	      writer.value( key , ftr.x1.*ptr );
	      
	      writer.level( "RAWLWR", "TRANS" );
	      writer.value( key , ftr.x1_p10.*ptr );
	      
	      writer.level( "RAWUPR", "TRANS" );
	      writer.value( key , ftr.x1_p90.*ptr );
	      
	      writer.level( "DIFF", "TRANS" );
	      writer.value( key , ftr.x1_diff.*ptr );
	      
	      writer.level( "DIFFLWR", "TRANS" );
	      writer.value( key , ftr.x1_diff_p10.*ptr );
	      
	      writer.level( "DIFFUPR", "TRANS" );
	      writer.value( key , ftr.x1_diff_p90.*ptr );
	      
	    }
	  writer.unlevel( "TRANS" );
	  writer.unlevel( "F" );
	}
           

      //
      // 128 Hz signal features
      //

      // calculate & aggregae epoch-wise features
      
      calc_128Hz_stats( *d , Fs_orig[s] , selected_epochs , epochs128, epoch2samples128, n128, &ftr , epoch_output );
      
      
      // output features
      if ( output )
	{
	  writer.level( 128 , "F" );
	  for (const auto& [key, ptr] : feature_map )
	    {
	      writer.level( "RAW", "TRANS" );
	      writer.value( key , ftr.x128.*ptr );
	      
	      writer.level( "DIFF", "TRANS" );
	      writer.value( key , ftr.x128_diff.*ptr );
	      
	      writer.level( "RAWLWR", "TRANS" );
	      writer.value( key , ftr.x128_p10.*ptr );
	      
	      writer.level( "DIFFLWR", "TRANS" );
	      writer.value( key , ftr.x128_diff_p10.*ptr );
	      
	      writer.level( "RAWUPR", "TRANS" );
	      writer.value( key , ftr.x128_p90.*ptr );
	      
	      writer.level( "DIFFUPR", "TRANS" );
	      writer.value( key , ftr.x128_diff_p90.*ptr );
	      
	    }
	  writer.unlevel( "TRANS" );
	  writer.unlevel( "F" );
	}

      
      //
      // predict
      //

      if ( make_predictions )
	{
	   
	  ctypes_pred_t prediction = predict( ftr );
      
	  for (int i=0; i<clslist.size(); i++)
	    {
	      writer.level( clslist[i] , "CTYPE" );
	      writer.value( "PP" , prediction.posteriors[i] );
	    }
	  writer.unlevel( "CTYPE" );
	}
      
      // next signal
      continue;
    }
  
  writer.unlevel( globals::signal_strat );
  
}


void ctypes_t::attach_model( const std::string & model_path ,  const std::string & model_lib )
{

  const std::string s = Helper::expand( model_path ) + "/" + model_lib;
  
  // already loaded?
  if ( lgbm_model_loaded == s ) 
    return;

  lgbm_model_loaded = s;
  	      
  const std::string f_mod = s + ".mod" ;
  const std::string f_ftr = s + ".ftr" ;
  const std::string f_cls = s + ".cls" ;

  if ( ! Helper::fileExists( f_mod ) ) Helper::halt( f_mod + " can not be opened" );
  if ( ! Helper::fileExists( f_ftr ) ) Helper::halt( f_ftr + " can not be opened" );
  if ( ! Helper::fileExists( f_cls ) ) Helper::halt( f_cls + " can not be opened" );

  // model
  
  lgbm.load_model( f_mod );

  // ftrs
  
  std::ifstream I1( f_ftr, std::ios::out );  
  varlist.clear();
  int n;  
  I1 >> n;
  for (int i=0;i<n;i++)
    {
      std::string x;
      I1 >> x;
      varlist.push_back( x );
    }
  I1.close();

  // cls
  std::ifstream I2( f_cls, std::ios::out );  
  clslist.clear();
  I2 >> n;
  for (int i=0;i<n;i++)
    {
      std::string x;
      I2 >> x;
      clslist.push_back( x );
    }
  I2.close();

  // parse ftrs
  var_root.clear();
  var_hz.clear();
  var_trans.clear();
  for (int i=0; i<varlist.size(); i++)
    {
      // expected either 3 or 4 _ delim tokens
      std::vector<std::string> tok = Helper::parse( varlist[i] , "_" );
      if ( tok.size() == 3 )
	{
	  var_root.push_back( tok[0] );
	  var_hz.push_back( tok[1] );
	  var_trans.push_back( tok[2] );
	}
      else if ( tok.size() == 4 )
	{
	  var_root.push_back( tok[0] + "_" + tok[1] ); // some roots have a single _
	  var_hz.push_back( tok[2] );
	  var_trans.push_back( tok[3] );
	}
      else
	Helper::halt( "bad format for " + f_ftr + " line " + varlist[i] );      
    }
  
}


// ---------- helpers ----------

bool ctypes_t::ends_with(const std::string& s, const char* suf) {
  const size_t n = s.size(), m = std::char_traits<char>::length(suf);
  return n >= m && s.compare(n - m, m, suf) == 0;
}

bool ctypes_t::contains(const std::string& s, const char* tok) {
  return s.find(tok) != std::string::npos;
}

const ctypes_specific_ftrs_t& ctypes_t::select_block(
    const ctypes_ftrs_t& f,
    const std::string& hz,        // "1" or "128"
    const std::string& trans      // e.g. "RAW", "RAWLWR", "DIFFLWR", "DIFFUPR"
						     ) {
  const bool is128 = (hz == "128");
  const bool isdiff = contains(trans, "DIFF");
  const bool is_lwr = ends_with(trans, "LWR");
  const bool is_upr = ends_with(trans, "UPR");
  
  // RAW vs DIFF
  if (!is128 && !isdiff) {
    if (is_lwr) return f.x1_p10;
    if (is_upr) return f.x1_p90;
    return f.x1;
  }
  if (!is128 && isdiff) {
    if (is_lwr) return f.x1_diff_p10;
    if (is_upr) return f.x1_diff_p90;
    return f.x1_diff;
  }
  if (is128 && !isdiff) {
    if (is_lwr) return f.x128_p10;
    if (is_upr) return f.x128_p90;
    return f.x128;
  }
  // is128 && isdiff
  if (is_lwr) return f.x128_diff_p10;
  if (is_upr) return f.x128_diff_p90;
  return f.x128_diff;
}

double ctypes_t::get_feature_value(
				   const ctypes_ftrs_t& ftrs,
				   const std::string& root,   // e.g. "ZCR" or "spectral_entropy" etc.
				   const std::string& hz,     // "1" or "128"
				   const std::string& trans   // "RAW", "DIFFLWR", ...
				   )
{
  auto it = feature_map.find(root);

  if (it == feature_map.end())
    Helper::halt( "Unknown root feature: " + root );

  const ctypes_specific_ftrs_t& blk = select_block(ftrs, hz, trans);

  return blk.*(it->second);
}


ctypes_pred_t ctypes_t::predict( const ctypes_ftrs_t & ftr ) {

  ctypes_pred_t p;

  const int m = varlist.size();
    
  // build single row ftr vector
  Eigen::VectorXd X = Eigen::VectorXd::Zero( m );

  for (int i=0; i<m; i++)
    {
      const std::string & root = var_root[i];
      const std::string & hz = var_hz[i];
      const std::string & trans = var_trans[i];
      
      // note: currently missing n_epochs, num_rem, num_n3 etc which do not have hz/trans attr.
      X[i] = get_feature_value(ftr, root, hz, trans );
      
    }
  
  // predict
  p.posteriors = lgbm.predict1( X );

  return p;  
  
}



void ctypes_t::calc_1Hz_stats( const std::vector<double> & x , const double Fs_orig,
			       const std::vector<int> & epochs , 				 
			       const std::vector<std::pair<int,int> > & epoch2samples,
			       const size_t n1,
			       ctypes_ftrs_t * ftr ,
			       bool epoch_output )
{

  const int Fs = 1;
  
  // make new channels  
  std::vector<double> x1, a1;    
  make_1s( x , Fs_orig , &x1, &a1 );
  // first-diff
  std::vector<double> d1 = make_diff( x1 );
  
  //
  // normalize
  //

  const double win = 0.01;

  normalize( &x1 );
  MiscMath::winsorize( &x1, win );

  normalize( &d1 );
  MiscMath::winsorize( &d1, win );

  normalize( &a1 );
  MiscMath::winsorize( &a1, win );

  //
  // compute features per epoch (e.g. 300-seconds, shifted)
  //

  std::vector<ctypes_ftrs_t> aggr;

  if ( epoch_output )
    writer.level( 1 , "F" );
  
  for (int e=0; e<epoch2samples.size(); e++)
    {

      const int s1 = epoch2samples[e].first;
      const int s2 = epoch2samples[e].second;
      
      std::vector<double> ex1( s2 - s1 );
      std::vector<double> ed1( s2 - s1 );
      std::vector<double> ea1( s2 - s1 );      
      
      for (int p=s1; p<s2; p++)
	{
	  const int idx = p - s1;   // 0..len-1
	  //	  std::cout << " p = " << p << " " << s1 << " " << s2 << " " << idx << "  "<< ex1.size() << " " << x1.size() << "\n";
	  ex1[idx] = x1[p];
	  ed1[idx] = d1[p];
	  ea1[idx] = a1[p];	  
	}
      
      ctypes_ftrs_t f;
 
      f.x1      = calc_specific_stats( ex1 , 1 );
      f.x1_diff = calc_specific_stats( ed1 , 1 );
      f.a1      = calc_specific_stats( ea1 , 1 );
 
      // store
      aggr.push_back( f );
      
      // output epoch-level stats?
      if ( epoch_output )
	{	  
	  writer.epoch( epochs[e] );
	  for (const auto& [key, ptr] : feature_map )
	    {
	      writer.level( "RAW", "TRANS" );
	      writer.value( key , f.x1.*ptr );

	      writer.level( "ENV", "TRANS" );
	      writer.value( key , f.a1.*ptr );
	      
	      writer.level( "DIFF", "TRANS" );
	      writer.value( key , f.x1_diff.*ptr );
	    }
	  writer.unlevel( "TRANS" );
	}
    }

  if ( epoch_output )
    {
      writer.unepoch();
      writer.unlevel( "F" );
    }
  
  //
  // aggregate over epochs -> median + 10/90th percentiles
  //

  aggregate1( aggr , ftr );
  
}


void ctypes_t::calc_128Hz_stats( const std::vector<double> & x , const double Fs_orig,
				 const std::set<int> & selected_epochs ,
				 const std::vector<int> & epochs , 				 
				 const std::vector<std::pair<int,int> > & epoch2samples,
				 const size_t n128, 
				 ctypes_ftrs_t * ftr ,
				 bool epoch_output )
{

  
  // *** 128 Hz stats are done on *epoch-wise* - only on N2, N3 and REM epochs
  
  // resample to 128 Hz

  const int Fs = 128;
  
  std::vector<double> x128;
  if ( Fs_orig == Fs )
    {
      x128 = x;
      x128.resize( n128 );
    }
  else
    {
      x128 = dsptools::resample( &x , Fs_orig , Fs , SRC_SINC_FASTEST );
      x128.resize( n128 );
    }

  std::vector<double> d128 = make_diff( x128 ); 

  //
  // normalize
  //

  const double win = 0.01;

  normalize( &x128 );
  MiscMath::winsorize( &x128, win );

  normalize( &d128 );
  MiscMath::winsorize( &d128, win );  
 
  
  //
  // compute features epoch-wise
  //
    
  std::vector<ctypes_ftrs_t> aggr;

  if ( epoch_output )
    writer.level( 128 , "F" );
  
  for (int e=0; e<epoch2samples.size(); e++)
    {

      // process this epoch?
      if ( selected_epochs.find( epochs[e] ) == selected_epochs.end() )
	continue;
      
      const int s1 = epoch2samples[e].first;

      const int s2 = epoch2samples[e].second;
      
      std::vector<double> ex128( s2 - s1 );
      std::vector<double> ed128( s2 - s1 );
      
      for (int p=s1; p<s2; p++)
	{
	  const int idx = p - s1;   // 0..len-1
	  ex128[idx] = x128[p];
	  ed128[idx] = d128[p];	  
	}
      
      ctypes_ftrs_t f;
 
      f.x128      = calc_specific_stats( ex128 , 128 );
      f.x128_diff = calc_specific_stats( ed128 , 128 );

      // store
      aggr.push_back( f );
      
      // output epoch-level stats?
      if ( epoch_output )
	{	  
	  writer.epoch( epochs[e] );
	  for (const auto& [key, ptr] : feature_map )
	    {
	      writer.level( "RAW", "TRANS" );
	      writer.value( key , f.x128.*ptr );

	      writer.level( "DIFF", "TRANS" );
	      writer.value( key , f.x128_diff.*ptr );
	    }
	  writer.unlevel( "TRANS" );
	}
    }
  
  if ( epoch_output )
    {
      writer.unepoch();
      writer.unlevel( "F" );
    }
  
  //
  // aggregate over epochs -> median + 10/90th percentiles
  //

  aggregate128( aggr , ftr );

}




bool ctypes_t::is_nan(double x) { return std::isnan(x); }

// Nearest-rank style using index floor(p*(n-1)) after sorting by nth_element.
// Mutates v order.
double ctypes_t::percentile_inplace(std::vector<double> & v, double p)
{
  if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
  if (p <= 0.0) p = 0.0;
  if (p >= 1.0) p = 1.0;
  
  const size_t n = v.size();
  const size_t k = static_cast<size_t>(std::floor(p * (n - 1)));

  std::nth_element(v.begin(), v.begin() + k, v.end());
  return v[k];
}



void ctypes_t::compute_p10_med_p90(std::vector<double> & vals,
				   double & out_p10,
				   double & out_med,
				   double & out_p90)
{
  if (vals.empty())
    {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    out_p10 = nan;
    out_med = nan;
    out_p90 = nan;
    return;
  }
  
  // Each call re-partitions; that's fine.
  out_p10 = percentile_inplace(vals, 0.10);
  out_med = percentile_inplace(vals, 0.50);
  out_p90 = percentile_inplace(vals, 0.90);
}



void ctypes_t::aggregate1(const std::vector<ctypes_ftrs_t> & aggr,
                         ctypes_ftrs_t * ftr)
{
  if (!ftr) return;

  // Zero outputs
  ftr->x1      = {};
  ftr->a1      = {};
  ftr->x1_diff = {};

  ftr->x1_p10      = {};
  ftr->a1_p10      = {};
  ftr->x1_diff_p10 = {};

  ftr->x1_p90      = {};
  ftr->a1_p90      = {};
  ftr->x1_diff_p90 = {};

  if (aggr.empty()) return;

  auto aggregate_specific =
    [&](ctypes_specific_ftrs_t & med,
        ctypes_specific_ftrs_t & p10,
        ctypes_specific_ftrs_t & p90,
        const auto & pick)
  {
    std::vector<double> vals;
    vals.reserve(aggr.size());

    auto do_field = [&](auto field_ptr, double & out_p10, double & out_med, double & out_p90)
    {
      vals.clear();
      for (const auto & a : aggr)
      {
        const double x = (pick(a).*field_ptr);
        if (!is_nan(x)) vals.push_back(x);
      }
      compute_p10_med_p90(vals, out_p10, out_med, out_p90);
    };

    do_field(&ctypes_specific_ftrs_t::flatline_frac,          p10.flatline_frac,          med.flatline_frac,          p90.flatline_frac);
    do_field(&ctypes_specific_ftrs_t::clip_frac,              p10.clip_frac,              med.clip_frac,              p90.clip_frac);
    do_field(&ctypes_specific_ftrs_t::unique_frac,            p10.unique_frac,            med.unique_frac,            p90.unique_frac);
    do_field(&ctypes_specific_ftrs_t::most_common_value_frac, p10.most_common_value_frac, med.most_common_value_frac, p90.most_common_value_frac);
    do_field(&ctypes_specific_ftrs_t::transition_rate,        p10.transition_rate,        med.transition_rate,        p90.transition_rate);

    do_field(&ctypes_specific_ftrs_t::logH1,                  p10.logH1,                  med.logH1,                  p90.logH1);
    do_field(&ctypes_specific_ftrs_t::H2,                     p10.H2,                     med.H2,                     p90.H2);
    do_field(&ctypes_specific_ftrs_t::H3,                     p10.H3,                     med.H3,                     p90.H3);
    do_field(&ctypes_specific_ftrs_t::line_length,            p10.line_length,            med.line_length,            p90.line_length);
    do_field(&ctypes_specific_ftrs_t::kurtosis,               p10.kurtosis,               med.kurtosis,               p90.kurtosis);
    do_field(&ctypes_specific_ftrs_t::skewness,               p10.skewness,               med.skewness,               p90.skewness);
    do_field(&ctypes_specific_ftrs_t::ZCR,                    p10.ZCR,                    med.ZCR,                    p90.ZCR);

    do_field(&ctypes_specific_ftrs_t::spectral_centroid,      p10.spectral_centroid,      med.spectral_centroid,      p90.spectral_centroid);
    do_field(&ctypes_specific_ftrs_t::spectral_edge,          p10.spectral_edge,          med.spectral_edge,          p90.spectral_edge);
    do_field(&ctypes_specific_ftrs_t::spectral_bandwidth,     p10.spectral_bandwidth,     med.spectral_bandwidth,     p90.spectral_bandwidth);
    do_field(&ctypes_specific_ftrs_t::spectral_flatness,      p10.spectral_flatness,      med.spectral_flatness,      p90.spectral_flatness);
    do_field(&ctypes_specific_ftrs_t::spectral_entropy,       p10.spectral_entropy,       med.spectral_entropy,       p90.spectral_entropy);
    do_field(&ctypes_specific_ftrs_t::spectral_lowpower,      p10.spectral_lowpower,      med.spectral_lowpower,      p90.spectral_lowpower);
    do_field(&ctypes_specific_ftrs_t::spectral_highpower,     p10.spectral_highpower,     med.spectral_highpower,     p90.spectral_highpower);

    do_field(&ctypes_specific_ftrs_t::acf1,                   p10.acf1,                   med.acf1,                   p90.acf1);
    do_field(&ctypes_specific_ftrs_t::acf_1s,                 p10.acf_1s,                 med.acf_1s,                 p90.acf_1s);
    do_field(&ctypes_specific_ftrs_t::acf_decay,              p10.acf_decay,              med.acf_decay,              p90.acf_decay);
    do_field(&ctypes_specific_ftrs_t::acf_peak,               p10.acf_peak,               med.acf_peak,               p90.acf_peak);
    do_field(&ctypes_specific_ftrs_t::acf_min,                p10.acf_min,                med.acf_min,                p90.acf_min);
  };

  // across epochs: epoch-level values live in a.x1 / a.x1_diff 

  aggregate_specific(
    ftr->x1, ftr->x1_p10, ftr->x1_p90,
    [](const ctypes_ftrs_t & a) -> const ctypes_specific_ftrs_t & { return a.x1; }
  );

  aggregate_specific(
    ftr->a1, ftr->a1_p10, ftr->a1_p90,
    [](const ctypes_ftrs_t & a) -> const ctypes_specific_ftrs_t & { return a.a1; }
  );
  
  aggregate_specific(
    ftr->x1_diff, ftr->x1_diff_p10, ftr->x1_diff_p90,
    [](const ctypes_ftrs_t & a) -> const ctypes_specific_ftrs_t & { return a.x1_diff; }
  );

}


void ctypes_t::aggregate128(const std::vector<ctypes_ftrs_t> & aggr,
                         ctypes_ftrs_t * ftr)
{
  if (!ftr) return;

  // Zero outputs
  ftr->x128      = {};
  ftr->x128_diff = {};

  ftr->x128_p10      = {};
  ftr->x128_diff_p10 = {};

  ftr->x128_p90      = {};
  ftr->x128_diff_p90 = {};

  if (aggr.empty()) return;

  auto aggregate_specific =
    [&](ctypes_specific_ftrs_t & med,
        ctypes_specific_ftrs_t & p10,
        ctypes_specific_ftrs_t & p90,
        const auto & pick)
  {
    std::vector<double> vals;
    vals.reserve(aggr.size());

    auto do_field = [&](auto field_ptr, double & out_p10, double & out_med, double & out_p90)
    {
      vals.clear();
      for (const auto & a : aggr)
      {
        const double x = (pick(a).*field_ptr);
        if (!is_nan(x)) vals.push_back(x);
      }
      compute_p10_med_p90(vals, out_p10, out_med, out_p90);
    };

    do_field(&ctypes_specific_ftrs_t::flatline_frac,          p10.flatline_frac,          med.flatline_frac,          p90.flatline_frac);
    do_field(&ctypes_specific_ftrs_t::clip_frac,              p10.clip_frac,              med.clip_frac,              p90.clip_frac);
    do_field(&ctypes_specific_ftrs_t::unique_frac,            p10.unique_frac,            med.unique_frac,            p90.unique_frac);
    do_field(&ctypes_specific_ftrs_t::most_common_value_frac, p10.most_common_value_frac, med.most_common_value_frac, p90.most_common_value_frac);
    do_field(&ctypes_specific_ftrs_t::transition_rate,        p10.transition_rate,        med.transition_rate,        p90.transition_rate);

    do_field(&ctypes_specific_ftrs_t::logH1,                  p10.logH1,                  med.logH1,                  p90.logH1);
    do_field(&ctypes_specific_ftrs_t::H2,                     p10.H2,                     med.H2,                     p90.H2);
    do_field(&ctypes_specific_ftrs_t::H3,                     p10.H3,                     med.H3,                     p90.H3);
    do_field(&ctypes_specific_ftrs_t::line_length,            p10.line_length,            med.line_length,            p90.line_length);
    do_field(&ctypes_specific_ftrs_t::kurtosis,               p10.kurtosis,               med.kurtosis,               p90.kurtosis);
    do_field(&ctypes_specific_ftrs_t::skewness,               p10.skewness,               med.skewness,               p90.skewness);
    do_field(&ctypes_specific_ftrs_t::ZCR,                    p10.ZCR,                    med.ZCR,                    p90.ZCR);

    do_field(&ctypes_specific_ftrs_t::spectral_centroid,      p10.spectral_centroid,      med.spectral_centroid,      p90.spectral_centroid);
    do_field(&ctypes_specific_ftrs_t::spectral_edge,          p10.spectral_edge,          med.spectral_edge,          p90.spectral_edge);
    do_field(&ctypes_specific_ftrs_t::spectral_bandwidth,     p10.spectral_bandwidth,     med.spectral_bandwidth,     p90.spectral_bandwidth);
    do_field(&ctypes_specific_ftrs_t::spectral_flatness,      p10.spectral_flatness,      med.spectral_flatness,      p90.spectral_flatness);
    do_field(&ctypes_specific_ftrs_t::spectral_entropy,       p10.spectral_entropy,       med.spectral_entropy,       p90.spectral_entropy);
    do_field(&ctypes_specific_ftrs_t::spectral_lowpower,      p10.spectral_lowpower,      med.spectral_lowpower,      p90.spectral_lowpower);
    do_field(&ctypes_specific_ftrs_t::spectral_highpower,     p10.spectral_highpower,     med.spectral_highpower,     p90.spectral_highpower);

    do_field(&ctypes_specific_ftrs_t::acf1,                   p10.acf1,                   med.acf1,                   p90.acf1);
    do_field(&ctypes_specific_ftrs_t::acf_1s,                 p10.acf_1s,                 med.acf_1s,                 p90.acf_1s);
    do_field(&ctypes_specific_ftrs_t::acf_decay,              p10.acf_decay,              med.acf_decay,              p90.acf_decay);
    do_field(&ctypes_specific_ftrs_t::acf_peak,               p10.acf_peak,               med.acf_peak,               p90.acf_peak);
    do_field(&ctypes_specific_ftrs_t::acf_min,                p10.acf_min,                med.acf_min,                p90.acf_min);
  };

  // across epochs: epoch-level values live in a.x128 / a.x128_diff 
  aggregate_specific(
    ftr->x128, ftr->x128_p10, ftr->x128_p90,
    [](const ctypes_ftrs_t & a) -> const ctypes_specific_ftrs_t & { return a.x128; }
  );

  aggregate_specific(
    ftr->x128_diff, ftr->x128_diff_p10, ftr->x128_diff_p90,
    [](const ctypes_ftrs_t & a) -> const ctypes_specific_ftrs_t & { return a.x128_diff; }
  );

}




ctypes_specific_ftrs_t ctypes_t::calc_specific_stats( const std::vector<double> & x , int Fs )
{

  // Fs is now either 1 or 128 Hz
  
  const bool slow_sr = Fs < 32 ; 
  
  ctypes_specific_ftrs_t out;

  //
  // 0) basic
  //

  out.flatline_frac = flatline_fraction( x );

  out.clip_frac = clip_fraction( x );

  out.unique_frac = unique_frac_q( x );

  out.most_common_value_frac = most_common_value_frac_q( x );

  out.transition_rate = transition_rate_q( x );
  
  //
  // 1) time-domain stats
  //

  out.line_length = line_length_norm_sd( x );

  const double m = mad(x);
  const double eps = (m > 0.0) ? 1e-9 * m : 0.0;
  out.ZCR = zero_cross_rate( x , eps );
  
  // hjorth
  MiscMath::hjorth( &x , &out.logH1 , &out.H2, &out.H3 );
  constexpr double EPS = 1e-12;
  out.logH1 = (out.logH1 > 0.0) ? std::log(out.logH1) : std::log(EPS);

  double mn = mean( x ) ;
  double sd = sd_sample( x ) ;
  
  out.skewness = fabs( MiscMath::skewness( x , mn , sd ) ); // |skew|
  out.kurtosis = MiscMath::kurtosis( x , mn );

  
  //
  // 2) spectral summaries (fixed 
  //
  
  spectral_stats_t spectral_stats( x , Fs );

  out.spectral_centroid = spectral_stats.centroid;
  out.spectral_edge = spectral_stats.edge90;
  out.spectral_bandwidth = spectral_stats.bandwidth;
  out.spectral_flatness = spectral_stats.flatness;
  out.spectral_entropy = spectral_stats.entropy;
  out.spectral_highpower = spectral_stats.high;
  out.spectral_lowpower = spectral_stats.low;
 
  
  //
  // 3) ACF features
  //

  if ( Fs <= 0 || x.size() < 4) {
    out.acf1 = out.acf_1s = out.acf_decay = out.acf_peak = out.acf_min = std::numeric_limits<double>::quiet_NaN();
  } else {
    
    // lags we care about - in both cases 640 samples, but either 5sec (5*128Hz=640) or ~10 mins (640*1Hz)
    const int lag_1s = Fs;
    const int max_sec = slow_sr ? 640 : 5 ; 
    const int lag_max = std::max(lag_1s, max_sec * Fs);
    
    // compute ACF
    acf_t A(x, lag_max);
    const std::vector<double> r = A.acf(); // r[0..lag_max], normalized by r[0]
    
    if (r.size() < 2) {
      out.acf1 = out.acf_1s = out.acf_decay = out.acf_peak = out.acf_min = std::numeric_limits<double>::quiet_NaN();
    } else {
      
      // 1) acf(1)
      out.acf1 = r.size() > 1 ? r[1] : 0.0;
      
      // 2) acf at 1 second
      if ( slow_sr )
	out.acf_1s = std::numeric_limits<double>::quiet_NaN();
      else
	out.acf_1s = (lag_1s < (int)r.size()) ? r[lag_1s] : r.back();
      
      // 3) decay time: first lag where ACF < 1/e
      const double thr = std::exp(-1.0); // ~0.3679
      int decay_lag = -1;
      for (int lag = 1; lag < (int)r.size(); ++lag) {
	if (r[lag] < thr) { decay_lag = lag; break; }
      }

      // If it never drops below threshold within lag_max, set to lag_max/Fs (capped)
      if (decay_lag < 0) decay_lag = (int)r.size() - 1;
      out.acf_decay = (double)decay_lag / (double)Fs;
      
      // 4) peak ACF in a physiologic lag range (exclude tiny lags)
      // For PSG channel ID: 0.5–5 s is a good default:
      //  - captures resp rhythm (2–6 s)
      //  - captures heart rhythm (~0.6–1.2 s) depending on HR
      const int lag_lo = std::max(1, (int)std::lround(0.5 * Fs));
      const int lag_hi = std::min((int)r.size() - 1, max_sec * Fs);

      double pk = -1.0;
      if (lag_lo <= lag_hi) {
	for (int lag = lag_lo; lag <= lag_hi; ++lag) {
	  if (r[lag] > pk) pk = r[lag];
	}
      }
      out.acf_peak = (pk >= 0.0 ? pk : 0.0);
      
      // 5) min ACF for short lags (artifact/ringing)
      if ( slow_sr )
	out.acf_min = std::numeric_limits<double>::quiet_NaN();
      else
	{
	  double min_val = std::numeric_limits<double>::infinity();
	  
	  const int lag2_lo = std::max(1, (int)std::round(0.1 * Fs));
	  const int lag2_hi = std::min((int)x.size() - 1, (int)std::round(1.0 * Fs));
	  
	  if (lag2_lo < lag2_hi)
	    {
	      for (int lag = lag2_lo; lag <= lag2_hi; ++lag) {
		if (std::isfinite(r[lag]) && r[lag] < min_val)
		  min_val = r[lag];
	      }
	      
	      // If something went wrong, fall back safely
	      if (!std::isfinite(min_val))
		out.acf_min = std::numeric_limits<double>::quiet_NaN();
	      else
		out.acf_min = min_val;
	    }
	}
            
    }
    
  }

  return out;
}


// Robust normalization: x := (x - median) / (1.4826 * MAD)
void ctypes_t::normalize( std::vector<double> * x )
{

  // Ignores non-finite values (NaN/Inf) when estimating center/scale
  // Leaves non-finite values unchanged
  // If MAD ~ 0, falls back to IQR-based scale; if still ~0, mean-abs-dev; else no-op/zero-center.

  if (x == nullptr) return;
  if (x->empty()) return;
  
  // Collect finite values for robust stats
  std::vector<double> v;
  v.reserve(x->size());
  for (double a : *x) {
    if (std::isfinite(a)) v.push_back(a);
  }
  if (v.empty()) return;
  
  // Median (modifies v)
  const double med = median_inplace(v);
  
  // MAD: median(|xi - med|)
  std::vector<double> dev;
  dev.reserve(v.size());
  for (double a : v) dev.push_back(std::fabs(a - med));
  double mad = median_inplace(dev);
  
  auto safe_pos = [](double s) {
    return std::isfinite(s) && s > 0.0;
  };
  
  double scale = std::numeric_limits<double>::quiet_NaN();
  
  // Primary robust scale (Gaussian-consistent)
  if (safe_pos(mad)) {
    scale = 1.4826 * mad;
  } else {
    // Fallback 1: IQR / 1.349 (Gaussian-consistent)
    // Compute Q1/Q3 via nth_element on a fresh copy
    std::vector<double> q = v;
    const size_t n = q.size();
    
    auto quantile_inplace = [&](std::vector<double>& w, double p) -> double {
      if (w.empty()) return std::numeric_limits<double>::quiet_NaN();
      if (p <= 0.0) return *std::min_element(w.begin(), w.end());
      if (p >= 1.0) return *std::max_element(w.begin(), w.end());
      
      const double idx = p * (w.size() - 1);
      const size_t k = static_cast<size_t>(idx);
      std::nth_element(w.begin(), w.begin() + k, w.end());
      const double a = w[k];
      
      const size_t k2 = std::min(k + 1, w.size() - 1);
      std::nth_element(w.begin(), w.begin() + k2, w.end());
      const double b = w[k2];
      
      const double frac = idx - static_cast<double>(k);
      return a + (b - a) * frac;
    };
    
    double q1 = std::numeric_limits<double>::quiet_NaN();
    double q3 = std::numeric_limits<double>::quiet_NaN();
    {
      std::vector<double> w = q;
      q1 = quantile_inplace(w, 0.25);
    }
    {
      std::vector<double> w = q;
      q3 = quantile_inplace(w, 0.75);
    }
    
    const double iqr = q3 - q1;
    if (safe_pos(iqr)) {
      scale = iqr / 1.349;
    } else {
      // Fallback 2: mean absolute deviation from median (not as robust, but better than 0)
      double sum = 0.0;
      for (double a : v) sum += std::fabs(a - med);
      const double mad_mean = sum / static_cast<double>(v.size());
      if (safe_pos(mad_mean)) scale = mad_mean;
    }
  }

  // If still no usable scale, just center (or no-op if med not finite)
  if (!safe_pos(scale)) {
    if (!std::isfinite(med)) return;
    for (double &a : *x) {
      if (std::isfinite(a)) a -= med;
    }
    return;
  }
  
  // Apply normalization
  const double inv = 1.0 / scale;
  for (double &a : *x) {
    if (std::isfinite(a)) a = (a - med) * inv;
  }
}


std::set<int> ctypes_t::select( std::set<int> & s , const int n )
{
  if (n <= 0 || s.empty()) return {};
  if (n >= static_cast<int>(s.size())) return s;

  // Copy set into vector for indexing/shuffling
  std::vector<int> v(s.begin(), s.end());

  // Random engine (seed once per call; acceptable here)
  static std::mt19937 rng{std::random_device{}()};

  // Shuffle and take first n
  std::shuffle(v.begin(), v.end(), rng);

  return std::set<int>(v.begin(), v.begin() + n);

}

#endif // LGBM

