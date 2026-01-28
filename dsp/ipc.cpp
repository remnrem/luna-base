
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

#include "dsp/ipc.h"

#include <vector>
#include "edf/edf.h"
#include "edf/slice.h"
#include "param.h"

#include "helper/logger.h"
#include "db/db.h"

extern writer_t writer;
extern logger_t logger;

void dsptools::ipc( edf_t & edf , param_t & param )
{

  const std::string hilbert_mag_prefix = "_ht_mag";
  const std::string hilbert_phase_prefix = "_ht_phase";
  
  
  // signals
  std::string s1, s2;
  
  bool all_by_all = true;
  
  // take sig or sig1 / sig2
  if ( param.has( "sig1" ) && param.has( "sig2" ) )
    {
      s1 = param.value( "sig1" );
      s2 = param.value( "sig2" );
      all_by_all = false;
    }
  else
    s1 = s2 = param.value( "sig" );
  
  // check signals
  const bool NO_ANNOTS = true;
  signal_list_t signals1 = edf.header.signal_list( s1 , NO_ANNOTS );
  signal_list_t signals2 = edf.header.signal_list( s2 , NO_ANNOTS );
  
  const int ns1 = signals1.size();
  const int ns2 = signals2.size();

  if ( ns1 == 0 || ns2 == 0 )
    {
      logger << "  no pairwise signal comparisons to perform\n";
      return ;
    }

  int np = all_by_all ? ns1*(ns1-1) : ns1 * ns2;


  //
  // check all signals have hilbert phase & amp ( sig_ht_mag and sig_ht_phase ) 
  //

  int sr = 0;
  for (int s=0; s<ns1; s++)
    {
      std::string s_mag   = signals1.label(s) + hilbert_mag_prefix;
      std::string s_phase = signals1.label(s) + hilbert_phase_prefix;      
      if ( ! edf.header.has_signal( s_mag ) ) Helper::halt( "could not find " + s_mag );
      if ( ! edf.header.has_signal( s_phase ) ) Helper::halt( "could not find " + s_phase );
      int slot = edf.header.signal( s_mag );
      if ( sr == 0 )
	sr = edf.header.sampling_freq( slot );
      else if (  edf.header.sampling_freq( slot ) != sr )
	Helper::halt( "requires uniform sampling rate across signals" );      
    }

  for (int s=0; s<ns2; s++)
    {
      std::string s_mag   = signals2.label(s) + hilbert_mag_prefix;
      std::string s_phase = signals2.label(s) + hilbert_phase_prefix;      
      if ( ! edf.header.has_signal( s_mag ) ) Helper::halt( "could not find " + s_mag );
      if ( ! edf.header.has_signal( s_phase ) ) Helper::halt( "could not find " + s_phase );
      int slot = edf.header.signal( s_mag );
      if ( sr == 0 )
	sr = edf.header.sampling_freq( slot );
      else if (  edf.header.sampling_freq( slot ) != sr )
	Helper::halt( "requires uniform sampling rate across signals" );      
    }
  

  
  //
  // for each seed signal
  //   -- look at each epoch
  //   -- bundle other channels
  //   -- compute IPC stats
  //   -- accumulate channel
  //   -- add channel
  //   -- report mean of per-epoch stats
  

  for (int s=0; s<ns1; s++)
    {
      std::string s_mag   = signals1.label(s) + hilbert_mag_prefix;
      std::string s_phase = signals1.label(s) + hilbert_phase_prefix;

        const int ne = edf.timeline.calc_epochs_contig();
	logger << "  iterating over " << ne << " contig-basad epochs\n";

	
      
      
    }
  
  
  //
  // Pre-compute 
  //
  
  ipc_t ipc;

  
}


// Key IPC computations (generic sampling rate), assuming you already have:
//  - bandpass_filter(x, sr, f_lo, f_hi, ...)
//  - analytic_signal(x_filtered) -> phase(t), amp(t)  [Hilbert-based]
// This code only shows the core IPC math + summaries.
//
// IPC definition (samplewise):
//   dphi(t) = wrap_to_pi( phi_seed(t) - phi_tgt(t) )
//   IPC(t)  = cos(dphi(t))                 // signed instantaneous phase coherence
//
// Optional amplitude weighting + gating:
//   w0(t) = min(amp_seed(t), amp_tgt(t))
//   w(t)  = w0(t) if w0(t) >= thr else 0   // gate low-amp times (phase unreliable)
//   IPCw(t) = w(t) * cos(dphi(t))
//
// Summaries:
//   mean_ipc          = mean(IPC(t)) over used samples
//   mean_ipc_weighted = sum(w(t)*IPC(t))/sum(w(t))
//   plv               = |sum(w(t)*exp(i*dphi(t)))| / sum(w(t))   (weighted PLV)
//   mean_phase        = arg( sum(w(t)*exp(i*dphi(t))) )          (circular mean lag)
//   frac_inphase      = fraction(|dphi(t)| < pi/6) among used samples
//
// Notes:
// - sr is arbitrary. Filtering band edges are in Hz.
// - dooes filtering/Hilbert on-the-fly per channel

#include <vector>
#include <cmath>
#include <complex>
#include <limits>
#include <algorithm>



// ------------------------------ IPC core ------------------------------


ipc_output_t ipc_t::compute_ipc(const ipc_phaseamp_t & seed , 
				const ipc_phaseamp_t & tgt,
				double sr,
				const ipc_param_t & P,
				bool return_timeseries )
{

  ipc_output_t out;

  const size_t N = seed.amp.size();

  if (N == 0 || tgt.amp.size() != N) return out;

  out.summary.n_total = N;
  
  // Build weights w0(t)=min(amp_seed, amp_tgt) or 1.0
  std::vector<double> w0(N, 1.0);
  if (P.amplitude_weighting) {
    w0.resize(N);
    for (size_t t = 0; t < N; ++t) w0[t] = std::min(seed.amp[t], tgt.amp[t]);
  }

  // Determine gating threshold
  double thr = 0.0;
  if (P.amplitude_weighting && P.gate_low_amp) {
    thr = P.gate_use_quantile ? quantile(w0, P.gate_quantile) : P.gate_abs;
    if (!finite(thr)) thr = 0.0;
  }
  
  // Compute per-sample IPC quantities (optionally returned)
  if (return_timeseries) {
    out.ipc.resize(N);
    out.ipcw.resize(N);
    out.dphi.resize(N);
    out.w.resize(N);
  }

  // Summaries (over “used” samples): apply edge drop in seconds
  const size_t edge = (size_t)std::llround(P.edge_drop_sec * sr);
  const size_t lo = std::min(edge, N);
  const size_t hi = (N > edge ? N - edge : 0);
  
  double sum_ipc = 0.0;
  double sum_w = 0.0;
  double sum_w_ipc = 0.0;

  size_t n_used = 0;
  size_t n_inphase = 0;
  
  std::complex<double> sum_e(0.0, 0.0); // for PLV + mean_phase

  for (size_t t = 0; t < N; ++t) {
    // Phase difference and IPC
    double d = wrap_to_pi(seed.phase[t] - tgt.phase[t]);
    double c = std::cos(d);

    // Weight
    double w = 1.0;
    if (P.amplitude_weighting) {
      w = w0[t];
      if (P.gate_low_amp && w < thr) w = 0.0;
    }

    if (return_timeseries) {
      out.dphi[t] = d;
      out.ipc[t]  = c;
      out.w[t]    = w;
      out.ipcw[t] = w * c;
    }
    
    // Summaries: only use interior + finite + (if weighted) w>0
    if (t < lo || t >= hi) continue;
    if (!finite(d) || !finite(c)) continue;
    if (P.amplitude_weighting && w <= 0.0) continue;

    ++n_used;
    sum_ipc += c;

    if (P.amplitude_weighting) {
      sum_w += w;
      sum_w_ipc += w * c;
      sum_e += std::polar(w, d);          // w * exp(i*dphi)
    } else {
      sum_e += std::polar(1.0, d);
    }

    if (std::fabs(d) < M_PI / 6.0) ++n_inphase;
  }

  out.summary.n_used = n_used;
  
  if (n_used > 0) {
    out.summary.mean_ipc = sum_ipc / (double)n_used;
    out.summary.frac_inphase = (double)n_inphase / (double)n_used;
    
    // PLV and mean phase lag from sum_e
    if (P.amplitude_weighting) {
      if (sum_w > 0.0) {
        out.summary.plv = std::abs(sum_e) / sum_w;
        out.summary.mean_ipc_weighted = sum_w_ipc / sum_w;
        out.summary.mean_phase = std::atan2(sum_e.imag(), sum_e.real());
      }
    } else {
      out.summary.plv = std::abs(sum_e) / (double)n_used;
      out.summary.mean_phase = std::atan2(sum_e.imag(), sum_e.real());
    }
  }
  
  return out;
}


// ------------------------------ Batch: seeds vs targets ------------------------------


ipc_batch_result_t ipc_t::compute_ipc_seed_to_set(const std::vector<ipc_phaseamp_t>& signals,
						  const std::vector<int>& s1,
						  const std::vector<int>& s2,
						  double sr,
						  const ipc_param_t & P,
						  ipc_derived_mode_t mode)
{

  ipc_batch_result_t out;
  
  for (int si : s1) {
    const auto& seed = signals[si];
    const size_t N = seed.amp.size();
    
    std::vector<double> acc(N, 0.0);
    std::vector<double> accw(N, 0.0);
    size_t n_tgt_used = 0;
    
    for (int tj : s2) {
      if (tj == si) continue;
      const auto& tgt = signals[tj];
      if (tgt.amp.size() != N) continue;
      
      const bool want_ts = (mode == ipc_derived_mode_t::PER_PAIR_IPC ||
                            mode == ipc_derived_mode_t::PER_PAIR_IPCW ||
                            mode == ipc_derived_mode_t::PER_SEED_MEAN_IPC ||
                            mode == ipc_derived_mode_t::PER_SEED_MEAN_IPCW);
      
      ipc_output_t R = compute_ipc(seed, tgt, sr, P, want_ts);
      
      // store summary row
      ipc_pair_summary_row_t row;
      row.seed_idx = si;
      row.tgt_idx = tj;
      row.summary = R.summary;
      out.summaries.push_back(row);
      
      // derived channels
      if (mode == ipc_derived_mode_t::PER_PAIR_IPC) {
        out.derived.push_back(R.ipc);        
      } else if (mode == ipc_derived_mode_t::PER_PAIR_IPCW) {
        out.derived.push_back(R.ipcw);
      } else if (mode == ipc_derived_mode_t::PER_SEED_MEAN_IPC) {
        for (size_t t = 0; t < N; ++t) acc[t] += R.ipc[t];
        ++n_tgt_used;
      } else if (mode == ipc_derived_mode_t::PER_SEED_MEAN_IPCW) {
        for (size_t t = 0; t < N; ++t) { acc[t] += R.ipcw[t]; accw[t] += R.w[t]; }
        ++n_tgt_used;
      }
    }

    if (mode == ipc_derived_mode_t::PER_SEED_MEAN_IPC && n_tgt_used > 0) {
      std::vector<double> v(N);
      for (size_t t = 0; t < N; ++t) v[t] = acc[t] / (double)n_tgt_used;
      out.derived.push_back(std::move(v));      
    } else if (mode == ipc_derived_mode_t::PER_SEED_MEAN_IPCW && n_tgt_used > 0) {
      std::vector<double> v(N);
      for (size_t t = 0; t < N; ++t) v[t] = (accw[t] > 0.0 ? acc[t] / accw[t] : 0.0); // or NA sentinel
      out.derived.push_back(std::move(v));
    }
  }

  return out;
}




//
// helpers
//

double ipc_t::wrap_to_pi(double x) { return std::atan2(std::sin(x), std::cos(x)); }

bool ipc_t::finite(double x) { return std::isfinite(x); }

double ipc_t::quantile(std::vector<double> v, double q) {
  v.erase(std::remove_if(v.begin(), v.end(), [](double x){return !finite(x);}), v.end());
  if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
  std::sort(v.begin(), v.end());
  double pos = q * (v.size() - 1);
  size_t lo = (size_t)std::floor(pos);
  size_t hi = (size_t)std::ceil(pos);
  double frac = pos - lo;
  return v[lo] * (1.0 - frac) + v[hi] * frac;
}

