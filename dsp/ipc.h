
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

#ifndef __IPC_H__
#define __IPC_H__

#include <vector>
#include <string>
#include <cmath>
#include <complex>
#include <algorithm>
#include <numeric>

struct edf_t;
struct param_t;

namespace dsptools
{

  void ipc( edf_t & , param_t & );

}



enum class ipc_derived_mode_t {
  NO_IPC,
  PER_PAIR_IPC,        // add IPC(t) per (seed,tgt)
  PER_PAIR_IPCW,       // add IPCw(t) per (seed,tgt)
  PER_SEED_MEAN_IPC,   // add mean IPC(t) across targets per seed
  PER_SEED_MEAN_IPCW   // add weighted mean across targets per seed
};

struct ipc_param_t {
  
  // freq-band
  double f_lo = -1.0;
  double f_hi = -1.0;
  
  // Amplitude weighting and gating
  bool amplitude_weighting = true;
  bool gate_low_amp = true;
  
  // Gating threshold mode
  bool gate_use_quantile = true;
  double gate_quantile = 0.30;   // e.g. discard bottom 30% of min-amp
  double gate_abs = 0.0;         // used if gate_use_quantile=false
  
  // Edge drop to avoid filter/Hilbert transients (in seconds)
  double edge_drop_sec = 0.0;    // e.g. 10–30 sec for slow bands
};


struct ipc_stats_t {
  size_t n_total = 0;
  size_t n_used = 0;  
  double mean_ipc = std::numeric_limits<double>::quiet_NaN();
  double mean_ipc_weighted = std::numeric_limits<double>::quiet_NaN();  
  double plv = std::numeric_limits<double>::quiet_NaN();         // weighted
  double mean_phase = std::numeric_limits<double>::quiet_NaN();  // [-pi,pi]
  double frac_inphase = std::numeric_limits<double>::quiet_NaN();// |dphi|<pi/6
};


// batch
struct ipc_pair_summary_row_t { 
  int seed_idx;
  int tgt_idx;
  ipc_stats_t summary;
};


struct ipc_batch_result_t { 
  std::vector<ipc_pair_summary_row_t> summaries;  

  // Optional: derived channels
  // (a) per pair IPC(t) channel (can explode in number)
  // (b) per seed aggregate IPC_mean(t) over all targets (preferred)
  std::vector<std::vector<double>> derived; // one vector per new channel
};

struct ipc_output_t {
  std::vector<double> ipc;    // cos(dphi)
  std::vector<double> ipcw;   // w*cos(dphi) (0 if gated)  <- output
  std::vector<double> dphi;   // wrapped phase diff
  std::vector<double> w;      // weights used (0 if gated) <- output
  ipc_stats_t summary;
};


struct ipc_phaseamp_t {
  std::vector<double> phase;
  std::vector<double> amp;
};


struct ipc_t
{
  
  static ipc_output_t compute_ipc(const ipc_phaseamp_t & seed,
				  const ipc_phaseamp_t & tgt,
				  double sr,
				  const ipc_param_t & param,
				  bool return_timeseries = false) ;
  
  static ipc_batch_result_t compute_ipc_seed_to_set(const std::vector<ipc_phaseamp_t>& signals,
						    const std::vector<int>& s1,
						    const std::vector<int>& s2,
						    double sr,
						    const ipc_param_t & P,
						    ipc_derived_mode_t mode);
  

  // helpers
  static double wrap_to_pi(double x);
  static bool finite(double x);
  static double quantile(std::vector<double> v, double q);

  struct circ_stats_t {
    double mean_phase;  // in (-pi, pi]
    double R;           // mean resultant length [0, 1]
  };
  
  static circ_stats_t circular_mean(const std::vector<double>& theta);
};



#endif
