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

#ifndef __LUNA_DPP_SPEC_H__
#define __LUNA_DPP_SPEC_H__

#include <string>
#include <vector>
#include <set>
#include <map>

// DPP's feature-spec grammar, modeled on pops_spec_t/pops_specs_t
// (pops/spec.h) but generalized for: (a) an explicit window-length
// dimension (crossed per feature block via windows=), (b) two-channel
// connectivity features (channel pairs, comma-separated, following the
// same convention pops_spec_t already uses for its COH feature), each
// channel slot optionally carrying its own named filter band, and (c) a
// simpler, flat feature list -- no SELECT/DROP or level-1/level-2 derived
// features, since DPP stage 2 has no SVD/smoothing/final-column-selection
// concept yet (every specified feature is always computed and output).

enum dpp_feature_t
  {
    DPP_PSD ,       // canonical-band log power (delta, theta, alpha, sigma,
                    // beta) -- natural log, unconditional, matching POPS's
                    // own POPS_BANDS feature convention
    DPP_SLOPE ,     // spectral slope/aperiodic exponent
    DPP_HJORTH ,    // activity, mobility, complexity
    DPP_SKEW ,
    DPP_KURTOSIS ,
    DPP_MSE ,       // sample entropy
    DPP_ENVELOPE ,  // filter-Hilbert magnitude summary (mean/SD/CV)
    DPP_PLV ,       // phase-locking value between two (channel,band) inputs
    DPP_COH ,       // coherence between two channels
    DPP_PSI         // phase slope index between two channels
  };

// a named bandpass definition, used both to pre-filter for ENVELOPE/PLV
// (via the causal padded-filter path, stats/dpp-filter.h) and, for
// COH/PSI, purely as a convenient named (lwr,upr) range for post-hoc
// band-summarizing an already-computed cross-spectrum (no filtering
// applied in that case)
struct dpp_filter_t
{
  std::string name;
  double lwr, upr;
  double ripple; // Kaiser ripple (default 0.02)
  double tw;     // transition width, Hz (default 1)
};

struct dpp_channel_t
{
  std::string ch;
  bool has_prefilter;
  double prefilter_lwr, prefilter_upr;
  dpp_channel_t() : has_prefilter(false), prefilter_lwr(0), prefilter_upr(0) { }
};

struct dpp_spec_t
{
  std::string block;
  dpp_feature_t ftr;

  // 1 channel for single-channel features, 2 for connectivity features
  std::vector<std::string> ch;

  // parallel to ch: named filter band per channel slot ("" = none/raw)
  std::vector<std::string> band;

  // this spec's own window length (seconds) -- already resolved from
  // windows= (block-level or global default) at read() time; one
  // dpp_spec_t per crossed window length
  double window_sec;

  std::map<std::string,std::string> arg;

  bool has( const std::string & key ) const;
  double narg( const std::string & key ) const;

  // number of output columns for this feature type
  int cols() const;

  // e.g. "PSD.C3", "PLV.C3-C4.sigma" -- everything except the trailing
  // ".w<N>.V<k>" (window length and column index), which the caller
  // appends once window_sec/column-index are known
  std::string label_root() const;
};

struct dpp_specs_t
{
  dpp_specs_t() : default_window_sec(0) , default_step_sec(0) { }

  void read( const std::string & f );
  void init();

  // zero-config default: PSD (canonical bands) + SLOPE + HJORTH,
  // independently for each channel in 'channels', at a single window/step
  void init_default( const std::vector<std::string> & channels );

  // CLI-inline shorthand: extend/override init_default()'s spec without a
  // file -- windows=, filters=name:lo-hi,..., features=PSD,HJORTH,...,
  // prefilter=lo-hi (applied to every channel already declared)
  void apply_inline_overrides( const std::string & windows_arg ,
				const std::string & filters_arg ,
				const std::string & features_arg ,
				const std::string & prefilter_arg );

  static std::map<std::string,dpp_feature_t> lab2ftr;
  static std::map<dpp_feature_t,std::string> ftr2lab;

  std::map<std::string,dpp_channel_t> chs;
  std::map<std::string,dpp_filter_t> filters;
  std::vector<dpp_spec_t> specs;

  double default_window_sec;
  double default_step_sec;

  bool loaded() const { return specs.size() != 0; }

  // does 'ch' exist (as a declared channel)?
  bool has_channel( const std::string & ch ) const { return chs.find( ch ) != chs.end(); }

  // does 'name' exist (as a declared filter band)?
  bool has_filter( const std::string & name ) const { return filters.find( name ) != filters.end(); }

};

#endif
