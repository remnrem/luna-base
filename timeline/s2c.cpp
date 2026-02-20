
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

#include "timeline/timeline.h"
#include "param.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "db/db.h"
#include "helper/logger.h"
#include "miscmath/miscmath.h"

#include <vector>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <random>
#include <limits>
#include <stdexcept>


extern writer_t writer;
extern logger_t logger;



class s2a2_t {

public:

  struct s2a2_param_t {
    s2a2_param_t() = default;

    explicit s2a2_param_t(const param_t& param)
    {

    // Signal settings

    sig_is_seg = !param.has("all-by-all");

    // Annotations 

    add_waves = param.has( "waves" );
    
    if (add_waves && param.empty("waves")) {
      Helper::halt("no annotation label specified waves={label}");
    }
    
    wave_label = add_waves ? param.value("waves") : "";
    
    add_ch_inst_label = param.has("add-channel-inst-label")
      ? param.yesno("add-channel-inst-label")
      : false;

    add_ch_class_label = param.has("add-channel-class-label")
      ? param.yesno("add-channel-class-label")
      : false;

    add_halfwaves = false;
    if (param.has("half-waves")) {
      const std::string k = "half-waves";
      if (param.empty(k)) {
        add_halfwaves = true;
      } else {
        const std::string v = param.value(k);
        const bool bool_like =
          !v.empty() &&
          (v[0] == '0' || v[0] == '1' ||
           v[0] == 'y' || v[0] == 'Y' ||
           v[0] == 'n' || v[0] == 'N' ||
           v[0] == 't' || v[0] == 'T' ||
           v[0] == 'f' || v[0] == 'F');
        if (bool_like) {
          add_halfwaves = param.yesno(k);
        } else {
          add_halfwaves = true;
          halfwaves_label = v;
        }
      }
    }
    add_peak_points = param.has("peak-points");
    if (add_peak_points) {
      const std::string k = "peak-points";
      if (param.empty(k)) {
        Helper::halt("no annotation label specified " + k + "={label}");
      }
      peak_points_label = param.value(k);
    }

    add_wave_bins = param.has("waves-bins");
    wave_bins_label = add_wave_bins ? param.value("waves-bins") : "";
    if (add_wave_bins && param.empty("waves-bins")) {
      Helper::halt("no annotation label specified waves-bins={label}");
    }

    
    // default: segment on pos2neg zero-crossings
    pos2neg = param.has("pos2neg") ? param.yesno("pos2neg") : true;
    dir = pos2neg ? cross_dir_t::POS2NEG : cross_dir_t::NEG2POS;

    // Selection criteria
    sel_tmin = param.has("t-min");
    sel_tmax = param.has("t-max");
    th_tmin = sel_tmin ? param.requires_dbl("t-min") : 0.0;
    th_tmax = sel_tmax ? param.requires_dbl("t-max") : 0.0;

    // Map overall duration constraints to cycle ticks (used by cycle detection).
    // t-min/t-max are specified in seconds; convert to tp ticks.
    if (sel_tmin && th_tmin > 0.0) {
      min_cycle_ticks = static_cast<uint64_t>(llround(th_tmin / globals::tp_duration));
    }
    if (sel_tmax && th_tmax > 0.0) {
      max_cycle_ticks = static_cast<uint64_t>(llround(th_tmax / globals::tp_duration));
    }

    sel_tmin_neg = param.has("t-min-neg");
    sel_tmax_neg = param.has("t-max-neg");
    th_tmin_neg = sel_tmin_neg ? param.requires_dbl("t-min-neg") : 0.0;
    th_tmax_neg = sel_tmax_neg ? param.requires_dbl("t-max-neg") : 0.0;

    sel_tmin_pos = param.has("t-min-pos");
    sel_tmax_pos = param.has("t-max-pos");
    th_tmin_pos = sel_tmin_pos ? param.requires_dbl("t-min-pos") : 0.0;
    th_tmax_pos = sel_tmax_pos ? param.requires_dbl("t-max-pos") : 0.0;

    // if not otherwise specified, use half of overall min/max as default
    if (!param.has("no-halfwave-t")) {
      if (sel_tmin && !sel_tmin_neg) {
        th_tmin_neg = th_tmin / 2.0;
        sel_tmin_neg = true;
      }
      if (sel_tmax && !sel_tmax_neg) {
        th_tmax_neg = th_tmax / 2.0;
        sel_tmax_neg = true;
      }
      if (sel_tmin && !sel_tmin_pos) {
        th_tmin_pos = th_tmin / 2.0;
        sel_tmin_pos = true;
      }
      if (sel_tmax && !sel_tmax_pos) {
        th_tmax_pos = th_tmax / 2.0;
        sel_tmax_pos = true;
      }
    }

    // Map half-wave duration constraints to ticks
    if (sel_tmin_neg && th_tmin_neg > 0.0) {
      min_neg_ticks = static_cast<uint64_t>(llround(th_tmin_neg / globals::tp_duration));
    }
    if (sel_tmax_neg && th_tmax_neg > 0.0) {
      max_neg_ticks = static_cast<uint64_t>(llround(th_tmax_neg / globals::tp_duration));
    }
    if (sel_tmin_pos && th_tmin_pos > 0.0) {
      min_pos_ticks = static_cast<uint64_t>(llround(th_tmin_pos / globals::tp_duration));
    }
    if (sel_tmax_pos && th_tmax_pos > 0.0) {
      max_pos_ticks = static_cast<uint64_t>(llround(th_tmax_pos / globals::tp_duration));
    }

    // percentile-based amplitude , e.g. mag=20 means in top 20% of average MAG value
    // compared to all other waves
    sel_mag = param.has("mag-percentile");

    th_mag = sel_mag ? param.requires_dbl("mag-percentile") : 0.0;
    if (sel_mag && (th_mag <= 0.0 || th_mag > 1.0)) {
      Helper::halt("mag-percentile must be between 0 and 1");
    }

    // norm values and only take is above th_magz
    sel_magz = param.has("mag-z");
    th_magz = sel_magz ? param.requires_dbl("mag-z") : 0.0;
    use_mag = sel_mag || sel_magz;

    // bootstrap settings for mean bins
    // Respect explicit bootstrap=0/F; otherwise allow bootstrap-n/ci to imply bootstrap.
    if (param.has("bootstrap")) {
      do_bootstrap = param.yesno("bootstrap", true, true);
    } else {
      do_bootstrap = param.has("bootstrap-n") || param.has("bootstrap-ci");
    }
    bootstrap_n = param.has("bootstrap-n") ? param.requires_int("bootstrap-n") : 1000;
    if (bootstrap_n < 1) bootstrap_n = 1;
    bootstrap_ci = param.has("bootstrap-ci") ? param.requires_dbl("bootstrap-ci") : 0.95;
    if (bootstrap_ci <= 0.0 || bootstrap_ci >= 1.0) {
      Helper::halt("bootstrap-ci must be between 0 and 1");
    }

    // lag settings
    lag_window_s = param.has("lag-window") ? param.requires_dbl("lag-window") : 0.0;
    if (lag_window_s < 0.0) lag_window_s = 0.0;
    lag_use_abs = param.has("lag-abs") ? param.yesno("lag-abs") : false;

    // per-cycle metric output
    emit_cycle_metrics = param.has("emit-per-cycle") ?
      param.yesno("emit-per-cycle") : false;

    // phase outputs
    emit_ph_grid = param.has("emit-ph-grid") ? param.yesno("emit-ph-grid", true, true) : true;
    emit_ph_amp = param.has("emit-ph-amp") ?
      param.yesno("emit-ph-amp", false, true) : false;
    amp_bins = param.has("amp-bins") ? param.requires_int("amp-bins") : 10;
    if (amp_bins < 2) amp_bins = 2;

    // time-domain stats around seg anchor
    emit_time_domain = param.has("time-domain") ? param.yesno("time-domain", false, true) : false;
    time_window_s = param.has("time-window") ? param.requires_dbl("time-window") : 100.0;
    if (time_window_s < 0.0) time_window_s = 0.0;
    time_bin_s = param.has("time-bin") ? param.requires_dbl("time-bin") : 1.0;
    if (time_bin_s <= 0.0) time_bin_s = 1.0;
    time_min_n = param.has("time-min-n") ? param.requires_int("time-min-n") : 1;
    if (time_min_n < 1) time_min_n = 1;
    time_lock = param.has("time-lock") ? param.value("time-lock") : "pos";
    if (time_lock != "pos" && time_lock != "neg") {
      Helper::halt("time-lock must be one of: pos, neg");
    }

    emit_td_grid = param.has("emit-td-grid") ? param.yesno("emit-td-grid", true, true) : true;
    emit_td_summary = param.has("emit-td-summary") ? param.yesno("emit-td-summary", true, true) : true;

    // optional outputs
    emit_se = param.has("emit-se") ? param.yesno("emit-se", false, true) : false;
    emit_mad = param.has("emit-mad") ? param.yesno("emit-mad", false, true) : false;
    std::string emit_asym_key = "emit-asym";
    if (param.has("emit-asymm")) emit_asym_key = "emit-asymm";
    else if (param.has("emit-asymmetry")) emit_asym_key = "emit-asymmetry";
    emit_asym_metrics = param.has(emit_asym_key) ?
      param.yesno(emit_asym_key, false, true) : false;

    emit_seed_summary = param.has("emit-seed") ? param.yesno("emit-seed", true, true) : true;
    emit_sig_summary = param.has("emit-sig") ? param.yesno("emit-sig", true, true) : true;

  }


  // definitions;
  bool sig_is_seg;
  std::string wave_label;
  bool add_ch_inst_label;
  bool add_ch_class_label;
  bool pos2neg;
 
  bool sel_tmin;
  bool sel_tmax;
  double th_tmin;
  double th_tmax;
  
  bool sel_tmin_neg;
  bool sel_tmax_neg;
  double th_tmin_neg;
  double th_tmax_neg;

  bool sel_tmin_pos;
  bool sel_tmax_pos;
  double th_tmin_pos;
  double th_tmax_pos;

  bool sel_mag;
  double th_mag;
  bool sel_magz;
  double th_magz;
  bool use_mag;

  // bootstrap for mean bins
  bool do_bootstrap = false;
  int bootstrap_n = 1000;
  double bootstrap_ci = 0.95;

  // lag calculation
  double lag_window_s = 0.0;   // 0 => use full cycle [t0,t1]
  bool lag_use_abs = false;

  // per-cycle metrics output
  bool emit_cycle_metrics = false;

  // phase outputs
  bool emit_ph_grid = true;
  bool emit_ph_amp = false;
  int amp_bins = 10;

  // time-domain stats output
  bool emit_time_domain = false;
  double time_window_s = 100.0;
  double time_bin_s = 1.0;
  int time_min_n = 1;
  std::string time_lock = "pos";
  bool emit_td_grid = true;
  bool emit_td_summary = true;
  bool emit_se = false;
  bool emit_mad = false;
  bool emit_asym_metrics = false;

  // summary output controls
  bool emit_seed_summary = true;
  bool emit_sig_summary = true;

  // annotation outputs
  bool add_waves = false;
  bool add_halfwaves = false;
  std::string halfwaves_label;
  bool add_peak_points = false;
  std::string peak_points_label;
  bool add_wave_bins = false;
  std::string wave_bins_label;

  // cycle detection
  // baseline for crossings
  bool use_epoch_median_zero = false;
  double zero = 0.0;

  // crossing direction
  enum class cross_dir_t { POS2NEG, NEG2POS } dir = cross_dir_t::POS2NEG;

  // sample rate
  double sr_hz;
  
  // hysteresis
  bool hysteresis = false;
  double h = 0.0;                 // absolute hysteresis half-width
  double h_frac_mad = 0.0;        // if >0 and h==0 => h = h_frac_mad * MAD(epoch)

  // duration constraints (in tp ticks)
  uint64_t min_cycle_ticks = 0;   // 0 disables
  uint64_t max_cycle_ticks = 0;   // 0 disables
  uint64_t min_neg_ticks = 0;     // 0 disables
  uint64_t max_neg_ticks = 0;     // 0 disables
  uint64_t min_pos_ticks = 0;     // 0 disables
  uint64_t max_pos_ticks = 0;     // 0 disables

  // debounce / cleanup
  uint64_t min_sep_ticks = 0;     // minimum time between crossings; 0 => derive from min_cycle_ticks/4 (heuristic)
  
  
};

struct cycle_bound_t {
  int i0 = -1;   // index into sample arrays (crossing between i0-1 and i0)
  int i1 = -1;   // next crossing index
  // sub-sample crossing times
  uint64_t t0 = 0;
  uint64_t t1 = 0;
  uint64_t t_mid = 0;
  bool has_mid = false;

  // peak markers within cycle
  int i_pos = -1;
  int i_neg = -1;
  uint64_t t_pos = 0;
  uint64_t t_neg = 0;
  double v_pos = std::numeric_limits<double>::quiet_NaN();
  double v_neg = std::numeric_limits<double>::quiet_NaN();

  // derived metrics
  double rel_pos = std::numeric_limits<double>::quiet_NaN();
  double rel_neg = std::numeric_limits<double>::quiet_NaN();
  double rel_i_pos = std::numeric_limits<double>::quiet_NaN();
  double rel_i_neg = std::numeric_limits<double>::quiet_NaN();
  double dt_pos_s = std::numeric_limits<double>::quiet_NaN();
  double dt_neg_s = std::numeric_limits<double>::quiet_NaN();
  double pos_slope = std::numeric_limits<double>::quiet_NaN();
  double neg_slope = std::numeric_limits<double>::quiet_NaN();
  double pos_slope_norm = std::numeric_limits<double>::quiet_NaN();
  double neg_slope_norm = std::numeric_limits<double>::quiet_NaN();
};

struct s2a2_out_t {
  std::vector<cycle_bound_t> cycles;
  // [channel][cycle][bin]
  std::vector<std::vector<std::vector<double>>> bins;
  int nbins = 0;
};

//  Implements S2A2
static s2a2_out_t s2a2_proc(
  const Eigen::MatrixXd&,
  const std::vector<uint64_t>* tp,
  int idx,
  const std::vector<int>& chs_idx,
  const s2a2_param_t& par,
  const std::string& seg_label,
  const std::vector<std::string>& sig_labels
);
 private:
  static inline bool finite(double x);
  static double epoch_median(std::vector<double> v);
  static double epoch_mad(const std::vector<double>& x, double med);
  static double median_dbl(std::vector<double> v);
  static double mad_dbl(const std::vector<double>& x, double med);
  static double sd_dbl(const std::vector<double>& x);
  static double mean_dbl(const std::vector<double>& x);
  static uint64_t interp_cross_time(uint64_t tA, uint64_t tB, double a, double b);
  static double half_width_s(const std::vector<double>& sig,
                             const std::vector<uint64_t>& tp,
                             int i0,
                             int i1,
                             int i_peak,
                             double level);
  struct window_metrics_t {
    double rms = std::numeric_limits<double>::quiet_NaN();
    double p2p = std::numeric_limits<double>::quiet_NaN();
    double duty = std::numeric_limits<double>::quiet_NaN();
    double max_slope = std::numeric_limits<double>::quiet_NaN();
  };
  static window_metrics_t window_metrics(const std::vector<double>& sig,
                                         const std::vector<uint64_t>& tp,
                                         uint64_t t_center,
                                         double window_s,
                                         double zero,
                                         int sr_hz);
  static bool window_range_ok(const std::vector<uint64_t>& tp,
                              int sr_hz,
                              uint64_t t_lo,
                              uint64_t t_hi,
                              size_t* i0,
                              size_t* i1);
  static std::vector<double> zscore_finite(const std::vector<double>& x);
  static std::vector<cycle_bound_t> filter_cycles_by_mag(
    const std::vector<double>& seed,
    const s2a2_param_t& par_work,
    const std::vector<cycle_bound_t>& cycles,
    int* mag_excl,
    std::vector<double>* mag_vals
  );
  static uint64_t median_dt_ticks(const std::vector<uint64_t>& tp);
  static double interp_value_at_tp_ld(const std::vector<uint64_t>& tp,
                                      const std::vector<double>& sig,
                                      long double t);
  static std::vector<double> bin_cycle_4pt(const std::vector<uint64_t>& tp,
                                           const std::vector<double>& sig,
                                           uint64_t t_start,
                                           uint64_t t_mid1,
                                           uint64_t t_mid2,
                                           uint64_t t_end,
                                           int nbins);
 public:
  static uint64_t phase_to_time4pt(uint64_t t_start,
                                   uint64_t t_mid1,
                                   uint64_t t_mid2,
                                   uint64_t t_end,
                                   long double ph);
 private:
  static std::vector<std::vector<std::vector<double>>> step4_piecewise_bins(
    const Eigen::MatrixXd& X,
    const std::vector<uint64_t>& tp,
    const std::vector<int>& chs_idx,
    const std::vector<cycle_bound_t>& cycles,
    const s2a2_param_t& par_work,
    int nbins
  );
  static std::vector<double> mean_bins_for_signal(
    const std::vector<uint64_t>& tp,
    const std::vector<double>& sig,
    const std::vector<cycle_bound_t>& cycles,
    const s2a2_param_t& par_work,
    int nbins
  );
  static int find_max_index_in_window(const std::vector<uint64_t>& tp,
                                      const std::vector<double>& sig,
                                      uint64_t t_lo,
                                      uint64_t t_hi,
                                      bool use_abs);
  static int find_min_index_in_window(const std::vector<uint64_t>& tp,
                                      const std::vector<double>& sig,
                                      uint64_t t_lo,
                                      uint64_t t_hi);
  static std::pair<int, double> crosscorr_lag_bins(const std::vector<double>& a,
                                                   const std::vector<double>& b);
  static int find_nearest_local_extremum(const std::vector<uint64_t>& tp,
                                         const std::vector<double>& sig,
                                         uint64_t t_lo,
                                         uint64_t t_hi,
                                         uint64_t t_target,
                                         bool find_max,
                                         bool use_abs);
  static double fit_sincos_amplitude(const std::vector<double>& y);
  static int argmax_index(const std::vector<double>& v);
  static bool refine_peak_parabolic(const std::vector<double>& sig,
                                    const std::vector<uint64_t>& tp,
                                    int i0,
                                    int i1,
                                    int i_peak,
                                    bool is_max,
                                    uint64_t& t_peak_out,
                                    double& v_peak_out);
  static bool find_opposite_cross_time(const std::vector<double>& seed,
                                       const std::vector<uint64_t>& tp,
                                       int i0,
                                       int i1,
                                       const s2a2_param_t& par,
                                       double zero,
                                       double h,
                                       uint64_t& t_cross);
  static void step1_extract_seed(const Eigen::MatrixXd& X,
                                 const std::vector<uint64_t>* tp,
                                 int idx,
                                 const std::vector<int>& chs_idx,
                                 const s2a2_param_t& par_in,
                                 std::vector<double>& seed,
                                 double& zero,
                                 double& h,
                                 uint64_t& dt_med,
                                 s2a2_param_t& par_work);
  static std::vector<cycle_bound_t> step2_detect_cycles(const std::vector<double>& seed,
                                                       const std::vector<uint64_t>& tp,
                                                       const s2a2_param_t& par,
                                                       double zero,
                                                       double h,
                                                       size_t* n_putative,
                                                       size_t* n_after_duration);
  static void step3_mark_peaks(const std::vector<double>& seed,
                               const std::vector<uint64_t>& tp,
                               std::vector<cycle_bound_t>& cycles);
  static void step3_derive_metrics(const std::vector<double>& seed,
                                   const std::vector<uint64_t>& tp,
                                   const s2a2_param_t& par_work,
                                   double zero,
                                   double h,
                                   std::vector<cycle_bound_t>& cycles);
};


//  Driver for S2C

void timeline_t::signal2cycle( const param_t & param )
{
  

  
  //
  // signal(s) to use: assume narrow-band inputs, pre-filtered
  //
  
  std::string signal_label = param.requires( "sig" );
  std::string segnal_label = param.has( "seed" ) ? param.value( "seed" ) : signal_label;
  
  signal_list_t signals = edf->header.signal_list(signal_label);
  signal_list_t segnals = edf->header.signal_list(segnal_label);

  logger << "  generating S2A2 annotations for signal(s) ["
	 << signal_label << "] and segmenting seed-signal(s) ["
	 << segnal_label << "]\n";
  
    
  if (segnals.size() == 0 || signals.size() == 0 ) {
    return;
  }

  const int nsig = static_cast<int>(signals.size());

  const int nseg = static_cast<int>(segnals.size());

  //
  // get all data
  //
  
  std::set<std::string> allsigs_set;
  for (int s = 0; s < nsig; ++s) allsigs_set.insert(signals.label(s));
  for (int s = 0; s < nseg; ++s) allsigs_set.insert(segnals.label(s));
  std::string allsigs_label = Helper::stringize(allsigs_set);
  signal_list_t allsigs = edf->header.signal_list(allsigs_label);
    
  // signals data
  eigen_matslice_t mslice(*edf, allsigs, wholetrace());
  const Eigen::MatrixXd& X = mslice.data_ref();
  
  // time points
  const std::vector<uint64_t>* tp = mslice.ptimepoints();
    
  // parameters
  s2a2_t::s2a2_param_t par(param);

  // log parameters once per run (grouped)
  logger << "  S2C params\n";
  logger << "    signals   : sig_is_seed=" << (par.sig_is_seg ? "T" : "F") << "\n";
  logger << "    annots    : waves=" << (par.add_waves ? par.wave_label : ".")
         << " add-channel-inst-label=" << (par.add_ch_inst_label ? "T" : "F")
         << " add-channel-class-label=" << (par.add_ch_class_label ? "T" : "F")
         << " half-waves=" << (par.add_halfwaves ? (par.halfwaves_label.empty() ? "T" : par.halfwaves_label) : ".")
         << " peak-points=" << (par.add_peak_points ? par.peak_points_label : ".")
         << " waves-bins=" << (par.add_wave_bins ? par.wave_bins_label : ".") << "\n";
  logger << "    crossings : pos2neg=" << (par.pos2neg ? "T" : "F") << "\n";
  logger << "    durations : t-min=" << (par.sel_tmin ? Helper::dbl2str(par.th_tmin) : ".")
         << " t-max=" << (par.sel_tmax ? Helper::dbl2str(par.th_tmax) : ".")
         << " t-min-neg=" << (par.sel_tmin_neg ? Helper::dbl2str(par.th_tmin_neg) : ".")
         << " t-max-neg=" << (par.sel_tmax_neg ? Helper::dbl2str(par.th_tmax_neg) : ".")
         << " t-min-pos=" << (par.sel_tmin_pos ? Helper::dbl2str(par.th_tmin_pos) : ".")
         << " t-max-pos=" << (par.sel_tmax_pos ? Helper::dbl2str(par.th_tmax_pos) : ".") << "\n";
  logger << "    magnitude : mag-percentile=" << (par.sel_mag ? Helper::dbl2str(par.th_mag) : ".")
         << " mag-z=" << (par.sel_magz ? Helper::dbl2str(par.th_magz) : ".") << "\n";
  logger << "    bootstrap : bootstrap=" << (par.do_bootstrap ? "T" : "F")
         << " bootstrap-n=" << par.bootstrap_n
         << " bootstrap-ci=" << Helper::dbl2str(par.bootstrap_ci) << "\n";
  logger << "    lag       : lag-window=" << Helper::dbl2str(par.lag_window_s)
         << " lag-abs=" << (par.lag_use_abs ? "T" : "F") << "\n";
  logger << "    output    : emit-per-cycle=" << (par.emit_cycle_metrics ? "T" : "F")
         << " emit-seed=" << (par.emit_seed_summary ? "T" : "F")
         << " emit-sig=" << (par.emit_sig_summary ? "T" : "F")
         << " emit-ph-grid=" << (par.emit_ph_grid ? "T" : "F")
         << " emit-ph-amp=" << (par.emit_ph_amp ? "T" : "F")
         << " amp-bins=" << par.amp_bins
         << " emit-se=" << (par.emit_se ? "T" : "F")
         << " emit-mad=" << (par.emit_mad ? "T" : "F")
         << " emit-asym=" << (par.emit_asym_metrics ? "T" : "F") << "\n";
  logger << "    time-dom  : time-domain=" << (par.emit_time_domain ? "T" : "F")
         << " time-window=" << Helper::dbl2str(par.time_window_s)
         << " time-bin=" << Helper::dbl2str(par.time_bin_s)
         << " time-min-n=" << par.time_min_n
         << " time-lock=" << par.time_lock
         << " emit-td-grid=" << (par.emit_td_grid ? "T" : "F")
         << " emit-td-summary=" << (par.emit_td_summary ? "T" : "F") << "\n";
  
  // iterate over segmenting-signals
  for (int s = 0; s < nseg; ++s) {
    const std::string seed_label = segnals.label(s);
    int seed_idx = allsigs.find(seed_label);
    if (seed_idx == -1) Helper::halt("internal error, cannot index seed " + seed_label);
    
    // get all other signals to be compared w.r.t. this segmenting signal
    std::vector<int> chs_idx;

    const bool use_only_seg =
      par.sig_is_seg && nsig == 1 && signal_label == segnal_label;

    if (use_only_seg) 
      chs_idx.push_back(seed_idx);
    else
      {
	for (int s2 = 0; s2 < nsig; ++s2) {
	  int idx = allsigs.find(signals.label(s2));
	  if (idx == -1) Helper::halt("internal error, cannot index sig");
	  chs_idx.push_back(idx);
	}
      }
    
    // we now have:
    //   seed_idx - one channel to make segments
    //   chs_idx[] - one or more channels to summarize w.r.t. those
    
    par.sr_hz = edf->header.sampling_freq(segnals(s));
    
    std::vector<std::string> sig_labels;
    sig_labels.reserve(chs_idx.size());
    for (int k = 0; k < static_cast<int>(chs_idx.size()); ++k) {
      int col = chs_idx[k];
      sig_labels.push_back(col >= 0 && col < allsigs.size() ? allsigs.label(col) : ".");
    }
    
    // process
    s2a2_t::s2a2_out_t res =
      s2a2_t::s2a2_proc(X, tp, seed_idx, chs_idx, par, seed_label, sig_labels);

    // annotations
    annot_t * a_full = nullptr;
    if ( par.add_waves )
      a_full = par.add_ch_class_label ?
	edf->annotations->add(par.wave_label + "_" + seed_label) :
	edf->annotations->add(par.wave_label);

    annot_t * a_pos = nullptr;
    annot_t * a_neg = nullptr;
    if (par.add_halfwaves) {
      const std::string base = par.halfwaves_label.empty() ? par.wave_label : par.halfwaves_label;
      a_pos = par.add_ch_class_label ?
        edf->annotations->add(base + "_POS_" + seed_label) :
        edf->annotations->add(base + "_POS");
      a_neg = par.add_ch_class_label ?
        edf->annotations->add(base + "_NEG_" + seed_label) :
        edf->annotations->add(base + "_NEG");
    }

    annot_t * a_pos_peak = nullptr;
    annot_t * a_neg_peak = nullptr;
    if (par.add_peak_points) {
      a_pos_peak = par.add_ch_class_label ?
        edf->annotations->add(par.peak_points_label + "_POS_PEAK_" + seed_label) :
        edf->annotations->add(par.peak_points_label + "_POS_PEAK");
      a_neg_peak = par.add_ch_class_label ?
        edf->annotations->add(par.peak_points_label + "_NEG_PEAK_" + seed_label) :
        edf->annotations->add(par.peak_points_label + "_NEG_PEAK");
    }

    annot_t * a_bins = nullptr;
    if (par.add_wave_bins) {
      a_bins = par.add_ch_class_label ?
        edf->annotations->add(par.wave_bins_label + "_" + seed_label) :
        edf->annotations->add(par.wave_bins_label);
    }

    // add annotations

    const bool neg_first = (par.dir == s2a2_t::s2a2_param_t::cross_dir_t::POS2NEG);

    for (const auto& c : res.cycles) {
      if (c.t1 <= c.t0) continue;
      interval_t full_iv(c.t0, c.t1);

      if ( par.add_waves ) 
	a_full->add(par.add_ch_inst_label ? seed_label : ".", full_iv, seed_label);
      
      if (par.add_halfwaves && c.has_mid && c.t_mid > c.t0 && c.t_mid < c.t1) {
        interval_t first_iv(c.t0, c.t_mid);
        interval_t second_iv(c.t_mid, c.t1);
        if (neg_first) {
          a_neg->add(par.add_ch_inst_label ? seed_label : "NEG", first_iv, seed_label);
          a_pos->add(par.add_ch_inst_label ? seed_label : "POS", second_iv, seed_label);
        } else {
          a_pos->add(par.add_ch_inst_label ? seed_label : "POS", first_iv, seed_label);
          a_neg->add(par.add_ch_inst_label ? seed_label : "NEG", second_iv, seed_label);
        }
      }
      if (par.add_peak_points) {
        if (c.i_pos >= 0) {
          a_pos_peak->add(par.add_ch_inst_label ? seed_label : "POS",
                          interval_t(c.t_pos, c.t_pos),
                          seed_label);
        }
        if (c.i_neg >= 0) {
          a_neg_peak->add(par.add_ch_inst_label ? seed_label : "NEG",
                          interval_t(c.t_neg, c.t_neg),
                          seed_label);
        }
      }

      if (par.add_wave_bins && a_bins != nullptr) {
        const uint64_t t0 = c.t0;
        const uint64_t t1 = c.t1;
        uint64_t t_mid1 = 0;
        uint64_t t_mid2 = 0;
        if (par.dir == s2a2_t::s2a2_param_t::cross_dir_t::POS2NEG) {
          t_mid1 = c.t_neg;
          t_mid2 = c.t_pos;
        } else {
          t_mid1 = c.t_pos;
          t_mid2 = c.t_neg;
        }
        if (t1 > t0 && t0 < t_mid1 && t_mid1 < t_mid2 && t_mid2 < t1) {
          for (int b = 0; b < 12; ++b) {
            long double f0 = static_cast<long double>(b) / 12.0L;
            long double f1 = static_cast<long double>(b + 1) / 12.0L;
            uint64_t tb0 = s2a2_t::phase_to_time4pt(t0, t_mid1, t_mid2, t1, f0);
            uint64_t tb1 = (b == 11) ? t1 : s2a2_t::phase_to_time4pt(t0, t_mid1, t_mid2, t1, f1);
            if (tb1 <= tb0) continue;
            interval_t bin_iv(tb0, tb1);
            char lab[4];
            snprintf(lab, sizeof(lab), "B%02d", b + 1);
            a_bins->add(lab, bin_iv, seed_label);
          }
        }
      }
    }

  } // next segmenting-signal
  
}




inline bool s2a2_t::finite(double x) { return std::isfinite(x); }

double s2a2_t::epoch_median(std::vector<double> v)
{
  v.erase(std::remove_if(v.begin(), v.end(), [](double x) { return !finite(x); }), v.end());
  if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
  const size_t n = v.size();
  const size_t k = n / 2;
  std::nth_element(v.begin(), v.begin() + k, v.end());
  double med = v[k];
  if (n % 2 == 0) {
    std::nth_element(v.begin(), v.begin() + k - 1, v.end());
    med = 0.5 * (med + v[k - 1]);
  }
  return med;
}

// robust MAD about median
double s2a2_t::epoch_mad(const std::vector<double>& x, double med)
{
  std::vector<double> d;
  d.reserve(x.size());
  for (double v : x) {
    if (finite(v)) d.push_back(std::fabs(v - med));
  }
  return epoch_median(std::move(d));
}

double s2a2_t::median_dbl(std::vector<double> v)
{
  v.erase(std::remove_if(v.begin(), v.end(), [](double x) { return !finite(x); }), v.end());
  if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
  const size_t n = v.size();
  const size_t k = n / 2;
  std::nth_element(v.begin(), v.begin() + k, v.end());
  double med = v[k];
  if (n % 2 == 0) {
    std::nth_element(v.begin(), v.begin() + k - 1, v.end());
    med = 0.5 * (med + v[k - 1]);
  }
  return med;
}

double s2a2_t::mad_dbl(const std::vector<double>& x, double med)
{
  std::vector<double> d;
  d.reserve(x.size());
  for (double v : x) {
    if (finite(v)) d.push_back(std::fabs(v - med));
  }
  return median_dbl(std::move(d));
}

double s2a2_t::sd_dbl(const std::vector<double>& x)
{
  double mu = 0.0;
  size_t n = 0;
  for (double v : x) {
    if (!finite(v)) continue;
    mu += v;
    ++n;
  }
  if (n < 2) return std::numeric_limits<double>::quiet_NaN();
  mu /= static_cast<double>(n);
  double var = 0.0;
  for (double v : x) {
    if (!finite(v)) continue;
    double d = v - mu;
    var += d * d;
  }
  var /= static_cast<double>(n - 1);
  return std::sqrt(var);
}

double s2a2_t::mean_dbl(const std::vector<double>& x)
{
  double mu = 0.0;
  size_t n = 0;
  for (double v : x) {
    if (!finite(v)) continue;
    mu += v;
    ++n;
  }
  if (n == 0) return std::numeric_limits<double>::quiet_NaN();
  return mu / static_cast<double>(n);
}

double s2a2_t::half_width_s(
  const std::vector<double>& sig,
  const std::vector<uint64_t>& tp,
  int i0,
  int i1,
  int i_peak,
  double level
)
{
  if (i0 < 0 || i1 < 0 || i_peak < 0) return std::numeric_limits<double>::quiet_NaN();
  if (i0 >= i1) return std::numeric_limits<double>::quiet_NaN();
  if (i_peak < i0 || i_peak > i1) return std::numeric_limits<double>::quiet_NaN();

  double t_left = std::numeric_limits<double>::quiet_NaN();
  double t_right = std::numeric_limits<double>::quiet_NaN();

  for (int i = i_peak; i > i0; --i) {
    double a = sig[static_cast<size_t>(i - 1)] - level;
    double b = sig[static_cast<size_t>(i)] - level;
    if (!finite(a) || !finite(b)) continue;
    if (a == 0.0) { t_left = tp[static_cast<size_t>(i - 1)] * globals::tp_duration; break; }
    if ((a < 0.0 && b > 0.0) || (a > 0.0 && b < 0.0) || b == 0.0) {
      uint64_t t = interp_cross_time(tp[static_cast<size_t>(i - 1)],
                                     tp[static_cast<size_t>(i)],
                                     a,
                                     b);
      t_left = static_cast<double>(t) * globals::tp_duration;
      break;
    }
  }

  for (int i = i_peak; i < i1; ++i) {
    double a = sig[static_cast<size_t>(i)] - level;
    double b = sig[static_cast<size_t>(i + 1)] - level;
    if (!finite(a) || !finite(b)) continue;
    if (b == 0.0) { t_right = tp[static_cast<size_t>(i + 1)] * globals::tp_duration; break; }
    if ((a < 0.0 && b > 0.0) || (a > 0.0 && b < 0.0) || a == 0.0) {
      uint64_t t = interp_cross_time(tp[static_cast<size_t>(i)],
                                     tp[static_cast<size_t>(i + 1)],
                                     a,
                                     b);
      t_right = static_cast<double>(t) * globals::tp_duration;
      break;
    }
  }

  if (!finite(t_left) || !finite(t_right)) return std::numeric_limits<double>::quiet_NaN();
  if (t_right < t_left) return std::numeric_limits<double>::quiet_NaN();
  return t_right - t_left;
}

bool s2a2_t::window_range_ok(
  const std::vector<uint64_t>& tp,
  int sr_hz,
  uint64_t t_lo,
  uint64_t t_hi,
  size_t* i0,
  size_t* i1
)
{
  if (tp.empty()) return false;
  if (t_hi < t_lo) return false;
  auto it0 = std::lower_bound(tp.begin(), tp.end(), t_lo);
  auto it1 = std::upper_bound(tp.begin(), tp.end(), t_hi);
  if (it0 == tp.end()) return false;
  size_t j0 = static_cast<size_t>(it0 - tp.begin());
  size_t j1 = (it1 == tp.begin()) ? 0 : static_cast<size_t>((it1 - tp.begin()) - 1);
  if (j0 >= tp.size() || j1 >= tp.size() || j0 > j1) return false;
  if (tp[j0] > t_lo) return false;
  if (tp[j1] < t_hi) return false;
  if (sr_hz > 0 && j1 > j0) {
    for (size_t i = j0 + 1; i <= j1; ++i) {
      if (timeline_t::discontinuity(tp, sr_hz, static_cast<int>(i - 1), static_cast<int>(i))) {
        return false;
      }
    }
  }
  if (i0) *i0 = j0;
  if (i1) *i1 = j1;
  return true;
}

s2a2_t::window_metrics_t s2a2_t::window_metrics(
  const std::vector<double>& sig,
  const std::vector<uint64_t>& tp,
  uint64_t t_center,
  double window_s,
  double zero,
  int sr_hz
)
{
  window_metrics_t out;
  if (sig.empty() || tp.empty() || sig.size() != tp.size()) return out;
  if (window_s <= 0.0) return out;

  const double half = 0.5 * window_s;
  const uint64_t t_lo = (t_center > static_cast<uint64_t>(llround(half / globals::tp_duration)))
    ? t_center - static_cast<uint64_t>(llround(half / globals::tp_duration))
    : 0;
  const uint64_t t_hi = t_center + static_cast<uint64_t>(llround(half / globals::tp_duration));

  size_t i0 = 0;
  size_t i1 = 0;
  if (!window_range_ok(tp, sr_hz, t_lo, t_hi, &i0, &i1)) return out;

  double vmin = std::numeric_limits<double>::infinity();
  double vmax = -std::numeric_limits<double>::infinity();
  double ss = 0.0;
  size_t n = 0;
  size_t npos = 0;
  double max_slope = 0.0;
  bool has_slope = false;

  for (size_t i = i0; i <= i1; ++i) {
    double v = sig[i];
    if (!finite(v)) continue;
    if (v < vmin) vmin = v;
    if (v > vmax) vmax = v;
    ss += v * v;
    ++n;
    if (v > zero) ++npos;
    if (i > i0) {
      double v0 = sig[i - 1];
      if (finite(v0)) {
        uint64_t dt_ticks = tp[i] > tp[i - 1] ? (tp[i] - tp[i - 1]) : 0;
        if (dt_ticks > 0) {
          double dt = dt_ticks * globals::tp_duration;
          double s = std::fabs((v - v0) / dt);
          if (!has_slope || s > max_slope) {
            max_slope = s;
            has_slope = true;
          }
        }
      }
    }
  }

  if (n > 0) {
    out.rms = std::sqrt(ss / static_cast<double>(n));
    out.duty = static_cast<double>(npos) / static_cast<double>(n);
    if (finite(vmin) && finite(vmax)) out.p2p = vmax - vmin;
    if (has_slope) out.max_slope = max_slope;
  }

  return out;
}

std::vector<double> s2a2_t::zscore_finite(const std::vector<double>& x)
{
  double mu = 0.0;
  size_t n = 0;
  for (double v : x) {
    if (!finite(v)) continue;
    mu += v;
    ++n;
  }
  if (n < 2) return x;
  mu /= static_cast<double>(n);
  double var = 0.0;
  for (double v : x) {
    if (!finite(v)) continue;
    double d = v - mu;
    var += d * d;
  }
  var /= static_cast<double>(n - 1);
  double sd = std::sqrt(var);
  if (!finite(sd) || sd == 0.0) return x;
  std::vector<double> z(x.size(), std::numeric_limits<double>::quiet_NaN());
  for (size_t i = 0; i < x.size(); ++i) {
    if (finite(x[i])) z[i] = (x[i] - mu) / sd;
  }
  return z;
}

std::vector<s2a2_t::cycle_bound_t> s2a2_t::filter_cycles_by_mag(
  const std::vector<double>& seed,
  const s2a2_param_t& par_work,
  const std::vector<cycle_bound_t>& cycles,
  int* mag_excl,
  std::vector<double>* mag_vals
)
{
  if (mag_excl) *mag_excl = 0;
  if (!par_work.use_mag || cycles.empty()) return cycles;

  std::vector<double> mags;
  mags.reserve(cycles.size());
  for (const auto& c : cycles) {
    if (c.i0 < 0 || c.i1 < 0 ||
        c.i0 >= static_cast<int>(seed.size()) ||
        c.i1 > static_cast<int>(seed.size()) ||
        c.i1 <= c.i0) {
      mags.push_back(std::numeric_limits<double>::quiet_NaN());
      continue;
    }
    double ss = 0.0;
    size_t n = 0;
    for (int i = c.i0; i < c.i1; ++i) {
      double v = seed[static_cast<size_t>(i)];
      if (!finite(v)) continue;
      ss += v * v;
      ++n;
    }
    if (n == 0) {
      mags.push_back(std::numeric_limits<double>::quiet_NaN());
    } else {
      mags.push_back(std::sqrt(ss / static_cast<double>(n)));
    }
  }

  if (par_work.sel_magz) mags = zscore_finite(mags);
  if (mag_vals) *mag_vals = mags;

  double percentile = 0.0;
  if (par_work.sel_mag) {
    std::vector<double> finite_mags;
    finite_mags.reserve(mags.size());
    for (double v : mags) {
      if (finite(v)) finite_mags.push_back(v);
    }
    if (finite_mags.size() > 1) {
      percentile = MiscMath::percentile(finite_mags, par_work.th_mag);
    }
  }

  std::vector<cycle_bound_t> out;
  out.reserve(cycles.size());
  for (size_t i = 0; i < cycles.size(); ++i) {
    bool okay = finite(mags[i]);
    if (okay && par_work.sel_mag && mags[i] < percentile) okay = false;
    if (okay && par_work.sel_magz && mags[i] < par_work.th_magz) okay = false;
    if (okay) {
      out.push_back(cycles[i]);
    } else if (mag_excl) {
      ++(*mag_excl);
    }
  }

  return out;
}

// Robust median dt from tp
uint64_t s2a2_t::median_dt_ticks(const std::vector<uint64_t>& tp)
{
  if (tp.size() < 2) return 0;
  std::vector<uint64_t> d;
  d.reserve(tp.size() - 1);
  for (size_t i = 1; i < tp.size(); ++i) {
    if (tp[i] > tp[i - 1]) d.push_back(tp[i] - tp[i - 1]);
  }
  if (d.empty()) return 0;
  const size_t k = d.size() / 2;
  std::nth_element(d.begin(), d.begin() + k, d.end());
  uint64_t med = d[k];
  if (d.size() % 2 == 0) {
    std::nth_element(d.begin(), d.begin() + k - 1, d.end());
    med = (med + d[k - 1]) / 2;
  }
  return med;
}

// Compute sub-sample crossing time using linear interpolation in amplitude,
// then linear interpolation in tp. Returns tp at crossing.
uint64_t s2a2_t::interp_cross_time(uint64_t tA, uint64_t tB, double a, double b)
{
  // a and b are (s - zero) at endpoints; crossing when == 0
  // frac = a / (a - b) in [0,1] if true crossing
  double denom = (a - b);
  if (denom == 0.0) return tA;
  double frac = a / denom;
  if (!finite(frac)) return tA;
  if (frac < 0.0) frac = 0.0;
  if (frac > 1.0) frac = 1.0;
  long double t = static_cast<long double>(tA)
    + static_cast<long double>(frac) * static_cast<long double>(tB - tA);
  if (t < 0) t = 0;
  if (t > (long double)std::numeric_limits<uint64_t>::max())
    t = (long double)std::numeric_limits<uint64_t>::max();
  return static_cast<uint64_t>(llround(static_cast<double>(t)));
}

// Linear interpolation of signal value at arbitrary (non-integer) time tick.
double s2a2_t::interp_value_at_tp_ld(
  const std::vector<uint64_t>& tp,
  const std::vector<double>& sig,
  long double t
)
{
  const size_t n = tp.size();
  if (n == 0 || sig.size() != n) return std::numeric_limits<double>::quiet_NaN();
  if (t <= static_cast<long double>(tp.front())) return sig.front();
  if (t >= static_cast<long double>(tp.back())) return sig.back();

  auto it = std::lower_bound(
    tp.begin(),
    tp.end(),
    t,
    [](uint64_t a, long double val) { return static_cast<long double>(a) < val; }
  );
  size_t i1 = static_cast<size_t>(it - tp.begin());
  if (i1 == 0) return sig.front();
  size_t i0 = i1 - 1;

  long double t0 = static_cast<long double>(tp[i0]);
  long double t1 = static_cast<long double>(tp[i1]);
  if (t1 <= t0) return sig[i0];

  long double frac = (t - t0) / (t1 - t0);
  if (frac < 0.0L) frac = 0.0L;
  if (frac > 1.0L) frac = 1.0L;
  return sig[i0] + static_cast<double>(frac) * (sig[i1] - sig[i0]);
}

// Build a 4-point piecewise map (start crossing -> peak -> trough -> end crossing)
// and return a fixed-bin (e.g., 101) resample via linear interpolation.
std::vector<double> s2a2_t::bin_cycle_4pt(
  const std::vector<uint64_t>& tp,
  const std::vector<double>& sig,
  uint64_t t_start,
  uint64_t t_mid1,
  uint64_t t_mid2,
  uint64_t t_end,
  int nbins = 101
)
{
  std::vector<double> out;
  if (nbins <= 1) return out;
  if (tp.empty() || sig.size() != tp.size()) return out;
  if (!(t_start < t_mid1 && t_mid1 < t_mid2 && t_mid2 < t_end)) return out;

  out.resize(static_cast<size_t>(nbins), std::numeric_limits<double>::quiet_NaN());

  const long double p0 = 0.0L;
  // For a sine wave with zero-crossings at start/end, trough at 0.25 and peak at 0.75
  const long double p1 = 0.25L;
  const long double p2 = 0.75L;
  const long double p3 = 1.0L;

  for (int b = 0; b < nbins; ++b) {
    long double ph = static_cast<long double>(b) / static_cast<long double>(nbins - 1);
    long double t = static_cast<long double>(t_start);

    if (ph <= p1) {
      long double frac = (ph - p0) / (p1 - p0);
      t = static_cast<long double>(t_start) +
          frac * static_cast<long double>(t_mid1 - t_start);
    } else if (ph <= p2) {
      long double frac = (ph - p1) / (p2 - p1);
      t = static_cast<long double>(t_mid1) +
          frac * static_cast<long double>(t_mid2 - t_mid1);
    } else {
      long double frac = (ph - p2) / (p3 - p2);
      t = static_cast<long double>(t_mid2) +
          frac * static_cast<long double>(t_end - t_mid2);
    }

    out[static_cast<size_t>(b)] = interp_value_at_tp_ld(tp, sig, t);
  }

  return out;
}

uint64_t s2a2_t::phase_to_time4pt(
  uint64_t t_start,
  uint64_t t_mid1,
  uint64_t t_mid2,
  uint64_t t_end,
  long double ph
)
{
  if (!(t_start < t_mid1 && t_mid1 < t_mid2 && t_mid2 < t_end)) return t_start;
  if (ph < 0.0L) ph = 0.0L;
  if (ph > 1.0L) ph = 1.0L;

  const long double p0 = 0.0L;
  const long double p1 = 0.25L;
  const long double p2 = 0.75L;
  const long double p3 = 1.0L;

  long double t = static_cast<long double>(t_start);
  if (ph <= p1) {
    long double frac = (ph - p0) / (p1 - p0);
    t = static_cast<long double>(t_start) +
        frac * static_cast<long double>(t_mid1 - t_start);
  } else if (ph <= p2) {
    long double frac = (ph - p1) / (p2 - p1);
    t = static_cast<long double>(t_mid1) +
        frac * static_cast<long double>(t_mid2 - t_mid1);
  } else {
    long double frac = (ph - p2) / (p3 - p2);
    t = static_cast<long double>(t_mid2) +
        frac * static_cast<long double>(t_end - t_mid2);
  }

  if (t < 0) t = 0;
  if (t > (long double)std::numeric_limits<uint64_t>::max())
    t = (long double)std::numeric_limits<uint64_t>::max();
  return static_cast<uint64_t>(llround(static_cast<double>(t)));
}

// ---- STEP 4: apply 4-point piecewise map to other signals ------------------
std::vector<std::vector<std::vector<double>>> s2a2_t::step4_piecewise_bins(
  const Eigen::MatrixXd& X,
  const std::vector<uint64_t>& tp,
  const std::vector<int>& chs_idx,
  const std::vector<cycle_bound_t>& cycles,
  const s2a2_param_t& par_work,
  int nbins = 101
)
{
  std::vector<std::vector<std::vector<double>>> out;
  if (chs_idx.empty() || cycles.empty() || nbins <= 1) return out;

  out.resize(chs_idx.size());

  // Extract each channel signal and bin per cycle
  for (size_t c = 0; c < chs_idx.size(); ++c) {
    int col = chs_idx[c];
    if (col < 0 || col >= X.cols()) continue;

    std::vector<double> sig(tp.size());
    for (size_t t = 0; t < tp.size(); ++t) {
      sig[t] = X(static_cast<int>(t), col);
    }

    out[c].reserve(cycles.size());
    for (const auto& cyc : cycles) {
      uint64_t t_start = cyc.t0;
      uint64_t t_end = cyc.t1;

      // Order midpoints by cycle definition
      uint64_t t_mid1 = 0;
      uint64_t t_mid2 = 0;
      if (par_work.dir == s2a2_param_t::cross_dir_t::POS2NEG) {
        // POS->NEG crossing cycle: trough then peak
        t_mid1 = cyc.t_neg;
        t_mid2 = cyc.t_pos;
      } else { // NEG2POS
        // NEG->POS crossing cycle: peak then trough
        t_mid1 = cyc.t_pos;
        t_mid2 = cyc.t_neg;
      }

      out[c].push_back(
        bin_cycle_4pt(tp, sig, t_start, t_mid1, t_mid2, t_end, nbins)
      );
    }
  }

  return out;
}

std::vector<double> s2a2_t::mean_bins_for_signal(
  const std::vector<uint64_t>& tp,
  const std::vector<double>& sig,
  const std::vector<cycle_bound_t>& cycles,
  const s2a2_param_t& par_work,
  int nbins
)
{
  if (nbins <= 1 || cycles.empty() || sig.size() != tp.size()) {
    return std::vector<double>();
  }

  std::vector<double> mean(static_cast<size_t>(nbins), 0.0);
  std::vector<size_t> nvalid(static_cast<size_t>(nbins), 0);

  for (const auto& cyc : cycles) {
    uint64_t t_start = cyc.t0;
    uint64_t t_end = cyc.t1;
    uint64_t t_mid1 = 0;
    uint64_t t_mid2 = 0;
    if (par_work.dir == s2a2_param_t::cross_dir_t::POS2NEG) {
      t_mid1 = cyc.t_neg;
      t_mid2 = cyc.t_pos;
    } else {
      t_mid1 = cyc.t_pos;
      t_mid2 = cyc.t_neg;
    }
    std::vector<double> cyc_bins = bin_cycle_4pt(tp, sig, t_start, t_mid1, t_mid2, t_end, nbins);
    if (static_cast<int>(cyc_bins.size()) != nbins) continue;
    for (int b = 0; b < nbins; ++b) {
      double v = cyc_bins[static_cast<size_t>(b)];
      if (!finite(v)) continue;
      mean[static_cast<size_t>(b)] += v;
      nvalid[static_cast<size_t>(b)] += 1;
    }
  }

  for (int b = 0; b < nbins; ++b) {
    if (nvalid[static_cast<size_t>(b)] > 0) {
      mean[static_cast<size_t>(b)] /= static_cast<double>(nvalid[static_cast<size_t>(b)]);
    } else {
      mean[static_cast<size_t>(b)] = std::numeric_limits<double>::quiet_NaN();
    }
  }

  return mean;
}

int s2a2_t::find_max_index_in_window(
  const std::vector<uint64_t>& tp,
  const std::vector<double>& sig,
  uint64_t t_lo,
  uint64_t t_hi,
  bool use_abs
)
{
  if (tp.empty() || sig.size() != tp.size()) return -1;
  if (t_hi < t_lo) std::swap(t_lo, t_hi);
  auto it0 = std::lower_bound(tp.begin(), tp.end(), t_lo);
  auto it1 = std::upper_bound(tp.begin(), tp.end(), t_hi);
  if (it0 == tp.end()) return -1;
  size_t i0 = static_cast<size_t>(it0 - tp.begin());
  size_t i1 = (it1 == tp.begin()) ? 0 : static_cast<size_t>((it1 - tp.begin()) - 1);
  if (i0 >= sig.size() || i1 >= sig.size() || i0 > i1) return -1;

  int idx = -1;
  double vmax = -std::numeric_limits<double>::infinity();
  for (size_t i = i0; i <= i1; ++i) {
    double v = sig[i];
    if (!finite(v)) continue;
    double z = use_abs ? std::fabs(v) : v;
    if (z > vmax) { vmax = z; idx = static_cast<int>(i); }
  }
  return idx;
}

int s2a2_t::find_min_index_in_window(
  const std::vector<uint64_t>& tp,
  const std::vector<double>& sig,
  uint64_t t_lo,
  uint64_t t_hi
)
{
  if (tp.empty() || sig.size() != tp.size()) return -1;
  if (t_hi < t_lo) std::swap(t_lo, t_hi);
  auto it0 = std::lower_bound(tp.begin(), tp.end(), t_lo);
  auto it1 = std::upper_bound(tp.begin(), tp.end(), t_hi);
  if (it0 == tp.end()) return -1;
  size_t i0 = static_cast<size_t>(it0 - tp.begin());
  size_t i1 = (it1 == tp.begin()) ? 0 : static_cast<size_t>((it1 - tp.begin()) - 1);
  if (i0 >= sig.size() || i1 >= sig.size() || i0 > i1) return -1;

  int idx = -1;
  double vmin = std::numeric_limits<double>::infinity();
  for (size_t i = i0; i <= i1; ++i) {
    double v = sig[i];
    if (!finite(v)) continue;
    if (v < vmin) { vmin = v; idx = static_cast<int>(i); }
  }
  return idx;
}

std::pair<int, double> s2a2_t::crosscorr_lag_bins(
  const std::vector<double>& a,
  const std::vector<double>& b
)
{
  if (a.empty() || b.empty() || a.size() != b.size()) {
    return std::make_pair(0, std::numeric_limits<double>::quiet_NaN());
  }
  const int n = static_cast<int>(a.size());
  int best_shift = 0;
  double best_r = -std::numeric_limits<double>::infinity();

  for (int shift = 0; shift < n; ++shift) {
    double sxa = 0.0, sxb = 0.0, sxx = 0.0, syy = 0.0, sxy = 0.0;
    size_t m = 0;
    for (int i = 0; i < n; ++i) {
      int j = i + shift;
      if (j >= n) j -= n;
      double x = a[static_cast<size_t>(i)];
      double y = b[static_cast<size_t>(j)];
      if (!finite(x) || !finite(y)) continue;
      sxa += x;
      sxb += y;
      ++m;
    }
    if (m < 3) continue;
    double mean_x = sxa / static_cast<double>(m);
    double mean_y = sxb / static_cast<double>(m);
    for (int i = 0; i < n; ++i) {
      int j = i + shift;
      if (j >= n) j -= n;
      double x = a[static_cast<size_t>(i)];
      double y = b[static_cast<size_t>(j)];
      if (!finite(x) || !finite(y)) continue;
      double dx = x - mean_x;
      double dy = y - mean_y;
      sxx += dx * dx;
      syy += dy * dy;
      sxy += dx * dy;
    }
    if (sxx <= 0.0 || syy <= 0.0) continue;
    double r = sxy / std::sqrt(sxx * syy);
    if (r > best_r) {
      best_r = r;
      best_shift = shift;
    }
  }

  if (!finite(best_r)) return std::make_pair(0, best_r);
  if (best_shift > n / 2) best_shift -= n;
  return std::make_pair(best_shift, best_r);
}

int s2a2_t::find_nearest_local_extremum(
  const std::vector<uint64_t>& tp,
  const std::vector<double>& sig,
  uint64_t t_lo,
  uint64_t t_hi,
  uint64_t t_target,
  bool find_max,
  bool use_abs
)
{
  if (tp.empty() || sig.size() != tp.size()) return -1;
  if (t_hi < t_lo) std::swap(t_lo, t_hi);
  auto it0 = std::lower_bound(tp.begin(), tp.end(), t_lo);
  auto it1 = std::upper_bound(tp.begin(), tp.end(), t_hi);
  if (it0 == tp.end()) return -1;
  size_t i0 = static_cast<size_t>(it0 - tp.begin());
  size_t i1 = (it1 == tp.begin()) ? 0 : static_cast<size_t>((it1 - tp.begin()) - 1);
  if (i0 >= sig.size() || i1 >= sig.size() || i0 > i1) return -1;
  if (i1 - i0 < 2) return static_cast<int>(i0);

  int best_idx = -1;
  uint64_t best_dt = std::numeric_limits<uint64_t>::max();

  for (size_t i = i0 + 1; i + 1 <= i1; ++i) {
    double v0 = sig[i - 1];
    double v1 = sig[i];
    double v2 = sig[i + 1];
    if (!finite(v0) || !finite(v1) || !finite(v2)) continue;
    if (use_abs) { v0 = std::fabs(v0); v1 = std::fabs(v1); v2 = std::fabs(v2); }
    bool is_ext = find_max ? (v1 >= v0 && v1 >= v2) : (v1 <= v0 && v1 <= v2);
    if (!is_ext) continue;

    uint64_t t = tp[i];
    uint64_t dt = (t >= t_target) ? (t - t_target) : (t_target - t);
    if (dt < best_dt) {
      best_dt = dt;
      best_idx = static_cast<int>(i);
    }
  }

  if (best_idx >= 0) return best_idx;

  // Fallback to global extremum within window if no local extrema found.
  if (find_max) return find_max_index_in_window(tp, sig, t_lo, t_hi, use_abs);
  return find_min_index_in_window(tp, sig, t_lo, t_hi);
}

double s2a2_t::fit_sincos_amplitude(const std::vector<double>& y)
{
  const size_t n = y.size();
  if (n < 4) return std::numeric_limits<double>::quiet_NaN();

  long double S00 = 0.0L, S01 = 0.0L, S02 = 0.0L;
  long double S11 = 0.0L, S12 = 0.0L, S22 = 0.0L;
  long double b0 = 0.0L, b1 = 0.0L, b2 = 0.0L;

  for (size_t i = 0; i < n; ++i) {
    double yi = y[i];
    if (!finite(yi)) continue;
    long double tau = static_cast<long double>(i) / static_cast<long double>(n - 1);
    long double ang = 2.0L * M_PI * tau;
    long double s = std::sin(ang);
    long double c = std::cos(ang);

    S00 += 1.0L;
    S01 += s;
    S02 += c;
    S11 += s * s;
    S12 += s * c;
    S22 += c * c;
    b0 += yi;
    b1 += yi * s;
    b2 += yi * c;
  }

  long double D =
    S00 * (S11 * S22 - S12 * S12) -
    S01 * (S01 * S22 - S12 * S02) +
    S02 * (S01 * S12 - S11 * S02);
  if (D == 0.0L) return std::numeric_limits<double>::quiet_NaN();

  long double D1 =
    S00 * (b1 * S22 - S12 * b2) -
    b0 * (S01 * S22 - S12 * S02) +
    S02 * (S01 * b2 - b1 * S02);
  long double D2 =
    S00 * (S11 * b2 - b1 * S12) -
    S01 * (S01 * b2 - b1 * S02) +
    b0 * (S01 * S12 - S11 * S02);

  long double b = D1 / D;
  long double c = D2 / D;
  return std::sqrt(static_cast<double>(b * b + c * c));
}

int s2a2_t::argmax_index(const std::vector<double>& v)
{
  int idx = -1;
  double vmax = -std::numeric_limits<double>::infinity();
  for (size_t i = 0; i < v.size(); ++i) {
    if (!finite(v[i])) continue;
    if (v[i] > vmax) { vmax = v[i]; idx = static_cast<int>(i); }
  }
  return idx;
}

bool s2a2_t::refine_peak_parabolic(
  const std::vector<double>& sig,
  const std::vector<uint64_t>& tp,
  int i0,
  int i1,
  int i_peak,
  bool is_max,
  uint64_t& t_peak_out,
  double& v_peak_out
)
{
  if (i_peak < 0 || static_cast<size_t>(i_peak) >= sig.size() ||
      static_cast<size_t>(i_peak) >= tp.size()) {
    return false;
  }
  t_peak_out = tp[static_cast<size_t>(i_peak)];
  v_peak_out = sig[static_cast<size_t>(i_peak)];
  if (!finite(v_peak_out)) return false;

  if (i0 < 0 || i1 < 0 || i_peak <= i0 || i_peak >= i1) return false;
  if (i_peak - 1 < i0 || i_peak + 1 > i1) return false;
  if (i_peak <= 0 || static_cast<size_t>(i_peak + 1) >= sig.size() ||
      static_cast<size_t>(i_peak + 1) >= tp.size()) {
    return false;
  }

  const double y1 = sig[static_cast<size_t>(i_peak - 1)];
  const double y2 = sig[static_cast<size_t>(i_peak)];
  const double y3 = sig[static_cast<size_t>(i_peak + 1)];
  if (!finite(y1) || !finite(y2) || !finite(y3)) return false;

  const long double x1 = static_cast<long double>(tp[static_cast<size_t>(i_peak - 1)]);
  const long double x2 = static_cast<long double>(tp[static_cast<size_t>(i_peak)]);
  const long double x3 = static_cast<long double>(tp[static_cast<size_t>(i_peak + 1)]);
  if (!(x1 < x2 && x2 < x3)) return false;

  // Fit y(z) = A z^2 + B z + C around z = x - x2 for numerical stability.
  const long double z1 = x1 - x2;
  const long double z3 = x3 - x2;
  const long double u1 = static_cast<long double>(y1 - y2);
  const long double u3 = static_cast<long double>(y3 - y2);
  const long double den = z1 * z3 * (z1 - z3);
  if (den == 0.0L) return false;

  const long double A = (u1 * z3 - u3 * z1) / den;
  const long double B = (u3 * z1 * z1 - u1 * z3 * z3) / den;
  if (!std::isfinite(static_cast<double>(A)) ||
      !std::isfinite(static_cast<double>(B)) ||
      A == 0.0L) {
    return false;
  }
  if (is_max && A >= 0.0L) return false;
  if (!is_max && A <= 0.0L) return false;

  const long double z_v = -B / (2.0L * A);
  if (!std::isfinite(static_cast<double>(z_v))) return false;
  if (z_v < z1 || z_v > z3) return false;

  const long double x_v = x2 + z_v;
  if (x_v < x1 || x_v > x3) return false;

  const long double y_v = static_cast<long double>(y2) + A * z_v * z_v + B * z_v;
  if (!std::isfinite(static_cast<double>(y_v))) return false;

  t_peak_out = static_cast<uint64_t>(llround(static_cast<double>(x_v)));
  v_peak_out = static_cast<double>(y_v);
  return true;
}

// Find the opposite-direction crossing within (i0,i1] and return its interpolated time.
bool s2a2_t::find_opposite_cross_time(
  const std::vector<double>& seed,
  const std::vector<uint64_t>& tp,
  int i0,
  int i1,
  const s2a2_param_t& par,
  double zero,
  double h,
  uint64_t& t_cross
)
{
  auto above = [&](double v) { return v > (zero + h); };
  auto below = [&](double v) { return v < (zero - h); };

  bool armed = false;

  for (int i = i0 + 1; i <= i1; ++i) {
    double a = seed[i-1];
    double b = seed[i];
    if (!finite(a) || !finite(b)) { armed = false; continue; }
    if (!(tp[static_cast<size_t>(i)] > tp[static_cast<size_t>(i - 1)])) {
      armed = false;
      continue;
    }

    bool cross = false;
    if (!par.hysteresis) {
      if (par.dir == s2a2_param_t::cross_dir_t::POS2NEG)
        cross = (a < zero && b >= zero); // opposite: NEG2POS
      else
        cross = (a > zero && b <= zero); // opposite: POS2NEG
    } else {
      if (par.dir == s2a2_param_t::cross_dir_t::POS2NEG) {
        if (!armed && below(a)) armed = true;
        if (armed && above(b)) cross = true;
      } else { // NEG2POS
        if (!armed && above(a)) armed = true;
        if (armed && below(b)) cross = true;
      }
    }

    if (!cross) continue;

    double a0 = a - zero;
    double b0 = b - zero;
    t_cross = interp_cross_time(tp[static_cast<size_t>(i - 1)], tp[static_cast<size_t>(i)], a0, b0);
    return true;
  }

  return false;
}


// ---- STEP 1: validate + extract seed + establish baseline/hysteresis --------
void s2a2_t::step1_extract_seed(
  const Eigen::MatrixXd& X,
  const std::vector<uint64_t>* tp,
  int idx,
  const std::vector<int>& chs_idx,
  const s2a2_param_t& par_in,
  // outputs
  std::vector<double>& seed,          // length T
  double& zero,                       // baseline crossing level
  double& h,                          // hysteresis half-width (may be 0)
  uint64_t& dt_med,                   // median dt in ticks
  s2a2_param_t& par_work              // working copy (with derived params)
)
{
  if (!tp) throw std::runtime_error("tp is null");
  const size_t T = tp->size();
  if (T == 0) throw std::runtime_error("tp empty");
  if (static_cast<size_t>(X.rows()) != T) throw std::runtime_error("X.rows != tp.size");
  if (idx < 0 || idx >= X.cols()) throw std::runtime_error("idx out of range");

  par_work = par_in;

  // Extract seed samples (time-major)
  seed.resize(T);
  for (size_t t = 0; t < T; ++t) seed[t] = X(static_cast<int>(t), idx);

  // Median dt in ticks (diagnostic; also used for default min_sep)
  dt_med = median_dt_ticks(*tp);

  // Baseline (zero)
  if (par_work.use_epoch_median_zero) {
    zero = epoch_median(seed);
    if (!finite(zero)) zero = par_work.zero;
  } else {
    zero = par_work.zero;
  }

  // Hysteresis half-width
  h = par_work.h;
  if (par_work.hysteresis && h <= 0.0 && par_work.h_frac_mad > 0.0) {
    double med = (par_work.use_epoch_median_zero ? zero : epoch_median(seed));
    if (!finite(med)) med = zero;
    double mad = epoch_mad(seed, med);
    if (finite(mad)) h = par_work.h_frac_mad * mad;
  }
  if (!par_work.hysteresis) h = 0.0;

  // Default min_sep_ticks heuristic if not given
  if (par_work.min_sep_ticks == 0) {
    if (par_work.min_cycle_ticks > 0) {
      par_work.min_sep_ticks = par_work.min_cycle_ticks / 4;
    } else if (dt_med > 0) {
      par_work.min_sep_ticks = 2 * dt_med; // at least 2 samples
    }
  }
}

// ---- STEP 2: detect crossings + build cycles --------------------------------
std::vector<s2a2_t::cycle_bound_t> s2a2_t::step2_detect_cycles(
  const std::vector<double>& seed,
  const std::vector<uint64_t>& tp,
  const s2a2_param_t& par,
  double zero,
  double h,
  size_t* n_putative,
  size_t* n_after_duration
)
{
  const int T = static_cast<int>(seed.size());
  std::vector<int> xc;
  xc.reserve(static_cast<size_t>(std::max(8, T / 10)));

  std::vector<cycle_bound_t> cycles;
  size_t n_put = 0;

  auto append_cycles = [&](const std::vector<int>& xc_local) {
    if (xc_local.size() < 2) return;
    for (size_t k = 1; k < xc_local.size(); ++k) {
      int i0 = xc_local[k - 1];
      int i1 = xc_local[k];
      ++n_put;

      // duration constraints using tp at crossing indices
      uint64_t t0 = tp[static_cast<size_t>(i0)];
      uint64_t t1 = tp[static_cast<size_t>(i1)];
      if (t1 <= t0) continue;
      uint64_t dur = t1 - t0;

      if (par.min_cycle_ticks > 0 && dur < par.min_cycle_ticks) continue;
      if (par.max_cycle_ticks > 0 && dur > par.max_cycle_ticks) continue;

      // Optional: refine crossing times sub-sample using linear interpolation about zero.
      // This uses samples at (i-1,i). i0 and i1 are the right endpoint indices.
      uint64_t t0s = t0, t1s = t1;
      {
        double a0 = seed[i0 - 1] - zero;
        double b0 = seed[i0] - zero;
        t0s = interp_cross_time(tp[static_cast<size_t>(i0 - 1)],
                                tp[static_cast<size_t>(i0)],
                                a0,
                                b0);
        double a1 = seed[i1 - 1] - zero;
        double b1 = seed[i1] - zero;
        t1s = interp_cross_time(tp[static_cast<size_t>(i1 - 1)],
                                tp[static_cast<size_t>(i1)],
                                a1,
                                b1);
        if (t1s <= t0s) { t0s = t0; t1s = t1; } // fallback
      }

      cycle_bound_t cb;
      cb.i0 = i0;
      cb.i1 = i1;
      cb.t0 = t0s;
      cb.t1 = t1s;

      uint64_t t_mid = 0;
      bool has_mid = find_opposite_cross_time(seed, tp, i0, i1, par, zero, h, t_mid);

      // Optional: half-wave duration constraints (neg/pos) inside each cycle
      if (par.min_neg_ticks > 0 || par.max_neg_ticks > 0 ||
          par.min_pos_ticks > 0 || par.max_pos_ticks > 0) {
        if (!has_mid) continue; // cannot determine half-wave durations

        uint64_t neg_ticks = 0;
        uint64_t pos_ticks = 0;

        if (par.dir == s2a2_param_t::cross_dir_t::POS2NEG) {
          if (t_mid <= t0s || t1s <= t_mid) continue;
          neg_ticks = t_mid - t0s;
          pos_ticks = t1s - t_mid;
        } else { // NEG2POS
          if (t_mid <= t0s || t1s <= t_mid) continue;
          pos_ticks = t_mid - t0s;
          neg_ticks = t1s - t_mid;
        }

        if (par.min_neg_ticks > 0 && neg_ticks < par.min_neg_ticks) continue;
        if (par.max_neg_ticks > 0 && neg_ticks > par.max_neg_ticks) continue;
        if (par.min_pos_ticks > 0 && pos_ticks < par.min_pos_ticks) continue;
        if (par.max_pos_ticks > 0 && pos_ticks > par.max_pos_ticks) continue;
      }

      cb.t_mid = t_mid;
      cb.has_mid = has_mid;
      cycles.push_back(cb);
    }
  };

  auto above = [&](double v) { return v > (zero + h); };
  auto below = [&](double v) { return v < (zero - h); };

  // For hysteresis: require a "armed" state
  bool armed = false;
  uint64_t last_cross_t = 0;

  for (int i = 1; i < T; ++i) {
    double a = seed[i-1];
    double b = seed[i];
    if (!finite(a) || !finite(b)) { armed = false; continue; }
    const uint64_t ti = tp[static_cast<size_t>(i)];
    const uint64_t tim1 = tp[static_cast<size_t>(i - 1)];
    if (!(ti > tim1)) {
      armed = false;
      continue;
    }
    if (timeline_t::discontinuity(tp, static_cast<int>(par.sr_hz), i - 1, i)) {
      armed = false;
      last_cross_t = 0;
      append_cycles(xc);
      xc.clear();
      continue;
    }

    bool cross = false;

    if (!par.hysteresis) {
      if (par.dir == s2a2_param_t::cross_dir_t::POS2NEG)
        cross = (a > zero && b <= zero);
      else
        cross = (a < zero && b >= zero);
    } else {
      if (par.dir == s2a2_param_t::cross_dir_t::POS2NEG) {
        if (!armed && above(a)) armed = true;
        if (armed && below(b)) cross = true;
      } else { // NEG2POS
        if (!armed && below(a)) armed = true;
        if (armed && above(b)) cross = true;
      }
    }

    if (!cross) continue;

    // Debounce based on time separation
    uint64_t tc = tp[static_cast<size_t>(i)];
    if (!xc.empty() && par.min_sep_ticks > 0) {
      if (tc <= last_cross_t + par.min_sep_ticks) {
        continue;
      }
    }

    xc.push_back(i);
    last_cross_t = tc;
    armed = false; // require re-arming after a crossing when hysteresis is on
  }

  append_cycles(xc);
  if (n_putative) *n_putative = n_put;
  if (n_after_duration) *n_after_duration = cycles.size();

  return cycles;
}

// ---- STEP 3: mark cycles with pos/neg peaks --------------------------------
void s2a2_t::step3_mark_peaks(
  const std::vector<double>& seed,
  const std::vector<uint64_t>& tp,
  std::vector<cycle_bound_t>& cycles
)
{
  const int T = static_cast<int>(seed.size());
  for (auto& c : cycles) {
    if (c.i0 < 0 || c.i1 < 0 || c.i0 >= T || c.i1 >= T || c.i0 > c.i1) {
      continue;
    }

    double vmax = -std::numeric_limits<double>::infinity();
    double vmin = std::numeric_limits<double>::infinity();

    for (int i = c.i0; i <= c.i1; ++i) {
      double v = seed[static_cast<size_t>(i)];
      if (!finite(v)) continue;
      if (v > vmax) { vmax = v; c.i_pos = i; }
      if (v < vmin) { vmin = v; c.i_neg = i; }
    }

    if (c.i_pos >= 0) {
      c.v_pos = seed[static_cast<size_t>(c.i_pos)];
      c.t_pos = tp[static_cast<size_t>(c.i_pos)];
      uint64_t t_ref = c.t_pos;
      double v_ref = c.v_pos;
      if (refine_peak_parabolic(seed, tp, c.i0, c.i1, c.i_pos, true, t_ref, v_ref)) {
        c.t_pos = t_ref;
        c.v_pos = v_ref;
      }
    }
    if (c.i_neg >= 0) {
      c.v_neg = seed[static_cast<size_t>(c.i_neg)];
      c.t_neg = tp[static_cast<size_t>(c.i_neg)];
      uint64_t t_ref = c.t_neg;
      double v_ref = c.v_neg;
      if (refine_peak_parabolic(seed, tp, c.i0, c.i1, c.i_neg, false, t_ref, v_ref)) {
        c.t_neg = t_ref;
        c.v_neg = v_ref;
      }
    }
  }
}

// ---- STEP 3b: derive cycle metrics -----------------------------------------
void s2a2_t::step3_derive_metrics(
  const std::vector<double>& seed,
  const std::vector<uint64_t>& tp,
  const s2a2_param_t& par_work,
  double zero,
  double h,
  std::vector<cycle_bound_t>& cycles
)
{
  (void)seed;
  const size_t n = tp.size();
  for (auto& c : cycles) {
    if (c.i0 < 0 || c.i1 < 0 || c.i0 >= static_cast<int>(n) || c.i1 >= static_cast<int>(n) ||
        c.i0 >= c.i1) {
      continue;
    }

    if (c.i1 > c.i0) {
      c.rel_i_pos = (static_cast<double>(c.i_pos) - static_cast<double>(c.i0)) /
                    (static_cast<double>(c.i1) - static_cast<double>(c.i0));
      c.rel_i_neg = (static_cast<double>(c.i_neg) - static_cast<double>(c.i0)) /
                    (static_cast<double>(c.i1) - static_cast<double>(c.i0));
    }
    if (c.t1 > c.t0) {
      c.rel_pos = (static_cast<double>(c.t_pos) - static_cast<double>(c.t0)) /
                  (static_cast<double>(c.t1) - static_cast<double>(c.t0));
      c.rel_neg = (static_cast<double>(c.t_neg) - static_cast<double>(c.t0)) /
                  (static_cast<double>(c.t1) - static_cast<double>(c.t0));
    }

    if (c.t_pos >= c.t0) c.dt_pos_s = (c.t_pos - c.t0) * globals::tp_duration;
    if (c.t_neg >= c.t0) c.dt_neg_s = (c.t_neg - c.t0) * globals::tp_duration;

    uint64_t t_mid = 0;
    bool has_mid = find_opposite_cross_time(seed, tp, c.i0, c.i1, par_work, zero, h, t_mid);
    if (!has_mid) continue;

    if (par_work.dir == s2a2_param_t::cross_dir_t::POS2NEG) {
      double dtp = (c.t_pos > t_mid) ? (c.t_pos - t_mid) * globals::tp_duration : 0.0;
      double dtn = (c.t_neg > c.t0) ? (c.t_neg - c.t0) * globals::tp_duration : 0.0;
      if (dtp > 0.0) {
        c.pos_slope = (c.v_pos - zero) / dtp;
        c.pos_slope_norm = 1.0 / dtp;
      }
      if (dtn > 0.0) {
        c.neg_slope = (c.v_neg - zero) / dtn;
        c.neg_slope_norm = 1.0 / dtn;
      }
    } else { // NEG2POS
      double dtp = (c.t_pos > c.t0) ? (c.t_pos - c.t0) * globals::tp_duration : 0.0;
      double dtn = (c.t_neg > t_mid) ? (c.t_neg - t_mid) * globals::tp_duration : 0.0;
      if (dtp > 0.0) {
        c.pos_slope = (c.v_pos - zero) / dtp;
        c.pos_slope_norm = 1.0 / dtp;
      }
      if (dtn > 0.0) {
        c.neg_slope = (c.v_neg - zero) / dtn;
        c.neg_slope_norm = 1.0 / dtn;
      }
    }
  }
}




s2a2_t::s2a2_out_t s2a2_t::s2a2_proc(
  const Eigen::MatrixXd& X,
  const std::vector<uint64_t>* tp,
  int idx,
  const std::vector<int>& chs_idx,
  const s2a2_param_t& par,
  const std::string& seg_label,
  const std::vector<std::string>& sig_labels )
{
  
  s2a2_out_t res;
  
  if (tp == nullptr) Helper::halt("internal error");
  if (tp->size() != static_cast<size_t>(X.rows())) Helper::halt("internal error");

  // start output for this seed
  writer.level(seg_label, "SEED");
  
  std::vector<double> seed;
  double zero = 0.0, h = 0.0;
  uint64_t dt_med = 0;
  s2a2_param_t par_work;

  step1_extract_seed(X, tp, idx, chs_idx, par, seed, zero, h, dt_med, par_work);
  
  size_t n_putative = 0;
  size_t n_after_dur = 0;
  std::vector<cycle_bound_t> cycles =
    step2_detect_cycles(seed, *tp, par_work, zero, h, &n_putative, &n_after_dur);

  step3_mark_peaks(seed, *tp, cycles);

  step3_derive_metrics(seed, *tp, par_work, zero, h, cycles);

  // Keep pre-magnitude-filter cycles for seed-level summaries so SEED and SEED,BIN are consistent.
  const std::vector<cycle_bound_t> cycles_seed = cycles;

  {
    std::vector<double> v_tpos_tneg;
    std::vector<double> v_tpos;
    std::vector<double> v_tneg;
    std::vector<double> v_ttot;
    std::vector<double> v_zc_neg;
    std::vector<double> v_neg_zc;
    std::vector<double> v_zc_pos;
    std::vector<double> v_pos_zc;
    std::vector<double> v_amp_pos;
    std::vector<double> v_amp_neg;
    std::vector<double> v_amp_p2p;
    std::vector<double> v_slope_pos;
    std::vector<double> v_slope_neg;
    std::vector<double> v_sharp_pos;
    std::vector<double> v_sharp_neg;
    std::vector<double> v_amp_asym;
    std::vector<double> v_slope_asym;
    std::vector<double> v_sharp_asym;
    const bool emit_cycle = par_work.emit_cycle_metrics;

    for (size_t ci = 0; ci < cycles.size(); ++ci) {
      const auto& c = cycles[ci];
      if (c.i0 < 0 || c.i1 < 0 || c.i0 >= c.i1) continue;
      uint64_t t_mid = 0;
      if (!find_opposite_cross_time(seed, *tp, c.i0, c.i1, par_work, zero, h, t_mid)) {
        continue;
      }

      double Tpos = std::numeric_limits<double>::quiet_NaN();
      double Tneg = std::numeric_limits<double>::quiet_NaN();
      if (par_work.dir == s2a2_param_t::cross_dir_t::POS2NEG) {
        if (t_mid > c.t0) Tneg = (t_mid - c.t0) * globals::tp_duration;
        if (c.t1 > t_mid) Tpos = (c.t1 - t_mid) * globals::tp_duration;
      } else {
        if (t_mid > c.t0) Tpos = (t_mid - c.t0) * globals::tp_duration;
        if (c.t1 > t_mid) Tneg = (c.t1 - t_mid) * globals::tp_duration;
      }

      double ratio = std::numeric_limits<double>::quiet_NaN();
      if (finite(Tpos) && finite(Tneg) && Tneg > 0.0) ratio = Tpos / Tneg;

      double zc_neg = std::numeric_limits<double>::quiet_NaN();
      double neg_zc = std::numeric_limits<double>::quiet_NaN();
      double zc_pos = std::numeric_limits<double>::quiet_NaN();
      double pos_zc = std::numeric_limits<double>::quiet_NaN();
      if (par_work.dir == s2a2_param_t::cross_dir_t::POS2NEG) {
        if (c.t_neg > c.t0) zc_neg = (c.t_neg - c.t0) * globals::tp_duration;
        if (t_mid > c.t_neg) neg_zc = (t_mid - c.t_neg) * globals::tp_duration;
        if (c.t_pos > t_mid) zc_pos = (c.t_pos - t_mid) * globals::tp_duration;
        if (c.t1 > c.t_pos) pos_zc = (c.t1 - c.t_pos) * globals::tp_duration;
      } else {
        if (c.t_pos > c.t0) zc_pos = (c.t_pos - c.t0) * globals::tp_duration;
        if (t_mid > c.t_pos) pos_zc = (t_mid - c.t_pos) * globals::tp_duration;
        if (c.t_neg > t_mid) zc_neg = (c.t_neg - t_mid) * globals::tp_duration;
        if (c.t1 > c.t_neg) neg_zc = (c.t1 - c.t_neg) * globals::tp_duration;
      }

      double Apos = c.v_pos - zero;
      double Aneg = c.v_neg - zero;
      double absAneg = std::fabs(Aneg);
      double Ap2p = std::numeric_limits<double>::quiet_NaN();
      if (finite(c.v_pos) && finite(c.v_neg)) Ap2p = c.v_pos - c.v_neg;

      double slope_pos = std::numeric_limits<double>::quiet_NaN();
      if (finite(c.pos_slope)) slope_pos = std::fabs(c.pos_slope);
      double slope_neg = std::numeric_limits<double>::quiet_NaN();
      if (finite(c.neg_slope)) slope_neg = std::fabs(c.neg_slope);

      double sharp_pos = std::numeric_limits<double>::quiet_NaN();
      double sharp_neg = std::numeric_limits<double>::quiet_NaN();
      if (c.i_pos >= 0 && c.i_neg >= 0) {
        double level_pos = zero + 0.5 * (c.v_pos - zero);
        double level_neg = zero + 0.5 * (c.v_neg - zero);
        double wpos = half_width_s(seed, *tp, c.i0, c.i1, c.i_pos, level_pos);
        double wneg = half_width_s(seed, *tp, c.i0, c.i1, c.i_neg, level_neg);
        if (finite(wpos) && wpos > 0.0) sharp_pos = std::fabs(c.v_pos - zero) / wpos;
        if (finite(wneg) && wneg > 0.0) sharp_neg = std::fabs(c.v_neg - zero) / wneg;
      }

      double amp_asym = std::numeric_limits<double>::quiet_NaN();
      if (finite(Apos) && finite(absAneg)) {
        double denom = Apos + absAneg;
        if (denom != 0.0) amp_asym = (Apos - absAneg) / denom;
      }

      double slope_asym = std::numeric_limits<double>::quiet_NaN();
      if (finite(slope_pos) && finite(slope_neg) && slope_neg > 0.0) {
        slope_asym = slope_pos / slope_neg;
      }

      double sharp_asym = std::numeric_limits<double>::quiet_NaN();
      if (finite(sharp_pos) && finite(sharp_neg)) {
        double denom = sharp_pos + sharp_neg;
        if (denom != 0.0) sharp_asym = (sharp_pos - sharp_neg) / denom;
      }

      v_tpos_tneg.push_back(ratio);
      v_tpos.push_back(Tpos);
      v_tneg.push_back(Tneg);
      if (finite(Tpos) && finite(Tneg)) v_ttot.push_back(Tpos + Tneg);
      v_zc_neg.push_back(zc_neg);
      v_neg_zc.push_back(neg_zc);
      v_zc_pos.push_back(zc_pos);
      v_pos_zc.push_back(pos_zc);
      v_amp_pos.push_back(Apos);
      v_amp_neg.push_back(absAneg);
      v_amp_p2p.push_back(Ap2p);
      v_slope_pos.push_back(slope_pos);
      v_slope_neg.push_back(slope_neg);
      v_sharp_pos.push_back(sharp_pos);
      v_sharp_neg.push_back(sharp_neg);
      v_amp_asym.push_back(amp_asym);
      v_slope_asym.push_back(slope_asym);
      v_sharp_asym.push_back(sharp_asym);

      if (emit_cycle) {
        writer.level( static_cast<int>(ci), globals::count_strat );
        writer.value("T0_S", c.t0 * globals::tp_duration);
        writer.value("T1_S", c.t1 * globals::tp_duration);
        writer.value("DUR_POS", Tpos);
        writer.value("DUR_NEG", Tneg);
        writer.value("DUR_ZC_NEG", zc_neg);
        writer.value("DUR_NEG_ZC", neg_zc);
        writer.value("DUR_ZC_POS", zc_pos);
        writer.value("DUR_POS_ZC", pos_zc);
        writer.value("AMP_POS", Apos);
        writer.value("AMP_NEG", absAneg);
        writer.value("AMP_P2P", Ap2p);
        writer.value("SLOPE_POS", slope_pos);
        writer.value("SLOPE_NEG", slope_neg);
        writer.value("SHARP_POS", sharp_pos);
        writer.value("SHARP_NEG", sharp_neg);
        if (par_work.emit_asym_metrics) {
          writer.value("DUR_RATIO", ratio);
          writer.value("AMP_ASYM", amp_asym);
          writer.value("SLOPE_ASYM", slope_asym);
          writer.value("SHARP_ASYM", sharp_asym);
        }
        writer.unlevel( globals::count_strat );
      }
    }

    if (par_work.emit_seed_summary) {
      writer.value("DUR_POS", mean_dbl(v_tpos));
      writer.value("DUR_POS_MD", median_dbl(v_tpos));
      writer.value("DUR_NEG", mean_dbl(v_tneg));
      writer.value("DUR_NEG_MD", median_dbl(v_tneg));
      writer.value("DUR", mean_dbl(v_ttot));
      writer.value("DUR_MD", median_dbl(v_ttot));
      writer.value("DUR_ZC_NEG", mean_dbl(v_zc_neg));
      writer.value("DUR_ZC_NEG_MD", median_dbl(v_zc_neg));
      writer.value("DUR_NEG_ZC", mean_dbl(v_neg_zc));
      writer.value("DUR_NEG_ZC_MD", median_dbl(v_neg_zc));
      writer.value("DUR_ZC_POS", mean_dbl(v_zc_pos));
      writer.value("DUR_ZC_POS_MD", median_dbl(v_zc_pos));
      writer.value("DUR_POS_ZC", mean_dbl(v_pos_zc));
      writer.value("DUR_POS_ZC_MD", median_dbl(v_pos_zc));
      writer.value("AMP_POS", mean_dbl(v_amp_pos));
      writer.value("AMP_POS_MD", median_dbl(v_amp_pos));
      writer.value("AMP_NEG", mean_dbl(v_amp_neg));
      writer.value("AMP_NEG_MD", median_dbl(v_amp_neg));
      writer.value("AMP_P2P", mean_dbl(v_amp_p2p));
      writer.value("AMP_P2P_MD", median_dbl(v_amp_p2p));
      writer.value("SLOPE_POS", mean_dbl(v_slope_pos));
      writer.value("SLOPE_POS_MD", median_dbl(v_slope_pos));
      writer.value("SLOPE_NEG", mean_dbl(v_slope_neg));
      writer.value("SLOPE_NEG_MD", median_dbl(v_slope_neg));
      writer.value("SHARP_POS", mean_dbl(v_sharp_pos));
      writer.value("SHARP_POS_MD", median_dbl(v_sharp_pos));
      writer.value("SHARP_NEG", mean_dbl(v_sharp_neg));
      writer.value("SHARP_NEG_MD", median_dbl(v_sharp_neg));
      if (par_work.emit_asym_metrics) {
        writer.value("DUR_RATIO", mean_dbl(v_tpos_tneg));
        writer.value("DUR_RATIO_MD", median_dbl(v_tpos_tneg));
        writer.value("AMP_ASYM", mean_dbl(v_amp_asym));
        writer.value("AMP_ASYM_MD", median_dbl(v_amp_asym));
        writer.value("SLOPE_ASYM", mean_dbl(v_slope_asym));
        writer.value("SLOPE_ASYM_MD", median_dbl(v_slope_asym));
        writer.value("SHARP_ASYM", mean_dbl(v_sharp_asym));
        writer.value("SHARP_ASYM_MD", median_dbl(v_sharp_asym));
      }
      if (par_work.emit_mad) {
        writer.value("DUR_POS_MAD", mad_dbl(v_tpos, median_dbl(v_tpos)));
        writer.value("DUR_NEG_MAD", mad_dbl(v_tneg, median_dbl(v_tneg)));
        writer.value("DUR_MAD", mad_dbl(v_ttot, median_dbl(v_ttot)));
        writer.value("DUR_ZC_NEG_MAD", mad_dbl(v_zc_neg, median_dbl(v_zc_neg)));
        writer.value("DUR_NEG_ZC_MAD", mad_dbl(v_neg_zc, median_dbl(v_neg_zc)));
        writer.value("DUR_ZC_POS_MAD", mad_dbl(v_zc_pos, median_dbl(v_zc_pos)));
        writer.value("DUR_POS_ZC_MAD", mad_dbl(v_pos_zc, median_dbl(v_pos_zc)));
        writer.value("AMP_POS_MAD", mad_dbl(v_amp_pos, median_dbl(v_amp_pos)));
        writer.value("AMP_NEG_MAD", mad_dbl(v_amp_neg, median_dbl(v_amp_neg)));
        writer.value("AMP_P2P_MAD", mad_dbl(v_amp_p2p, median_dbl(v_amp_p2p)));
        writer.value("SLOPE_POS_MAD", mad_dbl(v_slope_pos, median_dbl(v_slope_pos)));
        writer.value("SLOPE_NEG_MAD", mad_dbl(v_slope_neg, median_dbl(v_slope_neg)));
        writer.value("SHARP_POS_MAD", mad_dbl(v_sharp_pos, median_dbl(v_sharp_pos)));
        writer.value("SHARP_NEG_MAD", mad_dbl(v_sharp_neg, median_dbl(v_sharp_neg)));
        if (par_work.emit_asym_metrics) {
          writer.value("DUR_RATIO_MAD", mad_dbl(v_tpos_tneg, median_dbl(v_tpos_tneg)));
          writer.value("AMP_ASYM_MAD", mad_dbl(v_amp_asym, median_dbl(v_amp_asym)));
          writer.value("SLOPE_ASYM_MAD", mad_dbl(v_slope_asym, median_dbl(v_slope_asym)));
          writer.value("SHARP_ASYM_MAD", mad_dbl(v_sharp_asym, median_dbl(v_sharp_asym)));
        }
      }
    }

  }

  int mag_excl = 0;
  const size_t pre_mag_cycles = cycles.size();
  std::vector<double> mag_vals;
  cycles = filter_cycles_by_mag(seed, par_work, cycles, &mag_excl, &mag_vals);
  if (par_work.use_mag) {
    logger << "  seed=" << seg_label
           << " cycles: putative=" << n_putative
           << " after_dur=" << n_after_dur
           << " after_mag=" << cycles.size()
           << "\n";
  }
 
  if (par_work.emit_seed_summary) {
    writer.value("N_PRE", static_cast<int>(pre_mag_cycles));
    writer.value("N_POST", static_cast<int>(cycles.size()));
  }
  
  res.nbins = 101;

  res.bins = step4_piecewise_bins(X, *tp, chs_idx, cycles, par_work, res.nbins);

  const int phase_bins12 = 12;
  std::vector<double> seed_ph12_mean(static_cast<size_t>(phase_bins12),
                                      std::numeric_limits<double>::quiet_NaN());
  if (par_work.emit_seed_summary && !cycles_seed.empty()) {
    auto seed_interval_mean = [&](uint64_t t_lo, uint64_t t_hi) -> double {
      if (t_hi <= t_lo || tp->empty() || seed.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
      }
      auto it0 = std::lower_bound(tp->begin(), tp->end(), t_lo);
      auto it1 = std::lower_bound(tp->begin(), tp->end(), t_hi);
      double ssum = 0.0;
      size_t n = 0;
      for (auto it = it0; it != it1; ++it) {
        size_t i = static_cast<size_t>(it - tp->begin());
        if (i >= seed.size()) break;
        double v = seed[i];
        if (!finite(v)) continue;
        ssum += v;
        ++n;
      }
      if (n > 0) return ssum / static_cast<double>(n);
      long double t_mid = 0.5L * (static_cast<long double>(t_lo) + static_cast<long double>(t_hi));
      double v_mid = interp_value_at_tp_ld(*tp, seed, t_mid);
      return finite(v_mid) ? v_mid : std::numeric_limits<double>::quiet_NaN();
    };

    std::vector<double> ph12_sum(static_cast<size_t>(phase_bins12), 0.0);
    std::vector<size_t> ph12_n(static_cast<size_t>(phase_bins12), 0);

    for (const auto& cyc : cycles_seed) {
      uint64_t t0 = cyc.t0;
      uint64_t t1 = cyc.t1;
      uint64_t t_mid1 = 0;
      uint64_t t_mid2 = 0;
      if (par_work.dir == s2a2_param_t::cross_dir_t::POS2NEG) {
        t_mid1 = cyc.t_neg;
        t_mid2 = cyc.t_pos;
      } else {
        t_mid1 = cyc.t_pos;
        t_mid2 = cyc.t_neg;
      }
      if (!(t1 > t0 && t0 < t_mid1 && t_mid1 < t_mid2 && t_mid2 < t1)) continue;

      for (int b = 0; b < phase_bins12; ++b) {
        long double f0 = static_cast<long double>(b) / 12.0L;
        long double f1 = static_cast<long double>(b + 1) / 12.0L;
        uint64_t tb0 = phase_to_time4pt(t0, t_mid1, t_mid2, t1, f0);
        uint64_t tb1 = (b == phase_bins12 - 1)
          ? t1
          : phase_to_time4pt(t0, t_mid1, t_mid2, t1, f1);
        if (tb1 <= tb0) continue;

        double p = seed_interval_mean(tb0, tb1);
        if (!finite(p)) continue;
        ph12_sum[static_cast<size_t>(b)] += p;
        ph12_n[static_cast<size_t>(b)] += 1;
      }
    }

    for (int b = 0; b < phase_bins12; ++b) {
      if (ph12_n[static_cast<size_t>(b)] > 0) {
        seed_ph12_mean[static_cast<size_t>(b)] =
          ph12_sum[static_cast<size_t>(b)] /
          static_cast<double>(ph12_n[static_cast<size_t>(b)]);
      }
    }

    for (int b = 0; b < phase_bins12; ++b) {
      const std::string bin_label =
        std::string("B") + (b < 9 ? "0" : "") + Helper::int2str(b + 1);
      writer.level(bin_label, "BIN");
      writer.value("MEAN", seed_ph12_mean[static_cast<size_t>(b)]);
      writer.unlevel("BIN");
    }
  }

  std::vector<std::vector<double>> mean_bins;

  std::vector<double> mean_seed;

  if (!res.bins.empty() && res.nbins > 1) {
    mean_bins.assign(res.bins.size(), std::vector<double>(static_cast<size_t>(res.nbins),
                                                         std::numeric_limits<double>::quiet_NaN()));
    int seed_col = -1;
    if (!seg_label.empty()) {
      for (size_t c = 0; c < sig_labels.size(); ++c) {
        if (sig_labels[c] == seg_label) { seed_col = static_cast<int>(c); break; }
      }
    }

    // First pass: compute mean waveform per channel
    for (size_t c = 0; c < res.bins.size(); ++c) {
      std::vector<double> mean(static_cast<size_t>(res.nbins), 0.0);
      std::vector<size_t> nvalid(static_cast<size_t>(res.nbins), 0);

      for (const auto& cyc : res.bins[c]) {
        if (static_cast<int>(cyc.size()) != res.nbins) continue;
        for (int b = 0; b < res.nbins; ++b) {
          double v = cyc[static_cast<size_t>(b)];
          if (!finite(v)) continue;
          mean[static_cast<size_t>(b)] += v;
          nvalid[static_cast<size_t>(b)] += 1;
        }
      }

      for (int b = 0; b < res.nbins; ++b) {
        if (nvalid[static_cast<size_t>(b)] > 0) {
          mean[static_cast<size_t>(b)] /= static_cast<double>(nvalid[static_cast<size_t>(b)]);
        } else {
          mean[static_cast<size_t>(b)] = std::numeric_limits<double>::quiet_NaN();
        }
      }
      mean_bins[c] = mean;
    }

    if (seed_col >= 0 && static_cast<size_t>(seed_col) < mean_bins.size()) {
      mean_seed = mean_bins[static_cast<size_t>(seed_col)];
    } else {
      mean_seed = mean_bins_for_signal(*tp, seed, cycles, par_work, res.nbins);
    }

    // Second pass: compute stats and per-channel outputs
    for (size_t c = 0; c < res.bins.size(); ++c) {
      const std::vector<double>& mean = mean_bins[c];
      // summary stats
      double D = std::numeric_limits<double>::quiet_NaN();
      {
        double vmin = std::numeric_limits<double>::infinity();
        double vmax = -std::numeric_limits<double>::infinity();
        for (double v : mean) {
          if (!finite(v)) continue;
          if (v < vmin) vmin = v;
          if (v > vmax) vmax = v;
        }
        if (finite(vmin) && finite(vmax)) D = vmax - vmin;
      }

      int idx_max = argmax_index(mean);
      double tau_max_deg = std::numeric_limits<double>::quiet_NaN();
      if (idx_max >= 0 && res.nbins > 1) {
        tau_max_deg = 360.0 * static_cast<double>(idx_max) /
                      static_cast<double>(res.nbins - 1);
      }

      double M = fit_sincos_amplitude(mean);

      // per-cycle variability for D and tau_max
      std::vector<double> Dk;
      std::vector<double> tauk_deg;
      Dk.reserve(res.bins[c].size());
      tauk_deg.reserve(res.bins[c].size());

      for (const auto& cyc : res.bins[c]) {
        if (static_cast<int>(cyc.size()) != res.nbins) continue;
        double vmin = std::numeric_limits<double>::infinity();
        double vmax = -std::numeric_limits<double>::infinity();
        for (double v : cyc) {
          if (!finite(v)) continue;
          if (v < vmin) vmin = v;
          if (v > vmax) vmax = v;
        }
        if (finite(vmin) && finite(vmax)) Dk.push_back(vmax - vmin);

        int im = argmax_index(cyc);
        if (im >= 0 && res.nbins > 1) {
          double tau_deg = 360.0 * static_cast<double>(im) /
                           static_cast<double>(res.nbins - 1);
          tauk_deg.push_back(tau_deg);
        }
      }

      double D_sd = sd_dbl(Dk);

      // circular dispersion from per-cycle tau max
      double circ_disp = std::numeric_limits<double>::quiet_NaN();
      double circ_mean_deg = std::numeric_limits<double>::quiet_NaN();
      double circ_R = std::numeric_limits<double>::quiet_NaN();
      if (!tauk_deg.empty()) {
        double sx = 0.0;
        double sy = 0.0;
        for (double deg : tauk_deg) {
          double rad = deg * M_PI / 180.0;
          sx += std::cos(rad);
          sy += std::sin(rad);
        }
        double R = std::sqrt(sx * sx + sy * sy) /
                   static_cast<double>(tauk_deg.size());
        circ_R = R;
        circ_disp = 1.0 - R;
        double mean_rad = std::atan2(sy, sx);
        double mean_deg = mean_rad * 180.0 / M_PI;
        if (mean_deg < 0.0) mean_deg += 360.0;
        circ_mean_deg = mean_deg;
      }

      // per-cycle absolute-time lag to seed peaks
      std::vector<double> dt_pos_s;
      std::vector<double> dt_neg_s;
      std::vector<double> td_pos_p2p;
      std::vector<double> td_neg_p2p;
      std::vector<const cycle_bound_t*> td_cycles;
      std::vector<double> sig;
      bool has_sig = false;
      {
        // build signal vector for this channel
        int col = (c < chs_idx.size()) ? chs_idx[c] : -1;
        if (col >= 0 && col < X.cols()) {
          sig.resize(tp->size());
          for (size_t t = 0; t < tp->size(); ++t) sig[t] = X(static_cast<int>(t), col);
          has_sig = true;
          const bool td_pos = par_work.emit_time_domain && (par_work.time_lock == "pos");
          const bool td_neg = par_work.emit_time_domain && (par_work.time_lock == "neg");
          if (par_work.emit_time_domain && par_work.time_window_s > 0.0) {
            td_cycles.reserve(cycles.size());
          }
          uint64_t half_ticks = 0;
          if (par_work.time_window_s > 0.0) {
            half_ticks = static_cast<uint64_t>(
              llround(0.5 * par_work.time_window_s / globals::tp_duration));
          }
          uint64_t w_ticks = 0;
          if (par_work.lag_window_s > 0.0) {
            w_ticks = static_cast<uint64_t>(llround(par_work.lag_window_s / globals::tp_duration));
          }
          for (const auto& cyc : cycles) {
            uint64_t t0 = cyc.t0;
            uint64_t t1 = cyc.t1;
            uint64_t tpos = cyc.t_pos;
            uint64_t tneg = cyc.t_neg;

            uint64_t lo = t0;
            uint64_t hi = t1;
            if (w_ticks > 0) {
              lo = (tpos > w_ticks ? tpos - w_ticks : 0);
              hi = tpos + w_ticks;
            }
            int imax = find_nearest_local_extremum(
              *tp, sig, lo, hi, tpos, true, par_work.lag_use_abs
            );
            if (imax >= 0) {
              double dt = (static_cast<double>((*tp)[static_cast<size_t>(imax)]) -
                           static_cast<double>(tpos)) * globals::tp_duration;
              dt_pos_s.push_back(dt);
            }
            if (td_pos && par_work.time_window_s > 0.0) {
              uint64_t t_lo = (tpos > half_ticks) ? (tpos - half_ticks) : 0;
              uint64_t t_hi = tpos + half_ticks;
              if (window_range_ok(*tp, static_cast<int>(par_work.sr_hz), t_lo, t_hi, nullptr, nullptr)) {
                td_cycles.push_back(&cyc);
                window_metrics_t wm = window_metrics(sig, *tp, tpos, par_work.time_window_s, zero,
                                                    static_cast<int>(par_work.sr_hz));
              if (finite(wm.p2p)) td_pos_p2p.push_back(wm.p2p);
              }
            }

            lo = t0;
            hi = t1;
            if (w_ticks > 0) {
              lo = (tneg > w_ticks ? tneg - w_ticks : 0);
              hi = tneg + w_ticks;
            }
            int imin = find_nearest_local_extremum(
              *tp, sig, lo, hi, tneg, false, false
            );
            if (imin >= 0) {
              double dt = (static_cast<double>((*tp)[static_cast<size_t>(imin)]) -
                           static_cast<double>(tneg)) * globals::tp_duration;
              dt_neg_s.push_back(dt);
            }
            if (td_neg && par_work.time_window_s > 0.0) {
              uint64_t t_lo = (tneg > half_ticks) ? (tneg - half_ticks) : 0;
              uint64_t t_hi = tneg + half_ticks;
              if (window_range_ok(*tp, static_cast<int>(par_work.sr_hz), t_lo, t_hi, nullptr, nullptr)) {
                td_cycles.push_back(&cyc);
                window_metrics_t wm = window_metrics(sig, *tp, tneg, par_work.time_window_s, zero,
                                                    static_cast<int>(par_work.sr_hz));
              if (finite(wm.p2p)) td_neg_p2p.push_back(wm.p2p);
              }
            }
          }
        }
      }

      std::vector<double> ph12_mean(static_cast<size_t>(phase_bins12),
                                     std::numeric_limits<double>::quiet_NaN());
      if (has_sig && !cycles.empty()) {
        auto interval_mean = [&](uint64_t t_lo, uint64_t t_hi) -> double {
          if (t_hi <= t_lo || tp->empty() || sig.empty()) {
            return std::numeric_limits<double>::quiet_NaN();
          }
          auto it0 = std::lower_bound(tp->begin(), tp->end(), t_lo);
          auto it1 = std::lower_bound(tp->begin(), tp->end(), t_hi);
          double ssum = 0.0;
          size_t n = 0;
          for (auto it = it0; it != it1; ++it) {
            size_t i = static_cast<size_t>(it - tp->begin());
            if (i >= sig.size()) break;
            double v = sig[i];
            if (!finite(v)) continue;
            ssum += v;
            ++n;
          }
          if (n > 0) return ssum / static_cast<double>(n);
          long double t_mid = 0.5L * (static_cast<long double>(t_lo) + static_cast<long double>(t_hi));
          double v_mid = interp_value_at_tp_ld(*tp, sig, t_mid);
          return finite(v_mid) ? v_mid : std::numeric_limits<double>::quiet_NaN();
        };

        std::vector<double> ph12_sum(static_cast<size_t>(phase_bins12), 0.0);
        std::vector<size_t> ph12_n(static_cast<size_t>(phase_bins12), 0);
        for (const auto& cyc : cycles) {
          uint64_t t0 = cyc.t0;
          uint64_t t1 = cyc.t1;
          uint64_t t_mid1 = 0;
          uint64_t t_mid2 = 0;
          if (par_work.dir == s2a2_param_t::cross_dir_t::POS2NEG) {
            t_mid1 = cyc.t_neg;
            t_mid2 = cyc.t_pos;
          } else {
            t_mid1 = cyc.t_pos;
            t_mid2 = cyc.t_neg;
          }
          if (!(t1 > t0 && t0 < t_mid1 && t_mid1 < t_mid2 && t_mid2 < t1)) continue;

          for (int b = 0; b < phase_bins12; ++b) {
            long double f0 = static_cast<long double>(b) / 12.0L;
            long double f1 = static_cast<long double>(b + 1) / 12.0L;
            uint64_t tb0 = phase_to_time4pt(t0, t_mid1, t_mid2, t1, f0);
            uint64_t tb1 = (b == phase_bins12 - 1)
              ? t1
              : phase_to_time4pt(t0, t_mid1, t_mid2, t1, f1);
            if (tb1 <= tb0) continue;

            double p = interval_mean(tb0, tb1);
            if (!finite(p)) continue;
            ph12_sum[static_cast<size_t>(b)] += p;
            ph12_n[static_cast<size_t>(b)] += 1;
          }
        }
        for (int b = 0; b < phase_bins12; ++b) {
          if (ph12_n[static_cast<size_t>(b)] > 0) {
            ph12_mean[static_cast<size_t>(b)] =
              ph12_sum[static_cast<size_t>(b)] /
              static_cast<double>(ph12_n[static_cast<size_t>(b)]);
          }
        }
      }

      double dt_pos_med = median_dbl(dt_pos_s);
      double dt_pos_mad = mad_dbl(dt_pos_s, dt_pos_med);
      double dt_pos_mean = mean_dbl(dt_pos_s);
      double dt_neg_med = median_dbl(dt_neg_s);
      double dt_neg_mad = mad_dbl(dt_neg_s, dt_neg_med);
      double dt_neg_mean = mean_dbl(dt_neg_s);
      const bool use_pos = true;
      const bool use_neg = true;
      const bool td_pos = par_work.emit_time_domain && (par_work.time_lock == "pos");
      const bool td_neg = par_work.emit_time_domain && (par_work.time_lock == "neg");
      double td_pos_p2p_mean = mean_dbl(td_pos_p2p);
      double td_pos_p2p_med = median_dbl(td_pos_p2p);
      double td_pos_p2p_mad = mad_dbl(td_pos_p2p, td_pos_p2p_med);
      double td_neg_p2p_mean = mean_dbl(td_neg_p2p);
      double td_neg_p2p_med = median_dbl(td_neg_p2p);
      double td_neg_p2p_mad = mad_dbl(td_neg_p2p, td_neg_p2p_med);

      const std::string seg_lab = seg_label.empty() ? "." : seg_label;
      const std::string sig_lab =
        (c < sig_labels.size() && !sig_labels[c].empty()) ? sig_labels[c] : ".";
      int cc_lag_bins = 0;
      double cc_r = std::numeric_limits<double>::quiet_NaN();
      if (!mean_seed.empty() && c < mean_bins.size()) {
        std::pair<int, double> cc = crosscorr_lag_bins(mean_seed, mean_bins[c]);
        cc_lag_bins = cc.first;
        cc_r = cc.second;
      }
      double cc_lag_s = std::numeric_limits<double>::quiet_NaN();
      if (!cycles.empty()) {
        std::vector<double> cyc_len_s;
        cyc_len_s.reserve(cycles.size());
        for (const auto& cyc : cycles) {
          if (cyc.t1 > cyc.t0) {
            cyc_len_s.push_back((cyc.t1 - cyc.t0) * globals::tp_duration);
          }
        }
        double med_len = median_dbl(cyc_len_s);
        if (finite(med_len) && res.nbins > 1) {
          cc_lag_s = med_len * static_cast<double>(cc_lag_bins) /
                     static_cast<double>(res.nbins - 1);
        }
      }

      if (par_work.emit_sig_summary) {
        writer.level(sig_lab, "SIG");
        writer.value("D", D);
        writer.value("TAU_MAX_DEG", tau_max_deg);
        double M_over_D = std::numeric_limits<double>::quiet_NaN();
        if (finite(D) && D != 0.0 && finite(M)) M_over_D = M / D;
        writer.value("SINFIT_OVER_P2P", M_over_D);
        writer.value("D_SD", D_sd);
        writer.value("CIRC_MEAN_DEG", circ_mean_deg);
        writer.value("R", circ_R);
        writer.value("CC_LAG", cc_lag_s);
        writer.value("CC_R", cc_r);

        if (par_work.emit_time_domain && par_work.time_window_s > 0.0 && par_work.emit_td_summary) {
          const size_t td_cycles_total = cycles.size();
          const size_t td_cycles_used = td_cycles.size();
          const double td_cycles_used_pct = (td_cycles_total > 0)
            ? static_cast<double>(td_cycles_used) / static_cast<double>(td_cycles_total)
            : std::numeric_limits<double>::quiet_NaN();
          writer.value("TD_USED_PCT", td_cycles_used_pct);
        }

        if (use_pos || use_neg) {
          const bool use_pos_metrics = (par_work.time_lock == "pos");
          const double dt_mean = use_pos_metrics ? dt_pos_mean : dt_neg_mean;
          const double dt_md = use_pos_metrics ? dt_pos_med : dt_neg_med;
          const double dt_mad = use_pos_metrics ? dt_pos_mad : dt_neg_mad;

          writer.value("DT", dt_mean);
          writer.value("DT_MD", dt_md);
          if (par_work.emit_mad) {
            writer.value("DT_MAD", dt_mad);
          }

          if (par_work.emit_td_summary && ((use_pos_metrics && td_pos) || (!use_pos_metrics && td_neg))) {
            const double td_p2p_mean = use_pos_metrics ? td_pos_p2p_mean : td_neg_p2p_mean;
            const double td_p2p_md = use_pos_metrics ? td_pos_p2p_med : td_neg_p2p_med;
            const double td_p2p_mad = use_pos_metrics ? td_pos_p2p_mad : td_neg_p2p_mad;

            writer.value("TD_P2P", td_p2p_mean);
            writer.value("TD_P2P_MD", td_p2p_md);
            if (par_work.emit_mad) {
              writer.value("TD_P2P_MAD", td_p2p_mad);
            }
          }
        }

        for (int b = 0; b < phase_bins12; ++b) {
          const std::string bin_label =
            std::string("B") + (b < 9 ? "0" : "") + Helper::int2str(b + 1);
          writer.level(bin_label, "BIN");
          writer.value("MEAN", ph12_mean[static_cast<size_t>(b)]);
          writer.unlevel("BIN");
        }

        writer.unlevel("SIG");
      }

      if (par_work.emit_time_domain && par_work.time_window_s > 0.0 && par_work.emit_td_grid) {
        const double half = 0.5 * par_work.time_window_s;
        const int nbins_t = static_cast<int>(std::floor(par_work.time_window_s / par_work.time_bin_s)) + 1;
        if (nbins_t > 0) {
        auto emit_time_grid = [&](const std::vector<const cycle_bound_t*>& tcycles) {
            if (tcycles.empty()) return;
            std::vector<double> sum(static_cast<size_t>(nbins_t), 0.0);
            std::vector<size_t> cnt(static_cast<size_t>(nbins_t), 0);
            std::vector<double> flat;
            flat.reserve(static_cast<size_t>(nbins_t) * tcycles.size());

            for (const auto* cyc : tcycles) {
              uint64_t ta = (par_work.time_lock == "pos") ? cyc->t_pos : cyc->t_neg;
              double t0s = static_cast<double>(ta) * globals::tp_duration - half;
              double t1s = static_cast<double>(ta) * globals::tp_duration + half;
              if (t1s < 0.0) continue;
              uint64_t t_lo = t0s > 0.0 ? static_cast<uint64_t>(llround(t0s / globals::tp_duration)) : 0;
              uint64_t t_hi = static_cast<uint64_t>(llround(t1s / globals::tp_duration));
              size_t i0 = 0;
              size_t i1 = 0;
              if (!window_range_ok(*tp, static_cast<int>(par_work.sr_hz), t_lo, t_hi, &i0, &i1)) {
                continue;
              }
              if (i0 >= sig.size() || i1 >= sig.size() || i0 > i1) continue;
              for (int b = 0; b < nbins_t; ++b) {
                double dt = -half + static_cast<double>(b) * par_work.time_bin_s;
                double tsec = static_cast<double>(ta) * globals::tp_duration + dt;
                if (tsec < 0.0) continue;
                long double tt = tsec / globals::tp_duration;
                double v = interp_value_at_tp_ld(*tp, sig, tt);
                if (!finite(v)) continue;
                sum[static_cast<size_t>(b)] += v;
                cnt[static_cast<size_t>(b)] += 1;
                flat.push_back(v);
              }
            }

            writer.level(sig_lab, "SIG");
            std::mt19937 rng(1337);
            for (int b = 0; b < nbins_t; ++b) {
              double mean = (cnt[static_cast<size_t>(b)] >= static_cast<size_t>(par_work.time_min_n))
                ? sum[static_cast<size_t>(b)] / static_cast<double>(cnt[static_cast<size_t>(b)])
                : std::numeric_limits<double>::quiet_NaN();
              double sec = -half + static_cast<double>(b) * par_work.time_bin_s;
              writer.level(Helper::dbl2str(sec, 3), "SEC");
              writer.value("MEAN", mean);
              if (par_work.do_bootstrap) {
                std::vector<double> vals;
                vals.reserve(tcycles.size());
                for (const auto* cyc : tcycles) {
                  uint64_t ta = (par_work.time_lock == "pos") ? cyc->t_pos : cyc->t_neg;
                  double tsec = static_cast<double>(ta) * globals::tp_duration + sec;
                  if (tsec < 0.0) continue;
                  long double tt = tsec / globals::tp_duration;
                  double v = interp_value_at_tp_ld(*tp, sig, tt);
                  if (finite(v)) vals.push_back(v);
                }

                double se = std::numeric_limits<double>::quiet_NaN();
                double ci_lo = std::numeric_limits<double>::quiet_NaN();
                double ci_hi = std::numeric_limits<double>::quiet_NaN();

                if (vals.size() >= static_cast<size_t>(std::max(2, par_work.time_min_n))) {
                  std::uniform_int_distribution<size_t> dist(0, vals.size() - 1);
                  std::vector<double> boots;
                  boots.reserve(static_cast<size_t>(par_work.bootstrap_n));
                  for (int bi = 0; bi < par_work.bootstrap_n; ++bi) {
                    double sumv = 0.0;
                    for (size_t k = 0; k < vals.size(); ++k) {
                      sumv += vals[dist(rng)];
                    }
                    boots.push_back(sumv / static_cast<double>(vals.size()));
                  }
                  double mu = 0.0;
                  for (double v : boots) mu += v;
                  mu /= static_cast<double>(boots.size());
                  double var = 0.0;
                  for (double v : boots) {
                    double d = v - mu;
                    var += d * d;
                  }
                  var /= static_cast<double>(boots.size());
                  se = std::sqrt(var);

                  const double alpha = 0.5 * (1.0 - par_work.bootstrap_ci);
                  std::sort(boots.begin(), boots.end());
                  size_t lo_idx = static_cast<size_t>(std::floor(alpha * (boots.size() - 1)));
                  size_t hi_idx = static_cast<size_t>(std::floor((1.0 - alpha) * (boots.size() - 1)));
                  ci_lo = boots[lo_idx];
                  ci_hi = boots[hi_idx];
                }

                if (par_work.emit_se) {
                  writer.value("SE", se);
                }
                writer.value("CI_LO", ci_lo);
                writer.value("CI_HI", ci_hi);
              }
              writer.unlevel("SEC");
            }
            if (flat.size() >= 2) {
              std::vector<double> edges;
              edges.reserve(static_cast<size_t>(par_work.amp_bins + 1));
              for (int k = 0; k <= par_work.amp_bins; ++k) {
                double p = static_cast<double>(k) / static_cast<double>(par_work.amp_bins);
                edges.push_back(MiscMath::percentile(flat, p));
              }
              std::vector<std::vector<size_t>> counts(
                static_cast<size_t>(nbins_t),
                std::vector<size_t>(static_cast<size_t>(par_work.amp_bins), 0)
              );
              std::vector<size_t> totals(static_cast<size_t>(nbins_t), 0);

              for (const auto* cyc : tcycles) {
                uint64_t ta = (par_work.time_lock == "pos") ? cyc->t_pos : cyc->t_neg;
                for (int b = 0; b < nbins_t; ++b) {
                  double dt = -half + static_cast<double>(b) * par_work.time_bin_s;
                  double tsec = static_cast<double>(ta) * globals::tp_duration + dt;
                  if (tsec < 0.0) continue;
                  long double tt = tsec / globals::tp_duration;
                  double v = interp_value_at_tp_ld(*tp, sig, tt);
                  if (!finite(v)) continue;
                  int bin = par_work.amp_bins - 1;
                  for (int k = 0; k < par_work.amp_bins; ++k) {
                    double lo = edges[static_cast<size_t>(k)];
                    double hi = edges[static_cast<size_t>(k + 1)];
                    if (k == par_work.amp_bins - 1) {
                      if (v >= lo && v <= hi) { bin = k; break; }
                    } else {
                      if (v >= lo && v < hi) { bin = k; break; }
                    }
                  }
                  counts[static_cast<size_t>(b)][static_cast<size_t>(bin)] += 1;
                  totals[static_cast<size_t>(b)] += 1;
                }
              }
              for (int b = 0; b < nbins_t; ++b) {
                double sec = -half + static_cast<double>(b) * par_work.time_bin_s;
                writer.level(Helper::dbl2str(sec, 3), "SEC");
                for (int k = 0; k < par_work.amp_bins; ++k) {
                  double z = (totals[static_cast<size_t>(b)] >= static_cast<size_t>(par_work.time_min_n))
                    ? static_cast<double>(counts[static_cast<size_t>(b)][static_cast<size_t>(k)]) /
                      static_cast<double>(totals[static_cast<size_t>(b)])
                    : std::numeric_limits<double>::quiet_NaN();
                  writer.level(Helper::int2str(k), "AMP");
                  writer.value("DENS", z);
                  writer.unlevel("AMP");
                }
                writer.unlevel("SEC");
              }
            }
            writer.unlevel("SIG");
          };

          if (has_sig) {
            emit_time_grid(td_cycles);
          }
      }
      }

      if (par_work.emit_ph_amp) {
        std::vector<double> flat;
        flat.reserve(res.bins[c].size() * static_cast<size_t>(res.nbins));
        for (const auto& cyc : res.bins[c]) {
          if (static_cast<int>(cyc.size()) != res.nbins) continue;
          for (double v : cyc) {
            if (finite(v)) flat.push_back(v);
          }
        }
        if (flat.size() >= 2) {
          std::vector<double> edges;
          edges.reserve(static_cast<size_t>(par_work.amp_bins + 1));
          for (int k = 0; k <= par_work.amp_bins; ++k) {
            double p = static_cast<double>(k) / static_cast<double>(par_work.amp_bins);
            edges.push_back(MiscMath::percentile(flat, p));
          }

          const std::string sig_lab =
            (c < sig_labels.size() && !sig_labels[c].empty()) ? sig_labels[c] : ".";

          writer.level(sig_lab, "SIG");
          for (int b = 0; b < res.nbins; ++b) {
            std::vector<size_t> counts(static_cast<size_t>(par_work.amp_bins), 0);
            size_t total = 0;
            for (const auto& cyc : res.bins[c]) {
              if (static_cast<int>(cyc.size()) != res.nbins) continue;
              double v = cyc[static_cast<size_t>(b)];
              if (!finite(v)) continue;
              int bin = par_work.amp_bins - 1;
              for (int k = 0; k < par_work.amp_bins; ++k) {
                double lo = edges[static_cast<size_t>(k)];
                double hi = edges[static_cast<size_t>(k + 1)];
                if (k == par_work.amp_bins - 1) {
                  if (v >= lo && v <= hi) { bin = k; break; }
                } else {
                  if (v >= lo && v < hi) { bin = k; break; }
                }
              }
              counts[static_cast<size_t>(bin)] += 1;
              ++total;
            }
            writer.level(Helper::int2str(b), "PH");
            for (int k = 0; k < par_work.amp_bins; ++k) {
              double z = (total > 0)
                ? static_cast<double>(counts[static_cast<size_t>(k)]) / static_cast<double>(total)
                : std::numeric_limits<double>::quiet_NaN();
              writer.level(Helper::int2str(k), "AMP");
              writer.value("DENS", z);
              writer.unlevel("AMP");
            }
            writer.unlevel("PH");
          }
          writer.unlevel("SIG");
        }
      }
    }
  }

  if (par_work.emit_ph_grid && !mean_bins.empty()) {
    const std::string seg_lab = seg_label.empty() ? "." : seg_label;
    for (size_t c = 0; c < mean_bins.size(); ++c) {
      const std::string sig_lab =
        (c < sig_labels.size() && !sig_labels[c].empty()) ? sig_labels[c] : ".";

      writer.level(sig_lab, "SIG");
      std::mt19937 rng(1337);
      for (int b = 0; b < res.nbins; ++b) {
        writer.level(Helper::int2str(b), "PH");
        writer.value("MEAN", mean_bins[c][static_cast<size_t>(b)]);

        if (par_work.do_bootstrap) {
          std::vector<double> vals;
          vals.reserve(res.bins[c].size());
          for (const auto& cyc : res.bins[c]) {
            if (static_cast<int>(cyc.size()) != res.nbins) continue;
            double v = cyc[static_cast<size_t>(b)];
            if (finite(v)) vals.push_back(v);
          }

          double se = std::numeric_limits<double>::quiet_NaN();
          double ci_lo = std::numeric_limits<double>::quiet_NaN();
          double ci_hi = std::numeric_limits<double>::quiet_NaN();

          if (vals.size() >= 2) {
            std::uniform_int_distribution<size_t> dist(0, vals.size() - 1);
            std::vector<double> boots;
            boots.reserve(static_cast<size_t>(par_work.bootstrap_n));
            for (int bi = 0; bi < par_work.bootstrap_n; ++bi) {
              double sum = 0.0;
              for (size_t k = 0; k < vals.size(); ++k) {
                sum += vals[dist(rng)];
              }
              boots.push_back(sum / static_cast<double>(vals.size()));
            }
            double mu = 0.0;
            for (double v : boots) mu += v;
            mu /= static_cast<double>(boots.size());
            double var = 0.0;
            for (double v : boots) {
              double d = v - mu;
              var += d * d;
            }
            var /= static_cast<double>(boots.size());
            se = std::sqrt(var);

            const double alpha = 0.5 * (1.0 - par_work.bootstrap_ci);
            std::sort(boots.begin(), boots.end());
            size_t lo_idx = static_cast<size_t>(std::floor(alpha * (boots.size() - 1)));
            size_t hi_idx = static_cast<size_t>(std::floor((1.0 - alpha) * (boots.size() - 1)));
            ci_lo = boots[lo_idx];
            ci_hi = boots[hi_idx];
          }

          if (par_work.emit_se) {
            writer.value("SE", se);
          }
          writer.value("CI_LO", ci_lo);
          writer.value("CI_HI", ci_hi);
        }

        writer.unlevel("PH");
      }
      writer.unlevel("SIG");
      
    }
  }

  // all done for this seed
  writer.unlevel("SEED");
  
  res.cycles = std::move(cycles);
  return res;
}
