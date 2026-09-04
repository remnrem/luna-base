//    --------------------------------------------------------------------
//
//    This file is part of Luna.
//
//    --------------------------------------------------------------------

#include "stats/dpp-vector.h"
#include "stats/dpp.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "miscmath/miscmath.h"
#include "param.h"
#include "stats/dpp-fit.h"
#include "stats/dpp-io.h"
#include "stats/dpp-spec.h"
#include "timeline/timeline.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <map>
#include <set>

extern logger_t logger;

namespace {

std::vector<std::string> feature_set(const param_t &param) {
  std::vector<std::string> f;
  if (param.has("vector-features"))
    f = param.strvector("vector-features");
  else
    f = {"EMBEDDING", "HYPNOGRAM", "EMBEDDING_GEOMETRY",
         "EMBEDDING_DYNAMICS"};

  for (int i = 0; i < (int)f.size(); i++)
    for (int j = 0; j < (int)f[i].size(); j++)
      f[i][j] = (char)std::toupper((unsigned char)f[i][j]);
  return f;
}

bool has_feature(const std::vector<std::string> &f, const std::string &x) {
  return std::find(f.begin(), f.end(), x) != f.end();
}

std::string time_mode(const param_t &param) {
  std::string x =
      param.has("vector-time") ? param.value("vector-time") : "RELATIVE";
  for (char &c : x)
    c = (char)std::toupper((unsigned char)c);
  if (x != "RELATIVE" && x != "EDF" && x != "ONSET")
    Helper::halt("vector-time= must be RELATIVE, ONSET, or EDF");
  return x;
}

void add_label(std::vector<std::string> *l, const std::string &s) {
  l->push_back(s);
}

double safe_cosine(const std::vector<double> &a, const std::vector<double> &b) {
  double ab = 0, aa = 0, bb = 0;
  const int n = std::min(a.size(), b.size());
  for (int i = 0; i < n; i++)
    if (std::isfinite(a[i]) && std::isfinite(b[i])) {
      ab += a[i] * b[i];
      aa += a[i] * a[i];
      bb += b[i] * b[i];
    }
  if (aa <= 0 || bb <= 0)
    return std::numeric_limits<double>::quiet_NaN();
  return ab / std::sqrt(aa * bb);
}

double safe_distance(const std::vector<double> &a,
                     const std::vector<double> &b) {
  double s = 0;
  int n = 0;
  const int m = std::min(a.size(), b.size());
  for (int i = 0; i < m; i++)
    if (std::isfinite(a[i]) && std::isfinite(b[i])) {
      const double d = a[i] - b[i];
      s += d * d;
      ++n;
    }
  return n > 0 ? std::sqrt(s) : std::numeric_limits<double>::quiet_NaN();
}

double norm(const std::vector<double> &x) {
  double s = 0;
  int n = 0;
  for (int i = 0; i < (int)x.size(); i++)
    if (std::isfinite(x[i])) {
      s += x[i] * x[i];
      ++n;
    }
  return n > 0 ? std::sqrt(s) : std::numeric_limits<double>::quiet_NaN();
}

// Find the current epoch using the standard epoch geometry.  This keeps
// vector DPP independent of the vector sampling rate, including one
// sample per long EDF record.
int epoch_at(edf_t &edf, const uint64_t tp) {
  if (!edf.timeline.epoched())
    return -1;
  const uint64_t elen = edf.timeline.epoch_len_tp_uint64_t();
  const uint64_t einc = edf.timeline.epoch_increment_tp();
  const int ne = edf.timeline.num_total_epochs();
  if (elen == 0 || einc == 0 || ne == 0)
    return -1;
  return MiscMath::position2leftepoch(tp, elen, einc, ne);
}

} // namespace

bool dpp_vector::enabled(const param_t &param) {
  return !(param.has("classic") && param.yesno("classic")) &&
         !param.has("dpp-classic-capture");
}

dpp_vector::layout_t dpp_vector::layout(const int embedding_dim,
                                        const param_t &param) {
  const std::vector<std::string> f = feature_set(param);
  layout_t l;
  l.embedding_dim = embedding_dim;
  l.raw = has_feature(f, "EMBEDDING");
  l.context = has_feature(f, "HYPNOGRAM");
  l.geom = has_feature(f, "EMBEDDING_GEOMETRY");
  l.dyn = has_feature(f, "EMBEDDING_DYNAMICS");
  int o = 0;
  l.raw_offset = l.raw ? o : -1;
  if (l.raw)
    o += embedding_dim;
  l.context_offset = l.context ? o : -1;
  if (l.context)
    o += 10;
  l.geom_offset = l.geom ? o : -1;
  if (l.geom)
    o += 5;
  l.dyn_offset = l.dyn ? o : -1;
  if (l.dyn)
    o += 4;
  return l;
}

std::vector<std::string> dpp_vector::labels(const int n_channels,
                                            const param_t &param) {
  const std::vector<std::string> f = feature_set(param);
  const std::string tm = time_mode(param);
  std::vector<std::string> l;

  if (has_feature(f, "EMBEDDING"))
    for (int i = 0; i < n_channels; i++)
      add_label(&l, "VEC.EMBEDDING.V" + Helper::int2str(i + 1));

  if (has_feature(f, "HYPNOGRAM")) {
    add_label(&l, tm == "EDF"     ? "VEC.HYPNOGRAM.EDF_TIME_FRAC"
                  : tm == "ONSET" ? "VEC.HYPNOGRAM.ONSET_TIME_H"
                                  : "VEC.HYPNOGRAM.RETAINED_FRAC");
    add_label(&l, "VEC.HYPNOGRAM.STAGE_W");
    add_label(&l, "VEC.HYPNOGRAM.STAGE_N1");
    add_label(&l, "VEC.HYPNOGRAM.STAGE_N2");
    add_label(&l, "VEC.HYPNOGRAM.STAGE_N3");
    add_label(&l, "VEC.HYPNOGRAM.STAGE_R");
    add_label(&l, "VEC.HYPNOGRAM.STAGE_UNKNOWN");
    add_label(&l, "VEC.HYPNOGRAM.CYCLE");
    add_label(&l, "VEC.HYPNOGRAM.CYCLE_PHASE");
    add_label(&l, "VEC.HYPNOGRAM.VALID_DIM_FRAC");
  }

  if (has_feature(f, "EMBEDDING_GEOMETRY")) {
    add_label(&l, "VEC.EMBEDDING_GEOMETRY.NORM");
    add_label(&l, "VEC.EMBEDDING_GEOMETRY.BASELINE_DIST");
    add_label(&l, "VEC.EMBEDDING_GEOMETRY.BASELINE_COS");
    add_label(&l, "VEC.EMBEDDING_GEOMETRY.PREV_DIST");
    add_label(&l, "VEC.EMBEDDING_GEOMETRY.PREV_COS");
  }

  if (has_feature(f, "EMBEDDING_DYNAMICS")) {
    add_label(&l, "VEC.EMBEDDING_DYNAMICS.VELOCITY");
    add_label(&l, "VEC.EMBEDDING_DYNAMICS.ACCELERATION");
    add_label(&l, "VEC.EMBEDDING_DYNAMICS.NORM_DELTA");
    add_label(&l, "VEC.EMBEDDING_DYNAMICS.TURN_ANGLE");
  }

  if (has_feature(f, "DSP")) {
    dpp_specs_t specs;
    if (param.has("spec"))
      specs.read(param.value("spec"));
    else {
      const std::vector<std::string> channels = param.strvector("dsp-sig");
      if (channels.empty())
        Helper::halt("vector-features=DSP requires spec= or dsp-sig=");
      specs.init_default(channels);
      specs.apply_inline_overrides(
          param.has("windows") ? param.value("windows") : "",
          param.has("filters") ? param.value("filters") : "",
          param.has("features") ? param.value("features") : "",
          param.has("prefilter") ? param.value("prefilter") : "");
    }
    int n = 0;
    for (const dpp_spec_t &s : specs.specs)
      n += s.cols();
    for (int j = 0; j < n; j++)
      add_label(&l, "DPP.DSP.F" + Helper::int2str(j + 1));
  }

  return l;
}

bool dpp_vector::run(edf_t &edf, param_t &param) {
  if (!enabled(param))
    return false;

  signal_list_t signals = edf.header.signal_list(param.requires("sig"), true);
  if (signals.size() == 0) {
    logger << "  *** no signals selected for DPP vector mode\n";
    return true;
  }

  const int nc = signals.size();
  std::vector<double> fs(nc);
  std::vector<std::vector<double>> x(nc);
  std::vector<std::vector<uint64_t>> tp(nc);

  if (!param.has("hypno-context") || param.yesno("hypno-context")) {
    edf.timeline.ensure_epoched();

    const bool has_hypno = edf.timeline.epoch_annotation("W") ||
                           edf.timeline.epoch_annotation("N1") ||
                           edf.timeline.epoch_annotation("N2") ||
                           edf.timeline.epoch_annotation("N3") ||
                           edf.timeline.epoch_annotation("R");
    if (!has_hypno)
      Helper::halt("DPP vector mode requires HYPNO output when "
                   "hypno-context=T");
  }

  double common_fs = 0;
  for (int c = 0; c < nc; c++) {
    fs[c] = edf.header.sampling_freq(signals(c));
    if (!std::isfinite(fs[c]) || fs[c] <= 0)
      Helper::halt("DPP vector mode requires positive sampling rates");
    if (c == 0)
      common_fs = fs[c];
    else if (std::fabs(fs[c] - common_fs) > 1e-9 * std::max(1.0, common_fs))
      Helper::halt("DPP vector mode requires all selected channels to have the "
                   "same sampling rate");

    // Pull the complete continuous trace and apply any epoch MASK below at
    // the vector-row level; this avoids requiring EDF+D/RESTRUCTURE.
    slice_t w(edf, signals(c), edf.timeline.wholetrace(true));
    x[c] = *w.pdata();
    tp[c] = *w.ptimepoints();
    if (x[c].size() != tp[c].size())
      Helper::halt("DPP vector mode found invalid signal time points");
  }

  for (int c = 1; c < nc; c++)
    if (tp[c] != tp[0])
      Helper::halt("DPP vector mode requires time-aligned selected channels");

  // Index in the unmasked vector sequence.  It lets derived features detect
  // a masked gap after the retained rows have been compacted.
  std::vector<int> source_index(tp[0].size());
  for (int r = 0; r < (int)source_index.size(); r++)
    source_index[r] = r;

  // Vector observations are low-rate samples, so a normal Luna epoch MASK
  // cannot be delegated to slice_t/wholetrace(): those operate on retained
  // EDF records and otherwise require RE (RESTRUCTURE). Honor the current
  // epoch mask explicitly here, before deriving context, geometry, or
  // dynamics. This keeps the source EDF continuous and also prevents
  // prev/next-derived features from crossing a masked region.
  if (edf.timeline.is_epoch_mask_set()) {
    std::vector<int> keep;
    keep.reserve(tp[0].size());
    for (int r = 0; r < (int)tp[0].size(); r++) {
      const int epoch = epoch_at(edf, tp[0][r]);
      if (epoch >= 0 && !edf.timeline.masked_epoch(epoch))
        keep.push_back(r);
    }

    const int n_before = tp[0].size();
    for (int c = 0; c < nc; c++) {
      std::vector<double> xx;
      std::vector<uint64_t> tt;
      xx.reserve(keep.size());
      tt.reserve(keep.size());
      for (int i : keep) {
        xx.push_back(x[c][i]);
        tt.push_back(tp[c][i]);
      }
      x[c].swap(xx);
      tp[c].swap(tt);
    }

    std::vector<int> kept_source_index;
    kept_source_index.reserve(keep.size());
    for (int i : keep)
      kept_source_index.push_back(source_index[i]);
    source_index.swap(kept_source_index);

    logger << "  DPP vector mode: epoch MASK retained " << keep.size() << " of "
           << n_before << " observations\n";
    if (keep.empty()) {
      Helper::halt(
          "DPP vector mode: no observations remain after the epoch MASK for " +
          edf.id);
    }
  }

  const int nr = tp[0].size();
  if (nr == 0) {
    logger << "  *** no vector observations for DPP\n";
    return true;
  }

  // Cache cycle extents in epoch coordinates so a low-rate vector sample
  // (including one sample per long EDF record) can still receive a stable
  // within-cycle phase value without assuming that vector and hypnogram
  // sampling rates match.
  std::vector<int> epoch_cycle;
  std::map<int, std::pair<int, int>> cycle_extent;
  if ((!param.has("hypno-context") || param.yesno("hypno-context")) &&
      edf.timeline.epoched()) {
    const int ne = edf.timeline.num_total_epochs();
    epoch_cycle.assign(ne, 0);
    for (int e = 0; e < ne; e++)
      for (int k = 1; k <= 8; k++)
        if (edf.timeline.epoch_annotation("_NREMC_" + Helper::int2str(k), e)) {
          epoch_cycle[e] = k;
          break;
        }
    for (int e = 0; e < ne; e++)
      if (epoch_cycle[e]) {
        if (cycle_extent.find(epoch_cycle[e]) == cycle_extent.end())
          cycle_extent[epoch_cycle[e]] = std::make_pair(e, e);
        else {
          cycle_extent[epoch_cycle[e]].first =
              std::min(cycle_extent[epoch_cycle[e]].first, e);
          cycle_extent[epoch_cycle[e]].second =
              std::max(cycle_extent[epoch_cycle[e]].second, e);
        }
      }
  }

  std::vector<std::string> flabels = labels(nc, param);
  const std::vector<std::string> fset = feature_set(param);
  const std::string tm = time_mode(param);
  const double time_origin =
      tm == "ONSET" && !tp[0].empty() ? tp[0][0] / (double)globals::tp_1sec : 0;
  const bool raw = has_feature(fset, "EMBEDDING");
  const bool context = has_feature(fset, "HYPNOGRAM");
  const bool geom = has_feature(fset, "EMBEDDING_GEOMETRY");
  const bool dyn = has_feature(fset, "EMBEDDING_DYNAMICS");
  const bool dsp = has_feature(fset, "DSP");

  std::vector<double> baseline(nc, 0);
  std::vector<int> baseline_n(nc, 0);
  for (int r = 0; r < nr; r++)
    for (int c = 0; c < nc; c++)
      if (std::isfinite(x[c][r])) {
        baseline[c] += x[c][r];
        ++baseline_n[c];
      }
  for (int c = 0; c < nc; c++)
    if (baseline_n[c])
      baseline[c] /= baseline_n[c];

  dpp_matrix_t mat;
  mat.id = edf.id;
  mat.time_sec.resize(nr);
  mat.X.resize(nr);

  for (int r = 0; r < nr; r++) {
    const double t = tp[0][r] / (double)globals::tp_1sec;
    mat.time_sec[r] = t;
    std::vector<double> cur(nc), prev(nc), next(nc);
    for (int c = 0; c < nc; c++) {
      cur[c] = x[c][r];
      const bool prev_row = r > 0 && source_index[r] == source_index[r - 1] + 1;
      const bool next_row =
          r + 1 < nr && source_index[r + 1] == source_index[r] + 1;
      prev[c] =
          prev_row ? x[c][r - 1] : std::numeric_limits<double>::quiet_NaN();
      next[c] =
          next_row ? x[c][r + 1] : std::numeric_limits<double>::quiet_NaN();
    }

    std::vector<double> row;
    if (raw)
      row.insert(row.end(), cur.begin(), cur.end());

    int epoch = epoch_at(edf, tp[0][r]);
    std::string stage = "UNKNOWN";
    int cycle = 0;
    if (epoch >= 0) {
      if (edf.timeline.epoch_annotation("W", epoch))
        stage = "W";
      else if (edf.timeline.epoch_annotation("N1", epoch))
        stage = "N1";
      else if (edf.timeline.epoch_annotation("N2", epoch))
        stage = "N2";
      else if (edf.timeline.epoch_annotation("N3", epoch))
        stage = "N3";
      else if (edf.timeline.epoch_annotation("R", epoch))
        stage = "R";
      for (int k = 1; k <= 8; k++)
        if (edf.timeline.epoch_annotation("_NREMC_" + Helper::int2str(k),
                                          epoch)) {
          cycle = k;
          break;
        }
    }

    double cycle_phase = std::numeric_limits<double>::quiet_NaN();
    if (cycle > 0 && epoch >= 0 && epoch < (int)epoch_cycle.size()) {
      const std::pair<int, int> ex = cycle_extent[cycle];
      cycle_phase = ex.second > ex.first ? (double)(epoch - ex.first) /
                                               (double)(ex.second - ex.first)
                                         : 0;
      if (cycle_phase < 0)
        cycle_phase = 0;
      if (cycle_phase > 1)
        cycle_phase = 1;
    }

    if (context) {
      double frac = nr > 1 ? (double)r / (double)(nr - 1) : 0;
      if (tm == "EDF") {
        const double edf_sec = edf.header.nr * edf.header.record_duration;
        frac = edf_sec > 0 ? t / edf_sec : 0;
        if (frac < 0)
          frac = 0;
        if (frac > 1)
          frac = 1;
      } else if (tm == "ONSET")
        frac = (t - time_origin) / 3600.0;
      row.push_back(frac);
      row.push_back(stage == "W");
      row.push_back(stage == "N1");
      row.push_back(stage == "N2");
      row.push_back(stage == "N3");
      row.push_back(stage == "R");
      row.push_back(stage == "UNKNOWN");
      row.push_back(cycle);
      row.push_back(cycle_phase);
      int nv = 0;
      for (int c = 0; c < nc; c++)
        if (std::isfinite(cur[c]))
          ++nv;
      row.push_back(nc > 0 ? (double)nv / nc : 0);
    }

    const double cur_norm = norm(cur);
    const double base_dist = safe_distance(cur, baseline);
    const double base_cos = safe_cosine(cur, baseline);
    const double prev_dist = safe_distance(cur, prev);
    const double prev_cos = safe_cosine(cur, prev);
    if (geom) {
      row.push_back(cur_norm);
      row.push_back(base_dist);
      row.push_back(base_cos);
      row.push_back(prev_dist);
      row.push_back(prev_cos);
    }

    const bool prev_row = r > 0 && source_index[r] == source_index[r - 1] + 1;
    const bool next_row =
        r + 1 < nr && source_index[r + 1] == source_index[r] + 1;
    const double dt =
        prev_row ? (tp[0][r] - tp[0][r - 1]) / (double)globals::tp_1sec : 0;
    const double dt1 =
        next_row ? (tp[0][r + 1] - tp[0][r]) / (double)globals::tp_1sec : dt;
    const double vel =
        dt > 0 ? prev_dist / dt : std::numeric_limits<double>::quiet_NaN();
    const double next_dist = safe_distance(next, cur);
    const double vel_next =
        dt1 > 0 ? next_dist / dt1 : std::numeric_limits<double>::quiet_NaN();
    double turn = std::numeric_limits<double>::quiet_NaN();
    if (prev_row && next_row) {
      std::vector<double> v1(nc), v2(nc);
      for (int c = 0; c < nc; c++) {
        v1[c] = cur[c] - prev[c];
        v2[c] = next[c] - cur[c];
      }
      const double co = safe_cosine(v1, v2);
      if (std::isfinite(co))
        turn = std::acos(std::max(-1.0, std::min(1.0, co)));
    }
    if (dyn) {
      row.push_back(vel);
      row.push_back(std::isfinite(vel) && std::isfinite(vel_next)
                        ? std::fabs(vel_next - vel) / std::max(dt1, 1e-9)
                        : std::numeric_limits<double>::quiet_NaN());
      row.push_back(std::isfinite(cur_norm) && std::isfinite(norm(prev))
                        ? cur_norm - norm(prev)
                        : std::numeric_limits<double>::quiet_NaN());
      row.push_back(turn);
    }

    mat.X[r] = row;
  }

  if (dsp) {
    dpp_matrix_t classic;
    if (!dpp_classic::extract(edf, param, mat.time_sec, &classic))
      Helper::halt("DPP vector mode: failed to extract DSP features");
    if (classic.time_sec.size() != mat.time_sec.size())
      Helper::halt("DPP vector mode DSP/vector row-count mismatch; use the "
                   "same step= as the vector sampling interval");
    for (int r = 0; r < nr; r++) {
      if (std::fabs(classic.time_sec[r] - mat.time_sec[r]) > 1e-6)
        Helper::halt("DPP vector mode DSP/vector time-grid mismatch at row " +
                     Helper::int2str(r) + "; use the same step= as the vector "
                     "sampling interval");
      mat.X[r].insert(mat.X[r].end(), classic.X[r].begin(), classic.X[r].end());
    }
    logger << "  DPP vector mode: appended "
           << (classic.X.empty() ? 0 : classic.X[0].size())
           << " DSP feature columns on the aligned vector grid\n";
  }

  if ((int)mat.X[0].size() != (int)flabels.size())
    Helper::halt("internal DPP vector feature-label mismatch");

  if (param.has("data"))
    dpp_io::save(param.value("data"), mat, (int)flabels.size(), false);

  logger << "  DPP vector mode: " << nc << " channels, " << common_fs << " Hz, "
         << nr << " observations, " << flabels.size() << " features\n";

  if (param.has("model")) {
    dpp_specs_t dummy_specs;
    dpp_fit::apply(edf, param, dummy_specs, mat, NULL);
  }
  return true;
}
