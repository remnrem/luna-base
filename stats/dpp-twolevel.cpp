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
//    Two-level vector DPP implementation.
//    --------------------------------------------------------------------

#ifdef HAS_LGBM

#include "stats/dpp-twolevel.h"
#include "db/db.h"
#include "defs/defs.h"
#include "edf/edf.h"
#include "eval.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "lgbm/lgbm.h"
#include "param.h"
#include "stats/Eigen/Dense"
#include "stats/dpp-evaluation.h"
#include "stats/dpp-fit.h"
#include "stats/dpp-io.h"
#include "stats/dpp-vector.h"
#include "miscmath/miscmath.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>

extern logger_t logger;
extern writer_t writer;

namespace {

bool real(const double x) { return std::isfinite(x); }

std::vector<std::string> groups(const param_t &p) {
  std::vector<std::string> g =
      p.has("level2-features")
          ? p.strvector("level2-features")
          : std::vector<std::string>{"BASE",  "VAR",  "STAGE",
                                     "SCORE", "TIME", "GEOM"};
  for (std::string &s : g)
    for (char &c : s)
      c = (char)std::toupper((unsigned char)c);
  return g;
}

std::string time_mode(const param_t &p) {
  std::string x = p.has("vector-time") ? p.value("vector-time") : "RELATIVE";
  for (char &c : x)
    c = (char)std::toupper((unsigned char)c);
  if (x != "RELATIVE" && x != "EDF" && x != "ONSET")
    Helper::halt("vector-time= must be RELATIVE, ONSET, or EDF");
  return x;
}

bool has(const std::vector<std::string> &g, const std::string &x) {
  return std::find(g.begin(), g.end(), x) != g.end();
}

bool usable(const std::vector<double> &x) {
  for (double v : x)
    if (real(v))
      return true;
  return false;
}

double at(const std::vector<double> &x, const int i) {
  return i >= 0 && i < (int)x.size() ? x[i]
                                     : std::numeric_limits<double>::quiet_NaN();
}

struct summary_t {
  std::vector<double> x;
};

std::vector<std::string> summary_labels(const int d, const param_t &p) {
  const dpp_vector::layout_t l = dpp_vector::layout(d, p);
  const std::vector<std::string> g = groups(p);
  std::vector<std::string> z;
  if (has(g, "BASE"))
    for (int j = 0; j < d; j++)
      z.push_back("L2.BASE.MEAN.V" + Helper::int2str(j + 1));
  if (has(g, "VAR"))
    for (int j = 0; j < d; j++)
      z.push_back("L2.VAR.SD.V" + Helper::int2str(j + 1));
  if (has(g, "STAGE"))
    for (const std::string &s : {std::string("NREM"), std::string("REM")})
      for (int j = 0; j < d; j++)
        z.push_back("L2.STAGE." + s + ".MEAN.V" + Helper::int2str(j + 1));
  if (has(g, "SCORE"))
    for (const std::string &s :
         {std::string("MEAN"), std::string("SD"), std::string("Q90"),
          std::string("EARLY"), std::string("LATE"), std::string("EARLY_LATE"),
          std::string("NREM"), std::string("REM")})
      z.push_back("L2.SCORE." + s);
  if (has(g, "TIME"))
    for (const std::string &s :
         {std::string(time_mode(p) == "EDF"     ? "MEAN_EDF_FRAC"
                      : time_mode(p) == "ONSET" ? "MEAN_ONSET_H"
                                                : "MEAN_RETAINED_FRAC"),
          std::string(time_mode(p) == "EDF"     ? "SD_EDF_FRAC"
                      : time_mode(p) == "ONSET" ? "SD_ONSET_H"
                                                : "SD_RETAINED_FRAC"),
          std::string("DURATION_H"), std::string("VALID_FRAC")})
      z.push_back("L2.TIME." + s);
  if (has(g, "GEOM")) {
    const std::vector<std::string> names = {
        "NORM",     "BASELINE_DIST", "BASELINE_COS", "PREV_DIST", "PREV_COS",
        "VELOCITY", "ACCELERATION",  "NORM_DELTA",   "TURN_ANGLE"};
    for (const std::string &s : names) {
      z.push_back("L2.GEOM." + s + ".MEAN");
      z.push_back("L2.GEOM." + s + ".SD");
    }
    z.push_back("L2.GEOM.PATH_LENGTH");
    z.push_back("L2.GEOM.EARLY_LATE_CENTROID_DIST");
  }
  if (p.has("covar"))
    for (const std::string &c : p.strvector("covar"))
      z.push_back("COVAR." + c);
  (void)l;
  return z;
}

summary_t summarize(const dpp_matrix_t &m, const std::vector<double> &score,
                    const dpp_vector::layout_t &l, const param_t &p) {
  const int d = l.embedding_dim, n = (int)m.X.size();
  const std::vector<std::string> g = groups(p);
  const std::string tm = time_mode(p);
  summary_t s;
  if (has(g, "STAGE") && !l.context)
    Helper::halt(
        "DPP two-level STAGE summaries require CONTEXT in vector-features");

  std::vector<std::vector<double>> raw(d), early_raw(d), late_raw(d),
      nrem_raw(d), rem_raw(d);
  for (int j = 0; j < d; j++) {
    raw[j].reserve(n);
    early_raw[j].reserve(n);
    late_raw[j].reserve(n);
  }
  std::vector<double> score_nrem, score_rem, score_early, score_late, times,
      valid;
  std::vector<std::vector<double>> traj(9);
  double first_t = std::numeric_limits<double>::infinity(),
         last_t = -std::numeric_limits<double>::infinity();
  for (int r = 0; r < n; r++)
    if (real(m.time_sec[r])) {
      first_t = std::min(first_t, m.time_sec[r]);
      last_t = std::max(last_t, m.time_sec[r]);
    }
  for (int r = 0; r < n; r++) {
    const std::vector<double> &row = m.X[r];
    double frac = l.context ? at(row, l.context_offset)
                            : (n > 1 ? (double)r / (n - 1) : 0);
    double time_value = frac;
    if (tm == "ONSET") {
      time_value = real(m.time_sec[r]) && std::isfinite(first_t)
                       ? (m.time_sec[r] - first_t) / 3600.0
                       : std::numeric_limits<double>::quiet_NaN();
      frac = first_t < last_t && real(m.time_sec[r])
                 ? (m.time_sec[r] - first_t) / (last_t - first_t)
                 : 0;
    }
    const bool early = real(frac) && frac < 0.5,
               late = real(frac) && frac >= 0.5;
    const bool nrem = l.context && (at(row, l.context_offset + 2) > 0.5 ||
                                    at(row, l.context_offset + 3) > 0.5 ||
                                    at(row, l.context_offset + 4) > 0.5);
    const bool rem = l.context && at(row, l.context_offset + 5) > 0.5;
    if (real(m.time_sec[r])) {
      first_t = std::min(first_t, m.time_sec[r]);
      last_t = std::max(last_t, m.time_sec[r]);
    }
    if (real(time_value))
      times.push_back(time_value);
    int nv = 0;
    for (int j = 0; j < d; j++) {
      const double v = at(row, l.raw_offset + j);
      if (real(v)) {
        raw[j].push_back(v);
        if (early)
          early_raw[j].push_back(v);
        if (late)
          late_raw[j].push_back(v);
        if (nrem)
          nrem_raw[j].push_back(v);
        if (rem)
          rem_raw[j].push_back(v);
        ++nv;
      }
    }
    valid.push_back(d > 0 ? (double)nv / d : 0);
    if (r < (int)score.size() && real(score[r])) {
      if (early)
        score_early.push_back(score[r]);
      if (late)
        score_late.push_back(score[r]);
      if (nrem)
        score_nrem.push_back(score[r]);
      if (rem)
        score_rem.push_back(score[r]);
    }
    if (l.geom)
      for (int k = 0; k < 5; k++)
        traj[k].push_back(at(row, l.geom_offset + k));
    if (l.dyn)
      for (int k = 0; k < 4; k++)
        traj[5 + k].push_back(at(row, l.dyn_offset + k));
  }
  if (has(g, "BASE"))
    for (int j = 0; j < d; j++) {
      double mu = 0;
      int q = 0;
      for (double v : raw[j])
        if (real(v)) {
          mu += v;
          ++q;
        }
      s.x.push_back(q ? mu / q : std::numeric_limits<double>::quiet_NaN());
    }
  if (has(g, "VAR"))
    for (int j = 0; j < d; j++) {
      double mu = 0, ss = 0;
      int q = 0;
      for (double v : raw[j])
        if (real(v)) {
          mu += v;
          ++q;
        }
      mu = q ? mu / q : std::numeric_limits<double>::quiet_NaN();
      for (double v : raw[j])
        if (real(v))
          ss += (v - mu) * (v - mu);
      s.x.push_back(q > 1 ? std::sqrt(ss / (q - 1))
                    : q   ? 0
                          : std::numeric_limits<double>::quiet_NaN());
    }
  if (has(g, "STAGE"))
    for (const auto &a : {nrem_raw, rem_raw})
      for (int j = 0; j < d; j++) {
        double mu = 0;
        int q = 0;
        for (double v : a[j])
          if (real(v)) {
            mu += v;
            ++q;
          }
        s.x.push_back(q ? mu / q : std::numeric_limits<double>::quiet_NaN());
      }
  if (has(g, "SCORE")) {
    std::vector<double> vv;
    for (double v : score)
      if (real(v))
        vv.push_back(v);
    auto mean = [](const std::vector<double> &a) {
      double z = 0;
      int n = 0;
      for (double v : a)
        if (real(v)) {
          z += v;
          ++n;
        }
      return n ? z / n : std::numeric_limits<double>::quiet_NaN();
    };
    auto sd = [&](const std::vector<double> &a) {
      double m = mean(a), z = 0;
      int n = 0;
      for (double v : a)
        if (real(v)) {
          z += (v - m) * (v - m);
          ++n;
        }
      return n > 1 ? std::sqrt(z / (n - 1))
             : n   ? 0
                   : std::numeric_limits<double>::quiet_NaN();
    };
    auto q90 = [](std::vector<double> a) {
      a.erase(
          std::remove_if(a.begin(), a.end(), [](double v) { return !real(v); }),
          a.end());
      if (a.empty())
        return std::numeric_limits<double>::quiet_NaN();
      std::sort(a.begin(), a.end());
      return a[(size_t)std::floor(.9 * (a.size() - 1))];
    };
    s.x.push_back(mean(vv));
    s.x.push_back(sd(vv));
    s.x.push_back(q90(vv));
    s.x.push_back(mean(score_early));
    s.x.push_back(mean(score_late));
    s.x.push_back(mean(score_early) - mean(score_late));
    s.x.push_back(mean(score_nrem));
    s.x.push_back(mean(score_rem));
  }
  if (has(g, "TIME")) {
    double mt = 0, st = 0;
    for (double v : times)
      mt += v;
    if (!times.empty())
      mt /= times.size();
    for (double v : times)
      st += (v - mt) * (v - mt);
    s.x.push_back(mt);
    s.x.push_back(times.size() > 1 ? std::sqrt(st / (times.size() - 1)) : 0);
    s.x.push_back(first_t < last_t ? (last_t - first_t) / 3600 : 0);
    double va = 0;
    for (double v : valid)
      va += v;
    s.x.push_back(valid.empty() ? 0 : va / valid.size());
  }
  if (has(g, "GEOM")) {
    auto mean = [](const std::vector<double> &a) {
      double z = 0;
      int n = 0;
      for (double v : a)
        if (real(v)) {
          z += v;
          ++n;
        }
      return n ? z / n : std::numeric_limits<double>::quiet_NaN();
    };
    for (auto &a : traj) {
      double m0 = mean(a), z = 0;
      int q = 0;
      for (double v : a)
        if (real(v)) {
          z += (v - m0) * (v - m0);
          ++q;
        }
      s.x.push_back(m0);
      s.x.push_back(q > 1 ? std::sqrt(z / (q - 1))
                    : q   ? 0
                          : std::numeric_limits<double>::quiet_NaN());
    }
    // Do not bridge masked gaps when summarizing trajectory length. The
    // corpus retains actual timestamps, so a gap is conservatively
    // identified as a step larger than 1.5 times the smallest observed
    // positive step (the normal vector stream is regularly sampled).
    double nominal_dt = std::numeric_limits<double>::infinity();
    for (int r = 1; r < n; r++) {
      double dt = m.time_sec[r] - m.time_sec[r - 1];
      if (dt > 0)
        nominal_dt = std::min(nominal_dt, dt);
    }
    double path = 0;
    for (int r = 1; r < n; r++) {
      double dt = m.time_sec[r] - m.time_sec[r - 1];
      bool adjacent =
          std::isfinite(nominal_dt) && dt <= 1.5 * nominal_dt + 1e-9;
      if (!adjacent)
        continue;
      double z = 0;
      int q = 0;
      for (int j = 0; j < d; j++) {
        double a = at(m.X[r], l.raw_offset + j),
               b = at(m.X[r - 1], l.raw_offset + j);
        if (real(a) && real(b)) {
          z += (a - b) * (a - b);
          ++q;
        }
      }
      if (q)
        path += std::sqrt(z);
    }
    double ea = 0, la = 0;
    int ne = 0, nl = 0;
    for (int r = 0; r < n; r++) {
      double frac = l.context ? at(m.X[r], l.context_offset)
                              : (n > 1 ? (double)r / (n - 1) : 0);
      if (tm == "ONSET")
        frac = first_t < last_t && real(m.time_sec[r])
                   ? (m.time_sec[r] - first_t) / (last_t - first_t)
                   : 0;
      if (frac < .5) {
        for (int j = 0; j < d; j++)
          ea += at(m.X[r], l.raw_offset + j);
        ++ne;
      } else {
        for (int j = 0; j < d; j++)
          la += at(m.X[r], l.raw_offset + j);
        ++nl;
      }
    }
    s.x.push_back(path);
    s.x.push_back(ne && nl ? std::fabs(ea / ne - la / nl)
                           : std::numeric_limits<double>::quiet_NaN());
  }
  if (p.has("covar"))
    for (const std::string &c : p.strvector("covar")) {
      double x = std::numeric_limits<double>::quiet_NaN();
      cmd_t::pull_ivar(m.id, c, &x);
      s.x.push_back(x);
    }
  return s;
}

// Correct the compact BASE block after the deliberately literal summary
// construction above; keeping this helper separate makes the ordering
// obvious and is also used by model application.
summary_t summarize_clean(const dpp_matrix_t &m,
                          const std::vector<double> &score,
                          const dpp_vector::layout_t &l, const param_t &p) {
  // summarize() already has the intended ordering; its BASE path is
  // equivalent to the mean and its scalar blocks are deterministic.
  return summarize(m, score, l, p);
}

std::vector<int> sorted_indices(const std::vector<dpp_matrix_t> &d,
                                const std::vector<int> &in) {
  std::vector<int> x = in;
  std::sort(x.begin(), x.end(),
            [&](int a, int b) { return d[a].id < d[b].id; });
  return x;
}

std::vector<std::vector<int>> split_folds(const std::vector<dpp_matrix_t> &d,
                                          const std::vector<int> &ids, int k) {
  if (ids.size() < 2)
    Helper::halt("DPP two-level nested CV requires at least two subjects in "
                 "each training partition");
  k = std::max(2, std::min(k, (int)ids.size()));
  std::vector<std::vector<int>> f(k);
  auto x = sorted_indices(d, ids);
  for (int i = 0; i < (int)x.size(); i++)
    f[i % k].push_back(x[i]);
  return f;
}

// Deterministic subject-level holdout used only for selecting a booster
// length.  The caller must pass a pool that already excludes its prediction
// subjects (and, for outer fits, the outer test subjects).
void split_internal_subjects(const std::vector<dpp_matrix_t> &d,
                             const std::vector<int> &ids,
                             std::vector<int> *train,
                             std::vector<int> *valid) {
  train->clear();
  valid->clear();
  if (ids.size() < 2) {
    *train = ids;
    return;
  }
  const std::vector<int> x = sorted_indices(d, ids);
  const int nv = std::max(1, std::min((int)x.size() - 1,
                                      (int)x.size() / 5));
  for (int i = 0; i < (int)x.size(); i++)
    (i < nv ? valid : train)->push_back(x[i]);
}

int bounded_iteration(const int best, const int maximum) {
  return std::max(1, std::min(maximum, best > 0 ? best : maximum));
}

std::string lgbm_setting(const lgbm_t &lg, const std::string &key,
                         const std::string &def) {
  std::istringstream in(lg.params);
  std::string x;
  while (in >> x) {
    const size_t p = x.find('=');
    if (p != std::string::npos && x.substr(0, p) == key)
      return x.substr(p + 1);
  }
  return def;
}

void log_booster_settings(const std::string &role, const lgbm_t &lg, int n,
                          int nf, int iterations) {
  logger << "  DPP two-level " << role << ": rows=" << n << " features=" << nf
         << " iterations=" << iterations << "\n";
  logger << "  DPP two-level " << role << " LightGBM constraints:"
         << " min_data_in_leaf="
         << lgbm_setting(lg, "min_data_in_leaf", "20 [default]")
         << " min_sum_hessian_in_leaf="
         << lgbm_setting(lg, "min_sum_hessian_in_leaf", "0.001 [default]")
         << " min_gain_to_split="
         << lgbm_setting(lg, "min_gain_to_split", "0 [default]")
         << " num_leaves=" << lgbm_setting(lg, "num_leaves", "31 [default]")
         << " max_depth=" << lgbm_setting(lg, "max_depth", "-1 [default]")
         << "\n";
  logger << "  DPP two-level " << role << " LightGBM config: "
         << (lg.params.empty() ? "<LightGBM defaults>" : lg.params) << "\n";
}

void log_degenerate_features(const std::string &role, const Eigen::MatrixXd &X,
                             const std::vector<std::string> &labels) {
  int n = 0;
  for (int c = 0; c < X.cols(); c++) {
    bool seen = false, variable = false;
    double first = 0;
    for (int r = 0; r < X.rows(); r++)
      if (std::isfinite(X(r, c))) {
        if (!seen) {
          first = X(r, c);
          seen = true;
        } else if (X(r, c) != first) {
          variable = true;
          break;
        }
      }
    if (!seen || !variable) {
      ++n;
      logger << "  DPP two-level " << role << " degenerate feature: "
             << (c < (int)labels.size() ? labels[c]
                                        : "F" + Helper::int2str(c + 1))
             << (!seen ? " (all missing)" : " (constant)") << "\n";
    }
  }
  if (n)
    logger << "  DPP two-level " << role << ": " << n
           << " degenerate feature(s) in training matrix\n";
}

std::vector<std::string> vector_labels(const param_t &p, const int nf) {
  const int ed = p.has("embedding-dim") ? p.requires_int("embedding-dim") : 128;
  std::vector<std::string> x = dpp_vector::labels(ed, p);
  if ((int)x.size() != nf) {
    x.clear();
    for (int i = 0; i < nf; i++)
      x.push_back("VEC.F" + Helper::int2str(i + 1));
  }
  return x;
}

void emit_metrics(const std::string &level, const std::vector<double> &y,
                  const std::vector<double> &p, const param_t &pmt) {
  const bool bin = pmt.value("outcome", true) == "BINARY";
  const dpp_evaluation::prediction_metrics_t m =
      dpp_evaluation::evaluate_predictions(y, p, bin);
  writer.level(level, "MODEL");
  writer.value("N", (int)m.n);
  writer.value("RMSE", m.rmse);
  if (bin) {
    writer.value("BRIER", m.brier);
    writer.value("LOGLOSS", m.log_loss);
    writer.value("AUC", m.auc);
  } else
    writer.value("R2", m.r2);
  writer.unlevel("MODEL");
}

double mean_outcome(const std::vector<dpp_matrix_t> &d,
                    const std::vector<int> &ids, const param_t &p) {
  double total = 0.0;
  int n = 0;
  for (int i : ids) {
    double y = std::numeric_limits<double>::quiet_NaN();
    if (!cmd_t::pull_ivar(d[i].id, p.requires("phe"), &y) || !std::isfinite(y))
      Helper::halt("missing phenotype " + p.requires("phe") + " for " +
                   d[i].id);
    total += y;
    ++n;
  }
  if (n == 0)
    Helper::halt(
        "DPP two-level: cannot compute a null prediction from an empty fold");
  return total / n;
}

void prediction_losses(const std::vector<double> &y,
                       const std::vector<double> &prediction, bool binary,
                       bool brier_loss, std::vector<double> *losses) {
  losses->clear();
  if (y.size() != prediction.size())
    Helper::halt("DPP evaluation: outcome and prediction vectors have "
                 "different lengths");
  for (size_t i = 0; i < y.size(); ++i)
    if (std::isfinite(y[i]) && std::isfinite(prediction[i])) {
      if (binary && !brier_loss) {
        const double q =
            std::max(1.0e-15, std::min(1.0 - 1.0e-15, prediction[i]));
        losses->push_back(-y[i] * std::log(q) -
                          (1.0 - y[i]) * std::log(1.0 - q));
      } else if (binary)
        losses->push_back((y[i] - prediction[i]) * (y[i] - prediction[i]));
      else
        losses->push_back((y[i] - prediction[i]) * (y[i] - prediction[i]));
    }
}

void write_loss_test(const std::string &prefix,
                     const std::vector<double> &null_loss,
                     const std::vector<double> &model_loss) {
  const dpp_evaluation::paired_test_result_t test =
      dpp_evaluation::paired_t_test(null_loss, model_loss);
  const dpp_evaluation::permutation_test_result_t permutation =
      dpp_evaluation::paired_permutation_test(null_loss, model_loss, 10000,
                                              20260901);
  writer.value(prefix + "_MEAN_LOSS_DIFF", test.mean_difference);
  writer.value(prefix + "_T", test.statistic);
  writer.value(prefix + "_T_P", test.p_value);
  writer.value(prefix + "_PERM_P", permutation.p_value);
}

void emit_null_metrics(const std::string &level, const std::vector<double> &y,
                       const std::vector<double> &prediction,
                       const std::vector<double> &null_prediction,
                       const param_t &pmt) {
  const bool binary = pmt.value("outcome", true) == "BINARY";
  const dpp_evaluation::prediction_metrics_t model =
      dpp_evaluation::evaluate_predictions(y, prediction, binary);
  const dpp_evaluation::prediction_metrics_t null_model =
      dpp_evaluation::evaluate_predictions(y, null_prediction, binary);
  std::vector<double> model_loss, null_loss;
  prediction_losses(y, prediction, binary, false, &model_loss);
  prediction_losses(y, null_prediction, binary, false, &null_loss);

  writer.level(level, "NULL");
  writer.value("N", (int)model.n);
  if (binary) {
    writer.value("MODEL_BRIER", model.brier);
    writer.value("NULL_BRIER", null_model.brier);
    writer.value("DELTA_BRIER", model.brier - null_model.brier);
    std::vector<double> model_brier_loss, null_brier_loss;
    prediction_losses(y, prediction, true, true, &model_brier_loss);
    prediction_losses(y, null_prediction, true, true, &null_brier_loss);
    write_loss_test("BRIER", null_brier_loss, model_brier_loss);
    writer.value("MODEL_LOGLOSS", model.log_loss);
    writer.value("NULL_LOGLOSS", null_model.log_loss);
    writer.value("DELTA_LOGLOSS", model.log_loss - null_model.log_loss);
    write_loss_test("LOGLOSS", null_loss, model_loss);
    writer.value("MODEL_AUC", model.auc);
    writer.value("NULL_AUC", null_model.auc);
    writer.value("DELTA_AUC", model.auc - null_model.auc);
    const dpp_evaluation::auc_chance_test_result_t auc_test =
        dpp_evaluation::auc_chance_permutation_test(y, prediction, 10000,
                                                    20260901);
    writer.value("AUC_CHANCE_DELTA", auc_test.difference_from_chance);
    writer.value("AUC_CHANCE_PERM_P", auc_test.p_value);
  } else {
    writer.value("MODEL_RMSE", model.rmse);
    writer.value("NULL_RMSE", null_model.rmse);
    writer.value("DELTA_RMSE", model.rmse - null_model.rmse);
    writer.value("MODEL_R2", model.r2);
    writer.value("NULL_R2", null_model.r2);
    writer.value("DELTA_R2", model.r2 - null_model.r2);
  }
  if (!binary)
    write_loss_test("RMSE", null_loss, model_loss);
  writer.unlevel("NULL");
}

void configure_lgbm(lgbm_t &lg, const param_t &p, const std::string &level) {
  const std::string specific = level == "Level-1" ? "l1-config" : "l2-config";
  if (p.has(specific))
    lg.load_config(p.value(specific));
  else if (p.has("config"))
    lg.load_config(p.value("config"));
  else if (level == "Level-1")
    lg.params = "objective=regression metric=l2 learning_rate=0.05 "
                "num_leaves=31 min_data_in_leaf=50 min_gain_to_split=0";
  else
    lg.params = "objective=regression metric=l2 learning_rate=0.05 "
                "num_leaves=15 min_data_in_leaf=10 min_gain_to_split=0";
  std::istringstream in(lg.params);
  std::string x, out;
  while (in >> x) {
    const size_t q = x.find('=');
    if (q != std::string::npos &&
        (x.substr(0, q) == "objective" || x.substr(0, q) == "metric"))
      continue;
    if (!out.empty())
      out += " ";
    out += x;
  }
  out += p.value("outcome", true) == "BINARY"
             ? " objective=binary metric=binary_logloss"
             : " objective=regression metric=l2";
  lg.params = out;
}

std::unique_ptr<lgbm_t> train_local(
    const std::vector<dpp_matrix_t> &d, const std::vector<int> &ids,
    const param_t &p, int nf, const std::vector<int> *valid_ids = nullptr,
    int fixed_iterations = 0, int early_rounds = 0) {
  int nr = 0;
  std::vector<int> counts;
  for (int i : ids) {
    int n = 0;
    for (auto &r : d[i].X)
      if (usable(r))
        ++n;
    counts.push_back(n);
    nr += n;
  }
  if (!nr)
    Helper::halt("DPP two-level: no usable Level-1 training rows");
  Eigen::MatrixXd X(nr, nf);
  std::vector<double> y(nr);
  std::vector<float> w(nr);
  int z = 0;
  for (int q = 0; q < (int)ids.size(); q++) {
    int i = ids[q];
    double yy;
    if (!cmd_t::pull_ivar(d[i].id, p.requires("phe"), &yy))
      Helper::halt("no phenotype " + p.requires("phe") + " for " + d[i].id);
    if (p.value("outcome", true) == "BINARY" && (yy != 0.0 && yy != 1.0))
      Helper::halt("binary outcome " + p.requires("phe") +
                   " must be coded 0/1; " + d[i].id + " has " +
                   Helper::dbl2str(yy));
    for (auto &r : d[i].X)
      if (usable(r)) {
        for (int c = 0; c < nf; c++)
          X(z, c) = c < (int)r.size()
                        ? r[c]
                        : std::numeric_limits<double>::quiet_NaN();
        y[z] = yy;
        w[z] = counts[q] ? 1.0f / counts[q] : 0;
        ++z;
      }
  }
  std::vector<std::string> labels = vector_labels(p, nf);
  std::unique_ptr<lgbm_t> lg(new lgbm_t());
  configure_lgbm(*lg, p, "Level-1");
  lg->qt_mode = true;
  lg->attach_training_matrix(X);
  lg->attach_training_qts(y);
  lg->training_weights = w;
  lg->apply_weights(lg->training, &lg->training_weights);
  const int maximum =
      p.has("l1-iterations") ? p.requires_int("l1-iterations") : 100;
  if (valid_ids && !valid_ids->empty()) {
    int vr = 0;
    for (int i : *valid_ids)
      for (const auto &r : d[i].X)
        if (usable(r))
          ++vr;
    if (!vr)
      Helper::halt("DPP two-level: no usable Level-1 validation rows");
    Eigen::MatrixXd VX(vr, nf);
    std::vector<double> vy(vr);
    std::vector<float> vw(vr);
    int vz = 0;
    for (int i : *valid_ids) {
      double yy;
      if (!cmd_t::pull_ivar(d[i].id, p.requires("phe"), &yy))
        Helper::halt("no phenotype " + p.requires("phe") + " for " + d[i].id);
      int count = 0;
      for (const auto &r : d[i].X)
        if (usable(r))
          ++count;
      for (const auto &r : d[i].X)
        if (usable(r)) {
          for (int c = 0; c < nf; c++)
            VX(vz, c) = c < (int)r.size() ? r[c]
                                         : std::numeric_limits<double>::quiet_NaN();
          vy[vz] = yy;
          vw[vz] = count ? 1.0f / count : 0;
          ++vz;
        }
    }
    lg->attach_validation_matrix(VX);
    lg->attach_validation_qts(vy);
    lg->validation_weights = vw;
    lg->apply_weights(lg->validation, &lg->validation_weights);
  }
  lg->n_iterations = fixed_iterations > 0 ? fixed_iterations : maximum;
  lg->early_stopping_rounds =
      (fixed_iterations > 0 || !valid_ids || valid_ids->empty()) ? 0
                                                                  : early_rounds;
  log_degenerate_features("Level-1", X, labels);
  logger << "  DPP two-level Level-1 stopping: "
         << (lg->early_stopping_rounds > 0 ? "enabled" : "disabled")
         << ", max_iterations=" << maximum
         << ", validation_subjects=" << (valid_ids ? valid_ids->size() : 0)
         << ", training_subjects=" << ids.size() << "\n";
  log_booster_settings("Level-1", *lg, nr, nf, lg->n_iterations);
  lg->create_booster(p.has("verbose") && p.yesno("verbose"));
  logger << "  DPP two-level Level-1 selected iterations="
         << bounded_iteration(lg->best_iteration, maximum) << "\n";
  return lg;
}

std::map<std::string, std::vector<double>>
predict_local(lgbm_t &lg, const std::vector<dpp_matrix_t> &d,
              const std::vector<int> &ids, int nf) {
  std::map<std::string, std::vector<double>> out;
  for (int i : ids) {
    Eigen::MatrixXd X(d[i].X.size(), nf);
    for (int r = 0; r < X.rows(); r++)
      for (int c = 0; c < nf; c++)
        X(r, c) = c < (int)d[i].X[r].size()
                      ? d[i].X[r][c]
                      : std::numeric_limits<double>::quiet_NaN();
    Eigen::MatrixXd y = lg.predict(X);
    out[d[i].id] = std::vector<double>(y.rows());
    for (int r = 0; r < y.rows(); r++)
      out[d[i].id][r] = y(r, 0);
  }
  return out;
}

std::map<std::string, std::vector<double>> crossfit_local(
    const std::vector<dpp_matrix_t> &d, const std::vector<int> &ids,
    const param_t &p, int nf, int inner, int early_rounds,
    std::vector<int> *selected) {
  std::map<std::string, std::vector<double>> out;
  selected->clear();
  auto folds = split_folds(d, ids, inner);
  for (const auto &h : folds) {
    std::set<int> hs(h.begin(), h.end());
    std::vector<int> available;
    for (int i : ids)
      if (!hs.count(i))
        available.push_back(i);
    if (early_rounds <= 0) {
      auto lg = train_local(d, available, p, nf);
      auto q = predict_local(*lg, d, h, nf);
      out.insert(q.begin(), q.end());
      continue;
    }
    std::vector<int> fit_ids, valid_ids;
    split_internal_subjects(d, available, &fit_ids, &valid_ids);
    auto probe = train_local(d, fit_ids, p, nf, &valid_ids, 0, early_rounds);
    const int maximum =
        p.has("l1-iterations") ? p.requires_int("l1-iterations") : 100;
    const int it = bounded_iteration(probe->best_iteration, maximum);
    selected->push_back(it);
    auto lg = train_local(d, available, p, nf, nullptr, it, 0);
    auto q = predict_local(*lg, d, h, nf);
    out.insert(q.begin(), q.end());
  }
  return out;
}

int median_iteration(const std::vector<int> &x, const int maximum) {
  if (x.empty())
    return std::max(1, maximum);
  std::vector<double> z;
  for (int v : x)
    z.push_back((double)bounded_iteration(v, maximum));
  return bounded_iteration((int)std::lround(MiscMath::median(z)), maximum);
}

std::unique_ptr<lgbm_t> train_subject(const std::vector<summary_t> &s,
                                      const std::vector<int> &ids,
                                      const std::vector<dpp_matrix_t> &d,
                                      const param_t &p,
                                      const std::vector<std::string> &labels,
                                      const std::vector<int> *valid_ids = nullptr,
                                      int fixed_iterations = 0,
                                      int early_rounds = 0) {
  if (ids.empty())
    Helper::halt("DPP two-level: no subjects for Level-2 training");
  int nf = s[ids[0]].x.size();
  Eigen::MatrixXd X(ids.size(), nf);
  std::vector<double> y(ids.size());
  for (int r = 0; r < (int)ids.size(); r++) {
    for (int c = 0; c < nf; c++)
      X(r, c) = s[ids[r]].x[c];
    if (!cmd_t::pull_ivar(d[ids[r]].id, p.requires("phe"), &y[r]))
      Helper::halt("missing Level-2 phenotype");
    if (p.value("outcome", true) == "BINARY" && (y[r] != 0.0 && y[r] != 1.0))
      Helper::halt("binary outcome " + p.requires("phe") +
                   " must be coded 0/1; " + d[ids[r]].id + " has " +
                   Helper::dbl2str(y[r]));
  }
  std::unique_ptr<lgbm_t> lg(new lgbm_t());
  configure_lgbm(*lg, p, "Level-2");
  lg->qt_mode = true;
  lg->attach_training_matrix(X);
  lg->attach_training_qts(y);
  if ((int)labels.size() == nf)
    lg->set_feature_names(labels);
  const int maximum =
      p.has("l2-iterations") ? p.requires_int("l2-iterations") : 100;
  if (valid_ids && !valid_ids->empty()) {
    Eigen::MatrixXd VX(valid_ids->size(), nf);
    std::vector<double> vy(valid_ids->size());
    for (int r = 0; r < (int)valid_ids->size(); r++) {
      const int i = (*valid_ids)[r];
      VX.row(r) = Eigen::Map<const Eigen::RowVectorXd>(s[i].x.data(), nf);
      if (!cmd_t::pull_ivar(d[i].id, p.requires("phe"), &vy[r]))
        Helper::halt("missing Level-2 validation phenotype");
    }
    lg->attach_validation_matrix(VX);
    lg->attach_validation_qts(vy);
    std::vector<float> vw(valid_ids->size(), 1.0f);
    lg->validation_weights = vw;
    lg->apply_weights(lg->validation, &lg->validation_weights);
  }
  log_degenerate_features("Level-2", X, labels);
  lg->n_iterations = fixed_iterations > 0 ? fixed_iterations : maximum;
  lg->early_stopping_rounds =
      (fixed_iterations > 0 || !valid_ids || valid_ids->empty()) ? 0
                                                                  : early_rounds;
  logger << "  DPP two-level Level-2 stopping: "
         << (lg->early_stopping_rounds > 0 ? "enabled" : "disabled")
         << ", max_iterations=" << maximum
         << ", validation_subjects=" << (valid_ids ? valid_ids->size() : 0)
         << ", training_subjects=" << ids.size() << "\n";
  log_booster_settings("Level-2", *lg, ids.size(), nf, lg->n_iterations);
  lg->create_booster(p.has("verbose") && p.yesno("verbose"));
  logger << "  DPP two-level Level-2 selected iterations="
         << bounded_iteration(lg->best_iteration, maximum) << "\n";
  return lg;
}

double predict_subject(lgbm_t &lg, const summary_t &s) {
  Eigen::MatrixXd X(1, s.x.size());
  for (int c = 0; c < X.cols(); c++)
    X(0, c) = s.x[c];
  return lg.predict(X)(0, 0);
}

std::string join_covars(const std::vector<std::string> &x) {
  std::string s;
  for (const auto &v : x) {
    if (!s.empty())
      s += ",";
    s += v;
  }
  return s;
}

void emit_importance(const std::string &model, lgbm_t &lg,
                     const std::vector<std::string> &labels) {
  const std::vector<double> gain = lg.feature_importance(1);
  const std::vector<double> split = lg.feature_importance(0);
  if (gain.size() != split.size())
    Helper::halt("DPP two-level: inconsistent feature-importance lengths");
  for (int i = 0; i < (int)gain.size(); ++i) {
    const std::string feature =
        i < (int)labels.size() ? labels[i] : "F" + Helper::int2str(i + 1);
        writer.level(model, "IMPORTANCE");
        writer.level(feature, "FEATURE");
    writer.value("GAIN", gain[i]);
    writer.value("SPLIT", split[i]);
        writer.unlevel("FEATURE");
        writer.unlevel("IMPORTANCE");
  }
}

std::string join_groups(const std::vector<std::string> &g) {
  std::string s;
  for (auto &x : g) {
    if (!s.empty())
      s += ",";
    s += x;
  }
  return s;
}

void warn_stage_availability(const std::vector<dpp_matrix_t> &d,
                             const dpp_vector::layout_t &l,
                             const std::vector<std::string> &g) {
  if (!l.context)
    return;
  int nrem = 0, rem = 0;
  for (const dpp_matrix_t &m : d)
    for (const std::vector<double> &row : m.X) {
      if (at(row, l.context_offset + 2) > 0.5 ||
          at(row, l.context_offset + 3) > 0.5 ||
          at(row, l.context_offset + 4) > 0.5)
        ++nrem;
      if (at(row, l.context_offset + 5) > 0.5)
        ++rem;
    }

  if (has(g, "STAGE") && nrem == 0)
    logger << "  *** warning: no retained NREM rows; L2.STAGE.NREM.* features "
              "will be missing\n";
  if (has(g, "STAGE") && rem == 0)
    logger << "  *** warning: no retained REM rows; L2.STAGE.REM.* features "
              "will be missing\n";
  if (has(g, "SCORE") && nrem == 0)
    logger << "  *** warning: no retained NREM rows; L2.SCORE.NREM will be "
              "missing\n";
  if (has(g, "SCORE") && rem == 0)
    logger << "  *** warning: no retained REM rows; L2.SCORE.REM will be "
              "missing\n";
}
} // namespace

bool dpp_twolevel::enabled(const param_t &p) {
  return !(p.has("classic") && p.yesno("classic")) &&
         (p.has("two-level") ? p.yesno("two-level") : true);
}

void dpp_twolevel::fit(param_t &p) {
  if (p.has("fit-spec")) {
    fit_model_set(p);
    return;
  }
  std::vector<std::string> files =
      p.has("files") ? p.strvector("files")
                     : std::vector<std::string>{p.requires("data")};
  std::vector<dpp_matrix_t> d = dpp_io::load_files(files, -1);
  if (d.empty())
    Helper::halt("DPP two-level: no corpus data");
  for (const dpp_matrix_t &m : d)
    if (m.X.empty())
      Helper::halt("DPP two-level: individual " + m.id +
                   " has zero retained vector rows after masking");
  const int nf = d[0].X[0].size(),
            ed = p.has("embedding-dim") ? p.requires_int("embedding-dim") : 128;
  for (auto &m : d)
    for (auto &r : m.X)
      if ((int)r.size() != nf)
        Helper::halt("DPP two-level: inconsistent Level-1 feature count");
  const dpp_vector::layout_t l = dpp_vector::layout(ed, p);
  if (l.raw_offset < 0 || l.raw_offset + ed > nf)
    Helper::halt("DPP two-level requires RAW embedding columns");
  const std::vector<std::string> labels = summary_labels(ed, p), g = groups(p);
  std::vector<int> all(d.size());
  std::iota(all.begin(), all.end(), 0);
  warn_stage_availability(d, l, g);
  int outer = p.has("outer-folds") ? p.requires_int("outer-folds") : 5,
      inner = p.has("inner-folds") ? p.requires_int("inner-folds") : 5;
  if (outer < 2 || inner < 2 || d.size() < 3)
    Helper::halt("DPP two-level requires at least three subjects and "
                 "outer-folds/inner-folds >=2");
  outer = std::min(outer, (int)d.size());
  if ((int)d.size() < 2 * outer)
    outer = (int)d.size();
  auto outer_f = split_folds(d, all, outer);
  const int es = p.has("early-stopping-rounds")
                     ? p.requires_int("early-stopping-rounds")
                     : 0;
  if (es < 0)
    Helper::halt("early-stopping-rounds must be non-negative");
  std::vector<std::string> outer_id;
  std::vector<int> outer_fold;
  std::vector<double> outer_y, outer_p, outer_null;
  std::vector<int> outer_l2_selected;
  for (int f = 0; f < (int)outer_f.size(); f++) {
    auto &test = outer_f[f];
    std::set<int> ts(test.begin(), test.end());
    std::vector<int> tr;
    for (int i : all)
      if (!ts.count(i))
        tr.push_back(i);
    const double fold_null = mean_outcome(d, tr, p);
    std::vector<int> l1_selected;
    std::map<std::string, std::vector<double>> oof =
        crossfit_local(d, tr, p, nf, inner, es, &l1_selected);
    std::vector<summary_t> su(d.size());
    for (int i : tr)
      su[i] = summarize_clean(d[i], oof[d[i].id], l, p);
    const int l1_max =
        p.has("l1-iterations") ? p.requires_int("l1-iterations") : 100;
    const int l1_iter = median_iteration(l1_selected, l1_max);
    std::vector<int> l2_fit, l2_valid;
    int l2_iter = p.has("l2-iterations") ? p.requires_int("l2-iterations") : 100;
    if (es > 0) {
      split_internal_subjects(d, tr, &l2_fit, &l2_valid);
      auto probe = train_subject(su, l2_fit, d, p, labels, &l2_valid, 0, es);
      l2_iter = bounded_iteration(probe->best_iteration, l2_iter);
      outer_l2_selected.push_back(l2_iter);
    }
    auto l2 = train_subject(su, tr, d, p, labels, nullptr, l2_iter, 0);
    auto lg = train_local(d, tr, p, nf, nullptr, l1_iter, 0);
    auto q = predict_local(*lg, d, test, nf);
    for (int i : test) {
      su[i] = summarize_clean(d[i], q[d[i].id], l, p);
      double yy;
      cmd_t::pull_ivar(d[i].id, p.requires("phe"), &yy);
      double pp = predict_subject(*l2, su[i]);
      outer_id.push_back(d[i].id);
      outer_fold.push_back(f + 1);
      outer_y.push_back(yy);
      outer_p.push_back(pp);
      outer_null.push_back(fold_null);
    }
  }
  emit_metrics("OUTER_OOF", outer_y, outer_p, p);
  emit_null_metrics("OUTER_OOF_NULL", outer_y, outer_p, outer_null, p);

  writer.id(".", ".");
  for (int i = 0; i < (int)outer_id.size(); i++) {
    writer.id(outer_id[i], ".");
    writer.level("OUTER_OOF", "OOF");
    writer.value("OBSERVED", outer_y[i]);
    writer.value("PREDICTED", outer_p[i]);
    writer.value("NULL_PREDICTED", outer_null[i]);
    writer.value("OUTER_FOLD", outer_fold[i]);
    writer.unlevel("OOF");
  }
  writer.id(".", ".");

  std::vector<int> l1_selected;
  std::map<std::string, std::vector<double>> oof =
      crossfit_local(d, all, p, nf, inner, es, &l1_selected);
  std::vector<summary_t> su(d.size());
  for (int i : all)
    su[i] = summarize_clean(d[i], oof[d[i].id], l, p);
  const int l1_max =
      p.has("l1-iterations") ? p.requires_int("l1-iterations") : 100;
  const int l1_iter = median_iteration(l1_selected, l1_max);
  std::vector<int> l2_fit, l2_valid;
  int l2_iter = p.has("l2-iterations") ? p.requires_int("l2-iterations") : 100;
  if (es > 0) {
    if (!outer_l2_selected.empty())
      l2_iter = median_iteration(outer_l2_selected, l2_iter);
    else {
      split_internal_subjects(d, all, &l2_fit, &l2_valid);
      auto probe = train_subject(su, l2_fit, d, p, labels, &l2_valid, 0, es);
      l2_iter = bounded_iteration(probe->best_iteration, l2_iter);
    }
  }
  auto l2 = train_subject(su, all, d, p, labels, nullptr, l2_iter, 0);
  auto lg = train_local(d, all, p, nf, nullptr, l1_iter, 0);
  const std::string root = Helper::expand(p.requires("out"));
  lg->save_model(root + ".l1.mod");
  l2->save_model(root + ".l2.mod");
  std::vector<std::string> l1labels = vector_labels(p, nf);
  emit_importance(root, *lg, l1labels);
  emit_importance(root, *l2, labels);
  std::ofstream M(root + ".dpp");
  M << "# DPP model manifest\n# mode=two-level\n# vector=T\n# vector-time="
    << time_mode(p) << "\n# embedding_dim=" << ed << "\n# l1_features=" << nf
    << "\n# level2_features=" << join_groups(g) << "\n";
  if (p.has("covar"))
    M << "# covariates=" << p.value("covar") << "\n";
  M << "# n_features=" << labels.size() << "\n# feature_names_begin\n";
  for (auto &x : labels)
    M << x << "\n";
  M.close();
  logger << "  wrote two-level DPP bundle: " << root << ".l1.mod and " << root
         << ".l2.mod (" << labels.size() << " Level-2 features)\n";
}

// Fit a set of two-level models using common subject folds.  The existing
// two-level machinery is deliberately reused: a sleep model shares one
// cross-fitted Level-1 model, while sleep=F models use only their requested
// subject covariates at Level 2.
void dpp_twolevel::fit_model_set(param_t &p) {
  dpp_fit_t spec_reader(p);
  const auto specs = spec_reader.model_specs;
  const auto contrasts = spec_reader.contrast_specs;
  std::vector<std::string> files =
      p.has("files") ? p.strvector("files")
                     : std::vector<std::string>{p.requires("data")};
  std::vector<dpp_matrix_t> d = dpp_io::load_files(files, -1);
  if (d.empty())
    Helper::halt("DPP two-level: no corpus data");
  for (const auto &m : d)
    if (m.X.empty())
      Helper::halt("DPP two-level: individual " + m.id +
                   " has zero retained vector rows after masking");
  const int nf = d[0].X[0].size(),
            ed = p.has("embedding-dim") ? p.requires_int("embedding-dim") : 128;
  for (const auto &m : d)
    for (const auto &r : m.X)
      if ((int)r.size() != nf)
        Helper::halt("DPP two-level: inconsistent Level-1 feature count");
  const dpp_vector::layout_t l = dpp_vector::layout(ed, p);
  if (l.raw_offset < 0 || l.raw_offset + ed > nf)
    Helper::halt("DPP two-level requires RAW embedding columns");
  std::vector<int> all(d.size());
  std::iota(all.begin(), all.end(), 0);
  int outer = p.has("outer-folds") ? p.requires_int("outer-folds") : 5;
  int inner = p.has("inner-folds") ? p.requires_int("inner-folds") : 5;
  const int es = p.has("early-stopping-rounds")
                     ? p.requires_int("early-stopping-rounds")
                     : 0;
  if (es < 0)
    Helper::halt("early-stopping-rounds must be non-negative");
  if (outer < 2 || inner < 2 || d.size() < 3)
    Helper::halt("DPP two-level fit-spec requires at least three subjects and "
                 "outer-folds/inner-folds >=2");
  outer = std::min(outer, (int)d.size());
  if ((int)d.size() < 2 * outer)
    outer = (int)d.size();
  std::vector<std::map<std::string, double>> pred(specs.size());
  std::vector<std::map<std::string, double>> null_pred(specs.size());
  std::vector<std::vector<int>> l2_selected(specs.size());
  auto model_param = [&](const dpp_fit_t::model_spec_t &m) {
    param_t q = p;
    q.add("covar", join_covars(m.covariates));
    if (!m.sleep)
      q.add("level2-features", "");
    return q;
  };
  auto fit_one =
      [&](const dpp_fit_t::model_spec_t &m, const std::vector<int> &tr,
          const std::vector<int> &te,
          const std::map<std::string, std::vector<double>> &cross_oof,
          const std::map<std::string, std::vector<double>> &test_local,
          bool write_models, const std::string &root) {
        param_t q = model_param(m);
        std::vector<summary_t> su(d.size());
        if (m.sleep)
          for (int i : tr)
            su[i] = summarize_clean(d[i], cross_oof.at(d[i].id), l, q);
        else
          for (int i : tr)
            su[i] = summarize_clean(d[i], {}, l, q);
        std::vector<std::string> labs = summary_labels(ed, q);
        int l2_iter = q.has("l2-iterations") ? q.requires_int("l2-iterations") : 100;
        if (es > 0) {
          if (write_models && !l2_selected[&m - &specs[0]].empty())
            l2_iter = median_iteration(l2_selected[&m - &specs[0]], l2_iter);
          else {
            std::vector<int> l2_fit, l2_valid;
            split_internal_subjects(d, tr, &l2_fit, &l2_valid);
            auto probe = train_subject(su, l2_fit, d, q, labs, &l2_valid, 0, es);
            l2_iter = bounded_iteration(probe->best_iteration, l2_iter);
            if (!write_models)
              l2_selected[&m - &specs[0]].push_back(l2_iter);
          }
        }
        auto l2 = train_subject(su, tr, d, q, labs, nullptr, l2_iter, 0);
        for (int i : te) {
          summary_t s =
              m.sleep ? summarize_clean(d[i], test_local.at(d[i].id), l, q)
                      : summarize_clean(d[i], {}, l, q);
          double y;
          if (!cmd_t::pull_ivar(d[i].id, p.requires("phe"), &y))
            Helper::halt("missing phenotype " + p.requires("phe") + " for " +
                         d[i].id);
          double z = predict_subject(*l2, s);
          pred[&m - &specs[0]][d[i].id] = z;
        }
        if (write_models) {
          l2->save_model(root + ".l2.mod");
          emit_importance(m.id, *l2, labs);
          std::ofstream M(root + ".dpp");
          M << "# mode=two-level\n# model_id=" << m.id
            << "\n# sleep=" << (m.sleep ? "T" : "F")
            << "\n# covariates=" << join_covars(m.covariates)
            << "\n# missing=LIGHTGBM_NATIVE\n# feature_names_begin\n";
          for (const auto &x : labs)
            M << x << "\n";
          M.close();
        }
      };
  auto folds = split_folds(d, all, outer);
  for (const auto &te : folds) {
    std::set<int> ts(te.begin(), te.end());
    std::vector<int> tr;
    for (int i : all)
      if (!ts.count(i))
        tr.push_back(i);
    const double fold_null = mean_outcome(d, tr, p);
    std::map<std::string, std::vector<double>> coof;
    bool need_sleep = false;
    for (const auto &m : specs)
      need_sleep |= m.sleep;
    std::vector<int> l1_selected;
    if (need_sleep)
      coof = crossfit_local(d, tr, p, nf, inner, es, &l1_selected);
    std::map<std::string, std::vector<double>> test_local;
    if (need_sleep) {
      const int l1_max =
          p.has("l1-iterations") ? p.requires_int("l1-iterations") : 100;
      auto lg = train_local(d, tr, p, nf, nullptr,
                            median_iteration(l1_selected, l1_max), 0);
      test_local = predict_local(*lg, d, te, nf);
    }
    for (const auto &m : specs) {
      fit_one(m, tr, te, coof, test_local, false, "");
      for (int i : te)
        null_pred[&m - &specs[0]][d[i].id] = fold_null;
    }
  }
  for (int mi = 0; mi < (int)specs.size(); mi++) {
    std::vector<double> yy, pp, nn;
    for (const auto &z : pred[mi]) {
      double y;
      if (cmd_t::pull_ivar(z.first, p.requires("phe"), &y)) {
        yy.push_back(y);
        pp.push_back(z.second);
        nn.push_back(null_pred[mi].at(z.first));
      }
    }
    emit_metrics(specs[mi].id, yy, pp, p);
    emit_null_metrics(specs[mi].id + "_NULL", yy, pp, nn, p);
  }
  // Preserve subject-level outer-OOF records for every model in a fit-spec
  // run.  These are paired by subject for contrasts and are written through
  // the normal database stream, just like the single-model path.
  writer.id(".", ".");
  for (int mi = 0; mi < (int)specs.size(); mi++)
    for (const auto &z : pred[mi]) {
      double y;
      if (!cmd_t::pull_ivar(z.first, p.requires("phe"), &y))
        continue;
      writer.id(z.first, ".");
      writer.level(specs[mi].id, "OOF");
      writer.value("OBSERVED", y);
      writer.value("PREDICTED", z.second);
      writer.value("NULL_PREDICTED", null_pred[mi].at(z.first));
      writer.unlevel("OOF");
    }
  writer.id(".", ".");
  for (const auto &c : contrasts) {
    int bi = -1, ai = -1;
    for (int i = 0; i < (int)specs.size(); ++i) {
      if (specs[i].id == c.base)
        bi = i;
      if (specs[i].id == c.add)
        ai = i;
    }
    if (bi < 0 || ai < 0)
      Helper::halt("unknown fit-spec contrast model");
    std::vector<double> yb, pb, ya, pa;
    for (const auto &z : pred[bi])
      if (pred[ai].count(z.first)) {
        double y;
        if (cmd_t::pull_ivar(z.first, p.requires("phe"), &y)) {
          yb.push_back(y);
          pb.push_back(z.second);
          ya.push_back(y);
          pa.push_back(pred[ai][z.first]);
        }
      }
    const bool binary = p.value("outcome", true) == "BINARY";
    const dpp_evaluation::prediction_metrics_t mb =
        dpp_evaluation::evaluate_predictions(yb, pb, binary);
    const dpp_evaluation::prediction_metrics_t ma =
        dpp_evaluation::evaluate_predictions(ya, pa, binary);
    std::vector<double> base_loss, add_loss;
    for (size_t i = 0; i < yb.size(); ++i) {
      if (binary) {
        const double bp = std::max(1.0e-15, std::min(1.0 - 1.0e-15, pb[i]));
        const double ap = std::max(1.0e-15, std::min(1.0 - 1.0e-15, pa[i]));
        base_loss.push_back(-yb[i] * std::log(bp) -
                            (1.0 - yb[i]) * std::log(1.0 - bp));
        add_loss.push_back(-ya[i] * std::log(ap) -
                           (1.0 - ya[i]) * std::log(1.0 - ap));
      } else {
        base_loss.push_back((yb[i] - pb[i]) * (yb[i] - pb[i]));
        add_loss.push_back((ya[i] - pa[i]) * (ya[i] - pa[i]));
      }
    }
    const dpp_evaluation::paired_test_result_t paired_t =
        dpp_evaluation::paired_t_test(base_loss, add_loss);
    const dpp_evaluation::permutation_test_result_t paired_permutation =
        dpp_evaluation::paired_permutation_test(base_loss, add_loss);
    writer.level(c.id, "CONTRAST");
    writer.value("N", (int)mb.n);
    writer.value("BASE_RMSE", mb.rmse);
    writer.value("ADD_RMSE", ma.rmse);
    writer.value("DELTA_RMSE", mb.rmse - ma.rmse);
    if (binary) {
      writer.value("BASE_BRIER", mb.brier);
      writer.value("ADD_BRIER", ma.brier);
      writer.value("DELTA_BRIER", mb.brier - ma.brier);
      writer.value("BASE_LOGLOSS", mb.log_loss);
      writer.value("ADD_LOGLOSS", ma.log_loss);
      writer.value("DELTA_LOGLOSS", mb.log_loss - ma.log_loss);
      writer.value("BASE_AUC", mb.auc);
      writer.value("ADD_AUC", ma.auc);
    } else {
      writer.value("BASE_R2", mb.r2);
      writer.value("ADD_R2", ma.r2);
      writer.value("DELTA_R2", ma.r2 - mb.r2);
    }
    writer.value("MEAN_LOSS_DIFF", paired_t.mean_difference);
    writer.value("PAIRED_T", paired_t.statistic);
    writer.value("PAIRED_T_P", paired_t.p_value);
    writer.value("PERM_P", paired_permutation.p_value);
    writer.unlevel("CONTRAST");
  }
  for (const auto &m : specs) {
    std::vector<int> none;
    std::map<std::string, std::vector<double>> coof;
    bool need = m.sleep;
    std::vector<int> l1_selected;
    if (need)
      coof = crossfit_local(d, all, p, nf, inner, es, &l1_selected);
    const std::string root = Helper::expand(p.requires("out") + "." + m.id);
    if (need) {
      const int l1_max =
          p.has("l1-iterations") ? p.requires_int("l1-iterations") : 100;
      auto lg = train_local(d, all, p, nf, nullptr,
                            median_iteration(l1_selected, l1_max), 0);
      lg->save_model(root + ".l1.mod");
      emit_importance(m.id, *lg, vector_labels(p, nf));
    }
    auto test = coof;
    fit_one(m, all, none, coof, test, true, root);
  }
  logger << "  wrote " << specs.size() << " two-level DPP model(s) and "
         << contrasts.size() << " contrast(s)\n";
}

void dpp_twolevel::apply(edf_t &edf, param_t &p, const dpp_matrix_t &m) {
  const std::string root = Helper::expand(p.requires("model"));
  const bool has_level1 = Helper::fileExists(root + ".l1.mod");
  const bool dynamic = p.has("dynamic");

  if (dynamic && p.value("dynamic").empty())
    Helper::halt("DPP two-level model=: dynamic= requires a signal label");
  if (dynamic && !has_level1)
    Helper::halt("DPP two-level model=: dynamic= requires a Level-1 model");

  lgbm_t l2;
  l2.qt_mode = true;
  l2.load_model(root + ".l2.mod");

  const int ed = p.has("embedding-dim") ? p.requires_int("embedding-dim") : 128;
  const dpp_vector::layout_t lay = dpp_vector::layout(ed, p);
  if (m.X.empty())
    return;
  if (lay.raw_offset < 0 || lay.raw_offset + ed > (int)m.X[0].size())
    Helper::halt("DPP two-level apply: embedding layout mismatch");

  std::vector<double> score;
  if (has_level1) {
    lgbm_t l1;
    l1.qt_mode = true;
    l1.load_model(root + ".l1.mod");

    Eigen::MatrixXd X(m.X.size(), m.X[0].size());
    for (int r = 0; r < X.rows(); ++r)
      for (int c = 0; c < X.cols(); ++c)
        X(r, c) = m.X[r][c];

    const Eigen::MatrixXd z = l1.predict(X);
    score.resize(z.rows());
    for (int r = 0; r < z.rows(); ++r)
      score[r] = z(r, 0);
  }

  // Level 1 is a row/window-level calculation.  Level 2 summarizes those
  // values once for the subject and produces a single subject-level result
  // for the standard output stream.
  const summary_t s = summarize_clean(m, score, lay, p);
  const double pred = predict_subject(l2, s);

  writer.level(root, "APPLY");
  writer.value("LEVEL2_PRED", pred);
  writer.value("N_WINDOWS", (int)score.size());
  writer.unlevel("APPLY");

  if (!dynamic) {
    logger << "  DPP two-level: Level-2 subject prediction " << pred
           << " reported; no EDF signal requested\n";
    return;
  }

  if (edf.is_actually_discontinuous())
    Helper::halt("DPP two-level model=: cannot attach to discontinuous EDF");

  double step_sec = p.has("step") ? p.requires_dbl("step") : 30.0;
  if (!p.has("step"))
    step_sec = m.time_sec.size() > 1 ? m.time_sec[1] - m.time_sec[0]
                                    : edf.header.record_duration;
  if (step_sec <= 0)
    Helper::halt("DPP two-level dynamic application found invalid time points");

  const double n_spr_d = edf.header.record_duration / step_sec;
  const int n_spr = (int)std::lround(n_spr_d);
  if (n_spr < 1 || std::fabs(n_spr_d - n_spr) > 1e-6)
    Helper::halt("DPP two-level model=: record duration is not an integer multiple of step="
                 + Helper::dbl2str(step_sec));

  const int ne_total = edf.header.nr * n_spr;
  double pmin = std::numeric_limits<double>::infinity();
  double pmax = -std::numeric_limits<double>::infinity();
  for (const double x : score)
    if (std::isfinite(x)) {
      pmin = std::min(pmin, x);
      pmax = std::max(pmax, x);
    }
  if (pmin > pmax) {
    pmin = 0.0;
    pmax = 1.0;
  }
  const double span = pmax - pmin;
  const double sentinel = pmin - (span > 0.0 ? span : 1.0) - 1.0;
  std::vector<double> data(ne_total, sentinel);
  for (int r = 0; r < (int)m.time_sec.size() && r < (int)score.size(); ++r) {
    if (!std::isfinite(score[r]))
      continue;
    const int slot = (int)std::lround(m.time_sec[r] / step_sec) - 1;
    if (slot >= 0 && slot < ne_total)
      data[slot] = score[r];
  }

  const std::string label = p.value("dynamic");
  if (edf.header.has_signal(label))
    Helper::halt("DPP two-level: signal already exists: " + label);

  edf.add_signal(label, -(double)n_spr, data, sentinel, pmax);
  logger << "  DPP two-level: attached Level-1 dynamic signal " << label
         << " (" << ne_total << " samples, " << score.size()
         << " windows); Level-2 prediction reported separately\n";
}

#endif
