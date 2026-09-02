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

#ifndef __LUNA_DPP_FIT_H__
#define __LUNA_DPP_FIT_H__

#include "lgbm/lgbm.h"
#include "stats/dpp-io.h"

#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

struct param_t;
struct edf_t;
struct dpp_specs_t;

// Stage 3: simple pooled (non-stage-conditioned) LightGBM regression --
// train a single booster from a concatenated DPP data= corpus against a
// person-level phenotype (--dpp-fit, cohort-level, no attached EDF), then
// re-apply it to a recording as a new signal (DPP model=, per-EDF).
//
// Reuses, unchanged: lgbm_t (lgbm/lgbm.h), dpp_io::load_files (stage 2),
// cmd_t::pull_ivar (phenotype lookup), dpp_specs_t (feature-label
// derivation), edf_t::add_signal. Deliberately no .ranges/winsorization
// (LightGBM's tree splits are scale-invariant; POPS needs ranges only for
// its SVD step, which DPP has no equivalent of -- see the implementation
// plan's Stage 3 SS3), no stage-conditioning, no grouped CV/OOF (later
// stages). See the implementation plan's "Stage 3" section for the full
// design and reused-primitive citations.

struct dpp_fit_t {
  dpp_fit_t(param_t &param);
  ~dpp_fit_t(); // explicit, out-of-line -- see the definition's comment
                // (dpp-fit.cpp)

  // top-level driver
  void fit();

  // individual steps, public for direct unit testing (mirrors the
  // dpp_specs_t::read()-style direct-struct-testing pattern already used
  // for stage 2's spec grammar)
  void load_corpus();
  void build_feature_labels();
  void attach_phenotypes();
  void attach_covariates();
  void load_hypno_corpus(); // stage 4, only if stage_conditioned
  void flatten_and_split();
  void apply_row_weights();
  void train_booster();        // pooled (stage 3) path
  void train_stage_boosters(); // stage-conditioned (stage 4) path
  void save_bundle();
  void parse_fit_spec();
  void fit_model_set();

  // writes <root>.importance (pooled) / <root>.<stage>.importance
  // (stage-conditioned): full_labels paired with lg's gain- and
  // split-based feature importance, sorted by gain descending. Called
  // from save_bundle(), after each booster's save_model().
  void write_importance_table(lgbm_t &lg, const std::string &file);

  // stage-conditioned only: <root>.importance, gain/split summed across
  // all stage boosters (there's no single shared model to report otherwise)
  void write_aggregate_importance_table(const std::string &file);

  // stage 5: individual-level K-fold CV/OOF, purely an evaluation layer --
  // the final saved bundle (above) is always trained on the entire corpus,
  // regardless of whether cross_validate() ran. assign_folds() builds
  // fold_assignment (sorted-round-robin default, or fold-file= override);
  // cross_validate() repeatedly re-drives flatten_and_split()/
  // train_booster()/train_stage_boosters() with validation_ids set to each
  // fold in turn, collecting one out-of-fold prediction per usable row into
  // oof_rows, then restores validation_ids and calls write_oof().
  void assign_folds();
  void cross_validate();
  void write_oof();

  // subject-level OOF summary: the evaluation that actually matters, since
  // the phenotype is per-subject, not per-row -- row-level r/RMSE (above)
  // mixes in within-subject autocorrelation and unequal per-subject row
  // counts. aggregate_oof_by_subject() means each subject's rows down to
  // one (y_true, mean y_pred) pair (oof_subjects); write_oof_subject()
  // computes Pearson/Spearman/RMSE over those pairs and writes them
  void aggregate_oof_by_subject();
  void write_oof_subject();

  // stage 4: weighted regression of true phenotype on one stage's raw
  // (uncalibrated) validation predictions, restricted to validation rows
  // where that stage's PP_s > 0.5 (GLM has no weighted-least-squares path
  // -- this hard-threshold selection is the pragmatic alternative, see the
  // implementation plan's Stage 4 SS3). Returns {a_s,b_s} = {1,0} if no
  // validation set, too few qualifying rows, or GLM itself reports an
  // invalid/degenerate fit (GLM::beta() is declared but has no
  // implementation anywhere in the tree -- use GLM::display() instead,
  // which is implemented and does the same job)
  std::pair<double, double> calibrate(int stage_idx);

  param_t &param;

  std::vector<dpp_matrix_t> individuals;
  std::vector<std::string> full_labels; // per-column labels, mat.X order
  int n_features;
  std::string phe_label;
  // Optional subject-level covariates requested with covar=.  Loaded and
  // retained for model-set fitting; missing values remain NaN for LightGBM.
  std::vector<std::string> covariate_labels;
  std::set<std::string> validation_ids;
  std::map<std::string, double>
      phenotype; // individual ID -> outcome, filled by attach_phenotypes()
  std::map<std::string, std::map<std::string, double>>
      covariate; // ID -> label -> value (NaN = missing)

  struct model_spec_t {
    std::string id;
    std::vector<std::string> covariates;
    bool sleep;
  };
  struct contrast_spec_t {
    std::string id;
    std::string base;
    std::string add;
  };
  bool model_set_mode;
  std::vector<model_spec_t> model_specs;
  std::vector<contrast_spec_t> contrast_specs;

  Eigen::MatrixXd Xtrain, Xvalid;
  std::vector<double> ytrain, yvalid;

  // per-individual block starts (row index within Xtrain/Xvalid where that
  // individual's rows begin) and 1/n_i weight table, keyed the same way --
  // mirrors pops_t::attach_indiv_weights (pops/pops.cpp). Used by the
  // pooled (stage 3) path's apply_row_weights()/add_block_weights()
  std::vector<uint64_t> istart_train, istart_valid;
  std::map<uint64_t, float> wtable_train, wtable_valid;

  // stage 5: per-row (individual ID, time_sec) aligned with Xvalid, filled
  // by flatten_and_split() alongside Xvalid -- needed only to attribute an
  // out-of-fold prediction back to a specific row in cross_validate()'s
  // .oof output; Xtrain has no equivalent since it's never written out
  std::vector<std::pair<std::string, double>> valid_row_key;
  std::vector<std::pair<std::string, double>> train_row_key;

  lgbm_t lgbm; // pooled (stage 3) booster

  // stage 4: presence of hypno=/hypno-files= (checked in the constructor)
  bool stage_conditioned;
  bool vector_mode;
  bool two_level;
  bool hypno_three_state;
  std::vector<std::string> stage_labels; // "W","R","N1","N2","N3" or 3-state
  std::vector<dpp_matrix_t>
      hypno_individuals; // loaded from hypno=/hypno-files=

  // per-row hypnodensity, aligned row-for-row with Xtrain/Xvalid (built by
  // flatten_and_split() alongside them, from the same row selection/order)
  Eigen::MatrixXd Htrain, Hvalid;

  // per-row (1/n_i) individual weight, aligned with Xtrain/Xvalid -- unlike
  // wtable_train/istart_train above (block-constant, for add_block_weights),
  // stage 4 needs this expanded per-row since it's multiplied by a per-row
  // PP_s value that varies within an individual's own rows
  std::vector<float> indiv_weight_train, indiv_weight_valid;

  std::vector<lgbm_t> stage_lgbm; // one per stage, parallel to stage_labels
  std::vector<double> calib_a,
      calib_b; // one per stage, parallel to stage_labels

  std::string out_root;

  // stage 5: K-fold CV/OOF (evaluation-only; see cross_validate() above)
  int n_folds; // 0 = disabled
  bool save_fold_models;
  std::map<std::string, int> fold_assignment; // individual ID -> fold index

  struct oof_row_t {
    std::string id;
    double time_sec;
    double y_true;
    double y_pred;
    int fold;
  };
  std::vector<oof_row_t> oof_rows;

  struct oof_subject_t {
    std::string id;
    double y_true;
    double y_pred_mean;
    int n_rows;
  };
  std::vector<oof_subject_t> oof_subjects;

  // stage 5c: CV-derived calibration/iteration-count for the *final*,
  // full-corpus bundle -- reuses signal the K-fold loop already produces
  // rather than holding out any fresh data from the final fit (which
  // would otherwise mean some subjects get in-sample-fit predictions and
  // others genuinely out-of-sample ones when the same bundle is later
  // applied back across the whole training cohort -- rejected for exactly
  // that reason). See cross_validate()/fit() for how these get populated
  // and consumed.
  int early_stopping_rounds; // 0 = disabled; only meaningful when a fold's own
                             // validation set is attached (a no-op for the
                             // final, no-validation fit)

  // set by cross_validate() (median across folds), consumed by
  // train_booster()/train_stage_boosters() in place of iterations=/default
  // -- 0 / empty = unset, meaning "use iterations=/default as before"
  int forced_n_iterations;
  std::vector<int> forced_stage_n_iterations;

  // set by cross_validate() (mean across folds' own calibrate() calls),
  // applied by fit() to override the final fit's own calibrate() calls
  // (which always return {1,0}, since the final fit has no validation set)
  std::vector<double> cv_calib_a, cv_calib_b;
};

namespace dpp_fit {

// per-EDF apply path: load bundle, validate the caller's feature labels
// against the manifest, predict, reshape into an EDF signal (all-or-
// nothing -- halts on any mismatch rather than silently degrading).
// 'hmat' is the per-row hypnodensity (stats/dpp.cpp's hypno_mat, row-
// aligned with 'mat'), non-NULL whenever the caller supplied hypno-
// prefix=; required (halts if NULL) for a stage-conditioned bundle,
// ignored for a pooled one.
void apply(edf_t &edf, param_t &param, const dpp_specs_t &specs,
           const dpp_matrix_t &mat, const dpp_matrix_t *hmat);

// shared by dpp_fit_t::build_feature_labels() and dpp_fit::apply(): the
// full, ordered, per-column feature-label list for a spec (same suffix
// convention -- V/V1..Vk -- the interactive writer already emits in
// stats/dpp.cpp, so bundle feature names always match destrat output)
std::vector<std::string> feature_labels(const dpp_specs_t &specs);

} // namespace dpp_fit

#endif
#endif
