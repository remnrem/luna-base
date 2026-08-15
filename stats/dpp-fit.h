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

#include <string>
#include <vector>
#include <set>
#include <map>

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

struct dpp_fit_t
{
  dpp_fit_t( param_t & param );

  // top-level driver
  void fit();

  // individual steps, public for direct unit testing (mirrors the
  // dpp_specs_t::read()-style direct-struct-testing pattern already used
  // for stage 2's spec grammar)
  void load_corpus();
  void build_feature_labels();
  void attach_phenotypes();
  void flatten_and_split();
  void apply_row_weights();
  void train_booster();
  void save_bundle();

  param_t & param;

  std::vector<dpp_matrix_t> individuals;
  std::vector<std::string> full_labels;   // per-column labels, mat.X order
  int n_features;
  std::string phe_label;
  std::set<std::string> validation_ids;
  std::map<std::string,double> phenotype;   // individual ID -> outcome, filled by attach_phenotypes()

  Eigen::MatrixXd Xtrain, Xvalid;
  std::vector<double> ytrain, yvalid;

  // per-individual block starts (row index within Xtrain/Xvalid where that
  // individual's rows begin) and 1/n_i weight table, keyed the same way --
  // mirrors pops_t::attach_indiv_weights (pops/pops.cpp)
  std::vector<uint64_t> istart_train, istart_valid;
  std::map<uint64_t,float> wtable_train, wtable_valid;

  lgbm_t lgbm;

  std::string out_root;
};

namespace dpp_fit {

  // per-EDF apply path: load bundle, validate the caller's feature labels
  // against the manifest, predict, reshape into an EDF signal (all-or-
  // nothing -- halts on any mismatch rather than silently degrading)
  void apply( edf_t & edf , param_t & param , const dpp_specs_t & specs , const dpp_matrix_t & mat );

  // shared by dpp_fit_t::build_feature_labels() and dpp_fit::apply(): the
  // full, ordered, per-column feature-label list for a spec (same suffix
  // convention -- V/V1..Vk -- the interactive writer already emits in
  // stats/dpp.cpp, so bundle feature names always match destrat output)
  std::vector<std::string> feature_labels( const dpp_specs_t & specs );

}

#endif
#endif
