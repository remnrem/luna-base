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

// POPS-SFM: in-process SleepFM embedding features for POPS level-1.
//
// Given a comma-delimited candidate channel list (1-10 channels; a recording
// may legitimately supply only a subset -- SleepFM's channel-agnostic BAS
// branch tolerates this via its own channel_mask), this runs the SleepFM
// ONNX model once over the whole recording (one forward pass per 300s -- or
// 'step' -- window, non-overlapping by default) and derives four per-epoch
// views from that single pass:
//
//   tokens       6 x 128 raw per-5s sequence tokens for this 30s epoch
//                (already fully self-attended across the whole 300s window)
//   epoch-pool   128-d mean of just this epoch's 6 tokens (free, no extra
//                inference)
//   window-pool  128-d SleepFM's own learned pooled_embedding for the whole
//                300s window (identical for every epoch sharing a window;
//                not a simple mean of the tokens -- a separate learned
//                attention-pooling layer)
//   position     2 columns: seconds of real past/future context available
//                to this epoch within its host window (0 for epochs fully
//                centered; useful for the model to learn to discount
//                boundary epochs, since SleepFM itself has no notion of a
//                POPS epoch)
//
// An epoch only gets values if all 6 of its token slots fall inside a
// single complete window; otherwise its row is left as NaN (consistent with
// POPS-CODA's existing convention of NaN-filling flagged/unavailable rows).

#if defined(HAS_LGBM) && defined(HAS_ORT)

#ifndef __LUNA_POPS_SFM_H__
#define __LUNA_POPS_SFM_H__

#include "stats/Eigen/Dense"
#include <string>
#include <vector>

struct edf_t;
struct signal_list_t;

struct pops_sfm_result_t {
  Eigen::MatrixXd tokens;       // ne x (tokens_per_epoch * 128), NaN where unavailable
  Eigen::MatrixXd epoch_pool;   // ne x 128
  Eigen::MatrixXd window_pool;  // ne x 128
  Eigen::MatrixXd position;     // ne x 2  (past_context_sec, future_context_sec)
};

namespace pops_sfm_t
{
  // Computes all four SleepFM-derived views in a single pass. Result row i
  // corresponds to physical recording epoch E[i] (i.e. edf.timeline.epoch(
  // E[i])), *not* row index i -- callers must pass the same E used to index
  // X1, since level1() may call this before any bad-row pruning re-slots
  // rows away from a trivial 0..ne-1 epoch numbering. Epochs are assumed to
  // be 30 seconds long (POPS' standard epoch length).
  //   signals  already-resolved channels (1-10) to use -- callers must
  //            resolve names/aliases themselves (e.g. via the same
  //            CH-declaration + pops_opt_t::aliases path the general
  //            channel loop in level1() already uses) rather than passing
  //            raw candidate names for this function to re-resolve; this
  //            struct only reads from already-matched signal slots
  //   modality BAS, RESP, EKG or EMG (channel-count cap + SleepFM branch)
  //   path,lib SleepFM ONNX model location, as for the ORT command
  //   step_sec window step in seconds (300 = non-overlapping, block-aligned)
  pops_sfm_result_t run( edf_t & edf ,
                         const signal_list_t & signals ,
                         const std::string & modality ,
                         const std::string & path ,
                         const std::string & lib ,
                         int step_sec ,
                         const std::vector<int> & E );
}

#endif
#endif
