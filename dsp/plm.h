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

#ifndef __LUNA_PLM_H__
#define __LUNA_PLM_H__

//
// LM : leg-movement / periodic-leg-movement detection following the
// WASM 2016 (Ferri et al.) event grammar, with two interchangeable
// low-level EMG event detectors:
//
//   method=wasm     : literal fixed micro-volt excursions above a dynamic
//                     resting baseline
//   method=adaptive : robust local-baseline + local-robust-scale detector
//                     (scale-relative), sharing the identical downstream
//                     WASM LM/CLM/bilateral/PLM classification
//
// The processing is deliberately split into separable stages so the
// low-level detector is independent of the WASM event grammar:
//
//   per-channel LM detection            (signal domain)
//     -> bilateral L/R combination      (MODE=LR only)  [pure grammar]
//     -> CLM classification                             [pure grammar]
//     -> PLM sequence detection                         [pure grammar]
//     -> annotations + summaries + verbose event output
//
// The three [pure grammar] steps operate only on the abstract event
// structures below and are directly unit-testable (see tests/tests.cpp).
//

#include <cstdint>
#include <string>
#include <vector>

struct edf_t;
struct param_t;

namespace plm {

  // --- low-level detector -------------------------------------------------

  enum detector_t { DET_WASM , DET_ADAPTIVE };

  // laterality of an event/component
  enum side_t { SIDE_L , SIDE_R , SIDE_B , SIDE_LEG };

  const char * side_str( side_t s );

  // --- parameters ---------------------------------------------------------

  struct plm_param_t {

    plm_param_t() { set_defaults(); }
    void set_defaults();
    void parse( param_t & param );

    detector_t detector;         // DET_WASM | DET_ADAPTIVE

    // literal-WASM amplitude thresholds (micro-volts)
    double onset;                // onset excursion (default 8 uV)
    double offset;               // offset excursion (default 2 uV)
    double morph;                // morphology excursion (== offset, 2 uV)

    // adaptive (robust standardized) thresholds
    double k_on;                 // onset in robust-scale units (default 6)
    double k_off;                // offset in robust-scale units (default 2)
    double k_morph;              // morphology in robust-scale units (== k_off)

    // shared event-state-machine timing (seconds)
    double offset_dur;           // continuous below-offset confirmation (0.5)
    double min_dur;              // minimum LM duration (0.5)
    double morph_win;            // morphology sliding window (0.5)
    double max_clm_dur;          // CLM upper duration bound (10)

    // bilateral combination
    double bilateral_gap;        // max L/R gap to connect (strict <) (0.5)
    int    bilateral_max_components; // (4)
    double bilateral_max_dur;    // combined-duration CLM ceiling (15)

    // PLM periodicity (seconds)
    double imi_min;              // (10)
    double imi_max;              // (90)
    int    min_seq;              // consecutive CLMs for a sequence (4)

    // dynamic iterative baseline estimator
    double baseline_window;      // (60 s)
    double baseline_step;        // grid step (1 s)
    int    baseline_iter;        // exclusion/recompute rounds (3)
    double baseline_exclude_k;   // provisional-excursion cut (3)
    double baseline_exclude_pad; // dilation (0.5 s)

    // QC
    double qc_base_high;         // resting baseline > this uV -> QC flag (16)

    // annotation prefix
    std::string prefix;          // (LM)

    // optional pre-filtering (reuses Luna FIR)
    bool   do_filter;
    double f_lwr;                // high-pass edge (Hz) if set
    double f_upr;                // low-pass edge (Hz) if set

    // verbose event-level output
    bool   verbose;
  };

  // --- event structures ---------------------------------------------------

  // one detected event on one physical input signal
  struct lm_component_t {
    int comp_id;                 // 1-based
    side_t side;                 // SIDE_L / SIDE_R / SIDE_LEG
    std::string sig;             // channel label

    uint64_t onset_tp, offset_tp;
    double onset_sec;            // event onset in seconds from EDF start
    double dur;                  // seconds (real elapsed)

    double base_onset, scale_onset; // B(t),S(t) at onset
    double on_thr, off_thr;      // representative thresholds at onset

    double peak;                 // max rectified amplitude
    double peak_exc;             // max (y - B)
    double peak_z;               // max (y - B)/S
    double morph_value;          // max 0.5s window-median excursion
    bool   morph_ok;

    bool   clm_component;        // 0.5 <= dur <= max_clm_dur

    std::string stage, state;    // stage at onset
    bool   qc_base_high;
  };

  // final evaluation-unit LM (monolateral component or bilateral grouping)
  struct lm_event_t {
    int evt_id;                  // 1-based
    side_t side;                 // L / R / B / LEG

    uint64_t onset_tp, offset_tp;
    double onset_sec;
    double dur;                  // seconds

    std::vector<int> comp_ids;   // constituent component ids
    int n_comp, n_l, n_r;
    bool all_comps_clm;          // every constituent component is CLM

    bool clm;
    bool plm;

    int seq;                     // sequence id (0 = none)
    int seq_pos;                 // 1-based position within its sequence
    int seq_n;                   // sequence length

    double imi_prev, imi_next;   // seconds; NaN if not measurable

    std::string stage, state;
    double peak_exc_max, peak_z_max;
    bool   qc_base_high;
  };

  struct plm_sequence_t {
    int seq_id;                  // 1-based
    uint64_t onset_tp, offset_tp;
    double dur;                  // seconds
    int n;                       // number of PLMs
    std::vector<double> imis;    // within-sequence IMIs (seconds)
    int n_l, n_r, n_b, n_leg;
  };

  // --- pure WASM 2016 grammar (unit-testable) -----------------------------

  // MODE=LR : connected-component bilateral combination of L/R components.
  // Components are re-sorted by onset on entry; returns events sorted by
  // onset with side/n_comp/all_comps_clm populated (clm not yet set).
  std::vector<lm_event_t> combine_bilateral( std::vector<lm_component_t> & comps ,
                                             const plm_param_t & P );

  // MODE=LEG : each component becomes one final event.
  std::vector<lm_event_t> combine_single( std::vector<lm_component_t> & comps ,
                                          const plm_param_t & P );

  // classify each final event as CLM / non-CLM (WASM rules).
  void classify_clm( std::vector<lm_event_t> & events , const plm_param_t & P );

  // detect PLM sequences over the chronological final-event stream; sets
  // plm / seq / seq_pos / seq_n / imi_prev / imi_next on events and returns
  // the sequences. Events must already be CLM-classified.
  std::vector<plm_sequence_t> detect_sequences( std::vector<lm_event_t> & events ,
                                                const plm_param_t & P );

}

// --- command driver -------------------------------------------------------

struct leg_movements_t {
  leg_movements_t( edf_t & edf , param_t & param );
};

#endif
