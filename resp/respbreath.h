
//    --------------------------------------------------------------------
//
//    This file is part of Luna.
//
//    LUNA is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    LUNA is distributed in the hope that it will be useful,
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
// RESPBREATH : Respiratory breath segmentation
//
// Detects individual breaths from respiratory-like PSG channels (nasal
// cannula, thermistor, effort belt) and emits breath-level timing
// annotations.  The primary goal is robust timing for downstream
// phase-locking / coupling analyses with EEG.
//

#ifndef __RESPBREATH_H__
#define __RESPBREATH_H__

#include <vector>
#include <string>
#include <cstdint>

#include "intervals/intervals.h"

struct edf_t;
struct param_t;


// ---------------------------------------------------------------------------
// One detected breath event
// ---------------------------------------------------------------------------

struct breath_event_t {

  // Timing (seconds from recording start)
  double start_sec;         // breath onset (preceding trough)
  double peak_sec;          // inspiratory peak (main fiducial)
  double end_sec;           // breath end (following trough)

  // Derived durations (seconds)
  double t_insp;            // peak - start
  double t_exp;             // end  - peak
  double t_tot;             // end  - start

  // Amplitude
  double amp_insp;          // peak - preceding trough
  double amp_exp;           // peak - following trough
  double amp_sym;           // peak - mean(adjacent troughs)

  // Quality / metadata
  double confidence;        // [0, 1]
  bool   low_conf;          // convenience flag
  bool   near_artifact;     // adjacent to artifact interval
  bool   fused;             // came from multi-channel consensus
  int    source_ch;         // index into channels vector (0-based)
  int    n_supporting_ch;   // how many channels agreed on this breath

  // Absolute timepoints (for annotation)
  uint64_t tp_start;
  uint64_t tp_peak;
  uint64_t tp_end;

  // Segment / breath index bookkeeping
  int seg_idx;
  int breath_idx;           // within the full recording

  breath_event_t()
    : start_sec(0), peak_sec(0), end_sec(0),
      t_insp(0), t_exp(0), t_tot(0),
      amp_insp(0), amp_exp(0), amp_sym(0),
      confidence(0), low_conf(true), near_artifact(false),
      fused(false), source_ch(0), n_supporting_ch(1),
      tp_start(0), tp_peak(0), tp_end(0),
      seg_idx(0), breath_idx(0)
  { }
};


// ---------------------------------------------------------------------------
// One artifact / unusable interval
// ---------------------------------------------------------------------------

struct breath_artifact_t {
  double   start_sec;
  double   end_sec;
  uint64_t tp_start;
  uint64_t tp_end;
  int      source_ch;   // -1 = all-channel / fused artifact
  std::string reason;   // short diagnostic string

  breath_artifact_t()
    : start_sec(0), end_sec(0), tp_start(0), tp_end(0),
      source_ch(-1), reason("unknown")
  { }
};


// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

struct respbreath_t {
  respbreath_t( edf_t & edf , param_t & param );
};

#endif // __RESPBREATH_H__
