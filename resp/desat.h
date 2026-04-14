
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

#ifndef __DESAT_H__
#define __DESAT_H__

#include <vector>
#include <string>
#include <deque>

#include "intervals/intervals.h"

struct edf_t;
struct param_t;


// One detected desaturation event
struct desat_event_t {
  double start_sec;   // seconds from recording start
  double stop_sec;
  double duration;    // seconds
  double nadir;       // minimum SpO2 (%) during event
  double baseline;    // local baseline at event onset
  double drop;        // baseline - nadir (magnitude)
  int    seg;         // segment index (0-based)
};


// Per-segment summary
struct desat_seg_results_t {
  int    n_desats;
  double t_valid;     // seconds of valid (non-artifact) signal
  double t_artifact;  // seconds of artifact
  double t_desat;     // seconds spent in detected desats
  double mean_spo2;   // mean SpO2 of valid samples
  double pct_lt90;    // % time (valid) below 90
  double pct_lt88;
  double pct_lt85;
  double pct_lt80;
  std::vector<desat_event_t> events;
};


// Entry point
struct desat_t {
  desat_t( edf_t & edf , param_t & param );
};

#endif
