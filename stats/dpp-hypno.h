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

#ifndef __LUNA_DPP_HYPNO_H__
#define __LUNA_DPP_HYPNO_H__

#include <cstdint>
#include <string>
#include <vector>

struct edf_t;

// Reads already-attached POPS hypnodensity signals (PP_W, PP_R, PP_N1,
// PP_N2, PP_N3, or PP_W/PP_R/PP_NR for three-state -- exact strings from
// posterior_labels(), pops/posteriors.cpp) as ordinary EDF signals. No
// hypnogram_t/staging dependency: this is a generic signal point-lookup,
// the same in spirit as how stats/dpp.cpp reads any other channel, just
// for a single instantaneous nearest-timepoint value per call rather than
// a windowed feature. Deliberately does not reuse stats/dpp.cpp's
// dpp_trace_t/tr.sr: hypnodensity channels are sub-1Hz (e.g. one value per
// 30s epoch), and tr.sr's (int) cast of a fractional Fs truncates to 0,
// corrupting the spectral/window arithmetic that assumes sr>=1 -- not
// needed here anyway, since no windowing happens, just a nearest-sample
// lookup (see the implementation plan's Stage 4 design).
//
// Precondition: a prior `POPS ... posterior-channels=` run already
// attached these channels. Not triggered or computed by DPP itself.

namespace dpp_hypno {

// canonical stage label lists, matching pops/posteriors.cpp's
// posterior_labels() exactly
std::vector<std::string> stage_labels(bool three_state);

struct lookup_t {
  // halts if any required PP_<label> channel is missing from the EDF
  lookup_t(edf_t &edf, const std::string &prefix, bool three_state);

  // nearest-timepoint PP_s value for each stage, in labels order
  std::vector<double> at(uint64_t t_tp) const;

  std::vector<std::string> labels;

private:
  struct trace_t {
    std::vector<double> raw;
    std::vector<uint64_t> tp;
  };
  std::vector<trace_t> traces; // parallel to labels
};

} // namespace dpp_hypno

#endif
