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

#include "stats/dpp-hypno.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "helper/helper.h"
#include "timeline/timeline.h"

#include <algorithm>
#include <limits>

namespace {

// nearest-sample lookup: no shared/exported version of this exists to
// call instead (SIGDYN's and stats/dpp.cpp's own are each private,
// file-local helpers too -- same idiom repeated here rather than adding
// a new shared-header dependency for a 12-line function)
int nearest_sample(const std::vector<uint64_t> &tp, uint64_t t) {
  if (tp.empty())
    return -1;
  std::vector<uint64_t>::const_iterator it =
      std::lower_bound(tp.begin(), tp.end(), t);
  if (it == tp.begin())
    return 0;
  if (it == tp.end())
    return (int)tp.size() - 1;
  const int hi = (int)(it - tp.begin());
  const int lo = hi - 1;
  const uint64_t dh = *it > t ? *it - t : t - *it;
  const uint64_t dl = t > tp[lo] ? t - tp[lo] : tp[lo] - t;
  return dh < dl ? hi : lo;
}

} // namespace

std::vector<std::string> dpp_hypno::stage_labels(bool three_state) {
  if (three_state)
    return {"W", "R", "NR"};
  return {"W", "R", "N1", "N2", "N3"};
}

dpp_hypno::lookup_t::lookup_t(edf_t &edf, const std::string &prefix,
                              bool three_state) {
  labels = stage_labels(three_state);
  traces.resize(labels.size());

  for (size_t i = 0; i < labels.size(); i++) {
    const std::string ch = prefix + "_" + labels[i];

    if (!edf.header.has_signal(ch))
      Helper::halt("hypnodensity channel " + ch +
                   " not present in the attached EDF (run POPS ... "
                   "posterior-channels= first)");

    signal_list_t sl = edf.header.signal_list(ch);
    if (sl.size() != 1)
      Helper::halt("could not resolve hypnodensity channel " + ch);

    slice_t wtrace(edf, sl(0), edf.timeline.wholetrace());
    traces[i].raw = *wtrace.pdata();
    traces[i].tp = *wtrace.ptimepoints();
  }
}

std::vector<double> dpp_hypno::lookup_t::at(uint64_t t_tp) const {
  std::vector<double> vals(labels.size(),
                           std::numeric_limits<double>::quiet_NaN());
  for (size_t i = 0; i < labels.size(); i++) {
    const int idx = nearest_sample(traces[i].tp, t_tp);
    if (idx >= 0)
      vals[i] = traces[i].raw[idx];
  }
  return vals;
}
