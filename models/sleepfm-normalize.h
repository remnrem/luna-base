
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

// Shared SleepFM input normalization for ORT and POPS.
#ifndef LUNA_SLEEPFM_NORMALIZE_H
#define LUNA_SLEEPFM_NORMALIZE_H

#include "edf/edf.h"
#include "edf/slice.h"
#include "dsp/resample.h"
#include "helper/helper.h"

#include <cmath>
#include <vector>

namespace sleepfm_preprocessing
{

struct normalization_t
{
  double mean;
  double sd;

  explicit normalization_t(const std::vector<double> &samples)
  {
    if (samples.empty()) Helper::halt("SleepFM input contains no samples");
    const double origin = samples.front();
    double sum = 0;
    for (double v : samples)
      {
        if (!std::isfinite(v)) Helper::halt("SleepFM input contains non-finite samples");
        sum += v - origin;
      }
    mean = origin + sum / samples.size();
    // Two passes avoid cancellation in E[x^2] - E[x]^2 for near-constant channels.
    double sum2 = 0;
    for (double v : samples)
      {
        const double d = v - mean;
        sum2 += d * d;
      }
    sd = std::sqrt(sum2 / samples.size()); // population SD, as in numpy.std
    if (!std::isfinite(mean) || !std::isfinite(sd))
      Helper::halt("SleepFM normalization statistics are non-finite");
  }

  void apply(std::vector<double> *samples) const
  {
    for (double &v : *samples)
      {
        const double centered = v - mean;
        v = sd > 0 ? centered / sd : centered;
        if (!std::isfinite(v)) Helper::halt("SleepFM preprocessing produced non-finite data");
      }
  }
};

inline std::vector<normalization_t> recording_normalization(
    edf_t &edf, const signal_list_t &signals, int sample_rate)
{
  // Use all available recording samples, including the tail outside complete
  // model windows. Never derive this interval from epochs or the stride offset.
  // Recomputing on a later POPS pass therefore uses the same population; keeping
  // these statistics local avoids stale caches after EDF/channel modifications.
  const interval_t whole(0, edf.timeline.last_time_point_tp + 1);
  std::vector<normalization_t> result;
  result.reserve(signals.size());
  for (int c = 0; c < signals.size(); ++c)
    {
      slice_t sl(edf, signals(c), whole);
      const std::vector<double> &raw = *sl.pdata();
      if (raw.empty()) Helper::halt("SleepFM input contains no samples");
      for (double v : raw)
        if (!std::isfinite(v)) Helper::halt("SleepFM input contains non-finite samples");
      // Match upstream's normalization population: the whole resampled channel.
      // Retain Luna's existing sinc resampler and window extraction behavior.
      result.emplace_back(dsptools::resample(&raw,
          edf.header.sampling_freq(signals(c)), sample_rate, SRC_SINC_BEST_QUALITY));
    }
  return result;
}

} // namespace sleepfm_preprocessing
#endif
