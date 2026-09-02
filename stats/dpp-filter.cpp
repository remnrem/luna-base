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

#include "stats/dpp-filter.h"

#include "dsp/fir.h"
#include "timeline/timeline.h"

#include <algorithm>

double dpp_filters::pad_seconds(const dpp_filter_t &filt, int sr) {
  fir_t fir;
  int windowLength = 0;
  double beta = 0;
  fir.calculateKaiserParams(filt.ripple, filt.tw, (double)sr, &windowLength,
                            &beta);
  return 1.5 * (windowLength / 2.0) / (double)sr;
}

namespace {

// genuinely causal bandpass: design the same Kaiser-window coefficients
// dsptools::apply_fir() would (reused directly, not reimplemented), but
// apply them via fir_impl_t::getOutputSample() in plain sample order --
// no burn-in/output-relabeling trick, so output[i] depends only on
// input[0..i]. delayLine starts zero-initialized (fir_impl_t's own
// constructor), giving a settling transient over the first ~(taps-1)/2
// samples rather than ever reading ahead. See stats/dpp-filter.h's top
// comment for why dsptools::apply_fir() itself isn't used here.
std::vector<double> causal_bandpass(const std::vector<double> &x, int sr,
                                    double ripple, double tw, double lwr,
                                    double upr) {
  std::vector<double> fc =
      dsptools::design_bandpass_fir(ripple, tw, (double)sr, lwr, upr);
  fir_impl_t fir(fc);
  std::vector<double> out(x.size());
  for (int i = 0; i < (int)x.size(); i++)
    out[i] = fir.getOutputSample(x[i]);
  return out;
}

} // namespace

std::vector<double>
dpp_filters::prefilter_trace(const std::vector<double> &x,
                             const std::vector<uint64_t> &tp, int sr,
                             double lwr, double upr) {
  const int n = (int)x.size();

  // default: raw pass-through, only overwritten below for segments long
  // enough to filter safely (also correctly handles n==0)
  std::vector<double> out(x);

  // minimum segment length worth attempting to filter, from the same
  // Kaiser design this filter itself uses -- shorter contiguous runs are
  // left unfiltered rather than risking a degenerate/out-of-bounds FIR
  // burn-in on a near-empty segment
  fir_t fir;
  int windowLength = 0;
  double beta = 0;
  fir.calculateKaiserParams(0.02, 1.0, (double)sr, &windowLength, &beta);
  const int min_len = windowLength + 1;

  // filter each contiguous (gap-free) run independently, so filter state
  // never carries across a genuine discontinuity (e.g. from MASK+RE) --
  // otherwise samples just after a gap would inherit filter history from
  // before it, contaminating any window entirely on the post-gap side
  // even though that window itself passes get_window()'s own
  // discontinuity check
  int seg_start = 0;
  for (int i = 1; i <= n; i++) {
    const bool boundary =
        (i == n) || timeline_t::discontinuity(tp, sr, i - 1, i);
    if (!boundary)
      continue;

    const int seg_len = i - seg_start;
    if (seg_len >= min_len) {
      std::vector<double> seg(x.begin() + seg_start, x.begin() + i);
      std::vector<double> filtered =
          causal_bandpass(seg, sr, 0.02, 1.0, lwr, upr);
      std::copy(filtered.begin(), filtered.end(), out.begin() + seg_start);
    }

    seg_start = i;
  }

  return out;
}

std::vector<double> dpp_filters::apply_band(const std::vector<double> &x,
                                            int sr, const dpp_filter_t &filt) {
  return causal_bandpass(x, sr, filt.ripple, filt.tw, filt.lwr, filt.upr);
}
