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

#ifndef __LUNA_DPP_FILTER_H__
#define __LUNA_DPP_FILTER_H__

#include <vector>
#include <cstdint>

#include "stats/dpp-spec.h"

// Thin glue over the existing FIR kernel (dsptools::apply_fir(), dsp/fir.h)
// and its own Kaiser settling-length calculator (fir_t::calculateKaiserParams,
// called directly, not reimplemented). No new filtering algorithm here --
// just (a) a helper to size a causal run-in pad for a named band, (b) a
// one-shot whole-trace prefilter (no padding needed -- see stats/dpp.cpp),
// and (c) applying a named band to an already-extracted (padded) window
// buffer. Window extraction, padding, and gap-checking all live in
// stats/dpp.cpp, which owns the trailing-buffer/window loop; this file only
// does the actual filtering.

namespace dpp_filters {

  // causal run-in length (seconds) for this band, derived directly from
  // fir_t::calculateKaiserParams()'s own filter-length calculation for the
  // requested ripple/transition-width -- proportional to the real filter
  // length, not a fixed guess
  double pad_seconds( const dpp_filter_t & filt , int sr );

  // one-time bandpass (default ripple/tw), no padding within a contiguous
  // run -- a filter transient confined to the first/last fraction of a
  // second of a multi-hour recording is negligible. 'tp' is required and
  // used to detect internal discontinuities (e.g. from MASK+RE): each
  // contiguous run is filtered independently, so filter state never
  // carries across a genuine gap in the recording (a window entirely on
  // one side of a gap would otherwise still pick up a filtered value
  // contaminated by history from the other side, even though the window
  // itself passes get_window()'s own discontinuity check)
  std::vector<double> prefilter_trace( const std::vector<double> & x ,
				       const std::vector<uint64_t> & tp ,
				       int sr , double lwr , double upr );

  // apply a named band to an already-extracted buffer (the caller is
  // responsible for including whatever causal run-in padding it wants
  // filtered along with the reporting window, and for dropping/using that
  // padding afterward as appropriate for the feature being computed)
  std::vector<double> apply_band( const std::vector<double> & x , int sr ,
				  const dpp_filter_t & filt );

}

#endif
