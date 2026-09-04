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

#ifndef __LUNA_DPP_H__
#define __LUNA_DPP_H__

#include <vector>

struct edf_t;
struct param_t;
struct dpp_matrix_t;

namespace dsptools {
// DPP stage 2: generic multiscale feature-extraction over trailing,
// causal, fixed-length windows -- see stats/dpp-spec.h / dpp-filter.h /
// dpp-io.h for the spec grammar, filtering, and binary-corpus I/O this
// orchestrates. No train=/model= modes yet (stage 3).
void dpp(edf_t &edf, param_t &param);
} // namespace dsptools

namespace dpp_classic {
bool extract(edf_t &, param_t &, const std::vector<double> &, dpp_matrix_t *);
} // namespace dpp_classic

#endif
