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
//
//    Two-level vector DPP: cross-fitted local model plus subject summaries.
//    --------------------------------------------------------------------

#ifndef __LUNA_DPP_TWOLEVEL_H__
#define __LUNA_DPP_TWOLEVEL_H__

#include <string>

struct edf_t;
struct param_t;
struct dpp_matrix_t;

namespace dpp_twolevel {
bool enabled(const param_t &);
void fit(param_t &);
void fit_model_set(param_t &);
void apply(edf_t &, param_t &, const dpp_matrix_t &);
} // namespace dpp_twolevel

#endif
