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
//    Common ONNX Runtime helpers.
//    --------------------------------------------------------------------

#ifndef LUNA_ORT_COMMON_H
#define LUNA_ORT_COMMON_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#ifdef HAS_ORT
#include <onnxruntime/core/session/onnxruntime_cxx_api.h>

namespace ort_common
{

Ort::AllocatedStringPtr named_input(Ort::Session &, Ort::AllocatorWithDefaultOptions &,
                                    const std::string &);
Ort::AllocatedStringPtr named_output(Ort::Session &, Ort::AllocatorWithDefaultOptions &,
                                     const std::string &);
size_t input_index(Ort::Session &, const std::string &);
size_t output_index(Ort::Session &, const std::string &);
std::vector<int64_t> tensor_shape(const Ort::TypeInfo &);
void check_shape(const std::vector<int64_t> &, const std::vector<int64_t> &, const std::string &);

} // namespace ort_common
#endif

#endif
