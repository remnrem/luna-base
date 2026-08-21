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

#include "models/ort-common.h"

#ifdef HAS_ORT

#include "helper/helper.h"
#include <stdexcept>

namespace ort_common {

Ort::AllocatedStringPtr named_input(Ort::Session &s, Ort::AllocatorWithDefaultOptions &a, const std::string &want) {
  for (size_t i = 0; i < s.GetInputCount(); ++i) {
    auto n = s.GetInputNameAllocated(i, a);
    if (n && want == n.get()) return n;
  }
  throw std::runtime_error("ORT input not found: " + want);
}

Ort::AllocatedStringPtr named_output(Ort::Session &s, Ort::AllocatorWithDefaultOptions &a, const std::string &want) {
  for (size_t i = 0; i < s.GetOutputCount(); ++i) {
    auto n = s.GetOutputNameAllocated(i, a);
    if (n && want == n.get()) return n;
  }
  throw std::runtime_error("ORT output not found: " + want);
}

size_t input_index(Ort::Session &s, const std::string &want) {
  Ort::AllocatorWithDefaultOptions a;
  for (size_t i = 0; i < s.GetInputCount(); ++i) {
    auto n = s.GetInputNameAllocated(i, a);
    if (n && want == n.get()) return i;
  }
  throw std::runtime_error("ORT input not found: " + want);
}

size_t output_index(Ort::Session &s, const std::string &want) {
  Ort::AllocatorWithDefaultOptions a;
  for (size_t i = 0; i < s.GetOutputCount(); ++i) {
    auto n = s.GetOutputNameAllocated(i, a);
    if (n && want == n.get()) return i;
  }
  throw std::runtime_error("ORT output not found: " + want);
}

std::vector<int64_t> tensor_shape(const Ort::TypeInfo &t) {
  return t.GetTensorTypeAndShapeInfo().GetShape();
}

void check_shape(const std::vector<int64_t> &actual,
                 const std::vector<int64_t> &expected,
                 const std::string &what) {
  if (actual.size() != expected.size()) Helper::halt("ORT " + what + " rank does not match expected shape");
  for (size_t i = 0; i < actual.size(); ++i)
    if (expected[i] >= 0 && actual[i] >= 0 && actual[i] != expected[i])
      Helper::halt("ORT " + what + " shape does not match expected shape");
}

}

#endif
