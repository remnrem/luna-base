#ifndef LUNA_ORT_COMMON_H
#define LUNA_ORT_COMMON_H

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>

#ifdef HAS_ORT
#include <onnxruntime/core/session/onnxruntime_cxx_api.h>

namespace ort_common {

Ort::AllocatedStringPtr named_input(Ort::Session &, Ort::AllocatorWithDefaultOptions &, const std::string &);
Ort::AllocatedStringPtr named_output(Ort::Session &, Ort::AllocatorWithDefaultOptions &, const std::string &);
size_t input_index(Ort::Session &, const std::string &);
size_t output_index(Ort::Session &, const std::string &);
std::vector<int64_t> tensor_shape(const Ort::TypeInfo &);
void check_shape(const std::vector<int64_t> &, const std::vector<int64_t> &, const std::string &);

}
#endif

#endif
