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

// SleepFM-specific ONNX Runtime execution for signal windows.

#include "edf/edf.h"
#include "edf/slice.h"
#include "param.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "db/db.h"
#include "models/ort-common.h"

#ifdef HAS_ORT
#include <onnxruntime/core/session/onnxruntime_cxx_api.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

extern logger_t logger;
extern writer_t writer;

namespace {

static std::string basename(const std::string &p) {
  const size_t n = p.find_last_of("/\\");
  return n == std::string::npos ? p : p.substr(n + 1);
}

static std::string normalize_model_path(const std::string &path) {
  const std::string p = Helper::expand(path);
  const std::string suffix = ".onnx";
  if (p.size() >= suffix.size() && p.compare(p.size() - suffix.size(), suffix.size(), suffix) == 0) return p;
  return p + suffix;
}

static std::string join_path(const std::string &path, const std::string &root) {
  if (path.empty() || path == ".") return root;
  if (path.back() == '/' || path.back() == '\\') return path + root;
  return path + "/" + root;
}

static std::vector<float> resample(const std::vector<double> &x, double from, double to, int n) {
  std::vector<float> y(n, 0);
  if (x.empty()) return y;
  for (int i = 0; i < n; ++i) {
    const double at = i * from / to;
    const int j = std::min<int>(static_cast<int>(at), x.size() - 1);
    const int k = std::min<int>(j + 1, x.size() - 1);
    const double a = x[j], b = x[k], f = at - j;
    y[i] = std::isfinite(a) && std::isfinite(b)
      ? static_cast<float>(a + f * (b - a))
      : std::numeric_limits<float>::quiet_NaN();
  }
  return y;
}

static void standardize(std::vector<float> *x) {
  double sum = 0, sum2 = 0;
  int n = 0;
  for (float v : *x) {
    if (!std::isfinite(v)) Helper::halt("SleepFM input contains non-finite samples");
    sum += v;
    sum2 += v * v;
    ++n;
  }
  if (!n) Helper::halt("SleepFM input contains no samples");
  const double mean = sum / n;
  const double var = std::max(0.0, sum2 / n - mean * mean);
  const double sd = std::sqrt(var);
  for (float &v : *x) {
    v -= static_cast<float>(mean);
    if (sd > 0) v /= static_cast<float>(sd);
    if (!std::isfinite(v)) Helper::halt("SleepFM preprocessing produced non-finite data");
  }
}

}

void proc_ort(edf_t &edf, param_t &param) {
#ifndef HAS_ORT
  (void)edf;
  (void)param;
  Helper::halt("ORT support was not compiled in; rebuild with ORT=1");
#else
  constexpr int sample_rate = 128;
  constexpr int window_seconds = 300;
  constexpr int window_samples = sample_rate * window_seconds;
  constexpr int embedding_dim = 128;
  constexpr int sequence_length = 60;

  const std::string model_root = join_path(
    param.has("path") ? Helper::expand(param.value("path")) : ".",
    param.has("lib") ? Helper::expand(param.value("lib")) : "sleepfm");
  const std::string model = normalize_model_path(model_root);
  std::ifstream model_file(model);
  if (!model_file) Helper::halt("SleepFM ONNX model missing: " + model);

  std::string modality = param.has("modality") ? param.value("modality") : "";
  if (modality == "auto" || modality.empty())
    Helper::halt("SleepFM ORT requires modality=BAS, RESP, EKG, or EMG; channel aliases are resolved by Luna");
  const int max_channels = modality == "BAS" ? 10 : modality == "RESP" ? 7 : modality == "EKG" ? 2 : modality == "EMG" ? 4 : 0;
  if (!max_channels) Helper::halt("unknown SleepFM modality: " + modality);

  signal_list_t signals = edf.header.signal_list(param.requires("sig"), true);
  if (!signals.size()) Helper::halt("SleepFM ORT found no input signals");
  if (signals.size() > max_channels) Helper::halt("SleepFM ORT excessive channel count for modality " + modality);
  if (edf.timeline.last_time_point_tp < window_seconds * globals::tp_1sec)
    Helper::halt("SleepFM ORT requires at least one complete five-minute window");

  const uint64_t window_ticks = window_seconds * globals::tp_1sec;
  const int requested_step = param.has("step") ? param.requires_int("step") : window_seconds;
  if (requested_step <= 0) Helper::halt("SleepFM ORT step must be positive");
  const uint64_t step_ticks = requested_step * globals::tp_1sec;
  const int n_windows = static_cast<int>((edf.timeline.last_time_point_tp - window_ticks) / step_ticks) + 1;
  logger << "  ORT SleepFM: model=" << model << ", modality=" << modality
         << ", channels=" << signals.size() << "/" << max_channels << "\n"
         << "    preprocessing: target rate=" << sample_rate << " Hz; ";
  bool needs_resampling = false;
  for (int c = 0; c < signals.size(); ++c) {
    if (std::fabs(edf.header.sampling_freq(signals(c)) - sample_rate) > 1e-9)
      needs_resampling = true;
  }
  if (needs_resampling) {
    logger << "resampling input signal(s) to target rate; ";
  } else {
    logger << "input signal rate already matches target (no resampling needed); ";
  }
  logger << "standardized within each window\n"
         << "    windowing: " << n_windows << " complete " << window_seconds
         << "s windows, step=" << requested_step << "s\n"
         << "    input: " << window_samples << " samples/channel; unused slots zero-padded\n";

  Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "luna-sleepfm");
  Ort::SessionOptions opts;
  opts.SetIntraOpNumThreads(param.has("threads") ? param.requires_int("threads") : 1);
  opts.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_BASIC);
  Ort::Session session(env, model.c_str(), opts);
  Ort::AllocatorWithDefaultOptions alloc;

  if (session.GetInputCount() < 2 || session.GetOutputCount() < 2)
    Helper::halt("SleepFM ONNX model must have signal, channel_mask, pooled_embedding, and sequence_embedding tensors");
  auto inname = ort_common::named_input(session, alloc, "signal");
  auto maskname = ort_common::named_input(session, alloc, "channel_mask");
  auto pooled = ort_common::named_output(session, alloc, "pooled_embedding");
  auto sequence = ort_common::named_output(session, alloc, "sequence_embedding");
  const size_t signal_idx = ort_common::input_index(session, "signal");
  const size_t mask_idx = ort_common::input_index(session, "channel_mask");
  const size_t pooled_idx = ort_common::output_index(session, "pooled_embedding");
  const size_t sequence_idx = ort_common::output_index(session, "sequence_embedding");
  ort_common::check_shape(ort_common::tensor_shape(session.GetInputTypeInfo(signal_idx)), {1, -1, window_samples}, "SleepFM signal input");
  ort_common::check_shape(ort_common::tensor_shape(session.GetInputTypeInfo(mask_idx)), {1, -1}, "SleepFM channel-mask input");
  ort_common::check_shape(ort_common::tensor_shape(session.GetOutputTypeInfo(pooled_idx)), {1, embedding_dim}, "SleepFM pooled output");
  ort_common::check_shape(ort_common::tensor_shape(session.GetOutputTypeInfo(sequence_idx)), {1, sequence_length, embedding_dim}, "SleepFM sequence output");
  if (session.GetInputTypeInfo(signal_idx).GetTensorTypeAndShapeInfo().GetElementType() != ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT ||
      session.GetInputTypeInfo(mask_idx).GetTensorTypeAndShapeInfo().GetElementType() != ONNX_TENSOR_ELEMENT_DATA_TYPE_BOOL)
    Helper::halt("SleepFM ONNX input dtypes must be float32 signal and bool channel_mask");

  const std::string requested = param.has("output") ? param.value("output") : "pooled_embedding";
  if (requested != "pooled_embedding" && requested != "sequence_embedding")
    Helper::halt("SleepFM output must be pooled_embedding or sequence_embedding");
  const char *selected = requested == "pooled_embedding" ? pooled.get() : sequence.get();
  const bool add_channels = param.has("add-channels");
  const bool no_output = param.has("no-output");
  if (add_channels && requested != "sequence_embedding")
    Helper::halt("SleepFM add-channels is available only for sequence_embedding; pooled embeddings are not emitted as EDF channels");
  const std::string channel_root = add_channels && !param.value("add-channels").empty()
    ? param.value("add-channels") : modality;
  const int output_interval = requested == "sequence_embedding" ? 5 : window_seconds;
  logger << "    output: " << requested;
  if (requested == "sequence_embedding") logger << " (60 contextual 5s tokens/window, 128 values/token)";
  else logger << " (128 values/window)";
  if (no_output) logger << ", database output suppressed";
  logger << "\n";
  int output_samples_per_record = 0;
  std::vector<std::vector<double> > channel_data;
  if (add_channels) {
    if (edf.is_actually_discontinuous())
      Helper::halt("SleepFM add-channels requires a continuous EDF");
    if (requested_step != window_seconds)
      Helper::halt("SleepFM add-channels requires step=300 so outputs map to regular EDF channels");
    const double n_spr_d = edf.header.record_duration / output_interval;
    output_samples_per_record = static_cast<int>(std::lround(n_spr_d));
    if (output_samples_per_record < 1 || std::fabs(n_spr_d - output_samples_per_record) > 1e-6)
      Helper::halt("SleepFM add-channels requires an EDF record size compatible with " +
                   Helper::int2str(output_interval) + "-second outputs; use RECORD-SIZE dur=" +
                   Helper::int2str(output_interval) + " continue");
    const size_t n_output = static_cast<size_t>(edf.header.nr) * output_samples_per_record;
    channel_data.assign(embedding_dim, std::vector<double>(n_output, -9.0));
    logger << "  SleepFM: will add " << embedding_dim << " EDF channel(s) with root '"
           << channel_root << "' at " << output_interval << "-second resolution\n";
  }

  const uint64_t win = window_seconds * globals::tp_1sec;
  const uint64_t step = requested_step * globals::tp_1sec;
  Ort::MemoryInfo mem = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
  int w = 0;
  for (uint64_t start = 0; start + win <= edf.timeline.last_time_point_tp; start += step, ++w) {
    interval_t iv(start, start + win);
    std::vector<float> x(static_cast<size_t>(max_channels) * window_samples, 0);
    std::unique_ptr<bool[]> mask(new bool[max_channels]);
    for (int c = 0; c < max_channels; ++c) mask[c] = true;
    for (int c = 0; c < signals.size(); ++c) {
      slice_t sl(edf, signals(c), iv);
      std::vector<float> z = resample(*sl.pdata(), edf.header.sampling_freq(signals(c)), sample_rate, window_samples);
      standardize(&z);
      std::copy(z.begin(), z.end(), x.begin() + static_cast<size_t>(c) * window_samples);
      mask[c] = false;
    }
    std::array<int64_t, 3> xs{1, max_channels, window_samples};
    std::array<int64_t, 2> ms{1, max_channels};
    auto xt = Ort::Value::CreateTensor<float>(mem, x.data(), x.size(), xs.data(), 3);
    auto mt = Ort::Value::CreateTensor<bool>(mem, mask.get(), max_channels, ms.data(), 2);
    const char *ins[] = {inname.get(), maskname.get()};
    const char *outs[] = {selected};
    Ort::Value ivals[] = {std::move(xt), std::move(mt)};
    auto out = session.Run(Ort::RunOptions{nullptr}, ins, ivals, 2, outs, 1);
    auto ti = out[0].GetTensorTypeAndShapeInfo();
    if (ti.GetElementType() != ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT || ti.GetElementCount() == 0)
      Helper::halt("SleepFM ONNX output is not a nonempty float tensor");
    const float *v = out[0].GetTensorData<float>();
    const size_t count = ti.GetElementCount();
    std::string labels;
    for (int c = 0; c < signals.size(); ++c) { if (c) labels += ","; labels += signals.label(c); }
    if (!no_output) {
      writer.level(w + 1, "WIN");
      writer.value("START", iv.start_sec()); writer.value("END", iv.stop_sec());
      writer.value("MODEL", basename(model)); writer.value("MODALITY", modality);
      writer.value("CHANNELS", signals.size()); writer.value("CHANNEL_LIST", labels);
      writer.value("SR", sample_rate); writer.value("OUTPUT", requested); writer.value("N", static_cast<int>(count));
    }
    const bool seq = requested == "sequence_embedding";
    for (size_t i = 0; i < count; ++i) {
      if (!no_output) {
        writer.level(static_cast<int>(i + 1), "OUT");
        writer.value("TIME", iv.start_sec() + (seq ? (i / embedding_dim) * 5 : 0));
        writer.value("VALUE", v[i]);
      }
      if (add_channels) {
        const int token = seq ? static_cast<int>(i / embedding_dim) : 0;
        const double output_time = iv.start_sec() + token * output_interval;
        const long long sample = std::llround(output_time / output_interval);
        const size_t offset = static_cast<size_t>(sample);
        if (i % embedding_dim < channel_data.size() && offset < channel_data[i % embedding_dim].size())
          channel_data[i % embedding_dim][offset] = v[i];
      }
    }
    if (!no_output) writer.unlevel("OUT");
  }
  if (!no_output) writer.unlevel("WIN");
  if (add_channels) {
    int added = 0;
    for (int i = 0; i < embedding_dim; ++i) {
      const std::string label = channel_root + "_" + (i + 1 < 10 ? "00" : i + 1 < 100 ? "0" : "") + Helper::int2str(i + 1);
      if (edf.header.has_signal(label)) {
        logger << "  SleepFM: channel " << label << " already exists, skipping\n";
        continue;
      }
      edf.add_signal(label, -output_samples_per_record, channel_data[i]);
      ++added;
    }
    logger << "  SleepFM: added " << added << " channel(s) to EDF with root '" << channel_root << "'\n";
  }
  logger << "  ORT: evaluated " << w << " SleepFM window(s), model=" << basename(model) << ", modality=" << modality
         << (no_output ? ", database output suppressed" : "") << "\n";
#endif
}
