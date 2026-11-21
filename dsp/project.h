// Fixed numerical projection, independent resource/kernel API.
// Part of Luna; GPL-3.0-or-later (see LICENSE.txt).
#ifndef LUNA_PROJECT_H
#define LUNA_PROJECT_H
#include "stats/Eigen/Dense"
#include "helper/json.h"
#include <map>
#include <string>
#include <vector>
struct edf_t;
struct param_t;
namespace project {
struct paths_t {
  std::string matrix, channels, parcels, metadata, manifest, bundle;
};
struct operator_t {
  Eigen::MatrixXd weights;
  std::vector<std::string> inputs, outputs, parcel_ids;
  nlohmann::json metadata;
  std::map<std::string,std::string> hashes, files;
  std::string name, input_unit, output_unit, reference;
  bool average_reference = false, embedded_projector = false;
};
std::string sha256(const std::string & bytes);
std::string normalize(const std::string & label);
double voltage_scale(const std::string & unit);
operator_t load(paths_t paths);
Eigen::MatrixXd multiply(const Eigen::MatrixXd & weights, const Eigen::MatrixXd & inputs);
void run(edf_t & edf, param_t & param);
}
#endif
