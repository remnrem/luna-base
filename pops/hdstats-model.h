// Read-only representation of a structured HDSTATS-PC V3 model.
//
// This deliberately has no dependency on the compiled HDSTATS registry or on
// Eigen.  It is therefore the single model-file reader suitable for clients
// such as Lunascope as well as Luna's own PCA projection code.
#ifndef __HDSTATS_MODEL_H__
#define __HDSTATS_MODEL_H__

#include <string>
#include <vector>

namespace hdstats_model {

struct feature_t {
  std::string id;
  double mean = 0.0;
  double sd = 0.0;
  std::vector<double> loadings;
};

struct panel_t {
  std::string id;
  int n_input = 0;
  int n_fit = 0;
  double total_ss = 0.0;
  std::vector<std::string> declared_features;
  std::vector<feature_t> features;
  std::vector<double> singular_values;
};

// V1 JSON models are self-describing: unlike V3, their component registry is
// carried by the model instead of being supplied by Luna's compiled manifest.
struct component_t {
  std::string id;
  std::string panel_id;
  int panel_pc = 0;                 // one-based column within panel loadings
  int orientation = 1;
  std::string label;
  std::string description;
};

// Optional post-fit annotation records. entity_type is PANEL, PC, or FEATURE.
// Correlation endpoint order is retained by the reader.
struct description_t {
  std::string entity_type;
  std::string entity_id;
  std::string text;
};

struct external_association_t {
  std::string entity_type;
  std::string entity_id;
  std::string factor;
  double score = 0.0;
};

struct correlation_t {
  std::string first;
  std::string second;
  double value = 0.0;
};

struct model_t {
  std::string model_id;
  std::string fingerprint;
  std::string specification_json;   // canonical embedded V1 spec; empty for V3
  std::string component_id_style = "legacy"; // legacy PANEL/1 or R-safe PANEL_1
  std::vector<panel_t> panels;
  std::vector<component_t> components;
  std::vector<description_t> descriptions;
  std::vector<external_association_t> external_associations;
  std::vector<correlation_t> pc_correlations;
  std::vector<correlation_t> feature_correlations;

  const panel_t * find_panel( const std::string & id ) const;
  const feature_t * find_feature( const std::string & panel_id,
                                  const std::string & feature_id ) const;
  const description_t * find_description( const std::string & entity_type,
                                          const std::string & entity_id ) const;
};

// Reads a legacy __HDSTATS_PC_V3 model (with its optional trailer) or a
// self-contained JSON HDSTATS-PC V1 model.
//
//   __HDSTATS_PC_ANNOTATIONS_V1
//   DESC       <DOMAIN|PC|FEATURE> <id> <text>
//   EXT_ASSOC  <DOMAIN|PC|FEATURE> <id> <factor> <signed-score>
//   PC_CORR    <pc-id-1> <pc-id-2> <correlation>
//   FTR_CORR   <feature-id-1> <feature-id-2> <correlation>
//   END_HDSTATS_PC_ANNOTATIONS
//
// The trailer is optional.  Correlation endpoint order is retained, but
// self-pairs and duplicate unordered pairs are rejected.  On failure false is
// returned and, if supplied, error receives a concise diagnostic.
bool read( const std::string & file, model_t * model, std::string * error = nullptr );

// Read a post-fit descriptor and overlay its reviewed component semantics on
// a model. The model reference must match exactly; failure leaves the model
// unchanged.
bool read_descriptor( const std::string & file, model_t * model,
                      std::string * error = nullptr );

}

#endif
