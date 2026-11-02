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
  double weight = 1.0;              // applied to standardized values before PCA
  std::vector<double> loadings;
};

struct panel_t {
  std::string id;
  int n_input = 0;
  int n_fit = 0;
  double total_ss = 0.0;
  std::vector<std::string> declared_features;
  std::vector<feature_t> features;
  std::vector<double> singular_values; // unrotated basis, including in JSON V2/V3
  std::vector<double> score_sds;       // final raw score SDs; populated for JSON
  std::string rotation_method = "none";
  std::vector<std::vector<double>> rotation; // basis rows -> output columns
  // One-based axes selected before rotation; these index rotation rows for V2/V3.
  // Identity for ordinary models without pre-rotation pc-select.
  std::vector<int> source_components;
  // JSON V3 only: the rotated basis width (rows of rotation/length of
  // singular_values), which can exceed the final output width n_components
  // (rotation columns/score_sds length) once rotated discovery PCs have been
  // selected post-rotation. Zero for V1/V2, where it equals n_components.
  int n_basis = 0;
  // JSON V3 only: one-based indices into the n_basis rotated discovery PCs
  // that were kept as this panel's final output columns. Empty for V1/V2.
  std::vector<int> selected_components;
};

// JSON models are self-describing: panel/component coordinates are derived
// from the embedded effective specification, never from a compiled PC registry.
struct component_t {
  std::string id;
  std::string panel_id;
  int panel_pc = 0;                 // one-based column within panel loadings
  int source_panel_pc = 0;          // original axis; 0 for rotated mixtures
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
  std::string specification_json;   // canonical embedded v2 spec, for V1/V2/V3 alike
  std::string component_id_style = "legacy"; // legacy PANEL/1 or R-safe PANEL_1
  // True when the embedded v2 contract requests scores divided by their
  // training-sample SD. Legacy V3 models retain their historical raw scores.
  bool score_unit_variance = false;
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
// self-contained JSON V1 (PCA) or V2 (rotated) model embedding a v2 specification.
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
//
// read() above applies this automatically for any JSON model: a model file
// PANEL.model.json (or PANEL.json) picks up a sibling PANEL.desc.json in the
// same directory, if one exists, and the overlaid orientation is what every
// consumer of the resulting model_t uses, including real projection
// (hdstats_pc::project_components/project_panels) -- not display only. A
// missing sibling is not an error: every component keeps its default
// orientation of +1. A malformed sibling is an error: read() fails rather
// than silently projecting with the wrong sign.
bool read_descriptor( const std::string & file, model_t * model,
                      std::string * error = nullptr );

}

#endif
