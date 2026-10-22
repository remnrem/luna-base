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

// Optional, append-only annotation records.  entity_type is DOMAIN, PC, or
// FEATURE.  PC_CORR and FTR_CORR preserve the order written in the file.
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
  std::vector<panel_t> panels;
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

// Reads a __HDSTATS_PC_V3 model and its optional trailer:
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

}

#endif
