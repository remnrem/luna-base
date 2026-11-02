// HDSTATS level-1 feature files and structured level-2 PCA models.
#ifndef __HDSTATS_PC_H__
#define __HDSTATS_PC_H__

#include <map>
#include <ostream>
#include <string>
#include <vector>

struct param_t;

namespace hdstats_pc {

typedef std::map<std::string,double> feature_map_t;

struct panel_scores_t {
  std::string panel;
  std::vector<double> scores;
  int missing = 0;
};

struct feature_contribution_t {
  std::string feature;
  double value;
  double raw_value;
};

// A component coordinate from a fitted model. Semantic labels/orientations are
// optional post-fit descriptor data, not part of the fitting specification.
struct component_score_t {
  std::string id;
  std::string panel;
  int panel_pc;
  int orientation;
  std::string label;
  bool compatible;
  int missing;
  double score;
  std::vector<feature_contribution_t> contributions;
};

// Stable key in DESTRAT's combined wide-label form:
// VARIABLE.FAC_LVL.FAC_LVL, from a FAC/LVL;FAC/LVL stratum.
std::string key( const std::string & variable, const std::string & stratum );

// Level 1: write one tab-delimited key/value line for an individual.
void write_features( const std::string & directory,
                     const std::string & id,
                     const feature_map_t & features );

// Legacy V2 single-model projection. Structured V3 projection uses
// project_panels() so one unavailable panel does not suppress the others.
bool project( const std::string & model,
              const feature_map_t & features,
              std::vector<double> * scores,
              int * missing );

// Structured models project each panel independently. Legacy V2 models are
// returned as a single GENERIC panel.
std::vector<panel_scores_t> project_panels( const std::string & model,
                                           const feature_map_t & features );

// Projection of self-contained JSON models. IDs are PANEL_N coordinates until
// optional post-fit descriptor data provides reviewed semantics.
std::vector<component_score_t> project_components( const std::string & model,
                                                    const feature_map_t & features,
                                                    double min_loading = 0.10 );

// Write a model's PC registry as a TSV table. Each component row includes
// oriented loadings meeting the minimum absolute threshold.
void describe( const std::string & model, std::ostream & out, double min_loading = 0.10,
               const std::string & descriptor = "" );

// Level 2 fitting requires spec=<JSON>. The v2 layout owns panel membership
// and weights; nc=N optionally overrides every panel's default count.
// pc-select=<JSON> retains reviewed exploratory axes before rotation.
// rotate=T (default) performs deterministic group-quartimax within that subspace;
// rotate=F preserves ordinary PCA. JSON V2 embeds final loadings, score SDs, and
// rotation provenance; projection has no spec/selection/rotation-time dependency.
void fit( param_t & param );

}

// HYPNO-PC deliberately uses a distinct on-disk format from HDSTATS-PC.
// The numerical PCA implementation is shared, but the two model families
// must never be interchangeable.
namespace hypno_pc {

typedef std::map<std::string,double> feature_map_t;

void write_features( const std::string & directory,
                     const std::string & id,
                     const feature_map_t & features );

bool project( const std::string & model,
              const feature_map_t & features,
              std::vector<double> * scores,
              int * missing );

// Level 2: luna --hypno-pc data=all.hypc model=reference.txt nc=5
void fit( param_t & param );

}

#endif
