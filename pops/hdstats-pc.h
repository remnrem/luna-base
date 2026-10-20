// HDSTATS level-1 feature files and structured level-2 PCA models.
#ifndef __HDSTATS_PC_H__
#define __HDSTATS_PC_H__

#include <map>
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

// Stable key for a stdout-style VARIABLE and FAC/LVL;FAC/LVL stratum.
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

// Level 2: luna --hdstats-pc --options data=all.txt model=reference.txt nc=2
// with optional nc-<panel>=N overrides. Fitting uses hard-coded C++ panels.
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
