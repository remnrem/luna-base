// HDSTATS level-1 feature files and level-2 PCA model (V2 schema).
#ifndef __HDSTATS_PC_H__
#define __HDSTATS_PC_H__

#include <map>
#include <string>
#include <vector>

struct param_t;

namespace hdstats_pc {

typedef std::map<std::string,double> feature_map_t;

// Level 1: write one tab-delimited key/value line for an individual.
void write_features( const std::string & directory,
                     const std::string & id,
                     const feature_map_t & features );

// Level 1 projection against a fixed, text model. Returns false on a schema
// mismatch; missing contains the number of required absent/non-finite values.
bool project( const std::string & model,
              const feature_map_t & features,
              std::vector<double> * scores,
              int * missing );

// Level 2: luna --hdstats-pc data=all.txt model=reference.txt nc=5
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
