//    --------------------------------------------------------------------
//
//    This file is part of Luna.
//
//    --------------------------------------------------------------------

#ifndef __LUNA_DPP_VECTOR_H__
#define __LUNA_DPP_VECTOR_H__

#include <string>
#include <vector>

struct edf_t;
struct param_t;

// Separate low-rate/vector path for DPP.  Input channels are already
// window-level observations (e.g. SleepFM tokens), so this path does not
// apply raw-signal DSP or assume Fs >= 1 Hz.
namespace dpp_vector
{
  struct layout_t
  {
    int embedding_dim;
    int raw_offset;
    int context_offset;
    int geom_offset;
    int dyn_offset;
    bool raw;
    bool context;
    bool geom;
    bool dyn;
  };

  bool enabled( const param_t & );

  // Column layout shared by vector extraction and the two-level subject
  // summarizer.  The corpus format is intentionally numeric-only, so the
  // fit side reconstructs this layout from embedding-dim= and
  // vector-features=.
  layout_t layout( int embedding_dim , const param_t & );

  // Extract vector rows, optionally write a DPP corpus, and optionally apply
  // a fitted DPP bundle.  Returns true when vector mode was selected.
  bool run( edf_t & , param_t & );

  // Generic labels used by the cohort fit/apply manifest.  The vector input
  // channel names are deliberately not baked into the binary corpus format;
  // channel ordering is validated at extraction time and the feature count
  // is validated by the fitted bundle.
  std::vector<std::string> labels( int n_channels , const param_t & );
}

#endif
