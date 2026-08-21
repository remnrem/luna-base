//    --------------------------------------------------------------------
//    Two-level vector DPP: cross-fitted local model plus subject summaries.
//    --------------------------------------------------------------------

#ifndef __LUNA_DPP_TWOLEVEL_H__
#define __LUNA_DPP_TWOLEVEL_H__

#include <string>

struct edf_t;
struct param_t;
struct dpp_matrix_t;

namespace dpp_twolevel
{
  bool enabled( const param_t & );
  void fit( param_t & );
  void apply( edf_t &, param_t &, const dpp_matrix_t & );
}

#endif
