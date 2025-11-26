#ifndef SC_FLUCTANAL
#define SC_FLUCTANAL
#include <math.h>
#include <string.h>
#include "c22_stats.h"
#include "CO_AutoCorr.h"

#ifdef __cplusplus
extern "C" {
#endif
  
extern double SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1(const double y[], const int size);
extern double SC_FluctAnal_2_50_1_logi_prop_r1(const double y[], const int size, const char how[]);
extern double SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1(const double y[], const int size);

#ifdef __cplusplus
} /* extern "C" */
#endif


#endif
