#ifndef SB_COARSEGRAIN_H
#define SB_COARSEGRAIN_H
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "c22_stats.h"
#include "helper_functions.h"

extern "C" { 
extern void sb_coarsegrain(const double y[], const int size, const char how[], const int num_groups, int labels[]);
}

#endif
