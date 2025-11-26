//
//  splinefit.h
//  C_polished
//
//  Created by Carl Henning Lubba on 27/09/2018.
//  Copyright © 2018 Carl Henning Lubba. All rights reserved.
//

#ifndef splinefit_h
#define splinefit_h

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

extern int splinefit(const double *y, const int size, double *yOut);

#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* splinefit_h */
