//
//  butterworth.h
//  
//
//  Created by Carl Henning Lubba on 23/09/2018.
//

#ifndef butterworth_h
#define butterworth_h

#include <stdio.h>


#ifdef __cplusplus
extern "C" {
#endif

  extern void butterworthFilter(const double y[], const int size, const int nPoles, const double W, double out[]);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* butterworth_h */
