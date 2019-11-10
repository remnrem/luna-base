/**
 * @file libICA.h
 * 
 * Main FastICA functions.
 */

#ifndef LIBICA_H_
#define LIBICA_H_

#define MAX_ITERATIONS   1000
#define TOLERANCE        0.0001
#define LIBICA_ALPHA     1

#include "stats/matrix.h"

void fastICA( Data::Matrix<double> & X , int compc, Data::Matrix<double> & K , Data::Matrix<double> & W , Data::Matrix<double> & A , Data::Matrix<double> & S );


#endif /*LIBICA_H_*/
