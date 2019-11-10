
/**
 * @file matrix.h
 * 
 * Matrix/vector operations.
 */

#ifndef __ICA_MATRIX_H__
#define __ICA_MATRIX_H__

#include "stats/matrix.h"

#include <cfloat>
#define SCAL_EPSILON	DBL_EPSILON

typedef double scal;

// these resize too
void mat_zeroize(Data::Matrix<double> & , const int rows = 0 , const int cols = 0);
void vect_zeroize(Data::Vector<double> & , const int cols = 0 );

void vect_apply_fx(Data::Vector<double> & v, scal (*fx)(double,double), scal par);
void mat_apply_fx(Data::Matrix<double> & M, scal (*fx)(double,double), scal par);
void mat_mean_rows(Data::Matrix<double> & M, Data::Vector<double> & v);
scal mat_max_diag(Data::Matrix<double> & M );
scal mat_max_abs_diag(Data::Matrix<double> & M );
void mat_diag(Data::Vector<double> & v, Data::Matrix<double> & R);
void mat_transpose(Data::Matrix<double> & M, Data::Matrix<double> & R);
void mat_inverse(Data::Matrix<double> & M, Data::Matrix<double> & R);
void mat_sub(Data::Matrix<double> & A, Data::Matrix<double> & B, Data::Matrix<double> & R);
void mat_mult(Data::Matrix<double> & A, Data::Matrix<double> & B, Data::Matrix<double> & R);
void mat_center( Data::Matrix<double> & M,  Data::Vector<double> & means );
void mat_decenter(Data::Matrix<double> & M, Data::Vector<double> & means);

#endif /*MATRIX_H_*/
