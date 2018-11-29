/**
 * @file matrix.h
 * 
 * Matrix/vector operations.
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include <float.h>
#define SCAL_EPSILON	DBL_EPSILON

typedef double **mat;
typedef double *vect;
typedef double scal;

mat mat_create(int rows, int cols);
void mat_delete(mat M, int rows, int cols);
vect vect_create(int n);
void vect_delete(vect v);
void mat_zeroize(mat M, int rows, int cols);
void vect_zeroize(vect v, int n);

void mat_copy(mat M, int rows, int cols, mat Md);
void vect_apply_fx(vect v, int n, scal (*fx)(), scal par);
void mat_apply_fx(mat M, int rows, int cols, scal (*fx)(), scal par);
void mat_mean_rows(mat M, int rows, int cols, vect v);
scal mat_max_diag(mat M, int rows, int cols);
void mat_diag(vect v, int n, mat R);
void mat_transpose(mat M, int rows, int cols, mat R);
void mat_inverse(mat M, int dim, mat R);
void mat_sub(mat A, mat B, int rows, int cols, mat R);
void mat_mult(mat A, int rows_A, int cols_A, mat B, int rows_B, int cols_B, mat R);
void mat_center(mat M, int rows, int cols, vect means);
void mat_decenter(mat M, int rows, int cols, vect means);

#endif /*MATRIX_H_*/
