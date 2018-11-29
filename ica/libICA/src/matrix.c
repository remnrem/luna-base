/**
 * @file matrix.c
 * 
 * Matrix/vector operations.
 */

#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

#define min(a,b) ((a<b) ? a : b)
#define abs(x) ((x)<0 ? -(x) : (x))

/**
 * Creates a matrix of given size.
 */
mat mat_create(int rows, int cols)
{
	mat M; int i;

	M = (mat) malloc(rows * sizeof(vect));
	if (M == NULL) {
		perror("Error allocating memory!");
		exit(-1);
	}
	for (i=0; i<rows; i++) {
		M[i] = (vect) malloc(cols * sizeof(scal));
		if (M[i] == NULL) {
			perror("Error allocating memory!");
			exit(-1);
		}
	}
	
	return M;
}

/**
 * Deletes matrix M.
 */
void mat_delete(mat M, int rows, int cols)
{
	int i;

	for (i=0; i<rows; i++)
		free(M[i]);
	free(M);
	M = NULL;
}

/**
 * Fills the whole matrix with 0's.
 */
void mat_zeroize(mat M, int rows, int cols)
{
	int i,j;

	for(i=0; i<rows; i++)
		for(j=0; j<cols; j++)
			M[i][j] = 0;
}

/**
 * Creates a vector of given size.
 */
vect vect_create(int n)
{
	vect v;

	v = (vect) malloc(n * sizeof(scal));
	if (v == NULL) {
		perror("Error allocating memory!");
		exit(-1);
	}	
	return v;
}

/**
 * Deletes vector v
 */
void vect_delete(vect v)
{
	free(v); v = NULL;
}

/**
 * Fills the whole vector with 0's.
 */
void vect_zeroize(vect v, int n)
{
	int i;

	for (i=0; i<n; i++)
		v[i] = 0;
}

/**
 * Applyies function fx() with parameter par on
 * every vector element.
 */
void vect_apply_fx(vect v, int n, scal (*fx)(), scal par)
{
	int i;

	for (i=0; i<n; i++)
		v[i] = (*fx)(v[i], par);
}

/**
 * Copyies mat M to Md.
 */
void mat_copy(mat M, int rows, int cols, mat Md)
{
	int i, j;

	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
			Md[i][j] = M[i][j];
}

/**
 * Applyies function fx() with parameter par on
 * every matrix element.
 */
void mat_apply_fx(mat M, int rows, int cols, scal (*fx)(), scal par)
{
	int i,j;

	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
			M[i][j] = (*fx)(M[i][j], par);
}

/**
 * Computes the mean of every row.
 */
void mat_mean_rows(mat M, int rows, int cols, vect v)
{
	int i,j;
	scal sum;

	for (i=0; i< rows; i++) {
		sum = 0;
		for (j=0; j<cols; j++)
			sum += M[i][j];
		v[i] = sum / cols;
	}
}

/**
 * Returns the maximal element on the diagonal
 * of the matrix M.
 */
scal mat_max_diag(mat M, int rows, int cols)
{
	int i; scal max;

	max = M[0][0];
	for (i=1; i<rows; i++)
		if (M[i][i] > max)
			max = M[i][i];
	return max;
}

/**
 * Creates a diagonal matrix from vector v.
 */
void mat_diag(vect v, int n, mat R)
{
	int i;

	mat_zeroize(R, n, n);
	for (i=0; i<n; i++)
		R[i][i] = v[i];
}

/**
 * Transponse matrix M.
 */
void mat_transpose(mat M, int rows, int cols, mat R)
{
	int i,j;

	for (i=0; i<rows; i++)
		for(j=0; j<cols; j++)
			R[j][i] = M[i][j]; 
}

/**
 * Compute inverse matrix. R = M^(-1).
 */
void mat_inverse(mat M, int dim, mat R)
{
	mat T; int i, j, k, maxrow;
	scal temp;

	// Create M|E	
	T = mat_create(dim, 2*dim);

	for (i=0; i<dim; i++)
		for (j=0; j<dim; j++)
			T[i][j] = M[i][j];
	for (i=0; i<dim; i++)
		for (j=dim; j<2*dim; j++)
			if (i == j - dim)
				T[i][j] = 1;
			else
				T[i][j] = 0;

	// Gauss-Jordan elimination
	for (i=0; i<dim; i++) {
		maxrow = i;
		for (j=i+1; j<dim; j++)
			if (abs(T[j][i]) > abs(T[maxrow][i]))
				maxrow = j;
		for (j=0; j<2*dim; j++) {
			temp = T[i][j];
			T[i][j] = T[maxrow][j];
			T[maxrow][j] = temp; 
		}
		if (abs(T[i][i]) <= SCAL_EPSILON) {
			perror("Inversion error: Singular matrix!");
			exit(-1);
		}
		for (j=i+1; j<dim; j++) {
			temp = T[j][i] / T[i][i];
			for (k=i; k<2*dim; k++)
				T[j][k] -= T[i][k] * temp;
		}
	}
	for (i=dim-1; i>=0; i--) {
		temp  = T[i][i];
		for (j=0; j<i; j++)
			for (k=2*dim-1; k>=i; k--)
				T[j][k] -=  T[i][k] * T[j][i] / temp;
		T[i][i] /= temp;
		for (j=dim; j<2*dim; j++) {
			T[i][j] /= temp;
			R[i][j-dim] = T[i][j];
		}
	}

	mat_delete(T, dim, 2*dim);
}

/**
 * Matrix subtraction. R = A - B.
 */
void mat_sub(mat A, mat B, int rows, int cols, mat R)
{
	int i,j;

	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
			R[i][j] = A[i][j] - B[i][j];
}

/**
 * Matrix multiplication. R = A * B.
 */
void mat_mult(mat A, int rows_A, int cols_A, mat B, int rows_B, int cols_B, mat R)
{
	int i,j,k;

	mat_zeroize(R, rows_A, cols_B);
	for(i=0; i<rows_A; i++ )
		for(j=0; j<cols_B; j++)
			for(k=0; k<cols_A; k++)
				R[i][j] += A[i][k] * B[k][j];
}

/**
 * Centers mat M. (Subtracts the mean from every column)
 */
void mat_center(mat M, int rows, int cols, vect means)
{
	int i,j;

	vect_zeroize(means, cols);

	for (i=0; i<rows; i++)
		for(j=0; j<cols; j++)
			means[j] += M[i][j];		
	for (i=0; i<cols; i++)
		means[i] /= rows; 
	for (i=0; i<rows; i++)
		for(j=0; j<cols; j++)
			M[i][j] -= means[j];	
}

/**
 * De-center mat M. (Adds the mean to every column)
 */
void mat_decenter(mat M, int rows, int cols, vect means)
{
	int i,j;

	for(i=0; i<rows; i++)
		for(j=0; j<cols; j++)
			M[i][j] += means[j]; 
}
