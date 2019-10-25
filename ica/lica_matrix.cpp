/**
 * @file matrix.c
 * 
 * Matrix/vector operations.
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "ica/lica_matrix.h"
#include "helper/helper.h"


#define min(a,b) ((a<b) ? a : b)
#define abs(x) ((x)<0 ? -(x) : (x))



void mat_zeroize( Data::Matrix<double> & M , const int rows , const int cols )
{

  if ( rows != 0 || cols != 0 )
    M.resize( rows , cols );
  
  const int mr = M.dim1();
  const int mc = M.dim2();

  for(int i=0; i<mr; i++)
    for(int j=0; j<mc; j++)
      M[i][j] = 0;
}

void vect_zeroize( Data::Vector<double> & v , const int cols )
{
  if ( cols != 0 ) v.resize( cols );
  const int n = v.size();
  for (int i=0; i<n; i++) v[i] = 0;
}


/**
 * Applies function fx() with parameter par on
 * every vector element.
 */

void vect_apply_fx( Data::Vector<double> & v, scal (*fx)(double,double), scal par)
{
  const int n = v.size();
  for (int i=0; i<n; i++)
    v[i] = (*fx)(v[i], par);
}

/**
 * Applies function fx() with parameter par on
 * every matrix element.
 */
void mat_apply_fx( Data::Matrix<double> & M, scal (*fx)(double,double), scal par)
{
  const int rows = M.dim1();
  const int cols = M.dim2();
  
  for (int i=0; i<rows; i++)
    for (int j=0; j<cols; j++)
      M[i][j] = (*fx)(M[i][j], par);
}

/**
 * Computes the mean of every row.
 */
void mat_mean_rows( Data::Matrix<double> & M, Data::Vector<double> & v)
{
  const int rows = M.dim1();
  const int cols = M.dim2();
  
  scal sum;
  
  for (int i=0; i< rows; i++) {
    sum = 0;
    for (int j=0; j<cols; j++)
      sum += M[i][j];
    v[i] = sum / cols;
  }
}

/**
 * Returns the maximal element on the diagonal
 * of the matrix M.
 */
scal mat_max_diag( Data::Matrix<double> & M )
{

  const int rows = M.dim1();
  scal max = M[0][0];
  for (int i=1; i<rows; i++)
    if (M[i][i] > max)
      max = M[i][i];
  return max;
}

/**
 * Returns the maximal absolute element on the diagonal
 * of the matrix M.
 */
scal mat_max_abs_diag(Data::Matrix<double> & M )
{
  const int rows = M.dim1();
  scal max = fabs( M[0][0] ) ;
  for (int i=1; i<rows; i++)
    if ( fabs( M[i][i] ) > max)
      max = fabs( M[i][i] );
  return max;
}

/**
 * Creates a diagonal matrix from vector v.
 */
void mat_diag( Data::Vector<double> & v , Data::Matrix<double> & R )
{
  const int n = v.size();
  mat_zeroize( R );
  for (int i=0; i<n; i++)
    R[i][i] = v[i];
}

/**
 * Transponse matrix M.
 */
void mat_transpose( Data::Matrix<double> & M , Data::Matrix<double> & R )
{
  const int rows = M.dim1();
  const int cols = M.dim2();
  for (int i=0; i<rows; i++)
    for(int j=0; j<cols; j++)
      R[j][i] = M[i][j]; 
}

/**
 * Compute inverse matrix. R = M^(-1).
 */

void mat_inverse( Data::Matrix<double> & M,  Data::Matrix<double> & R )
{

  const int dim = M.dim1();

  int maxrow = 0;

  scal temp;

  // Create M|E	
  Data::Matrix<double> T( dim , 2*dim );
  //  T = mat_create(dim, 2*dim);
  
  for (int i=0; i<dim; i++)
    for (int j=0; j<dim; j++)
      T[i][j] = M[i][j];
  for (int i=0; i<dim; i++)
    for (int j=dim; j<2*dim; j++)
      if (i == j - dim)
	T[i][j] = 1;
      else
	T[i][j] = 0;
  
  // Gauss-Jordan elimination
  for (int i=0; i<dim; i++) {
    maxrow = i;
    for (int j=i+1; j<dim; j++)
      if (abs(T[j][i]) > abs(T[maxrow][i]))
	maxrow = j;
    for (int j=0; j<2*dim; j++) {
      temp = T[i][j];
      T[i][j] = T[maxrow][j];
      T[maxrow][j] = temp; 
    }
    if (abs(T[i][i]) <= SCAL_EPSILON) {
      Helper::halt( "ica_t:: inversion error, singular matrix");      
    }
    for (int j=i+1; j<dim; j++) {
      temp = T[j][i] / T[i][i];
      for (int k=i; k<2*dim; k++)
	T[j][k] -= T[i][k] * temp;
    }
  }
  for (int i=dim-1; i>=0; i--) {
    temp  = T[i][i];
    for (int j=0; j<i; j++)
      for (int k=2*dim-1; k>=i; k--)
	T[j][k] -=  T[i][k] * T[j][i] / temp;
    T[i][i] /= temp;
    for (int j=dim; j<2*dim; j++) {
      T[i][j] /= temp;
      R[i][j-dim] = T[i][j];
    }
  }
    
}


/**
 * Matrix subtraction. R = A - B.
 */
void mat_sub( Data::Matrix<double> & A, Data::Matrix<double> & B,  Data::Matrix<double> & R)
{
  const int rows = A.dim1();
  const int cols = A.dim2();
  for (int i=0; i<rows; i++)
    for (int j=0; j<cols; j++)
      R[i][j] = A[i][j] - B[i][j];
}

/**
 * Matrix multiplication. R = A * B.
 */
void mat_mult( Data::Matrix<double> & A, Data::Matrix<double> & B, Data::Matrix<double> &  R )
{
  const int rows_A = A.dim1();
  const int cols_A = A.dim2();
  const int cols_B = B.dim2();

  mat_zeroize(R , rows_A, cols_B);

  for(int i=0; i<rows_A; i++ )
    for(int j=0; j<cols_B; j++)
      for(int k=0; k<cols_A; k++)
	R[i][j] += A[i][k] * B[k][j];
}


/**
 * Centers mat M. (Subtracts the mean from every column)
 */

void mat_center( Data::Matrix<double> & M , Data::Vector<double> & means)
{

  int rows = M.dim1();
  int cols = M.dim2();
  
  vect_zeroize( means, cols );
  
  for (int i=0; i<rows; i++)
    for(int j=0; j<cols; j++)
      means[j] += M[i][j];		
  for (int i=0; i<cols; i++)
    means[i] /= rows; 
  for (int i=0; i<rows; i++)
    for(int j=0; j<cols; j++)
      M[i][j] -= means[j];	
}

/**
 * De-center mat M. (Adds the mean to every column)
 */
void mat_decenter( Data::Matrix<double> & M , Data::Vector<double> &  means )
{
  const int rows = M.dim1();
  const int cols = M.dim2();
  
  for(int i=0; i<rows; i++)
    for(int j=0; j<cols; j++)
      M[i][j] += means[j]; 
}
