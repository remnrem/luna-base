
//    --------------------------------------------------------------------
//
//    This file is part of Luna.
//
//    LUNA is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    Luna is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with Luna. If not, see <http://www.gnu.org/licenses/>.
//
//    Please see LICENSE.txt for more details.
//
//    --------------------------------------------------------------------

#include "stats/eigen_ops.h"
#include "stats/statistics.h"

#include "miscmath/crandom.h"

void eigen_ops::random_normal( Eigen::MatrixXd & M )
{
  const int rows = M.rows();
  const int cols = M.cols();
  for (int r = 0 ; r < rows ; r++ )
    for (int c = 0 ; c < cols ; c++)
      M(r,c) = Statistics::ltqnorm( CRandom::rand() );
}

void eigen_ops::scale( Eigen::MatrixXd & M , bool normalize )
{
  const int N = M.rows();

  Eigen::Array<double, 1, Eigen::Dynamic> means = M.colwise().mean();
  
  if ( normalize )
    {
      Eigen::Array<double, 1, Eigen::Dynamic> sds = ((M.array().rowwise() - means ).square().colwise().sum()/(N-1)).sqrt();
      M.array().rowwise() -= means;
      M.array().rowwise() /= sds;
    }
  else
    {
      M.array().rowwise() -= means;
    }
}


/* 

void eigen_ops::zeroize( Eigen::MatrixXd & M , const int rows , const int cols )
{
  M = Eigen::MatrixXd::Zero( rows , cols );
}

void eigen_ops::random_normal( Eigen::MatrixXd & M , const int rows , const int cols )
{
  M = Eigen::MatrixXd::Random( rows , cols );
  
  
}

void eigen_ops::vect_zeroize( Eigen::ArrayXd & v , const int cols )
{
  v = Eigen::ArrayXd::Zero( cols );
}


// apply function fx() with parameter param, to each array element

void eigen_ops::apply_fx( Eigen::ArrayXd & v, double (*fx)(double,double), double param)
{
  const int n = v.size();

  for (int i=0; i<n; i++)
    v[i] = (*fx)(v[i], param);

}

// apply function fx() with parameter param, to each matrix element

void eigen_ops::apply_fx( Eigen::MatrixXd & M, double (*fx)(double,double), double param)
{
  const int rows = M.rows();
  const int cols = M.cols();  
  for (int i=0; i<rows; i++)
    for (int j=0; j<cols; j++)
      M[i][j] = (*fx)(M[i][j], param);
}


// row means
void eigen_ops::mat_mean_rows( Eigen::MatrixXd & M, Eigen::ArrayXd & v)
{

  const int rows = M.dim1();
  const int cols = M.dim2();
  
  double sum;
  
  for (int i=0; i< rows; i++) {
    sum = 0;
    for (int j=0; j<cols; j++)
      sum += M[i][j];
    v[i] = sum / cols;
  }

}

 // * Returns the maximal element on the diagonal
 // * of the matrix M.

double eigen_ops::mat_max_diag( Eigen::MatrixXd & M )
{

  const int rows = M.dim1();
  double max = M[0][0];
  for (int i=1; i<rows; i++)
    if (M[i][i] > max)
      max = M[i][i];
  return max;
}

 // Returns the maximal absolute element on the diagonal
 // of the matrix M.
 
double eigen_ops::mat_max_abs_diag(Eigen::MatrixXd & M )
{
  const int rows = M.dim1();
  double max = fabs( M[0][0] ) ;
  for (int i=1; i<rows; i++)
    if ( fabs( M[i][i] ) > max)
      max = fabs( M[i][i] );
  return max;
}

// Creates a diagonal matrix from vector v.
 
void eigen_ops::mat_diag( Eigen::ArrayXd & v , Eigen::MatrixXd & R )
{
  const int n = v.size();
  mat_zeroize( R );
  for (int i=0; i<n; i++)
    R[i][i] = v[i];
}

// Transponse matrix M.

void eigen_ops::mat_transpose( Eigen::MatrixXd & M , Eigen::MatrixXd & R )
{
  const int rows = M.dim1();
  const int cols = M.dim2();
  for (int i=0; i<rows; i++)
    for(int j=0; j<cols; j++)
      R[j][i] = M[i][j]; 
}

// Centers mat M. (Subtracts the mean from every column)

void eigen_ops::mat_center( Eigen::MatrixXd & M )
{
  //
}

void eigen_ops::mat_center( Eigen::MatrixXd & M , Eigen::ArrayXd & means )
{

  // int rows = M.dim1();
  // int cols = M.dim2();
  
  // vect_zeroize( means, cols );
  
  // for (int i=0; i<rows; i++)
  //   for(int j=0; j<cols; j++)
  //     means[j] += M[i][j];		
  // for (int i=0; i<cols; i++)
  //   means[i] /= rows; 
  // for (int i=0; i<rows; i++)
  //   for(int j=0; j<cols; j++)
  //     M[i][j] -= means[j];	
}


void eigen_op::mat_decenter( Eigen::MatrixXd & M , Eigen::ArrayXd &  means )
{

  const int rows = M.dim1();
  const int cols = M.dim2();
  
  for(int i=0; i<rows; i++)
    for(int j=0; j<cols; j++)
      M[i][j] += means[j]; 
}


*/
