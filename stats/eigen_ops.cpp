
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
#include "miscmath/miscmath.h"
#include "stats/statistics.h"
#include <vector>
#include "miscmath/crandom.h"


// nb. using Eigen:::Ref<>
// for a writable reference:    Eigen::Ref<Eigen::VectorXd> 
// for a const ref:             const Eigen::Ref<const Eigen::VectorXd> & 

std::vector<double> eigen_ops::copy_vector( const Eigen::VectorXd & e ) 
{
  std::vector<double> v( &e[0] , e.data() + e.size() );
  return v;
}

std::vector<double> eigen_ops::copy_array( const Eigen::ArrayXd & e ) 
{
  std::vector<double> v( &e[0] , e.data() + e.size() );
  return v;
}

Eigen::ArrayXd eigen_ops::copy_array( const std::vector<double> & e ) 
{
  Eigen::ArrayXd v = Eigen::ArrayXd::Zero( e.size() );
  for ( int i=0; i<e.size(); i++) v[i] = e[i];
  return v;
}

void eigen_ops::random_normal( Eigen::MatrixXd & M )
{
  const int rows = M.rows();
  const int cols = M.cols();
  for (int r = 0 ; r < rows ; r++ )
    for (int c = 0 ; c < cols ; c++)
      M(r,c) = Statistics::ltqnorm( CRandom::rand() );
}

void eigen_ops::scale( Eigen::Ref<Eigen::MatrixXd> M , bool normalize )
{
  const int N = M.rows();

  Eigen::Array<double, 1, Eigen::Dynamic> means = M.colwise().mean();

  if ( normalize )
    {
      Eigen::Array<double, 1, Eigen::Dynamic> sds = ((M.array().rowwise() - means ).square().colwise().sum()/(N-1)).sqrt();
      M.array().rowwise() -= means;
      // for (int i=0;i<sds.size();i++) 
      // 	if ( sds[i] == 0 ) sds[i] = 1;
      M.array().rowwise() /= sds;
    }
  else
    {
      M.array().rowwise() -= means;
    }
}


void eigen_ops::robust_scale( Eigen::Ref<Eigen::MatrixXd> m , double w )
{
  // 1) winsorize at +/- w 
  // 2) 
  const int rows = m.rows();
  const int cols = m.cols();
  for (int c=0;c<cols;c++)
    {
      std::vector<double> v = copy_vector( m.col(c) );
      double median = MiscMath::median( v );
      double iqr = MiscMath::iqr( v );
      
      double robust_sd = 0.7413 * iqr;

     
      // winsorize?

      if ( w > 0 )
	{
	  double lwr = MiscMath::percentile( v , w );
	  double upr = MiscMath::percentile( v , 1-w );
	  if ( lwr >= upr ) 
	    Helper::halt( "cannot robust_scale().. pls fix me" );

	  for (int i=0; i<rows; i++)
	    {	      
	      if      ( m(i,c) < lwr ) m(i,c) = lwr;
	      else if ( m(i,c) > upr ) m(i,c) = upr;
	    }
	}

      // median / IQR normalize
      for (int i=0; i<rows; i++)	
	m(i,c) = ( m(i,c) - median ) / robust_sd;
      
    }

  // hmm... a bad idea, or unnecessary?
  
  // finally, also scale by mean/variance too, just to ensure correct 
  // overall scale

  scale( m , true );

}



double eigen_ops::sdev( const Eigen::VectorXd & x )
{
  const int N = x.size();
  double mean = x.mean();
  return sqrt( ( x.array() - mean ).square().sum()/(N-1) );
}



// unit_scale, specifying min/max and truncating at 0 and 1
Eigen::VectorXd eigen_ops::unit_scale( const Eigen::VectorXd & x , double xmin , double xmax )
{
  const int n = x.size();
  if ( n == 0 ) return x;
  if ( xmin >= xmax ) return x;
  
  Eigen::VectorXd r( n );
  for (int i=0;i<n;i++) 
    {
      if ( x[i] <= xmin ) r[i] = 0;
      else if ( x[i] >= xmax ) r[i] = 1;
      else r[i] = ( x[i] - xmin ) / ( xmax - xmin );
    }
  return r;

}


Eigen::VectorXd eigen_ops::unit_scale( const Eigen::VectorXd & x )
{

  const int n = x.size();
  if ( n == 0 ) return x;

  double xmin = x[0] , xmax = x[0];
  for (int i=0;i<n;i++)
    {
      if ( x[i] < xmin ) xmin = x[i];
      else if ( x[i] > xmax ) xmax = x[i];
    }  

  if ( xmin == xmax ) return x;

  Eigen::VectorXd r( n );
  for (int i=0;i<n;i++) r[i] = ( x[i] - xmin ) / ( xmax - xmin );
  return r;
}



// // apply function fx() with parameter param, to each matrix element

// void eigen_ops::apply_fx( Eigen::MatrixXd & M, double (*fx)(double,double), double param)
// {
//   const int rows = M.rows();
//   const int cols = M.cols();  
//   for (int i=0; i<rows; i++)
//     for (int j=0; j<cols; j++)
//       M[i][j] = (*fx)(M[i][j], param);
// }


// // row means
// void eigen_ops::mat_mean_rows( Eigen::MatrixXd & M, Eigen::ArrayXd & v)
// {

//   const int rows = M.dim1();
//   const int cols = M.dim2();
  
//   double sum;
  
//   for (int i=0; i< rows; i++) {
//     sum = 0;
//     for (int j=0; j<cols; j++)
//       sum += M[i][j];
//     v[i] = sum / cols;
//   }

// }

//  // * Returns the maximal element on the diagonal
//  // * of the matrix M.

// double eigen_ops::mat_max_diag( Eigen::MatrixXd & M )
// {

//   const int rows = M.dim1();
//   double max = M[0][0];
//   for (int i=1; i<rows; i++)
//     if (M[i][i] > max)
//       max = M[i][i];
//   return max;
// }

//  // Returns the maximal absolute element on the diagonal
//  // of the matrix M.
 
// double eigen_ops::mat_max_abs_diag(Eigen::MatrixXd & M )
// {
//   const int rows = M.dim1();
//   double max = fabs( M[0][0] ) ;
//   for (int i=1; i<rows; i++)
//     if ( fabs( M[i][i] ) > max)
//       max = fabs( M[i][i] );
//   return max;
// }

// // Creates a diagonal matrix from vector v.
 
// void eigen_ops::mat_diag( Eigen::ArrayXd & v , Eigen::MatrixXd & R )
// {
//   const int n = v.size();
//   mat_zeroize( R );
//   for (int i=0; i<n; i++)
//     R[i][i] = v[i];
// }

// // Transponse matrix M.

// void eigen_ops::mat_transpose( Eigen::MatrixXd & M , Eigen::MatrixXd & R )
// {
//   const int rows = M.dim1();
//   const int cols = M.dim2();
//   for (int i=0; i<rows; i++)
//     for(int j=0; j<cols; j++)
//       R[j][i] = M[i][j]; 
// }

// // Centers mat M. (Subtracts the mean from every column)

// void eigen_ops::mat_center( Eigen::MatrixXd & M )
// {
//   //
// }

// void eigen_ops::mat_center( Eigen::MatrixXd & M , Eigen::ArrayXd & means )
// {

//   // int rows = M.dim1();
//   // int cols = M.dim2();
  
//   // vect_zeroize( means, cols );
  
//   // for (int i=0; i<rows; i++)
//   //   for(int j=0; j<cols; j++)
//   //     means[j] += M[i][j];		
//   // for (int i=0; i<cols; i++)
//   //   means[i] /= rows; 
//   // for (int i=0; i<rows; i++)
//   //   for(int j=0; j<cols; j++)
//   //     M[i][j] -= means[j];	
// }


// void eigen_op::mat_decenter( Eigen::MatrixXd & M , Eigen::ArrayXd &  means )
// {

//   const int rows = M.dim1();
//   const int cols = M.dim2();
  
//   for(int i=0; i<rows; i++)
//     for(int j=0; j<cols; j++)
//       M[i][j] += means[j]; 
// }



