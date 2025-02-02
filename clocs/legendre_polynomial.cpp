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

// This library downloaded from
// https://people.sc.fsu.edu/~jburkardt/cpp_src/legendre_polynomial/legendre_polynomial.html

// The computer code and data files described and made available on
// this web page are distributed under the GNU LGPL license.

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iomanip>
#include <iostream>

#include "dsp/r8lib.h"

using namespace std;

#include "legendre_polynomial.h"



Eigen::VectorXd legendre( const int N , double x )
{

  // note: not sure why this is structured this way... need to check
  // originally, this did not free allocated memory from
  // pm_polynomial_value()

  // returns N+1 vector 
  const int NP1 = N+1;
  
  Eigen::VectorXd L = Eigen::VectorXd::Zero( NP1 );
  
  for (int M=0;M<=N;M++)
    {

      // std::vector<double> ans( NP1 );
      // double * pl = &(ans)[0];      
      // pl = pm_polynomial_value( 1 , N , M , &x );
      // // return last elemnt
      // pl += N;
      // L[M] = *pl;
      
      // new version
      double * pl = pm_polynomial_value( 1 , N , M , &x );
      
      // return last elemnt
      double * pl_orig = pl;
      pl += N;
      L[M] = *pl;
      
      // clean-up
      if ( pl_orig != NULL )
	delete [] pl_orig;
    }
  return L;  
  
}


// Redundant:: original version that did not free up allocated memory
// std::vector<std::vector<double> > legendre( const int N , const std::vector<double> & x )
// {
//   // returns N+1 by x.size() matrix 
//   const int MM = x.size();
//   const int NP1 = N+1;
//   double * px = (double*)&(x[0]);
//   std::vector<std::vector<double> > L( NP1 );
//   for (int i=0;i<NP1;i++) L[i].resize( MM );
//   for (int M=0;M<=N;M++)
//     {
//       std::vector<double> ans( MM * NP1 );
//       double * pl = &(ans)[0];      
//       pl = pm_polynomial_value( MM , N , M , px );
//       // skip to last set
//       pl += MM*N;
//       for (int j=0;j< MM ; j++) 
// 	L[M][j] = *(pl++);
      
//     }
//   return L;  
// }



std::vector<Eigen::MatrixXd> legendre( const int N , const Eigen::MatrixXd & D )
{
  
  // return N vector, each element mirrors the original input matrix
  
  const int nr = D.rows();
  const int nc = D.cols();

  std::vector<Eigen::MatrixXd> R(N);
  for (int n=0;n<N;n++) 
    R[n] = Eigen::MatrixXd::Zero( nr,nc );
  
  // for each datapoint
  for (int n=1;n<=N;n++)
    for (int r=0;r<nr;r++)
      for (int c=0;c<nc;c++)
	{
	  Eigen::VectorXd y = legendre( n , D(r,c) );
	  R[n-1](r,c) = y[0];
	}
  
  return R;  
}




//****************************************************************************80

char digit_to_ch ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    DIGIT_TO_CH returns the base 10 digit character corresponding to a digit.
//
//  Example:
//
//     I     C
//   -----  ---
//     0    '0'
//     1    '1'
//   ...    ...
//     9    '9'
//    10    '*'
//   -83    '*'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the digit, which should be between 0 and 9.
//
//    Output, char DIGIT_TO_CH, the appropriate character '0'
//    through '9' or '*'.
//
{
  char c;

  if ( 0 <= i && i <= 9 )
  {
    c = '0' + i;
  }
  else
  {
    c = '*';
  }

  return c;
}
//****************************************************************************80

// int i4_log_10 ( int i )

// //****************************************************************************80
// //
// //  Purpose:
// //
// //    I4_LOG_10 returns the whole part of the logarithm base 10 of an I4.
// //
// //  Discussion:
// //
// //    It should be the case that 10^I4_LOG_10(I) <= |I| < 10^(I4_LOG_10(I)+1).
// //    (except for I = 0).
// //
// //    The number of decimal digits in I is I4_LOG_10(I) + 1.
// //
// //  Example:
// //
// //        I    I4_LOG_10(I)
// //
// //        0     0
// //        1     0
// //        2     0
// //
// //        9     0
// //       10     1
// //       11     1
// //
// //       99     1
// //      100     2
// //      101     2
// //
// //      999     2
// //     1000     3
// //     1001     3
// //
// //     9999     3
// //    10000     4
// //    10001     4
// //
// //  Licensing:
// //
// //    This code is distributed under the GNU LGPL license.
// //
// //  Modified:
// //
// //    17 June 2003
// //
// //  Author:
// //
// //    John Burkardt
// //
// //  Parameters:
// //
// //    Input, int I, the integer.
// //
// //    Output, int I4_LOG_10, the whole part of the logarithm of abs ( I ).
// //
// {
//   int ten_pow;
//   int value;

//   i = abs ( i );

//   ten_pow = 10;
//   value = 0;

//   while ( ten_pow <= i )
//   {
//     ten_pow = ten_pow * 10;
//     value = value + 1;
//   }

//   return value;
// }
// //****************************************************************************80

// int i4_max ( int i1, int i2 )

// //****************************************************************************80
// //
// //  Purpose:
// //
// //    I4_MAX returns the maximum of two I4's.
// //
// //  Licensing:
// //
// //    This code is distributed under the GNU LGPL license.
// //
// //  Modified:
// //
// //    13 October 1998
// //
// //  Author:
// //
// //    John Burkardt
// //
// //  Parameters:
// //
// //    Input, int I1, I2, are two integers to be compared.
// //
// //    Output, int I4_MAX, the larger of I1 and I2.
// //
// {
//   int value;

//   if ( i2 < i1 )
//   {
//     value = i1;
//   }
//   else
//   {
//     value = i2;
//   }
//   return value;
// }
// //****************************************************************************80

// int i4_min ( int i1, int i2 )

// //****************************************************************************80
// //
// //  Purpose:
// //
// //    I4_MIN returns the minimum of two I4's.
// //
// //  Licensing:
// //
// //    This code is distributed under the GNU LGPL license.
// //
// //  Modified:
// //
// //    13 October 1998
// //
// //  Author:
// //
// //    John Burkardt
// //
// //  Parameters:
// //
// //    Input, int I1, I2, two integers to be compared.
// //
// //    Output, int I4_MIN, the smaller of I1 and I2.
// //
// {
//   int value;

//   if ( i1 < i2 )
//   {
//     value = i1;
//   }
//   else
//   {
//     value = i2;
//   }
//   return value;
// }

//****************************************************************************80

string i4_to_s ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_S converts an I4 to a string.
//
//  Example:
//
//    INTVAL  S
//
//         1  1
//        -1  -1
//         0  0
//      1952  1952
//    123456  123456
//   1234567  1234567
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, an integer to be converted.
//
//    Output, string I4_TO_S, the representation of the integer.
//
{
  int digit;
  int j;
  int length;
  int ten_power;
  string s;
  char s_char[80];
  static double ten = 10.0;

  length = i4_log_10 ( i );

  ten_power = ( int ) ( pow ( ten, length ) );

  if ( i < 0 )
  {
    length = length + 1;
  }
//
//  Add one position for the trailing null.
//
  length = length + 1;

  if ( i == 0 )
  {
    s_char[0] = '0';
    s_char[1] = '\0';
    s = string ( s_char );
    return s;
  }
//
//  Now take care of the sign.
//
  j = 0;
  if ( i < 0 )
  {
    s_char[j] = '-';
    j = j + 1;
    i = abs ( i );
  }
//
//  Find the leading digit of I, strip it off, and stick it into the string.
//
  while ( 0 < ten_power )
  {
    digit = i / ten_power;
    s_char[j] = digit_to_ch ( digit );
    j = j + 1;
    i = i - digit * ten_power;
    ten_power = ten_power / 10;
  }
//
//  Tack on the trailing NULL.
//
  s_char[j] = '\0';
  j = j + 1;

  s = string ( s_char );

  return s;
}
//****************************************************************************80

void imtqlx ( int n, double d[], double e[], double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    IMTQLX diagonalizes a symmetric tridiagonal matrix.
//
//  Discussion:
//
//    This routine is a slightly modified version of the EISPACK routine to 
//    perform the implicit QL algorithm on a symmetric tridiagonal matrix. 
//
//    The authors thank the authors of EISPACK for permission to use this
//    routine. 
//
//    It has been modified to produce the product Q' * Z, where Z is an input 
//    vector and Q is the orthogonal matrix diagonalizing the input matrix.  
//    The changes consist (essentialy) of applying the orthogonal transformations
//    directly to Z as they are generated.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//    Roger Martin, James Wilkinson,
//    The Implicit QL Algorithm,
//    Numerische Mathematik,
//    Volume 12, Number 5, December 1968, pages 377-383.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input/output, double D(N), the diagonal entries of the matrix.
//    On output, the information in D has been overwritten.
//
//    Input/output, double E(N), the subdiagonal entries of the 
//    matrix, in entries E(1) through E(N-1).  On output, the information in
//    E has been overwritten.
//
//    Input/output, double Z(N).  On input, a vector.  On output,
//    the value of Q' * Z, where Q is the matrix that diagonalizes the
//    input symmetric tridiagonal matrix.
//
{
  double b;
  double c;
  double f;
  double g;
  int i;
  int ii;
  int itn = 30;
  int j;
  int k;
  int l;
  int m;
  int mml;
  double p;
  double prec;
  double r;
  double s;

  prec = r8_epsilon ( );

  if ( n == 1 )
  {
    return;
  }

  e[n-1] = 0.0;

  for ( l = 1; l <= n; l++ )
  {
    j = 0;
    for ( ; ; )
    {
      for ( m = l; m <= n; m++ )
      {
        if ( m == n )
        {
          break;
        }

        if ( fabs ( e[m-1] ) <= prec * ( fabs ( d[m-1] ) + fabs ( d[m] ) ) )
        {
          break;
        }
      }
      p = d[l-1];
      if ( m == l )
      {
        break;
      }
      if ( itn <= j )
      {
        cout << "\n";
        cout << "IMTQLX - Fatal error!\n";
        cout << "  Iteration limit exceeded\n";
        exit ( 1 );
      }
      j = j + 1;
      g = ( d[l] - p ) / ( 2.0 * e[l-1] );
      r =  sqrt ( g * g + 1.0 );
      g = d[m-1] - p + e[l-1] / ( g + fabs ( r ) * r8_sign ( g ) );
      s = 1.0;
      c = 1.0;
      p = 0.0;
      mml = m - l;

      for ( ii = 1; ii <= mml; ii++ )
      {
        i = m - ii;
        f = s * e[i-1];
        b = c * e[i-1];

        if ( fabs ( g ) <= fabs ( f ) )
        {
          c = g / f;
          r =  sqrt ( c * c + 1.0 );
          e[i] = f * r;
          s = 1.0 / r;
          c = c * s;
        }
        else
        {
          s = f / g;
          r =  sqrt ( s * s + 1.0 );
          e[i] = g * r;
          c = 1.0 / r;
          s = s * c;
        }
        g = d[i] - p;
        r = ( d[i-1] - g ) * s + 2.0 * c * b;
        p = s * r;
        d[i] = g + p;
        g = c * r - b;
        f = z[i];
        z[i] = s * z[i-1] + c * f;
        z[i-1] = c * z[i-1] - s * f;
      }
      d[l-1] = d[l-1] - p;
      e[l-1] = g;
      e[m-1] = 0.0;
    }
  }
//
//  Sorting.
//
  for ( ii = 2; ii <= m; ii++ )
  {
    i = ii - 1;
    k = i;
    p = d[i-1];

    for ( j = ii; j <= n; j++ )
    {
      if ( d[j-1] < p )
      {
         k = j;
         p = d[j-1];
      }
    }

    if ( k != i )
    {
      d[k-1] = d[i-1];
      d[i-1] = p;
      p = z[i-1];
      z[i-1] = z[k-1];
      z[k-1] = p;
    }
  }
  return;
}
//****************************************************************************80

double *p_exponential_product ( int p, double b )

//****************************************************************************80
//
//  Purpose:
//
//    P_EXPONENTIAL_PRODUCT: exponential products for P(n,x).
//
//  Discussion:
//
//    Let P(n,x) represent the Legendre polynomial of degree n.  
//
//    For polynomial chaos applications, it is of interest to know the
//    value of the integrals of products of exp(B*X) with every possible pair
//    of basis functions.  That is, we'd like to form
//
//      Tij = Integral ( -1.0 <= X <= +1.0 ) exp(B*X) * P(I,X) * P(J,X) dx
//
//    We will estimate these integrals using Gauss-Legendre quadrature.
//    Because of the exponential factor exp(B*X), the quadrature will not 
//    be exact.
//
//    However, when B = 0, the quadrature is exact, and moreoever, the
//    table will be the identity matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P, the maximum degree of the polyonomial 
//    factors.  0 <= P.
//
//    Input, double B, the coefficient of X in the exponential factor.
//
//    Output, double P_EXPONENTIAL_PRODUCT[(P+1)*(P+1)], the table of integrals.  
//
{
  double *h_table;
  int i;
  int j;
  int k;
  int order;
  double *table;
  double *w_table;
  double x;
  double *x_table;

  table = new double[(p+1)*(p+1)];

  for ( j = 0; j <= p; j++ )
  {
    for ( i = 0; i <= p; i++ )
    {
      table[i+j*(p+1)] = 0.0;
    }
  }

  order = ( 3 * p + 4 ) / 2;

  x_table = new double[order];
  w_table = new double[order];

  p_quadrature_rule ( order, x_table, w_table );

  for ( k = 0; k < order; k++ )
  {
    x = x_table[k];
    h_table = p_polynomial_value ( 1, p, x_table + k );
//
//  The following formula is an outer product in H_TABLE.
//
    for ( j = 0; j <= p; j++ )
    {
      for ( i = 0; i <= p; i++ )
      {
        table[i+j*(p+1)] = table[i+j*(p+1)] 
          + w_table[k] * exp ( b * x ) * h_table[i] * h_table[j];
      }
    }
    delete [] h_table;
  }

  delete [] w_table;
  delete [] x_table;

  return table;
}
//****************************************************************************80

double p_integral ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    P_INTEGRAL evaluates a monomial integral associated with P(n,x).
//
//  Discussion:
//
//    The integral:
//
//      integral ( -1 <= x < +1 ) x^n dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the exponent.
//    0 <= N.
//
//    Output, double P_INTEGRAL, the value of the integral.
//
{
  double value;

  if ( ( n % 2 ) == 1 )
  {
    value = 0.0;
  }
  else
  {
    value = 2.0 / ( double ) ( n + 1 );
  }
  return value;
}
//****************************************************************************80

double *p_polynomial_coefficients ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    P_POLYNOMIAL_COEFFICIENTS: coefficients of Legendre polynomial P(n,x).
//
//  Discussion:
//
//     1
//     0     1
//    -1/2   0      3/2
//     0    -3/2    0     5/2
//     3/8   0    -30/8   0     35/8
//     0    15/8    0   -70/8    0     63/8
//    -5/16  0    105/16  0   -315/16   0    231/16
//     0   -35/16   0   315/16   0   -693/16   0    429/16
//
//     1.00000
//     0.00000  1.00000
//    -0.50000  0.00000  1.50000
//     0.00000 -1.50000  0.00000  2.5000
//     0.37500  0.00000 -3.75000  0.00000  4.37500
//     0.00000  1.87500  0.00000 -8.75000  0.00000  7.87500
//    -0.31250  0.00000  6.56250  0.00000 -19.6875  0.00000  14.4375
//     0.00000 -2.1875   0.00000  19.6875  0.00000 -43.3215  0.00000  26.8125
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996.
//
//  Parameters:
//
//    Input, int N, the highest order polynomial to evaluate.
//    Note that polynomials 0 through N will be evaluated.
//
//    Output, double P_POLYNOMIAL_COEFFICIENTS[(N+1)*(N+1)], the coefficients of 
//    the Legendre polynomials of degree 0 through N.
//
{
  double *c;
  int i;
  int j;
  double t;

  if ( n < 0 )
  {
    return NULL;
  }

  c = new double[(n+1)*(n+1)];

  for ( i = 0; i <= n; i++ )
  {
    for ( j = 0; j <= n; j++ )
    {
      c[i+j*(n+1)] = 0.0;
    }
  }

  c[0+0*(n+1)] = 1.0;

  if ( 0 < n )
  {
    c[1+1*(n+1)] = 1.0;
  }

  for ( i = 2; i <= n; i++ )
  {
    for ( j = 0; j <= i-2; j++ )
    {
      c[i+j*(n+1)] =          
          ( double ) (   - i + 1 ) * c[i-2+j*(n+1)] / ( double ) i;
    }
    for ( j = 1; j <= i; j++ )
    {
      c[i+j*(n+1)] = c[i+j*(n+1)] 
        + ( double ) ( i + i - 1 ) * c[i-1+(j-1)*(n+1)] / ( double ) i;
    }
  }

  return c;
}
//****************************************************************************80

double *p_polynomial_prime ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P_POLYNOMIAL_PRIME evaluates the derivative of Legendre polynomials P(n,x).
//
//  Discussion:
//
//    P(0,X) = 1
//    P(1,X) = X
//    P(N,X) = ( (2*N-1)*X*P(N-1,X)-(N-1)*P(N-2,X) ) / N
//
//    P'(0,X) = 0
//    P'(1,X) = 1
//    P'(N,X) = ( (2*N-1)*(P(N-1,X)+X*P'(N-1,X)-(N-1)*P'(N-2,X) ) / N
//
//    Thanks to Dimitriy Morozov for pointing out a memory leak caused by
//    not deleting the work array V before return, 19 March 2013.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 March 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996.
//
//  Parameters:
//
//    Input, int M, the number of evaluation points.
//
//    Input, int N, the highest order polynomial to evaluate.
//    Note that polynomials 0 through N will be evaluated.
//
//    Input, double X[M], the evaluation points.
//
//    Output, double P_POLYNOMIAL_PRIME[M*(N+1)], the values of the derivatives 
//    of the Legendre polynomials of order 0 through N at the points.
//
{
  int i;
  int j;
  double *v;
  double *vp;

  if ( n < 0 )
  {
    return NULL;
  }

  vp = new double[m*(n+1)];

  for ( i = 0; i < m; i++ )
  {
    vp[i+0*m] = 0.0;
  }

  if ( n < 1 )
  {
    return vp;
  }

  v = new double[m*(n+1)];

  for ( i = 0; i < m; i++ )
  {
    v[i+0*m] = 1.0;
  }

  for ( i = 0; i < m; i++ )
  {
    v[i+1*m] = x[i];
    vp[i+1*m] = 1.0;
  }

  for ( j = 2; j <= n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      v[i+j*m] = ( ( double ) ( 2 * j - 1 ) * x[i] * v[i+(j-1)*m]   
                 - ( double ) (     j - 1 ) *        v[i+(j-2)*m] ) 
                 / ( double ) (     j     );
 
      vp[i+j*m] = ( ( double ) ( 2 * j - 1 ) * ( v[i+(j-1)*m] + x[i] * vp[i+(j-1)*m] ) 
                  - ( double ) (     j - 1 ) *   vp[i+(j-2)*m]               ) 
                  / ( double ) (     j     );
    }
  }
 
  delete [] v;

  return vp;
}
//****************************************************************************80

double *p_polynomial_prime2 ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P_POLYNOMIAL_PRIME2: second derivative of Legendre polynomials P(n,x).
//
//  Discussion:
//
//    P(0,X) = 1
//    P(1,X) = X
//    P(N,X) = ( (2*N-1)*X*P(N-1,X)-(N-1)*P(N-2,X) ) / N
//
//    P'(0,X) = 0
//    P'(1,X) = 1
//    P'(N,X) = ( (2*N-1)*(P(N-1,X)+X*P'(N-1,X)-(N-1)*P'(N-2,X) ) / N
//
//    P"(0,X) = 0
//    P"(1,X) = 0
//    P"(N,X) = ( (2*N-1)*(2*P'(N-1,X)+X*P"(N-1,X)-(N-1)*P"(N-2,X) ) / N
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 May 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996.
//
//  Parameters:
//
//    Input, int M, the number of evaluation points.
//
//    Input, int N, the highest order polynomial to evaluate.
//    Note that polynomials 0 through N will be evaluated.
//
//    Input, double X[M], the evaluation points.
//
//    Output, double P_POLYNOMIAL_PRIME2[M*(N+1)], the second derivative
//    of the Legendre polynomials of order 0 through N at the points.
//
{
  int i;
  int j;
  double *v;
  double *vp;
  double *vpp;

  if ( n < 0 )
  {
    return NULL;
  }

  vpp = new double[m*(n+1)];

  for ( i = 0; i < m; i++ )
  {
    vpp[i+0*m] = 0.0;
  }

  if ( n < 1 )
  {
    return vpp;
  }

  v = new double[m*(n+1)];
  vp = new double[m*(n+1)];

  for ( i = 0; i < m; i++ )
  {
    v[i+0*m] = 1.0;
    vp[i+0*m] = 0.0;
  }

  for ( i = 0; i < m; i++ )
  {
    v[i+1*m] = x[i];
    vp[i+1*m] = 1.0;
    vpp[i+1*m] = 0.0;
  }

  for ( j = 2; j <= n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      v[i+j*m] = 
        ( ( double ) ( 2 * j - 1 ) * x[i] * v[i+(j-1)*m]   
        - ( double ) (     j - 1 ) *        v[i+(j-2)*m] ) 
        / ( double ) (     j     );
 
      vp[i+j*m] = 
        ( ( double ) ( 2 * j - 1 ) * ( v[i+(j-1)*m] + x[i] * vp[i+(j-1)*m] ) 
        - ( double ) (     j - 1 ) *   vp[i+(j-2)*m]               ) 
        / ( double ) (     j     );

      vpp[i+j*m] = 
        ( ( double ) ( 2 * j - 1 ) * ( 2.0 * vp[i+(j-1)*m] + x[i] * vpp[i+(j-1)*m] ) 
        - ( double ) (     j - 1 ) *   vpp[i+(j-2)*m]               ) 
        / ( double ) (     j     );
    }
  }
 
  delete [] v;
  delete [] vp;

  return vpp;
}


//****************************************************************************80

double *p_polynomial_value ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P_POLYNOMIAL_VALUE evaluates the Legendre polynomials P(n,x).
//
//  Discussion:
//
//    P(n,1) = 1.
//    P(n,-1) = (-1)^N.
//    | P(n,x) | <= 1 in [-1,1].
//
//    The N zeroes of P(n,x) are the abscissas used for Gauss-Legendre
//    quadrature of the integral of a function F(X) with weight function 1
//    over the interval [-1,1].
//
//    The Legendre polynomials are orthogonal under the inner product defined
//    as integration from -1 to 1:
//
//      Integral ( -1 <= X <= 1 ) P(I,X) * P(J,X) dX 
//        = 0 if I =/= J
//        = 2 / ( 2*I+1 ) if I = J.
//
//    Except for P(0,X), the integral of P(I,X) from -1 to 1 is 0.
//
//    A function F(X) defined on [-1,1] may be approximated by the series
//      C0*P(0,x) + C1*P(1,x) + ... + CN*P(n,x)
//    where
//      C(I) = (2*I+1)/(2) * Integral ( -1 <= X <= 1 ) F(X) P(I,x) dx.
//
//    The formula is:
//
//      P(n,x) = (1/2^N) * sum ( 0 <= M <= N/2 ) C(N,M) C(2N-2M,N) X^(N-2*M)
//
//  Differential equation:
//
//    (1-X*X) * P(n,x)'' - 2 * X * P(n,x)' + N * (N+1) = 0
//
//  First terms:
//
//    P( 0,x) =      1
//    P( 1,x) =      1 X
//    P( 2,x) = (    3 X^2 -       1)/2
//    P( 3,x) = (    5 X^3 -     3 X)/2
//    P( 4,x) = (   35 X^4 -    30 X^2 +     3)/8
//    P( 5,x) = (   63 X^5 -    70 X^3 +    15 X)/8
//    P( 6,x) = (  231 X^6 -   315 X^4 +   105 X^2 -     5)/16
//    P( 7,x) = (  429 X^7 -   693 X^5 +   315 X^3 -    35 X)/16
//    P( 8,x) = ( 6435 X^8 - 12012 X^6 +  6930 X^4 -  1260 X^2 +   35)/128
//    P( 9,x) = (12155 X^9 - 25740 X^7 + 18018 X^5 -  4620 X^3 +  315 X)/128
//    P(10,x) = (46189 X^10-109395 X^8 + 90090 X^6 - 30030 X^4 + 3465 X^2-63)/256
//
//  Recursion:
//
//    P(0,x) = 1
//    P(1,x) = x
//    P(n,x) = ( (2*n-1)*x*P(n-1,x)-(n-1)*P(n-2,x) ) / n
//
//    P'(0,x) = 0
//    P'(1,x) = 1
//    P'(N,x) = ( (2*N-1)*(P(N-1,x)+X*P'(N-1,x)-(N-1)*P'(N-2,x) ) / N
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996.
//
//  Parameters:
//
//    Input, int M, the number of evaluation points.
//
//    Input, int N, the highest order polynomial to evaluate.
//    Note that polynomials 0 through N will be evaluated.
//
//    Input, double X[M], the evaluation points.
//
//    Output, double P_POLYNOMIAL_VALUE[M*(N+1)], the values of the Legendre
//    polynomials of order 0 through N.
//
{
  int i;
  int j;
  double *v;

  if ( n < 0 )
  {
    return NULL;
  }

  v = new double[m*(n+1)];

  for ( i = 0; i < m; i++ )
  {
    v[i+0*m] = 1.0;
  }

  if ( n < 1 )
  {
    return v;
  }

  for ( i = 0; i < m; i++ )
  {
    v[i+1*m] = x[i];
  }
 
  for ( j = 2; j <= n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      v[i+j*m] = ( ( double ) ( 2 * j - 1 ) * x[i] * v[i+(j-1)*m]   
                 - ( double ) (     j - 1 ) *        v[i+(j-2)*m] ) 
                 / ( double ) (     j     );
    }
  }
 
  return v;
}
//****************************************************************************80

void p_polynomial_values ( int &n_data, int &n, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    P_POLYNOMIAL_VALUES: selected values of the Legendre polynomials P(n,x).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the function.
//
//    Output, double &X, the point where the function is evaluated.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 22

  static double fx_vec[N_MAX] = {
      0.1000000000000000E+01,
      0.2500000000000000E+00,
     -0.4062500000000000E+00,
     -0.3359375000000000E+00,
      0.1577148437500000E+00,
      0.3397216796875000E+00,
      0.2427673339843750E-01,
     -0.2799186706542969E+00,
     -0.1524540185928345E+00,
      0.1768244206905365E+00,
      0.2212002165615559E+00,
      0.0000000000000000E+00,
     -0.1475000000000000E+00,
     -0.2800000000000000E+00,
     -0.3825000000000000E+00,
     -0.4400000000000000E+00,
     -0.4375000000000000E+00,
     -0.3600000000000000E+00,
     -0.1925000000000000E+00,
      0.8000000000000000E-01,
      0.4725000000000000E+00,
      0.1000000000000000E+01 };

  static int n_vec[N_MAX] = {
     0,  1,  2,
     3,  4,  5,
     6,  7,  8,
     9, 10,  3,
     3,  3,  3,
     3,  3,  3,
     3,  3,  3,
     3 };

  static double x_vec[N_MAX] = {
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.00E+00,
     0.10E+00,
     0.20E+00,
     0.30E+00,
     0.40E+00,
     0.50E+00,
     0.60E+00,
     0.70E+00,
     0.80E+00,
     0.90E+00,
     1.00E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double *p_polynomial_zeros ( int nt )

//****************************************************************************80
//
//  Purpose:
//
//    P_POLYNOMIAL_ZEROS: zeros of Legendre function P(n,x).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NT, the order of the rule.
//
//    Output, double P_POLYNOMIAL_ZEROS[NT], the zeros.
//
{
  double *bj;
  int i;
  double *t;
  double *wts;
  
  t = new double[nt];

  for ( i = 0; i < nt; i++ )
  {
    t[i] = 0.0;
  }

  bj = new double[nt];
  for ( i = 0; i < nt; i++ )
  {
    bj[i] = ( double ) ( ( i + 1 ) * ( i + 1 ) ) 
          / ( double ) ( 4 *  ( i + 1 ) * ( i + 1 ) - 1 );
    bj[i] = sqrt ( bj[i] );
  }

  wts = new double[nt];
  wts[0] = sqrt ( 2.0 );
  for ( i = 1; i < nt; i++ )
  {
    wts[i] = 0.0;
  }

  imtqlx ( nt, t, bj, wts );

  delete [] bj;
  delete [] wts;

  return t;
}
//****************************************************************************80

double *p_power_product ( int p, int e )

//****************************************************************************80
//
//  Purpose:
//
//    P_POWER_PRODUCT: power products for Legendre polynomial P(n,x).
//
//  Discussion:
//
//    Let P(n,x) represent the Legendre polynomial of degree n.  
//
//    For polynomial chaos applications, it is of interest to know the
//    value of the integrals of products of X with every possible pair
//    of basis functions.  That is, we'd like to form
//
//      Tij = Integral ( -1.0 <= X <= +1.0 ) X^E * P(I,x) * P(J,x) dx
//
//    We will estimate these integrals using Gauss-Legendre quadrature.
//
//    When E is 0, the computed table should be the identity matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P, the maximum degree of the polyonomial 
//    factors.  0 <= P.
//
//    Input, int E, the exponent of X in the integrand.
//    0 <= E.
//
//    Output, double P_POWER_PRODUCT[(P+1)*(P+1)], the table of integrals.  
//
{
  double *h_table;
  int i;
  int j;
  int k;
  int order;
  double *table;
  double *w_table;
  double x;
  double *x_table;

  table = new double[(p+1)*(p+1)];

  for ( j = 0; j <= p; j++ )
  {
    for ( i = 0; i <= p; i++ )
    {
      table[i+j*(p+1)] = 0.0;
    }
  }

  order = p + 1 + ( ( e + 1 ) / 2 );

  x_table = new double[order];
  w_table = new double[order];

  p_quadrature_rule ( order, x_table, w_table );

  for ( k = 0; k < order; k++ )
  {
    x = x_table[k];
    h_table = p_polynomial_value ( 1, p, x_table + k );
//
//  The following formula is an outer product in H_TABLE.
//
    if ( e == 0 )
    {
      for ( i = 0; i <= p; i++ )
      {
        for ( j = 0; j <= p; j++ )
        {
          table[i+j*(p+1)] = table[i+j*(p+1)] + w_table[k] * h_table[i] * h_table[j];
        }
      }
    }
    else
    {
      for ( i = 0; i <= p; i++ )
      {
        for ( j = 0; j <= p; j++ )
        {
          table[i+j*(p+1)] = table[i+j*(p+1)] 
            + w_table[k] * pow ( x, e ) * h_table[i] * h_table[j];
        }
      }
    }
    delete [] h_table;
  }

  delete [] w_table;
  delete [] x_table;

  return table;
}
//****************************************************************************80

void p_quadrature_rule ( int nt, double t[], double wts[] )

//****************************************************************************80
//
//  Purpose:
//
//    P_QUADRATURE_RULE: quadrature for Legendre function P(n,x).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NT, the order of the rule.
//
//    Output, double T[NT], WTS[NT], the points and weights
//    of the rule.
//
{
  double *bj;
  int i;
  
  for ( i = 0; i < nt; i++ )
  {
    t[i] = 0.0;
  }

  bj = new double[nt];
  for ( i = 0; i < nt; i++ )
  {
    bj[i] = ( double ) ( ( i + 1 ) * ( i + 1 ) ) 
          / ( double ) ( 4 *  ( i + 1 ) * ( i + 1 ) - 1 );
    bj[i] = sqrt ( bj[i] );
  }

  wts[0] = sqrt ( 2.0 );
  for ( i = 1; i < nt; i++ )
  {
    wts[i] = 0.0;
  }

  imtqlx ( nt, t, bj, wts );

  for ( i = 0; i < nt; i++ )
  {
    wts[i] = pow ( wts[i], 2 );
  }
  delete [] bj;

  return;
}
//****************************************************************************80

double *pm_polynomial_value ( int mm, int n, int m, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    PM_POLYNOMIAL_VALUE evaluates the Legendre polynomials Pm(n,m,x).
//
//  Differential equation:
//
//    (1-X*X) * Y'' - 2 * X * Y + ( N (N+1) - (M*M/(1-X*X)) * Y = 0
//
//  First terms:
//
//    M = 0  ( = Legendre polynomials of first kind P(N,X) )
//
//    Pm(0,0,x) =    1
//    Pm(1,0,x) =    1 X
//    Pm(2,0,x) = (  3 X^2 -   1)/2
//    Pm(3,0,x) = (  5 X^3 -   3 X)/2
//    Pm(4,0,x) = ( 35 X^4 -  30 X^2 +   3)/8
//    Pm(5,0,x) = ( 63 X^5 -  70 X^3 +  15 X)/8
//    Pm(6,0,x) = (231 X^6 - 315 X^4 + 105 X^2 -  5)/16
//    Pm(7,0,x) = (429 X^7 - 693 X^5 + 315 X^3 - 35 X)/16
//
//    M = 1
//
//    Pm(0,1,x) =   0
//    Pm(1,1,x) =   1 * SQRT(1-X^2)
//    Pm(2,1,x) =   3 * SQRT(1-X^2) * X
//    Pm(3,1,x) = 1.5 * SQRT(1-X^2) * (5*X^2-1)
//    Pm(4,1,x) = 2.5 * SQRT(1-X^2) * (7*X^3-3*X)
//
//    M = 2
//
//    Pm(0,2,x) =   0
//    Pm(1,2,x) =   0
//    Pm(2,2,x) =   3 * (1-X^2)
//    Pm(3,2,x) =  15 * (1-X^2) * X
//    Pm(4,2,x) = 7.5 * (1-X^2) * (7*X^2-1)
//
//    M = 3
//
//    Pm(0,3,x) =   0
//    Pm(1,3,x) =   0
//    Pm(2,3,x) =   0
//    Pm(3,3,x) =  15 * (1-X^2)^1.5
//    Pm(4,3,x) = 105 * (1-X^2)^1.5 * X
//
//    M = 4
//
//    Pm(0,4,x) =   0
//    Pm(1,4,x) =   0
//    Pm(2,4,x) =   0
//    Pm(3,4,x) =   0
//    Pm(4,4,x) = 105 * (1-X^2)^2
//
//  Recursion:
//
//    if N < M:
//      Pm(N,M,x) = 0
//    if N = M:
//      Pm(N,M,x) = (2*M-1)!! * (1-X*X)^(M/2) where N!! means the product of
//      all the odd integers less than or equal to N.
//    if N = M+1:
//      Pm(N,M,x) = X*(2*M+1)*Pm(M,M,x)
//    if M+1 < N:
//      Pm(N,M,x) = ( X*(2*N-1)*Pm(N-1,M,x) - (N+M-1)*Pm(N-2,M,x) )/(N-M)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input, int MM, the number of evaluation points.
//
//    Input, int N, the maximum first index of the Legendre
//    function, which must be at least 0.
//
//    Input, int M, the second index of the Legendre function,
//    which must be at least 0, and no greater than N.
//
//    Input, double X[MM], the point at which the function is to be
//    evaluated.
//
//    Output, double PM_POLYNOMIAL_VALUE[MM*(N+1)], the function values.
//
{
  double fact;
  int i;
  int j;
  int k;
  double *v;

  v = new double[mm*(n+1)];

  for ( j = 0; j < n + 1; j++ )
  {
    for ( i = 0; i < mm; i++ )
    {
      v[i+j*mm] = 0.0;
    }
  }

//
//  J = M is the first nonzero function.
//
  if ( m <= n )
  {
    for ( i = 0; i < mm; i++ )
    {
      v[i+m*mm] = 1.0;
    }

    fact = 1.0;
    for ( k = 0; k < m; k++ )
    {
      for ( i = 0; i < mm; i++ )
      {
        v[i+m*mm] = - v[i+m*mm] * fact * sqrt ( 1.0 - x[i] * x[i] );
      }
      fact = fact + 2.0;
    }
  }

//
//  J = M + 1 is the second nonzero function.
//
  if ( m + 1 <= n )
  {
    for ( i = 0; i < mm; i++ )
    {
      v[i+(m+1)*mm] = x[i] * ( double ) ( 2 * m + 1 ) * v[i+m*mm];
    }
  }

//
//  Now we use a three term recurrence.
//
  for ( j = m + 2; j <= n; j++ )
  {
    for ( i = 0; i < mm; i++ )
    {
      v[i+j*mm] = ( ( double ) ( 2 * j     - 1 ) * x[i] * v[i+(j-1)*mm]
                  + ( double ) (   - j - m + 1 ) *        v[i+(j-2)*mm] )
                  / ( double ) (     j - m     );
    }
  }

  return v;
}

//****************************************************************************80

void pm_polynomial_values ( int &n_data, int &n, int &m, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    PM_POLYNOMIAL_VALUES: selected values of Legendre polynomials Pm(n,m,x).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, int &M, double &X,
//    the arguments of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
      0.0000000000000000E+00,
     -0.5000000000000000E+00,
      0.0000000000000000E+00,
      0.3750000000000000E+00,
      0.0000000000000000E+00,
     -0.8660254037844386E+00,
     -0.1299038105676658E+01,
     -0.3247595264191645E+00,
      0.1353164693413185E+01,
     -0.2800000000000000E+00,
      0.1175755076535925E+01,
      0.2880000000000000E+01,
     -0.1410906091843111E+02,
     -0.3955078125000000E+01,
     -0.9997558593750000E+01,
      0.8265311444100484E+02,
      0.2024442836815152E+02,
     -0.4237997531890869E+03,
      0.1638320624828339E+04,
     -0.2025687389227225E+05  };

  static int m_vec[N_MAX] = {
    0, 0, 0, 0,
    0, 1, 1, 1,
    1, 0, 1, 2,
    3, 2, 2, 3,
    3, 4, 4, 5 };

  static int n_vec[N_MAX] = {
    1,  2,  3,  4,
    5,  1,  2,  3,
    4,  3,  3,  3,
    3,  4,  5,  6,
    7,  8,  9, 10 };

  static double x_vec[N_MAX] = {
     0.00E+00,
     0.00E+00,
     0.00E+00,
     0.00E+00,
     0.00E+00,
     0.50E+00,
     0.50E+00,
     0.50E+00,
     0.50E+00,
     0.20E+00,
     0.20E+00,
     0.20E+00,
     0.20E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    m = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    m = m_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double *pmn_polynomial_value ( int mm, int n, int m, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    PMN_POLYNOMIAL_VALUE: normalized Legendre polynomial Pmn(n,m,x).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input, int MM, the number of evaluation points.
//
//    Input, int N, the maximum first index of the Legendre
//    function, which must be at least 0.
//
//    Input, int M, the second index of the Legendre function,
//    which must be at least 0, and no greater than N.
//
//    Input, double X[MM], the evaluation points.
//
//    Output, double PMN_POLYNOMIAL_VALUE[MM*(N+1)], the function values.
//
{
  double factor;
  int i;
  int j;
  double *v;

  v = pm_polynomial_value ( mm, n, m, x );
//
//  Normalization.
//
  for ( j = m; j <= n; j++ )
  {
    factor = sqrt ( ( ( double ) ( 2 * j + 1 ) * r8_factorial ( j - m ) ) 
      / ( 2.0 * r8_factorial ( j + m ) ) );
    for ( i = 0; i < mm; i++ )
    {
      v[i+j*mm] = v[i+j*mm] * factor;
    }
  }

  return v;
}
//****************************************************************************80

void pmn_polynomial_values ( int &n_data, int &n, int &m, double &x, 
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    PMN_POLYNOMIAL_VALUES: selected values of normalized Legendre polynomial Pmn(n,m,x).
//
//  Discussion:
//
//    In Mathematica, the unnormalized function can be evaluated by:
//
//      LegendreP [ n, m, x ]
//
//    The function is normalized by dividing by
//
//      sqrt ( 2 * ( n + m )! / ( 2 * n + 1 ) / ( n - m )! )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, int &M, double &X,
//    the arguments of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 21

  static double fx_vec[N_MAX] = {
    0.7071067811865475E+00, 
    0.6123724356957945E+00, 
   -0.7500000000000000E+00, 
   -0.1976423537605237E+00, 
   -0.8385254915624211E+00, 
    0.7261843774138907E+00, 
   -0.8184875533567997E+00, 
   -0.1753901900050285E+00, 
    0.9606516343087123E+00, 
   -0.6792832849776299E+00, 
   -0.6131941618102092E+00, 
    0.6418623720763665E+00, 
    0.4716705890038619E+00, 
   -0.1018924927466445E+01, 
    0.6239615396237876E+00, 
    0.2107022704608181E+00, 
    0.8256314721961969E+00, 
   -0.3982651281554632E+00, 
   -0.7040399320721435E+00, 
    0.1034723155272289E+01, 
   -0.5667412129155530E+00 };

  static int m_vec[N_MAX] = {
    0, 0, 1, 0,
    1, 2, 0, 1,
    2, 3, 0, 1,
    2, 3, 4, 0,
    1, 2, 3, 4,
    5 };

  static int n_vec[N_MAX] = {
    0,  1,  1,  2,
    2,  2,  3,  3,
    3,  3,  4,  4,
    4,  4,  4,  5,
    5,  5,  5,  5,
    5 };

  static double x_vec[N_MAX] = {
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    m = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    m = m_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double *pmns_polynomial_value ( int mm, int n, int m, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    PMNS_POLYNOMIAL_VALUE: sphere-normalized Legendre polynomial Pmn(n,m,x).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input, int MM, the number of evaluation points.
//
//    Input, int N, the maximum first index of the Legendre
//    function, which must be at least 0.
//
//    Input, int M, the second index of the Legendre function,
//    which must be at least 0, and no greater than N.
//
//    Input, double X[MM], the evaluation points.
//
//    Output, double PMNS_POLYNOMIAL_VALUE[MM*(N+1)], the function values.
//
{
  double factor;
  int i;
  int j;
  const double pi = 3.141592653589793;
  double *v;

  v = pm_polynomial_value ( mm, n, m, x );
//
//  Normalization.
//
  for ( j = m; j <= n; j++ )
  {
    factor = sqrt ( ( ( double ) ( 2 * j + 1 ) * r8_factorial ( j - m ) ) 
      / ( 4.0 * pi * r8_factorial ( j + m ) ) );
    for ( i = 0; i < mm; i++ )
    {
      v[i+j*mm] = v[i+j*mm] * factor;
    }
  }

  return v;
}
//****************************************************************************80

void pmns_polynomial_values ( int &n_data, int &n, int &m, double &x, 
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    PMNS_POLYNOMIAL_VALUES: selected values of sphere-normalized Legendre polynomial Pmns(n,m,x).
//
//  Discussion:
//
//    In Mathematica, the unnormalized function can be evaluated by:
//
//      LegendreP [ n, m, x ]
//
//    The function is normalized for the sphere by dividing by
//
//      sqrt ( 4 * pi * ( n + m )! / ( 2 * n + 1 ) / ( n - m )! )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, int &M, double &X,
//    the arguments of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 21

  static double fx_vec[N_MAX] = {
     0.2820947917738781,
     0.2443012559514600,
    -0.2992067103010745,
    -0.07884789131313000,
    -0.3345232717786446,
     0.2897056515173922,
    -0.3265292910163510,
    -0.06997056236064664,
     0.3832445536624809,
    -0.2709948227475519,
    -0.2446290772414100,
     0.2560660384200185,
     0.1881693403754876,
    -0.4064922341213279,
     0.2489246395003027,
     0.08405804426339821,
     0.3293793022891428,
    -0.1588847984307093,
    -0.2808712959945307,
     0.4127948151484925,
    -0.2260970318780046 };

  static int m_vec[N_MAX] = {
    0, 0, 1, 0,
    1, 2, 0, 1,
    2, 3, 0, 1,
    2, 3, 4, 0,
    1, 2, 3, 4,
    5 };

  static int n_vec[N_MAX] = {
    0,  1,  1,  2,
    2,  2,  3,  3,
    3,  3,  4,  4,
    4,  4,  4,  5,
    5,  5,  5,  5,
    5 };

  static double x_vec[N_MAX] = {
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    m = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    m = m_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double *pn_pair_product ( int p )

//****************************************************************************80
//
//  Purpose:
//
//    PN_PAIR_PRODUCT: pair products for normalized Legendre polynomial Pn(n,x).
//
//  Discussion:
//
//    Let P(n,x) represent the Legendre polynomial of degree n.  
//
//    To check orthonormality, we compute
//
//      Tij = Integral ( -1.0 <= X <= +1.0 ) Pn(i,x) * Pn(j,x) dx
//
//    We will estimate these integrals using Gauss-Legendre quadrature.
//
//    The computed table should be the identity matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P, the maximum degree of the polyonomial 
//    factors.  0 <= P.
//
//    Input, int E, the exponent of X in the integrand.
//    0 <= E.
//
//    Output, double PN_PAIR_PRODUCT[(P+1)*(P+1)], the table of integrals.  
//
{
  double *h_table;
  int i;
  int j;
  int k;
  int order;
  double *table;
  double *w_table;
  double x;
  double *x_table;

  table = new double[(p+1)*(p+1)];

  for ( j = 0; j <= p; j++ )
  {
    for ( i = 0; i <= p; i++ )
    {
      table[i+j*(p+1)] = 0.0;
    }
  }

  order = p + 1;

  x_table = new double[order];
  w_table = new double[order];

  p_quadrature_rule ( order, x_table, w_table );

  for ( k = 0; k < order; k++ )
  {
    x = x_table[k];
    h_table = pn_polynomial_value ( 1, p, x_table + k );

    for ( i = 0; i <= p; i++ )
    {
      for ( j = 0; j <= p; j++ )
      {
        table[i+j*(p+1)] = table[i+j*(p+1)] 
          + w_table[k] * h_table[i] * h_table[j];
      }
    }
    delete [] h_table;
  }

  delete [] w_table;
  delete [] x_table;

  return table;
}
//****************************************************************************80

double *pn_polynomial_coefficients ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    PN_POLYNOMIAL_COEFFICIENTS: coefficients of normalized Legendre Pn(n,x).
//
//  Discussion:
//
//    Pn(n,x) = P(n,x) * sqrt ( (2n+1)/2 )
//
//          1       x       x^2     x^3     x^4      x^5    x^6     x^7
//
//    0   0.707
//    1   0.000   1.224
//    2  -0.790   0.000   2.371
//    3   0.000  -2.806   0.000   4.677
//    4   0.795   0.000  -7.954   0.000   9.280
//    5   0.000   4.397   0.000 -20.520   0.000   18.468
//    6  -0.796   0.000  16.731   0.000 -50.193    0.000  36.808
//    7   0.000  -5.990   0.000  53.916   0.000 -118.616   0.000  73.429 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996.
//
//  Parameters:
//
//    Input, int N, the highest order polynomial to evaluate.
//    Note that polynomials 0 through N will be evaluated.
//
//    Output, double PN_POLYNOMIAL_COEFFICIENTS[(N+1)*(N+1)], the coefficients of 
//    the normalized Legendre polynomials of degree 0 through N.
//
{
  double *c;
  int i;
  int j;
  double t;

  if ( n < 0 )
  {
    return NULL;
  }
//
//  Compute P(i,x) coefficients.
//
  c = new double[(n+1)*(n+1)];

  for ( i = 0; i <= n; i++ )
  {
    for ( j = 0; j <= n; j++ )
    {
      c[i+j*(n+1)] = 0.0;
    }
  }

  c[0+0*(n+1)] = 1.0;

  if ( 0 < n )
  {
    c[1+1*(n+1)] = 1.0;
  }

  for ( i = 2; i <= n; i++ )
  {
    for ( j = 0; j <= i-2; j++ )
    {
      c[i+j*(n+1)] =          
          ( double ) (   - i + 1 ) * c[i-2+j*(n+1)] / ( double ) i;
    }
    for ( j = 1; j <= i; j++ )
    {
      c[i+j*(n+1)] = c[i+j*(n+1)] 
        + ( double ) ( i + i - 1 ) * c[i-1+(j-1)*(n+1)] / ( double ) i;
    }
  }
//
//  Normalize them.
//
  for ( i = 0; i <= n; i++ )
  {
    t = sqrt ( ( double ) ( 2 * i + 1 ) / 2.0 );
    for ( j = 0; j <= i; j++ )
    {
      c[i+j*(n+1)] = c[i+j*(n+1)] * t;
    }
  }

  return c;
}
//****************************************************************************80

double *pn_polynomial_value ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    PN_POLYNOMIAL_VALUE evaluates the normalized Legendre polynomials Pn(n,x).
//
//  Discussion:
//
//    The normalized Legendre polynomials are orthonormal under the inner product 
//    defined as integration from -1 to 1:
//
//      Integral ( -1 <= x <= +1 ) Pn(i,x) * Pn(j,x) dx = delta(i,j)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996.
//
//  Parameters:
//
//    Input, int M, the number of evaluation points.
//
//    Input, int N, the highest order polynomial to evaluate.
//    Note that polynomials 0 through N will be evaluated.
//
//    Input, double X[M], the evaluation points.
//
//    Output, double PN_POLYNOMIAL_VALUE[M*(N+1)], the values of the Legendre
//    polynomials of order 0 through N.
//
{
  int i;
  int j;
  double norm;
  double *v;

  v = p_polynomial_value ( m, n, x );

  for ( j = 0; j <= n; j++ )
  {
    norm = sqrt ( 2 / ( double ) ( 2 * j + 1 ) );
    for ( i = 0; i < m; i++ )
    {
      v[i+j*m] = v[i+j*m] / norm;
    }
  }
  return v;
}
//****************************************************************************80

void pn_polynomial_values ( int &n_data, int &n, double &x, 
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    PN_POLYNOMIAL_VALUES: selected values of the normalized Legendre polynomials.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 March 2016
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the function.
//
//    Output, double &X, the point where the function is evaluated.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 22

  static double fx_vec[N_MAX] = {
    0.7071067811865475, 
    0.3061862178478972, 
   -0.642337649721702, 
   -0.6284815141846855, 
    0.3345637065282053, 
    0.7967179601799685, 
    0.06189376866246124, 
   -0.766588850921089, 
   -0.4444760242953344, 
    0.5450094674858101, 
    0.7167706229835538, 
    0.0000000000000000, 
   -0.2759472322745781, 
   -0.5238320341483518, 
   -0.7155919752205163, 
   -0.823164625090267, 
   -0.8184875533567997, 
   -0.6734983296193094, 
   -0.360134523476992, 
    0.1496662954709581, 
    0.8839665576253438, 
    1.870828693386971 };

  static int n_vec[N_MAX] = {
     0,  1,  2,
     3,  4,  5,
     6,  7,  8,
     9, 10,  3,
     3,  3,  3,
     3,  3,  3,
     3,  3,  3,
     3 };

  static double x_vec[N_MAX] = {
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.00E+00,
     0.10E+00,
     0.20E+00,
     0.30E+00,
     0.40E+00,
     0.50E+00,
     0.60E+00,
     0.70E+00,
     0.80E+00,
     0.90E+00,
     1.00E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}



// //****************************************************************************80

// double r8_epsilon ( )

// //****************************************************************************80
// //
// //  Purpose:
// //
// //    R8_EPSILON returns the R8 roundoff unit.
// //
// //  Discussion:
// //
// //    The roundoff unit is a number R which is a power of 2 with the
// //    property that, to the precision of the computer's arithmetic,
// //      1 < 1 + R
// //    but
// //      1 = ( 1 + R / 2 )
// //
// //  Licensing:
// //
// //    This code is distributed under the GNU LGPL license.
// //
// //  Modified:
// //
// //    01 September 2012
// //
// //  Author:
// //
// //    John Burkardt
// //
// //  Parameters:
// //
// //    Output, double R8_EPSILON, the R8 round-off unit.
// //
// {
//   const double value = 2.220446049250313E-016;

//   return value;
// }
// //****************************************************************************80

// double r8_factorial ( int n )

// //****************************************************************************80
// //
// //  Purpose:
// //
// //    R8_FACTORIAL computes the factorial of N.
// //
// //  Discussion:
// //
// //    factorial ( N ) = product ( 1 <= I <= N ) I
// //
// //  Licensing:
// //
// //    This code is distributed under the GNU LGPL license.
// //
// //  Modified:
// //
// //    16 January 1999
// //
// //  Author:
// //
// //    John Burkardt
// //
// //  Parameters:
// //
// //    Input, int N, the argument of the factorial function.
// //    If N is less than 1, the function value is returned as 1.
// //
// //    Output, double R8_FACTORIAL, the factorial of N.
// //
// {
//   int i;
//   double value;

//   value = 1.0;

//   for ( i = 1; i <= n; i++ )
//   {
//     value = value * ( double ) ( i );
//   }

//   return value;
// }
// //****************************************************************************80

// double r8_sign ( double x )

// //****************************************************************************80
// //
// //  Purpose:
// //
// //    R8_SIGN returns the sign of an R8.
// //
// //  Licensing:
// //
// //    This code is distributed under the GNU LGPL license.
// //
// //  Modified:
// //
// //    18 October 2004
// //
// //  Author:
// //
// //    John Burkardt
// //
// //  Parameters:
// //
// //    Input, double X, the number whose sign is desired.
// //
// //    Output, double R8_SIGN, the sign of X.
// //
// {
//   double value;

//   if ( x < 0.0 )
//   {
//     value = -1.0;
//   }
//   else
//   {
//     value = 1.0;
//   }
//   return value;
// }
// //****************************************************************************80

// void r8mat_print ( int m, int n, double a[], string title )

// //****************************************************************************80
// //
// //  Purpose:
// //
// //    R8MAT_PRINT prints an R8MAT.
// //
// //  Discussion:
// //
// //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
// //    in column-major order.
// //
// //    Entry A(I,J) is stored as A[I+J*M]
// //
// //  Licensing:
// //
// //    This code is distributed under the GNU LGPL license.
// //
// //  Modified:
// //
// //    10 September 2009
// //
// //  Author:
// //
// //    John Burkardt
// //
// //  Parameters:
// //
// //    Input, int M, the number of rows in A.
// //
// //    Input, int N, the number of columns in A.
// //
// //    Input, double A[M*N], the M by N matrix.
// //
// //    Input, string TITLE, a title.
// //
// {
//   r8mat_print_some ( m, n, a, 1, 1, m, n, title );

//   return;
// }
// //****************************************************************************80

// void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
//   int jhi, string title )

// //****************************************************************************80
// //
// //  Purpose:
// //
// //    R8MAT_PRINT_SOME prints some of an R8MAT.
// //
// //  Discussion:
// //
// //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
// //    in column-major order.
// //
// //  Licensing:
// //
// //    This code is distributed under the GNU LGPL license.
// //
// //  Modified:
// //
// //    20 August 2010
// //
// //  Author:
// //
// //    John Burkardt
// //
// //  Parameters:
// //
// //    Input, int M, the number of rows of the matrix.
// //    M must be positive.
// //
// //    Input, int N, the number of columns of the matrix.
// //    N must be positive.
// //
// //    Input, double A[M*N], the matrix.
// //
// //    Input, int ILO, JLO, IHI, JHI, designate the first row and
// //    column, and the last row and column to be printed.
// //
// //    Input, string TITLE, a title.
// //
// {
// # define INCX 5

//   int i;
//   int i2hi;
//   int i2lo;
//   int j;
//   int j2hi;
//   int j2lo;

//   cout << "\n";
//   cout << title << "\n";

//   if ( m <= 0 || n <= 0 )
//   {
//     cout << "\n";
//     cout << "  (None)\n";
//     return;
//   }
// //
// //  Print the columns of the matrix, in strips of 5.
// //
//   for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
//   {
//     j2hi = j2lo + INCX - 1;
//     j2hi = i4_min ( j2hi, n );
//     j2hi = i4_min ( j2hi, jhi );

//     cout << "\n";
// //
// //  For each column J in the current range...
// //
// //  Write the header.
// //
//     cout << "  Col:    ";
//     for ( j = j2lo; j <= j2hi; j++ )
//     {
//       cout << setw(7) << j - 1 << "       ";
//     }
//     cout << "\n";
//     cout << "  Row\n";
//     cout << "\n";
// //
// //  Determine the range of the rows in this strip.
// //
//     i2lo = i4_max ( ilo, 1 );
//     i2hi = i4_min ( ihi, m );

//     for ( i = i2lo; i <= i2hi; i++ )
//     {
// //
// //  Print out (up to) 5 entries in row I, that lie in the current strip.
// //
//       cout << setw(5) << i - 1 << ": ";
//       for ( j = j2lo; j <= j2hi; j++ )
//       {
//         cout << setw(12) << a[i-1+(j-1)*m] << "  ";
//       }
//       cout << "\n";
//     }
//   }

//   return;
// # undef INCX
// }
// //****************************************************************************80

// double r8vec_dot_product ( int n, double a1[], double a2[] )

// //****************************************************************************80
// //
// //  Purpose:
// //
// //    R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.
// //
// //  Discussion:
// //
// //    An R8VEC is a vector of R8's.
// //
// //  Licensing:
// //
// //    This code is distributed under the GNU LGPL license.
// //
// //  Modified:
// //
// //    03 July 2005
// //
// //  Author:
// //
// //    John Burkardt
// //
// //  Parameters:
// //
// //    Input, int N, the number of entries in the vectors.
// //
// //    Input, double A1[N], A2[N], the two vectors to be considered.
// //
// //    Output, double R8VEC_DOT_PRODUCT, the dot product of the vectors.
// //
// {
//   int i;
//   double value;

//   value = 0.0;
//   for ( i = 0; i < n; i++ )
//   {
//     value = value + a1[i] * a2[i];
//   }
//   return value;
// }
// //****************************************************************************80

// double *r8vec_linspace_new ( int n, double a_first, double a_last )

// //****************************************************************************80
// //
// //  Purpose:
// //
// //    R8VEC_LINSPACE_NEW creates a vector of linearly spaced values.
// //
// //  Discussion:
// //
// //    An R8VEC is a vector of R8's.
// //
// //    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
// //
// //    In other words, the interval is divided into N-1 even subintervals,
// //    and the endpoints of intervals are used as the points.
// //
// //  Licensing:
// //
// //    This code is distributed under the GNU LGPL license.
// //
// //  Modified:
// //
// //    29 March 2011
// //
// //  Author:
// //
// //    John Burkardt
// //
// //  Parameters:
// //
// //    Input, int N, the number of entries in the vector.
// //
// //    Input, double A_FIRST, A_LAST, the first and last entries.
// //
// //    Output, double R8VEC_LINSPACE_NEW[N], a vector of linearly spaced data.
// //
// {
//   double *a;
//   int i;

//   a = new double[n];

//   if ( n == 1 )
//   {
//     a[0] = ( a_first + a_last ) / 2.0;
//   }
//   else
//   {
//     for ( i = 0; i < n; i++ )
//     {
//       a[i] = ( ( double ) ( n - 1 - i ) * a_first 
//              + ( double ) (         i ) * a_last ) 
//              / ( double ) ( n - 1     );
//     }
//   }
//   return a;
// }
// //****************************************************************************80

// void r8vec_print ( int n, double a[], string title )

// //****************************************************************************80
// //
// //  Purpose:
// //
// //    R8VEC_PRINT prints an R8VEC.
// //
// //  Discussion:
// //
// //    An R8VEC is a vector of R8's.
// //
// //  Licensing:
// //
// //    This code is distributed under the GNU LGPL license.
// //
// //  Modified:
// //
// //    16 August 2004
// //
// //  Author:
// //
// //    John Burkardt
// //
// //  Parameters:
// //
// //    Input, int N, the number of components of the vector.
// //
// //    Input, double A[N], the vector to be printed.
// //
// //    Input, string TITLE, a title.
// //
// {
//   int i;

//   cout << "\n";
//   cout << title << "\n";
//   cout << "\n";
//   for ( i = 0; i < n; i++ )
//   {
//     cout << "  " << setw(8)  << i
//          << ": " << setw(14) << a[i]  << "\n";
//   }

//   return;
// }
// //****************************************************************************80

// void r8vec2_print ( int n, double a1[], double a2[], string title )

// //****************************************************************************80
// //
// //  Purpose:
// //
// //    R8VEC2_PRINT prints an R8VEC2.
// //
// //  Discussion:
// //
// //    An R8VEC2 is a dataset consisting of N pairs of real values, stored
// //    as two separate vectors A1 and A2.
// //
// //  Licensing:
// //
// //    This code is distributed under the GNU LGPL license.
// //
// //  Modified:
// //
// //    19 November 2002
// //
// //  Author:
// //
// //    John Burkardt
// //
// //  Parameters:
// //
// //    Input, int N, the number of components of the vector.
// //
// //    Input, double A1[N], double A2[N], the vectors to be printed.
// //
// //    Input, string TITLE, a title.
// //
// {
//   int i;

//   cout << "\n";
//   cout << title << "\n";
//   cout << "\n";
//   for ( i = 0; i <= n - 1; i++ )
//   {
//     cout << setw(6)  << i
//          << ": " << setw(14) << a1[i]
//          << "  " << setw(14) << a2[i] << "\n";
//   }

//   return;
// }
// //****************************************************************************80

// void timestamp ( )

// //****************************************************************************80
// //
// //  Purpose:
// //
// //    TIMESTAMP prints the current YMDHMS date as a time stamp.
// //
// //  Example:
// //
// //    31 May 2001 09:45:54 AM
// //
// //  Licensing:
// //
// //    This code is distributed under the GNU LGPL license.
// //
// //  Modified:
// //
// //    08 July 2009
// //
// //  Author:
// //
// //    John Burkardt
// //
// //  Parameters:
// //
// //    None
// //
// {
// # define TIME_SIZE 40

//   static char time_buffer[TIME_SIZE];
//   const struct std::tm *tm_ptr;
//   size_t len;
//   std::time_t now;

//   now = std::time ( NULL );
//   tm_ptr = std::localtime ( &now );

//   len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

//   std::cout << time_buffer << "\n";

//   return;
// # undef TIME_SIZE
// }
