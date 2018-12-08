
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

#include "statistics.h"
#include "helper/helper.h"
#include "matrix.h"
#include "dcdflib.h"
#include "ipmpar.h"

#include <iostream>
#include <cmath>
#include <algorithm>


#ifndef M_2PI
#define M_2PI 6.283185307179586476925286766559/* 2*pi */
#endif

#ifndef LOG_BOOL
#define LOG_BOOL false
#endif

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI 0.918938533204672741780329736406/* log(sqrt(2*pi)) */
#endif

#ifndef M_LN_SQRT_PId2
#define M_LN_SQRT_PId2 0.225791352644727432363097614947/* log(sqrt(pi/2)) */
#endif

const double EPS = 1.0e-6;

double Statistics::sum( const Data::Vector<double> & a )
{
  double s = 0;
  for (int i=0; i<a.size(); i++) s += a[i];
  return s;
}

double Statistics::sum_squares( const Data::Vector<double> & a )
{
  double s = 0;
  for (int i=0; i<a.size(); i++) s += a[i] * a[i];
  return s;
}

Data::Vector<double> Statistics::row_sums( const Data::Matrix<double> & a )
{
  Data::Vector<double> r( a.dim1() );
  for (int i=0;i<a.dim1();i++)
    for (int j=0;j<a.dim2();j++)
      r[i] += a(i,j);   // avoid making temporary vector from row
  return r;
}

Data::Vector<double> Statistics::col_sums( const Data::Matrix<double> & a)
{
  Data::Vector<double> r( a.dim2() );
  for (int i=0;i<a.dim2();i++)
    r[i] = sum( a.col(i) ); // use fast column reference  
  return r;
}


Data::Vector<double> Statistics::mean_center_cols( const Data::Matrix<double> & d )
{
  // and return the vector of original means

  Data::Vector<double> means = mean( d );
  
  const int nr = d.dim1();  
  const int nc = d.dim2();
  
  for (int c=0;c<nc;c++)
    d.col(c).inplace_add( -means[c] );
  
  return means;
}

Data::Vector<double> Statistics::mean( const Data::Matrix<double> & d )
{
  Data::Vector<double> m( d.dim2() );

  for (int j=0; j<d.dim2(); j++)
    {
      for (int i=0; i<d.dim1(); i++)
	m[j] += d(i,j);
      m[j] /= d.dim1();
    }
  return m;
}

Data::Vector<double> Statistics::variance( const Data::Matrix<double> & d )
{
  return variance( d , mean(d) );
}

Data::Vector<double> Statistics::variance( const Data::Matrix<double> & d , const Data::Vector<double> & u )
{
  Data::Vector<double> v( d.dim2() ); 
  Data::Matrix<double> s = covariance_matrix( d , u , d , u );
  for (int i=0; i<d.dim2(); i++) v(i) = s(i,i);
  return v;
}

Data::Matrix<double> Statistics::covariance_matrix( const Data::Matrix<double> & d )
{
  return covariance_matrix( d , mean(d) );
}

Data::Matrix<double> Statistics::covariance_matrix( const Data::Matrix<double> & d , const Data::Vector<double> & u )
{
  return covariance_matrix( d , u , d , u );  
}  



std::vector<double> Statistics::as_vector( const Data::Vector<double> & d )
{
  std::vector<double> v( d.size() ) ;
  for (int i=0; i<d.size(); i++) v[i] = d[i];
  return v;
}

Data::Matrix<double> Statistics::covariance_matrix( const Data::Matrix<double> & x , 
						      const Data::Vector<double> & u ,
						      const Data::Matrix<double> & y , 
						      const Data::Vector<double> & v )
{
  
  // calculate Sxy e.g. lower quadrant of partitioned covariance matrix
  // note -- ignores symmetry -- easy speedup to add 
  if ( x.dim1() != y.dim1() ) Helper::halt("internal error, unequal row numbers in covariance_matrix()"); 
  const int n = x.dim1();
  Data::Matrix<double> s( x.dim2() , y.dim2() ) ;
  for (int i=0; i<x.dim2(); i++)
    for (int j=0; j<y.dim2(); j++)
      {
	const double mx = u[i];
	const double my = v[j];
	for (int k=0;k<n;k++)
	  s(i,j) += ( x(k,i) - mx ) * ( y(k,j) - my ); 
	s(i,j) /= n-1;
	//	s(j,i) = s(i,j);
      }
  return s;
}

Data::Matrix<double> Statistics::covariance_matrix( const Data::Matrix<double> & x , const Data::Matrix<double> & y )
{
  return covariance_matrix( x , mean(x) , y , mean(y) );
}

Data::Matrix<double> Statistics::transpose( const Data::Matrix<double> & d )
{
  const int row = d.dim1();
  const int col = d.dim2();
  Data::Matrix<double> r( col, row );
  for (int i = 0; i < row; i++)
    for (int j = 0; j < col; j++)
      r(j,i) = d(i,j);
  return r;
}

Data::Matrix<double> Statistics::inverse( const Data::Matrix<double> & u_orig, bool * flag )
{
  
  const double eps = 1e-24; 

  Data::Matrix<double> u = u_orig;
  
  if ( u.dim1() == 0 || u.dim1() != u.dim2() ) 
    Helper::halt("cannot inverted non-square matrix");
  int n = u.dim1();
  
  Data::Vector<double> w(n);  
  Data::Matrix<double> v(n,n);

  if ( flag ) 
    *flag = Statistics::svdcmp( u, w, v ); 
  else 
    Statistics::svdcmp( u, w, v ); 
  
  // Look for singular values

  double wmax = 0;
  for (int i=0; i<n; i++)
    wmax = w[i] > wmax ? w[i] : wmax;
  double wmin = wmax * eps;
  for (int i=0; i<n; i++)
    w[i] = w[i] < wmin ? 0 : 1/w[i];
  
  // u w t(v)
  // row U * 1/w
  
  // results matrix
  Data::Matrix<double> r(n,n);

  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      u(i,j) *= w[j];
  
  // [nxn].[t(v)] 
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      for (int k=0; k<n; k++)
	r(i,j) += u(i,k) * v(j,k);
    
  return r;
}


Data::Matrix<double> Statistics::matrix_sqrt( const Data::Matrix<double> & u_orig )
{
  

  Data::Matrix<double> u = u_orig;
  
  // Using SVD, square root is U . sqrt(D) . V_T
  //  msqrt <- function(m) { m <- svd(m); m$u %*% sqrt(diag(m$d)) %*% t(m$v) }  

  const double eps = 1e-12;
  
  int n = u.dim1();

  Data::Vector<double> d(n);
  Data::Matrix<double> v(n,n);

  Statistics::svdcmp(u,d,v);

  // Take square root of diagonal values                                                                                                                                                                                                                                                                          
  for (int i=0; i<n; i++)
    d[i] = sqrt(d[i]);

  // Multiplication to reconstruct original                                                                                                                                                                                                                                                                             
  Data::Matrix<double> r(n,n);
  Data::Matrix<double> r2(n,n);

  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      r(i,j) = u(i,j) * d[j];
  
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      for (int k=0; k<n; k++)
 	r2(i,j) += r(i,k) * v(j,k);
  
  return r2;
  
}

bool Statistics::svdcmp( Data::Matrix<double> & a, Data::Vector<double> & w , Data::Matrix<double> & v )
{
  
  bool flag;
  
  int i,its,j,jj,k,l,nm;
  
  double anorm,c,f,g,h,s,scale,x,y,z;
  
  double volatile temp;
  
  int m=a.dim1();
  if ( m == 0 ) Helper::halt("Internal problem in SVD function (no observations left?)");

  int n = a.dim2();

  std::vector<double> rv1( n );
  
  g=scale=anorm=0.0;
  for (i=0;i<n;i++) {
    l=i+2;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i < m) {
      for (k=i;k<m;k++) scale += fabs(a(k,i));
      if (scale != 0.0) {
	for (k=i;k<m;k++) {
	  a[k][i] /= scale;
	  s += a(k,i)*a(k,i);
	}
	f=a(i,i);
	g = -Statistics::SIGN(sqrt(s),f);
	h=f*g-s;
	a(i,i)=f-g;
	for (j=l-1;j<n;j++) {
	  for (s=0.0,k=i;k<m;k++) s += a(k,i)*a(k,j);
	  f=s/h;
	  for (k=i;k<m;k++) a(k,j) += f*a(k,i);
	}
	for (k=i;k<m;k++) a[k][i] *= scale;
      }
    }
    w[i]=scale *g;
    g=s=scale=0.0;
    if (i+1 <= m && i+1 != n) {
      for (k=l-1;k<n;k++) scale += fabs(a(i,k));
      if (scale != 0.0) {
	for (k=l-1;k<n;k++) {
	  a(i,k) /= scale;
	  s += a(i,k)*a(i,k);
	}
	f=a(i,l-1);
	g = -Statistics::SIGN(sqrt(s),f);
	h=f*g-s;
	a(i,l-1)=f-g;
	for (k=l-1;k<n;k++) rv1[k]=a(i,k)/h;
	for (j=l-1;j<m;j++) {
	  for (s=0.0,k=l-1;k<n;k++) s += a(j,k)*a(i,k);
	  for (k=l-1;k<n;k++) a(j,k) += s*rv1[k];
	}
	for (k=l-1;k<n;k++) a(i,k) *= scale;
      }
    }
    anorm=FNMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  for (i=n-1;i>=0;i--) {
    if (i < n-1) {
      if (g != 0.0) {
	for (j=l;j<n;j++)
	  v[j][i]=(a[i][j]/a[i][l])/g;
	for (j=l;j<n;j++) {
	  for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
	  for (k=l;k<n;k++) v[k][j] += s*v[k][i];
	}
      }
      for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=FNMIN(m,n)-1;i>=0;i--) {
    l=i+1;
    g=w[i];
    for (j=l;j<n;j++) a[i][j]=0.0;
    if (g != 0.0) {
      g=1.0/g;
      for (j=l;j<n;j++) {
	for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
	f=(s/a[i][i])*g;
	for (k=i;k<m;k++) a[k][j] += f*a[k][i];
      }
      for (j=i;j<m;j++) a[j][i] *= g;
    } else for (j=i;j<m;j++) a[j][i]=0.0;
    ++a[i][i];
  }
  for (k=n-1;k>=0;k--) {
    for (its=0;its<30;its++) {
      flag=true;
      for (l=k;l>=0;l--) {
	nm=l-1;
	temp=fabs(rv1[l])+anorm;
	if (temp == anorm) {
	  flag=false;
	  break;
	}
	temp=fabs(w[nm])+anorm;
	if (temp == anorm) break;
      }
      if (flag) {
	c=0.0;
	s=1.0;
	for (i=l;i<k+1;i++) {
	  f=s*rv1[i];
	  rv1[i]=c*rv1[i];
	  temp = fabs(f)+anorm;
	  if (temp == anorm) break;
	  g=w[i];
	  h=pythag(f,g);
	  w[i]=h;
	  h=1.0/h;
	  c=g*h;
	  s = -f*h;
	  for (j=0;j<m;j++) {
	    y=a[j][nm];
	    z=a[j][i];
	    a[j][nm]=y*c+z*s;
	    a[j][i]=z*c-y*s;
	  }
	}
      }
      z=w[k];
      if (l == k) {
	if (z < 0.0) {
	  w[k] = -z;
	  for (j=0;j<n;j++) v[j][k] = -v[j][k];
	}
	break;
      }

      if (its == 29) 
	{
	  Helper::warn("cannot converge SVD, perhaps due to multi-colinearity"); 
	  return false;
	}

      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=Statistics::pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+Statistics::SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
	i=j+1;
	g=rv1[i];
	y=w[i];
	h=s*g;
	g=c*g;
	z=pythag(f,h);
	rv1[j]=z;
	c=f/z;
	s=h/z;
	f=x*c+g*s;
	g=g*c-x*s;
	h=y*s;
	y *= c;
	for (jj=0;jj<n;jj++) {
	  x=v[jj][j];
	  z=v[jj][i];
	  v[jj][j]=x*c+z*s;
	  v[jj][i]=z*c-x*s;
	}
	z=pythag(f,h);
	w[j]=z;
	if (z) {
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f=c*g+s*y;
	x=c*y-s*g;
	for (jj=0;jj<m;jj++) {
	  y=a[jj][j];
	  z=a[jj][i];
	  a[jj][j]=y*c+z*s;
	  a[jj][i]=z*c-y*s;
	}
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  return true;
}

void Statistics::svbksb( Data::Matrix<double> & u , Data::Vector<double> & w , Data::Matrix<double> & v, Data::Vector<double> & b, Data::Vector<double> & x)
{
  
  int m = u.dim1();
  int n = u.dim2();

  Data::Vector<double> tmp(n);
  
  for (int j=0; j<n; j++) 
    {
      double s = 0.0;
      if ( w[j] != 0.0) 
	{
	  for (int i=0; i<m; i++) s += u(i,j)*b[i];
	  s /= w[j];
	}
      tmp[j] = s;
    }
  
  for (int j=0; j<n; j++) 
    {
      double s = 0.0;
      for (int jj=0; jj<n; jj++) 
	s += v(j,jj) * tmp[jj]; 
      x[j]=s;
    }
}


double Statistics::pythag(const double a, const double b)
{
  double absa,absb;
 
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

inline double SQR(double a)
{
  return a*a;
}





std::vector<double> Statistics::canonical_correlation( const Data::Matrix<double> & x , const Data::Matrix<double> & y , double * pv )
{
  
  // 1. Partitioned covariance matrix  
  //    S_XX  S_YX
  //    S_XY  S_YY
              
  const int nx = x.dim2();
  const int ny = y.dim2();
  const int ne = nx < ny ? nx : ny ;             
  
  // total covariance matrix
  
  if ( x.dim1() != y.dim1() ) Helper::halt("different number of individuals on left and right hand of canonical correlation");
  int nind = x.dim1(); // total N individuals

  Data::Matrix<double> I11 = covariance_matrix( x, x );
  Data::Matrix<double> I12 = covariance_matrix( x, y );
  Data::Matrix<double> I21 = covariance_matrix( y, x );
  Data::Matrix<double> I22 = covariance_matrix( y, y );
  
  Data::Matrix<double> I11b( nx,nx );
  Data::Matrix<double> I22b( ny,ny );

  // 2. Calculate the p x p matrix M1 = inv(sqrt(sig11)) %*% sig12 %*% inv(sig22) %*% sig21 %*% inv(sqrt(sig11))

  bool flag = true;
  I11 = Statistics::matrix_sqrt( I11 );
  I11 = Statistics::inverse(I11,&flag);
  if ( ! flag ) Helper::warn( "could not invert matrix in canonical_correlation()" );
  I22 = Statistics::inverse(I22,&flag);
  if ( ! flag ) Helper::warn( "could not invert matrix in canonical_correlation()" );
  I22b = Statistics::matrix_sqrt( I22b ); // For Step 4b
  I22b = Statistics::inverse( I22b , &flag );
  if ( ! flag ) Helper::warn( "could not invert matrix in canonical_correlation()" );
  I11b = Statistics::inverse( I11b, &flag );
  if ( ! flag ) Helper::warn( "could not invert matrix in canonical_correlation()" );
  
  
  Data::Matrix<double> M1 = Statistics::matrix_multiply( 
                              Statistics::matrix_multiply( 
                              Statistics::matrix_multiply( 
                              Statistics::matrix_multiply( I11 , I12 ) , I22 ) , I21 ) , I11 );

  // Compute and sort eigen values only 

  bool okay = true;
  std::vector<double> sorted_eigenvalues = Statistics::as_vector( Statistics::eigenvalues( M1 , & okay ) );
  std::sort( sorted_eigenvalues.begin(), sorted_eigenvalues.end(), std::greater<double>() );

  // Display largest canonical correlation and its position

  // Use Bartlett's test to get p-values for each canonical correlation

  if ( pv ) *pv = Statistics::bartlett( nind, nx , ny , sorted_eigenvalues ); 
  
  return sorted_eigenvalues;
  
}


double Statistics::bartlett(const int N, 
			    const int p, 
			    const int q, 
			    const std::vector<double> & eigen )
{
  int p2 = p < q ? p : q; // Number of canonical correlations  
  double prod_eigen=1.0;
  for (int j=0; j<p2; j++)
    prod_eigen *= (1-eigen[j]);      
  double chisq = -1*(N - 1 - 0.5*(p+q+1)) * log(prod_eigen);      
  return chi2_prob( chisq, p*q );
}


double Statistics::chi2_prob(double x, double df)
{

  if ( ! Helper::realnum(x) ) return -9;

  double p, q;
  int st = 0;      // error variable
  int w = 1;      // function variable
  double bnd = 1; // boundary function

  // NCP is set to 0
  cdfchi(&w,&p,&q,&x,&df,&st,&bnd);

  // Check status
  if (st != 0 ) return -9;

  // Return p-value
  return q;
  
}


double Statistics::noncentral_chi2_prob( double x , double df , double ncp )
{
  int w = 1;
  double bnd = 1;
  int st = 0;
  double p, q;  
  cdfchn(&w,&p,&q,&x,&df,&ncp,&st,&bnd);
  return q;
}


// Inverse normal distribution

/*
 * Lower tail quantile for standard normal distribution function.
 *
 * This function returns an approximation of the inverse cumulative
 * standard normal distribution function.  I.e., given P, it returns
 * an approximation to the X satisfying P = Pr{Z <= X} where Z is a
 * random variable from the standard normal distribution.
 *
 * The algorithm uses a minimax approximation by rational functions
 * and the result has a relative error whose absolute value is less
 * than 1.15e-9.
 *
 * Author:      Peter J. Acklam
 * Time-stamp:  2002-06-09 18:45:44 +0200
 * E-mail:      jacklam@math.uio.no
 * WWW URL:     http://www.math.uio.no/~jacklam
 *
 * C implementation adapted from Peter's Perl version
 */


/* Coefficients in rational approximations. */

static const double a[] =
  {
    -3.969683028665376e+01,
    2.209460984245205e+02,
    -2.759285104469687e+02,
    1.383577518672690e+02,
    -3.066479806614716e+01,
     2.506628277459239e+00
  };

static const double b[] =
  {
    -5.447609879822406e+01,
    1.615858368580409e+02,
    -1.556989798598866e+02,
    6.680131188771972e+01,
    -1.328068155288572e+01
  };

static const double c[] =
  {
    -7.784894002430293e-03,
    -3.223964580411365e-01,
    -2.400758277161838e+00,
    -2.549732539343734e+00,
    4.374664141464968e+00,
     2.938163982698783e+00
  };

static const double d[] =
  {
    7.784695709041462e-03,
    3.224671290700398e-01,
    2.445134137142996e+00,
    3.754408661907416e+00
  };

#define LOW 0.02425
#define HIGH 0.97575

double Statistics::ltqnorm( double p )
{
  
  double q, r;
  
  if (p < 0 || p > 1)
    {
      return 0.0;
    }
  else if (p == 0)
    {
      return -HUGE_VAL /* minus "infinity" */;
    }
  else if (p == 1)
    {
      return HUGE_VAL /* "infinity" */;
    }
  else if (p < LOW)
    {
      /* Rational approximation for lower region */
      q = sqrt(-2*log(p));
      return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
        ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    }
  else if (p > HIGH)
    {
      /* Rational approximation for upper region */
      q  = sqrt(-2*log(1-p));
      return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
        ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    }
  else
    {
      /* Rational approximation for central region */
      q = p - 0.5;
      r = q*q;
      return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
        (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
    }
}



double Statistics::normden(double scr, double mean, double var)
{
  return (1.0 / sqrt(2*M_PI*var)) * exp( - ( ((scr-mean)*(scr-mean)) / (2*var) ) ) ;
}

double Statistics::t_prob(double T, double df)
{

  if ( ! Helper::realnum(T) ) return -9; 

  T = fabs(T);
  
  double p, q;
  int st = 0;      // error variable
  int w = 1;       // function variable
  double bnd = 1;  // boundary function
  
  // NCP is set to 0
  cdft(&w,&p,&q,&T,&df,&st,&bnd);
  
  // Check status
  if (st != 0 ) return -9;
  
  // Return two-sided p-value
  return 2*q;
  
}

  
//   // Sort evectors. Rows must be ordered according to cancor value (highest first)
  
//   Data::Matrix<double> sorted_eigenvectors = eigen.z;
  
//   std::vector<int> order_eigenvalues(nx);

//   for (int i=0; i<nx; i++)
//     {
 
//      // Determine position of the vector associated with the ith cancor

//       for (int j=0; j<n1; j++)
// 	{

// 	  if ( eigen.d[j] == sorted_eigenvalues_gene1[i] )               
// 	    {
// 	      if (i==0)        
// 		{
// 		  order_eigenvalues_gene1[i]=j;
// 		  break;
// 		}
// 	      else
// 		{
// 		  if (j!=order_eigenvalues_gene1[i-1])
//                     {
//                       order_eigenvalues_gene1[i]=j;
//                       break;
//                     }
// 		}
// 	    }
// 	}
//     }

//   for (int i=0; i<n1; i++)
//     {
//       sorted_eigenvectors_gene1[i] = eigen.z[order_eigenvalues_gene1[i]];
//     }

//   //   cout << "Eigenvector matrix - unsorted:\n";
//   // display(gene1_eigen.z);
//   //cout << "Eigenvector matrix - sorted:\n";
//   //display(sorted_eigenvectors_gene1);


//   ////////////////////////////////////////////////////////
//   // Step 4b. Calculate the q x q eigenvectors of M2 (f). These are
//   // required to compute the coefficients used to build the p
//   // canonical variates b[k] for gene2 (see below). The first p are
//   // given by: f[k] = (1/sqrt(eigen[k])) * inv_sqrt_I22 %*% I21 %*%
//   // inv_sqrt_sig11 %*% e[k] for (k in 1:p) { e.vectors.gene2[,k] =
//   // (1/sqrt(e.values[k])) * inv.sqrt.sig22 %*% sig21 %*%
//   // inv.sqrt.sig11 %*% e.vectors.gene1[,k] }
           
//   matrix_t M2;

//   multMatrix(I22b, I21, tmp);
//   multMatrix(tmp, I11b, M2);
//   multMatrix(M2, I12, tmp);
//   multMatrix(tmp, I22b, M2);
//   Eigen gene2_eigen = eigenvectors(M2);

//   //cout << "Eigenvalues Gene 2 - unsorted:\n";

//   //display(gene2_eigen.d);
 
//   // Sort evalues for gene2
//   vector<double> sorted_eigenvalues_gene2 = gene2_eigen.d;
//   sort(sorted_eigenvalues_gene2.begin(),sorted_eigenvalues_gene2.end(),greater<double>());

//   // Sort eigenvectors for gene2
//   matrix_t sorted_eigenvectors_gene2 = gene2_eigen.z;
//   vector<int> order_eigenvalues_gene2(ny);

//   for (int i=0; i<ny; i++)
//     {
//       // Determine position of the vector associated with the ith cancor
//       for (int j=0; j<n2; j++)
// 	{
// 	  if (gene2_eigen.d[j]==sorted_eigenvalues_gene2[i])
// 	    {
// 	      if (i==0)
// 		{
// 		  order_eigenvalues_gene2[i]=j;
// 		  break;
// 		}
// 	      else
// 		{
// 		  if (j!=order_eigenvalues_gene2[i-1])
// 		    {
// 		      order_eigenvalues_gene2[i]=j;
// 		      break;
// 		    }
// 		}
// 	    }
// 	}
//     }

//   for (int i=0; i<n2; i++)
//     {
//       sorted_eigenvectors_gene2[i] = gene2_eigen.z[order_eigenvalues_gene2[i]];
//     }

//   //cout << "Eigenvector matrix Gene 2 - unsorted:\n";
//   //display(gene2_eigen.z);

//   //cout << "Eigenvector matrix Gene 2 - sorted:\n";
//   //display(sorted_eigenvectors_gene2);

//   //exit(0);


//   //////////////////////////////////////////////////////////////////////////////////
//   // Step 5 - Calculate the gene1 (pxp) and gene2 (pxq) coefficients
//   // used to create the canonical variates associated with the p
//   // canonical correlations

//   transposeMatrix(gene1_eigen.z);
//   transposeMatrix(gene2_eigen.z);

//   matrix_t coeff_gene1;
//   matrix_t coeff_gene2;

//   multMatrix(gene1_eigen.z, I11, coeff_gene1);
//   multMatrix(gene2_eigen.z, I22b, coeff_gene2);

//   //cout << "Coefficients for Gene 1:\n";
//   //display(coeff_gene1);

//   //cout << "Coefficients for Gene 2:\n";
//   //display(coeff_gene2);

//   //exit(0);

//   ///////////////////////////////////////////////////////////////////////
//   // Step 6 - Compute the gene1 and gene2 canonical variates
//   // associated with the highest canonical correlation NOTE: the
//   // original variables of data need to have the mean subtracted  first!
//   // Otherwise, the resulting correlation between variate.gene1 and
//   // variate.gene1 != estimated cancor.

//   // For each individual, eg compos.gene1 =
//   // evector.gene1[1]*SNP1.gene1 + evector.gene1[2]*SNP2.gene1 + ...

//   /////////////////////////////////
//   // Consider each SNP in gene1

//   vector<double> gene1(nind);
            
//   for (int j=0; j<n1; j++)
//     {

//       CSNP * ps = pSNP[j];


//       ///////////////////////////
//       // Iterate over individuals

//       for (int i=0; i< P.n ; i++)
// 	{

// 	  // Only need to look at one perm set
// 	  bool a1 = ps->one[i];
// 	  bool a2 = ps->two[i];

// 	  if ( a1 )
// 	    {
// 	      if ( a2 ) // 11 homozygote
// 		{
// 		  gene1[i] += (1 - mean[j]) * coeff_gene1[order_eigenvalues_gene1[cancor1_pos]][j];
// 		}
// 	      else      // 12 
// 		{
// 		  gene1[i] += (0 - mean[j]) * coeff_gene1[order_eigenvalues_gene1[cancor1_pos]][j];
// 		}
// 	    }
// 	  else
// 	    {
// 	      if ( a2 )      // 21
// 		{
// 		  gene1[i] += (0 - mean[j]) * coeff_gene1[order_eigenvalues_gene1[cancor1_pos]][j];
// 		}
// 	      else           // 22 homozygote
// 		{
// 		  gene1[i] += (-1 - mean[j]) * coeff_gene1[order_eigenvalues_gene1[cancor1_pos]][j];
// 		}
// 	    }

// 	} // Next individual

//     } // Next SNP in gene1

//   /////////////////////////////////
//   // Consider each SNP in gene2
//   vector<double> gene2(P.n);
//   int cur_snp = -1;            
//   for (int j=n1; j<n1+n2; j++)
//     {

//       cur_snp++;
//       CSNP * ps = pSNP[j];


                
//       // Iterate over individuals

//       for (int i=0; i<P.n; i++)
// 	{
                    
// 	  // Only need to look at one perm set
// 	  bool a1 = ps->one[i];
// 	  bool a2 = ps->two[i];

// 	  if ( a1 )
// 	    {
// 	      if ( a2 ) // 11 homozygote
// 		{
// 		  gene2[i] += (1 - mean[j]) * coeff_gene2[order_eigenvalues_gene2[cancor1_pos]][cur_snp];
// 		}
// 	      else      // 12
// 		{
// 		  gene2[i] += (0 - mean[j]) * coeff_gene2[order_eigenvalues_gene2[cancor1_pos]][cur_snp];
// 		}
// 	    }
// 	  else
// 	    {
// 	      if ( a2 )      // 21
// 		{
// 		  gene2[i] += (0 - mean[j]) * coeff_gene2[order_eigenvalues_gene2[cancor1_pos]][cur_snp];
// 		}
// 	      else           // 22 homozygote
// 		{
// 		  gene2[i] += (-1 - mean[j]) * coeff_gene2[order_eigenvalues_gene2[cancor1_pos]][cur_snp];
// 		}
// 	    }
                    
// 	} // Next individual
                
//     } // Next SNP in gene2


//   // Store gene1.variate and gene2.variate in the multiple_covariates field of P.sample
//   // TO DO: NEED TO CHECK IF FIELDS ARE EMPTY FIRST!

//   for (int i=0; i<P.n; i++)
//     {
//       P.sample[i]->clist.resize(2);
//       P.sample[i]->clist[0] = gene1[i];
//       P.sample[i]->clist[1] = gene2[i];
//     }



  
//  return cc;  
//}


Data::Vector<double> Statistics::eigenvalues( Data::Matrix<double> & a , bool * okay )
{
  *okay = true;
  // 'a' should be a square, symmetric matrix
  int n=a.dim1();
  Data::Vector<double> e(n);
  Data::Vector<double> d(n);
  if ( ! tred2(a,d,e) ) *okay = false; 
  if ( ! tqli(d,e) ) *okay = false;
  return d;
}

// Householder method to reduce real, symmetric matrix
// to tridiagonal form
// Modified to return only eigenvalues.
bool Statistics::tred2( Data::Matrix<double> & a , 
			Data::Vector<double> & d,
			Data::Vector<double> & e)
{
  int l,k,j,i;
  double scale,hh,h,g,f;

  int n=d.dim1();
  for (i=n-1;i>0;i--) {
    l=i-1;
    h=scale=0.0;
    if (l > 0) {
      for (k=0;k<l+1;k++)
	scale += fabs(a[i][k]);
      if (scale == 0.0)
	e[i]=a[i][l];
      else {
	for (k=0;k<l+1;k++) {
	  a[i][k] /= scale;
	  h += a[i][k]*a[i][k];
	}
	f=a[i][l];
	g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
	e[i]=scale*g;
	h -= f*g;
	a[i][l]=f-g;
	f=0.0;
	for (j=0;j<l+1;j++) {
	  // Next statement can be omitted if eigenvectors not wanted
// 	  a[j][i]=a[i][j]/h;
	  g=0.0;
	  for (k=0;k<j+1;k++)
	    g += a[j][k]*a[i][k];
	  for (k=j+1;k<l+1;k++)
	    g += a[k][j]*a[i][k];
	  e[j]=g/h;
	  f += e[j]*a[i][j];
	}
	hh=f/(h+h);
	for (j=0;j<l+1;j++) {
	  f=a[i][j];
	  e[j]=g=e[j]-hh*f;
	  for (k=0;k<j+1;k++)
	    a[j][k] -= (f*e[k]+g*a[i][k]);
	}
      }
    } else
      e[i]=a[i][l];
    d[i]=h;
  }
  // Next statement can be omitted if eigenvectors not wanted
//   d[0]=0.0;
  e[0]=0.0;
  // Contents of this loop can be omitted if eigenvectors not
  //	wanted except for statement d[i]=a[i][i];
  for (i=0;i<n;i++) {
//     l=i;
//     if (d[i] != 0.0) {
//       for (j=0;j<l;j++) {
// 	g=0.0;
// 	for (k=0;k<l;k++)
// 	  g += a[i][k]*a[k][j];
// 	for (k=0;k<l;k++)
// 	  a[k][j] -= g*a[k][i];
//       }
//     }
    d[i]=a[i][i];
//     a[i][i]=1.0;
//     for (j=0;j<l;j++) a[j][i]=a[i][j]=0.0;
  }
  return true;
}

// Modified to return only eigenvalues.
bool Statistics::tqli( Data::Vector<double> & d, Data::Vector<double> & e )
{
  const int MAXIT = 60 ; // hmm, was 30 but SKAT had issues with large genes
  // -- I believe these routines have some issues with convergenece when 
  //    using newer compilers -- should add in a EPS difference test when 
  //    testing for convergence here.

  int m,l,iter,i,k;
  double s,r,p,g,f,dd,c,b;
  double volatile temp;
  int n=d.dim1();
  for (i=1;i<n;i++) e[i-1]=e[i];
  e[n-1]=0.0;
  for (l=0;l<n;l++) {
    iter=0;
    do {
      for (m=l;m<n-1;m++) {
	dd=fabs(d[m])+fabs(d[m+1]);
	temp=fabs(e[m])+dd;
	if (temp == dd) break;
      }
      if (m != l) {
	if (iter++ == MAXIT ) { Helper::warn( "convergence problem in tqli()" ); return false; } 
	g=(d[l+1]-d[l])/(2.0*e[l]);
	r=pythag(g,1.0);
	g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
	s=c=1.0;
	p=0.0;
	for (i=m-1;i>=l;i--) {
	  f=s*e[i];
	  b=c*e[i];
	  e[i+1]=(r=pythag(f,g));
	  if (r == 0.0) {
	    d[i+1] -= p;
	    e[m]=0.0;
	    break;
	  }
	  s=f/r;
	  c=g/r;
	  g=d[i+1]-p;
	  r=(d[i]-g)*s+2.0*c*b;
	  d[i+1]=g+(p=s*r);
	  g=c*r-b;
	  // Next loop can be omitted if eigenvectors not wanted
	  /* for (k=0;k<n;k++) {
	     f=z[k][i+1];
	     z[k][i+1]=s*z[k][i]+c*f;
	     z[k][i]=c*z[k][i]-s*f;
	     } */
	}
	if (r == 0.0 && i >= l) continue;
	d[l] -= p;
	e[l]=g;
	e[m]=0.0;
      }
    } while (m != l);
  }
  return true;
}



Statistics::Eigen Statistics::eigenvectors( Data::Matrix<double> & a , bool * okay )
{

  *okay = true;

  // 'a' should be a square, symmetric matrix
  int n=a.dim1();
  
  Statistics::Eigen E(n);
  
  Data::Vector<double> e(n);
  if ( ! Statistics::EV_tred2( a, E.d, e) ) *okay = false;
  if ( ! Statistics::EV_tqli( E.d, e, a) ) *okay = false;
  E.z = a;
  return E;
}


// Householder method to reduce real, symmetric matrix
// to tridiagonal form
// Modified to return both eigenvalues and eigenvectors

bool Statistics::EV_tred2( Data::Matrix<double> & a ,
			   Data::Vector<double> & d , 
			   Data::Vector<double> & e )
{
  int l,k,j,i;
 
  double scale,hh,h,g,f;
  
  int n=d.dim1();
  for (i=n-1;i>0;i--) {
    l=i-1;
    h=scale=0.0;
    if (l > 0) {
      for (k=0;k<l+1;k++)
        scale += fabs(a[i][k]);
      if (scale == 0.0)
        e[i]=a[i][l];
      else {
        for (k=0;k<l+1;k++) {
          a[i][k] /= scale;
          h += a[i][k]*a[i][k];
        }
        f=a[i][l];
        g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
        e[i]=scale*g;
        h -= f*g;
        a[i][l]=f-g;
        f=0.0;
        for (j=0;j<l+1;j++) {
          a[j][i]=a[i][j]/h;
          g=0.0;
          for (k=0;k<j+1;k++)
            g += a[j][k]*a[i][k];
          for (k=j+1;k<l+1;k++)
            g += a[k][j]*a[i][k];
          e[j]=g/h;
          f += e[j]*a[i][j];
        }
        hh=f/(h+h);
        for (j=0;j<l+1;j++) {
          f=a[i][j];
          e[j]=g=e[j]-hh*f;
          for (k=0;k<j+1;k++)
            a[j][k] -= (f*e[k]+g*a[i][k]);
        }
      }
    } else
      e[i]=a[i][l];
    d[i]=h;
  }

  d[0]=0.0;
  e[0]=0.0;

  for (i=0;i<n;i++) {
    l=i;
    if (d[i] != 0.0) {
      for (j=0;j<l;j++) {
        g=0.0;
        for (k=0;k<l;k++)
          g += a[i][k]*a[k][j];
        for (k=0;k<l;k++)
          a[k][j] -= g*a[k][i];
      }
    }
    d[i]=a[i][i];
    a[i][i]=1.0;
    for (j=0;j<l;j++) a[j][i]=a[i][j]=0.0;
  }
  return true;
}


bool Statistics::EV_tqli( Data::Vector<double> & d , 
			  Data::Vector<double> & e , 
			  Data::Matrix<double> & z )
{
  int m,l,iter,i,k;
  double s,r,p,g,f,dd,c,b;
  
  int n=d.dim1();
  for (i=1;i<n;i++) e[i-1]=e[i];
  e[n-1]=0.0;
  for (l=0;l<n;l++) {
    iter=0;
    do {
      for (m=l;m<n-1;m++) {
        dd=fabs(d[m])+fabs(d[m+1]);
        if (fabs(e[m])+dd == dd) break;
      }
      if (m != l) {
        if (iter++ == 30) { Helper::warn("convergence issue in EVtqli()"); return false; } 
        g=(d[l+1]-d[l])/(2.0*e[l]);
        r=Statistics::pythag(g,1.0);
        g=d[m]-d[l]+e[l]/(g+Statistics::SIGN(r,g));
        s=c=1.0;
        p=0.0;
        for (i=m-1;i>=l;i--) {
          f=s*e[i];
          b=c*e[i];
          e[i+1]=(r=Statistics::pythag(f,g));
          if (r == 0.0) {
            d[i+1] -= p;
            e[m]=0.0;
            break;
          }
          s=f/r;
          c=g/r;
          g=d[i+1]-p;
          r=(d[i]-g)*s+2.0*c*b;
          d[i+1]=g+(p=s*r);
          g=c*r-b;

          for (k=0;k<n;k++) {
            f=z[k][i+1];
            z[k][i+1]=s*z[k][i]+c*f;
            z[k][i]=c*z[k][i]-s*f;
          }
        }
        if (r == 0.0 && i >= l) continue;
        d[l] -= p;
        e[l]=g;
        e[m]=0.0;
      }
    } while (m != l);
  }
  return true;
}


Data::Matrix<double> Statistics::matrix_multiply( const Data::Matrix<double> & a, const Data::Matrix<double> & b )
{
  //  int ar ac x br bc
  if ( a.dim2() != b.dim1() ) Helper::halt("non-conformable matrix multiplication requested");     

  const int nrow = a.dim1();
  const int ncol = b.dim2();
  const int nk = a.dim2();
  Data::Matrix<double> r(nrow,ncol);
  for (int i=0;i<nrow;i++)
    for(int j=0; j<ncol;j++)
      for (int k=0; k<nk; k++)
	r(i,j) += a(i,k) * b(k,j); 	
  return r;
}

Data::Vector<double> Statistics::matrix_multiply( const Data::Matrix<double> & a , const Data::Vector<double> & b)
{
  if ( a.dim2() != b.dim1() ) Helper::halt("non-conformable matrix multiplication requested");     
  Data::Vector<double> r( a.dim1() );
  const int nrow = a.dim1();
  const int nk = a.dim2();
  for (int i=0;i<nrow;i++)
    for (int k=0; k<nk; k++)
      r(i) += a(i,k) * b(k); 	
  return r;
}

Data::Vector<double> Statistics::matrix_multiply( const Data::Vector<double> & a , const Data::Matrix<double> & b )
{
  if ( a.dim1() != b.dim1() ) Helper::halt("non-conformable matrix multiplication requested");     
  Data::Vector<double> r( b.dim2() );
  const int nrow = b.dim2();
  const int nk = a.dim1();
  for (int i=0;i<nrow;i++)
    for (int k=0; k<nk; k++)
      r(i) += a(k) * b(k,i);
  return r;
}

double Statistics::matrix_inner_product( const Data::Vector<double> & a , const Data::Vector<double> & b )
{
  if ( a.dim1() != b.dim1() ) 
    {
      Helper::warn("internal error: non-comformable inner-product" );
      return 0;
    }
  double r;
  for (int i = 0 ; i < a.dim1() ; i++ )
    r += a[i] * b[i];
  return r;
}

Data::Matrix<double> Statistics::matrix_outer_product( const Data::Vector<double> & a , const Data::Vector<double> & b )
{
  Data::Matrix<double> r(a.dim1() , b.dim1() );
  for (int i = 0 ; i < r.dim1() ; i++ )
    for (int j = 0 ; j < r.dim2() ; j++ )
      r[i][j] = a[i] * b[j];
  return r;
}



long unsigned int Statistics::factorial(int n)
{
  long unsigned int z = 1;
  for (int i = 1 ; i <= n ; i++)
    z *= i; 
  return z;
}

long unsigned int Statistics::combin(int n, int k)
{
  if (k>n) return 0;
  long double z = 1;
  int r = k;
  if ( k > (n-k) ) r = n-k;
  for (int i = 0 ; i <= r-1 ; i++)
    z *= (long double)(n-i) / (long double)(r-i);
  return (long unsigned int)z;
}

long double factorial(int x) {
  int i;
  long double result = 1;
  for (i = 2; i <= x; i++)
    result *= i;
  return result;
}


/**
 * AUTHOR
 *   Catherine Loader, catherine\research.bell-labs.com.
 *   October 23, 2000.
 *
 *  Merge in to R:
 * Copyright (C) 2000, The R Core Development Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
 *
 *
 * DESCRIPTION
 *
 *   To compute the binomial probability, call dbinom(x,n,p).
 *   This checks for argument validity, and calls dbinom_raw().
 *
 *   dbinom_raw() does the actual computation; note this is called by
 *   other functions in addition to dbinom()).
 *     (1) dbinom_raw() has both p and q arguments, when one may be represented
 *         more accurately than the other (in particular, in df()).
 *     (2) dbinom_raw() does NOT check that inputs x and n are integers. This
 *         should be done in the calling function, where necessary.
 *     (3) Also does not check for 0 <= p <= 1 and 0 <= q <= 1 or NaN's.
 *         Do this in the calling function.
 */

inline double rd0(const bool give_log) {
  return 0.0;
}

inline double rd1(const bool give_log) {
  return give_log ? 0.0 : 1.0;
}

inline double rdexp(const double x, const bool give_log) {
  return give_log ? x : exp(x);
}

inline double rdfexp(const double f, const double x, const bool give_log) {
  return give_log ? -0.5*log(f) + (x) : exp(x) / sqrt(f);
}
/// Natural log of the gamma function
///
/// \param x
///
/// from
///
/// NIST Guide to Available Math Software.
/// Source for module GAMLN from package CMLIB.
/// Retrieved from TIBER on Wed Apr 29 17:30:20 1998.
// =====================================================================
//    WRITTEN BY D. E. AMOS, SEPTEMBER, 1977.
//
//    REFERENCES
//        SAND-77-1518
//
//        COMPUTER APPROXIMATIONS BY J.F.HART, ET.AL., SIAM SERIES IN
//        APPLIED MATHEMATICS, WILEY, 1968, P.135-136.
//
//        NBS HANDBOOK OF MATHEMATICAL FUNCTIONS, AMS 55, BY
//        M. ABRAMOWITZ AND I.A. STEGUN, DECEMBER. 1955, P.257.
//
//    ABSTRACT
//        GAMLN COMPUTES THE NATURAL LOG OF THE GAMMA FUNCTION FOR
//        X.GT.0. A RATIONAL CHEBYSHEV APPROXIMATION IS USED ON
//        8.LT.X.LT.1000., THE ASYMPTOTIC EXPANSION FOR X.GE.1000. AND
//        A RATIONAL CHEBYSHEV APPROXIMATION ON 2.LT.X.LT.3. FOR
//        0.LT.X.LT.8. AND X NON-INTEGRAL, FORWARD OR BACKWARD
//        RECURSION FILLS IN THE INTERVALS  0.LT.X.LT.2 AND
//        3.LT.X.LT.8. FOR X=1.,2.,...,100., GAMLN IS SET TO
//        NATURAL LOGS OF FACTORIALS.
//
//    DESCRIPTION OF ARGUMENTS
//
//        INPUT
//          X      - X.GT.0
//
//        OUTPUT
//          GAMLN  - NATURAL LOG OF THE GAMMA FUNCTION AT X
//
//    ERROR CONDITIONS
//        IMPROPER INPUT ARGUMENT - A FATAL ERROR
//
double Statistics::gamln(const double x) {
  static double xlim1 = 8.0;
  static double xlim2 = 1e3;
  static double rtwpil = .918938533204673;
  static double p[5] =
    { 7.66345188e-4, -5.9409561052e-4,
      7.936431104845e-4, -.00277777775657725,
      .0833333333333169
    };
  static double q[2] =
    { -.00277777777777778, .0833333333333333 }
  ;
  static double pcoe[9] =
    { .00297378664481017,
      .0092381945590276, .109311595671044, .398067131020357,
      2.15994312846059, 6.33806799938727,
      20.7824725317921, 36.0367725300248, 62.0038380071273
    }
  ;
  static double qcoe[4] =
    { 1.0, -8.90601665949746,
      9.82252110471399, 62.003838007127
    };
  static double gln[100] =
    { 0., 0., .693147180559945,
      1.79175946922806, 3.17805383034795,
      4.78749174278205, 6.5792512120101,
      8.52516136106541, 10.6046029027453,
      12.8018274800815, 15.1044125730755,
      17.5023078458739, 19.9872144956619,
      22.5521638531234, 25.1912211827387,
      27.8992713838409, 30.6718601060807,
      33.5050734501369, 36.3954452080331,
      39.3398841871995, 42.3356164607535,
      45.3801388984769, 48.4711813518352,
      51.6066755677644, 54.7847293981123,
      58.0036052229805, 61.261701761002,
      64.5575386270063, 67.8897431371815,
      71.257038967168, 4.6582363488302,
      78.0922235533153, 81.557959456115,
      85.0544670175815, 88.5808275421977,
      92.1361756036871, 95.7196945421432,
      99.3306124547874, 102.968198614514,
      106.631760260643, 110.320639714757,
      114.034211781462, 117.771881399745,
      121.533081515439, 125.317271149357,
      129.123933639127, 132.952575035616,
      136.802722637326, 140.673923648234,
      144.565743946345, 148.477766951773,
      152.409592584497, 156.360836303079,
      160.331128216631, 164.320112263195,
      168.327445448428, 172.352797139163,
      176.395848406997, 180.456291417544,
      184.533828861449, 188.628173423672,
      192.739047287845, 196.86618167289,
      201.009316399282, 205.168199482641,
      209.342586752537, 213.532241494563,
      217.736934113954, 221.95644181913,
      226.190548323728, 230.439043565777,
      234.701723442818, 238.978389561834,
      243.268849002983, 247.572914096187,
      251.890402209723, 256.22113555001,
      260.564940971863, 264.921649798553,
      269.29109765102, 273.673124285694,
      278.067573440366, 282.47429268763,
      286.893133295427, 291.32395009427,
      295.766601350761, 300.220948647014,
      304.686856765669, 309.164193580147,
      313.652829949879, 318.152639620209,
      322.663499126726, 327.185287703775,
      331.717887196928, 336.261181979198,
      340.815058870799, 345.379407062267,
      349.95411804077, 354.539085519441,
      359.134205369575
    };

  /* System generated locals */
  long int i__1;
  double ret_val = 0.0;

  /* Local variables */
  static double dgam;
  static long int i__;
  static double t, dx, px, qx, rx, xx;
  static long int ndx, nxm;
  static double sum, rxx;

  if (x <= 0.0) {
	   goto L90;
	  } else {
	    goto L5;
	  }
	 L5:
	  ndx = static_cast<long>(x);
	  t = x - static_cast<double>(ndx);
	  if (t == 0.0) {
	    goto L51;
	  }
	  dx = xlim1 - x;
	  if (dx < 0.0) {
	    goto L40;
	  }
	  /*     RATIONAL CHEBYSHEV APPROXIMATION ON 2.LT.X.LT.3 FOR GAMMA(X) */
	  nxm = ndx - 2;
	  px = pcoe[0];
	  for (i__ = 2; i__ <= 9; ++i__) {
	    /* L10: */
	    px = t * px + pcoe[i__ - 1];
	  }
	  qx = qcoe[0];
	  for (i__ = 2; i__ <= 4; ++i__) {
	    /* L15: */
	    qx = t * qx + qcoe[i__ - 1];
	  }
	  dgam = px / qx;
	  if (nxm > 0) {
	    goto L22;
	  }
	  if (nxm == 0) {
	    goto L25;
	  }
	  /*     BACKWARD RECURSION FOR 0.LT.X.LT.2 */
	  dgam /= t + 1.0;
	  if (nxm == -1) {
	    goto L25;
	  }
	  dgam /= t;
	  ret_val = log (dgam);
	  return ret_val;
	  /*     FORWARD RECURSION FOR 3.LT.X.LT.8 */
	 L22:
	  xx = t + 2.0;
	  i__1 = nxm;
	  for (i__ = 1; i__ <= i__1; ++i__) {
	    dgam *= xx;
	    /* L24: */
	    xx += 1.0;
	  }
	 L25:
	  ret_val = log (dgam);
	  return ret_val;
	  /*     X.GT.XLIM1 */
	 L40:
	  rx = 1.0 / x;
	  rxx = rx * rx;
	  if (x - xlim2 < 0.0) {
	    goto L41;
	  }
	  px = q[0] * rxx + q[1];
	  ret_val = px * rx + (x - 0.5) * log (x) - x + rtwpil;
	  return ret_val;
	  /*     X.LT.XLIM2 */
	 L41:
	  px = p[0];
	  sum = (x - 0.5) * log (x) - x;
	  for (i__ = 2; i__ <= 5; ++i__) {
	    px = px * rxx + p[i__ - 1];
	    /* L42: */
	  }
	  ret_val = px * rx + sum + rtwpil;
	  return ret_val;
	  /*     TABLE LOOK UP FOR INTEGER ARGUMENTS LESS THAN OR EQUAL 100. */
	 L51:
	  if (ndx > 100) {
	    goto L40;
	  }
	  ret_val = gln[ndx - 1];
	  return ret_val;
	 L90:
	  return ret_val;
	}


// for stirlerr()
const double S0 = 0.083333333333333333333;       /* 1/12 */
const double S1 = 0.00277777777777777777778;     /* 1/360 */
const double S2 = 0.00079365079365079365079365;  /* 1/1260 */
const double S3 = 0.000595238095238095238095238; /* 1/1680 */
const double S4 = 0.0008417508417508417508417508; /* 1/1188 */

/*
 *  AUTHOR
 *    Catherine Loader, catherine\research.bell-labs.com.
 *    October 23, 2000.
 *
 *  Merge in to R:
 * Copyright (C) 2000, The R Core Development Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
 *
 *
 *  DESCRIPTION
 *
 *    Computes the log of the error term in Stirling's formula.
 *      For n > 15, uses the series 1/12n - 1/360n^3 + ...
 *      For n <=15, integers or half-integers, uses stored values.
 *      For other n < 15, uses lgamma directly (don't use this to
 *        write lgamma!)
 *
 * Merge in to R:
 * Copyright (C) 2000, The R Core Development Team
 * R has lgammafn, and lgamma is not part of ISO C
 *
 */
/* stirlerr(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n ) */
double stirlerr(double n) {
  /*
    error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
  */
  const double sferr_halves[31] = {
    0.0,       /* n=0 - wrong, place holder only */
    0.1534264097200273452913848,        /* 0.5 */
    0.0810614667953272582196702,        /* 1.0 */
    0.0548141210519176538961390,        /* 1.5 */
    0.0413406959554092940938221,        /* 2.0 */
    0.03316287351993628748511048,       /* 2.5 */
    0.02767792568499833914878929,       /* 3.0 */
    0.02374616365629749597132920,       /* 3.5 */
    0.02079067210376509311152277,       /* 4.0 */
    0.01848845053267318523077934,       /* 4.5 */
    0.01664469118982119216319487,       /* 5.0 */
    0.01513497322191737887351255,       /* 5.5 */
    0.01387612882307074799874573,       /* 6.0 */
    0.01281046524292022692424986,       /* 6.5 */
    0.01189670994589177009505572,       /* 7.0 */
    0.01110455975820691732662991,       /* 7.5 */
    0.010411265261972096497478567,       /* 8.0 */
    0.009799416126158803298389475,       /* 8.5 */
    0.009255462182712732917728637,       /* 9.0 */
    0.008768700134139385462952823,       /* 9.5 */
    0.008330563433362871256469318,       /* 10.0 */
    0.007934114564314020547248100,       /* 10.5 */
    0.007573675487951840794972024,       /* 11.0 */
    0.007244554301320383179543912,       /* 11.5 */
    0.006942840107209529865664152,       /* 12.0 */
    0.006665247032707682442354394,       /* 12.5 */
    0.006408994188004207068439631,       /* 13.0 */
    0.006171712263039457647532867,       /* 13.5 */
    0.005951370112758847735624416,       /* 14.0 */
    0.005746216513010115682023589,       /* 14.5 */
    0.005554733551962801371038690  /* 15.0 */
  };

  double nn;
  if (n <= 15.0) {
    nn = n + n;
    if (nn == (int)nn)
      return (sferr_halves[(int)nn]);
    return (Statistics::gamln(n + 1.0) - (n + 0.5)*log(n) + n - M_LN_SQRT_2PI);
  }
  nn = n * n;
  if (n > 500)
    return ((S0 -S1 / nn) / n);
  if (n > 80)
    return ((S0 -(S1 - S2 / nn) / nn) / n);
  if (n > 35)
    return ((S0 -(S1 - (S2 - S3 / nn) / nn) / nn) / n);
  /* 15 < n <= 35 : */
  return ((S0 -(S1 - (S2 - (S3 - S4 / nn) / nn) / nn) / nn) / n);
}


/*
 *  AUTHOR
 * Catherine Loader, catherine\research.bell-labs.com.
 * October 23, 2000.
 *
 *  Merge in to R:
 * Copyright (C) 2000, The R Core Development Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
 *
 *
 *  DESCRIPTION
 * Evaluates the "deviance part"
 * bd0(x,M) :=  M * D0(x/M) = M*[ x/M * log(x/M) + 1 - (x/M) ] =
 *    =  x * log(x/M) + M - x
 * where M = E[X] = n*p (or = lambda), for   x, M > 0
 *
 * in a manner that should be stable (with small relative error)
 * for all x and np. In particular for x/np close to 1, direct
 * evaluation fails, and evaluation is based on the Taylor series
 * of log((1+v)/(1-v)) with v = (x-np)/(x+np).
 */
double bd0(double x, double np) {
  double ej, s, s1, v;
  int j;
  if (fabs(x - np) < 0.1*(x + np)) {
    v = (x - np) / (x + np);
    s = (x - np) * v; /* s using v -- change by MM */
    ej = 2 * x * v;
    v = v * v;
    for (j = 1; ; j++) { /* Taylor series */
      ej *= v;
      s1 = s + ej / ((j << 1) + 1);

      if (s1 == s) /* last term was effectively 0 */
        return (s1);

      s = s1;
    }
  }
  /* else:  | x - np |  is not too small */
  return (x*log(x / np) + np - x);
}





double Statistics::dbinom_raw(const double k, const double n, const double p) {
  double f, lc;
  const double q = 1 - p;
  if (p == 0)
    return ((k == 0) ? rd1(LOG_BOOL) : rd0(LOG_BOOL));
  if (q == 0)
    return ((k == n) ? rd1(LOG_BOOL) : rd0(LOG_BOOL));
  if (k == 0) {
    if (n == 0)
      return rd1(LOG_BOOL);
    lc = (p < 0.1) ? -bd0(n, n * q) - n * p : n * log(q);
    return rdexp(lc, LOG_BOOL);
  }
  if (k == n) {
    lc = (q < 0.1) ? -bd0(n, n * p) - n * q : n * log(p);
    return rdexp(lc, LOG_BOOL);
  }
  if (k < 0 || k > n)
    return rd0(LOG_BOOL);
  lc = stirlerr(n) - stirlerr(k) - stirlerr(n - k) - bd0(k, n * p) - bd0(n - k, n * q);
  f = (M_2PI * k * (n - k)) / n;
  return rdfexp(f, lc, LOG_BOOL);
}


/// Density of the binomial distribution.
///
/// \f[ f(k) = {n \choose k}p^k (1-p)^{n-k} \f]
///
/// Returns probability of getting k successes in n binomial trials
/// with a probability p of success on each trial, if give_log == false.
/// If give_log == true, returns the natural logarithm of the probability.
///
/// \param k          Number of successes
/// \param n          Number of trials
/// \param p          Probability of success on each trial
/// \param give_log   FALSE
///
/// \f[ \mbox{E}(x) = np \f]
/// \f[ \mbox{Var}(x) = np(1-p) \f]
///
/// The implementation uses dbinom_raw from R
///
double Statistics::dbinom(const int k, const int n, double p ) {
  return dbinom_raw(k, n, p);
}




double Statistics::gammln(double xx)
{
  
  //  Returns the value ln[Γ(xx)] for xx > 0. 
  
  // Internal arithmetic will be done in double precision, a nicety
  // that you can omit if five-figure accuracy is good enough.
  
  static double cof[6]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5}; 
  
  int j;
  double y=xx, x=xx; 
  double tmp=x+5.5; 
  tmp -= (x+0.5)*log(tmp); 
  double ser=1.000000000190015; 
  for (j=0;j<=5;j++) ser += cof[j]/++y; 
  return -tmp+log(2.5066282746310005*ser/x);
}

double Statistics::factrl(int n)
{  
  static int ntop=4;
  static double a[33]={1.0,1.0,2.0,6.0,24.0}; // Fill in table only as required.
  int j;
  if (n < 0) { std::cerr << "exit1\n"; }
  if (n > 32 ) return exp(gammln(n+1));
  while (ntop<n) { 
    j=ntop++;
    a[ntop]=a[j]*ntop;
  }
  return a[n];
}

double Statistics::bico(int n, int k)
{  
  // binomial coeff
  return floor( 0.5 + exp( factln(n) -factln(k) -factln(n-k) ) );
}

double Statistics::factln(int n)
{
  static double a[101];
  if (n <= 1) return 0.0;
  if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0));
  else return gammln(n+1.0);
}

double Statistics::beta( double x , double a1 , double a2 )
{
  Helper::halt("not implemented beta()");
  return 0;
}



//double Statistics::integrate_trapezoidal(double a, double b, double (*f)(double x,void*), void*d, double eps)
double Statistics::integrate_old(double a, double b, double (*f)(double x,void*,bool*), bool*okay, void*d, double eps)
{
  // trapezoidal rule
  const double ZEPS = 1e-10;
  double old = update_integral(a, b, f, d, okay, 0.0, 0), result;
  int round = 1;
  *okay = true;

  while (1)
    {
      //std::cout << "round = " << round << "\n";
      result = update_integral(a, b, f, d, okay, old, round++);
      if ( ! *okay ) return 0;
      if ( fabs(result-old) < eps*(fabs(result)+fabs(old))+ZEPS)
	return result;
      old = result;
    }
}


double Statistics::integrate(double a, double b, double (*f)(double x,void*,bool*), bool *okay, void*d,  double eps)
{

  // simpson's rule
  const double ZEPS = 1e-10;
  double old = update_integral(a, b, f, d, okay, 0.0, 0), result;
  double sold = old, sresult;
  int round = 1;
  *okay = true;

  while (1)
    {
      //std::cout << "round = " << round << "\n";
      result = update_integral(a, b, f, d, okay, old, round++);
      if ( ! *okay ) return 0;
      sresult = (4.0 * result - old) / 3.0;

      std::cout << "res = " << result << " " << sresult << " " << sold << " " << fabs( sresult - sold ) << eps*(fabs(sresult)+fabs(sold))+ZEPS ;

      if ( fabs(sresult-sold) <  eps*(fabs(sresult)+fabs(sold))+ZEPS )
	return sresult;
      old = result; sold=sresult;
    }
}


double Statistics::update_integral(double a, double b,
				   double (*f)(double x,void*,bool*), void * d, bool * okay, 
				   double previous, int round)
{
  double h, sum;
  int i, n = 1 << (round - 1);
  if (round == 0)
    return 0.5 * ((*f)(a,d,okay) + (*f)(b,d,okay)) * (b - a);
  sum = previous * n / (b - a);
  h = (b - a) / (2 * n);
  for (int i = 1; i < 2 * n; i += 2)
    sum += (*f)(a + i*h,d,okay);
  return sum * h;
}


// double Statistics::trapzd( double a, double b, double (*f)(double x,void*,bool*),void * d, bool *, int )
// {
//   //
// }



void Statistics::polint( double * xa, 
			 double * ya, 
			 int n, 
			 double x, 
			 double *y, 
			 double *dy )
{
  // returns values 'y' and error estimate 'dy'
  int m,ns=1;
  double den,dift,ho,hp,w;

  double dif=fabs(x-xa[1]);
  Data::Vector<double> c(n+1);
  Data::Vector<double> d(n+1);
  
  for (int i=1;i<=n;i++) 
    {
      if ( (dift=fabs(x-xa[i])) < dif)
	{
	  ns=i;
	  dif=dift;
	}
      c[i]=ya[i];
      d[i]=ya[i];
    }

  *y=ya[ns--];

  for (int m=1;m<n;m++)
    {
      for (int i=1;i<=n-m;i++)
	{
	  ho=xa[i]-x;
	  hp=xa[i+m]-x;
	  w=c[i+1]-d[i];
	  if ( (den=ho-hp) == 0.0) Helper::halt("error in polint");
	  den=w/den;
	  d[i]=hp*den;
	  c[i]=ho*den;
	}
      *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
    }
}

double Statistics::qromo( double a, double b, double (*f)(double x,void*,bool*), bool * okay , void*d , int method , double eps )
{

  const double EPS = eps;
  const int JMAX = 14;
  const int JMAXP = JMAX + 1;
  const int K = 5;
  double previous = 0;
  
  double ss = 0;
  double dss = 0;

  Data::Vector<double> h( JMAXP+1);
  Data::Vector<double> s( JMAXP );
  
  // note -- uses 1-based indexing here...

  h[1] = 1.0;
  
  for (int j=1;j<=JMAX;j++)
    {
      switch ( method ) 
	{
	case 1 : s[j]=midpnt(a,b,f,previous,d,okay,j); break;
	case 2 : s[j]=midpnt(a,b,f,previous,d,okay,j); break;
	case 3 : s[j]=midsql(a,b,f,previous,d,okay,j); break;
	// case 4 : s[j]=midsqu(a,b,f,previous,d,okay,j); break;
	// case 5 : s[j]=midexp(a,b,f,previous,d,okay,j); break;  
	}
      
      if ( ! *okay ) return 0;

      if (j >= K) 
	{
	  polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
	  if ( fabs(dss) <= EPS*fabs(ss) ) return ss;
	}
      h[j+1]=h[j]/9.0;      

      previous = s[j];
      
    }
  *okay = false; // too many steps
  return 0;
}
  

double Statistics::midpnt( double aa, double bb, double (*f)(double x,void*,bool*), double previous, void * d, bool * okay , int n )
{
  double b = 1/aa;
  double a = 1/bb;
  if ( n == 1 ) 
    return (b-a) * f( 0.5 * (a+b) , d , okay ) ;
  
  int it=1;
  for (int j=1;j<n-1;j++) it *= 3;
  double tnm = it;
  double del = (b-a)/(3.0*tnm);
  double ddel = del+del;
  double x = a + 0.5 * del;
  double sum = 0;
  for (int j=0; j<it;j++)
    {
      sum += f(x,d,okay);
      x += ddel;
      sum += f(x,d,okay);
      x += del;
    }
  return (previous+(b-a)*sum/tnm)/3.0;
}
  


double Statistics::midsql( double aa, double bb, double (*f)(double x,void*,bool*), double previous, void * d, bool * okay , int n )
{
  double b = sqrt(bb-aa);
  double a = 0.0;

  if ( n == 1 ) 
    return (b-a) * f( 0.5 * (a+b) , d , okay ) ;

#define FUNC(x) (2.0*(x)*(*funk)(aa+(x)*(x)))

  int it=1;
  for (int j=1;j<n-1;j++) it *= 3;
  double tnm = it;
  double del = (b-a)/(3.0*tnm);
  double ddel = del+del;
  double x = a + 0.5 * del;
  double sum = 0;
  for (int j=0; j<it;j++)
    {
      sum += 2 * x * f(aa + x*x  ,d,okay);
      x += ddel;
      sum += 2 * x * f(aa + x*x  ,d,okay);
      x += del;
    }
  return (previous+(b-a)*sum/tnm)/3.0;
}
  








double Statistics::qsimp( double a, double b, double (*f)(double x,void*,bool*), bool * okay , void * d , double eps )
{
  const int JMAX = 15;
  const double EPS = eps;
  double s, st, ost=0.0, os=0.0; 
  for (int j=0; j<JMAX; j++)
    {
      st = midpnt(a,b,f,os,d,okay,j);
      s=(9.0*st-ost)/8.0;
      if ( j > 5 ) 
	{
	  if ( fabs( s-os ) < EPS * fabs(os) || 
	       ( s == 0.0 && os == 0.0 ) ) return s;
	  os = s;
	  ost = st;
	}
    }
  *okay = false;
  return 0;
}



Data::Matrix<double> Statistics::cholesky( const Data::Matrix<double> & b )
{

  // return Cholesky decomposition in upper-triangular matrix

  if ( b.dim1() != b.dim2() ) Helper::halt("cholesky of non-square matrix requested");

  const int n = b.dim1();

  Data::Matrix<double> a = b;
  
  if ( n == 0 ) Helper::halt("cholesky: 0-element matrix");

  for (int i=0;i<n;i++)
    for (int j=i;j<n;j++)
      {
	double sum = a(i,j);	
	for ( int k = i-1; k>=0; k--) 
	  sum -= a(i,k) * a(j,k);
	if ( i == j ) 
	  {
	    if ( sum <= 0.0 ) Helper::halt("cholesky failed");
	    a(i,i) = sqrt(sum);
	  }
	else
	  {
	    a(j,i) = sum/a(i,i);
	    a(i,j) = 0.0;
	  }
      }
  return a;
}


bool Statistics::qchisq( double q , double df , double * x )
{
  
  if ( ! Helper::realnum(q) ) return false;
  else if ( q>=1 ) return 0;
  
  double p = 1 - q;
  int st = 0;      // error variable
  int w = 2;       // function variable
  double bnd = 1;  // boundary function
  
  // NCP is set to 0
  cdfchi(&w,&p,&q,x,&df,&st,&bnd);
  if (st != 0 ) return false;
  return true; 
}


double Statistics::dchisq( double x , double df )
{
  // return dgamma(x, df / 2., 2., give_log);
  return -1;
}


bool Statistics::t_test( double u1, double s1, int n1 ,
			 double u2 , double s2 , int n2 ,
			 double * pvalue , double * p_lower , double * p_upper )
{
  // Welch's t-test, unequal samples, unequal variances
  if ( n1 < 2 || n2 < 2 ) return false;
  if ( s1 <= 0 || s2 <= 0 ) return false;
  if ( pvalue == NULL ) return false;
  double x1 = s1 / (double) n1;
  double x2 = s2 / (double) n2;
  double t = ( u1 - u2 ) / sqrt( x1 + x2 );
  double df = ( ( x1 + x2 ) * ( x1 + x2 ) ) / ( x1*x1/(double)(n1-1) + x2*x2/(double)(n2-1)  );
  *pvalue = t_prob( t , df );
  if ( p_lower ) *p_lower = u1 < u2 ? *pvalue * 0.5 : 1.0 ;
  if ( p_upper ) *p_upper = u1 > u2 ? *pvalue * 0.5 : 1.0 ;
  return true;
}



double Statistics::correlation( const std::vector<double> & x , const std::vector<double> & y )
{
  // basic correlation 

  // r = cov(1,2) / sqrt( var(1).var(2) )                                                                                                                 
  double X = 0;
  double X2 = 0;
  double Y = 0;
  double Y2 = 0;
  double XY = 0;

  const int n = x.size();
  if ( y.size() != n ) Helper::halt("error in correl()");

  for (int i=0; i<n; i++)
    {
      X += x[i];
      Y += y[i];
      XY += x[i]*y[i];
      // Sum squares	   
      X2 += x[i]*x[i];
      Y2 += y[i]*y[i];
    }
  
  X /= (double)n;
  X2 /= (double)n;
  Y /= (double)n;
  Y2 /= (double)n;
  XY /= (double)n;

  double var1 = X2 - X*X;
  double var2 = Y2 - Y*Y;
  double cov12 = XY - X*Y;

  double r = cov12 / (sqrt(var1)*sqrt(var2));
  
  return r;
}


double Statistics::mean( const Data::Vector<double> & x )
{
  double s = 0;
  const int n = x.size();
  for (int i=0;i<n;i++) s += x[i];
  s /= (double)n;
  return s;  
}


double Statistics::variance( const Data::Vector<double> & x ) 
{ 
  double m = mean(x);
  double ss = 0;
  const int n = x.size();
  for (int i=0;i<n;i++) ss += ( x[i] - m ) * ( x[i] - m );
  ss /= (double)n;
  return ss;
}
