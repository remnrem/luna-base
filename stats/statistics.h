
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

#ifndef __STATISTICS_H__
#define __STATISTICS_H__

#include "matrix.h"
#include "fisher.h"

#include <vector>

namespace Statistics { 
  
  template<class T> inline const T SQR(const T a) {return a*a;}
  template<class T> inline const T FNMAX(const T &a, const T &b) {return b > a ? (b) : (a);}
  template<class T> inline const T FNMIN(const T &a, const T &b) {return b < a ? (b) : (a);}
  template<class T> inline const T SIGN(const T &a, const T &b) {return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}
  template<class T> inline void SWAP(T &a, T &b) {T dum=a; a=b; b=dum;}

  std::vector<double> as_vector( const Data::Vector<double> & );
  
  // Singular Value Decomposition
  bool svdcmp( Data::Matrix<double> & , Data::Vector<double> & , Data::Matrix<double> & );
  void svbksb( Data::Matrix<double> & , Data::Vector<double> & , Data::Matrix<double> & , Data::Vector<double> & , Data::Vector<double> & );

  // 1-dimensional numerical integration
  double integrate(double a, double b, double (*f)(double x,void*,bool*), bool * , void*d = NULL , double eps = 10e-15 );  
  double integrate_old(double a, double b, double (*f)(double x,void*,bool*), bool*,void*d, double eps);
  double update_integral(double a, double b, double (*f)(double x,void*,bool*),void * d, bool *, double previous, int round);

  // NR 1-dimensional integration
  
  double qsimp( double a, double b, double (*f)(double x,void*,bool*), bool * okay , void*d = NULL , double eps = 10e-6 );
  
  //  double trapzd( double a, double b, double (*f)(double x,void*,bool*),void * d, bool *, int );

  // (Open) Romberg integration
  
  void polint( double * xa , 
	       double * ya , 
	       int n , 
	       double x, double *y, double *dy );

  double qromo( double a, double b, double (*f)(double x,void*,bool*), bool * okay , void*d = NULL , int method = 1 , double eps = 10e-6 );
  
  double midpnt( double aa, double bb, double (*f)(double x,void*,bool*), double previous, void * d, bool * okay , int n ); //method1 
  double midinf( double aa, double bb, double (*f)(double x,void*,bool*), double previous, void * d, bool * okay , int n ); //2
  double midsql( double aa, double bb, double (*f)(double x,void*,bool*), double previous, void * d, bool * okay , int n ); //3
  double midsqu( double aa, double bb, double (*f)(double x,void*,bool*), double previous, void * d, bool * okay , int n ); //4
  double midexp( double aa, double bb, double (*f)(double x,void*,bool*), double previous, void * d, bool * okay , int n ); //5


  
  
  long unsigned int factorial(int n);
  long unsigned int combin(int n, int k);
  
  double dbinom( const int k , const int n , double p );
  double chi2_prob( double x, double df );
  bool   qchisq( double p , double df , double * );
  double dchisq( double x , double df );
  double noncentral_chi2_prob( double x, double df , double  );
  double t_prob( double x, double df );
  double ltqnorm( double p );
  double normden( double x , double m = 0 , double sd = 1 );
  
  bool t_test( double u1, double sd1, int n1 ,
	       double u2 , double s2 , int n2 ,
	       double * pvalue , double * p_lower = NULL , double * p_upper = NULL );


  // Beta distribution density probability function
  double beta( double , double , double );

  // mainly helper functions for the above
  double factln(int);
  double factrl(int);
  double gamln(double);
  double gammln(double);
  double factln(int n);
  double bico(int n, int k);
  double dbinom_raw( const double , const double , const double);
  double pythag( const double , const double );

  // matrix operations

  Data::Matrix<double> matrix_sqrt( const Data::Matrix<double> & );
  Data::Matrix<double> transpose( const Data::Matrix<double> & );
  Data::Matrix<double> inverse( const Data::Matrix<double> & , bool * flag = NULL ); // SVD-based

  Data::Matrix<double> matrix_multiply( const Data::Matrix<double> & , const Data::Matrix<double> & );
  Data::Vector<double> matrix_multiply( const Data::Matrix<double> & , const Data::Vector<double> & );  
  Data::Vector<double> matrix_multiply( const Data::Vector<double> & , const Data::Matrix<double> & );
  double matrix_inner_product( const Data::Vector<double> & , const Data::Vector<double> & );
  Data::Matrix<double> matrix_outer_product( const Data::Vector<double> & , const Data::Vector<double> & );  

  // Mean, variance and covariance

  double sum( const Data::Vector<double> & );
  double sum_squares( const Data::Vector<double> & );
  Data::Vector<double> row_sums( const Data::Matrix<double> & );
  Data::Vector<double> col_sums( const Data::Matrix<double> & );

  double mean( const Data::Vector<double> & );
  double variance( const Data::Vector<double> & );

  Data::Vector<double> mean( const Data::Matrix<double> & );
  Data::Vector<double> variance( const Data::Matrix<double> & );
  Data::Vector<double> variance( const Data::Matrix<double> & , const Data::Vector<double> & );
  
  Data::Vector<double> mean_center_cols( const Data::Matrix<double> & d );
    
  Data::Matrix<double> covariance_matrix( const Data::Matrix<double> & );
  Data::Matrix<double> covariance_matrix( const Data::Matrix<double> & , const Data::Vector<double> & );

  Data::Matrix<double> covariance_matrix( const Data::Matrix<double> & , const Data::Matrix<double> & );
  Data::Matrix<double> covariance_matrix( const Data::Matrix<double> & , const Data::Vector<double> & , 
					    const Data::Matrix<double> & , const Data::Vector<double> & );
  
  std::vector<double> canonical_correlation( const Data::Matrix<double> & , const Data::Matrix<double> & , double * pv = NULL );     

  Data::Matrix<double> cholesky( const Data::Matrix<double> & );
  
  double correlation( const std::vector<double> & a , const std::vector<double> & b );
  
  double bartlett(const int N, 
		  const int p, 
		  const int q, 
		  const std::vector<double> & eigen );
  

  struct Eigen
  {
    Eigen(const int n) : d(n), z(n,n) { }     
    Data::Vector<double> d; // eigenvalues
    Data::Matrix<double> z; // eigenvectors
  };
  
  Statistics::Eigen eigenvectors( Data::Matrix<double> & , bool * okay );
  
  bool EV_tqli( Data::Vector<double> & d , 
		Data::Vector<double> & e , 
		Data::Matrix<double> & z );

  bool EV_tred2( Data::Matrix<double> & a , 
		 Data::Vector<double> & d , 
		 Data::Vector<double> & e );
  

  Data::Vector<double> eigenvalues( Data::Matrix<double> & , bool * okay );  
  bool tqli( Data::Vector<double> & d , 
	     Data::Vector<double> & e );
	     
  
  bool tred2( Data::Matrix<double> & a , 
	      Data::Vector<double> & d , 
	      Data::Vector<double> & e );
  

}

  

#endif
