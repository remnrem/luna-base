
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

#ifndef __LEGENDRE_H__
#define __LEGENDRE_H__

#include <vector>
#include "stats/matrix.h"

std::vector<std::vector<double> > legendre( const int , const std::vector<double> & );
std::vector<double> legendre( const int N , double x );
std::vector<Data::Matrix<double> > legendre( const int , const Data::Matrix<double> & );


char digit_to_ch ( int i );
/* int i4_log_10 ( int i ); */
/* int i4_max ( int i1, int i2 ); */
/* int i4_min ( int i1, int i2 ); */
std::string i4_to_s ( int i );
void imtqlx ( int n, double d[], double e[], double z[] );
double *p_exponential_product ( int p, double b );
double p_integral ( int n );
double *p_polynomial_coefficients ( int n );
double *p_polynomial_prime ( int m, int n, double x[] );
double *p_polynomial_prime2 ( int m, int n, double x[] );
double *p_polynomial_value ( int m, int n, double x[] );
void p_polynomial_values ( int &n_data, int &n, double &x, double &fx );
double *p_polynomial_zeros ( int nt );
double *p_power_product ( int p, int e );
void p_quadrature_rule ( int nt, double t[], double wts[] );
double *pm_polynomial_value ( int mm, int n, int m, double x[] );
void pm_polynomial_values ( int &n_data, int &n, int &m, double &x, double &fx );
double *pmn_polynomial_value ( int mm, int n, int m, double x[] );
void pmn_polynomial_values ( int &n_data, int &n, int &m, double &x, 
  double &fx );
double *pmns_polynomial_value ( int mm, int n, int m, double x[] );
void pmns_polynomial_values ( int &n_data, int &n, int &m, double &x, 
  double &fx );
double *pn_pair_product ( int p );
double *pn_polynomial_coefficients ( int n );
double *pn_polynomial_value ( int m, int n, double x[] );
void pn_polynomial_values ( int &n_data, int &n, double &x, 
  double &fx );

/* double r8_epsilon ( ); */
/* double r8_factorial ( int n ); */
/* double r8_sign ( double x ); */
/* void r8mat_print ( int m, int n, double a[], std::string title ); */
/* void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, */
/* 			int jhi, std::string title ); */
/* double r8vec_dot_product ( int n, double a1[], double a2[] ); */
/* double *r8vec_linspace_new ( int n, double a_first, double a_last ); */
/* void r8vec_print ( int n, double a[], std::string title ); */
/* void r8vec2_print ( int n, double a1[], double a2[], std::string title ); */

#endif
