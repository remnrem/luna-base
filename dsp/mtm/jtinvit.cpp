
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

#include "mtm.h"

#include <cmath>

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (double)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (double)min(a,b)
#define dmax(a,b) (double)max(a,b)



int mtm::jtinvit_(int *nm, int *n, double *d, double *e, double *e2, 
		  int *m, double *w, int *ind,double *z, int *ierr, double *rv1, double *rv2, 
		  double *rv3, double *rv4, double *rv6)
{
  /* Initialized data */
  
  static double machep = 1.25e-15;
  
  /* System generated locals */
  int z_dim1, z_offset, i__1, i__2, i__3;
  double d__1, d__2;
  
  /* Builtin functions */
  //    double sqrt();
  
  /* Local variables */
  static double norm;
  static int i, j, p, q, r, s;
  static double u, v, order;
  static int group;
  static double x0, x1;
  static int ii, jj, ip;
  static double uk, xu;
  static int tag, its;
  static double eps2, eps3, eps4;
  
  static double rtem;

  
  /*     this subroutine is a translation of the inverse iteration tech- */
  /*     nique in the algol procedure tristurm by peters and wilkinson. */
  /*     handbook for auto. comp., vol.ii-linear algebra, 418-439(1971). */
  
  /*     this subroutine finds those eigenvectors of a tridiagonal */
  /*     symmetric matrix corresponding to specified eigenvalues, */
  /*     using inverse iteration. */
  
  /*     on input: */
  
  /*        nm must be set to the row dimension of two-dimensional */
  /*          array parameters as declared in the calling program */
  /*          dimension statement; */
  
  /*        n is the order of the matrix; */
  
  /*        d contains the diagonal elements of the input matrix; */
  
  /*        e contains the subdiagonal elements of the input matrix */
  /*          in its last n-1 positions.  e(1) is arbitrary; */
  
  /*        e2 contains the squares of the corresponding elements of e, */
  /*          with zeros corresponding to negligible elements of e. */
  /*          e(i) is considered negligible if it is not larger than */
  /*          the product of the relative machine precision and the sum */
  /*          of the magnitudes of d(i) and d(i-1).  e2(1) must contain */
  /*          0.0d0 if the eigenvalues are in ascending order, or 2.0d0 */
  /*          if the eigenvalues are in descending order.  if  bisect, */
  /*          tridib, or  imtqlv  has been used to find the eigenvalues, */
  /*          their output e2 array is exactly what is expected here; */
  
  /*        m is the number of specified eigenvalues; */
  
  /*        w contains the m eigenvalues in ascending or descending order; 
   */
  
  /*        ind contains in its first m positions the submatrix indices */
  /*          associated with the corresponding eigenvalues in w -- */
  /*          1 for eigenvalues belonging to the first submatrix from */
  /*          the top, 2 for those belonging to the second submatrix, etc. 
   */
  
  /*     on output: */
  
  /*        all input arrays are unaltered; */
  
  /*        z contains the associated set of orthonormal eigenvectors. */
  /*          any vector which fails to converge is set to zero; */
  
  /*        ierr is set to */
  /*          zero       for normal return, */
  /*          -r         if the eigenvector corresponding to the r-th */
  /*                     eigenvalue fails to converge in 5 iterations; */
  
  /*        rv1, rv2, rv3, rv4, and rv6 are temporary storage arrays. */
  
  /*     questions and comments should be directed to b. s. garbow, */
  /*     applied mathematics division, argonne national laboratory */
  
  /*     ------------------------------------------------------------------ 
   */
  
  /*     :::::::::: machep is a machine dependent parameter specifying */
  /*                the relative precision of floating point arithmetic. */
  /*                machep = 16.0d0**(-13) for long form arithmetic */
  /*                on s360 :::::::::: */
  /*  for f_floating dec fortran */
  /*      data machep/1.1d-16/ */
  /*  for g_floating dec fortran */
  /* Parameter adjustments */
  --rv6;
  --rv4;
  --rv3;
  --rv2;
  --rv1;
  --e2;
  --e;
  --d;
  z_dim1 = *nm;
  z_offset = z_dim1 + 1;
  z -= z_offset;
  --ind;
  --w;
  
  /* Function Body */
  
  *ierr = 0;
  if (*m == 0) {
    goto L1001;
  }
  tag = 0;
  order = 1. - e2[1];
  q = 0;
  /*     :::::::::: establish and process next submatrix :::::::::: */
 L100:
  p = q + 1;
  
  i__1 = *n;
  for (q = p; q <= i__1; ++q) {
    if (q == *n) {
      goto L140;
    }
    if (e2[q + 1] == 0.) {
      goto L140;
    }
    /* L120: */
  }
  /*     :::::::::: find vectors by inverse iteration :::::::::: */
 L140:
  ++tag;
  s = 0;
  
  i__1 = *m;
  for (r = 1; r <= i__1; ++r) {
    if (ind[r] != tag) {
      goto L920;
    }
    its = 1;
    x1 = w[r];
    if (s != 0) {
      goto L510;
    }
    /*     :::::::::: check for isolated root :::::::::: */
    xu = 1.;
    if (p != q) {
      goto L490;
    }
    rv6[p] = 1.;
    goto L870;
  L490:
    norm = (d__1 = d[p], abs(d__1));
    ip = p + 1;
    
    i__2 = q;
    for (i = ip; i <= i__2; ++i) {
      /* L500: */
      norm = norm + (d__1 = d[i], abs(d__1)) + (d__2 = e[i], abs(d__2));
    }
    /*     :::::::::: eps2 is the criterion for grouping, */
    /*                eps3 replaces zero pivots and equal */
    /*                roots are modified by eps3, */
    /*                eps4 is taken very small to avoid overflow :::::::::
		      : */
    eps2 = norm * .001;
    eps3 = machep * norm;
    uk = (double) (q - p + 1);
    eps4 = uk * eps3;
    uk = eps4 / sqrt(uk);
    s = p;
  L505:
    group = 0;
    goto L520;
    /*     :::::::::: look for close or coincident roots :::::::::: */
  L510:
    if ((d__1 = x1 - x0, abs(d__1)) >= eps2) {
      goto L505;
    }
    ++group;
    if (order * (x1 - x0) <= 0.) {
      x1 = x0 + order * eps3;
    }
    /*     :::::::::: elimination with interchanges and */
    /*                initialization of vector :::::::::: */
  L520:
    v = 0.;
    
    i__2 = q;
    for (i = p; i <= i__2; ++i) {
      rv6[i] = uk;
      if (i == p) {
	goto L560;
      }
      if ((d__1 = e[i], abs(d__1)) < abs(u)) {
	goto L540;
      }
      /*     :::::::::: warning -- a divide check may occur here if */
      /*                e2 array has not been specified correctly ::::::
			:::: */
      xu = u / e[i];
      rv4[i] = xu;
      rv1[i - 1] = e[i];
      rv2[i - 1] = d[i] - x1;
      rv3[i - 1] = 0.;
      if (i != q) {
	rv3[i - 1] = e[i + 1];
      }
      u = v - xu * rv2[i - 1];
      v = -xu * rv3[i - 1];
      goto L580;
    L540:
      xu = e[i] / u;
      rv4[i] = xu;
      rv1[i - 1] = u;
      rv2[i - 1] = v;
      rv3[i - 1] = 0.;
    L560:
      u = d[i] - x1 - xu * v;
      if (i != q) {
	v = e[i + 1];
      }
    L580:
      ;
    }
    
    if (u == 0.) {
      u = eps3;
    }
    rv1[q] = u;
    rv2[q] = 0.;
    rv3[q] = 0.;
    /*     :::::::::: back substitution */
    /*                for i=q step -1 until p do -- :::::::::: */
  L600:
    i__2 = q;
    for (ii = p; ii <= i__2; ++ii) {
      i = p + q - ii;
      rtem = rv6[i] - u * rv2[i] - v * rv3[i];
      rv6[i] = (rtem) / rv1[i];
      v = u;
      u = rv6[i];
      /* L620: */
    }
    /*     :::::::::: orthogonalize with respect to previous */
    /*                members of group :::::::::: */
    if (group == 0) {
      goto L700;
    }
    j = r;
    
    i__2 = group;
    for (jj = 1; jj <= i__2; ++jj) {
    L630:
      --j;
      if (ind[j] != tag) {
	goto L630;
      }
      xu = 0.;
      
      i__3 = q;
      for (i = p; i <= i__3; ++i) {
	/* L640: */
	xu += rv6[i] * z[i + j * z_dim1];
      }
      
      i__3 = q;
      for (i = p; i <= i__3; ++i) {
	/* L660: */
	rv6[i] -= xu * z[i + j * z_dim1];
      }
      
      /* L680: */
    }
    
  L700:
    norm = 0.;
    
    i__2 = q;
    for (i = p; i <= i__2; ++i) {
      /* L720: */
      norm += (d__1 = rv6[i], abs(d__1));
    }
    
    if (norm >= 1.) {
      goto L840;
    }
    /*     :::::::::: forward substitution :::::::::: */
    if (its == 5) {
      goto L830;
    }
    if (norm != 0.) {
      goto L740;
    }
    rv6[s] = eps4;
    ++s;
    if (s > q) {
      s = p;
    }
    goto L780;
  L740:
    xu = eps4 / norm;
    
    i__2 = q;
    for (i = p; i <= i__2; ++i) {
      /* L760: */
      rv6[i] *= xu;
    }
    /*     :::::::::: elimination operations on next vector */
    /*                iterate :::::::::: */
  L780:
    i__2 = q;
    for (i = ip; i <= i__2; ++i) {
      u = rv6[i];
      /*     :::::::::: if rv1(i-1) .eq. e(i), a row interchange */
      /*                was performed earlier in the */
      /*                triangularization process :::::::::: */
      if (rv1[i - 1] != e[i]) {
	goto L800;
      }
      u = rv6[i - 1];
      rv6[i - 1] = rv6[i];
    L800:
      rv6[i] = u - rv4[i] * rv6[i - 1];
      /* L820: */
    }
    
    ++its;
    goto L600;
    /*     :::::::::: set error -- non-converged eigenvector :::::::::: */
  L830:
    *ierr = -r;
    xu = 0.;
    goto L870;
    /*     :::::::::: normalize so that sum of squares is */
    /*                1 and expand to full order :::::::::: */
  L840:
    u = 0.;
    
    i__2 = q;
    for (i = p; i <= i__2; ++i) {
      /* L860: */
      /* Computing 2nd power */
      d__1 = rv6[i];
      u += d__1 * d__1;
    }
    
    xu = 1. / sqrt(u);
    
  L870:
    i__2 = *n;
    for (i = 1; i <= i__2; ++i) {
      /* L880: */
      z[i + r * z_dim1] = 0.;
    }
    
    i__2 = q;
    for (i = p; i <= i__2; ++i) {
      /* L900: */
      z[i + r * z_dim1] = rv6[i] * xu;
    }
    
    x0 = x1;
  L920:
    ;
  }
  
  if (q < *n) {
    goto L100;
  }
 L1001:
  return 0;
  /*     :::::::::: last card of tinvit :::::::::: */
} /* tinvit_ */

