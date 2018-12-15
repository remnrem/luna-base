
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




int mtm::jtridib_(int *n, double *eps1, double *d,
		  double *e, double *e2, double *lb,
		  double *ub, int *m11, int *m, double *w, int *ind, int *ierr, 
		  double *rv4, double *rv5)
{
  /* Initialized data */
  
  static double machep = 1.25e-15;
  
  /* System generated locals */
  int i__1, i__2;
  double d__1, d__2, d__3;
  
  /* Local variables */
  static int i, j, k, l, p, q, r, s;
  static double u, v;
  static int m1, m2;
  static double t1, t2, x0, x1;
  static int m22, ii;
  static double xu;
  static int isturm, tag;
  
  
  
  /*     this subroutine is a translation of the algol procedure bisect, */
  /*     num. math. 9, 386-393(1967) by barth, martin, and wilkinson. */
  /*     handbook for auto. comp., vol.ii-linear algebra, 249-256(1971). */
  
  /*     this subroutine finds those eigenvalues of a tridiagonal */
  /*     symmetric matrix between specified boundary indices, */
  /*     using bisection. */
  
  /*     on input: */
  
  /*        n is the order of the matrix; */
  
  /*        eps1 is an absolute error tolerance for the computed */
  /*          eigenvalues.  if the input eps1 is non-positive, */
  /*          it is reset for each submatrix to a default value, */
  /*          namely, minus the product of the relative machine */
  /*          precision and the 1-norm of the submatrix; */
  
  /*        d contains the diagonal elements of the input matrix; */
  
  /*        e contains the subdiagonal elements of the input matrix */
  /*          in its last n-1 positions.  e(1) is arbitrary; */
  
  /*        e2 contains the squares of the corresponding elements of e. */
  /*          e2(1) is arbitrary; */
  
  /*        m11 specifies the lower boundary index for the desired */
  /*          eigenvalues; */
  
  /*        m specifies the number of eigenvalues desired.  the upper */
  /*          boundary index m22 is then obtained as m22=m11+m-1. */
  
  /*     on output: */
  
  /*        eps1 is unaltered unless it has been reset to its */
  /*          (last) default value; */
  
  /*        d and e are unaltered; */
  
  /*        elements of e2, corresponding to elements of e regarded */
  /*          as negligible, have been replaced by zero causing the */
  /*          matrix to split into a direct sum of submatrices. */
  /*          e2(1) is also set to zero; */
  
  /*        lb and ub define an interval containing exactly the desired */
  /*          eigenvalues; */
  
  /*        w contains, in its first m positions, the eigenvalues */
  /*          between indices m11 and m22 in ascending order; */
  
  /*        ind contains in its first m positions the submatrix indices */
  /*          associated with the corresponding eigenvalues in w -- */
  /*          1 for eigenvalues belonging to the first submatrix from */
  /*          the top, 2 for those belonging to the second submatrix, etc.; 
   */
  
  /*        ierr is set to */
  /*          zero       for normal return, */
  /*          3*n+1      if multiple eigenvalues at index m11 make */
  /*                     unique selection impossible, */
  /*          3*n+2      if multiple eigenvalues at index m22 make */
  /*                     unique selection impossible; */
  
  /*        rv4 and rv5 are temporary storage arrays. */
  
  /*     note that subroutine tql1, imtql1, or tqlrat is generally faster */
  /*     than tridib, if more than n/4 eigenvalues are to be found. */
  
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
  --rv5;
  --rv4;
  --e2;
  --e;
  --d;
  --ind;
  --w;
  
  /* Function Body */
  
  *ierr = 0;
  tag = 0;
  xu = d[1];
  x0 = d[1];
  u = 0.;
  /*     :::::::::: look for small sub-diagonal entries and determine an */
  /*                interval containing all the eigenvalues :::::::::: */
  i__1 = *n;
  for (i = 1; i <= i__1; ++i) {
    x1 = u;
    u = 0.;
    if (i != *n) {
      u = (d__1 = e[i + 1], abs(d__1));
    }
    /* Computing MIN */
    d__1 = d[i] - (x1 + u);
    xu = min(d__1,xu);
    /* Computing MAX */
    d__1 = d[i] + (x1 + u);
    x0 = max(d__1,x0);
    if (i == 1) {
      goto L20;
    }
    if ((d__1 = e[i], abs(d__1)) > machep * ((d__2 = d[i], abs(d__2)) + (
									 d__3 = d[i - 1], abs(d__3)))) {
      goto L40;
    }
  L20:
    e2[i] = 0.;
  L40:
    ;
  }
  
  /* Computing MAX */
  d__1 = abs(xu), d__2 = abs(x0);
  x1 = max(d__1,d__2) * machep * (double) (*n);
  xu -= x1;
  t1 = xu;
  x0 += x1;
  t2 = x0;
  /*     :::::::::: determine an interval containing exactly */
  /*                the desired eigenvalues :::::::::: */
  p = 1;
  q = *n;
  m1 = *m11 - 1;
  if (m1 == 0) {
    goto L75;
  }
  isturm = 1;
 L50:
  v = x1;
  x1 = xu + (x0 - xu) * .5;
  if (x1 == v) {
    goto L980;
  }
  goto L320;
 L60:
  if ((i__1 = s - m1) < 0) {
    goto L65;
  } else if (i__1 == 0) {
    goto L73;
  } else {
    goto L70;
  }
 L65:
  xu = x1;
  goto L50;
 L70:
  x0 = x1;
  goto L50;
 L73:
  xu = x1;
  t1 = x1;
 L75:
  m22 = m1 + *m;
  if (m22 == *n) {
    goto L90;
  }
  x0 = t2;
  isturm = 2;
  goto L50;
 L80:
  if ((i__1 = s - m22) < 0) {
    goto L65;
  } else if (i__1 == 0) {
    goto L85;
  } else {
    goto L70;
  }
 L85:
  t2 = x1;
 L90:
  q = 0;
  r = 0;
  /*     :::::::::: establish and process next submatrix, refining */
  /*                interval by the gerschgorin bounds :::::::::: */
 L100:
  if (r == *m) {
    goto L1001;
  }
  ++tag;
  p = q + 1;
  xu = d[p];
  x0 = d[p];
  u = 0.;
  
  i__1 = *n;
  for (q = p; q <= i__1; ++q) {
    x1 = u;
    u = 0.;
    v = 0.;
    if (q == *n) {
      goto L110;
    }
    u = (d__1 = e[q + 1], abs(d__1));
    v = e2[q + 1];
  L110:
    /* Computing MIN */
    d__1 = d[q] - (x1 + u);
    xu = min(d__1,xu);
    /* Computing MAX */
    d__1 = d[q] + (x1 + u);
    x0 = max(d__1,x0);
    if (v == 0.) {
      goto L140;
    }
    /* L120: */
  }
  
 L140:
  /* Computing MAX */
  d__1 = abs(xu), d__2 = abs(x0);
  x1 = max(d__1,d__2) * machep;
  if (*eps1 <= 0.) {
    *eps1 = -x1;
  }
  if (p != q) {
    goto L180;
  }
  /*     :::::::::: check for isolated root within interval :::::::::: */
  if (t1 > d[p] || d[p] >= t2) {
    goto L940;
  }
  m1 = p;
  m2 = p;
  rv5[p] = d[p];
  goto L900;
 L180:
  x1 *= (double) (q - p + 1);
  /* Computing MAX */
  d__1 = t1, d__2 = xu - x1;
  *lb = max(d__1,d__2);
  /* Computing MIN */
  d__1 = t2, d__2 = x0 + x1;
  *ub = min(d__1,d__2);
  x1 = *lb;
  isturm = 3;
  goto L320;
 L200:
  m1 = s + 1;
  x1 = *ub;
  isturm = 4;
  goto L320;
 L220:
  m2 = s;
  if (m1 > m2) {
    goto L940;
  }
  /*     :::::::::: find roots by bisection :::::::::: */
  x0 = *ub;
  isturm = 5;
  
  i__1 = m2;
  for (i = m1; i <= i__1; ++i) {
    rv5[i] = *ub;
    rv4[i] = *lb;
    /* L240: */
  }
  /*     :::::::::: loop for k-th eigenvalue */
  /*                for k=m2 step -1 until m1 do -- */
  /*                (-do- not used to legalize -computed go to-) :::::::::: 
   */
  k = m2;
 L250:
  xu = *lb;
  /*     :::::::::: for i=k step -1 until m1 do -- :::::::::: */
  i__1 = k;
  for (ii = m1; ii <= i__1; ++ii) {
    i = m1 + k - ii;
    if (xu >= rv4[i]) {
      goto L260;
    }
    xu = rv4[i];
    goto L280;
  L260:
    ;
  }
  
 L280:
  if (x0 > rv5[k]) {
    x0 = rv5[k];
  }
  /*     :::::::::: next bisection step :::::::::: */
 L300:
  x1 = (xu + x0) * .5;
  if (x0 - xu <= machep * 2. * (abs(xu) + abs(x0)) + abs(*eps1)) {
    goto L420;
  }
  /*     :::::::::: in-line procedure for sturm sequence :::::::::: */
 L320:
  s = p - 1;
  u = 1.;
  
  i__1 = q;
  for (i = p; i <= i__1; ++i) {
    if (u != 0.) {
      goto L325;
    }
    v = (d__1 = e[i], abs(d__1)) / machep;
    if (e2[i] == 0.) {
      v = 0.;
    }
    goto L330;
  L325:
    v = e2[i] / u;
  L330:
    u = d[i] - x1 - v;
    if (u < 0.) {
      ++s;
    }
    /* L340: */
  }
  
  switch ((int)isturm) {
  case 1:  goto L60;
  case 2:  goto L80;
  case 3:  goto L200;
  case 4:  goto L220;
  case 5:  goto L360;
  }
  /*     :::::::::: refine intervals :::::::::: */
 L360:
  if (s >= k) {
    goto L400;
  }
  xu = x1;
  if (s >= m1) {
    goto L380;
  }
  rv4[m1] = x1;
  goto L300;
 L380:
  rv4[s + 1] = x1;
  if (rv5[s] > x1) {
    rv5[s] = x1;
  }
  goto L300;
 L400:
  x0 = x1;
  goto L300;
  /*     :::::::::: k-th eigenvalue found :::::::::: */
 L420:
  rv5[k] = x1;
  --k;
  if (k >= m1) {
    goto L250;
  }
  /*     :::::::::: order eigenvalues tagged with their */
  /*                submatrix associations :::::::::: */
 L900:
  s = r;
  r = r + m2 - m1 + 1;
  j = 1;
  k = m1;
  
  i__1 = r;
  for (l = 1; l <= i__1; ++l) {
    if (j > s) {
      goto L910;
    }
    if (k > m2) {
      goto L940;
    }
    if (rv5[k] >= w[l]) {
      goto L915;
    }
    
    i__2 = s;
    for (ii = j; ii <= i__2; ++ii) {
      i = l + s - ii;
      w[i + 1] = w[i];
      ind[i + 1] = ind[i];
      /* L905: */
    }
    
  L910:
    w[l] = rv5[k];
    ind[l] = tag;
    ++k;
    goto L920;
  L915:
    ++j;
  L920:
    ;
  }
  
 L940:
    if (q < *n) {
      goto L100;
    }
    goto L1001;
    /*     :::::::::: set error -- interval cannot be found containing */
    /*                exactly the desired eigenvalues :::::::::: */
 L980:
    *ierr = *n * 3 + isturm;
 L1001:
    *lb = t1;
    *ub = t2;
    return 0;
    /*     :::::::::: last card of tridib :::::::::: */
} /* tridib_ */

