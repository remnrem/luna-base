
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

// All mtm functions are adapted code from: Lees, J. M. and J. Park
// (1995): Multiple-taper spectral analysis: A stand-alone
// C-subroutine: Computers & Geology: 21, 199-236.

#include "mtm.h"

#include <iostream>

#include <cstdio>
#include <cmath>

#define PI 3.14159265358979

#define ABS(a) ((a) < (0) ? -(a) : (a))

#include "nrutil.h"

#define DIAG1 0
#define MAX(a,b) ((a) >= (b) ? (a) : (b))


// int jtridib_(int *n, double *eps1, double *d, double *e, double *e2,
// 	     double *lb, double *ub, int *m11, int *m, double *w, int *ind, int *ierr, 
// 	     double *rv4, double *rv5);

// int jtinvit_(int *nm, int *n, double *d, double *e, double *e2, 
// 	     int *m, double *w, int *ind,double *z, int *ierr, double *rv1, double *rv2, 
// 	      double *rv3, double *rv4, double *rv6);



int  mtm::multitap(int num_points, int nwin, double *lam, double npi, double *tapers, double *tapsum)
{
  /*
   * get the multitaper slepian functions: 
   num_points = number of points in data stream
   nwin = number of windows 
   lam= vector of eigenvalues 
   npi = order of
   slepian functions 
   tapsum = sum of each taper, saved for use in adaptive weighting  
   tapers =  matrix of slepian tapers, packed in a 1D double array
  */
  
  int             i, j, k, kk;
  double         *z, ww, cs, ai, an, eps, rlu, rlb, aa;
  double          dfac, drat, gamma, bh, tapsq, TWOPI, DPI;
  double         *diag, *offdiag, *offsq;
  //char           *k1[4];
  
  char            name[81];
  double         *scratch1, *scratch2, *scratch3, *scratch4, *scratch6;
  
  
  /* need to initialize iwflag = 0 */
  double          anpi;
  double         *ell;
  int             key, nbin, npad;
  int            *ip;
  double         *evecs;
  //	double         *zee;
  
  long            len;
  int             ierr;
  int             m11;
  DPI = (double) PI;
  TWOPI = (double) 2 *DPI;
  
  anpi = npi;
  an = (double) (num_points);
  ww = (double) (anpi) / an;	/* this corresponds to P&W's W value  */
  cs = cos(TWOPI * ww);
  
  //  std::cerr << "P&W W value = " << ww << "\n";

  ell = dvector((long) 0, (long) nwin);
  
  diag = dvector((long) 0, (long) num_points);
  offdiag = dvector((long) 0, (long) num_points);
  offsq = dvector((long) 0, (long) num_points);
  
  scratch1 = dvector((long) 0, (long) num_points);
  scratch2 = dvector((long) 0, (long) num_points);
  scratch3 = dvector((long) 0, (long) num_points);
  scratch4 = dvector((long) 0, (long) num_points);
  scratch6 = dvector((long) 0, (long) num_points);
  


  /* make the diagonal elements of the tridiag matrix  */
  
  for (i = 0; i < num_points; i++) {
    ai = (double) (i);
    diag[i] = -cs * (((an - 1.) / 2. - ai)) * (((an - 1.) / 2. - ai));
    offdiag[i] = -ai * (an - ai) / 2.;
    offsq[i] = offdiag[i] * offdiag[i];
  }
  
  eps = 1.0e-13;
  m11 = 1;
  

  ip = ivector((long) 0, (long) nwin);

  /* call the eispac routines to invert the tridiagonal system */

  jtridib_(&num_points, &eps, diag, offdiag, offsq, &rlb, &rlu, &m11, &nwin, lam,
	   ip, &ierr, scratch1, scratch2);
#if DIAG1
  fprintf(stderr, "ierr=%d rlb=%.8f rlu=%.8f\n", ierr, rlb, rlu);
  
  fprintf(stderr, "eigenvalues for the eigentapers\n");
  
  for (k = 0; k < nwin; k++)
    fprintf(stderr, "%.20f ", lam[k]);
  
#endif
  

  len = num_points * nwin;
  evecs = dvector((long) 0, (long) len);
  
  
  jtinvit_(&num_points, &num_points, diag, offdiag, offsq, &nwin, lam, ip, evecs, &ierr,
	   scratch1, scratch2, scratch3, scratch4, scratch6);
  


  
  free_dvector(scratch1, (long) 0, (long) num_points);
  free_dvector(scratch2, (long) 0, (long) num_points);
  free_dvector(scratch3, (long) 0, (long) num_points);
  free_dvector(scratch4, (long) 0, (long) num_points);
  free_dvector(scratch6, (long) 0, (long) num_points);
  


  /*
   * we calculate the eigenvalues of the dirichlet-kernel problem i.e.
   * the bandwidth retention factors from slepian 1978 asymptotic
   * formula, gotten from thomson 1982 eq 2.5 supplemented by the
   * asymptotic formula for k near 2n from slepian 1978 eq 61 more
   * precise values of these parameters, perhaps useful in adaptive
   * spectral estimation, can be calculated explicitly using the
   * rayleigh-quotient formulas in thomson (1982) and park et al (1987)
   * 
   */
  
  dfac = (double) an *DPI * ww;
  drat = (double) 8. *dfac;
  
  
  dfac = (double) 4. *sqrt(DPI * dfac) * exp((double) (-2.0) * dfac);
  
  
  for (k = 0; k < nwin; k++) {
    lam[k] = (double) 1.0 - (double) dfac;
    dfac = dfac * drat / (double) (k + 1);
    
    
    
    /* fails as k -> 2n */
  }
  
  
  gamma = log((double) 8. * an * sin((double) 2. * DPI * ww)) + (double) 0.5772156649;
  


  for (k = 0; k < nwin; k++) {
    bh = -2. * DPI * (an * ww - (double) (k) /
		      (double) 2. - (double) .25) / gamma;
    ell[k] = (double) 1. / ((double) 1. + exp(DPI * (double) bh));
    
  }
  
  for (i = 0; i < nwin; i++)
    lam[i] = MAX(ell[i], lam[i]);
  

  //   normalize the eigentapers to preserve power for a white process
  //   i.e. they have rms value unity
  //   tapsum is the average of the eigentaper, should be near zero for
  //   antisymmetric tapers

  for (k = 0; k < nwin; k++) {
    kk = (k) * num_points;
    tapsum[k] = 0.;
    tapsq = 0.;
    for (i = 0; i < num_points; i++) {
      aa = evecs[i + kk];
      tapers[i + kk] = aa;
      tapsum[k] = tapsum[k] + aa;
      tapsq = tapsq + aa * aa;
    }
    aa = sqrt(tapsq / (double) num_points);
    tapsum[k] = tapsum[k] / aa;
    
    for (i = 0; i < num_points; i++) {
      tapers[i + kk] = tapers[i + kk] / aa;
      
    }
  }
  
  
  /* Free Memory */
  
  
  free_dvector(ell, (long) 0, (long) nwin);
  free_dvector(diag, (long) 0, (long) num_points);
  free_dvector(offdiag, (long) 0, (long) num_points);
  free_dvector(offsq, (long) 0, (long) num_points);
  free_ivector(ip, (long) 0, (long) nwin);
  
  free_dvector(evecs, (long) 0, (long) len);
  
  
  return 1;
}

