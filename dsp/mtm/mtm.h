
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

#ifndef __MTM_H__
#define __MTM_H__

struct edf_t;
struct param_t;

#include <vector>
#include <string>
#include <stdint.h>

#include "stats/Eigen/Dense"
#include "fftw/fftwrap.h"


namespace mtm
{  
  void wrapper( edf_t & edf , param_t & param );  
}


struct mtm_t
{
  
  mtm_t( const double npi = 3 , const int nwin = 5 );
  
  // pre-compute tapers once for fixed size segment
  void store_tapers( const int seg_size );
  
  // do actual MT (optionally. passing pre-computed tapers in mt_tapers)
  void apply( const std::vector<double> * , const int fs ,
	      const int seg_size , const int seg_step ,
	      bool verbose = false , mtm_t * mt_tapers = NULL );
  
  //
  // MTM parameters
  //
  
  // nw: number of pi-prolate functions (default = 3), usually something like:  3, 3.5, 4  
  double npi; 
  
  // t: the number of summing windows, nwin (default = 2*npi-1), e.g. 3*2-1 = 5 , i.e. drop last taper
  int nwin;
  
  // kind of analysis
  int kind;  //  1  hires
             //  2  adwait  (default)
             // naive periodogram
  
  int inorm; // default: 4 = 1/(N.Fs) weighting
  
  // p. 335, Percival and Walden, choose 
  // npi=2,3,4 some small integer
  // W = npi/(num_points*dt);
  // or num_points*W = npi/dt ;
  // K < 2*num_points*W*dt
  //  nwin = 0...K-1

  // dB output?
  bool dB;

  // hold results here
  std::vector<double> f;
  std::vector<double> spec; // average over segments
  std::vector<double> raw_spec; // average over segments (not dB)

  std::vector<std::vector<double> > espec; // spectrogram
  std::vector<std::vector<double> > raw_espec; // spectrogram
  
  // mean center / detrend
  bool opt_remove_mean;
  bool opt_remove_trend;


  // restrict to only some segments?
  std::vector<bool> restrict;
  
  //
  // Core functions
  //

  // get the multitaper slepian functions:
  // num_points = number of points in data stream
  // nwin = number of windows
  // npi = order of slepian functions
  // generates...
  // lam = vector of eigenvalues  
  // tapsum = sum of each taper, saved for use in adaptive weighting
  // tapers =  matrix of slepian tapers, packed in a 1D double array
  
  void generate_tapers( int num_points, int nwin, double npi );

  Eigen::VectorXd lam;
  Eigen::VectorXd tapsum;
  Eigen::MatrixXd tapers;
  
  void  do_mtap_spec(real_FFT * , double *data, int npoints, int kind,
 		     int nwin, double npi, int inorm, double dt,
 		     double *ospec,
		     int klen , 
		     double *dof = NULL , double *Fvalues = NULL );
  
  int adwait(double *sqr_spec,
	     double *dcf,
	     double *el,
	     int nwin,
	     int num_freq,
	     double *ares,
	     double *degf,
	     double avar);


  //
  // Helpers (mtm-helpers.cpp)
  //
  
  //  void get_F_values(double *sr, double *si, int nf, int nwin, double *Fvalue, double *b);

  int hires(double *sqr_spec,  double *el, int nwin, int num_freq, double *ares);

  int jtinvit_(int *nm, int *n, double *d, double *e, double *e2, 
	       int *m, double *w, int *ind,double *z, int *ierr, double *rv1, double *rv2, 
	       double *rv3, double *rv4, double *rv6);
  
  int jtridib_(int *n, double *eps1, double *d,
	       double *e, double *e2, double *lb,
	       double *ub, int *m11, int *m, double *w, int *ind, int *ierr, 
	       double *rv4, double *rv5);
  int get_pow_2(int inum);  
  double remove_mean(double x[], int lx);
  void rm_lintrend(double *x,  double *y, int length, double a, double b);
  void get_abfit(double *x, double *y, int length, double *slope, double *intercept);
  void rm_lin_sig_trend(double *y, int n, double dt);

  
};




#endif
