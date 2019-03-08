
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

struct mtm_t
{
  
  mtm_t( const double npi = 3 , const int nwin = 5 );
  
  void apply( const std::vector<double> * , const int fs , 
	      std::vector<double> * write_tapers = NULL , 
	      std::vector<double> * write_tapsum = NULL , 
	      std::vector<double> * write_lambda = NULL , 
	      const std::vector<double> * read_tapers = NULL , 
	      const std::vector<double> * read_tapsum = NULL , 
	      const std::vector<double> * read_lambda = NULL );

  // MTM parameters
  
  // number of pi-prolate functions (default = 3), usually something like:  3, 3.5, 4  
  double npi; 
  
  // the number of summing windows, nwin (default = 2*npi-1), e.g. 3*2-1 = 5 , i.e. drop last taper
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

  //  void bin( double w , double , double );

  // hold results here
  std::vector<double> f;
  std::vector<double> spec;
  
  //  void smooth( double w , double );
  
  // binned spectra  
/*   std::vector<double> bfa, bfb; */
/*   std::vector<double> bspec; */
  
  // output 
  bool display_tapers;
};


namespace mtm
{  

  void wrapper( edf_t & edf , param_t & param );
  
  int adwait(double *sqr_spec,
	     double *dcf,
	     double *el,
	     int nwin,
	     int num_freq,
	     double *ares,
	     double *degf,
	     double avar);

  void get_F_values(double *sr, double *si, int nf, int nwin, double *Fvalue, double *b);

  int hires(double *sqr_spec,  double *el, int nwin, int num_freq, double *ares);

  void dfour1(double data[], unsigned long nn, int isign);
  
  void jfour1(double data[], unsigned long nn, int isign);

  void jrealft(double data[], unsigned long n, int isign);

  int jtinvit_(int *nm, int *n, double *d, double *e, double *e2, 
	       int *m, double *w, int *ind,double *z, int *ierr, double *rv1, double *rv2, 
	       double *rv3, double *rv4, double *rv6);

  int jtridib_(int *n, double *eps1, double *d,
	       double *e, double *e2, double *lb,
	       double *ub, int *m11, int *m, double *w, int *ind, int *ierr, 
	       double *rv4, double *rv5);

  
  // get the multitaper slepian functions:
  // num_points = number of points in data stream
  // nwin = number of windows
  // lam = vector of eigenvalues
  // npi = order of slepian functions
  // tapsum = sum of each taper, saved for use in adaptive weighting
  // tapers =  matrix of slepian tapers, packed in a 1D double array
  
  int  multitap(int num_points, int nwin, double *lam, double npi, double *tapers, double *tapsum);
  
  //    series = input time series
  //    inum   = length of time series
  //    klength = number of elements in power spectrum (a power of 2)
  //    amp = returned power spectrum

  void  mt_get_spec(double *series, int inum, int klength, double *amp);
  
  // data = (type double) input time series
  // npoints = number of points in data

  // kind = flag for choosing hires or adaptive weighting coefficients
  // nwin = number of taper windows to calculate
  // npi = order of the slepian functions
  // inorm = flag for choice of normalization
  // dt = sampling interval (time)
  
  // ospec = output spctrum
  // dof = degrees of freedom at each frequency
  // Fvalues = Ftest value at each frequency estimate
  // klen = number of frequecies calculated (power of 2)

  void  do_mtap_spec(double *data, int npoints, int kind,
		     int nwin, double npi, int inorm, double dt,
		     double *ospec, double *dof, double *Fvalues, int klen, 
		     bool display_tapers , 
		     std::vector<double> * write_tapers = NULL, std::vector<double> * write_tapsum = NULL, std::vector<double> * write_lambda = NULL , 
		     const std::vector<double> * read_tapers = NULL , const std::vector<double> * read_tapsum = NULL, const std::vector<double> * read_lambda = NULL );
  
  
  // NR utilities
  
  void nrerror( const std::string & );
  double *vector(long nl, long nh);
  int *ivector(long nl, long nh);
  unsigned char *cvector(long nl, long nh);
  unsigned long *lvector(long nl, long nh);
  double *dvector(long nl, long nh);
  double **matrix(long nrl, long nrh, long ncl, long nch);
  double **dmatrix(long nrl, long nrh, long ncl, long nch);
  int **imatrix(long nrl, long nrh, long ncl, long nch);
  double **submatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch,
		     long newrl, long newcl);
  double **convert_matrix(double *a, long nrl, long nrh, long ncl, long nch);
  double ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
  void free_vector(double *v, long nl, long nh);
  void free_ivector(int *v, long nl, long nh);
  void free_cvector(unsigned char *v, long nl, long nh);
  void free_lvector(unsigned long *v, long nl, long nh);
  void free_dvector(double *v, long nl, long nh);
  void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);
  void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
  void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
  void free_submatrix(double **b, long nrl, long nrh, long ncl, long nch);
  void free_convert_matrix(double **b, long nrl, long nrh, long ncl, long nch);
  void free_f3tensor(double ***t, long nrl, long nrh, long ncl, long nch,long ndl, long ndh);

  // signal proc utilities

  void find_max_min(void *p, int n, double *max, double *min, int  flag);

  int get_pow_2(int inum);

  void  Db_scale(double *spec1, double *spec2, int num_freqs);

  void  Log_scale(double *spec1, double *spec2, int num_freqs);

  void Scale_Trace2(double *spec1, int num1, double *spec2, int num2, double *spec3);

  double  scale_trace_RMS(double x[], int lx);

  double  remove_mean(double x[], int lx);

  void complex_array( double input[], int inlength,double  output[], int olength);

  void zero_pad(double  output[], int start , int olength);

  void copy_trace(double a[], double b[], int num);

  void window_trace( double input[], double output[], int start, int num);

  void print_array( double array[], int inum);

  void rm_lintrend(double *x,  double *y, int length, double a, double b);

  void get_abfit(double *x, double *y, int length, double *slope, double *intercept);

  void rm_lin_sig_trend(double *y, int n, double dt, double *slope, double *cept);

  void get_indies(double t1,double t2,double dt,double tref, int numtot, int *ibeg, int *inum);

  void get_indtim(double *t1,double *t2,double dt,double tref, int numtot, int ibeg, int inum);

  double get_taper(int itype,int n, int k, double percent);
  
  int  apply_taper(double x[],int lx,int itype,double tsv);

  double get_cos_taper(int n, int k, double percent);

  void  smooth_fft(double *data, int npoints, double dt, double *naive_spec, int klen, double fWidth);
  
}


#endif
