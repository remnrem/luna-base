
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
#include <cstdlib>

#include <string>

#if defined(__MACH__)
#include <stdlib.h>
#else 
#include <malloc.h>
#endif

#define PI 3.141592654
#define ABS(a) ((a) < (0) ? -(a) : (a))

#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#define MIN(A, B) ((A) < (B) ? (A) : (B))


#include "nrutil.h"


//void jrealft(double data[], unsigned long n, int isign);

void mtm::find_max_min(void *p, int n, double *max, double *min, int  flag)
 {

   int i;
   int * pi = (int*)p;
   double *pf=(double*)p;
  
  if(flag){
    (*max) = (*min) = *pf;
    for(i=1; i<n; i++){
      *max = MAX(*max, *(pf+i));
      *min = MIN(*min, *(pf+i));
    }
  }
  else{
    (*max) = (*min) = *pi;
    for(i=1; i<n; i++){
      *max = MAX(*max, *(pi+i));
      *min = MIN(*min, *(pi+i));
    }
  }
}

/**************************************************************/

int mtm::get_pow_2(int inum)
{
  int             j, klength;
  /* find smallest power of 2 that encompasses the data */
  
  for (j = 1; pow((double) 2, (double) j) < inum; j++);
  return klength = pow((double) 2, (double) j);
}



void  mtm::Db_scale(double *spec1, double *spec2, int num_freqs)
{
  /* make a copy of spec2 onto spec1 with the Db scale  */
  int             i;
  for (i = 0; i < num_freqs; i++) {
    if (spec2[i] <= 0.0) {
      fprintf(stderr, "negative or zero spectrum: %d\n", i);
      fprintf(stderr, "%g \n", spec2[i]);
      exit(0);
    }
    spec1[i] = 10. * log10((double) spec2[i]);
  }
}

void  mtm::Log_scale(double *spec1, double *spec2, int num_freqs)
{
  /* make a copy of spec2 onto spec1 with the Db scale  */
  int             i;
  for (i = 0; i < num_freqs; i++) {
    if (spec2[i] <= 0.0) {spec1[i]=0.0;
      fprintf(stderr, "negative or zero spectrum: %d %g \n", i, spec2[i]);
    }
    spec1[i] =  log10( spec2[i]);
  }
}



void mtm::Scale_Trace2(double *spec1, int num1, double *spec2, int num2, double *spec3)
{
  /*
   * this routine scales one vector to the other for plotting purposes
   */
  double           diff1, diff2, max1, min1, max2, min2;
  int             double_or_int = 1;
  int             i;
  find_max_min(spec1, num1, &max1, &min1, double_or_int);
  find_max_min(spec2, num2, &max2, &min2, double_or_int);
  diff1 = max2 - min2;
  diff2 = max1 - min1;
  /* fprintf(stderr, "Scale2: max1=%e, min1=%e ,max2=%e ,min2=%e\n", max1, min1, max2, min2); */
  for (i = 0; i < num2; i++) {
    spec3[i] = ((spec2[i] - min2) / diff1) * diff2 + min1;
  }
  
}



double  mtm::scale_trace_RMS(double x[], int lx)
{
  
  int k;
  double mean;
  double std;
  mean = 0.;
  if(lx < 2 ) return mean;
  
  for( k=0; k<lx; k++)
    {
      mean += x[k];
    }
  
  mean = mean/ (double)lx;
  
  for( k=0; k<lx; k++)
    {
      x[k] = x[k] - mean;
    }
  
  std = 0.;
  for( k=0; k<lx; k++)
    {
      std += x[k] * x[k];
    }
  std = sqrt(std)/(lx-1);
  
  for( k=0; k<lx; k++)
    {
      x[k] = x[k] / std ;
    }
  
  /*  fprintf(stderr," %e %e \n", mean, std);*/
  return  mean;
}
/*********************************************************/

double  mtm::remove_mean(double x[], int lx)
{
  
  double mean = 0.;
  if( lx < 2 ) return mean;
  
  for( int k=0; k<lx; k++)
    {
      mean += x[k];
    }
  
  mean = mean/ (double)lx;
  
  for( int k=0; k<lx; k++)
    {
      x[k] = x[k] - mean;
    }

  return  mean;
}


void mtm::complex_array( double input[], int inlength,double  output[], int olength)
{
  int i,k,j;
  
  for(i=0; i< inlength ;  i++)
    {
      k = 2*i;
      j = k+1;
      
      if(i>olength) break;
      output[k] = input[i];
      
      output[j] = 0.0;
    }
  
}

void mtm::zero_pad(double  output[], int start , int olength)
{
  int i;
  for( i= start ; i< olength; i++) 
    {   
      output[i] = 0.0; 
    }
}


void mtm::copy_trace(double a[], double b[], int num)
{
  /* copies vector a ONTO vector b  */
  int i;
  for(i=0; i< num ;  i++)
    {
      b[i] = a[i];
    } 
}


void mtm::window_trace( double input[], double output[], int start, int num)
{
  int i;
  for(i=0; i< num ;  i++)
    {
      output[i] = input[i+start];
    }
}


void mtm::print_array( double array[], int inum)
{
  int i;
  for(i=0; i<inum; i++)
    {
      printf("%d %f\n", i, array[i]);
    }
}


void mtm::rm_lintrend(double *x,  double *y, int length, double a, double b)
{   /* program removes a linear trend of the data in vector y
       returns a de-trended vector y  */
  int i;
  /* fprintf(stderr, "in rm_lintrend %f %f\n",a,b); */
  for (i = 0; i < length; i++)
    
    { y[i] = y[i] - x[i]*a - b; }
  
  /*  fprintf(stderr, "done in rm trend....\n"); */
}


void mtm::get_abfit(double *x, double *y, int length, double *slope, double *intercept)
{   
  
  double s=0.0, sx=0.0, sy=0.0, sxx=0.0, sxy=0.0;
  double del;
  int i;
  
  for (i = 0; i < length; i++)    
    {
      
      sx += x[i];
      sy += y[i];
      sxx += x[i] * x[i];
      sxy += x[i]*y[i];
    } 
  
  s = length;
  
  del = s*sxx - sx*sx;
  
  if(del != 0.0)
    {
      *intercept = (sxx*sy - sx*sxy)/del;
      *slope = (s*sxy - sx*sy)/del;
    }
  
}


void mtm::rm_lin_sig_trend(double *y, int n, double dt, double *slope, double *cept)
{   /* program removes a linear trend of a time series vector y
       returns a de-trended vector y  */
  int i;
  double *x;
  double a, b;
  
  /* create an x vector of time */
  
  fprintf(stderr, "starting rm_lin_sig_trend....\n");
  
  x = (double *)malloc( ((n)*sizeof(double)) ); 
  for(i=0; i<n; i++)
    { 
      x[i] = i*dt;
    }
  
  /*  find the line of the data  */
  get_abfit(x, y, n, &a, &b);
  
  /*   remove the trend from the data  */
  rm_lintrend(x, y, n, a, b);
  fprintf(stderr, "fixing slope and cept....\n");
  
#if 0
  *slope = a;
  *cept  = b;
  
  
  fprintf(stderr, "done in rm_lin_sig_trend %f %f \n", *slope, *cept);
#endif
  free(x);
  
}



void mtm::get_indies(double t1,double t2,double dt,double tref, int numtot, int *ibeg, int *inum)
{
  
  *inum = (int)((t2 - t1)/dt )+1 ;  
  
  *ibeg = (int)((t1 - tref)/dt);
  
  if( *ibeg < 0 ) *ibeg = 0;
  if( (*ibeg + *inum)  >  numtot ) *inum = numtot- *ibeg;
  
}


void mtm::get_indtim(double *t1,double *t2,double dt,double tref, int numtot, int ibeg, int inum)
{
  
  *t1 = (double)ibeg*dt+tref;
  *t2 = (inum-1)*dt+ *t1;
  
}


double mtm::get_taper(int itype,int n, int k, double percent)
{

  /*
    c-this function generates a single sample of a data window.
    c-itype=1(rectangular), 2(tapered rectangular), 3(triangular),
    c-      4(hanning), 5(hamming), or 6(blackman).
    c-      (note:  tapered rectangular has cosine-tapered 10% ends.)
    c-n=size (total no. samples) of window.
    c-k=sample number within window, from 0 through n-1.
    c-  (if k is outside this range, spwndo is set to 0.)
  */
  int l;
  double vwin;
  vwin = 0.0;
  if(itype < 1 || itype > 6) return vwin;
  if(k<0 || k > n) return vwin;
  vwin = 1.0;
  switch(itype)
    {
    case 1:
      return vwin;
      break;
    case 2:
      l=(n-2)*percent;
      if(k<=l) vwin=0.5*(1.0-cos(k*PI/(l+1)));
      if(k>=n-l-2) vwin=0.5*(1.0-cos((n-k-1)*PI/(l+1)));
      return vwin;
      break;
    case 3:
      vwin=1.0-ABS(1.0-2*k/(n-1.0));
      return vwin;
      break;
    case 4:
      vwin=0.5*(1.0-cos(2*k*PI/(n-1))); 
      return vwin;    
      break;
    case 5:
      vwin=0.54-0.46*cos(2*k*PI/(n-1));
      return vwin;
      break;
    case 6:
      vwin=0.42-0.5*cos(2*k*PI/(n-1))+0.08*cos(4*k*PI/(n-1));
      return vwin;
      break;
    }
  return vwin;
}


int  mtm::apply_taper(double x[],int lx,int itype,double tsv)
{
  int k, ierror;
  double w;
  ierror=1;
  if(itype<1  || itype > 6) return ierror;
  tsv=0.0;
  for( k=0; k<=lx; k++)
    {
      w=get_taper(itype,lx+1,k, 0.05);
      x[k]=x[k]*w;
      tsv=tsv+w*w;
    }
  ierror=0;
  return ierror;
}

double mtm::get_cos_taper(int n, int k, double percent)
{
  /* n = number of data points in window. k = index of data point  */
  int l;
  double vwin;
  vwin = 0.0;
  
  if(k<0 || k > n) return vwin;
  vwin = 1.0;
  
  l=(n-2)*percent;
  if(k<=l) vwin=0.5*(1.0-cos(k*PI/(l+1)));
  if(k>=n-l-2) vwin=0.5*(1.0-cos((n-k-1)*PI/(l+1)));
  
  return vwin;
}


void  mtm::smooth_fft(double *data, int npoints, double dt, double *naive_spec, int klen, double fWidth)
{
  int             isign = 1;
  double          *dtemp, df, freqwin, tem;
  int             num_freqs;
  int             i, k, j;
  num_freqs = 1 + klen / 2;
  
  dtemp = vector(0, klen);
  
  copy_trace(data, dtemp, npoints);
  
  
  zero_pad(dtemp, npoints, klen);
  jrealft(dtemp - 1, (unsigned long) klen, isign);
   
  for (i = 1; i < num_freqs - 1; i++) {
    naive_spec[i] = (SQR(dtemp[2 * i + 1]) + SQR(dtemp[2 * i]));
    
  }
  naive_spec[0] = SQR(fabs(dtemp[0]));
  naive_spec[num_freqs - 1] = SQR(fabs(dtemp[1]));
  
  df = 2 * (0.5 / dt) / klen;
  freqwin = (int) (fWidth / df) / 2;
  
#if 1
  
  /* smooth the periodogram   
     fprintf(stderr, "smooth the periodogram 4, freqwin=%d\n", freqwin);
  */
  

  for (i = 0; i < num_freqs; i++) {
    tem = 0.0;
    k = 0;
    for (j = i - freqwin; j <= i + freqwin; j++) {
      if (j > 0 && j < num_freqs - 1) {
	tem += naive_spec[j];
	k++;
      }
    }
    
    if (k > 0) {
      naive_spec[i] = tem / (double) k;
    } else
      naive_spec[i] = naive_spec[i]; 
    
  }
  
#endif
  free_vector(dtemp, 0 , klen);
}

