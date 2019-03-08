
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
#include "helper/helper.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

#define ABS(a) ((a) < (0) ? -(a) : (a))

#include "nrutil.h"

#define perr(x,y)  (fprintf(stderr, x , y))
#define prbl (fprintf(stderr,"\n"))


void  mtm::mt_get_spec(double *series, int inum, int klength, double *amp)
{
  
  /*    series = input time series
	inum   = length of time series
	klength = number of elements in power spectrum (a power of 2)
	amp = returned power spectrum
  */

  int             i, j, isign = 1;
  
  unsigned long   nn;

  double           tsv;
    
  nn = klength;
  
  /* copy amp onto series and apply zero padding to  klength */
  
  for (i = 0; i < inum; i++) { amp[i] = series[i];  }

  zero_pad(amp, inum, klength);
  
  
  /*  Fast Fourier Transform Routine:  here we are using the Numerical Recipes
      routine jrealft which returns the fft in the 1-D input array
      packed as pairs of real numbers.
      The jrealft routine requires the input array to start at index=1
      so we must decrement the index of amp
  */
 
  jrealft(amp-1, nn, isign);
  
}

void  mtm::do_mtap_spec( double *data, 
			 int npoints, 
			 int kind,
			 int nwin, 
			 double npi, 
			 int inorm, 
			 double dt,
			 double *ospec, 
			 double *dof, 
			 double *Fvalues, 
			 int klen, 
			 bool display_tapers , 
			 std::vector<double> * write_tapers , 
			 std::vector<double> * write_tapsum , 
			 std::vector<double> * write_lambda , 
			 const std::vector<double> * read_tapers , 
			 const std::vector<double> * read_tapsum , 
			 const std::vector<double> * read_lambda )

{

  /*
    data = (type double) input time series
    npoints = number of points in data
    kind = flag for choosing hires or adaptive weighting coefficients

    nwin = number of taper windows to calculate

    npi = order of the slepian functions

    inorm = flag for choice of normalization

    dt = sampling interval (time)

    ospec = output spctrum

    dof = degrees of freedom at each frequency

    Fvalues = Ftest value at each frequency estimate

    klen = number of frequecies calculated (power of 2)
    
*/

  
  int             i, j, k;
  double         *lambda, *tapers;
  long            len, longlen;
  double          *xt;
  
  int             logg;
  int             nn;
  double          *b;
  int             iwin, kk;
  
  /*************/
  double          anrm, norm;
  double            *ReSpec, *ImSpec;
  double         *sqr_spec,  *amu;
  double          *amp, *fv;
  double          avamp, temp, sqramp;
  double          sum, *tapsum;
  /************/
  int num_freqs;
  int len_taps, num_freq_tap;
  
  double         *dcf, *degf, avar;
  int             n1, n2, kf;
  int             flag;
  int             one = 1;
  
  double tem1, tem2;
  
  /* lambda = vector of eigenvalues   
     tapsum = sum of each taper, saved for use in adaptive weighting  
     tapers =  matrix of slepian tapers, packed in a 1D double array    
  */
  
  lambda = dvector((long)0, (long)nwin);
  tapsum=dvector((long)0,(long)nwin);
  len_taps = npoints * nwin;
  tapers = dvector((long)0,(long) len_taps);
  
  num_freqs = 1+klen/2;
  num_freq_tap = num_freqs*nwin;
  
  //  std::cerr << "dets = " << klen << " " << num_freqs << "\n";

    
  //
  // calculate or read in Slepian tapers
  //
  
  if ( read_tapers && read_tapsum && read_lambda ) 
    {
      if ( read_tapers->size() != len_taps ) Helper::halt( "internal error, wrong saved taper length" );
      if ( read_tapsum->size() != nwin ) Helper::halt( "internal error, wrong saved taper length" );
      if ( read_lambda->size() != nwin ) Helper::halt( "internal error, wrong saved taper length" );
      for (int i=0;i<len_taps;i++) tapers[i] = (*read_tapers)[i];
      for (int i=0;i<nwin;i++) tapsum[i] = (*read_tapsum)[i];
      for (int i=0;i<nwin;i++) lambda[i] = (*read_lambda)[i];
    }
  else  // calculate tapers
    {
      k = multitap(npoints, nwin, lambda,  npi, tapers, tapsum);
    }
  

  
  //
  // write tapers?
  //

  if ( write_tapers ) 
    {
      write_tapers->resize( len_taps );
      for (int i=0;i<len_taps;i++) (*write_tapers)[i] = tapers[i];
    }

  if ( write_tapsum ) 
    {
      write_tapsum->resize( nwin );
      for (int i=0;i<nwin;i++) (*write_tapsum)[i] = tapsum[i];
    }

  if ( write_lambda ) 
    {
      write_lambda->resize( nwin );
      for (int i=0;i<nwin;i++) (*write_lambda)[i] = lambda[i];
    }

  // display tapers
  if ( display_tapers ) 
    {  
      for(i=0; i<npoints; i++)
	{
	  std::cout << "MTM" << "t"<<i;
	  for(j=0; j<nwin; j++) std::cout << "\t" << tapers[i+j*npoints];
	  std::cout << "\n";
	}
        
      for(j=0; j<nwin; j++) 
	{
	  std::cout << "LAMBDA " << j+1 << "\t" << lambda[j] << "\n";
	}
      
    }


  /* choose normalization based on inorm flag  */
  
  anrm = 1.;
  
  switch (inorm) {
  case 0:
    anrm = 1.;
    break;    
  case 1:
    anrm = npoints;
    break;
  case 2:
    anrm = 1 / dt;
    break;
  case 3:
    anrm = sqrt((double) npoints);
    break;
  case 4:
    //anrm = sqrt(npoints) * sqrt(1.0/dt);  // 1/(NFs)    
    anrm = sqrt(npoints/dt);
    break;
  default:
    anrm = 1.;
    break;
  }
  
  
  /* apply the taper in the loop.  do this nwin times  */
  
  b = vector((long)0, (long)npoints);
  amu = dvector((long)0,(long) num_freqs);
  sqr_spec = dvector((long)0,(long) num_freq_tap);
  ReSpec = dvector((long)0,(long) num_freq_tap);
  ImSpec = dvector((long)0,(long) num_freq_tap);
  
  
  for (iwin = 0; iwin < nwin; iwin++) {
    kk = iwin * npoints;
    kf = iwin * num_freqs;
    
    for (j = 0; j < npoints; j++)
      b[j] = data[j] * tapers[kk + j];   /*  application of  iwin-th taper   */
    
    amp = vector((long)0,(long) klen);
    
    mt_get_spec(b, npoints, klen, amp);  /* calculate the eigenspectrum */
    
    sum = 0.0;

    
    /* get spectrum from real fourier transform    */
    
    norm = 1.0/(anrm*anrm);
    
    for(i=1; i<num_freqs-1; i++){
      if(2*i+1 > klen) Helper::halt( "mtm_t error in index");
      if(i+kf > num_freq_tap ) Helper::halt( "mtm_t error in index");
      
      sqramp = SQR(amp[2*i+1])+SQR(amp[2*i]);
      
      ReSpec[i+kf] = amp[2*i];
      ImSpec[i+kf] = amp[2*i+1];
      
      sqr_spec[i+kf] =    norm*(sqramp);
      
      sum += sqr_spec[i+kf];
    }
    sqr_spec[0+kf] = norm*SQR(fabs(amp[0]));
    sqr_spec[num_freqs-1+kf] = norm*SQR(fabs(amp[1]));
    
    ReSpec[0+kf] = amp[0];
    ImSpec[0+kf] = 0.0;
    
    ReSpec[num_freqs-1+kf] = amp[1];
    ImSpec[num_freqs-1+kf] = 0.0;
    
    sum += sqr_spec[0+kf] + sqr_spec[num_freqs-1+kf];
    
    if(num_freqs-1+kf>num_freq_tap ) Helper::halt( "mtm_t error in index");
    
    temp = sum / (double) num_freqs;
    if (temp > 0.0)
      avamp = sqrt(temp) / anrm;
    else {
      avamp = 0.0;
      /* fprintf(stderr," avamp = 0.0! \n"); */ 
    }
    
    
    free_vector(amp,(long) 0,(long) klen);
    
  }
  
  free_vector(b, (long)0, (long)npoints);
  fv = vector((long)0,(long) num_freqs);


  
  //
  // Hi-res or adaptive weighting for spectra
  //
  
  
  switch (kind) 
    {
      //
      // hi-res
      //

    case 1:
      
      hires(sqr_spec,  lambda, nwin, num_freqs, amu);
      get_F_values(ReSpec, ImSpec, num_freqs, nwin, fv, tapsum);
      
      for (i = 0; i < num_freqs; i++) 
	{
	  ospec[i] = amu[i];
	  dof[i] = nwin-1;
	  Fvalues[i] = fv[i];
	}
      
      

      break;

      //
      // adaptive weighting for spectra
      //
      
    case 2:
      
      /* get avar = variance*/
      
      n1 = 0;
      n2 = npoints;
      
      avar = 0.0;
      
      for (i = n1; i < n2; i++)
	avar += (data[i]) * (data[i]);
      
      
      switch (inorm) {
      case 0:
	avar = avar / npoints;
	break;
	
      case 1:
	avar = avar / (npoints * npoints);
	break;
	
      case 2:
	avar = avar * dt * dt;
	break;
	
      case 3:	
	avar = avar / npoints;
	break;
	
      case 4:
	avar = avar / ( npoints / dt ) ;
      default:
	break;
      }
      
      
      dcf = dvector((long)0,(long) num_freq_tap);
      degf = dvector((long)0,(long) num_freqs);
      
      adwait(sqr_spec, dcf, lambda, nwin, num_freqs, amu, degf, avar);
      
      get_F_values(ReSpec, ImSpec, num_freqs, nwin, fv, tapsum);
    
      /* rap up   */
      
      for (i = 0; i < num_freqs; i++) {
	ospec[i] =amu[i];
	dof[i] = degf[i];
	Fvalues[i] = fv[i];
      }
    
      
      free_dvector(dcf,(long)0,(long) num_freq_tap);
      free_dvector(degf,(long)0,(long) num_freqs);
      free_vector(fv,(long)0,(long) num_freqs);
      
      
      break;
    }
  

  //
  // Free up memory, return...
  //
  
  free_dvector(amu,(long)0,(long) num_freqs);  
  
  free_dvector(sqr_spec, (long)0,(long) num_freq_tap);

  free_dvector(ReSpec, (long)0,(long) num_freq_tap);
  
  free_dvector(ImSpec, (long)0,(long) num_freq_tap);
  
  free_dvector(lambda,(long) 0,(long) nwin);
  
  free_dvector(tapers,(long) 0, (long)len_taps);

  free_dvector(tapsum,(long) 0, (long)nwin);
  
}
