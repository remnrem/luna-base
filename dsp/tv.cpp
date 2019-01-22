
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

#include "tv.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "eval.h"

#include "helper/helper.h"
#include "helper/logger.h"

#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>


void dsptools::tv( edf_t & edf , param_t & param )
{
  
  int lambda = param.requires_dbl( "lamdba" );
  if ( lambda < 0 ) Helper::halt( "lambda must be >= 0" ) ; 
  
  std::string signal_label = param.requires( "signal" );
  signal_list_t signals = edf.header.signal_list( signal_label );    
  const int ns = signals.size();
  
  // for each signal

  for (int s=0;s<ns;s++)
    {
      
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
      
      // pull entire trace (assumes continuous/contiguous data)
      
      interval_t interval = edf.timeline.wholetrace();
      
      slice_t slice1( edf , signals(s) , interval );    
      const std::vector<double> * d = slice1.pdata();
      
      // denoise
      std::vector<double> denoised = TV1D_denoise_copy( *d , lambda );
      
      // update 
      edf.update_signal( signals(s) , &denoised );

    }

  // all done
}


// Here we just provide a wrapper around the core function written by Laurent Condat, as described below:

/*
 #  File            : condat_fast_tv.c 
 #
 #  Version History : 1.0, Feb. 2012
 #
 #  Author          : Laurent Condat, PhD, CNRS research fellow in France.
 #
 #  Description     : This file contains an implementation in the C language
 #                    of algorithms described in the research paper:
 #	
 #                    L. Condat, "A Direct Algorithm for 1D Total Variation
 #                    Denoising", preprint hal-00675043, 2012.
 #
 #                    This implementation comes with no warranty: due to the
 #                    limited number of tests performed, there may remain
 #                    bugs. In case the functions would not do what they are
 #                    supposed to do, please email the author (contact info
 #                    to be found on the web).
 #
 #                    If you use this code or parts of it for any purpose,
 #                    the author asks you to cite the paper above or, in 
 #                    that event, its published version. Please email him if 
 #                    the proposed algorithms were useful for one of your 
 #                    projects, or for any comment or suggestion.
 #
 #  Usage rights    : Copyright Laurent Condat.
 #                    This file is distributed under the terms of the CeCILL
 #                    licence (compatible with the GNU GPL), which can be
 #                    found at the URL "http://www.cecill.info".
 #
 #  This software is governed by the CeCILL license under French law and
 #  abiding by the rules of distribution of free software. You can  use,
 #  modify and or redistribute the software under the terms of the CeCILL
 #  license as circulated by CEA, CNRS and INRIA at the following URL :
 #  "http://www.cecill.info".
 #
 #  As a counterpart to the access to the source code and rights to copy,
 #  modify and redistribute granted by the license, users are provided only
 #  with a limited warranty  and the software's author,  the holder of the
 #  economic rights,  and the successive licensors  have only  limited
 #  liability.
 #
 #  In this respect, the user's attention is drawn to the risks associated
 #  with loading,  using,  modifying and/or developing or reproducing the
 #  software by the user in light of its specific status of free software,
 #  that may mean  that it is complicated to manipulate,  and  that  also
 #  therefore means  that it is reserved for developers  and  experienced
 #  professionals having in-depth computer knowledge. Users are therefore
 #  encouraged to load and test the software's suitability as regards their
 #  requirements in conditions enabling the security of their systems and/or
 #  data to be ensured and,  more generally, to use and operate it in the
 #  same conditions as regards security.
 #
 #  The fact that you are presently reading this means that you have had
 #  knowledge of the CeCILL license and that you accept its terms.
*/




//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>


/* 
   This function implements the 1D total variation denoising 
   algorithm described in the paper referenced above. 
   If output=input, the process is performed in place. Else, 
   the values of input are left unchanged. 
   lambda must be nonnegative. lambda=0 is admissible and 
   yields output[k]=input[k] for all k. 
   If width<=0, nothing is done. 
*/


std::vector<double> dsptools::TV1D_denoise_copy(const std::vector<double> & input, const double lambda) 
{
  std::vector<double> copy = input;
  TV1D_denoise( copy , lambda );
  return copy;
}


void dsptools::TV1D_denoise(std::vector<double> & input, const double lambda) 
{
  
  // do this so that operation is performed in place
  std::vector<double> & output = input;

  const int width = input.size();
  
  if (width>0) 
    {				/*to avoid invalid memory access to input[0]*/
      int k=0, k0=0;			/*k: current sample location, k0: beginning of current segment*/
      double umin=lambda, umax=-lambda;	/*u is the dual variable*/
      double vmin=input[0]-lambda, vmax=input[0]+lambda;	/*bounds for the segment's value*/
      int kplus=0, kminus=0; 	/*last positions where umax=-lambda, umin=lambda, respectively*/
      const double twolambda=2.0*lambda;	/*auxiliary variable*/
      const double minlambda=-lambda;		/*auxiliary variable*/
      for (;;) {				/*simple loop, the exit test is inside*/
	while (k==width-1) {	/*we use the right boundary condition*/
	  if (umin<0.0) {			/*vmin is too high -> negative jump necessary*/
	    do output[k0++]=vmin; while (k0<=kminus);
	    umax=(vmin=input[kminus=k=k0])+(umin=lambda)-vmax;
	  } else if (umax>0.0) {	/*vmax is too low -> positive jump necessary*/
	    do output[k0++]=vmax; while (k0<=kplus);
	    umin=(vmax=input[kplus=k=k0])+(umax=minlambda)-vmin;
	  } else {
	    vmin+=umin/(k-k0+1); 
	    do output[k0++]=vmin; while(k0<=k); 
	    return;
	  }
	}
	if ((umin+=input[k+1]-vmin)<minlambda) {		/*negative jump necessary*/
	  do output[k0++]=vmin; while (k0<=kminus);
	  vmax=(vmin=input[kplus=kminus=k=k0])+twolambda;
	  umin=lambda; umax=minlambda;
	} else if ((umax+=input[k+1]-vmax)>lambda) {	/*positive jump necessary*/
	  do output[k0++]=vmax; while (k0<=kplus);
	  vmin=(vmax=input[kplus=kminus=k=k0])-twolambda;
	  umin=lambda; umax=minlambda;
	} else { 	/*no jump necessary, we continue*/
	  k++;
	  if (umin>=lambda) {		/*update of vmin*/
	    vmin+=(umin-lambda)/((kminus=k)-k0+1);
	    umin=lambda;
	  } 
	  if (umax<=minlambda) {	/*update of vmax*/
	    vmax+=(umax+lambda)/((kplus=k)-k0+1);
	    umax=minlambda;
	  } 	
	}
      }
    }
}


/* 
   This function implements the fused lasso signal approximator,
   described in the paper referenced above. 
   If output=input, the process is performed in place. Else, the 
   values of input are left unchanged. 
   lambda and mu must be nonnegative. If mu=0, the function does 
   the same as TV1D_denoise.
   If width<=0, nothing is done. 
*/

void dsptools::fused_lasso(double* input, double* output, const int width, const double lambda, const double mu) {
  if (width>0) {				/*to avoid invalid memory access to input[0]*/
    int k=0, k0=0;			/*k: current sample location, k0: beginning of current segment*/
    double umin=lambda, umax=-lambda;	/*u is the dual variable*/
    double vmin=input[0]-lambda, vmax=input[0]+lambda;	/*bounds for the segment's value*/
    int kplus=0, kminus=0; 	/*last positions where umax=-lambda, umin=lambda, respectively*/
    const double twolambda=2.0*lambda;	/*auxiliary variable*/
    const double minlambda=-lambda;		/*auxiliary variable*/
    for (;;) {				/*simple loop, the exit test is inside*/
      while (k==width-1) {	/*we use the right boundary condition*/
	if (umin<0.0) {			/*vmin is too high -> negative jump necessary*/
	  vmin = vmin>mu ? vmin-mu : vmin<-mu ? vmin+mu : 0.0;  
	  do output[k0++]=vmin; while (k0<=kminus);
	  umax=(vmin=input[kminus=k=k0])+(umin=lambda)-vmax;
	} else if (umax>0.0) {	/*vmax is too low -> positive jump necessary*/
	  vmax = vmax>mu ? vmax-mu : vmax<-mu ? vmax+mu : 0.0;
	  do output[k0++]=vmax; while (k0<=kplus);
	  umin=(vmax=input[kplus=k=k0])+(umax=minlambda)-vmin;
	} else {
	  vmin+=umin/(k-k0+1);
	  vmin = vmin>mu ? vmin-mu : vmin<-mu ? vmin+mu : 0.0;  
	  do output[k0++]=vmin; while(k0<=k); 
	  return;
	}
      }
      if ((umin+=input[k+1]-vmin)<minlambda) {		/*negative jump necessary*/
	vmin = vmin>mu ? vmin-mu : vmin<-mu ? vmin+mu : 0.0; 
	do output[k0++]=vmin; while (k0<=kminus);
	vmax=(vmin=input[kplus=kminus=k=k0])+twolambda;
	umin=lambda; umax=minlambda;
      } else if ((umax+=input[k+1]-vmax)>lambda) {	/*positive jump necessary*/
	vmax = vmax>mu ? vmax-mu : vmax<-mu ? vmax+mu : 0.0;
	do output[k0++]=vmax; while (k0<=kplus);
	vmin=(vmax=input[kplus=kminus=k=k0])-twolambda;
	umin=lambda; umax=minlambda;
      } else { 	/*no jump necessary, we continue*/
	k++;
	if (umin>=lambda) {		/*update of vmin*/
	  vmin+=(umin-lambda)/((kminus=k)-k0+1);
	  umin=lambda;
	} 
	if (umax<=minlambda) {	/*update of vmax*/
	  vmax+=(umax+lambda)/((kplus=k)-k0+1);
	  umax=minlambda;
	} 	
      }
    }
  }
}




