
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


// MTM functions are adapted code from: Lees, J. M. and J. Park
// (1995): Multiple-taper spectral analysis: A stand-alone
// C-subroutine: Computers & Geology: 21, 199-236.


#include "mtm.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "eval.h"
#include "fftw/fftwrap.h"

#include "db/db.h"
#include "helper/helper.h"
#include "helper/logger.h"

extern writer_t writer;
extern logger_t logger; 


mtm_t::mtm_t( const double npi , const int nwin ) : npi(npi) , nwin(nwin) 
{

  // by default, set to use 'adaptive weights' (2)
  kind = 2 ; 
  
  // set to use 1/(N.Fs) weights (4)
  inorm = 4 ; 
  
}


void mtm_t::store_tapers( const int seg_size ) 
{
  
  // Vector of eigenvalues
  lam    = Eigen::VectorXd::Zero( nwin );
  
  // Sum of each taper (used in adaptive weighting)
  tapsum = Eigen::VectorXd::Zero( nwin );

  // samples x tapers
  tapers = Eigen::MatrixXd::Zero( seg_size , nwin );
  
  // calculate Slepian tapers                                                                                                                                                 
  generate_tapers( seg_size, nwin, npi );
  
}


void mtm_t::apply( const std::vector<double> * d , const int fs , 
		   const int seg_size , const int seg_step , bool verbose , mtm_t * precomputed )
{

  const bool allsegs = restrict.size() == 0 ;
  
  const double dt = 1.0/(double)fs;

  const int total_npoints = d->size();
  
  // spectral window

  const int npoints = seg_size;
  
  const double fWidth =  npi/((double)npoints*dt);
  
  const int K = (int) 2*npoints*fWidth*dt;
  
  const double nyquist = 0.5/dt;
  
  const int klen = mtm_t::get_pow_2( npoints );
  
  const double df = 2*nyquist/klen;  
  
  const int nfreqs = 1+klen/2;
  
  int k = 1;

  int n_segs = 0;
  for (int p=0; p<total_npoints; p += seg_step )
    {
      if ( p + seg_size > total_npoints ) break;      
      ++n_segs;
    }

  if ( !allsegs )
    if ( restrict.size() != n_segs )
      Helper::halt("internal error in mtm w/ comp-segs implied" );

  int n_segs_actual = 0;
  if ( !allsegs )
    {
      for (int p=0; p<n_segs; p++)
	if ( ! restrict[p] )
	  ++n_segs_actual;
    }
  else
    n_segs_actual = n_segs;
  
				   
  const double spectral_resolution = ( 2 * npi ) / ( seg_size / (double)fs );

  if ( verbose ) 
    logger << "  assuming all channels have the same sample rate of " << fs << "Hz:\n"
	   << "    time half-bandwidth (nw) = " << npi << "\n"
	   << "    number of tapers         = " << nwin << "\n"
	   << "    spectral resolution      = " << spectral_resolution << "Hz\n"      
	   << "    segment duration         = " << seg_size / (double)fs << "s\n"
	   << "    segment step             = " << seg_step / (double)fs << "s\n"
	   << "    FFT size                 = " << klen << "\n"
	   << "    number of segments       = " << n_segs << "\n";

  if ( ! allsegs )
    logger << "    computed segments        = " << n_segs_actual << "\n"; 
  
  logger << "    adjustment               = "
	 << ( opt_remove_trend ? "detrend" : ( opt_remove_mean ? "constant" : "none" ) ) << "\n";
  
  //
  // Generate and store tapers
  //
  
  // Vector of eigenvalues
  lam    = Eigen::VectorXd::Zero( nwin );
  
  // Sum of each taper (used in adaptive weighting)
  tapsum = Eigen::VectorXd::Zero( nwin );

  // samples x tapers
  tapers = Eigen::MatrixXd::Zero( npoints , nwin );


  //
  // calculate (or attached pre-computed) Slepian tapers
  //

  if ( precomputed == NULL ) 
    generate_tapers( npoints, nwin, npi );
  else
    {
      lam = precomputed->lam;
      tapsum = precomputed->tapsum;
      tapers = precomputed->tapers;
    }
  
  //
  // Spectrogram output (for all n_segs whether all computed or not)
  //
  
  espec.resize( n_segs );
  raw_espec.resize( n_segs );

  //
  // Initiate FFT (no window, as tapers applied to data beforehand)
  //

  // use next pow2 FFT by default:
  real_FFT fftseg( seg_size , klen , fs , WINDOW_NONE );
  //  real_FFT fftseg( seg_size , seg_size , fs , WINDOW_NONE );
  
  //
  // Iterate over segments
  //

  int sn = 0; // count of segment  (whether processed or no)
  bool done_first = false;
  
  for ( int p = 0; p < total_npoints ; p += seg_step )
    {

      // all done?
      if ( p + seg_size > total_npoints ) break;

      // skip this segment?
      if ( ( ! allsegs ) && restrict[sn] ) { ++sn; continue; }
      
      // need to copy segment (i.e. if detrending)
      std::vector<double> segment( npoints );  // == seg_size
      
      int p2 = p;
      for (int j=0; j<seg_size; j++)
	segment[j] = (*d)[p2++];
      
      double * psegment = &(segment)[0];

      //
      // remove mean or detrend?
      //
      
      if ( opt_remove_mean ) 
	{
	  double m = mtm_t::remove_mean( psegment, npoints );
	}
      else if ( opt_remove_trend )
	{
	  rm_lin_sig_trend( psegment , npoints , dt );
	}

      
      // allocate storage for spectra
      espec[sn].resize( klen ,  0 );  
      
      // track raw spectra too 
      raw_espec[sn].resize( klen ,  0 );  

      // do actual MTM analysis
      do_mtap_spec( &fftseg,
		    &(segment)[0],
		    npoints ,
		    kind, nwin, npi, inorm, dt,
		    &(raw_espec)[sn][0],  klen );      
      
      //
      // shrink to positive spectrum (already scled x2)?
      // 
      
      raw_espec[sn].resize( nfreqs );
      espec[sn].resize( nfreqs );
      
      if ( ! done_first )
	f.resize( nfreqs , 0 );
      
      for (int i = 0; i < nfreqs; i++)
	{

	  if ( ! done_first )
	    f[i] = df*i;
	  
	  // report dB?
	  if ( dB )
	    espec[sn][i] = 10 * log10( raw_espec[sn][i] );
	  else
	    espec[sn][i] = raw_espec[sn][i] ;	  
	}  
      
      
      //
      // Next segment
      //

      done_first = true;
      
      ++sn;

    }

  //
  // Compute average spectrum (using only actually computed segments)
  //

  spec.resize( nfreqs , 0 );
  raw_spec.resize( nfreqs , 0 );
  
  for (int f=0; f<nfreqs; f++)
    {
      for (int i=0; i<n_segs; i++)
	{
	  if ( allsegs || ! restrict[i] )
	    {
	      spec[f] += espec[i][f];
	      raw_spec[f] += raw_espec[i][f];
	    }
	}
      spec[f] /= (double)n_segs_actual;      
      raw_spec[f] /= (double)n_segs_actual;
    }
 
}




// -------------------------------------------------------------------------------
//
// Generate tapers
//
// -------------------------------------------------------------------------------


#define PI 3.14159265358979
#define ABS(a) ((a) < (0) ? -(a) : (a))
#define DIAG1 0
#define MAX(a,b) ((a) >= (b) ? (a) : (b))


void  mtm_t::generate_tapers( int num_points, int nwin, double npi )
{

  //
  // Get the multitaper slepian functions: 

  // num_points = number of points in data stream
  // nwin = number of windows 
  // npi = order of slepian functions 
  
  // modifies members: (already sized before this is called)
  //  lam    = vector of eigenvalues 
  //  tapsum = sum of each taper, saved for use in adaptive weighting  
  //  tapers =  matrix of slepian tapers, packed in a 1D double array

  /* need to initialize iwflag = 0 */

  int             key, nbin, npad;
  long            len;
  int             ierr;

  const double DPI = (double) PI;
  const double TWOPI = (double) 2 * DPI;
  
  const double anpi = npi;
  const double an = (double) (num_points);
  const double ww = (double) (anpi) / an;	/* this corresponds to P&W's W value  */
  const double cs = cos( TWOPI * ww );
  
  Eigen::VectorXd ell = Eigen::VectorXd::Zero( nwin );
  Eigen::VectorXd diag = Eigen::VectorXd::Zero( num_points );
  Eigen::VectorXd offdiag = Eigen::VectorXd::Zero( num_points );
  Eigen::VectorXd offsq = Eigen::VectorXd::Zero( num_points );

  Eigen::VectorXd scratch1 = Eigen::VectorXd::Zero( num_points );
  Eigen::VectorXd scratch2 = Eigen::VectorXd::Zero( num_points );
  Eigen::VectorXd scratch3 = Eigen::VectorXd::Zero( num_points );
  Eigen::VectorXd scratch4 = Eigen::VectorXd::Zero( num_points );
  Eigen::VectorXd scratch6 = Eigen::VectorXd::Zero( num_points );
  
  //
  //  make the diagonal elements of the tridiag matrix 
  //

  for (int i=0; i < num_points; i++)
    {
      double ai = (double) (i);
      diag[i] = -cs * (((an - 1.) / 2. - ai)) * (((an - 1.) / 2. - ai));
      offdiag[i] = -ai * (an - ai) / 2.;
      offsq[i] = offdiag[i] * offdiag[i];
    }
  
  double eps = 1.0e-13;

  int m11 = 1;
  
  std::vector<int> ip( nwin , 0 );

  // call the eispac routines to invert the tridiagonal system 

  double rlb , rlu;
  
  jtridib_( &num_points,
	    &eps,
	    diag.data(),
	    offdiag.data(),
	    offsq.data(),
	    &rlb,
	    &rlu,
	    &m11,
	    &nwin,
	    lam.data(),
	    &(ip[0]),
	    &ierr,
	    scratch1.data(),
	    scratch2.data());
  
  len = num_points * nwin;

  Eigen::VectorXd evecs = Eigen::VectorXd::Zero( len );
  
  jtinvit_( &num_points,
	    &num_points,
	    diag.data(),
	    offdiag.data(),
	    offsq.data(),
	    &nwin,
	    lam.data(),
	    &(ip[0]),
	    evecs.data(),
	    &ierr,
	    scratch1.data(), scratch2.data(), scratch3.data(), scratch4.data(), scratch6.data() );
  

  // we calculate the eigenvalues of the dirichlet-kernel problem i.e.
  // the bandwidth retention factors from slepian 1978 asymptotic
  // formula, gotten from thomson 1982 eq 2.5 supplemented by the
  // asymptotic formula for k near 2n from slepian 1978 eq 61 more
  // precise values of these parameters, perhaps useful in adaptive
  // spectral estimation, can be calculated explicitly using the
  // rayleigh-quotient formulas in thomson (1982) and park et al (1987)
  
  
  double dfac = (double) an * DPI * ww;
  double drat = (double) 8. * dfac;
    
  dfac = (double) 4. * sqrt(DPI * dfac) * exp( (double) (-2.0) * dfac);
  
  
  for (int k = 0; k < nwin; k++)
    {
      lam[k] = (double) 1.0 - (double) dfac;
      dfac = dfac * drat / (double) (k + 1);
      /* fails as k -> 2n */
    }
  
  
  const double gamma = log( (double) 8. * an * sin((double) 2. * DPI * ww)) + (double) 0.5772156649;
  
  for (int k = 0; k < nwin; k++)
    {
      const double bh = -2. * DPI * (an * ww - (double) (k) /
				     (double) 2. - (double) .25) / gamma;
      ell[k] = (double) 1. / ((double) 1. + exp(DPI * (double) bh));      
    }
  
  for (int i = 0; i < nwin; i++)
    lam[i] = MAX( ell[i], lam[i] );
  
  
  // Normalize the eigentapers to preserve power for a white process
  //   i.e. they have rms value unity tapsum is the average of the
  //   eigentaper, should be near zero for antisymmetric tapers

  for (int k = 0; k < nwin; k++) 
    {
      const int kk = (k) * num_points;
      tapsum[k] = 0.;
      double tapsq = 0.;
      for (int i = 0; i < num_points; i++) 
	{
	  double aa = evecs[i + kk];
	  tapers(i , k ) = aa;
	  tapsum[k] = tapsum[k] + aa;
	  tapsq = tapsq + aa * aa;
	}
      
      double aa = sqrt( tapsq / (double) num_points);
      
      tapsum[k] = tapsum[k] / aa;
      
      for (int i = 0; i < num_points; i++) 
	tapers( i ,k ) = tapers( i , k ) / aa;
      
    } // next taper
    
}


// ----------------------------------------------------------------------------
//
// Core MTM function
//
// ----------------------------------------------------------------------------


void  mtm_t::do_mtap_spec( real_FFT * fftseg ,
			   double *data, 
			   int npoints, 
			   int kind,
			   int nwin, 
			   double npi, 
			   int inorm, 
			   double dt,
			   double *ospec, 
			   int klen ,
			   double *dof, 
			   double *Fvalues )

{
  
  // data = (type double) input time series
  // npoints = number of points in data
  // kind = flag for choosing hires or adaptive weighting coefficients
  // klen = number of frequecies calculated (power of 2)
  // nwin = number of taper windows to calculate
  // npi = order of the slepian functions  
  // inorm = flag for choice of normalization
  // dt = sampling interval (time)
  // ospec = output spectrum
  
  // len_taps = npoints * nwin;
  
  // number of frequencies
  int num_freqs = 1+klen/2;

  int num_freq_tap = num_freqs * nwin;
  

  //
  // Set up normalization 
  //
  
  double anrm = 1.;
  
  switch (inorm)
    {
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
      anrm = sqrt(npoints/dt);
      break;
    default:
      anrm = 1.;
      break;
    }
  
  
  //
  // Apply the taper in the loop, nwin times
  //

  std::vector<double> b( npoints , 0 );
  std::vector<double> amu( num_freqs , 0 );
  std::vector<double> sqr_spec( num_freq_tap , 0 );
  // std::vector<double> ReSpec( num_freq_tap , 0 );
  // std::vector<double> ImSpec( num_freq_tap , 0 );
  
  for (int iwin = 0; iwin < nwin; iwin++)
    {
      const int kk = iwin * npoints;
      const int kf = iwin * num_freqs;
  
      for (int j = 0; j < npoints; j++)
	b[j] = data[j] * tapers( j , iwin );  /*  application of  iwin-th taper   */
  
      std::vector<double> amp( klen , 0 );
  
      //
      // do FFT
      //

      fftseg->apply( &(b)[0] , npoints );
      
      //
      // get spectrum from real fourier transform: populate sqr_spec
      //
      
      for (int f=0; f<fftseg->cutoff; f++)	
	sqr_spec[ f + kf ] = fftseg->X[f];
	  
      if( num_freqs-1+kf > num_freq_tap )
	Helper::halt( "mtm_t error in index");
      
      // next taper
    }

  
  std::vector<double> fv( num_freqs , 0 );
  
  //
  // Weighting of spectra ('hi-res' or 'adaptive')
  //
  
  switch ( kind ) 
    {
      
    case 1: // 'hi-res'
      
      hires( &(sqr_spec)[0],  lam.data() , nwin, num_freqs, &(amu)[0] );

      // nb. disabled this part: requires extracting real and imaginary components above
      // if ( Fvalues != NULL )
      // 	get_F_values(ReSpec, ImSpec, num_freqs, nwin, fv, tapsum);
      
      for (int i = 0; i < num_freqs; i++) 
	{
	  ospec[i] = amu[i];
	  // if ( Fvalues != NULL )
	  //   {
	  //     dof[i] = nwin-1;
	  //     Fvalues[i] = fv[i];
	  //   }
	}
      
      break;
      
    case 2: // 'adaptive'
      
      // get avar = variance
      
      int n1 = 0;
      int n2 = npoints;    
      double avar = 0.0;
      
       for (int i = n1; i < n2; i++)
	 avar += (data[i]) * (data[i]);
       
       switch ( inorm )
	 {
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

       std::vector<double> dcf( num_freq_tap , 0 );
       std::vector<double> degf( num_freqs , 0 );
       
       adwait( &(sqr_spec)[0], &(dcf)[0], lam.data(), nwin, num_freqs, &(amu)[0], &(degf)[0], avar );
       
        // if ( Fvalues != NULL ) 
  	//  get_F_values(ReSpec, ImSpec, num_freqs, nwin, fv, tapsum);
       
       for (int i = 0; i < num_freqs; i++)
	 {
	   ospec[i] = amu[i];
	   // if ( Fvalues != NULL )
	   //   {
	   //     dof[i] = degf[i];
	   //     Fvalues[i] = fv[i];
	   //   }
	 }
    
       break;
     }
  
}


