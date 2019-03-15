
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

#include "coherence.h"

#include "resample.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "eval.h"
#include "db/db.h"
#include "fftw/fftwrap.h"
#include "helper/logger.h"
#include "helper/helper.h"

// Coherence code from : https://www.physionet.org/physiotools/wfdb/psd/coherence.c

extern writer_t writer;

extern logger_t logger;

#include <cmath>


void dsptools::coherence( edf_t & edf , param_t & param , bool legacy )
{
  
  std::string signal_label = param.requires( "sig" );
  
  signal_list_t signals = edf.header.signal_list( signal_label );  
  
  const int ns = signals.size();
  
  int sr = param.has( "sr" ) ? param.requires_int( "sr" ) : 0 ;
  
  writer.var( "COH" , "Spectral coherence (0..1)" );
  
  bool show_spectrum = param.has( "spectrum" );
  
  double upper_freq = param.has("max") ? param.requires_dbl( "max" ) : 20 ; 

  // adjust all SRs now
  if ( sr )
    {
      for (int s=0;s<ns;s++)
	{
	  
	  if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
	  
	  if ( edf.header.sampling_freq( signals(s) ) != sr ) 
	    {
	      logger << "resampling channel " << signals.label(s) 
		     << " from " << edf.header.sampling_freq( signals(s) )
		     << " to " << sr << "\n";
	      resample_channel( edf, signals(s) , sr );
	    }
	}
    }
  
  
  //
  // Epochs or whole signal?
  //

  bool epoched = edf.timeline.epoched() && param.has("epoch") ;
  

  //
  // Ensure similar sampling rates
  //
  
  for (int i=0;i<ns-1;i++)
    {
      
      if ( edf.header.is_annotation_channel( signals(i) ) ) continue;
      
      for (int j=i+1;j<ns;j++)
	{
	  
	  if ( edf.header.is_annotation_channel( signals(j) ) ) continue;
	  
	  if ( epoched ) 
	    {
	      
	      int epoch = edf.timeline.first_epoch();      

	      const int sr1 = edf.header.sampling_freq( signals(i) );
	      const int sr2 = edf.header.sampling_freq( signals(j) );	      
	      if ( sr1 != sr2 ) Helper::halt( "'COH epoch' requires similiar sampling rates (or specify, e.g., sr=200)" );

	      
	      // stratify output by SIGNALS
	      writer.level( signals.label(i) + "_x_" + signals.label(j) , "CHS" );

	      while ( 1 ) 
		{

		  int epoch = edf.timeline.next_epoch();      
		  if ( epoch == -1 ) break;
		  interval_t interval = edf.timeline.epoch( epoch );

// 		  slice_t slice1( edf , signals(i) , interval );
// 		  slice_t slice2( edf , signals(j) , interval );
		  
// 		  const std::vector<double> * d1 = slice1.pdata();
// 		  const std::vector<double> * d2 = slice2.pdata();

		  coh_t coh = dsptools::coherence( edf , signals(i) , signals(j) , interval , legacy );
		  
		  //
		  // Summarize into bands
		  //

		  double coh_slow = 0 , 
		    coh_delta = 0 , 
		    coh_theta = 0 , 
		    coh_alpha = 0 , 
		    coh_sigma = 0 , 
		    coh_beta = 0; 
		  
		  double n_slow = 0 , 
		    n_delta = 0 , 
		    n_theta = 0 , 
		    n_alpha = 0 , 
		    n_sigma = 0 , 
		    n_beta = 0; 
		  
		  int sz = coh.frq.size();

		  for (int k=0;k<sz;k++)
		    {
		      
		      if ( coh.frq[k] >= globals::freq_band[ SLOW ].first 
			   && coh.frq[k] < globals::freq_band[ SLOW ].second ) 
			{			  
			  coh_slow += coh.coh[k];
			  n_slow++;
			}

		      if ( coh.frq[k] >= globals::freq_band[ DELTA ].first 
			   && coh.frq[k] < globals::freq_band[ DELTA ].second ) 
			{
			  coh_delta += coh.coh[k];
			  n_delta++;
			}
		      
		      if ( coh.frq[k] >= globals::freq_band[ THETA ].first 
			   && coh.frq[k] < globals::freq_band[ THETA ].second ) 
			{
			  coh_theta += coh.coh[k];
			  n_theta++;
			}

		      if ( coh.frq[k] >= globals::freq_band[ ALPHA ].first 
			   && coh.frq[k] < globals::freq_band[ ALPHA ].second ) 
			{
			  coh_alpha += coh.coh[k];
			  n_alpha++;
			}

		      if ( coh.frq[k] >= globals::freq_band[ SIGMA ].first 
			   && coh.frq[k] < globals::freq_band[ SIGMA ].second ) 
			{
			  coh_sigma += coh.coh[k];
			  n_sigma++;
			}

		      if ( coh.frq[k] >= globals::freq_band[ BETA ].first 
			   && coh.frq[k] < globals::freq_band[ BETA ].second ) 
			{
			  coh_beta += coh.coh[k];
			  n_beta++;
			}
		      
		    }
		  
		  if ( n_slow  ) coh_slow /= (double)n_slow;
		  if ( n_delta ) coh_delta /= (double)n_delta; 
		  if ( n_theta ) coh_theta /= (double)n_theta;
		  if ( n_alpha ) coh_alpha /= (double)n_alpha; 
		  if ( n_sigma ) coh_sigma /= (double)n_sigma;
		  if ( n_beta ) coh_beta /= (double)n_beta;

		  writer.epoch( edf.timeline.display_epoch( epoch ) );		  				

		  writer.level( globals::band( SLOW ) , globals::band_strat );
		  writer.value( "COH" , coh_slow );

		  writer.level( globals::band( DELTA ) , globals::band_strat );
		  writer.value( "COH" , coh_delta );

		  writer.level( globals::band( THETA ) , globals::band_strat );
		  writer.value( "COH" , coh_theta );

		  writer.level( globals::band( ALPHA ) , globals::band_strat );
		  writer.value( "COH" , coh_alpha );

		  writer.level( globals::band( SIGMA ) , globals::band_strat );
		  writer.value( "COH" , coh_sigma );

		  writer.level( globals::band( BETA ) , globals::band_strat );
		  writer.value( "COH" , coh_beta );

		  writer.unlevel( globals::band_strat );
		  

		  if ( show_spectrum )
		    {
		      
		      for (int k=0;k<sz;k++)
			{
			  if ( upper_freq < 0 || coh.frq[k] <= upper_freq )
			    {
			      writer.level( coh.frq[k] , globals::freq_strat );			  
			      writer.value( "COH" , coh.coh[k] );
			      writer.value( "COH" , coh.coh[k] );
			      writer.value( "CSPEC" , coh.cross_spectrum[k] );
			      writer.value( "ASPEC1" , coh.auto_spectrum1[k] );
			      writer.value( "ASPEC2" , coh.auto_spectrum2[k] );
			      
			      if ( ! legacy )
				{
				  writer.value( "CSPEC.N1" , coh.cross_norm1[k] );
				  writer.value( "CSPEC.N2" , coh.cross_norm2[k] );
				}
			    }
			}
		      
		      writer.unlevel( globals::freq_strat );
		      
		    } // end of epoch-spectrum level output
		  
		} // next epoch
	      
	      writer.unlevel( "CHS" );
	      writer.unepoch();
	      
	    }
	  else
	    {
	      
	      //
	      // Coherence for entire signal
	      //

	      coh_t coh = coherence( edf , signals(i) , signals(j) , edf.timeline.wholetrace() , legacy  );
	      
	      int sz = coh.frq.size();
	      
	      writer.level( signals.label(i) + "_x_" + signals.label(j) , "CHS" );
	      
	      //
	      // Band-level summaries
	      //

	      double coh_slow = 0 , 
		coh_delta = 0 , 
		coh_theta = 0 , 
		coh_alpha = 0 , 
		coh_sigma = 0 , 
		coh_beta = 0; 
		  
	      double n_slow = 0 , 
		n_delta = 0 , 
		n_theta = 0 , 
		n_alpha = 0 , 
		n_sigma = 0 , 
		n_beta = 0; 
		  
	      for (int k=0;k<sz;k++)
		{
		  
		  if ( coh.frq[k] >= globals::freq_band[ SLOW ].first 
		       && coh.frq[k] < globals::freq_band[ SLOW ].second ) 
		    {			  
		      coh_slow += coh.coh[k];
		      n_slow++;
		    }
		  
		  if ( coh.frq[k] >= globals::freq_band[ DELTA ].first 
		       && coh.frq[k] < globals::freq_band[ DELTA ].second ) 
		    {
		      coh_delta += coh.coh[k];
		      n_delta++;
		    }
		  
		  if ( coh.frq[k] >= globals::freq_band[ THETA ].first 
		       && coh.frq[k] < globals::freq_band[ THETA ].second ) 
		    {
		      coh_theta += coh.coh[k];
		      n_theta++;
		    }
		  
		  if ( coh.frq[k] >= globals::freq_band[ ALPHA ].first 
		       && coh.frq[k] < globals::freq_band[ ALPHA ].second ) 
		    {
		      coh_alpha += coh.coh[k];
		      n_alpha++;
		    }
		  
		  if ( coh.frq[k] >= globals::freq_band[ SIGMA ].first 
		       && coh.frq[k] < globals::freq_band[ SIGMA ].second ) 
		    {
		      coh_sigma += coh.coh[k];
		      n_sigma++;
		    }
		  
		  if ( coh.frq[k] >= globals::freq_band[ BETA ].first 
		       && coh.frq[k] < globals::freq_band[ BETA ].second ) 
		    {
		      coh_beta += coh.coh[k];
		      n_beta++;
		    }
		  
		}
	      
	      if ( n_slow  ) coh_slow /= (double)n_slow;
	      if ( n_delta ) coh_delta /= (double)n_delta; 
	      if ( n_theta ) coh_theta /= (double)n_theta;
	      if ( n_alpha ) coh_alpha /= (double)n_alpha; 
	      if ( n_sigma ) coh_sigma /= (double)n_sigma;
	      if ( n_beta ) coh_beta /= (double)n_beta;

	      writer.level( globals::band( SLOW ) , globals::band_strat );
	      writer.value( "COH" , coh_slow );

	      writer.level( globals::band( DELTA ) , globals::band_strat );
	      writer.value( "COH" , coh_delta );
	      
	      writer.level( globals::band( THETA ) , globals::band_strat );
	      writer.value( "COH" , coh_theta );

	      writer.level( globals::band( ALPHA ) , globals::band_strat );
	      writer.value( "COH" , coh_alpha );

	      writer.level( globals::band( SIGMA ) , globals::band_strat );
	      writer.value( "COH" , coh_sigma );

	      writer.level( globals::band( BETA ) , globals::band_strat );
	      writer.value( "COH" , coh_beta );

	      writer.unlevel( globals::band_strat );
		  

	      //
	      // Frequency-bin level output
	      //

	      for (int k=0;k<sz;k++)
		{
		  if ( upper_freq < 0 || coh.frq[k] <= upper_freq )
		    {
		      writer.level( coh.frq[k] , globals::freq_strat );
		      writer.value( "COH" , coh.coh[k] );
		      writer.value( "CSPEC" , coh.cross_spectrum[k] );
		      writer.value( "ASPEC1" , coh.auto_spectrum1[k] );
		      writer.value( "ASPEC2" , coh.auto_spectrum2[k] );

		      if ( ! legacy )
			{
			  writer.value( "CSPEC.N1" , coh.cross_norm1[k] );
			  writer.value( "CSPEC.N2" , coh.cross_norm2[k] );
			}

		    }		  
		}

	      writer.unlevel( globals::freq_strat );
	      writer.unlevel( "CHS" );

	    }
	}
    }
  
}


coh_t dsptools::legacy_coherence( const std::vector<double> * s1 , 
				  const std::vector<double> * s2 , 
				  double sampfreq )
{

  int pps = 1024;  // segment size
  // int pps = 512 ;  // segment size
  // int pps = 256 ;  // segment size
  
  // scaling factors (default = 1.0)
  double sfx = 1.0;
  double sfy = 1.0;

  if ( pps < 2 || pps > 32768 ) Helper::halt( "coherence::pps must be 2 < x < 32768" );
  
  /* points per Fourier transform segment (a power of 2) */
  int npfft;      
  
  for (npfft = 2; npfft < pps; npfft <<= 1)
    ;
  
  // Buffers
  
  // Inputs ( size 'npfft')
  // Outputs ( size 'npfft/2 + 1')

  const int szin = npfft;

  // npfft will always be even;
  // so include DC and Nyquist frequency
  const int szout = npfft/2 + 1;
  
  //std::cout << "npfft , sz in, out = " << npfft << " " << szin << " " << szout << " " << pps << "\n";

  // outputs
  //  double *gxx, *gyy, *gxyre, *gxyim, *phi, *weight;

  std::vector<double> gxx( szout );
  std::vector<double> gyy( szout );

  std::vector<double> gxyre( szout );
  std::vector<double> gxyim( szout );

  std::vector<double> phi( szout );


  
  //
  // Set ip
  //
  
  const int ntot = s1->size();

  int pos = 0;

  int nnn = pps; // number of points per segment

  //  Number of FFT outputs (half of the number of inputs).
  //int nd2 = npfft/2;  // note,this should have had a +1, given even nfft
  int nd2 = szout;
  

  // Compute Hann window.
  std::vector<double> weight = MiscMath::hann_window( szin );

  int nffts = 0;
  
  while ( 1 ) 
    {
      
      int nloaded = 0;

      // clear xx and yy
      std::vector<double> xx( szin , 0  );  // npfft (power of 2)
      std::vector<double> yy( szin , 0 ); // npfft (power of 2)
      

      for (int p = pos ; p < ntot ; p++ )
	{
	  if ( nloaded == nnn ) break;
	  xx[nloaded] = (*s1)[p];
	  yy[nloaded] = (*s2)[p];
	  ++nloaded;
	}
      
      // advance (a half-window)
      pos += nnn/2;
      
      if ( pos >= ntot || nloaded == 0 ) break;
      

      // Detrend and zero-mean xx[] and yy[].
      const bool detrend = true ;
      if ( detrend ) 
	{
	  coh_lremv( xx , nloaded );
	  coh_lremv( yy , nloaded );
	}

      // Apply Hann window
      for (int i = 0; i < nloaded; i++) 
	{
	  xx[i] *= weight[i];
	  yy[i] *= weight[i];
	}


      // Compute forward FFT
      coh_fft842(0, npfft, &(xx[0]), &(yy[0]) );
      
      // Compute auto- and cross-spectra. 
      gxx[0] += 4.0 * xx[0] * xx[0];
      gyy[0] += 4.0 * yy[0] * yy[0];
      gxyre[0] += 2.0 * xx[0] * yy[0];
      gxyim[0] = 0.0;
      
      for (int i = 1; i < nd2; i++) 
	{
	  double xi = xx[i];
	  double xj = xx[npfft-i];
	  double yi = yy[i];
	  double yj = yy[npfft-i];
	  
	  gxx[i] += (xi+xj)*(xi+xj) + (yi-yj)*(yi-yj);
	  gyy[i] += (yi+yj)*(yi+yj) + (xi-xj)*(xi-xj);
	  gxyre[i] += xi*yj + xj*yi;
	  gxyim[i] += xj*xj + yj*yj - xi*xi - yi*yi;
      }
      
      ++nffts;
      
    }


  //
  // Compile results
  //

  coh_t res(nd2);

  if (nffts == 0) return res;

  // Sample interval (seconds).
  double dt = 1.0/sampfreq;

  // Frequency interval (Hz). 
  double df = 1.0/(dt*npfft);

  // Normalize estimates. 
  double temp1 = sfx * dt / (4.0 * nnn * nffts);
  double temp2 = sfy * dt / (4.0 * nnn * nffts);
  double sf = sqrt(fabs(sfx*sfy));
  double temp3 = sf  * dt / (2.0 * nnn * nffts);
  double temp4 = sf  * dt / (4.0 * nnn * nffts);

  
  // final calculations
  
  for (int i = 0; i < nd2; i++) 
    {
      gxx[i] *= temp1;
      gyy[i] *= temp2;
      gxyre[i] *= temp3;
      gxyim[i] *= temp4;
      
      // Compute and print magnitude squared coherence (dimensionless), and
      // cross- and auto-spectra (in dB). 

      phi[i] = gxyre[i]*gxyre[i] + gxyim[i]*gxyim[i];
      
      //      std::cout << "phi " << i << " " << phi[i] << " " << gxyre[i] << " " << gxyim[i] << " " << gxx[i] << " " << gyy[i] << "\n";

      if (gxx[i] == 0.0 || gyy[i] == 0.0) res.coh[i] = 1.0;
      else res.coh[i] = phi[i] / (gxx[i]*gyy[i]);
      
      // truncate?
      res.cross_spectrum[i] = phi[i] > 1.0e-10 ? 5.0*log10(phi[i]) : -50.0 ;
      res.auto_spectrum1[i] = gxx[i] > 1.0e-10 ? 10.0*log10(gxx[i]) : -100.0 ;
      res.auto_spectrum2[i] = gyy[i] > 1.0e-10 ? 10.0*log10(gyy[i]) : -100.0 ;
      
      // frq
      res.frq[i] = df*i;
      
    }
  
  return res;
}
  

coh_t dsptools::coherence( edf_t & edf , const int signal1 , const int signal2 , const interval_t & interval , bool legacy )
{  
  
  //
  // (re)check sampling rate is equal
  //
  
  const int sr1 = edf.header.sampling_freq( signal1 );
  const int sr2 = edf.header.sampling_freq( signal2 );
  const int ns = sr1 < sr2 ? sr1 : sr2 ;  
  if ( sr1 != sr2 ) 
    {      
      if ( sr1 != ns ) resample_channel( edf , signal1 , ns );
      if ( sr2 != ns ) resample_channel( edf , signal2 , ns );
    }
  
  //
  // extract signals
  //
  
  slice_t slice1( edf , signal1 , interval );    
  slice_t slice2( edf , signal2 , interval );    
  
  const std::vector<double> * d1 = slice1.pdata();
  const std::vector<double> * d2 = slice2.pdata();

  if ( d1->size() != d2->size() ) Helper::halt( "internal error, signals different length in coherence()");

  //
  // compute coherence 
  //
  
  if ( legacy ) 
    return legacy_coherence( d1 , d2 , ns );
  else
    {
      
      const int Fs = ns;
      const double segment_sec = 5; 
      const double overlap_sec = 0;
      const bool average_adj = false;
      const bool detrend = false;

      coherence_t coherence( *d1, *d2 , Fs , 
			     segment_sec , overlap_sec , 
			     WINDOW_HANN , average_adj , detrend );

      return coherence.res;
    }

}


void dsptools::coh_lremv( std::vector<double> & x , const int n )
{
  
  // use n and not x.size(), as x may be zero-padded

  // DC component of data
  double dc = 0.0;
  
  // slope of data 
  double slope = 0;

  for (int i = 0; i < n; i++) 
    {
      dc += x[i];
      slope += x[i]*(i+1);
    }
  dc /= (double)n;
  
  slope *= 12.0/(n*(n*(double)n-1.0));
  slope -= 6.0*dc/(n-1.0);

  double fln = dc - 0.5*(n+1.0)*slope;

  for (int i = 0; i < n; i++)
    x[i] -= (i+1)*slope + fln;
  
}

void dsptools::coh_r2tx(int nthpo, double* cr0, double* cr1, double * ci0, double * ci1)
{
  int i;
  double temp;
  for (i = 0; i < nthpo; i += 2) {
    temp = cr0[i] + cr1[i]; cr1[i] = cr0[i] - cr1[i]; cr0[i] = temp;
    temp = ci0[i] + ci1[i]; ci1[i] = ci0[i] - ci1[i]; ci0[i] = temp;
  }

}

void dsptools::coh_r4tx(int nthpo, 
			double * cr0, 
			double * cr1, 
			double * cr2, 
			double * cr3, 
			double * ci0, 
			double * ci1, 
			double * ci2, 
			double * ci3)
{
  int i;
  double i1, i2, i3, i4, r1, r2, r3, r4;

  for (i = 0; i < nthpo; i += 4) {
    r1 = cr0[i] + cr2[i];
    r2 = cr0[i] - cr2[i];
    r3 = cr1[i] + cr3[i];
    r4 = cr1[i] - cr3[i];
    i1 = ci0[i] + ci2[i];
    i2 = ci0[i] - ci2[i];
    i3 = ci1[i] + ci3[i];
    i4 = ci1[i] - ci3[i];
    cr0[i] = r1 + r3;
    ci0[i] = i1 + i3;
    cr1[i] = r1 - r3;
    ci1[i] = i1 - i3;
    cr2[i] = r2 - i4;
    ci2[i] = i2 + r4;
    cr3[i] = r2 + i4;
    ci3[i] = i2 - r4;
  }
  
}


void dsptools::coh_r8tx(int nx, int nthpo, int length, 
			double * cr0, double * cr1, double * cr2, double * cr3, double * cr4, double * cr5, 
			double * cr6, double * cr7, double * ci0, double * ci1,
			double * ci2, double * ci3, double * ci4, double * ci5, double * ci6, double * ci7)
{

  double scale = 2.0*M_PI/length, arg, tr, ti;
  double c1, c2, c3, c4, c5, c6, c7;
  double s1, s2, s3, s4, s5, s6, s7;
  double ar0, ar1, ar2, ar3, ar4, ar5, ar6, ar7;
  double ai0, ai1, ai2, ai3, ai4, ai5, ai6, ai7;
  double br0, br1, br2, br3, br4, br5, br6, br7;
  double bi0, bi1, bi2, bi3, bi4, bi5, bi6, bi7;
  int j, k;

  for (j = 0; j < nx; j++) {
    arg = j*scale;
    c1 = cos(arg);
    s1 = sin(arg);
    c2 = c1*c1 - s1*s1;
    s2 = 2.0*c1*s1;
    c3 = c1*c2 - s1*s2;
    s3 = c2*s1 + s2*c1;
    c4 = c2*c2 - s2*s2;
    s4 = 2.0*c2*s2;
    c5 = c2*c3 - s2*s3;
    s5 = c3*s2 + s3*c2;
    c6 = c3*c3 - s3*s3;
    s6 = 2.0*c3*s3;
    c7 = c3*c4 - s3*s4;
    s7 = c4*s3 + s4*c3;
    for (k = j; k < nthpo; k += length) {
      ar0 = cr0[k] + cr4[k];      ar4 = cr0[k] - cr4[k];
      ar1 = cr1[k] + cr5[k];      ar5 = cr1[k] - cr5[k];
      ar2 = cr2[k] + cr6[k];      ar6 = cr2[k] - cr6[k];
      ar3 = cr3[k] + cr7[k];      ar7 = cr3[k] - cr7[k];

      ai0 = ci0[k] + ci4[k];      ai4 = ci0[k] - ci4[k];
      ai1 = ci1[k] + ci5[k];      ai5 = ci1[k] - ci5[k];
      ai2 = ci2[k] + ci6[k];      ai6 = ci2[k] - ci6[k];
      ai3 = ci3[k] + ci7[k];      ai7 = ci3[k] - ci7[k];

      br0 = ar0 + ar2;            br2 = ar0 - ar2;
      br1 = ar1 + ar3;            br3 = ar1 - ar3;

      br4 = ar4 - ai6;            br6 = ar4 + ai6;
      br5 = ar5 - ai7;            br7 = ar5 + ai7;

      bi0 = ai0 + ai2;            bi2 = ai0 - ai2;
      bi1 = ai1 + ai3;            bi3 = ai1 - ai3;

      bi4 = ai4 + ar6;            bi6 = ai4 - ar6;
      bi5 = ai5 + ar7;            bi7 = ai5 - ar7;


      cr0[k] = br0 + br1;
      ci0[k] = bi0 + bi1;
      if (j > 0) {
	cr1[k] = c4*(br0-br1) - s4*(bi0-bi1);
	ci1[k] = c4*(bi0-bi1) + s4*(br0-br1);
	cr2[k] = c2*(br2-bi3) - s2*(bi2+br3);
	ci2[k] = c2*(bi2+br3) + s2*(br2-bi3);
	cr3[k] = c6*(br2+bi3) - s6*(bi2-br3);
	ci3[k] = c6*(bi2-br3) + s6*(br2+bi3);
	tr = M_SQRT1_2*(br5-bi5);
	ti = M_SQRT1_2*(br5+bi5);
	cr4[k] = c1*(br4+tr) - s1*(bi4+ti);
	ci4[k] = c1*(bi4+ti) + s1*(br4+tr);
	cr5[k] = c5*(br4-tr) - s5*(bi4-ti);
	ci5[k] = c5*(bi4-ti) + s5*(br4-tr);
	tr = -M_SQRT1_2*(br7+bi7);
	ti =  M_SQRT1_2*(br7-bi7);
	cr6[k] = c3*(br6+tr) - s3*(bi6+ti);
	ci6[k] = c3*(bi6+ti) + s3*(br6+tr);
	cr7[k] = c7*(br6-tr) - s7*(bi6-ti);
	ci7[k] = c7*(bi6-ti) + s7*(br6-tr);
      }
      else {
	cr1[k] = br0 - br1;
	ci1[k] = bi0 - bi1;
	cr2[k] = br2 - bi3;
	ci2[k] = bi2 + br3;
	cr3[k] = br2 + bi3;
	ci3[k] = bi2 - br3;
	tr = M_SQRT1_2*(br5 - bi5);
	ti = M_SQRT1_2*(br5 + bi5);
	cr4[k] = br4 + tr;
	ci4[k] = bi4 + ti;
	cr5[k] = br4 - tr;
	ci5[k] = bi4 - ti;
	tr = -M_SQRT1_2*(br7 + bi7);
	ti =  M_SQRT1_2*(br7 - bi7);
	cr6[k] = br6 + tr;
	ci6[k] = bi6 + ti;
	cr7[k] = br6 - tr;
	ci7[k] = bi6 - ti;
      }
    }
  }


}


void dsptools::coh_fft842(int in, int n, double * x, double * y)
{

  double temp;
  
  int i, j, ij, j1, j2, j3, j4, j5, j6, j7, j8, j9, j10, j11, j12, j13, j14,
    ji, l[15], nt, nx, n2pow, n8pow;
  
  for (n2pow = nt = 1; n2pow <= 15 && n > nt; n2pow++)
    nt <<= 1;
  
  n2pow--;
  if (n != nt) {
    (void)fprintf(stderr, "fft842: %d is not a power of 2\n", n);
    exit(2);
  }

  n8pow = n2pow/3;
  if (in == 0) {
    for (i = 0; i < n; i++)
      y[i] = -y[i];
  }
  
  /* Do radix 8 passes, if any. */
  for (i = 1; i <= n8pow; i++) {
    nx = 1 << (n2pow - 3*i);
    coh_r8tx(nx, n, 8*nx,
	     x, x+nx, x+2*nx, x+3*nx, x+4*nx, x+5*nx, x+6*nx, x+7*nx,
	     y, y+nx, y+2*nx, y+3*nx, y+4*nx, y+5*nx, y+6*nx, y+7*nx);
  }
  
  /* Do final radix 2 or radix 4 pass. */
  switch (n2pow - 3*n8pow) {
  case 0:   break;
  case 1:   coh_r2tx(n, x, x+1, y, y+1); break;
  case 2:   coh_r4tx(n, x, x+1, x+2, x+3, y, y+1, y+2, y+3); break;
  }
  
  for (j = 0; j < 15; j++) {
    if (j <= n2pow) l[j] = 1 << (n2pow - j);
    else l[j] = 1;
  }   
  ij = 0;
  for (j1 = 0; j1 < l[14]; j1++)
    for (j2 = j1; j2 < l[13]; j2 += l[14])
      for (j3 = j2; j3 < l[12]; j3 += l[13])
	for (j4 = j3; j4 < l[11]; j4 += l[12])
	  for (j5 = j4; j5 < l[10]; j5 += l[11])
	    for (j6 = j5; j6 < l[9]; j6 += l[10])
	      for (j7 = j6; j7 < l[8]; j7 += l[9])
		for (j8 = j7; j8 < l[7]; j8 += l[8])
		  for (j9 = j8; j9 < l[6]; j9 += l[7])
		    for (j10 = j9; j10 < l[5]; j10 += l[6])
		      for (j11 = j10; j11 < l[4]; j11 += l[5])
			for (j12 = j11; j12 < l[3]; j12 += l[4])
			  for (j13 = j12; j13 < l[2]; j13 += l[3])
			    for (j14 = j13; j14 < l[1]; j14 += l[2])
			      for (ji = j14; ji < l[0]; ji += l[1]) {
				if (ij < ji) {
				  temp = x[ij]; x[ij] = x[ji]; x[ji] = temp;
				  temp = y[ij]; y[ij] = y[ji]; y[ji] = temp;
				}
				ij++;
			      }
  if (in == 0) {
    for (i = 0; i < n; i++)
      y[i] = -y[i];
  }
  else {
    for (i = 0; i < n; i++) {
      x[i] /= (double)n;
      y[i] /= (double)n;
    }
  }

}


