
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

#include "hilbert.h"

#include "fftw/fftwrap.h"
#include "miscmath/crandom.h"
#include "miscmath/miscmath.h"
#include "dsp/fir.h"
#include "defs/defs.h"

#include <iostream>
#include <cmath>




hilbert_t::hilbert_t( const std::vector<double> & d ) : input(d)
{
  // this mode assumes we've already BPF the input
  proc();
}


hilbert_t::hilbert_t( const std::vector<double> & d , const int sr , double lwr , double upr , double ripple , double tw )
{
  // include band-pass filter
  input = dsptools::apply_fir( d , sr , fir_t::BAND_PASS , ripple , tw , lwr , upr );
  
  proc();
}
      
void hilbert_t::proc()
{

  int n = input.size();
  
  // 1) take FFT
  FFT fft( n , 1 , FFT_FORWARD );
  fft.apply( input );
  std::vector<dcomp> f = fft.transform();
  if ( f.size() != n ) Helper::halt( "internal error in hilbert()" );

  // 2) Adjusted postive/negative frequencies
  
  int pos_idx = floor(n/2.0) + ( n % 2 ) - 1;
  int neg_idx = ceil(n/2.0) + ( ! ( n % 2 ) );
  
  for (int i = 1 ; i <= pos_idx ; i++ ) f[i] *= 2;
  for (int i = neg_idx ; i < n ; i++ ) f[i] = 0;

  // equivalent to rotating Fourier coefficients by computing the
  // iAsin(2pft) component, i.e., the phase quadrature) positive
  // frequencies are rotated counter-clockwise; negative frequencies
  // are rotated clockwise)

  // f(posF) = f(posF) + -1i*complexf(posF);
  // f(negF) = f(negF) +  1i*complexf(negF);
  
  
  // 3) Inverse FFT of the rotated coefficients 
  FFT ifft( n , 1 , FFT_INVERSE );
  ifft.apply( f );
  std::vector<dcomp> ht = ifft.scaled_transform();
  
  if ( ht.size() != n ) Helper::halt( "problem in hilbert()" );

  // 4) Store phase, magnitude

  ph.resize( n );
  mag.resize( n );

  for(int i=0;i<n;i++)
    {
      double a = std::real( ht[i] ) ;
      double b = std::imag( ht[i] ) ;
      ph[i] = atan2( b , a );
      mag[i] = sqrt( a*a + b*b );     
    }
}

const std::vector<double> * hilbert_t::phase() const
{
  return & ph;
}

const std::vector<double> * hilbert_t::magnitude() const
{
  return & mag;
}

const std::vector<double> * hilbert_t::signal() const
{
  // band-pass filtered version
  return & input;
}

std::vector<double> hilbert_t::instantaneous_frequency( double Fs ) const
{
  std::vector<double> angles = ph; 
  unwrap( &angles );

  const int nm1 = angles.size() - 1;
  std::vector<double> f( nm1 );
  for (int i=0;i<nm1;i++)
    f[i] = Fs / ( 2.0 * M_PI ) * ( angles[i+1] - angles[i] ) ;
  return f;
        
}


void hilbert_t::bin( double p , int bs , std::vector<int> * acc ) const
{
  int a = floor( MiscMath::as_angle_0_pos2neg( p ) );
  int b = a / bs;
  
  if( b < 0 || b >= acc->size() ) 
    {      
      std::cerr << "p, a,b " << p << " " << MiscMath::as_angle_0_pos2neg( p ) << " " << a << " " << b << " " << acc->size() << "\n";
      Helper::halt( "internal error in hilbert_t::bin() " );
    }
  (*acc)[b]++;
}



itpc_t::itpc_t( const int ne , const int nbins )
{
  if ( 360 % nbins ) Helper::halt( "number of bins must imply integer number of degrees per bin" );    
  phase.resize( ne , 0 );
  event_included.resize( ne , false );
  phasebin.resize( nbins );
}

itpc_t hilbert_t::phase_events( const std::vector<int> & e , const std::vector<bool> * mask , 
				const int nreps , 
				const int sr , 
				const double epoch_sec 
				) const
{

  // given a set of 'events' (i.e. spindles) in 'e' defined by sample-point
 
  // if we are given a 'mask', only consider events that fall within this mask 
  // when permuting, mask stays w/ hilbert phase (i.e. typically indicates whether a SO is present)
  // but events get permuted

  // if nreps > 0 then use permutation to get empirical distribution of key metrics

  // return itpc_t, which contains the phase for each event, 
  //  - the ITPC statistic
  //  - counts of how many times 
  // calculate 'inter-trial phase clustering' metric (i.e. measure of
  // consistency

  const int n = e.size();

  const int nbins = 18; // 18 * 20-degree bins
  const int binsize = 360 / nbins;
 
  itpc_t itpc( n , nbins );

  // the maximum expected sample-point (+1)
  const int mx = ph.size();
  
  if ( mask != NULL && mask->size() != mx ) 
    Helper::halt( "internal error in hilbert_t::phase_events()" );

  //
  // shuffle only /within/ epoch, if epochs are defined?  If so, for
  // each event, we need to find the start and stop sp for that epoch,
  // and shuffle only within those bounds (with wrapping)
  //
  
  int es = sr * epoch_sec;

  // from e[], also get eoffset[], which is the relative within-epoch position
  std::vector<int> eoffset( n ) ;

  if ( es ) 
    {
      for (int i=0; i<n; i++) 
	eoffset[i] = e[i] % es;
    }

  
  //
  // Observed events
  //
  
  int counted = 0;
  
  std::vector<int> pbacc( nbins , 0 ); // phase-bin accumulator 

  dcomp s(0,0);

  for (int i=0;i<n;i++)
    {      

      // check range
      if ( e[i] >= mx || e[i] < 0 ) 
	Helper::halt( "problem requesting value outside range in hilbert()" );

      // if hilbert transform to be included? 
      if ( mask == NULL || (*mask)[ e[i] ] ) 
	{
	  
	  // track phase
	  itpc.phase[i] = ph[ e[i] ];      

	  // phase bin
	  bin( itpc.phase[i] , binsize , &pbacc );
	  
	  // accumulate ITPC
	  s += exp( dcomp(0,itpc.phase[i]) );

	  // and record that this event was included
	  ++counted;
	  
	  itpc.event_included[i] = true; 

	}      

    }


  // normalise ITPC and return
  s /= double(counted);
  
  // record ITPC etc
  itpc.itpc.obs = abs( s );
  itpc.ninc.obs  = counted;
  itpc.pv.obs    = exp( -counted * itpc.itpc.obs * itpc.itpc.obs );   
  itpc.sig.obs   = itpc.pv.obs < 0.05;
  itpc.angle.obs = MiscMath::as_angle_0_pos2neg( arg( s ) );
  for (int b=0; b < nbins; b++) itpc.phasebin[b].obs = pbacc[b];

  if ( nreps == 0 ) return itpc;


  //
  // Permutations
  //

  // max shuffle, either within epoch, or whole signal

  const int maxshuffle = es ? es : mx ;
  
  for (int r = 0 ; r < nreps ; r++ ) 
    {
      
      // get permutation shift
      int pp = CRandom::rand( maxshuffle );
      
      std::vector<int> pbacc( nbins , 0 ); // phase-bin accumulator 
      
      int counted = 0 ; 
      
      dcomp s(0,0);
      
      for (int i=0;i<n;i++)
	{      
	  
	  // get permuted position, and wrap if needed
	  int pei = e[i] + pp ; 
	  
	  // check for wrapping
	  if ( es == 0 ) // whole-signal shuffle
	    {
	      if ( pei >= maxshuffle ) pei -= maxshuffle;
	    }
	  else // within-epoch shuffle
	    {
	      if ( eoffset[i] + pp >= maxshuffle ) pei -= maxshuffle;
	    }

	  // if hilbert transform to be included? 

	  if ( mask == NULL || (*mask)[ pei ] ) 
	    {
	      
	      // track phase
	      itpc.phase[i] = ph[ pei ];      
	  
	      // phase bin
	      bin( itpc.phase[i] , binsize , &pbacc );

	      // accumulate ITPC
	      s += exp( dcomp(0,itpc.phase[i]) );

	      // and record that this event was included
	      ++counted;
	      
	      itpc.event_included[i] = true; 
	      
	    }      
	  
	}
      
      // normalise ITPC and return
      s /= double(counted);
  
      // record ITPC etc
      itpc.itpc.perm.push_back( abs( s ) );
      itpc.ninc.perm.push_back( counted );
      double pv = exp( -counted * itpc.itpc.obs * itpc.itpc.obs ) ;
      itpc.pv.perm.push_back( pv );
      itpc.sig.perm.push_back( pv < 0.05 );
      itpc.angle.perm.push_back( MiscMath::as_angle_0_pos2neg( arg( s ) ) );
      for (int b=0; b < nbins; b++) itpc.phasebin[b].perm.push_back ( pbacc[b] );

    }

  // get empirical p-values
  itpc.itpc.calc_stats();
  itpc.ninc.calc_stats();
  itpc.sig.calc_stats();
  for (int b=0; b < nbins; b++) itpc.phasebin[b].calc_stats();

  
  //
  // All done
  //
  
  return itpc;

}



void hilbert_t::unwrap( std::vector<double> * p ) const 
{
  
  // http://homepages.cae.wisc.edu/~brodskye/mr/phaseunwrap/unwrap.c

  const int n = p->size();
  
  std::vector<double> dp( n );
  std::vector<double> dps( n );
  std::vector<double> dp_corr( n );
  std::vector<double> cumsum( n );

  // default tol in matlab
  double cutoff = M_PI;   
  
  // incremental phase variation 
  // MATLAB: dp = diff(p, 1, 1);
  for (int j = 0; j < n-1; j++)
    dp[j] = (*p)[j+1] - (*p)[j];
  
  // equivalent phase variation in [-pi, pi]
  // MATLAB: dps = mod(dp+dp,2*pi) - pi;
  for (int j = 0; j < n-1; j++)
    dps[j] = (dp[j]+M_PI) - floor((dp[j]+M_PI) / (2*M_PI))*(2*M_PI) - M_PI;
  
  // preserve variation sign for +pi vs. -pi
  // MATLAB: dps(dps==pi & dp>0,:) = pi;
  for (int j = 0; j < n-1; j++)
    if ((dps[j] == -M_PI) && (dp[j] > 0))
      dps[j] = M_PI;

  // incremental phase correction
  // MATLAB: dp_corr = dps - dp;
  for (int j = 0; j < n-1; j++)
    dp_corr[j] = dps[j] - dp[j];
      
  // Ignore correction when incremental variation is smaller than cutoff
  // MATLAB: dp_corr(abs(dp)<cutoff,:) = 0;
  for (int j = 0; j < n-1; j++)
    if (fabs(dp[j]) < cutoff)
      dp_corr[j] = 0;

  // Find cumulative sum of deltas
  // MATLAB: cumsum = cumsum(dp_corr, 1);
  cumsum[0] = dp_corr[0];
  for (int j = 1; j < n-1; j++)
    cumsum[j] = cumsum[j-1] + dp_corr[j];

  // Integrate corrections and add to P to produce smoothed phase values
  // MATLAB: p(2:m,:) = p(2:m,:) + cumsum(dp_corr,1);
  for (int j = 1; j < n; j++)
    (*p)[j] += cumsum[j-1];
}
