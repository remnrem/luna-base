
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


#include "dsp/dfa.h"
#include "fftw/fftwrap.h"
#include "helper/helper.h"
#include "db/db.h"

#include "edf/edf.h"
#include "edf/slice.h"

extern logger_t logger;
extern writer_t writer;

// Direct implementation of Fourier domain DFA
// Nolte et al (2019) Scientific Reports

void dsptools::dfa_wrapper( edf_t & edf , param_t & param )
{

  // get signals

  std::string signal_label = param.requires( "sig" );
  const bool no_annotations = true;
  signal_list_t signals = edf.header.signal_list( signal_label , no_annotations );
  const int ns = signals.size();

  // parameters

  const int wn = param.has( "n" ) ? param.requires_int( "n" ) : 100;
  const double wmin = param.has( "min" ) ? param.requires_dbl( "min" ) : 0.1;
  const double wmax = param.has( "max" ) ? param.requires_dbl( "max" ) : 10;

  const double fmin = param.requires_dbl( "f-lwr" );
  const double fmax = param.requires_dbl( "f-upr" );
  const double ripple = param.has( "ripple" ) ? param.requires_dbl( "ripple" ) : 0.02;
  const double tw = param.has( "tw" ) ? param.requires_dbl( "tw" ) : 0.5;

  const bool by_epoch = param.yesno( "epoch" );
  
  //
  // iterate over each signal
  //
  
  for (int s=0; s<ns; s++)
    {
      logger << "  processing " << signals.label(s)
	     << " for "
	     << fmin << " - " << fmax << " Hz\n";

      //
      // start iterating over epochs
      //
      
      int ne = by_epoch ? edf.timeline.first_epoch() : 0 ;
      if ( by_epoch && ne == 0 ) return;
      
      //
      // Set up DFA params
      //

      const double Fs = edf.header.sampling_freq( signals(s) );       
      dfa_t dfa;
      dfa.set_windows( Fs , wmin );
      dfa.filter_hilbert( fmin, fmax, ripple, tw );

      //
      // track epoch level stats
      //

      while ( 1 )
	{

	  int epoch = by_epoch ? edf.timeline.next_epoch() : 0 ; 

          if ( epoch == -1 ) break;

	  if ( by_epoch )
            writer.epoch( edf.timeline.display_epoch( epoch ) );
	  
	  // epoch or whole trace
          interval_t interval = by_epoch ? edf.timeline.epoch( epoch ) : edf.timeline.wholetrace();
	  
          slice_t slice( edf , signals(s) , interval );
	  
	  const std::vector<double> * d = slice.pdata();
	  
	  // perform DFA

	  dfa.proc( d );	        
	  
	  // output

	  const int nw = dfa.w.size();

	  for (int i=0; i<nw; i++)
	    {
	      writer.level( dfa.t[i] , globals::sec_strat );
	      writer.value( "FLUCT" , dfa.fluctuations[i] );
	      writer.value( "SLOPE" , dfa.slopes[i] );
	    }
	  writer.unlevel( globals::sec_strat );

	  // if whole-trace mode, break now
	  if ( ! by_epoch ) break;

	}

      if ( by_epoch )
	writer.unepoch();
	  
      // next signal
    }
  
}

dfa_t::dfa_t()
{
  flwr = -1;
  fupr = -1;
  ripple = -1;
  tw = -1;
}


void dfa_t::set_windows( double sr1 , double l, int mexp , int c )
{
  sr = sr1;
  if ( c < 2 ) Helper::halt( "bad DFA values" );
  if ( mexp < 2 ) Helper::halt( "bad DFA values" ); 
  if ( l <= 0 ) Helper::halt( "bad wmin and wmax values" ); 

  w.resize( c );
  t.resize( c );
  
  // add on exponential scale up to ^2 (0.1 - 10s)
      
  for (int i=0; i<c; i++)
    {
      t[ i ] = l * pow( 10 , i / (double)(c-1) * mexp );
      w[ i ] = sr * t[ i ];
    }
}

void dfa_t::proc( const std::vector<double> * d )
{
  
  // d - original data vector  
  
  // number of samples
  const int n = d->size();
  
  // number of windows
  const int nw = w.size();
  
  const bool boxcar = true;

  //
  // step 1 : absolute amplitude from Hilbert transform
  //

  std::vector<double> d0 = *d;
  
  if ( flwr > 0 && fupr > flwr )
    {
      // filter-Hilbert
      hilbert_t hilbert( d0 , sr , flwr , fupr , ripple , tw );      
      d0 = *(hilbert.magnitude());
      for (int i=0; i<d0.size(); i++) d0[i] = abs( d0[i] );
    }
  

  //
  // step 2 : Fourier-based DFA on this signal
  //

  //
  // remove mean from signal
  //
  
  MiscMath::centre( &d0 );

  // nx - size of FFT , depends on odd/even

  //
  // FFT
  //

  int index_start = 0;
  
  FFT fftseg( n , n , sr , FFT_FORWARD , WINDOW_NONE );  
  fftseg.apply( &(d0[0]) , n );
  
  int nx = fftseg.cutoff;  
  std::vector<std::complex<double> > t = fftseg.transform();  

  const int np = nx - 1;
  
  Eigen::ArrayXd p = Eigen::ArrayXd::Zero( np );  
  for (int i=1; i<nx; i++)
    p[i-1] = 2 * pow( std::abs( t[i] ) , 2 );
  // half last element?
  if ( d->size() % 2 == 0 ) 
    p[np-1] /= 2;
  
  //
  // boxcar method implementation
  //

  Eigen::ArrayXd ff = Eigen::ArrayXd::Zero( np );
  for (int i=0; i<np; i++) ff[i] = i+1;
  
  Eigen::ArrayXd g1 = sin( M_PI * ff / (double)n );

  fluctuations.clear();
  slopes.clear();
  
  for (int k=0; k<nw; k++)
    {

      const double wl = w[k];
	    
      Eigen::ArrayXd hsin = sin( M_PI * ff * wl/(double)n ) ;
      Eigen::ArrayXd hcos = cos( M_PI * ff * wl/(double)n ) ;

      Eigen::ArrayXd hx = ( 1 - hsin / ( wl * g1 ) );
      Eigen::ArrayXd h = hx / ( 2 * g1 );
      Eigen::ArrayXd h2 = h.pow(2);
      double F2 = ( h2 * p ).sum();      
      fluctuations.push_back( sqrt(F2) / (double)n );

      Eigen::ArrayXd hy = -hx * ( hcos * M_PI * ff/(double)n - hsin/wl ) / ( wl*g1 );
      Eigen::ArrayXd h3 = hy / ( 4 * g1.pow(2) );
      slopes.push_back(  ( h3 * p ).sum() / F2 * wl );
      
    }

}


