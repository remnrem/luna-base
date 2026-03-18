
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
#include "param.h"
#include "fftw/fftwrap.h"
#include "helper/helper.h"
#include "db/db.h"

#include "edf/edf.h"
#include "edf/slice.h"

#include <cmath>
#include <set>

extern logger_t logger;
extern writer_t writer;

// Direct implementation of Fourier domain DFA
// Nolte et al (2019) Scientific Reports

namespace {

bool dfa_linear_fit( const std::vector<double> & t ,
                     const std::vector<double> & fluctuations ,
                     const double t_lwr ,
                     const double t_upr ,
                     double * alpha ,
                     double * r2 )
{
  std::vector<double> x, y;

  for (int i = 0; i < t.size(); ++i)
    {
      if ( t[i] <= 0 || fluctuations[i] <= 0 ) continue;
      if ( t[i] < t_lwr || t[i] > t_upr ) continue;
      x.push_back( std::log( t[i] ) );
      y.push_back( std::log( fluctuations[i] ) );
    }

  const int n = x.size();
  if ( n < 2 ) return false;

  double mx = 0 , my = 0;
  for (int i = 0; i < n; ++i)
    {
      mx += x[i];
      my += y[i];
    }
  mx /= n;
  my /= n;

  double sxx = 0 , sxy = 0 , syy = 0;
  for (int i = 0; i < n; ++i)
    {
      const double dx = x[i] - mx;
      const double dy = y[i] - my;
      sxx += dx * dx;
      sxy += dx * dy;
      syy += dy * dy;
    }

  if ( sxx <= 0 ) return false;

  const double b = sxy / sxx;
  const double a = my - b * mx;

  double sse = 0;
  for (int i = 0; i < n; ++i)
    {
      const double yhat = a + b * x[i];
      const double err = y[i] - yhat;
      sse += err * err;
    }

  *alpha = b;
  *r2 = syy > 0 ? 1.0 - ( sse / syy ) : 1.0;
  return true;
}

}

void dsptools::dfa_wrapper( edf_t & edf , param_t & param )
{

  // get signals

  std::string signal_label = param.requires( "sig" );
  const bool no_annotations = true;
  signal_list_t signals = edf.header.signal_list( signal_label , no_annotations );
  const int ns = signals.size();

  // parameters

  const int wn = param.has( "n" ) ? param.requires_int( "n" ) : 100;
  const bool explicit_time_grid = param.has( "min" ) || param.has( "max" );
  const bool classical_grid = explicit_time_grid ? false : param.value( "grid" , true ) != "TIME";
  const int    scale = param.has( "m" ) ? param.requires_int( "m" ) : 2;
  const double wmin = param.has( "min" ) ? param.requires_dbl( "min" ) : 0.1;
  const double wmax = param.has( "max" ) ? param.requires_dbl( "max" ) :
                      wmin * pow( 10.0 , scale );
  const double alpha_lwr = param.has( "alpha-lwr" ) ? param.requires_dbl( "alpha-lwr" ) : 0;
  const double alpha_upr = param.has( "alpha-upr" ) ? param.requires_dbl( "alpha-upr" ) : 1e12;
 
  const bool narrowband = param.has( "f-lwr" );   
  const double fmin = narrowband ? param.requires_dbl( "f-lwr" ) : -1 ; 
  const double fmax = narrowband ? param.requires_dbl( "f-upr" ) : -1 ;
  const double ripple = param.has( "ripple" ) ? param.requires_dbl( "ripple" ) : 0.01;
  const double tw = param.has( "tw" ) ? param.requires_dbl( "tw" ) : 1;
  const bool envelope = param.has( "envelope" ) ? param.yesno( "envelope" ) : true;
  const bool by_epoch = param.yesno( "epoch" );
  
  logger << "  DFA parameters\n"
	 << "     n (points) = " << wn << "\n"
	 << "     grid       = " << ( classical_grid ? "classical" : "time" ) << "\n";

  if ( classical_grid )
    logger << "     points     = " << wn << " between 4 and T/4 samples\n";
  else
    logger << "     min        = " << wmin << " sec\n"
	   << "     max        = " << wmax << " sec\n"
	   << "     m (scale)  = " << scale << "\n";

  logger
	 << "     alpha-lwr  = " << alpha_lwr << " sec\n"
	 << "     alpha-upr  = " << alpha_upr << " sec\n"
	 << "     epoch      = " << ( by_epoch ? "T" : "F" ) << "\n";
  
  if ( narrowband )
    logger << "  applying narrowband filter:\n"
	   << "    f-lwr    = " << fmin << "\n"
	   << "    f-upr    = " << fmax << "\n"
	   << "    ripple   = " << ripple << "\n"
	   << "    tw       = " << tw << "\n"
	   << "    envelope = " << ( envelope ? "T" : "F" ) << "\n";
  
  
  //
  // iterate over each signal
  //

  logger << "  processing:";
  
  for (int s=0; s<ns; s++)
    {

      logger << " " << signals.label(s);

      writer.level( signals.label(s) , globals::signal_strat );

      
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
      dfa.filter_hilbert( fmin, fmax, ripple, tw , envelope );
      
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

	  if ( classical_grid )
	    dfa.set_windows_classical( Fs , d->size() , wn );
	  else
	    dfa.set_windows_seconds( Fs , wmin , wmax , wn );

	  dfa.proc( d );	        
	  
	  // output

	  const int nw = dfa.w.size();
	  double alpha = 0 , r2 = 0;
	  if ( dfa_linear_fit( dfa.t , dfa.fluctuations , alpha_lwr , alpha_upr , &alpha , &r2 ) )
	    {
	      writer.value( "ALPHA" , alpha );
	      writer.value( "R2" , r2 );
	    }

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

  writer.unlevel( globals::signal_strat );

  logger << "\n";
}

dfa_t::dfa_t()
{
  flwr = -1;
  fupr = -1;
  ripple = -1;
  tw = -1;
}


void dfa_t::set_windows( double sr1 , double l , int mexp , int c )
{
  set_windows_seconds( sr1 , l , l * pow( 10.0 , mexp ) , c );
}


void dfa_t::set_windows_seconds( double sr1 , double lwr_sec , double upr_sec , int c )
{
  sr = sr1;
  if ( c < 2 ) Helper::halt( "bad DFA values" );
  if ( lwr_sec <= 0 || upr_sec <= 0 || upr_sec <= lwr_sec ) Helper::halt( "bad DFA min/max values" );

  w.resize( c );
  t.resize( c );
  
  // log-spaced timescale grid in seconds
      
  for (int i=0; i<c; i++)
    {
      t[ i ] = lwr_sec * pow( upr_sec / lwr_sec , i / (double)( c - 1 ) );
      w[ i ] = sr * t[ i ];
    }
}


void dfa_t::set_windows_classical( double sr1 , int np , int c )
{
  sr = sr1;
  if ( c < 2 ) Helper::halt( "bad DFA values" );
  if ( np < 16 ) Helper::halt( "bad DFA trace length for classical grid" );

  const int lwr = 4;
  const int upr = np / 4;

  if ( upr < lwr ) Helper::halt( "bad DFA trace length for classical grid" );

  std::set<int> windows;

  if ( upr == lwr )
    windows.insert( lwr );
  else
    for (int i=0; i<c; ++i)
      {
        const double pos = i / (double)( c - 1 );
        const double wl = lwr * pow( upr / (double)lwr , pos );
        windows.insert( (int)std::floor( wl + 0.5 ) );
      }

  if ( windows.size() < 2 )
    {
      windows.insert( lwr );
      windows.insert( upr );
    }

  if ( windows.size() < 2 ) Helper::halt( "bad DFA classical scale grid" );

  w.clear();
  t.clear();
  w.reserve( windows.size() );
  t.reserve( windows.size() );

  for (std::set<int>::const_iterator ii = windows.begin(); ii != windows.end(); ++ii)
    {
      w.push_back( *ii );
      t.push_back( *ii / sr );
    }
}

void dfa_t::proc( const std::vector<double> * d )
{

  //
  // step 0 : initialize
  //
  
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

      // use either the envelope 
      if ( envelope ) 
	d0 = *(hilbert.magnitude());
      else // or the filtered signal
	d0 = *(hilbert.signal());      
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
