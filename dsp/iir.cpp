
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

#include "dsp/iir.h"
#include "dsp/filter.h"

#include "helper/helper.h"
#include "helper/logger.h"
#include "fftw/fftwrap.h"
#include "db/db.h"
#include "eval.h"

#include "edf/slice.h"
#include "edf/edf.h"
#include "dsp/ngaus.h"

extern logger_t logger;

extern writer_t writer;


Eigen::MatrixXd dsptools::butterworth( const Eigen::MatrixXd & X , int order , int fs, double f1, double f2 )
{

  const int r = X.rows();
  const int c = X.cols();
  Eigen::MatrixXd F = Eigen::MatrixXd::Zero( r , c );  
  
  for (int j=0; j<c; j++)
    {
      iir_t iir;
      iir.init( BUTTERWORTH_BANDPASS , order, fs , f1, f2 );
      F.col(j) = iir.apply( X.col(j) );
    }

  return F;
  
}

void dsptools::apply_iir( edf_t & edf , param_t & param )
{

  
  if ( param.has( "butterworth" ) == param.has( "chebyshev" ) )
    Helper::halt( "IIR requires either butterworth or chebyshev" );

  const bool butterworth = param.has( "butterworth" );
  const bool low_pass = param.has( "lowpass" );
  const bool high_pass = param.has( "highpass" );
  const bool band_pass = param.has( "bandpass" );
  const bool band_stop = param.has( "bandstop" );

  if ( low_pass + high_pass + band_pass + band_stop != 1 )
    Helper::halt( "IIR requires one of lowpass, lowpass, bandpass or bandstop" );
  
  const std::vector<double> p;

  const std::vector<double> p0 = butterworth ? param.dblvector( "butterworth" ) : param.dblvector( "chebyshev" );

  if ( butterworth )
    {
      if ( p0.size() != 1 )
	Helper::halt( "expecting butterworth=<order>" );
    }
  else
    {
      if ( p0.size() != 2 )
        Helper::halt( "expecting chebyshev=<order>,<eps>" );
    }
  
  const int order = p0[0];

  const double ceps = butterworth ? 0 : p0[1];

  //
  // get signals
  //

  const bool no_annotations = true;
  
  const std::string signal_label = param.requires( "sig" );
  
  signal_list_t signals = edf.header.signal_list( signal_label , no_annotations );
  
  const int ns = signals.size();

  for (int s=0; s<ns; s++)
    {
      logger << "  filtering " << signals.label(s) << " with " << order << "-order " << ( butterworth ? "Butterworth" : "Chebyshev" ) << " ";
      if ( ! butterworth ) logger << "(eps=" << ceps << ") ";
      if ( low_pass ) logger << "low-pass";
      else if ( high_pass ) logger << "high-pass";
      else if ( band_pass ) logger << "band-pass";
      else logger << "band-stop" ;
      logger << " IIR filter\n";
      
      const int fs = edf.header.sampling_freq( signals(s) );
      
      iir_t iir;
      
      iir_type_t iir_type;
      
      if ( butterworth )
	{
	  if ( low_pass )
	    {
	      const	std::vector<double> p = param.dblvector( "lowpass" );
	      if ( p.size() != 1 ) Helper::halt( "expecting lowpass=<frq>" );
	      iir.init( BUTTERWORTH_LOWPASS , order, fs , p[0] );
	    }
	  else if ( high_pass )
	    {
	      const	std::vector<double> p = param.dblvector( "highpass" );
	      if ( p.size() != 1 ) Helper::halt( "expecting highpass=<frq>" );
	      iir.init( BUTTERWORTH_HIGHPASS , order, fs , p[0] );
	    }
	  else if ( band_pass )
	    {
	      const	std::vector<double> p = param.dblvector( "bandpass" );
	      if ( p.size() != 2 ) Helper::halt( "expecting lowpass=<frq>,<frq>" );
	      iir.init( BUTTERWORTH_BANDPASS , order, fs , p[0] , p[1] );
	    }
	  else if ( band_stop )
	    {
	      const	std::vector<double> p = param.dblvector( "bandstop" );
	      if ( p.size() != 2 ) Helper::halt( "expecting bandstops=<frq>,<frq>" );
	      iir.init( BUTTERWORTH_BANDSTOP , order, fs , p[0] , p[1] );
	    }
	}
      else
	{
	  if ( low_pass )
	    {
	      const	std::vector<double> p = param.dblvector( "lowpass" );
	      if ( p.size() != 1 ) Helper::halt( "expecting lowpass=<frq>" );
	      iir.init( CHEBYSHEV_LOWPASS , order, ceps, fs , p[0] );
	    }
	  else if ( high_pass )
	    {
	      const	std::vector<double> p = param.dblvector( "highpass" );
	      if ( p.size() != 1 ) Helper::halt( "expecting highpass=<frq>" );
	      iir.init( CHEBYSHEV_HIGHPASS , order, ceps, fs , p[0]);
	    }
	  else if ( band_pass )
	    {
	      const	std::vector<double> p = param.dblvector( "bandpass" );
	      if ( p.size() != 2 ) Helper::halt( "expecting bandpass=<frq>,<frq>" );
	      iir.init( CHEBYSHEV_BANDPASS , order, ceps, fs , p[0] , p[1] );
	    }
	  else if ( band_stop )
	    {
	      const	std::vector<double> p = param.dblvector( "bandstop" );
	      if ( p.size() != 2 ) Helper::halt( "expecting bandstop=<frq>,<frq>" );
	      iir.init( CHEBYSHEV_BANDSTOP , order, ceps, fs , p[0] , p[1] );
	    }
	}    


      //
      // Get signal, filter and update
      //
      
      interval_t interval = edf.timeline.wholetrace();
      
      slice_t slice( edf , signals(s) , interval );
      
      const std::vector<double> * d = slice.pdata();
      
      std::vector<double> filtered = iir.apply( *d );
      
      edf.update_signal( signals(s) , &filtered );

      // next signal
    }

  // all done
}


// implements Butterworth and Chebyshev IIR filters

iir_t::iir_t()
{
  bwlp = NULL;
  bwhp = NULL;
  bwbp = NULL;
  bwbs = NULL;
  chelp = NULL;
  chehp = NULL;
  chebp = NULL;
  chebs = NULL;  
}

void iir_t::init( iir_type_t type , int order , double p1, double p2 , double p3 , double p4 )
{

  //                               p1  p2  p3   p4
  // type == BUTTERWORTH_    order sr  f1 (f2)
  // type == CHEBYSHEV       order eps sr  f1  (f2)

  // 
  // CHEBandPass* create_che_band_pass_filter(int order, FTR_PRECISION epsilon, FTR_PRECISION sampling_frequency, FTR_PRECISION lower_half_power_frequency, FTR_PRECISION upper_half_power_frequency);

  if ( type == BUTTERWORTH_LOWPASS )
    bwlp = create_bw_low_pass_filter( order , p1 , p2 );
  else if ( type == BUTTERWORTH_HIGHPASS )
    bwhp = create_bw_high_pass_filter( order , p1 , p2 );
  else if ( type == BUTTERWORTH_BANDPASS )
    bwbp = create_bw_band_pass_filter( order , p1 , p2 , p3 );
  else if ( type == BUTTERWORTH_BANDSTOP )
    bwbs = create_bw_band_stop_filter( order , p1 , p2 , p3 );
  else if ( type == CHEBYSHEV_LOWPASS )
    chelp = create_che_low_pass_filter( order , p1 , p2 , p3 );
  else if ( type == CHEBYSHEV_HIGHPASS )
    chehp = create_che_high_pass_filter( order , p1 , p2 , p3 );
  else if ( type == CHEBYSHEV_BANDPASS )
    chebp = create_che_band_pass_filter( order , p1 , p2 , p3 , p4 );
  else if ( type == CHEBYSHEV_BANDSTOP )
    chebs = create_che_band_stop_filter( order , p1 , p2 , p3 , p4 );


  
}


iir_t::~iir_t( )
{
  if ( bwlp != NULL ) free_bw_low_pass( bwlp );
  if ( bwhp != NULL ) free_bw_high_pass( bwhp );
  if ( bwbp != NULL ) free_bw_band_pass( bwbp );
  if ( bwbs != NULL ) free_bw_band_stop( bwbs );

  if ( chelp != NULL ) free_che_low_pass( chelp );
  if ( chehp != NULL ) free_che_high_pass( chehp );
  if ( chebp != NULL ) free_che_band_pass( chebp );
  if ( chebs != NULL ) free_che_band_stop( chebs );

}


std::vector<double> iir_t::apply( const std::vector<double> & x )
{
  const int n = x.size();
  std::vector<double> y( n , 0 );
  
  if ( bwlp != NULL )
    {
      for (int i=0; i<n; i++) y[i] = bw_low_pass( bwlp , x[i] );
      return y;
    }

  if ( bwhp != NULL )
    {
      for (int i=0; i<n; i++) y[i] = bw_high_pass( bwhp , x[i] );
      return y;
    }
    
  if ( bwbp != NULL )
    {
      for (int i=0; i<n; i++) y[i] = bw_band_pass( bwbp , x[i] );
      return y;
    }
    
  if ( bwbs != NULL )
    {
      for (int i=0; i<n; i++) y[i] = bw_band_stop( bwbs , x[i] );
      return y;
    }

  
  if ( chelp != NULL )
    {
      for (int i=0; i<n; i++) y[i] = che_low_pass( chelp , x[i] );
      return y;
    }

  if ( chehp != NULL )
    {
      for (int i=0; i<n; i++) y[i] = che_high_pass( chehp , x[i] );
      return y;
    }
    
  if ( chebp != NULL )
    {
      for (int i=0; i<n; i++) y[i] = che_band_pass( chebp , x[i] );
      return y;
    }
    
  if ( chebs != NULL )
    {
      for (int i=0; i<n; i++) y[i] = che_band_stop( chebs , x[i] );
      return y;
    }
  
  return y;
  
}



Eigen::VectorXd iir_t::apply( const Eigen::VectorXd & x )
{
  const int n = x.size();

  Eigen::VectorXd y = Eigen::VectorXd::Zero( n );
  
  if ( bwbp != NULL )
    {
      for (int i=0; i<n; i++) y[i] = bw_band_pass( bwbp , x[i] );
      return y;
    }
  else
    Helper::halt( "internal Eigen BWBP error" );
    
  return y;
  
}
