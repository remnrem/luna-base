
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

#include "dsp/gc.h"

#include "param.h"
#include "edf/slice.h"
#include "edf/edf.h"
#include <iostream>
#include <vector>

#include "stats/eigen_ops.h"
#include "db/db.h"
#include "helper/logger.h"

std::map<int,std::map<int,double> > gc_t::y2x_sum;
std::map<int,std::map<int,double> > gc_t::x2y_sum;
std::map<int,std::map<int,std::map<double,double> > > gc_t::tf_x2y_sum;
std::map<int,std::map<int,std::map<double,double> > > gc_t::tf_y2x_sum;
int gc_t::ne;

extern logger_t logger;
extern writer_t writer;

void gc_wrapper( edf_t & edf , param_t & param )
{

  //
  // Get signals
  //

  const bool no_annotations = true;

  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) , no_annotations );
  
  if ( signals.size() < 2 ) return;
  
  const int ns = signals.size();
  

  //
  // Check sample rates
  //

  std::vector<double> Fs = edf.header.sampling_freq( signals );

  int sr = Fs[0];
  for (int s=1;s<ns;s++)
    if ( Fs[s] != sr )
      Helper::halt( "all sampling rates must be similar for PSI" );

  
  //
  // Analysis parameters
  //
  
  double timewin_ms = param.requires_dbl( "w" );
  double order_ms   = param.requires_dbl( "order" );
  
  logger << "  given sample rate of " << sr << "Hz:\n";
  logger << "  window = " << timewin_ms << " ( " << round( ( timewin_ms / (double)1000.0 ) * (double)sr ) << " sample points)\n";
  logger << "  model order = " << order_ms << " ( " << round( ( order_ms / (double)1000.0 ) * (double)sr ) << " sample points)\n";

  //
  // BIC for optimal model order
  //
  
  int compute_bic = param.has( "bic" ) ? param.requires_int( "bic" ) : 0 ;
  

  //
  // Frequency
  //
  
  bool has_frqs = param.has( "f-log" ) || param.has( "f" );
  bool take_logs = param.has( "f-log" );

  std::vector<double> frqs;
  if ( has_frqs ) 
    {
      frqs = take_logs ? param.dblvector( "f-log" ) : param.dblvector( "f" )  ;
      if ( frqs.size() != 3 ) Helper::halt( "expecting f=lwr,upr,n or f-log=lwr,upr,n" );
      
      frqs = take_logs ? 
	MiscMath::logspace( frqs[0] , frqs[1] , frqs[2] ) :
	MiscMath::linspace( frqs[0] , frqs[1] , frqs[2] );
	}


  //
  // Clear tracker
  //

  gc_t::init();


  //
  // Get data, epoch by epoch
  //

  gc_t::ne = edf.timeline.first_epoch();
  
  bool first = true;

  while ( 1 ) 
    {
      
      int epoch = edf.timeline.next_epoch();
      
      if ( epoch == -1 ) break;

      writer.epoch( edf.timeline.display_epoch( epoch ) );

      interval_t interval = edf.timeline.epoch( epoch );

      eigen_matslice_t mslice( edf , signals , interval );

      const Eigen::MatrixXd & X = mslice.data_ref();

      if ( first )
	{
	  const int timewin = round( ( timewin_ms / (double)1000.0 ) * (double)sr  );
	  const int Nr = X.rows() / timewin;
	  logger << "  split each epoch into " << Nr << " non-overlapping windows\n";
	  first = false;
	}

      // fit Granger prediction

      gc_t gc( X , signals , sr , 
	       timewin_ms , order_ms , 	       
	       has_frqs ? &frqs : NULL , 
	       compute_bic );
      
    }

  writer.unepoch();

  
  //
  // Report final averages (non-epoch level)
  //

  gc_t::report( signals );

}



gc_t::gc_t( const Eigen::MatrixXd & X , 
	    const signal_list_t & signals , 
	    int sr , 
	    double timewin_ms , 
	    double order_ms ,
	    const std::vector<double> * frqs ,
	    int compute_bic ,
	    const bool outputs
	    ) 
{
  
  const int timewin = round( ( timewin_ms / (double)1000.0 ) * (double)sr  );
  const int order = round( ( order_ms / (double)1000.0 ) * (double)sr );

  if ( timewin == 0 || order == 0 ) 
    Helper::halt( "bad w, overlap and/or order parameters" );

  // pairs of channels
  
  const int ns = signals.size();  
  
  const int np = X.rows();

  // how many times does timewin fit into np
  
  const int Nr = np / timewin; 
  
  // this is for a single epoch (e.g. 30 seconds)
  // perform sliding window (of timewin samples, each shifted timewin_inc)
  // within each window, extract all signals; fit univariate model for each individually
  // then do all pairs, and calculate GC values
  
  
  //
  // copy data, which will be modified 
  //
  
  Eigen::MatrixXd Z = X;
  
  //
  // detrend and scale (within blocks) 
  //
  
  int pnt = 0;

  for (int block=0; block < Nr; block++)
    {      
      eigen_ops::detrend( Z.block( pnt , 0 , timewin , ns )  );      
      eigen_ops::scale( Z.block( pnt , 0 , timewin , ns ) , true , true );            
      pnt += timewin;
    }


  //
  // check for stationarity (in this epoch)
  //

  
  // ( to do )
  
  
  //
  // for each signal: univariate autoregression
  //
  
  std::vector<double> uniE( ns );

  for (int s=0; s<ns; s++)
    {	  
      armorf_t arm1( Z.col(s) , Nr , timewin , order );
      uniE[s] = arm1.E(0,0) ; 
    }
  
  //
  // Now consider all pairs 
  //
  
  Eigen::MatrixXd ZZ( np , 2 );  
  
  for (int s1=0; s1<ns; s1++)
    {

      if ( outputs ) 
	writer.level( signals.label(s1) , globals::signal1_strat );

      ZZ.col(0) = Z.col(s1);
      
      for (int s2=s1+1; s2<ns; s2++)
	{

	  if ( outputs ) 
	    writer.level( signals.label(s2) , globals::signal2_strat );
	  
	  ZZ.col(1) = Z.col(s2);
	  
	  // fit AR models : two univariate autoregressive models, one bivariate autoregressive model
	  
	  // fit AR models (armorf_t is a straight C/C++ implementation 
	  // if armorf.m from the BSMART toolbox
	  	      	      
	  armorf_t arm12( ZZ  , Nr , timewin , order );
	  
	  // time-domain causal estimate
	  
	  y2x = log( uniE[s1] / arm12.E(0,0) );
	  x2y = log( uniE[s2] / arm12.E(1,1) );
	  
	  //
	  // BIC for optimal model order at each time point?
	  //
	  
	  if ( compute_bic )
	    {
	      
	      int idx = -1;
	      bic = 0; 
	      
	      for (int o=1; o<=compute_bic; o++)
		{		  
		  
		  armorf_t arm12( ZZ  , Nr , timewin , o );
		  
		  // compute BIC 
		  // bic(timei,bici) = log(det(E)) + (log(length(tempdata))*bici*2^2)/length(tempdata);
		  
		  double this_bic =  log( arm12.E.determinant() ) + ( log( np ) * o * 4 ) / (double)np ;
		  
		  if ( o == 1 || this_bic < bic )
		    {
		      idx = o;
		      bic = this_bic;
		    }
		}

	      if ( outputs )     
		writer.value( "BIC" , bic );
	    }
	  
	  
	  //
	  // Frequency-dependent GC
	  //
	  
	  if ( frqs != NULL )
	    {
	      
	      // code below is adapted from bsmart toolbox function pwcausal.m
	      // corrected covariance
	      
	      double eyx = arm12.E(1,1) - arm12.E(0,1) * arm12.E(0,1) / arm12.E( 0,0 );
	      double exy = arm12.E(0,0) - arm12.E(1,0) * arm12.E(1,0) / arm12.E( 1,1 );

	      const int N = arm12.E.rows(); // i.e. 2 

	      Eigen::MatrixXcd Ec = arm12.E;

	      for (int fi=0; fi<frqs->size(); fi++)
		{
		  
		  double f = (*frqs)[fi];
		  	      
		  // transfer matrix
		  Eigen::MatrixXcd H = Eigen::MatrixXd::Identity( N , N );
	      
		  for (int m=1; m<= order; m++)
		    H += (arm12.coeff.block( 0 , (m-1)*N , N , N ).array() 
			  * exp( dcomp( 0 , -m * 2 * M_PI * f / (double)sr ) )).matrix();
		  
       		  Eigen::MatrixXcd Hi = H.inverse();
		  		  
		  Eigen::MatrixXcd S = H.lu().solve( Ec ) * Hi.adjoint() ;
		  S /= (double)sr;
		  
		  tf_x2y[f] = log( abs( S(1,1) ) / abs( S(1,1) - ( Hi(1,0) * exy * std::conj( Hi(1,0) ) ) / (double)sr ) );
		  tf_y2x[f] = log( abs( S(0,0) ) / abs( S(0,0) - ( Hi(0,1) * eyx * std::conj( Hi(0,1) ) ) / (double)sr ) );
		  
		  // next frequency 
		}
	      
	    }
	  
	  //
	  // Outputs
	  //
	  
	  // time-domain

	  if ( outputs )
	    {
	      writer.value( "Y2X" , y2x );  
	      writer.value( "X2Y" , x2y );
	    }
	  
	  y2x_sum[ s1 ][ s2 ] += y2x;  
	  x2y_sum[ s1 ][ s2 ] += x2y;  

	  // frequency-stratified
	  
	  if ( frqs != NULL )
	    {
  
	      std::map<double,double>::const_iterator ff = tf_y2x.begin();
	      while ( ff != tf_y2x.end() )
		{
		  if ( outputs ) {
		    writer.level( ff->first , globals::freq_strat );
		    writer.value( "Y2X" , ff->second );
		  }
		  tf_y2x_sum[ s1 ][ s2 ][ ff->first ] += ff->second;
		  ++ff;
		}
	      
	      ff = tf_x2y.begin();
	      while ( ff != tf_x2y.end() )
		{      
		  if ( outputs ) {
		    writer.level( ff->first , globals::freq_strat );
		    writer.value( "X2Y" , ff->second );
		  }
		  tf_x2y_sum[ s1 ][ s2 ][ ff->first ] += ff->second;
		  ++ff;
		}
	      if ( outputs ) 
		writer.unlevel( globals::freq_strat );

	    }	  

	} // next s2
      if ( outputs ) 
	writer.unlevel( globals::signal2_strat );  	      
    } // next s1 
  if ( outputs ) 
    writer.unlevel( globals::signal1_strat );  	      
}


void gc_t::report( const signal_list_t & signals )
{
  
  // s1
  std::map<int,std::map<int,double> >::const_iterator ss1 = y2x_sum.begin();
  while ( ss1 != y2x_sum.end() )
    {
      const int s1 = ss1->first;
      
      writer.level( signals.label( s1 ) , globals::signal1_strat );

      const std::map<int,double> & sum2 = ss1->second;

      std::map<int,double>::const_iterator ss2 = sum2.begin();
      while( ss2 != sum2.end() )
	{
	  const int s2 = ss2->first;

	  writer.level( signals.label( s2 ) , globals::signal2_strat );

	  // time-domain outputs 
	  writer.value( "Y2X" , y2x_sum[ s1 ][ s2 ] / (double)ne ); 
	  writer.value( "X2Y" , x2y_sum[ s1 ][ s2 ] / (double)ne ); 
	  
	  // frequency domain outputs
	  const std::map<double,double> & fres = tf_x2y_sum[ s1 ][ s2 ];
	  
	  std::map<double,double>::const_iterator ff = fres.begin();
	  while ( ff != fres.end() )
	    {
	      writer.level( ff->first , globals::freq_strat );
	      writer.value( "X2Y" , ff->second / (double)ne ); 
	      ++ff;
	    }
	  
          const std::map<double,double> & fres2 = tf_y2x_sum[ s1 ][ s2 ];
	  ff = fres2.begin();
	  while ( ff != fres2.end() )
	    {
	      writer.level( ff->first , globals::freq_strat );
	      writer.value( "Y2X" , ff->second / (double)ne ); 
	      ++ff;
	    }
	  writer.unlevel( globals::freq_strat );
	 
	  // next channel pair
	  ++ss2;
	}
      writer.unlevel( globals::signal2_strat );

      ++ss1;
    }
  writer.unlevel( globals::signal1_strat );
}



