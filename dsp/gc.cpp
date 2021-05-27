
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

#include "edf/slice.h"
#include "edf/edf.h"

#include <iostream>
#include "stats/Eigen/Cholesky"
#include "stats/Eigen/SVD"
#include <vector>

#include "stats/eigen_ops.h"

#include "db/db.h"
#include "helper/logger.h"

extern logger_t logger;
extern writer_t writer;

void gc_wrapper( edf_t & edf , param_t & param )
{

  //
  // Get signals
  //

  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) );

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
  
  bool has_frqs = param.has( "f" );

  std::vector<double> frqs;
  if ( has_frqs ) 
    {
      frqs = param.dblvector( "f" );
      if ( frqs.size() != 3 ) Helper::halt( "expecting f=lwr,upr,n" );
      frqs = MiscMath::logspace( frqs[0] , frqs[1] , frqs[2] );
    }


  //
  // Get data, epoch by epoch
  //

  const int ne = edf.timeline.first_epoch();
  
  bool first = true;

  while ( 1 ) 
    {

      int epoch = edf.timeline.next_epoch();
      
      if ( epoch == -1 ) break;

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

      // outputs

      writer.epoch( edf.timeline.display_epoch( epoch ) );
      
      gc.report();

    }
  writer.unepoch();

}


gc_t::gc_t( const Eigen::MatrixXd & X , 
	    const signal_list_t & signals , 
	    int sr , 
	    double timewin_ms , 
	    double order_ms ,
	    const std::vector<double> * frqs ,
	    int compute_bic
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

  //  std::cout << "Z1\n" << Z << "\n\n";

  for (int block=0; block < Nr; block++)
    {      
      eigen_ops::detrend( Z.block( pnt , 0 , timewin , ns )  );      
      eigen_ops::scale( Z.block( pnt , 0 , timewin , ns ) , true );            
      pnt += timewin;
    }

  //  std::cout << "Z2\n" << Z << "\n\n";

  //
  // check for stationarity (in this epoch)
  //

  
  // ( to do )
  
  
  //
  // for each signal: univariate autoregression
  //
  
  std::vector<double> uniE( ns );
  //  std::cout << "Nr " << Nr << " " << timewin << " " << Nr * timewin << " " << np << " " << Z.rows() << " " << order << "\n";

  for (int s=0; s<ns; s++)
    {	  
      
      armorf_t arm1( Z.col(s) , Nr , timewin , order );
      uniE[s] = arm1.E(0,0) ; 
      //      std::cout << "EE" <<  arm1.E(0,0) << "\n";
    }
  
  //
  // Now consider all pairs 
  //
  
  Eigen::MatrixXd ZZ( np , 2 );  
  
  for (int s1=0; s1<ns; s1++)
    {
      
      writer.level( signals.label(s1) , globals::signal1_strat );

      ZZ.col(0) = Z.col(s1);
      
      for (int s2=s1+1; s2<ns; s2++)
	{

	  writer.level( signals.label(s2) , globals::signal2_strat );
	  
	  ZZ.col(1) = Z.col(s2);
	  
	  // fit AR models : two univariate autoregressive models, one bivariate autoregressive model
	  
	  // % fit AR models (model estimation from bsmart toolbox)
	  // 	  [Ax,Ex] = armorf(tempdata(1,:),EEG.trials,timewin_points,order_points);
	  // [Ay,Ey] = armorf(tempdata(2,:),EEG.trials,timewin_points,order_points);
	  // [Axy,E] = armorf(tempdata     ,EEG.trials,timewin_points,order_points);
	      	      
	  armorf_t arm12( ZZ  , Nr , timewin , order );
	  
	  // % time-domain causal estimate
	  // 	  y2x(timei)=log(Ex/E(1,1));
	  // x2y(timei)=log(Ey/E(2,2));
	  
	  y2x = log( uniE[s1] / arm12.E(0,0) );
	  x2y = log( uniE[s2] / arm12.E(1,1) );
	  
	  //
	  // BIC for optimal model order at each time point?
	  //
	  
	  if ( compute_bic )
	    {
	      
	      int idx = -1;
	      double min = 0; 

	      for (int o=1; o<=compute_bic; o++)
		{		  
		  // [Axy,E] = armorf(tempdata,EEG.trials,timewin_points,bici);
		  armorf_t arm12( ZZ  , Nr , timewin , o );
		  
		  // compute BIC 
		  // bic(timei,bici) = log(det(E)) + (log(length(tempdata))*bici*2^2)/length(tempdata);

		  double bic =  log( arm12.E.determinant() ) + ( log( np ) * o * 4 ) / (double)np ;
		  
		  std::cout << "bic = " << bic << "\n";

		  if ( o == 1 || bic < min )
		    {
		      idx = o;
		      min = bic;
		    }
		 
		  
		}
	      
	      // report best model order
	      writer.value( "BIC" , idx );
	      
	    }
	  
	  //
	  // Frequency-dependent GC
	  //
	  
	  if ( frqs != NULL )
	    {
	      
	      // % code below is adapted from bsmart toolbox function pwcausal.m
	      // % corrected covariance
	      //    eyx = E(2,2) - E(1,2)^2/E(1,1);
	      //    exy = E(1,1) - E(2,1)^2/E(2,2);
	      //    N = size(E,1);
	      
	      double eyx = arm12.E(1,1) - arm12.E(0,1) * arm12.E(0,1) / arm12.E( 0,0 );
	      double exy = arm12.E(0,0) - arm12.E(1,0) * arm12.E(1,0) / arm12.E( 1,1 );

	      const int N = arm12.E.rows(); // i.e. 2 

	      Eigen::MatrixXcd Ec = arm12.E;

	      // for fi=1:length(frequencies)

	      for (int fi=0; fi<frqs->size(); fi++)
		{
		  
		  double f = (*frqs)[fi];
		  	      
		  // transfer matrix
		  // H = eye(N);

		  Eigen::MatrixXcd H = Eigen::MatrixXd::Identity( N , N );

		  // for m = 1:order_points
		  //     H = H + Axy(:,(m-1)*N+1:m*N)*exp(-1i*m*2*pi*frequencies(fi)/EEG.srate);
		  // end
	      
		  for (int m=1; m<= order; m++)
		    H += (arm12.coeff.block( 0 , (m-1)*N , N , N ).array() * exp( dcomp( 0 , -m * 2 * M_PI * f / (double)sr ) )).matrix();
		  
		  // 	  Hi = inv(H);
       		  Eigen::MatrixXcd Hi = H.inverse();
		  		  
		  //       S  = H\E*Hi'/EEG.srate;
		  Eigen::MatrixXcd S = H.lu().solve( Ec ) * Hi.adjoint() ;
		  S /= (double)sr;
		  
		  // granger prediction per frequency
		  //  tf_granger(1,fi,timei) = log( abs(S(2,2))/abs(S(2,2)-(Hi(2,1)*exy*conj(Hi(2,1)))/EEG.srate) );
		  //  tf_granger(2,fi,timei) = log( abs(S(1,1))/abs(S(1,1)-(Hi(1,2)*eyx*conj(Hi(1,2)))/EEG.srate) );
		 
		  tf_x2y[f] = log( abs( S(1,1) ) / abs( S(1,1) - ( Hi(1,0) * exy * std::conj( Hi(1,0) ) ) / (double)sr ) );
		  tf_y2x[f] = log( abs( S(0,0) ) / abs( S(0,0) - ( Hi(0,1) * eyx * std::conj( Hi(0,1) ) ) / (double)sr ) );
		  
		  // next frequency 
		}
	      
	    }
	  
	} // next s2
      writer.unlevel( globals::signal2_strat );  	      
    } // next s1 
  writer.unlevel( globals::signal1_strat );  	      
}




void gc_t::report()
{
  
  writer.value( "Y2X" , y2x );

  writer.value( "X2Y" , x2y );

  std::map<double,double>::const_iterator ff = tf_y2x.begin();
  while ( ff != tf_y2x.end() )
    {      
      writer.level( ff->first , globals::freq_strat );
      writer.value( "Y2X" , ff->second );
      ++ff;
    }

  ff = tf_x2y.begin();
  while ( ff != tf_x2y.end() )
    {      
      writer.level( ff->first , globals::freq_strat );
      writer.value( "X2Y" , ff->second );
      ++ff;
    }
  writer.unlevel( globals::freq_strat );

}
