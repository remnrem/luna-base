
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

#include "dsp/peri.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "db/db.h"
#include "helper/logger.h"
#include "helper/helper.h"
#include "cwt/cwt.h"
#include "dsp/wrappers.h"
#include "dsp/xcorr.h"
#include "dsp/gc.h"

extern writer_t writer;
extern logger_t logger;


peri_param_t::peri_param_t( param_t & param )
{

  // time centering
  time0 = param.has( "c" ) ? param.requires_dbl( "c" ) : 0 ;

  // left/right clipping
  time_left = time_right = 0 ; 

  if ( param.has( "w" ) )
    time_left = time_right = param.requires_dbl( "w" );
  else if ( param.has( "l" ) )
    time_left = param.requires_dbl( "l" );
  else if ( param.has( "r" ) )
    time_right = param.requires_dbl( "r" );
  
  // CWT
  cwt_do = param.has( "cwt" );
  cwt_f.clear();
  
  if ( cwt_do )
    {
      // F
      std::vector<double> f = param.dblvector( "cwt" );
      if ( f.size() != 3 ) Helper::halt( "expecting cwt=start,stop,inc" );
      for (double ff = f[0] ; ff <= f[1] ; ff += f[2] ) cwt_f.push_back(ff);

      // FWHM - pick default for now
      cwt_fwhm.resize( cwt_f.size() );
      for (int i=0; i<cwt_fwhm.size(); i++)
	cwt_fwhm[i] = CWT::pick_fwhm( cwt_f[i] );
      
      // 20 seconds default
      cwt_timelength = 20;
      
    }

  // XCORR

  xcorr_w_sec = 0;
  xcorr_c_sec = 0;

      
}

void dsptools::peri( edf_t & edf , param_t & param )
{
  
  if ( ! edf.timeline.epoched() ) 
    Helper::halt( "data must be epoched" );
  
  //
  // get signals
  //
  
  std::string signal_label = param.requires( "sig" );
      
  const bool no_annotations = true;

  signal_list_t signals = edf.header.signal_list( signal_label , no_annotations );

  const int ns = signals.size();

  if ( ns == 0 )
    {
      logger << "  *** none of the requested signals found... bailing\n";
      return;
    }

  //
  // sample rates
  //

  std::vector<double> FsV = edf.header.sampling_freq( signals );
  
  const int Fs = FsV[0];
  for (int i=0; i<signals.size(); i++)
    if ( FsV[i] != Fs )
      Helper::halt( "unequal sampling frequencies" );
  
  
  //
  // iterate over epochs
  //
  
  int ne = edf.timeline.first_epoch();
  
  std::vector<Eigen::MatrixXd> X( ne );
  
  int cnt_epoch = 0 ;
  int epoch_size = -1;
  
  while ( 1 )
    {
      
      int epoch = edf.timeline.next_epoch();
      
      if ( epoch == -1 ) break;
      
      interval_t interval = edf.timeline.epoch( epoch );
      
      // get data      
      eigen_matslice_t mslice( edf , signals , interval );

      const Eigen::MatrixXd & S = mslice.data_ref();
      
      if ( epoch_size == -1 )
	epoch_size = S.rows();
      else if ( epoch_size != S.rows() )
	Helper::halt( "all epochs must be a similar duration" );
      
      X[ cnt_epoch ] = S;
      
      ++cnt_epoch;
      
    } // next epoch

  // get peri param

  peri_param_t pp( param );
  
  // perform peri-event analyses
  
  peri_t peri( X , pp, signals , Fs );
  
  // outputs
  
}

peri_t::peri_t( const std::vector<Eigen::MatrixXd> & X ,
		const peri_param_t & pp , 
		const signal_list_t & signals , 
		const int fs )
{
  
  // epoch-by-epoch data
  // epoch(ne) x sample-point(fs*epoch-dur) x channel(ns) 
  //   each epoch must have a similar size - or else be zero-padded 
  
  const int ne = X.size();

  if ( ne == 0 ) return;

  const int np = X[0].rows();
  const int ns = X[0].cols();
  const int nt = ne * np;
  
  //
  // define the time-track
  //
  
  // *assumption* that every epoch will have a similar structure
  //  (this was enforced in the caller above)
    
  std::vector<double> ts( np );
  const double tdur = np * 1.0/(double)fs ; 
  const double t1 = 1.0/(double)fs;
  for (int t=0; t<np; t++) ts[t] = t * t1 ; 

  //
  // 0) epoch-wise means for normalization (over points)
  //

  Eigen::MatrixXd Emean = Eigen::MatrixXd::Zero( ne , ns );
  Eigen::MatrixXd Esd   = Eigen::MatrixXd::Zero( ne , ns );
  
  for (int e=0; e<ne; e++)
    {
      for (int s=0; s<ns; s++)
	{
	  Emean(e,s) = X[e].col(s).mean();
	}
    }
  
  //
  // 1) means, medians, min/max and SDs
  //

  Eigen::VectorXd r_mean   = Eigen::VectorXd::Zero( np );
  Eigen::VectorXd r_median = Eigen::VectorXd::Zero( np );
  Eigen::VectorXd r_min    = Eigen::VectorXd::Zero( np );
  Eigen::VectorXd r_max    = Eigen::VectorXd::Zero( np );
  Eigen::VectorXd r_sd     = Eigen::VectorXd::Zero( np );

  // for each signal
  for (int s=0; s<ns; s++)
    {
      
      writer.level( signals.label(s) , globals::signal_strat );
      
      // for each time-point, collate all epochs
      for (int p=0; p<np; p++)
	{
	  writer.level( ts[p] , "SEC" );

	  // get mean-corrected vals
	  std::vector<double> x( ne );
	  for (int e=0; e<ne; e++)
	    x[e] = X[e](p,s) - Emean(e,s);


	  // means
	  r_mean[p] = MiscMath::mean( x );
	  r_median[p] = MiscMath::median( x );
	  MiscMath::minmax( x , &(r_min[p]) , &(r_max[p]) );
	  r_sd[p] = MiscMath::sdev( x , r_mean[p] );

	  writer.value( "MEAN" , r_mean[p] );
	  writer.value( "MEDIAN" , r_median[p] );
	  writer.value( "MIN" , r_min[p] );
	  writer.value( "MAX" , r_max[p] );
	  writer.value( "SD" , r_sd[p] );
	  
	 	  
	  
	}
      writer.unlevel( "SEC" );

    }
  writer.unlevel( globals::signal_strat );


  //
  // XCORR
  //

  if ( pp.xcorr_do )
    {

      int mxlag = fs * pp.xcorr_w_sec;
      int cent = fs * pp.xcorr_c_sec;

      // for each pair of signals
      for (int s1=0; s1<ns-1; s1++)
	{
	  writer.level( signals.label(s1) , "CH1" );
	  for (int s2=s1+1; s2<ns; s2++)
	    {
	      writer.level( signals.label(s2) , "CH2" );	      

	      double delay = 0 ;
	      
	      // for each epoch
	      for (int e=0; e<ne; e++)
		{
		  
		  std::vector<double> x1( np );
		  std::vector<double> x2( np );
		  for (int p=0; p<np; p++)
		    {
		      x1[p] = X[e](p,s1);
		      x2[p] = X[e](p,s2);
		    }

		  xcorr_t xc( x1, x2 , mxlag , cent );
		  		  
		  delay += xc.lags[xc.mx] / (double)fs ; 
		}

	      // mean D over epochs
	      writer.value( "D" , delay / (double)ne );
	    }
	  writer.unlevel( "CH2" );
	}
      writer.unlevel( "CH1" );      
    }


  //
  // GC
  //
  
  // go across one window at a time, place
  // all epochs in a single dataframe,
  // and then do a call of gc2_t()
  //  i.e. we have ne timeseries, each
  //       of which is a window of the same position
  //       from across different epochs

  if ( 1  ) //  pp.gc_do
    {

      int order = 2;
      
      // for each pair of signals
      for (int s1=0; s1<ns-1; s1++)
	{
	  writer.level( signals.label(s1) , "CH1" );
	  for (int s2=s1+1; s2<ns; s2++)
	    {
	      writer.level( signals.label(s2) , "CH2" );	      

	      // duration of window, seconds
	      double gc_win_sec = 0.5;

	      // duration of window, sample points
	      int gc_win = fs * gc_win_sec;
	      
	      // number of windows (non-overlapping)
	      int nw = np - gc_win + 1 ; // last possible window slot

	      std::cout << " doing " << nw << " windows  per epoch, each of gc_win = " << gc_win << " samples, np = " << np << "\n";
	      
	      // for each window, each to 'nw'
	      for (int w=0; w<nw; w++)
		{

		  double delay = 0;
		  // average over epochs

		  //std::cout << " going " << w << " to " << w+gc_win-1 << "\n";
		  for (int e=0; e<ne; e++)
		   {
		     
		     std::vector<double> xx( gc_win ) , yy( gc_win );
		     for (int p=0;p<gc_win; p++)
		       {
			 xx[p] = X[e](w+p,s1);
			 yy[p] = X[e](w+p,s2);
		       }

		     if ( e == 22 ) { 
		       std::cout << "\n\n epoch " << e << " window " << w << "\n";
		       for (int i=0; i<10; i++)
			 std::cout << w+i << " = " << xx[i] << "\t" << yy[i] << "\n";
		     }
		     
		     xcorr_t xc( xx, yy );
		     delay += xc.lags[xc.mx] / (double)fs ;

		   }

		  delay /= (double)ne;
		  
		  std::cout << "xc " << s1 << " " << s2 << " w = " << w << " " << delay << "\n";

		  // S.col(0).segment( e * gc_win , gc_win ) = X[e].col(s1).segment( w1, gc_win );
		  // S.col(1).segment( e * gc_win , gc_win ) = X[e].col(s2).segment( w1, gc_win );

		  
		  //     //		  Eigen::MatrixXd S = Eigen::MatrixXd::Zero( ne * gc_win , 2 );
		  
		  // // build up across epochs

		  // for (int e=0; e<ne; e++)
		  //   {
		  //     S.col(0).segment( e * gc_win , gc_win ) = X[e].col(s1).segment( w1, gc_win );
		  //     S.col(1).segment( e * gc_win , gc_win ) = X[e].col(s2).segment( w1, gc_win );
		  //   }

		  		  
		  // // run GC		  

		  // armorf_t arm1( S.col(0) , ne , gc_win , order );
		  // const double uniE_s1 = arm1.E(0,0) ;

		  // armorf_t arm2( S.col(1) , ne , gc_win , order );
		  // const double uniE_s2 = arm2.E(0,0) ;


		  // armorf_t arm12( S  , ne , gc_win , order );

		  // double y2x = log( uniE_s1 / arm12.E(0,0) );
		  // double x2y = log( uniE_s2 / arm12.E(1,1) );

		     
		    
		  // some output : writer.value( "D" , delay / (double)ne );

		} // next window
	      
	    } // next signal2
	  writer.unlevel( "CH2" );
	} // next signal1
      writer.unlevel( "CH1" );      
    }

  
  
  //
  // CWT 
  //

  if ( pp.cwt_do )
    {

      const int nf = pp.cwt_f.size();
      
      logger << "  performing CWT for " << ns << " signals and " << pp.cwt_f.size() << " frequencies...\n";

      // for each signal
      for (int s=0; s<ns; s++)
	{
	  
	  writer.level( signals.label(s) , globals::signal_strat );

	  // store CWT coefs
	  std::vector<std::vector<double> > grand_mag(nf);
	  for (int fi=0;fi<nf;fi++)
	    grand_mag[fi].resize(np,0);

	  
	  // for each epoch
	  for (int e=0; e<ne; e++)
	    {
	      
	      // get mean-corrected vals
	      std::vector<double> x( np ); 
	      for (int p=0; p<np; p++)
		x[p] = X[e](p,s);
	      
	      // for each frequency
	      for (int fi=0; fi<nf; fi++)
		{

		  std::vector<double> mag, phase;
		  
		  const bool pp_wrapped = true;
		  
		  dsptools::alt_run_cwt( x , fs ,
					 pp.cwt_f[fi] , pp.cwt_fwhm[fi] ,
					 pp.cwt_timelength , pp_wrapped ,
					 &mag, NULL ); // last = &phase
		  
		  if ( mag.size() != np )
		    Helper::halt( "internal error in CWT alignment" );
		  
		  for (int p=0;p<np;p++)
		    grand_mag[fi][p] += mag[p];
		  
		} // next frequency

	    } // next epoch

	  // outputs

	  for (int fi=0; fi<nf; fi++)
	    {
	      writer.level( pp.cwt_f[fi] , globals::freq_strat );

	      for (int p=0; p<np; p++)
		{
		  writer.level( ts[p] , "SEC" );		  
		  writer.value( "CWT" , grand_mag[fi][p] / (double)ne );
		} // next time point
	      writer.unlevel( "SEC" );
	      
	    } // next freq
	  writer.unlevel( globals::freq_strat );
	  
	} // next signal
      writer.unlevel( globals::signal_strat );
    }


  
  
}

 
