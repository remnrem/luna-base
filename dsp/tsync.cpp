
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


#include "dsp/tsync.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "helper/helper.h"

#include "db/db.h"
#include "helper/logger.h"

extern writer_t writer;
extern logger_t logger;

void dsptools::tsync( edf_t & edf , param_t & param )
{

  // get signals: assume that X_hilbert_phase and X_hilbert_mag exist

  std::string signal_label = param.requires( "sig" );
  if ( signal_label == "*" ) Helper::halt( "must specify an explicit sig list (roots)" );
  std::vector<std::string> sig_roots = param.strvector( "sig" );
    
  const bool no_annotations = true;

  // window (in seconds )

  const double w_sec = param.has( "w" ) ? param.requires_dbl( "w" ) : 0.25;

  const bool verbose = param.has( "verbose" );

  const bool epoch_output = param.has( "epoch" ); 
  
  // xcorr versus assume signals are from HT?

  bool ht_analysis =  param.has( "ht" ) ;
  
  bool xcorr_analysis = ! ht_analysis;


  if ( xcorr_analysis )
    {
      
      signal_list_t signals = edf.header.signal_list( signal_label , no_annotations );
      
      if ( signals.size() == 0 )
	Helper::halt( "problem locating signals" );
      
      std::vector<double> FsV = edf.header.sampling_freq( signals );
      
      const int Fs = FsV[0];
      for (int i=0; i<signals.size(); i++)
	if ( FsV[i] != Fs )
	  Helper::halt( "unequal sampling frequencies" );
      
      const int ns = signals.size();
      
      int cnt_epoch = 0 ;
      
      //
      // iterate over epochs
      //
      
      int ne = edf.timeline.first_epoch();
      
      std::map<int,std::map<int,std::map<int,double> > > xcorr;
      std::map<int,std::map<int,int> > delay;
      
      while ( 1 )
	{
	  
	  int epoch = edf.timeline.next_epoch();
	  
	  if ( epoch == -1 ) break;
	  
	  interval_t interval = edf.timeline.epoch( epoch );
	  
	  ++cnt_epoch;

	  if ( epoch_output )
	    writer.epoch( edf.timeline.display_epoch( epoch ) );		
	  

	  // get data
	  
	  eigen_matslice_t mslice( edf , signals , interval );
	  
	  const Eigen::MatrixXd & X = mslice.data_ref();
	  
	  const int rows = X.rows();
	  const int cols = X.cols();
	  
	  tsync_t tsync( X , Fs , (int)( w_sec * Fs ) );

	  
	  // epoch-level tracking/reporting
	  
	  for (int s1=0; s1<ns; s1++)
	    {
	      
	      if ( epoch_output )
		writer.level( signals.label(s1) , globals::signal1_strat );
	      
	      for (int s2=s1+1; s2<ns; s2++)
		{

		  if ( epoch_output )
		    writer.level( signals.label(s2) , globals::signal2_strat );

		  // track delay 
		  delay[s1][s2] += tsync.delay[s1][s2];

		  // oputput?
		  if ( epoch_output)
		    writer.value( "D" , tsync.delay[s1][s2] / (double)Fs );

		  // track Xcorr
		  std::map<int,double>::const_iterator ff = tsync.xcorr[s1][s2].begin();
		  while( ff != tsync.xcorr[s1][s2].end() )
		    {
		      xcorr[s1][s2][ ff->first ] += ff->second;
		      
		      if ( verbose )
			{
			  // level ff->first --> ff->second
			}

		      ++ff;
		    }
		}
	      if ( epoch_output ) writer.unlevel( globals::signal2_strat );	  
	    }
	  if ( epoch_output ) writer.unlevel( globals::signal1_strat );
	  
	  // next epoch 
	} 
      if ( epoch_output ) writer.unepoch();

      
      //
      // Report
      //
      

      for (int s1=0; s1<ns; s1++)
	{

	  writer.level( signals.label(s1) , globals::signal1_strat );

	  for (int s2=s1+1; s2<ns; s2++)
	    {

	      writer.level( signals.label(s2) , globals::signal2_strat );

	      writer.value( "S" , ( delay[ s1 ][ s2 ] / (double)cnt_epoch ) / (double)Fs );
	      
	      if ( verbose )
		{
		  std::map<int,double>::const_iterator ff = xcorr[s1][s2].begin();
		  
		  while( ff != xcorr[s1][s2].end() )
		    {
		      // delay
		      writer.level( ff->first , "D" );
		      writer.value( "XCORR" , ff->second / (double)cnt_epoch );		      
		      ++ff;
		    }
		  writer.unlevel( "D" );
		}	      
	    }
	  writer.unlevel( globals::signal2_strat );
	}
      writer.unlevel( globals::signal1_strat );

    }


  //
  // HT analysis 
  //

  
  
  if ( ht_analysis )
    {
      // phase
      signal_label = "";
      for (int s=0; s<sig_roots.size(); s++) signal_label += ( s != 0 ? "," : "" ) + sig_roots[s] + "_ht_ph" ; 
      signal_list_t signals_phase = edf.header.signal_list( signal_label , no_annotations );
      
      // magnitude
      signal_label = "";
      for (int s=0; s<sig_roots.size(); s++) signal_label += ( s != 0 ? "," : "" ) + sig_roots[s] + "_ht_mag" ; 
      signal_list_t signals_mag = edf.header.signal_list( signal_label , no_annotations );
      
      if ( signals_phase.size() != signals_mag.size() || signals_phase.size() == 0 )
	Helper::halt( "problem locating signals with associated _ht_ph and _ht_mag components" );
      
      std::vector<double> FsP = edf.header.sampling_freq( signals_phase );
      std::vector<double> FsM = edf.header.sampling_freq( signals_mag );
      
      const int Fs = FsP[0];
      for (int i=0; i<signals_phase.size(); i++)
	if ( FsP[i] != Fs || FsM[i] != Fs ) Helper::halt( "unequal sampling frequencies" );
      
      const int ns = signals_phase.size();
      
      int cnt_epoch = 0 ;
      
      //
      // iterate over epochs
      //
      
      int ne = edf.timeline.first_epoch();
      
      std::map<int,std::map<int,std::map<int,double> > > ph_diff , ph_lock;
      
      while ( 1 )
	{
	  
	  int epoch = edf.timeline.next_epoch();
	  
	  if ( epoch == -1 ) break;
	  
	  interval_t interval = edf.timeline.epoch( epoch );
	  
	  ++cnt_epoch;
	  
	  eigen_matslice_t mslice_phase( edf , signals_phase , interval );
	  eigen_matslice_t mslice_mag( edf , signals_mag , interval );
	  
	  const Eigen::MatrixXd & P = mslice_phase.data_ref();
	  const Eigen::MatrixXd & M = mslice_mag.data_ref();
	  
	  const int rows = P.rows();
	  const int cols = M.cols();
	  
	  if ( M.rows() != rows || M.cols() != cols )
	    Helper::halt( "problem locating signals with associated _ht_ph and _ht_mag components" ); 
	  
	  tsync_t tsync( P , M , Fs , (int)(0.25 * Fs ) );
	  
	  for (int s1=0; s1<ns; s1++)
	    for (int s2=0; s2<ns; s2++)
	      {
		std::map<int,double>::const_iterator ff = tsync.phase_diff[s1][s2].begin();
		while( ff != tsync.phase_diff[s1][s2].end() )
		  {
		    ph_diff[s1][s2][ ff->first ] += ff->second;
		    ph_lock[s1][s2][ ff->first ] += tsync.phase_lock[s1][s2][ ff->first ];	  

		    // std::cout << "E\t" << epoch << "\t" << s1 << "\t" << s2 << "\t" << ff->first << "\t" 
		    // 	      << ff->second << "\t" << tsync.phase_lock[s1][s2][ ff->first ] << "\n";
		    
		    ++ff;
		  }
	      }
	  
	}
         
  
      //
      // Report
      //
      
      for (int s1=0; s1<ns; s1++)
	{

	  writer.level( signals_phase.label(s1) , globals::signal1_strat );

	  for (int s2=0; s2<ns; s2++)
	    {

	      writer.level( signals_phase.label(s2) , globals::signal2_strat );

	      //	      writer.value( delay[ s1 ][ s2 ] / (double)cnt_epoch );

	      if ( verbose )
		{
		  std::map<int,double>::const_iterator ff = ph_lock[s1][s2].begin();
		  while( ff != ph_lock[s1][s2].end() )
		    {
		      // delay
		      writer.level( ff->first , "D" );
		      writer.value( "PHL" , ff->second / (double)cnt_epoch );		      
		      ++ff;
		    }
		  writer.unlevel( "D" );
		}	      
	    }
	  writer.unlevel( globals::signal2_strat );
	}
      writer.unlevel( globals::signal1_strat );
    }
   
}  



tsync_t::tsync_t( const Eigen::MatrixXd & P ,
		  const Eigen::MatrixXd & A ,
		  double sr ,
		  int w )
{
  
  
  // number of signals
  
  const int ns = P.cols();
  const int np = P.rows();
  
  // consider each pair of signals
  
  for (int s1=0; s1<ns; s1++)
    for (int s2=s1+1; s2<ns; s2++)
      {

	const Eigen::VectorXd & P1 = P.col(s1);
	const Eigen::VectorXd & P2 = P.col(s2);
	
	const Eigen::VectorXd & A1 = A.col(s1);
	const Eigen::VectorXd & A2 = A.col(s2);
	
	// pre-calculate exponents

	std::vector<dcomp> E1( np ) , E2( np );
	
	for (int p = 0 ; p < np ; p++)
	  {
	    E1[p] = exp( dcomp ( 0 , P1[p] ) );
	    E2[p] = exp( dcomp ( 0 , P2[p] ) );
	  }
	
	// look from p = w to p < np-w

	for (int idx = -w ; idx <= w ; idx++ )
	  {	    
	    double sumst0 = 0;
	    double sumst1 = 0 , sumst2 = 0;
	    for (int p = w ; p < np - w ; p++)
	      {
		// angular difference
		sumst0 +=  pdiff( P1[p+idx] , P2[p] ) ;

		// phase locking
		double t = abs( ( E1[p+idx] + E2[p] ) / 2.0 );
		sumst1 += t;
		
		// phase-locking weighted by signal amplitude  
		sumst2 += ( A1[p+idx] + A2[p] ) * t ;
		
	      }

	    sumst0 /= (double)np;
	    sumst1 /= (double)np;
	    sumst2 /= (double)np;
	    
	    phase_diff[ s1 ][ s2 ][ idx ] = sumst0;
	    phase_lock[ s1 ][ s2 ][ idx ] = sumst1;
	    //std::cout << "s1,s2,idx,st = " << s1 << " " << s2 << " " << idx << " " << sumst0 << " " << sumst1 << " " << sumst2  << "\n";
	  }
		
      }
  
}




tsync_t::tsync_t( const Eigen::MatrixXd & X ,
		  double sr ,
		  int w )

{
  xcorr.clear();
  delay.clear();

    // number of signals
  
  const int ns = X.cols();
  const int np = X.rows();
  
  // consider each pair of signals
  
  for (int s1=0; s1<ns; s1++)
    for (int s2=s1+1; s2<ns; s2++)
      {
	
	// for first signal, extract elements minus flanking 'w' ones
	const int npp = np - 2 * w; 
	const Eigen::VectorXd & xX1 = X.col(s1).segment( w , npp ) ;


	// second signal will be shifted from -w to +w samples
	// (extract done below), so need the full signal here

	const Eigen::VectorXd & X2 = X.col(s2);
	
	// look from p = w to p < np-w
	
	double xr_max = 0;

	int xr_max_idx = -w;
	
	for (int idx = -w ; idx <= w ; idx++ )
	  {	    
	    
	    double xr = xX1.dot( X2.segment( w + idx , npp ) );
	    
	    if ( xr > xr_max )
	      {
		xr_max = xr;
		xr_max_idx = idx;
	      }

	    // scale by N
	    xcorr[ s1 ][ s2 ][ idx ] = xr / (double)np;
	    
	  }
	
	// estimated delay in samples
	delay[ s1 ][ s2 ] = xr_max_idx;

	
      }
  
}

