
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

#include "dsp/standardize.h"

#include "param.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "stats/eigen_ops.h"


void dsptools::standardize( edf_t & edf , param_t & param )
{

  const bool silent = param.yesno( "silent", false, true );
  				   
  const bool by_epoch = param.has( "epoch" );

  // only do simple de-meaning? (perhaps per epoch)
  const bool simple_demean = param.has( "simple-demean" ); 
  
  // ( X - median ) / ( IQR ) 
  const bool iqr_norm = param.has( "IQR" );
  
  // center (based on median in first round)
  const bool center = param.has( "center" ) ? param.yesno( "center" ) : ( ! simple_demean ) ;
  
  // scale (based on IQR robust estimate of SD)
  const bool scale = param.has( "scale" ) ? param.yesno( "scale" ) : ( ! simple_demean )  ; 
  
  // also winsorize? (default= no)
  const bool winsor = param.has( "winsor" );
  const double wt = winsor ? param.requires_dbl( "winsor" ) : 0 ; 
  
  // after Winsorization, (non-robust) norm a second time (i.e. to ensure mean of 0, SD of 1)
  // default = no
  const bool second_norm = param.has( "second-norm" ) ? param.yesno( "second-norm" ) : false ;

  if ( ! ( center || scale || winsor || simple_demean ) )
    {
      if ( ! silent ) 
	logger << "  nothing to do, leavning standardization\n";
      return;
    }

  if ( simple_demean && ( center || scale || second_norm || winsor || iqr_norm ) )
    Helper::halt( "cannot combine simple-demean with center or scale or second-norm or winsor or IQR" );
  
  if ( winsor && ! ( center || scale || simple_demean ) )
    if ( ! silent ) 
      logger << "  only winsorizing signals, not performing initial standardization\n";
   
  // get signals
  signal_list_t signals = edf.header.signal_list( param.value( "sig" ) );
  edf.header.drop_annots_from_signal_list( &signals );
  const int ns = signals.size();


  // do by epoch?
  if ( by_epoch ) edf.timeline.ensure_epoched();

  if ( ! silent ) {
    
    if ( by_epoch ) logger << "  iterating over epochs\n";
    else logger << "  correcting for entire signal\n";
    
    if ( iqr_norm )
      logger << "  IQR-based standardization of " << ns << " signals\n";	    
    else
      {
	logger << "  robust standardization of " << ns << " signals";
	if ( winsor > 0 ) logger << ", winsorizing at " << wt;
	logger << "\n";
      }
  }
  
  //
  // get data (whole signal) 
  //

  eigen_matslice_t mslice( edf , signals , edf.timeline.wholetrace() );
  
  // i.e. we will modify this
  Eigen::MatrixXd & X = mslice.nonconst_data_ref();
  const int rows = X.rows();
  const int cols = X.cols();
  
  
  //
  // iterate over each epoch / or do whole signal in one go?
  //

  int ne = by_epoch ? edf.timeline.first_epoch() : 1 ;
  
  // for sample updates
  int r = 0;

  while ( 1 )
    {
      
      // next epoch
      int epoch = by_epoch ? edf.timeline.next_epoch() : 1 ; 
      
      // all done?
      if ( epoch == -1 ) break;
	  
      // epoch or whole trace?
      interval_t interval = by_epoch ? edf.timeline.epoch( epoch ) : edf.timeline.wholetrace();
      
      // get data (yes, dupes effort...)
      eigen_matslice_t mslice( edf , signals , interval );
      Eigen::MatrixXd & T = mslice.nonconst_data_ref();      

      // process
      if ( simple_demean )	
	eigen_ops::scale( T , true , false ); 
      else if ( iqr_norm )
	eigen_ops::IQR_norm( T );
      else
	eigen_ops::robust_scale( T , center , scale , wt , second_norm );
      
      // update X
      const int trows = T.rows();      
      for (int i=0; i<trows; i++)
	{
	  for (int j=0; j<cols; j++)
	    X(r,j) = T(i,j);
	  // next row of X
	  ++r;
	}
      
      // done?
      if ( ! by_epoch ) break;
      
      // next epoch
    }
      
  // update all signals
  
  for (int s=0; s<ns; s++)
    {
      std::vector<double> nsig = eigen_ops::copy_vector( X.col(s) );
      edf.update_signal( signals(s) , &nsig );
    }
   
}


void dsptools::rolling_standardize( edf_t & edf , param_t & param )
{
  // window size ( seconds )
  const double w = param.requires_dbl( "w" );

  if ( w < 1 )
    Helper::halt( "w must be at least 1 second" );
  
  // get signals
  signal_list_t signals = edf.header.signal_list( param.value( "sig" ) );
  edf.header.drop_annots_from_signal_list( &signals );
  const int ns = signals.size();
      
  // whole trace?
  interval_t interval = edf.timeline.wholetrace();
  
  // get data (yes, dupes effort...)
  eigen_matslice_t mslice( edf , signals , interval );
  Eigen::MatrixXd & T = mslice.nonconst_data_ref();      
  
  // update all signals  
  for (int s=0; s<ns; s++)
    {
      const int sr = edf.header.sampling_freq( signals(s) );
      Eigen::VectorXd Z = eigen_ops::rolling_norm( T.col(s) , sr * w );
      std::vector<double> nsig = eigen_ops::copy_vector( Z );
      edf.update_signal( signals(s) , &nsig );
    }
   
}

