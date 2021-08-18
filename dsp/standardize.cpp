
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
#include "edf/edf.h"
#include "edf/slice.h"
#include "stats/eigen_ops.h"


void dsptools::standardize( edf_t & edf , param_t & param )
{

  const bool by_epoch = param.has( "epoch" );
  
  // center (based on median in first round)
  const bool center = param.has( "center" ) ? param.yesno( "center" ) : true ; 
  
  // scale (based on IQR robust estimate of SD)
  const bool scale = param.has( "scale" ) ? param.yesno( "scale" ) : true ; 
  
  // also winsorize? (default= no)
  const bool winsor = param.has( "winsor" );
  const double wt = winsor ? param.requires_dbl( "winsor" ) : 0 ; 
  
  // after Winsorization, (non-robust) norm a second time (i.e. to ensure mean of 0, SD of 1)
  // default = no
  const bool second_norm = param.has( "second-norm" ) ? param.yesno( "second-norm" ) : false ; 

  // get signals
  signal_list_t signals = edf.header.signal_list( param.value( "sig" ) );
  edf.header.drop_annots_from_signal_list( &signals );
  const int ns = signals.size();


  // do by epoch?
  if ( by_epoch ) edf.timeline.ensure_epoched();

  if ( by_epoch ) logger << "  iterating over epochs\n";
  else logger << "  correcting for entire signal\n";

  logger << "  robust standardization of " << ns << " signals";
  if ( winsor > 0 ) logger << ", winsorizing at " << wt;
  logger << "\n";

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
      
      // get data (yes, dupes)
      eigen_matslice_t mslice( edf , signals , interval );
      Eigen::MatrixXd & T = mslice.nonconst_data_ref();      

      // process
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

