
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
  
  // initially, perform on whole signal
  //  (can add epoch option later)
  
  // get signals

  signal_list_t signals = edf.header.signal_list( param.value( "sig" ) );
  edf.header.drop_annots_from_signal_list( &signals );
  const int ns = signals.size();

  // get data 
  
  eigen_matslice_t mslice( edf , signals , edf.timeline.wholetrace() );

  // i.e. we will modify this
  Eigen::MatrixXd & X = mslice.nonconst_data_ref();

  const int rows = X.rows();
  const int cols = X.cols();

  bool winsor = param.has( "winsor" );
  double wt = winsor ? param.requires_dbl( "winsor" ) : 0 ; 
  bool second_norm = ! param.has( "no-second-norm" );

  logger << "  robust standardization of " << ns << " signals";
  if ( winsor > 0 ) logger << ", winsorizing at " << wt;
  logger << "\n";

  //
  // do work
  //
  
  eigen_ops::robust_scale( X , wt , second_norm );

  //
  // place back 
  //

  for (int s=0; s<ns; s++)
    {
      std::vector<double> nsig = eigen_ops::copy_vector( X.col(s) );
      edf.update_signal( signals(s) , &nsig );
    }
 
  
}

