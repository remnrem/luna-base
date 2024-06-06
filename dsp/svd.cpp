
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


#include "dsp/svd.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "stats/Eigen/Dense"

#include "helper/helper.h"
#include "helper/logger.h"
#include "db/db.h"

extern logger_t logger;
extern writer_t writer;

#include "stats/eigen_ops.h"

void dsptools::svd_wrapper( edf_t & edf , param_t & param )
{

  //
  // Fetch all data signals from 'sig'
  //
  
  const std::string signal_label = param.requires( "sig" );

  const bool no_annotations = true;
  
  signal_list_t signals = edf.header.signal_list( signal_label , no_annotations );  
  
  const int ns = signals.size();

  
  //
  // Options
  //
    
  const std::string component_tag = param.has( "tag" ) ? param.value( "tag" ) : "U_";
  
  const bool do_not_add_channels = param.has( "no-new-channels" );
  
  // Number of components 'nc' parameter
  int nc = param.has( "nc" ) ? param.requires_int( "nc" ) : ns ;

  const bool norm_chs = param.has( "norm" ) ? param.yesno( "norm" ) : false; 

  const bool do_winsor = param.has( "winsor" ) ;
  
  const double winsor_q = do_winsor ? param.requires_dbl( "winsor" ) : -1;

  if ( do_winsor && ( winsor_q > 0.5 || winsor_q < 0 ) )
    Helper::halt( "winsor must be between 0 and 0.5" );

  logger << "  extracting " << nc << " components from " << ns << " channels\n";
  if ( norm_chs ) logger << "  standardizing each channel to unit variance\n";
  if ( do_winsor ) logger << "  winsorizing at " << winsor_q * 100 << "%\n";

  //
  // Check sample rates
  //
  
  if ( ns < 2 ) return;
  
  const int sr = edf.header.sampling_freq( signals(0) );
  
  for (int i=1;i<ns;i++)
    {      
      if ( edf.header.sampling_freq( signals(i) ) != sr ) 
	Helper::halt( "all signals must have similar SR for SVD" );
    }


  //
  // Fetch sample matrix
  //
  
  eigen_matslice_t mslice( edf , signals , edf.timeline.wholetrace() );
  
  Eigen::MatrixXd & X = mslice.nonconst_data_ref();

  const int rows = X.rows();
  const int cols = X.cols();
  
  //
  // prep (mean-center cols, optionally normalize and winsorize) 
  //
  
  // mean-center (within each individual)
  eigen_ops::robust_scale( X, true , // center
			   norm_chs , // normalize
			   winsor_q ); 
  
  //
  // SVD
  //
  
  Eigen::BDCSVD<Eigen::MatrixXd> svd( X , Eigen::ComputeThinU | Eigen::ComputeThinV );
  Eigen::MatrixXd U = svd.matrixU();
  Eigen::MatrixXd V1 = svd.matrixV();
  Eigen::VectorXd W1 = svd.singularValues();

  double wsum = W1.sum();
  double cum = 0;
  // output 
  for (int i = 0 ; i < W1.size(); i++ )
    {
      writer.level( i+1 , "C" );
      writer.value( "W" ,W1[i] );
      writer.value( "VE" ,W1[i] / wsum );
      cum += W1[i];
      writer.value( "CVE" ,cum / wsum );
      writer.value( "INC" , (int)( i < nc ) );
    }
  writer.unlevel( "C" );

  // V
  for (int i=0;i<V1.rows(); i++)
    {
      writer.level( signals.label(i) , "FTR" ); 
      for (int j=0;j<nc; j++)
      {
	writer.level( j+1, "C" );
	writer.value( "V" , V1(i,j) );
      }
      writer.unlevel( "C" );
    }
  writer.unlevel( "FTR" );
  
  //
  // Add new signals
  //
  
  if ( ! do_not_add_channels ) 
    {
      logger << "  adding " << nc << " new signals to EDF:";
      
      for (int c=0;c<nc;c++)
	{
	  // scale by W
	  U.col(c) = U.col(c).array() * W1[c];

	  // add as channel
	  std::vector<double> copy( rows );
	  Eigen::VectorXd::Map( &copy[0], rows ) = U.col(c);
	  logger << " " << component_tag + Helper::int2str( c+1 ) ;
	  edf.add_signal( component_tag + Helper::int2str( c+1 ) , sr , copy );	  
	}
      logger << "\n";
    }
  
}

