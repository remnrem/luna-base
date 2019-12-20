
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

#include "ica-wrapper.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "eval.h"
#include "db/db.h"

extern writer_t writer;

extern logger_t logger;

void dsptools::ica_wrapper( edf_t & edf , param_t & param )
{

  std::string signal_label = param.requires( "sig" );

  std::string component_tag = param.has( "tag" ) ? param.value( "tag" ) : "ICA";

  bool write_S_matrix = param.has( "file" );
  
  std::string S_matrix_fileroot = write_S_matrix ? param.value( "file" ) : "xxx";

  bool original_signals = param.has( "original-signals" );

  bool do_not_add_channels = param.has( "no-new-channels" );

  //
  // Fetch all data signals (specified via sig=,  which defaults to all signals if not given)
  //

  const bool no_annotations = true;
  
  signal_list_t signals = edf.header.signal_list( signal_label , no_annotations );  

  const int ns = signals.size();


  //
  // Check sample rates
  //
  
  if ( ns < 2 ) return;
  
  const int sr = edf.header.sampling_freq( signals(0) );
  
  for (int i=1;i<ns;i++)
    {      
      if ( edf.header.sampling_freq( signals(i) ) != sr ) 
	Helper::halt( "all signals must have similar SR for ICA" );
    }


  //
  // Fetch sample matrix
  //
  
  matslice_t mslice( edf , signals , edf.timeline.wholetrace() );
  
  // nb. fastICA() this will modify X/mslice...
  Data::Matrix<double> & X = mslice.nonconst_data_ref();
    

  //
  // Number of components
  //

  int compc = param.has( "compc" ) ? param.requires_int( "compc" ) : ns ;


  //
  // ICA
  //
  
  ica_t ica( X , compc );


  //
  // Add new signals
  //

  if ( ! do_not_add_channels ) 
    {
      logger << "  adding " << compc << " new signals to EDF\n";
    
      for (int c=0;c<compc;c++)
	edf.add_signal( component_tag + "_" + Helper::int2str( c+1 ) , sr , *ica.S.col(c).data_pointer() );
      
    }
 
  
  
  //
  // Output mixing matrix, etc
  //

  const int rows = X.dim1();

  const int cols = X.dim2();


  // K : cols x compc
  // A : compc x compc
  // W : compc x compc
  
  for (int i=0;i<compc;i++)
    {
      writer.level( i+1 , "IC1" );
      for (int j=0;j<compc;j++)
	{
	  writer.level( j+1 , "IC2" );	  
	  writer.value( "A" , ica.A[i][j] );
	  writer.value( "W" , ica.W[i][j] );
	}
      writer.unlevel( "IC2" );
    }
  writer.unlevel( "IC1" );


  for (int i=0;i<cols;i++)
    {
      writer.level( signals.label(i) , globals::signal_strat );
      for (int j=0;j<compc;j++)
	{
	  writer.level( j+1 , "IC" );	  
	  writer.value( "K" , ica.K[i][j] );
	}
      writer.unlevel( "IC" );
    }
  writer.unlevel( globals::signal_strat );


  //
  // File-based output
  //
  
  if ( write_S_matrix )
    {
      
      std::string froot = S_matrix_fileroot + "_" ;
      
      std::ofstream S( (froot + "S.txt").c_str() , std::ios::out );
      for (int j=0;j<compc;j++) S << ( j ? "\t" : "" ) << "S" << j+1;
      S << "\n";  
      for (int i=0;i<rows;i++)
	{
	  for (int j=0;j<compc;j++) S << ( j ? "\t" : "" ) << ica.S[i][j] ;      
	  S << "\n";
	}
      S.close();

      if ( original_signals ) 
	{
	  std::ofstream F( (froot + "X.txt").c_str() , std::ios::out );
	  for (int j=0;j<cols;j++) F << ( j ? "\t" : "" ) << "X" << j+1;
	  F << "\n";
	  for (int i=0;i<rows;i++)
	    {
	      for (int j=0;j<cols;j++) F << ( j ? "\t" : "" ) << X[i][j];
	      F << "\n";
	    }
	  F.close();
	}

      //
      // other matrices
      //

      // K : cols x compc
      // A : compc x compc
      // W : compc x compc
      // S : as original data      
      
      std::ofstream K( (froot + "K.txt").c_str() , std::ios::out );
      for (int i=0;i<cols;i++)
	{
	  for (int j=0;j<compc;j++) K << ( j ? "\t" : "" ) << ica.K[i][j];
	  K << "\n";
	}
      K.close();
      
      std::ofstream W( (froot + "W.txt").c_str() , std::ios::out );
      for (int i=0;i<compc;i++)
	{
	  for (int j=0;j<compc;j++) W << ( j ? "\t" : "" ) << ica.W[i][j];
	  W << "\n";
	}
      W.close();
            
      std::ofstream A( (froot + "A.txt").c_str() , std::ios::out );
      for (int i=0;i<compc;i++)
	{
	  for (int j=0;j<compc;j++) A << ( j ? "\t" : "" ) << ica.A[i][j];
	  A << "\n";
	}
      A.close();

    }


}
