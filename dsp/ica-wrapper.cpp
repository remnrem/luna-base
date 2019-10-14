
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
  
  const bool no_annotations = true;

  signal_list_t signals = edf.header.signal_list( signal_label , no_annotations );  

  const int ns = signals.size();

  std::cout << "ns = " << ns << "\n";
  
  // Assuming multiple signals, all with similar sample rates
  if ( ns < 2 ) return;
  
  const int sr = edf.header.sampling_freq( signals(0) );
  for (int i=1;i<ns;i++)
    {      
      if ( edf.header.sampling_freq( signals(i) ) != sr ) 
	Helper::halt( "all signals must have similar SR for ICA" );
    }
  
  // Fetch sample matrix
  mslice_t mslice( edf , signals , edf.timeline.wholetrace() );
  
  const std::vector<double> * data = mslice.channel[0]->pdata();

  int rows = data->size();
  int cols = ns;
  mat pX = mat_create( rows , cols );
  
  for (int j=0;j<ns;j++)
    {
      const std::vector<double> * data = mslice.channel[j]->pdata();      
      for (int i=0;i<rows;i++) pX[i][j] = (*data)[i];
    }
  
  //
  // Number of components
  //

  int compc = param.has( "compc" ) ? param.requires_int( "compc" ) : ns ;

  logger << "  running with " << compc << " components\n"; 
  
  //
  // ICA
  //

  ica_t ica( pX , rows , cols , compc );
  
  logger << "  finished ICA\n";
  
  //
  // Output
  //

  std::string froot = "ica_";
  
  std::ofstream S( (froot + "S.txt").c_str() , std::ios::out );
  for (int j=0;j<compc;j++) S << ( j ? "\t" : "" ) << "S" << j+1;
  S << "\n";  
  for (int i=0;i<rows;i++)
    {
      for (int j=0;j<compc;j++) S << ( j ? "\t" : "" ) << ica.S[i][j] ;      
      S << "\n";
    }
  S.close();

  std::ofstream X( (froot + "X.txt").c_str() , std::ios::out );
  for (int j=0;j<cols;j++) X << ( j ? "\t" : "" ) << "X" << j+1;
  X << "\n";
  for (int i=0;i<rows;i++)
    {
      for (int j=0;j<cols;j++) X << ( j ? "\t" : "" ) << pX[i][j];
      X << "\n";
    }
  X.close();
  

  // other matrices
  
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
