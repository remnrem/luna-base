
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

  std::cout << "h0\n";

  std::string signal_label = param.requires( "sig" );
  
  signal_list_t signals = edf.header.signal_list( signal_label );  


  const int ns = signals.size();

  std::cout << "ns = " << ns << "\n";
  
  // Assuming multiple signals, all with similar sample rates
  if ( ns < 2 ) return;
  
  const int sr = edf.header.sampling_freq( signals(0) );
  logger << "st = " << sr << "\n";
  for (int i=1;i<ns;i++)
    {
      logger << "s2 = " <<  edf.header.sampling_freq( signals(i) ) << "\n";
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
      for (int i=0;i<cols;i++) pX[i][j] = (*data)[i];
    }
  
  //
  // Number of components
  //

  int compc = param.has( "compc" ) ? param.requires_int( "compc" ) : ns ;

  std::cout << "h1\n";
  
  //
  // ICA
  //

  ica_t ica( pX , rows , cols , compc );

  std::cout << "h2\n";
  
  //
  // Output
  //

  for (int i=0;i<rows;i++)
    {
      std::cout << i ;
      for (int j=0;i<cols;j++) std::cout << "\t" << pX[i][j] << "\t";
      for (int j=0;i<compc;j++) std::cout << "\t" << ica.S[i][j] ;      
    }

  // other matrices
  
  // K : cols x compc
  // A : compc x compc
  // W : compc x compc
  // S : as original data

  logger << "K\n";
  for (int i=0;i<cols;i++)
    {
      for (int j=0;j<compc;j++) logger << "\t" << ica.K[i][j];
      logger << "\n\n";
    }

  logger << "W\n";
  for (int i=0;i<compc;i++)
    {
      for (int j=0;j<compc;j++) logger << "\t" << ica.W[i][j];
      logger << "\n\n";
    }
  
  logger << "A\n";
  for (int i=0;i<compc;i++)
    {
      for (int j=0;j<compc;j++) logger << "\t" << ica.A[i][j];
      logger << "\n\n";
    }


}
