
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


#include "ica.h"

#include "edf/edf.h"
#include "helper/helper.h"
#include "main.h"
#include "db/db.h"

extern writer_t writer;

// wrapper around libICA function

bool ica_t::proc( const std::vector<std::vector<double> > & X , int compc )
{
  int rows = X.size();
  if ( rows == 0 ) return false;
  int cols = X[0].size();
  
  mat pX = mat_create( rows , cols );
  
  mat pW = mat_create(compc, compc);
  mat pA = mat_create(compc, compc);
  mat pK = mat_create(cols, compc);
  mat pS = mat_create(rows, cols);     

  // compute ICA
  fastICA(pX, rows, cols, compc, pK, pW, pA, pS);
  
  copy( W, pW, compc, compc );
  copy( A, pA, compc, compc );
  copy( K, pA, cols , compc );
  copy( S, pA, rows , cols );
  
  return true;
}


bool ica_t::proc( mat pX , int rows , int cols , int compc )
{  
  std::cout << "in here.\n";
  double ** pW = mat_create(compc, compc);
  double ** pA = mat_create(compc, compc);
  double ** pK = mat_create(cols, compc);
  double ** pS = mat_create(rows, cols);     

  // mean center
  double * means = vect_create( cols );
  mat_center( pX, rows, cols, means );

  // compute ICA
  fastICA(pX, rows, cols, compc, pK, pW, pA, pS);

  copy( W, pW, compc, compc );
  copy( A, pA, compc, compc );
  copy( K, pK, cols , compc );
  copy( S, pS, rows , cols );

  mat_delete( pW , compc , compc );
  mat_delete( pA , compc , compc );
  mat_delete( pK , cols  , compc );
  mat_delete( pS , rows  , cols );
  
  return true;
}





 void ica_wrapper( edf_t & edf , param_t & param )
 {

   std::string signal_label = param.requires( "signal" );
   
   signal_list_t all_signals = edf.header.signal_list( signal_label );  

   // only keep data-channels
   signal_list_t signals;
   
   for (int s=0;s< all_signals.size();s++)
     if ( edf.header.is_data_channel( all_signals(s) ) ) 
       signals.add( all_signals(s) , all_signals.label(s) );

   const int ns = signals.size();
   
   // Assuming multiple signals, all with similar sample rates
   if ( ns < 2 ) return;
   
   const int sr = edf.header.sampling_freq( signals(0) );
   std::cout << "ns = " << ns << " st = " << sr << "\n";
   
   for (int i=1;i<ns;i++)
     {
       std::cout << "sr2 = " << edf.header.sampling_freq( signals(i) )  << "\t"
		 << edf.header.label[ signals(i) ] << "\n";
       
       if ( edf.header.sampling_freq( signals(i) ) != sr ) 
	 Helper::halt( "all signals must have similar SR for ICA" );
     }
   
   // Fetch sample matrix
   mslice_t mslice( edf , signals , edf.timeline.wholetrace() );

   const std::vector<double> * data = mslice.channel[0]->pdata();

   int rows = data->size();
   int cols = ns;
   
   mat pX = mat_create( rows , cols );
   
   for (int j=0;j<cols;j++)
     {
       const std::vector<double> * data = mslice.channel[j]->pdata();      
       for (int i=0;i<rows;i++) pX[i][j] = (*data)[i];
     }

   //
   // Number of components
   //

   int compc = param.has( "compc" ) ? param.requires_int( "compc" ) : ns ;
   
   //
   // ICA
   //

   ica_t ica( pX , rows , cols , compc );


   //
   // Output
   //

   for (int i=0;i<rows;i++)
     {
       std::cout << i ;
       for (int j=0;j<cols;j++) std::cout << "\t" << pX[i][j] << "\t";
       for (int j=0;j<compc;j++) std::cout << "\t" << ica.S[i][j] ;      
       std::cout << "\n";
     }

   // other matrices

   // K : cols x compc
   // A : compc x compc
   // W : compc x compc
   // S : as original data

   std::cout << "K\n";
   for (int i=0;i<cols;i++)
     {
       for (int j=0;j<compc;j++) std::cout << "\t" << ica.K[i][j];
       std::cout << "\n\n";
     }

   std::cout << "W\n";
   for (int i=0;i<compc;i++)
     {
       for (int j=0;j<compc;j++) std::cout << "\t" << ica.W[i][j];
       std::cout << "\n\n";
     }

   std::cout << "A\n";
   for (int i=0;i<compc;i++)
     {
       for (int j=0;j<compc;j++) std::cout << "\t" << ica.A[i][j];
       std::cout << "\n\n";
     }


 }
