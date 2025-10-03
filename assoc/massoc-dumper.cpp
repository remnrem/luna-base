
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

#ifdef HAS_LGBM

#include "assoc/massoc.h"

#include "param.h"
#include "edf/edf.h"
#include "edf/slice.h"


// by epoch, dump out raw binary files for MASSOC inputs
// constraints - epochs must be the same size
// signals are concatenated in same epoch-row, so do no need to have the same SR

void massoc_t::massoc_dumper( edf_t & edf , param_t & param ) 
{
  
  const std::string signal_label = param.requires( "sig" );
  
  signal_list_t signals = edf.header.signal_list( signal_label );  
  
  const int ns = signals.size();
  
  if ( ns == 0 ) return ;

  std::vector<double> Fs = edf.header.sampling_freq( signals );

  int totsp = 0 ;
  for (int s=0;s<ns; s++)
    totsp +=  Fs[s];
  
  const bool to_stdout = param.has( "dump" ) ;

  const std::string class_id = param.requires( "id" );

  
  //
  // Iterate over epochs 
  //

  int ne = edf.timeline.first_epoch();
  
  // track epoch size (samples)

  int es = 0;
  
  int rcnt = -1;
  
  Eigen::MatrixXd X;
  
  // IDs to track
  
  std::vector<std::string> IID, ID, EID;
  
  while ( 1 ) 
    {
      
      int epoch = edf.timeline.next_epoch();      
       
      if ( epoch == -1 ) break;
      
      interval_t interval = edf.timeline.epoch( epoch );
      
      ++rcnt;
      
      int ccnt = -1;
      
      std::vector<double> ex;
      
      for ( int s=0; s<ns; s++ )
	{
	  
	  if ( edf.header.is_annotation_channel( signals(s) ) ) continue;	  
	  
	  slice_t slice( edf , signals(s) , interval );
	  
	  std::vector<double> * d = slice.nonconst_pdata();
	  
	  for (int p=0;p<d->size();p++)
	    ex.push_back( (*d)[p] );
	  
	} // next signal
      

      
      // check
      
      if ( es == 0 ) 
	{
	  es = ex.size();
	  
	  // rows = epochs/events
	  // cols = time-points X channels
	  X = Eigen::MatrixXd::Zero( ne , es );
	}
      else
	{
	  if ( es != ex.size() )
	    Helper::halt( "epochs must be of similar size" );
	}
      

      // save into X
      for (int p=0; p<es; p++)
	X(rcnt , p ) = ex[p];
      
      // also dump to console?
      if ( to_stdout ) 
	{
	  std::cout << rcnt << "\t" << ccnt;
	  for (int p=0; p<es; p++)
	    std::cout << "\t" << ex[p];
	  std::cout << "\n";
	}
      
      
      //
      // IDs
      //
      
      IID.push_back( edf.id );
      ID.push_back( class_id );
      EID.push_back( Helper::int2str( rcnt ) );
      
    } // next epoch
  
  
  //
  // Dump for MASSOC?
  //

  // write a binary file of rows = feature vectors
  // each has an ID, which can either be unique to the row or to the individual
  
  const std::string filename = param.requires( "file" );
 
  const int nrow = X.rows();
  const int ncol = X.cols();
  
  logger << "  constructed feature matrix, " << nrow << " observations by " << ncol << " features\n";
  
  
#ifdef HAS_LGBM

  // MASSOC IDS
  //   IID    individual ID (for concatentating files) 
  //   ID     event type, e.g. NOCS, SOCS, etc
  //   EID    event ID, e.g. 1,2,3,4,... within class of ID
  
  
  // column IDs: simply 1, 2, 3, etc
  std::vector<std::string> COLID( ncol );
  for (int i=0; i<ncol; i++)
    COLID[i] = Helper::int2str( i+1 );
  
  
  // save
  
  massoc_t massoc( IID, ID, EID, COLID, X, filename );
  
#else
  Helper::halt( "LGBM support not compiled in" );
#endif
  

}


#endif
