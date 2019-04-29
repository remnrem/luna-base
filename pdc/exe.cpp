
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

#include "pdc.h"

#include "helper/helper.h"
#include "helper/logger.h"
#include "db/db.h"
#include "eval.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "dsp/resample.h"
#include "stats/cluster.h"

#include <vector>
#include <map>


extern logger_t logger;

extern writer_t writer;

void pdc_t::similarity_matrix( edf_t & edf , param_t & param )
{


  
  // ExE uni mat=output-root  

  //
  // For multiple signals, create a single matrix based on all channels 
  // *unless* uni is specified, in which case write clusters separately
  // for each channel (i.e. output stratified by channel) 
  // -- do not allow mat and uni to be specified together/
  //
  
  bool univariate = param.has( "uni" );
  

  //
  // write cluster solution to normal output stream
  // write to a separate txt file, not an output-db
  //

  
  bool write_matrix = param.has( "mat" );
  
  if ( univariate && write_matrix ) 
    Helper::halt( "cannot specify uni and mat together" );
  
  std::string outfile = "";
  if ( write_matrix ) 
    outfile = param.requires( "mat" ) + "-" + edf.id + ".mat" ;

  
  
  //
  // Signals and sample-rate
  //
  
  std::string signal_label = param.requires( "sig" );   
  
  signal_list_t signals = edf.header.signal_list( signal_label );  

  const int ns = signals.size();

  // desired
  int sr = param.has( "sr" ) ? param.requires_int( "sr" ) : -1;

  // actual
  std::vector<double> Fs = edf.header.sampling_freq( signals );

  // resampling?
  if ( sr != -1 ) 
    {
      for (int s=0;s<ns;s++)
	{

	  if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
	  
	  if ( edf.header.sampling_freq( signals(s) ) != sr ) 
	    dsptools::resample_channel( edf, signals(s) , sr );
	  
	}  
    }

  
  //
  // Requires data to be epoched
  //

  if ( !edf.timeline.epoched() )
    Helper::halt( "ExE requires epoched data" );
  

  //
  // Set entropy values?
  //

  if ( param.has( "entropy" ) ) 
    {
      if ( univariate ) Helper::halt( "cannot specify uni and entropy options together" );

      // automatically set m and t
      entropy_heuristic_wrapper( param );
    }
  else 
    {
      m = param.has( "m" ) ? param.requires_int( "m" ) : 5 ;
      t = param.has( "t" ) ? param.requires_int( "t" ) : 1 ;
    }
  
  //
  // Outer loop (if in univariate mode)
  //
  
  for (int s0=0; s0<ns; s0++)
    {     

      
      //
      // Reset obs()
      //
      
      clear();

      
      //
      // Add channels
      //
      
      if ( univariate )
	{
	  // just add one channel
	  add_channel( signals.label(s0) );
	  // and stratify output by channel
	  writer.level( signals.label(s0) , globals::signal_strat );
	}
      else
	{
	  // add all channels (and we'll quite afterwards)
	  for ( int s=0; s<ns; s++ )
	    add_channel( signals.label(s) );
	}
      
      //
      // Get list of included epochs
      //
     
      int ne = edf.timeline.first_epoch();
 
      std::vector<int> epochs;
      
      while ( 1 ) 
	{
	  int epoch = edf.timeline.next_epoch();      
	  
	  if ( epoch == -1 ) break;
	  
	  interval_t interval = edf.timeline.epoch( epoch );
	  
	  // record for this epoch
	  
	  pdc_obs_t ob( ns );
	  
	  ob.id = epoch;
	  ob.label = ".";
	  
	  epochs.push_back( edf.timeline.display_epoch( epoch ) ) ; 
	  
	  
	  //
	  // get signal(s)	  
	  //
	  
	  if ( univariate )
	    {
	      if ( edf.header.is_annotation_channel( signals(s0) ) ) continue;

	      slice_t slice( edf , signals(s0) , interval );
	      
	      std::vector<double> * d = slice.nonconst_pdata();
	      
	      ob.ch[0] = true;
	      ob.ts[0] = *d;
	    }
	  else
	    {
	      for ( int s=0; s<ns; s++ )
		{
		  
		  // only consider data tracks
		  
		  if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
		  
		  slice_t slice( edf , signals(s) , interval );
		  
		  std::vector<double> * d = slice.nonconst_pdata();
		  
		  ob.ch[s] = true;
		  ob.ts[s] = *d;
		  
		}
	    }

	  //
	  // add data 
	  //
	  
	  add( ob );
          
	  //
	  // next epoch
	  //
	}

      
      //
      // Encode time-series
      //
      
      encode_ts();
      
      //
      // Calculate distance matrix
      //
      
      Data::Matrix<double> D = all_by_all();
      
      
      if ( D.dim1() != ne ) 
	Helper::halt( "internal error in pdc_t::similarity_matrix()" );
      
      if ( write_matrix )
	{
	  std::ofstream OUT1( outfile.c_str() , std::ios::out );
	  
	  for (int i=0;i<ne;i++)
	    {
	      for (int j=0;j<ne;j++) OUT1 << ( j ? "\t" : "" ) << D[i][j];
	      OUT1 << "\n";
	    }
	  OUT1.close();
	  
	  logger << " output distance matrix for " << ne
		 << " epochs (" << ns << " signals) to " << outfile << "\n";
	}
      
      
      //
      // Cluster
      //
      
      cluster_t cluster;
      
      int preK = param.has( "k" ) ? param.requires_int( "k" ) : 0 ;
      int maxS = param.has( "mx" ) ? param.requires_int( "mx" ) : 0 ;

      cluster_solution_t sol = cluster.build( D , preK , maxS );
      
      if ( sol.best.size() != ne ) Helper::halt( "internal error in ExE" );
      
      for (int e=0; e<ne; e++ )
	{
	  writer.epoch( edf.timeline.display_epoch( e ) );
	  writer.value( "CL" , sol.best[e] );
	}
      writer.unepoch();
     
      //
      // in multivariate-mode, we are all done after one pass
      //
      
      if ( ! univariate ) break;
     
      
      //
      // next channel in univariate mode
      //
    }

  if ( univariate ) 
    writer.unlevel( globals::signal_strat );

}
