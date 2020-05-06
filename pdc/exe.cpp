
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
  
  
  //
  // ensure we've cleared any prior PDC obs
  //

  clear();
  
  
  // ExE uni cat mat=output-root  

  //
  // Three modes:
  //
  //   for E epochs and K channels
  //
  //       uni       K different channel-specific ExE matrices
  //       [default] 1 combined multi-channel ExE matrix, w/ distance based on K channels
  //       cat       1 combined KExKE matrix, i.e. concatenating channels/epochs 
   
  //d
  // For multiple signals, create a single matrix based on all channels 
  // *unliess* uni is specified, in which case write clusters separately
  // for each channel (i.e. output stratified by channel) 
  // -- do not allow mat and uni to be specified together/
  //
  
  bool univariate = param.has( "uni" );

  //
  // Alternatively, concatenate multiple signals to get, for N channels,
  // an NE x NE matrix, i.e. to jointly cluster channels and epochs together
  //

  bool ne_by_ne = param.has( "cat" );

  if ( ne_by_ne && univariate )
    Helper::halt( "cannot specify both uni and cat" );

    
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

  const bool no_annotations = true;

  signal_list_t signals = edf.header.signal_list( signal_label , no_annotations );
  
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
  // Requires data to be epoched, ubnless 'cat' mode and no epochs,se whole signal
  //

  bool use_whole_trace = false;
  
  if ( ! edf.timeline.epoched() )
    {
      if ( ne_by_ne )
	{
	  // if not epoched but in cat mode, set one single epoch
	  // spanning the whole night, i.e. to het a KxK matrix
	  // of only channels
	  use_whole_trace = true;
	  logger << "  clustering channels only, not epochs\n";
	}
      else
	Helper::halt( "ExE requires epoched data" );
    }

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


  logger << "  PDC-encoding with m=" <<m << ", t=" << t << "\n";  
  
  //
  // Outer loop (if in univariate mode || cat mode , otherwise exits after one pass)
  //
  
  for (int s0=0; s0<ns; s0++)
    {     

      //
      // Only consider data channels
      //

      if ( univariate || ne_by_ne )
	{
	  if ( edf.header.is_annotation_channel( signals(s0) ) ) continue;
	}
      
      //
      // Reset obs() unless accumulating over all channels (cat mode)
      //
      
      if ( ! ne_by_ne ) 
	clear();
      
      
      //
      // Add channels
      //
      
      if ( univariate )
	{
	  logger << "  calculating epoch-by-epoch distances for " << signals.label(s0) << "\n";	  
	  // just add one channel
	  add_channel( signals.label(s0) );
	  
	  // and stratify output by channel
	  if ( univariate )
	    writer.level( signals.label(s0) , globals::signal_strat );
	}
      else if ( ne_by_ne )
	{
	  logger << "  calculating epoch-by-epoch distances for " << signals.label(s0) << "\n";
	  // only have a single, dummy 'channe'
	  if ( obs.size() == 0 )
	    add_channel( "_dummy" );
	}
      else
	{
	  // add all channels (and we'll quit afterwards)
	  logger << "  calculating epoch-by-epoch distances for " ;
	  for ( int s=0; s<ns; s++ )
	    {	      
	      add_channel( signals.label(s) );
	      logger << signals.label(s) << " " ;
	    }
	  logger << "\n";

	}
      
      //
      // Get list of included epochs
      //
     
      int ne = use_whole_trace ? 1 : edf.timeline.first_epoch();
 
      std::vector<int> epochs;
      
      while ( 1 ) 
	{

	  int epoch = use_whole_trace ? -9 : edf.timeline.next_epoch();      

	  if ( epoch == -1 ) break;

	  interval_t interval = use_whole_trace ? edf.timeline.wholetrace() : edf.timeline.epoch( epoch );

	  // record for this epoch/signal
	  
	  pdc_obs_t ob( ns );

	  ob.id = use_whole_trace ? "0" : Helper::int2str( edf.timeline.display_epoch( epoch ) ); 

	  ob.label = signals.label(s0);
	  
	  epochs.push_back( use_whole_trace ? 0 : edf.timeline.display_epoch( epoch ) ) ; 
		  
	  //
	  // get signal(s)	  
	  //
	  
	  if ( univariate || ne_by_ne )
	    {
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
	  // If we read the entire trace, we can exit epoch loop now
	  //

	  if ( use_whole_trace ) break;
	  
	  //
	  // next epoch
	  //
	}

      //
      // Encode (all unencoded) time-series
      //
      
      encode_ts();
      
      //
      // Purge unnecessary TS after encoding as PD (helps w/ cat mode)
      //
      
      purge_ts();


      //
      // For either univariate or multivariat mode, we have all epochs and so calculate
      // solution
      //

      if ( ! ne_by_ne )
	{
	  
	  exe_calc_matrix_and_cluster( edf , param , write_matrix , outfile );	  
	  
	  //
	  // in multivariate-mode, we are all done after one pass
	  //
	  
	  if ( ! univariate ) break;
	  
	}
      
      //
      // next channel in univariate || NExNE mode
      //

    }


  //
  // Final call after we've built all NExNE obs 
  // 

  if ( ne_by_ne )
    exe_calc_matrix_and_cluster( edf , param , write_matrix , outfile );

  //
  // else all done
  //
    
  if ( univariate ) 
    writer.unlevel( globals::signal_strat );

}




void pdc_t::exe_calc_matrix_and_cluster( edf_t & edf , param_t & param ,
					 const bool write_matrix , const std::string & outfile )
{

  const bool ne_by_ne = param.has( "cat" );
  
  const int nobs = obs.size();
    
  const bool use_whole_trace = ne_by_ne && ! edf.timeline.epoched() ; 

  
  //
  // Calculate distance matrix
  //
  
  
  Data::Matrix<double> D = all_by_all();
  
  
  if ( D.dim1() != nobs ) 
    Helper::halt( "internal error in pdc_t::similarity_matrix()" );
  
  if ( write_matrix )
    {

      std::ofstream OUT1( outfile.c_str() , std::ios::out );
      for (int i=0;i<nobs;i++)
	{
	  for (int j=0;j<nobs;j++) OUT1 << ( j ? "\t" : "" ) << D[i][j];
	  OUT1 << "\n";
	}
      OUT1.close();

      logger << "  output distance matrix for " << nobs
	     << " observations to " << outfile << "\n";

      // channel/epoch labels?
      if ( ne_by_ne )
	{	  
	  std::ofstream OUT1( ( outfile+".idx").c_str() , std::ios::out );

	  OUT1 << "ID\tCH";

	  if ( ! use_whole_trace )
	    OUT1 << "\tE";

	  OUT1 << "\n";

	  for (int i=0;i<nobs;i++)
	    {
	      OUT1 << edf.id << "\t"		 
		   << obs[i].label ;
	      if ( ! use_whole_trace )
		OUT1 << "\t" << obs[i].id ;
	      OUT1 << "\n";
	    }
	  
	  OUT1.close();
	  
	  logger << "  paired with channel/epoch index: " << outfile << ".idx\n";
	  
	}
      
    }

  
  //
  // Cluster
  //
  
  if ( ne_by_ne )
    logger << "  clustering channels/epochs...\n";
  else
    logger << "  clustering epochs...\n";

  
  cluster_t cluster;

  // max. number of clusters (stopping rule)
  int preK = param.has( "k" ) ? param.requires_int( "k" ) : 0 ;

  // constraint on maximum size of each cluster (0 = no constraint)
  int maxS = param.has( "mx" ) ? param.requires_int( "mx" ) : 0 ;

  //
  // do we want to cluster?
  //

  if ( preK == 0 && maxS == 0 ) return;
    
  //
  // get cluster solution
  //

  cluster_solution_t sol = cluster.build( D , preK , maxS );
  
  if ( sol.best.size() != nobs )
    Helper::halt( "internal error in ExE" );


  //
  // Report output
  //
  
  for (int e=0; e<nobs; e++ )
    {
      
      // [ option: uni ] (report K different ExE]
      // 1 C4
      // 2 C4
      // 3 C4

      // [ default, multivariate ]  (reports single ExE based on K channels)
      // 1 C3,C4
      // 2 C3,C4
      // 3 C3.C4

      // option: concatenate (reports single KExKE)
      // 1 C3
      // 2 C3
      // 3 C3
      // 1 C4
      // 2 C4
      // 3 C4

      // if reporting in CAT mode, here we need to specify the CH level
      // for each observation
      
      if ( ne_by_ne )
	writer.level( obs[e].label , globals::signal_strat );

      //
      // Stratify by epoch?
      //

      if ( ! use_whole_trace )
	writer.epoch( edf.timeline.display_epoch( e ) );
      
      //
      // cluster 'label' is the exemplar epoch for that cluster
      //

      if ( use_whole_trace )
	writer.value( "CL" , sol.exemplars[ sol.best[e] ] ); // use channel # 
      else
	writer.value( "CL" , edf.timeline.display_epoch( sol.exemplars[ sol.best[e] ] ) );
      
    }

  //
  // tidy up
  //

  if ( ne_by_ne )
    writer.unlevel( globals::signal_strat );

  if ( ! use_whole_trace )
    writer.unepoch();
  
}

