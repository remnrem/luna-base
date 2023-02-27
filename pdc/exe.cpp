
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
  // *unless* uni is specified, in which case write clusters separately
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
    outfile = param.requires( "mat" ) ;

  
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
  // Requires data to be epoched, unless 'cat' mode and no epochs,se whole signal
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
      if ( univariate )
	Helper::halt( "cannot specify uni and entropy options together" );

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
  // Find representative epochs
  //
  
  const int find_representatives = param.has( "representative" ) ?
    ( param.empty("representative" ) ? 1 : param.requires_int( "representative" ) ) : 0 ; 

  // by defauly, use median similarity 
  const bool use_median = ! param.has( "sum" );

  if ( find_representatives )
    {

      // 1) get epoch most similar to all other epochs
      // 2) get epoch most dissimilar to this epoch
      // 3) for next j cases, find epoch min( all_i D_i + 
      
      int mini1 = 0;

      double mind = 0;

      for (int i=0;i<nobs;i++)
	{
	  double d1 = 0;
	  if ( use_median )
	    {
	      const std::vector<double> * dd = D.col(i).data_pointer();
	      
	      //const Vector<T> * col_pointer( const int c ) const { return &data[c]; }
	      d1 = MiscMath::median( *dd );
	    }
	  else // sum
	    for (int j=0;j<nobs;j++)
	      if ( j != i ) d1 += D[i][j];
	  
	  if ( i == 0 || d1 < mind ) 
	    {
	      mind = d1;
	      mini1 = i;
	    }
	}
      
      std::vector<int> es;
      std::set<int> ess;
      es.push_back( mini1 );
      ess.insert( mini1 );
      
      // 2) find least similar to this
      int miniN = 0;
      mind = 0;
      for (int i=0;i<nobs;i++)
	{
	  if ( D[i][mini1] > mind )
	    {
	      mind = D[i][mini1] ;
	      miniN = i;
	    }
	}
      
      //      std::cout	<< " eDIS = " << obs[miniN].id << "\n";
      es.push_back( miniN );
      ess.insert( miniN );
      
      // now find next epoch most similar to all others, but not similar
      // to existing picks

      for (int k=0; k<find_representatives; k++)
	{
	  int minik = 0;
	  
	  double mind = 0;
	  int nn = 0;
	  
	  for (int i=0;i<nobs;i++)
	    {
	      if ( ess.find( i ) != ess.end() )
		continue;
	      
	      // mean distance to all other points
	      double d1 = 0;
	      for (int j=0;j<nobs;j++)
		if ( j != i )
		  d1 += D[i][j];
	      d1 /= (double)( nobs-1 );
	      
	      // adjust for distance to closest existing
	      // points (i.e. subtract these values)
	      // so less likely to select if close to
	      // another existing point
	      double d2 = 0;
	      for ( int q=0; q<es.size(); q++)
		//if ( q != 1 )  // skip 'worst case' pick here initially
		  if ( q == 0 || d2 < D[i][es[q]] )
		    d2 = D[i][es[q]];
	      							 
	      d1 /= d2;
	      
	      if ( i == 0 || d1 < mind )
		{
		  mind = d1;
		  minik = i;
		}
	      
	    }

	  //		std::cout << " found next point rep " << minik << " at " << mind << "\n";
	      es.push_back( minik );
	      ess.insert( minik );
	    }

      //
      // now, for the k+2 epochs, assign every other epoch to closest
      //

      int q = es.size();
      std::vector<int> asgn( nobs , 0 );
      std::vector<std::vector<int> > q2a( q );
      
      // 0 = most likely
      // 1 = least likely
      // 2 ... others

      for (int i=0; i<nobs; i++)
	for (int j=1; j<q; j++)
	  {
	    if ( D[i][es[j]] < D[i][ es[ asgn[i] ] ] )
	      {
		//std::cout << " setting " << i << " -> " << j << "\n";
		asgn[i] = j;		
	      }
	  }
     
      std::map<int,int> ec;
      for (int i=0; i<nobs; i++)
	{
	  ec[ asgn[i] ]++;
	  q2a[ asgn[i] ].push_back( i );
	}

      //
      // Now go back and find best representative for this whole set 
      // (and edit es[] ) 
      //
	   
      for (int j=0; j<q; j++)
	{
	  const int nq = ec[ j ];
	  const std::vector<int> & e = q2a[ j ];
	  if ( e.size() != nq ) Helper::halt( "internal error" );

	  // requires this cluster to have at least 3 elements
	  // otherwise, keep prototype as is
	  if ( nq >= 3 )
	    {
	      double d1 = 0;
	      int didx = 0;
	      for (int p=0; p<nq; p++)
		{
		  std::vector<double> dst;
		  for (int q=0; q<nq; q++)
		    if ( p!=q )
		      {
			//std::cout << " adding p,q " << e[p] << " " << e[q] << " " << D[e[p]][e[q]] << "\n";
			dst.push_back( D[e[p]][e[q]] );
		      }
		  
		  double dd = MiscMath::median( dst );
		  //std::cout << "  DD = " << e[p] << " " << dd << "\n";
		  if ( p == 0 || dd < d1 )
		    {
		      d1 = dd;
		      didx = e[p];
		    }
		}
	      
	      // std::cout << " final min median= " << d1 << "\n";
	      
	      // std::cout << " grp j = " << j << "\n"
	      // 	    << " updating " << es[ j ] << " tp " << didx << "\n";
	      
	      // update es[]
	      es[ j ] = didx;
	    }
	}
      
      //
      // outputs
      //

      for (int j=0; j<q; j++)
	{
	  writer.level( j+1 , "K" );
	  writer.value( "E" , obs[ es[j] ].id );
	  writer.value( "N" , ec[ j ] );
	}
      writer.unlevel( "K" );

      for (int i=0; i<nobs; i++)
	{
	  int e = 0;
	  if ( ! Helper::str2int( obs[ i ].id , &e ) )
	    Helper::halt( "internal error in exe-rep" );
	  
	  writer.epoch( e );
	  writer.value( "K" , 1+ asgn[i] );
	  writer.value( "KE" , obs[ es[ asgn[i] ] ].id );
	}
      writer.unepoch();
      
      // for (int i=0;i<nobs;i++)
      //   {
      //     for (int j=0;j<nobs;j++) OUT1 << ( j ? "\t" : "" ) << D[i][j];

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

