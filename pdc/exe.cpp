
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
  // Signals and sample-rate
  //
  
  std::string signal_label = param.requires( "sig" );   

  const bool no_annotations = true;

  signal_list_t signals = edf.header.signal_list( signal_label , no_annotations );
  
  const int ns = signals.size();


  //
  // write cluster solution to normal output stream
  // write to a separate txt file, not an output-db
  //

  
  bool write_matrix = param.has( "mat" );
  
  if ( univariate && write_matrix && ns != 1 ) 
    Helper::halt( "cannot specify uni and mat together if >1 signal" );
  
  std::string outfile = "";
  if ( write_matrix ) 
    outfile = param.requires( "mat" ) ;


  
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
       
  
  if ( find_representatives )
    {
      
      // nobs = number of epochs
      // 

      // track splits 0 (all), 1, 2, etc
      std::map<int,std::set<int> > splits;
      std::vector<int> split( nobs , 0 );
      
      // initialize: all in one
      for (int i=0;i<nobs;i++)
	splits[0].insert( i );
      // std::cout << " splits.size() " << splits.size() << "\n";
      // std::cout << " ss " << splits[0].size() << "\n";

      // code for next splits
      int scnt = 0;
      
      // start loop
      while ( 1 )
	{
	  

	  //
	  // Within each split, get an Otsu threshold
	  //  - if we return 0 freqs for all, then we're done 

	  bool zero = true;

	  std::map<int,double > th;
	  //std::cout <<"---------------------------------------------------------\n";
	  std::map<int,std::set<int> >::const_iterator ss = splits.begin();
          while ( ss != splits.end() )
	    {
	      std::set<int>::const_iterator pp = ss->second.begin();
	      std::vector<double> x;
              while ( pp != ss->second.end() )
                {
                  std::set<int>::const_iterator qq = pp;
		  qq++;
		  while ( qq != ss->second.end() )
		    {
		      x.push_back( D[*pp][*qq] );
		      //std::cout << " intootsu  " << D[*pp][*qq] << "\n";
		      ++qq;
		    }
		  ++pp;
		}

	      double empf = 0;
	      std::map<double,double> fvals;
	      th[ ss->first ] = MiscMath::threshold2( x , &empf , 0 , &fvals );
	      if ( empf > 1e-8 ) zero = false; 
	      
	      //std::cout << "\n\notsu " << th[ ss->first ] << " " << empf << "\n";
	      ++ss;
	    }

	  //
	  // No splits?
	  //

	  if ( zero ) break;
	  

	  //
	  // Based on the Otsu split, derive estimate of neighbors, based on below th
	  // distances only - i.e. we want lots 
	  //
	       
	  std::vector<double> w( nobs , 0 );
	  
	  for (int i=0;i<nobs;i++)
	    {
	      const double & st = th[ split[i] ] ;

	      const std::set<int> & s = splits.find( split[i] )->second;
	      std::set<int>::const_iterator ss = s.begin();
	      while ( ss != s.end() )
		{
		  
		  const double & dst = D[i][*ss] ;
		  //std::cout << " st dst " << st << " " << dst << " " << w[i] << "\n";
		  if ( dst <= st ) 
		    w[i]++;
		  ++ss;		  
		}
	    }

	   // std::cout << "\n";
	   // for (int i=0; i<nobs; i++)
	   //   std::cout << " pre " << i << "\t" << split[i] << "\t" << w[i] << "\n";
	   // std::cout << "\n";


	  
	  // 2) Within each split, find the largest distance
	  
	  
	  std::map<int,std::pair<int,int> > didx;
	  std::map<int,double > d1;
	  
	  ss = splits.begin();
	  while ( ss != splits.end() )
	    {
	      //std::cout << " ss->first looking " << ss->second.size() << "\n";
	      
	      std::set<int>::const_iterator pp = ss->second.begin();
	      
	      double d0 = 0;
	      int ix = 0, jx = 0;
	      while ( pp != ss->second.end() )
		{
		  std::set<int>::const_iterator qq = pp;
		  ++qq; // skip self
		  
		  while ( qq != ss->second.end() )
		    {

		      // raw distance for this pair
		      double dst = D[*pp][*qq];

		      // weighting: smallest 'w' of the pair
		      double w1 = w[ *pp ] < w[*qq ] ? w[*pp] : w[*qq ] ; 

		      // statistic
		      
		      double st = dst * w1 ;

		      //std::cout << " dst w1 st " << dst <<  " " << w1 << " " << st << "\n";
		      
		      if ( st > d0 )
			{
			  d0 = st;
			  ix = *pp; jx = *qq;			 
			}
		      ++qq;
		    }
		  ++pp;
		}

	      //std::cout << "  putative split " << ss->first << " = " << ix << " " << jx << " " << d0 << "\n";
	      // save largest distance for this split
	      didx[ss->first] = std::make_pair( ix, jx );
	      d1[ss->first] = d0;

	      ++ss;
	    }


	  //
	  // Determine which split to make
	  //
	  
	  //std::cout << "\n\ncurrently " << splits.size() << " splits\n";

	  double mind = 0;
	  int mins = 0;

	  //	  std::cout << " d1 sz = " << d1.size() << "\n";
	  
	  std::map<int,double>::const_iterator uu = d1.begin();
	  while ( uu != d1.end() )
	    {
	      //std::cout << " uu-> split " << uu->first << " " << uu->second << " d = " << mind << "\n";
	      if ( uu->second > mind )
		{
		  mind = uu->second;
		  mins = uu->first; 
		}
	      ++uu;
	    }
	  
	  //std::cout << " will split split " << mins << "\n";
	  int ix = didx[mins].first;
	  int jx = didx[mins].second;

	  // the two new labels to add
	  // 0, 1+2, 3+4,
          int lab1 = ++scnt;
	  int lab2 = ++scnt;
	  //std::cout << " splitting " << mins << " --> " << lab1 << " + " << lab2 << "\n";
	  
	  //std::cout << " D(ix,jx) = [ " << ix << " , " << jx << " ]  " << D[ix][jx] << " " << w[ix] << " " << w[jx] << "\n";
	  // make the split
	  const std::set<int> & s1 = splits[ mins ];
	  std::set<int>::const_iterator pp = s1.begin();
	  while ( pp != s1.end() )
	    {
	  //std::cout << " ix jx " << ix << " " << jx << "\n";
	      // does 'p' go w/ i or j in this split?
	      if ( D[*pp][ix] < D[*pp][jx] )
		{
	          split[ *pp ] = lab1;
		  splits[ lab1 ].insert( *pp );
		  //	  std::cout << " adding " << lab1 << "\n";
		}
	      else
		{
		  split[ *pp ] = lab2;
                  splits[ lab2 ].insert( *pp );
		  // std::cout << " adding " << lab2 << "\n";
	         }       
	      
	      ++pp;
	    }

	  // finally, remove the partent split					
	  splits.erase( splits.find( mins ) );
	       
	  // all done, count splits
	  std::map<int,int> cnts;
	  for (int i=0; i<nobs; i++)
	    cnts[split[i]]++;

	  // print:

	  // std::cout << "\n\n";
	  // if ( 0 )
	  //   {
	  //     for (int i=0; i<nobs; i++)
	  // 	std::cout <<"  end " << i << "\t" << split[i] << "\t" << w[i] << "\n";
	  //     std::cout << "\n\n";
	  //   }
	  
	  // std::map<int,int>::const_iterator cc = cnts.begin();
	  // while ( cc != cnts.end() )
	  //   {
	  //     std::cout << " splitN " << cc->first << " = " << cc->second << "\n";
	  //     ++cc;
	  //   }

	  //std::cout <<" FIN " << cnts.size() << " " << find_representatives << "\n";
	  if ( cnts.size() >= find_representatives ) break;
	}
    
      //
      // given splits, now find the most representative epoch of each 
      //

      
      std::map<int,int> repe;
           
      //
      // Now go back and find best representative for this whole set 
      //

      std::map<int,std::set<int> >::const_iterator ss = splits.begin();
      while ( ss != splits.end() )
       {
	 const int s = ss->first;
	 
	 const std::set<int> & epochs = ss->second;

	 // set representative as the first
	 repe[ s ] = *(epochs.begin());

	 // but otherwise get median 
	 if ( epochs.size() > 1 )
	   {
	     
	     double d1 = 0;                                                                                                                  
	     int didx = 0;                                                                                                                   
	     
	     std::set<int>::const_iterator pp = epochs.begin();
	     while ( pp != epochs.end() )
	       {

		 // candidate rep for this split
		 const int c = *pp;

		 // get median dst to all others in this split
		 std::vector<double> ds;
		 std::set<int>::const_iterator qq = epochs.begin();
		 while ( qq != epochs.end() )
		   {
		     if ( pp != qq )
		       ds.push_back( D[*pp][*qq] );
		     ++qq;
		   }
		 double medd = MiscMath::median ( ds );

		 // got a new representative?
		 if ( pp == epochs.begin() || medd < d1 )
		   {
		     d1 = medd;
		     didx = *pp;
		   }
		 
		 ++pp;
	       }

	     // update repe
	     repe[ s ] = didx;
	   }
	 
	 // next split
	 ++ss;
       }

    
      //
      // outputs
      //

      const int q = splits.size();
      int k = 0;
      std::map<int,int> track;
      std::map<int,int>::const_iterator kk = repe.begin();
      while ( kk != repe.end() )
	{
	  writer.level( ++k , "K" );
	  track[ kk->first ] = k; // for output below
	  writer.value( "E" , obs[ kk->second ].id );                                                                                              
	  writer.value( "N" , (int)splits[ kk->first ].size() );
	  ++kk;
	}
      writer.unlevel( "K" );
      
      
      for (int i=0; i<nobs; i++)
       	{
	  // actual epoch label
	  int e = 0;
       	  if ( ! Helper::str2int( obs[ i ].id , &e ) )
       	    Helper::halt( "internal error in exe-rep" );

	  // class in split[i]
	  // representative in repe[ split[i] ] 
	  const std::string rep = obs[ repe[ split[i] ] ].id ;      

	  writer.epoch( e );	  
	  writer.value( "K" , track[ split[i] ] );
	  writer.value( "KE" , rep );
	}
      writer.unepoch();      
      
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

