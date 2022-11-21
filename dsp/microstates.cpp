
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

#include "dsp/microstates.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "db/db.h"
#include "helper/logger.h"
#include "stats/kmeans.h"
#include "dsp/lzw.h"
#include "dsp/mse.h"
#include "dsp/tv.h"
#include "miscmath/crandom.h"
#include "timeline/cache.h"

#include <limits>
#include <algorithm>

extern writer_t writer;
extern logger_t logger;

#include "stats/eigen_ops.h"

std::vector<char> ms_prototypes_t::ms_labels;

bool ms_prototypes_t::ul_groups;

std::map<int,int> ms_prototypes_t::ul_mapping;

//
// TODO:  to add annot from states;
//        to add peaks to a cache_t from states
//        to specify labels directly for a prototype file (i.e. if diff order)
//        option to make aggregare EDF out of prototypes (i.e. to subsequently cluster prototypes)


void dsptools::microstates( edf_t & edf , param_t & param )
{

  //
  // This function can be called in one of several modes
  //
  // Single-EDF mode: find peaks, segment, backfit, smooth then calculate stats
  //
  // Multi-sample mode:  1) find peaks, aggregating into a single EDF (e.g. peaks.edf)  [ 'peaks' ]
  //                     2) read peaks.edf, segment and save prototype file [ assume a single EDF given, not sample-list ]  
  //                     3) read prototype file(s), then backfit, smooth, caclulate stats for all the sample-list
  //

  bool multi_peaks = param.has( "peaks" );
  
  bool multi_segment = param.has( "segment" ); // but will point to a single EDF (i.e. containing peaks from all EDFs)
  
  bool multi_backfit = param.has( "backfit" );

  // if none of the above three optinos, assume single-sample mode
  bool single_sample = ! ( multi_peaks || multi_segment || multi_backfit );

  bool skip_peaks = param.has( "all-points" );

  bool epoch = param.has( "epoch" );
  if ( epoch && ! multi_backfit )
    Helper::halt( "can only specify epoch when running in backfit mode" );
  if ( epoch && param.has( "dump-gfp" ) )
    Helper::halt( "cannot specify epoch and dump-gfp (to dump sample-level GFP) together" );
  if ( epoch && param.has( "write-states" ) )
    Helper::halt( "cannot specify epoch and write-states (to dump state order) together" );
  
  int run_kmers = param.has( "kmers" ) ? param.requires_int( "kmers" ) : 0 ;
  
  bool add_annot = param.has( "add-annot" );
  const std::string annot_tag = add_annot ? param.value( "add-annot" ) : "" ;

  bool add_sig = param.has( "add-sig" );
  const std::string sig_tag = add_sig ? param.value( "add-sig" ) : "" ;
  
  bool add_corrs = param.has( "add-spc-sig" );
  const std::string sig_corr_tag = add_corrs ? param.value( "add-spc-sig" ) : "" ; 

  bool add_conf = param.has( "add-conf" ); // only applies for 'add-spc-sig'
  const std::string sig_conf = add_conf ? param.value( "add-conf" ) : "" ;
  
  bool save_transitions = param.has( "cache" );
  const std::string cache_name = save_transitions ? param.value( "cache" ) : "" ; 
  
  if ( ( add_annot || add_sig || add_corrs ) && epoch ) 
    Helper::halt( "cannot use add-annot or add-sig in epoch mode" );

  //
  // Group-specific (two groups only: upper/lower)
  // microstates (i.e. pairs 'A' and 'a')
  //

  ms_prototypes_t::ul_groups = param.yesno( "grouped" ) ;
  
  // internally, code using (e.g.)
  // '1' for 'A/a', '2' for 'B/b' etc
  // nb. assumes letter encoding for states
  // std::map<char,int> state2group; 
  // std::map<int,char> grp2state;

  //
  // Channels
  //
  
  const bool no_annotations = true;  
  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) , no_annotations );  

  const int ns = signals.size();

  //
  // Check for equal sample rates
  //
  
  if ( ns < 2 ) return;
  
  int sr = edf.header.sampling_freq( signals(0) );
  
  for (int i=1;i<ns;i++)
    {      
      if ( edf.header.sampling_freq( signals(i) ) != sr )
	Helper::halt( "all signals must have similar SR for MS" );      
    }
  

  //
  // If in epoch mode, only read prototypes once
  //

  ms_prototypes_t prior_prototypes;
    
  if ( multi_backfit && epoch )
    {
      const std::string filename = param.value( "backfit" ); 
      prior_prototypes.read( filename );
      
      // make UL mapping?
      if ( ms_prototypes_t::ul_groups )
	prior_prototypes.make_ul_map();
      
      // check channels line up      
      if ( signals.size() != prior_prototypes.C )
	Helper::halt( "number of channels in " + filename + " does not match current signal selection" );

      for (int s=0; s<signals.size(); s++)
	if ( ! Helper::iequals( signals.label(s) , prior_prototypes.chs[s] ) )
	  Helper::halt( signals.label(s) + " does not match " + prior_prototypes.chs[s]  + " for slot " + Helper::int2str( s+1) );      
    }

  
  //
  // Fetch sample matrix: typically whole trace, but if in backfit mode
  // could be epoch level
  //

  int ne = edf.timeline.first_epoch();

  while ( 1 )
    {  
      
      interval_t interval;
      
      if ( epoch )
	{
	  int e = edf.timeline.next_epoch();      
       	  if ( e == -1 ) break;      
	  interval = edf.timeline.epoch( e );
	  writer.epoch( edf.timeline.display_epoch( e ) );
	  logger << " -- processing epoch " <<  edf.timeline.display_epoch( e ) << " (" << e+1 << " of " << ne << ")\n";
	}
      else
	interval = edf.timeline.wholetrace();

      //
      // Get actual signals
      //

      matslice_t mslice( edf , signals , interval );
      
      const Data::Matrix<double> & X = mslice.data_ref();
  
      //
      // Set up microstate class
      //
          
      microstates_t mstates( param , edf.id , sr );      
      

      //
      // Find peaks?
      //
      
      std::vector<int> peaks;
      
      if ( ( single_sample && ! skip_peaks ) || multi_peaks )
	{
	  logger << "  find GPF peaks\n";
	  peaks = mstates.find_peaks( X , signals );      
	  
	  // in multi-sample mode, just save peak (data) and move on
	  if ( multi_peaks )
	    {
	      microstates_t::aggregate2edf( X , signals, peaks , sr ,
					    param.requires_dbl( "pmin" ) ,
					    param.requires_dbl( "pmax" ) , 
					    param.value( "peaks" ) );
	      return;
	    }
	}
      
  
      //
      // Segmentatation based on peaks
      //
      
      ms_prototypes_t prototypes;
      
      if ( single_sample || multi_segment )
	{
	  logger << "  segmenting peaks to microstates\n";
	  
	  // optional 'canonical' file to use in assigning labels
          const bool has_canonicals = param.has( "canonical" ) ;	  
	  std::string canonical_file = has_canonicals ? param.value( "canonical" ) : "" ;


	  // nb, if peaks is empty, just takes all rows
	  // of X;  i.e. if coming from an aggregate peak EDF
	  
	  prototypes = mstates.segment( X , signals , peaks , 
					has_canonicals ? &canonical_file : NULL );

	  // re-assign labels based on a canonical set? 
	  // (based on spatial correlations) 
	  if ( param.has( "canonical" ) )
	    {
	      const std::string filename = Helper::expand( param.value( "canonical" ) );
	      
	    }
	  
	  
	  // In multi-sample mode, just save prototypes and move on
	  
	  if ( multi_segment )
	    {
	      const std::string filename = Helper::expand( param.value( "segment" ) );
	      prototypes.write( filename );
	      return;
	    }
	}
      

      //
      // Or read prototypes from a previous segmentation, and check that channels match
      //
      
      if ( multi_backfit )
	{
	 if ( epoch )
	   {
	     // already read... no need to re-read
	     // ms_labels will stay as is
	     prototypes = prior_prototypes;
	   }
	 else
	   {
	     const std::string filename = param.value( "backfit" );
	     prototypes.read( filename );

	     // make UL mapping?
	     if ( ms_prototypes_t::ul_groups )
	       prototypes.make_ul_map();
	     
	     // check channels line up 
	     if ( signals.size() != prototypes.C )
	       Helper::halt( "number of channels in " + filename + " does not match current signal selection" );
	     for (int s=0; s<signals.size(); s++)
	       if ( ! Helper::iequals( signals.label(s) , prototypes.chs[s] ) )
		 Helper::halt( signals.label(s) + " does not match " 
			       + prototypes.chs[s]  + " for slot " + Helper::int2str( s+1) );
	   }
       }

      
      //
      // Backfitting (w/ optional smoothing via TV)
      //
      
      const bool store_GMD = true; 
      
      logger << "  back-fitting solution to all time points\n";
      
      ms_backfit_t bf = mstates.backfit( Statistics::transpose( X ) , 
					 microstates_t::eig2mat( prototypes.A ) , 
					 param.has( "lambda" ) ? param.requires_dbl( "lambda" ) : 0 , 
					 store_GMD );
      
      //
      // Set ambiguous points to missing?
      //

      if ( param.has( "ambig" ) )
	{
	  std::vector<double> pr = param.dblvector( "ambig" );
	  //if ( ! ( pr.size() == 1 || pr.size() == 2 ) ) Helper::halt( "ambig expects two args" );
	  if ( pr.size() != 2 ) Helper::halt( "ambig expects two args" );

	  // e.g. 1.1 --> 10% increase in best over second best
	  const double th1 = pr[0] ;

	  // absolute value of max SPC (e.g. needs to be at least 0.5) 
	  const double th2 = pr[1] ;

	  // NOT IMPLEMENTED YET
	  // if ambiguous intervals < this threshold AND they are flanked by 
	  // two similar states, then fill-in this ambiguous region to be the
	  // same as the flanking regions
	  
	  // default = 20 msec (4 samples)
	  // const double fillin_msec = pr.size() == 2 ? pr[1] : 20 ;
	  // const int fillin_samples = round( fillin_msec * sr/1000.0 );
	  
	  bf.determine_ambiguity( th1 , th2 ); //, fillin_samples ); 
	}

  
      //
      // Smoothing / rejection of small intervals
      //
      
      ms_backfit_t smoothed = bf;
  
      if ( param.has( "min-msec" ) ) 
	{
	  double minTime_msec = param.has( "min-msec" ) ? param.requires_dbl( "min-msec" ) : 20 ; 
	  
	  int minTime_samples = round( minTime_msec * sr/1000.0 );
	  
	  logger << "  smoothing: rejecting segments <= " << minTime_msec << " msec\n";
	  
	  smoothed = mstates.smooth_reject( bf , minTime_samples );
	  
	} 
      else
	logger << "  no minimum segment duration specified (e.g. min-msec=20)\n";

      
      
      //
      // Final stats
      //
      
      logger << "  getting final microstate statistics\n";
      
      ms_stats_t stats = mstates.stats( Statistics::transpose( X ) , 
					microstates_t::eig2mat( prototypes.A ) , 
					smoothed.best() ,
					smoothed.best_unreduced() 
					);
      

      //
      // Verbose dumping of GFP and L point by point? (nb. this only works in whole-trace mode)
      //

      std::string dump_file = param.has( "gfp" ) ? param.value( "gfp" ) : "" ;
      
      if ( dump_file != "" )
	{

	  logger << "  dumping GFP and states to " << dump_file << "\n";
	  std::ofstream O1( Helper::expand( dump_file ).c_str() , std::ios::out );
	  
           
	  // get TP...
	  // for now, just raw stats
	  
	  const int N = stats.GFP.size();
	  std::vector<int> states = smoothed.best();
	  if ( states.size() != N ) Helper::halt( "hmmm" );
	  for (int i=0;i<N;i++)
	    O1 << ( states[i] != -1 ? ms_prototypes_t::ms_labels[ states[i] ] : '?' ) << "\t"
	       << stats.GFP[i] << "\n";
	  O1.close();
	}
    
      
      //
      // Add new annotations and/or channels
      //
      
      if ( add_annot || add_sig || add_corrs || save_transitions )
	{
	  
	  std::vector<int> states = smoothed.best();
	  const int N = states.size();

	  //
	  // Add annotations
	  //
	  
	  if ( add_annot )
	    {
	      
	      const std::vector<uint64_t> * tp = mslice.ptimepoints();
	      if ( tp->size() != N ) 
		Helper::halt( "internal error in add-annot" );
	      
	      // for each class (+ ambig)
	      std::vector<annot_t*> k2a( prototypes.K + 1 );
	      std::vector<std::string> k2l( prototypes.K + 1 );

	      for (int k=0; k<prototypes.K + 1 ; k++)
		{
		  std::string s="?";
		  if ( k < prototypes.K )
		    s[0] = ms_prototypes_t::ms_labels[k];
                  
		  const std::string alab = annot_tag + s;                  
		  k2l[k] = alab;
		  logger << "  adding annotation " << alab << "\n";		  
		  annot_t * a = edf.timeline.annotations.add( alab );
		  a->description = "MS " + alab;
		  k2a[k] = a;
		}
	
	      // iterate over points
	      int curr = states[0];
	      int idx = 0;
	      	      
	      for (int i=1;i<N; i++)
		{
		  // point marks end of a state?
		  if ( states[i] != curr || i == N-1 )
		    {
		      // i.e. annotation is right up to the start of the next one
		      // so end == start of new internal in end+1 encoding
		      interval_t interval( (*tp)[idx] , (*tp)[i] );
		      
		      // add annot to the appropriate annotation (decoding -1 to K for ambig states)
		      k2a[ curr == -1 ? prototypes.K : curr  ]->add( k2l[ curr == -1 ? prototypes.K : curr ] , interval , "." );

		      // reset new state
		      curr = states[i];
		      idx = i;
		    }
		}	      
	    }
	  
	  //
	  // Add as signals
	  //
	  
	  if ( add_sig || add_corrs ) 
	    {
	      
	      // fetch spatial correlations (from GMDs)  [ add-spc-sig ]
	      Data::Matrix<double> SPC;	      
	      
	      // optionally, add 'confidence' score of max versus next best
	      std::vector<double> conf;

	      if ( add_corrs )
		{		  
		  SPC = bf.GMD; // K x N matrix (i.e. original, pre-smoothing)

		  for (int k=0; k< prototypes.K; k++)		    
		    for (int i=0; i<N; i++) 
		      SPC(k,i) = 1 - ( SPC(k,i) * SPC(k,i) ) / 2.0; 		  
		  
		  if ( add_conf ) 
		    {
		      conf.resize( N );
		      for (int i=0; i<N; i++)
			{ 
			  double best = 0;
			  double next = 0;
			  for (int k=0; k<prototypes.K; k++)
			    {
			      if ( SPC(k,i) > best )
				{
				  next = best;
				  best = SPC(k,i);				  
				}
			      else if ( SPC(k,i) > next )
				next = SPC(k,i); 
			    }			  
			  
			  // set to 100 as max
			  conf[i] = next > 0 ? best / next : 100; 
			  if ( conf[i] > 100 ) conf[i] = 100;
			}
		      
		      logger << "  adding confidence channel " << sig_conf << "\n";
		      edf.add_signal( sig_conf , sr , conf );

		    }
		  
		}

	      //
	      // for SPC and binary flags (new channels), only do this for
	      // unambiguous states currently
	      //

	      logger << "  adding" ;
	      if ( add_corrs ) logger << " spatial correlations";
	      if ( add_sig ) logger << ( add_corrs ? "," : "" ) << " assigned states";
	      logger << " as channels (" << sr << "Hz) :";

	      for (int k=0; k<prototypes.K; k++)
		{
		  std::string s="?";
		  s[0] = ms_prototypes_t::ms_labels[k];
		  const std::string clab = sig_tag + s;		  
		  const std::string clab_corrs = sig_corr_tag + s;		  
		  
		  if ( add_sig ) logger << " " << clab;
		  if ( add_corrs ) logger << " " << clab_corrs;
		  
		  std::vector<double> dat( N );
		  
		  //
		  // add spatial correlations as new channel in EDF
		  //

		  if ( add_corrs )
		    {
		      for (int i=0; i<N; i++) 
			dat[i] = SPC(k,i);
		      
		      edf.add_signal( clab_corrs , sr , dat );

		    }

		  //
		  // add binary 0/1 for assigned state in EDF
		  //
		  
		  if ( add_sig )
		    {
		      for (int i=0; i<N; i++) 
			dat[i] = states[i] == k ? 1 : 0 ;		      
		      edf.add_signal( clab , sr , dat );
		    }
		  
		}
	      logger << "\n";
	    }
	  
	  
	  //
	  // Save transitions
	  //
	  
	  if ( save_transitions )
	    {
	      
	      cache_t<int> * cache = edf.timeline.cache.find_int( cache_name );
	      //cache_t<int> * cache = edf.timeline.cache.find_int( cache_name );
	      
	      int idx = 0;
	      
	      // save first sample-point in a new state
	      std::vector<int> transitions;
	      for (int i=1; i<N; i++) 
		if ( states[i] != states[i-1] ) 
		  {
		    transitions.push_back( i );
		  }
	      
	      // cache sample-points: transitions 
	      cache->add( ckey_t( "points" , writer.faclvl() ) , transitions );
	    }
	  	
	}

  
      //
      // Output all final stats
      //
      // By class 'K'
      //
      
      std::map<int,std::pair<int,double> > cnts = microstates_t::counts( smoothed.best() );
      
      for (int k=0; k < prototypes.K ; k++)
	{
	  std::string s="?";
	  s[0] = ms_prototypes_t::ms_labels[k];
	  writer.level( s , "K" );
	  writer.value( "GFP" , stats.m_gfp[k] );
	  
	  writer.value( "COV" , stats.m_cov[k] ); // denom = all time-points
	  writer.value( "COV2" , stats.m_cov_unambig[k] ); // denom = all unambiguous time-points	  

	  writer.value( "OCC" , stats.m_occ[k] );
	  writer.value( "OCC2" , stats.m_occ_unambig[k] );
	  		  	  
	  writer.value( "DUR" , stats.m_dur[k] );

	  writer.value( "WGT" , stats.m_wcov[k] ); // mean SPC for all samples (i.e. probabilistic coverage)
	  writer.value( "SPC" , stats.m_spc[k] );  // mean SPC for K=k assigned maps
	  
	  writer.value( "GEV" , stats.m_gev[k] );      
	  writer.value( "N" , cnts[k].first );
	  //writer.value( "F" , cnts[k].second );
	}
      writer.unlevel( "K" );
      
    
  
      //
      // State transition probabilities: up to K+1 for ambig
      //
      
      for (int k=0; k < prototypes.K + 1 ; k++)
	{
	  std::string s1="?";
	  if ( k < prototypes.K ) 
	    s1[0] = ms_prototypes_t::ms_labels[k];
	  
	  writer.level( s1 , "PRE" );
	  for (int k2=0; k2 < prototypes.K + 1 ; k2++)
	    {
	      if ( k != k2 )
		{
		  std::string s2="?";
		  if ( k2 < prototypes.K ) 
		    s2[0] = ms_prototypes_t::ms_labels[k2];
		  
		  writer.level( s2 , "POST" );
		  writer.value( "P" , stats.tr(k,k2) );
		}
	    }
	}
      writer.unlevel( "PRE" );
      writer.unlevel( "POST" );
      
      
      // Overall stats: complexity
      
      writer.value( "LZW"     , stats.lwz_states );

      writer.value( "GEV"     , stats.GEV_tot );

      // Entropy: skip for now
      
      // for (int m=1;m<=8;m++)
      // 	writer.value( "SE" + Helper::int2str(m)  , stats.samplen[m] );

      
      // kmer stats: single obs, so no group differences here
    
      if ( run_kmers )
	{

	  //
	  // output by obserevd sequence (OBS and OBS / GROUP ) 
	  //

	  std::map<std::string,double>::const_iterator pp = stats.kmers.basic.pval.begin();
	  while ( pp != stats.kmers.basic.pval.end() )
	    {
	      writer.level( (int)pp->first.size() , "L" );
	      writer.level( pp->first , "S" );
	      	      
	      writer.value( "NG"  , stats.kmers.equiv_set_size[ pp->first ] );
	      writer.value( "SG"  , stats.kmers.obs2equiv[ pp->first ] );
	      
	      writer.value( "OBS" , stats.kmers.basic.obs[ pp->first ] );
	      writer.value( "EXP" , stats.kmers.basic.exp[ pp->first ] );
	      writer.value( "P" , pp->second );
	      writer.value( "Z" , stats.kmers.basic.zscr[ pp->first ] );

	      if ( stats.kmers.equiv_set_size[ pp->first ] > 1 )
		{
		  writer.value( "W_OBS" , stats.kmers.equiv.obs[ pp->first ] );
		  writer.value( "W_EXP" , stats.kmers.equiv.exp[ pp->first ] );
		  writer.value( "W_P" , stats.kmers.equiv.pval[ pp->first ] );
		  writer.value( "W_Z" , stats.kmers.equiv.zscr[ pp->first ] );
		}
	      
	      ++pp;
	    }  
	  
	  writer.unlevel( "S" );
	  writer.unlevel( "L" );

	  //
	  // Output by EQ group 
	  //
	  
	  pp = stats.kmers.group.pval.begin();
	  while ( pp != stats.kmers.group.pval.end() )
	    {
	      writer.level( (int)pp->first.size() , "L" );	      
	      writer.level( pp->first , "SG" );
	      writer.value( "NG"  , stats.kmers.equiv_set_size[ pp->first ] );	      
	      writer.value( "OBS" , stats.kmers.group.obs[ pp->first ] );
	      writer.value( "EXP" , stats.kmers.group.exp[ pp->first ] );
	      writer.value( "P" , pp->second );
	      writer.value( "Z" , stats.kmers.group.zscr[ pp->first ] );	      
	      ++pp;
	    }  
	  writer.unlevel( "SG" );
	  writer.unlevel( "L" );
	
	}
      

      //
      // either next epoch or all done
      //
      
      if ( ! epoch ) break;
    }

  if ( epoch ) writer.unepoch();

  
  //
  // All done for MS analysis
  //
  
}


microstates_t::microstates_t( param_t & param , const std::string & subj_id_, const int sr_ ) 
{

  //
  // Microstate analysis parameters
  //

  sr = sr_;

  subj_id = subj_id_;
  
  multi_peaks = param.has( "peaks" );
  
  multi_segment = param.has( "segment" ); // but will point to a single EDF (i.e. containing peaks from all EDFs)
  
  multi_backfit = param.has( "backfit" );

  // if none of the above three optinos, assume single-sample mode
  single_sample = ! ( multi_peaks || multi_segment || multi_backfit );

  if ( (int)multi_peaks + (int)multi_segment + (int)multi_backfit > 1 )
    Helper::halt( "cannot specify more than one of: peaks, segment and backfit" );

  // range of classes, if segmenting
  if ( single_sample || multi_segment ) 
    {
      if ( ! param.has( "k" ) ) Helper::halt( "requires k to be specified" );
      ks = param.intvector( "k" );
    }
  
  dump_file = param.has( "dump-gfp" ) ? param.value( "dump-gfp" ) : "";
  
  standardize = param.has( "standardize" );
  
  verbose = param.has( "verbose" );

  // write sequence to file? 
  statesfile = param.has( "write-states" ) ? param.value( "write-states" ) : "" ;
  
  // all points (i.e. just just GFP peaks)
  skip_peaks = param.has( "all-points" );
  
  // reject peaks > T times std(GFP) above the mean GFP if (T>0)
  gfp_max_threshold = param.has( "gfp-max" ) ? param.requires_dbl( "gfp-max" ) : 0 ;

  // reject peaks < T times std(GFP) below the mean GFP (if T>0)
  gfp_min_threshold = param.has( "gfp-min" ) ? param.requires_dbl( "gfp-min" ) : 0 ; 

  // reject peaks if skewness is > T SD from the mean kurtosis of all GFP peaks (if T>0)
  gfp_kurt_threshold = param.has( "gfp-kurt" ) ? param.requires_dbl( "gfp-kurt" ) : 0 ;

  // if > 0 , select (randomly) only this many peaks per observation
  restrict_npeaks = param.has( "npeaks" ) ? param.requires_int( "npeaks" ) : 0; 

  // not imlpemented yet
  // select only peaks this far apart (to ensure distinct peaks)
  min_peak_dist = 0;    

  if ( param.has( "kmers" ) )
    {
      std::vector<int> k = param.intvector( "kmers" );
      if ( k.size() != 3 && k.size() != 4 ) Helper::halt( "expecting 3 or 4 args for kmers=min,max,nreps(,w)" );
      kmers_min = k[0];
      kmers_max = k[1];
      kmers_nreps = k[2];
      kmers_w = k.size() == 4 ? k[3] : 0 ;  
    }
  else
    {
      kmers_nreps = 0;
    }

  
}



std::vector<int> microstates_t::find_peaks( const Data::Matrix<double> & X , 
					    const signal_list_t & signals )
{
  
  //
  // Global field power 
  //

  // nb. X not transposed here
  const int np = X.dim1();
  const int nc = X.dim2();
  
  logger << "  calculating GFP for sample\n";

  Data::Vector<double> GFP( np );

  for (int i=0; i<np; i++)
    {
      // get time-points across channels
      Data::Vector<double> p = X.row( i );    
      // get SD of raw data
      GFP[i] = sqrt( Statistics::variance( p , 0 ) ); // use N denom      
    }
  
  //
  // Find peaks in GFP
  //

  std::vector<int> peak_idx;
  int n_peaks = 0;
  
  for (int i=1; i<(np-1); i++)
    {
      if ( GFP[i] > GFP[i-1] && GFP[i] > GFP[i+1] ) 
	{	  
	  peak_idx.push_back(i);
	  ++n_peaks;
	}
    }
      

  //
  // GFP thresholding: round 1) min/max values
  //

  if ( gfp_max_threshold > 0 || gfp_min_threshold > 0 )
    {
      
      Data::Vector<double> peak_gfp( n_peaks );
      
      for (int r=0; r<n_peaks; r++)
	peak_gfp[r] = GFP[ peak_idx[r] ];
      
      double mean = Statistics::mean( peak_gfp );
      double sd = sqrt( Statistics::variance( peak_gfp , 1 ) ); // use N-1 denom
      double th_max = mean + gfp_max_threshold * sd;
      double th_min = mean - gfp_min_threshold * sd;
      
      bool has_max_th = gfp_max_threshold > 0 ; 
      bool has_min_th = gfp_min_threshold > 0 ;
      int cnt_min = 0 , cnt_max = 0;
      std::vector<int> peak_idx2;
      for (int r=0; r<n_peaks; r++)
	{
	  bool okay = true;
	  if ( has_max_th && GFP[ peak_idx[r] ] > th_max ) { ++cnt_max; okay = false; } 
	  if ( has_min_th && GFP[ peak_idx[r] ] < th_min ) { ++cnt_min; okay = false; }
	  if ( okay ) 
	    peak_idx2.push_back( peak_idx[r] );
	}
      
      logger << "  given mean GFP of " << mean << ", applying threshold to require:\n";
      if ( has_max_th ) logger << "  - GFP < " << th_max << " [ mean(GFP) + " << gfp_max_threshold << " * SD(GFP) ]\n";
      if ( has_min_th ) logger << "  - GFP > " << th_min << " [ mean(GFP) - " << gfp_min_threshold << " * SD(GFP) ]\n";
      logger << "  keeping " << peak_idx2.size() << " of " << n_peaks << " peaks";
      if ( has_max_th ) logger << ", dropping " << cnt_max << " for gfp-max";
      if ( has_min_th ) logger << ", dropping " << cnt_min << " for gfp-min";
      logger << "\n";
      
      writer.value( "GFP_MEAN" , mean );
      writer.value( "GFP_SD" , sd );
      
      writer.value( "NP0" , (int)peak_idx.size() );

      n_peaks = peak_idx2.size();
      peak_idx = peak_idx2;
      
      
    }

  
  //
  // Round 2: exclude peaks based on abberrant skewness/mean
  //
  
  //  bool has_skew_th = gfp_skew_threshold > 0;
  bool has_kurt_th = gfp_kurt_threshold > 0;
  
  if ( has_kurt_th )
    {
      
      std::vector<double> k( n_peaks );
      
      for (int i=0; i<n_peaks; i++)
	{
	  
	  const int idx = peak_idx[i];

	  // get time-points across channels
	  Data::Vector<double> p = X.row( idx );
	  
	  const int n = p.size();

	  // should be 0, but just in case
	  double mn = Statistics::mean( p );
	  for (int j=0; j<n; j++)
	    p[j] -= mn;

	  // get kurtosis (assumes mean = 0)
	  double numer = 0 , denom = 0;
	  for (int j=0; j<n; j++)
	    {
	      numer += pow( p[j] , 4 );
	      denom += pow( p[j] , 2 );
	    }
	  
	  numer /= (double)n;
	  denom /= (double)n;
	  denom *= denom;
	  
	  k[i] = numer / denom - 3.0;

	} // next peak
      
      // get mean/SD of kurtosis
      
      double mean = Statistics::mean( k );
      double sd = sqrt( Statistics::variance( k  , 1 ) ); // use N-1 denom
      double th_kurt = mean + gfp_kurt_threshold * sd;

      std::vector<int> peak_idx2;
      for (int r=0; r<n_peaks; r++)
	if ( k[r] <= th_kurt ) peak_idx2.push_back( peak_idx[r] );

      logger << "  applying GFP kurtosis threshold, to require "
	     << ": kurtosis < mean(kurtosis) + " << gfp_kurt_threshold << " * SD(kurtosis)\n";
      logger << "  dropping " << peak_idx.size() - peak_idx2.size() << " GFP peaks, to leave " << peak_idx2.size() << "\n";

      n_peaks = peak_idx2.size();
      peak_idx = peak_idx2;
      
    }


  //
  // Round 3: Only select (at most) N peaks?
  //

  if ( restrict_npeaks > 0 )
    {
      if ( n_peaks > restrict_npeaks )
	{
	  std::vector<int> a( restrict_npeaks );
	  CRandom::random_draw( a );
	  std::vector<int> peak_idx2 = peak_idx;
	  peak_idx.clear();
	  
	  for (int r=0; r<restrict_npeaks; r++)
	    peak_idx.push_back( peak_idx2[ a[r] ] );
	  n_peaks = peak_idx.size();
	  logger << "  randomly selected " << restrict_npeaks << " of " << peak_idx2.size() << " peaks\n";
	}
    }

  //
  // final # GFP peaks 
  //

  writer.value( "NP" , (int)peak_idx.size() );

  //
  // Output full GFP
  //
  
  if ( verbose )
    {
      for (int i=0; i<peak_idx.size(); i++)
	{
	  writer.level( peak_idx[i] , "SP" );
	  writer.value( "GFP" , GFP[ peak_idx[i] ] );
	}
      writer.unlevel( "SP" );
    }

  //
  // All done
  //

  logger << "  extracted " << n_peaks << " peaks from " << np << " samples ("
	 << round( 100 * ( n_peaks / (double)np ) ) << "%)\n";  
  
  return peak_idx;
  
}

ms_prototypes_t microstates_t::segment( const Data::Matrix<double> & X , 
					const signal_list_t & signals ,
					const std::vector<int> & peak_idx , 
					const std::string * canonical_file )
{

  
  //
  // Copy subset of data (GFP peaks) prior to clustering?x
  //

  bool has_peak_list = peak_idx.size() > 0 ;
  
  Data::Matrix<double> Z( has_peak_list ? peak_idx.size() : X.dim1() , X.dim2() );
  
  if ( has_peak_list )
    {
      const int n_peaks = peak_idx.size();
      const int nc = X.dim2();
      for (int r=0; r<n_peaks; r++)
	for (int c=0; c<nc; c++)
	  Z( r , c ) = X( peak_idx[r] , c ) ;
    }
  else
    Z = X;

    
  //
  // Standardize values?
  //

  if ( standardize )
    Statistics::standardize( Z );

  //
  // Optionally, dump input prior to clustering? ifnore for now;  dump_file
  // will dump GFP instead (see below) 
  //

  if ( 0 || dump_file != "" )
    {
      logger << "  dumping raw matrix to " << dump_file << "\n";
      std::ofstream O1( dump_file.c_str() , std::ios::out );      
      O1 << Z.dump();       
      O1.close();      
    }
    

  //
  // Modified K-Means clustering for microstates
  //

  modkmeans_t kmeans( ks , false , 10 , 1000 , 1e-6 , verbose );

  modkmeans_all_out_t results = kmeans.fit( Z );


  //
  // optimal K selected
  //
  
  writer.value( "OPT_K" , results.K );

  //
  // Set default labels: 1, 2, 3, etc (i.e. not 'canonical')
  //

  ms_prototypes_t::ms_labels.resize( results.K , '?' );
  for (int i=0; i<results.K; i++)
    ms_prototypes_t::ms_labels[i] = (char)(49 + i);
  


  //
  // Normalize A results
  //
  
  eigen_ops::scale( results.A , true , true );


  //
  // Or, overwrite with 'best guess' labels given a canonical file;
  // NOTE: this only changes the 'optimal' labels (and not same K in the
  // more verbose output);   if states are written with 'sol' they will 
  // have these assigned labels
  //

  ms_prototypes_t prototypes( signals , results.A ) ;
  
  if ( canonical_file != NULL ) 
    prototypes.map_to_canonicals( *canonical_file ) ;

  //
  // Maps
  //

  const int C = Z.dim2();
  const int N = Z.dim1();

  //
  // Optimal prototype maps
  //

  for (int i=0; i<C; i++)
    {
      writer.level( signals.label(i) , globals::signal_strat );

      // optimal solution
      for (int j=0; j<results.K; j++)
	{	  
	  std::string s="?";
	  s[0] = ms_prototypes_t::ms_labels[j];
	  writer.level( s , "K" );
	  writer.value( "A" , results.A(i,j) );
	}
      writer.unlevel( "K" );	    
      
    }
  writer.unlevel( globals::signal_strat );

  
  //
  // Correlations between maps (optimal solution)
  //

  for (int k1=0; k1<results.K; k1++)
    {
      std::string s="?";
      s[0] = ms_prototypes_t::ms_labels[k1];
      writer.level( s , "K1" );

      for (int k2=0; k2<results.K; k2++)
	{
	  std::string s2="?";
	  s2[0] = ms_prototypes_t::ms_labels[k2];
	  writer.level( s2 , "K2" );      
      
	  writer.value( "SPC" , ms_prototypes_t::spatial_correlation( results.A.col(k1) , results.A.col(k2) ) );
	}
      writer.unlevel( "K2" );
    }
  writer.unlevel( "K1" );
  

  //
  // All prototype maps (will include optimal A)
  //
  
  for (int ki=0; ki<ks.size(); ki++)
    {
      const int K = ks[ki];

      writer.level( K , "KN" );
      
      for (int i=0; i<C; i++)
	{
	  writer.level( signals.label(i) , globals::signal_strat );
	  
	  for (int j=0; j<K; j++)
	    {
	      std::string s="?";	      
	      s[0] = ms_prototypes_t::ms_labels[j];
	      writer.level( s , "K" );
	      writer.value( "A" , results.kres[K].A(i,j) );
	    }
	  writer.unlevel( "K" );
	  
	}
      writer.unlevel( globals::signal_strat );
    }
  writer.unlevel( "KN" );

  
  //
  // Detailed fit outputs (over all K considered --> NK)
  //

  for (int ki=0; ki<ks.size(); ki++)
    {
      const int K = ks[ki];      
      writer.level( K , "NK" );
      writer.value( "R2" , results.kres[K].R2 );
      writer.value( "MSE" , results.kres[K].MSE );
      writer.value( "SIG2" , results.kres[K].sig2 );
      writer.value( "SIG2_MCV" , results.kres[K].sig2_modk_mcv );
    }
  writer.unlevel( "NK" );


  //
  // Save prototypes
  //

  return prototypes;
  
}


ms_backfit_t microstates_t::backfit( const Data::Matrix<double> & X_ ,
				     const Data::Matrix<double> & A_ ,
				     const double lambda , 
				     bool return_GMD )
{
  
  Data::Matrix<double> X = X_ ;
  Data::Matrix<double> A = A_ ; 
  
  // X will be C x N
  // A will be C x K
  
  // polarity invariant back-fitting

  const int C = A.dim1();
  const int K = A.dim2();
  const int N = X.dim2(); // assumes X is already transposed as C x N 

  //
  // GMD: global map dissimilarity
  //
  
  // Assumes average reference already set  
  // Normalise EEG and maps (average reference and gfp = 1 for EEG)
  
  //
  // Global field power 
  //

  Data::Vector<double> GFP( N );
  Data::Vector<double> avg( N );
  for (int j=0; j<N; j++)
    {
      // get time-points across channels
      const Data::Vector<double> & p = X.col( j );
      GFP[j] = sqrt( Statistics::variance( p , 0 ) ); // use N denom
      avg[j] = Statistics::mean( p );
    }

  Data::Vector<double> GFP_A( K );
  Data::Vector<double> avg_A( N );
  for (int j=0; j<K; j++)
    {
      // get time-points across channels
      const Data::Vector<double> & p = A.col( j );
      GFP_A[j] = sqrt( Statistics::variance( p , 0 ) ); // use N 
      avg_A[j] = Statistics::mean( p );
    }

  //
  // Normalize each
  //

  for (int i=0; i<C; i++)
    for (int j=0; j<N; j++)
      X(i,j) = ( X(i,j) - avg[j] ) / GFP[j];
  
  for (int i=0; i<C; i++)
    for (int j=0; j<K; j++)
      A(i,j) = ( A(i,j) - avg_A[j] ) / GFP_A[j];
  
  
  // Global map dissilarity

  Data::Matrix<double> GMD( K , N );

  for (int k=0; k<K; k++)
    {
      Data::Matrix<double> XX = X;
      // for each time-point
      for (int j=0;j<N;j++)
	{
	  double t = 0 , t2 = 0;
	  for (int i=0;i<C;i++)
	    {
	      t += ( X(i,j) - A(i,k) ) * ( X(i,j) - A(i,k) );
	      t2 += ( X(i,j) + A(i,k) ) * ( X(i,j) + A(i,k) ) ;
	    }
	  t = sqrt( t / (double)C );
	  t2 = sqrt( t2 / (double)C );

	  // pick smallest distance to ensure polarity invariance
	  GMD(k,j) = t < t2 ? t : t2 ;
	}
    }

  //
  // Smooth GMDs? 
  //
  
  if ( lambda > 0 )     
    {
      logger << "  applying total-variation denoiser on GMDs, lambda = " << lambda << "\n";
      for (int k=0; k<K; k++)
	{
	  std::vector<double> row = *GMD.row(k).data_pointer();
	  dsptools::TV1D_denoise( row , lambda );
	  for (int i=0;i<N;i++) GMD(k,i) = row[i];
	}
      
    }


  
  //
  // Get matching labels, i.e. assign states (min. GMD)
  //
  
  ms_backfit_t bf(N);
  
  for (int j=0;j<N;j++)
    {
      // add all labels/GMDs which will be sorted by add()
      for (int k=0;k<K;k++)
	{

	  // default: add label 'as is' 
	  if ( ! ms_prototypes_t::ul_groups )
	    {
	      bf.labels[j].add( k , GMD(k,j) );
	    }
	  else
	    {
	      // special case: in u/l grouping, merge 'A' and 'a' to both 'A'
	      // the labels[] ms_assignments_t struct can handle having >1
	      // version of the same label (i.e. will pick the
	      
	      bf.labels[j].add( ms_prototypes_t::ul_reduction(k) , k, GMD(k,j) );
	      
	    }
	}
      
      // pick the best for each time point
      bf.labels[j].set_picks();
      
    }
  
  //
  // Optionally, store GMD for smoothing
  //

  if ( return_GMD )
    bf.GMD = GMD;
  
  return bf;
    
}
  
			     



ms_backfit_t microstates_t::smooth_reject( const ms_backfit_t & sol , 
					   int minTime )
{
  
  const int N = sol.labels.size();

  if ( N == 0 )
    Helper::halt( "solution not populated in smooth_reject()" );

  // make a working copy, which will be editted;
  // no need to populate GMD (esp. as labels may change in any case)
  
  ms_backfit_t bf(N);
  bf.labels = sol.labels;
  bf.ambiguous = sol.ambiguous;
  
  for (int k=1; k <= minTime; k++)
    {
      // track changes
      std::vector<int> cruns( N , k );      
      int iter_num = 0;
      while ( 1 && iter_num < 1000 ) 
	{

	  // should not happen, but set upper limit (1000 iterations) 
	  // for any edge cases...

	  ++iter_num;

	  //
	  // count remaining short unambiguous segments
	  //

	  int sum_cruns = 0;
	  for (int c=0; c<cruns.size(); c++)
	    if ( (!bf.ambiguous[c]) && cruns[c] <= k ) ++sum_cruns;
	  if ( sum_cruns == 0 ) break;
	  
	  //
	  // perform run-length encoding (RLE)
	  //

	  ms_rle_t runs = rle( bf.best() );

	  //
	  // for any short segment that is flanked by ambiguous segments (or start/stop + an ambig segment)
	  // this will never get expanded, and so set this to be amiguous too.
	  //
	  
	  int cnt = 0;
	  for (int r=0; r<runs.c.size(); r++)
	    {
	      const bool start = r == 0 ;
	      const bool last  = r == runs.c.size() - 1 ;
	      const bool unsalvageable = ( start || runs.d[r-1] == -1 ) && ( last || runs.d[r+1] == -1 ) ; 

	      const bool too_short = runs.c[r] <= k;
	      const bool ambig = runs.d[r] == -1 ;

	      for (int j=0; j<runs.c[r]; j++)
		{		
		  		  
		  // shift if any (non-ambiguous) segment is too short
		  
		  if ( too_short && ! ambig ) 
		    {
		      // check that this segment is salvageable 
		      // i.e. not already flanked by ambiguous segments
		      //  or the start/stop		      
		      
		      if ( ! unsalvageable ) 
			bf.labels[cnt].shift();
		      else
			bf.ambiguous[cnt] = true;
		    }

		  cruns[cnt] = runs.c[r];
		  ++cnt;
		}
	    }
	}
    }
    
  return bf;
}


ms_backfit_t microstates_t::smooth_windowed( const ms_backfit_t & labels ,
					     const Eigen::MatrixXd & X_ ,
					     const Eigen::MatrixXd & A_ ,
					     int smooth_width ,	
					     double smooth_weight ,
					     int max_iterations ,
					     double threshold )
{

  Helper::halt( "microstates_t::smooth_windowed() not yet implemented" );
  
  Data::Matrix<double> X = eig2mat( X_ );  // C x N  EEG 
  Data::Matrix<double> A = eig2mat( A_ );  // C x K  prototypes
  const int C = X.dim1();
  const int N = X.dim2();
  const int K = A.dim2();
  ms_backfit_t bf(N);
  return bf;
  
}


ms_rle_t microstates_t::rle( const std::vector<int> & x )
{

  ms_rle_t ret;

  int ind = 0;

  ret.d.push_back( x[0] );
  ret.c.push_back( 1 );

  const int n = x.size();
  
  for (int i=1; i<n; i++)
    {
      if ( x[i-1] == x[i] )
	++ret.c[ind];
      else
	{
	  ++ind;
	  ret.d.push_back( x[i] );
	  ret.c.push_back( 1 );
	}
    }

  return ret;
}


ms_stats_t microstates_t::stats( const Data::Matrix<double> & X_ ,
				 const Data::Matrix<double> & A_ ,
				 const std::vector<int> & L , 
				 const std::vector<int> & L2 )
{

  // assignments are in 'L'
  // under default run, L2 == L
  // under U/l mapping, then L2 reflects the 'original' (i.e. A or a)
  // whereas L reflects the reduced assignment (e.g 'A' irrespective or whether
  // 'A' or 'a' was best fit)

  // when calculating the mean spatial correlation, etc, for a reduced state, we
  // need to look up the original/actual mapping (i.e. best) - i.e. anything
  // that involves the lookup into the GMD() table

  // but when doing the RLE, or looking at coverage stats, we should use
  // the standard/reduced L
  
  // note that K will still be the full set (e.g. 8) if using a reduced map
  //   A B C D a b c d 
  // but only the first ones will be populated w/ >0 coverage

  // note: is okay if one group does not have a given state:
  //    A B C D a b c e
  // i.e. this will imply the following reduced groups (K=5 effectively) 
  //    A B C D E
  
  
  ms_stats_t stats;
  
  Data::Matrix<double> X = X_;
  Data::Matrix<double> A = A_;
  
  const int C = X.dim1();
  const int N = X.dim2();
  const int K = A.dim2();
  
  //
  // Normalize X and A (by mean / set GFP = 1 )
  // (same code as backfit()
  //

  Data::Vector<double> GFP( N );
  Data::Vector<double> GFP_minus1( N ); // also get w/ N-1 denom for comparability 
  Data::Vector<double> avg( N );
  for (int j=0; j<N; j++)
    {
      // get time-points across channels
      const Data::Vector<double> & p = X.col( j );
      GFP[j] = sqrt( Statistics::variance( p , 0 ) ); // use N denomx
      GFP_minus1[j] = sqrt( Statistics::variance( p , 1 ) ); // use N-1 denomx (to get same output as Matlab implementation)
      avg[j] = Statistics::mean( p );
    }

  Data::Vector<double> GFP_A( K );
  Data::Vector<double> avg_A( N );
  for (int j=0; j<K; j++)
    {
      const Data::Vector<double> & p = A.col( j );
      GFP_A[j] = sqrt( Statistics::variance( p , 0 ) ); // use N denom
      avg_A[j] = Statistics::mean( p );
    }

  for (int i=0; i<C; i++)
    for (int j=0; j<N; j++)
      X(i,j) = ( X(i,j) - avg[j] ) / GFP[j];
  
  for (int i=0; i<C; i++)
    for (int j=0; j<K; j++)
      A(i,j) = ( A(i,j) - avg_A[j] ) / GFP_A[j];


  //
  // GMD Global map dissilarity  (K x N) 
  //
  
  Data::Matrix<double> GMD( K , N );

  for (int k=0; k<K; k++)
    {
      Data::Matrix<double> XX = X;
      // for each time-point
      for (int j=0;j<N;j++)
	{
	  double t = 0 , t2 = 0;
	  for (int i=0;i<C;i++)
	    {
	      t += ( X(i,j) - A(i,k) ) * ( X(i,j) - A(i,k) );
	      t2 += ( X(i,j) + A(i,k) ) * ( X(i,j) + A(i,k) ) ;
	    }
	  t = sqrt( t / (double)C );
	  t2 = sqrt( t2 / (double)C );

	  // pick smallest distance to ensure polarity invariance
	  GMD(k,j) = t < t2 ? t : t2 ;
	}
    }

  // Spatial correlations
  
  Data::Matrix<double> SpatCorr( K , N );
  for (int i=0;i<K;i++)
    for (int j=0;j<N;j++)
      SpatCorr(i,j) = 1 - ( GMD(i,j) * GMD(i,j) ) / 2.0;
  
  // Total GEV (nb. based on un-normalized X_)..
  
  Data::Vector<double> var = Statistics::sdev( X_ ,Statistics::mean(  X_ ) );
  double denom = 0;
  for (int j=0;j<N;j++)
    {
      var[j] *= var[j];
      denom += var[j];
    }

  // nb. uses L2 for original assigment 
  //   (only has effect with U/l mapping)
  stats.GEV_tot = 0;
  for (int j=0;j<N;j++)
    if ( L2[j] != -1 ) // for unambiguous assignments...
      stats.GEV_tot += SpatCorr( L2[j] , j ) * var[j] ;
  stats.GEV_tot /= denom;
  
  // transition probabilities: 
  //   (U/l mapping note: here use L not L2)
  
  ms_rle_t runs = rle( L );
  

  // copy GFP (per point -- needed?)
  
  stats.GFP = GFP;

  // Get spatial correlations
  
  stats.SpatCorr = SpatCorr;
  
  // means: note: these are for "ALL" maps
  stats.m_gfp.resize(K);
  stats.m_dur.resize(K);
  stats.m_occ.resize(K);
  stats.m_occ_unambig.resize(K);
  stats.m_cov.resize(K);
  stats.m_cov_unambig.resize(K);
  stats.m_wcov.resize(K);
  stats.m_gev.resize(K);
  stats.m_spc.resize(K);

  // track ambiguous points
  //  U/l mapping note: does not matter if we use L or L2 here

  int ambig = 0;
  for (int j=0;j<N;j++)
    if ( L[j] == -1 ) ++ambig;
  
  writer.value( "AMBIG" , ambig / (double)N );


  //
  // stats for each class
  //

  
  for (int k=0; k<K; k++)
    {
      
      // Mean GFP
      // U/l note:  here use L[] 
      // nb. uses /n-1 denom. here only to match ML Toolbox 
      
      std::vector<double> gfp_k;
      for (int j=0;j<N;j++)
	if ( L[j] == k )
	  gfp_k.push_back( GFP_minus1[j] );
      stats.m_gfp[k] = MiscMath::mean( gfp_k );
      
      // occur/duration
      // U/l note: uses L[] implicitly (rle_t above)
      std::vector<double> times;
      for (int i=0; i<runs.d.size(); i++)
	if ( runs.d[i] == k ) times.push_back( runs.c[i] * ( 1000.0 / sr ) );
      
      stats.m_occ[k] = times.size() / (double) N * sr;
      stats.m_dur[k] = MiscMath::mean( times );
      stats.m_cov[k] = ( stats.m_occ[k] * stats.m_dur[k] ) / 1000.0;
      
      // versions w/ denominator as only unambiguously assigned points
      stats.m_occ_unambig[k] = stats.m_occ[k] / ( 1.0 - ambig / (double)N );
      stats.m_cov_unambig[k] = stats.m_cov[k] / ( 1.0 - ambig / (double)N );
      
      // weighted coverage (mean spatial correlation across /all/ points)
      // U/l note: implicitly, this is being calculated for all K=4+4=8 maps
      //  i.e. is map specific, can be easily averaged / summed in the output
      stats.m_wcov[k] = MiscMath::mean( *SpatCorr.row(k).data_pointer() );

      
      // mean spatial correl (of assigned points) 
      // U/l note: this uses L2[] mapping to get the value
      //  but pools the output for the reduced state 
      std::vector<double> spc_k;
      for (int j=0;j<N;j++)
        if ( L[j] == k ) // matches on reduced form L[] 
          spc_k.push_back( SpatCorr( L2[j], j ) ); // but uses actual/original L2[]
      stats.m_spc[k] = MiscMath::mean( spc_k );
      
      // GEV
      // U/l mapping note: as above
      double numer = 0;
      double denom = 0;
      for (int j=0;j<N;j++)
	{
	  if ( L[j] == k ) // nb.  uses L[] for putput, but L2[] for lookup
	    numer += ( SpatCorr( L2[j], j ) * GFP[j] ) * ( SpatCorr( L2[j], j ) * GFP[j] ) ;
	  denom += GFP[j] * GFP[j];
	}
      
      stats.m_gev[k] = numer / denom;
      
      // next K
    }
  

  //
  // transition probs: these all based on reduced L mapping
  //
  
  // runs.d contains sequence of states
  
  const int seqlen = runs.d.size();
  
  // count ambig (-1) as the class 'K+1' here 

  stats.tr.resize(K+1,K+1);
  Data::Vector<double> row(K+1);
  
  for (int s = 0 ; s < seqlen - 1 ; s++)
    {
      const int klabel1 = runs.d[s] != -1 ? runs.d[s] : K ;
      const int klabel2 = runs.d[s+1] != -1 ? runs.d[s+1] : K ;
      ++stats.tr( klabel1 , klabel2 );
      ++row( klabel1 );
    }

  for (int i=0;i<K+1;i++)
    for (int j=0;j<K+1;j++)
      if ( i != j && row[i] > 0 ) 
	stats.tr(i,j) /= row[i];

  //
  // LZW complexity & SE
  //
  
  lzw_t lzw( runs.d , &stats.lwz_states );

  // mse_t se;
  // for (int m=1; m<=8; m++)
  //   stats.samplen[m] = se.sampen( runs.d , m );
  
  //
  // dump sequences to file?
  //


  if ( statesfile != ""  )
    {      
      logger << "  writing sequence order to " << statesfile << "\n";

      // only write unambiguous state sequences
      // as we are splicing out ambiguous segments, we need to manually check 
      // that we do not have a duplicate sequence here...
      // U/l mapping: this is implicitly based on L[], not L2[]
      std::ofstream OUT1( Helper::expand( statesfile ).c_str() , std::ios::out );
      OUT1 << subj_id << "\t";

      char last = '?';
      for (int i=0; i<runs.d.size(); i++) 
	if ( runs.d[i] != -1 )
	  {
	    char curr = ms_prototypes_t::ms_labels[ runs.d[i] ];
	    if ( curr != last ) 
	      {
		OUT1 << curr;
		last = curr; 
	      }
	  }
      OUT1 << "\n";
      OUT1.close();      
    }

  
  //
  // k-mer distributions, optionally (in single-obs mode)
  //
  
  if ( kmers_nreps )    
    stats.kmers.run( runs.d , kmers_min , kmers_max , kmers_nreps , kmers_w );
     
  return stats;
}



void ms_kmer_t::run( const std::map<std::string,std::vector<int> > & lall , int k1 , int k2 , int nreps , int w, 
		     const std::map<std::string,int> * grp , bool verbose )
{
  std::map<std::string,std::string> sall;
  std::map<std::string,std::vector<int> >::const_iterator ii = lall.begin();
  while ( ii != lall.end() )
    {
      const std::vector<int> & l = ii->second;
      std::string & s = sall[ ii->first ];
      const int n = l.size();
      
      std::vector<char> newseq;
      char last = '?';
      for (int i=0; i<n; i++)
        if ( l[i]  != -1 ) 
	  {
            char curr = ms_prototypes_t::ms_labels[ l[i] ];
            if ( curr != last )
              {
                newseq.push_back( curr );
                last = curr;
              }
	  }
            
      // do not allow for ambiguous labels here
      // as we are splicing out '?' codes, we need to manually make sure we don't 
      // have additional repeats  A?A --> A  not AA

      const int n2 = newseq.size();
      s = std::string( n2 , '?' );
      for (int i=0; i<n2; i++)
	s[i] = newseq[i];
      
      ++ii;
    }
  run( sall , k1 , k2 , nreps , w, grp , verbose );
}

void ms_kmer_t::run( const std::map<std::string,std::string> & sall , int k1 , int k2 , int nreps , int w, 
		     const std::map<std::string,int> * grp , bool verbose )
{

  // single obs mode?
  bool single_obs = sall.size() == 1;

  // group/phenotype contrast?
  // must be coded 0 / 1  (any other value == missing ) 
  bool grp_contrast = grp != NULL;  
  
  //
  // ensure bounded if using this naive algorithm 
  //
  
  if ( k1 < 2 ) k1 = 2;
  if ( k2 > 8 ) k2 = 8;


  //
  // Observed data: basic counts (pooling across all individuals)
  //

  std::map<std::string,std::string>::const_iterator ii = sall.begin();
  while ( ii != sall.end() )
    {
      // sequences for this individual
      const std::string & s = ii->second;
      const int n = s.size();

      // phenotype? (expecting 0/1)
      int phe = -1;
      if ( grp != NULL && grp->find( ii->first ) != grp->end() ) 
	phe = grp->find( ii->first )->second;
      
      // count observed sequences (pools across all individuals)
      for (int i=0; i<n; i++)
	for (int j=k1; j<=k2; j++)
	  if ( i+j < n )
	    {
	      const std::string ss = s.substr( i , j ) ;
	      ++basic.obs[ ss ];
	      if      ( phe == 0 ) ++basic_controls.obs[ ss ];
	      else if ( phe == 1 ) ++basic_cases.obs[ ss ];
	    }
      
      // next individual
      ++ii;
    }

  

  //
  // Observed data: C/C diffs in the overall counts
  //
  
  if ( grp != NULL )
    {
      std::map<std::string,double>::const_iterator cc = basic.obs.begin();
      while ( cc != basic.obs.end() )
	{
	  basic_diffs.obs[ cc->first ] = basic_cases.obs[ cc->first ] - basic_controls.obs[ cc->first ];
	  ++cc;
	}
    }
  

  //
  // Observed data: form equivalance groups, and get both sums (group) and within-group diffs (equiv)
  //

  // for each unique sequence, find the equivalance group:
  //  - the other sequences with the same numbers of each state, and also with no
  //    matching contiguous states (i.e. like the original sequences)

  // label each equivalance group by the first sorted value then
  // permute data to calculate the sum counts (i.e. to get a
  // group-normalized relative frequency )

  // observed -->  equiv group key.  e.g.  BCA --> ABC
  obs2equiv.clear();

  // equiv group key --> all members e.g.  ABC -->  ABC, ACB, BAC, BCA, CAB, CBA
  equivs.clear();

  // for each observed sequence, find the set of equivalent sequences  
  std::map<std::string,double>::const_iterator cc = basic.obs.begin();
  while ( cc != basic.obs.end() )
    {
      std::string str = cc->first;
      std::string key = first_permute( str );
      obs2equiv[ str ] = key;
      ++cc;
    }

  std::map<std::string,std::string>::const_iterator ee = obs2equiv.begin();
  while ( ee != obs2equiv.end() )
    {
      // based just on the keys, get the corresponding equiv. group
      equivs[ ee->second ] = permute( ee->second );

      // for convenienc of output, track size of equiv set directly
      equiv_set_size[ ee->first ] = equivs[ ee->second ].size();
     
      ++ee;
    }


  //
  // Group statistics (by e-label)
  //
  
  std::map<std::string,std::set<std::string> >::const_iterator ll = equivs.begin();
  while ( ll != equivs.end() )
    {
      // 'group' = denominator for : obs / obs-group 
      const std::set<std::string> & st = ll->second;
      std::set<std::string>::const_iterator qq = st.begin();
      while ( qq != st.end() )
	{	  
	  group.obs[ ll->first ] += basic.obs[ *qq ];
	  if ( grp != NULL )
	    {
	      group_cases.obs[ ll->first ] += basic_cases.obs[ *qq ];
	      group_controls.obs[ ll->first ] += basic_controls.obs[ *qq ];	      
	    }
	  ++qq;	  
	}      

      // C/C difference stat
      if ( grp != NULL )
	group_diffs.obs[ ll->first ] = group_cases.obs[ ll->first ] - group_controls.obs[ ll->first ];
      
      ++ll;
    }


  //
  // count up for OBS - GROUP (within group diffs)
  //

  ee = obs2equiv.begin();
  while ( ee != obs2equiv.end() )
    {
     
      // within-group enrichment: express as proportion of total (group.obs is set to sum of equiv groups above)      
      const std::string ess = obs2equiv[ ee->first ] ;

      equiv.obs[ ee->first ] = basic.obs[ ee->first ] / group.obs[ ess ];
      
      if ( grp != NULL )
	{
	  equiv_cases.obs[ ee->first ] = basic_cases.obs[ ee->first ] / group_cases.obs[ ess ];
	  equiv_controls.obs[ ee->first ] = basic_controls.obs[ ee->first ] / group_controls.obs[ ess ];	  
	  // diffs based on difference in these two relative freqs	  
	  equiv_diffs.obs[ ee->first ] = equiv_cases.obs[ ee->first ] - equiv_controls.obs[ ee->first ];
	}
	      
      ++ee;
    }

  
  if ( verbose ) 
    {
      logger << "  kmers: considering length " << k1 << " to " << k2 << "\n";
      
      logger << "  kmers: for " << basic.obs.size() << " sequences, "
	     << equivs.size() << " equivalence groups\n";
      
      logger << "  kmers: running " << nreps << " replicates, ";

      if ( w == 0 ) logger << "using global picks\n";
      else logger << "using local picks (w=" <<w << ")\n";

    }

  
  //
  // Begin permutations 
  //

  for (int r = 0 ; r < nreps ; r++)
    {
      
      if ( verbose ) 
	{
	  if ( r == 0 ) logger << "  ";
	  logger << ".";
	  if ( r % 10 == 9 ) logger << ". " << r+1 << " of " << nreps << " replicates done\n  ";
	}


      //
      // Permutation stratified by individual, although all counting is
      // pooled across individuals
      //
      
      // rather than standard CRandom::random_draw(), this ensures
      // no (*) similar states are contiguous
      // ... actually, there may be one or two contiguous sequences
      // at the end, but in the big picture, this should not matter 

      std::map<std::string,int> stat_basic, stat_basic_cases, stat_basic_controls;
      std::map<std::string,int> stat_group, stat_group_cases, stat_group_controls;
      std::map<std::string,double> stat_equiv, stat_equiv_cases, stat_equiv_controls;
      
      //
      // Iterate over individuals, accumulating statistics (i.e. summed over individuals)
      //

      std::map<std::string,std::string>::const_iterator ii = sall.begin();
      while ( ii != sall.end() )
	{      

	  //
	  // Get the optional phenotype
	  //
	  
	  int phe = -1;
	  if ( grp != NULL && grp->find( ii->first ) != grp->end() )
	    phe = grp->find( ii->first )->second;

	  //
	  // Get a permuted sequence (for this individual only)
	  //
	  
	  std::string ps = modified_random_draw( ii->second , w );	  

	  const int n = ps.size();
	  
	  //
	  // basic kmer count for this individual/replicate
	  //
	  
	  std::map<std::string,int> perm; // temporary per counts for this individual
	  
	  for (int i=0; i<n; i++)
	    for (int j=k1; j<=k2; j++)
	      if ( i+j < n )
		{
		  const std::string ss = ps.substr( i , j );

		  ++perm[ ss ]; // temporary (for equiv stats)
		  ++stat_basic[ ss ]; // main aggregators
		  if      ( phe == 0 ) ++stat_basic_controls[ ss ];
		  else if ( phe == 1 ) ++stat_basic_cases[ ss ];
		}
	  
	  
	  //
	  // next individual
	  //
	  
	  ++ii;
	  
	}
      
      
      //
      // Group sum statistics under the null
      //
      
      std::map<std::string,std::set<std::string> >::const_iterator pp = equivs.begin();
      while ( pp != equivs.end() )
	{
	  int sum = 0 , sum_cases = 0 , sum_controls = 0;
	  const std::set<std::string>  & eqs = pp->second;
	  std::set<std::string>::const_iterator ee = eqs.begin();
	  while ( ee != eqs.end() )
	    {	      
	      sum += stat_basic[ *ee ];
	      if ( grp != NULL )
		{
		  sum_cases += stat_basic_cases[ *ee ];
		  sum_controls += stat_basic_controls[ *ee ]; 
		}
	      ++ee;
	    }
	  
	  // sum across individuals; e-groups are labelled by pp->first, the group key
	  stat_group[ pp->first ] = sum;
	  if ( grp != NULL )
	    {
	      stat_group_controls[ pp->first ] = sum_controls;
	      stat_group_cases[ pp->first ] = sum_cases;
	    }
	  
	  ++pp;
	}

      
      //
      // Track permuted statistics for this replicate
      //

      std::map<std::string,double>::const_iterator ss = basic.obs.begin();
      while ( ss != basic.obs.end() )
	{
	  
	  basic.perm[ ss->first ].push_back( stat_basic[ ss->first ] );

	  if ( grp != NULL )
	    {
	      basic_cases.perm[ ss->first ].push_back( stat_basic_cases[ ss->first ] );
	      basic_controls.perm[ ss->first ].push_back( stat_basic_controls[ ss->first ] );
	      basic_diffs.perm[ ss->first ].push_back( stat_basic_cases[ ss->first ] - stat_basic_controls[ ss->first ] );
	    }

	  //
	  // Group relative diffs ( equiv = base / group )
	  //
	  
	  const std::string & ess = obs2equiv[ ss->first ];
	  
	  double estat = stat_group[ ess ] > 0 ? stat_basic[ ss->first ]  / (double)stat_group[ ess ] : 0 ;
	  equiv.perm[ ss->first ].push_back( estat );
	  
	  if ( grp != NULL )
	    {
	      
	      double estat_cases = stat_group_cases[ ess ] > 0 ? stat_basic_cases[ ss->first ]  / (double)stat_group_cases[ ess ] : 0 ;
	      double estat_controls = stat_group_controls[ ess ] > 0 ? stat_basic_controls[ ss->first ]  / (double)stat_group_controls[ ess ] : 0 ;
	      double estat_diffs  = estat_cases - estat_controls;

	      equiv_cases.perm[ ss->first ].push_back( estat_cases );
	      equiv_controls.perm[ ss->first ].push_back( estat_controls );
	      equiv_diffs.perm[ ss->first ].push_back( estat_diffs );

	    }

	  // next obs sequence
	  ++ss;
	}


      //
      // Group sums (keyed on equiv-group label, not the obs label) 
      //
      
      std::map<std::string,std::set<std::string> >::const_iterator ee = equivs.begin();
      while ( ee != equivs.end() )
        {
	  
	  group.perm[ ee->first ].push_back( stat_group[ ee->first ] );
	  
	  if ( grp != NULL )
	    {
	      group_cases.perm[ ee->first ].push_back( stat_group_cases[ ee->first ] );
	      group_controls.perm[ ee->first ].push_back( stat_group_controls[ ee->first ] );
	      group_diffs.perm[ ee->first ].push_back( stat_group_cases[ ee->first ] - stat_group_controls[ ee->first ] );
	    }
	  
	  // next obs equiv-group 
	  ++ee;
	}

      
      //
      // Next replicate
      //
	  
    } 

  
  //
  // All replicates complete: we have obs and perm[] from which we can get all stats to report;
  // do separately for each sequence
  //
	  
  std::map<std::string,double>::const_iterator oo = basic.obs.begin();
  while ( oo != basic.obs.end() )
    {
      
      const std::string & ss = oo->first ;

      const std::vector<double> & pp_basic = basic.perm[ ss ] ;
      const std::vector<double> & pp_basic_cases = basic_cases.perm[ ss ] ;
      const std::vector<double> & pp_basic_controls = basic_controls.perm[ ss ] ;
      const std::vector<double> & pp_basic_diffs = basic_diffs.perm[ ss ] ;

      // expected value
      basic.exp[ ss ] = MiscMath::mean( pp_basic );
      
      // Z score
      basic.zscr[ ss ] = ( basic.obs[ ss ] - basic.exp[ ss ] ) / MiscMath::sdev( pp_basic , basic.exp[ ss ]  );
      
      // Empirical P
      int pv = 0;
      for (int r=0; r<nreps; r++) if ( pp_basic[r] >= basic.obs[ ss ] ) ++pv;
      basic.pval[ ss ] = ( 1 + pv ) / (double)( 1 + nreps );

      // phenotype-based
      if ( grp != NULL )
	{
	  basic_cases.exp[ ss ] = MiscMath::mean( pp_basic_cases );
	  basic_controls.exp[ ss ] = MiscMath::mean( pp_basic_controls );
	  basic_diffs.exp[ ss ] = MiscMath::mean( pp_basic_diffs );      

	  basic_cases.zscr[ ss ] = ( basic_cases.obs[ ss ] - basic_cases.exp[ ss ] ) / MiscMath::sdev( pp_basic_cases , basic_cases.exp[ ss ]  );
	  basic_controls.zscr[ ss ] = ( basic_controls.obs[ ss ] - basic_controls.exp[ ss ] ) / MiscMath::sdev( pp_basic_controls , basic_controls.exp[ ss ]  );
	  basic_diffs.zscr[ ss ] = ( basic_diffs.obs[ ss ] - basic_diffs.exp[ ss ] ) / MiscMath::sdev( pp_basic_diffs , basic_diffs.exp[ ss ]  );

	  // C/C only
	  int pv = 0;
	  for (int r=0; r<nreps; r++) if ( pp_basic_cases[r] >= basic_cases.obs[ ss ] ) ++pv;
	  basic_cases.pval[ ss ] = ( 1 + pv ) / (double)( 1 + nreps );
	  
	  pv = 0;
	  for (int r=0; r<nreps; r++) if ( pp_basic_controls[r] >= basic_controls.obs[ ss ] ) ++pv;
	  basic_controls.pval[ ss ] = ( 1 + pv ) / (double)( 1 + nreps );
	}

      //
      // Same, for equiv-group stats
      //

      const std::vector<double> & pp_equiv = equiv.perm[ ss ] ;
      const std::vector<double> & pp_equiv_cases = equiv_cases.perm[ ss ] ;
      const std::vector<double> & pp_equiv_controls = equiv_controls.perm[ ss ] ;
      const std::vector<double> & pp_equiv_diffs = equiv_diffs.perm[ ss ] ;

      // expected value
      equiv.exp[ ss ] = MiscMath::mean( pp_equiv );
      
      // Z score
      equiv.zscr[ ss ] = ( equiv.obs[ ss ] - equiv.exp[ ss ] ) / MiscMath::sdev( pp_equiv , equiv.exp[ ss ]  );
      
      // Empirical P
      pv = 0;
      for (int r=0; r<nreps; r++) if ( pp_equiv[r] >= equiv.obs[ ss ] ) ++pv;
      equiv.pval[ ss ] = ( 1 + pv ) / (double)( 1 + nreps );

      // phenotype-based
      if ( grp != NULL )
	{
	  equiv_cases.exp[ ss ] = MiscMath::mean( pp_equiv_cases );
	  equiv_controls.exp[ ss ] = MiscMath::mean( pp_equiv_controls );
	  equiv_diffs.exp[ ss ] = MiscMath::mean( pp_equiv_diffs );      

	  equiv_cases.zscr[ ss ] = ( equiv_cases.obs[ ss ] - equiv_cases.exp[ ss ] ) / MiscMath::sdev( pp_equiv_cases , equiv_cases.exp[ ss ]  );
	  equiv_controls.zscr[ ss ] = ( equiv_controls.obs[ ss ] - equiv_controls.exp[ ss ] ) / MiscMath::sdev( pp_equiv_controls , equiv_controls.exp[ ss ]  );
	  equiv_diffs.zscr[ ss ] = ( equiv_diffs.obs[ ss ] - equiv_diffs.exp[ ss ] ) / MiscMath::sdev( pp_equiv_diffs , equiv_diffs.exp[ ss ]  );

	  // C/C only
	  int pv = 0;
	  for (int r=0; r<nreps; r++) if ( pp_equiv_cases[r] >= equiv_cases.obs[ ss ] ) ++pv;
	  equiv_cases.pval[ ss ] = ( 1 + pv ) / (double)( 1 + nreps );
	  
	  pv = 0;
	  for (int r=0; r<nreps; r++) if ( pp_equiv_controls[r] >= equiv_controls.obs[ ss ] ) ++pv;
	  equiv_controls.pval[ ss ] = ( 1 + pv ) / (double)( 1 + nreps );
	}

      
      //
      // next sequence
      //
      
      ++oo;
    }



  //
  // Finally, report for GROUP statistics
  //

  oo = group.obs.begin();
  while ( oo != group.obs.end() )
    {
      
      const std::string & ss = oo->first ;

      const std::vector<double> & pp_group = group.perm[ ss ] ;
      const std::vector<double> & pp_group_cases = group_cases.perm[ ss ] ;
      const std::vector<double> & pp_group_controls = group_controls.perm[ ss ] ;
      const std::vector<double> & pp_group_diffs = group_diffs.perm[ ss ] ;
      
      // expected value
      group.exp[ ss ] = MiscMath::mean( pp_group );
      
      // Z score
      group.zscr[ ss ] = ( group.obs[ ss ] - group.exp[ ss ] ) / MiscMath::sdev( pp_group , group.exp[ ss ]  );
      
      // Empirical P
      int pv = 0;
      for (int r=0; r<nreps; r++) if ( pp_group[r] >= group.obs[ ss ] ) ++pv;
      group.pval[ ss ] = ( 1 + pv ) / (double)( 1 + nreps );

      // phenotype-based
      if ( grp != NULL )
	{
	  group_cases.exp[ ss ] = MiscMath::mean( pp_group_cases );
	  group_controls.exp[ ss ] = MiscMath::mean( pp_group_controls );
	  group_diffs.exp[ ss ] = MiscMath::mean( pp_group_diffs );      

	  group_cases.zscr[ ss ] = ( group_cases.obs[ ss ] - group_cases.exp[ ss ] ) / MiscMath::sdev( pp_group_cases , group_cases.exp[ ss ]  );
	  group_controls.zscr[ ss ] = ( group_controls.obs[ ss ] - group_controls.exp[ ss ] ) / MiscMath::sdev( pp_group_controls , group_controls.exp[ ss ]  );
	  group_diffs.zscr[ ss ] = ( group_diffs.obs[ ss ] - group_diffs.exp[ ss ] ) / MiscMath::sdev( pp_group_diffs , group_diffs.exp[ ss ]  );

	  // C/C only
	  int pv = 0;
	  for (int r=0; r<nreps; r++) if ( pp_group_cases[r] >= group_cases.obs[ ss ] ) ++pv;
	  group_cases.pval[ ss ] = ( 1 + pv ) / (double)( 1 + nreps );
	  
	  pv = 0;
	  for (int r=0; r<nreps; r++) if ( pp_group_controls[r] >= group_controls.obs[ ss ] ) ++pv;
	  group_controls.pval[ ss ] = ( 1 + pv ) / (double)( 1 + nreps );
	}

      // next e-group
      ++oo;
    }

  //
  // Scale OBS and EXP rto be frequencies (within L, kmer size)
  //

  std::map<int,double> sum;

  // basic, OBS
  std::map<std::string,double>::iterator nn = basic.obs.begin();
  while ( nn != basic.obs.end() ) { sum[ nn->first.size() ] += nn->second; ++nn; }
  nn = basic.obs.begin();
  while ( nn != basic.obs.end() ) { nn->second /= sum[ nn->first.size() ] ; ++nn; }
  
  // basic, EXP
  sum.clear(); nn = basic.exp.begin();
  while ( nn != basic.exp.end() ) { sum[ nn->first.size() ] += nn->second; ++nn; }
  nn = basic.exp.begin();
  while ( nn != basic.exp.end() ) { nn->second /= sum[ nn->first.size() ] ; ++nn; }

  if ( grp != NULL )
    {
      // basic, cases, OBS
      sum.clear(); nn = basic_cases.obs.begin();
      while ( nn != basic_cases.obs.end() ) { sum[ nn->first.size() ] += nn->second; ++nn; }
      nn = basic_cases.obs.begin();
      while ( nn != basic_cases.obs.end() ) { nn->second /= sum[ nn->first.size() ] ; ++nn; }
      
      // basic, cases, EXP
      sum.clear(); nn = basic_cases.exp.begin();
      while ( nn != basic_cases.exp.end() ) { sum[ nn->first.size() ] += nn->second; ++nn; }
      nn = basic_cases.exp.begin();
      while ( nn != basic_cases.exp.end() ) { nn->second /= sum[ nn->first.size() ] ; ++nn; }

      // basic, controls, OBS
      sum.clear(); nn = basic_controls.obs.begin();
      while ( nn != basic_controls.obs.end() ) { sum[ nn->first.size() ] += nn->second; ++nn; }
      nn = basic_controls.obs.begin();
      while ( nn != basic_controls.obs.end() ) { nn->second /= sum[ nn->first.size() ] ; ++nn; }
      
      // basic, controls, EXP
      sum.clear(); nn = basic_controls.exp.begin();
      while ( nn != basic_controls.exp.end() ) { sum[ nn->first.size() ] += nn->second; ++nn; }
      nn = basic_controls.exp.begin();
      while ( nn != basic_controls.exp.end() ) { nn->second /= sum[ nn->first.size() ] ; ++nn; }
      
    }

  // ----- GROUP STATS -----

  // group, OBS
  sum.clear(); nn = group.obs.begin();
  while ( nn != group.obs.end() ) { sum[ nn->first.size() ] += nn->second; ++nn; }
  nn = group.obs.begin();
  while ( nn != group.obs.end() ) { nn->second /= sum[ nn->first.size() ] ; ++nn; }
  
  // group, EXP
  sum.clear(); nn = group.exp.begin();
  while ( nn != group.exp.end() ) { sum[ nn->first.size() ] += nn->second; ++nn; }
  nn = group.exp.begin();
  while ( nn != group.exp.end() ) { nn->second /= sum[ nn->first.size() ] ; ++nn; }

  if ( grp != NULL )
    {
      // group, cases, OBS
      sum.clear(); nn = group_cases.obs.begin();
      while ( nn != group_cases.obs.end() ) { sum[ nn->first.size() ] += nn->second; ++nn; }
      nn = group_cases.obs.begin();
      while ( nn != group_cases.obs.end() ) { nn->second /= sum[ nn->first.size() ] ; ++nn; }
      
      // group, cases, EXP
      sum.clear(); nn = group_cases.exp.begin();
      while ( nn != group_cases.exp.end() ) { sum[ nn->first.size() ] += nn->second; ++nn; }
      nn = group_cases.exp.begin();
      while ( nn != group_cases.exp.end() ) { nn->second /= sum[ nn->first.size() ] ; ++nn; }

      // group, controls, OBS
      sum.clear(); nn = group_controls.obs.begin();
      while ( nn != group_controls.obs.end() ) { sum[ nn->first.size() ] += nn->second; ++nn; }
      nn = group_controls.obs.begin();
      while ( nn != group_controls.obs.end() ) { nn->second /= sum[ nn->first.size() ] ; ++nn; }
      
      // group, controls, EXP
      sum.clear(); nn = group_controls.exp.begin();
      while ( nn != group_controls.exp.end() ) { sum[ nn->first.size() ] += nn->second; ++nn; }
      nn = group_controls.exp.begin();
      while ( nn != group_controls.exp.end() ) { nn->second /= sum[ nn->first.size() ] ; ++nn; }
      
    }


 
  //
  // All done...
  //
  
  
}



std::string ms_kmer_t::first_permute( std::string str )
{
  const int n = str.size();
  std::sort( str.begin() , str.end() );
  do {
    // do not include permutations with similar contiguous blocks
    // as originals will never have this feature
    bool okay = true;
    for (int i=1;i<n;i++)
      if ( str[i-1] == str[i] )
	{ okay = false; break; }
    if ( okay ) return str;
  }
  while ( next_permutation( str.begin() , str.end() ) );
  // should not happen, i.e. str should always be a possible first key
  // if no similar contiguous states
  Helper::halt( "invalid sequence given" );
  return "";

}

std::set<std::string> ms_kmer_t::permute( std::string str )
{
  std::set<std::string> perms;
  if ( str.size() == 0) return perms;
  const int n = str.size();
  std::sort( str.begin() , str.end() );
  do {
    // do not include permutations with similar contiguous blocks
    // as originals will never have this feature
    bool okay = true;
    for (int i=1;i<n;i++)
      if ( str[i-1] == str[i] )
	{ okay = false; break; }
    if ( okay ) perms.insert( str );    
  }
  while ( next_permutation( str.begin() , str.end() ) );
  return perms;
}
 

  
void ms_prototypes_t::write( const std::string & filename )
{
  logger << "  writing " << K << "-class prototypes to " << filename << "\n";

  // now have a header
  std::ofstream O1( filename.c_str() , std::ios::out );
  O1 << "CH";
  for (int k=0;k<K;k++) O1 << "\t" << ms_labels[k];
  O1 << "\n";
  
  // channels
  for (int c=0; c<C; c++)
    {
      O1 << chs[c];
      for (int k=0;k<K;k++)
	O1 << "\t" << A(c,k);
      O1 << "\n";
    }
  O1.close();
}



void ms_prototypes_t::read( const std::string & f )
{
  
  const std::string & filename = Helper::expand( f );

  if ( ! Helper::fileExists( filename ) )
    Helper::halt( "could not find " + filename );
  
  A.resize(0,0);
  chs.clear();
  C = 0;

  std::vector<double> t;
  
  std::ifstream IN1( filename.c_str() , std::ios::in );
  
  // get header row
  std::string line;
  Helper::safe_getline( IN1 , line );
  if ( line == "" || IN1.eof() ) 
    Helper::halt( "bad format for " + filename );

  std::vector<std::string> tok = Helper::parse( line );
  if ( tok.size() < 3 )
    Helper::halt( "problem reading prototypes from " + filename + "\n fewer than 2 classes\n" + line );
  
  if ( tok[0] != "CH" )
    Helper::halt( "expecting first column to be 'CH' in " + filename );
  

  // set class labels

  K = tok.size() - 1 ;
  
  logger << "  found " << K << " classes:" ;
  ms_labels.resize( K );
  for (int k=0; k<K; k++)
    {
      if ( tok[ k+1 ].size() != 1 ) 
	Helper::halt( "state label cannot be >1 char : " + tok[k+1] );
      ms_labels[ k ] = tok[ k+1 ][ 0 ];
      logger << " " << ms_labels[ k ] ;
    }
  logger << "\n";
  

  while ( ! IN1.eof() )
    {
      std::string line;
      Helper::safe_getline( IN1 , line );
      if ( IN1.eof() ) break;
      if ( line == "" ) break;
      
      std::vector<std::string> tok = Helper::parse( line );

      if ( tok.size() != K + 1 ) 
	Helper::halt( "problem reading prototypes (bad column number) from " + filename );
      
      for (int i=1;i<tok.size();i++)
	{
	  // nb. effective +1 encoding in file, due to channel col 0 
	  double x;
	  
	  if ( ! Helper::str2dbl( tok[i] , &x ) )
	    Helper::halt( "problem reading prototypes from " 
			  + filename + "\n in coversion to numeric: " 
			  + tok[i] + "\n" + line );
	  t.push_back(x);
	}
      
      chs.push_back( tok[0] );
      
      ++C;
    }
  
  IN1.close();
  
  if ( K == 0 || C == 0 )
    Helper::halt( "problem reading prototypes from " + filename + ": K or C == 0" );
  
  if ( t.size() != K * C )
    Helper::halt( "problem reading prototypes from " + filename + ": KC != # data points" );
  
  // store MAP in A
  A.resize( C , K );
  int tc = 0;
  for (int c=0; c<C; c++)
    for (int k=0;k<K;k++)
      A(c,k) = t[tc++];
  
  logger << "  read " << K << "-class prototypes for " << C << " channels from " << filename << "\n";

}


//
// two-group (U/l) reductions?
//

void ms_prototypes_t::make_ul_map()
{
  // to populate case-invariant mapping:
  // static std::vectormap<int,int> ul_mapping;
  
  // pool UPPER and lower variants
  // to the first instance of either A or a
  //  (in output the string will be always
  //   changed to uppercase anyway)
  //  but need to handle case of
  //   A B C D   a b c e
  //  i.e. output is (case/group invariant)
  //    A B C D E  

  // construct targets (first instance of
  //  the class), to make sure this points
  //  to the correct ms_labels[] (ignoring case)
  
  std::map<char,int> lab2num;
  
  for (int k=0; k<K; k++)
    {      
      char U = toupper( ms_labels[ k ] );      
      if ( lab2num.find( U ) == lab2num.end() )
	lab2num[ U ] = k;
    }

  // assign homes
  ul_mapping.clear();
  
  for (int k=0; k<K; k++)
    {
      ul_mapping[ k ] = lab2num[ toupper( ms_labels[ k ] ) ]; 
      //      std::cout << " assigning " << k << " --> " << ul_mapping[ k ]  << "  --- " << ms_labels[ k ] << "\n";
    }
  
}

int ms_prototypes_t::ul_reduction(const int l)
{
  // should always be present
  if ( ul_mapping.find( l ) == ul_mapping.end() )
    Helper::halt( "internal error in ul-mapping lookups" );
  return ul_mapping[ l ];
}


void microstates_t::aggregate2edf( const Data::Matrix<double> & X ,
				   const signal_list_t & signals ,
				   const std::vector<int> & peak_idx , 
				   const int srate ,
				   const double pmin ,
				   const double pmax , 
				   const std::string & edfname )
{

  if ( ! Helper::file_extension( edfname , "edf" ) )
    Helper::halt( "peaks file should have an .edf extension:" + edfname );

  bool exists = Helper::fileExists( edfname );

  // Nb:: currentyly, just look at whether it exists or no... i.e.
  //  as we might write from multiple processes... and so user has to delete file
  //  if they want to overwrite
  
  //
  // Extract peaks
  //

  bool has_peak_list = peak_idx.size() > 0 ;

  Data::Matrix<double> Z( has_peak_list ? peak_idx.size() : X.dim1() , X.dim2() );

  if ( has_peak_list )
    {
      const int n_peaks = peak_idx.size();
      const int nc = X.dim2();
      for (int r=0; r<n_peaks; r++)
	for (int c=0; c<nc; c++)
          Z( r , c ) = X( peak_idx[r] , c ) ;
    }
  else
    Z = X;
  
  
  
  //
  // key variables
  //

  const int N = Z.dim1();  // i.e. not transposed here
  const int C = Z.dim2(); 

  // so we can do easy appends, set n_samples = 1 for each channel, i.e. one sample per EDF record
  //  we then just change NR in the header when we add stuff;   we can set the SR (rec_dur) to 1, it
  //  will not matter
  
  //
  // write a new EDF header and first set of records
  //
  
  if ( ! exists ) 
    {
      
      logger << "  writing an aggregated GPF-peak EDF to " << edfname << "\n";
      
      edf_t agg_edf;
      
      const int recdur = 1;
      const int nr = N;
      const int ns = C;
      
      //
      // Set header
      //
      
      agg_edf.header.version = "0";
      agg_edf.header.patient_id = "GFP";
      agg_edf.header.recording_info = ".";
      agg_edf.header.startdate = ".";
      agg_edf.header.starttime = ".";
      agg_edf.header.nbytes_header = 256 + ns * 256;
      agg_edf.header.ns = 0;        // these will be added by add_signal()
      agg_edf.header.ns_all = ns;   // check this... should only matter for EDF access, so okay... 
      agg_edf.header.nr = agg_edf.header.nr_all = nr;   // likewise, value of nr_all should not matter, but set anyway
      agg_edf.header.record_duration = 1;    // arbitrary for an aggregate GPF-peak EDF
      agg_edf.header.record_duration_tp = agg_edf.header.record_duration * globals::tp_1sec; 

      
      //
      // create a (continuous) timeline  [ TODO: see about SEDF+ for masked EDF ] 
      //
      
      agg_edf.set_edf();
      agg_edf.set_continuous();
      agg_edf.timeline.init_timeline();

      //
      // resize data[][], by adding empty records (one per SEDF record == EDF epoch )
      //

      for (int r=0;r<nr;r++)
	{
	  edf_record_t record( &agg_edf ); 
	  agg_edf.records.insert( std::map<int,edf_record_t>::value_type( r , record ) );
	}

      //
      // add signals (this populates channel-specific 
      //

      for (int c=0;c<C;c++)
	{	  
	  // -1 implies 1 sample per record, i.e. if positive would be the SR
	  // for that signal, but this allows slower channels, by directly specifying
	  // the samples per record, (rather than samples per second)
	  
	  const std::vector<double> * p = Z.col(c).data_pointer();
	  
	  // nb. as we will add new records, it is our responsibility to set pmin and pmax
	  // such that they are reasonable values for all (future-aggregated) signals too
	  // [ otherwise, the newly appended data will be clipped at these values ]
	  
	  agg_edf.add_signal( signals.label(c) , -1  , *p , pmin , pmax );

	}
      
      //
      // Save this SEDF
      //
      
      bool saved = agg_edf.write( edfname );
      
      if ( ! saved ) Helper::halt( "problem trying to write " + edfname );
      
    }
  else
    {

      logger << "  appending to an existing aggregated GPF-peak EDF " << edfname << "\n";
 
      // data[rec][channels][samples]

      std::vector<std::vector<std::vector<double> > > data(N);

      // here, each record contaions only a single sample (n_samples[s] == 1), i.e. 1 GFP peak
      // i.e. so we can add new records without having to change n_samples[] etc )
      //      to the aggregate EDF

      // # records --> # of GFP-peaks (N)
      // NS        --> C
      // samples per record --> 1 


      for (int i=0;i<N;i++)
	{
	  data[i].resize(C);
	  for (int j=0;j<C;j++)
	    {
	      data[i][j].resize(1);
	      data[i][j][0] = Z(i,j);
	    }
	}

      //
      // And channel names
      //

      std::vector<std::string> channels( C );
      for (int c=0; c<C; c++) channels[c] = signals.label(c);
      
      //
      // edf_t utility function to append this strucrture to an EDF on disk
      //

      if ( ! edf_t::append( edfname , channels , data ) )
	Helper::halt( "problem appending new data" );
      
    }

}


char ms_kmer_t::pick( const std::map<char,int> & urns , char skip )
{
  int tot = 0;
  std::map<char,int>::const_iterator ii = urns.begin();
  std::vector<int> counts;
  std::vector<char> labels;
  while ( ii != urns.end() )
    {
      if ( ii->first != skip || ii->second == 0 )
	{
	  tot += ii->second;
	  counts.push_back( ii->second );
	  labels.push_back( ii->first );
	  //	  std::cout << " urn " << ii->first << " = " <<  ii->second  << "\n";
	}
      ++ii;
    }

  // if only one class left, and that was skipped previously, then we will not
  // be able to satisfy the constraint...  for now, here we have to return
  // the skipped value...

  if ( tot == 0 ) return skip;
  
  int rn = CRandom::rand( tot );
  int p = 0;
  //  std::cout << " rn / tot = " << rn << " " << tot << "\n";
  while ( 1 )
    {
      //std::cout << "p, rn, cnt = " << p << " " << rn << " " << counts[p] << "\n";
      if ( rn < counts[p] ) break;
      rn -= counts[p];
      ++p;
    }
  //  std::cout << " PICK " << p << " " << labels[p] << "\n";
  return labels[p];
}


std::string ms_kmer_t::modified_random_draw( const std::string & l , const int w )
{

  // permute 's' chars but in such a way that similar states (i.e. values of l, 0, 1, 2, ...)
  // are not contiguous
  // first element: pick any element (e.g. of 5)
  // next elemnt: pick from the remaining 5-1 elements, i.e. not the previously selected
  // repeat, until all done

  // two modes:
  //  global, then 'urns' is only made once, based on all states in the sequence
  //  local (W), then 'urns' is based on a window of W states 
  //    this is updated based on the W+1 state as we shift each

  // A B A C D E A C B A E C B
  // X
  // -----------

  //  X
  // -----------

  const bool global_picks = w == 0 ; 
   
  std::map<char,int> urns;

  // ovreall sequence length
  const int n = l.size();

  // urn total (but cap at total seq. length)
  const int nurn = global_picks ? n : ( w > n ? n : w )  ; 
  
  // fill urn (from start)
  for (int i=0;i<nurn;i++)
    ++urns[l[i]];
  
  // return string
  std::string p(n,'.');
  
  // initiate (from first of set, i.e. == urn contents)
  p[0] = l[ CRandom::rand( nurn ) ];
  
  char last = p[0];
  --urns[ last ];

  int lc = w;
  
  for (int i=1;i<n;i++)
    {

      // if ( 1 )
      // 	{
      // 	  std::cout << " SLOT = " << i << "\n";
      // 	  std::map<char,int>::const_iterator uu = urns.begin();
      // 	  while ( uu != urns.end() )
      // 	    {
      // 	      std::cout << uu->first << " = " << uu->second << "\n";
      // 	      ++uu;
      // 	    }
      // 	}

      // attempt to skip last pick
      char c = pick( urns , last );
      
      // typically, we expect to be able to pick one 
      // that is not equal to the last pick:
      if ( c != last ) 
	{
	  p[i] = c;
	  --urns[ c ];
	  last = c;
	}
      else // however, this might not always be the case...
	{
	   // std::cout << "\n\nfound issue:\n";
	   // std::cout << "pick i = " << i << " " << c << "\n";
	   // std::cout << "current:\n" << p.substr(0,i) << "?\n";

	  // i.e. if being forced to pick something the same as the prior, 
	  //   ABCABA[A]
	  //   012345 6
	  // we want to randomly find a prior spot where the 'A' can fit (i.e. between B and C)
	  // and insert it there, and move everything else forward by one...
	  //  nb; this issue should neve happen with only two classes: i.e. as sequence stats
	  //      don't make sense there in any case... 
	  
	  // pick a random starting slot between 0 and i-2
	  //  i.e. could go in slot 4 above; (as slot 5 (i-1) == i candidate) 
	  
	  int search = CRandom::rand( i - 2 );

	  // should not need this... but in case of weird
	  // edge cases, makes sure we bail rather than get stuck in 
	  // an infinite loop

	  bool looped = false;

	  while ( 1 ) 
	    {
	      //std::cout << " searching = " << search << "\n";
	      
	      // valid position?
	      bool okay = search == 0 ? p[0] != c : ( p[search-1] != c &&  p[search] != c ) ;

	      if ( ! okay ) 
		{ 
		  ++search;
		  if ( search == i-1 ) 
		    {
		      if ( looped ) 
			Helper::halt( "internal problem in modified_random_draw()... cannot create sequence" );
		      search = 0 ; 
		      looped = true;
		    }
		  continue;
		}
	      
	      // otherwise, shift everything else /after/ search forward by 1, and then insert at search spot
	      
	      for (int s = i ; s != search; s-- )
		p[s] = p[s-1];
	      
	      p[search] = c;
	      
	      // std::cout << "placed\n"
	      // 		<< p.substr( 0 , i+1 ) << "\n";
	      break;
		  
	    }
	  
	  // still need to drop urn count:
	  --urns[ c ];
	  
	  // and the last will still be 'c' (i.e. as this was also the prior pick)
	  last = c;
	  
	}
      
      // sanity check, can remove
      if ( urns[p[i]] < 0 ) Helper::halt( "error!" );      

      // in local mode, we need to fill up the urn
      if ( ! global_picks )
	{
	  if ( lc < n ) 
	    ++urns[l[lc++]];
	}
      
      // next pick
      
    }

  // sanity check... can remove 
  std::map<char,int>::const_iterator ll = urns.begin();
  while ( ll != urns.end() )
    {
      if ( ll->second != 0 ) Helper::halt( "bad" );
      ++ll;
    }
  
  // all done
  return p;
    
}




Data::Matrix<double> microstates_t::eig2mat( const Eigen::MatrixXd & E )
  {
    const int rows = E.rows();
    const int cols = E.cols();
    Data::Matrix<double> M( rows , cols );
    for (int r=0;r<rows;r++)
      for (int c=0;c<cols;c++)
	M(r,c) = E(r,c);
    return M;
  }

Eigen::MatrixXd microstates_t::mat2eig( const Data::Matrix<double> & M )
{
  const int rows = M.dim1();
  const int cols = M.dim2();
  Eigen::MatrixXd E( rows , cols );
  for (int r=0;r<rows;r++)
    for (int c=0;c<cols;c++)
      E(r,c) = M(r,c);
  return E;
}

Eigen::MatrixXd microstates_t::mat2eig_tr( const Data::Matrix<double> & M )
{
  const int rows = M.dim1();
  const int cols = M.dim2();
  Eigen::MatrixXd E( cols , rows );
  for (int r=0;r<rows;r++)
    for (int c=0;c<cols;c++)
      E(c,r) = M(r,c);
 return E;
}


void ms_prototypes_t::map_to_canonicals( const std::string & filename )
{
  // spatial correlation
  // expect prototype file with headers that are labels
  // take current maps A[][] and for each figure out the closest canonical state 
  // (based on canonicals) 
  // canonical file need not have the same channels

  std::string filename2 = Helper::expand( filename );
  if ( ! Helper::fileExists( filename2 ) ) 
    Helper::halt( "could not find canonical prototype file " + filename2 );
  
  // tmp store
 
  std::map<std::string,std::map<char,double> > ch2k2a;
  

  //
  // get header
  //

  std::ifstream IN1( filename2.c_str() , std::ios::in );

  std::string line;
  Helper::safe_getline( IN1 , line );
  if ( IN1.eof() || line == "" ) Helper::halt( "invalid header for " + filename2 );
  std::vector<std::string> tok = Helper::parse( line , "\t " );  
  if ( tok.size() < 2 ) Helper::halt( "bad format for " + filename2 );
  const int nmaps = tok.size() - 1 ; 
  
  if ( tok[0] != "CH" ) 
    Helper::halt( "column 1 should have header 'CH'" );
  
  std::vector<std::string> classes;

  for (int i=1; i<tok.size(); i++)
    {
      if ( tok[i].size() != 1 ) 
	Helper::halt( tok[i] + " -- state labels can only be single characters," + filename2 );
      classes.push_back( tok[i] );
    }

  // canonical map names are tok[]  (starting from X=1 to X <= nmaps)  

  while ( ! IN1.eof() )
    {
      std::string line;
      Helper::safe_getline( IN1 , line );
      if ( IN1.eof() ) break;
      if ( line == "" ) continue;
      std::vector<std::string> tok = Helper::parse( line , "\t " );
      if ( tok.size() != nmaps + 1 ) Helper::halt( "bad ... " );
      
      const std::string & channel = tok[0];

      for (int i=1; i<tok.size(); i++)
	{
	  double a;
	  if ( ! Helper::str2dbl( tok[i] , &a ) )
	    Helper::halt( "problem reading value: " + tok[i] );
	  ch2k2a[ channel ][ classes[ i-1 ][0] ] = a;
	}
      
    }
  IN1.close();

  // 
  // Create matrix in same order as A[][] 
  //

  const int canK = classes.size();
  
  Eigen::MatrixXd CA = Eigen::MatrixXd::Zero( C , canK );

  for (int c=0; c<C; c++)
    {
      // was this channel specified?
      if ( ch2k2a.find( chs[ c ] ) == ch2k2a.end() )
	Helper::halt( "could not find channel " + chs[c] + " in " + filename2 );
      
      for (int k=0; k<canK; k++)
	CA(c,k) = ch2k2a[ chs[c] ][ classes[k][0] ];

    }

  //
  // Now we have CA with canK 'canonical' prototypes
  // and our estimated/observed prototypes with K maps
  // and these are oriented for the same channels
  
  
  //
  // Find each observed map, find best match in canA
  //

  
  //
  // Normalize both maprs 
  //

  Eigen::MatrixXd ZA = A;
  eigen_ops::scale( ZA , true , true );  
  eigen_ops::scale( CA , true , true );

  Eigen::MatrixXd R( K , canK );

  for (int k=0; k<K; k++)
    for (int k2=0; k2<canK; k2++)
      R(k,k2) = spatial_correlation( ZA.col(k) , CA.col(k2) );
  
  std::cout << R << "\n";

  
  
  // For each 
  
  //
  // Re-order A to 
  //

  // Eigen::MatrixXd::Index maxIndex[2];
  // VectorXf maxVal(2);
  // for(int i=0;i<2;++i)
  //   maxVal(i) = mat.row(i).maxCoeff( &maxIndex[i] );


  // Eigen::MatrixXd A; // C x K                                                                                                                  

}


double ms_prototypes_t::spatial_correlation( const Eigen::VectorXd & M1, const Eigen::VectorXd & M2 , bool * flip )
{
  
  // calculate spatial correlation between two maps 
  // optionally, return 'flip==T' if we need to flip a map (for viz)
  
  const int nc = M1.size();
  
  if ( M2.size() != nc ) 
    Helper::halt( "internal error in spatial_correlation() : different channel N" );
  
  // Global Map Dissimilarity (GMD)
  
  double t = 0 , t2 = 0;
  for (int i=0; i<nc; i++)
    {
      t +=  ( M1[i] - M2[i] ) * ( M1[i] - M2[i] ) ;
      t2 += ( M1[i] + M2[i] ) * ( M1[i] + M2[i] ) ;
    }
  
  t = sqrt( t / (double)nc );
  t2 = sqrt( t2 / (double)nc );

  // pick smallest distance to ensure polarity invariance
  double gmd = t < t2 ? t : t2 ;

  if ( flip != NULL ) *flip = t2 < t ; 
  
  // return spatial correlation
  return 1 - ( gmd * gmd ) / 2.0; 
  
}


// confidence threshold 
void ms_backfit_t::determine_ambiguity( double conf , double th2 ) 
{
  const int K = GMD.dim1();
  const int N = GMD.dim2();

  ambiguous.resize( N , false );

  // note: lots of redundancy here... streamline these calcs...

  // get SPC
  Data::Matrix<double> SPC = GMD;
  for (int k=0; k<K; k++)
    for (int i=0; i<N; i++)
      SPC(k,i) = 1 - ( SPC(k,i) * SPC(k,i) ) / 2.0;
  
  int cnt = 0; 

  // get CONF
  for (int i=0; i<N; i++)
    { 
      double best = 0;
      double next = 0;
      for (int k=0; k<K; k++)
	{
	  if ( SPC(k,i) > best )
	    {
	      next = best;
	      best = SPC(k,i);				  
	    }
	  else if ( SPC(k,i) > next )
	    next = SPC(k,i); 
	}			  
      
      // set to 100 as max
      double c = next > 0 ? best / next : 100; 
      
      // mark as ambiguous
      if ( c < conf || best < th2 ) 
	{ ambiguous[i] = true; ++cnt; } 
    }

  // // 
  // // any fill-ins?
  // // 

  // if ( fillin_samples )
  //   {
  //     int prior = -1;
  //     for (int i=1; i<N-1; i++)
  // 	{
  // 	  if ( ! ambiguous[i] ) 
  // 	    prior = 

  // 	  if ( ambiguous[i] ) 
	    
  // 	}

  //   }

  logger << "  set " << Helper::dbl2str_fixed( 100 * cnt / (double)N , 2 ) << "% points as ambiguous\n";
  
}
  




ms_cmp_maps_t::ms_cmp_maps_t( const std::map<std::string,std::map<std::string,std::map<std::string,double> > > & data ,
			      const Eigen::MatrixXd * fixed , 
			      const std::vector<std::string> * fixed_chs , 
			      const std::map<std::string,int> & phe ,
			      const int nreps ,
			      const double p )
{

  // if fixed is non-NULL, then do a C/C comparison against this fixed map
  // otherwise, compare all concordant pairs / all discordant pairs
  
  // d: id -> k -> ch -> a 
  
  std::vector<int> g;
  
  // track channels and K's (i.e. to check that all is the same)
  std::set<std::string> chs;
  std::set<std::string> ks;
  std::vector<std::string> ids;
  
  std::map<std::string,int>::const_iterator ii = phe.begin();
  while ( ii != phe.end() )
    {
      // no data?
      if ( data.find( ii->first ) == data.end() )
	{ ++ii; continue; } 

      // missing phenotye? (not 0 or 1)
      if ( ii->second != 0 && ii->second != 1 ) { ++ii; continue; }

      // track 0 or 1 phenotype
      g.push_back( ii->second );

      ids.push_back( ii->first );
      
      // make matrix
      // track channels
      const std::map<std::string,std::map<std::string,double> > & dat = data.find( ii->first )->second;
      
      // track # of classes
      std::map<std::string,std::map<std::string,double> >::const_iterator kk = dat.begin();
      while ( kk != dat.end() )
	{
	  ks.insert( kk->first ); // track class
	  const std::map<std::string,double> & dat2 = kk->second;
	  std::map<std::string,double>::const_iterator cc = dat2.begin();
	  while ( cc != dat2.end() )
	    {
	      chs.insert( cc->first );
	      ++cc;
	    }
	  ++kk;
	}
      
      ++ii;
    }

  //
  // size of the problem 
  //

  const int ni = g.size();  
  const int nc = chs.size();
  const int nk = ks.size();

  //
  // the actual data
  //

  std::vector<Eigen::MatrixXd> m( ni );

  for (int i=0; i<ni; i++)
    {
      Eigen::MatrixXd X = Eigen::MatrixXd::Zero( nc , nk );

      const std::map<std::string,std::map<std::string,double> > & dat = data.find( ids[i] )->second;      
      if ( dat.size() != nk ) Helper::halt( "uneven data1") ;

      // classes = columns
      int c = 0;
      std::set<std::string>::const_iterator kk = ks.begin();
      while ( kk != ks.end() )
	{
	  
	  const std::map<std::string,double> & dat2 = dat.find( *kk )->second;	  
	  if ( dat2.size() != nc ) Helper::halt( "uneven data2") ;

	  
	  int r = 0;
	  std::set<std::string>::const_iterator cc = chs.begin();
	  while ( cc != chs.end() )
	    {
	      X( r , c ) = dat2.find( *cc )->second;
	      // next row/channel
	      ++r;
	      ++cc;
	    }
	  // next column/class
	  ++kk;
	  ++c; 
	}

      // store data for this individual...
      m[i] = X;
      
      // next individual

    }
  
  //
  // Now all data should be nicely loaded in and squared off
  //


  //
  // Individual vs. template comparisons? 
  //

  if ( fixed != NULL )
    {
      
      logger << "  comparing individuals to a fixed template\n";
     
      Eigen::MatrixXd FX;
      
      const int nt = fixed->cols();
      if ( nt < nk ) Helper::halt( "template contains fewer prototypes than data" );
      if ( fixed->rows() != nc ) 
	{
	  logger << "  template = " << fixed->rows() << " channels\n"
		 << "  data     = " << nc << " channels\n";
	  Helper::halt( "fixed map has bad # of channels" );
	}
      FX = Eigen::MatrixXd::Zero( nc , nt );
      std::map<std::string,int> ch2row;
      for (int r=0; r<nc; r++) ch2row[ (*fixed_chs)[r] ] = r;
      std::set<std::string>::const_iterator cc = chs.begin();
      int rn = 0;
      while ( cc != chs.end() )
	{
	  if ( ch2row.find( *cc ) == ch2row.end() )
	    Helper::halt( "mismatch of channel labels" );	  
	  for (int j=0; j<nt; j++)
	    FX(rn,j) = (*fixed)(ch2row[*cc],j);
	  ++rn;
	  ++cc;
	}
      
      Eigen::VectorXd R = Eigen::VectorXd::Zero( ni );
      
      // store best matches?
      // after reading template, (with prototypes::read())
      // the static ms_prototypes_t::ms_labels[] will contain
      // the T canonical labels
      // nb, we allow for the template to have more classes than
      // the data;  we use a special comparison function cmp_map_templates()
      // that returns the optimal match returnning a subset of nk of the total nt
      
      std::vector<std::vector<int> > best(ni);
      for (int i=0;i<ni; i++) best[i].resize( nk );
      
      // always brute-force ; special comparison as we
      // allow FX to have more classes than 'm'
      for (int i=0;i<ni; i++)
	R[i] = cmp_maps_template( m[i] , FX , p, &(best[i]) ) ;


      //
      // Statistics on R (distance matrix)
      //
      
      // Original data
      
      std::vector<int> perm( ni );
      for (int i=0;i<ni;i++) perm[i] = i;
      
      // | case-F - control-F |
      double within[2]; // also , mean 11, and mean 00 groups
      double het_between = het_template_statistic( g , perm , R , within );
  
      // Permutations      
      int r_het = 0 , r_w_cas = 0 , r_w_con = 0;
      for (int p = 0; p < nreps ; p++ )
	{

	  // shuffle
	  CRandom::random_draw( perm );
	  
	  // within group stats
	  double pwithin[2];
	  double phet = het_template_statistic( g , perm , R , pwithin );
	  if ( phet >= het_between ) ++r_het;      
	  
	  // n.b sign of comparisons:: testing if /less distant/
	  if ( pwithin[0] <= within[0] ) ++r_w_cas;
	  if ( pwithin[1] <= within[1] ) ++r_w_con;
	}
      
      double p_het = ( r_het + 1 ) / (double)( nreps + 1 );
      double p_cas = ( r_w_cas + 1 ) / (double)( nreps + 1 );
      double p_con = ( r_w_con + 1 ) / (double)( nreps + 1 );  
      
      logger << "  case-template (more similar?)                  : " << within[0] << " p = " << p_cas   << "\n";
      logger << "  control-template (more similar?)               : " << within[1] << " p = " << p_con   << "\n";
      logger << "  | case-template - control-template | distance  : " << het_between << " p = " << p_het   << "\n";
      
      writer.id( "." , "." );
      writer.level( "ALL" , "SUMM" );
      writer.value( "P_TEMPLATE_HET" , p_het );
      writer.value( "P_TEMPLATE_CAS" , p_cas );
      writer.value( "P_TEMPLATE_CON" , p_con );  
      writer.unlevel( "SUMM" );

      for (int i=0;i<ni; i++)
	{
	  writer.id( ids[i] , "." );
	  writer.value( "S" , R[i] );

	  // matches w/ the canonical template
	  for (int k=0;k<nk;k++)
	    {
	      writer.level( k , "SLOT" );
	      writer.value( "T" , std::string( 1, ms_prototypes_t::ms_labels[ best[i][k] ] ) );
	    }
	  writer.unlevel( "SLOT" );
	}
      
      // all done now
      return;
    }

  
  
  //
  // Otherwise, report pairwise comparisons
  //
  
  logger << "  creating individual-by-individual global distance matrix\n";
    
  Eigen::MatrixXd R = Eigen::MatrixXd::Zero( ni, ni );

  for (int i=0;i<ni; i++)
    for (int j=i+1;j<ni;j++)
      R(i,j) = R(j,i) = cmp_maps_bf( m[i] , m[j] , p ) ;
  
  //
  // Original data
  //

  std::vector<int> perm( ni );
  for (int i=0;i<ni;i++) perm[i] = i;

  // individual-level summary statistics
  Eigen::VectorXd ires( ni );

  // concordant pairs / discordant pairs
  double pobs = statistic( g , perm , R , &ires );

  // | case-pairs - control-pairs |
  double within[2]; // also , mean 11, and mean 00 groups
  double het_between = het_statistic( g , perm , R , within );
  
  
  //
  // Permutations
  //

  int r = 0 ; 
  int r_het = 0 , r_w_cas = 0 , r_w_con = 0;
  
  for (int p = 0; p < nreps ; p++ )
    {
      // shuffle
      CRandom::random_draw( perm );

      // pairwise stats
      double pst = statistic( g , perm , R , NULL );
      if ( pst >= pobs ) ++r;      

      // within group stats
      double pwithin[2];
      double phet = het_statistic( g , perm , R , pwithin );
      if ( phet >= het_between ) ++r_het;      

      // n.b signs here - testing whether phenotypic group is /less distant/
      // than expected by chance (to eachother)
      if ( pwithin[0] <= within[0] ) ++r_w_cas;
      if ( pwithin[1] <= within[1] ) ++r_w_con;
    }
  
  double p_pairs = ( r + 1 ) / (double)( nreps + 1 );
  double p_het = ( r_het + 1 ) / (double)( nreps + 1 );
  double p_cas = ( r_w_cas + 1 ) / (double)( nreps + 1 );
  double p_con = ( r_w_con + 1 ) / (double)( nreps + 1 );  
  
  logger << "  within-case (more similar?)                : " << within[0] << " p = " << p_cas   << "\n";
  logger << "  within-control (more similar?)             : " << within[1] << " p = " << p_con   << "\n";
  logger << "  | within-case - within-control | distance  : " << het_between << " p = " << p_het   << "\n";
  logger << "  concordant / discordant pair distance      : " << pobs  << " p = " << p_pairs << "\n";
  
  writer.id( "." , "." );
  writer.level( "ALL" , "SUMM" );
  writer.value( "P_PAIRS" , p_pairs );
  writer.value( "P_HET" , p_het );
  writer.value( "P_CAS" , p_cas );
  writer.value( "P_CON" , p_con );  
  writer.unlevel( "SUMM" );

  for (int i=0;i<ni; i++)
    {
      writer.id( ids[i] , "." );
      writer.value( "S" , ires[i] );
    }
  
}

//
// Get statistic (C/C and also for each person) given R
//

double ms_cmp_maps_t::statistic( const std::vector<int> & phe ,
				 const std::vector<int> & perm ,
				 const Eigen::MatrixXd & R ,
				 Eigen::VectorXd * ires )
{

  // if ires != NULL, return summary statistic for mean similarity of each
  // person to all others

  if ( ires != NULL )
    {
      *ires = R.array().colwise().sum();
      *ires /= (double)(R.rows() - 1);
    }

  
  // C/C permutation
  // only count discordant C/C pairs
  // return sum

  double st = 0;
  double st_within = 0;
  
  const int ni = R.rows();

  int n_disc = 0 , n_conc = 0;
  
  for (int i=0; i<ni; i++)
    {
      for (int j=0; j<ni; j++)
	{
      
	  // if only looking at discordant pairs, we don't
	  // need to worry abount count self matches (although
	  // those should be set to 0.0 in any case)
	  
	  if ( phe[ perm[ i ] ] != phe[ perm[ j ] ]  )
	    {
	      st += R( i , j );
	      ++n_disc;
	    }
	  else
	    {
	      st_within += R( i , j );
	      ++n_conc;
	    }
	}
    }
  
  // discordant pairs / concordant pairs distance
  return (st/(double)n_disc )  / ( st_within/(double)n_conc) ;
  
}


double ms_cmp_maps_t::het_statistic( const std::vector<int> & phe ,
				     const std::vector<int> & perm ,
				     const Eigen::MatrixXd & R ,
				     double *within )
{

  
  // this only looks at case-case and control-control pairs
  //  (versus the other statistic is based on case-control versus other pairs)
  
  double st = 0;
  double st_cas = 0 , st_con = 0;
  int n_cas_pairs = 0 , n_con_pairs = 0;  

  const int ni = R.rows();
  
  for (int i=0; i<ni; i++)
    for (int j=0; j<ni; j++)
      {

	// here, ignore discordant pairs
	if ( phe[ perm[ i ] ] != phe[ perm[ j ] ]  ) continue;
	
	if ( phe[ perm[ i ] ] == 1 )
	  {
	    st_cas += R( i , j );
	    ++n_cas_pairs;
	  }
	else
	  {
	    st_con += R( i , j );
            ++n_con_pairs;
	  }
      }

  // mean within-case distance
  within[0] = st_cas / (double)n_cas_pairs;

  // mean within-control distance
  within[1] = st_con / (double)n_con_pairs;
  
  // | case-case - control-control |
  return fabs( within[0] - within[1] );
  
  
}



double ms_cmp_maps_t::het_template_statistic( const std::vector<int> & phe ,
					      const std::vector<int> & perm ,
					      const Eigen::VectorXd & R ,					      
					      double *within )

{

  
  // this only looks at case-case and control-control pairs
  //  (versus the other statistic is based on case-control versus other pairs)
  
  double st = 0;
  double st_cas = 0 , st_con = 0;
  int n_cas = 0 , n_con = 0;  

  const int ni = R.rows();
  
  for (int i=0; i<ni; i++)
    {
      if ( phe[ perm[ i ] ] == 1 )
	{
	  st_cas += R[i] ;
	  ++n_cas;
	}
      else
	{
	  st_con += R[i];
	  ++n_con;
	}
    }
  
  // mean case-template distance
  within[0] = st_cas / (double)n_cas;
  
  // mean control-template distance
  within[1] = st_con / (double)n_con;
  
  // | case-case - control-control |
  return fabs( within[0] - within[1] );
  
  
}






//
// Get spatial correlation given two maps, and find the best join: greedy
//

// double ms_cmp_maps_t::cmp_maps( const Eigen::MatrixXd & A , const Eigen::MatrixXd & B )
// {

//   // here, we can assume that A and B are of the same dimension
//   // we can also assume that maps are normalized
  
//   const int nk = A.cols();

//   // get matrix of spatial correlations
  
//   Eigen::MatrixXd R = Eigen::MatrixXd::Zero( nk , nk );
//   for (int i=0; i<nk; i++)
//     for (int j=0; j<nk; j++)
//       R(i,j) = ms_prototypes_t::spatial_correlation( A.col(i) , B.col(j) );
  

//   // find best match, greedily
  
//   double res = 0;
//   for (int k=0; k<nk; k++)
//     {
//       Eigen::Index r0, c0;      
//       R.maxCoeff(&r0,&c0);
//       res += R(r0,c0);
//       R.row(r0) = Eigen::VectorXd::Zero( nk );
//       R.col(c0) = Eigen::VectorXd::Zero( nk );      
//     }

//   return res / (double)nk;
  
// }



//
// Get spatial correlation given two maps, and find the best join: brute-force
//

double ms_cmp_maps_t::cmp_maps_bf( const Eigen::MatrixXd & A , const Eigen::MatrixXd & B , double p )
{

  const int nk = A.cols();
  Eigen::MatrixXd R = Eigen::MatrixXd::Zero( nk , nk );
  for (int i=0; i<nk; i++)
    for (int j=0; j<nk; j++)
      R(i,j) = ms_prototypes_t::spatial_correlation( A.col(i) , B.col(j) );
  
  // find best match, brute-force over all combinations
  //  for K pairs , find all possible matches

  // keep person 'A' fixed (1, 2, 3, ..., K)
  // then make all possible permutations of 'B' vector

  // metric = minimum of sum (1-R)^p / K 
  
  std::vector<int> kb( nk );
  for (int i=0; i<nk; i++) kb[i] = i;

  double min_res = 999;
  
  do {
    double res = 0;
    for (int k=0; k<nk; k++)
      res += pow( 1.0 - R(k,kb[k]) , p );
    if ( res < min_res )
      min_res = res;
    
  } while ( std::next_permutation( kb.begin() , kb.end() ) );
    
  return min_res / (double)nk;
  
}


//
// Get spatial correlation given a map and a 'template' set;
// Using brute-force enumeration of all possibilities, here we allow the template (T)
// to have more columns than map "A" and we select the best
//

double ms_cmp_maps_t::cmp_maps_template( const Eigen::MatrixXd & A ,
					 const Eigen::MatrixXd & T ,
					 double p ,
					 std::vector<int> * best )
{

  const int nk = A.cols();
  const int nt = T.cols(); 

  // nk x nt correlation matrix:
  Eigen::MatrixXd R = Eigen::MatrixXd::Zero( nk , nt );
  for (int i=0; i<nk; i++)
    for (int j=0; j<nt; j++)
      R(i,j) = ms_prototypes_t::spatial_correlation( A.col(i) , T.col(j) );
  
  // find best match, brute-force over all combinations
  //  for K pairs , find all possible matches

  // keep person 'A' fixed (1, 2, 3, ..., K)
  // then make all possible permutations of 'B' vector
  //  by permuting all 'nt' labels, but only selecting the first 'nk'
  // involves redundant work, but should not be an issue really

  // metric as above
  
  std::vector<int> kt( nt );
  for (int i=0; i<nt; i++) kt[i] = i;
  
  double min_res = 999;
  
  do {
    double res = 0;    

    // nb, here only looking at the first nk of nt
    for (int k=0; k<nk; k++) 
      res += pow( 1.0 - R(k,kt[k]) , p );
    
    if ( res < min_res )
      {
	min_res = res;
	if ( best != NULL ) *best = kt; 
      }
    
  } while ( std::next_permutation( kt.begin() , kt.end() ) );
  
  // reduce 'best' down to the first nk elements from nt
  if ( best != NULL )
    best->resize( nk );
  
  return min_res / (double)nk;
  
}


std::vector<char> ms_cmp_maps_t::label_maps( const ms_prototypes_t & T , 
					     const std::vector<char> & Tl ,
					     ms_prototypes_t * A , 
					     const std::vector<char> & Al ,
					     double minr ,
					     double p,
					     bool verbose )

{

  //
  // same logic as cmp_maps_template() above
  // but also flip polarity of A to match best-fit T
  //
  
  // get number of maps from a) template T , b) candidate map set A
  const int nt = T.K;
  const int nk = A->K;

  if ( nk > nt )
    Helper::halt( "template must have same or larger number of maps than the target" );
  
  std::vector<char> r(nk);
  std::vector<int> best(nt);
  std::vector<double> spatialr(nk);
  std::vector<bool> flip(nk);
  
  // nk x nt correlation matrix:
  Eigen::MatrixXd R = Eigen::MatrixXd::Zero( nk , nt );
  for (int i=0; i<nk; i++)
    for (int j=0; j<nt; j++)
      R(i,j) = ms_prototypes_t::spatial_correlation( A->A.col(i) , T.A.col(j) );
  
  // find best match, brute-force over all combinations
  //  for K pairs , find all possible matches
  
  // keep person 'A' fixed (1, 2, 3, ..., K)
  // then make all possible permutations of 'B' vector
  //  by permuting all 'nt' labels, but only selecting the first 'nk'
  // involves redundant work, but should not be an issue really

  if ( verbose )
    {
      std::cout << "R\n" << R << "\n\n";
    }
  
  std::vector<int> kt( nt );
  for (int i=0; i<nt; i++) kt[i] = i;
    
  double min_res = 999;
  
  do {

    double res = 0;
    
    // nb, here only looking at the first nk of nt    
    for (int k=0; k<nk; k++) 
      res += pow( 1 - R(k,kt[k]) , p );
            
    // verbose?
    if ( verbose )
      {
	double mr = 1;
	for (int k=0; k<nk; k++)
	  {
	    if ( mr > R(k,kt[k]) ) mr = R(k,kt[k]) ;
	    std::cout << "  "
		      << Al[k] << ">" << Tl[ kt[k] ] 
		      << " " << R(k,kt[k])
		      << ""; 
	  }
	std::cout << " S = " << res << "\n";
      }

    // best?
    if ( res < min_res )
      {

	min_res = res;

	best = kt;
	
	for (int k=0; k<nk; k++)
	  {
	    // recalc to get polarity 
	    //  spatialr[k] = R(k,kt[k]);

	    bool f;
	    spatialr[k] = ms_prototypes_t::spatial_correlation( A->A.col(k) , T.A.col(kt[k]) , &f );
	    flip[k] = f;
	  }
      }
    
  } while ( std::next_permutation( kt.begin() , kt.end() ) );


  bool all_optimal = true;
  bool all_mapped = true;
  
  // copy first 'nk' elements from T to A
  for (int k=0; k<nk; k++)
    {
      bool okay = spatialr[k] >= minr;

      char newlab = okay ? Tl[ best[k] ] : '?';
      
      logger << "   mapping [" << Al[k] << "] --> template [" << Tl[ best[k] ] << "] with R = " << spatialr[k];      
      if ( flip[k] ) logger << " [flip]";
      else logger << "       ";
      if ( ! okay ) logger << " *** below threshold corr. -- assigning '?'";
      logger << "\n";

      if ( ! okay ) all_mapped = false;
      
      // is this an 'optimal' match (i.e. no other template with a higher SPC?
      bool optimal = true;
      for (int j=0; j<nt; j++)
	{
	  if ( j == k ) continue;
	  if ( R(k,kt[j]) > spatialr[k] )
	    {
	      optimal = false;
	      break;
	    }
	}

      writer.level( std::string( 1, Tl[ best[k] ] ) , "KT" );
      writer.value( "K1" , std::string( 1, Al[k] ) );
      writer.value( "SPC" , spatialr[k] );
      writer.value( "FLIP" , flip[k] );
      writer.value( "MAPPED" , okay );
      writer.value( "OPTIMAL" , optimal );

      if ( ! optimal ) all_optimal = false;
      
      // do flip?
      if ( flip[k] )
	A->A.col(k) = -1 * A->A.col(k).array();
      
      // copy label
      r[k] = newlab;

    }

  writer.unlevel( "KT" );

  writer.value( "S" , min_res );
  writer.value( "OPTIMAL" , all_optimal );
  writer.value( "MAPPED" , all_mapped );

  return r;
  
}



