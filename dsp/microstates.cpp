
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
#include "miscmath/crandom.h"
#include <limits>

extern writer_t writer;
extern logger_t logger;


//
// TODO:  to add annot from states;
//        to add peaks to a cache_t from states
//        to specify labels directly for a prototype file (i.e. if diff order)
//        option to make aggregare EDF out of prototypes (i.e. to subsequently cluster prototypes)

// static member to track header of aggrgated EDF in multi-sample mode
bool microstates_t::wrote_header = false;

void dsptools::microstates( edf_t & edf , param_t & param )
{

  //
  // This function can be called in one of several modes
  //
  // Single-EDF mode: find peaks, segmentat, backfit, smooth then calculate stats
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
  if ( epoch && param.has( "gfp" ) )
    Helper::halt( "cannot specify epoch and gfp (to dump sample-level GFP) together" );
  if ( epoch && param.has( "write-states" ) )
    Helper::halt( "cannot specify epoch and write-states (to dump state order) together" );
  
  int run_kmers = param.has( "kmers" ) ? param.requires_int( "kmers" ) : 0 ;
  
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
      const std::string filename = Helper::expand( param.value( "backfit" ) );
      prior_prototypes.read( filename );
      
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
	  logger << "  segmenting peaks to microstrates\n";
	  
	  // nb, if peaks is empty, just takes all rows
	  // of X;  i.e. if coming from an aggregate peak EDF
	  
	  prototypes = mstates.segment( X , signals , peaks );
	  
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
	     prototypes = prior_prototypes;
	   }
	 else
	   {
	     const std::string filename = Helper::expand( param.value( "backfit" ) );
	     prototypes.read( filename );
	     
	     // check channels line up 
	     if ( signals.size() != prototypes.C )
	       Helper::halt( "number of channels in " + filename + " does not match current signal selection" );
	     for (int s=0; s<signals.size(); s++)
	       if ( ! Helper::iequals( signals.label(s) , prototypes.chs[s] ) )
		 Helper::halt( signals.label(s) + " does not match " + prototypes.chs[s]  + " for slot " + Helper::int2str( s+1) );
	   }
       }

      
      //
      // Backfitting
      //
      
      const bool store_GMD = true; // not sure this is needed...    
      
      logger << "  back-fitting solution to all time points\n";
      
      ms_backfit_t bf = mstates.backfit( Statistics::transpose( X ) , microstates_t::eig2mat( prototypes.A ) , store_GMD );
      
  
      //
      // Smoothing / rejection of small intervals
      //
      
      // smooth_reject takes minTime in sample points
      
      double minTime_msec = param.has( "min-msec" ) ? param.requires_dbl( "min-msec" ) : 20 ; 
      
      int minTime_samples = round( minTime_msec * sr/1000.0 );
      
      logger << "  smoothing: rejecting segments <= " << minTime_msec << " msec\n";
      
      ms_backfit_t smoothed = mstates.smooth_reject( bf , minTime_samples );
      
  
      //
      // Final stats
      //
      
      ms_stats_t stats = mstates.stats( Statistics::transpose( X ) , microstates_t::eig2mat( prototypes.A ) , smoothed.best() );
  


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
	    O1 << (char)(states[i] + 65 ) << "\t"
	       << stats.GFP[i] << "\n";
	  O1.close();
	}
  

      //
      // Output all final stats
      //

      
      //
      // By class 'K'
      //
      
      std::map<int,std::pair<int,double> > cnts = microstates_t::counts( smoothed.best() );
      
      for (int k=0; k < prototypes.K ; k++)
	{
	  std::string s="?";
	  s[0] = (char)(65+k);
	  writer.level( s , "K" );
	  writer.value( "GFP" , stats.m_gfp[k] );
	  writer.value( "OCC" , stats.m_occ[k] );
	  writer.value( "DUR" , stats.m_dur[k] );
	  writer.value( "COV" , stats.m_cov[k] );
	  writer.value( "SPC" , stats.m_spc[k] );
	  writer.value( "GEV" , stats.m_gev[k] );      
	  writer.value( "N" , cnts[k].first );
	  writer.value( "F" , cnts[k].second );
	}
      writer.unlevel( "K" );
      

  
      //
      // State transition probabilities
      //
      
      for (int k=0; k < prototypes.K ; k++)
	{
	  std::string s1="?";
	  s1[0] = (char)(65+k);
	  writer.level( s1 , "PRE" );
	  for (int k2=0; k2 < prototypes.K ; k2++)
	    {
	      if ( k != k2 )
		{
		  std::string s2="?";
		  s2[0] = (char)(65+k2);
		  writer.level( s2 , "POST" );
		  writer.value( "P" , stats.tr(k,k2) );
		}
	    }
	}
      writer.unlevel( "PRE" );
      writer.unlevel( "POST" );
      
      
      // Overall stats: complexity
      
      writer.value( "LZW"     , stats.lwz_states );
      writer.value( "LZW_ALL" , stats.lwz_points );
      writer.value( "GEV"     , stats.GEV_tot );
      
      // kmer stats: single obs, so no group differences here
    
      if ( run_kmers )
	{
	  std::map<std::string,double>::const_iterator pp = stats.kmers.basic.pval.begin();
	  while ( pp != stats.kmers.basic.pval.end() )
	    {
	      writer.level( pp->first , "KMER" );
	      writer.level( (int)pp->first.size() , "L" );
	      
	      writer.value( "EQ"   , stats.kmers.obs2equiv[ pp->first ] );
	      writer.value( "EQN"  , stats.kmers.equiv_set_size[ pp->first ] );
	      
	      writer.value( "O_OBS" , stats.kmers.basic.obs[ pp->first ] );
	      writer.value( "O_EXP" , stats.kmers.basic.exp[ pp->first ] );
	      writer.value( "O_P" , pp->second );
	      writer.value( "O_Z" , stats.kmers.basic.zscr[ pp->first ] );

	      if ( stats.kmers.equiv_set_size[ pp->first ] > 1 )
		{
		  writer.value( "M_OBS" , stats.kmers.equiv.obs[ pp->first ] );
		  writer.value( "M_EXP" , stats.kmers.equiv.exp[ pp->first ] );
		  writer.value( "M_P" , stats.kmers.equiv.pval[ pp->first ] );
		  writer.value( "M_Z" , stats.kmers.equiv.zscr[ pp->first ] );
		}
	      
	      ++pp;
	    }  
	  writer.unlevel( "L" );
	  writer.unlevel( "KMER" );
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
  
  // reject peaks > T times std(GFP) above the mean GFP if T>0 
  gfp_threshold = param.has( "gfp-th" ) ? param.requires_dbl( "gfp-th" ) : 0 ;
  
  // if > 0 , select (randomly) only this many peaks per observation
  restrict_npeaks = param.has( "npeaks" ) ? param.requires_int( "npeaks" ) : 0; 

  // not imlpemented yet
  // select only peaks this far apart (to ensure distinct peaks)
  min_peak_dist = 0;    

  if ( param.has( "kmers" ) )
    {
      std::vector<int> k = param.intvector( "kmers" );
      if ( k.size() != 3 ) Helper::halt( "expecting 3 args for kmers=min,max,nreps" );
      kmers_min = k[0];
      kmers_max = k[1];
      kmers_nreps = k[2];
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
  // GFP threshold
  //

  if ( gfp_threshold > 0 )
    {
      Data::Vector<double> peak_gfp( n_peaks );
      for (int r=0; r<n_peaks; r++)
	peak_gfp[r] = GFP[ peak_idx[r] ];
      
      double mean = Statistics::mean( peak_gfp );
      double sd = sqrt( Statistics::variance( peak_gfp , 1 ) ); // use N-1 denom
      double th = mean + gfp_threshold * sd;

      std::vector<int> peak_idx2;
      for (int r=0; r<n_peaks; r++)
	if ( GFP[ peak_idx[r] ] <= th )
	  peak_idx2.push_back( peak_idx[r] );

      logger << "  applying GFP threshold mean + " << gfp_threshold << "SDs, keeping " << peak_idx2.size() << " of " << n_peaks << " peaks\n";
      
      n_peaks = peak_idx2.size();
      peak_idx = peak_idx2;
      
    }

  //
  // Only select (at most) N peaks?
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
					const std::vector<int> & peak_idx )
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

  // modkmeans_t( const std::vector<int> & ks ,
  // 	       const bool normalize = false ,
  // 	       const int nreps = 10 ,
  // 	       const int max_iterations = 1000 ,
  // 	       const double threshold = 1e-6 ,
  // 	       const bool verbose = false )

  modkmeans_all_out_t results = kmeans.fit( Z );

  //
  // optimal K selected
  //
  
  writer.value( "OPT_K" , results.K );


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
	  s[0] = (char)(65+j);
	  writer.level( s , "K" );
	  writer.value( "A" , results.A(i,j) );
	}
      writer.unlevel( "K" );	    
      
    }
  writer.unlevel( globals::signal_strat );

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
	      s[0] = (char)(65+j);
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
      writer.value( "MSE" , results.kres[K].MSE );
      writer.value( "R2" , results.kres[K].R2 );
      writer.value( "MSE" , results.kres[K].MSE );
      writer.value( "SIG2" , results.kres[K].sig2 );
      writer.value( "SIG2_MCV" , results.kres[K].sig2_modk_mcv );
    }
  writer.unlevel( "NK" );

  //
  // Save prototypes
  //


  ms_prototypes_t prototypes( signals , results.A ) ;
  return prototypes;
  
}


ms_backfit_t microstates_t::backfit( const Data::Matrix<double> & X_ ,
				     const Data::Matrix<double> & A_ ,
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
  // Standardize EEG first?  hmm check polarity here, etc
  //
  
  // if ( 0 )
  //   Statistics::standardize( X );


  //
  // GMD: global map dissimilarity
  //
  
  // Assumes average reference already set  
  // Normalise EEG and maps (average reference and gfp = 1 for EEG)

  //X = X ./ repmat(std(X,1), C, 1); % already have average reference  
  //A = (A - repmat(mean(A,1), C, 1)) ./ repmat(std(A,1), C, 1);
  
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

  // GMD = nan(K,N*T);
  // for k = 1:K
  //   GMD(k,:) = sqrt(mean( (X - repmat(A(:,k),1,N*T)).^2 ));
  // end

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


  // get matching labels (min. GMD)

  ms_backfit_t bf(N);

  for (int j=0;j<N;j++)
    {
      // add all labels/GMDs which will be sorted by add()
      for (int k=0;k<K;k++)
	bf.labels[j].add( k , GMD(k,j) );

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
  // no need to populate GMD (esp. of labels may change in any case)
  
  ms_backfit_t bf(N);
  bf.labels = sol.labels;
    
  for (int k=1; k <= minTime; k++)
    {
      // track changes
      std::vector<int> cruns( N , k );      
      while ( 1 )
	{
	  int sum_cruns = 0;
	  for (int c=0; c<cruns.size(); c++)
	    if ( cruns[c] <= k ) ++sum_cruns;
	  if ( sum_cruns == 0 ) break;
	  //	  std::cout << " iter " << k << "\t" << sum_cruns << "\n";
	  
	  ms_rle_t runs = rle( bf.best() );
	  
	  int cnt = 0;
	  for (int r=0; r<runs.c.size(); r++)
	    for (int j=0; j<runs.c[r]; j++)
	      {
		// shift if segment is too short
		if ( runs.c[r] <= k ) bf.labels[cnt].shift();
		cruns[cnt] = runs.c[r];
		++cnt;
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

  Helper::halt( "not yet implemented" );

  // TEMP..
  Data::Matrix<double> X = eig2mat( X_ );  // C x N  EEG 
  Data::Matrix<double> A = eig2mat( A_ );  // C x K  prototypes

  const int C = X.dim1();
  const int N = X.dim2();
  const int K = A.dim2();
  
  ms_backfit_t bf(N);

  // % [L,sig2,R2,MSE,ind] = window_smoothing(X,A,opts)
  // %  Implementation of the Segmentation Smoothing Algorithm, as described in
  // %  Table II of [1]. Smoothes using the interval t-b to t+b excluding t.
  // %  Note, that temporary allocation of labels (denoted with Lambda in [1])
  // %  is not necessary in this implementation, and steps 3 and 6 are therefore
  // %  left out.

  // const = sum(sum(X.^2));
  // const double const1 = Statistics::sum( Statistics::sum_squares( X ) );


  // double sig2_old = 0;
  // double sig2 = std::numeric_limits<double>::max();
  // xxx
  
  

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
				 const std::vector<int> & L )
{
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

  Data::Matrix<double> SpatCorr( K , N );
  for (int i=0;i<K;i++)
    for (int j=0;j<N;j++)
      SpatCorr(i,j) = 1 - ( GMD(i,j) * GMD(i,j) ) / 2.0;

  // Total GEV (nb. bsed on un-normalized X_)..
  
  Data::Vector<double> var = Statistics::sdev( X_ ,Statistics::mean(  X_ ) );
  double denom = 0;
  for (int j=0;j<N;j++)
    {
      var[j] *= var[j];
      denom += var[j];
    }

  stats.GEV_tot = 0;
  for (int j=0;j<N;j++)
    stats.GEV_tot += SpatCorr( L[j] , j ) * var[j] ;
  stats.GEV_tot /= denom;
  
  //   .Gfp        - Global field power
  //   .Occurence  - Occurence of a microstate per s
  //   .Duration   - Average duration of a microstate
  //   .Coverage   - % of time occupied by a microstate
  //   .GEV        - Global Explained Variance of microstate
  //   .MspatCorr  - Spatial Correlation between template maps and microstates
  //   .TP         - transition probabilities

  // transition probabilities

  ms_rle_t runs = rle( L );

  // copy GFP (per point -- needed?)
  
  stats.GFP = GFP;
  
  // means
  stats.m_gfp.resize(K);
  stats.m_dur.resize(K);
  stats.m_occ.resize(K);
  stats.m_cov.resize(K);
  stats.m_gev.resize(K);
  stats.m_spc.resize(K);

  for (int k=0; k<K; k++)
    {

      // Mean GFP (nb. uses /n-1 version here only to match ML Toolbox
      // but we will swap to /n
      
      std::vector<double> gfp_k;
      for (int j=0;j<N;j++)
	if ( L[j] == k )
	  gfp_k.push_back( GFP_minus1[j] );
      stats.m_gfp[k] = MiscMath::mean( gfp_k );
      
      // occur/duration
      std::vector<double> times;
      for (int i=0; i<runs.d.size(); i++)
	if ( runs.d[i] == k ) times.push_back( runs.c[i] * ( 1000.0 / sr ) );

      stats.m_occ[k] = times.size() / (double) N * sr;
      stats.m_dur[k] = MiscMath::mean( times );
      stats.m_cov[k] = ( stats.m_occ[k] * stats.m_dur[k] ) / 1000.0;

      // mean spatial correl
      std::vector<double> spc_k;
      for (int j=0;j<N;j++)
        if ( L[j] == k )
          spc_k.push_back( SpatCorr( L[j], j ) );
      stats.m_spc[k] = MiscMath::mean( spc_k );
      
      // GEV

      //GEV(trial,k) = sum( (GFP(trial,L(trial,:)==k) .* MspatCorrTMP(k,L(trial,:)==k)).^2) ./ sum(GFP(trial,:).^2);

      double numer = 0;
      double denom = 0;
      for (int j=0;j<N;j++)
	{
	  if ( L[j] == k )
	    numer += ( SpatCorr( L[j], j ) * GFP[j] ) * ( SpatCorr( L[j], j ) * GFP[j] ) ;
	  denom += GFP[j] * GFP[j];
	}
      
      stats.m_gev[k] = numer / denom;

      // next K
    }
  
  
  //
  // transition probs
  //

  // runs.d contains sequence of states

  const int seqlen = runs.d.size();

  stats.tr.resize(K,K);
  Data::Vector<double> row(K);
  
  for (int s = 0 ; s < seqlen - 1 ; s++)
    {
      ++stats.tr( runs.d[s] , runs.d[s+1] );
      ++row( runs.d[s] );
    }

  for (int i=0;i<K;i++)
    for (int j=0;j<K;j++)
      if ( i != j ) stats.tr(i,j) /= row[i];


  //
  // LZW complexity
  //

  lzw_t lzw( L , &stats.lwz_points );

  lzw_t lzw2( runs.d , &stats.lwz_states );


  //
  // dump sequences to file?
  //

  if ( statesfile != ""  )
    {      
      logger << "  writing sequence order to " << statesfile << "\n";
      const int n = runs.d.size();      
      // Encode as A, B, C, ...
      std::string s = std::string( n , '?' );
      for (int i=0; i<n; i++)
        s[ i ] = (char)(65 + runs.d[i] );
      std::ofstream OUT1( Helper::expand( statesfile ).c_str() , std::ios::out );
      OUT1 << subj_id << "\t" << s << "\n";
      OUT1.close();      
    }

  
  //
  // k-mer distributions, optionally (in single-obs mode)
  //

  if ( kmers_nreps )
    stats.kmers.run( runs.d , kmers_min , kmers_max , kmers_nreps );
   
  return stats;
}



void ms_kmer_t::run( const std::map<std::string,std::vector<int> > & lall , int k1 , int k2 , int nreps ,
		     const std::map<std::string,int> * grp )
{
  std::map<std::string,std::string> sall;
  std::map<std::string,std::vector<int> >::const_iterator ii = lall.begin();
  while ( ii != lall.end() )
    {
      const std::vector<int> & l = ii->second;
      std::string & s = sall[ ii->first ];
      const int n = l.size();
      // Encode as A, B, C, ...                                                                                                                                              
      s = std::string( n , '?' );
      for (int i=0; i<n; i++)
        s[ i ] = (char)(65 + l[i] );
      ++ii;
    }
  run( sall , k1 , k2 , nreps , grp );
}

void ms_kmer_t::run( const std::map<std::string,std::string> & sall , int k1 , int k2 , int nreps ,
		     const std::map<std::string,int> * grp )
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
  // Observed data: form equivalance groups
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

  // count up for each equivalence group
  std::map<std::string,std::string>::const_iterator ee = obs2equiv.begin();

  while ( ee != obs2equiv.end() )
    {
      // based just on the keys, get the corresponding equiv. group
      equivs[ ee->second ] = permute( ee->second );

      // for convenienc of output, track size of equiv set directly
      equiv_set_size[ ee->first ] = equivs[ ee->second ].size();

      // denominator for : obs / obs-group 
      const std::set<std::string> & st = equivs[ ee->second ];
      std::set<std::string>::const_iterator qq = st.begin();
      while ( qq != st.end() )
	{	  
	  equiv.obs[ ee->first ] += basic.obs[ *qq ];
	  if ( grp != NULL )
	    {
	      equiv_cases.obs[ ee->first ] += basic_cases.obs[ *qq ];
	      equiv_controls.obs[ ee->first ] += basic_controls.obs[ *qq ];	      
	    }
	  ++qq;
	}

      // express as proportion of total (equiv.obs is set to sum of equiv groups above)

      equiv.obs[ ee->first ] = basic.obs[ ee->first ] / equiv.obs[ ee->first ];
      
      if ( grp != NULL )
	{
	  equiv_cases.obs[ ee->first ] = basic_cases.obs[ ee->first ] / equiv_cases.obs[ ee->first ];
	  equiv_controls.obs[ ee->first ] = basic_controls.obs[ ee->first ] / equiv_controls.obs[ ee->first ];
	  // diffs based on difference in these two relative freqs
	  equiv_diffs.obs[ ee->first ] = equiv_cases.obs[ ee->first ] - equiv_controls.obs[ ee->first ];
	}
	      
      ++ee;
    }
  
  logger << "  kmers: considering length " << k1 << " to " << k2 << "\n";
  
  logger << "  kmers: for " << basic.obs.size() << " sequences, "
	 << equivs.size() << " equivalence groups\n";
	
  

  
  //
  // Begin permutations 
  //
    

  logger << "  kmers: running " << nreps << " replicates...\n";
  

  for (int r = 0 ; r < nreps ; r++)
    {
      
      //
      // Permutation stratified by individual, although all counting is
      // pooled across individuals
      //
      
      // rather than standard CRandom::random_draw(), this ensures
      // no (*) similar states are contiguous
      // ... actually, there may be one or two contiguous sequences
      // at the end, but in the big picture, this should not matter 

      std::map<std::string,int> stat_basic, stat_basic_cases, stat_basic_controls;
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
	  
	  std::string ps = modified_random_draw( ii->second );	  

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
      // Normalize equiv group stats (based on total-numer / total-denom , i.e. where total is across all people)
      //

      //
      // equivalance group stats: denominators
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
	  
	  // sum across individuals
	  stat_equiv[ pp->first ] = sum;
	  if ( grp != NULL )
	    {
	      stat_equiv_controls[ pp->first ] = sum_cases;
	      stat_equiv_cases[ pp->first ] = sum_controls;
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
	  // Normalize equiv group stats (based on total-numer / total-denom , i.e. where total is across all people)
	  //

	  
	  const std::string & ess = obs2equiv[ ss->first ];
	  
	  double estat = stat_equiv[ ess ] > 0 ? stat_basic[ ss->first ]  / (double)stat_equiv[ ess ] : 0 ;
	  equiv.perm[ ss->first ].push_back( estat );
	  
	  if ( grp != NULL )
	    {
	      
	      double estat_cases = stat_equiv_cases[ ess ] > 0 ? stat_basic_cases[ ss->first ]  / (double)stat_equiv_cases[ ess ] : 0 ;
	      double estat_controls = stat_equiv_controls[ ess ] > 0 ? stat_basic_controls[ ss->first ]  / (double)stat_equiv_controls[ ess ] : 0 ;
	      double estat_diffs  = estat_cases - estat_controls;

	      equiv_cases.perm[ ss->first ].push_back( estat_cases );
	      equiv_controls.perm[ ss->first ].push_back( estat_controls );
	      equiv_diffs.perm[ ss->first ].push_back( estat_diffs );

	    }

	  // next obs sequence
	  ++ss;
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
  std::ofstream O1( filename.c_str() , std::ios::out );
  for (int c=0; c<C; c++)
    {
      O1 << chs[c];
      for (int k=0;k<K;k++)
	O1 << "\t" << A(c,k);
      O1 << "\n";
    }
  O1.close();
}

void ms_prototypes_t::read( const std::string & filename )
{
  if ( ! Helper::fileExists( filename ) )
    Helper::halt( "could not find " + filename );

  A.resize(0,0);
  chs.clear();
  bool first = true;
  std::vector<double> t;

  C = 0;
  
  std::ifstream IN1( filename.c_str() , std::ios::in );
  while ( ! IN1.eof() )
    {
      std::string line;
      Helper::safe_getline( IN1 , line );
      if ( IN1.eof() ) break;
      if ( line == "" ) break;
      std::vector<std::string> tok = Helper::parse( line );
      if ( tok.size() < 2 ) 
	Helper::halt( "problem reading prototypes from " + filename + "\n fewer than 2 tokens\n" + line );

      if ( first )
	{
	  K = tok.size() - 1;
	  first = false;
	  logger << "  assuming " << K << " classes in " << filename << "\n";
	}
      else if ( tok.size() - 1 != K )
	Helper::halt( "problem reading prototypes from " + filename + "\n col-1 != K\n" + line );

      for (int i=1;i<tok.size();i++)
	{
	  double x;
	  if ( ! Helper::str2dbl( tok[i] , &x ) )
	    Helper::halt( "problem reading prototypes from " + filename + "\n in coversion to numeric: " + tok[i] + "\n" + line );
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

  A.resize( C , K );
  int tc = 0;
  for (int c=0; c<C; c++)
    for (int k=0;k<K;k++)
      A(c,k) = t[tc++];
  
  logger << "  read " << K << "-class prototypes for " << C << " channels from " << filename << "\n";
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
  
  // possibilities:
  // 1) ! exists                     --> write header & records
  // 2)   exists && ! wrote_header   --> overwrite, as above
  // 3)   exists &&   wrote_header   --> append to end of existing EDF

  //
  // this case should not happen
  //
  
  if ( ( !exists )  && microstates_t::wrote_header )
    Helper::halt( "problem writing EDF " + edfname );


  //
  // Exract peaks
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


std::string ms_kmer_t::modified_random_draw( const std::string & l )
{

  // permute 's' chars but in such a way that similar states (i.e. values of l, 0, 1, 2, ...)
  // are not contiguous
  // first element: pick any element (e.g. of 5)
  // next elemnt: pick from the remaining 5-1 elements, i.e. not the previously selected
  // repeat, until all done

  std::map<char,int> urns;
  const int n = l.size();
  for (int i=0;i<n;i++) ++urns[l[i]];

  // return string
  std::string p(n,'?');
  
  // initiate
  p[0] = l[ CRandom::rand( n ) ];
  
  char last = p[0];
  --urns[ last ];

  for (int i=1;i<n;i++)
    {
      // skip last pick
      p[i] = pick( urns , last );
      --urns[ p[i] ];
      last = p[i];
      
      // sanity check, can remove
      if ( urns[p[i]] < 0 ) Helper::halt( "error!" );      
      
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

Eigen::Matrix3d microstates_t::mat2eig( const Data::Matrix<double> & M )
{
  const int rows = M.dim1();
  const int cols = M.dim2();
  Eigen::MatrixXd E( rows , cols );
  for (int r=0;r<rows;r++)
    for (int c=0;c<cols;c++)
      E(r,c) = M(r,c);
  return E;
}

Eigen::Matrix3d microstates_t::mat2eig_tr( const Data::Matrix<double> & M )
{
  const int rows = M.dim1();
  const int cols = M.dim2();
  Eigen::MatrixXd E( cols , rows );
  for (int r=0;r<rows;r++)
    for (int c=0;c<cols;c++)
      E(c,r) = M(r,c);
 return E;
}
