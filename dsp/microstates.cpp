
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

  
  //
  // Which channels
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
  // Fetch sample matrix
  //
  
  matslice_t mslice( edf , signals , edf.timeline.wholetrace() );

  const Data::Matrix<double> & X = mslice.data_ref();
  

  //
  // Perform microstate segmentation 
  //
  
  microstates_t mstates( param , sr );      

  //
  // Find peaks?
  //

  std::vector<int> peaks;
  
  if ( single_sample || multi_peaks )
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
      const std::string filename = Helper::expand( param.value( "backfit" ) );
      prototypes.read( filename );

      // check channels line up 
      if ( signals.size() != prototypes.C )
	Helper::halt( "number of channels in " + filename + " does not match current signal selection" );
      for (int s=0; s<signals.size(); s++)
	if ( ! Helper::iequals( signals.label(s) , prototypes.chs[s] ) )
	  Helper::halt( signals.label(s) + " does not match " + prototypes.chs[s]  + " for slot " + Helper::int2str( s+1) );
    }


  //
  // Backfitting
  //

  const bool store_GMD = true; // not sure this is needed...    

  ms_backfit_t bf = mstates.backfit( Statistics::transpose( X )  , prototypes.A , store_GMD );
  
  
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
  
  ms_stats_t stats = mstates.stats( Statistics::transpose( X ) , prototypes.A , smoothed.best() );


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


  // State transition probabilities
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

  // State frequencies


  
  std::map<int,std::pair<int,double> >::const_iterator cc = cnts.begin();
  while ( cc != cnts.end() )
    {
      std::cout << cc->first << "\t"
		<< cc->second.first << "\t"
		<< cc->second.second << "\n";
      ++cc;
    }
  
  
  
}


microstates_t::microstates_t( param_t & param , const int sr_ )
{

  //
  // Microstate analysis parameters
  //

  sr = sr_;
  
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
  
  dump_file = param.has( "dump" ) ? param.value( "dump" ) : "";
  
  standardize = param.has( "standardize" );
  
  verbose = param.has( "verbose" );
  
  // reject peaks > T times std(GFP) above the mean GFP if T>0 
  gfp_threshold = param.has( "gfp-th" ) ? param.requires_dbl( "gfp-th" ) : 0 ;
  
  // if > 0 , select (randomly) only this many peaks per observation
  restrict_npeaks = param.has( "npeaks" ) ? param.requires_int( "npeaks" ) : 0; 

  // not imlpemented yet
  // select only peaks this far apart (to ensure distinct peaks)
  min_peak_dist = 0;    
  
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
  // Optionally, dump input prior to clustering?
  //

  if ( dump_file != "" )
    {
      logger << "  dumping raw matrix to " << dump_file << "\n";
      std::ofstream O1( dump_file.c_str() , std::ios::out );      
      O1 << Z.dump();       
      O1.close();      
    }
    

  //
  // Modified K-Means clustering for microstates
  //

  modkmeans_t kmeans( ks );

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


  Data::Matrix<double> X = X_;
  Data::Matrix<double> A = A_;
  
  // X will be C x N
  // A will be C x K
  
  // polarity invariant back-fitting

  const int C = A.dim1();
  const int K = A.dim2();
  const int N = X.dim2(); // assumes X is already transposed as C x N 
  
  //
  // Standardize EEG first?  hmm check polarity here, etc
  //
  
  if ( 0 )
    Statistics::standardize( X );

  // GMD: global map dissimilarity
  
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


ms_backfit_t smooth_windowed( const ms_backfit_t & labels ,
			      const Data::Matrix<double> & X_ ,
			      const Data::Matrix<double> & A_ ,
			      int smooth_width ,	
			      double smooth_weight ,
			      int max_iterations ,
			      double threshold )
{

  Helper::halt( "not yet implemented" );
  
  Data::Matrix<double> X = X_;  // C x N  EEG 
  Data::Matrix<double> A = A_;  // C x K  prototypes

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

  // Total GEV (nb. bsed on un-normalized X_)

  Data::Vector<double> var = Statistics::sdev( X_ ,Statistics::mean( X_ ) );
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
  // k-mer distributions
  //

  ms_kmer_t kmers( runs.d , 2 , 6 );

   
  return stats;
}



ms_kmer_t::ms_kmer_t( const std::vector<int> & l , int k1 , int k2 )
{

  // ensure bounded if using this naive algorithm
  if ( k1 < 2 ) k1 = 2;
  if ( k2 > 10 ) k2 = 8;

  // Encode as A, B, C, ...
  const int n = l.size();
  s = std::string( n , '?' );
  for (int i=0; i<n; i++) s[i] = (char)(65 + l[i] ); 
       
  for (int i=0; i<n; i++)
    for (int j=k1; j<=k2; j++)
      if ( i+j < n ) ++obs[ s.substr( i , j ) ];

  
  // std::map<std::string,int>::const_iterator cc = obs.begin();
  // while ( cc != obs.end() )
  //   {
  //     std::cout << cc->first << "\t" << cc->second << "\n";
      
  //     std::vector<std::string> perms = permute( cc->first );
  //     std::cout << "  " << perms.size() << " perms: ";
  //     for (int p=0; p<perms.size(); p++)
  // 	std::cout << " " << perms[p] ;
  //     std::cout << "\n";
      
  //     ++cc;
  //   }

  //  std::cout << "found " << obs.size() << " distinct kmers\n";

  // https://math.stackexchange.com/questions/220521/expected-number-of-times-random-substring-occurs-inside-of-larger-random-string
  // https://www.nature.com/articles/s41598-018-33433-8

  // compare to max of each permutation group
  //  i.e. if ABCD  then all permutations of A,B,C and D

  
  
}
 
// Function to print permutations of string str,
// out is used to store permutations one by one
std::vector<std::string> ms_kmer_t::permute( std::string str )
{
  std::vector<std::string> perms;
  if ( str.size() == 0) return perms;
  std::sort( str.begin() , str.end() );
  do { perms.push_back( str ); }
  while ( next_permutation( str.begin() , str.end() ) );
  return perms;
}
 


void ms_kmer_t::shuffle( const int rep )
{


  if ( obs.size() == 0 ) return;

  // use 'rep' random replicates, fully shuffling the sequences, to
  // generate the null distribution of kmers expected by chance and
  // compare to what we observed

  for (int r = 0 ; r < rep ; r++)
    {

      std::vector<int> a;
    }
  
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

  A.clear();
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
      
      logger << "  writing an aggregated GPF-peak SEDF to " << edfname << "\n";
      
      edf_t agg_edf;
      
      const int recdur = 1;
      const int nr = N;
      const int ns = C;
      
      //
      // Set header
      //
      
      agg_edf.header.version = "0";
      agg_edf.header.patient_id = "GFP-peaks";
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

