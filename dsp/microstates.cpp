
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

#include <limits>

extern writer_t writer;
extern logger_t logger;

void dsptools::microstates( edf_t & edf , param_t & param )
{

  const bool no_annotations = true;
  
  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) , no_annotations );  

  const int ns = signals.size();

  // Check sample rates
  
  if ( ns < 2 ) return;
  
  const int sr = edf.header.sampling_freq( signals(0) );
  
  for (int i=1;i<ns;i++)
    {      
      if ( edf.header.sampling_freq( signals(i) ) != sr ) 
	Helper::halt( "all signals must have similar SR for MICROSTATES" );
    }


  
  //
  // Fetch sample matrix
  //
  
  matslice_t mslice( edf , signals , edf.timeline.wholetrace() );

  const Data::Matrix<double> & X = mslice.data_ref();
  

  //
  // Segment, or load an existing solution from disk?
  //

  bool segment = false;

  
  
  //
  // Perform basic microstate analysis
  //

  if ( 0 )
    {
      microstates_t mstates( param );
      
      mstates.segment( X , signals );
    }
  
  
  //
  // Get solutions
  //

  if ( 1 )
    {
      
      Data::Matrix<double> A(105,4); 
      std::ifstream IN1( "sol.4" , std::ios::in );
      std::vector<std::string> channels;
      int row = 0;
      while ( !IN1.eof() )
	{
	  std::string ch;
	  double a1,a2,a3,a4;

	  IN1 >> ch;
	  std::cout << "reading " << ch << "\n";
	  if ( ch == "" || IN1.eof() ) continue;
	  IN1 >> a1 >> a2 >> a3 >> a4;
	  channels.push_back( ch );
	  A(row,0) = a1;
	  A(row,1) = a2;
	  A(row,2) = a3;
	  A(row,3) = a4;	  
	  ++row;
	}
      IN1.close();

      // check that signals match

      if ( ns != A.dim1() ) Helper::halt( "different number of signals" );
      for (int s=0; s<ns; s++)
	if ( ! Helper::iequals( signals.label(s) , channels[s] ) )
	  Helper::halt( "signals do not match" );
    
      // assumes a transposed X

      microstates_t mstates( param );

      std::cout << "abnout to BF\n";

      bool store_GMD = true;
      
      ms_backfit_t bf = mstates.backfit( Statistics::transpose( X )  , A , store_GMD );
      

      //
      // Smooth
      //

      // smooth_reject takes minTime in sample points

      double minTime_msec = param.has( "min-msec" ) ? param.requires_dbl( "min-msec" ) : 20 ; 
      int minTime_samples = round( minTime_msec * sr/1000.0 );

      // use minimum segment length smoothing

      logger << "  smoothing: rejecting segments <= " << minTime_msec << " msec\n";
      ms_backfit_t smoothed = mstates.smooth_reject( bf , minTime_samples );


      //
      // print stats post smoothing
      //

      std::map<int,std::pair<int,double> > cnts = microstates_t::counts( smoothed.best() );
      
      std::map<int,std::pair<int,double> >::const_iterator cc = cnts.begin();
      while ( cc != cnts.end() )
	{
	  std::cout << cc->first << "\t"
		    << cc->second.first << "\t"
		    << cc->second.second << "\n";
	  ++cc;
        }


      //
      // Final stats
      //

      ms_stats_t stats = mstates.stats( Statistics::transpose( X ), A, smoothed.best(), sr );

    }
  
  
  
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

  std::cout << GMD.print( "GMD" , 3 , 10 ) << "\n";  

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


  if ( 1 )
    {
      std::cout << "N = " << N << " " << bf.labels.size() << "\n";
      
      std::map<int,std::pair<int,double> > cnts = microstates_t::counts( bf.best() );
      
      std::map<int,std::pair<int,double> >::const_iterator cc = cnts.begin();
      while ( cc != cnts.end() )
	{
      std::cout << cc->first << "\t"
		<< cc->second.first << "\t"
		<< cc->second.second << "\n";
      ++cc;
	}
    }
  
  return bf;
    
}
  
			     




microstates_t::microstates_t( param_t & param )
{
  //
  // Microstate analysis parameters
  //
  
  ks = param.intvector( "k" );
  
  dump_file = param.has( "dump" ) ? param.value( "dump" ) : "";
  
  take_abs = false;
  
  standardize = false;

  verbose = true;

  gfp_threshold = 0; // reject peaks > T times std(GFP) above the mean GFP if T>0 

  npeaks = 0; // if > 0 , select (randomly) only this many peaks per observation

  // not imlpemented yet
  min_peak_dist = 0; // e.g. select only peaks this far apart (to ensure distinct peaks)
  
  
}


void microstates_t::segment( const Data::Matrix<double> & X , 
			     const signal_list_t & signals )
{

  //
  // Copy inputs
  //

  Data::Matrix<double> Z = X;
 

  //
  // Standardize values?
  //

  if ( standardize )
    Statistics::standardize( Z );
  
  //
  // Global field power 
  //

  const int np = X.dim1();
  const int nc = X.dim2();

  logger << "  calculating GFP for sample\n";

  Data::Vector<double> GFP( np );

  for (int i=0; i<np; i++)
    {

      // get time-points across channels
      Data::Vector<double> p = Z.row( i );
    
      // get SD of raw data
      GFP[i] = sqrt( Statistics::variance( p , 0 ) ); // use N denom
      
    }
  
  //
  // Find peaks in GFP
  //

  gfp_threshold = -1 ;  // tmp KLUDGE
  
  bool find_peaks = gfp_threshold >= 0 ;

  std::vector<bool> peak( np , false );
  std::vector<int> peak_idx;
  int n_peaks = 0;

  if ( find_peaks )
    {
      for (int i=1; i<(np-1); i++)
	{
	  if ( GFP[i] > GFP[i-1] && GFP[i] > GFP[i+1] ) 
	    {
	      peak[i] = true;
	      peak_idx.push_back(i);
	      ++n_peaks;
	    }
	}
    }
  else // just copy all data 
    {
      for (int i=0; i<np; i++)
	{
	  peak[i] = true;
	  peak_idx.push_back(i);
	  ++n_peaks;
	}
    }
  
      
  //
  // Output GFP
  //
  
  if ( find_peaks && verbose )
    {
      for (int i=0; i<sol.size(); i++)
	{
	  writer.level( peak_idx[i] , "SP" );
	      writer.value( "GFP" , GFP[i] );
	}
      writer.unlevel( "SP" );
    }
  
  //
  // Copy subset of data for clustering (nb. by default, ignores signal polarity for clustring)
  //
  
  Data::Matrix<double> P( n_peaks , nc );
  
  for (int r=0; r<n_peaks; r++)
    for (int c=0; c<nc; c++)
      P( r , c ) = take_abs ? fabs( Z( peak_idx[r] , c ) ) : Z( peak_idx[r] , c ) ;
  
  
  logger << "  extracted " << n_peaks << " peaks from " << np << " samples ("
	 << round( 100 * ( n_peaks / (double)np ) ) << "%)\n";
  

  
  //
  // Optionally, dump input prior to clustering (to check K-means)
  //

  if ( dump_file != "" )
    {
      logger << "  dumping raw matrix to " << dump_file << "\n";
      std::ofstream O1( dump_file.c_str() , std::ios::out );      
      O1 << P.dump();       
      O1.close();      
    }


  //
  // Modified K-Means
  //

  modkmeans_t kmeans( ks );

  modkmeans_all_out_t results = kmeans.fit( P );

  std::cout << "got results back here\n";

  
  //
  // summarize class frequencies
  //


  const int C = P.dim2();
  const int N = P.dim1();
  

  //
  // optimal K selected
  //
  
  writer.value( "KN" , results.K );


  //
  // Maps
  //
  
  for (int i=0; i<C; i++)
    {
      writer.level( signals.label(i) , globals::signal_strat );

      // optimal solution
      for (int j=0; j<results.K; j++)
	{
	  writer.level( j+1 , "K" );
	  writer.value( "A" , results.A(i,j) );
	}
      writer.unlevel( "K" );	    
    }
  writer.unlevel( globals::signal_strat );
  
  
  
  //
  // Assignments
  //

  std::map<int,int> cnts;
  for (int i=0; i< results.L.size(); i++)
    cnts[ results.L[i] ]++;
  
  // optimal solution                                                                                                                                                               
  for (int j=0; j<results.K; j++)
    {
      writer.level( j+1 , "K" );
      writer.value( "N" , cnts[ j ] );
      writer.value( "F" , cnts[ j ] / (double)N );
    }
  writer.unlevel( "K" );
  
  
  //
  // Detailed outputs
  //

  for (int ki=0; ki<ks.size(); ki++)
    {
      const int K = ks[ki];      
      writer.level( K , "KN" );
      writer.value( "MSE" , results.kres[K].MSE );
      writer.value( "R2" , results.kres[K].R2 );
      writer.value( "MSE" , results.kres[K].MSE );
      writer.value( "SIG2" , results.kres[K].sig2 );
      writer.value( "SIG2_MCV" , results.kres[K].sig2_modk_mcv );
    }
  writer.unlevel( "KN" );
  
  return;
  
  //
  // K-means clustering on peaks
  //

  for (int ki = 0 ; ki < ks.size(); ki ++ )
    {

      //
      // set number of clusters
      //
      
      const int k = ks[ki];
      
      writer.level( k , "NK" );

      //
      // do clustering 
      //

      kmeans_t kmeans;
      std::vector<int> sol;

      // returns channels x classes matrix of means
      // sol contains observation -> class assignments
      
      Data::Matrix<double> means = kmeans.kmeans( P , k , &sol );

      // Get unit scaled vesion of means (for plotting)

      Data::Matrix<double> means01 = Statistics::unit_scale_cols( means );

      // summarize class frequencies
      std::map<int,int> cnts;
      for (int i=0; i<sol.size(); i++) cnts[ sol[i] ]++;
           

      //
      // Output means
      //
      
      for (int i=0; i<k; i++)
	{
	  writer.level( i , "KI" );
	  for (int s=0;s<signals.size();s++)	    
	    {
	      writer.level( signals.label(s) , globals::signal_strat );
	      writer.value( "M" , means( s , i ) );
	      writer.value( "M01" , means01( s , i ) );
	    }
	  writer.unlevel( globals::signal_strat );
	}
      writer.unlevel( "KI" );


      //
      // Output solution
      //

      for (int i=0; i<sol.size(); i++)
	{
	  writer.level( peak_idx[i] , "SP" );
	  writer.value( "S" , sol[i] );
	  writer.value( "GFP" , GFP[i] );
	}
      writer.unlevel( "SP" );

      //
      // Summary counts      
      //

      std::map<int,int>::const_iterator cc = cnts.begin();
      while ( cc != cnts.end() )
	{
	  writer.level( cc->first , "KI" );
	  writer.value( "N" , cc->second );
	  writer.value( "PCT" , cc->second / (double)n_peaks );
	  
	  ++cc;
	}
      writer.unlevel( "KI" );


    } // next K 
  
  writer.unlevel( "K" );
   
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
				 const std::vector<int> & L ,
				 const int sr )
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


  std::cout << X.print( "X NORM" , 10 , 5 ) << "\n";

  std::cout << A.print( "A NORM" , 10 , 4 ) << "\n";
  
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

  std::cout << GMD.print( "GMD" , 3 , 10 ) << "\n";  

  // Spatial correlation (prototypes/templates and microstates)
  std::cout << " K, N , C = " << K << " " << N << " " << C << "\n";
  
  Data::Matrix<double> SpatCorr( K , N );
  for (int i=0;i<K;i++)
    for (int j=0;j<N;j++)
      SpatCorr(i,j) = 1 - ( GMD(i,j) * GMD(i,j) ) / 2.0;

  std::cout << SpatCorr.print( "SpatCorr" , 4 , 5 ) << "\n";
  
  // Total GEV (nb. bsed on un-normalized X_)

  Data::Vector<double> var = Statistics::sdev( X_ ,Statistics::mean( X_ ) );
  double denom = 0;
  for (int j=0;j<N;j++)
    {
      var[j] *= var[j];
      denom += var[j];
    }

  double GEV_tot = 0;
  for (int j=0;j<N;j++)
    GEV_tot += SpatCorr( L[j] , j ) * var[j] ;
  GEV_tot /= denom;

  std::cout << "GEV_tot = " << GEV_tot << "\n";
  std::cout << "denom = " << denom << "\n";
  
  //   .Gfp        - Global field power
  //   .Occurence  - Occurence of a microstate per s
  //   .Duration   - Average duration of a microstate
  //   .Coverage   - % of time occupied by a microstate
  //   .GEV        - Global Explained Variance of microstate
  //   .MspatCorr  - Spatial Correlation between template maps and microstates
  //   .TP         - transition probabilities
  //   .seq        - Sequence of reoccurrence of MS ((trials x)  time).
  //   .msFirstTF  - First occurence of a microstate (similar to use like seq)
  //                 ((trials x)  time)

  // transition probabilities

  ms_rle_t runs = rle( L );

  //  stats.GFP = GFP;
  
  // means
  Data::Vector<double> m_gfp(K);
  Data::Vector<double> m_dur(K);
  Data::Vector<double> m_occ(K);
  Data::Vector<double> m_cov(K);
  Data::Vector<double> m_gev(K);
  Data::Vector<double> m_spc(K);

  for (int k=0; k<K; k++)
    {

      // Mean GFP (nb. uses /n-1 version here only to match ML Toolbox
      // but we will swap to /n

      std::vector<double> gfp_k;
      for (int j=0;j<N;j++)
	if ( L[j] == k )
	  gfp_k.push_back( GFP_minus1[j] );
      m_gfp[k] = MiscMath::mean( gfp_k );
      
      // occur/duration
      std::vector<double> times;
      for (int i=0; i<runs.d.size(); i++)
	if ( runs.d[i] == k ) times.push_back( runs.c[i] * ( 1000.0 / sr ) );

      m_occ[k] = times.size() / (double) N * sr;
      m_dur[k] = MiscMath::mean( times );
      m_cov[k] = ( m_occ[k] * m_dur[k] ) / 1000.0;

      // mean spatial correl
      std::vector<double> spc_k;
      for (int j=0;j<N;j++)
        if ( L[j] == k )
          spc_k.push_back( SpatCorr( L[j], j ) );
      m_spc[k] = MiscMath::mean( spc_k );
      
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
      
      m_gev[k] = numer / denom;

      // next K
    }
  
  std::cout << m_gfp.print( "m_gfp" ) << "\n";
  std::cout << m_occ.print( "m_occ" ) << "\n";
  std::cout << m_dur.print( "m_dur" ) << "\n";
  std::cout << m_cov.print( "m_cov" ) << "\n";

  std::cout << m_spc.print( "m_spc" ) << "\n";
  std::cout << m_gev.print( "m_gev" ) << "\n";
  
  //
  // transition probs
  //

  // runs.d contains sequence of states

  const int seqlen = runs.d.size();

  Data::Matrix<double> tr(K,K);
  Data::Vector<double> row(K);
  
  for (int s = 0 ; s < seqlen - 1 ; s++)
    {
      ++tr( runs.d[s] , runs.d[s+1] );
      ++row( runs.d[s] );
    }

  for (int i=0;i<K;i++)
    for (int j=0;j<K;j++)
      if ( i != j ) tr(i,j) /= row[i];

  std::cout << tr.print( "TR" ) << "\n";

  //
  // LZW complexity
  //

  double ratio = 0;

  lzw_t lzw( L , &ratio );
  std::cout << "ratio = " << ratio << "\n";


  lzw_t lzw2( runs.d , &ratio );
  std::cout << "ratio (2) = " << ratio << "\n";

  return stats;
}


