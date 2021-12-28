
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

#include "dsp/tclst.h"
#include "timeline/timeline.h"
#include "stats/statistics.h"
#include "db/db.h"
#include "helper/logger.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "dsp/hilbert.h"
#include "stats/eigen_ops.h"
#include "stats/kmeans.h"

extern writer_t writer;
extern logger_t logger;

void dsptools::tclst( edf_t & edf , param_t & param )
{
  
  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) );
  
  if ( signals.size() < 1 ) return;
  
  const int ns = signals.size();
  
  // check sample rates: must be similar across all 
  
  std::vector<double> Fs = edf.header.sampling_freq( signals );
  
  for (int s=1;s<Fs.size();s++) 
    if ( Fs[s] != Fs[0] ) 
      Helper::halt( "sample rates must be similar across signals for TCLST" );

  //
  // complex
  //

  const bool use_complex_dist = param.has( "complex" );
  
  bool use_amp = param.has( "amp" );
  bool use_phase = param.has( "phase" );
  
  if ( use_amp && use_complex_dist ) Helper::halt( "can only specify complex OR amp" );
  if ( use_phase && use_complex_dist ) Helper::halt( "can only specify complex OR phase" );

  // default:
  if ( ! ( use_complex_dist || use_amp || use_phase ) ) { use_amp = use_phase = true; }

  //
  // Verbose report for one SW? (1-based syntax)
  //

  const int verbose_interval = param.has( "report" ) ? param.requires_int( "report" ) - 1 : -1; 
  
  //
  // clustering
  //

  int k1 = 0;
  int k2 = 0;
  
  if ( param.has( "k" ) )
    k1 = k2 = param.requires_int( "k" );
  else
    {
      if ( param.has( "k1" ) ) k1 = param.requires_int( "k1" );
      if ( param.has( "k2" ) ) k2 = param.requires_int( "k2" );      
    }

  if ( k1 > k2 || k1 < 0 || k2 < 0 )
    Helper::halt( "bad specification of k/k1/k2" );

  //
  // hierarchical complete-linkage clustering
  // (if -1, use silhouette, and replace hcK w/ selected value)
  //

  int hcK = param.has( "hc" ) ? param.requires_int( "hc" ) : 0 ;
  

  //
  // specify a seed channel, i.e. against which the phase of over channels will be measured
  //

  const std::string seed = param.requires( "seed" );
  if ( signals.find( seed ) == -1 )
    Helper::halt( "seed " + seed + " not specified in sig" );  
  int seed_n = -1;
  for (int s=0; s<ns; s++)
    if ( Helper::iequals( signals.label(s) , seed ) )
      seed_n = s;
  if ( seed_n == -1 )
    Helper::halt( "internal error in tclst_t sig/seed selection" );
  
  
  //
  // get time-points: expecting a cache in sample-point units (i.e. 
  // to ensure same number of points per window; i.e. this is also why
  // all signals must have the same SR
  //
  
  //
  // get window
  //
  
  const double half_window1 = param.has( "w" ) ? param.requires_dbl( "w" ) : param.requires_dbl( "w1" );
  const double half_window2 = param.has( "w" ) ? param.requires_dbl( "w" ) : param.requires_dbl( "w2" );

  if ( half_window1 < 0 || half_window2 < 0 || half_window1 + half_window2 == 0 )
    Helper::halt( "invalid windows (w, w1 and/or w2) " );

  const int half_points1 = half_window1 * Fs[0] ;
  const int half_points2 = half_window2 * Fs[0] ;
  const int points = 1 + half_points1 + half_points2;

  logger << "  using a window of " << half_window1 + half_window2 << " seconds, "
	 << points << " sample points\n";
  
  //
  // filter-hilbert
  //
  
  const bool filter = param.has( "f-lwr" );
  
  const double f_lwr = filter ? param.requires_dbl( "f-lwr" ) : 0 ;
  const double f_upr = filter ? param.requires_dbl( "f-upr" ) : 100 ;
  if ( f_lwr >= f_upr ) Helper::halt( "f-lwr must be lower than f-upr" );
  const double fir_ripple = param.has( "ripple" ) ? param.requires_dbl( "ripple" ) : 0.01;
  const double fir_tw = param.has( "tw" ) ? param.requires_dbl( "tw" ) : 0.5;
  
  // build 
  std::vector<double> t(points,0);
  
  const double inc = 1.0/Fs[0];
  const int midp = half_points1;
  
  t[ midp ] = 0;
  int idx = midp - 1;
  int cnt = 0;

  while ( idx >= 0 )
    {      
      t[ idx ] = - inc * cnt;
      --idx;
      --cnt;
    }

  idx = midp + 1;
  cnt = 0;
  while ( cnt != half_points2 )
    {
      t[ idx ] = inc * cnt;
      ++idx;
      ++cnt;
    }
  
  // for (int w=0; w<points; w++)
  //   std::cout << " win = " << t[w] << "\n";
  
  //
  // Get seed sample-points from cache
  //
  
  std::string cache_name = param.requires( "cache" );
  
  if ( ! edf.timeline.cache.has_int( cache_name ) )
    Helper::halt( "cache not found for this individual: " + cache_name );
  
  cache_t<int> * cache = edf.timeline.cache.find_int( cache_name );
  
  //
  // Get any/all keys associated with the 'points' intrnal name cache
  //
  
  std::set<ckey_t> ckeys = cache->keys( "points" );
  
  std::set<ckey_t>::const_iterator cc = ckeys.begin();
  
  std::vector<int> starts, stops;
  
  while ( cc != ckeys.end() )
    {
      
      int scnt = 0;
      
      std::vector<int> cx = cache->fetch( *cc );
      
      //
      // number of intervals 
      //
      
      const int ni = cx.size();
      
      logger << "  found " << ni << " intervals in the cache " << cache_name << "\n";
      
      //
      // Get time-line
      //

      if ( signals.size() == 0 ) Helper::halt( "no signals specified" );
      
      slice_t slice( edf , signals(0) , edf.timeline.wholetrace() );
      const std::vector<uint64_t> * tp = slice.ptimepoints();                                                     
      
      //
      // Get intervals
      //
      
      for ( int i=0; i<cx.size(); i++)
	{
	  
	  int lower = cx[i] - half_points1;
	  int upper = cx[i] + half_points2;
	  
	  // interval out-of-range
	  if ( lower < 0 || upper >= tp->size() ) 
	    continue;
	  
	  // check for discontinuities
	  if ( edf.timeline.discontinuity( *tp , Fs[0] , lower , upper ) )
	    continue;
	  
	  // add to the list
	  starts.push_back( lower );
	  stops.push_back( upper );
	  
	}
      

      
      //
      // build up data stores
      //
      
      // peaks/intervals (ni) * time-points (points) * signals (ns) 
      //                                     complex | non-complex
      std::vector<Eigen::MatrixXd> X(ni); // real    | signal amplitudes (filtered)
      std::vector<Eigen::MatrixXd> P(ni); // imag    | phases
      std::vector<Eigen::MatrixXd> P2(ni); // store phases for output, if in complex mode
      std::vector<Eigen::MatrixXd> ZP2(ni); // store phases for output, if in complex mode

      for (int i=0; i<ni; i++)
	{
	  X[i] = Eigen::MatrixXd::Zero( points , ns );
	  P[i] = Eigen::MatrixXd::Zero( points , ns );
	  if ( use_complex_dist )
	    {
	      P2[i] = Eigen::MatrixXd::Zero( points , ns );
	      ZP2[i] = Eigen::MatrixXd::Zero( points , ns );
	    }
	}
      
      //
      // Get each signal, filter, and populate X and P
      //

      std::vector<std::string> chs;
      
      for (int s=0; s< ns; s++)
	{
	  logger << "  pulling data for " << signals.label(s) << "\n";
	  chs.push_back( signals.label(s) );
	  
	  // get all data
	  slice_t slice( edf , signals(s) , edf.timeline.wholetrace() );
	  const std::vector<double> * d = slice.pdata();	  

	  if ( ! filter ) Helper::halt( "requires f-lwr, f-upr" );

	  // filter-hilbert
	  hilbert_t ht( *d , Fs[s] , f_lwr , f_upr , fir_ripple , fir_tw , use_complex_dist );


	  //
	  // -----
	  //

	  if ( use_complex_dist )
	    {
	      
	      const std::vector<dcomp> cmp = ht.get_complex();

	      // also store phase for outputs
	      const std::vector<double> * phase = ht.phase();
	      
	      // populate each interval
	      //   X = real, P = imag
              for (int i=0; i<ni; i++)
                {
                  int s1 = starts[i];
                  for (int p=0; p<points; p++)
                    {
                      X[i](p,s) = std::real( cmp[s1] );
                      P[i](p,s) = std::imag( cmp[s1] );
		      P2[i](p,s) = (*phase)[s1];
		      ZP2[i](p,s) = (*phase)[s1];
                      ++s1;
                    }
                }	      
	    }
	  else
	    {
	      // get filtered signal
	      const std::vector<double> * signal = ht.signal();
	      
	      // get phase angles	  
	      const std::vector<double> * phase = ht.phase();
	      
	      // populate each interval
	      for (int i=0; i<ni; i++)
		{	      	      
		  int s1 = starts[i];
		  for (int p=0; p<points; p++)
		    {		  
		      //std::cout << "det " << p << " " << (*signal)[s1] << " " << (*phase)[s1] << "\n";
		      P[i](p,s) = (*phase)[s1];
		      X[i](p,s) = (*signal)[s1];
		      ++s1;
		    }	      
		}
	    }
	 
	  
	  // next signal
	}

      //
      // Save originals if verbose-reporting for one interval?
      //

      if ( verbose_interval >= ni )
	Helper::halt( "bad report=interval specified" );
      
      Eigen::MatrixXd P0, X0;
      if ( verbose_interval >= 0 )
	{
	  P0 = P[verbose_interval];
	  X0 = X[verbose_interval];	  
	}
    
      //
      // Normalize phases by the seed channel (i.e. relative phase,
      // shortest distance wrapping angles)
      //


      // phase : 0      = positive peak       
      //         +pi/2  = pos-2-neg ZC
      //         +/- pi = negative peak
      //         -pi/2  = neg-2-pos ZC

      // +ve phase difference: A-B --> A comes first 
      // -ve phase difference: A-B --> B comes first
      
      if ( ! use_complex_dist )
	{

	  //
	  // Normalize amplitudes by mean of seed at peak
	  //

	  double amp_mean = 0;
	  for (int i=0; i<ni; i++)
	    amp_mean += X[i](half_points1,seed_n);
	  amp_mean /= (double)ni;
	  
	  for (int i=0; i<ni; i++)
	    X[i] /= amp_mean;

	}

      
      //
      // Normalize phases 
      //
      
      std::vector<double> mean_phase( ni , 0 );
      std::vector<double> mean_phase_normed( ni , 0 );
      
      for (int i=0; i<ni; i++)
	{

	  // nb. given mode, point to correct lot
	  Eigen::MatrixXd & m = use_complex_dist ? P2[i] : P[i];
	  
	  const int r = m.rows();
	  const int midp = half_points1;
	  
	  std::vector<double> midp_diff( ns );
	  
	  // get difference between seed and each channel at mid-points
	  for ( int s=0; s<ns; s++)
	    {
	      
	      midp_diff[s] =
		MiscMath::deg2rad( MiscMath::angle_difference( MiscMath::rad2deg( M_PI + m(midp,seed_n) ) ,
							       MiscMath::rad2deg( M_PI + m(midp,s) ) ) );
	      
	      mean_phase[s] += m(midp,s);
	      mean_phase_normed[s] += midp_diff[s];
	      
	    }
	  
	  //
	  // unwrap each channel
	  //
	  
	  for ( int s=0; s<ns; s++)
	    {
	      std::vector<double> pp = eigen_ops::copy_array( m.col(s) ) ;
	      hilbert_t::unwrap( &pp );
	      
	      // mid-point for this channel
	      const double midph = pp[midp];
	      
	      // normalize by midpoint for this channel
	      //  but also scale relative to seed channel ( midp_diff[] )
	      for (int p=0; p<r; p++)
		pp[p] = pp[p] - midph + midp_diff[s] ;	      
	      
	      m.col(s) = eigen_ops::copy_array( pp );
	    }
	  
	  //                  XXX
	  // SEED   1    1.5   2   2.5  3    0   
	  // CH     1.1  1.6   2.1 2.6  3.2  0.2 
	  
	  //        -1   -0.5  0   0.5  1    PI         
	  //                   0.1 0.6  1.2  
	  
	  // next interval
	}
      

      // make ZPH for complex mode (output)

      if ( use_complex_dist )
	{
	  for (int i=0;i<ni;i++)
	    {
	      Eigen::VectorXd seed = ZP2[i].col(seed_n);
	      for (int s=0;s<ns;s++)
		for (int p=0;p<points;p++)
		  {
		    ZP2[i](p,s) = MiscMath::deg2rad( MiscMath::angle_difference( MiscMath::rad2deg( M_PI + ZP2[i](p,seed_n) ) ,
										 MiscMath::rad2deg( M_PI + ZP2[i](p,s) ) ) ) ;
		  }
	    }
	}
      
      // report phase means (at seed points)
      
      for ( int s=0; s<ns; s++)
	{
	  writer.level( signals.label(s) , globals::signal_strat );
	  writer.value( "PH" , mean_phase[s] / (double) ni );
	  writer.value( "ZPH" , mean_phase_normed[s] / (double) ni );
	}
      writer.unlevel( globals::signal_strat );
      
      
      //
      // verbose report?
      //

      if ( verbose_interval >= 0 )
	{
	  const Eigen::MatrixXd & XX = X[ verbose_interval ];
	  const Eigen::MatrixXd & PP = P[ verbose_interval ];
	  	  
	  for (int s=0; s<ns; s++)
	    for (int p=0; p<points; p++)
	      std::cout << chs[s] << "\t"
			<< p - half_points1 << "\t"
			<< X0(p,s) << "\t"
			<< XX(p,s) << "\t"
			<< P0(p,s) << "\t"
			<< PP(p,s) << "\n";
	}
      
      //
      // initiate distance calculation and clustering for this set of intervals
      //

      
      tclst_t tc( use_complex_dist || use_amp ? &X : NULL ,
		  use_complex_dist || use_phase ? &P : NULL ,
		  chs , t , k1 , k2 ,
		  hcK,
		  use_complex_dist );
     
      //
      // output overall feature means
      //

      int c = 0;
      for (int s=0; s<ns; s++)
	{
	  writer.level( signals.label(s) , globals::signal_strat );
	  
	  for (int p=0; p<points; p++)
	    {
	      writer.level( p - half_points1 , globals::sample_strat );
	      if ( use_amp ) writer.value( "A" , tc.tm[ c++ ] );
	      if ( use_phase ) writer.value( "P" , tc.tm[ c++ ] );      
	    }	      
	  writer.unlevel( globals::sample_strat );
	}
      writer.unlevel( globals::signal_strat );	      

      
      //
      // k-means features
      //
      
      Data::Matrix<double> Pa, ZPa;
      if ( use_complex_dist )
	{
	  Pa.resize( ni , ns * points );
	  ZPa.resize( ni , ns * points );
	  for (int i=0; i<ni; i++)
	    {
	      int c = 0;
	      for (int s=0; s<ns; s++)
		for (int p=0; p<points; p++)
		  {
		    Pa(i,c) = P2[i](p,s);
		    ZPa(i,c) = ZP2[i](p,s);
		    ++c;
		  }
	      
	    }	  	  
	}
      
      for (int kn=k1; kn<=k2; kn++)
	{
	  writer.level( kn , "KN" );
	  
	  const Data::Matrix<double> & km = tc.kmeans[kn];	  
	  
	  std::map<int,std::vector<double> > pclmeans, zpclmeans;
	  if ( use_complex_dist )
	    {
	      pclmeans = Statistics::group_means( Pa , tc.ksol[kn] );
	      zpclmeans = Statistics::group_means( ZPa , tc.ksol[kn] );
	    }
	  
	  for (int k=0; k<kn; k++)
	    {
	      writer.level( k+1, globals::cluster_strat );
	      
	      int c = 0, c2 = 0;
	      for (int s=0; s<ns; s++)
		{
		  writer.level( signals.label(s) , globals::signal_strat );
		  
		  for (int p=0; p<points; p++)
		    {
		      writer.level( p - half_points1 , globals::sample_strat );
		      if ( use_complex_dist || use_amp ) writer.value(  use_complex_dist ? "REAL" : "A"  , km( c++ , k ) );
		      if ( use_complex_dist || use_phase ) writer.value(  use_complex_dist ? "IMAG" : "P"  , km( c++ , k ) );
		      if ( use_complex_dist ) writer.value( "PH" , pclmeans[k][c2++] );
		      if ( use_complex_dist ) writer.value( "ZPH" , zpclmeans[k][c2++] );
		    }	      
		  writer.unlevel( globals::sample_strat );
		}
	      writer.unlevel( globals::signal_strat );	      
	    }
	  writer.unlevel( globals::cluster_strat );
	}

      writer.unlevel( "KN" );


      //
      // hierarchical clustering group means
      //

      hcK = tc.sol.k;

      std::map<int,std::vector<double> > pclmeans, zpclmeans;
      if ( use_complex_dist )
	{
	  pclmeans = Statistics::group_means( Pa , tc.sol.best );
	  zpclmeans = Statistics::group_means( ZPa , tc.sol.best );
	}
      
      for (int k=0; k<hcK ; k++)
	{
	  writer.level( k+1, globals::cluster_strat );
	  
	  int c = 0, c2 = 0;
	  for (int s=0; s<ns; s++)
	    {
	      writer.level( signals.label(s) , globals::signal_strat );
	      
	      for (int p=0; p<points; p++)
		{
		  writer.level( p - half_points1 , globals::sample_strat );
		  if ( use_complex_dist || use_amp ) writer.value( use_complex_dist ? "REAL" : "A" , tc.clmeans[k][ c++ ] );
		  if ( use_complex_dist || use_phase ) writer.value( use_complex_dist ? "IMAG" : "P" , tc.clmeans[k][ c++ ]);
		  if ( use_complex_dist ) writer.value( "PH" , pclmeans[k][ c2++ ] );
		  if ( use_complex_dist ) writer.value( "ZPH" , zpclmeans[k][ c2++ ] );		  
		}	      
	      writer.unlevel( globals::sample_strat );
	    }
	  writer.unlevel( globals::signal_strat );	      
	}
      writer.unlevel( globals::cluster_strat );
      
      
      //
      // Class frequencies/assignments
      //

      if ( k1 )
	{
	  for (int k=k1; k<=k2; k++)
	    {
	      writer.level( k , "KN" );

	      std::map<int,int> c;
	      for (int i=0;i<ni;i++)
		{
		  writer.level( i+1 , globals::count_strat );
		  writer.value( "C" , tc.ksol[k][i] );
		  c[ tc.ksol[k][i] ]++;
		}
	      writer.unlevel( globals::count_strat );
	      	      
	      for (int kk=0; kk<k; kk++)
		{
		  writer.level( kk+1 , globals::cluster_strat );
		  writer.value( "P" , c[ kk ] / (double)ni );
		  writer.value( "CNT" , c[ kk ] );
		}
	      writer.unlevel( globals::cluster_strat );
	      
	    }
	  writer.unlevel( "KN" );
	}

      if ( hcK )
	{
	  std::map<int,int> c;
	  for (int i=0;i<ni;i++)
	    {
	      writer.level( i+1 , globals::count_strat );
	      writer.value( "C" , tc.sol.best[i] );
	      c[ tc.sol.best[i] ]++;
	    }
	  writer.unlevel( globals::count_strat ); 
	  
	  for (int kk=0; kk<hcK; kk++)
	    {
	      writer.level( kk+1 , globals::cluster_strat );
	      writer.value( "P" , c[ kk ] / (double)ni );
	      writer.value( "CNT" , c[ kk ] );
	    }
	  writer.unlevel( globals::cluster_strat );	  
	}
      
      
      //
      // Output distance matrix of intervals, if computed
      //

      if ( hcK && param.has( "distance" ) )
	{
	  if ( param.value( "distance" ) == "" )
	    Helper::halt( "no distance=filename given" );
	  
	  const std::string filename = Helper::expand( param.value( "distance" ) );
	  std::ofstream O1( filename.c_str() , std::ios::out );
	  for (int i=0;i<ni; i++)
	    {
	      for (int j=0; j<ni; j++)
		O1 << ( j == 0 ? "" : "\t" ) << tc.D(i,j) ;
	      O1 << "\n";
	    }
	  O1.close();
	}

      //
      // If in complex mode, output phases
      //

      if ( use_complex_dist )
	{
	  // P2[i](p,s) contains the transformed phase information
	  // get means by K 
	  if ( hcK )
	    {
	      
	    }

	  // kmeans
	  if ( k1 )
	    {

	    }
	}
      
      
      //
      // Next cache seed
      //
      
      ++cc;
      
    }
  
}

tclst_t::tclst_t( const std::vector<Eigen::MatrixXd> * X ,
		  const std::vector<Eigen::MatrixXd> * P ,
		  const std::vector<std::string> & chs ,
		  const std::vector<double> & t ,
		  const int k1, const int k2,
		  const int hcK , 
		  const bool use_complex_dist )
{

  const bool hasX = X != NULL;
  const bool hasP = P != NULL;
  if ( ! ( hasX || hasP ) )
    Helper::halt( "bad call of tclst_t" );

  if ( use_complex_dist && ! ( hasX  ||  hasP ) )
    Helper::halt( "bad call of tclst_t" );
  
  // number of intervals
  n = hasX ? X->size() : P->size();
  if ( hasX && hasP &&  P->size() != X->size() )
    Helper::halt( "internal problem in tclst_t() " );

  // number of signals
  const int ns = hasX ? (*X)[0].cols() : (*P)[0].cols();
  
  // number of sample points in each interval
  const int np = hasX ? (*X)[0].rows() : (*P)[0].rows();
  
  logger << "  time-locked clustering for " << n << " " << np << "-point intervals, based on " << ns << " channels\n";
  
  // distance matrix
  D.resize( n , n );
  
  // for k-means clustering
  Data::Matrix<double> XPa( n , ( (int)hasX+(int)hasP ) * ns * np );

  for (int i=0; i<n; i++)
    {
      int c = 0;
      for (int s=0; s<ns; s++)
	for (int p=0; p<np; p++)
	  {
	    if ( hasX ) XPa(i,c++) = (*X)[i](p,s);
	    if ( hasP ) XPa(i,c++) = (*P)[i](p,s);
	  }
    }
  

  //
  // hierarchical clustering on D
  //
  
  if ( hcK )
    {

      //
      // Make distance matrix
      //
      
      for (int i=0; i<n; i++)
	for (int j=0; j<n; j++)
	  {
	    if ( i == j )
	      {
		D(i,j) = 0;
		continue;
	      }
	    
	    if ( i > j )
	      {
		D(i,j) = D(j,i);
		continue;
	      }
	    
	    const Eigen::MatrixXd * xi = hasX ? &(*X)[i] : NULL ;
	    const Eigen::MatrixXd * pi = hasP ? &(*P)[i] : NULL ;
	    
	    const Eigen::MatrixXd * xj = hasX ? &(*X)[j] : NULL ;
	    const Eigen::MatrixXd * pj = hasP ? &(*P)[j] : NULL ;
	    
	    if ( use_complex_dist )
	      {
		double d = 0;
		for (int s=0; s<ns; s++)
		  for (int p=0; p<np; p++)
		    d += pow( ( (*xi)(p,s) - (*xj)(p,s) ) , 2 ) + pow( ( (*pi)(p,s) - (*pj)(p,s) ) , 2 );
		d = sqrt( d );
		D(i,j) = d;
	      }
	    else
	      {
		
		// signals
		double dx = 0;
		
		// phases
		double px = 0; 
		
		for (int s=0; s<ns; s++)
		  for (int p=0; p<np; p++)
		    {
		      if ( hasX ) 
			dx += pow( ( (*xi)(p,s) - (*xj)(p,s) ) , 2 );
		      if ( hasP )
			px += pow( ( (*pi)(p,s) - (*pj)(p,s) ) , 2 );
		    }
		
		if ( hasX )
		  dx = sqrt( dx );
		if ( hasP )
		  px = sqrt( px );
		
		// take (weighted?) sum of these two measures?
		D(i,j) = dx + px;
	      }
	    
	  }

      // scale and combine ?
      if ( ! use_complex_dist )
	{
	  // const double dn = D.sum() / (n*n);
	  // conat double d2n = D2.sum() / (n*n);
	  // std::cerr << " scalibg " << dn <<" " << d2n << "\n";
	  // D /= dn;
	  // D2 /= d2n;
	  // D = D + D2;
	  // D /= 2;
	  // std::cout << " now " << D.sum() /	(n*n) << "\n";
	  
	}
      
      //
      // do clustering
      //

      cluster_t cluster;
      
      // constraint on maximum size of each cluster (0 = no constraint)
      int maxS = 0;
      
      sol = cluster.build( D , hcK , maxS );

      // cluster means

      clmeans = Statistics::group_means( XPa , sol.best );

    }
  

  //
  // Total means
  //
      
  tm = Statistics::mean( XPa ).extract();
  
  //
  // k-means clustering on original data 
  //

  if ( k1 ) 
    {
      for (int k=k1; k<=k2; k++)
	{
	  kmeans_t km;      
	  std::vector<int> ksol1;      
	  Data::Matrix<double> m = km.kmeans( XPa , k , &ksol1 );
	  
	  varexp[ k ] = km.between / ( km.within + km.between );	  
	  //	  logger << "  fit k-means, K = " << k << ", " << varexp[ k ] << " variance explained\n";
	  
	  // save assignments
	  ksol[ k ] = ksol1;
	  
	  // save means
	  kmeans[ k ] = m;
	  
	}
      
    }
       
  
}

