
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

#include "dsp/tlock.h"
#include "timeline/timeline.h"
#include "stats/statistics.h"
#include "db/db.h"
#include "helper/logger.h"
#include "edf/edf.h"
#include "edf/slice.h"

extern writer_t writer;
extern logger_t logger;


void dsptools::tlock( edf_t & edf , param_t & param )
{
  
  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) );

  if ( signals.size() < 1 ) return;

  const int ns = signals.size();

  // check sample rates: must be similar across all 

  std::vector<double> Fs = edf.header.sampling_freq( signals );
  

  for (int s=1;s<Fs.size();s++) 
    if ( Fs[s] != Fs[0] ) 
      Helper::halt( "sample rates must be similar across signals for TLOCK" );

  //
  // Options for how to handle inputs:  take logs before means, or treat as 
  // angles (circular values);  here, rather than calculate a circular mean, or ITPC
  // etc, for now output a binned histogram of counts (column nornalized) by angle bins
  // e.g. 18 bins of 20 degrees
  //

  bool take_log = param.has( "tolog" ) ;

  int angle_bins = param.has( "phase" ) ? param.requires_int( "phase" ) : 0 ;

  if ( take_log && angle_bins ) Helper::halt( "cannot specify both tolog and phase" );

  //
  // get time-points: expecting a cache in sample-point units (i.e. 
  // to ensure same number of points per window; i.e. this is also why
  // all signals must have the same SR
  //
  
  
  //
  // get window
  //
  
  double half_window = param.requires_dbl( "w" );
  if ( half_window <= 0 ) Helper::halt( "w must be a positive number" );

  int half_points = half_window * Fs[0] ;
  int points = 1 + 2 * half_points;
  
  std::vector<double> t;
  for ( double w = -half_window ; w <= half_window ; w += 1.0/Fs[0] ) 
    t.push_back(w);

  if ( t.size() != points ) 
    Helper::halt( "internal error constructing window" );

  //
  // Normalisation: by default, 20% + 20% of window, i.e. skip middle 60%
  //
  
  int norm_points = points * ( param.has( "np" ) ? param.requires_dbl( "np" ) : 0.2 );
  
  if ( norm_points < 0 || norm_points > points / 2 ) 
    Helper::halt( "expecting np between 0 and 0.5" );
  
  //
  // Get seed sample-points from cache
  //

  std::string cache_name = param.requires( "cache" );

  if ( ! edf.timeline.cache.has_int( cache_name ) )
    Helper::halt( "cache not found for this individual: " + cache_name );

  cache_t<int> * cache = edf.timeline.cache.find_int( cache_name );


  cache->dump();

  //
  // Get any/all keys associated with the 'peak' intrnal name cache
  //

  std::set<ckey_t> ckeys = cache->keys( "peak" );

  std::set<ckey_t>::const_iterator cc = ckeys.begin();

  while ( cc != ckeys.end() )
    {
      
      std::vector<int> cx = cache->fetch( *cc );

      //
      // add output stratifiers based on this key
      //

      std::map<std::string,std::string>::const_iterator ss = cc->stratum.begin();
      while ( ss != cc->stratum.end() )
	{
	  writer.level( ss->second , "S_" + ss->first );
	  ++ss;
	}

      //
      // perform separately for each signal
      //

      for (int s=0; s<ns; s++)
	{
	  
	  writer.level( signals.label(s) , globals::signal_strat );
	  
	  // get data and TP information 
	  
	  slice_t slice( edf , signals(s) , edf.timeline.wholetrace() );
	  
	  const std::vector<double> * d = slice.pdata();
	  
	  const std::vector<uint64_t> * tp = slice.ptimepoints();
	  
	  // build up TLOCK matrix
	  
	  tlock_t tlock(t, norm_points );
      
	  // loop time-points to sync/lock on
	  
	  for ( int i=0; i<cx.size(); i++)
	    {
	  
	      int lower = cx[i] - half_points;
	      int upper = cx[i] + half_points;
	      
	      // interval out-of-range
	      if ( lower < 0 || upper > d->size() ) 
		continue;
	      
	      // TODO: check for discontunuities (EDF+)	  
	      //static bool discontinuity( const std::vector<uint64_t> & t , int sr, int sp1, int sp2 );
	      

	      tlock.add( d , lower , upper , take_log , angle_bins );
	      
	      // next sample-point
	    }
	  
	  

	  //
	  // Report either as angles, or as means
	  //

	  if ( angle_bins ) 
	    {
	      Data::Matrix<double> angbin = tlock.angles();
	      
	      std::cout << angbin.dim1() << " "<< angbin.dim2() << "\n";
	      std::cout << angle_bins << " " << t.size() << "\n";

	      if ( angbin.dim1() != angle_bins || angbin.dim2() != t.size() )
		Helper::halt( "internal error in tlock_t()" );

	      for (int i=0; i<angle_bins; i++)
		{
		  writer.level( i+1 , "PHASE" );
		  for (int j=0; j<t.size(); j++)
		    {
		      writer.level( t[j] , "SEC" );		      
		      writer.value( "M" , angbin(i,j) );
		    }
		}
	      writer.unlevel( "PHASE" );
	      writer.unlevel( "SEC" );
	      
	    }
	  
	 
	  //
	  // report summaries: for regulat values, get the average
	  // and normalize (by edges)
	  //
	  
	  else
	    {
	      Data::Vector<double> means = tlock.average();
	  
	      if ( means.size() != t.size() ) 
		Helper::halt( "internal error in tlock_t()" );
	      
	      for (int i=0; i<means.size(); i++) 
		{
		  writer.level( t[i] , "SEC" );
		  writer.value( "M" , means[i] );
		}

	      writer.unlevel( "SEC" );
	    }


	} // next signal
      
      writer.unlevel( globals::signal_strat );

      //
      // clear key output stratifiers
      // (clear now versus after last loop, as the strata may differ between keys)
      //
      

      ss = cc->stratum.begin();
      while ( ss != cc->stratum.end() )
	{
	  writer.unlevel( "S_" + ss->first );
	  ++ss;
	}


      //
      // Next strataum
      //

      ++cc;
    }
  
  

}


tlock_t::tlock_t( const std::vector<double> & t , const int norm_points ) 
  : t(t) , norm_points(norm_points)
{
  
}

void tlock_t::add( const std::vector<double> * x , 
		   const int lower , const int higher , 
		   const bool take_log , 
		   const int angle_bins )
{

  if ( higher - lower + 1 != t.size() )
    Helper::halt( "internal error");

  Data::Vector<double> d( t.size() ) ;


  //
  // is this an angle? (build histogram)
  //

  if ( angle_bins )
    {
      // assume angle as 0..360
      double wa = 360.0 / (double)angle_bins;

      int j = 0 ;

      for ( int i = lower ; i <= higher ; i++ )
	{
	  // expect radians, -PI .. +PI
	  double deg = ( M_PI + (*x)[i] ) * 180.0 / M_PI;
	  if ( deg < 0 || deg > 360 ) Helper::halt( "value not a valid angle" );
	  int ia = deg / wa;
	  if ( ia == angle_bins ) ia = 0;
	  // store index
	  d[j++] = ia;
	}
      
      // at start, need to set up X to be [ angle bins x epochs ] 
      
      if ( X.dim1() == 0 )
	X.resize( d.size() , angle_bins  );
     
      // accumlate counts
      for (int j=0;j<d.size(); j++)
	X(j,d[j])++;
            
      
    }
  
  //
  // or else a regular value (take means, after normalizing )
  //
  
  else
    {

      int j = 0;
      
      if ( take_log ) 
	for ( int i = lower ; i <= higher ; i++ )
	  d[j++] = log( (*x)[i] );
      else
	for ( int i = lower ; i <= higher ; i++ )
	  d[j++] = (*x)[i] ;
      
      if ( X.dim1() == 0 )
	{
	  X.resize( d.size() , 1 );
	  for (int j=0;j<d.size(); j++)
	    X(j,0) = d[j];
	}
      else
	X.add_col ( d ) ;
      
    }

  
}


Data::Vector<double> tlock_t::average( ) const
{  
  // return row means (transpose + col means ) 
  Data::Matrix<double> Xt = Statistics::transpose( X );
  Data::Vector<double> means = Statistics::mean( Xt );

  // normalize by window edges (defaukt 20% either side)
  if ( norm_points )
    {
      double norm = 0; 

      const int n = means.size() ; // should be 'points'

      for (int i=0; i<norm_points; i++)
	{
	  norm += means[i];
	  norm += means[n-(i+1)];      
	}
      norm /= 2.0 * norm_points;
      for (int i=0; i<n; i++) means[i] /= norm;
    }
  return means;
} 


Data::Matrix<double> tlock_t::angles() const 
{
  // [ time-points/rows X bins/cols ] 
  // normalize each time points 
  Data::Matrix<double> C = Statistics::transpose( X );
  // make each row sum to 1.0
  Data::Vector<double> S = Statistics::col_sums( C );
  
  for (int i=0;i<C.dim1();i++)
    for(int j=0;j<C.dim2();j++)
      C(i,j) /= S[j];
  
  return C;
}
