
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

#include "dsp/peaks.h"
#include "param.h"
#include "timeline/cache.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "db/db.h"
#include "defs/defs.h"
#include "miscmath/miscmath.h"

extern writer_t writer;
extern logger_t logger;

void dsptools::peaks( edf_t & edf , param_t & param )
{
  // write peaks to a cache (e.g. for TLOCK)
  const bool to_cache = param.has( "cache" );
  const std::string cache_name = to_cache ? param.requires( "cache" ) : "";
  cache_t<int> * cache = to_cache ? edf.timeline.cache.find_int( cache_name ) : NULL ;
  if ( to_cache ) logger << "  writing peaks to cache " << cache_name << "\n" ;
  
  // write peaks to an annot (e.g. for GED)
  const bool to_annot = param.has( "annot" );  
  const std::string annot_name = to_annot ? param.requires( "annot" ) : "";
  annot_t * annot = to_annot ? edf.annotations->add( annot_name ) : NULL ; 
  uint64_t w_tp = 0; // window around each point
  if ( param.has( "w" ) )
    w_tp = globals::tp_1sec * param.requires_dbl( "w" );
  if ( to_annot ) logger << "  writing peaks to annotation " << annot_name << ", +/- " << param.value( "w" ) << "sec\n";
  
       
  // signals
  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) );  
  const int ns = signals.size();
  
  // options
  bool by_epoch = param.has( "epoch" );  
  bool min = param.has( "min" ) || param.has( "min-only" );
  bool max = param.has( "min-only" ) ? false : true;
  int th_clipped = param.has( "clipped" ) ? param.requires_int( "clipped" ) : 3 ;
  int percentile = param.has( "percentile" ) ? param.requires_int( "percentile" ) : 0 ; 


  // nb.  because of how we track sample-points below, (for now) we cannot
  // allow epochs that are not exactly contiguous (i.e. no overlap, no gaps)
  
  if ( by_epoch && ! edf.timeline.exactly_contiguous_epochs() ) 
    Helper::halt( "can only have exactly contiguous epochs in PEAKS currently" );

  //
  // iterate over signals
  //

  for (int s=0; s<ns; s++)
    {

      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;

      const int sr = edf.header.sampling_freq( signals(s) );

      int ne = by_epoch ? edf.timeline.first_epoch() : -9;

      writer.level( signals.label(s) , globals::signal_strat );

      std::map<std::string,std::string> faclvl = writer.faclvl();
      
      int curr_point = -1;

      // get data
      
      while ( 1 ) 
	{
	  
	  int epoch = by_epoch ? edf.timeline.next_epoch() : -9;

	  if ( epoch == -1 ) break;

	  std::map<std::string,std::string> faclvl1 = faclvl;

	  if ( by_epoch )
	    {
	      writer.epoch( edf.timeline.display_epoch( epoch ) );
	      faclvl1[ "E" ] = Helper::int2str( edf.timeline.display_epoch( epoch ) );
	    }

	  interval_t interval = by_epoch ? edf.timeline.epoch( epoch ) : edf.timeline.wholetrace();
	  
	  slice_t slice( edf , signals(s) , interval );

	  const std::vector<double> * d = slice.pdata();

	  const std::vector<uint64_t> * tp = slice.ptimepoints();
	  
	  // need sample points w.r.t. current EDF structure (i.e. not time-points per se)
	  // but based on the current epoch

	  std::vector<int> sp;

	  if ( by_epoch )
	    {
	      const int n = d->size();
	      sp.resize( n );
	      for (int i=0; i<n; i++) 
		sp[i] = ++curr_point;	      
	      // i.e. assumes contiguous epochs
	      // but this way, the sample-points for this epoch
	      // will map back onto the whole-trace sample-point
	      // list (assuming no RE events are called)
	    }

	  // find peaks
	  peaks_t peaks;
	  peaks.min = min;
	  peaks.max = max;
	  peaks.th_clipped = th_clipped;
	  peaks.percentile = percentile;
	  peaks.detect( d , by_epoch ? &sp : NULL );	  
	  
	  // cache: nb. 'points' is the variable name that TLOCK looks for
	  if ( to_cache )
	    cache->add( ckey_t( "points" , faclvl1 ) , peaks.pk );	  

	  // annots
	  if ( to_annot )
	    {
	      const int np = peaks.pk.size();
	      for (int p=0; p<np; p++)
		{
		  interval_t interval( (*tp)[ peaks.pk[p] ] , (*tp)[ peaks.pk[p] ] );
		  interval.expand( w_tp );
		  annot->add( peaks.ismin[p] ? "-ve" : "+ve" , interval , signals.label(s) );
		}
	    }

	    
	  // all done for this channel if we are in whole-signal mode
	  if ( ! by_epoch ) break;
	  
	} // next interval
      
      if ( by_epoch )
	writer.unepoch();

      // next signal
    }
  
  writer.unlevel( globals::signal_strat );
        
}
      

void peaks_t::detect( const std::vector<double> * x , const std::vector<int> * p_sp )
{

  pk.clear();
  values.clear();
  ismin.clear();

  const int n = x->size();

  // set sample-points
  std::vector<int> sp;
  if ( p_sp != NULL )
    {
      sp = *p_sp;
      if ( sp.size() != n ) Helper::halt( "internal error in PEAKS" );
    }
  else
    {
      sp.resize(n);
      for (int i=0; i<n; i++) sp[i] = i;
    }

  //
  // mask clipped points?
  //

  std::vector<bool> clipped;
  double overall_min, overall_max;

  if ( th_clipped )
    {
      // get signal min/max empirically
      MiscMath::minmax( *x , &overall_min , &overall_max );
      
      clipped.resize( n , false );
      
      if ( max )
	{
	  for (int i=0; i<n; i++)
	    {
	      
	      if ( fabs( overall_max - (*x)[i] ) > EPS ) continue;
	      
	      int j;
	      for ( j = i + 1 ; j < n; j++)
		{
		  if ( fabs( overall_max - (*x)[j] ) > EPS ) 
		    break; 
		}
	      	      
	      // 1 2 3 4
              // X X X .

	      if ( j - i >= th_clipped )
		for (int k=i; k<j; k++) clipped[k] = true;
	      
	      // advance (point 'j' is not clipped, so we can set to j, next
	      // iteration will be j+1
	      i = j;

	    }

	}
      
      if ( min )
	{
	  for (int i=0; i<n; i++)
	    {	      
	      if ( fabs( (*x)[i] - overall_min ) > EPS ) continue;
	      int j;
	      for ( j = i + 1 ; j < n; j++)
		{
		  if ( fabs( (*x)[j] - overall_min ) > EPS ) 
		    break; 
		}	      
	      if ( j - i >= th_clipped )
		for (int k=i; k<j; k++) clipped[k] = true;	      
	      i=j;
	    }	  
	}
      
      // count total clipping
      int n_clipped = 0;
      for (int i=0; i<n; i++) 
	if ( clipped[i] ) ++n_clipped;
            
    }
  

  //
  // end points not included, so 0 means 'not a peak'
  //

  int last = 0; 
  
  for (int i=1; i<(n-1); i++)
    {

      // if at/by a clipped signal part, we can skip
      if ( clipped[i] || clipped[i-1] || clipped[i+1] ) continue;
      
      // otherwise, search for simple peaks (any peak)

      if ( max ) 
	{
	  if ( (*x)[i] > (*x)[i-1] && (*x)[i] > (*x)[i+1] )
	    {
	      pk.push_back(sp[i]);
	      values.push_back( (*x)[i] );
	      ismin.push_back( false );
	    }	  
	}
      
      if ( min ) 
	{
	  if ( (*x)[i] < (*x)[i-1] && (*x)[i] < (*x)[i+1] )
	    {
	      pk.push_back(sp[i]);
	      values.push_back( (*x)[i] );
	      ismin.push_back( true );
	    }	  
	}

      // next point
    }
  
  
  //
  // Percentile thresholds?  Do separately for max and min 
  //
    
  // get absolute values
  
  if ( percentile > 0 && percentile < 100 ) 
    {
      std::vector<double> maxval, minval;
      std::vector<int> maxpk, minpk;
      for (int i=0; i<values.size(); i++)
	{
	  if ( ismin[i] )
	    {
	      minpk.push_back( pk[i] );
	      minval.push_back( -values[i] ); // nb, sign flipped
	    }
	  else 
	    {
	      maxpk.push_back( pk[i] );
	      maxval.push_back( values[i] );
	    }
	}
      
      pk.clear();
      values.clear();
      ismin.clear();
      
      if ( max ) 
	{
	  // get value X = top N% and set to 0/1 if below/above X                                                                                                                                  
	  double threshold = MiscMath::percentile( maxval , 1.0 - percentile / 100.0 ) ;
	  std::vector<double> maxval2;
	  std::vector<int> maxpk2;
	  for (int i=0;i<maxval.size(); i++)
	    {
	      if ( maxval[i] >= threshold ) 
		{
		  maxval2.push_back( maxval[i] );
		  maxpk2.push_back( maxpk[i] );
		}
	    }
	  
	  // copy back into originals
	  for (int i=0; i<maxpk2.size(); i++)
	    {
	      values.push_back( maxval2[i] );
	      pk.push_back( maxpk2[i] );
	      ismin.push_back( false );
	    }
	}
      
      
      // do the same for minima (taking the smalltest peaks)
      
      if ( min ) 
	{
	  
	  // get value X = top N% and set to 0/1 if below/above X
	  // nb. sign flipped here

	  double threshold = MiscMath::percentile( minval , 1.0 - percentile / 100.0 ) ;

	  std::vector<double> minval2;
	  std::vector<int> minpk2;

	  for (int i=0;i<minval.size(); i++)
	    {
	      if ( minval[i] >= threshold ) // nb sign flipped
		{
		  minval2.push_back( minval[i] );
		  minpk2.push_back( minpk[i] );
		}
	    }
	  
	  // copy back into originals
	  for (int i=0; i<minpk2.size(); i++)
	    {
	      values.push_back( -minval2[i] ); // nb flip sign back
	      pk.push_back( minpk2[i] );
	      ismin.push_back( false );
	    }
	}
      
    } // end of percentile code


  //
  // all done
  //

}
 
