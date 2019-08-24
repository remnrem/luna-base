
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

#include "fiplot.h"

#include "miscmath/miscmath.h"
#include "cwt/cwt.h"
#include "fftw/fftwrap.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "miscmath/miscmath.h"

#include <string>
#include "db/db.h"
#include "eval.h"

#include "helper/helper.h"
#include "helper/logger.h"

extern writer_t writer;

extern logger_t logger;

void fiplot_wrapper( edf_t & edf , const param_t & param , const std::vector<double> * raw , const int * sr )
{

  //
  // Parameters
  //

  //
  // Time (or cycle) bins
  //
  
  double t_lwr = param.has( "t-lwr" ) ? param.requires_dbl( "t-lwr" ) : 0.1 ;
  double t_upr = param.has( "t-upr" ) ? param.requires_dbl( "t-upr" ) : 4 ;
  double t_inc = param.has( "t-inc" ) ? param.requires_dbl( "t-inc" ) : 0.1  ;

  //
  // Scale by cycles instead of seconds?
  //

  const bool cycles = param.has( "by-cycles" ) || param.has( "c-lwr" ); 

  if ( param.has( "c-lwr" ) ) t_lwr = param.requires_dbl( "c-lwr" );
  if ( param.has( "c-upr" ) ) t_upr = param.requires_dbl( "c-upr" );
  if ( param.has( "c-inc" ) ) t_inc = param.requires_dbl( "c-inc" );

  //
  // Frequencies
  //

  // use either f-log {# of steps} OR f-inc {increment size}

  bool logspace = param.has("f-log"); // now interpret f-inc as the number of steps to have

  const double f_lwr = param.has( "f-lwr" ) ? param.requires_dbl( "f-lwr" ) : 1  ;
  const double f_upr = param.has( "f-upr" ) ? param.requires_dbl( "f-upr" ) : 20 ;
  const double f_inc = param.has( "f-log" ) ? param.requires_dbl( "f-log" ) : 
    ( param.has( "f-inc" ) ? param.requires_dbl( "f-inc" ) : 1 )  ;
  const int num_cyc = param.has( "cycles" ) ? param.requires_int( "cycles" ) : 7 ;

  
  //
  // Threshold, i.e. only consider signals above th times the mean
  //
  
  const double th = param.has( "th" ) ? param.requires_dbl( "th" ) : 0 ; 
  
  const bool normalize = param.has( "norm" );
  
  const bool logit = param.has( "log" );

  //
  // Input data 
  //

  //
  // Use EDF signals
  //
  
  if ( raw == NULL )
    {
      
      std::string signal_label = param.requires( "sig" );
  
      signal_list_t signals = edf.header.signal_list( signal_label );
      
      const int ns = signals.size();
      
      std::vector<double> Fs = edf.header.sampling_freq( signals );  
      
      interval_t interval = edf.timeline.wholetrace(); 
      
      for (int s = 0 ; s < ns ; s++ ) 
	{
	  
	  // Only consider raw signal channels
	  
	  if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
	  
	  // Output stratifier
	  
	  writer.level( signals.label(s) , globals::signal_strat );
	  
	  // Pull all data
	  
	  slice_t slice( edf , signals(s) , interval );
	  
	  const std::vector<double> * d = slice.pdata();
	  
	  const std::vector<uint64_t> * tp = slice.ptimepoints();
	  
	  // ~no data?
	  if ( d->size() <= 2 ) continue;
	  
	  // do plot 
	  
	  fiplot_t fp( *d , tp , Fs[s] , 
		       th , normalize , logit , 
		       t_lwr, t_upr, t_inc , cycles , 
		       f_lwr, f_upr, f_inc , num_cyc , logspace );
	  
	  
	  //
	  // Next signal
	  //
	  
	  writer.unlevel( globals::signal_strat );
      
	} 
      
    }

  //
  // F/I plot for a raw signal
  //

  else
    {

      if ( sr == NULL ) Helper::halt( "no SR specified for raw signal" );
      
      // ~no data?                                                                                                                                                    
      if ( raw->size() <= 2 ) Helper::halt( "no signal" );
      
      writer.level( "RAW" , globals::signal_strat );
      
      double dt = 1 / (double)(*sr);
      
      std::vector<uint64_t> tp( raw->size() );

      for (int i=0;i<raw->size();i++) tp[i] = i * dt;

      // do plot                                                                                                                                                      
      
      fiplot_t fp( *raw , &tp , *sr , 
		   th, normalize , logit , 
		   t_lwr, t_upr, t_inc , cycles ,
		   f_lwr, f_upr, f_inc , num_cyc , logspace );


      writer.unlevel( globals::signal_strat );

    }
      
}


void fiplot_t::set_f( double lwr , double upr , double inc , bool logspace , int num_cyc ) 
{

  frqs.clear();
  
  f_lwr = lwr;
  f_upr = upr;
  f_inc = inc; // # of inc if logspace == T
  num_cycles = num_cyc;

  if ( ! logspace )
    {
      for (double f = f_lwr ; f <= f_upr ; f += f_inc ) frqs.push_back( f );
    }
  else
    {
      frqs = MiscMath::logspace( f_lwr , f_upr , f_inc );
    }
}


void fiplot_t::proc( const std::vector<double> & x , const std::vector<uint64_t> * tp , const int fs )
{
  
  // time points (secs)
  nt = ( t_upr - t_lwr ) / t_inc;
  
  // frq points (Hz)
  nf = frqs.size();

  for ( int fi = 0 ; fi < nf ; fi++ )
    {
      
      double f = frqs[fi];

      logger << "  assessing " << f << " Hz ...";

      writer.level( f , globals::freq_strat );
      
      // get CWT 

      std::vector<double> c = cwt( x , fs , f , num_cycles );

      //        for (int j=0;j<c.size();j++)
//       for (int j=0;j<1000;j++)
// 	  {
//   	  std::cout << "x\t" 
//   		    << f << "\t"
//   		    << j << "\t"
//   		    << c[j] << "\n";
//   	}
      
      // get intervals
      
      fibin_t r = intervalize( c , tp , fs , t_lwr , t_upr , t_inc , cycles , f ); 
      
      // report
      std::map<int,fipair_t>::const_iterator rr = r.r.begin();
      while ( rr != r.r.end() )
	{
	  writer.level( r.t[ rr->first ] , "TBIN" );
	  writer.value ( "FIP" , rr->second.w );
	  //writer.value ( "FIPN" , rr->second.n );
	  ++rr;
	}
      writer.unlevel( "TBIN" );
    } 
  writer.unlevel( globals::freq_strat );
}


fibin_t fiplot_t::intervalize( const std::vector<double> & x_ , 
			       const std::vector<uint64_t> * tp , 
			       const int fs , 
			       const double t_lwr , 
			       const double t_upr , 
			       const double t_inc , 
			       const bool plot_by_cycles ,
			       const double fc )
{

  fibin_t r;
  
  //
  // number of time/cycle points
  //
  
  int nt = 1 + ( t_upr - t_lwr ) / t_inc ;
  
  //
  // time-bin increment; dt2 defined to find discontinuities
  //
  
  const double dt = 1.0/(double)fs;
  
  const uint64_t dt2 = globals::tp_1sec * ( dt * 1.5 );
  
  //
  // point duration to time-point: sample-points are inclusive so +1
  // double time = ( j - i + 1 ) * dt  
  //


  

  //
  // Scale?
  //

  // copy, as we need to rescale
  std::vector<double> x = x_;

  const int n = x.size();

  if ( logit ) 
    {
      for (int i=1;i<n;i++)
	x[i] = log( x[i] );
    }

    
  //
  // Rescale measure to 0..1 (i.e. to normalize across different Fc
  //

  if ( normalize ) 
    {

      double min = x[0] , max = x[0];
      
      for (int i=1;i<n;i++)
	{
	  if      ( x[i] < min ) min = x[i];
	  else if ( x[i] > max ) max = x[i];
	}
      
      if ( max == min ) Helper::halt( "flat signal" );
      
      const double rng = max - min;
      
      for (int i=0;i<n;i++) 
	x[i] = ( x[i] - min ) / rng ; 
    }


  //
  // Take only above threshold points?
  //

  double yt = 0;
  
  if ( th > 0 ) 
    {
      yt = th * MiscMath::mean( x );
      logger << " setting " << th << "x threshold to " << yt << "\n";
    }

  //
  // rule of thumb: only consider intervals that are at least twice
  // as long as the longest time span considered e.g. 5 seconds
  // requires a full 30-second epoch, etc
  //
  
  uint64_t required_tp = (uint64_t)2 * t_upr * globals::tp_1sec ; 
  
  //
  // scan for discontinuities
  //
  
  int first_idx = -1;
  
  //
  // T means the next point is not adjacent (or last point)
  //

  std::vector<bool> disc( n , false ); 
  
  //
  // track how of the total time interval gets excluded (i.e. if intervals longer than epochs)
  //
  
  double included_seconds = 0;
  
  double all_seconds = 0;
  
  for (int i=0;i<n;i++) 
    {
      
      //
      // track leftmost point of contiguous time
      //
      
      if ( first_idx == -1 ) first_idx = i;
      
      //
      // discontinuity?
      //

      if ( i == n-1 ) disc[i] = true;
      else if ( (*tp)[i+1] - (*tp)[i] > dt2 ) { disc[i] = true; } 
      
      //
      // was this segment long enough?
      //
      
      if ( disc[i] ) 
	{
	  
	  uint64_t length = i == n-1  // last point?
	    ? (*tp)[i] - (*tp)[first_idx] + ( (*tp)[1] - (*tp)[0] ) // add extra point
	    : (*tp)[i+1] - (*tp)[first_idx];
	  
	  // duration of segment
	  all_seconds += length * globals::tp_duration ; 
	  
	  
	  if ( length < required_tp ) 
	    {	 
	      // skip whole region if not long enough 
	      for (int j=first_idx;j<=i;j++) disc[j] = true;
	    }
	  else
	    included_seconds +=  length * globals::tp_duration ;
	  
	  // reset leftmost point of segment 
	  first_idx = -1;
	}
    }
  
  logger << "  including " << ( included_seconds / all_seconds ) * 100 << "% of " << all_seconds << " seconds\n";
  

  //
  // Initialise time-bins
  //
  
  r.t.clear();
  
  for (int t=0;t<nt;t++) 
    {
      // score as the midpoint of each time-bin      
      r.t.push_back( t_lwr + t * t_inc + 0.5 * t_inc );
      r.r[ t ] = fipair_t(0,0);  
    }


  //
  // for each point, see how far we can go before we find a lower
  // point (or a discontinuity);   i.e. there may be boundary issues if 
  // a very fragmented time-series, as we will underestimate the true length
  // of any intervals that start... but not sure an easy workaround at this
  // point
  //
  
  // fipoint_t (i,h,t)  
  
  std::set<fipoint_t> pts;
  
  for (int i=0;i<n;i++)
    {

      //
      // min threshold at which we ignore
      //

      if ( x[i] <= yt ) continue;
      
      //
      // at a discontinuity?
      //

      if ( disc[i] ) continue;


      for (int j=i+1;j<n;j++)
	{

	  // did we meet the returning down-slope? 
	  if ( ! disc[j] )
	    {
	      if ( x[j] < x[i] ) 
		{
		  fipoint_t f( i , j , x[i] ) ;
		  pts.insert( f ) ;
		  break;
		}
	    }
	  else
	    {
	      // or, we hit end/discontinuity?
	      
	      fipoint_t f( i , j , x[i] ) ;
	      pts.insert( f ) ;
	      //i = j+1; // skip ahead
	      break;	      
	    }
	}
    }
  
  logger << " decomposed signal into " << pts.size() << " elements\n";
  
  if ( 0 )
    {
      
      std::cout << "decomposed signal into " << pts.size() << " elements\n";
      std::set<fipoint_t>::const_iterator ff = pts.begin();
      while ( ff != pts.end() ) 
	{
	  std::cout << ff->t << "\t"
		    << ff->i << " - " << ff->j << "\t"
		    << ff->h << "\n";
	  ++ff;
	}
    }

  
  //
  // Build the frequency/interval map
  //
  
  // track how many 'height' values for each point actually are represented in a
  // all intervals spanning that position (should be ~all)

  std::vector<double> used( n , 0 );

  //
  // go through intervals, starting with the longest (i.e. based on
  // fipoint_t::operator<() ) and add up the interval map
  //
  
  std::set<fipoint_t>::const_iterator ff = pts.begin();
  while ( ff != pts.end() ) 
    {
      
      //
      // time duration ( seconds )
      //
      
      double t = (double)( ff->j - ff->i + 1 ) * dt;
      
      //
      // or in cycles?
      //

      if ( plot_by_cycles ) t *= fc;
      
      
      //
      // Inside range?
      //
      
      if ( t >= t_lwr && t < t_upr )
	{
	  
	  int bin = ( t - t_lwr ) / t_inc; 
	  
// 	  double tbin = t_lwr + bin * t_inc + 0.5 * t_inc;  	  
// 	  std::cout  << "bin = " << nt << "\t" << bin << " " << tbin << "\n";
	  
	  //
	  // add value to bin
	  //
	  
	  double add = 0;
	  double amt = 0;
	  
	  // for each spanned sample-point

	  for (int i= ff->i ; i <= ff->j ; i++ ) 

	    {
	      // the additional height beyond what is already accounted for 
	      // by longer spanning intervals

	      double part = ff->h - used[i];
	      
	      // keep track for the F/I bin
	      add += part;
	      
	      // and for this sample point, for the next (shorter) spanning interval
	      used[i] += part;	      
	      
	      ++amt;
	    }
	  
	  // store in F/I bin
	  r.r[ bin ].w += add;
	  r.r[ bin ].n += amt;
	  
	}

      // bin       
      ++ff;
    }

  // at end, x should equal used
  //    for (int i=0;i<n;i++) std::cout << "x/u " << fc << "\t" << x[i] << "\t" << used[i] << "\n";  
  
  // all done
  
  return r;
}


std::vector<double> fiplot_t::cwt( const std::vector<double> & x , const int fs , const double fc , const int num_cycles )
{
  std::vector<double> r;

  CWT cwt;
  cwt.set_sampling_rate( fs );
  cwt.add_wavelet( fc , num_cycles );  
  cwt.load( &x );  
  cwt.run();
  
  const std::vector<double> & results = cwt.results(0);

  return results;

  // OR... do some smoothing?

   const double moving_window_sec = 0.1;
   int window_points = moving_window_sec * fs;
   if ( window_points % 2 == 0 ) ++window_points;
	  
   const std::vector<double> averaged 
     = MiscMath::moving_average( results , window_points );
  
  return averaged;


  
}
