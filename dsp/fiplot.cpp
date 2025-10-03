
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
#include "dsp/emd.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "miscmath/miscmath.h"

#include <string>
#include "db/db.h"
#include "param.h"

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

  double f_lwr = param.has( "f-lwr" ) ? param.requires_dbl( "f-lwr" ) : 1  ;
  double f_upr = param.has( "f-upr" ) ? param.requires_dbl( "f-upr" ) : 20 ;
  double f_inc = param.has( "f-log" ) ? param.requires_dbl( "f-log" ) : 
    ( param.has( "f-inc" ) ? param.requires_dbl( "f-inc" ) : 1 )  ;
  int num_cyc = param.has( "cycles" ) ? param.requires_int( "cycles" ) : 7 ;

  //
  // Or, we are just looking at a raw signal (i.e. already we have peaks) 
  //

  if ( param.yesno( "envelope" ) )
    {
      // i.e. instead of doing CWT for various bands, just take the
      // signal as is and intervaliz / make a single "interval plot"
      // IP .. i.e. no "F" component now; just take the envelope usin HT

      // set f to negative
      f_lwr = f_upr = f_inc = -1;
      num_cyc = 0;
    }

  
  //
  // Threshold, i.e. only consider signals above th times the mean
  //
  
  const double th = param.has( "th" ) ? param.requires_dbl( "th" ) : 0 ; 
  
  const bool normalize = param.has( "norm" );
  
  const bool logit = param.has( "log" );

  //
  // Output options
  //

  const bool verbose = param.has( "verbose" );
  
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
		       f_lwr, f_upr, f_inc , num_cyc , logspace , 
		       verbose );
	  
	  
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
		   f_lwr, f_upr, f_inc , num_cyc , logspace ,
		   verbose );


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

  // special case: do not do CWT
  
  if ( num_cycles == 0 )
    {      
      frqs.push_back( -1 ); // min ENV
      frqs.push_back( -2 ); // max ENV
      frqs.push_back( -3 ); // original
      return;
    }
  
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
  nt = ( t_upr - t_lwr ) / t_inc + 1; 
  
  // frq points (Hz)
  nf = frqs.size();

  // if not CWT, but getting envelope of a single signal, do here once
  std::vector<double> mine, maxe;
  if ( num_cycles == 0 )
    {
      std::vector<double> ignore_mean_env = emd_t::envelope_mean( x , false , &mine, &maxe );
    }
  
  for ( int fi = 0 ; fi < nf ; fi++ )
    {
      
      double f = frqs[fi];

      std::vector<double> c;

      // CWT

      if ( f > 0 )
	{
	  logger << "  assessing " << f << " Hz ...";
	  
	  writer.level( f , globals::freq_strat );
	  
	  c = cwt( x , fs , f , num_cycles );
	  
	  if ( verbose ) 
	    for (int i=0;i<c.size();i++)
	      std::cout << "CWT\t" << c[i] << "\t" << x[i] << "\n";
	}
      else
	{
	  // get envelope -- fudge, should make
	  // things clearer...
	  if ( fi == 0 ) // min
	    c = mine;
	  else if ( fi == 1 ) // max
	    c = maxe;
	  else // orig-raw
	    {
	      // get 0 .. 1 scale
	      double min, max;
	      MiscMath::minmax( x , &min, &max );
	      min *= 1.01; max *= 1.01;
	      double rng = max - min;
	      c = x;
	      for (int i=0; i<c.size(); i++) c[i] = ( c[i] - min ) / rng ;
	    }
	  
	  // std::cout << "fi = " << fi << "\n";
	  // for (int i=20000; i<20400; i++)
	  //   std::cout << x[i] << "\t" << c[i] << "\n";

	  // indicate w/ F == -1 and +1 for min/max env, and 0 for orig
	  if ( fi == 0 )
	    writer.level( -1 , globals::freq_strat );
	  else if ( fi == 1 ) 
	    writer.level( +1 , globals::freq_strat );
	  else
	    writer.level( 0 , globals::freq_strat );
	}
      
      // get intervals

      fibin_t r = intervalize( c , tp , fs , t_lwr , t_upr , t_inc , cycles , f ); 
      
      // report

      // FIP : sum-FIP / number of seconds
      // ZIP : sum-FIP normalized to sum to 1.0 across the row (F)

      const double tsec = c.size() / (double)fs;
      
      double fsum = 0;
      std::map<int,fipair_t>::const_iterator rr = r.r.begin();
      while ( rr != r.r.end() )
	{
	  fsum += rr->second.w;
	  ++rr;
	}
      
      rr = r.r.begin();
      while ( rr != r.r.end() )
	{
	  writer.level( r.t[ rr->first ] , "TBIN" );
	  writer.value ( "FIP" , rr->second.w / tsec );
	  writer.value ( "ZIP" , rr->second.w / fsum );
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
  // number of time/cycle points (intervals/bins)
  //
  
  int nt = ( t_upr - t_lwr ) / t_inc ;
  
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
  // Rescale measure (for each frequency)
  //

  if ( normalize ) 
    {

      
      //x = MiscMath::Z( x );
      
      // scaling: 0..1 (i.e. to normalize across different Fc)
      
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
      double m0 = MiscMath::mean( x );
      yt = th * m0;
      logger << " setting " << th << "x threshold to " << yt
	     << " ( = " << th << " * mean of " << m0 << ")\n";
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
  // track how of the total time interval gets excluded (i.e. if
  // intervals longer than epochs)
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
      // discontinuity? (or last point)
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
  

  if ( verbose ) 
    logger << "  including " << ( included_seconds / all_seconds ) * 100
	   << "% of " << all_seconds << " seconds\n";
  

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

  if ( verbose )
    {
      int n_disc = 0 , n_above_th = 0;
      for (int i=0;i<n;i++)
	{ 
	  if ( disc[i] ) ++n_disc;
	  if ( x[i] >= yt ) ++n_above_th;
	  //std::cerr << "dets = " << x[i] << "\n";
	}
      
      logger << " of " << n << " points, " << n_disc << " discordancies, "
	     << n_above_th << " (" << (n_above_th/double(n))*100.0 << "%) above threshold\n";
    }


  //
  // fipoint_t (i,h,t)  
  //
  
  std::set<fipoint_t> pts;
  
  for (int i=0;i<n;i++)
    {

      //
      // min threshold at which we ignore
      //
      
      if ( x[i] < yt ) continue;
      
      //
      // at a discontinuity?
      //

      if ( disc[i] ) continue;
      
      bool rising_slope = false;
	
      for (int j=i+1;j<n;j++)
	{
	  
	  // hit a discontinuity? (or last point)
	  // still need to add, so that any higher
	  // points that are within this longer, lower
	  // truncated interval are still baseline-corrected
	  
	  if ( disc[j] )
	    {
	      if ( rising_slope )
                {
		  // true --> implies a truncated interval
		  // (is not added to summ stats below)
                  fipoint_t f( i , j , x[i] , true ) ;
                  pts.insert( f ) ;
                }	      
	      break;
	    }


	  //
	  // do we see any rising slope here?
	  //
	  
	  if ( ( ! rising_slope ) && x[j] >= x[i] )
	    rising_slope = true;
	  
	  // did we meet the returning down-slope? 
	  // need to have spanned at least one other point
	  // i.e. else on a down-slope
	  
	  if ( x[j] < x[i] ) 
	    {
	      if ( rising_slope )
		{
		  fipoint_t f( i , j , x[i] ) ;
		  pts.insert( f ) ;
		}
	      break;
	    }
	}
    }

  // last point: ignore
  // fipoint_t f( n-1 , n-1 , x[n-1] ) ;
  // pts.insert( f ) ;
  
  if ( verbose )
    {
      
      logger << "decomposed signal into " << pts.size() << " elements\n";

      std::set<fipoint_t>::const_iterator ff = pts.begin();
      while ( ff != pts.end() ) 
	{
	  std::cout << "EL\t" << ff->t << "\t"
		    << ff->i << " - " << ff->j << "\t"
		    << ff->h << "\n";
	  ++ff;
	}
    }
  
  
  //
  // Build the frequency/interval map
  //
  
  // track how many 'height' values for each point are represented across
  // all intervals spanning that position

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

      //      std::cout << " considering " << t << " " << ff->i << "\t" << ff->j << "\t" << ff->h << "\n";
      
      //
      // or in cycles?
      //

      if ( plot_by_cycles ) t *= fc;

      //
      // As testing in descending order, if this interval is lower, then we're all done
      //
      
      if ( t < t_lwr ) break;
      
      //
      // If interval is longer, *and/or* if it is truncated, we still need to adjust 'used' so
      // that subsequent added intervals are appropriately baseline-adjusted
      //

      // check both for floating point misses?
      int bin = ( t - t_lwr ) / t_inc;
      
      if ( ff->trunc || t >= t_upr || bin >= nt )
	{

	  // for each spanned sample-point
	  for (int i = ff->i ; i <= ff->j ; i++ ) 
	    {
	      // this is as below, i.e.
	      //   part = ff->h - used[i]
	      //   used += part
	      // i.e. implies:
	      // used = ff->h, if we don't need to track
	      used[i] = ff->h; 
	    }
	}
	   
      //
      // ... otherwise, implies inside the range, so add to a bin
      //
      
      else
	{
	 
	  
	  // TMP
	  // double tbin = t_lwr + bin * t_inc + 0.5 * t_inc;  	  
	  // std::cout  << "adding bin = " << nt << "\t" << bin << " " << tbin << "\n";
	  // TMP
	  
	  //
	  // add value to bin
	  //
	  
	  double add = 0;
	  double amt = 0;
	  
	  // for each spanned sample-point

	  for (int i = ff->i ; i <= ff->j ; i++ ) 
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

	  //std::cout << " int " << tbin << "\t" << ff->i << "\t" << ff->j << "\t" << ff->h << "\t" << add << "\n";
	  // store in F/I bin
	  r.r[ bin ].w += add;
	  r.r[ bin ].n += amt;
	  
	}
      

      // bin       
      ++ff;
    }

  // at end, x should equal used
  //    for (int i=0;i<n;i++) std::cout << "x/u " << fc << "\t" << x[i] << "\t" << used[i] << "\n";  
  //for (int i=0;i<n;i++) std::cout << "x/u " << x[i] << "\t" << used[i] << "\n";  
  
  // all done
  
  return r;
}


std::vector<double> fiplot_t::cwt( const std::vector<double> & x , const int fs , const double fc , const int num_cycles )
{
  std::vector<double> r;

  CWT cwt;
  cwt.set_sampling_rate( fs );
  //cwt.add_wavelet( fc , num_cycles );  
  cwt.alt_add_wavelet( fc , CWT::pick_fwhm( fc ) , 10 );
  logger << "  cwt: fc = " << fc << " FWHM = " << CWT::pick_fwhm( fc ) << "\n";
  cwt.load( &x );  
  cwt.run();

  
  const std::vector<double> & results = cwt.results(0);

  return results;

  // OR... do some smoothing?
  // const double moving_window_sec = 0.1;
  // int window_points = moving_window_sec * fs;
  // if ( window_points % 2 == 0 ) ++window_points;  
  // const std::vector<double> averaged 
  //   = MiscMath::moving_average( results , window_points );  
  // return averaged;
  
}
