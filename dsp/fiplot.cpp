
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

void fiplot_wrapper( edf_t & edf , const param_t & param )
{

  std::string signal_label = param.requires( "sig" );
  
  signal_list_t signals = edf.header.signal_list( signal_label );

  const int ns = signals.size();
  
  std::vector<double> Fs = edf.header.sampling_freq( signals );  

  interval_t interval = edf.timeline.wholetrace(); 
  
  //
  // Time (or cycle) bins
  //
  
  double t_lwr = param.has( "t-lwr" ) ? param.requires_dbl( "t-lwr" ) : 0.1 ;
  double t_upr = param.has( "t-upr" ) ? param.requires_dbl( "t-upr" ) : 4 ;
  double t_inc = param.has( "t-inc" ) ? param.requires_dbl( "t-inc" ) : 0.1  ;

  const bool cycles = param.has( "cycles" ) || param.has( "c-lwr" ); 

  if ( param.has( "c-lwr" ) ) t_lwr = param.requires_dbl( "c-lwr" );
  if ( param.has( "c-upr" ) ) t_upr = param.requires_dbl( "c-upr" );
  if ( param.has( "c-inc" ) ) t_inc = param.requires_dbl( "c-inc" );

  //
  // frequencies
  //

  // use either f-log {# of steps} OR f-inc {increment size}

  bool logspace = param.has("f-log"); // now interpret f-inc as the number of steps to have

  const double f_lwr = param.has( "f-lwr" ) ? param.requires_dbl( "f-lwr" ) : 1  ;
  const double f_upr = param.has( "f-upr" ) ? param.requires_dbl( "f-upr" ) : 20 ;
  const double f_inc = param.has( "f-log" ) ? param.requires_dbl( "f-log" ) : 
    ( param.has( "f-inc" ) ? param.requires_dbl( "f-inc" ) : 1 )  ;
  



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
		   t_lwr, t_upr, t_inc , cycles , 
	 	   f_lwr, f_upr, f_inc , logspace );
      
      
      //
      // Next signal
      //
      
      writer.unlevel( globals::signal_strat );
      
    } 

}


void fiplot_t::set_f( double lwr , double upr , double inc , bool logspace ) 
{
  
  frqs.clear();
  
  f_lwr = lwr;
  f_upr = upr;
  f_inc = inc; // # of inc if logspace == T
    
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
      
      // get CWT (12 cycles, fixed for now...)
      std::vector<double> c = cwt( x , fs , f , 7 );
      
      // get intervals
      fibin_t r = intervalize( c , tp , fs , t_lwr , t_upr , t_inc , cycles , f ); 
      
      // report
      std::map<double,fipair_t>::const_iterator rr = r.r.begin();
      while ( rr != r.r.end() )
	{
	  writer.level( rr->first , "TBIN" );
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
			       const bool cycles ,
			       const double fc )
{

  const double EPS = 0;
  
  fibin_t r;

  // number of time/cycle points
  int nt = ( t_upr - t_lwr ) / t_inc ;
  
  const double dt = 1.0/(double)fs;
  const uint64_t dt2 = globals::tp_1sec * ( dt * 1.5 );
  
  // point duration to time-point: sample-points are inclusive so +1
  // double time = ( j - i + 1 ) * dt  
  
  const int n = x_.size();

  // norm by min/max of 0/1
  // and also track discontinuities
  double min = x_[0];
  double max = x_[0];  

  // rule of thumb:
  // only consider intervals that are at least six times as long as the longest time span
  // considered
  // e.g. 5 seconds requires a full 30-second epoch, etc

  uint64_t required_tp = (uint64_t)6 * t_upr * globals::tp_1sec ; 
  
  int first_idx = -1;
  
  // T means the next point is not adjacent (or last point)
  std::vector<bool> disc( n , false ); 
  
  // track how much gets excluded
  double included_seconds = 0;
  double all_seconds = 0;

  for (int i=0;i<n;i++) 
    {
      
      if      ( x_[i] < min ) min = x_[i];
      else if ( x_[i] > max ) max = x_[i]; 
      
      if ( first_idx == -1 ) first_idx = i;
      
      // discontinuities
      if ( i == n-1 ) disc[i] = true;
      else if ( (*tp)[i+1] - (*tp)[i] > dt2 ) { disc[i] = true; } 
      
      // was this previous segment long enough anyway? 
      if ( disc[i] ) 
	{
	  
	  uint64_t length = i == n-1  // last point?
	    ? (*tp)[i] - (*tp)[first_idx] + ( (*tp)[1] - (*tp)[0] ) // add extra point
	    : (*tp)[i+1] - (*tp)[first_idx];
	  
	  all_seconds += length * globals::tp_duration ; 

	  //logger << "seg length " << length * globals::tp_duration << "\t";

	  //	  logger << required_tp * globals::tp_duration << "\n";
	  if ( length < required_tp ) 
	    {
	      //logger << " not long enough...\n";
	      for (int j=first_idx;j<=i;j++) disc[j] = true;
	    }
	  else
	    included_seconds +=  length * globals::tp_duration ;
	  first_idx = -1;
	}
    }
  
  logger << "  including " << ( included_seconds / all_seconds ) * 100 << "% of " << all_seconds << " seconds\n";
  
  // normalize
  std::vector<double> x( n , 0 );
  for (int i=0;i<n;i++) 
    {
      std::cout << "x\t" << i << "\t" << fc << "\t" << x_[i] ;
      //      x[i] = ( x_[i] - min ) / ( max - min );
      x[i] = x_[i];
      std::cout << "\t" << x[i] << "\n";
    }

  // init. each time-bin
  for (int t=0;t<nt;t++) 
    {
      // score as midpoint of each time-bin
      r.r[ t_lwr + t * t_inc + 0.5 * t_inc ] = fipair_t(0,0);  
    }
  

  
  // for each point, see how far we can go before we find a lower
  // point; store that time value, plus position plus height unless we
  // come to a discontinuity

  // fipoint_t (i,h,t)  

  
  std::set<fipoint_t> pts;
  for (int i=0;i<n;i++)
    {
      
      // min threshold at which we ignore
      if ( x[i] <= EPS ) continue;

      // for now, simply skip if we hit a discontinuity
      if ( disc[i] ) continue;

      for (int j=i+1;j<n;j++)
	{

	  // stop: original
	  // if ( x[j] <= x[i] ) 
	  if ( x[j] < x[i] ) 
	    {
	      // original
	      // fipoint_t f( i , j-1 , x[i] ) ;
	      
	      // new: go right up to last point
	      // may want to halve difference is a lot less than 
	      fipoint_t f( i , j , x[i] ) ;
	      pts.insert( f ) ;
	      break;
	    }
	  
	  // skip if next point is too far, but do not save
	  // might as well update to this point too

	  if ( disc[j] ) 
	    {
	      i = j+1;
	      break;	      
	    }
	}
    }
    
  // last point: ignore
  // fipoint_t f( n-1 , n-1 , x[n-1] ) ;
  // pts.insert( f ) ;

  
  
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


  std::vector<double> used( n , 0 );
  
  std::set<fipoint_t>::const_iterator ff = pts.begin();
  while ( ff != pts.end() ) 
    {
      
      // time duration ( seconds )
      double t = (double)( ff->j - ff->i + 1 ) * dt;
      
      // or in cycles?
      if ( cycles ) t *= fc;
      
      
//       logger << "t = " << t << "\t"
// 		<< ff->i << ".." << ff->j << " " << dt << "\t";
      
      // outside of range? 
      // if ( t < t_lwr ) logger << "below\n";
      // else if ( t >= t_upr )  logger << "above\n";

      if ( t >= t_lwr && t < t_upr )
	{
	  
	  int bin = ( t - t_lwr ) / t_inc; 
	  
	  double tbin = t_lwr + bin * t_inc + 0.5 * t_inc;  
	  
	  //	  logger << "bin = " << bin << " " << tbin << "\n";
	  
	  // add value to bin
	  double add = 0;
	  double amt = 0;
	  for (int i= ff->i ; i <= ff->j ; i++ ) 
	    {
	      double part = ff->h - used[i];
	      //logger << "  i " << i << " " << x[i] << " " << used[i] << " " << part << "\n";
	      add += part;
	      used[i] += part;	      
	      ++amt;
	    }
	  
	  // store 
	  r.r[ tbin ].w += add;
	  r.r[ tbin ].n += amt;
	  
	}

      // bin       
      ++ff;
    }
  
  // at end, x should equal used
  //for (int i=0;i<n;i++) std::cout << "x/u " << x[i] << "\t" << used[i] << "\n";  
  
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
  const double moving_window_sec = 0.1;
  int window_points = moving_window_sec * fs;
  if ( window_points % 2 == 0 ) ++window_points;
	  
  const std::vector<double> averaged 
    = MiscMath::moving_average( results , window_points );
  
  return averaged;
  
}
