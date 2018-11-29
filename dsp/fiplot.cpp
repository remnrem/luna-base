#include "fiplot.h"
#include "../miscmath/miscmath.h"
#include "../cwt/cwt.h"
#include "../fftw/fftwrap.h"

extern writer_t writer;


void fiplot_wrapper( edf_t & edf , const param_t & param )
{

  std::string signal_label = param.requires( "signal" );
  
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



void fiplot_t::proc( const std::vector<double> & x , const std::vector<uint64_t> * tp , const int fs )
{

  // time points (secs)
  nt = ( t_upr - t_lwr ) / t_inc;
  
  // frq points (Hz)
  nf = frqs.size();

  for ( int fi = 0 ; fi < nf ; fi++ )
    {
      
      double f = frqs[fi];

      writer.level( f , globals::freq_strat );
      
      // get CWT (7 cycles, fixed for now...)
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

	  //std::cerr << "seg length " << length * globals::tp_duration << "\t";

	  //	  std::cerr << required_tp * globals::tp_duration << "\n";
	  if ( length < required_tp ) 
	    {
	      //std::cerr << " not long enough...\n";
	      for (int j=first_idx;j<=i;j++) disc[j] = true;
	    }
	  else
	    included_seconds +=  length * globals::tp_duration ;
	  first_idx = -1;
	}
    }
  
  std::cerr << "including " << ( included_seconds / all_seconds ) * 100 << "% of " << all_seconds << "\n";
  
  // normalize
  std::vector<double> x( n , 0 );
  for (int i=0;i<n;i++) 
    x[i] = ( x_[i] - min ) / ( max - min );
  
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
      
      
//       std::cerr << "t = " << t << "\t"
// 		<< ff->i << ".." << ff->j << " " << dt << "\t";
      
      // outside of range? 
      // if ( t < t_lwr ) std::cerr << "below\n";
      // else if ( t >= t_upr )  std::cerr << "above\n";

      if ( t >= t_lwr && t < t_upr )
	{
	  
	  int bin = ( t - t_lwr ) / t_inc; 
	  
	  double tbin = t_lwr + bin * t_inc + 0.5 * t_inc;  
	  
	  //	  std::cerr << "bin = " << bin << " " << tbin << "\n";
	  
	  // add value to bin
	  double add = 0;
	  double amt = 0;
	  for (int i= ff->i ; i <= ff->j ; i++ ) 
	    {
	      double part = ff->h - used[i];
	      //std::cerr << "  i " << i << " " << x[i] << " " << used[i] << " " << part << "\n";
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
  //  for (int i=0;i<n;i++) std::cout << "x/u " << x[i] << "\t" << used[i] << "\n";  
  
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
