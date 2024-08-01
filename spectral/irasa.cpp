
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


#include "spectral/irasa.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "dsp/resample.h"
#include "fftw/fftwrap.h"
#include "miscmath/qdynam.h"

#include "helper/helper.h"
#include "db/db.h"

extern writer_t writer;
extern logger_t logger;


void irasa_wrapper( edf_t & edf , param_t & param )
{

  //
  // Get signals
  //

  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) );
  
  if ( signals.size() == 0 ) return;

  const int ns = signals.size();
  
  std::vector<double> Fs = edf.header.sampling_freq( signals );
  
  //
  // Analysis parameters
  //

  const bool silent = param.has( "silent" );

  const double h_min = param.has( "h-min" ) ? param.requires_dbl( "h-min" ) : 1.05;
  const double h_max = param.has( "h-max" ) ? param.requires_dbl( "h-max" ) : 1.95;
  const int    h_cnt = param.has( "h-steps" ) ? param.requires_dbl( "h-steps" ) : 19;
  
  const double f_lwr = param.has( "min" ) ? param.requires_dbl( "min" ) : 1 ;
  const double f_upr = param.has( "max" ) ? param.requires_dbl( "max" ) : 30 ;
  
  const bool average_adj = false;
  const double segment_sec = param.has( "segment-sec" ) ? param.requires_dbl( "segment-sec" ) : 4 ;
  const double overlap_sec = param.has( "segment-overlap" ) ?  param.requires_dbl( "segment-overlap" ) : 2 ;

  const bool calc_dynamics = param.has( "dynam" );
  
  const bool logout = param.has( "dB" );
  const bool epoch_lvl_output = param.has( "epoch" );
  
  window_function_t window_function = WINDOW_HAMMING;	   
  if      ( param.has( "no-window" ) ) window_function = WINDOW_NONE;
  else if ( param.has( "hann" ) ) window_function = WINDOW_HANN;
  else if ( param.has( "hamming" ) ) window_function = WINDOW_HAMMING;
  else if ( param.has( "tukey50" ) ) window_function = WINDOW_TUKEY50;
  
  const bool segment_median = ! param.yesno( "segment-mean" );
  const bool epoch_median = ! param.yesno( "epoch-mean" ); 
  
  const int converter = param.has( "fast" ) ? SRC_LINEAR : SRC_SINC_FASTEST ;

  std::vector<double> slope_range(2);
  slope_range[0] = f_lwr;
  slope_range[1] = f_upr;
  const double slope_outlier = 2 ; 

  const double fmin = f_lwr / h_max;
  const double fmax = f_upr * h_max;

  

  logger << "  specified frequency range is " << f_lwr << " - " << f_upr << " Hz\n";
  logger << "  full evaluated frequency range given h_max = " << h_max
	 << " is " << fmin << " - " << fmax << " Hz\n";

  bool problem = false;
  for (int s=0; s<ns; s++)
    {
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;
      if ( fmax > Fs[s] / 2.0 )
	{
	  logger << "  for " << signals.label(s) << ", Nyquist = " << Fs[s] / 2.0
		 << " Hz is less than implied upper evaluated of " << h_max << " * " << f_upr << " = " << fmax << " Hz\n";
	  problem = true;
	}
      
    }

  if ( problem )
    logger << "  *** warning *** evaluated frequency range exceeds Nyquist for one or more signals\n";
  

  //
  // Caching
  //

  const bool cache_data = param.has( "cache" );

  const std::string cache_name = cache_data ? param.requires( "cache" ) : "" ;

  const bool cache_epochs = param.has( "cache-epochs" );

  cache_t<double> * cache = NULL ;

  if ( cache_data )
    cache = edf.timeline.cache.find_num( cache_name );

  //
  // Iterate over signals
  //

  for (int s=0; s<ns; s++)
    {

      //
      // Skip non-data channels
      //

      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;

      writer.level( signals.label(s) , globals::signal_strat );

      //
      // Get data
      //

      slice_t slice( edf , signals(s) , edf.timeline.wholetrace() );
      
      const std::vector<double> * d = slice.pdata();

      const int ne = edf.timeline.first_epoch();
	    
      //
      // analysis 
      //

      irasa_t irasa( edf , param, *d , Fs[s] , edf.timeline.epoch_length(), ne, h_min, h_max, h_cnt , f_lwr, f_upr ,
		     segment_sec , overlap_sec , converter , epoch_lvl_output , logout , slope_range , slope_outlier ,
		     window_function , segment_median , epoch_median , cache , cache_epochs , silent , calc_dynamics );
      

      //
      // output
      //
      
      for (int f=0; f<irasa.n; f++)
	{
	  writer.level( irasa.frq[f] , globals::freq_strat );

	  if ( ! silent ) 
	    {
	      if ( logout )	    
		writer.value( "LOGF" , log( irasa.frq[f] ) );
	      
	      writer.value( "APER" , irasa.aperiodic[f] );
	      writer.value( "PER" , irasa.periodic[f] );	      
	    }
	  
	  if ( cache_data )
	    {
	      cache->add( ckey_t( "APER" , writer.faclvl() ) , irasa.aperiodic[f] );
	      cache->add( ckey_t( "PER" , writer.faclvl() ) , irasa.periodic[f] );
	    }
	  
	}
      writer.unlevel( globals::freq_strat );

      //
      // spectral slope?
      //
      
      bool okay = spectral_slope_helper( irasa.aperiodic_raw , 
					 irasa.frq ,
					 slope_range ,
					 slope_outlier ,
					 true );
      
      // next signal
    }

  writer.unlevel( globals::signal_strat );

}


irasa_t::irasa_t( edf_t & edf ,
		  param_t & param , 
		  const std::vector<double> & d ,
		  const int sr ,
		  const double epoch_sec,
		  const int ne, 
		  const double h_min ,
		  const double h_max ,
		  const int h_cnt ,
		  const double f_lwr,
		  const double f_upr ,
		  const double segment_sec , 
		  const double overlap_sec ,
		  const int converter, 
		  const bool epoch_lvl_output ,
		  const bool logout ,
		  const std::vector<double> & slope_range , 
		  const double slope_outlier ,
		  const int window_function ,
		  const bool segment_median ,
		  const bool epoch_median , 
		  cache_t<double> * cache , 
		  const bool cache_epochs , 
		  const bool silent ,
		  const bool calc_dynamics )
{
  
  const double h_inc = ( h_max - h_min ) / (double)(h_cnt-1);
  
  const int orig_epoch_smps = sr * epoch_sec;

  const int segment_points = segment_sec * sr;
  
  const int noverlap_points  = overlap_sec * sr;

   
  //
  // Get resampled versions of channels
  //
    
  std::vector<std::vector<double> > up, down;
  std::vector<int> up_epoch_smps, down_epoch_smps;
  
  for (int hi=0; hi<h_cnt; hi++)
    {
      const double h = h_min + hi * h_inc;	  
      //logger << "  creating resampled signals for " << signals.label(s) << " h = " << h << "\n";
      
      up.push_back( dsptools::resample( &d , sr , sr * h , converter ) );
      down.push_back( dsptools::resample( &d , sr , sr / h , converter ) );
      
      const int up_smps = up[ up.size() - 1 ].size();
      const int down_smps = down[ down.size() - 1 ].size();
      up_epoch_smps.push_back( up_smps / ne );
      down_epoch_smps.push_back( down_smps / ne );	  
    }      
  

   
  //
  // track epoch level stats, to get mean/median at the end
  //
  
  std::vector<std::vector<double> > apers, apers_raw, pers;

  //
  // dynamics?
  //

  qdynam_t qd;
  if ( calc_dynamics )
    qd.init( edf , param ); 
    
  //
  // Process epoch-wise
  //

  edf.timeline.first_epoch();


  
  for (int ec = 0; ec < ne ; ec++)
    {
      
      int epoch = edf.timeline.next_epoch();
      
      if ( epoch == -1 )
	Helper::halt( "internal error in irasa_t() - we've lost track of epoch counts" );

      // get original 
      std::vector<double> x( orig_epoch_smps );
      for (int i=0; i<orig_epoch_smps; i++)
	x[i] = d[ ec * orig_epoch_smps + i ] ;
      
      MiscMath::centre( x );
      
      const int total_points = orig_epoch_smps;
      
      // implied number of segments                                                                                                                                        
      const int noverlap_segments = floor( ( total_points - noverlap_points)
					   / (double)( segment_points - noverlap_points ) );

      PWELCH pwelch( x ,
		     sr, 
		     segment_sec ,
		     noverlap_segments ,
		     (window_function_t)window_function ,
		     segment_median );      
      
      std::vector<std::vector<double> > updowns( h_cnt );
      
      //
      // Up/down-sampled versions
      //
      
      for (int hi=0; hi<h_cnt; hi++)
	{
	  const double h = h_min + hi * h_inc;

	  const std::vector<double> & hup = up[ hi ];
	  const std::vector<double> & hdown = down[ hi ];
	      
	  const int up_smps = up_epoch_smps[ hi ];
	  const int down_smps = down_epoch_smps[ hi ];
	  
	  std::vector<double> up1( up_smps );
	  for (int i=0; i<up_smps; i++)
	    up1[i] = hup[ ec * up_smps + i ];
	  
	  std::vector<double> down1( down_smps );
	  for (int i=0; i<down_smps; i++)
	    down1[i] = hdown[ ec * down_smps + i ];
	  
	  //
	  // up
	  //

	  
	  const int up_noverlap_segments = floor( ( up_smps - noverlap_points )
						  / (double)( segment_points - noverlap_points ) );
	  
	  MiscMath::centre( up1 );
	  
	  PWELCH up_pwelch( up1 ,
			    sr, 
			    segment_sec ,
			    noverlap_segments ,
			    (window_function_t)window_function ,
			    segment_median );
	  
	      
	  //
	  // down
	  //

	  const int down_noverlap_segments = floor( ( down_smps - noverlap_points )
						    / (double)( segment_points - noverlap_points ) );
	  
	  MiscMath::centre( down1 );
	  
	  PWELCH down_pwelch( down1 ,
			      sr, 
			      segment_sec ,
			      noverlap_segments ,
			      (window_function_t)window_function ,
			      segment_median );
	  
	  

	  //
	  // collate geometric means (for freq range only)
	  //
	  
	  std::vector<double> ud;
	  for (int i=0; i<pwelch.psd.size() ; i++)
	    ud.push_back( sqrt( up_pwelch.psd[i] * down_pwelch.psd[i] ) );
	  updowns[ hi ] = ud;
	  
	}

      //
      // track size of FFT
      //

      if ( frq.size() == 0  )
	{
	  frq.clear();

	  for (int i=0; i<pwelch.psd.size() ; i++)
	    if ( pwelch.freq[ i ] >= f_lwr && pwelch.freq[ i ] <= f_upr )
	      frq.push_back( pwelch.freq[ i ] );

	  n = frq.size();
	  
	  periodic.resize( n );
	  aperiodic.resize( n );
	  aperiodic_raw.resize( n );

	  apers.resize( n );
	  apers_raw.resize( n );
	  pers.resize( n );
	}
      
      //
      // take median for each frequency
      //

      int cnt = 0;

      // verbose, epoch level output?
      if ( epoch_lvl_output || cache_epochs )
	writer.epoch( edf.timeline.display_epoch( epoch ) );

      // get main statistics for this epoch
      std::vector<double> aper_spectrum, aper_frq;

      for (int i=0; i<pwelch.psd.size() ; i++)
	{
	  if ( pwelch.freq[ i ] >= f_lwr && pwelch.freq[ i ] <= f_upr )
	    {
	      
	      std::vector<double> du( h_cnt );
	      for (int hi=0; hi<h_cnt; hi++) du[hi] = updowns[hi][i];

	      double aper = MiscMath::median( du , true ) ;
	      double per = pwelch.psd[ i ] - aper ;
	      
	      const bool okay = aper > 0 &&  pwelch.psd[ i ] > 0 ;

	      double log_aper, log_per;
	      
	      if ( ( calc_dynamics || logout ) && okay )
		{
		  log_aper = 10 * log10( aper );
                  log_per = 10 * log10( pwelch.psd[ i ] ) - log_aper;
		}
	      
	      // verbose, epoch level output? (or pass to qdynam_t?)
	      if ( epoch_lvl_output || cache_epochs || calc_dynamics )
		{		            
		  
		  writer.level( pwelch.freq[ i ] , globals::freq_strat );
		  
		  if ( epoch_lvl_output ) 
		    {
		      
		      // for epoch-level slope (below) [ always raw PSD ]
		      aper_frq.push_back( pwelch.freq[ i ] );
		      aper_spectrum.push_back( aper );
		  
		      if ( ! silent ) 
			{
			  if ( logout )
			    {
			      if ( okay )
				{
				  writer.value( "PER" , log_per );
				  writer.value( "APER" , log_aper );			    
				}
			    }
			  else
			    {
			      writer.value( "PER" , per );
			      writer.value( "APER" , aper );
			    }
			}
		    }

		  //
		  // dynamics?
		  //

		  if ( calc_dynamics && okay )
		    {
		      const int e = edf.timeline.display_epoch( epoch ) - 1;
		      qd.add( writer.faclvl_notime() , "APER" , e  , log_aper );
		      qd.add( writer.faclvl_notime() , "PER" , e  , log_per );
		    }
		  
		  //
		  // add epoch level data to cache 
		  //

		  if ( cache_epochs ) 
		    {
		      cache->add( ckey_t( "APER" , writer.faclvl() ) , aper );
		      cache->add( ckey_t( "PER" , writer.faclvl() ) , per );
		    }

		}
	    	      

	      //
	      // track for average over all 
	      //
	      
	      if ( logout )
		{
		  if ( okay )
		    {
		      apers[ cnt ].push_back( log_aper );
		      apers_raw[ cnt ].push_back( aper );
		      pers[ cnt ].push_back( log_per );
		    }
		}
	      else
		{
		  apers[ cnt ].push_back( aper ); 
		  apers_raw[ cnt ].push_back( aper );
		  pers[ cnt ].push_back( per ) ;
		}
	      
	      ++cnt;
	    }
	}
      
      
      if ( epoch_lvl_output || cache_epochs || calc_dynamics )
	writer.unlevel( globals::freq_strat );
      

      //
      // spectral slope
      //
      
      if ( epoch_lvl_output ) 
	{
	  bool okay = spectral_slope_helper( aper_spectrum , 
					     aper_frq , 
					     slope_range ,
					     slope_outlier ,
					     true );	  
	}

            
    }
  
  if ( epoch_lvl_output )
    writer.unepoch();
  
  //
  // average
  //

  for (int i=0; i<n; i++)
    {
      periodic[i] = epoch_median ? MiscMath::median( pers[i] ) : MiscMath::mean( pers[i] );
      aperiodic[i] = epoch_median ? MiscMath::median( apers[i] ) : MiscMath::mean( apers[i] );
      aperiodic_raw[i] = epoch_median ? MiscMath::median( apers_raw[i] ) : MiscMath::mean( apers_raw[i] );
    }

  //
  // dynamics?
  //

  if ( calc_dynamics )
    qd.proc_all();
  
}




