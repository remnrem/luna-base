
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

#include "mtm.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "eval.h"
#include "fftw/fftwrap.h"

#include "db/db.h"
#include "helper/helper.h"
#include "helper/logger.h"

extern writer_t writer;
extern logger_t logger; 


void mtm::wrapper( edf_t & edf , param_t & param )
{
  
  std::string signal_label = param.requires( "sig" );

  signal_list_t signals = edf.header.signal_list( signal_label );    

  std::vector<double> Fs = edf.header.sampling_freq( signals );
  
  const int ns = signals.size();

  // nb. for efficiency, MTM uses its own segmentation of signals, rather than the
  // general 'epoch' function;  it also gives an average spectrum at the end,
  // but this can provide epoch level output too
  
  // nb. this assumes a continuous EDF

  const double segment_size_sec = param.has( "segment-sec" ) ? param.requires_dbl( "segment-sec" ) : 30;
  const double segment_step_sec = param.has( "segment-inc" ) ? param.requires_dbl( "segment-inc" ) : segment_size_sec ;

  //
  // report epoch-level ?
  //
  
  bool epoch_level_output = param.has( "epoch" );

  bool display_tapers = param.has( "dump-tapers" );

  bool mean_center = param.has( "mean-center" );
  bool remove_linear_trend = param.has( "detrend" );
  if ( mean_center && remove_linear_trend )
    Helper::halt( "cannot specify both mean-center and detrend" );
  
  //
  // MTM parameters (tw or nw)
  //
  
  double npi = 3;
  if ( param.has( "nw" ) ) npi = param.requires_dbl( "nw" );
  else if ( param.has( "tw" ) ) npi = param.requires_dbl( "tw" );
  
  int nwin = param.has( "t" ) ? param.requires_int( "t" ) : 2*floor(npi)-1 ;

  //
  // Required minimum SR to attempt MTM
  //

  int min_sr = param.has( "sr" ) ? param.requires_int( "sr" ) : 0 ; 

  //
  // Reporting full spectrum? (default 0.5 to 25 Hz)
  //

  double min_f = param.has( "min" ) ? param.requires_dbl( "min" ) : 0.5; 
  double max_f = param.has( "max" ) ? param.requires_dbl( "max" ) : 25;  

  //
  // Spectral slope
  //
  
  const bool spectral_slope = param.has( "slope" );
  const std::vector<double> slope_range = param.dblvector( "slope" );
  const bool spectral_slope_show_epoch = param.has( "epoch-slope" ) || param.has( "slope-epoch" );
  
  if ( spectral_slope )
    {
      if ( slope_range.size() != 2 ||
           slope_range[0] >= slope_range[1] ||
           slope_range[0] <= 0 ||
           slope_range[1] <= 0 )
	Helper::halt( "expecting slope=lwr,upr" );
    }

  // outlier threshold to remove individual PSD points when calculating a single slope  
  const double slope_outlier = param.has( "slope-th" ) ? param.requires_dbl( "slope-th" ) : 3 ;

  // threshold to reove epochs when summarizing slopes over all epochs  
  const double slope_th2     = param.has( "slope-th2" ) ? param.requires_dbl( "slope-th2" ) : 3 ;

    
  // output
  
  bool dB = param.has( "dB" );

  //
  // start analysis
  //

  interval_t interval = edf.timeline.wholetrace();

  //
  // Get each signal
  //
  
  for (int s = 0 ; s < ns; s++ )
    {
      
      //
      // only consider data tracks
      //
      
      if ( edf.header.is_annotation_channel( signals(s) ) )
	continue;
      
      //
      // Min. required SR?
      //

      if ( min_sr && Fs[s] < min_sr )
	continue;
      
      //
      // Stratify output by channel
      //
      
      writer.level( signals.label(s) , globals::signal_strat );
      
      
      //
      // Get data
      //
      
      slice_t slice( edf , signals(s) , interval );
      
      const std::vector<double> * d = slice.pdata();	   

      //
      // Step size in sample-points
      //
      
      const int segment_size = Fs[s] * segment_size_sec;
      const int segment_step = Fs[s] * segment_step_sec;
      
      //
      // call MTM
      //
      
      mtm_t mtm( npi , nwin );
      
      mtm.dB = dB;
      mtm.opt_remove_mean = mean_center;
      mtm.opt_remove_trend = remove_linear_trend;
      
      // s==0 means only give verbose output on first channel
      mtm.apply( d , Fs[s] , segment_size , segment_step , s == 0 );
      
      if ( s == 0 ) logger << "  processed channel(s):";
      logger << " " << signals.label(s) ;
      
      //
      // Output: tapers?
      //

      if ( display_tapers )
	{	  
	  for( int i=0; i< mtm.tapers.rows(); i++)
	    {
	      writer.level( i+1 , "SP" );
	      for(int j=0; j< mtm.tapers.cols(); j++)
		{
		  writer.level( j+1 , "TAPER" );
		  writer.value( "W" , mtm.tapers(i,j) );
		}
	      writer.unlevel( "TAPER" );
	    }
	  writer.unlevel( "SP" );
      
	  for(int j=0; j< mtm.lam.size(); j++) 
	    {
	      writer.level( j+1 , "TAPER" );
	      writer.value( "LAMBDA" , mtm.lam[j] );	  
	    }
	  writer.unlevel( "TAPER" );

	}

      
      //
      // Output: averaged spectrum
      //

      for ( int i = 0 ; i < mtm.f.size() ; i++ ) 
	{
	  if ( mtm.f[i] >= min_f && mtm.f[i] <= max_f ) 
	    {
	      writer.level( mtm.f[i] , globals::freq_strat  );
	      writer.value( "MTM" , mtm.spec[i] );
	    }
	}
      writer.unlevel( globals::freq_strat );


      //
      // display spectral slope?
      //
      
      // if dB mode, then based on the mean of db-scaled power (where
      // we need to convert back to raw); otherwise, just use raw

      if ( spectral_slope ) 
	{		   
	  
	  std::vector<double> xx = dB ? mtm.spec : mtm.raw_spec;
	  
	  if ( dB ) 
	    for (int i=0; i<xx.size(); i++) 
	      xx[i] = pow( 10 , xx[i]/10.0) ;	      
	  
	  spectral_slope_helper( xx , 
				 mtm.f , 
				 slope_range , 
				 slope_outlier );
	}
      
      
      
      //
      // segment-wise output? (do not call 'epoch' as this is typically different
      // i.e. not using epoch encoding etc
      //

      // store spectral slope per epoch for this channel?                                                                                                               
      std::vector<double> slopes;
      
      if ( epoch_level_output || spectral_slope )
	{
	  const int nsegs = mtm.espec.size();
	  
	  if ( epoch_level_output ) 
	    {
	      for ( int j = 0 ; j < nsegs ; j++)
		{
		  writer.level( j+1 , "SEG" );	  
		  
		  for ( int i = 0 ; i < mtm.f.size() ; i++ ) 
		    {
		      if ( mtm.f[i] >= min_f && mtm.f[i] <= max_f ) 
			{
			  writer.level( mtm.f[i] , globals::freq_strat  );
			  writer.value( "MTM" , mtm.espec[j][i] );
			}
		    }
		  writer.unlevel( globals::freq_strat );	      
		}
	    }
	  
	  // epoch level spectral slope? (based on the raw power )
	  
	  if ( spectral_slope )
	    {		  
	      
	      for ( int j = 0 ; j < nsegs ; j++)
		{
		  double es1 = 0 ;
		  
		  bool okay = spectral_slope_helper( mtm.raw_espec[j] , 
						     mtm.f ,						     
						     slope_range ,
						     slope_outlier ,
						     spectral_slope_show_epoch ,
						     &es1 );
		  
		  if ( okay ) slopes.push_back( es1 );
		}
	    }
	  
	  
	}
      
      if ( epoch_level_output ) 
	writer.unlevel( "SEG" );
      
      

      //
      // spectral slope based on distribution of epoch-level slopes?
      //
      
      if ( spectral_slope ) 
	{
	  if ( slopes.size() > 2 )
	    {
	      std::vector<double> s2 = MiscMath::outliers( &slopes , slope_th2 );
	      double s_mean = MiscMath::mean( s2 );
	      double s_med  = MiscMath::median( s2 );
	      double s_sd   = MiscMath::sdev( s2 , s_mean );
	      writer.value( "SPEC_SLOPE_MN" , s_mean );
	      writer.value( "SPEC_SLOPE_MD" , s_med );
	      writer.value( "SPEC_SLOPE_SD" , s_sd );
	    }
	}

      
    } // next signal

  logger << "\n";
  
  writer.unlevel( globals::signal_strat );
  
    
}



