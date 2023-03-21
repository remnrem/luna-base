
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

#include "stats/Eigen/Dense"
#include "stats/eigen_ops.h"

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
  // create new signals?
  // prefix_CH_N ... where N = 1,2,3, that correspond to Fs in range
  //
  
  const bool new_sigs = param.has( "add" ) ;
  
  std::string new_sig_prefix = new_sigs ? param.value( "add" ) : "" ; 

  
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
  // Start/stop times?
  //

  const bool   restrict_start = param.has( "start" );
  const bool   restrict_stop = param.has( "stop" );
  const double restrict_start_sec = param.has( "start" ) ? param.requires_dbl( "start" ) : 0 ;
  const double restrict_stop_sec = param.has( "stop" ) ? param.requires_dbl( "stop" ) : 0 ;
  
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
      const uint64_t delta_tp = globals::tp_1sec / Fs[s] ;

      //
      // Get time points (and flags for segments that span discontinuities)
      //  - also, indicate whether all should be computed (at segment/eppch level)
      //  - this is for moonlight MTM interactive viewer mainly
      
      const std::vector<uint64_t> * tp = slice.ptimepoints();	   
      const int np = tp->size();
      
      std::vector<double> start, stop;
      std::vector<int> start_sp, stop_sp; // original signal encoding
      
      // when 'add' option, count number of segments spanning this sample point
      std::vector<int> addn( np , 0 ); 

      // get 
      int nf = -1;
      Eigen::MatrixXd addX;

      std::vector<bool> disc, restrict;
      
      int p = 0;
      int nn = 0;
      int actual = 0;

      // actual segment size (may differ from requested due to sample rates)
      // here, check against this)
      const double segment_sec2 = segment_size / (double)Fs[s];

      while ( 1 ) {
	if ( p + segment_size > np ) break;
	
	double start_sec = (*tp)[p] * globals::tp_duration;
	double stop_sec = ( (*tp)[ p + segment_size - 1 ] + delta_tp ) * globals::tp_duration; // '1past'
	double implied_sec = stop_sec - start_sec;
	
	start_sp.push_back( p );
	stop_sp.push_back( p + segment_size - 1 ) ; // 'last point in seg'
	start.push_back( start_sec );
	stop.push_back( stop_sec );
	disc.push_back( fabs( implied_sec - segment_sec2 ) > 0.0001 );

	bool okay = true;
	if ( restrict_start && start_sec < restrict_start_sec ) okay = false;
	if ( restrict_stop && stop_sec > restrict_stop_sec ) okay = false;
	restrict.push_back( ! okay );
	if ( okay ) ++actual;
	++nn;
	// std::cout << "seg " << nn << "\t" << p << "\t" << start_sec << "\t" << stop_sec << "\t"
	// 	  << " sz " << implied_sec << " " << segment_sec2 << " " 
	// 	  << ( fabs( implied_sec - segment_sec2 ) > 0.001 )
	//   	  << "\trestrict=" << ! okay << "\n";
	
	// next segment
	p += segment_step;
	
      }


      if ( actual == 0 )
	{
	  logger << "  *** no segments to process, leaving MTM...\n";
	  return;
	}

      //
      // call MTM
      //
      
      mtm_t mtm( npi , nwin );
      
      mtm.dB = dB;
      mtm.opt_remove_mean = mean_center;
      mtm.opt_remove_trend = remove_linear_trend;

      if ( restrict_start || restrict_stop )
	mtm.restrict = restrict;
	  
      // s==0 means only give verbose output on first channel
      mtm.apply( d , Fs[s] , segment_size , segment_step , s == 0 );
      
      if ( s == 0 ) logger << "  processed channel(s):";
      logger << " " << signals.label(s) ;
      

      //
      // count freq bins if not already done? could add only areduced set perhaps?
      //

      if ( new_sigs && nf == -1 ) 
	{	  
	  nf = 0;
	  for ( int i = 0 ; i < mtm.f.size() ; i++ )
	    if ( mtm.f[i] >= min_f && mtm.f[i] <= max_f )
	      ++nf;
	  addX = Eigen::MatrixXd::Zero( np , nf );	  
	}
      
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
      
      if ( epoch_level_output || spectral_slope || new_sigs )
	{
	  const int nsegs = mtm.espec.size();

	  if ( nsegs != start.size() )
	    Helper::halt( "internal error in MTM timing:" + Helper::int2str( nsegs ) + " vs " + Helper::int2str( (int)start.size() ) );
	  
	  if ( epoch_level_output || new_sigs ) 
	    {
	      for ( int j = 0 ; j < nsegs ; j++)
		{
		  
		  if ( ! restrict[j] ) 
		    {
		      
		      //
		      // add main output
		      //
		      
		      if ( epoch_level_output ) 
			{
			  writer.level( j+1 , "SEG" );	  
			  writer.value( "START" , start[j] );
			  writer.value( "STOP" , stop[j] );
			  writer.value( "DISC" , (int)disc[j] );
			  
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
		      
		      //
		      // make new signals?		      
		      //
		      
		      if ( new_sigs ) 
			{
			  int s1 = start_sp[j];
			  int s2 = stop_sp[j];
			  //std::cout << " s1 s2 = " << s1 << " .. " << s2 << "\n";
			  for (int p=s1; p<=s2; p++)
			    {
			      addn[p]++;
			      
			      int fidx = 0;
			      for ( int i = 0 ; i < mtm.f.size() ; i++ )
				if ( mtm.f[i] >= min_f && mtm.f[i] <= max_f )
				  {
				    addX(p,fidx++) += mtm.espec[j][i];
				    //std::cout << " adding seg << " << j << " sample " << p << " freq " << fidx-1 <<"=" << mtm.f[i] << " = " << mtm.espec[j][i] << "\n";
				  }
			    }
			  
			}
		    }
		}
	    }
	  
	  // epoch level spectral slope? (based on the raw power )
	  
	  if ( spectral_slope )
	    {		  
	      
	      for ( int j = 0 ; j < nsegs ; j++)
		{

		  if ( ! restrict[j] )
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

      
      //
      // add new signals?
      //

      if ( new_sigs )
	{

	  int fidx = 0;
	  for ( int i = 0 ; i < mtm.f.size() ; i++ )
	    if ( mtm.f[i] >= min_f && mtm.f[i] <= max_f )
	      {
		
		const std::string new_sig_label = new_sig_prefix + "_" + signals.label(s) + "_" + Helper::int2str( fidx + 1 );
		
		std::vector<double> dat = eigen_ops::copy_array( addX.col( fidx ) );
		
		// normalize
		for (int p=0; p<np; p++) dat[p] /= (double)addn[p];
		
		if ( dat.size() != np ) Helper::halt( "internal error in MTM 'add'" );
		
		logger << "  adding new signal " << new_sig_label << " ( MTM @ " << mtm.f[i] << " Hz )\n";
		
		edf.add_signal( new_sig_label , Fs[s] , dat );
		
		++fidx;
	      }
	}
      
      
    } // next signal

  logger << "\n";
  
  writer.unlevel( globals::signal_strat );
  
    
}



