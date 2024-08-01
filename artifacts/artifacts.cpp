
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

#include "artifacts.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "annot/annot.h"
#include "pdc/pdc.h"
#include "miscmath/qdynam.h"

#include "helper/helper.h"
#include "helper/logger.h"
#include "eval.h"
#include "db/db.h"

#include "dsp/resample.h"
#include "dsp/mse.h"
#include "dsp/lzw.h"

#include "miscmath/miscmath.h"
#include "fftw/fftwrap.h"

#include <iostream>
#include <fstream>

extern writer_t writer;

extern logger_t logger;


annot_t * buckelmuller_artifact_detection( edf_t & edf , 
					   param_t & param , 
					   const std::string & signal_label , 
					   const double delta_threshold , 
					   const double beta_threshold , 
					   const double delta_lwr , 
					   const double delta_upr ,
					   const double beta_lwr , 
					   const double beta_upr ,
					   const std::string & filename )
{


  //
  // parameters
  //
  
  bool set_mask = ! param.has( "no-mask" );
  bool verbose = param.has("verbose") || param.has( "epoch" );
  
  //
  // attach signal
  //
  
  signal_list_t signals = edf.header.signal_list( signal_label );  

  const int ns = signals.size();

  //
  // Sampling frequencies
  //

  std::vector<double> Fs = edf.header.sampling_freq( signals );


  //
  // Point to first epoch (assume 30 seconds, but could be different)
  //
  
  int ne = edf.timeline.first_epoch();

 
  //
  // Store per-epoch power
  //

  std::vector<std::vector<double> > delta(ns);
  std::vector<std::vector<double> > beta(ns);

    
  //
  // for each each epoch 
  //

  int cnt = 0;
  std::vector<int> track_epochs;
  
  while ( 1 ) 
    {
      
      int epoch = edf.timeline.next_epoch();      
       
      if ( epoch == -1 ) break;
      
      track_epochs.push_back( epoch );
      
      //
      // Get data 
      //

      interval_t interval = edf.timeline.epoch( epoch );

      for ( int s=0; s<ns; s++ )
	{

	  //
	  // only consider data tracks
	  //
	  
	  if ( edf.header.is_annotation_channel( signals(s) ) ) continue;	  
	  
	  slice_t slice( edf , signals(s) , interval );
	  
	  std::vector<double> * d = slice.nonconst_pdata();

	  const std::vector<uint64_t> * tp = slice.ptimepoints();
	  
	  //
	  // Mean-centre 30-second window
	  //

	  MiscMath::centre( d );
            
	  //
	  // Apply PWELCH to this epoch
	  //

	  // aim to get 10 windows of 4 seconds in 30sec epoch
	  
	  int noverlap_segments = 10;         
	  int segment_size_sec = 4;
	  
	  PWELCH pwelch( *d , Fs[s] , segment_size_sec , noverlap_segments );
      
	  // track power bands     
	  
	  delta[s].push_back( pwelch.psdsum( DELTA ) );
	  beta[s].push_back( pwelch.psdsum( beta_lwr  , beta_upr ) );
	  

	} // next signal

    } // next epoch



  
  //
  // Report for each signal
  //
  
  std::vector<std::vector<double> > delta_average(ns);
  std::vector<std::vector<double> > beta_average(ns);
  
  for (int s=0;s<ns;s++)
    {

      //
      // only consider data tracks
      //
      
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;

      //
      // Output stratifier
      //

      writer.level( signals.label(s) , globals::signal_strat );


      //
      // Make running averages
      //
      
      delta_average[s] = MiscMath::moving_average( delta[s] , 15 );
      beta_average[s]  = MiscMath::moving_average( beta[s]  , 15 );
      
      int total = 0, altered = 0;

      for (int e=0;e<ne;e++)
	{
	  
	  double dfac = delta[s][e] / delta_average[s][e];
	  double bfac = beta[s][e] / beta_average[s][e];
	  
	  bool dmask = dfac > delta_threshold ;
	  bool bmask = bfac > beta_threshold ;
	  bool mask  = dmask || bmask;
	  
	  if ( verbose ) 
	    {
	      
	      writer.epoch( edf.timeline.display_epoch( track_epochs[e] ) );

	      writer.var( "DELTA" , "Delta power" );
	      writer.var( "DELTA_AVG" , "Local average delta power" );
	      writer.var( "DELTA_FAC" , "Relative delta power factor" );
	      writer.var( "BETA" , "Beta power" );
	      writer.var( "BETA_AVG" , "Local average beta power" );
	      writer.var( "BETA_FAC" , "Relative beta power factor" );
	      writer.var( "DELTA_MASK" , "Masked based on delta power" );
	      writer.var( "BETA_MASK" , "Masked based on beta power" );
	      writer.var( "MASK" , "Masked" );		

	      writer.value( "DELTA" , delta[s][e] );
	      writer.value( "DELTA_AVG" , delta_average[s][e] );
	      writer.value( "DELTA_FAC" , dfac );

	      writer.value( "BETA" , beta[s][e] );
	      writer.value( "BETA_AVG" , beta_average[s][e] );
	      writer.value( "BETA_FAC" , bfac );
	      
	      writer.value( "DELTA_MASK" , dmask );
	      writer.value( "BETA_MASK" , bmask );
	      writer.value( "MASK" , mask );

	    }
	  
	  //
	  // Mask this epoch?
	  //
	  
	  if ( set_mask && mask )
	    {
	      if ( !edf.timeline.masked(e) ) ++altered;
	      edf.timeline.set_epoch_mask( e );
	      ++total;
	    }


	  writer.unepoch();
	  
	}

      if ( set_mask )
	logger << " masked " << total << " of " << ne << " epochs, altering " << altered << "\n";
      
      // # masked (i.e. actually), # masked, # total
      // altered , 

      writer.var( "FLAGGED_EPOCHS" , "Number of epochs failing Buckelmueller" );
      writer.var( "ALTERED_EPOCHS" , "Number of epochs actually masked" ); 
      writer.var( "TOTAL_EPOCHS" , "Number of epochs tested" );
      
      writer.value( "FLAGGED_EPOCHS" , total );
      writer.value( "ALTERED_EPOCHS" , altered );
      writer.value( "TOTAL_EPOCHS" , ne );

      writer.unlevel( globals::signal_strat );
      
    } // next signal    
  
  

  //
  // For now, do not return any annot_t
  //
  
  return NULL;


  //
  // In future, we may routinely expand all functions to return annotations
  //
  
  annot_t * a = edf.timeline.annotations.add( "Buckelmuller" );
  
  if ( a == NULL ) std::cout << "is null;\n";

  a->description = "Buckelmuller et al (2006) automatic artifact detection" ;

  for (int s=0;s<ns;s++)
    {
      //
      // only consider data tracks
      //
      
      if ( edf.header.is_annotation_channel( signals(s) ) ) 
	continue;
      
      for (int e=0;e<ne;e++)
	{
	  double dfac = delta[s][e] / delta_average[s][e];
	  double bfac = beta[s][e] / beta_average[s][e];
	  
	  bool reject =  dfac > delta_threshold || bfac > beta_threshold;
	  
	  if ( reject ) 
	    a->add( "buckelmuller:" + signals.label(s)  , edf.timeline.epoch(e) , signals.label(s) );
	}
      
    }
  
  return a;
  
}



//
// SIGSTATS
//

void  rms_per_epoch( edf_t & edf , param_t & param )
{
  
  if ( param.has( "th" ) || param.has( "chep" ) || param.has( "cstats" ) || param.has( "astats" ) || param.has( "mask" ) ) 
    Helper::halt( "use CHEP-MASK to find channel/epoch outliers: SIGSTATS now only reports epoch-level/individual-level statistics" );

  // Hjorth parameters: H1, H2, H3
  // Second-order Hjorth
  // Optional: RMS, % clipped signals
  // Optional: permutation entropy
  // Optional: fractal dimensions
  
  std::string signal_label = param.requires( "sig" );  

  bool verbose = param.has( "verbose" ) || param.has( "epoch") ;

  bool calc_rms = param.has( "rms" );

  bool calc_clipped = param.has( "clipped" );

  bool calc_flat = param.has( "flat" );

  bool calc_maxxed = param.has( "max" );

  int required_sr = param.has( "sr-over" ) ? param.requires_int( "sr-over" ) : 0 ; 

  bool calc_pfd = param.has( "pfd" ) ;

  bool calc_dynamics = param.has( "dynam" );
  
  bool calc_pe = param.has( "pe" ) || param.has( "pe-m" ) || param.has( "pe-t" );
  std::vector<int> pe_m = { 3,4,5,6,7 };
  if ( param.has( "pe-m" ) ) pe_m = param.intvector( "pe-m" ) ;
  int pe_t = param.has( "pe-t" ) ? param.requires_int( "pe-t" ) : 1 ; 
  
  bool calc_hjorth2 = param.has( "hjorth2" );

  double hjorth2_win = param.has( "hjorth2-win" ) ? param.requires_dbl( "hjorth2-win" ) : 1 ;

  double hjorth2_inc = param.has( "hjorth2-inc" ) ? param.requires_dbl( "hjorth2-inc" ) : hjorth2_win ;
  if ( hjorth2_inc < 0 ) hjorth2_inc = 0;


  double flat_eps = 1e-6;
  if ( calc_flat )
    {
      if ( ! param.empty( "flat" ) ) flat_eps = param.requires_dbl( "flat" );
      logger << "  epsilon for flat signals: |X[i]-X[i-1]| < " << flat_eps << "\n";
    }
  
  double max_value = 0;
  if ( calc_maxxed )
    {
      max_value = param.requires_dbl( "max" );
      logger << "  reporting max proportion, |X| > " << max_value << "\n";
    }
  
  //
  // Optionally calculate turning rate
  //
  
  // bool turning_rate = param.has( "tr" )  || param.has( "tr-epoch" ) || param.has( "tr-d" ) || param.has( "tr-smooth" );

  // double tr_epoch_sec = 1.0;
  // int    tr_d = 4;
  // int    tr_epoch_smooth = 30; // +1 is added afterwards
  
  // if ( turning_rate ) 
  //   {
  //     if ( param.has( "tr-epoch" ) ) tr_epoch_sec = param.requires_dbl( "tr-epoch" );
  //     if ( param.has( "tr-d" ) ) tr_d = param.requires_int( "tr-d" );
  //     if ( param.has( "tr-smooth" ) ) tr_epoch_smooth = param.requires_int( "tr-smooth" );
  //     logger << " calculating turning rate: d="<< tr_d << " for " << tr_epoch_sec << "sec epochs, smoothed over " << tr_epoch_smooth << " epochs\n";
  //   }

  
  //
  // Attach signals
  //
  
  signal_list_t signals = edf.header.signal_list( signal_label );  

  const int ns_all = signals.size();

  
  // data channels (& optional required sampling-rate)
  std::vector<int> sdata;
  for (int s=0;s<ns_all;s++)
    if ( ! edf.header.is_annotation_channel( signals(s) ) )
      if ( required_sr == 0 || edf.header.sampling_freq( signals( s ) ) >= required_sr )	
	sdata.push_back(s);
  
  int ns = 0;
  for (int s=0;s<ns_all;s++)
    if ( ! edf.header.is_annotation_channel( signals(s) ) )
      if ( required_sr == 0 || edf.header.sampling_freq( signals( s ) ) >= required_sr ) 
	++ns;
  
  if ( ns == 0 ) return;

  //
  // Store per-epoch statistics
  //

  std::vector<int>    n( ns , 0 );

  std::vector<double> rms( ns , 0 );
  std::vector<double> clipped( ns , 0 );
  std::vector<double> flat( ns , 0 );
  std::vector<double> maxxed( ns , 0 );
  std::vector<double> mean_activity( ns , 0 );
  std::vector<double> mean_mobility( ns , 0 );
  std::vector<double> mean_complexity( ns , 0 );
  //  std::vector<double> mean_turning_rate( ns , 0 ); // nb. this is based on turning-rate epoch size

  //
  // dynamics
  //

  qdynam_t qd;

  if ( calc_dynamics )
    qd.init( edf , param ) ;

  
  //
  // Point to first epoch 
  //
  
  int ne = edf.timeline.first_epoch();

  if ( ne == 0 ) return;
  
  
  //
  // For each signal
  //

  int si = -1;

  for (int s=0;s<ns_all;s++)
    {

      //
      // only consider data tracks
      //
      
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;

      //
      // and any SR requirements
      //

      if ( required_sr != 0 && edf.header.sampling_freq( signals( s ) ) < required_sr ) continue;

      //
      // what's left here should match the stores defined above
      //

      ++si;

      //
      // output stratifier (only needed at this stage if verbose, epoch-level output will
      // also be written)
      //

      if ( verbose || calc_dynamics ) 
	writer.level( signals.label(s) , globals::signal_strat );

      //
      // reset to first epoch
      //
      
      edf.timeline.first_epoch();


      //
      // Get sampling rate
      //
    
      int sr = edf.header.sampling_freq( signals( s ) );
      
      //
      // for each each epoch 
      //
      
      while ( 1 ) 
	{
	  
	  //
	  // Get next epoch
	  //
	  
	  int epoch = edf.timeline.next_epoch();
	  
	  if ( epoch == -1 ) break;
	  
	  interval_t interval = edf.timeline.epoch( epoch );
	  	  
	  slice_t slice( edf , signals(s) , interval );
	  
	  std::vector<double> * d = slice.nonconst_pdata();

	  //
	  // get clipped, flat and/or maxxed points (each is a proportion of points in the epoch)
	  //

	  double c = calc_clipped ? MiscMath::clipped( *d ) : 0 ;

	  double f = calc_flat ? MiscMath::flat( *d , flat_eps ) : 0 ;
	  
	  double m = calc_maxxed ? MiscMath::max( *d , max_value ) : 0 ; 


	  //
	  // Mean-centre 30-second window, calculate RMS
	  //

	  MiscMath::centre( d );
	  
	  double x = calc_rms ? MiscMath::rms( *d ) : 0 ;
	  	  

	  //
	  // Permutation entropy
	  //

	  std::vector<double> pe( pe_m.size() );

	  if ( calc_pe )
	    {
	      for ( int p=0; p<pe_m.size(); p++)
		{
		  int sum1 = 1;
		  std::vector<double> pd = pdc_t::calc_pd( *d , pe_m[p] , pe_t , &sum1 );
		  pe[p] = pdc_t::permutation_entropy( pd );
		}
	    }

	  //
	  // Fractal dimension
	  //

	  double pfd = 0;
	  if ( calc_pfd )
	    pfd = MiscMath::petrosian_FD( *d );
	  
	  //
	  // Hjorth parameters
	  //

	  double activity = 0 , mobility = 0 , complexity = 0;
	  
	  MiscMath::hjorth( d , &activity , &mobility , &complexity );


	  //
	  // 'Second-order' Hjorth
	  //

	  double hjorth2[9];

	  if ( calc_hjorth2 && sr >= 50 )
	    {
	      MiscMath::hjorth2( d , &(hjorth2)[0] , hjorth2_win * sr , hjorth2_inc * sr );
	    }		    


	  //
	  // Turning rate
	  //
	  
	  // double turning_rate_mean = 0;

	  // if ( turning_rate )
	  //   {

	  //     std::vector<double> subepoch_tr;
	      
	  //     turning_rate_mean = MiscMath::turning_rate( d , sr, tr_epoch_sec , tr_d , &subepoch_tr );
	      
	  //     for (int i=0;i<subepoch_tr.size();i++)
	  // 	e_tr[s].push_back( subepoch_tr[i] );
	  //   }

	  //
	  // Store for dynamics
	  //

	  if ( calc_dynamics )
	    {
	      const int e = edf.timeline.display_epoch( epoch ) - 1;
	      qd.add( writer.faclvl_notime() , "H1" , e  , log1p( activity ) );
	      qd.add( writer.faclvl_notime() , "H2" , e  , mobility );
	      qd.add( writer.faclvl_notime() , "H3" , e  , complexity );	      
	    }
	  
	  //
	  // Verbose output
	  //
	  
	  if ( verbose )
	    {

	      writer.epoch( edf.timeline.display_epoch( epoch ) );

	      //
	      // Report calculated values
	      //
	      
	      writer.value( "H1" , activity );
	      writer.value( "H2" , mobility );
	      writer.value( "H3" , complexity );
	      
	      if ( calc_hjorth2 && sr >= 50 )
		{
		  writer.value( "H1H1" , hjorth2[0] );
		  writer.value( "H1H2" , hjorth2[1] );
		  writer.value( "H1H3" , hjorth2[2] );

		  writer.value( "H2H1" , hjorth2[3] );
		  writer.value( "H2H2" , hjorth2[4] );
		  writer.value( "H2H3" , hjorth2[5] );
		  
		  writer.value( "H3H1" , hjorth2[6] );
		  writer.value( "H3H2" , hjorth2[7] );
		  writer.value( "H3H3" , hjorth2[8] );
		}
	      
	      
	      if ( calc_rms )
		writer.value( "RMS" , x );

	      if ( calc_pe )
		for (int p=0;p<pe_m.size();p++)
		  writer.value( "PE" + Helper::int2str( pe_m[p] ) , pe[p] ) ;

	      if ( calc_pfd )
		writer.value( "PFD" , pfd );

	      if ( calc_clipped )
		writer.value( "CLIP" , c );

	      if ( calc_flat )
		writer.value( "FLAT" , f );

	      if ( calc_maxxed )
		writer.value( "MAX" , m );
	      
	      // if ( turning_rate ) 
	      // 	writer.value( "TR" , turning_rate_mean );
	    }


	  //
	  // Tot up for individual-level means
	  //

	  if ( calc_rms )
	    rms[si] += x;

	  if ( calc_clipped )
	    clipped[si] += c;
	  
	  if ( calc_flat )
	    flat[si] += f;

	  if ( calc_maxxed )
	    maxxed[si] += m;
	  
	  mean_activity[si] += activity;
	  mean_mobility[si] += mobility;
	  mean_complexity[si] += complexity;

	  n[si]   += 1;
		  
	  //
	  // Next epoch
	  //
	  
	} 

      if ( verbose ) 
	writer.unepoch();
      
      //
      // Next signal
      //
      
    }

  if ( verbose || calc_dynamics )
    writer.unlevel( globals::signal_strat );


  //
  // dynamics
  //

  if ( calc_dynamics )
    qd.proc_all();

      
  //
  // Turning rate sub-epoch level reporting (including smoothing over sub-epochs)
  //

  // if ( turning_rate )
  //   {
  //     Helper::halt( "turning-rate not implemented" );
  //     Helper::halt( "need to fix signal indexing, ns vs ns_all, si, as above");
  //     // i.e. to skip annotation channels
      
  //     for (int s=0;s<ns;s++)
  // 	{
  // 	  int sr = edf.header.sampling_freq( signals( s ) );
	  
  // 	  // how many units (in # of sub-epoch units);  +1 means includes self
  // 	  int winsize = 1 + tr_epoch_smooth / tr_epoch_sec ; 

  // 	  //	  logger << "sz = " << e_tr[s].size() << " " << winsize << "\n";
  // 	  e_tr[s] = MiscMath::moving_average( e_tr[s] , winsize );

  // 	  // output
  // 	  writer.level( signals.label(s) , globals::signal_strat );
  // 	  for (int i=0;i<e_tr[s].size();i++)
  // 	    {
  // 	      writer.level( i+1 , "SUBEPOCH" );
  // 	      writer.value( "TR" , e_tr[s][i] );
  // 	    }
  // 	  writer.unlevel( "SUBEPOCH" );
  // 	}
  //     writer.unlevel( globals::signal_strat );
  //   }
  

   
  //
  // Individual level summary
  //
  
  for (int si=0;si<ns;si++)
    {
      writer.level( signals.label(sdata[si]) , globals::signal_strat );
      
      // mean Hjorth parameters
      writer.value( "H1"   , mean_activity[si] / (double)n[si] );
      writer.value( "H2"   , mean_mobility[si] / (double)n[si] );
      writer.value( "H3"   , mean_complexity[si] / (double)n[si] );

      if ( calc_clipped )
	writer.value( "CLIP" , clipped[si] / (double)n[si] );

      if ( calc_flat )
	writer.value( "FLAT" , flat[si] / (double)n[si] );

      if ( calc_maxxed )
	writer.value( "MAX" , maxxed[si] / (double)n[si] );

      if ( calc_rms )
	writer.value( "RMS"   , rms[si] / (double)n[si] );
      
    }

  writer.unlevel( globals::signal_strat );

}




//
// CHEP-MASK: (1) fixed values only (no iterative/Hjorth procedures)
//

void  chep_mask_fixed( edf_t & edf , param_t & param )
{

  // set CHEP mask based on :
  // - % clipped signals
  // - flat signals
  // - max values
  
  // nb. this command only ever sets the masks; we do not unmasked good values
  bool calc_clipped = param.has( "clipped" );
  bool calc_flat = param.has( "flat" );
  bool calc_maxxed = param.has( "max" );
  bool calc_minmax = param.has( "min-max" );

  // nothing to do...
  if ( !( calc_clipped || calc_flat || calc_maxxed || calc_minmax ) ) return;
  
  // e.g. exclude EPOCH is more than 5% of points are clipped  
  double clip_threshold = calc_clipped ? param.requires_dbl( "clipped" ) : 0.05 ;
  if ( calc_clipped )
    logger << "  flagging epochs with " << clip_threshold << " proportion X[i] == max(X) or min(X)\n";
  
  double flat_threshold = 0.05;
  double flat_eps = 1e-6;
  if ( calc_flat )
    {
      std::vector<double> x = param.dblvector( "flat" );
      if ( x.size() != 1 && x.size() != 2 ) Helper::halt( "flat requires 1 or 2 param: flat=<prop>,<eps>" );
      flat_threshold = x[0];
      if ( x.size() == 2 ) flat_eps = x[1];
      logger << "  flagging epochs with " << flat_threshold << " proportion |X[i]-X[i-1]| < " << flat_eps << "\n";
    }
  
  double max_threshold = 0.05;
  double max_value = 0;
  if ( calc_maxxed )
    {
      std::vector<double> x = param.dblvector( "max" );
      if ( x.size() != 2 ) Helper::halt( "max requires 2 params: max=<value>,<prop>" );
      max_value = x[0];
      max_threshold = x[1];
      logger << "  flagging epochs with " << max_threshold << " proportion |X| > " << max_value << "\n";
    }
  
  // reject if MAX is *below* a set threshold
  const double minmax_threshold = calc_minmax ? param.requires_dbl( "min-max" ) : 0 ;
  if ( calc_minmax )
    {
      if ( minmax_threshold <= 0 )  Helper::halt( "expecting min-max to be > 0" );
      logger << "  flagging epochs with a max |X| less than " << minmax_threshold << "\n";
    }

  
  //
  // Attach signals
  //

  std::string signal_label = param.requires( "sig" );  
  
  signal_list_t signals = edf.header.signal_list( signal_label );  

  // count channels
  const int ns = signals.size();  
  
  if ( ns == 0 ) return;
  
  //
  // Point to first epoch 
  //
  
  const int ne = edf.timeline.first_epoch();

  if ( ne == 0 ) return;
  
  //
  // Track what we remove, all channels/epochs
  //

  int count_all = 0 , count_unmasked = 0 , count_masked = 0;
  
  //
  // For each signal
  //

  for (int s=0; s<ns; s++)
    {

      //
      // only consider data tracks
      //
      
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;

      //
      // reset to first epoch
      //
      
      edf.timeline.first_epoch();

      //
      // Get sampling rate
      //
    
      int sr = edf.header.sampling_freq( signals( s ) );
      
      //
      // Track what we remove
      //

      int cnt_clp = 0 , cnt_flt = 0 , cnt_max = 0 , cnt_minmax = 0 , cnt_any = 0;
      

      //
      // for each each epoch 
      //
      
      while ( 1 ) 
	{
	  
	  ++count_all;

	  //
	  // Get next epoch, which respects the epoch-level mask and CHEP mask
	  //
	  
	  int epoch = edf.timeline.next_epoch();
	  
	  if ( epoch == -1 ) break;
		  
	  if ( edf.timeline.masked( epoch , signals.label(s) ) ) continue;

	  ++count_unmasked;

	  //
	  // Get data
	  //
	  
	  interval_t interval = edf.timeline.epoch( epoch );
	  	  
	  slice_t slice( edf , signals(s) , interval );
	  
	  std::vector<double> * d = slice.nonconst_pdata();

	  //
	  // clipped, flat and/or maxxed points (each is a proportion of points in the epoch)
	  //

	  double c = calc_clipped ? MiscMath::clipped( *d ) : 0 ;

	  double f = calc_flat ? MiscMath::flat( *d , flat_eps ) : 0 ;

	  double m = calc_maxxed ? MiscMath::max( *d , max_value ) : 0 ; 
	  
	  // actual max(|X|)
	  double mxval = 0;
	  if ( calc_minmax )
	    {
	      double mx = 0 , mn = 0;
	      MiscMath::minmax( *d , &mn , &mx );
	      mxval = fabs(mx) > fabs(mn) ? fabs(mx) : fabs(mn) ;
	    }
	  
	  //
	  // Mask?
	  //

	  bool set_mask = false;
		  
	  // For clipping/flat/max, use a fixed threshold
	  if ( calc_clipped && c > clip_threshold ) 
	    {
	      set_mask = true;
	      ++cnt_clp;
	    }
	  	  
	  if ( calc_flat && f > flat_threshold ) 
	    {
	      set_mask = true;
	      ++cnt_flt;
	    }
	  
	  if ( calc_maxxed && m > max_threshold ) 
	    {	      
	      set_mask = true;
	      ++cnt_max;
	    }	  
	  
	  if ( calc_minmax && mxval < minmax_threshold )
	    {
	      set_mask = true;
              ++cnt_minmax;      
	    }
	  
	  if ( set_mask ) 
	    {	      
	      edf.timeline.set_chep_mask( epoch , signals.label(s) );	      
	      ++count_masked;
	      ++cnt_any;
	    }
	  
	  //
	  // Next epoch
	  //

	}

      //
      // Report signal level stats
      //

      
      logger << "  for " << signals.label(s) << ", clipped: " << cnt_clp 
	     << " flat: " << cnt_flt
	     << " max: " << cnt_max
	     << " min-max: " << cnt_minmax << "\n";

      // writer.level( signals.label(s) , globals::signal_strat );
      
      // // channel-level output about # of epochs being masked
      // if ( calc_clipped) writer.value( "CLP" , cnt_clp );
      // if ( calc_flat) writer.value( "FLT" , cnt_flt );
      // if ( calc_maxxed) writer.value( "MAX" , cnt_max );

      // // if >1 criteria, report intersection 
      // if ( calc_clipped + calc_flat + calc_maxxed > 1 ) 
      // 	writer.value( "ANY" , cnt_any );
	
      //
      // Next signal
      //
      
    } 

  //  writer.unlevel( globals::signal_strat );
  
  
  logger << "  masked " << count_masked << " epoch/channel pairs of "
	 << count_unmasked << " previously unmasked (" << count_all << " in total)\n";
  
}



//
// CHEP-MASK: (2) statistical Hjorth-based outlier detection
//

void  chep_mask( edf_t & edf , param_t & param )
{
  
  // Hjorth parameters: H1, H2, H3
  // respects any existing CHEP mask, and/or epoch-level mask
  
  std::string signal_label = param.requires( "sig" );  

  // within-channel, across epochs:  ep-th
  // within-epoch, across channels: ch-th
  // across-channels, across-epoochs: chep-th

  // by default, respects existing CHEP mask;  if in form ep-th0 instead, 
  // then any existing chep mask is ignore;  within type of outlier detection, 
  // iterative runs respect any mask (obviously).   Specifically, we
  //   make a copy of the CHEP
  //   clear CHEP
  //   run the command as is
  //   merge the orignial CHEP back into the new CHEP
  
  // all can be applied iterativly;  all respect the existing masks
  //  use CHEP epochs=X  
  //      CHEP drop-channels 
  //   to drop channels and/or set the epoch mask

  //
  // allow for iterative outlier detection, i.e. with multiple
  // comma-delimited thresholds
  //
  
  bool ep_ignore = param.has( "ep-th0" );
  bool ch_ignore = param.has( "ch-th0" );
  bool chep_ignore = param.has( "chep-th0" );

  if ( ep_ignore   && param.has( "ep-th" ) ) Helper::halt( "cannot specify both ep-th and ep-th0" );
  if ( ch_ignore   && param.has( "ch-th" ) ) Helper::halt( "cannot specify both ch-th and ch-th0" );
  if ( chep_ignore && param.has( "chep-th" ) ) Helper::halt( "cannot specify both chep-th and chep-th0" );

  std::vector<double> ep_th, ch_th, chep_th;

  if ( param.has( ep_ignore ? "ep-th0" : "ep-th" ) ) ep_th = param.dblvector( ep_ignore ? "ep-th0" : "ep-th"  );
  if ( param.has( ch_ignore ? "ch-th0" : "ch-th" ) ) ch_th = param.dblvector( ch_ignore ? "ch-th0" : "ch-th" );
  if ( param.has( chep_ignore ? "chep-th0" : "chep-th" ) ) chep_th = param.dblvector( chep_ignore ? "chep-th0" : "chep-th" );

  //
  // Attach signals
  //
  
  signal_list_t signals = edf.header.signal_list( signal_label );  

  // all channels
  const int ns_all = signals.size();
  
  // data channels (slot number)
  std::vector<int> sdata;
  for (int s=0;s<ns_all;s++)
    if ( ! edf.header.is_annotation_channel( signals(s) ) )
      sdata.push_back( signals(s) );
  
  int ns = sdata.size();
  
  if ( ns == 0 ) return;

  //
  // Track epoch-level statistics
  //
  
  std::vector<std::vector<double> > e_act(ns);
  std::vector<std::vector<double> > e_mob(ns);
  std::vector<std::vector<double> > e_cmp(ns);
  std::vector<std::vector<int> > e_epoch(ns);
  

  //
  // Point to first epoch 
  //
  
  int ne = edf.timeline.first_epoch();

  if ( ne == 0 ) return;
    
  //
  // For each signal
  //


  for (int si=0; si < ns; si++ )
    {

      const std::string clabel = edf.header.label[ sdata[si] ];

      //
      // reset to first epoch
      //
      
      edf.timeline.first_epoch();

      //
      // Get sampling rate
      //
      
      int sr = edf.header.sampling_freq( sdata[si] );

      //
      // for each each epoch 
      //
      
      while ( 1 ) 
	{
	  
	  //
	  // Get next epoch, which respects the epoch-level mask
	  //
	  
	  int epoch = edf.timeline.next_epoch();
	  
	  if ( epoch == -1 ) break;
	  
	  //
	  // Is already CHEP-masked?
	  //

	  if ( edf.timeline.masked( epoch , clabel ) )
	    {
	      // just insert zeros, as this will be skipped over downstream 
	      // in any case;  this way, we have all channels with the same number
	      // of epochs;   nb.  we can make e_epoch a single vector, as all e_epoch[si] should
	      // now be identical
	      e_act[si].push_back( -9 );
	      e_mob[si].push_back( -9 );
	      e_cmp[si].push_back( -9 );
	      e_epoch[si].push_back( epoch );
	      continue;
	    }

	  //
	  // Get data
	  //
	  
	  interval_t interval = edf.timeline.epoch( epoch );
	  
	  slice_t slice( edf , sdata[si] , interval );
	  
	  std::vector<double> * d = slice.nonconst_pdata();
	  

	  //
	  // Mean-centre 30-second window
	  //

	  MiscMath::centre( d );
	  
	  //
	  // Hjorth parameters
	  //

	  double activity = 0 , mobility = 0 , complexity = 0;
	  
	  MiscMath::hjorth( d , &activity , &mobility , &complexity );

	  //
	  // Track all epoch/channel level Hjorth values for outlier detection
	  //

	  e_act[si].push_back( activity );
	  e_mob[si].push_back( mobility );
	  e_cmp[si].push_back( complexity );

	  // track unmasked epochs (which may vary by channel)
	  e_epoch[si].push_back( epoch );
	  	  
	  //
	  // Next epoch
	  //

	} 

      //
      // Next signal
      //

    }

  
  //
  // Apply statistical masks to Hjorth paramters: e_act, e_mob and e_cmp values 
  //
  // 1) within-channel, between-epoch masking (ep-th)
  // 2) within-epoch, between channel masking (ch-th) 
  // 3) between-channel, between-epoch masking (chep-th)
  //
  // Apply iterative masks, in the above order

  // if at any step we are ignoring the prior mask, copy/clear/merge in using this:

  std::map<int,std::set<std::string> > chep_copy;

  if ( ep_th.size() > 0 ) 
    logger << "  within-channel/between-epoch outlier detection, ep-th" << ( ep_ignore ? "0" : "" ) << " = " 
	   << param.value( ep_ignore ? "ep-th0" : "ep-th" ) << "\n";
  

  //
  // 1) ep-th masks
  //


  if ( ep_ignore ) 
    {
      logger << "   (ignoring existing CHEP mask)\n"; 
      chep_copy = edf.timeline.make_chep_copy();
      edf.timeline.clear_chep_mask();
    }
  
  int total = 0;
  
  for (int o=0 ; o<ep_th.size(); o++ )
    {
      
      int total_this_iteration = 0;
      
      //
      // for each channel separately
      //
      
      for (int si = 0; si < ns; si++)
	{
	  
	  const std::string clabel = edf.header.label[ sdata[si] ];

	  int cnt_act = 0 , cnt_mob = 0 , cnt_cmp = 0;
	  
	  const int nepochs = e_epoch[si].size();
      
	  if ( nepochs < 3 ) continue;
      
	  std::vector<double> act_act;
	  std::vector<double> act_mob;
	  std::vector<double> act_cmp;
	  std::vector<int> act_epoch;

	  // calculate current mean for RMS and Hjorth parameters
	  // (clipping based on a single fixed threshold, not statistically)

	  for (int j=0;j<nepochs;j++)
	    {
	      if ( ! edf.timeline.masked( e_epoch[si][j] , clabel ) )
		{
		  act_act.push_back( e_act[si][j] );
		  act_mob.push_back( e_mob[si][j] );
		  act_cmp.push_back( e_cmp[si][j] );
		  act_epoch.push_back( e_epoch[si][j] );
		}	      
	    }
	
	  // actual number of currently included epochs
	  const int ne = act_epoch.size();
	  
	  // if down to too few epochs, skip
	  if ( ne < 3 ) continue;

	  const double mean_act = MiscMath::mean( act_act );
	  const double mean_mob = MiscMath::mean( act_mob );
	  const double mean_cmp = MiscMath::mean( act_cmp );
	  
	  const double sd_act = MiscMath::sdev( act_act , mean_act ); 
	  const double sd_mob = MiscMath::sdev( act_mob , mean_mob ); 
	  const double sd_cmp = MiscMath::sdev( act_cmp , mean_cmp ); 

	  const double this_th = ep_th[o];
	  
	  const double lwr_act = mean_act - this_th * sd_act;
	  const double lwr_mob = mean_mob - this_th * sd_mob;
	  const double lwr_cmp = mean_cmp - this_th * sd_cmp;

	  const double upr_act = mean_act + this_th * sd_act;
	  const double upr_mob = mean_mob + this_th * sd_mob;
	  const double upr_cmp = mean_cmp + this_th * sd_cmp;

	  for (int ei=0; ei<ne; ei++)
	    {
	      
	      bool set_mask = false;
		  		  
	      if ( act_act[ei] < lwr_act || act_act[ei] > upr_act ) 
		{
		  set_mask = true;
		  cnt_act++;
		}
		      
	      if ( act_mob[ei] < lwr_mob || act_mob[ei] > upr_mob ) 
		{
		  set_mask = true;
		  cnt_mob++;
		}
		      
	      if ( act_cmp[ei] < lwr_cmp || act_cmp[ei] > upr_cmp )
		{
		  set_mask = true;
		  cnt_cmp++;
		}
	      
	      if ( set_mask ) 
		{
		  const int e = act_epoch[ei];
		  edf.timeline.set_chep_mask( e , clabel ) ;
		  ++total_this_iteration;
		  ++total;
		}
	      
	      
	    } // next epoch	  
	  	  
	  
	} // next channel
      
      	  //
	  // report stats for this iteration 
	  //
      
      logger << "   iteration " << o+1 << ": removed " 
	     << total_this_iteration 
	     << " channel/epoch pairs this iteration (" << total << " in total)\n";
      
    
    } // next iteration

  
  if ( ep_ignore )
    {
      // add original mask back in
      edf.timeline.merge_chep_mask( chep_copy );
    }

  //
  // 2) ch-th masks
  //
  
  if ( ch_th.size() > 0 ) 
    logger << "  between-channel/within-epoch outlier detection, ch-th" << ( ch_ignore ? "0" : "" ) << " = " 
	   << param.value( ch_ignore ? "ch-th0" : "ch-th" ) << "\n";
  
  total = 0;
  
  if ( ch_ignore ) 
    {
      logger << "   (ignoring existing CHEP mask)\n"; 
      chep_copy = edf.timeline.make_chep_copy();
      edf.timeline.clear_chep_mask();
    }

  for (int o=0 ; o<ch_th.size(); o++ )
    {
      
      int total_this_iteration = 0;
      
      //
      // for each epoch separately
      //
      
      for (int ei=0; ei<e_epoch[0].size(); ei++) 
	{

	  // e_epoch is redundant... will be same across all channels... so just pick the first
	  
	  int epoch = e_epoch[0][ei];
	  
	  // track stats across comparable channels for this epoch (*ee)
	  
	  int cnt_act = 0 , cnt_mob = 0 , cnt_cmp = 0 ;
	  
	  std::vector<double> act_act;
	  std::vector<double> act_mob;
	  std::vector<double> act_cmp;
	  std::vector<std::string> act_clabel; // which channels actually used?
	  
	  // for each channel (not already masked)
	  for (int si = 0; si < ns; si++)
	    {

	      const std::string clabel = edf.header.label[ sdata[si] ];

	      if ( ! edf.timeline.masked( epoch , clabel ) )
		{
		  act_act.push_back( e_act[si][ ei ] );
		  act_mob.push_back( e_mob[si][ ei ] );
		  act_cmp.push_back( e_cmp[si][ ei ] );
		  act_clabel.push_back( clabel );
		}	      
	    }
	
	  // actual number of currently included channels
	  const int act_ns = act_clabel.size();
	  
	  // if down to too few signals, skip
	  if ( act_ns < 3 ) continue;

	  const double mean_act = MiscMath::mean( act_act );
	  const double mean_mob = MiscMath::mean( act_mob );
	  const double mean_cmp = MiscMath::mean( act_cmp );
	  
	  const double sd_act = MiscMath::sdev( act_act , mean_act ); 
	  const double sd_mob = MiscMath::sdev( act_mob , mean_mob ); 
	  const double sd_cmp = MiscMath::sdev( act_cmp , mean_cmp ); 

	  const double this_th = ch_th[o];
	  
	  const double lwr_act = mean_act - this_th * sd_act;
	  const double lwr_mob = mean_mob - this_th * sd_mob;
	  const double lwr_cmp = mean_cmp - this_th * sd_cmp;

	  const double upr_act = mean_act + this_th * sd_act;
	  const double upr_mob = mean_mob + this_th * sd_mob;
	  const double upr_cmp = mean_cmp + this_th * sd_cmp;
	  
	  for (int si=0; si < act_ns; si++)
	    {

	      bool set_mask = false;

	      if ( act_act[si] < lwr_act || act_act[si] > upr_act ) 
		{
		  set_mask = true;
		  cnt_act++;
		}
		      
	      if ( act_mob[si] < lwr_mob || act_mob[si] > upr_mob ) 
		{
		  set_mask = true;
		  cnt_mob++;
		}
		      
	      if ( act_cmp[si] < lwr_cmp || act_cmp[si] > upr_cmp )
		{
		  set_mask = true;
		  cnt_cmp++;
		}
	      
	      if ( set_mask ) 
		{
		  edf.timeline.set_chep_mask( epoch , act_clabel[si] ) ;
		  ++total_this_iteration;
		  ++total;
		}
	      
	    } // next channel
	  
	} // next channel
      
      	  //
	  // report stats for this iteration 
	  //
      
      logger << "   iteration " << o+1 << ": removed " 
	     << total_this_iteration 
	     << " channel/epoch pairs this iteration (" << total << " in total)\n";
      
      
    } // next iteration


  if ( ch_ignore )
    {
      // add original mask back in
      edf.timeline.merge_chep_mask( chep_copy );
    }

  
  //
  // 3) chep-th masks
  //
  
  if ( chep_th.size() > 0 ) 
    logger << "  between-channel/between-epoch outlier detection, chep-th" << ( chep_ignore ? "0" : "" ) << " = " 
	   << param.value( chep_ignore ? "chep-th0" : "chep-th" ) << "\n";
  
  total = 0 ;

  if ( chep_ignore ) 
    {
      logger << "   (ignoring existing CHEP mask)\n"; 
      chep_copy = edf.timeline.make_chep_copy();
      edf.timeline.clear_chep_mask();
    }

  for (int o=0 ; o<chep_th.size(); o++)
    {
      
      // standardize over all epochs and channels
      double h1_mean = 0 , h2_mean = 0 , h3_mean = 0;
      double h1_sd = 0 , h2_sd = 0 , h3_sd = 0;
      uint64_t cnt = 0;
	  
      for (int si=0;si<ns;si++)
	{
	  const std::string clabel = edf.header.label[ sdata[si] ];

	  for (int ei=0; ei< e_epoch[si].size(); ei++)
	    if ( ! edf.timeline.masked( e_epoch[si][ei] , clabel ) )
	      {
		h1_mean += e_act[si][ei];
		h2_mean += e_mob[si][ei];
		h3_mean += e_cmp[si][ei];
		++cnt;
	      }
	}

      h1_mean /= (double)cnt;
      h2_mean /= (double)cnt;
      h3_mean /= (double)cnt;
      
      for (int si=0;si<ns;si++)
	{
	  const std::string clabel = edf.header.label[ sdata[si] ];

	  for (int ei=0; ei< e_epoch[si].size(); ei++)
	    if ( ! edf.timeline.masked( e_epoch[si][ei] , clabel ) )
	      {
		h1_sd += ( e_act[si][ei] - h1_mean ) * ( e_act[si][ei] - h1_mean ) ;
		h2_sd += ( e_mob[si][ei] - h2_mean ) * ( e_mob[si][ei] - h2_mean ) ;
		h3_sd += ( e_cmp[si][ei] - h3_mean ) * ( e_cmp[si][ei] - h3_mean ) ;
	      }
	}

      
      h1_sd = sqrt( h1_sd / (double)(cnt-1) );
      h2_sd = sqrt( h2_sd / (double)(cnt-1) );
      h3_sd = sqrt( h3_sd / (double)(cnt-1) );

      double h1_lwr = h1_mean - chep_th[o] * h1_sd;
      double h2_lwr = h2_mean - chep_th[o] * h2_sd;
      double h3_lwr = h3_mean - chep_th[o] * h3_sd;
      
      double h1_upr = h1_mean + chep_th[o] * h1_sd;
      double h2_upr = h2_mean + chep_th[o] * h2_sd;
      double h3_upr = h3_mean + chep_th[o] * h3_sd;
	  
      // mask
      uint64_t masked = 0;
      for (int si=0;si<ns;si++)
	{
	  const std::string clabel = edf.header.label[ sdata[si] ];
	  for (int ei=0; ei< e_epoch[si].size(); ei++)
	    if ( ! edf.timeline.masked( e_epoch[si][ei] , clabel ) )
	      {
		bool set_mask = false; 
		if      ( e_act[si][ei] < h1_lwr || e_act[si][ei] > h1_upr ) set_mask = true;
		else if ( e_mob[si][ei] < h2_lwr || e_mob[si][ei] > h2_upr ) set_mask = true;
		else if ( e_cmp[si][ei] < h3_lwr || e_cmp[si][ei] > h3_upr ) set_mask = true;
		
		if ( set_mask ) 
		  {
		    edf.timeline.set_chep_mask( e_epoch[si][ei] , clabel );
		    ++masked;
		  }	      
	      }	    
	}

      // all done for this iteration
      logger << "  masked " << masked << " CHEPs of " << cnt << " unmasked CHEPs ("
	     << 100*(masked/(double(ns*ne))) << "%), from " << ne*ns << " total CHEPs, "
	     << "on iteration " << o+1 << "\n";
      
      
      //
      // Next iteration
      //
    }


  if ( chep_ignore )
    {
      // add original mask back in 
      edf.timeline.merge_chep_mask( chep_copy );
    }
  
  //
  // All done for CHEP-MASK based on Hjorth parameters
  //
}








void  mse_per_epoch( edf_t & edf , param_t & param )
{
  

  //
  // Calculate MSE per epoch and average
  //

  //
  // MSE parameters
  //

  int    m  = param.has("m") ? param.requires_int("m") : 2;  

  double r  = param.has("r") ? param.requires_dbl("r") : 0.15;
  
  std::vector<int> scale;
  if ( param.has("s") ) 
    {
      scale = param.intvector("s");
      if ( scale.size() != 3 ) Helper::halt( "mse s=lwr,upr,inc" );
    }
  else
    {
      scale.push_back(1);
      scale.push_back(10);
      scale.push_back(2);
    }

  //
  // Output
  //
    
  bool verbose = param.has( "verbose" ) ;

  
  //
  // Attach signal(s)
  //

  std::string signal_label = param.requires( "sig" );  
  
  signal_list_t signals = edf.header.signal_list( signal_label );  

  const int ns = signals.size();

  
  
  //
  // For each signal  
  //
  
  for (int s=0;s<ns;s++)
    {

      //
      // only consider data tracks
      //
      
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;


      logger << " estimating MSE for " << signals.label(s) << "\n";
      
      //
      // output stratifier
      //
      
      writer.level( signals.label(s) , globals::signal_strat );	  
      
      //
      // to track overall mean over epochs: scale -> vector of per-epoch MSEs
      //
      
      std::map<int,std::vector<double> > all_mses;


      //
      // Point to first epoch (assume 30 seconds, but could be different)
      //
      
      int ne = edf.timeline.first_epoch();
      
      if ( ne == 0 ) return;
      

      //
      // for each each epoch 
      //

      while ( 1 ) 
	{
	  
	  //
	  // Get next epoch
	  //
	  
	  int epoch = edf.timeline.next_epoch();
	  
	  if ( epoch == -1 ) break;
	  
	  interval_t interval = edf.timeline.epoch( epoch );
	  	  
	  //
	  // get data
	  //

	  slice_t slice( edf , signals(s) , interval );
	  
	  std::vector<double> * d = slice.nonconst_pdata();

	  
	  //
	  // Mean-centre 30-second window, calculate RMS
	  //

	  mse_t mse( scale[0] , scale[1] , scale[2] , m , r );
	  
	  std::map<int,double> mses = mse.calc( *d );
	  
	  //
	  // track
	  //
	  

	  if ( verbose )
	    writer.epoch( edf.timeline.display_epoch( epoch ) );
	  
	  std::map<int,double>::const_iterator ii = mses.begin();
	  while ( ii != mses.end() )
	    {
	      const int & scale = ii->first ; 
	      
	      all_mses[  scale ].push_back( ii->second ); 
	      
	      // verbose output?

	      if ( verbose )
		{
		  writer.level( scale , "SCALE" );
		  writer.value( "MSE" , ii->second );      	  
		}
	      
	      // next scale
	      ++ii;
	    }
	  
	  if ( verbose )
	    writer.unlevel( "SCALE" );
	      
	} // next epoch
  
      if ( verbose ) 
	writer.unepoch();
      
      
      // overall means
      
      std::map<int,std::vector<double> >::const_iterator ii = all_mses.begin();
      while ( ii != all_mses.end() )
	{
	  double mse = 0;
	  const int & scale = ii->first;
	  const std::vector<double> & x = ii->second;
	  for (int i=0;i<x.size();i++) mse += x[i];
	  mse /= (double)x.size();
	  
	  writer.level( scale , "SCALE" );
	  writer.value( "MSE" , mse );
	  
	  ++ii;
	}
      writer.unlevel( "SCALE" );
      
    } // next signal
  
  writer.unlevel( globals::signal_strat );
  
}




void lzw_per_epoch( edf_t & edf , param_t & param )
{
  
  //
  // Calculate LZW per epoch and average
  //

  //
  // LZW parameters
  //

  int nbins   = param.has( "nbins" ) ? param.requires_int( "nbins" ) : 20 ;
  int nsmooth = param.has( "nsmooth" ) ? param.requires_int( "nsmooth" ) : 1 ;

  bool epoched = param.has( "epoch" );

     
  //
  // Attach signal(s)
  //

  std::string signal_label = param.requires( "sig" );  
  
  signal_list_t signals = edf.header.signal_list( signal_label );  

  const int ns = signals.size();


  //
  // Point to first epoch (assume 30 seconds, but could be different)
  //
  
  int ne = edf.timeline.first_epoch();

  if ( ne == 0 ) return;
  
  
  //
  // For each signal  
  //
  
  for (int s=0;s<ns;s++)
    {
      
      //
      // only consider data tracks
      //
      
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;


      //
      // output stratifier
      //
      
      writer.level( signals.label(s) , globals::signal_strat );	  
      

      //
      // per-epoch, or whole-signal calculation?
      //
      
      if ( ! epoched ) 
	{

	  // get all data
	  interval_t interval = edf.timeline.wholetrace();
	  
	  slice_t slice( edf , s , interval );
	  const std::vector<double> * d = slice.pdata();
	  const int sz = d->size();

	  // designed for per-epoch data, but just use first 
	  // slot for entire signal
	  std::vector<std::vector<double> > track_lzw;
	  track_lzw.push_back( *d );
	  	  
	  // coarse-grain signal
	  coarse_t c( track_lzw , nbins , nsmooth );
	  
	  // compress	  
	  lzw_t lzw( c );
      
	  // index	  
	  double index = lzw.size(0) / (double)track_lzw[0].size();

	  // output
	  writer.value( "LZW" , index );
	  	  
	}
     

      //
      // Epoch level analyses
      //
     
      if ( epoched )
	{
	  
	  
	  int ne = edf.timeline.first_epoch();
	  
	  if ( ne == 0 ) return;

	  // track all signals here

	  std::vector<std::vector<double> > track_lzw;
	  std::vector<int> track_e;

	  //
	  // for each each epoch 
	  //
	  
	  while ( 1 ) 
	    {
	      
	      // Get next epoch
	  
	      int epoch = edf.timeline.next_epoch();
	  
	      if ( epoch == -1 ) break;
	  
	      interval_t interval = edf.timeline.epoch( epoch );
	      
	      // get data
	      
	      slice_t slice( edf , signals(s) , interval );
	      
	      std::vector<double> * d = slice.nonconst_pdata();
	      
	      // lzw_t class is designed for per-epoch data to be taken 
	      // all in one structure	      
	      track_lzw.push_back( *d );
	      track_e.push_back( epoch );
	      
	    } // next epoch	      
	      	      
	  // coarse-grain signal
	  coarse_t c( track_lzw , nbins , nsmooth );
	  
	  // compress	  
	  lzw_t lzw( c );
      
	  // index	  
	  for (int e=0; e<track_e.size(); e++)
	    {

	      double index = lzw.size(e) / (double)track_lzw[e].size();
	      
	      // output
	      writer.epoch( edf.timeline.display_epoch( track_e[e] ) );
	      writer.value( "LZW" , index );
	      writer.unepoch();
			     
	    }
	  
	}
      
      writer.unlevel( globals::signal_strat );
    } // next signal
  
}



void    spike_signal( edf_t & edf , int s1 , int s2 , double wgt , const std::string & ns )
{

  if ( s1 == s2 ) return;
  if ( edf.header.is_annotation_channel(s1) ) Helper::halt( "annotation channel specified for SPIKE" ) ;
  if ( edf.header.is_annotation_channel(s2) ) Helper::halt( "annotation channel specified for SPIKE" ) ;

  bool append_new_channel = ns != "";
  
  interval_t interval = edf.timeline.wholetrace();
  
  // Right now, need a similar sampling rate
  
  const int Fs1 = edf.header.sampling_freq( s1 );
  const int Fs2 = edf.header.sampling_freq( s2 );
  
  const std::string label1 = edf.header.label[s1];
  const std::string label2 = edf.header.label[s2];

  if ( Fs1 != Fs2 ) 
  {
    logger << "Note: resampling " << label2 << " to " << Fs1 << " to match " << label1 << "\n";
    dsptools::resample_channel( edf, s2 , Fs1 );
  }

  slice_t slice1( edf , s1 , interval );
  const std::vector<double> * d1 = slice1.pdata();
  const int sz1 = d1->size();
 
  slice_t slice2( edf , s2 , interval );
  const std::vector<double> * d2 = slice2.pdata();
  const int sz2 = d2->size();
  
  if ( sz1 != sz2 ) Helper::halt( "problem in SPIKE, unequal channel lengths" );
      
  // apply SPIKE
  std::vector<double> spiked( sz1 , 0);
  for (int i=0;i<sz1;i++)
    spiked[i] = (*d1)[i] + wgt * (*d2)[i];
  
  
  // either UPDATE
  if ( append_new_channel )
    {
      const std::string label = edf.header.label[s1] + "-spike-" + edf.header.label[s1] + "-wgt-" + Helper::dbl2str( wgt );      
      edf.add_signal( label , Fs1 , spiked ); 
    }
  else // ... else UPDATE
    edf.update_signal( s1 , &spiked );

}
