
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

#include "polarity.h"
#include "param.h"

#include "dsp/fir.h"
#include "fftw/fftwrap.h"
#include "dsp/hilbert.h"

#include "helper/helper.h"
#include "helper/logger.h"
#include "db/db.h"

#include "edf/edf.h"
#include "edf/slice.h"

#include <map>

extern writer_t writer;

extern logger_t logger;

void dsptools::polarity( edf_t & edf , const param_t & param )
{

  
  std::string signal_label = param.requires( "sig" );   
  
  signal_list_t signals = edf.header.signal_list( signal_label );  
  
  const int ns = signals.size();


  //
  // Parameters
  //
  
  // SD threshold for extracting from BPF signal  
  double th = param.has( "th" ) ? param.requires_dbl( "th" ) : 1.0;
  
  // extract whole ZC-to-ZC as long as peak meets threshold (default)
  bool zc2zc = param.has( "not-zc2zc" ) ? false : true ;
  
  // up calcate PSD up to 5 Hz
  double flim = param.has( "flim" ) ? param.requires_dbl( "flim" ) : 5.0 ;

  double f_lwr = param.has( "f-lwr" ) ? param.requires_dbl( "f-lwr" ) : 0.5 ;
  double f_upr = param.has( "f-upr" ) ? param.requires_dbl( "f-upr" ) : 4.0 ;
     
  bool mirror_mode = param.has( "not-mirror" ) ? false : true;
  
  bool double_up = param.has( "double" );

  if ( double_up ) mirror_mode = false;
  
  bool analyse_bpf_signal = param.has( "raw" ) ? false : true;

  // delta-mode?  i.e. downward/upward rather than pos/neg ?
  bool d_mode = param.has( "d-mode" );
  if ( d_mode ) 
    {
      mirror_mode = true;
      analyse_bpf_signal = false;
      th = 0;
      flim = 20; // calculate PSD up to 20 Hz
    }

  bool ht_mode = param.has( "ht" );


  logger << " running polarity checks, th=" << th << " for " << f_lwr << "-" << f_upr << "Hz\n";
  
  //
  // Check each signal (assuming it is EEG, filtered for NREM sleep)
  //
  
  for (int s = 0 ; s < ns ; s++ ) 
    {
      
      //
      // only consider data tracks
      //
      
      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;

      //
      // get signal data and meta-data
      //

      double sr = edf.header.sampling_freq( signals )[s];
      
      interval_t interval = edf.timeline.wholetrace();
      
      slice_t slice( edf , signals(s) , interval );
      
      const std::vector<double> * d = slice.pdata();
      
      const std::vector<uint64_t> * tp = slice.ptimepoints(); 
      
      //
      // perform polarity check
      //

      writer.level( signals.label(s) , globals::signal_strat );

      if ( ht_mode) 
	ht_polarity_check( *d , tp , sr , f_lwr , f_upr );
      else
	polarity_check(  *d , tp , sr , 
			 th, zc2zc , flim, f_lwr , f_upr , 
			 mirror_mode , double_up , analyse_bpf_signal , d_mode );
      

      writer.unlevel( globals::signal_strat );

    }
  
  
}



// Helper to extract regions up to zero-crossings given |signal| > |threshold|

std::vector<bool> dsptools::make_mask( const std::vector<double> & x , double th  )
{

  const int ns = x.size();

  std::vector<bool> masked( ns , true );

  for (int i=0;i<ns;i++) 
    {

      bool downpeak = x[i] < -th;
      bool uppeak   = x[i] >  th;
      bool peak = uppeak | downpeak; 
      
      if ( peak ) // track back/forward to next ZC
	{
	  int j = i;
	  while ( 1 ) 
	    {
	      if ( j == 0 ) break;
	      --j;
	      if ( ( downpeak && x[j] > 0 ) || ( uppeak && x[j] < 0 ) ) { ++j; break; }
	    }
	  for (;j!=i;j++) masked[j] = false; // unmask	      
	  // now forward
	  j=i;
	  while ( 1 ) 
	    {
	      ++j;
	      if ( j == ns ) { --j; break; } 
	      if ( ( downpeak && x[j] > 0 ) || ( uppeak && x[j] < 0 ) ) { --j; break; } 
	    }
	  for (;j!=i;j--) masked[j] = false; // unmask	      
	}
    }    
  return masked;
}


void dsptools::ht_polarity_check( const std::vector<double> & x , 
				  const std::vector<uint64_t> * tp , 
				  int fs , double f_lwr , double f_upr )
{

  double th = 2;
  
  double ripple = 0.01;
  double tw = 0.5;

  hilbert_t hilbert( x , fs , f_lwr , f_upr , ripple, tw );
  
  const std::vector<double> * phase = hilbert.phase();
  const std::vector<double> * mag = hilbert.magnitude();
  std::vector<double> frq = hilbert.instantaneous_frequency( fs );

  if ( phase->size() != frq.size()+1 ) Helper::halt( "internal problem in ht_polarity_check()" );
  if ( phase->size() != x.size() ) Helper::halt( "internal problem in ht_polarity_check()" );

  // for now just bin by integer value of phase
  std::map<int,int> cnt;
  std::map<int,double> val;
  std::map<int,double> eeg;

  //
  // only consider large peaks in BPF signal
  //
  
  const std::vector<double> signal = MiscMath::Z( *hilbert.signal() );

  std::vector<bool> masked = make_mask( signal , th );

  //  std::vector<double> Zmag = MiscMath::Z( *mag );

  const int n = x.size();
  const double lim = f_upr * 2 ;
  for (int i = 0 ; i < n-1 ; i++ ) 
    {
      int pi = round( (*phase)[i] );
      if ( frq[i] > 0 & frq[i] < lim & ! masked[i] )
	{
	  cnt[pi]++;
	  val[pi]+=frq[i];
	  eeg[pi]+=x[i];
	}
      std::cout << "zzz\t" << x[i] << "\t" << masked[i] << "\t" << (*phase)[i] << "\t" << frq[i] << "\n";
    }


  
  std::cout << "xxx";

  std::map<int,int>::const_iterator ii = cnt.begin();
  while ( ii != cnt.end() ) 
    {
      std::cout << "\t" << val[ ii->first ] / (double) ii->second ;
      ++ii;
    }

  ii = cnt.begin();
  while ( ii != cnt.end() ) 
    {
      std::cout << "\t" << eeg[ ii->first ] / (double) ii->second ;
      ++ii;
    }

  ii = cnt.begin();
  while ( ii != cnt.end() ) 
    {
      std::cout << "\t" << ii->second ;
      ++ii;
    }

  std::cout << "\n";
  
  
}



void dsptools::polarity_check( const std::vector<double> & x0 , const std::vector<uint64_t> * tp , int fs , 
			       double th ,  // threshold for extracted filtered segments
			       bool zc2zc , // extract around peaks up to ZC's
			       double flim , // calculate PSD up to flim Hz only
			       double f_lwr , // lower BPF transition frequency
			       double f_upr , // upper 
			       bool mirror_mode ,   // mirror odd numbered up/down segments, i.e. make 'wave-like' signal
			       bool double_up, // instead of mirror-mode alternates, double/reflect each interval
			       bool analyse_bpf_signal , // analysis the BPF signal, not raw data
			       bool d_mode  // segment 'upward' and 'downward' rather than pos and neg			       
			       )
{

  if ( double_up && d_mode ) Helper::halt( "not implemented yet, d-mode and double" );

  // for PSD reporting ;; lower and increment;  flim is upper
  double flow = 0.25;
  double finc = 0.25;


  //
  // Trackers
  //

  std::vector<double> delta_activity, delta_mobility, delta_complexity;

  // for skewness tests
  std::vector<double> sigmean , sigmedian , sigdiff;

  // in d-mode, track mean up and down components too
  std::vector<double> avg_up_activity, avg_up_mobility, avg_up_complexity;
  std::vector<double> avg_down_activity, avg_down_mobility, avg_down_complexity;
  
  // assume PSD freq 1..20 Hz
  std::map<int,std::vector<double> > delta_psd;
  std::map<int,std::vector<double> > delta_relpsd;

  std::map<int,std::vector<double> > avg_up_psd;
  std::map<int,std::vector<double> > avg_down_psd;
  std::map<int,std::vector<double> > avg_up_relpsd;
  std::map<int,std::vector<double> > avg_down_relpsd;
  
  double avg_up_time = 0 , avg_down_time = 0;
  
  std::vector<double> x = x0;
  
  const int ns = x.size();
  const int epoch = 30 * (double)fs; 
  const int ne = ns / epoch;
  
  //
  // Band-pass filter and normalize signal
  //

  std::vector<double> ripple( 1, 0.01 ); // 40dB attenuation
  std::vector<double> tw( 1, 0.5 );
  
  std::vector<double> f = MiscMath::Z( dsptools::apply_fir( x , fs , fir_t::BAND_PASS ,
							    1 , // Kaiser
							    ripple , tw , // ripple/ TW
							    f_lwr , f_upr ) );

  if ( f.size() != ns ) Helper::halt( "problem in dsp::polarity()" );

  //
  // Adjust filtered and original signal by 'th' ammount, if only
  // taking, but only if not taking the entire signal
  //
  
  if ( ! zc2zc )
    {
      for (int i=0; i < ns; i++) 
	{      
	  if      ( f[i] > -th && f[i] < th ) x[i] = 0 ; 
	  else if ( f[i] >= th  ) { x[i] -= th; f[i] -= th; }
	  else if ( f[i] <= -th ) { x[i] += th; f[i] += th; }      
	}
    }

  //
  // Determine masked points, based on filtered signal
  //

  // if in ZC mode, start with all masked 'out'
  
  std::vector<bool> masked( ns , true );
  
  if ( zc2zc ) 
    {
      masked = make_mask( f , th );
    }
  else
    {
      for (int i=0;i<f.size();i++) 
	masked[i] = f[i] > -th && f[i] < th ; 
    }


  //
  // for each epoch
  //

  int c = 0;
  for (int e = 0 ; e < ne ; e++ ) 
    {
      
      std::vector<double> edata;
      std::vector<double> fdata;

      for (int i=0;i<epoch;i++)
	{
	  if ( ! masked[c] ) 
	    {
	      edata.push_back( analyse_bpf_signal ? f[c] : x[c] );
	      fdata.push_back( f[c] );
	    }
	  ++c;
	}


      
      std::vector<double> down, up;

      if ( d_mode ) 
	{
	  
	  int nxt = 1;
	  int direction = 0;

	  while ( 1 )
	    { 
	      if ( nxt == edata.size() ) break;
	      if ( edata[nxt] > edata[0] ) { direction = 1 ; break; }
	      if ( edata[nxt] < edata[0] ) { direction = -1 ; break; }
	      ++nxt;
	    }

	  // starts both with 0
	  if ( direction == 1 ) up.push_back( edata[0] );
	  else down.push_back( edata[0] );
	  
	  int up_mirror = 1;
	  int down_mirror = 1;

	  for (int i=1;i<edata.size();i++)
	    {
	      
	      int fdir = fdata[i] - fdata[i-1] > 0 ? + 1 : ( fdata[i] - fdata[i-1] < 0 ? -1 : 0 ) ;
	      
	      //double delta = edata[i] - edata[i-1] ;
	      //double value = delta;
	      double value = edata[i];

	      // change mode?   up to down
	      if ( direction == 1 && fdir == -1 )
		{
		  direction = -1;
		  //		  down_mirror *= -1;
		  down.push_back( value );
		  // if ( down.size() == 0 ) down.push_back( down_mirror * value );
		  // else down.push_back( down[ down.size() - 1 ] - down_mirror * value );
		}
	      else if ( direction == -1 && fdir == 1 )
		{
		  direction = +1;
		  //up_mirror *= -1;
		  // if ( up.size() == 0 ) up.push_back( up_mirror * delta );
		  // else up.push_back( up[ up.size() - 1 ] + up_mirror * value );
		  up.push_back( value );
		}	      
	      else if ( direction == 1 ) // continue in same direction
		{
		  up.push_back( value );
		  //up.push_back( up[ up.size() - 1 ] + up_mirror * value );
		}
	      else if ( direction == -1 ) // is still checked here, as could be flat  epoch (==0) 
		{
		  down.push_back( value );
		  //down.push_back( down[ down.size() - 1 ] - down_mirror * value );
		}

	      // std::cout << edata[i] << "\t" 
	      // 		<< fdata[i] << "\t"
	      // 		<< fdir << "\t"
	      // 		<< direction << "\t";
	      // 	//<< delta << "\t";
	      // if ( direction == 1 ) 
	      // 	std::cout <<  up[ up.size() - 1 ]  << "\t" << "" << "\n";
	      // else if ( direction == -1 ) 
	      // 	std::cout <<  "" << "\t" << down[ down.size() - 1 ] << "\n";
	      // else 
	      // 	std::cout << "."  << "\t" << "." << "\n";

	    }

	}
      else
	{

	  int mirror_up = 1;
	  int mirror_down = 1;

	  // double-up mode?
	  if ( double_up ) 
	    {
	      
	      std::vector<double> buffer;
	      
	      for (int i=0;i<edata.size();i++)
		{
		  bool is_up = edata[i] > 0 ;
		  bool is_down = edata[i] < 0;
		  
		  // state change?
		  if ( is_up && i != 0 && edata[i-1] <= 0 ) 
		    {
		      for (int j=0;j<buffer.size();j++) down.push_back( buffer[j] );
		      for (int j=0;j<buffer.size();j++) down.push_back( -buffer[j] );
		      buffer.clear();
		    }
		  else if ( is_down && i != 0 && edata[i-1] >= 0 ) 
		    {
		      for (int j=0;j<buffer.size();j++) up.push_back( buffer[j] );
		      for (int j=0;j<buffer.size();j++) up.push_back( -buffer[j] );
		      buffer.clear();
		    }
		   
		  buffer.push_back( edata[i] );
		  
		}
	      
	      //	      logger << " edata " << edata.size() << " " << up.size() << " " << down.size() << "\n";
	      
	    }
	  else
	    {
	      for (int i=0;i<edata.size();i++)
		{
		  bool is_up = edata[i] > 0 ;
		  bool is_down = edata[i] < 0;
		  
		  if ( mirror_mode ) 
		    {
		      // state change?
		      if ( is_up && up.size() > 0 && edata[i-1] <= 0 ) mirror_up *= -1;
		      if ( is_down && down.size() > 0 && edata[i-1] >= 0 ) mirror_down *= -1;
		      
		      if ( is_up ) up.push_back( edata[i] * mirror_up ); 
		      if ( is_down ) down.push_back( edata[i] * mirror_down ); 	      
		      
		    }
		  else
		    {
		      if      ( is_up ) up.push_back( edata[i] );
		      else if ( is_down ) down.push_back( edata[i] );
		    }
		}
	    }
	}
    	        

      // for (int u=0;u<up.size();u++) std::cout << "u\t" << up[u] << "\n";
      // for (int d=0;d<down.size();d++) std::cout << "d\t" << down[d] << "\n";
      
      if ( up.size() < 10 || down.size() < 10 ) continue; 
      
      // std::cout << " << EPOCH << \n\n";
      
      // for (int i=0;i<2000;i++) std::cout << "e " << edata[i] << "\n";

      // for (int i=0;i<2000;i++) std::cout << "t " << up[i] << "\n";
      // for (int i=0;i<2000;i++) std::cout << "b " << down[i] << "\n";

      
      //
      // Hjorth param
      //

      double up_activity = 0 , up_mobility = 0 , up_complexity = 0;
      MiscMath::hjorth( &up , &up_activity , &up_mobility , &up_complexity , ! globals::legacy_hjorth );
      
      double down_activity = 0 , down_mobility = 0 , down_complexity = 0;
      MiscMath::hjorth( &down , &down_activity , &down_mobility , &down_complexity , ! globals::legacy_hjorth );
            
      
      //
      // Mean/skewness tests
      //

      double epoch_mean = MiscMath::mean ( edata );
      double epoch_median = MiscMath::median( edata );
      sigmean.push_back( epoch_mean );
      sigmedian.push_back( epoch_median );
      sigdiff.push_back( epoch_mean - epoch_median );
      


      //
      // PSD at 1 Hz bins
      //
      
      double up_time = up.size() / (double)fs;
      double down_time = down.size() / (double)fs;

      avg_up_time += up_time;
      avg_down_time += down_time;


      bool do_fft = up_time > 4 && down_time > 4;
      std::map<int,double> up_psd;
      std::map<int,double> down_psd;
      
      std::map<int,double> up_relpsd;
      std::map<int,double> down_relpsd;

      
      // require atleast 4 seconds 
      if ( do_fft )
	{

	  // be sure to take same length of data from each
	  // as we don't need to process these beyond the FFT
	  // we can directly chop the vectors

	  int minx = up.size() < down.size() ? up.size() : down.size();	  
	  up.resize( minx );
	  down.resize( minx );

	  // for
	  const double segment_sec  = 4;
	  const int segment_points = segment_sec * fs;
	  const double overlap_sec = 2;
	  
	  const int total_points = up.size();
	  const int noverlap_points  = overlap_sec * fs;
           
	  // implied number of segments
	  int noverlap_segments = floor( ( total_points - noverlap_points) 
					 / (double)( segment_points - noverlap_points ) );
	  
	  PWELCH up_pwelch( up , fs , segment_sec , noverlap_segments );
	  PWELCH down_pwelch( down , fs , segment_sec , noverlap_segments );
	  
	  double up_tot_pow = 0 , down_tot_pow = 0;
	
	  // track power bands     
	  for (double i=flow;i<=flim;i+=finc)
	    {
	      int idx = 100 * i ; 
	      up_psd[idx] = up_pwelch.psdsum( i , i+1);
	      down_psd[idx] = down_pwelch.psdsum( i , i+1 );
	      up_tot_pow += up_psd[idx];
	      down_tot_pow += down_psd[idx];
	    }
	  
	  // relative power
	  for (double i=flow;i<=flim;i+=finc)
	    {
	      int idx = 100 * i ; 
	      up_relpsd[idx] = up_psd[idx] / up_tot_pow;
	      down_relpsd[idx] = down_psd[idx] / down_tot_pow;
	    }

	}
      


      //
      // Differences in summary stats
      //

      double d_activity = up_activity - down_activity;
      double d_mobility = up_mobility - down_mobility;
      double d_complexity = up_complexity - down_complexity;
      
      delta_activity.push_back( d_activity );
      delta_mobility.push_back( d_mobility );
      delta_complexity.push_back( d_complexity );

     
      if ( do_fft )
	{
	  for (double i=flow;i<=flim;i+=finc)
	    {
	      int idx = 100 * i ; 
	      double d_psd = log( up_psd[idx] ) - log( down_psd[idx] );
	      delta_psd[idx].push_back( d_psd );
	      
	      double d_relpsd = up_relpsd[idx]  - down_relpsd[idx] ;
	      delta_relpsd[idx].push_back( d_relpsd );
	      
	    }
	}


      if ( d_mode )
	{
	  avg_up_activity.push_back( up_activity );
	  avg_up_mobility.push_back( up_mobility );
	  avg_up_complexity.push_back( up_complexity );

	  avg_down_activity.push_back( down_activity );
	  avg_down_mobility.push_back( down_mobility );
	  avg_down_complexity.push_back( down_complexity );

	  if ( do_fft )
	    {
	      for (double i=flow;i<=flim;i+=finc)
		{
		  int idx = 100 * i ; 
		  avg_up_psd[idx].push_back( up_psd[idx] );
		  avg_down_psd[idx].push_back( down_psd[idx] );

		  avg_up_relpsd[idx].push_back( up_relpsd[idx] );
		  avg_down_relpsd[idx].push_back( down_relpsd[idx] );
		}
	    }
	  
	}
      
    }
  


  const int min_epochs_required = 10;
  const int n_h_epochs = delta_activity.size();

  writer.value("UP_TIME" , avg_up_time / (double)n_h_epochs );
  writer.value("DOWN_TIME" , avg_down_time / (double)n_h_epochs );

  if ( delta_complexity.size() > min_epochs_required )
    {
      
      double mean_sigmean = MiscMath::mean( sigmean );
      double mean_sigmedian = MiscMath::mean( sigmedian );
      double mean_sigdiff = MiscMath::mean( sigdiff );
      double sd_sigdiff = MiscMath::sdev( sigdiff );
      double t_sigdiff = mean_sigdiff / ( sd_sigdiff / sqrt( n_h_epochs ) );
      writer.value( "MN" , mean_sigmean );
      writer.value( "MD" , mean_sigmedian );
      writer.value( "T_DIFF" , t_sigdiff );
      
      double mean_delta_activity = MiscMath::mean( delta_activity );
      double sd_delta_activity = MiscMath::sdev( delta_activity );
      
      double mean_delta_mobility = MiscMath::mean( delta_mobility );
      double sd_delta_mobility = MiscMath::sdev( delta_mobility );
      
      double mean_delta_complexity = MiscMath::mean( delta_complexity );
      double sd_delta_complexity = MiscMath::sdev( delta_complexity );
      
      double t_activity = mean_delta_activity / ( sd_delta_activity / sqrt( n_h_epochs ) );
      double t_mobility = mean_delta_mobility / ( sd_delta_mobility / sqrt( n_h_epochs ) );
      double t_complexity = mean_delta_complexity / ( sd_delta_complexity / sqrt(n_h_epochs ) );
      
      writer.value( "T_H1" , t_activity );
      writer.value( "T_H2" , t_mobility );
      writer.value( "T_H3" , t_complexity );

      if ( d_mode ) 
	{

	  double mean_up_activity = MiscMath::mean( avg_up_activity );
	  double mean_up_mobility = MiscMath::mean( avg_up_mobility );
	  double mean_up_complexity = MiscMath::mean( avg_up_complexity );
	  
	  writer.value( "UP_H1" , mean_up_activity );
	  writer.value( "UP_H2" , mean_up_mobility );
	  writer.value( "UP_H3" , mean_up_complexity );

	  double mean_down_activity = MiscMath::mean( avg_down_activity );
	  double mean_down_mobility = MiscMath::mean( avg_down_mobility );
	  double mean_down_complexity = MiscMath::mean( avg_down_complexity );
	  
	  writer.value( "DOWN_H1" , mean_down_activity );
	  writer.value( "DOWN_H2" , mean_down_mobility );
	  writer.value( "DOWN_H3" , mean_down_complexity );
	 
	}

      writer.value( "N_H" , (int)delta_activity.size() );
      writer.value( "N_FFT" , (int)delta_psd[ 0.5 * 100 ].size() );
      
      
      if ( delta_psd[ 0.5 * 100 ].size() > min_epochs_required ) 
	{

	  for (double i=flow;i<=flim;i+=finc)
	    {
	      
	      int idx = 100 * i ; 

	      writer.level( i , globals::freq_strat );
	      
	      if ( 0 ) 
		{
		  double mean = MiscMath::mean( delta_psd[idx] );
		  double sd = MiscMath::sdev( delta_psd[idx] , mean );
		  double t = mean / ( sd / sqrt( delta_psd[idx].size() ) );	      
		  writer.value( "T_PSD" , t );
		}

	      double mean = MiscMath::mean( delta_relpsd[idx] );
	      double sd = MiscMath::sdev( delta_relpsd[idx] , mean );
	      double t = mean / ( sd / sqrt( delta_relpsd[idx].size() ) );
	      writer.value( "T_RELPSD" , t );
	      
	      if ( d_mode )
		{
	  
		  mean = MiscMath::mean( avg_up_psd[idx] );		  
		  writer.value( "UP_PSD" , mean );

		  mean = MiscMath::mean( avg_down_psd[idx] );
		  writer.value( "DOWN_PSD" , mean );

		  mean = MiscMath::mean( avg_up_relpsd[idx] );		  
		  writer.value( "UP_RELPSD" , mean );

		  mean = MiscMath::mean( avg_down_relpsd[idx] );
		  writer.value( "DOWN_RELPSD" , mean );
		}
	      
	    }
	  writer.unlevel( globals::freq_strat );
	}

    }
      
    
}
