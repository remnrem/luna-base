
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

#include "dsp/qc.h"

#include "dsp/fir.h"
#include "stats/eigen_ops.h"
#include "edf/edf.h"
#include "edf/slice.h"

#include "db/db.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "fftw/fftwrap.h"

extern logger_t logger;

extern writer_t writer;

dsptools::qc_t::qc_t( edf_t & edf1 , param_t & param ) : edf(edf1) 
{
  
  //
  // break-down by signal types
  //

  //
  // get signals
  //

  const bool no_annotations = true;
  
  const std::string resp_signal_labels = param.has( "resp" ) ? param.value( "resp" ) : "" ;
  const std::string spo2_signal_labels = param.has( "SpO2" ) ? param.value( "SpO2" ) : "" ;
  const std::string eeg_signal_labels = param.has( "eeg" ) ? param.value( "eeg" ) : "" ; 
  
  signal_list_t resp_signals = edf.header.signal_list( resp_signal_labels , no_annotations );
  signal_list_t spo2_signals = edf.header.signal_list( spo2_signal_labels , no_annotations );
  signal_list_t eeg_signals  = edf.header.signal_list( eeg_signal_labels , no_annotations );

  
  //
  // parameters
  //

  // SNR threshold for respiratoruy signals

  resp_th = param.has( "resp-snr-th" ) ? param.requires_dbl( "resp-snr-th" ) : 10 ; 
  resp_window_dur = param.has( "resp-win" ) ? param.requires_dbl( "resp-win" ) : 120 ;
  resp_window_inc = param.has( "resp-inc" ) ? param.requires_dbl( "resp-inc" ) : 10 ;

  // allow 10% of night to be 'bad' (noise-level 2 or worse) 
  resp_prop_th = param.has( "resp-th" ) ? param.requires_dbl( "resp-th" ) : 0.1;  
  
  resp_p1_lwr = param.has( "resp-p1-lwr" ) ? param.requires_dbl( "resp-p1-lwr" ) : 0.1 ;
  resp_p1_upr = param.has( "resp-p1-upr" ) ? param.requires_dbl( "resp-p1-upr" ) : 1 ;
  resp_p2_lwr = param.has( "resp-p2-lwr" ) ? param.requires_dbl( "resp-p2-lwr" ) : 1 ;
  resp_p2_upr = param.has( "resp-p2-upr" ) ? param.requires_dbl( "resp-p2-upr" ) : 10 ;
    
  resp_min_sr = param.has( "resp-min-sr" ) ? param.requires_int( "resp-min-sr" ) : 32 ;
  resp_epsilon = param.has( "resp-epsilon" ) ? param.requires_dbl( "resp-epsilon" ) : 1e-8 ; 

  //
  // eeg
  //
  
  eeg_window_dur = param.has( "eeg-win" ) ? param.requires_dbl( "eeg-win" ) : 30 ;
  eeg_window_inc = param.has( "eeg-inc" ) ? param.requires_dbl( "eeg-inc" ) : 30 ;
  
  eeg_eps = 1e-2;

  eeg_min_sr = 100; // 100 Hz min. sample rate
  
  eeg_min_amp_th = 5;   // uV
  eeg_max_amp_th = 500; // uV

  eeg_spectral_peakedness_th = 0;
  eeg_spectral_skewness_th = 0;

  eeg_h1_min = 0.1; 
  eeg_h1_max = 100;
  
  eeg_h2_min = 0.1;
  eeg_h2_max = 1.0;
  
  eeg_h3_min = 0.1;
  eeg_h3_max = 2.0;
  
  //
  // Generic options
  //
  
  // verbose epoch-level output
  
  by_epoch = param.has( "epoch" );
  
  
  
  //
  // Annotations & channels
  //

  resp_add_annot = param.has( "resp-add-annot" );
  if ( resp_add_annot ) resp_annot_label = param.value( "resp-add-annot" );

  resp_add_channel = param.has( "resp-add-channel" );
  if ( resp_add_channel ) resp_channel_label = param.value( "resp-add-channel" );

  
  //
  // Respiratory signals
  //
  
  do_resp( resp_signals );
  
  
  //
  // EEG 
  //

  do_eeg( eeg_signals );
  
  
}





void dsptools::qc_t::do_resp( signal_list_t & signals )
{
  
  const int ns = signals.size();

  if ( ns == 0 ) return;

  logger << "  checking " << ns << " respiratory channels:\n"
	 << "     window size (sec) = " << resp_window_dur << "\n"
	 << "     window step (sec) = " << resp_window_inc << "\n"
	 << "     signal range (Hz) = " << resp_p1_lwr << " - " << resp_p1_upr << "\n"
    	 << "     noise range (Hz)  = " << resp_p2_lwr << " - " << resp_p2_upr << "\n"
	 << "     SNR threshold     = " << resp_th << "\n"
	 << "     max noise prop.   = " << resp_prop_th << "\n"
	 << "     minimum SR (Hz)   = " << resp_min_sr << "\n";
  
  // force 120 epochs, sliding in 10 seconds
  edf.timeline.set_epoch( resp_window_dur , resp_window_inc );
  

  //
  // iterate over all signals       
  //
  
  for (int s=0; s<ns; s++)
    {

      //
      // just bail for now if SR is too low
      //

      const double Fs = edf.header.sampling_freq( signals(s) );
      
      if ( Fs < resp_min_sr )
	Helper::halt( signals.label(s) + " has a sample rate of "
		      + Helper::dbl2str( Fs ) + ", lower resp-min-sr="
		      + Helper::dbl2str( resp_min_sr ) );      


      //
      // start processing this signal
      //
      
      writer.level( signals.label(s) , globals::signal_strat );

      // for each (overlapping) epoch, track SNR as well as the start/stop sample-points
      // (i.e. use the latter to remake the original signal) 

      std::vector<double> snr;
      std::vector<double> P1;
      std::vector<bool> valid;
      std::vector<std::pair<int,int> > smps;
      
      const int ne = edf.timeline.first_epoch();

      logger << "  considering " << ne << " windows\n"; 
           
      while ( 1 ) 
	{
	  
	  int epoch = edf.timeline.next_epoch();

	  if ( epoch == -1 ) break;
	  
	  interval_t interval = edf.timeline.epoch( epoch );

	  // note: final true returns sample points
	  slice_t slice( edf , signals(s) , interval, 1, false , true );
	  
	  std::vector<double> * d = slice.nonconst_pdata();

	  const std::vector<int> * sp = slice.psmps();

	  const int n = d->size();
	  
	  FFT fftseg( n , n , Fs , FFT_FORWARD , WINDOW_NONE );
	  
	  fftseg.apply( &((*d)[0]) , n );
	  
	  double p1 = 0 , p2 = 0;

	  for (int f=0; f<fftseg.cutoff; f++)
	    {
	      if ( fftseg.frq[f] >= resp_p1_lwr && fftseg.frq[f] < resp_p1_upr )
		{ p1 += fftseg.X[f];  } 
	      else if ( fftseg.frq[f] >= resp_p2_lwr && fftseg.frq[f] < resp_p2_upr )
		{ p2 += fftseg.X[f];  } 
	      else if ( fftseg.frq[f] > resp_p2_upr ) break;
	    }
	  	  
	  // catch degenerate cases	  
	  const bool valid1 = ! ( fabs( p1 ) < resp_epsilon || fabs( p2 ) < resp_epsilon  );
	  
	  // get SNR (?how to handle degenerate cases)
	  const double snr1 = p1 / p2 ; 
	  
	  // sample span of this window (for reconstructing an original signal)
	  std::pair<int,int> smps1( (*sp)[0] , (*sp)[n-1] );

	  // track for this epoch
	  snr.push_back( snr1 );
	  P1.push_back( p1 );
	  valid.push_back( valid1 );
	  smps.push_back( smps1 );

	} // next epoch
      

      //
      // 90th-percentile SNR must be > 10x to consider the signal
      //
      
      const double snr90 = MiscMath::percentile( snr , 0.9 );
      
      bool valid_signal = snr90 > resp_th ;  // defaults to 10x
      
      
      //
      // Get P1 scaling factor (based on log-mean of over valid regions)
      //

      // calculate:  P1scale = 10^mean(log10(P1(StoN_x>minStoN)))
      //   i.e. mean of P1 for SNR-positive regions (using log-transform for robustness)

      int n1 = 0;
      double sum = 0;
      
      const int total_epochs = snr.size();
      int valid_epochs = 0;
      double snr1 = 0;

      for (int e=0; e<total_epochs; e++)
	{
	  if ( valid[e] )
	    {
	      ++valid_epochs;
	      if ( snr[e] > resp_th )
		{
		  sum += log10( P1[e] );
		  ++n1;
		}
	    }
	}
      
      const double P1scale = n1 > 0 ? pow( 10.0 , sum / (double)n1 ) : 1 ; 

      //
      // Make criteria (0 for invalid regions, will be flagged below) 
      //
      
      std::vector<double> criteria( total_epochs );

      for (int e=0; e<total_epochs; e++)
	criteria[e] = valid[e] ? ( P1[e] / P1scale ) * snr[e] : 0 ;
      
      
      //
      // Verbose outputs
      //

      if ( by_epoch )
	{
	  for (int e=0; e<total_epochs; e++)
	    {
	      writer.level( e+1 , "WIN" ) ;
	      writer.value( "P1" , P1[e] );
	      writer.value( "SNR" , snr[e] );
	      writer.value( "CRIT" , criteria[e] );
	    }
	  writer.unlevel( "WIN" );
	}
            
      //
      // Derive % of signal that is noise
      //

      // noisewav below is created as a continuous signal and is
      // assigned values 1, 2 or 3. 

      // 1 : criteria < 1/4    = 0.25
      // 2 : criteria < 0.5/4  = 0.125 
      // 3 : criteria < 0.24/4 = 0.0625 

      // each sample spanned by multiple windows: assign the worst (3 worse than 2 worse than 1)

      
      // pull original signal for entire trace jsut to get N and time-points
      slice_t slice( edf , signals(s) , edf.timeline.wholetrace());
      const std::vector<uint64_t> * tp = slice.ptimepoints();
      const int n = tp->size();

      // new noise waveform (initialize at 0) 
      std::vector<double> noisewav( n , 0 );

      // iterate over each window:
      std::vector<double> fac = { 1 / 4.0 , 0.5 / 4.0 , 0.25 / 4.0 };

      for (int f=0; f<3; f++) // <1, <2, <3	
	{
	  // set for this window
	  for (int e=0; e<total_epochs; e++)
	    {
	      if ( criteria[e] < fac[f] )
		{
		  const int start = smps[e].first;
		  const int stop  = smps[e].second;
		  for (int p=start; p<=stop; p++)
		    noisewav[p] = f+1;
		}
	    }
	}
      

      //
      // Estimate fraction of recording that is noise (nb: allow for floating point
      // issues, though this should be fine) 
      //

      double f1 = 0 , f2 = 0 , f3 = 0;
      for (int p=0; p<n; p++)
	{
	  if ( noisewav[p] >= .99 )
	    {
	      ++f1;
	      if ( noisewav[p] >= 1.99 )
		{
		  ++f2;
		  if ( noisewav[p] >= 2.99 )
		    {
		      ++f3;
		    }
		}
	    }
	}

      // as a proportion of all samples
      f1 /= (double)n;
      f2 /= (double)n;
      f3 /= (double)n;
      
      // take 'f2' as the noise estimate (by default, cannot have more
      // than 0.1 of record as 'bad' based on the F2 noise level
      
      const bool bad_signal = f2 > resp_prop_th; 

      
      //
      // Main outputs
      //

      writer.value( "BAD" , (int)(bad_signal) );
      writer.value( "N_VALID_WIN" , valid_epochs );
      writer.value( "P_VALID_WIN" , valid_epochs /(double)total_epochs );
      writer.value( "P_NOISE1" , f1 );
      writer.value( "P_NOISE2" , f2 );
      writer.value( "P_NOISE3" , f3 );

      
      //
      // add as a channel?
      //
      
      if ( resp_add_channel )
	{
	  const std::string lab = signals.label(s) + "_" + resp_channel_label ;
	  edf.add_signal( lab , Fs , noisewav );
	  logger << "  adding new QC signal " << lab << ", " << Fs << " Hz\n"; 
	}

      //
      // add as annotations?
      //
      
      if ( resp_add_annot )
	{
	  // new annotation
	  annot_t * a = edf.annotations->add( resp_annot_label );
	  
	  // get collapsed regions
	  std::set<interval_t> ints;
	  
	  bool in = false;
	  uint64_t start = 0LLU;
	  int ac = 0;

	  // Flag if noise 2 or more
	  const double fmin = 1.9;
	  
	  for (int p=0; p<n; p++)
	    {
	      
	      // end of a segment
	      if ( in && ( noisewav[p] >= fmin || p == n-1 ) )
		  {
		    a->add( "." ,
			    interval_t( start , (*tp)[p] ) ,
			    signals.label(s) );
		    ++ac;
		    in = false;
		  }
	      else if ( (!in) && noisewav[p] < fmin )
		{
		  start = (*tp)[p];
		  in = true;
		}
	    }
	      
	  logger << "  added " << ac << " " << resp_annot_label
		 << " annotations, marking likely artifact for "
		 << signals.label(s) << "\n";
	      
	}
      
    } // next signal
  
  writer.unlevel( globals::signal_strat );
  
}



void dsptools::qc_t::do_eeg( signal_list_t & signals )
{
  
  const int ns = signals.size();
  
  if ( ns == 0 ) return;
  
  logger << "  checking " << ns << " EEG channels:\n"
	 << "     window size (sec) = " << eeg_window_dur << "\n"
	 << "     window step (sec) = " << eeg_window_inc << "\n"
	 << "     minimum SR (Hz)   = " << eeg_min_sr << "\n";


  // epoch definition: by default 30+30 but can be changed
  edf.timeline.set_epoch( eeg_window_dur , eeg_window_inc );

  // Welch PSD parameters
  const double eeg_fft_seg_sec = 4;
  const double eeg_fft_inc_sec = 2;

  // bandpass filtering (of normed epoch)
  const double eeg_lwr_frq = 1;
  const double eeg_upr_frq = 35;
  const std::vector<double> eeg_ripple = { 0.01 , 0.01 };
  const std::vector<double> eeg_tw = { 0.5 , 5 };

  // freq range for peakedness calculation
  const double eeg_peak_minf = 2;
  const double eeg_peak_maxf = 28;
  const double eeg_peak_median_filter_n = 11;
  
  //
  // iterate over all signals       
  //
  
  for (int s=0; s<ns; s++)
    {

      //
      // just bail for now if SR is too low
      //

      const double Fs = edf.header.sampling_freq( signals(s) );

      if ( Fs < eeg_min_sr ) Helper::halt( "SR too low" );
      
      //
      // pull signal
      //
      
      writer.level( signals.label(s) , globals::signal_strat );
      
      // for each (overlapping) epoch, track SNR as well as the start/stop sample-points
      // (i.e. use the latter to remake the original signal) 

      // each window/epoch in SPs
      std::vector<std::pair<int,int> > smps;

      // is this window/epoch good/bad?
      std::vector<bool> valid;
      
      // track epoch-level values for statistical comparison
      std::vector<double> h1, h2, h3;
      
      const int ne = edf.timeline.first_epoch();
      
      logger << "  considering " << ne << " windows\n"; 
           
      while ( 1 ) 
	{
	  
	  int epoch = edf.timeline.next_epoch();

	  if ( epoch == -1 ) break;
	  
	  interval_t interval = edf.timeline.epoch( epoch );
	  
	  // note: final true returns sample points
	  slice_t slice( edf , signals(s) , interval, 1, false , true );
	  
	  std::vector<double> d = *slice.nonconst_pdata();
	  
	  const std::vector<int> * sp = slice.psmps();
	  
	  const int n = d.size();
	  
	  //
	  // max/min values and clipping ( each is proportion of samples ) 
	  //

	  // clipped signal
	  double c = MiscMath::clipped( d ) ;

	  // flat signal
	  double f = MiscMath::flat( d , eeg_eps ) ;

	  // max signal
	  double mx = MiscMath::max( d , eeg_max_amp_th ) ; 

	  // min signal
	  double mn = MiscMath::min( d , eeg_max_amp_th ) ; 
	    
	    
	   	  
	  //
	  // Filter and normalize (robust time-domain Z)
	  //

	  std::vector<double> flt = dsptools::apply_fir( d ,
							 Fs ,
							 fir_t::BAND_PASS ,
							 1 , // Kaiser window
							 eeg_ripple , eeg_tw ,
							 eeg_lwr_frq , eeg_upr_frq ,
							 0 ,  // order/ignored
							 fir_t::HAMMING , // ignored for KW
							 true // use FFT
							 );

	  //
	  // Robust norm 
	  //

	  // urgh... lazy but got to avoid all this pointless
	  // deep-copying, even if not critical for performance right
	  // now...

	  
	  Eigen::VectorXd zflt = Eigen::VectorXd::Map(flt.data(), flt.size());

	  eigen_ops::robust_scale( zflt , 
				   true, // center (median)
				   true, // normalize (IQR)
				   0 ,   // no winsorization
				   false ); // no second-round normalization needed
	  
	  std::vector<double> zd(zflt.data(), zflt.data() + zflt.size());
	  
	  
	  //
	  // Hjorth parameters
	  //

	  double activity = 0 , mobility = 0 , complexity = 0;
	  double zactivity = 0 , zmobility = 0 , zcomplexity = 0;

	  // use variance-based method

	  const bool USE_VARIANCE_METHOD = true;
	  
	  MiscMath::hjorth( &d , &activity , &mobility , &complexity , USE_VARIANCE_METHOD);
	  
          MiscMath::hjorth( &zd , &zactivity , &zmobility , &zcomplexity , USE_VARIANCE_METHOD );

	  //
	  // Time-domain skew and kurtosis 
	  //

	  double mean   = MiscMath::mean( zd ); // should be 0 / 1 approx...
 	  double sd     = MiscMath::sdev( zd , mean );
	  double skew   = MiscMath::skewness( zd , mean , sd );
	  double kurt   = MiscMath::kurtosis( zd , mean ) ; 
	  
	  
	  //
	  // Spectral metrics (PSD)
	  //
	  
	  const int total_points = zd.size();
	  const int segment_points = eeg_fft_seg_sec * Fs;
	  const int noverlap_points  = eeg_fft_inc_sec * Fs;
	  
	  // implied number of segments
	  int noverlap_segments = floor( ( total_points - noverlap_points) 
					 / (double)( segment_points - noverlap_points ) );
	  
	  PWELCH pwelch( zd , 
			 Fs , 
			 eeg_fft_seg_sec , 
			 noverlap_segments , 
			 WINDOW_TUKEY50 );
	  
	  
	  //
	  // power spectra 
	  //

	  // using bin_t, 1 means no binning [ do we need this? ] 
	  bin_t bin( eeg_lwr_frq, eeg_upr_frq , 1 );
	  bin.bin( pwelch.freq , pwelch.psd );	      
	  
	  //
	  // check for ~zero power values
	  //
	  
	  for ( int i = 0 ; i < bin.bfa.size() ; i++ )
	    {
	      if ( bin.bfb[i] > eeg_peak_maxf ) break;
	      if ( bin.bspec[i] <= 0 && bin.bfa[i] >= eeg_peak_minf ) 
		bin.bspec[i] = 1e-4 ; // set to -40dB as a fudge		   	       
	    }
	  

	  //
	  // Spectral peakedness
	  //

	 
	  // extract out
	  std::vector<double> frq;
	  std::vector<double> logged;

	  for ( int i = 0 ; i < bin.bfa.size() ; i++ )
            {
              if ( bin.bfb[i] > eeg_peak_maxf ) break;
              if ( bin.bfa[i] >= eeg_peak_minf )
		{
		  frq.push_back( bin.bfb[i] );
                  logged.push_back( 10*log10( bin.bspec[i] ) );
                }
            }
	   
	  // not a wide enough region for median smoothing even?
	  const bool no_peakedness = frq.size() < eeg_peak_median_filter_n * 1.5 ;

	  // SPK, KURT
	  double m1 = 0 , m2 = 0;
	  
	  // helper function
	  if ( ! no_peakedness ) 
	    psd_shape_metrics( frq ,
			       logged , 
			       eeg_peak_median_filter_n , 
			       &m1, &m2 ,
			       NULL, NULL, NULL );
	  

	  
	  //
	  // relative band powers
	  //

	  double pow_delta      = pwelch.psdsum( DELTA ) ;      
	  double pow_theta      = pwelch.psdsum( THETA ) ;      
	  double pow_alpha      = pwelch.psdsum( ALPHA ) ;      
	  double pow_sigma      = pwelch.psdsum( SIGMA ) ;      
	  double pow_beta       = pwelch.psdsum( BETA )  ;      
	  double pow_total      = pow_delta + pow_theta + pow_alpha + pow_sigma + pow_beta;
	    
	  // --> store relative powers
	  double rel_delta = pow_delta / pow_total;
	  double rel_theta = pow_theta / pow_total;
	  double rel_alpha = pow_alpha / pow_total;
	  double rel_sigma = pow_sigma / pow_total;
	  double rel_beta  = pow_beta  / pow_total;

	  // --> absolute, normalized (by Hz) values
	  pow_delta /= globals::band_width( DELTA );
	  pow_theta /= globals::band_width( THETA );
	  pow_alpha /= globals::band_width( ALPHA );
	  pow_sigma /= globals::band_width( SIGMA );
	  pow_beta  /= globals::band_width( BETA );

	  std::cout << edf.id << "\t"
		    << epoch << "\t"
		    << c << "\t"
		    << f << "\t"
		    << mx << "\t"
		    << mn << "\t"
		    << mean << "\t"
		    << sd << "\t"
		    << skew << "\t"
		    << kurt << "\t"
		    << activity << "\t"
		    << mobility << "\t"
		    << complexity << "\t"
		    << zactivity << "\t"
		    << zmobility << "\t"
		    << zcomplexity << "\t"
		    << m1 << "\t"
		    << m2 << "\t"
		    << rel_delta << "\t"
		    << rel_theta << "\t"
		    << rel_alpha << "\t"
		    << rel_sigma << "\t"
		    << rel_beta << "\t"
		    << pow_delta << "\t"
		    << pow_theta << "\t"
		    << pow_alpha << "\t"
		    << pow_sigma << "\t"
		    << pow_beta << "\n";
	  
	  
	  //
	  // track sample points
	  //
	  
	  std::pair<int,int> smps1( (*sp)[0] , (*sp)[n-1] );

	  //
	  // all done for this epoch
	  //


	  //
	  // Heuristic to flag whether valid or not
	  //

	  bool valid1 = true;
	  

	  //
	  // track for this epoch
	  //
	  
	  valid.push_back( valid1 );
	  
	  smps.push_back( smps1 );

	} // next epoch
      
      
    }
  
  writer.unlevel( globals::signal_strat );
  
}
