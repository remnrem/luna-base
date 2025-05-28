
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

  resp_p1_lwr = param.has( "resp-p1-lwr" ) ? param.requires_dbl( "resp-p1-lwr" ) : 0.1 ;
  resp_p1_upr = param.has( "resp-p1-upr" ) ? param.requires_dbl( "resp-p1-upr" ) : 1 ;
  resp_p2_lwr = param.has( "resp-p2-lwr" ) ? param.requires_dbl( "resp-p2-lwr" ) : 1 ;
  resp_p2_upr = param.has( "resp-p2-upr" ) ? param.requires_dbl( "resp-p2-upr" ) : 10 ;
    
  resp_min_sr = param.has( "resp-min-sr" ) ? param.requires_int( "resp-min-sr" ) : 32 ;
  resp_epsilon = param.has( "resp-epsilon" ) ? param.requires_dbl( "resp-epsilon" ) : 1e-8 ; 
  
  //
  // Generic options
  //
  
  by_epoch = param.has( "epoch" );

  //
  // Annotations & channels
  //

  resp_add_annot = param.has( "resp-add-annot" );
  if ( resp_add_annot ) resp_annot_label = param.value( "resp-add-annot" );

  resp_add_channel = param.has( "resp-add-channel" );
  if ( resp_add_channel ) resp_channel_label = param.value( "resp-add-channel" );
  


  //
  // Perform checks
  //
  
  
  //
  // Respiratory signals
  //
  
  do_resp( resp_signals );
  
  
  //
  // EEG 
  //

  
  
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
	 << "     minimum SR (Hz)   = " << resp_min_sr << "\n";
  
  // force 120 epochs, sliding in 10 seconds
  edf.timeline.set_epoch( resp_window_dur , resp_window_inc );
  
  
  // iterate over all signals       

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

	  //std::cout << snr1 << " " << p1 << " " << valid1 << "\n";
	  // track for this epoch
	  snr.push_back( snr1 );
	  P1.push_back( p1 );
	  valid.push_back( valid1 );
	  smps.push_back( smps1 );

	} // next epoch

            
      //
      // Mean SNR (over valid regions)
      //
      
      const int total_epochs = snr.size();
      int valid_epochs = 0;
      double snr1 = 0;
      for (int e=0; e<total_epochs; e++)
	{
	  if ( valid[e] )
	    {
	      ++valid_epochs;
	      snr1 += snr[e];
	    }
	}

      if ( valid_epochs > 1 )
	snr1 /= (double)valid_epochs;
      
      //
      // Derive % of signal that is noise
      //
      
      // minumum SNR is defined resp_th (by default, 10x)

      // calculate:  P1scale = 10^mean(log10(P1(StoN_x>minStoN)))
      //   i.e. mean of P1 for SNR-positive regions (using log-transform for robustness)
      
      std::vector<double> mm;
      for (int e=0; e<total_epochs; e++)
	{
	  if ( valid[e] && snr[e] > resp_th )
	    mm.push_back( log10( P1[e] ) );
	}
      
      const int nn = mm.size();
      const double P1scale = nn ? pow( 10 , MiscMath::mean( mm ) ) : 0 ; 
      
      // criteria = P1/P1scale.*StoN_x
      //   --> put back in P1
      //   --> 
      for (int e=0; e<total_epochs; e++)
	if ( valid[e] )
	  {
	    P1[e] = P1[e] / P1scale * snr[e] ;
	  }
	
      // EstFNoise2 = mean(criteria<(0.5/4/))
      // Estimated Noise percentage= (100*EstFNoise2)      
      //  ???

      double prop_noise = 0;
      for (int e=0; e<total_epochs; e++)
	if (  valid[e] )
	  prop_noise += P1[e] < 0.5/4 ; // ??
      
      prop_noise /= (double)valid_epochs;

      
      //
      // Main outputs
      //

      writer.value( "SNR" , snr1 );
      
      writer.value( "N_VALID" , valid_epochs );

      writer.value( "P_VALID" , valid_epochs /(double)total_epochs );
      
      writer.value( "P_NOISE" , prop_noise );

      
      //
      // Reconstruct the signal and/or add annotations?
      //
      
      if ( resp_add_channel || resp_add_annot )
	{

	  // pull original signal for entire trace jsut to get N and time-points
	  slice_t slice( edf , signals(s) , edf.timeline.wholetrace());	  
	  const std::vector<uint64_t> * tp = slice.ptimepoints();	  
	  const int n = tp->size();	
	  
	  // create the averaged signal
	  std::vector<double> x1(n,0);
	  std::vector<int> n1(n,0);
	  
	  for (int e=0; e<total_epochs; e++)
	    {
	      const int start = smps[e].first;
	      const int stop  = smps[e].second;	      
	      for (int p=start; p<=stop; p++)
		{
		  x1[p] += snr[e];
		  ++n1[p];
		}
	    }

	  // get mean over windows
	  // for invalid regions, this should be set to 1
	  for (int p=0; p<n; p++)
	    if (n1[p]>1) x1[p] /= (double)n1[p];

	  //
	  // add as a channel?
	  //

	  if ( resp_add_channel )
	    {
	      const std::string lab = signals.label(s) + "_" + resp_channel_label ;
	      edf.add_signal( lab , Fs , x1 );
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

	      // TEMP
	      resp_th = 500;
	      for (int p=0; p<n; p++)
		{
		  //std::cout << " in " << in << " " << x1[p] << "\n";
		  
		  // end of a segment
		  if ( in && ( x1[p] >= resp_th || p == n-1 ) )
		    {
		      a->add( "." ,
			      interval_t( start , (*tp)[p] ) ,
			      signals.label(s) );
		      ++ac;
		      in = false;
		    }
		  else if ( (!in) && x1[p] < resp_th )
		    {
		      start = (*tp)[p];
		      in = true;
		    }
		}
	      
	      logger << "  added " << ac << " " << resp_annot_label
		     << " annotations, marking likely artifact for "
		     << signals.label(s) << "\n";
	      
	    }
	}
      
    } // next signal
  
  writer.unlevel( globals::signal_strat );
  
}


