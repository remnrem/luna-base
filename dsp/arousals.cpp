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

#if 0

#include "dsp/arousals.h"
#include "param.h"
#include "defs/defs.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "stats/eigen_ops.h"
#include "stats/kmeans_eigen.h"
#include "miscmath/miscmath.h"
#include "fftw/fftwrap.h"
#include "miscmath/ghmm.h"
#include "dsp/spline.h"

#include "helper/logger.h"
#include "db/db.h"

#include <memory>

extern logger_t logger;
extern writer_t writer;

namespace {

bool arousal_verbose_annots = false;
std::string arousal_verbose_prefix = "arv_";

bool arousal_bool_token( const std::string & s )
{
  return s == "T" || s == "t" || s == "Y" || s == "y" ||
    s == "1" || s == "F" || s == "f" || s == "N" || s == "n" || s == "0";
}

bool arousal_true_token( const std::string & s )
{
  return s == "T" || s == "t" || s == "Y" || s == "y" || s == "1";
}

}

// expected inputs
//  High-pass ~0.3–0.5 Hz
//  Low-pass 30–35 Hz
//  Notch at 50/60 Hz if needed
// EMG:
//   Bandpass: 10–100 Hz

arousals_t::arousals_t( edf_t & edf , param_t & param )
{

  parent = &edf;
  
  // use EPOCH mechanism to add new features
  // requires EPOCH inc is 1 Hz
  //  default epoch size = 4 seconds, 0.5 second shift
  //  i.e. 8x overlap, 2 Hz feature sampling

  const double epoch_win = param.has( "win" ) ? param.requires_dbl( "win" ) : 4.0 ;
  const double epoch_inc = param.has( "inc" ) ? param.requires_dbl( "inc" ) : 0.5 ;
  if ( epoch_win < 1 || epoch_inc <= 0 || epoch_win < epoch_inc )
    Helper::halt( "invalid epoch win/inc values" );
  
  int ne = edf.timeline.set_epoch( epoch_win, epoch_inc ) ;
  
  logger << "  deriving features for " << ne << " "
	 << epoch_win << "s epochs, w/ "
	 << 100*(epoch_inc / epoch_win) << "% overlap\n";

  
  // to save the signal, the EDF record size must be consistent w/ the
  // new SR of the signal i.e. to have same number of samples per EDF
  // record, cannot have SR = 0.5 Hz with 1-second records, etc
  // perhaps for now, enfore that EDF record size is a multiple of 1
  // second, and always output at 1 Hz
  
  uint64_t d = edf.header.record_duration_tp % globals::tp_1sec;
  if ( d != 0 )
    Helper::halt( "must have EDF records that are a multiple of 1 second (use RECORD-SIZE)" );
 

  //
  // Options
  //
  
	  if ( param.has( "winsor" ) || param.has( "no-winsor" ) )
	    Helper::halt( "AROUSALS winsor/no-winsor options are not supported by the current heuristic" );
	  if ( param.has( "annot" ) )
	    Helper::halt( "AROUSALS annot option is not supported; annotation class names are fixed" );
	  if ( param.has( "prefix" ) )
	    Helper::halt( "AROUSALS prefix option is not supported; use add=<prefix> for derived channels" );
	  if ( param.has( "per-channel" ) )
	    Helper::halt( "AROUSALS per-channel option is not supported by the current heuristic" );


  // optional diagnostic annotation tracks for the intermediate masks and
  // candidate intervals used by the heuristic detector
  arousal_verbose_annots = false;
  arousal_verbose_prefix = "arv_";
  const std::string verbose_key = param.has( "verbose-annot" ) ? "verbose-annot" :
    ( param.has( "verbose-annots" ) ? "verbose-annots" : "" );
  if ( verbose_key != "" )
    {
      if ( param.empty( verbose_key ) )
	arousal_verbose_annots = true;
      else
	{
	  const std::string v = param.value( verbose_key );
	  if ( arousal_bool_token( v ) )
	    arousal_verbose_annots = arousal_true_token( v );
	  else
	    {
	      arousal_verbose_annots = true;
	      arousal_verbose_prefix = v;
	    }
	}
    }

  // add channels
  const bool add_chs = param.has( "add" );
  const std::string ch_prefix = ! param.empty( "add" ) ? param.value( "add" ) : "a_" ;

  // NREM/REM arousal heuristic detector thresholds & toggles
  arousal_heuristic_params_t hp;
  hp.epoch_inc = epoch_inc;
  hp.do_nrem = param.yesno( "nrem" , true , true );
  hp.do_rem  = param.yesno( "rem"  , true , true );
  hp.verbose = param.has( "verbose" );
  if ( param.has( "mode" ) )
    {
      const std::string mode = param.value( "mode" );
      if ( mode == "emg" || mode == "EMG" )
	hp.emg_mode = true;
      else if ( mode != "shift" && mode != "SHIFT" )
	Helper::halt( "mode must be shift or emg" );
    }
  if ( param.has( "shift-peak" )       ) hp.shift_peak       = param.requires_dbl( "shift-peak" );
  if ( param.has( "shift-hysteresis" ) ) hp.shift_hysteresis = param.requires_dbl( "shift-hysteresis" );
  if ( param.has( "sigma-veto" )        ) hp.sigma_veto_nrem   = param.requires_dbl( "sigma-veto" );
  if ( param.has( "sigma-veto2" )       ) hp.sigma_veto2_nrem  = param.requires_dbl( "sigma-veto2" );
  if ( param.has( "sigma-veto-rem" )    ) hp.sigma_veto_rem    = param.requires_dbl( "sigma-veto-rem" );
  if ( param.has( "sigma-veto2-rem" )   ) hp.sigma_veto2_rem   = param.requires_dbl( "sigma-veto2-rem" );
  if ( param.has( "th-emg-artifact" ) ) hp.th_emg_artifact = param.requires_dbl( "th-emg-artifact" );
  if ( param.has( "th-h3-artifact" )  ) hp.th_h3_artifact  = param.requires_dbl( "th-h3-artifact" );
  if ( param.has( "th-pwr-artifact" ) ) hp.th_pwr_artifact = param.requires_dbl( "th-pwr-artifact" );
  if ( param.has( "emg-median" )      ) hp.emg_median_sec  = param.requires_dbl( "emg-median" );
  if ( hp.emg_median_sec < 0 )
    Helper::halt( "emg-median must be 0 or greater" );
  if ( param.has( "min-dur" )     ) hp.min_dur       = param.requires_dbl( "min-dur" );
  if ( param.has( "max-dur" )     ) hp.max_dur       = param.requires_dbl( "max-dur" );
  if ( param.has( "long-dur" )    ) hp.long_dur      = param.requires_dbl( "long-dur" );
  if ( param.has( "arousal-dur" ) ) hp.arousal_dur   = param.requires_dbl( "arousal-dur" );
  if ( param.has( "merge-gap" )   ) hp.merge_gap_sec = param.requires_dbl( "merge-gap" );
  if ( param.has( "pre-sleep" )   ) hp.pre_sleep_sec = param.requires_dbl( "pre-sleep" );
  if ( param.has( "emg-rise-th" )      ) hp.emg_rise_th      = param.requires_dbl( "emg-rise-th" );
  if ( param.has( "emg-rise-min-dur" ) ) hp.emg_rise_min_dur = param.requires_dbl( "emg-rise-min-dur" );
  if ( param.has( "emg-rise-buffer" )  ) hp.emg_rise_buffer  = param.requires_dbl( "emg-rise-buffer" );

  const int wake_bridge = param.has( "wake-bridge" ) ? param.requires_int( "wake-bridge" ) : 1;
  if ( wake_bridge < 0 )
    Helper::halt( "wake-bridge must be 0 or greater" );

  //
  // Get signals
  //
  
  const bool has_eeg = param.has( "eeg" );
  const bool has_emg = param.has( "emg" );
  
  const std::string eeg_signal_label = has_eeg ? param.value( "eeg" ) : "" ;
  const std::string emg_signal_label = has_emg ? param.value( "emg" ) : "" ;
    
  const bool NO_ANNOTS = true;   
  signal_list_t eeg_signals = edf.header.signal_list( eeg_signal_label , NO_ANNOTS );
  signal_list_t emg_signals = edf.header.signal_list( emg_signal_label , NO_ANNOTS );  
  
  const int ns_eeg = eeg_signals.size();
  const int ns_emg = emg_signals.size();
  const int ns     = ns_eeg + ns_emg;
  if ( hp.emg_mode && ns_emg == 0 )
    Helper::halt( "mode=emg requires an EMG channel" );
  if ( hp.do_rem && ns_emg == 0 )
    Helper::halt( "AROUSALS REM detection requires an EMG channel; use rem=F to run NREM-only" );

  if ( ns_eeg == 0 && ! hp.emg_mode )
    {
      logger << "  no valid EEG signals detected... leaving AROUSALS\n";
      return;
    }
  
  logger << "  running AROUSALS for " << ns_eeg << " EEG and " << ns_emg << " EMG signals\n";
  
    

  //
  // Obtain sampling freqs (Hz) - these must be similar
  //
  
  std::vector<double> Fs_eeg = edf.header.sampling_freq( eeg_signals );
  std::vector<double> Fs_emg = edf.header.sampling_freq( emg_signals );

  if ( ns_eeg )
    {
      const double Fs0 = Fs_eeg[0];
      for (int i=0; i<ns_eeg; i++)
	{
	  if ( fabs( Fs_eeg[i] - Fs0 ) > 1e-4 )
	    Helper::halt( "all EEG signals must have similar sample rates" );
	}

      const double sr_rounded = std::round( Fs0 );
      if ( fabs( Fs0 - sr_rounded ) > 1e-4 )
	Helper::halt( "AROUSALS requires integer EEG sample rates" );

      sr = (int)sr_rounded;

      if ( sr < 60 )
	Helper::halt( "EEG sample rate too low" );
    }
  
  for (int i=1; i<ns_emg; i++)
    {
      if ( fabs( Fs_emg[i] - Fs_emg[0] ) > 1e-4 )
	Helper::halt( "all EMG signals must have similar sample rates" );      
    }
  
  //
  // Construct initial feature matrix, w/ time-track and NREM/REM annotation
  //
  
  // state (based on *start* of epoch)
  std::vector<int> state; // 0 = NR , 1 = R  , 2 = W
  std::vector<int> seq; // contig w/in each stage/state
  std::vector<double> sec; // time of each feature-window center
  
  // features:
  //    log-pwr / arousal-shift / emg-rms / [rel-sigma] / [h3 ]
  //  only first 3 ftrs used for hmm/detector

  // low-pwr --> artifact: broadband (0.5-35Hz) “movement/artifact” axis (optionally: 30-45 Hz HF-EEG)
  // arousal-shift --> EEG arousal
  // emg     --> EMG bursts
  // sigma (subtyping events)
  // h3    (subtyping)

  // build features
  Eigen::MatrixXd Xeeg, Xemg;
  build_ftr_matrix( edf , eeg_signals, emg_signals , Xeeg, Xemg, &state, &seq, &sec , wake_bridge , epoch_inc ) ;
  
  // process [ epoch x (3+2) ftrs ] 
  Eigen::MatrixXd Xftr = process_ftr_matrix( &Xeeg , &Xemg , state , sec , hp.emg_median_sec , hp.epoch_inc );

  // tidy up
  Xeeg.resize(0,0); Xemg.resize(0,0);
  
  // time-track
  std::vector<std::vector<std::vector<double> > > tt;

  // assemble into sleep-state-specific contigs
  std::vector<std::vector<std::vector<Eigen::VectorXd> > > X = assemble( Xftr, state, seq, sec , &tt );

  // extract key features for HMM detection
  //         0        1      2       3      4    
  //         totpow , shift  emg     sigma  h3
  // first three ftrs (of 5)
  // std::vector<int> ex = { 1,2 };
  // X = extract( X , ex );

  
  // optionally,  dump to stdout
  // dump( X , tt );

  std::map<std::string,std::set<interval_t> > anns = event_heuristic( X , tt , hp );
  
  std::map<std::string,std::set<interval_t> >::const_iterator aa = anns.begin();
  while ( aa != anns.end() )
    {
      annot_t * annot = parent->annotations->add( aa->first );
      const std::set<interval_t> & evts = aa->second;
      std::set<interval_t>::const_iterator ee = evts.begin();
      while ( ee != evts.end() )
	{
	  annot->add( "." , *ee , "." );
	  ++ee;
	}
      ++aa;
    }


  
  //
  // add channels
  //

  if ( add_chs && hp.do_nrem )
    add_channels( X[0], tt[0] , ch_prefix + "nrem_" );

  if ( add_chs && hp.do_rem )
    add_channels( X[1], tt[1] , ch_prefix + "rem_" );

}

int arousals_t::annotate( const int state ,
			  std::vector<std::vector<double> > & tt ,
			  std::vector<std::vector<Eigen::VectorXd> > & X , 
			  std::vector<std::vector<int> > & paths ,
			  std::vector<Eigen::MatrixXd> & posteriors ,
			  int st_arousal, int st_artifact, int st_baseline,
			  const std::string & aname )
{
  if ( X.size() == 0 || X[0].size() == 0 ) return 0;

  const int nstates = state == 0 ? 3 : 2 ;
  
  // summary stats per class (over all segments)
  const int nftr = X[0][0].size();
  const int nseq = X.size();
  
  std::vector<int> counts( nstates );
  Eigen::MatrixXd means = Eigen::MatrixXd::Zero( nstates , nftr );
  
  for (int i=0; i<nseq; i++)
    for (int e=0; e<X[i].size(); e++)
      {
	int s = paths[i][e];
	counts[ s ]++;
	for (int j=0; j<nftr; j++)
	  means( s, j ) += X[i][e][j];
      }
  
  int tot = 0;
  for (int s=0; s<nstates; s++)
    {
      tot += counts[s];
      for (int j=0; j<nftr; j++)
	means( s , j ) /= (double)counts[s];
    }
  
  for (int s=0; s<nstates; s++)
    std::cerr << "state " << s << " " << counts[s] << " " << counts[s] / (double)tot << "\n";
  
  for (int j=0; j<nftr; j++)
    {
      for (int s=0; s<nstates; s++)
	std::cerr << "\t" << means(s,j) ;
      std::cerr << "\n";
    }


  //
  // AASM rules
  //

  if ( st_artifact != -1 )
    {
      add_annot( st_artifact , paths , tt , aname+"_artifact" , "NREM" );
    }
  
  if ( st_arousal  != -1 )
    {
      add_annot( st_arousal , paths , tt , aname+"_arousal" , "NREM" );
    }
 
  int na = 0;
    
  return na;
}

void arousals_t::add_channels( const std::vector<std::vector<Eigen::VectorXd> > & X ,
			       const std::vector<std::vector<double> > & tt ,
			       const std::string & ch_prefix )
{
  // make channel w/ 0.5 Hz
  const int nr = parent->header.nr;
  const double rs = parent->header.record_duration;
  const int nsamples = 2 * rs;

  std::vector<std::string> chlab = { "totpwr", "shift", "emg", "sigma", "h3" };
    
  for (int cidx = 0; cidx<5; cidx++)
    {
      
      // make a new signal (2Hz)
      std::string label = ch_prefix;
      if ( ! label.empty() && label[ label.size() - 1 ] != '_' ) label += "_";
      label += chlab[cidx];
      parent->init_signal( label , 2 );
      
      const int slot = parent->header.signal( label );
      slice_t slice( *parent , slot ,  parent->timeline.wholetrace() );
      const std::vector<uint64_t> * tp = slice.ptimepoints();
      const int np = tp->size();
      if ( np != nr * nsamples )
	Helper::halt( "internal error in arousals_t::add_channels()" );
      
	      // get t-track of new signal in seconds
	      std::vector<double> sec( np );
	      for (int i=0; i<np; i++)
		sec[i] = (*tp)[i] / (double)globals::tp_1sec;
      
      // new data
      std::vector<double> xx( np , 0 );
      
	      // get original times, and values
	      for (int sq=0; sq<tt.size(); sq++)
		{
		  // within this sequence (contig), interpolate
		  const int ne = tt[sq].size();
		  if ( ne == 0 ) continue;

		  std::vector<double> xx0(ne);
		  for (int i=0; i<ne; i++)
		    xx0[i] = X[sq][i][cidx];

		  const double tmin = tt[sq][0];
		  const double tmax = tt[sq][ tt[sq].size() - 1 ];

		  tk::spline spline;
		  if ( ne > 2 )
		    spline.set_points( tt[sq] , xx0 );

		  // interpolation at this resolution
		  for (int i = 0 ; i < np; i++ )
		    {
		      if ( sec[i] > tmax ) break;
		      if ( sec[i] < tmin ) continue;

		      if ( ne == 1 )
			xx[i] = xx0[0];
		      else if ( ne == 2 )
			{
			  const double w = ( sec[i] - tmin ) / ( tmax - tmin );
			  xx[i] = xx0[0] + w * ( xx0[1] - xx0[0] );
			}
		      else
			xx[i] = spline( sec[i] );
		    }

		  // next contig
		}
      
      // add signal
      parent->update_signal( slot , &xx );
    } 
}
			       


void arousals_t::add_annot( int idx ,   			  
			    std::vector<std::vector<int> > & paths ,
			    std::vector<std::vector<double> > & tt ,			    
			    const std::string & class_label ,
			    const std::string & inst_label )
			    
{
  // add the annot.
  annot_t * a = parent->annotations->add( class_label );

  const int nseqs = paths.size();

  for (int sq = 0; sq < nseqs; sq++)
    {
      int prior = -1;
      const int nepochs = paths[sq].size();
      for (int i=0; i<nepochs; i++)
	{
	  // starting a new stretch?
	  if ( prior == -1 && paths[sq][i] == idx )
	    prior = i;
	  else if ( prior != -1 && paths[sq][i] != idx )
	    {
	      // ending an existing stretch
	      double start = tt[sq][prior];
	      double stop = tt[sq][i];
	      a->add( inst_label , interval_t( start * globals::tp_1sec , stop * globals::tp_1sec ) , "." );
	      prior = -1; 
	    }
	}
      // end segment
      if ( prior != -1 )
	{
	  double start = tt[sq][prior];
	  double stop = tt[sq][nepochs-1] + 0.5; // if using fixed 0.5s steps
	  a->add( inst_label , interval_t( start * globals::tp_1sec , stop * globals::tp_1sec ) , "." );
	}
    }
  
}



std::vector<std::vector<std::vector<Eigen::VectorXd> > > arousals_t::assemble( const Eigen::MatrixXd & Xftr ,
									       const std::vector<int> & state ,
									       const std::vector<int> & seq ,
									       const std::vector<double> & sec , 
									       std::vector<std::vector<std::vector<double> > > * tt )
{

  
  // assume , for NR and R separatetely:
  // return --> [ state (NR/R) ][ contig ][ epoch x feature ]
  // also, time-track tt for each contig (to make annots)

  // state x seq x epoch x ftr
  std::vector<std::vector<std::vector<Eigen::VectorXd> > > X;
  
  tt->clear();
  
  // enumerate state-specific sleep contigs [ state ][ contig ] -> epoch-list
  std::map<int,std::map<int,std::vector<int> > > contigs;
  
  const int n = state.size(); // num epochs
  const int nftr = Xftr.cols(); // shoiuld always be 2 + 1 ( + 2 ) = 5
   
  // extract and split by stage/contig

  for (int i=0; i<n; i++)
    {
      // skip wake
      if ( state[i] == 2 ) continue;
      // track which epochs
      contigs[ state[i] ][ seq[i] ].push_back( i );
    }

  // size up X
  X.resize( 2 ); // NR, R
  X[0].resize( contigs[0].size() );
  X[1].resize( contigs[1].size() );

  // size up tt
  tt->resize( 2 ); 
  (*tt)[0].resize( contigs[0].size() );
  (*tt)[1].resize( contigs[1].size() );
  
	  for (int st=0; st<2; st++)
	    {
	      int sq = 0;
	      for (std::map<int,std::vector<int> >::const_iterator cc = contigs[st].begin();
		   cc != contigs[st].end(); ++cc, ++sq)
		{
		  const std::vector<int> & ee = cc->second;
		  const int nce = ee.size();

		  // tt
		  (*tt)[st][sq].resize(nce);
		  for (int i=0; i<nce; i++)
		    (*tt)[st][sq][i] = sec[ ee[i] ] ;

		  // data
		  X[st][sq].resize( nce );
		  for (int i=0; i<nce; i++)
		    X[st][sq][i] = Xftr.row( ee[i] ) ;

		}
	    }

  // [st][contig][epoch x ftr]
  return X;
}



void arousals_t::dump( const std::vector<std::vector<std::vector<Eigen::VectorXd> > > & X ,
		       const std::vector<std::vector<std::vector<double> > > & tt ) const
{

  for (int st=0; st<2; st++)
    {
      const int nc = tt[st].size();
      for (int sq=0; sq<nc; sq++)
	{
	  const int ne = tt[st][sq].size();
	  //	  std::cout << " st=" << st << ", sq=" << sq << ", ne=" << ne << "\n";
	  
	  for (int i=0; i<ne; i++)
	    {
	      std::cout << " " << i << "\t" << tt[st][sq][i];
	      const Eigen::VectorXd & XX = X[st][sq][i];
	      for (int j=0; j<XX.size(); j++)
		std::cout << "\t" << XX(j);
	      std::cout << "\n";		
	    }
	}
    }
  

}



void arousals_t::build_ftr_matrix( edf_t & edf ,
				   const signal_list_t & eeg_signals ,
				   const signal_list_t & emg_signals ,
				   Eigen::MatrixXd & Xeeg,
				   Eigen::MatrixXd & Xemg,
				   std::vector<int> * state ,
				   std::vector<int> * sequence ,
				   std::vector<double> * sec ,
				   const int wake_bridge ,
				   const double epoch_inc
				   )
					      
{

  // state: annotate whether start of each window is NR, R or W (0,1,2)
  // seq: track contiguous, within-state sequences i.e. (0,0,0,1,1,1,2,2,2,etc)

  // three features
  //   artifact          | Neeg x 1 
  //   EEG arousal       | core      = Neeg x 3 [ (theta+alpha+beta) / delta ]
  //                     | subtyping = Neeg x 2 [ H3, rel-sigma ]
  //   EMG burst         | Nemg x 1 

  // total EEG features = 4 x Neeg
  // EMG = 1 x Nemg
  
  // if multiple channels, average before Z scoring to create single
  // final signals
  
  const int ns_eeg = eeg_signals.size();
  const int ns_emg = emg_signals.size();

  const int nftr_eeg = ns_eeg * 4;
  const int nftr_emg = ns_emg * 1;

  const int nftr = nftr_eeg + nftr_emg;
  
  state->clear();
  sequence->clear();
  sec->clear();
  
  int ne = edf.timeline.first_epoch();

  // allocate storage for features, per epoch
  Xeeg = Eigen::MatrixXd::Zero( ne , nftr_eeg );
  Xemg = Eigen::MatrixXd::Zero( ne , nftr_emg );

  std::unique_ptr<FFT> eeg_fft;
  if ( ns_eeg )
    {
      const int index_length = std::lround( sr * edf.timeline.epoch_length() );
      if ( index_length < 1 )
	Helper::halt( "invalid EEG FFT window length" );
      eeg_fft.reset( new FFT( index_length , index_length , sr , FFT_FORWARD , WINDOW_TUKEY50 ) );
    }

  // ensure staging is present in SleepStage
  edf.annotations->make_sleep_stage( edf.timeline );

  annot_t * staging = edf.annotations->find( "SleepStage" );

  if ( staging == NULL )
    Helper::halt( "no staging information present" );
  
  // track state-specific contigs
  int prior_state = 2;  // 0/1 = NR/R
  int seq_nr = -1;
  int seq_r = -1;
  bool have_prior_feature_sec = false;
  double prior_feature_sec = 0;

  // Optionally count the first N scored wake epochs after NR/R as the
  // preceding sleep state, so sleep-transition arousals can start in W.
  // SleepStage W annotations may be merged across many epochs, so use elapsed
  // time from the first W window rather than annotation identity.
  const double wake_bridge_sec = wake_bridge * globals::default_epoch_len;
  int last_sleep_state = 2;
  bool in_wake_run = false;
  int wake_bridge_state = 2;
  double wake_run_start_sec = 0;
  
  // get EMG thresholds [ whole-night ] 
  const double emg_SD_threshold = 9;  // 8-10 good range
  std::vector<double> emg_th( ns_emg );
  for (int s=0; s<ns_emg; s++)
    {
      // pull whole signal
      slice_t slice( edf , emg_signals(s) , edf.timeline.wholetrace() );
      const std::vector<double> * d = slice.pdata();
      const int n = d->size();
      // get median absolute deviation (MAD)
      const double median = MiscMath::median( *d );
      
      std::vector<double> dd( n );
      for (int i=0; i<n; i++) dd[i] = fabs( (*d)[i] - median );
      const double MAD = MiscMath::median( dd );
      const double sigma = 1.4826 * MAD;
      if ( sigma < 1e-6 ) Helper::halt( "low EMG amplitude" );
      emg_th[s] = sigma * emg_SD_threshold;
    }
  
  
  // iterate over epochs

  int eidx = -1;
  
  while ( 1 )
    {
      // next epoch
      int epoch = edf.timeline.next_epoch();

      // all done?
      if ( epoch == -1 ) break;

      ++eidx;
      
      interval_t interval = edf.timeline.epoch( epoch );
      
      // which stage does this start in?
      annot_map_t events = staging->extract( interval );

      int scored_st = 2; // wake/unknown
      bool scored_wake = false;
      if ( events.size() >= 1 )
	{
	  // get first annotaiton, i.e. for start, if >1
	  const instance_idx_t & idx = events.begin()->first;
	  if ( idx.id == "N1" || idx.id == "N2" || idx.id == "N3" || idx.id == "NREM4" || idx.id == "NR" )
	    scored_st = 0;
	  else if ( idx.id == "R" )
	    scored_st = 1;
	  else if ( idx.id == "W" )
	    scored_wake = true;
	}	  

      int st = scored_st;
      if ( scored_st == 0 || scored_st == 1 )
	{
	  last_sleep_state = scored_st;
	  in_wake_run = false;
	  wake_bridge_state = 2;
	  wake_run_start_sec = 0;
	}
      else if ( scored_wake && wake_bridge > 0 )
	{
	  if ( ! in_wake_run )
	    {
	      in_wake_run = true;
	      wake_bridge_state = last_sleep_state;
	      wake_run_start_sec = interval.start_sec();
	    }

	  if ( ( wake_bridge_state == 0 || wake_bridge_state == 1 ) &&
	       interval.start_sec() - wake_run_start_sec < wake_bridge_sec )
	    st = wake_bridge_state;
	}
      else
	{
	  last_sleep_state = 2;
	  in_wake_run = false;
	  wake_bridge_state = 2;
	  wake_run_start_sec = 0;
	}

      // Time-stamp each feature row at the center of its analysis window.
      const double feature_sec = interval.start_sec() + 0.5 * interval.duration_sec();
      const bool time_discontinuity =
	have_prior_feature_sec && feature_sec - prior_feature_sec > 1.5 * epoch_inc;

      // a new NR or R sequence?
      if ( st != prior_state || time_discontinuity )
	{
	  if ( st == 0 ) ++seq_nr;	    
	  else if ( st == 1 ) ++seq_r;
	}

      // for next time round
      prior_state = st;
      prior_feature_sec = feature_sec;
      have_prior_feature_sec = true;
      
      // track
      state->push_back( st );      
      if ( st == 0 ) sequence->push_back( seq_nr );
      else if ( st == 1 ) sequence->push_back( seq_r );
      else sequence->push_back( -1 ); // we don't care: W or ?
      sec->push_back( feature_sec );
      
      // if non-sleep, can skip
      if ( st == 2 ) continue;
      
	      if ( ns_eeg )
		{
		  // get EEGs
		  eigen_matslice_t mslice( edf , eeg_signals , interval );
		  const Eigen::MatrixXd & X = mslice.data_ref();

		  // derive features
			  Xeeg.row( eidx ) = calc_eeg_ftrs( X , *eeg_fft );
		}

      if ( ns_emg )
	{
	  eigen_matslice_t mslice_emg( edf , emg_signals , interval );
	  const Eigen::MatrixXd & X = mslice_emg.data_ref();
	  Xemg.row( eidx ) = calc_emg_ftrs( X , emg_th );	  
	}

    }
  
}



Eigen::VectorXd arousals_t::calc_emg_ftrs( const Eigen::MatrixXd & X , const std::vector<double> & thr )
{

  // thr is threshold for clipping EMG signal, base don median-of-median
  const int ns = X.cols();
  
  // currently, only a single ftr per EMG
  Eigen::VectorXd F = Eigen::VectorXd::Zero( 1 * ns ) ;
  const double eps = 1e-8;
  
  int fi = -1;
  for (int s=0; s<ns; s++)
    {
      // RMS of clipped signals (based on whole recording threshold)
      const Eigen::VectorXd & d = X.col(s);
      const int n = d.size();
      const double t = thr[s];
      double rms = 0;
      for (int i=0; i<n; i++)
	{
	  if ( d[i] < -t || d[i] > t ) rms += t*t;
	  else rms += d[i]*d[i];
	}
      rms /= (double)n;
      rms = sqrt( rms );
      F(++fi) = log( rms + eps );
    }

  return F;
}

Eigen::VectorXd arousals_t::calc_eeg_ftrs( const Eigen::MatrixXd & X , FFT & fftseg )
{

  // currently, get 4 ftrs from each EEG: 2 for detector, 2 for subtyping
  Eigen::VectorXd F = Eigen::VectorXd::Zero( 4 * X.cols() );

  // per channel

  const double eps = 1e-12;

  int fi = 0;
  
  for (int s=0; s<X.cols(); s++)
    {
      fftseg.apply( X.col(s).data() , X.rows() );

      double p_tot = 0;   // 4-30
      double p_delta = 0; // 1-4
      double p_theta = 0; // 4-8
      double p_alpha = 0; // 8-10
      double p_sigma = 0; // 10-16
      double p_beta = 0;  // 16-30
      double p_shift = 0; // log((theta+alpha+beta)/delta)

      double hjorth3 = 0;

      for (int f = 0; f < fftseg.cutoff; f++)
	{
	  double frq = fftseg.frq[f] ;
	  if ( frq < 1 || frq >= 30 ) continue;
	  if ( frq >= 4 ) p_tot += fftseg.X[f] ;
	  if      ( frq >= 16 ) p_beta  += fftseg.X[f];
	  else if ( frq >= 10 ) p_sigma += fftseg.X[f];
	  else if ( frq >= 8  ) p_alpha += fftseg.X[f];
	  else if ( frq >= 4  ) p_theta += fftseg.X[f];
	  else                  p_delta += fftseg.X[f];
	}

      // transform to log/relative power

      p_shift = log( ( p_theta + p_alpha + p_beta + eps ) / ( p_delta + eps ) );
      p_sigma = log( ( p_sigma + eps ) / ( p_tot + eps ) );
      p_tot = log( p_tot + eps );

	      // Hjorth
	      double h1, h2, h3;
	      bool okay = hjorth( X.col(s) , &h1, &h2, &h3 , true );
	      if ( ! okay ) h3 = 0;
      
      F(fi++) = p_tot;
      F(fi++) =	p_shift;
      F(fi++) =	p_sigma;
      F(fi++) =	h3;
    }
  
  return F;
}



bool arousals_t::hjorth( const Eigen::VectorXd & x ,
			 double * activity ,
			 double * mobility ,
			 double * complexity ,
			 const bool mean_center ) const
{

  const int n = x.size();
  if (n < 3) return false;  // Need at least 3 points for second derivative

  const Eigen::VectorXd dxV = x.tail(n - 1) - x.head(n - 1);
  const Eigen::VectorXd ddxV = dxV.tail(n - 2) - dxV.head(n - 2);

  const double eps = 1e-12;  // Small epsilon to avoid division by zero

  const double mx2 = (mean_center ? (x.array() - x.mean()).matrix().squaredNorm()
                                  : x.squaredNorm()) / double(n);
  const double mdx2 = dxV.squaredNorm() / double(n - 1);
  const double mddx2 = ddxV.squaredNorm() / double(n - 2);

  if (mx2 < eps || mdx2 < eps) return false;  // Avoid division by zero

  *activity = mx2;
  *mobility = sqrt(mdx2 / mx2);
  *complexity = sqrt((mddx2 * mx2) / (mdx2 * mdx2));

  if (!Helper::realnum(*activity)) return false;
  if (!Helper::realnum(*mobility)) return false;
  if (!Helper::realnum(*complexity)) return false;
  
  return true;
}


void arousals_t::init_kmeans_hmm( gaussian_hmm_t & hmm ,
				  const std::vector<std::vector<Eigen::VectorXd>> & sequences )
{
  
  // get single matrix
  Eigen::MatrixXd X = stack_sequences(sequences, hmm.dim()); // N_total x M
  
  // k-means w/ K = number of HMM states
  kmeans_result_t km = kmeans(X, hmm.n_states(), /*max_iters=*/50, /*tol=*/1e-4, /*seed=*/123);

  // HMM expects mu as M x N (columns = state means), kmeans gives K x M
  Eigen::MatrixXd Mu = km.centroids.transpose();  // (M x K)

  // Initialize means; covariances stay as current/default (identity)
  hmm.set_emission(Mu);
  
  // Rough π from first sample of each sequence
  Eigen::VectorXd pi = Eigen::VectorXd::Zero(hmm.n_states());
  int n_seq = static_cast<int>(sequences.size());
  int offset = 0;
  for (int s = 0; s < n_seq; ++s) {
    if (!sequences[s].empty()) {
      int k = km.labels(offset);  // label of first sample in sequence s
      pi(k) += 1.0;
      offset += sequences[s].size();
    }
  }
  if (pi.sum() > 0) pi /= pi.sum();
  hmm.set_initial(pi);
  
  // Rough A from label transitions in flat X
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(hmm.n_states(), hmm.n_states());
  for (int i = 0; i + 1 < km.labels.size(); ++i) {
    int k1 = km.labels(i);
    int k2 = km.labels(i+1);
    A(k1,k2) += 1.0;
  }

  // Normalize rows
  for (int i = 0; i < A.rows(); ++i) {
    double s = A.row(i).sum();
    if (s > 0) A.row(i) /= s;
  }
  
  hmm.set_transition(A);
    
}

std::vector<std::vector<std::vector<Eigen::VectorXd> > >
arousals_t::extract( const std::vector<std::vector<std::vector<Eigen::VectorXd> > > & X , const std::vector<int> & ex )
  
{
  // ugly...
  std::vector<std::vector<std::vector<Eigen::VectorXd> > > R = X;  
  for (int st=0; st<X.size(); st++)
    for (int sq=0; sq<X[st].size(); sq++)
      for (int i=0; i<X[st][sq].size(); i++)
	R[st][sq][i] = X[st][sq][i](ex);
    return R;
}


void arousals_t::map_states( const Eigen::MatrixXd & mu , int* st_arousal ,int * st_artifact , int* st_baseline , bool nrem ) const
{
  // for NREM:
  // pick artifact (movement) first
  // then pick arousal among the remaining two
  // leave a “no clear arousal/artifact” escape hatch if separation is weak.

  // for now assume single EEG, single EMG so alpha/shift/EMG
  // --> otherwise, make collapsed (mean) matrix
  
  // mu - rows    = metrics ; 
  //      cols    = states:

  //  std::cout << "MMMuuu\n" << mu << "\n";

  Eigen::VectorXd alpha = mu.row(0);
  Eigen::VectorXd shift = mu.row(1);
  Eigen::VectorXd emg = mu.row(2);
  
  // max EMG
  int max_emg;
  emg.maxCoeff(&max_emg);
  //  std::cout << "max_emg = " << max_emg << "\n";

  // second largest
  int min_emg;
  emg.minCoeff(&min_emg);

  int second_emg = 3 - max_emg - min_emg;
    
  double delta_emg = emg(max_emg) - emg(second_emg);

  // artifact requires:
  //  emg(max_emg) > 0.5
  //  delta(max_emg) > 0.5
  //  shift not very high, i.e. less than 0.3 from max shift
  int max_shift;
  shift.maxCoeff( &max_shift );
  double delta_shift = shift(max_shift) - shift(max_emg);

  bool has_artifact_class = true;
  if ( max_emg == min_emg ) has_artifact_class = false;
  if ( emg(max_emg) < 0.5 ) has_artifact_class = false;
  if ( delta_emg < 0.5 ) has_artifact_class = false;
  if ( delta_shift < 0.3 ) has_artifact_class = false;

  // if no reliable artifact: treat as two-class model

  // from non-artifact classes:
  //   heuristic: favor frequency shift, with some contribution from alpha and EMG

  // Sar​(k)=bk​+0.5ak​+0.3mk

  double s0 = shift(0) + 0.5 * alpha(0) + 0.3 * emg(0);
  double s1 = shift(1) + 0.5 * alpha(1) + 0.3 * emg(1);
  double s2 = shift(2) + 0.5 * alpha(2) + 0.3 * emg(2);
  
  int max_s = 0;  
  if ( s1 > s0 ) max_s = 1;
  if ( s2 > s1 && s2 > s0 ) max_s = 2;

  int min_s = 0;  
  if ( s1 < s0 ) min_s = 1;
  if ( s2 < s1 && s2 < s0 ) min_s = 2;

  // clver hack - make this whatever min/max are not, given 0/1/2 encoding
  int mid_s = 3 - min_s - max_s;

  // std::cout << " s0,s1,s2 = " << s0 << " " << s1 << " " << s2 << "\n";  
  // std::cout << "min, mid, max s= " << min_s << " " << mid_s << " " << max_s << "\n";
  
  // arousal - max_s (that is not the artifact class;
  *st_arousal = max_s;
  *st_baseline = min_s; // baseline state w/ min arousal

  // ... unless we need to remove the artifact class out
  
  if ( has_artifact_class )
    {
      // ensure not artifact class
      if ( max_emg == max_s )
	*st_arousal = mid_s;

      // baseline is the remaining class
      // clever hack:
      *st_baseline = 3 - *st_arousal - max_emg;
      
    }
      
  // are arousal and baseline meaningfully separated?
  double delta_shift2 = shift(*st_arousal) - shift(*st_baseline);
  double delta_emg2 = emg(*st_arousal) - emg(*st_baseline);

  // NREM: cortical arousal only, do not require EMG
  // REM: requires both EEG and EMG arousals, as per AASM
  bool has_arousal_state = nrem 
    ? delta_shift2 > 0.4
    : delta_shift2 > 0.4 && delta_emg2 > 0.5;

  // std::cout << "shift arousal criteria: " << shift(*st_arousal) << " - " << shift(*st_baseline) << " = " << delta_shift2 << "\n";
  // std::cout << "EMG arousal criteria: " << emg(*st_arousal) << " - " << emg(*st_baseline) << " = " << delta_emg2 << "\n";
   
  // std::cout << "has artifact " << has_artifact_class << " " 
  // 	    << ( has_artifact_class ? max_emg : -1 )
  // 	    << "\n";
  
  // std::cout << "has arousal " << has_arousal_state << " " 
  // 	    << ( has_arousal_state ? *st_arousal : -1 )
  // 	    << "\n";
    
  if ( has_artifact_class ) *st_artifact = max_emg;
  else *st_artifact = -1;
  if ( ! has_arousal_state ) *st_arousal = -1;
			       
  std::cout << "(artifact / arousal / baseline ) = " << *st_artifact << " / " << *st_arousal << " / " << *st_baseline << "\n";

  
}



Eigen::MatrixXd arousals_t::process_ftr_matrix( Eigen::MatrixXd * Xeeg , Eigen::MatrixXd * Xemg , const std::vector<int> & st , const std::vector<double> & sec , const double emg_median_sec , const double epoch_inc )
{
  // input:
  //   per-channel eeg: tot-pow, arousal-shift, [ sigma, h3 ]
  //   per-channel emg: rms
  //   log-scaling/relative power already done

  // output:
  //  subtract local baseline for each feature
  //  average over channels
  //  robust norm based on sleep epochs
  //  return 5 values in single matrix 

  if ( Xeeg->rows() != Xemg->rows() )
    Helper::halt( "internal error in arousals_t::process_ftr_matrix()" );
  if ( (int)sec.size() != Xeeg->rows() )
    Helper::halt( "internal time-track error in arousals_t::process_ftr_matrix()" );
  
  const int nch_eeg = Xeeg->cols() / 4;
  const int nch_emg = Xemg->cols();

  const int neeg = Xeeg->cols();

  // median filter (including *all* epochs, wake+sleep)
  // 120 seconds, converted to feature rows based on the current inc
  const double baseline_sec = 120.0;
  int baseline_n = epoch_inc > 0 ? (int)( baseline_sec / epoch_inc + 0.5 ) : 1;
  if ( baseline_n < 1 ) baseline_n = 1;

  // For the long baseline, bridge short masked gaps; otherwise tiny retained
  // islands would get effectively no useful baseline estimate.  Only split
  // the baseline filter across gaps larger than the baseline window itself.
  auto subtract_segmented_median = [&]( Eigen::MatrixXd * M )
    {
      for (int s=0; s<M->cols(); s++)
	{
	  int start = 0;
	  while ( start < M->rows() )
	    {
	      int stop = start;
	      while ( stop + 1 < M->rows() && sec[stop+1] - sec[stop] <= baseline_sec )
		++stop;

	      const int n = stop - start + 1;
	      const int local_baseline_n = std::min( baseline_n , n );
	      Eigen::VectorXd x = M->col(s).segment( start , n );
	      M->col(s).segment( start , n ) = x - eigen_ops::median_filter( x , local_baseline_n );
	      start = stop + 1;
	    }
	}
    };

  subtract_segmented_median( Xeeg );
  subtract_segmented_median( Xemg );
  
  // merge channels  : note column re-ordering/splice in detector ftrs
  //                   to first 3 cols, vs subtyping ftrs in last 2

  Eigen::MatrixXd X = Eigen::MatrixXd::Zero( Xeeg->rows() , 5 );
  for (int s=0; s<nch_eeg; s++)
    {
      X.col( 0 ) += Xeeg->col( s*4 + 0);
      X.col( 1 ) += Xeeg->col( s*4 + 1);
      X.col( 3 ) += Xeeg->col( s*4 + 2); // nb. is correct Xeeg[,2] -> X[,3]
      X.col( 4 ) += Xeeg->col( s*4 + 3); // nb. same, Xeeg[,3] -> X[,4]
    }
  
  for (int s=0; s<nch_emg; s++)
    X.col( 2 ) += Xemg->col( s ); // EMG goes into X[,2]

  if ( nch_eeg > 1 )
    {
      X.col( 0 ) = X.col( 0 ).array() / (double)nch_eeg;
      X.col( 1 ) = X.col( 1 ).array() / (double)nch_eeg;
      X.col( 3 ) = X.col( 3 ).array() / (double)nch_eeg;
      X.col( 4 ) = X.col( 4 ).array() / (double)nch_eeg;
    }
  
  // normalize (based on sleep epochs only)
  if ( nch_emg > 1 )
    X.col( 2 ) = X.col( 2 ).array() / (double)nch_emg;
  
  // now -> X = EEG-pwr EEG-shift EMG-rms [ EEG-spindle EEG-h3 ]
  //  a single, averaged value per channel


  // finally, per-stage (NREM vs REM) robust norm: REM has systematically
  // lower baseline EMG (atonia), so a single NREM+REM pooled median/MAD
  // would under-represent REM's own variability and mute genuine REM EMG
  // rises; normalize NREM and REM epochs against their own stage statistics
  for (int i=0; i<X.cols(); i++)
    {
      const Eigen::VectorXd z_nrem = robust_mad_norm( X.col(i) , st , 0 );
      const Eigen::VectorXd z_rem  = robust_mad_norm( X.col(i) , st , 1 );
      Eigen::VectorXd z( X.rows() );
      for (int e=0; e<X.rows(); e++)
	z(e) = ( st[e] == 1 ) ? z_rem(e) : z_nrem(e);
      X.col(i) = z;
    }

  // De-spike the normalized EMG feature before thresholding/clipping.  This
  // leaves EEG shift/totpow/sigma/H3 untouched.
  if ( nch_emg > 0 && emg_median_sec > 0 )
    {
      int n = epoch_inc > 0 ? (int)( emg_median_sec / epoch_inc + 0.5 ) : 1;
      if ( n < 1 ) n = 1;
      if ( n > 1 ) X.col(2) = eigen_ops::median_filter( X.col(2) , n );
    }

  // and clip; keep headroom above artifact thresholds (4-5 z) so diagnostics
  // and added channels can still grade artifact severity.
  const double zth = 12;
  for (int i=0; i<X.cols(); i++)
    X.col(i) = X.col(i).cwiseMax(-zth).cwiseMin(zth);
  
  // upweight the arousal-shift track
  //  X.col(1) = X.col(1).array() ;
  
  // all done, give back
  return X;
}





Eigen::VectorXd arousals_t::robust_mad_norm(const Eigen::VectorXd & x,
                                            const std::vector<int> & st,
                                            const int target_state)
{

  const int n = x.size();

  if ((int)st.size() != n)
    Helper::halt( "robust_mad_norm: x and st must have same length" );

  // collect values for the target stage only (0=NREM, 1=REM), so that e.g.
  // REM's own (typically lower) EMG baseline/variability is used to scale
  // REM values, rather than a stage-pooled scale
  std::vector<double> vals;
  vals.reserve(n);
  for (int i = 0; i < n; ++i)
    {
      if (st[i] == target_state && std::isfinite(x[i])) // this stage only, skip NaNs
	vals.push_back(x[i]);
    }

  // this stage may simply be absent from the recording (e.g. no REM at
  // all); the returned z-scores are unused for epochs outside this stage,
  // so return neutral (all-zero) values rather than halting
  if (vals.empty())
    return Eigen::VectorXd::Zero(n);

  const double median = MiscMath::median(vals);

  // MAD over this stage only
  std::vector<double> dev;
  dev.reserve(vals.size());
  for (double v : vals)
    dev.push_back(std::fabs(v - median));

  const double MAD = MiscMath::median(dev);

  double sigma = 1.4826 * MAD;
  if (sigma <= 0.0)
    sigma = 1.0;  // fallback: avoids division by zero

  // normalize all samples (wake + sleep) using this stage's median/sigma
  Eigen::VectorXd z(n);
  for (int i = 0; i < n; ++i)
    z[i] = (x[i] - median) / sigma;

  return z;
}


std::map<std::string,std::set<interval_t> > arousals_t::event_heuristic( const std::vector<std::vector<std::vector<Eigen::VectorXd> > > & X ,
									 const std::vector<std::vector<std::vector<double> > > & tt ,
									 const arousal_heuristic_params_t & p )
{
  // metrics 0 tot-pwr
  //         1 arousal-shift
  //         2 EMG
  //         3 sigma
  //         4 H3

  std::map<std::string,std::set<interval_t> > ret;

  const int idx_pwr = 0;
  const int idx_shift = 1;
  const int idx_emg = 2;
  const int idx_sigma = 3;
  const int idx_h3 = 4;

  // gap-merge tolerance in samples, based on the actual epoch grid spacing
  // (rather than assuming a fixed 2 Hz / 0.5s grid)
  const int max_gap_samples = std::max( 1 , (int)std::lround( p.merge_gap_sec / p.epoch_inc ) );

  // NREM (st=0) and REM (st=1); each stage is independently toggled,
  // thresholded and reported (writer stratifier SS = NR/R)
  for (int st=0; st<2; st++)
    {
      if ( st == 0 && ! p.do_nrem ) continue;
      if ( st == 1 && ! p.do_rem  ) continue;

      const std::string stg_lab = st == 0 ? "nrem" : "rem";
      const bool is_rem = st == 1;

      // spindle (sigma) veto: NREM guards against spindle contamination;
      // REM defaults are permissive (effectively disabled) since spindles
      // are a NREM phenomenon
      const double th_sigma_veto  = is_rem ? p.sigma_veto_rem  : p.sigma_veto_nrem;
      const double th_sigma_veto2 = is_rem ? p.sigma_veto2_rem : p.sigma_veto2_nrem;

      // check ftrs inside/outside of artifact regions
      Eigen::VectorXd ftr_art( 5 ), ftr_nonart( 5 );
      Eigen::VectorXd ftr_baseline( 5 ), ftr_arousal( 5 ), ftr_uarousal( 5 );
      ftr_art.setZero(); ftr_nonart.setZero();
      ftr_baseline.setZero(); ftr_arousal.setZero(); ftr_uarousal.setZero();
      int n_art = 0 , n_nonart = 0, n_baseline = 0 , n_arousal = 0, n_uarousal = 0;

      // event counts
      int cnt_evts = 0 , cnt_uevts = 0, cnt_long_evts = 0, cnt_arts = 0;

      // verbose candidate-pipeline counts; denominators are raw candidate
      // intervals unless a variable name explicitly states otherwise
      int vcand_raw = 0;
      int vcand_dur_ok = 0;
      int vcand_merged = 0;
      int vcand_final_dur_ok = 0;
      int vcand_presleep_ok = 0;
      int vcand_rem_emg_confirmed = 0;
      int vcand_rem_emg_rejected = 0;

      // verbose threshold/mask diagnostics.  Window-mask proportions use
      // VDIAG_N_WIN as denominator; peak-block proportions use
      // VDIAG_N_LOCAL_SHIFT_PEAK as denominator.
      int vdiag_win = 0;
      int vdiag_high_shift_peak = 0;
      int vdiag_high_shift_grow = 0;
      int vdiag_high_sigma_peak = 0;
      int vdiag_high_sigma_grow = 0;
      int vdiag_high_emg = 0;
      int vdiag_high_h3 = 0;
      int vdiag_high_pwr = 0;
      int vdiag_art_emg = 0;
      int vdiag_art_h3 = 0;
      int vdiag_art_pwr = 0;
      int vdiag_artifact = 0;
      int vdiag_artifact_block = 0;
      int vdiag_local_shift_peak = 0;
      int vdiag_peak_block_artifact = 0;
      int vdiag_peak_block_sigma = 0;
      int vdiag_peak_accepted = 0;

      // REM-only: cortical-shift candidates that failed EMG confirmation
      int cnt_rem_cortical_only = 0;
      int cnt_long_rem_cortical_only = 0;
      int cnt_rem_emg_confirmed_long = 0;

      // mean event durations
      double dur_major = 0 , dur_micro = 0, dur_long = 0;

      // Keep the original max-dur bucket unchanged; optionally carry
      // longer candidates through the same cleanup path for separate output.
      const bool do_long = p.long_dur > p.max_dur;
      const double candidate_max_dur = do_long ? p.long_dur : p.max_dur;

      // within each contig
      const int nc = X[st].size();

      for (int c=0; c<nc; c++)
	{
	  const std::vector<Eigen::VectorXd> & D = X[st][c];
	  const int ne = D.size();
	  const int idx_detect = p.emg_mode ? idx_emg : idx_shift;
	  const std::string detect_label = p.emg_mode ? "emg" : "shift";

	  auto add_verbose_interval = [&]( const std::string & label ,
					   const double t0 ,
					   const double t1 )
	    {
	      if ( ! arousal_verbose_annots ) return;
	      if ( t1 <= t0 ) return;
	      ret[ arousal_verbose_prefix + stg_lab + "_" + label ].insert(
	        interval_t( globals::tp_1sec * t0 , globals::tp_1sec * t1 ) );
	    };

	  auto add_verbose_mask = [&]( const std::string & label ,
				       const std::vector<bool> & mask )
	    {
	      if ( ! arousal_verbose_annots ) return;
	      const std::vector<std::pair<int,int> > ints = mask_to_intervals( mask );
	      for (int ii=0; ii<ints.size(); ii++)
		add_verbose_interval( label ,
				      tt[st][c][ ints[ii].first ] ,
				      tt[st][c][ ints[ii].second ] + p.epoch_inc );
	    };

	  auto add_verbose_events = [&]( const std::string & label ,
					 const std::vector<std::pair<int,int> > & ev )
	    {
	      if ( ! arousal_verbose_annots ) return;
	      for (int ii=0; ii<ev.size(); ii++)
		add_verbose_interval( label ,
				      tt[st][c][ ev[ii].first ] ,
				      tt[st][c][ ev[ii].second ] + p.epoch_inc );
	    };

	  std::vector<bool> high_detect_peak( ne , false );
	  std::vector<bool> high_detect_grow( ne , false );
	  std::vector<bool> high_sigma_peak( ne , false );
	  std::vector<bool> high_sigma_grow( ne , false );
	  std::vector<bool> high_emg( ne , false );
	  std::vector<bool> high_h3( ne , false );
	  std::vector<bool> high_pwr( ne , false );
	  std::vector<bool> art_emg( ne , false );
	  std::vector<bool> art_h3( ne , false );
	  std::vector<bool> art_pwr( ne , false );

	  // flag artifact annotations separately from the stricter mask used to
	  // block candidate peaks/growth.
	  std::vector<bool> artifact( ne , false );
	  std::vector<bool> artifact_block( ne , false );
	  for (int i=0; i<ne; i++)
	    {
	      const Eigen::VectorXd & ftr = D[i];

	      bool is_artifact = false;
	      bool is_artifact_block = false;

	      high_detect_peak[i]  = ftr(idx_detect) >= p.shift_peak;
	      high_detect_grow[i]  = ftr(idx_detect) >= p.shift_hysteresis;
	      high_sigma_peak[i] = ftr(idx_sigma) > th_sigma_veto;
	      high_sigma_grow[i] = ftr(idx_sigma) > th_sigma_veto2;
	      high_emg[i]        = ftr(idx_emg)   > p.th_emg_artifact;
	      high_h3[i]         = ftr(idx_h3)    > p.th_h3_artifact;
	      high_pwr[i]        = ftr(idx_pwr)   > p.th_pwr_artifact;

	      // Very high EMG
	      if ( ftr(idx_emg) > p.th_emg_artifact && ftr(idx_shift) < 0.5 )
		{
		  is_artifact = true;
		  art_emg[i] = true;
		  if ( ftr(idx_shift) < 0.0 )
		    is_artifact_block = true;
		}

	      // High H3 (very noisy)
	      if ( ftr(idx_h3) > p.th_h3_artifact )
		{
		  is_artifact = true;
		  art_h3[i] = true;
		  if ( ftr(idx_shift) < 0.0 )
		    is_artifact_block = true;
		}

	      // implausible broadband tot-pwr
	      if ( ftr(idx_pwr) > p.th_pwr_artifact && ftr(idx_shift) < 0.5 )
		{
		  is_artifact = true;
		  art_pwr[i] = true;
		  if ( ftr(idx_shift) < 0.0 )
		    is_artifact_block = true;
		}

	      artifact[i] = is_artifact;
	      artifact_block[i] = is_artifact_block;

	      ++vdiag_win;
	      if ( high_detect_peak[i] ) ++vdiag_high_shift_peak;
	      if ( high_detect_grow[i] ) ++vdiag_high_shift_grow;
	      if ( high_sigma_peak[i] ) ++vdiag_high_sigma_peak;
	      if ( high_sigma_grow[i] ) ++vdiag_high_sigma_grow;
	      if ( high_emg[i] ) ++vdiag_high_emg;
	      if ( high_h3[i] ) ++vdiag_high_h3;
	      if ( high_pwr[i] ) ++vdiag_high_pwr;
	      if ( art_emg[i] ) ++vdiag_art_emg;
	      if ( art_h3[i] ) ++vdiag_art_h3;
	      if ( art_pwr[i] ) ++vdiag_art_pwr;
	      if ( artifact[i] ) ++vdiag_artifact;
	      if ( artifact_block[i] ) ++vdiag_artifact_block;
	    }

	  add_verbose_mask( "high_" + detect_label + "_peak" , high_detect_peak );
	  add_verbose_mask( "high_" + detect_label + "_grow" , high_detect_grow );
	  add_verbose_mask( "high_sigma_peak_veto" , high_sigma_peak );
	  add_verbose_mask( "high_sigma_grow_veto" , high_sigma_grow );
	  add_verbose_mask( "high_emg" , high_emg );
	  add_verbose_mask( "high_h3" , high_h3 );
	  add_verbose_mask( "high_pwr" , high_pwr );
	  add_verbose_mask( "art_emg_lowshift" , art_emg );
	  add_verbose_mask( "art_h3" , art_h3 );
	  add_verbose_mask( "art_pwr_lowshift" , art_pwr );
	  add_verbose_mask( "artifact" , artifact );

	  //
	  // track means by art/not
	  //

	  for (int i=0; i<ne; i++)
	    {
	      if ( artifact[i] )
		{
		  ftr_art += D[i];
		  ++n_art;
		}
	      else
		{
		  ftr_nonart += D[i];
		  ++n_nonart;
		}
	    }

	  // get peaks in the selected detector feature (shift by default; EMG in
	  // diagnostic EMG-only mode).  Allow high sigma to veto shift-based NREM
	  // candidates only.

	  std::vector<int> pks;
	  for (int i=0; i<ne; i++)
	    {
	      const Eigen::VectorXd & ftr = D[i];

	      if ( ftr(idx_detect) < p.shift_peak ) continue;

	      if ( i != 0 )
		if ( ftr(idx_detect) <= D[i-1](idx_detect) )
		  continue;

	      if ( i != ne-1 )
		if ( ftr(idx_detect) < D[i+1](idx_detect) )
		  continue;

	      ++vdiag_local_shift_peak;

	      if ( ! p.emg_mode && artifact_block[i] )
		{
		  ++vdiag_peak_block_artifact;
		  continue;
		}

	      if ( ! p.emg_mode && ftr(idx_sigma) > th_sigma_veto )
		{
		  ++vdiag_peak_block_sigma;
		  continue;
		}

	      ++vdiag_peak_accepted;
	      pks.push_back(i);

	    }

	  if ( arousal_verbose_annots )
	    {
	      std::vector<bool> peak_mask( ne , false );
	      for (int i=0; i<pks.size(); i++) peak_mask[ pks[i] ] = true;
	      add_verbose_mask( detect_label + "_local_peak" , peak_mask );
	    }

	  // peaks -> events
	  std::vector<std::pair<int,int> > evts;
	  for (int pki=0; pki<pks.size(); pki++)
	    {
	      // walk back
	      int start = pks[pki];
	      while ( 1 ) {
		if ( start == 0 )
		  break;
		if ( ! p.emg_mode && artifact_block[start-1] )
		  break;
		if ( D[start-1](idx_detect) < p.shift_hysteresis )
		  break;
		if ( ! p.emg_mode && D[start-1](idx_sigma) > th_sigma_veto2 )
		  break;
		--start;
	      }

	      // walk forward
	      int stop = pks[pki];
	      while ( 1 ) {
		if ( stop == ne-1 )
		  break;
		if ( ! p.emg_mode && artifact_block[stop+1] )
		  break;
		if ( D[stop+1](idx_detect) < p.shift_hysteresis )
		  break;
		if ( ! p.emg_mode && D[stop+1](idx_sigma) > th_sigma_veto2 )
		  break;
		++stop;
	      }

	      evts.push_back( std::pair<int,int>(start,stop) );

	    }

	  add_verbose_events( "candidate_raw" , evts );
	  vcand_raw += evts.size();

	  // prune/merge events (duration on the actual epoch grid, not
	  // assuming a fixed 2 Hz grid)

	  std::vector<std::pair<int,int> > evts2;
	  for (int e=0; e<evts.size(); e++)
	    {
	      const std::pair<int,int> & evt = evts[e];
	      const double duration_sec = (evt.second - evt.first + 1) * p.epoch_inc;
	      if ( duration_sec >= p.min_dur && duration_sec <= candidate_max_dur ) evts2.push_back( evt );
	    }

	  add_verbose_events( "candidate_duration_ok" , evts2 );
	  vcand_dur_ok += evts2.size();

	  // merge shorter events that are near
	  evts = merge_events_with_gap_sorted( evts2 , max_gap_samples );

	  add_verbose_events( "candidate_merged" , evts );
	  vcand_merged += evts.size();


	  // final pruning
	  evts2.clear();
	  for (int e=0; e<evts.size(); e++)
	    {
	      const std::pair<int,int> & evt = evts[e];
	      const double duration_sec = (evt.second - evt.first + 1) * p.epoch_inc;
	      if ( duration_sec >= p.min_dur && duration_sec <= candidate_max_dur ) evts2.push_back( evt );
	    }

	  add_verbose_events( "candidate_final_duration_ok" , evts2 );
	  vcand_final_dur_ok += evts2.size();

	  // require at least 'pre-sleep' of stable sleep (i.e. base just on the contig)
	  evts.clear();
	  for (int e=0; e<evts2.size(); e++)
	    {
	      const std::pair<int,int> & evt = evts2[e];
	      int start_idx = evt.first;
	      double t_event_start = tt[st][c][start_idx];
	      double t_contig_start = tt[st][c][0];
	      double pre_sleep_sec = t_event_start - t_contig_start;
	      if (pre_sleep_sec >= p.pre_sleep_sec)
		evts.push_back(evt);
	    }

	  add_verbose_events( "candidate_presleep_ok" , evts );
	  vcand_presleep_ok += evts.size();

		  //
		  // REM only: require a concurrent chin-EMG amplitude rise (AASM) to
		  // confirm each cortical-shift candidate; candidates that fail this
		  // are dropped (not scored) but kept as a diagnostic annotation
		  //
		  if ( is_rem && ! p.emg_mode )
		    {
	      std::vector<bool> emg_high( ne , false );
	      for (int i=0; i<ne; i++)
		emg_high[i] = D[i](idx_emg) > p.emg_rise_th;

	      const std::vector<std::pair<int,int> > emg_runs = mask_to_intervals( emg_high );

	      add_verbose_mask( "rem_emg_rise" , emg_high );

	      std::vector<std::pair<double,double> > emg_windows;
	      for (int r=0; r<emg_runs.size(); r++)
		{
		  const double e0 = tt[st][c][ emg_runs[r].first ];
		  const double e1 = tt[st][c][ emg_runs[r].second ] + p.epoch_inc;
		  if ( e1 - e0 >= p.emg_rise_min_dur )
		    {
		      emg_windows.push_back( std::make_pair( e0 , e1 ) );
		      add_verbose_interval( "rem_emg_rise_duration_ok" , e0 , e1 );
		    }
		}

	      std::vector<std::pair<int,int> > evts_confirmed;
	      std::vector<std::pair<int,int> > evts_rejected;
	      for (int e=0; e<evts.size(); e++)
		{
		  const double t0 = tt[st][c][ evts[e].first ];
		  const double t1 = tt[st][c][ evts[e].second ] + p.epoch_inc;

		  bool confirmed = false;
		  for (int w=0; w<emg_windows.size(); w++)
		    if ( emg_windows[w].first <= t1 + p.emg_rise_buffer &&
			 emg_windows[w].second >= t0 - p.emg_rise_buffer )
		      { confirmed = true; break; }

		  if ( confirmed )
		    {
		      evts_confirmed.push_back( evts[e] );
		      const double dur = t1 - t0;
		      if ( do_long && dur > p.max_dur && dur <= p.long_dur )
			++cnt_rem_emg_confirmed_long;
		    }
		  else
		    {
		      evts_rejected.push_back( evts[e] );
		      const double dur = t1 - t0;
		      if ( do_long && dur > p.max_dur && dur <= p.long_dur )
			{
			  ret[ "long_rem_cortical_only" ].insert( interval_t( globals::tp_1sec * t0, globals::tp_1sec * t1 ) );
			  ++cnt_long_rem_cortical_only;
			}
		      else
			{
			  ret[ "rem_cortical_only" ].insert( interval_t( globals::tp_1sec * t0, globals::tp_1sec * t1 ) );
			  ++cnt_rem_cortical_only;
			}
		    }
		}
	      add_verbose_events( "rem_emg_confirmed" , evts_confirmed );
	      add_verbose_events( "rem_emg_rejected" , evts_rejected );
	      vcand_rem_emg_confirmed += evts_confirmed.size();
	      vcand_rem_emg_rejected += evts_rejected.size();
	      evts = evts_confirmed;
	    }

	  //
	  // convert to intervals
	  //
	  std::vector<std::pair<int,int>> arr_major, arr_micro;

	  for (int e=0; e<evts.size(); e++)
	    {
	      const double t0 = tt[st][c][evts[e].first];
	      const double t1 = tt[st][c][evts[e].second] + p.epoch_inc; // up to end

	      // annotate arousal (>= arousal-dur) versus micro-arousal
	      const double dur = t1 - t0;

	      if ( do_long && dur > p.max_dur && dur <= p.long_dur )
		{
		  ret[ "long_arousal_" + stg_lab ].insert( interval_t( globals::tp_1sec * t0, globals::tp_1sec * t1 ) );
		  dur_long += dur;
		  cnt_long_evts++;
		}
	      else if ( dur >= p.arousal_dur && dur <= p.max_dur )
		{
		  ret[ "arousal_" + stg_lab ].insert( interval_t( globals::tp_1sec * t0, globals::tp_1sec * t1 ) );
		  arr_major.push_back( evts[e] );
		  dur_major += dur;
		  cnt_evts++;
		}
	      else
		{
		  ret[ "micro_arousal_" + stg_lab ].insert( interval_t( globals::tp_1sec * t0, globals::tp_1sec * t1 ) );
		  arr_micro.push_back( evts[e] );
		  dur_micro += dur;
		  cnt_uevts++;
		}
	    }

	  // also track artifacts

	  std::vector<std::pair<int,int>> arts = mask_to_intervals( artifact );
	  for (int e=0; e<arts.size(); e++)
	    {
	      const double t0 = tt[st][c][arts[e].first];
	      const double t1 = tt[st][c][arts[e].second] + p.epoch_inc; // up to end
	      ret[ "art_" + stg_lab ].insert( interval_t( globals::tp_1sec * t0, globals::tp_1sec * t1 ) );
	    }

	  cnt_arts += arts.size();

	  //
	  // also track ftr means by non-art baseline vs non-art arousal
	  //

	  std::vector<int> arr( ne , 0 ); // 0 / 1 / 2 = baseline / micro-arousal / arousal
	  for (int e=0; e<arr_micro.size(); e++)
	    for (int i=arr_micro[e].first; i<=arr_micro[e].second; i++)
	      arr[i] = 1;
	  for (int e=0; e<arr_major.size(); e++)
	    for (int i=arr_major[e].first; i<=arr_major[e].second; i++)
	      arr[i] = 2;

	  for (int i=0; i<ne; i++)
	    {
	      if ( ! artifact[i] )
		{
		  if ( arr[i] == 2 )
		    {
		      ftr_arousal += D[i];
		      ++n_arousal;
		    }
		  else if ( arr[i] == 1 )
		    {
		      ftr_uarousal += D[i];
		      ++n_uarousal;
		    }
		  else
		    {
		      ftr_baseline += D[i];
		      ++n_baseline;
		    }
		}
	    }

	  // next contig
	}

      //
      // per-stage output
      //

      writer.level( is_rem ? "R" : "NR" , "SS" );

      writer.level( "artifact" , "CLS" );
      writer.value( "NE", n_art );
      if ( n_art > 0 )
	{
	  writer.value( "PWR", ftr_art(0) / (double)n_art );
	  writer.value( "SHIFT", ftr_art(1) / (double)n_art );
	  writer.value( "EMG", ftr_art(2) / (double)n_art );
	  writer.value( "SIGMA", ftr_art(3) / (double)n_art );
	  writer.value( "CMPLX", ftr_art(4) / (double)n_art );
	}

      writer.level( "non_artifact" , "CLS" );
      writer.value( "NE", n_nonart );
      if ( n_nonart > 0 )
	{
	  writer.value( "PWR", ftr_nonart(0) / (double)n_nonart );
	  writer.value( "SHIFT", ftr_nonart(1) / (double)n_nonart );
	  writer.value( "EMG", ftr_nonart(2) / (double)n_nonart );
	  writer.value( "SIGMA", ftr_nonart(3) / (double)n_nonart );
	  writer.value( "CMPLX", ftr_nonart(4) / (double)n_nonart );
	}

      writer.level( "arousal" , "CLS" );
      writer.value( "NE", n_arousal );
      if ( n_arousal > 0 )
	{
	  writer.value( "PWR", ftr_arousal(0) / (double)n_arousal );
	  writer.value( "SHIFT", ftr_arousal(1) / (double)n_arousal );
	  writer.value( "EMG", ftr_arousal(2) / (double)n_arousal );
	  writer.value( "SIGMA", ftr_arousal(3) / (double)n_arousal );
	  writer.value( "CMPLX", ftr_arousal(4) / (double)n_arousal );
	}

      writer.level( "micro_arousal" , "CLS" );
      writer.value( "NE", n_uarousal );
      if ( n_uarousal > 0 )
	{
	  writer.value( "PWR", ftr_uarousal(0) / (double)n_uarousal );
	  writer.value( "SHIFT", ftr_uarousal(1) / (double)n_uarousal );
	  writer.value( "EMG", ftr_uarousal(2) / (double)n_uarousal );
	  writer.value( "SIGMA", ftr_uarousal(3) / (double)n_uarousal );
	  writer.value( "CMPLX", ftr_uarousal(4) / (double)n_uarousal );
	}

      writer.level( "baseline" , "CLS" );
      writer.value( "NE", n_baseline );
      if ( n_baseline > 0 )
	{
	  writer.value( "PWR", ftr_baseline(0) / (double)n_baseline );
	  writer.value( "SHIFT", ftr_baseline(1) / (double)n_baseline );
	  writer.value( "EMG", ftr_baseline(2) / (double)n_baseline );
	  writer.value( "SIGMA", ftr_baseline(3) / (double)n_baseline );
	  writer.value( "CMPLX", ftr_baseline(4) / (double)n_baseline );
	}

      writer.unlevel( "CLS" );


      // arousal rate
      double tot_sec = 0;
      const int nc2 = X[st].size();
      for (int c=0; c<nc2; c++)
	{
	  const int ne = tt[st][c].size();
	  if ( ne == 0 ) continue;
	  double mint = tt[st][c][0];
	  double maxt = tt[st][c][ne-1];
	  tot_sec += maxt - mint + p.epoch_inc;
	}

      // total time
      writer.value( "MINS" , tot_sec / 60.0 );

      // report arousal counts
      writer.value( "N" , cnt_evts );
      writer.value( "AI" , tot_sec > 0 ? cnt_evts / ( tot_sec / 3600.0 ) : 0.0 );
      if ( cnt_evts > 0 ) writer.value( "DUR" , dur_major / (double)cnt_evts );

      // report micro-arousal counts
      writer.value( "N_MICRO" , cnt_uevts );
      writer.value( "AI_MICRO" , tot_sec > 0 ? cnt_uevts / ( tot_sec / 3600.0 ) : 0.0 );
      if ( cnt_uevts > 0 ) writer.value( "DUR_MICRO" , dur_micro / (double)cnt_uevts );

      // report long arousal-like events separately so the original arousal
      // set, index and feature summaries remain backward-compatible.
      writer.value( "N_LONG" , cnt_long_evts );
      writer.value( "AI_LONG" , tot_sec > 0 ? cnt_long_evts / ( tot_sec / 3600.0 ) : 0.0 );
      if ( cnt_long_evts > 0 ) writer.value( "DUR_LONG" , dur_long / (double)cnt_long_evts );

      // artifacts
      writer.value( "N_ART" , cnt_arts );
      writer.value( "AI_ART" , tot_sec > 0 ? cnt_arts / ( tot_sec / 3600.0 ) : 0.0 );

      if ( p.verbose )
	{
	  auto vdiag_win_prop = [&]( const int n ) {
	    return vdiag_win > 0 ? n / (double)vdiag_win : 0.0;
	  };
	  auto vdiag_peak_prop = [&]( const int n ) {
	    return vdiag_local_shift_peak > 0 ? n / (double)vdiag_local_shift_peak : 0.0;
	  };
	  auto vcand_prop = [&]( const int n ) {
	    return vcand_raw > 0 ? n / (double)vcand_raw : 0.0;
	  };

	  writer.value( "VDIAG_N_WIN" , vdiag_win );
	  writer.value( "VDIAG_N_HIGH_SHIFT_PEAK" , vdiag_high_shift_peak );
	  writer.value( "VDIAG_PROP_HIGH_SHIFT_PEAK" , vdiag_win_prop( vdiag_high_shift_peak ) );
	  writer.value( "VDIAG_N_HIGH_SHIFT_GROW" , vdiag_high_shift_grow );
	  writer.value( "VDIAG_PROP_HIGH_SHIFT_GROW" , vdiag_win_prop( vdiag_high_shift_grow ) );
	  writer.value( "VDIAG_N_HIGH_SIGMA_PEAK" , vdiag_high_sigma_peak );
	  writer.value( "VDIAG_PROP_HIGH_SIGMA_PEAK" , vdiag_win_prop( vdiag_high_sigma_peak ) );
	  writer.value( "VDIAG_N_HIGH_SIGMA_GROW" , vdiag_high_sigma_grow );
	  writer.value( "VDIAG_PROP_HIGH_SIGMA_GROW" , vdiag_win_prop( vdiag_high_sigma_grow ) );
	  writer.value( "VDIAG_N_HIGH_EMG" , vdiag_high_emg );
	  writer.value( "VDIAG_PROP_HIGH_EMG" , vdiag_win_prop( vdiag_high_emg ) );
	  writer.value( "VDIAG_N_HIGH_H3" , vdiag_high_h3 );
	  writer.value( "VDIAG_PROP_HIGH_H3" , vdiag_win_prop( vdiag_high_h3 ) );
	  writer.value( "VDIAG_N_HIGH_PWR" , vdiag_high_pwr );
	  writer.value( "VDIAG_PROP_HIGH_PWR" , vdiag_win_prop( vdiag_high_pwr ) );
	  writer.value( "VDIAG_N_ART_EMG_LOWSHIFT" , vdiag_art_emg );
	  writer.value( "VDIAG_PROP_ART_EMG_LOWSHIFT" , vdiag_win_prop( vdiag_art_emg ) );
	  writer.value( "VDIAG_N_ART_H3" , vdiag_art_h3 );
	  writer.value( "VDIAG_PROP_ART_H3" , vdiag_win_prop( vdiag_art_h3 ) );
	  writer.value( "VDIAG_N_ART_PWR_LOWSHIFT" , vdiag_art_pwr );
	  writer.value( "VDIAG_PROP_ART_PWR_LOWSHIFT" , vdiag_win_prop( vdiag_art_pwr ) );
	  writer.value( "VDIAG_N_ARTIFACT" , vdiag_artifact );
	  writer.value( "VDIAG_PROP_ARTIFACT" , vdiag_win_prop( vdiag_artifact ) );
	  writer.value( "VDIAG_N_ARTIFACT_BLOCK" , vdiag_artifact_block );
	  writer.value( "VDIAG_PROP_ARTIFACT_BLOCK" , vdiag_win_prop( vdiag_artifact_block ) );
	  writer.value( "VDIAG_N_LOCAL_SHIFT_PEAK" , vdiag_local_shift_peak );
	  writer.value( "VDIAG_N_PEAK_BLOCK_ARTIFACT" , vdiag_peak_block_artifact );
	  writer.value( "VDIAG_PROP_PEAK_BLOCK_ARTIFACT" , vdiag_peak_prop( vdiag_peak_block_artifact ) );
	  writer.value( "VDIAG_N_PEAK_BLOCK_SIGMA" , vdiag_peak_block_sigma );
	  writer.value( "VDIAG_PROP_PEAK_BLOCK_SIGMA" , vdiag_peak_prop( vdiag_peak_block_sigma ) );
	  writer.value( "VDIAG_N_PEAK_ACCEPTED" , vdiag_peak_accepted );
	  writer.value( "VDIAG_PROP_PEAK_ACCEPTED" , vdiag_peak_prop( vdiag_peak_accepted ) );

	  writer.value( "VCAND_N_RAW" , vcand_raw );
	  writer.value( "VCAND_N_DUR_OK" , vcand_dur_ok );
	  writer.value( "VCAND_PROP_DUR_OK" , vcand_prop( vcand_dur_ok ) );
	  writer.value( "VCAND_N_MERGED" , vcand_merged );
	  writer.value( "VCAND_PROP_MERGED" , vcand_prop( vcand_merged ) );
	  writer.value( "VCAND_N_FINAL_DUR_OK" , vcand_final_dur_ok );
	  writer.value( "VCAND_PROP_FINAL_DUR_OK" , vcand_prop( vcand_final_dur_ok ) );
	  writer.value( "VCAND_N_PRESLEEP_OK" , vcand_presleep_ok );
	  writer.value( "VCAND_PROP_PRESLEEP_OK" , vcand_prop( vcand_presleep_ok ) );

	  if ( is_rem )
	    {
	      writer.value( "VCAND_N_REM_EMG_CONFIRMED" , vcand_rem_emg_confirmed );
	      writer.value( "VCAND_PROP_REM_EMG_CONFIRMED" , vcand_prop( vcand_rem_emg_confirmed ) );
	      writer.value( "VCAND_N_REM_EMG_REJECTED" , vcand_rem_emg_rejected );
	      writer.value( "VCAND_PROP_REM_EMG_REJECTED" , vcand_prop( vcand_rem_emg_rejected ) );
	    }
	}

	      // REM-only: EMG-confirmation diagnostics
	      if ( is_rem )
		{
		  writer.value( "N_REM_CORTICAL_ONLY" , cnt_rem_cortical_only );
		  const int denom = vcand_rem_emg_confirmed + vcand_rem_emg_rejected;
		  if ( denom > 0 ) writer.value( "PCT_REM_EMG_CONFIRMED" , 100.0 * vcand_rem_emg_confirmed / (double)denom );
		  writer.value( "N_LONG_REM_CORTICAL_ONLY" , cnt_long_rem_cortical_only );
		  const int denom_long = cnt_rem_emg_confirmed_long + cnt_long_rem_cortical_only;
		  if ( denom_long > 0 ) writer.value( "PCT_LONG_REM_EMG_CONFIRMED" , 100.0 * cnt_rem_emg_confirmed_long / (double)denom_long );
		}

      writer.unlevel( "SS" );

      // next stage
    }

  return ret;
}


std::vector<std::pair<int,int>>
arousals_t::merge_events_with_gap_sorted(const std::vector<std::pair<int,int>> &events,
					 int max_gap)
{
  if (events.empty()) return {};
  
  std::vector<std::pair<int,int>> merged;
  int cur_start = events[0].first;
  int cur_end   = events[0].second;
  
  for (size_t i = 1; i < events.size(); ++i) {
    int s = events[i].first;
    int e = events[i].second;
    
    // closed intervals; merge if overlap or gap <= max_gap
    if (s <= cur_end + max_gap) {
      if (e > cur_end) cur_end = e;
    } else {
      merged.emplace_back(cur_start, cur_end);
      cur_start = s;
      cur_end   = e;
    }
  }
  
  merged.emplace_back(cur_start, cur_end);
  return merged;
}


#include <vector>
#include <utility>

// mask[i] == true => inside the interval
// events are closed intervals [start, stop]
std::vector<std::pair<int,int>>
arousals_t::mask_to_intervals(const std::vector<bool> &mask)
{
  std::vector<std::pair<int,int>> out;
  const int n = mask.size();
  if (n == 0) return out;
  
  int start = -1;
  
  for (int i = 0; i < n; ++i) {
    if (mask[i]) {
      if (start == -1)
	start = i;     // begin a new run
    }
    else {
      if (start != -1) {
	// end previous run just before i
	out.emplace_back(start, i - 1);
	start = -1;
      }
    }
  }
  
  // end run at end of vector
  if (start != -1)
    out.emplace_back(start, n - 1);
  
  return out;
}

#endif // 0: retained as the historical legacy AROUSALS implementation
