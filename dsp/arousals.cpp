
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

extern logger_t logger;
extern writer_t writer;

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
  //  epoch size = 1 second, 0.5 second shift
  //  will have 50% shift, i.e. 2 Hz signal effectively

  const double epoch_win = param.has( "win" ) ? param.requires_dbl( "win" ) : 1.0 ;
  const double epoch_inc = param.has( "inc" ) ? param.requires_dbl( "inc" ) : 0.5 ;
  if ( epoch_win < 1 || epoch_win < epoch_inc ) Helper::halt( "invalid epoch win/inc values" );
  
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
  
  const std::string new_sig_prefix = param.has( "prefix" ) ? param.value( "prefix" ) : "" ; 
  const bool write_channel_specific_metrics = param.has( "per-channel" ) ? param.yesno( "per-channel" ) : false;
  
  // handle outliers
  double winsor_th = param.has( "no-winsor" ) ? -9 : 0.005;
  if ( param.has( "winsor" ) ) winsor_th = param.requires_dbl( "winsor" ) ;
  if ( winsor_th > 0.2 ) Helper::halt( "winsor should be less than 0.2" );
  
  // add annotations (instance ID is 'REM' or 'NREM')
  const std::string aname = param.has( "annot" ) ? param.value( "annot" ) : "l" ;    

  // add channels
  const bool add_chs = param.has( "add" );
  const std::string ch_prefix = ! param.empty( "add" ) ? param.value( "add" ) : "a_" ; 
  
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

  if ( ns_eeg == 0 )
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
  
  sr = Fs_eeg[0];
  
  for (int i=0; i<ns_eeg; i++)
    {
      if ( fabs( Fs_eeg[i] - sr ) > 1e-4 )
	Helper::halt( "all EEG signals must have similar sample rates" );      
    }
  
  for (int i=1; i<ns_emg; i++)
    {
      if ( fabs( Fs_emg[i] - Fs_emg[0] ) > 1e-4 )
	Helper::halt( "all EMG signals must have similar sample rates" );      
    }

  if ( sr < 60 )
    Helper::halt( "EEG sample rate too low" );
  
  //
  // Construct initial feature matrix, w/ time-track and NREM/REM annotation
  //
  
  // state (based on *start* of epoch)
  std::vector<int> state; // 0 = NR , 1 = R  , 2 = W
  std::vector<int> seq; // contig w/in each stage/state
  std::vector<double> sec; // time of each epoch start
  
  // features:
  //    log-pwr / rel-beta / emg-rms / [rel-sigma] / [h3 ]
  //  only first 3 ftrs used for hmm/detector

  // low-pwr --> artifact: broadband (0.5-35Hz) “movement/artifact” axis (optionally: 30-45 Hz HF-EEG)
  // beta    --> EEG arousal
  // emg     --> EMG bursts
  // sigma (subtyping events)
  // h3    (subtyping)

  // build features
  Eigen::MatrixXd Xeeg, Xemg;
  build_ftr_matrix( edf , eeg_signals, emg_signals , Xeeg, Xemg, &state, &seq, &sec ) ;
  
  // process [ epoch x (3+2) ftrs ] 
  Eigen::MatrixXd Xftr = process_ftr_matrix( &Xeeg , &Xemg , state );

  // tidy up
  Xeeg.resize(0,0); Xemg.resize(0,0);
  
  // time-track
  std::vector<std::vector<std::vector<double> > > tt;

  // assemble into sleep-state-specific contigs
  std::vector<std::vector<std::vector<Eigen::VectorXd> > > X = assemble( Xftr, state, seq, sec , &tt );

  // extract key features for HMM detection
  //         0        1      2       3      4    
  //         totpow , beta   emg     sigma  h3
  // first three ftrs (of 5)
  // std::vector<int> ex = { 1,2 };
  // X = extract( X , ex );

  
  // optionally,  dump to stdout
  // dump( X , tt );

  std::map<std::string,std::set<interval_t> > anns = event_heuristic( X , tt );
  
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
  // add channels0
  //
  
  // NREM
  if ( add_chs )
    add_channels( X[0], tt[0] , ch_prefix ); 


  // all done (skip HMM part)

  return;

  // --------------------------------------------------------------------------------


  // not used
  // hmm( # states , # ftrs ) 
  gaussian_hmm_t nrem_hmm( 3, 2 );
  init_kmeans_hmm( nrem_hmm , X[0] );
  const int max_iters = 30;
  const double tol = 1e-4;
  nrem_hmm.train_multi( X[0], max_iters, tol);
  // output
  const int nrem_seq = X[0].size();
  std::vector<std::vector<int> > nrem_paths( nrem_seq );
  std::vector<double> nrem_loglik( nrem_seq );
  std::vector<Eigen::MatrixXd> nrem_posteriors( nrem_seq );
  
  for (int i=0; i<nrem_seq; i++)
    {
      nrem_hmm.viterbi( X[0][i], nrem_paths[i]);
      Eigen::MatrixXd gamma;
      double ll_seq = 0.0;      
      nrem_hmm.posteriors( X[0][i], gamma , &nrem_loglik[i] );     
      for (int p=0;p<nrem_paths[i].size();p++)
	{
	   std::cout << i << "\t"
	   	    << nrem_loglik[i] << "\t"
	   	    << p << "\t"
	   	    << nrem_paths[i][p] << "\t"
	   	    << gamma(p,0) << "\t"
	   	    << gamma(p,1) << "\t"
	   	    << gamma(p,2) << "\n";
	}
      
      
      nrem_posteriors[i] = gamma;
    }
  
  for (int k = 0; k < 3; ++k) {
    std::cerr << "state " << k << " mu = " << nrem_hmm.mu().col(k).transpose() << "\n";
    std::cerr << "state " << k << " cov:\n" << nrem_hmm.covariances()[k] << "\n";
  }

   
 
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

  std::vector<std::string> chlab = { "totpwr", "beta", "emg", "sigma", "h3" };
    
  for (int cidx = 0; cidx<5; cidx++)
    {
      
      // make a new signal (2Hz)
      std::string label = ch_prefix + "_" + chlab[cidx] ;
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
	sec[i] = (*tp)[i] / globals::tp_1sec;
      
      // new data
      std::vector<double> xx( np , 0 );
      
      // get original times, and values      
      for (int sq=0; sq<tt.size(); sq++)
	{
	  // within this sequence (contig), interpolate
	  const int ne = tt[sq].size();
	  std::vector<double> xx0(ne);      
	  for (int i=0; i<ne; i++)
	    xx0[i] = X[sq][i][cidx];
	  
	  // use spline to get values
	  tk::spline spline;
	  spline.set_points( tt[sq] , xx0 );
	  
	  const double tmin = tt[sq][0];
	  const double tmax = tt[sq][ tt[sq].size() - 1 ];
	  
	  // spline interpolation at this resolution
	  std::vector<double> interp;
	  for (int i = 0 ; i < np; i++ )
	    {
	      if ( sec[i] > tmax ) break;
	      if ( sec[i] < tmin ) continue;
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
    for (int sq=0; sq<contigs[st].size(); sq++)
      {	
	const std::vector<int> & ee = contigs[st][sq];
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
				   std::vector<double> * sec 
				   )
					      
{

  // state: annotate whether start of each window is NR, R or W (0,1,2)
  // seq: track contiguous, within-state sequences i.e. (0,0,0,1,1,1,2,2,2,etc)

  // three features
  //   artifact          | Neeg x 1 
  //   EEG arousal       | core      = Neeg x 3 [ beta / (alpha, theta) ]
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
  
  // ensure staging is present in SleepStage  
  edf.annotations->make_sleep_stage( edf.timeline );

  annot_t * staging = edf.annotations->find( "SleepStage" );

  if ( staging == NULL )
    Helper::halt( "no staging information present" );
  
  // track state-specific contigs
  int prior_state = 2;  // 0/1 = NR/R
  int seq_nr = -1;
  int seq_r = -1;
  
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
      
      int st = 2; // wake
      if ( events.size() >= 1 )
	{
	  // get first annotaiton, i.e. for start, if >1
	  const instance_idx_t & idx = events.begin()->first;
	  if ( idx.id == "N1" || idx.id == "N2" || idx.id == "N3" || idx.id == "NREM4" || idx.id == "NR" )
	    st = 0;
	  else if ( idx.id == "R" )
	    st = 1; 
	}	  

      // a new NR or R sequence?
      if ( st != prior_state ) 
	{
	  if ( st == 0 ) ++seq_nr;	    
	  else if ( st == 1 ) ++seq_r;
	}

      // for next time round
      prior_state = st;
      
      // track
      state->push_back( st );      
      if ( st == 0 ) sequence->push_back( seq_nr );
      else if ( st == 1 ) sequence->push_back( seq_r );
      else sequence->push_back( -1 ); // we don't care: W or ?
      sec->push_back( interval.start_sec() );
      
      // if non-sleep, can skip
      if ( st == 2 ) continue;
      
      // get EEGs
      eigen_matslice_t mslice( edf , eeg_signals , interval );
      const Eigen::MatrixXd & X = mslice.data_ref();

      // derive features
      Xeeg.row( eidx ) = calc_eeg_ftrs( X );

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

Eigen::VectorXd arousals_t::calc_eeg_ftrs( const Eigen::MatrixXd & X  )
{

  // currently, get 4 ftrs from each EEG: 2 for detector, 2 for subtyping
  Eigen::VectorXd F = Eigen::VectorXd::Zero( 4 * X.cols() );

  // FFT
  int index_length = X.rows();

  FFT fftseg( index_length , index_length , sr , FFT_FORWARD , WINDOW_TUKEY50 );
  
  // per channel

  const double eps = 1e-12;

  int fi = 0;
  
  for (int s=0; s<X.cols(); s++)
    {
      fftseg.apply( X.col(s).data() , index_length );

      double p_tot = 0;   // 4-30 
      double p_theta = 0; // 4-8
      double p_alpha = 0; // 8-10
      double p_sigma = 0; // 10-16
      double p_beta = 0;  // 16-30

      double hjorth3 = 0;      
      
      for (int f = 0; f < fftseg.cutoff; f++)
	{
	  double frq = fftseg.frq[f] ;	  
	  if ( frq < 4 || frq >= 30 ) continue;	  
	  p_tot += fftseg.X[f] ;
	  if      ( frq >= 16 ) p_beta  += fftseg.X[f];
	  else if ( frq >= 10 ) p_sigma += fftseg.X[f];
	  else if ( frq >= 8  ) p_alpha += fftseg.X[f];
	  else if ( frq >= 4  ) p_theta += fftseg.X[f];	  
	}
      
      // transform to log/relative power
      
      p_beta = log( p_beta / ( p_alpha + p_theta + eps ) );
      p_sigma = log( p_sigma / ( p_tot + eps ) );
      p_tot = log( p_tot + eps );

      // Hjorth
      double h1, h2, h3;
      bool okay = hjorth( X.col(s) , &h1, &h2, &h3 );
      if ( ! okay ) h3 = 0;
      
      F(fi++) = p_tot;
      F(fi++) =	p_beta;
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

  // for now assume single EEG, single EMG so  alpha beta EMG
  // --> otherwise, make collapsed (mean) matrix
  
  // mu - rows    = metrics ; 
  //      cols    = states:

  //  std::cout << "MMMuuu\n" << mu << "\n";

  Eigen::VectorXd alpha = mu.row(0);
  Eigen::VectorXd beta = mu.row(1);
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
  //  beta not very high (beta is not extremely high, i.e. less than 0.3 from max beta(
  int max_beta;
  beta.maxCoeff( &max_beta );
  double delta_beta = beta(max_beta) - beta(max_emg);

  bool has_artifact_class = true;
  if ( max_emg == min_emg ) has_artifact_class = false;
  if ( emg(max_emg) < 0.5 ) has_artifact_class = false;
  if ( delta_emg < 0.5 ) has_artifact_class = false;
  if ( delta_beta < 0.3 ) has_artifact_class = false;

  // if no reliable artifact: treat as two-class model

  // from non-artifact classes:
  //   heuristic: favor beta, with some contribution from alpha and EMG

  // Sar​(k)=bk​+0.5ak​+0.3mk

  double s0 = beta(0) + 0.5 * alpha(0) + 0.3 * emg(0);
  double s1 = beta(1) + 0.5 * alpha(1) + 0.3 * emg(1);
  double s2 = beta(2) + 0.5 * alpha(2) + 0.3 * emg(2);
  
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
  double delta_beta2 = beta(*st_arousal) - beta(*st_baseline);
  double delta_emg2 = emg(*st_arousal) - emg(*st_baseline);

  // NREM: cortical arousal only, do not require EMG
  // REM: requires both EEG and EMG arousals, as per AASM
  bool has_arousal_state = nrem 
    ? delta_beta2 > 0.4 
    : delta_beta2 > 0.4 && delta_emg2 > 0.5;

  // std::cout << "beta arousal criteria: " << beta(*st_arousal) << " - " << beta(*st_baseline) << " = " << delta_beta2 << "\n";
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



Eigen::MatrixXd arousals_t::process_ftr_matrix( Eigen::MatrixXd * Xeeg , Eigen::MatrixXd * Xemg , const std::vector<int> & st )
{
  // input:
  //   per-channel eeg: tot-pow, beta, [ sigma, h3 ]
  //   per-channel emg: rms
  //   log-scaling/relative power already done

  // output:
  //  subtract local baseline for each feature
  //  average over channels
  //  robust norm based on sleep epochs
  //  return 5 values in single matrix 

  if ( Xeeg->rows() != Xemg->rows() )
    Helper::halt( "internal error in arousals_t::process_ftr_matrix()" );
  
  const int nch_eeg = Xeeg->cols() / 4;
  const int nch_emg = Xemg->cols();

  const int neeg = Xeeg->cols();

  // median filter (including *all* epochs, wake+sleep)
  // 240 = two mins filter at 2Hz

  for (int s=0; s< Xeeg->cols(); s++)
    Xeeg->col(s) = Xeeg->col(s) - eigen_ops::median_filter( Xeeg->col(s) , 240 );
  
  for (int s=0; s< Xemg->cols(); s++)
    Xemg->col(s) = Xemg->col(s) - eigen_ops::median_filter( Xemg->col(s) , 240 );
  
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
  
  // now -> X = EEG-pwr EEG-beta EMG-rms [ EEG-spindle EEG-h3 ] 
  //  a single, averaged value per channel


  // finally, sleep-specific robust norm
  for (int i=0; i<X.cols(); i++)
    X.col(i) = robust_mad_norm( X.col(i) , st );

  // and clip
  const double zth = 6;
  for (int i=0; i<X.cols(); i++)
    X.col(i) = X.col(i).cwiseMax(-zth).cwiseMin(zth);
  
  // upweight the beta track
  //  X.col(1) = X.col(1).array() ;
  
  // all done, give back
  return X;
}





Eigen::VectorXd arousals_t::robust_mad_norm(const Eigen::VectorXd & x,
                                            const std::vector<int> & st)
{

  const int n = x.size();

  if ((int)st.size() != n)
    Helper::halt( "robust_mad_norm: x and st must have same length" );
  
  // collect sleep-only values (e.g. st != 2)
  std::vector<double> vals;
  vals.reserve(n);
  for (int i = 0; i < n; ++i)
    {
      if (st[i] != 2 && std::isfinite(x[i])) // sleep only, skip NaNs
	vals.push_back(x[i]);
    }
  
  if (vals.empty())
    Helper::halt( "robust_mad_norm: no sleep samples to estimate median/MAD");
  
  const double median = MiscMath::median(vals);
  
  // MAD over sleep-only
  std::vector<double> dev;
  dev.reserve(vals.size());
  for (double v : vals)
    dev.push_back(std::fabs(v - median));

  const double MAD = MiscMath::median(dev);

  double sigma = 1.4826 * MAD;
  if (sigma <= 0.0)
    sigma = 1.0;  // fallback: avoids division by zero
  
  // normalize all samples (wake + sleep) using sleep-based median/sigma
  Eigen::VectorXd z(n);
  for (int i = 0; i < n; ++i)
    z[i] = (x[i] - median) / sigma;
  
  return z;
}


std::map<std::string,std::set<interval_t> > arousals_t::event_heuristic( const std::vector<std::vector<std::vector<Eigen::VectorXd> > > & X ,
									 const std::vector<std::vector<std::vector<double> > > & tt )
{
  // metrics 0 tot-pwr
  //         1 beta
  //         2 EMG
  //         3 sigma
  //         4 H3

  std::map<std::string,std::set<interval_t> > ret;
  
  const int idx_pwr = 0;
  const int idx_beta = 1;
  const int idx_emg = 2;
  const int idx_sig = 3;
  const int idx_h3 = 4;
    
  // Only NREM for now
  int st = 0;

  std::string stg_lab = "nrem";
  if ( st == 1 ) stg_lab = "rem";
  
  // check ftrs inside/outside of artifact regions
  Eigen::VectorXd ftr_art( 5 ), ftr_nonart( 5 );
  Eigen::VectorXd ftr_baseline( 5 ), ftr_arousal( 5 ), ftr_uarousal( 5 );
  int n_art = 0 , n_nonart = 0, n_baseline = 0 , n_arousal = 0, n_uarousal = 0;

  // event counts
  int cnt_evts = 0 , cnt_uevts = 0, cnt_arts = 0;

  // mean event durations
  double dur_major = 0 , dur_micro = 0;
  
  // within each contig
  const int nc = X[st].size();

  for (int c=0; c<nc; c++)
    {
      const std::vector<Eigen::VectorXd> & D = X[st][c];
      const int ne = D.size();

      // flag artifact
      std::vector<bool> artifact( ne , false );
      int cnt1 = 0, cnt2=0, cnt3=0;
      for (int i=0; i<ne; i++)
	{
	  const Eigen::VectorXd & ftr = D[i];

	  bool is_artifact = false;
	  
	  // Very high EMG
	  if ( ftr(idx_emg) > 5 && ftr(idx_beta) < 0.5 )
	    { is_artifact = true; ++cnt1; } 

	  // High H3 (very noisy)
	  if ( ftr(idx_h3) > 4 )
	    { is_artifact = true; ++cnt2; } 

	  // implausible broadband tot-pwr
	  if ( ftr(idx_pwr) > 4 && ftr(idx_beta) < 0.5 )
	    { is_artifact = true; ++cnt3; } 
	  
	  // set
	  artifact[i] = is_artifact;
	  
	  // expand?
	  // if ( is_artifact )
	  //   {
	  //     if ( i != 0 ) artifact[i-1] = true;
	  //     if ( i != ne - 1 ) artifact[i+1] = true;	  
	  //   }

	}

      std::cout << " art1,2,3 = "
		<< cnt1 / (double)ne << " "
		<< cnt2 / (double)ne << " "
		<< cnt3 / (double)ne << "\n";
      
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


      const double th_beta_peak = 1.2;
      const double th_beta_hysteresis = 0.6;
      
      // get peaks: high beta, peak, not artifact
      std::vector<int> pks;
      for (int i=0; i<ne; i++)
	{
	  if ( artifact[i] ) continue;

	  const Eigen::VectorXd & ftr = D[i];
		    
	  if ( ftr(idx_beta) < th_beta_peak ) continue;

	  if ( i != 0 )
	    if ( ftr(idx_beta) <= D[i-1](idx_beta) )
	      continue;

	  if ( i != ne-1 )
	    if ( ftr(idx_beta ) < D[i+1](idx_beta) )
	      continue;

	  pks.push_back(i);
	  
	}

      
      // peaks -> events  
      std::vector<std::pair<int,int> > evts;
      for (int p=0; p<pks.size(); p++)
	{
	  // walk back
	  int start = pks[p];
	  while ( 1 ) {
	    if ( start == 0 )
	      break;	    
	    if ( artifact[start-1] )
	      break;
	    if ( D[start-1](idx_beta) < th_beta_hysteresis )
	      break;
	    --start;
	  }
	  
	  // walk forward
	  int stop = pks[p];
	  while ( 1 ) {
            if ( stop == ne-1 )
              break;
            if ( artifact[stop+1] )
              break;
            if ( D[stop+1](idx_beta) < th_beta_hysteresis )
              break;
            ++stop;
          }

	  evts.push_back( std::pair<int,int>(start,stop) );
	  
	}
      
      // prune/merge events
      // we now have event from start -> stop
      // assuming 2Hz signal, get duration
      // const double duration_sec = (stop - start + 1) * 0.5;
      
      std::vector<std::pair<int,int> > evts2;
      for (int e=0; e<evts.size(); e++)
	{
	  // ignore long events > 15, and v. short events (<2)
	  const std::pair<int,int> & evt = evts[e];
	  const double duration_sec = (evt.second - evt.first + 1) * 0.5;
	  if ( duration_sec >= 2 && duration_sec <= 15 ) evts2.push_back( evt );
	}
      
      // merge shorter events that are near
      const int max_gap = 5;  // < 2.5 s gap at 2 Hz (0.5 s/epoch)
      evts = merge_events_with_gap_sorted( evts2 , max_gap );
      
      
      // final prunig
      evts2.clear();
      for (int e=0; e<evts.size(); e++)
	{	  
          const std::pair<int,int> & evt = evts[e];
          const double duration_sec = (evt.second - evt.first + 1) * 0.5;	
	  if ( duration_sec >= 2 && duration_sec <= 15 ) evts2.push_back( evt );
        }

      // require at least 10s of stable sleep (i.e. base just on the contig)
      evts.clear();
      for (int e=0; e<evts2.size(); e++)
	{
	  const std::pair<int,int> & evt = evts2[e];
	  int start_idx = evt.first;
	  double t_event_start = tt[st][c][start_idx];
	  double t_contig_start = tt[st][c][0];
	  double pre_sleep_sec = t_event_start - t_contig_start;
	  if (pre_sleep_sec >= 10.0)
	    evts.push_back(evt);
	}
      
      //
      // convert to intervals
      //
      std::vector<std::pair<int,int>> arr_major, arr_micro;
      
      for (int e=0; e<evts.size(); e++)
	{
	  const double t0 = tt[st][c][evts[e].first];
	  const double t1 = tt[st][c][evts[e].second] + 0.5; // up to end

	  // annotate arousal (3+) versus micro-arousal (1.5-2) 
	  const double dur = t1 - t0;

	  if ( dur >= 3 )
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
          const double t1 = tt[st][c][arts[e].second] + 0.5; // up to end
          std::cout << " artifact = " << t0 << " " << t1 << " | " << t1 - t0 << "\n";
          ret[ "art_nrem" ].insert( interval_t( globals::tp_1sec * t0, globals::tp_1sec * t1 ) );
        }

      cnt_arts += arts.size();
      
      //
      // also track ftr means by non-art baseline vs non-art arousal
      //

      std::vector<int> arr( ne , 0 ); // 0 / 1 / 2 = baseline / micro-arousal / arousal 
      for (int e=0; e<arr_micro.size(); e++)
	for (int p=arr_micro[e].first; p<=arr_micro[e].second; p++)
	  arr[p] = 1;
      for (int e=0; e<arr_major.size(); e++)
	for (int p=arr_major[e].first; p<=arr_major[e].second; p++)
	  arr[p] = 2;
	  
      for (int i=0; i<ne; i++)
        {
          if ( ! artifact[i] )
            {
	      if ( arr[i] == 2 )
		{
		  ftr_arousal += D[i];
		  ++n_arousal;
		}
	      if ( arr[i] == 1 )
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
  // report beta per art/non-art
  //


  writer.level( "artifact" , "CLS" );
  writer.value( "NE", n_art );
  writer.value( "PWR", ftr_art(0) / (double)n_art );
  writer.value( "BETA", ftr_art(1) / (double)n_art );
  writer.value( "EMG", ftr_art(2) / (double)n_art );
  writer.value( "SIGMA", ftr_art(3) / (double)n_art );
  writer.value( "CMPLX", ftr_art(4) / (double)n_art );

  writer.level( "non_artifact" , "CLS" );
  writer.value( "NE", n_nonart );  
  writer.value( "PWR", ftr_nonart(0) / (double)n_nonart );
  writer.value( "BETA", ftr_nonart(1) / (double)n_nonart );
  writer.value( "EMG", ftr_nonart(2) / (double)n_nonart );
  writer.value( "SIGMA", ftr_nonart(3) / (double)n_nonart );
  writer.value( "CMPLX", ftr_nonart(4) / (double)n_nonart );

  writer.level( "arousal" , "CLS" );
  writer.value( "NE", n_arousal );
  writer.value( "PWR", ftr_arousal(0) / (double)n_arousal );
  writer.value( "BETA", ftr_arousal(1) / (double)n_arousal );
  writer.value( "EMG", ftr_arousal(2) / (double)n_arousal );
  writer.value( "SIGMA", ftr_arousal(3) / (double)n_arousal );
  writer.value( "CMPLX", ftr_arousal(4) / (double)n_arousal );

  writer.level( "micro_arousal" , "CLS" );
  writer.value( "NE", n_uarousal );
  writer.value( "PWR", ftr_uarousal(0) / (double)n_uarousal );
  writer.value( "BETA", ftr_uarousal(1) / (double)n_uarousal );
  writer.value( "EMG", ftr_uarousal(2) / (double)n_uarousal );
  writer.value( "SIGMA", ftr_uarousal(3) / (double)n_uarousal );
  writer.value( "CMPLX", ftr_uarousal(4) / (double)n_uarousal );
  
  writer.level( "baseline" , "CLS" );
  writer.value( "NE", n_baseline );
  writer.value( "PWR", ftr_baseline(0) / (double)n_baseline );
  writer.value( "BETA", ftr_baseline(1) / (double)n_baseline );
  writer.value( "EMG", ftr_baseline(2) / (double)n_baseline );
  writer.value( "SIGMA", ftr_baseline(3) / (double)n_baseline );
  writer.value( "CMPLX", ftr_baseline(4) / (double)n_baseline );

  writer.unlevel( "CLS" );


  // arousal rate
  double tot_sec = 0;
  for (int c=0; c<nc; c++)
    {
      const int ne = tt[st][c].size();
      if ( ne == 0 ) continue;
      double mint = tt[st][c][0];
      double maxt = tt[st][c][ne-1];
      tot_sec += maxt - mint + 0.5;      
    }

  // total time
  writer.value( "MINS" , tot_sec / 60.0 );

  // report arousal counts
  writer.value( "N" , cnt_evts );
  writer.value( "AI" , cnt_evts / ( tot_sec / 3600.0 ) );
  writer.value( "DUR" , dur_major / (double)cnt_evts );
  
  // report arousal counts
  writer.value( "N_MICRO" , cnt_uevts );
  writer.value( "AI_MICRO" , cnt_uevts / ( tot_sec / 3600.0 ) );
  writer.value( "DUR_MICRO" , dur_micro / (double)cnt_uevts );

  // artifacts
  writer.value( "N_ART" , cnt_arts );
  writer.value( "AI_ART" , cnt_arts / ( tot_sec / 3600.0 ) );

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

