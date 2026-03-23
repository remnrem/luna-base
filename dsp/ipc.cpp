
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

#include "dsp/ipc.h"

#include <vector>
#include <map>
#include "edf/edf.h"
#include "edf/slice.h"
#include "param.h"

#include "helper/logger.h"
#include "db/db.h"

extern writer_t writer;
extern logger_t logger;

namespace
{

ipc_stats_t average_ipc_stats( const std::vector<ipc_stats_t> & stats_in )
{
  ipc_stats_t out;
  out.plv = out.mean_ipc = out.mean_ipc_weighted = out.frac_inphase = 0.0;

  std::vector<double> phases;
  int n_mean_ipc = 0, n_mean_ipc_weighted = 0, n_plv = 0, n_frac_inphase = 0;

  for (int i=0; i<stats_in.size(); i++)
    {
      const ipc_stats_t & stat = stats_in[i];
      out.n_total += stat.n_total;
      out.n_used += stat.n_used;

      if ( std::isfinite( stat.mean_ipc ) )
        {
          out.mean_ipc += stat.mean_ipc;
          ++n_mean_ipc;
        }
      if ( std::isfinite( stat.mean_ipc_weighted ) )
        {
          out.mean_ipc_weighted += stat.mean_ipc_weighted;
          ++n_mean_ipc_weighted;
        }
      if ( std::isfinite( stat.plv ) )
        {
          out.plv += stat.plv;
          ++n_plv;
        }
      if ( std::isfinite( stat.frac_inphase ) )
        {
          out.frac_inphase += stat.frac_inphase;
          ++n_frac_inphase;
        }
      if ( std::isfinite( stat.mean_phase ) ) phases.push_back( stat.mean_phase );
    }

  out.mean_ipc = n_mean_ipc ? out.mean_ipc / (double)n_mean_ipc : std::numeric_limits<double>::quiet_NaN();
  out.mean_ipc_weighted = n_mean_ipc_weighted ? out.mean_ipc_weighted / (double)n_mean_ipc_weighted : std::numeric_limits<double>::quiet_NaN();
  out.plv = n_plv ? out.plv / (double)n_plv : std::numeric_limits<double>::quiet_NaN();
  out.frac_inphase = n_frac_inphase ? out.frac_inphase / (double)n_frac_inphase : std::numeric_limits<double>::quiet_NaN();
  out.mean_phase = ipc_t::circular_mean( phases ).mean_phase;

  return out;
}

void write_ipc_summary_row( const std::string & ch1 ,
                            const std::string & ch2 ,
                            const ipc_stats_t & stats ,
                            const bool flip_phase = false )
{
  writer.level( ch1 , globals::signal1_strat );
  writer.level( ch2 , globals::signal2_strat );
  writer.value( "N_TOT" , (int)stats.n_total );
  writer.value( "N_USED" , (int)stats.n_used );
  writer.value( "IPC" , stats.mean_ipc );
  writer.value( "WIPC" , stats.mean_ipc_weighted );
  writer.value( "PLV" , stats.plv );
  writer.value( "PHASE" , flip_phase && std::isfinite( stats.mean_phase ) ? -stats.mean_phase : stats.mean_phase );
  writer.value( "P_INPHASE" , stats.frac_inphase );
  writer.unlevel( globals::signal2_strat );
  writer.unlevel( globals::signal1_strat );
}

void write_ipc_summary_values( const ipc_stats_t & stats ,
                               const bool flip_phase = false )
{
  writer.value( "N_TOT" , (int)stats.n_total );
  writer.value( "N_USED" , (int)stats.n_used );
  writer.value( "IPC" , stats.mean_ipc );
  writer.value( "WIPC" , stats.mean_ipc_weighted );
  writer.value( "PLV" , stats.plv );
  writer.value( "PHASE" , flip_phase && std::isfinite( stats.mean_phase ) ? -stats.mean_phase : stats.mean_phase );
  writer.value( "P_INPHASE" , stats.frac_inphase );
}

void write_ipc_epoch_row( const std::string & ch1 ,
                          const std::string & ch2 ,
                          const int epoch ,
                          const ipc_stats_t & stats ,
                          const bool flip_phase = false )
{
  writer.epoch( epoch );
  writer.level( ch1 , globals::signal1_strat );
  writer.level( ch2 , globals::signal2_strat );
  write_ipc_summary_values( stats , flip_phase );
  writer.unlevel( globals::signal2_strat );
  writer.unlevel( globals::signal1_strat );
}

}

void dsptools::ipc( edf_t & edf , param_t & param )
{

  const std::string hilbert_mag_prefix = "_ht_mag";
  const std::string hilbert_phase_prefix = "_ht_ph";
  
  
  // signals
  std::string s1, s2;
  
  bool all_by_all = true;
  
  // take sig or sig1 / sig2
  if ( param.has( "sig1" ) && param.has( "sig2" ) )
    {
      s1 = param.value( "sig1" );
      s2 = param.value( "sig2" );
      all_by_all = false;
    }
  else
    s1 = s2 = param.value( "sig" );
  
  // check signals
  const bool NO_ANNOTS = true;
  signal_list_t signals1 = edf.header.signal_list( s1 , NO_ANNOTS );
  signal_list_t signals2 = edf.header.signal_list( s2 , NO_ANNOTS );
  
  const int ns1 = signals1.size();
  const int ns2 = signals2.size();

  if ( ns1 == 0 || ns2 == 0 )
    {
      logger << "  no pairwise signal comparisons to perform\n";
      return ;
    }

  int np = all_by_all ? ns1*(ns1-1) : ns1 * ns2;


  //
  // check all signals have hilbert phase & amp ( sig_ht_mag and sig_ht_phase ) 
  //

  std::vector<int> slots1_mag, slots1_phase;
  int sr = 0;
  for (int s=0; s<ns1; s++)
    {
      std::string s_mag   = signals1.label(s) + hilbert_mag_prefix;
      std::string s_phase = signals1.label(s) + hilbert_phase_prefix;      
      if ( ! edf.header.has_signal( s_mag ) ) Helper::halt( "could not find " + s_mag );
      if ( ! edf.header.has_signal( s_phase ) ) Helper::halt( "could not find " + s_phase );
      int slot = edf.header.signal( s_mag );
      int phase_slot = edf.header.signal( s_phase );
      slots1_mag.push_back( slot );
      slots1_phase.push_back( phase_slot );
      if ( sr == 0 )
	sr = edf.header.sampling_freq( slot );
      else if (  edf.header.sampling_freq( slot ) != sr )
	Helper::halt( "requires uniform sampling rate across signals" );      
      else if (  edf.header.sampling_freq( phase_slot ) != sr )
	Helper::halt( "requires uniform sampling rate across signals" );
    }


  // 
  std::vector<int> slots2_mag, slots2_phase;
  
  for (int s=0; s<ns2; s++)
    {
      std::string s_mag   = signals2.label(s) + hilbert_mag_prefix;
      std::string s_phase = signals2.label(s) + hilbert_phase_prefix;      
      if ( ! edf.header.has_signal( s_mag ) ) Helper::halt( "could not find " + s_mag );
      if ( ! edf.header.has_signal( s_phase ) ) Helper::halt( "could not find " + s_phase );

      int slot_mag = edf.header.signal( s_mag );
      int slot_phase = edf.header.signal( s_phase );

      // add
      slots2_mag.push_back( slot_mag );
      slots2_phase.push_back( slot_phase );
      
      if ( sr == 0 )
	sr = edf.header.sampling_freq( slot_mag );
      else if (  edf.header.sampling_freq( slot_mag ) != sr )
	Helper::halt( "requires uniform sampling rate across signals" );
      else if (  edf.header.sampling_freq( slot_phase ) != sr )
	Helper::halt( "requires uniform sampling rate across signals" );
      
    }

  //
  // IPC params
  //
  
  ipc_param_t ipc_param;
  if ( param.has( "no-weights" ) ) ipc_param.amplitude_weighting = false;
  if ( param.has( "no-gate" ) ) ipc_param.gate_low_amp = false;
  if ( param.has( "q" ) )
    {
      ipc_param.gate_use_quantile = true;
      ipc_param.gate_quantile = param.requires_dbl( "q" );
      if ( ipc_param.gate_quantile < 0.0 || ipc_param.gate_quantile > 1.0 )
        Helper::halt( "q must be between 0 and 1" );
    }
  if ( param.has( "th" ) )
    {
      ipc_param.gate_use_quantile = false;
      ipc_param.gate_abs = param.requires_dbl( "th" );
    }
  if ( param.has( "edge" ) )
    {
      ipc_param.edge_drop_sec = param.requires_dbl( "edge" );
      if ( ipc_param.edge_drop_sec < 0.0 )
        Helper::halt( "edge must be non-negative" );
    }

  const bool lag_output = param.has( "w" );
  const bool epoch_output = param.has( "epoch" );
  const int max_lag = lag_output ? (int)std::llround( param.requires_dbl( "w" ) * sr ) : 0;

  // channel(s) to add
  const bool add_channels = param.has( "add-channels" );
  const bool add_combined_channels = add_channels && param.value( "add-channels" ) == "seed" ;
  ipc_derived_mode_t ctype = ipc_derived_mode_t::NO_IPC;
  if ( add_channels ) ctype = add_combined_channels ?
			ipc_derived_mode_t::PER_SEED_MEAN_IPCW :
			ipc_derived_mode_t::PER_PAIR_IPCW ; 
  
  const std::string ch_prefix = param.has( "prefix" ) ?
    param.value( "prefix" ) : "IPC_" ; 

  const bool use_existing_epochs = epoch_output && edf.timeline.epoched();
  const int ne = epoch_output ? edf.timeline.first_epoch()
                              : edf.timeline.calc_epochs_contig();
  logger << "  IPC settings: weights="
         << ( ipc_param.amplitude_weighting ? "T" : "F" )
         << ", gate=" << ( ipc_param.gate_low_amp ? "T" : "F" );
  if ( ipc_param.gate_low_amp )
    logger << ", gate-mode=" << ( ipc_param.gate_use_quantile ? "quantile" : "absolute" )
           << ", gate-th=" << ( ipc_param.gate_use_quantile ? ipc_param.gate_quantile : ipc_param.gate_abs );
  logger << ", edge=" << ipc_param.edge_drop_sec << "s\n";
  logger << "  iterating over " << ne << " "
         << ( epoch_output ? ( use_existing_epochs ? "existing" : "default" )
                           : "contig-based" )
         << " epochs\n";

  if ( all_by_all && ! lag_output && ! add_channels )
    {
      std::vector<std::vector<std::vector<ipc_stats_t> > > pair_stats(
        ns1 , std::vector<std::vector<ipc_stats_t> >( ns1 ) );

      while ( 1 )
        {
          int epoch = edf.timeline.next_epoch();
          if ( epoch == -1 ) break;
          interval_t interval = edf.timeline.epoch( epoch );
          const int display_epoch = edf.timeline.display_epoch( epoch );

          std::vector<ipc_phaseamp_t> dat( ns1 );
          for (int s=0; s<ns1; s++)
            {
              slice_t slice_mag( edf , slots1_mag[s] , interval );
              slice_t slice_phase( edf , slots1_phase[s] , interval );
              dat[s].amp = *slice_mag.nonconst_pdata();
              dat[s].phase = *slice_phase.nonconst_pdata();
            }

          for (int s=0; s<ns1; s++)
            for (int s2=s+1; s2<ns1; s2++)
              {
                ipc_stats_t stats = ipc_t::compute_ipc( dat[s] , dat[s2] , sr , ipc_param , false ).summary;
                pair_stats[s][s2].push_back( stats );
                if ( epoch_output )
                  {
                    write_ipc_epoch_row( signals1.label(s) , signals1.label(s2) , display_epoch , stats , false );
                    write_ipc_epoch_row( signals1.label(s2) , signals1.label(s) , display_epoch , stats , true );
                  }
              }
        }

      if ( epoch_output ) writer.unepoch();

      for (int s=0; s<ns1; s++)
        for (int s2=s+1; s2<ns1; s2++)
          {
            ipc_stats_t stats = average_ipc_stats( pair_stats[s][s2] );
            write_ipc_summary_row( signals1.label(s) , signals1.label(s2) , stats , false );
            write_ipc_summary_row( signals1.label(s2) , signals1.label(s) , stats , true );
          }

      return;
    }
  
  //
  // for each seed signal
  //   -- look at each epoch
  //   -- bundle other channels
  //   -- compute IPC stats
  //   -- accumulate channel
  //   -- add channel
  //   -- report mean of per-epoch stats
  //
    
  for (int s=0; s<ns1; s++)
    {
      const std::string s_mag   = signals1.label(s) + hilbert_mag_prefix;
      const std::string s_phase = signals1.label(s) + hilbert_phase_prefix;

      const int seed_mag   = edf.header.signal( s_mag );
      const int seed_phase = edf.header.signal( s_phase );
      
      // track results (over epochs)
      std::vector<ipc_batch_result_t> results;
      std::vector<int> epoch_display;
      std::map<int,std::map<int,std::map<int,std::vector<ipc_stats_t> > > > lagged_results;

      // iterate over epochs      
      while ( 1 ) 
	{
	  
	  int epoch = edf.timeline.next_epoch();
	  if ( epoch == -1 ) break;	  
	  interval_t interval = edf.timeline.epoch( epoch );
          epoch_display.push_back( edf.timeline.display_epoch( epoch ) );
	  
	  // seed channel
	  slice_t slice_seed_mag( edf , seed_mag , interval );
	  slice_t slice_seed_phase( edf , seed_phase , interval );	    

	  std::vector<ipc_phaseamp_t> dat( 1 + ns2 );
	  ipc_phaseamp_t seed;
	  seed.amp   = *slice_seed_mag.nonconst_pdata();
	  seed.phase = *slice_seed_phase.nonconst_pdata();	  
	  dat[0] = seed;
			
	  // compile comparison channels
	  for (int s2=0; s2<ns2; s2++)
	    {
	      slice_t slice_mag( edf , slots2_mag[s2] , interval );
 	      slice_t slice_phase( edf , slots2_phase[s2] , interval );

	      ipc_phaseamp_t d;
	      d.amp   = *slice_mag.nonconst_pdata();
	      d.phase = *slice_phase.nonconst_pdata();
	      dat[1 + s2 ] = d;
	    }
	  
	  // IPC
	  ipc_t ipc;
	  
	  std::vector<int> idx1 = { 0 } ;
	  std::vector<int> idx2;
          idx2.reserve( ns2 );
	  for (int i=1; i<=ns2; i++) idx2.push_back( i ); 
	  
	  ipc_batch_result_t res = ipc.compute_ipc_seed_to_set( dat , idx1, idx2 , sr, 
								ipc_param , ctype );
	  
	  results.push_back( res );

          if ( lag_output )
            {
              for (int i=1; i<=ns2; i++)
                {
                  if ( all_by_all && signals1.label(s) == signals2.label(i-1) ) continue;
                  ipc_lag_output_t lagged = ipc.compute_ipc_lagged( dat[0] , dat[i] , sr ,
                                                                    ipc_param , max_lag );
                  for (int li=0; li<lagged.rows.size(); li++)
                    lagged_results[0][i][ lagged.rows[li].lag ].push_back( lagged.rows[li].summary );
                }
            }
	  
	  // next epoch
	}

      
      
      //
      // summarize output, iterating over pairs of channels
      //
      
      if ( epoch_output ) writer.unepoch();

      writer.level( signals1.label(s) , globals::signal1_strat );	  
      
      const std::vector<ipc_pair_summary_row_t> & p = results[0].summaries;

      for (int j=0; j<p.size(); j++)
	{

	  int seed_idx = p[j].seed_idx;
	  int tgt_idx = p[j].tgt_idx;
	  
	  // ignore seed-self stats
	  //	  if ( signals2( tgt_idx - 1 ) == signals1( s ) ) continue;
	  
	  // output stratified by CH2; note, as seed is slot 0, we need to -1 to align
	  // with signals2
	  writer.level( signals2.label( tgt_idx - 1 ) , globals::signal2_strat );
	  
	  // get average over epochs
	  std::vector<ipc_stats_t> epoch_stats;
	  const int n1 = results.size();
          epoch_stats.reserve( n1 );

	  for (int i=0; i < n1; i++)
	    epoch_stats.push_back( results[i].summaries[j].summary );

          if ( epoch_output )
            for (int i=0; i < n1; i++)
              write_ipc_epoch_row( signals1.label(s) ,
                                   signals2.label( tgt_idx - 1 ) ,
                                   epoch_display[i] ,
                                   results[i].summaries[j].summary );

          ipc_stats_t stats = average_ipc_stats( epoch_stats );
          write_ipc_summary_values( stats );

          if ( lag_output )
            {
              std::map<int,std::vector<ipc_stats_t> > & lag_profiles = lagged_results[seed_idx][tgt_idx];
              std::map<int,std::vector<ipc_stats_t> >::iterator ll = lag_profiles.begin();
              while ( ll != lag_profiles.end() )
                {
                  ipc_stats_t lag_stat = average_ipc_stats( ll->second );
                  writer.level( ll->first , "D" );
                  writer.value( "T" , ll->first / (double)sr );
                  writer.value( "N_TOT" , (int)lag_stat.n_total );
                  writer.value( "N_USED" , (int)lag_stat.n_used );
                  writer.value( "IPC" , lag_stat.mean_ipc );
                  writer.value( "WIPC" , lag_stat.mean_ipc_weighted );
                  writer.value( "PLV" , lag_stat.plv );
                  writer.value( "PHASE" , lag_stat.mean_phase );
                  writer.value( "P_INPHASE" , lag_stat.frac_inphase );
                  writer.unlevel( "D" );
                  ++ll;
                }
            }

          writer.unlevel( globals::signal2_strat );
	  
	}

      writer.unlevel( globals::signal1_strat );
      
      
      //
      // add channel?
      //
      
      if ( add_channels )
	{
	  // if add_combined_channels expect 1 channel per seed
	  // if not, expected 'ns2' channels per seed

	  int nch = results[0].derived.size();
	  int np = 0;
	  for (int j=0; j<results.size(); j++)
	    np += results[j].derived[0].size();
	      
	  for (int j=0; j<nch; j++)
	    {
	      std::string label = ch_prefix + signals1.label(s);
	      if ( ! add_combined_channels )
		label += "_" + signals2.label( j );

	      if ( signals1.label(s) == signals2.label( j ) )
		continue;
	      
	      std::vector<double> xx(np);
	      int idx = 0; 
	      for (int i=0; i<results.size(); i++)
		{
		  const std::vector<double> & xx1 = results[i].derived[j];
		  for (int k=0; k<xx1.size(); k++)
		    xx[idx++] = xx1[k];
		}
	      
	      logger << "  adding " << label << "\n";
	      // add channel	      
	      edf.add_signal( label , sr , xx );
	    }
	}
      
      
      // next seed
    }
  writer.unlevel( globals::signal1_strat );
}


// Key IPC computations (generic sampling rate), assuming you already have:
//  - bandpass_filter(x, sr, f_lo, f_hi, ...)
//  - analytic_signal(x_filtered) -> phase(t), amp(t)  [Hilbert-based]

//
// IPC definition (samplewise):
//   dphi(t) = wrap_to_pi( phi_seed(t) - phi_tgt(t) )
//   IPC(t)  = cos(dphi(t))   // signed instantaneous phase coherence
//
// Optional amplitude weighting + gating:
//   w0(t) = min(amp_seed(t), amp_tgt(t))
//   w(t)  = w0(t) if w0(t) >= thr else 0   // gate low-amp times (phase unreliable)
//   IPCw(t) = w(t) * cos(dphi(t))
//
// Summaries:
//   mean_ipc          = mean(IPC(t)) over used samples
//   mean_ipc_weighted = sum(w(t)*IPC(t))/sum(w(t))
//   plv               = |sum(w(t)*exp(i*dphi(t)))| / sum(w(t))   (weighted PLV)
//   mean_phase        = arg( sum(w(t)*exp(i*dphi(t))) )          (circular mean lag)
//   frac_inphase      = fraction(|dphi(t)| < pi/6) among used samples


#include <vector>
#include <cmath>
#include <complex>
#include <limits>
#include <algorithm>


// ------------------------------ IPC core ------------------------------


ipc_output_t ipc_t::compute_ipc(const ipc_phaseamp_t & seed , 
				const ipc_phaseamp_t & tgt,
				double sr,
				const ipc_param_t & P,
				bool return_timeseries )
{

  ipc_output_t out;

  const size_t N = seed.amp.size();

  if (N == 0 || tgt.amp.size() != N) return out;

  out.summary.n_total = N;
  
  // Build weights w0(t)=min(amp_seed, amp_tgt) or 1.0
  std::vector<double> w0(N, 1.0);
  if (P.amplitude_weighting) {
    w0.resize(N);
    for (size_t t = 0; t < N; ++t) w0[t] = std::min(seed.amp[t], tgt.amp[t]);
  }

  // Determine gating threshold
  double thr = 0.0;
  if (P.amplitude_weighting && P.gate_low_amp) {
    thr = P.gate_use_quantile ? quantile(w0, P.gate_quantile) : P.gate_abs;
    if (!finite(thr)) thr = 0.0;
  }
  
  // Compute per-sample IPC quantities (optionally returned)
  if (return_timeseries) {
    out.ipc.resize(N);
    out.ipcw.resize(N);
    out.dphi.resize(N);
    out.w.resize(N);
  }

  // Summaries (over “used” samples): apply edge drop in seconds
  const size_t edge = (size_t)std::llround(P.edge_drop_sec * sr);
  const size_t lo = std::min(edge, N);
  const size_t hi = (N > edge ? N - edge : 0);
  
  double sum_ipc = 0.0;
  double sum_w = 0.0;
  double sum_w_ipc = 0.0;

  size_t n_used = 0;
  size_t n_inphase = 0;
  
  std::complex<double> sum_e(0.0, 0.0); // for PLV + mean_phase

  for (size_t t = 0; t < N; ++t) {
    // Phase difference and IPC
    double d = wrap_to_pi(seed.phase[t] - tgt.phase[t]);
    double c = std::cos(d);

    // Weight
    double w = 1.0;
    if (P.amplitude_weighting) {
      w = w0[t];
      if (P.gate_low_amp && w < thr) w = 0.0;
    }

    if (return_timeseries) {
      out.dphi[t] = d;
      out.ipc[t]  = c;
      out.w[t]    = w;
      out.ipcw[t] = w * c;
    }
    
    // Summaries: only use interior + finite + (if weighted) w>0
    if (t < lo || t >= hi) continue;
    if (!finite(d) || !finite(c)) continue;
    if (P.amplitude_weighting && w <= 0.0) continue;

    ++n_used;
    sum_ipc += c;

    if (P.amplitude_weighting) {
      sum_w += w;
      sum_w_ipc += w * c;
      sum_e += std::polar(w, d);          // w * exp(i*dphi)
    } else {
      sum_e += std::polar(1.0, d);
    }

    if (std::fabs(d) < M_PI / 6.0) ++n_inphase;
  }

  out.summary.n_used = n_used;
  
  if (n_used > 0) {
    out.summary.mean_ipc = sum_ipc / (double)n_used;
    out.summary.frac_inphase = (double)n_inphase / (double)n_used;
    // PLV and mean phase lag from sum_e
    if (P.amplitude_weighting) {
      if (sum_w > 0.0) {
        out.summary.plv = std::abs(sum_e) / sum_w;
        out.summary.mean_ipc_weighted = sum_w_ipc / sum_w;
        out.summary.mean_phase = std::atan2(sum_e.imag(), sum_e.real());
      }
    } else {
      out.summary.plv = std::abs(sum_e) / (double)n_used;
      out.summary.mean_phase = std::atan2(sum_e.imag(), sum_e.real());
    }
  }
  
  return out;
}

ipc_lag_output_t ipc_t::compute_ipc_lagged(const ipc_phaseamp_t & seed,
					   const ipc_phaseamp_t & tgt,
					   double sr,
					   const ipc_param_t & P,
					   int max_lag)
{
  ipc_lag_output_t out;

  const size_t N = seed.amp.size();
  if ( N == 0 || tgt.amp.size() != N || seed.phase.size() != N || tgt.phase.size() != N )
    return out;

  if ( max_lag < 0 ) max_lag = -max_lag;

  for ( int lag = -max_lag ; lag <= max_lag ; ++lag )
    {
      const size_t seed_start = lag < 0 ? 0 : (size_t)lag;
      const size_t tgt_start = lag < 0 ? (size_t)(-lag) : 0;
      const size_t offset = seed_start > tgt_start ? seed_start : tgt_start;
      const size_t span = N > offset ? N - offset : 0;

      ipc_lag_row_t row;
      row.lag = lag;

      if ( span == 0 )
        {
          out.rows.push_back( row );
          continue;
        }

      ipc_phaseamp_t seed_shifted;
      ipc_phaseamp_t tgt_shifted;
      seed_shifted.phase.assign( seed.phase.begin() + seed_start , seed.phase.begin() + seed_start + span );
      seed_shifted.amp.assign( seed.amp.begin() + seed_start , seed.amp.begin() + seed_start + span );
      tgt_shifted.phase.assign( tgt.phase.begin() + tgt_start , tgt.phase.begin() + tgt_start + span );
      tgt_shifted.amp.assign( tgt.amp.begin() + tgt_start , tgt.amp.begin() + tgt_start + span );

      row.summary = compute_ipc( seed_shifted , tgt_shifted , sr , P , false ).summary;
      out.rows.push_back( row );
    }

  return out;
}


// ------------------------------ Batch: seeds vs targets ------------------------------


ipc_batch_result_t ipc_t::compute_ipc_seed_to_set(const std::vector<ipc_phaseamp_t>& signals,
						  const std::vector<int>& s1,
						  const std::vector<int>& s2,
						  double sr,
						  const ipc_param_t & P,
						  ipc_derived_mode_t mode)
{

  ipc_batch_result_t out;
  
  for (int si : s1) {
    const auto& seed = signals[si];
    const size_t N = seed.amp.size();
    
    std::vector<double> acc(N, 0.0);
    std::vector<double> accw(N, 0.0);
    size_t n_tgt_used = 0;
    
    for (int tj : s2) {
      if (tj == si) continue;
      const auto& tgt = signals[tj];
      if (tgt.amp.size() != N) continue;
      
      const bool want_ts = (mode == ipc_derived_mode_t::PER_PAIR_IPC ||
                            mode == ipc_derived_mode_t::PER_PAIR_IPCW ||
                            mode == ipc_derived_mode_t::PER_SEED_MEAN_IPC ||
                            mode == ipc_derived_mode_t::PER_SEED_MEAN_IPCW);
      
      ipc_output_t R = compute_ipc(seed, tgt, sr, P, want_ts);
      
      // store summary row
      ipc_pair_summary_row_t row;
      row.seed_idx = si;
      row.tgt_idx = tj;
      row.summary = R.summary;
      out.summaries.push_back(row);
      
      // derived channels
      if (mode == ipc_derived_mode_t::PER_PAIR_IPC) {
        out.derived.push_back(R.ipc);        
      } else if (mode == ipc_derived_mode_t::PER_PAIR_IPCW) {
        out.derived.push_back(R.ipcw);
      } else if (mode == ipc_derived_mode_t::PER_SEED_MEAN_IPC) {
        for (size_t t = 0; t < N; ++t) acc[t] += R.ipc[t];
        ++n_tgt_used;
      } else if (mode == ipc_derived_mode_t::PER_SEED_MEAN_IPCW) {
        for (size_t t = 0; t < N; ++t) { acc[t] += R.ipcw[t]; accw[t] += R.w[t]; }
        ++n_tgt_used;
      }
    }

    if (mode == ipc_derived_mode_t::PER_SEED_MEAN_IPC && n_tgt_used > 0) {
      std::vector<double> v(N);
      for (size_t t = 0; t < N; ++t) v[t] = acc[t] / (double)n_tgt_used;
      out.derived.push_back(std::move(v));      
    } else if (mode == ipc_derived_mode_t::PER_SEED_MEAN_IPCW && n_tgt_used > 0) {
      std::vector<double> v(N);
      for (size_t t = 0; t < N; ++t) v[t] = (accw[t] > 0.0 ? acc[t] / accw[t] : 0.0); // or NA sentinel
      out.derived.push_back(std::move(v));
    }
  }

  return out;
}




//
// helpers
//

double ipc_t::wrap_to_pi(double x) { return std::atan2(std::sin(x), std::cos(x)); }

bool ipc_t::finite(double x) { return std::isfinite(x); }

double ipc_t::quantile(std::vector<double> v, double q) {
  v.erase(std::remove_if(v.begin(), v.end(), [](double x){return !finite(x);}), v.end());
  if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
  std::sort(v.begin(), v.end());
  double pos = q * (v.size() - 1);
  size_t lo = (size_t)std::floor(pos);
  size_t hi = (size_t)std::ceil(pos);
  double frac = pos - lo;
  return v[lo] * (1.0 - frac) + v[hi] * frac;
}

ipc_t::circ_stats_t ipc_t::circular_mean(const std::vector<double>& theta)
{
  double s = 0.0;
  double c = 0.0;
  const size_t n = theta.size();
  
  if (n == 0)
    return { NAN, NAN };

  for (double x : theta) {
    s += std::sin(x);
    c += std::cos(x);
  }

  s /= n;
  c /= n;

  circ_stats_t out;
  out.mean_phase = std::atan2(s, c);
  out.R = std::sqrt(s*s + c*c);

  return out;
}
