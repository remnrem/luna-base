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
#include <stdexcept>
#include "edf/edf.h"
#include "edf/slice.h"
#include "param.h"

#include "helper/logger.h"
#include "db/db.h"

extern writer_t writer;
extern logger_t logger;

namespace
{

void validate_ipc_parameters( double sr, const ipc_param_t & p )
{
  if ( ! std::isfinite(sr) || sr <= 0 ||
       ! std::isfinite(p.gate_quantile) || p.gate_quantile < 0 || p.gate_quantile > 1 ||
       ! std::isfinite(p.gate_abs) || p.gate_abs < 0 ||
       ! std::isfinite(p.edge_drop_sec) || p.edge_drop_sec < 0 || ! std::isfinite(p.edge_drop_sec * sr) )
    throw std::invalid_argument("IPC requires a finite positive sampling rate and finite valid gating/edge parameters");
}

void validate_ipc_vectors( const ipc_phaseamp_t & seed, const ipc_phaseamp_t & target )
{
  if ( seed.amp.size() != seed.phase.size() || target.amp.size() != target.phase.size() ||
       seed.amp.size() != target.amp.size() )
    throw std::invalid_argument("IPC amplitude and phase vectors must have matching lengths");
}

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

void write_ipc_summary_row( const std::string & ch1 ,
                            const std::string & ch2 ,
                            const ipc_stats_t & stats ,
                            const bool flip_phase = false )
{
  writer.unepoch();
  writer.level( ch1 , globals::signal1_strat );
  writer.level( ch2 , globals::signal2_strat );
  write_ipc_summary_values(stats, flip_phase);
  writer.unlevel( globals::signal2_strat );
  writer.unlevel( globals::signal1_strat );
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
  writer.unepoch();
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



  //
  // check all signals have hilbert phase & amp ( sig_ht_mag and sig_ht_phase ) 
  //

  std::vector<int> slots1_mag, slots1_phase;
  double sr = 0;
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
      if ( ! std::isfinite(sr) || sr <= 0 || edf.header.sampling_freq( phase_slot ) != sr )
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
      if ( ! std::isfinite(sr) || sr <= 0 || edf.header.sampling_freq( slot_phase ) != sr )
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
      if ( ! std::isfinite(ipc_param.gate_quantile) || ipc_param.gate_quantile < 0.0 || ipc_param.gate_quantile > 1.0 )
        Helper::halt( "IPC q must be finite and between 0 and 1" );
    }
  if ( param.has( "th" ) )
    {
      ipc_param.gate_use_quantile = false;
      ipc_param.gate_abs = param.requires_dbl( "th" );
      if ( ! std::isfinite(ipc_param.gate_abs) || ipc_param.gate_abs < 0 ) Helper::halt("IPC th must be finite and nonnegative");
    }
  if ( param.has( "edge" ) )
    {
      ipc_param.edge_drop_sec = param.requires_dbl( "edge" );
      if ( ! std::isfinite(ipc_param.edge_drop_sec) || ipc_param.edge_drop_sec < 0.0 || ! std::isfinite(ipc_param.edge_drop_sec * sr) )
        Helper::halt( "IPC edge must be finite and nonnegative" );
    }

  const bool lag_output = param.has( "w" );
  const bool epoch_output = param.has( "epoch" );
  const double lag_seconds = lag_output ? param.requires_dbl("w") : 0;
  if ( ! std::isfinite(lag_seconds) || lag_seconds < 0 ||
       ! std::isfinite(lag_seconds * sr) || lag_seconds * sr >= std::numeric_limits<int>::max() - 1 )
    Helper::halt("IPC w must be finite, nonnegative, and representable as a sample lag");
  const int max_lag = std::llround(lag_seconds * sr);

  // channel(s) to add
  const bool add_channels = param.has( "add-channels" );
  const bool add_combined_channels = add_channels && param.value( "add-channels" ) == "seed" ;
  ipc_derived_mode_t ctype = ipc_derived_mode_t::NO_IPC;
  if ( add_channels ) ctype = add_combined_channels ?
			ipc_derived_mode_t::PER_SEED_MEAN_IPCW :
			ipc_derived_mode_t::PER_PAIR_IPCW ; 
  
  const std::string ch_prefix = param.has( "prefix" ) ?
    param.value( "prefix" ) : "IPC_" ; 

  // unique: emit each unordered pair once; (B,A) only differs from (A,B) by PHASE sign
  // (and, for lag output, by the sign of D), so drop the reverse when both are visited
  const bool unique_pairs = param.has( "unique" );
  if ( unique_pairs && add_combined_channels )
    Helper::halt( "IPC unique cannot be combined with add-channels=seed" );
  std::map<int,int> sigpos1; // EDF slot -> position in sig/sig1 list
  std::set<int> sigset2;
  for (int s=0; s<ns1; s++) sigpos1.insert( std::make_pair( signals1(s) , s ) );
  for (int s=0; s<ns2; s++) sigset2.insert( signals2(s) );
  auto keep_pair = [&]( const int s , const int t ) {
    const int a = signals1(s), b = signals2(t);
    if ( a == b ) return false;
    if ( ! unique_pairs ) return true;
    // reverse (b,a) is also visited: keep the orientation following sig/sig1 list order
    if ( sigpos1.count(b) && sigset2.count(a) ) return sigpos1[a] < sigpos1[b];
    return true;
  };

  // Build an explicit work list: IPC never consumes or replaces caller epoch state.
  std::vector<interval_t> intervals;
  std::vector<int> display_epochs;
  const bool use_existing_epochs = epoch_output && edf.timeline.epoched();
  if ( epoch_output )
    {
      // A local copy supplies default epoching without changing an unepoched caller.
      timeline_t local = edf.timeline;
      if ( ! local.epoched() ) local.first_epoch();
      for ( int e = 0; e < local.num_total_epochs(); ++e )
        if ( ! local.is_epoch_mask_set() || ! local.masked(e) )
          {
            intervals.push_back( local.epoch(e) );
            display_epochs.push_back( local.display_epoch(e) );
          }
    }
  else
    {
      // Merge adjacent retained records directly, including EDF+D record onsets.
      // Contig mode retains its historical mask-independent analysis scope.
      for ( int r = edf.timeline.first_record(); r != -1; r = edf.timeline.next_record(r) )
        {
          const uint64_t start = edf.timeline.timepoint(r);
          const uint64_t stop = start + edf.header.record_duration_tp;
          if ( ! intervals.empty() && intervals.back().stop == start ) intervals.back().stop = stop;
          else intervals.push_back( interval_t(start, stop) );
        }
    }
  logger << "  IPC settings: weights=" << (ipc_param.amplitude_weighting ? "T" : "F")
         << ", gate=" << (ipc_param.amplitude_weighting && ipc_param.gate_low_amp ? "T" : "F");
  if ( ! ipc_param.amplitude_weighting ) logger << " (no-weights also disables amplitude gating)";
  else if ( ipc_param.gate_low_amp )
    logger << ", gate-mode=" << (ipc_param.gate_use_quantile ? "quantile" : "absolute")
           << ", gate-th=" << (ipc_param.gate_use_quantile ? ipc_param.gate_quantile : ipc_param.gate_abs);
  logger << ", edge=" << ipc_param.edge_drop_sec << "s\n"
         << "  iterating over " << intervals.size() << " "
         << (epoch_output ? (use_existing_epochs ? "existing" : "default") : "contig-based")
         << " intervals; summaries equally weight usable interval statistics\n";
  if ( intervals.empty() )
    {
      writer.unepoch();
      logger << "  no usable IPC intervals; no results or derived channels\n";
      return;
    }

  // Map every retained sample once; EDF+D gaps occupy no synthetic buffer slots.
  std::vector<uint64_t> recording_tp;
  std::vector<size_t> interval_offset;
  if ( add_channels )
    {
      slice_t full(edf, slots1_mag.front(), edf.timeline.wholetrace());
      recording_tp = *full.ptimepoints();
      if ( recording_tp.size() != size_t(edf.header.nr) * edf.header.n_samples[slots1_mag.front()] )
        Helper::halt("IPC derived channels require complete retained recording samples");
      for ( const auto & interval : intervals )
        {
          if ( edf.timeline.valid_tps(interval) != interval.duration() )
            Helper::halt("IPC derived channels do not support epochs spanning EDF+D gaps");
          slice_t slice(edf, slots1_mag.front(), interval);
          const auto & tp = *slice.ptimepoints();
          auto begin = tp.empty() ? recording_tp.end() : std::lower_bound(recording_tp.begin(), recording_tp.end(), tp.front());
          const size_t offset = begin - recording_tp.begin();
          if ( tp.size() > recording_tp.size() - offset || ! std::equal(tp.begin(), tp.end(), begin) )
            Helper::halt("IPC interval samples do not align to the recording sample grid");
          interval_offset.push_back(offset);
        }
      // Preflight all names before adding any channels; same-identity targets are absent.
      std::set<std::string> labels;
      for ( int seed = 0; seed < ns1; ++seed )
        for ( int target = 0; target < ns2; ++target )
          if ( keep_pair( seed , target ) )
            {
              const std::string label = ch_prefix + signals1.label(seed) +
                (add_combined_channels ? "" : "_" + signals2.label(target));
              if ( add_combined_channels && labels.count(label) ) continue;
              if ( edf.header.has_signal(label) || ! labels.insert(label).second )
                Helper::halt("IPC derived channel label already exists: " + label);
            }
      logger << "  IPC derived channels: recording sample placement, zero fill for uncovered/invalid samples,"
             << " mean of finite epoch contributions at overlapping samples\n";
    }

  if ( all_by_all && ! lag_output && ! add_channels )
    {
      std::vector<std::vector<std::vector<ipc_stats_t> > > pair_stats(
        ns1 , std::vector<std::vector<ipc_stats_t> >( ns1 ) );

      for ( size_t ei = 0; ei < intervals.size(); ++ei )
        {
          const interval_t & interval = intervals[ei];
          const int display_epoch = epoch_output ? display_epochs[ei] : 0;

          std::vector<ipc_phaseamp_t> dat( ns1 );
          for (int s=0; s<ns1; s++)
            {
              slice_t slice_mag( edf , slots1_mag[s] , interval );
              slice_t slice_phase( edf , slots1_phase[s] , interval );
              dat[s].amp = *slice_mag.nonconst_pdata();
              dat[s].phase = *slice_phase.nonconst_pdata();
            }

          for ( const auto & signal : dat ) validate_ipc_vectors(dat.front(), signal);
          for (int s=0; s<ns1; s++)
            for (int s2=s+1; s2<ns1; s2++)
              {
                if ( signals1(s) == signals1(s2) ) continue;
                ipc_stats_t stats = ipc_t::compute_ipc( dat[s] , dat[s2] , sr , ipc_param , false ).summary;
                pair_stats[s][s2].push_back( stats );
                if ( epoch_output )
                  {
                    write_ipc_epoch_row( signals1.label(s) , signals1.label(s2) , display_epoch , stats , false );
                    if ( ! unique_pairs )
                      write_ipc_epoch_row( signals1.label(s2) , signals1.label(s) , display_epoch , stats , true );
                  }
              }
        }

      for (int s=0; s<ns1; s++)
        for (int s2=s+1; s2<ns1; s2++)
          {
            if ( signals1(s) == signals1(s2) ) continue;
            ipc_stats_t stats = average_ipc_stats( pair_stats[s][s2] );
            write_ipc_summary_row( signals1.label(s) , signals1.label(s2) , stats , false );
            if ( ! unique_pairs )
              write_ipc_summary_row( signals1.label(s2) , signals1.label(s) , stats , true );
          }

      return;
    }
  
  // Each seed traverses the same immutable interval list, independently.
  for ( int s = 0; s < ns1; ++s )
    {
      std::vector<int> targets;
      for ( int t = 0; t < ns2; ++t )
        if ( keep_pair( s , t ) ) targets.push_back(t);
      if ( targets.empty() ) continue;

      std::vector<ipc_batch_result_t> results;
      std::map<int,std::map<int,std::vector<ipc_stats_t> > > lagged_results;
      for ( size_t ei = 0; ei < intervals.size(); ++ei )
        {
          const interval_t & interval = intervals[ei];
          std::vector<ipc_phaseamp_t> dat(1 + targets.size());
          slice_t seed_mag(edf, slots1_mag[s], interval);
          slice_t seed_phase(edf, slots1_phase[s], interval);
          dat[0].amp = *seed_mag.pdata();
          dat[0].phase = *seed_phase.pdata();
          std::vector<int> idx2;
          for ( size_t j = 0; j < targets.size(); ++j )
            {
              slice_t mag(edf, slots2_mag[targets[j]], interval);
              slice_t phase(edf, slots2_phase[targets[j]], interval);
              dat[j+1].amp = *mag.pdata();
              dat[j+1].phase = *phase.pdata();
              idx2.push_back(j+1);
            }
          results.push_back( ipc_t::compute_ipc_seed_to_set(dat, {0}, idx2, sr, ipc_param, ctype) );
          if ( lag_output )
            for ( int j : idx2 )
              for ( const auto & row : ipc_t::compute_ipc_lagged(dat[0], dat[j], sr, ipc_param, max_lag).rows )
                lagged_results[j][row.lag].push_back(row.summary);
        }

      // Empty inputs have already returned; derive target labels from the filtered list.
      for ( size_t j = 0; j < targets.size(); ++j )
        {
          const std::string ch1 = signals1.label(s), ch2 = signals2.label(targets[j]);
          std::vector<ipc_stats_t> epoch_stats;
          for ( size_t ei = 0; ei < results.size(); ++ei )
            {
              const ipc_stats_t & stats = results[ei].summaries[j].summary;
              epoch_stats.push_back(stats);
              if ( epoch_output ) write_ipc_epoch_row(ch1, ch2, display_epochs[ei], stats);
            }
          write_ipc_summary_row(ch1, ch2, average_ipc_stats(epoch_stats));
          if ( lag_output )
            for ( const auto & lag : lagged_results[j+1] )
              {
                writer.unepoch();
                writer.level(ch1, globals::signal1_strat);
                writer.level(ch2, globals::signal2_strat);
                writer.level(lag.first, "D");
                writer.value("T", lag.first / sr);
                write_ipc_summary_values(average_ipc_stats(lag.second));
                writer.unlevel("D");
                writer.unlevel(globals::signal2_strat);
                writer.unlevel(globals::signal1_strat);
              }
        }

      if ( add_channels )
        {
          const size_t nch = results.front().derived.size();
          for ( size_t j = 0; j < nch; ++j )
            {
              std::string label = ch_prefix + signals1.label(s);
              if ( ! add_combined_channels ) label += "_" + signals2.label(targets[j]);
              std::vector<double> values(recording_tp.size(), 0.0);
              std::vector<size_t> counts(recording_tp.size(), 0);
              for ( size_t ei = 0; ei < results.size(); ++ei )
                {
                  const auto & source = results[ei].derived[j];
                  const size_t offset = interval_offset[ei];
                  if ( source.size() > values.size() - offset )
                    Helper::halt("IPC derived slice exceeds recording buffer");
                  for ( size_t k = 0; k < source.size(); ++k )
                    if ( std::isfinite(source[k]) )
                      {
                        values[offset+k] += source[k];
                        ++counts[offset+k];
                      }
                }
              for ( size_t k = 0; k < values.size(); ++k )
                if ( counts[k] ) values[k] /= counts[k];
              logger << "  adding " << label << "\n";
              edf.add_signal(label, sr, values);
            }
        }
    }
}


// Key IPC computations use precomputed Hilbert phase and magnitude vectors.

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


// ------------------------------ IPC core ------------------------------


ipc_output_t ipc_t::compute_ipc(const ipc_phaseamp_t & seed , 
				const ipc_phaseamp_t & tgt,
				double sr,
				const ipc_param_t & P,
				bool return_timeseries )
{

  ipc_output_t out;

  const size_t N = seed.amp.size();

  validate_ipc_parameters(sr, P);
  validate_ipc_vectors(seed, tgt);
  if (N == 0) return out;

  out.summary.n_total = N;
  
  // A sample requires finite phase AND amplitude inputs, even with no-weights.
  const double missing = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> w0(N, missing);
  for ( size_t t = 0; t < N; ++t )
    if ( finite(seed.amp[t]) && finite(tgt.amp[t]) &&
         finite(seed.phase[t]) && finite(tgt.phase[t]) && finite(seed.phase[t] - tgt.phase[t]) )
      w0[t] = P.amplitude_weighting ? std::min(seed.amp[t], tgt.amp[t]) : 1.0;

  // Determine gating threshold
  double thr = 0.0;
  if (P.amplitude_weighting && P.gate_low_amp) {
    thr = P.gate_use_quantile ? quantile(w0, P.gate_quantile) : P.gate_abs;
    if (!finite(thr)) thr = 0.0;
  }
  
  // Compute per-sample IPC quantities (optionally returned)
  if (return_timeseries) {
    out.ipc.assign(N, missing);
    out.ipcw.assign(N, missing);
    out.dphi.assign(N, missing);
    out.w.assign(N, 0.0);
  }

  // Summaries (over “used” samples): apply edge drop in seconds
  const double edge_samples = P.edge_drop_sec * sr;
  const size_t edge = edge_samples >= N ? N : (size_t)std::llround(edge_samples);
  const size_t lo = std::min(edge, N);
  const size_t hi = (N > edge ? N - edge : 0);
  
  double sum_ipc = 0.0;
  double sum_w = 0.0;
  double sum_w_ipc = 0.0;

  size_t n_used = 0;
  size_t n_inphase = 0;
  
  std::complex<double> sum_e(0.0, 0.0); // for PLV + mean_phase

  for (size_t t = 0; t < N; ++t) {
    if ( ! finite(w0[t]) ) continue;
    // Phase difference and IPC
    double d = wrap_to_pi(seed.phase[t] - tgt.phase[t]);
    double c = std::cos(d);
    if ( ! finite(d) || ! finite(c) ) continue;

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
  validate_ipc_parameters(sr, P);
  validate_ipc_vectors(seed, tgt);
  const int64_t magnitude = max_lag < 0 ? -int64_t(max_lag) : int64_t(max_lag);
  if ( magnitude >= std::numeric_limits<int>::max() )
    throw std::invalid_argument("IPC lag magnitude is not representable");
  if ( N == 0 ) return out;
  max_lag = int(magnitude);

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

// GCC 13 false positive: sized-delete range analysis misidentifies N*sizeof(double) as a pointer offset
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfree-nonheap-object"
ipc_batch_result_t ipc_t::compute_ipc_seed_to_set(const std::vector<ipc_phaseamp_t>& signals,
						  const std::vector<int>& s1,
						  const std::vector<int>& s2,
						  double sr,
						  const ipc_param_t & P,
						  ipc_derived_mode_t mode)
{

  validate_ipc_parameters(sr, P);
  for ( int i : s1 )
    if ( i < 0 || size_t(i) >= signals.size() ) throw std::invalid_argument("IPC seed index out of range");
  for ( int i : s2 )
    if ( i < 0 || size_t(i) >= signals.size() ) throw std::invalid_argument("IPC target index out of range");
  for ( int si : s1 )
    {
      validate_ipc_vectors(signals[si], signals[si]);
      for ( int tj : s2 ) validate_ipc_vectors(signals[si], signals[tj]);
    }
  ipc_batch_result_t out;
  
  for (int si : s1) {
    const auto& seed = signals[si];
    const size_t N = seed.amp.size();
    
    std::vector<double> acc(N, 0.0);
    std::vector<double> accw(N, 0.0);
    size_t n_tgt_used = 0;
    std::vector<size_t> n_valid(N, 0);
    
    for (int tj : s2) {
      if (tj == si) continue;
      const auto& tgt = signals[tj];
      
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
        for (size_t t = 0; t < N; ++t)
          if ( finite(R.ipc[t]) ) { acc[t] += R.ipc[t]; ++n_valid[t]; }
        ++n_tgt_used;
      } else if (mode == ipc_derived_mode_t::PER_SEED_MEAN_IPCW) {
        for (size_t t = 0; t < N; ++t)
          if ( finite(R.ipcw[t]) && finite(R.w[t]) )
            { acc[t] += R.ipcw[t]; accw[t] += R.w[t]; ++n_valid[t]; }
        ++n_tgt_used;
      }
    }

    if (mode == ipc_derived_mode_t::PER_SEED_MEAN_IPC && n_tgt_used > 0) {
      std::vector<double> v(N);
      for (size_t t = 0; t < N; ++t) v[t] = n_valid[t] ? acc[t] / (double)n_valid[t] : std::numeric_limits<double>::quiet_NaN();
      out.derived.push_back(std::move(v));      
    } else if (mode == ipc_derived_mode_t::PER_SEED_MEAN_IPCW && n_tgt_used > 0) {
      std::vector<double> v(N);
      for (size_t t = 0; t < N; ++t) v[t] = ! n_valid[t] ? std::numeric_limits<double>::quiet_NaN() : (accw[t] > 0.0 ? acc[t] / accw[t] : 0.0);
      out.derived.push_back(std::move(v));
    }
  }

  return out;
}
#pragma GCC diagnostic pop




//
// helpers
//

double ipc_t::wrap_to_pi(double x) { return std::atan2(std::sin(x), std::cos(x)); }

bool ipc_t::finite(double x) { return std::isfinite(x); }

double ipc_t::quantile(std::vector<double> v, double q) {
  if ( ! finite(q) || q < 0 || q > 1 ) throw std::invalid_argument("IPC quantile must be finite and in [0,1]");
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
