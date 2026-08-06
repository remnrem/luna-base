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

#include "dsp/arousals2.h"
#include "param.h"
#include "defs/defs.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "stats/eigen_ops.h"
#include "miscmath/miscmath.h"
#include "fftw/fftwrap.h"
#include "dsp/spline.h"

#include "helper/logger.h"
#include "db/db.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <numeric>

extern logger_t logger;
extern writer_t writer;

namespace {

bool arousals2_verbose_annots = false;
std::string arousals2_verbose_prefix = "arv_";

std::string arousals2_verbose_label( const std::string & label )
{
  // Keep the diagnostic tracks in one shared NREM/REM order so that
  // visual review does not interleave stage-specific copies of the same
  // pipeline step.
  static const std::vector<std::string> labels = {
    "eeg_fast_rise_pass",
    "eeg_tilt_rise_pass",
    "eeg_seed_supported",
    "eeg_fast_active",
    "emg_rise_seed",
    "delta_artifact",
    "eeg_artifact",
    "candidate_eeg_raw",
    "candidate_emg_raw",
    "candidate_raw",
    "candidate_duration_ok",
    "candidate_merged",
    "candidate_duration_capped",
    "candidate_presleep_ok",
    "candidate_artifact_ok",
    "candidate_level_refined",
    "candidate_final_duration_ok",
    "rem_emg_rise",
    "rem_emg_rise_duration_ok",
    "rem_emg_confirmed",
    "rem_emg_rejected"
  };

  for ( size_t i = 0; i < labels.size(); ++i )
    if ( labels[i] == label )
      {
        const std::string n = i < 10 ? "0" + std::to_string( i ) : std::to_string( i );
        return arousals2_verbose_prefix + n + "_" + label;
      }

  return arousals2_verbose_prefix + "99_" + label;
}

bool arousals2_bool_token( const std::string & s )
{
  return s == "T" || s == "F" || s == "t" || s == "f" ||
    s == "Y" || s == "N" || s == "y" || s == "n" ||
    s == "1" || s == "0" || s == "true" || s == "false" ||
    s == "TRUE" || s == "FALSE";
}

bool arousals2_true_token( const std::string & s )
{
  return s == "T" || s == "t" || s == "Y" || s == "y" ||
    s == "1" || s == "true" || s == "TRUE";
}

struct arousals2_manual_stage_t {
  std::set<interval_t> events;
  std::vector<interval_t> stage_intervals;
  double stage_sec = 0;
};

bool arousals2_stage_label( const std::string & id , std::string * stage )
{
  if ( id == "R" ) { *stage = "R"; return true; }
  if ( id == "NR" || id == "N1" || id == "N2" || id == "N3" || id == "NREM4" )
    { *stage = "NR"; return true; }
  return false;
}

uint64_t arousals2_overlap_tp( const interval_t & a , const interval_t & b )
{
  const uint64_t s = std::max( a.start , b.start );
  const uint64_t e = std::min( a.stop , b.stop );
  return e > s ? e - s : 0;
}

uint64_t arousals2_set_union_tp( const std::set<interval_t> & events )
{
  if ( events.empty() ) return 0;
  uint64_t start = events.begin()->start;
  uint64_t stop = events.begin()->stop;
  uint64_t total = 0;
  for ( std::set<interval_t>::const_iterator i = std::next( events.begin() );
        i != events.end(); ++i )
    {
      if ( i->start <= stop ) stop = std::max( stop , i->stop );
      else { total += stop - start; start = i->start; stop = i->stop; }
    }
  return total + stop - start;
}

std::map<std::string,arousals2_manual_stage_t>
arousals2_collect_manual( edf_t & edf , const std::string & label )
{
  std::map<std::string,arousals2_manual_stage_t> out;
  annot_t * manual = edf.annotations->find( label );
  if ( manual == NULL ) Helper::halt( "AROUSALS manual annotation not found: " + label );

  edf.annotations->make_sleep_stage( edf.timeline );
  annot_t * staging = edf.annotations->find( "SleepStage" );
  if ( staging == NULL ) Helper::halt( "AROUSALS manual mode requires sleep staging" );

  for ( annot_map_t::const_iterator i = staging->interval_events.begin();
        i != staging->interval_events.end(); ++i )
    {
      std::string stage;
      if ( arousals2_stage_label( i->first.id , &stage ) )
        {
          out[ stage ].stage_intervals.push_back( i->first.interval );
          out[ stage ].stage_sec += i->first.interval.duration_sec();
        }
    }

  for ( annot_map_t::const_iterator i = manual->interval_events.begin();
        i != manual->interval_events.end(); ++i )
    {
      std::string best_stage;
      uint64_t best_overlap = 0;
      for ( annot_map_t::const_iterator s = staging->interval_events.begin();
            s != staging->interval_events.end(); ++s )
        {
          std::string stage;
          if ( ! arousals2_stage_label( s->first.id , &stage ) ) continue;
          const uint64_t overlap = arousals2_overlap_tp( i->first.interval , s->first.interval );
          if ( overlap > best_overlap ) { best_overlap = overlap; best_stage = stage; }
        }
      if ( best_overlap > 0 ) out[ best_stage ].events.insert( i->first.interval );
    }
  return out;
}

struct arousals2_bin_metrics_t {
  uint64_t n = 0, tp = 0, fp = 0, tn = 0, fn = 0;
  double precision = 0, recall = 0, f1 = 0, kappa = 0;
};

arousals2_bin_metrics_t arousals2_bin_compare(
    const std::vector<interval_t> & stage_intervals,
    const std::set<interval_t> & automated,
    const std::set<interval_t> & manual )
{
  const uint64_t bin_tp = 30 * globals::tp_1sec;
  std::set<uint64_t> valid_bins, auto_bins, manual_bins;

  auto add_stage_bins = [ & ]( const interval_t & e )
    {
      if ( e.stop <= e.start ) return;
      const uint64_t b0 = e.start / bin_tp;
      const uint64_t b1 = ( e.stop - 1 ) / bin_tp;
      for ( uint64_t b = b0; b <= b1; ++b ) valid_bins.insert( b );
    };

  auto add_event_bin = [ & ]( const interval_t & e , std::set<uint64_t> * bins )
    {
      if ( e.stop <= e.start ) return;
      const uint64_t duration = e.stop - e.start;
      if ( duration <= 2 * bin_tp )
        {
          // Match the reference 30-s evaluation: short events contribute
          // only to the bin containing their midpoint.
          const uint64_t middle = e.start + std::min( duration , 6 * globals::tp_1sec ) / 2;
          bins->insert( middle / bin_tp );
        }
      else
        {
          // For unusually long events, mark only the interior bins, as in
          // the reference event-to-time-array implementation.
          const uint64_t b0 = e.start / bin_tp + 1;
          const uint64_t b1 = e.stop / bin_tp;
          for ( uint64_t b = b0; b < b1; ++b ) bins->insert( b );
        }
    };

  for ( const interval_t & e : stage_intervals ) add_stage_bins( e );
  for ( const interval_t & e : automated ) add_event_bin( e , &auto_bins );
  for ( const interval_t & e : manual ) add_event_bin( e , &manual_bins );

  arousals2_bin_metrics_t out;
  for ( const uint64_t b : valid_bins )
    {
      const bool a = auto_bins.find( b ) != auto_bins.end();
      const bool m = manual_bins.find( b ) != manual_bins.end();
      if ( a && m ) ++out.tp;
      else if ( a ) ++out.fp;
      else if ( m ) ++out.fn;
      else ++out.tn;
    }
  out.n = out.tp + out.fp + out.tn + out.fn;
  out.precision = out.tp + out.fp > 0 ? out.tp / (double)( out.tp + out.fp ) : 0;
  out.recall = out.tp + out.fn > 0 ? out.tp / (double)( out.tp + out.fn ) : 0;
  out.f1 = out.precision + out.recall > 0 ?
    2.0 * out.precision * out.recall / ( out.precision + out.recall ) : 0;
  if ( out.n > 0 )
    {
      const double p0 = ( out.tp + out.tn ) / (double)out.n;
      const double pe =
        ( ( out.tp + out.fn ) * ( out.tp + out.fp ) +
          ( out.fp + out.tn ) * ( out.fn + out.tn ) ) /
        ( (double)out.n * out.n );
      out.kappa = 1.0 - pe > 0 ? ( p0 - pe ) / ( 1.0 - pe ) : 0;
    }
  return out;
}

bool arousals2_manual_match( const interval_t & automated ,
                              const interval_t & manual ,
                              const uint64_t min_overlap_tp )
{
  const uint64_t overlap = arousals2_overlap_tp( automated , manual );
  if ( overlap < min_overlap_tp || overlap == 0 ) return false;
  return true;
}

void arousals2_write_manual( const std::map<std::string,arousals2_manual_stage_t> & manual ,
                             const std::map<std::string,std::set<interval_t> > & automated ,
                             const uint64_t min_overlap_tp ,
                             const bool eval_by_bin )
{
  for ( const std::string stage : { std::string( "NR" ) , std::string( "R" ) } )
    {
      const std::map<std::string,arousals2_manual_stage_t>::const_iterator m = manual.find( stage );
      const std::set<interval_t> empty;
      const std::set<interval_t> & me = m == manual.end() ? empty : m->second.events;
      const double stage_sec = m == manual.end() ? 0 : m->second.stage_sec;
      double duration = 0, duration_long = 0, duration_std = 0;
      int n_long = 0, n_std = 0;
      for ( std::set<interval_t>::const_iterator e = me.begin(); e != me.end(); ++e )
        {
          const double d = e->duration_sec();
          duration += d;
          if ( d > 15.0 ) { ++n_long; duration_long += d; }
          else { ++n_std; duration_std += d; }
        }

      writer.level( stage , "SS" );
      writer.value( "MINS" , stage_sec / 60.0 );
      // Canonical manual metrics: default = all events; explicit STD/LONG
      // fields split the same events at the 15-second boundary.
      writer.value( "N_MAN" , (int)me.size() );
      writer.value( "AI_MAN" , stage_sec > 0 ? me.size() / ( stage_sec / 3600.0 ) : 0.0 );
      if ( ! me.empty() ) writer.value( "DUR_MAN" , duration / (double)me.size() );
      writer.value( "N_MAN_STD" , n_std );
      writer.value( "AI_MAN_STD" , stage_sec > 0 ? n_std / ( stage_sec / 3600.0 ) : 0.0 );
      if ( n_std > 0 ) writer.value( "DUR_MAN_STD" , duration_std / (double)n_std );
      writer.value( "N_MAN_LONG" , n_long );
      writer.value( "AI_MAN_LONG" , stage_sec > 0 ? n_long / ( stage_sec / 3600.0 ) : 0.0 );
      if ( n_long > 0 ) writer.value( "DUR_MAN_LONG" , duration_long / (double)n_long );
      writer.level( "manual" , "CLS" );
      writer.value( "NE" , (int)me.size() );
      writer.value( "AI" , stage_sec > 0 ? me.size() / ( stage_sec / 3600.0 ) : 0.0 );
      if ( ! me.empty() ) writer.value( "DUR" , duration / (double)me.size() );
      writer.unlevel( "CLS" );
      writer.unlevel( "SS" );

      const std::string auto_label = stage == "NR" ? "arousal_nrem" : "arousal_rem";
      const std::map<std::string,std::set<interval_t> >::const_iterator a = automated.find( auto_label );
      const std::set<interval_t> & ae = a == automated.end() ? empty : a->second;
      struct candidate_t { uint64_t overlap; int ai; int mi; };
      std::vector<interval_t> av( ae.begin() , ae.end() ), mv( me.begin() , me.end() );
      std::vector<candidate_t> candidates;
      for ( int ai = 0; ai < av.size(); ++ai )
        for ( int mi = 0; mi < mv.size(); ++mi )
          if ( arousals2_manual_match( av[ai] , mv[mi] , min_overlap_tp ) )
            candidates.push_back( { arousals2_overlap_tp( av[ai] , mv[mi] ) , ai , mi } );
      std::sort( candidates.begin() , candidates.end() ,
                 []( const candidate_t & x , const candidate_t & y ) { return x.overlap > y.overlap; } );
      std::set<int> used_a, used_m;
      uint64_t overlap_tp = 0;
      uint64_t union_tp = 0;
      double prop_auto = 0, prop_manual = 0, start_diff = 0, start_diff_abs = 0, dur_diff = 0;
      int tp = 0;
      for ( const candidate_t & c : candidates )
        if ( used_a.insert( c.ai ).second && used_m.insert( c.mi ).second )
          {
            ++tp;
            overlap_tp += c.overlap;
            union_tp += ( av[c.ai].stop - av[c.ai].start ) +
              ( mv[c.mi].stop - mv[c.mi].start ) - c.overlap;
            prop_auto += c.overlap / (double)( av[c.ai].stop - av[c.ai].start );
            prop_manual += c.overlap / (double)( mv[c.mi].stop - mv[c.mi].start );
            const double sd = ( (double)av[c.ai].start - (double)mv[c.mi].start ) /
              (double)globals::tp_1sec;
            start_diff += sd;
            start_diff_abs += std::fabs( sd );
            dur_diff += ( (double)( av[c.ai].stop - av[c.ai].start ) -
                          (double)( mv[c.mi].stop - mv[c.mi].start ) ) /
              (double)globals::tp_1sec;
          }
      const int fp = (int)av.size() - tp;
      const int fn = (int)mv.size() - tp;
      const double precision = av.empty() ? 0 : tp / (double)av.size();
      const double recall = mv.empty() ? 0 : tp / (double)mv.size();
      const double f1 = precision + recall > 0 ? 2.0 * precision * recall / ( precision + recall ) : 0;
      writer.level( "1" , "COMP" );
      writer.level( stage , "SS" );
      writer.value( "N_AUTO" , (int)av.size() );
      writer.value( "N_MANUAL" , (int)mv.size() );
      writer.value( "TP" , tp );
      writer.value( "FP" , fp );
      writer.value( "FN" , fn );
      writer.value( "PRECISION" , precision );
      writer.value( "RECALL" , recall );
      writer.value( "F1" , f1 );
      if ( eval_by_bin )
        {
          const arousals2_bin_metrics_t bm = arousals2_bin_compare(
            m == manual.end() ? std::vector<interval_t>() : m->second.stage_intervals,
            ae , me );
          writer.value( "BIN_SEC" , 30 );
          writer.value( "BIN_N" , (int)bm.n );
          writer.value( "BIN_TP" , (int)bm.tp );
          writer.value( "BIN_FP" , (int)bm.fp );
          writer.value( "BIN_TN" , (int)bm.tn );
          writer.value( "BIN_FN" , (int)bm.fn );
          writer.value( "BIN_PRECISION" , bm.precision );
          writer.value( "BIN_RECALL" , bm.recall );
          writer.value( "BIN_F1" , bm.f1 );
          writer.value( "BIN_KAPPA" , bm.kappa );
        }
      writer.value( "OVERLAP_SEC" , overlap_tp / (double)globals::tp_1sec );
      std::set<interval_t> all_events = ae;
      all_events.insert( me.begin() , me.end() );
      const uint64_t all_union_tp = arousals2_set_union_tp( all_events );
      writer.value( "UNION_ALL_SEC" , all_union_tp / (double)globals::tp_1sec );
      writer.value( "OVERLAP_UNION_ALL" , all_union_tp > 0 ? overlap_tp / (double)all_union_tp : 0.0 );
      writer.value( "UNION_SEC" , union_tp / (double)globals::tp_1sec );
      writer.value( "OVERLAP_UNION" , union_tp > 0 ? overlap_tp / (double)union_tp : 0.0 );
      if ( tp > 0 )
        {
          writer.value( "P_AUTO" , prop_auto / tp );
          writer.value( "P_MANUAL" , prop_manual / tp );
          writer.value( "START_D" , start_diff / tp );
          writer.value( "START_DABS" , start_diff_abs / tp );
          writer.value( "DUR_D" , dur_diff / tp );
        }
      writer.unlevel( "SS" );
      writer.unlevel( "COMP" );
    }
}

}

arousals2_t::arousals2_t( edf_t & edf , param_t & param )
{
  parent = &edf;

  const bool have_manual = param.has( "manual" );
  const std::string manual_label = have_manual ? param.requires( "manual" ) : "";
  const bool eval_by_bin = param.has( "eval-by-bin" ) &&
    ( param.empty( "eval-by-bin" ) || param.yesno( "eval-by-bin" , true , true ) );
  const double manual_min_overlap = param.has( "manual-min-overlap" ) ?
    param.requires_dbl( "manual-min-overlap" ) : 0.0;
  if ( manual_min_overlap < 0 ) Helper::halt( "AROUSALS manual-min-overlap must be non-negative" );

  const bool NO_ANNOTS = true;
  const std::string eeg_signal_label = param.has( "eeg" ) ? param.value( "eeg" ) : "";
  const std::string emg_signal_label = param.has( "emg" ) ? param.value( "emg" ) : "";
  signal_list_t eeg_signals = edf.header.signal_list( eeg_signal_label , NO_ANNOTS );
  signal_list_t emg_signals = edf.header.signal_list( emg_signal_label , NO_ANNOTS );
  const int ns_eeg = eeg_signals.size();
  const int ns_emg = emg_signals.size();
  std::map<std::string,arousals2_manual_stage_t> manual;
  if ( have_manual ) manual = arousals2_collect_manual( edf , manual_label );
  if ( have_manual && ns_eeg == 0 && ns_emg == 0 )
    {
      arousals2_write_manual( manual , std::map<std::string,std::set<interval_t> >() ,
                              (uint64_t)std::llround( manual_min_overlap * globals::tp_1sec ) ,
                              eval_by_bin );
      logger << "  no valid EEG or EMG signals; summarized manual annotation " << manual_label << "\n";
      return;
    }

  const double epoch_win = param.has( "win" ) ? param.requires_dbl( "win" ) : 4.0;
  const double epoch_inc = param.has( "inc" ) ? param.requires_dbl( "inc" ) : 0.5;
  if ( epoch_win < 1 || epoch_inc <= 0 || epoch_win < epoch_inc )
    Helper::halt( "invalid epoch win/inc values" );

  int ne = edf.timeline.set_epoch( epoch_win , epoch_inc );
  logger << "  deriving arousal features for " << ne << " "
         << epoch_win << "s epochs, w/ "
         << 100*( epoch_inc / epoch_win ) << "% overlap\n";

  uint64_t d = edf.header.record_duration_tp % globals::tp_1sec;
  if ( d != 0 )
    Helper::halt( "must have EDF records that are a multiple of 1 second (use RECORD-SIZE)" );

  if ( param.has( "winsor" ) || param.has( "no-winsor" ) )
    Helper::halt( "winsor/no-winsor options are not supported for arousal detection" );
  if ( param.has( "prefix" ) )
    Helper::halt( "prefix option is not supported; use add=<prefix> for derived channels" );
  if ( param.has( "per-channel" ) )
    Helper::halt( "per-channel option is not supported for arousal detection" );
  if ( param.has( "sigma-veto" ) || param.has( "sigma-veto2" ) ||
       param.has( "sigma-veto-rem" ) || param.has( "sigma-veto2-rem" ) )
    Helper::halt( "sigma veto options are not used for arousal detection" );

  arousals2_verbose_annots = false;
  arousals2_verbose_prefix = "arv_";
  const std::string verbose_key = param.has( "verbose-annot" ) ? "verbose-annot" :
    ( param.has( "verbose-annots" ) ? "verbose-annots" : "" );
  if ( verbose_key != "" )
    {
      if ( param.empty( verbose_key ) )
        arousals2_verbose_annots = true;
      else
        {
          const std::string v = param.value( verbose_key );
          if ( arousals2_bool_token( v ) )
            arousals2_verbose_annots = arousals2_true_token( v );
          else
            {
              arousals2_verbose_annots = true;
              arousals2_verbose_prefix = v;
            }
        }
    }

  const bool add_chs = param.has( "add" );
  const std::string ch_prefix = ! param.empty( "add" ) ? param.value( "add" ) : "a_";
  const std::string event_root = ! param.empty( "annot" ) ? param.value( "annot" ) : "arousal";

  arousals2_params_t hp;
  hp.epoch_inc = epoch_inc;
  hp.post_sec = epoch_win;
  hp.do_nrem = param.yesno( "nrem" , true , true );
  hp.do_rem  = param.yesno( "rem"  , true , true );
  hp.verbose = param.has( "verbose" );
  hp.nrem_emg_only = param.yesno( "nrem-emg-only" , false , true );
  if ( param.has( "eeg-tilt-rise-th" ) ) hp.eeg_tilt_rise_th = param.requires_dbl( "eeg-tilt-rise-th" );
  if ( param.has( "eeg-fast-rise-th" ) ) hp.eeg_fast_rise_th = param.requires_dbl( "eeg-fast-rise-th" );

  // Global EEG thresholds apply to both stages; stage-specific options below
  // can then override them independently.
  hp.eeg_tilt_rise_th_nrem = hp.eeg_tilt_rise_th;
  hp.eeg_tilt_rise_th_rem  = hp.eeg_tilt_rise_th;
  hp.eeg_fast_rise_th_nrem = hp.eeg_fast_rise_th;
  hp.eeg_fast_rise_th_rem  = hp.eeg_fast_rise_th;
  if ( param.has( "eeg-tilt-rise-th-nrem" ) )
    hp.eeg_tilt_rise_th_nrem = param.requires_dbl( "eeg-tilt-rise-th-nrem" );
  if ( param.has( "eeg-tilt-rise-th-rem" ) )
    hp.eeg_tilt_rise_th_rem = param.requires_dbl( "eeg-tilt-rise-th-rem" );
  if ( param.has( "eeg-fast-rise-th-nrem" ) )
    hp.eeg_fast_rise_th_nrem = param.requires_dbl( "eeg-fast-rise-th-nrem" );
  if ( param.has( "eeg-fast-rise-th-rem" ) )
    hp.eeg_fast_rise_th_rem = param.requires_dbl( "eeg-fast-rise-th-rem" );
  if ( param.has( "eeg-active-th" )    ) hp.eeg_active_th    = param.requires_dbl( "eeg-active-th" );
  if ( param.has( "emg-rise-th" )      ) hp.emg_rise_th      = param.requires_dbl( "emg-rise-th" );
  if ( param.has( "emg-active-th" )    ) hp.emg_active_th    = param.requires_dbl( "emg-active-th" );
  if ( param.has( "active-gap" )       ) hp.active_gap_sec   = param.requires_dbl( "active-gap" );
  if ( param.has( "delta-artifact-th" )) hp.delta_artifact_th = param.requires_dbl( "delta-artifact-th" );
  if ( param.has( "artifact-frac" )    ) hp.artifact_frac    = param.requires_dbl( "artifact-frac" );
  if ( param.has( "eeg-artifact-th" )  ) hp.eeg_artifact_th  = param.requires_dbl( "eeg-artifact-th" );
  if ( param.has( "emg-median" )       ) hp.emg_median_sec   = param.requires_dbl( "emg-median" );
  if ( param.has( "pre-local" )        ) hp.pre_sec          = param.requires_dbl( "pre-local" );
  if ( param.has( "pre-gap" )          ) hp.pre_gap_sec      = param.requires_dbl( "pre-gap" );
  if ( param.has( "min-dur" )          ) hp.min_dur          = param.requires_dbl( "min-dur" );
  if ( param.has( "max-dur" )          ) hp.max_dur          = param.requires_dbl( "max-dur" );
  if ( param.has( "long-dur" )         ) hp.long_dur         = param.requires_dbl( "long-dur" );
  if ( param.has( "arousal-dur" )      ) hp.arousal_dur      = param.requires_dbl( "arousal-dur" );
  if ( param.has( "duration-level-frac" ) ) hp.duration_level_frac = param.requires_dbl( "duration-level-frac" );
  if ( param.has( "merge-gap" )        ) hp.merge_gap_sec    = param.requires_dbl( "merge-gap" );
  if ( param.has( "pre-sleep" )        ) hp.pre_sleep_sec    = param.requires_dbl( "pre-sleep" );
  if ( param.has( "emg-rise-min-dur" ) ) hp.emg_rise_min_dur = param.requires_dbl( "emg-rise-min-dur" );
  if ( param.has( "emg-rise-buffer" )  ) hp.emg_rise_buffer  = param.requires_dbl( "emg-rise-buffer" );

  if ( ! std::isfinite( hp.eeg_tilt_rise_th ) || ! std::isfinite( hp.eeg_fast_rise_th ) ||
       ! std::isfinite( hp.eeg_tilt_rise_th_nrem ) || ! std::isfinite( hp.eeg_tilt_rise_th_rem ) ||
       ! std::isfinite( hp.eeg_fast_rise_th_nrem ) || ! std::isfinite( hp.eeg_fast_rise_th_rem ) ||
       ! std::isfinite( hp.eeg_active_th ) || ! std::isfinite( hp.emg_rise_th ) ||
       ! std::isfinite( hp.emg_active_th ) || ! std::isfinite( hp.active_gap_sec ) ||
       ! std::isfinite( hp.delta_artifact_th ) || ! std::isfinite( hp.artifact_frac ) ||
       ! std::isfinite( hp.eeg_artifact_th ) ||
       ! std::isfinite( hp.emg_median_sec ) || ! std::isfinite( hp.pre_sec ) ||
       ! std::isfinite( hp.pre_gap_sec ) || ! std::isfinite( hp.min_dur ) ||
       ! std::isfinite( hp.max_dur ) || ! std::isfinite( hp.long_dur ) ||
       ! std::isfinite( hp.arousal_dur ) || ! std::isfinite( hp.duration_level_frac ) ||
       ! std::isfinite( hp.merge_gap_sec ) ||
       ! std::isfinite( hp.pre_sleep_sec ) || ! std::isfinite( hp.emg_rise_min_dur ) ||
       ! std::isfinite( hp.emg_rise_buffer ) )
    Helper::halt( "invalid non-finite arousal-detection parameter" );
  if ( hp.emg_median_sec < 0 || hp.pre_sec <= hp.pre_gap_sec ||
       hp.pre_gap_sec < 0 || hp.active_gap_sec < 0 )
    Helper::halt( "invalid arousal-detection local window or smoothing parameter" );
  if ( hp.min_dur <= 0 || hp.arousal_dur < hp.min_dur ||
       hp.max_dur < hp.arousal_dur || hp.long_dur <= hp.max_dur )
    Helper::halt( "invalid arousal-detection duration parameter ordering" );
  if ( hp.duration_level_frac <= 0 || hp.duration_level_frac > 1 )
    Helper::halt( "duration-level-frac must be > 0 and <= 1" );
  if ( hp.merge_gap_sec < 0 || hp.pre_sleep_sec < 0 ||
       hp.emg_rise_min_dur < 0 || hp.emg_rise_buffer < 0 )
    Helper::halt( "invalid arousal-detection duration/gap parameter" );
  if ( hp.artifact_frac < 0 || hp.artifact_frac > 1 )
    Helper::halt( "artifact-frac must be between 0 and 1" );
  if ( hp.eeg_artifact_th < 0 || hp.eeg_artifact_th > 1 )
    Helper::halt( "eeg-artifact-th must be between 0 and 1" );

  const int wake_bridge = param.has( "wake-bridge" ) ? param.requires_int( "wake-bridge" ) : 1;
  if ( wake_bridge < 0 )
    Helper::halt( "wake-bridge must be 0 or greater" );

  if ( ns_eeg == 0 && ns_emg == 0 )
    {
      logger << "  no valid EEG or EMG signals detected... leaving arousal detection\n";
      return;
    }
  if ( hp.do_rem && ns_emg == 0 )
    Helper::halt( "REM arousal detection requires an EMG channel; use rem=F to run NREM-only" );

  logger << "  running arousal detection for " << ns_eeg << " EEG and " << ns_emg << " EMG signals\n";

  std::vector<double> Fs_eeg = edf.header.sampling_freq( eeg_signals );
  std::vector<double> Fs_emg = edf.header.sampling_freq( emg_signals );

  if ( ns_eeg )
    {
      const double Fs0 = Fs_eeg[0];
      for (int i=0; i<ns_eeg; i++)
        if ( fabs( Fs_eeg[i] - Fs0 ) > 1e-4 )
          Helper::halt( "all EEG signals must have similar sample rates" );

      const double sr_rounded = std::round( Fs0 );
      if ( fabs( Fs0 - sr_rounded ) > 1e-4 )
        Helper::halt( "arousal detection requires integer EEG sample rates" );
      sr = (int)sr_rounded;
      if ( sr < 60 )
        Helper::halt( "EEG sample rate too low" );
    }

  for (int i=1; i<ns_emg; i++)
    if ( fabs( Fs_emg[i] - Fs_emg[0] ) > 1e-4 )
      Helper::halt( "all EMG signals must have similar sample rates" );

  std::vector<int> state;
  std::vector<int> seq;
  std::vector<double> sec;
  Eigen::MatrixXd Xeeg, Xemg;
  build_ftr_matrix( edf , eeg_signals , emg_signals , Xeeg , Xemg ,
                    &state , &seq , &sec , wake_bridge , epoch_inc );

  Eigen::MatrixXd Xftr = process_ftr_matrix( &Xeeg , &Xemg , state , seq , sec , hp );
  Xeeg.resize(0,0);
  Xemg.resize(0,0);

  std::vector<std::vector<std::vector<double> > > tt;
  std::vector<std::vector<std::vector<Eigen::VectorXd> > > X = assemble( Xftr , state , seq , sec , &tt );

  std::vector<std::vector<std::vector<double> > > eeg_artifact;
  // Compute the raw-EEG quality score before event scoring.  add_channels()
  // also owns the diagnostic-channel emission, so it is reused here even
  // when no derived channels were requested.
  add_channels( X , tt , add_chs ? ch_prefix : "" , eeg_signals , epoch_win , &eeg_artifact );

  std::map<std::string,std::set<interval_t> > anns = event_heuristic( X , tt , eeg_artifact , hp );
  if ( event_root != "arousal" )
    {
      std::map<std::string,std::set<interval_t> > renamed;
      for ( std::map<std::string,std::set<interval_t> >::const_iterator a = anns.begin(); a != anns.end(); ++a )
        {
          const std::string prefix = "arousal";
          const std::string label = a->first.compare( 0 , prefix.size() , prefix ) == 0
            ? event_root + a->first.substr( prefix.size() ) : a->first;
          renamed[ label ].insert( a->second.begin() , a->second.end() );
        }
      anns.swap( renamed );
    }
  if ( have_manual )
    arousals2_write_manual( manual , anns ,
                            (uint64_t)std::llround( manual_min_overlap * globals::tp_1sec ) ,
                            eval_by_bin );
  std::map<std::string,std::set<interval_t> >::const_iterator aa = anns.begin();
  while ( aa != anns.end() )
    {
      // Keep the ordinary output compact: arousal is the canonical
      // all-stage, all-duration event class.  Attach the stage-specific,
      // standard/long, artifact, and verbose pipeline classes only when
      // diagnostic annotations have explicitly been requested.
      if ( ! arousals2_verbose_annots && aa->first != event_root )
        {
          ++aa;
          continue;
        }
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

}

void arousals2_t::build_ftr_matrix( edf_t & edf ,
                                    const signal_list_t & eeg_signals ,
                                    const signal_list_t & emg_signals ,
                                    Eigen::MatrixXd & Xeeg ,
                                    Eigen::MatrixXd & Xemg ,
                                    std::vector<int> * state ,
                                    std::vector<int> * sequence ,
                                    std::vector<double> * sec ,
                                    const int wake_bridge ,
                                    const double epoch_inc )
{
  const int ns_eeg = eeg_signals.size();
  const int ns_emg = emg_signals.size();

  state->clear();
  sequence->clear();
  sec->clear();

  int ne = edf.timeline.first_epoch();
  Xeeg = Eigen::MatrixXd::Zero( ne , ns_eeg * 4 );
  Xemg = Eigen::MatrixXd::Zero( ne , ns_emg );

  int eeg_index_length = 0;
  std::unique_ptr<FFT> eeg_fft;
  if ( ns_eeg )
    {
      eeg_index_length = std::lround( sr * edf.timeline.epoch_length() );
      if ( eeg_index_length < 1 )
        Helper::halt( "invalid EEG FFT window length" );
      eeg_fft.reset( new FFT( eeg_index_length , eeg_index_length , sr , FFT_FORWARD , WINDOW_TUKEY50 ) );
    }

  edf.annotations->make_sleep_stage( edf.timeline );
  annot_t * staging = edf.annotations->find( "SleepStage" );
  if ( staging == NULL )
    Helper::halt( "no staging information present" );

  int prior_state = 2;
  int seq_nr = -1;
  int seq_r = -1;
  bool have_prior_feature_sec = false;
  double prior_feature_sec = 0;

  const double wake_bridge_sec = wake_bridge * globals::default_epoch_len;
  int last_sleep_state = 2;
  bool in_wake_run = false;
  int wake_bridge_state = 2;
  double wake_run_start_sec = 0;

  const double emg_SD_threshold = 9;
  std::vector<double> emg_th( ns_emg );
  for (int s=0; s<ns_emg; s++)
    {
      slice_t slice( edf , emg_signals(s) , edf.timeline.wholetrace() );
      const std::vector<double> * d = slice.pdata();
      const int n = d->size();
      const double median = MiscMath::median( *d );
      std::vector<double> dd( n );
      for (int i=0; i<n; i++) dd[i] = fabs( (*d)[i] - median );
      const double MAD = MiscMath::median( dd );
      const double sigma = 1.4826 * MAD;
      if ( sigma < 1e-6 ) Helper::halt( "low EMG amplitude" );
      emg_th[s] = sigma * emg_SD_threshold;
    }

  int eidx = -1;
  while ( 1 )
    {
      int epoch = edf.timeline.next_epoch();
      if ( epoch == -1 ) break;
      ++eidx;

      interval_t interval = edf.timeline.epoch( epoch );
      const double feature_sec = interval.start_sec() + 0.5 * interval.duration_sec();
      const uint64_t feature_tp = (uint64_t)std::round( feature_sec * globals::tp_1sec );
      annot_map_t events = staging->extract( interval_t( feature_tp , feature_tp + 1LLU ) );

      int scored_st = 2;
      bool scored_wake = false;
      if ( events.size() >= 1 )
        {
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
              wake_run_start_sec = feature_sec;
            }
          if ( ( wake_bridge_state == 0 || wake_bridge_state == 1 ) &&
               feature_sec - wake_run_start_sec < wake_bridge_sec )
            st = wake_bridge_state;
        }
      else
        {
          last_sleep_state = 2;
          in_wake_run = false;
          wake_bridge_state = 2;
          wake_run_start_sec = 0;
        }

      const bool time_discontinuity =
        have_prior_feature_sec && feature_sec - prior_feature_sec > 1.5 * epoch_inc;
      if ( st != prior_state || time_discontinuity )
        {
          if ( st == 0 ) ++seq_nr;
          else if ( st == 1 ) ++seq_r;
        }

      prior_state = st;
      prior_feature_sec = feature_sec;
      have_prior_feature_sec = true;

      state->push_back( st );
      if ( st == 0 ) sequence->push_back( seq_nr );
      else if ( st == 1 ) sequence->push_back( seq_r );
      else sequence->push_back( -1 );
      sec->push_back( feature_sec );

      if ( st == 2 ) continue;

      if ( ns_eeg )
        {
          eigen_matslice_t mslice( edf , eeg_signals , interval );
          if ( mslice.data_ref().rows() != eeg_index_length )
            Helper::halt( "internal EEG epoch length mismatch in arousal detection" );
          Xeeg.row( eidx ) = calc_eeg_ftrs( mslice.data_ref() , *eeg_fft );
        }

      if ( ns_emg )
        {
          eigen_matslice_t mslice_emg( edf , emg_signals , interval );
          Xemg.row( eidx ) = calc_emg_ftrs( mslice_emg.data_ref() , emg_th );
        }
    }
}

Eigen::VectorXd arousals2_t::calc_eeg_ftrs( const Eigen::MatrixXd & X , FFT & fftseg )
{
  Eigen::VectorXd F = Eigen::VectorXd::Zero( 4 * X.cols() );
  const double eps = 1e-12;
  int fi = 0;

  for (int s=0; s<X.cols(); s++)
    {
      fftseg.apply( X.col(s).data() , X.rows() );

      double p_delta = 0;
      double p_theta = 0;
      double p_alpha = 0;
      double p_beta = 0;

      for (int f=0; f<fftseg.cutoff; f++)
        {
          const double frq = fftseg.frq[f];
          if ( frq < 1 || frq >= 30 ) continue;
          if      ( frq >= 16 ) p_beta  += fftseg.X[f];
          else if ( frq >= 10 ) continue;
          else if ( frq >= 8  ) p_alpha += fftseg.X[f];
          else if ( frq >= 4  ) p_theta += fftseg.X[f];
          else                  p_delta += fftseg.X[f];
        }

      F(fi++) = log( p_delta + eps );
      F(fi++) = log( p_theta + eps );
      F(fi++) = log( p_alpha + eps );
      F(fi++) = log( p_beta + eps );
    }

  return F;
}

Eigen::VectorXd arousals2_t::calc_emg_ftrs( const Eigen::MatrixXd & X , const std::vector<double> & thr )
{
  const int ns = X.cols();
  Eigen::VectorXd F = Eigen::VectorXd::Zero( ns );
  const double eps = 1e-8;

  for (int s=0; s<ns; s++)
    {
      const Eigen::VectorXd & d = X.col(s);
      const int n = d.size();
      const double t = thr[s];
      double rms = 0;
      for (int i=0; i<n; i++)
        {
          if ( d[i] < -t || d[i] > t ) rms += t*t;
          else rms += d[i] * d[i];
        }
      rms /= (double)n;
      F(s) = log( sqrt( rms ) + eps );
    }

  return F;
}

Eigen::MatrixXd arousals2_t::process_ftr_matrix( Eigen::MatrixXd * Xeeg ,
                                                 Eigen::MatrixXd * Xemg ,
                                                 const std::vector<int> & st ,
                                                 const std::vector<int> & seq ,
                                                 const std::vector<double> & sec ,
                                                 const arousals2_params_t & p )
{
  if ( Xeeg->rows() != Xemg->rows() )
    Helper::halt( "internal error in arousals2_t::process_ftr_matrix()" );
  if ( (int)sec.size() != Xeeg->rows() )
    Helper::halt( "internal time-track error in arousals2_t::process_ftr_matrix()" );
  if ( (int)seq.size() != Xeeg->rows() )
    Helper::halt( "internal sequence-track error in arousals2_t::process_ftr_matrix()" );

  const int n = Xeeg->rows();
  const int nch_eeg = Xeeg->cols() / 4;
  const int nch_emg = Xemg->cols();

  const double baseline_sec = 120.0;
  int baseline_n = p.epoch_inc > 0 ? (int)( baseline_sec / p.epoch_inc + 0.5 ) : 1;
  if ( baseline_n < 1 ) baseline_n = 1;

  auto subtract_segmented_median = [&]( Eigen::MatrixXd * M )
    {
      for (int s=0; s<M->cols(); s++)
        {
          int start = 0;
          while ( start < M->rows() )
            {
              while ( start < M->rows() && st[start] == 2 )
                ++start;
              if ( start >= M->rows() ) break;

              int stop = start;
              while ( stop + 1 < M->rows() &&
                      st[stop+1] == st[start] &&
                      seq[stop+1] == seq[start] )
                ++stop;
              const int nn = stop - start + 1;
              const int local_baseline_n = std::min( baseline_n , nn );
              Eigen::VectorXd x = M->col(s).segment( start , nn );
              M->col(s).segment( start , nn ) = x - eigen_ops::median_filter( x , local_baseline_n );
              start = stop + 1;
            }
        }
    };

  subtract_segmented_median( Xeeg );
  subtract_segmented_median( Xemg );

  Eigen::MatrixXd X = Eigen::MatrixXd::Zero( n , NFTR );
  for (int s=0; s<nch_eeg; s++)
    {
      X.col( IDX_DELTA ) += Xeeg->col( s*4 + 0 );
      X.col( IDX_THETA ) += Xeeg->col( s*4 + 1 );
      X.col( IDX_ALPHA ) += Xeeg->col( s*4 + 2 );
      X.col( IDX_BETA  ) += Xeeg->col( s*4 + 3 );
    }
  if ( nch_eeg > 1 )
    for (int c=IDX_DELTA; c<=IDX_BETA; c++)
      X.col(c) = X.col(c).array() / (double)nch_eeg;

  for (int s=0; s<nch_emg; s++)
    X.col( IDX_EMG ) += Xemg->col(s);
  if ( nch_emg > 1 )
    X.col( IDX_EMG ) = X.col( IDX_EMG ).array() / (double)nch_emg;

  for (int i=0; i<=IDX_EMG; i++)
    {
      const Eigen::VectorXd z_nrem = robust_mad_norm( X.col(i) , st , 0 );
      const Eigen::VectorXd z_rem  = robust_mad_norm( X.col(i) , st , 1 );
      for (int r=0; r<n; r++)
        {
          if      ( st[r] == 0 ) X(r,i) = z_nrem[r];
          else if ( st[r] == 1 ) X(r,i) = z_rem[r];
          else X(r,i) = 0;
          if ( X(r,i) < -12 ) X(r,i) = -12;
          if ( X(r,i) >  12 ) X(r,i) =  12;
        }
    }

  if ( nch_emg && p.emg_median_sec > 0 )
    {
      int mw = (int)( p.emg_median_sec / p.epoch_inc + 0.5 );
      if ( mw > 1 ) X.col( IDX_EMG ) = eigen_ops::median_filter( X.col( IDX_EMG ) , mw );
    }

  std::vector<double> delta(n), theta(n), alpha(n), beta(n), emg(n);
  for (int i=0; i<n; i++)
    {
      delta[i] = X(i,IDX_DELTA);
      theta[i] = X(i,IDX_THETA);
      alpha[i] = X(i,IDX_ALPHA);
      beta[i] = X(i,IDX_BETA);
      emg[i] = X(i,IDX_EMG);
    }

  for (int i=0; i<n; i++)
    {
      if ( st[i] != 0 && st[i] != 1 ) continue;

      // "current activity" window, centered on i: [ i - half_post , i + half_post ].
      // Using a centered (rather than purely forward [i, i+post_sec]) window removes
      // the systematic lead: the rise now increases when activity actually rises at i,
      // not post_sec seconds beforehand.
      const double half_post = 0.5 * p.post_sec;
      int cur0 = i, cur1 = i;
      while ( cur0 > 0 && st[cur0-1] == st[i] && seq[cur0-1] == seq[i] && sec[i] - sec[cur0-1] <= half_post )
        --cur0;
      while ( cur1 + 1 < n && st[cur1+1] == st[i] && seq[cur1+1] == seq[i] && sec[cur1+1] - sec[i] <= half_post )
        ++cur1;

      // backward baseline, ending before the centered window (gap >= half_post so they do not overlap)
      const double base_gap = p.pre_gap_sec > half_post ? p.pre_gap_sec : half_post;
      int pre0 = i, pre1 = i;
      while ( pre0 > 0 && st[pre0-1] == st[i] && seq[pre0-1] == seq[i] && sec[i] - sec[pre0-1] <= p.pre_sec )
        --pre0;
      while ( pre1 > pre0 && sec[i] - sec[pre1] < base_gap )
        --pre1;

      if ( pre1 <= pre0 ) continue;

      const double rd = local_mean( delta , cur0 , cur1 ) - local_median( delta , pre0 , pre1 );
      const double rt = local_mean( theta , cur0 , cur1 ) - local_median( theta , pre0 , pre1 );
      const double ra = local_mean( alpha , cur0 , cur1 ) - local_median( alpha , pre0 , pre1 );
      const double rb = local_mean( beta  , cur0 , cur1 ) - local_median( beta  , pre0 , pre1 );
      const double re = local_mean( emg   , cur0 , cur1 ) - local_median( emg   , pre0 , pre1 );

      const double fast_rise = std::max( rt , std::max( ra , rb ) );
      const double tilt_rise = std::max( rt - rd , std::max( ra - rd , rb - rd ) );

      X(i,IDX_EEG_TILT_RISE) = tilt_rise;
      X(i,IDX_EEG_FAST_RISE) = fast_rise;
      X(i,IDX_EMG_RISE) = re;
      X(i,IDX_DELTA_ART) = rd >= p.delta_artifact_th ? rd : 0;
    }

  return X;
}

double arousals2_t::local_median( const std::vector<double> & x , int a , int b ) const
{
  if ( a < 0 ) a = 0;
  if ( b >= (int)x.size() ) b = x.size() - 1;
  if ( b < a ) return 0;
  std::vector<double> y;
  y.reserve( b - a + 1 );
  for (int i=a; i<=b; i++) y.push_back( x[i] );
  return MiscMath::median( y );
}

double arousals2_t::local_mean( const std::vector<double> & x , int a , int b ) const
{
  if ( a < 0 ) a = 0;
  if ( b >= (int)x.size() ) b = x.size() - 1;
  if ( b < a ) return 0;
  double s = 0;
  for (int i=a; i<=b; i++) s += x[i];
  return s / (double)( b - a + 1 );
}

Eigen::VectorXd arousals2_t::robust_mad_norm( const Eigen::VectorXd & x ,
                                              const std::vector<int> & st ,
                                              const int target_state )
{
  std::vector<double> vals;
  for (int i=0; i<x.size(); i++)
    if ( st[i] == target_state && std::isfinite( x[i] ) )
      vals.push_back( x[i] );
  Eigen::VectorXd z = Eigen::VectorXd::Zero( x.size() );
  if ( vals.size() < 5 ) return z;

  const double median = MiscMath::median( vals );
  std::vector<double> dev( vals.size() );
  for (int i=0; i<vals.size(); i++) dev[i] = fabs( vals[i] - median );
  double sigma = 1.4826 * MiscMath::median( dev );
  if ( sigma < 1e-8 ) sigma = 1.0;
  for (int i=0; i<x.size(); i++) z[i] = ( x[i] - median ) / sigma;
  return z;
}

std::vector<std::vector<std::vector<Eigen::VectorXd> > >
arousals2_t::assemble( const Eigen::MatrixXd & Xftr ,
                       const std::vector<int> & state ,
                       const std::vector<int> & seq ,
                       const std::vector<double> & sec ,
                       std::vector<std::vector<std::vector<double> > > * tt )
{
  std::vector<std::vector<std::vector<Eigen::VectorXd> > > X;
  tt->clear();

  std::map<int,std::map<int,std::vector<int> > > contigs;
  for (int i=0; i<state.size(); i++)
    {
      if ( state[i] == 2 ) continue;
      contigs[ state[i] ][ seq[i] ].push_back( i );
    }

  X.resize( 2 );
  X[0].resize( contigs[0].size() );
  X[1].resize( contigs[1].size() );
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
          (*tt)[st][sq].resize( nce );
          X[st][sq].resize( nce );
          for (int i=0; i<nce; i++)
            {
              (*tt)[st][sq][i] = sec[ ee[i] ];
              X[st][sq][i] = Xftr.row( ee[i] );
            }
        }
    }

  return X;
}

void arousals2_t::add_channels( const std::vector<std::vector<std::vector<Eigen::VectorXd> > > & X ,
                                const std::vector<std::vector<std::vector<double> > > & tt ,
                                const std::string & ch_prefix ,
                                const signal_list_t & eeg_signals ,
                                const double epoch_win ,
                                std::vector<std::vector<std::vector<double> > > * eeg_artifact_out )
{
  const int nr = parent->header.nr;
  const double rs = parent->header.record_duration;
  const int nsamples = 2 * rs;

  std::vector<std::string> chlab = {
    "delta", "theta", "alpha", "beta", "emg",
    "eeg_tilt_rise", "eeg_fast_rise", "emg_rise", "delta_art", "activation", "eeg_artifact"
  };

  // Diagnostic-only bounded score for two kinds of non-EEG-like activity:
  // narrow-band/harmonic trains and long plateaus with abrupt edges.  This
  // never enters event detection.  It uses one EEG channel at a time; channels
  // are pooled by max below.
  auto artifact_score = [&]( const std::vector<double> & x , const int fs )
    {
      if ( x.size() < 20 || fs < 10 ) return 0.0;
      for ( int i=0; i<x.size(); ++i )
        if ( ! std::isfinite( x[i] ) ) return 0.0;
      const int step = std::max( 1 , (int)std::ceil( x.size() / 64.0 ) );
      std::vector<double> z;
      for ( int i=0; i<x.size(); i+=step ) z.push_back( x[i] );
      const double med = MiscMath::median( z );
      std::vector<double> ad( z.size() );
      for ( int i=0; i<z.size(); ++i ) ad[i] = std::fabs( z[i] - med );
      const double mad = MiscMath::median( ad );
      if ( mad < 1e-12 ) return 0.0;
      for ( int i=0; i<z.size(); ++i ) z[i] = ( z[i] - med ) / ( 1.4826 * mad );

      auto bound = [&]( const double v )
        { return std::max( 0.0 , std::min( 1.0 , v ) ); };

      std::vector<double> dx( z.size() - 1 );
      for ( int i=1; i<z.size(); ++i ) dx[i-1] = std::fabs( z[i] - z[i-1] );
      const double dsum = std::accumulate( dx.begin() , dx.end() , 0.0 );
      double jump_concentration = 0;
      if ( dsum > 1e-12 )
        {
          std::vector<double> sorted = dx;
          std::sort( sorted.begin() , sorted.end() );
          const int cut = std::max( 1 , (int)std::floor( 0.05 * sorted.size() ) );
          double top = 0;
          for ( int i=sorted.size()-cut; i<sorted.size(); ++i ) top += sorted[i];
          jump_concentration = std::max( 0.0 , std::min( 1.0 , ( top / dsum - 0.12 ) / 0.45 ) );
        }

      // Plateau/edge branch: long low-slope runs plus concentrated jumps.
      const double flat_fraction = bound( ( std::count_if( dx.begin() , dx.end(),
                                      [&]( const double v ) { return v < 0.20; } ) /
                                      (double)dx.size() - 0.45 ) / 0.45 );
      const double plateau_score = bound( flat_fraction * ( 0.35 + 0.65 * jump_concentration ) );

      double mean = std::accumulate( dx.begin() , dx.end() , 0.0 ) / dx.size();
      double var = 0;
      for ( int i=0; i<dx.size(); ++i ) var += ( dx[i] - mean ) * ( dx[i] - mean );
      var /= dx.size();
      const int lag0 = std::max( 1 , (int)std::lround( 0.20 * fs / step ) );
      const int lag1 = std::min( (int)dx.size() / 3 , (int)std::lround( 1.50 * fs / step ) );
      double periodicity = 0;
      if ( var >= 1e-12 )
        for ( int lag=lag0; lag<=lag1; ++lag )
          {
            double ac = 0;
            const int n = dx.size() - lag;
            for ( int i=0; i<n; ++i ) ac += ( dx[i] - mean ) * ( dx[i+lag] - mean );
            periodicity = std::max( periodicity , std::fabs( ac / ( n * var ) ) );
          }
      periodicity = std::max( 0.0 , std::min( 1.0 , ( periodicity - 0.25 ) / 0.50 ) );

      // Harmonic branch: use the raw waveform, since smooth periodic artifacts
      // need not have sharp transitions.  A narrow spectral line is stronger
      // evidence when it has additional harmonics above the local background.
      std::vector<double> y( x.size() );
      double ymean = 0;
      for ( int i=0; i<x.size(); ++i ) ymean += x[i];
      ymean /= x.size();
      for ( int i=0; i<x.size(); ++i ) y[i] = x[i] - ymean;
      FFT fft( y.size() , y.size() , fs , FFT_FORWARD , WINDOW_HANN );
      fft.apply( y );
      std::vector<double> power;
      std::vector<double> freq;
      for ( int i=0; i<fft.cutoff; ++i )
        if ( fft.frq[i] >= 0.5 && fft.frq[i] <= 40.0 && std::isfinite( fft.X[i] ) )
          {
            freq.push_back( fft.frq[i] );
            power.push_back( std::max( 0.0 , fft.X[i] ) );
          }
      double harmonic_score = 0;
      if ( power.size() >= 8 )
        {
          std::vector<double> residual( power.size() , 0 );
          for ( int i=0; i<power.size(); ++i )
            {
              std::vector<double> local;
              for ( int j=std::max( 0 , i-5 ); j<=std::min( (int)power.size()-1 , i+5 ); ++j )
                if ( j != i ) local.push_back( power[j] );
              const double bg = MiscMath::median( local );
              residual[i] = std::log( ( power[i] + 1e-12 ) / ( bg + 1e-12 ) );
            }
          double best = 0;
          for ( int i=1; i+1<residual.size(); ++i )
            if ( residual[i] > residual[i-1] && residual[i] >= residual[i+1] )
              {
                const double prominence = bound( ( residual[i] - 0.8 ) / 2.0 );
                if ( prominence == 0 ) continue;
                int nh = 0;
                for ( int k=2; k<=5; ++k )
                  {
                    const double target = k * freq[i];
                    if ( target > 40 ) break;
                    int jbest = -1;
                    double dbest = 1e9;
                    for ( int j=1; j<residual.size(); ++j )
                      if ( std::fabs( freq[j] - target ) < dbest )
                        { dbest = std::fabs( freq[j] - target ); jbest = j; }
                    if ( jbest >= 0 && dbest <= 0.5 && residual[jbest] > 0.8 ) ++nh;
                  }
                best = std::max( best , prominence * ( nh >= 1 ? 1.0 : 0.35 ) );
              }
          harmonic_score = bound( best * ( 0.45 + 0.55 * periodicity ) );
        }

      // Piecewise-linearity branch: fit short local straight lines.  This
      // catches held/interpolated-looking segments even when there is no
      // sharp edge or narrow spectral peak.
      double straight_score = 0;
      const int nsub = 8;
      const int msub = z.size() / nsub;
      if ( msub >= 4 )
        {
          int good = 0, longest = 0, run = 0;
          for ( int s=0; s<nsub; ++s )
            {
              const int a = s * msub;
              double st = 0, sy = 0, sty = 0, stt = 0, syy = 0;
              for ( int j=0; j<msub; ++j )
                {
                  const double t = j;
                  const double v = z[a+j];
                  st += t; sy += v; sty += t*v; stt += t*t; syy += v*v;
                }
              const double den = msub * stt - st * st;
              const double slope = std::fabs( den ) > 1e-12 ? ( msub * sty - st * sy ) / den : 0;
              const double intercept = ( sy - slope * st ) / msub;
              double rss = 0, tss = 0;
              const double ybar = sy / msub;
              for ( int j=0; j<msub; ++j )
                {
                  const double e = z[a+j] - ( intercept + slope*j );
                  rss += e*e;
                  const double d = z[a+j] - ybar;
                  tss += d*d;
                }
              const double rms = std::sqrt( rss / msub );
              const double r2 = tss > 1e-12 ? std::max( 0.0 , 1.0 - rss / tss ) : 1.0;
              // Either a locally very good line, or a genuinely held/flat
              // segment with little within-segment variation.
              const bool is_straight = rms < 0.40 && ( r2 > 0.80 || tss / msub < 0.10 );
              if ( is_straight ) { ++good; ++run; longest = std::max( longest , run ); }
              else run = 0;
            }
          const double fraction = (double)good / nsub;
          const double persistence = (double)longest / nsub;
          const double evidence = bound( ( fraction - 0.25 ) / 0.50 ) *
            bound( ( persistence - 0.125 ) / 0.375 );
          straight_score = 1.0 - std::exp( -4.0 * evidence );
        }

      // The branches are deliberately independent: a smooth harmonic train,
      // an edge/plateau artifact, and a piecewise-linear artifact need not
      // share the same signature.
      const double score = std::max( straight_score , std::max( harmonic_score , plateau_score ) );
      return std::isfinite( score ) ? std::max( 0.0 , std::min( 1.0 , score ) ) : 0.0;
    };

  std::vector<std::vector<std::vector<double> > > eeg_artifact( tt.size() );
  for ( int st=0; st<tt.size(); ++st )
    {
      eeg_artifact[st].resize( tt[st].size() );
      for ( int sq=0; sq<tt[st].size(); ++sq )
        eeg_artifact[st][sq].assign( tt[st][sq].size() , 0 );
    }
  if ( eeg_signals.size() )
    {
      const int fs = (int)std::lround( parent->header.sampling_freq( eeg_signals(0) ) );
      std::vector<const std::vector<double> *> raw( eeg_signals.size() );
      std::vector<const std::vector<uint64_t> *> raw_tp( eeg_signals.size() );
      std::vector<slice_t *> slices;
      for ( int s=0; s<eeg_signals.size(); ++s )
        {
          slice_t * sl = new slice_t( *parent , eeg_signals(s) , parent->timeline.wholetrace() );
          slices.push_back( sl );
          raw[s] = sl->pdata();
          raw_tp[s] = sl->ptimepoints();
        }
      for ( int st=0; st<tt.size(); ++st )
        for ( int sq=0; sq<tt[st].size(); ++sq )
          {
            const int ne = tt[st][sq].size();
            for ( int i=0; i<ne; ++i )
              {
                const uint64_t center = (uint64_t)std::llround( tt[st][sq][i] * globals::tp_1sec );
                const uint64_t half = (uint64_t)std::llround( 0.5 * epoch_win * globals::tp_1sec );
                double pooled = 0;
                for ( int s=0; s<eeg_signals.size(); ++s )
                  {
                    const std::vector<uint64_t> & tpv = *raw_tp[s];
                    std::vector<uint64_t>::const_iterator lo = std::lower_bound( tpv.begin() , tpv.end() , center > half ? center-half : 0 );
                    std::vector<uint64_t>::const_iterator hi = std::upper_bound( tpv.begin() , tpv.end() , center+half );
                    if ( lo == hi ) continue;
                    std::vector<double> w( raw[s]->begin() + ( lo - tpv.begin() ) , raw[s]->begin() + ( hi - tpv.begin() ) );
                    const double score = artifact_score( w , fs );
                    if ( std::isfinite( score ) ) pooled = std::max( pooled , score );
                  }
                eeg_artifact[st][sq][i] = pooled;
              }
            // The diagnostic is intended to identify sustained contamination,
            // not isolated feature-window excursions.  Smooth the pooled
            // per-channel score over five 2-Hz points (about 2.5 seconds).
            // This reduces single blips while retaining artifacts lasting
            // several seconds, including events shorter than 10 seconds.
            const std::vector<double> raw_score = eeg_artifact[st][sq];
            for ( int i=0; i<ne; ++i )
              {
                const int lo = std::max( 0 , i - 2 );
                const int hi = std::min( ne - 1 , i + 2 );
                double sum = 0;
                for ( int j=lo; j<=hi; ++j ) sum += raw_score[j];
                eeg_artifact[st][sq][i] = std::max( 0.0 , std::min( 1.0 , sum / ( hi - lo + 1 ) ) );
              }
          }
      for ( int s=0; s<slices.size(); ++s ) delete slices[s];
    }

  if ( eeg_artifact_out != NULL ) *eeg_artifact_out = eeg_artifact;

  if ( ch_prefix.empty() ) return;

  for (int cidx=0; cidx<chlab.size(); cidx++)
    {
      std::string label = ch_prefix;
      if ( ! label.empty() && label[ label.size() - 1 ] != '_' ) label += "_";
      label += chlab[cidx];
      parent->init_signal( label , 2 );

      const int slot = parent->header.signal( label );
      slice_t slice( *parent , slot , parent->timeline.wholetrace() );
      const std::vector<uint64_t> * tp = slice.ptimepoints();
      const int np = tp->size();
      if ( np != nr * nsamples )
        Helper::halt( "internal error in arousals2_t::add_channels()" );

      std::vector<double> sec( np );
      for (int i=0; i<np; i++)
        sec[i] = (*tp)[i] / (double)globals::tp_1sec;

      std::vector<double> xx( np , 0 );
      for (int st=0; st<tt.size(); st++)
        for (int sq=0; sq<tt[st].size(); sq++)
          {
            const int ne = tt[st][sq].size();
            if ( ne == 0 ) continue;

            std::vector<double> xx0(ne);
            for (int i=0; i<ne; i++)
              {
                if ( cidx < NFTR )
                  xx0[i] = X[st][sq][i][cidx];
                else if ( cidx == NFTR + 0 )
                  xx0[i] = std::max( X[st][sq][i](IDX_THETA),
                                     std::max( X[st][sq][i](IDX_ALPHA), X[st][sq][i](IDX_BETA) ) );
                else
                  xx0[i] = eeg_artifact[st][sq][i];
                if ( ! std::isfinite( xx0[i] ) ) xx0[i] = 0;
              }

            const double tmin = tt[st][sq][0];
            const double tmax = tt[st][sq][ne-1];

            tk::spline spline;
            if ( ne > 2 ) spline.set_points( tt[st][sq] , xx0 );

            for (int i=0; i<np; i++)
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
          }

      for ( int i=0; i<np; ++i )
        if ( ! std::isfinite( xx[i] ) ) xx[i] = 0;
      parent->update_signal( slot , &xx );
    }
}

std::map<std::string,std::set<interval_t> >
arousals2_t::event_heuristic( const std::vector<std::vector<std::vector<Eigen::VectorXd> > > & X ,
                              const std::vector<std::vector<std::vector<double> > > & tt ,
                              const std::vector<std::vector<std::vector<double> > > & eeg_artifact ,
                              const arousals2_params_t & p )
{
  std::map<std::string,std::set<interval_t> > ret;

  const int max_gap_samples = std::max( 1 , (int)std::lround( p.merge_gap_sec / p.epoch_inc ) );
  const int active_gap_samples = std::max( 0 , (int)std::lround( p.active_gap_sec / p.epoch_inc ) );
  const bool do_long = p.long_dur > p.max_dur;
  const double candidate_max_dur = do_long ? p.long_dur : p.max_dur;

  for (int st=0; st<2; st++)
    {
      if ( st == 0 && ! p.do_nrem ) continue;
      if ( st == 1 && ! p.do_rem  ) continue;

      const bool is_rem = st == 1;
      const std::string stg_lab = st == 0 ? "nrem" : "rem";

      Eigen::VectorXd ftr_art( NFTR ), ftr_nonart( NFTR );
      Eigen::VectorXd ftr_baseline( NFTR ), ftr_arousal( NFTR );
      ftr_art.setZero(); ftr_nonart.setZero();
      ftr_baseline.setZero(); ftr_arousal.setZero();
      int n_art = 0, n_nonart = 0, n_baseline = 0, n_arousal = 0;

      int cnt_evts = 0, cnt_long_evts = 0, cnt_arts = 0;
      int cnt_rem_cortical_only = 0, cnt_long_rem_cortical_only = 0, cnt_rem_emg_confirmed_long = 0;
      double dur_major = 0, dur_long = 0;

      int vcand_raw = 0, vcand_eeg_raw = 0, vcand_emg_raw = 0, vcand_dur_ok = 0;
      int vcand_dur_too_short = 0;
      int vcand_merged = 0, vcand_final_dur_ok = 0, vcand_presleep_ok = 0;
      int vcand_final_dur_too_short = 0, vcand_final_dur_trimmed = 0;
      int vcand_artifact_block = 0, vcand_rem_emg_confirmed = 0, vcand_rem_emg_rejected = 0;
      int vdiag_win = 0, vdiag_eeg_seed = 0, vdiag_emg_seed = 0, vdiag_delta_art = 0;

      const int nc = X[st].size();
      for (int c=0; c<nc; c++)
        {
          const std::vector<Eigen::VectorXd> & D = X[st][c];
          const int ne = D.size();
          if ( ne == 0 ) continue;

	          auto add_verbose_interval = [&]( const std::string & label , const double t0 , const double t1 )
	            {
              if ( ! arousals2_verbose_annots || t1 <= t0 ) return;
              ret[ arousals2_verbose_label( label ) ].insert(
	                interval_t( globals::tp_1sec * t0 , globals::tp_1sec * t1 ) );
            };

	          auto event_interval = [&]( const std::pair<int,int> & evt )
	            {
	              const double half_inc = 0.5 * p.epoch_inc;
	              const double t0 = tt[st][c][ evt.first ] - half_inc;
	              const double t1 = tt[st][c][ evt.second ] + half_inc;
	              return std::make_pair( t0 , t1 );
            };

          auto add_verbose_mask = [&]( const std::string & label , const std::vector<bool> & mask )
            {
              if ( ! arousals2_verbose_annots ) return;
	              const std::vector<std::pair<int,int> > ints = mask_to_intervals( mask );
	              for (int ii=0; ii<ints.size(); ii++)
                {
                  const std::pair<double,double> tint = event_interval( ints[ii] );
                  add_verbose_interval( label , tint.first , tint.second );
                }
            };

          auto add_verbose_events = [&]( const std::string & label ,
                                          const std::vector<std::pair<int,int> > & ev )
            {
              if ( ! arousals2_verbose_annots ) return;
	              for (int ii=0; ii<ev.size(); ii++)
                {
                  const std::pair<double,double> tint = event_interval( ev[ii] );
                  add_verbose_interval( label , tint.first , tint.second );
                }
            };

          const int candidate_max_samples =
            std::max( 1 , (int)std::floor( candidate_max_dur / p.epoch_inc + 1e-9 ) );

          auto trim_to_candidate_max = [&]( std::pair<int,int> evt )
            {
              if ( evt.second - evt.first + 1 > candidate_max_samples )
                evt.second = std::min( ne - 1 , evt.first + candidate_max_samples - 1 );
              return evt;
            };

          auto require_seed_support = [&]( const std::vector<bool> & seed ,
                                           const double support_sec ,
                                           const double min_seed_sec )
            {
              std::vector<bool> supported( ne , false );
              const int support_samples = std::max( 1 , (int)std::lround( support_sec / p.epoch_inc ) );
              const int min_seed_samples = std::max( 1 , (int)std::ceil( min_seed_sec / p.epoch_inc - 1e-9 ) );

              for (int i=0; i<ne; i++)
                {
                  if ( ! seed[i] ) continue;
                  int nseed = 0;
                  const int stop = std::min( ne - 1 , i + support_samples - 1 );
                  for (int j=i; j<=stop; j++)
                    if ( seed[j] ) ++nseed;
                  supported[i] = nseed >= min_seed_samples;
                }

              return supported;
            };

          std::vector<bool> eeg_seed( ne , false );
          std::vector<bool> emg_seed( ne , false );
          std::vector<bool> eeg_tilt_rise_pass( ne , false );
          std::vector<bool> eeg_fast_rise_pass( ne , false );
          std::vector<bool> eeg_active( ne , false );
          std::vector<bool> emg_active( ne , false );
          std::vector<bool> artifact( ne , false );
          std::vector<bool> delta_artifact( ne , false );
          std::vector<bool> raw_eeg_artifact( ne , false );
          std::vector<bool> artifact_blocked( ne , false );

          for (int i=0; i<ne; i++)
            {
              const Eigen::VectorXd & f = D[i];
              const double eeg_tilt_rise_th = is_rem ? p.eeg_tilt_rise_th_rem : p.eeg_tilt_rise_th_nrem;
              const double eeg_fast_rise_th = is_rem ? p.eeg_fast_rise_th_rem : p.eeg_fast_rise_th_nrem;
              eeg_tilt_rise_pass[i] = f(IDX_EEG_TILT_RISE) >= eeg_tilt_rise_th;
              eeg_fast_rise_pass[i] = f(IDX_EEG_FAST_RISE) >= eeg_fast_rise_th;
              eeg_seed[i] = eeg_tilt_rise_pass[i] && eeg_fast_rise_pass[i];
              emg_seed[i] = f(IDX_EMG_RISE) >= p.emg_rise_th;
              eeg_active[i] = f(IDX_EEG_FAST_RISE) >= p.eeg_active_th &&
                f(IDX_EEG_TILT_RISE) > 0;
              emg_active[i] = f(IDX_EMG_RISE) >= p.emg_active_th;
              delta_artifact[i] = f(IDX_DELTA_ART) > 0;
              raw_eeg_artifact[i] =
                st < eeg_artifact.size() && c < eeg_artifact[st].size() &&
                i < eeg_artifact[st][c].size() &&
                eeg_artifact[st][c][i] >= p.eeg_artifact_th;
              artifact[i] = delta_artifact[i] || raw_eeg_artifact[i];

              ++vdiag_win;
              if ( emg_seed[i] ) ++vdiag_emg_seed;
              if ( artifact[i] ) ++vdiag_delta_art;
            }

          const std::vector<bool> eeg_seed_supported = require_seed_support( eeg_seed , 2.0 , 1.0 );

          for (int i=0; i<ne; i++)
            if ( eeg_seed_supported[i] ) ++vdiag_eeg_seed;

          add_verbose_mask( "eeg_seed_supported" , eeg_seed_supported );
          add_verbose_mask( "eeg_tilt_rise_pass" , eeg_tilt_rise_pass );
          add_verbose_mask( "eeg_fast_rise_pass" , eeg_fast_rise_pass );
          add_verbose_mask( "eeg_fast_active" , eeg_active );
          add_verbose_mask( "emg_rise_seed" , emg_seed );
          add_verbose_mask( "delta_artifact" , delta_artifact );
          add_verbose_mask( "eeg_artifact" , raw_eeg_artifact );

          for (int i=0; i<ne; i++)
            {
              if ( artifact[i] ) { ftr_art += D[i]; ++n_art; }
              else { ftr_nonart += D[i]; ++n_nonart; }
            }

          std::vector<double> delta(ne), theta(ne), alpha(ne), beta(ne);
          for (int i=0; i<ne; i++)
            {
              delta[i] = D[i](IDX_DELTA);
              theta[i] = D[i](IDX_THETA);
              alpha[i] = D[i](IDX_ALPHA);
              beta[i]  = D[i](IDX_BETA);
            }

          auto grow_events = [&]( const std::vector<bool> & seed ,
                                  const std::vector<bool> & active )
            {
              std::vector<std::pair<int,int> > ev;
              std::vector<bool> used( ne , false );
              const int onset_backpad_samples = std::max( 0 , (int)std::lround( 0.5 / p.epoch_inc ) );
              for (int i=0; i<ne; i++)
                {
                  if ( ! seed[i] || used[i] ) continue;
                  int start = i, stop = i;
                  int gaps = 0;
                  while ( start > 0 )
                    {
                      if ( i - ( start - 1 ) >= candidate_max_samples ) break;
                      if ( used[start-1] ) break;
                      if ( active[start-1] ) gaps = 0;
                      else if ( ++gaps > active_gap_samples ) break;
                      --start;
                    }
                  gaps = 0;
                  while ( stop + 1 < ne )
                    {
                      if ( ( stop + 1 ) - start + 1 > candidate_max_samples ) break;
                      if ( active[stop+1] ) gaps = 0;
                      else if ( ++gaps > active_gap_samples ) break;
                      ++stop;
                    }
                  while ( start < i && ! active[start] ) ++start;
                  while ( stop > i && ! active[stop] ) --stop;
                  if ( stop >= start )
                    {
                      const int support_start = start;
                      const int emit_start = std::max( start , i - onset_backpad_samples );
                      std::pair<int,int> evt = trim_to_candidate_max( std::make_pair( emit_start , stop ) );
                      for (int j=support_start; j<=stop; j++) used[j] = true;
                      ev.push_back( evt );
                    }
                }
              return ev;
            };

          // EEG candidate growth is anchored to the baseline used at the
          // initial supported seed.  Seed detection above remains based on
          // the per-epoch rise features, but once a seed starts an event we
          // do not let the backward baseline adapt as the event grows.
          //
          // For a seed at i, the fixed baseline is the preceding pre_sec
          // seconds, ending before the base gap.  At each candidate epoch j,
          // the current activity window remains centered on j; only the
          // baseline is fixed to i.
          auto grow_eeg_events = [&]( const std::vector<bool> & seed )
            {
              std::vector<std::pair<int,int> > ev;
              std::vector<bool> used( ne , false );
              const int onset_backpad_samples = std::max( 0 , (int)std::lround( 0.5 / p.epoch_inc ) );

              const double half_post = 0.5 * p.post_sec;
              const double base_gap = p.pre_gap_sec > half_post ? p.pre_gap_sec : half_post;

              auto baseline_bounds = [&]( const int anchor , int * pre0 , int * pre1 )
                {
                  *pre0 = anchor;
                  *pre1 = anchor;
                  while ( *pre0 > 0 &&
                          tt[st][c][anchor] - tt[st][c][*pre0-1] <= p.pre_sec )
                    --(*pre0);
                  while ( *pre1 > *pre0 && tt[st][c][anchor] - tt[st][c][*pre1] < base_gap )
                    --(*pre1);
                  return *pre1 > *pre0;
                };

              auto current_bounds = [&]( const int j , int * cur0 , int * cur1 )
                {
                  *cur0 = j;
                  *cur1 = j;
                  while ( *cur0 > 0 &&
                          tt[st][c][j] - tt[st][c][*cur0-1] <= half_post )
                    --(*cur0);
                  while ( *cur1 + 1 < ne &&
                          tt[st][c][*cur1+1] - tt[st][c][j] <= half_post )
                    ++(*cur1);
                };

              for (int i=0; i<ne; i++)
                {
                  if ( ! seed[i] || used[i] ) continue;

                  int pre0 = i, pre1 = i;
                  if ( ! baseline_bounds( i , &pre0 , &pre1 ) ) continue;

                  const double bd = local_median( delta , pre0 , pre1 );
                  const double bt = local_median( theta , pre0 , pre1 );
                  const double ba = local_median( alpha , pre0 , pre1 );
                  const double bb = local_median( beta  , pre0 , pre1 );

                  auto active = [&]( const int j )
                    {
                      int cur0 = j, cur1 = j;
                      current_bounds( j , &cur0 , &cur1 );
                      const double rd = local_mean( delta , cur0 , cur1 ) - bd;
                      const double rt = local_mean( theta , cur0 , cur1 ) - bt;
                      const double ra = local_mean( alpha , cur0 , cur1 ) - ba;
                      const double rb = local_mean( beta  , cur0 , cur1 ) - bb;
                      const double fast_rise = std::max( rt , std::max( ra , rb ) );
                      const double tilt_rise = std::max( rt - rd , std::max( ra - rd , rb - rd ) );
                      return fast_rise >= p.eeg_active_th && tilt_rise > 0;
                    };

                  int start = i, stop = i;
                  int gaps = 0;
                  while ( start > 0 )
                    {
                      if ( i - ( start - 1 ) >= candidate_max_samples ) break;
                      if ( used[start-1] ) break;
                      if ( active( start - 1 ) ) gaps = 0;
                      else if ( ++gaps > active_gap_samples ) break;
                      --start;
                    }
                  gaps = 0;
                  while ( stop + 1 < ne )
                    {
                      if ( ( stop + 1 ) - start + 1 > candidate_max_samples ) break;
                      if ( active( stop + 1 ) ) gaps = 0;
                      else if ( ++gaps > active_gap_samples ) break;
                      ++stop;
                    }
                  while ( start < i && ! active( start ) ) ++start;
                  while ( stop > i && ! active( stop ) ) --stop;
                  if ( stop >= start )
                    {
                      const int support_start = start;
                      const int emit_start = std::max( start , i - onset_backpad_samples );
                      std::pair<int,int> evt = trim_to_candidate_max( std::make_pair( emit_start , stop ) );
                      for (int j=support_start; j<=stop; j++) used[j] = true;
                      ev.push_back( evt );
                    }
                }
              return ev;
            };

          std::vector<std::pair<int,int> > eeg_evts = grow_eeg_events( eeg_seed_supported );
          std::vector<std::pair<int,int> > emg_evts;
          const bool need_emg_evts = p.verbose || arousals2_verbose_annots || ( ! is_rem && p.nrem_emg_only );
          if ( need_emg_evts )
            emg_evts = grow_events( emg_seed , emg_active );
          vcand_eeg_raw += eeg_evts.size();
          vcand_emg_raw += emg_evts.size();
          add_verbose_events( "candidate_eeg_raw" , eeg_evts );
          add_verbose_events( "candidate_emg_raw" , emg_evts );

          // Scored arousals are EEG-seeded by default. EMG-only NREM scoring is
          // available as an explicit opt-in for exploratory use; REM always
          // requires EEG candidates with independent EMG confirmation.
          std::vector<std::pair<int,int> > evts = eeg_evts;
          if ( ! is_rem && p.nrem_emg_only )
            evts.insert( evts.end() , emg_evts.begin() , emg_evts.end() );
          std::sort( evts.begin() , evts.end() );
          vcand_raw += evts.size();
          add_verbose_events( "candidate_raw" , evts );

          std::vector<std::pair<int,int> > evts2;
          for (int e=0; e<evts.size(); e++)
            {
              std::pair<int,int> evt = evts[e];
              const double duration_sec = ( evt.second - evt.first + 1 ) * p.epoch_inc;
              if ( duration_sec < p.min_dur )
                {
                  ++vcand_dur_too_short;
                  continue;
                }
	              evts2.push_back( evt );
            }
          vcand_dur_ok += evts2.size();
          add_verbose_events( "candidate_duration_ok" , evts2 );

          evts = merge_events_with_gap_sorted( evts2 , max_gap_samples );
          vcand_merged += evts.size();
          add_verbose_events( "candidate_merged" , evts );

          evts2.clear();
          for (int e=0; e<evts.size(); e++)
            {
              std::pair<int,int> evt = evts[e];
              const double duration_sec = ( evt.second - evt.first + 1 ) * p.epoch_inc;
              if ( duration_sec < p.min_dur )
                {
                  ++vcand_final_dur_too_short;
                  continue;
                }
              if ( duration_sec > candidate_max_dur )
                {
                  evt = trim_to_candidate_max( evt );
                  ++vcand_final_dur_trimmed;
                }
              evts2.push_back( evt );
          }
          vcand_final_dur_ok += evts2.size();
          add_verbose_events( "candidate_duration_capped" , evts2 );

          evts.clear();
          for (int e=0; e<evts2.size(); e++)
            {
              const double t_event_start = tt[st][c][ evts2[e].first ];
              const double t_contig_start = tt[st][c][0];
              if ( t_event_start - t_contig_start >= p.pre_sleep_sec )
                evts.push_back( evts2[e] );
            }
          vcand_presleep_ok += evts.size();
          add_verbose_events( "candidate_presleep_ok" , evts );

          evts2.clear();
          for (int e=0; e<evts.size(); e++)
            {
              int n_art_evt = 0;
              for (int i=evts[e].first; i<=evts[e].second; i++)
                if ( artifact[i] ) ++n_art_evt;
              const double frac = n_art_evt / (double)( evts[e].second - evts[e].first + 1 );
              // block only when a real fraction of the event is extreme-movement artifact;
              // the onset-epoch veto was dropped (the tilt-rise seed already rejects
              // movement-driven onsets, so it only deleted genuine K-complex-onset arousals)
              if ( frac >= p.artifact_frac )
                {
                  ++vcand_artifact_block;
                  for (int i=evts[e].first; i<=evts[e].second; i++)
                    artifact_blocked[i] = true;
                  const std::pair<double,double> tint = event_interval( evts[e] );
                  const double t0 = tint.first;
                  const double t1 = tint.second;
                  ret[ "suppressed_" + stg_lab ].insert( interval_t( globals::tp_1sec * t0 , globals::tp_1sec * t1 ) );
                }
              else
                evts2.push_back( evts[e] );
            }
          evts = evts2;
          add_verbose_events( "candidate_artifact_ok" , evts );

          if ( is_rem )
            {
              std::vector<bool> emg_high( ne , false );
              for (int i=0; i<ne; i++) emg_high[i] = emg_active[i] || emg_seed[i];
              const std::vector<std::pair<int,int> > emg_runs = mask_to_intervals( emg_high );
              add_verbose_mask( "rem_emg_rise" , emg_high );

              std::vector<std::pair<double,double> > emg_windows;
	              for (int r=0; r<emg_runs.size(); r++)
                {
                  const std::pair<double,double> eint = event_interval( emg_runs[r] );
                  const double e0 = eint.first;
                  const double e1 = eint.second;
                  if ( e1 - e0 >= p.emg_rise_min_dur )
                    {
                      emg_windows.push_back( std::make_pair( e0 , e1 ) );
                      add_verbose_interval( "rem_emg_rise_duration_ok" , e0 , e1 );
                    }
                }

              std::vector<std::pair<int,int> > confirmed, rejected;
	              for (int e=0; e<evts.size(); e++)
                {
                  const std::pair<double,double> tint = event_interval( evts[e] );
                  const double t0 = tint.first;
                  const double t1 = tint.second;
                  bool ok = false;
                  for (int w=0; w<emg_windows.size(); w++)
                    if ( emg_windows[w].first <= t1 + p.emg_rise_buffer &&
                         emg_windows[w].second >= t0 - p.emg_rise_buffer )
                      { ok = true; break; }
                  if ( ok ) confirmed.push_back( evts[e] );
                  else rejected.push_back( evts[e] );
                }
              vcand_rem_emg_confirmed += confirmed.size();
              vcand_rem_emg_rejected += rejected.size();
              add_verbose_events( "rem_emg_confirmed" , confirmed );
              add_verbose_events( "rem_emg_rejected" , rejected );

	              for (int e=0; e<rejected.size(); e++)
                {
                  const std::pair<double,double> tint = event_interval( rejected[e] );
                  const double t0 = tint.first;
                  const double t1 = tint.second;
                  const double dur = t1 - t0;
                  ret[ "suppressed_rem" ].insert( interval_t( globals::tp_1sec * t0 , globals::tp_1sec * t1 ) );
                  if ( do_long && dur > p.max_dur && dur <= p.long_dur )
                    ++cnt_long_rem_cortical_only;
                  else
                    ++cnt_rem_cortical_only;
                }
	              for (int e=0; e<confirmed.size(); e++)
                {
                  const std::pair<double,double> tint = event_interval( confirmed[e] );
                  const double t0 = tint.first;
                  const double t1 = tint.second;
                  const double dur = t1 - t0;
                  if ( do_long && dur > p.max_dur && dur <= p.long_dur )
                    ++cnt_rem_emg_confirmed_long;
                }
              evts = confirmed;
            }

          // Place onset/offset at a configurable fraction of each event's own peak fast-band LEVEL.  The rise
          // feature is used to *detect* the transition, but it decays mid-event as the
          // baseline adapts, so it is a poor ruler for extent.  The (centered, zero-phase)
          // band-power level stays elevated for the whole event, so its level crossings
          // coincide with the true onset and offset regardless of arousal strength.
          auto activation = [&]( int j ) {
            return std::max( D[j](IDX_THETA) , std::max( D[j](IDX_ALPHA) , D[j](IDX_BETA) ) );
          };
          for (int e=0; e<evts.size(); e++)
            {
              const int a = evts[e].first, b = evts[e].second;
              if ( b <= a ) continue;
              int pk = a;
              double pkv = activation( a );
              for (int j=a+1; j<=b; j++)
                {
                  const double v = activation( j );
                  if ( v > pkv ) { pkv = v; pk = j; }
                }
              if ( pkv <= 0 ) continue;
              const double level = p.duration_level_frac * pkv;
              // walk outward across the whole contig (not just the grown event, whose end was
              // set by the decaying rise) until the level falls below the configured fraction of peak
              int nf = pk, nl = pk;
              while ( nf > 0      && activation( nf-1 ) >= level ) --nf;
              while ( nl < ne - 1 && activation( nl+1 ) >= level ) ++nl;
              // cap to the same max the rest of the pipeline enforces, so the refined
              // event always lands in a valid duration bucket (avoids runaway extents)
              evts[e] = trim_to_candidate_max( std::make_pair( nf , nl ) );
            }

          // This is the final duration gate after level-based boundary
          // refinement.  Emit it separately so the diagnostic track matches
          // the events eligible for final standard/long annotation output.
          add_verbose_events( "candidate_level_refined" , evts );
          std::vector<std::pair<int,int> > final_duration_evts;
          for (int e=0; e<evts.size(); e++)
            {
              const std::pair<double,double> tint = event_interval( evts[e] );
              const double dur = tint.second - tint.first;
              const bool duration_ok =
                dur >= p.arousal_dur &&
                ( dur <= p.max_dur || ( do_long && dur <= p.long_dur ) );
              if ( duration_ok ) final_duration_evts.push_back( evts[e] );
            }
          add_verbose_events( "candidate_final_duration_ok" , final_duration_evts );
          evts = final_duration_evts;

          std::vector<std::pair<int,int> > arr_major;
	          for (int e=0; e<evts.size(); e++)
	            {
	              const std::pair<double,double> tint = event_interval( evts[e] );
	              const double t0 = tint.first;
	              const double t1 = tint.second;
	              const double dur = t1 - t0;
	      if ( do_long && dur > p.max_dur && dur <= p.long_dur )
		{
		  ret[ "arousal_long_" + stg_lab ].insert( interval_t( globals::tp_1sec * t0 , globals::tp_1sec * t1 ) );
		  ret[ "arousal_long" ].insert( interval_t( globals::tp_1sec * t0 , globals::tp_1sec * t1 ) );
		  ret[ "arousal_" + stg_lab ].insert( interval_t( globals::tp_1sec * t0 , globals::tp_1sec * t1 ) );
		  ret[ "arousal" ].insert( interval_t( globals::tp_1sec * t0 , globals::tp_1sec * t1 ) );
		  arr_major.push_back( evts[e] );
		  dur_long += dur;
		  ++cnt_long_evts;
		}
	      else if ( dur >= p.arousal_dur && dur <= p.max_dur )
		{
		  ret[ "arousal_std_" + stg_lab ].insert( interval_t( globals::tp_1sec * t0 , globals::tp_1sec * t1 ) );
		  ret[ "arousal_std" ].insert( interval_t( globals::tp_1sec * t0 , globals::tp_1sec * t1 ) );
		  ret[ "arousal_" + stg_lab ].insert( interval_t( globals::tp_1sec * t0 , globals::tp_1sec * t1 ) );
		  ret[ "arousal" ].insert( interval_t( globals::tp_1sec * t0 , globals::tp_1sec * t1 ) );
		  arr_major.push_back( evts[e] );
		  dur_major += dur;
		  ++cnt_evts;
		}
	      else if ( dur < p.arousal_dur )
		continue;
	      else
		Helper::halt( "internal arousal-detection duration bucket error" );
            }

          std::vector<std::pair<int,int> > arts = mask_to_intervals( artifact );
	          for (int e=0; e<arts.size(); e++)
	            {
	              const std::pair<double,double> tint = event_interval( arts[e] );
	              const double t0 = tint.first;
	              const double t1 = tint.second;
              ret[ "artifact_" + stg_lab ].insert( interval_t( globals::tp_1sec * t0 , globals::tp_1sec * t1 ) );
            }
          cnt_arts += arts.size();

          std::vector<int> arr( ne , 0 );
          for (int e=0; e<arr_major.size(); e++)
            for (int i=arr_major[e].first; i<=arr_major[e].second; i++) arr[i] = 2;

          for (int i=0; i<ne; i++)
            {
              if ( artifact[i] ) continue;
              if ( artifact_blocked[i] ) continue;
              if ( arr[i] == 2 ) { ftr_arousal += D[i]; ++n_arousal; }
              else { ftr_baseline += D[i]; ++n_baseline; }
            }
        }

      writer.level( is_rem ? "R" : "NR" , "SS" );

      auto write_cls = [&]( const std::string & cls , const Eigen::VectorXd & f , const int n )
        {
          writer.level( cls , "CLS" );
          writer.value( "NE" , n );
          if ( n > 0 )
            {
              writer.value( "DELTA" , f(IDX_DELTA) / (double)n );
              writer.value( "THETA" , f(IDX_THETA) / (double)n );
              writer.value( "ALPHA" , f(IDX_ALPHA) / (double)n );
              writer.value( "BETA" , f(IDX_BETA) / (double)n );
              writer.value( "EMG" , f(IDX_EMG) / (double)n );
              writer.value( "EEG_TILT_RISE" , f(IDX_EEG_TILT_RISE) / (double)n );
              writer.value( "EEG_FAST_RISE" , f(IDX_EEG_FAST_RISE) / (double)n );
              writer.value( "EMG_RISE" , f(IDX_EMG_RISE) / (double)n );
              writer.value( "DELTA_ART" , f(IDX_DELTA_ART) / (double)n );
            }
        };

      write_cls( "artifact" , ftr_art , n_art );
      write_cls( "non_artifact" , ftr_nonart , n_nonart );
      write_cls( "arousal" , ftr_arousal , n_arousal );
      write_cls( "baseline" , ftr_baseline , n_baseline );
      writer.unlevel( "CLS" );

      double tot_sec = 0;
      for (int c=0; c<X[st].size(); c++)
        {
          const int ne = tt[st][c].size();
          if ( ne == 0 ) continue;
          tot_sec += tt[st][c][ne-1] - tt[st][c][0] + p.epoch_inc;
        }

      writer.value( "MINS" , tot_sec / 60.0 );
      // Report final unique annotation intervals.  Candidate counters above
      // are retained as diagnostics, but endpoint refinement can make two
      // candidate records resolve to one final interval in the annotation set.
      auto final_count = [&]( const std::string & label ) {
        const std::map<std::string,std::set<interval_t> >::const_iterator it = ret.find( label );
        return it == ret.end() ? 0 : (int)it->second.size();
      };
      auto final_duration = [&]( const std::string & label ) {
        double d = 0;
        const std::map<std::string,std::set<interval_t> >::const_iterator it = ret.find( label );
        if ( it != ret.end() )
          for ( std::set<interval_t>::const_iterator e = it->second.begin(); e != it->second.end(); ++e )
            d += e->duration_sec();
        return d;
      };

      const int n_std_final = final_count( "arousal_std_" + stg_lab );
      const int n_long_final = final_count( "arousal_long_" + stg_lab );
      const double dur_std_final = final_duration( "arousal_std_" + stg_lab );
      const double dur_long_final = final_duration( "arousal_long_" + stg_lab );
      const int n_all_final = n_std_final + n_long_final;
      const double dur_all_final = dur_std_final + dur_long_final;

      // Canonical automated metrics: default = all events; explicit STD/LONG
      // fields split the same events at the 15-second boundary.
      writer.value( "N" , n_all_final );
      writer.value( "AI" , tot_sec > 0 ? n_all_final / ( tot_sec / 3600.0 ) : 0.0 );
      if ( n_all_final > 0 ) writer.value( "DUR" , dur_all_final / (double)n_all_final );
      writer.value( "N_STD" , n_std_final );
      writer.value( "AI_STD" , tot_sec > 0 ? n_std_final / ( tot_sec / 3600.0 ) : 0.0 );
      if ( n_std_final > 0 ) writer.value( "DUR_STD" , dur_std_final / (double)n_std_final );
      writer.value( "N_LONG" , n_long_final );
      writer.value( "AI_LONG" , tot_sec > 0 ? n_long_final / ( tot_sec / 3600.0 ) : 0.0 );
      if ( n_long_final > 0 ) writer.value( "DUR_LONG" , dur_long_final / (double)n_long_final );
      const int n_art_final = final_count( "artifact_" + stg_lab );
      const int n_art_block_final = final_count( "suppressed_" + stg_lab );
      writer.value( "N_ART" , n_art_final );
      writer.value( "AI_ART" , tot_sec > 0 ? n_art_final / ( tot_sec / 3600.0 ) : 0.0 );
      writer.value( "N_SUPPRESSED" , n_art_block_final );
      writer.value( "AI_SUPPRESSED" , tot_sec > 0 ? n_art_block_final / ( tot_sec / 3600.0 ) : 0.0 );

      if ( p.verbose )
        {
          auto vcand_prop = [&]( const int n ) { return vcand_raw > 0 ? n / (double)vcand_raw : 0.0; };
          auto vdiag_prop = [&]( const int n ) { return vdiag_win > 0 ? n / (double)vdiag_win : 0.0; };
          writer.value( "VDIAG_N_WIN" , vdiag_win );
          writer.value( "VDIAG_N_EEG_TILT_SEED" , vdiag_eeg_seed );
          writer.value( "VDIAG_PROP_EEG_TILT_SEED" , vdiag_prop( vdiag_eeg_seed ) );
          writer.value( "VDIAG_N_EMG_RISE_SEED" , vdiag_emg_seed );
          writer.value( "VDIAG_PROP_EMG_RISE_SEED" , vdiag_prop( vdiag_emg_seed ) );
          writer.value( "VDIAG_N_DELTA_ARTIFACT" , vdiag_delta_art );
          writer.value( "VDIAG_PROP_DELTA_ARTIFACT" , vdiag_prop( vdiag_delta_art ) );
          writer.value( "VCAND_N_RAW" , vcand_raw );
          writer.value( "VCAND_N_EEG_RAW" , vcand_eeg_raw );
          writer.value( "VCAND_N_EMG_RAW" , vcand_emg_raw );
          writer.value( "VCAND_N_DUR_OK" , vcand_dur_ok );
          writer.value( "VCAND_PROP_DUR_OK" , vcand_prop( vcand_dur_ok ) );
          writer.value( "VCAND_N_DUR_TOO_SHORT" , vcand_dur_too_short );
          writer.value( "VCAND_PROP_DUR_TOO_SHORT" , vcand_prop( vcand_dur_too_short ) );
          writer.value( "VCAND_N_MERGED" , vcand_merged );
          writer.value( "VCAND_PROP_MERGED" , vcand_prop( vcand_merged ) );
          writer.value( "VCAND_N_FINAL_DUR_OK" , vcand_final_dur_ok );
          writer.value( "VCAND_PROP_FINAL_DUR_OK" , vcand_prop( vcand_final_dur_ok ) );
          writer.value( "VCAND_N_FINAL_DUR_TOO_SHORT" , vcand_final_dur_too_short );
          writer.value( "VCAND_PROP_FINAL_DUR_TOO_SHORT" , vcand_prop( vcand_final_dur_too_short ) );
          writer.value( "VCAND_N_FINAL_DUR_TRIMMED" , vcand_final_dur_trimmed );
          writer.value( "VCAND_PROP_FINAL_DUR_TRIMMED" , vcand_prop( vcand_final_dur_trimmed ) );
          writer.value( "VCAND_N_PRESLEEP_OK" , vcand_presleep_ok );
          writer.value( "VCAND_PROP_PRESLEEP_OK" , vcand_prop( vcand_presleep_ok ) );
          writer.value( "VCAND_N_ARTIFACT_BLOCKED" , vcand_artifact_block );
          writer.value( "VCAND_PROP_ARTIFACT_BLOCKED" , vcand_prop( vcand_artifact_block ) );
          if ( is_rem )
            {
              writer.value( "VCAND_N_REM_EMG_CONFIRMED" , vcand_rem_emg_confirmed );
              writer.value( "VCAND_PROP_REM_EMG_CONFIRMED" , vcand_prop( vcand_rem_emg_confirmed ) );
              writer.value( "VCAND_N_REM_EMG_REJECTED" , vcand_rem_emg_rejected );
              writer.value( "VCAND_PROP_REM_EMG_REJECTED" , vcand_prop( vcand_rem_emg_rejected ) );
            }
        }

      if ( is_rem && p.verbose )
        {
          writer.value( "N_REM_CORTICAL_ONLY" , cnt_rem_cortical_only );
          const int denom = vcand_rem_emg_confirmed + vcand_rem_emg_rejected;
          if ( denom > 0 )
            writer.value( "PCT_REM_EMG_CONFIRMED" , 100.0 * vcand_rem_emg_confirmed / (double)denom );
          writer.value( "N_LONG_REM_CORTICAL_ONLY" , cnt_long_rem_cortical_only );
          const int denom_long = cnt_rem_emg_confirmed_long + cnt_long_rem_cortical_only;
          if ( denom_long > 0 )
            writer.value( "PCT_LONG_REM_EMG_CONFIRMED" , 100.0 * cnt_rem_emg_confirmed_long / (double)denom_long );
        }

      logger << "  " << ( is_rem ? "R" : "NR" )
             << ": standard=" << n_std_final
             << " long=" << n_long_final
             << " all=" << n_all_final
             << " artifact=" << n_art_final
             << " suppressed=" << n_art_block_final << "\n";

      writer.unlevel( "SS" );
    }

  return ret;
}

std::vector<std::pair<int,int> >
arousals2_t::merge_events_with_gap_sorted( const std::vector<std::pair<int,int> > & events ,
                                           const int max_gap )
{
  if ( events.empty() ) return std::vector<std::pair<int,int> >();
  std::vector<std::pair<int,int> > merged;
  int cur_start = events[0].first;
  int cur_end = events[0].second;
  for (int i=1; i<events.size(); i++)
    {
      const int s = events[i].first;
      const int e = events[i].second;
      if ( s <= cur_end + max_gap )
        {
          if ( e > cur_end ) cur_end = e;
        }
      else
        {
          merged.push_back( std::make_pair( cur_start , cur_end ) );
          cur_start = s;
          cur_end = e;
        }
    }
  merged.push_back( std::make_pair( cur_start , cur_end ) );
  return merged;
}

std::vector<std::pair<int,int> >
arousals2_t::mask_to_intervals( const std::vector<bool> & mask )
{
  std::vector<std::pair<int,int> > out;
  const int n = mask.size();
  int start = -1;
  for (int i=0; i<n; i++)
    {
      if ( mask[i] )
        {
          if ( start == -1 ) start = i;
        }
      else if ( start != -1 )
        {
          out.push_back( std::make_pair( start , i-1 ) );
          start = -1;
        }
    }
  if ( start != -1 )
    out.push_back( std::make_pair( start , n-1 ) );
  return out;
}
