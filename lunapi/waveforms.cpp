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

#include "luna.h"
#include "lunapi/waveforms.h"

#include "helper/helper.h"
#include "helper/logger.h"
#include "stats/catch22/catch22.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>

extern logger_t logger;

namespace {

struct event_seed_t {
  std::string annot;
  std::string instance;
  std::string annot_ch;
  std::string meta;
  interval_t annot_interval;
  double anchor_sec = std::numeric_limits<double>::quiet_NaN();
};

struct channel_cache_t {
  std::string label;
  std::string unit;
  int slot = -1;
  double sr = 0.0;
  double phys_min = 0.0;
  double phys_max = 0.0;
  uint64_t sample_step_tp = 0LLU;
  std::vector<uint64_t> timepoints;
  std::vector<double> times;
  std::vector<double> values;
};

double align_sec(const interval_t & interval, const std::string & align)
{
  if ( align == "start" ) return interval.start_sec();
  if ( align == "stop" ) return interval.stop_sec();
  return interval.mid_sec();
}

std::string normalize_align(const std::string & align)
{
  if ( align == "midpoint" ) return "mid";
  if ( align == "start" || align == "mid" || align == "stop" ) return align;
  Helper::halt( "align must be start, mid or stop" );
  return "mid";
}

std::vector<event_seed_t> collect_events(edf_t & edf, const std::vector<std::string> & annots,
                                         const std::string & align)
{
  std::vector<event_seed_t> events;
  for (int a = 0; a < annots.size(); ++a)
    {
      annot_t * annot = edf.annotations->find( annots[a] );
      if ( annot == NULL ) continue;
      annot_map_t::const_iterator ii = annot->interval_events.begin();
      while ( ii != annot->interval_events.end() )
        {
          const instance_idx_t & idx = ii->first;
          const instance_t * inst = ii->second;
          event_seed_t ev;
          ev.annot = annots[a];
          ev.instance = idx.id;
          ev.annot_ch = idx.ch_str;
          ev.meta = inst != NULL ? inst->print() : "";
          ev.annot_interval = idx.interval;
          ev.anchor_sec = align_sec( idx.interval , align );
          events.push_back( ev );
          ++ii;
        }
    }

  std::sort( events.begin() , events.end() ,
             []( const event_seed_t & a , const event_seed_t & b ) {
               if ( a.annot_interval < b.annot_interval ) return true;
               if ( b.annot_interval < a.annot_interval ) return false;
               if ( a.annot < b.annot ) return true;
               if ( a.annot > b.annot ) return false;
               if ( a.instance < b.instance ) return true;
               if ( a.instance > b.instance ) return false;
               return a.annot_ch < b.annot_ch;
             } );
  return events;
}

bool annotation_channel_matches(const std::string & annot_ch,
                                const std::vector<std::string> & channels)
{
  if ( annot_ch == "" || annot_ch == "." ) return true;
  for (int i = 0; i < channels.size(); ++i)
    if ( annot_ch == channels[i] ) return true;
  return false;
}

std::vector<channel_cache_t> build_channel_caches(edf_t & edf, const std::vector<std::string> & channels)
{
  std::vector<channel_cache_t> caches;
  for (int c = 0; c < channels.size(); ++c)
    {
      const int slot = edf.header.signal( channels[c] );
      if ( slot == -1 ) continue;
      if ( edf.header.is_annotation_channel( slot ) ) continue;
      channel_cache_t ch;
      ch.label = channels[c];
      ch.unit = edf.header.phys_dimension[ slot ];
      ch.slot = slot;
      ch.sr = edf.header.sampling_freq( slot );
      ch.phys_min = edf.header.physical_min[ slot ];
      ch.phys_max = edf.header.physical_max[ slot ];
      if ( ch.sr <= 0 ) continue;

      slice_t slice( edf , slot , edf.timeline.wholetrace() );
      const std::vector<double> * data = slice.pdata();
      const std::vector<uint64_t> * tp = slice.ptimepoints();
      if ( data == NULL || tp == NULL || data->size() == 0 || tp->size() != data->size() ) continue;

      if ( tp->size() >= 2 )
        {
          const uint64_t sr_step_tp =
            ch.sr > 0 ? (uint64_t)std::llround( globals::tp_1sec / ch.sr ) : 0LLU;
          const uint64_t sr_tol_tp =
            sr_step_tp > 0 ? std::max<uint64_t>( 1LLU , sr_step_tp / 20LLU ) : 0LLU;
          uint64_t step_tp = 0LLU;
          for (int i = 1; i < tp->size(); ++i)
            {
              if ( (*tp)[i] <= (*tp)[i - 1] )
                {
                  continue;
                }
              const uint64_t delta = (*tp)[i] - (*tp)[i - 1];
              if ( sr_step_tp != 0LLU )
                {
                  const uint64_t diff = delta > sr_step_tp ? delta - sr_step_tp : sr_step_tp - delta;
                  if ( diff > sr_tol_tp )
                    {
                      continue;
                    }
                }
              if ( step_tp == 0LLU )
                step_tp = delta;
            }
          if ( step_tp == 0LLU ) continue;
          ch.sample_step_tp = step_tp;
        }
      if ( ch.sample_step_tp == 0LLU ) continue;

      ch.values.assign( data->begin() , data->end() );
      ch.timepoints.assign( tp->begin() , tp->end() );
      ch.times.resize( tp->size() );
      for (int i = 0; i < tp->size(); ++i)
        ch.times[i] = (*tp)[i] * globals::tp_duration;
      caches.push_back( ch );
    }
  return caches;
}

int nearest_index(const std::vector<double> & times, const double target)
{
  std::vector<double>::const_iterator it = std::lower_bound( times.begin() , times.end() , target );
  if ( it == times.begin() ) return 0;
  if ( it == times.end() ) return times.size() - 1;
  const int right = it - times.begin();
  const int left = right - 1;
  return std::fabs( target - times[left] ) <= std::fabs( times[right] - target ) ? left : right;
}

std::vector<std::string> feature_names_internal(const waveform_feature_options_t & options)
{
  std::vector<std::string> names;
  if ( options.catch_enabled )
    {
      const int ncatch = options.catch24 ? 24 : 22;
      for (int i = 0; i < ncatch; ++i)
        names.push_back( catch22_t::short_name( i ) );
    }
  if ( options.basic_stats )
    {
      names.push_back( "duration_sec" );
      names.push_back( "range_raw" );
      names.push_back( "mean_raw" );
      names.push_back( "sd_raw" );
    }
  return names;
}

void add_basic_features(waveform_block_t * block)
{
  if ( block == NULL ) return;
  const int n = block->values.size();
  if ( n == 0 )
    {
      block->features["duration_sec"] = std::numeric_limits<double>::quiet_NaN();
      block->features["range_raw"] = std::numeric_limits<double>::quiet_NaN();
      block->features["mean_raw"] = std::numeric_limits<double>::quiet_NaN();
      block->features["sd_raw"] = std::numeric_limits<double>::quiet_NaN();
      return;
    }

  double sum = 0.0;
  double minv = block->values[0];
  double maxv = block->values[0];
  for (int i = 0; i < n; ++i)
    {
      const double v = block->values[i];
      if ( ! std::isfinite( v ) ) continue;
      sum += v;
      if ( v < minv ) minv = v;
      if ( v > maxv ) maxv = v;
    }
  const double mean = sum / (double)n;
  double ss = 0.0;
  for (int i = 0; i < n; ++i)
    {
      const double d = block->values[i] - mean;
      ss += d * d;
    }
  const double sd = n > 1 ? std::sqrt( ss / (double)( n - 1 ) ) : 0.0;
  block->features["duration_sec"] = std::max( 0.0 , block->data_stop_sec - block->data_start_sec );
  block->features["range_raw"] = maxv - minv;
  block->features["mean_raw"] = mean;
  block->features["sd_raw"] = sd;
}

} // namespace

namespace waveform_core
{

std::vector<std::string> feature_names(const waveform_feature_options_t & options)
{
  return feature_names_internal( options );
}

waveform_extract_result_t extract_fixed_window_waveforms(
  edf_t & edf ,
  const std::vector<std::string> & annots ,
  const std::vector<std::string> & channels ,
  const double pre_secs ,
  const double post_secs ,
  const std::string & align0 ,
  const std::string & require )
{
  const std::string align = normalize_align( align0 );
  if ( require != "full" )
    Helper::halt( "fixed-window waveform extraction currently requires require=full" );
  waveform_extract_result_t out;
  out.requested_annots = annots;
  out.requested_channels = channels;
  out.align = align;
  out.mode = "fixed";
  out.pre_secs = pre_secs;
  out.post_secs = post_secs;

  std::vector<event_seed_t> events = collect_events( edf , annots , align );
  out.total_events = events.size();
  std::vector<channel_cache_t> caches = build_channel_caches( edf , channels );
  if ( caches.size() == 0 ) return out;

  for (int e = 0; e < events.size(); ++e)
    {
      waveform_event_t wave;
      wave.annot = events[e].annot;
      wave.instance = events[e].instance;
      wave.annot_ch = events[e].annot_ch;
      wave.meta = events[e].meta;
      wave.annot_start_sec = events[e].annot_interval.start_sec();
      wave.annot_stop_sec = events[e].annot_interval.stop_sec();
      wave.anchor_sec = events[e].anchor_sec;

      for (int c = 0; c < caches.size(); ++c)
        {
          const channel_cache_t & ch = caches[c];
          if ( ch.times.size() < 2 || ch.timepoints.size() != ch.times.size() || ch.sample_step_tp == 0LLU )
            {
              ++out.dropped["too_short"];
              continue;
            }
          const uint64_t pre_tp = (uint64_t)std::llround( std::max( 0.0 , pre_secs ) * globals::tp_1sec );
          const uint64_t post_tp = (uint64_t)std::llround( std::max( 0.0 , post_secs ) * globals::tp_1sec );
          const int n_pre = (int)std::llround( (double)pre_tp / (double)ch.sample_step_tp );
          const int n_post = (int)std::llround( (double)post_tp / (double)ch.sample_step_tp );
          const int center = nearest_index( ch.times , events[e].anchor_sec );
          const int lo = center - n_pre;
          const int hi = center + n_post;
          if ( lo < 0 || hi >= ch.times.size() )
            {
              ++out.dropped["edge"];
              continue;
            }

          const double sample_dt = (double)ch.sample_step_tp * globals::tp_duration;
          waveform_block_t block;
          block.label = ch.label;
          block.unit = ch.unit;
          block.sr = ch.sr;
          block.sample_step_tp = ch.sample_step_tp;
          block.data_start_sec = ch.times[lo];
          block.data_stop_sec = ch.times[hi] + sample_dt;
          block.rel_time.resize( hi - lo + 1 );
          block.values.resize( hi - lo + 1 );
          bool okay = true;
          const uint64_t center_tp = ch.timepoints[center];
          for (int i = lo; i <= hi; ++i)
            {
              const int j = i - lo;
              const int64_t offset = (int64_t)( j - n_pre );
              const int64_t expected_tp = (int64_t)center_tp + offset * (int64_t)ch.sample_step_tp;
              if ( (int64_t)ch.timepoints[i] != expected_tp )
                {
                  okay = false;
                  break;
                }
              block.rel_time[j] = ( (double)ch.timepoints[i] - (double)center_tp ) * globals::tp_duration;
              block.values[j] = ch.values[i];
            }
          if ( ! okay )
            {
              ++out.dropped["gap"];
              continue;
            }

          if ( std::isnan( wave.wave_start_sec ) || block.data_start_sec < wave.wave_start_sec )
            wave.wave_start_sec = block.data_start_sec;
          if ( std::isnan( wave.wave_stop_sec ) || block.data_stop_sec > wave.wave_stop_sec )
            wave.wave_stop_sec = block.data_stop_sec;
          wave.blocks.push_back( block );
        }

      if ( wave.blocks.size() == 0 )
        ++out.dropped["empty"];
      else
        out.events.push_back( wave );
    }
  return out;
}


waveform_extract_result_t extract_annotation_window_waveforms(
  edf_t & edf ,
  const std::vector<std::string> & annots ,
  const std::vector<std::string> & channels ,
  const double flank_left_secs ,
  const double flank_right_secs ,
  const std::string & align0 ,
  const std::string & require ,
  const bool match_annot_channel )
{
  const std::string align = normalize_align( align0 );
  if ( require != "full" && require != "any" )
    Helper::halt( "require must be full or any" );

  waveform_extract_result_t out;
  out.requested_annots = annots;
  out.requested_channels = channels;
  out.align = align;
  out.mode = "annot_window";
  out.flank_left_secs = flank_left_secs;
  out.flank_right_secs = flank_right_secs;

  std::vector<event_seed_t> events = collect_events( edf , annots , align );
  out.total_events = events.size();
  std::vector<channel_cache_t> caches = build_channel_caches( edf , channels );
  if ( caches.size() == 0 ) return out;

  const uint64_t left_tp = (uint64_t)( std::max( 0.0 , flank_left_secs ) * globals::tp_1sec );
  const uint64_t right_tp = (uint64_t)( std::max( 0.0 , flank_right_secs ) * globals::tp_1sec );

  for (int e = 0; e < events.size(); ++e)
    {
      if ( match_annot_channel )
        {
          if ( events[e].annot_ch == "" || events[e].annot_ch == "." ||
               ! annotation_channel_matches( events[e].annot_ch , channels ) )
            {
              ++out.dropped["annot_ch"];
              continue;
            }
        }

      interval_t wave_iv = events[e].annot_interval;
      const uint64_t requested_tp = wave_iv.duration() + left_tp + right_tp;
      if ( left_tp ) wave_iv.expand_left( left_tp );
      if ( right_tp ) wave_iv.expand_right( right_tp );

      const uint64_t retained_tp = edf.timeline.valid_tps( wave_iv );
      if ( requested_tp == 0
           || ( require == "full"
                ? ( wave_iv.duration() != requested_tp || retained_tp != requested_tp )
                : retained_tp == 0LLU ) )
        {
          ++out.dropped["retention"];
          continue;
        }

      waveform_event_t wave;
      wave.annot = events[e].annot;
      wave.instance = events[e].instance;
      wave.annot_ch = events[e].annot_ch;
      wave.meta = events[e].meta;
      wave.annot_start_sec = events[e].annot_interval.start_sec();
      wave.annot_stop_sec = events[e].annot_interval.stop_sec();
      wave.anchor_sec = events[e].anchor_sec;
      wave.wave_start_sec = wave_iv.start_sec();
      wave.wave_stop_sec = wave_iv.stop_sec();

      bool okay = true;
      for (int c = 0; c < caches.size(); ++c)
        {
          const channel_cache_t & ch = caches[c];
          if ( match_annot_channel && ch.label != events[e].annot_ch )
            continue;
          slice_t slice( edf , ch.slot , wave_iv );
          const std::vector<double> * data = slice.pdata();
          const std::vector<uint64_t> * tp = slice.ptimepoints();
          if ( data == NULL || tp == NULL || data->size() == 0 || tp->size() != data->size() )
            {
              okay = false;
              break;
            }

          waveform_block_t block;
          block.label = ch.label;
          block.unit = ch.unit;
          block.sr = ch.sr;
          block.sample_step_tp = ch.sample_step_tp;
          block.phys_min = ch.phys_min;
          block.phys_max = ch.phys_max;
          block.data_start_sec = (*tp)[0] * globals::tp_duration;
          block.data_stop_sec = (*tp)[ tp->size() - 1 ] * globals::tp_duration + ( (double)ch.sample_step_tp * globals::tp_duration );
          block.rel_time.resize( tp->size() );
          block.values.resize( data->size() );
          for (int i = 0; i < tp->size(); ++i)
            {
              const double sec = (*tp)[i] * globals::tp_duration;
              block.rel_time[i] = sec - wave.anchor_sec;
              block.values[i] = (*data)[i];
            }
          wave.blocks.push_back( block );
        }

      const int expected_blocks = match_annot_channel ? 1 : (int)caches.size();
      if ( ! okay || (int)wave.blocks.size() != expected_blocks )
        {
          ++out.dropped["empty"];
          continue;
        }
      out.events.push_back( wave );
    }

  return out;
}

void compute_features(waveform_extract_result_t * result, const waveform_feature_options_t & options)
{
  if ( result == NULL ) return;
  result->feature_names = feature_names_internal( options );

  for (int e = 0; e < result->events.size(); ++e)
    {
      waveform_event_t & ev = result->events[e];
      for (int b = 0; b < ev.blocks.size(); ++b)
        {
          waveform_block_t & block = ev.blocks[b];
          block.features.clear();
          block.feature_qc = -1;
          if ( options.catch_enabled )
            {
              catch22_t c22( options.catch24 );
              block.feature_qc = c22.calc( block.values.data() , block.values.size() ) ? 0 : c22.qc;
              const int ncatch = options.catch24 ? 24 : 22;
              for (int i = 0; i < ncatch; ++i)
                {
                  const std::string name = catch22_t::short_name( i );
                  block.features[name] = block.feature_qc == 0 ? c22.stat( i ) : std::numeric_limits<double>::quiet_NaN();
                }
            }
          if ( options.basic_stats ) add_basic_features( &block );
        }
    }
}
} // namespace waveform_core


waveform_extract_result_t lunapi_inst_t::extract_event_waveforms(
  const std::vector<std::string> & annots ,
  const std::vector<std::string> & chs ,
  const double pre_secs ,
  const double post_secs ,
  const std::string & align ,
  const std::string & require ) const
{
  if ( require != "full" && require != "any" )
    Helper::halt( "require must be full or any" );
  lunapi_inst_t * self = const_cast<lunapi_inst_t*>( this );
  return waveform_core::extract_fixed_window_waveforms( self->edf , annots , chs , pre_secs , post_secs , align , require );
}


waveform_extract_result_t lunapi_inst_t::compute_waveform_features(
  const waveform_extract_result_t & input ,
  const bool catch24 ,
  const bool basic_stats ) const
{
  waveform_extract_result_t out = input;
  waveform_feature_options_t options;
  options.catch_enabled = true;
  options.catch24 = catch24;
  options.basic_stats = basic_stats;
  waveform_core::compute_features( &out , options );
  return out;
}


waveform_extract_result_t lunapi_inst_t::extract_event_waveforms_with_features(
  const std::vector<std::string> & annots ,
  const std::vector<std::string> & chs ,
  const double pre_secs ,
  const double post_secs ,
  const std::string & align ,
  const std::string & require ,
  const bool catch24 ,
  const bool basic_stats ) const
{
  waveform_extract_result_t out = extract_event_waveforms( annots , chs , pre_secs , post_secs , align , require );
  waveform_feature_options_t options;
  options.catch_enabled = true;
  options.catch24 = catch24;
  options.basic_stats = basic_stats;
  waveform_core::compute_features( &out , options );
  return out;
}
