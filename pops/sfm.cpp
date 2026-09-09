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

#include "pops/sfm.h"

#if defined(HAS_LGBM) && defined(HAS_ORT)

#include "edf/edf.h"
#include "edf/slice.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "models/ort-common.h"
#include "models/sleepfm-normalize.h"
#include "dsp/resample.h"

#include <onnxruntime/core/session/onnxruntime_cxx_api.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <memory>

extern logger_t logger;

namespace {

  // Path/window-resampling helpers parallel models/ort-sleepfm.cpp.
  // Whole-record channel normalization is shared in sleepfm-normalize.h.

  std::string normalize_model_path( const std::string & path )
  {
    const std::string p = Helper::expand( path );
    const std::string suffix = ".onnx";
    if ( p.size() >= suffix.size() && p.compare( p.size() - suffix.size(), suffix.size(), suffix ) == 0 ) return p;
    return p + suffix;
  }

  std::string join_path( const std::string & path , const std::string & root )
  {
    if ( path.empty() || path == "." ) return root;
    if ( path.back() == '/' || path.back() == '\\' ) return path + root;
    return path + "/" + root;
  }

  // proper sinc-based resampling (libsamplerate, via Luna's existing
  // dsptools::resample) rather than naive linear interpolation, so that
  // downsampling is correctly anti-aliased -- matching SleepFM's own
  // training-data preprocessing, which applies a low-pass filter before
  // downsampling whenever the source rate exceeds the 128 Hz target
  // (sleepfm/preprocessing/preprocessing.py: filter_signal(), butter(4,
  // ..., btype='low') + filtfilt, applied iff rate > resample_rate).
  // dsptools::resample() is a no-op when rates already match, and its
  // sinc converter anti-aliases on downsampling and is artifact-free on
  // upsampling -- matching that conditional-filter behavior without
  // needing to replicate it explicitly.
  std::vector<double> resample( const std::vector<double> & x , double from , double to , int n )
  {
    if ( x.empty() ) return std::vector<double>( n , 0 );
    std::vector<double> y = dsptools::resample( &x , from , to , SRC_SINC_BEST_QUALITY );
    y.resize( n , 0.0 ); // guarantee exact length, as dsptools::resample_channel() itself does
    return y; // normalize before conversion to the model's float32 input
  }

  int modality_max_channels( const std::string & modality )
  {
    if ( modality == "BAS" ) return 10;
    if ( modality == "RESP" ) return 7;
    if ( modality == "EKG" ) return 2;
    if ( modality == "EMG" ) return 4;
    Helper::halt( "unknown SFM modality: " + modality );
    return 0;
  }

}


pops_sfm_result_t pops_sfm_t::run( edf_t & edf ,
                                   const signal_list_t & signals ,
                                   const std::string & modality ,
                                   const std::string & path ,
                                   const std::string & lib ,
                                   int step_sec ,
                                   const std::vector<int> & E )
{

  const int ne = (int)E.size();

  constexpr int sample_rate = 128;
  constexpr int window_seconds = 300;
  constexpr int window_samples = sample_rate * window_seconds;
  constexpr int embedding_dim = 128;
  constexpr int sequence_length = 60; // 300s / 5s
  constexpr int token_seconds = 5;
  constexpr int tokens_per_epoch = 6; // 30s epoch / 5s token -- POPS' fixed epoch length

  const int max_channels = modality_max_channels( modality );

  if ( ! signals.size() ) Helper::halt( "SFM found no matching channels" );
  if ( signals.size() > max_channels ) Helper::halt( "SFM: too many channels given for modality " + modality );

  if ( edf.is_actually_discontinuous() )
    Helper::halt( "SFM requires a continuous EDF (RECORD-SIZE dur=30 continue)" );

  pops_sfm_result_t res;
  res.tokens      = Eigen::MatrixXd::Constant( ne , tokens_per_epoch * embedding_dim , std::numeric_limits<double>::quiet_NaN() );
  res.epoch_pool  = Eigen::MatrixXd::Constant( ne , embedding_dim , std::numeric_limits<double>::quiet_NaN() );
  res.window_pool = Eigen::MatrixXd::Constant( ne , embedding_dim , std::numeric_limits<double>::quiet_NaN() );
  res.position    = Eigen::MatrixXd::Constant( ne , 2 , std::numeric_limits<double>::quiet_NaN() );

  if ( step_sec <= 0 ) Helper::halt( "SFM step= must be positive" );

  // last_time_point_tp is the tp of the last representable instant
  // (inclusive), not one-past-the-end -- e.g. edf.timeline.epoch() for the
  // final epoch has .stop == last_time_point_tp + 1. All exclusive
  // end-of-recording arithmetic below must use total_tp, not
  // last_time_point_tp directly, or window boundaries anchored at the
  // recording's end land one tick short and silently exclude the final
  // epoch.
  const uint64_t total_tp = edf.timeline.last_time_point_tp + 1;

  if ( total_tp < (uint64_t)window_seconds * globals::tp_1sec )
    {
      logger << "  SFM: recording shorter than one 300s window, skipping\n";
      return res;
    }

  const std::string model = normalize_model_path( join_path( Helper::expand( path ) , lib ) );
  std::ifstream model_file( model );
  if ( ! model_file ) Helper::halt( "SFM ONNX model missing: " + model );

  Ort::Env env( ORT_LOGGING_LEVEL_WARNING , "luna-pops-sfm" );
  Ort::SessionOptions opts;
  opts.SetIntraOpNumThreads(1);
  opts.SetGraphOptimizationLevel( GraphOptimizationLevel::ORT_ENABLE_BASIC );
  Ort::Session session( env , model.c_str() , opts );
  Ort::AllocatorWithDefaultOptions alloc;

  if ( session.GetInputCount() < 2 || session.GetOutputCount() < 2 )
    Helper::halt( "SFM ONNX model must have signal, channel_mask, pooled_embedding, and sequence_embedding tensors" );

  auto inname   = ort_common::named_input( session , alloc , "signal" );
  auto maskname = ort_common::named_input( session , alloc , "channel_mask" );
  auto pooled   = ort_common::named_output( session , alloc , "pooled_embedding" );
  auto sequence = ort_common::named_output( session , alloc , "sequence_embedding" );
  const size_t signal_idx   = ort_common::input_index( session , "signal" );
  const size_t mask_idx     = ort_common::input_index( session , "channel_mask" );
  const size_t pooled_idx   = ort_common::output_index( session , "pooled_embedding" );
  const size_t sequence_idx = ort_common::output_index( session , "sequence_embedding" );
  ort_common::check_shape( ort_common::tensor_shape( session.GetInputTypeInfo(signal_idx) ) ,   {1,-1,window_samples} , "SFM signal input" );
  ort_common::check_shape( ort_common::tensor_shape( session.GetInputTypeInfo(mask_idx) ) ,      {1,-1} ,               "SFM channel-mask input" );
  ort_common::check_shape( ort_common::tensor_shape( session.GetOutputTypeInfo(pooled_idx) ) ,   {1,embedding_dim} ,    "SFM pooled output" );
  ort_common::check_shape( ort_common::tensor_shape( session.GetOutputTypeInfo(sequence_idx) ) , {1,sequence_length,embedding_dim} , "SFM sequence output" );

  const uint64_t win  = (uint64_t)window_seconds * globals::tp_1sec;
  const uint64_t step = (uint64_t)step_sec * globals::tp_1sec;
  Ort::MemoryInfo mem = Ort::MemoryInfo::CreateCpu( OrtArenaAllocator , OrtMemTypeDefault );

  // Phase-align the window grid to the current epoch grid. Normally epoch 0
  // starts at tp=0 (phase=0, i.e. the grid used before resolution=5
  // support). But POPS resolution=5 (hypnodensity.cpp) re-epochs the whole
  // EDF timeline once per stride via edf.timeline.set_epoch(30,30,offset_tp)
  // -- offset_tp = k*5s for k=0..5 -- and reruns level1() (hence this
  // function) fresh each time. Deriving the phase from the live timeline
  // (rather than requiring an explicit parameter) keeps every epoch boundary
  // aligned with a window boundary regardless of stride offset, since
  // 300s = 10 x 30s and offset_tp is always < 30s: without this, an epoch
  // straddling an (unshifted) window boundary would be dropped (left NaN)
  // once every ~10 epochs at every stride but k=0.
  const uint64_t phase = edf.timeline.num_epochs() > 0 ? ( edf.timeline.epoch(0).start % step ) : 0;

  // for each result row, the (start,stop) of its actual physical epoch
  // E[row] -- assumed 30s epochs. Kept in both integer tp-ticks (for exact,
  // rounding-free window-membership comparisons) and seconds (for value
  // computations, e.g. token offsets and position-in-window).
  std::vector<uint64_t> epoch_start_tp( ne ), epoch_stop_tp( ne );
  std::vector<double> epoch_start_sec( ne ), epoch_stop_sec( ne );
  for (int e=0; e<ne; e++)
    {
      interval_t iv = edf.timeline.epoch( E[e] );
      epoch_start_tp[e]  = iv.start;
      epoch_stop_tp[e]   = iv.stop;
      epoch_start_sec[e] = iv.start_sec();
      epoch_stop_sec[e]  = iv.stop_sec();
      if ( std::fabs( (epoch_stop_sec[e] - epoch_start_sec[e]) - tokens_per_epoch * token_seconds ) > 1e-6 )
	Helper::halt( "SFM requires standard 30-second POPS epochs" );
    }

  // Normalization always uses the original whole recording, independent of
  // the live epoch grid (including every resolution=5 stride).
  const auto normalization = sleepfm_preprocessing::recording_normalization( edf , signals , sample_rate );
  logger << "  SFM: standardized per channel over the whole recording\n";

  // runs one window starting at w_start_sec and assigns its epochs; if
  // 'only_fill_gaps' is set, only epochs not already assigned are touched
  // (used for the end-of-recording catch-up window below, so it never
  // clobbers a regular grid window's assignment, only fills what the grid
  // missed)
  auto run_window = [&]( uint64_t start , bool only_fill_gaps )
    {
      interval_t iv( start , start + win );
      const double w_start_sec = iv.start_sec();

      std::vector<float> x( (size_t)max_channels * window_samples , 0 );
      std::unique_ptr<bool[]> mask( new bool[max_channels] );
      for (int c=0; c<max_channels; c++) mask[c] = true;
      for (int c=0; c<signals.size(); c++)
	{
	  slice_t sl( edf , signals(c) , iv );
	  std::vector<double> z = resample( *sl.pdata() , edf.header.sampling_freq( signals(c) ) , sample_rate , window_samples );
	  normalization[c].apply( &z );
	  std::copy( z.begin() , z.end() , x.begin() + (size_t)c * window_samples );
	  mask[c] = false;
	}

      std::array<int64_t,3> xs{1,max_channels,window_samples};
      std::array<int64_t,2> ms{1,max_channels};
      auto xt = Ort::Value::CreateTensor<float>( mem , x.data() , x.size() , xs.data() , 3 );
      auto mt = Ort::Value::CreateTensor<bool>( mem , mask.get() , max_channels , ms.data() , 2 );
      const char * ins[]  = { inname.get() , maskname.get() };
      const char * outs[] = { pooled.get() , sequence.get() };
      Ort::Value ivals[] = { std::move(xt) , std::move(mt) };
      auto out = session.Run( Ort::RunOptions{nullptr} , ins , ivals , 2 , outs , 2 );

      const float * pooled_v   = out[0].GetTensorData<float>();
      const float * sequence_v = out[1].GetTensorData<float>();

      // which epochs are fully contained within this window?
      for (int e=0; e<ne; e++)
	{
	  if ( only_fill_gaps && ! std::isnan( res.window_pool(e,0) ) ) continue;

	  // exact integer-tick comparison -- avoids any floating-point
	  // rounding at the window boundary (e.g. the last epoch of a
	  // window/recording landing a hair outside a double-precision
	  // seconds comparison)
	  if ( epoch_start_tp[e] < iv.start || epoch_stop_tp[e] > iv.stop ) continue;

	  const int first_token = (int)std::lround( ( epoch_start_sec[e] - w_start_sec ) / token_seconds );
	  if ( first_token < 0 || first_token + tokens_per_epoch > sequence_length ) continue;

	  std::array<double,embedding_dim> pool_sum{};
	  pool_sum.fill( 0.0 );
	  for (int j=0; j<tokens_per_epoch; j++)
	    {
	      const float * tok = sequence_v + (size_t)(first_token + j) * embedding_dim;
	      for (int d=0; d<embedding_dim; d++)
		{
		  res.tokens(e, j*embedding_dim + d) = tok[d];
		  pool_sum[d] += tok[d];
		}
	    }
	  for (int d=0; d<embedding_dim; d++)
	    res.epoch_pool(e,d) = pool_sum[d] / tokens_per_epoch;

	  for (int d=0; d<embedding_dim; d++)
	    res.window_pool(e,d) = pooled_v[d];

	  res.position(e,0) = epoch_start_sec[e] - w_start_sec;                         // past context
	  res.position(e,1) = ( w_start_sec + window_seconds ) - epoch_stop_sec[e];     // future context
	}
    };

  int nw = 0;
  for ( uint64_t start = phase ; start + win <= total_tp ; start += step , ++nw )
    run_window( start , false );

  // catch-up window: covers whatever the regular grid left short at the
  // tail of the recording (duration not an exact multiple of 'step'),
  // so no epoch is ever left without SFM values. Safe to subtract here:
  // we already returned early above if total_tp < win.
  bool any_gap = false;
  for (int e=0; e<ne; e++) if ( std::isnan( res.window_pool(e,0) ) ) { any_gap = true; break; }
  if ( any_gap )
    {
      run_window( total_tp - win , true );
      ++nw;
    }

  logger << "  SFM: evaluated " << nw << " " << window_seconds << "s window(s) (step=" << step_sec
	 << "s), modality=" << modality << ", channels=" << signals.size() << "/" << max_channels << "\n";

  return res;
}

#endif
