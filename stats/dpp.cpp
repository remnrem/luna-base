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

#include "stats/dpp.h"

#include "stats/dpp-spec.h"
#include "stats/dpp-filter.h"
#include "stats/dpp-io.h"
#include "stats/dpp-fit.h"
#include "stats/dpp-hypno.h"
#include "stats/dpp-vector.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "timeline/timeline.h"
#include "param.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "miscmath/miscmath.h"
#include "dsp/hilbert.h"
#include "dsp/ipc.h"
#include "dsp/psi.h"
#include "dsp/mse.h"
#include "stats/catch22/catch22.h"
#include "fftw/fftwrap.h"
#include "fftw/cohfft.h"
#include "stats/matrix.h"
#include "db/db.h"
#include "defs/defs.h"

#include <vector>
#include <string>
#include <map>
#include <set>
#include <complex>
#include <chrono>
#include <cmath>
#include <limits>
#include <algorithm>
#include <memory>

extern writer_t writer;
extern logger_t logger;

// DPP stage 2 feature-extraction engine. Reuses existing computational
// primitives throughout (PWELCH, spectral_slope_helper, MiscMath::hjorth/
// skewness/kurtosis/clipped/flat/outliers, mse_t, hilbert_t, ipc_t,
// coherence_t, psi_t, timeline_t::discontinuity/interval_overlaps_masked_region,
// CHEP's masked()) -- see stats/dpp-spec.h and stats/dpp-filter.h for the
// spec grammar and the only new filtering glue. No hypnogram_t/NREM-cycle/
// hypnodensity reads anywhere here (out of scope for DPP itself; see the
// implementation plan's revised SS E/F).

namespace {

  // nearest-sample lookup: no shared/exported version of this exists to
  // call instead (SIGDYN's own is a private, file-local helper too)
  int nearest_sample( const std::vector<uint64_t> & tp , uint64_t t )
  {
    if ( tp.empty() ) return -1;
    std::vector<uint64_t>::const_iterator it = std::lower_bound( tp.begin() , tp.end() , t );
    if ( it == tp.begin() ) return 0;
    if ( it == tp.end() ) return (int)tp.size() - 1;
    const int hi = (int)( it - tp.begin() );
    const int lo = hi - 1;
    const uint64_t dh = *it > t ? *it - t : t - *it;
    const uint64_t dl = t > tp[lo] ? t - tp[lo] : tp[lo] - t;
    return dh < dl ? hi : lo;
  }

  struct dpp_trace_t
  {
    int slot;
    int sr;
    std::vector<double> raw;
    std::vector<double> effective; // == raw, or prefiltered if requested
    std::vector<uint64_t> tp;
  };

  // window-level QC gate, entirely from existing MiscMath primitives
  bool qc_bad( const std::vector<double> & x , double flat_th , double clip_th , double outlier_th )
  {
    if ( x.size() < 2 ) return true;
    if ( flat_th > 0 && MiscMath::flat( x ) > flat_th ) return true;
    if ( clip_th > 0 && MiscMath::clipped( x ) > clip_th ) return true;
    if ( outlier_th > 0 )
      {
	std::vector<bool> inc( x.size() , true );
	const int nout = MiscMath::outliers( &x , outlier_th , &inc );
	if ( (double)nout / (double)x.size() > 0.5 ) return true;
      }
    return false;
  }

  // MASK (interval-based) + CHEP (epoch/channel-based) check, both reused
  // directly; CHEP is skipped gracefully if the recording isn't epoched
  bool window_masked( edf_t & edf , uint64_t tp_start , uint64_t tp_stop , const std::string & ch )
  {
    interval_t iv( tp_start , tp_stop + 1 );
    if ( edf.timeline.interval_overlaps_masked_region( iv ) ) return true;

    if ( edf.timeline.epoched() )
      {
	const uint64_t elen = edf.timeline.epoch_len_tp_uint64_t();
	const uint64_t einc = edf.timeline.epoch_increment_tp();
	const int ne = edf.timeline.num_total_epochs();
	const int e1 = MiscMath::position2leftepoch( tp_start , elen , einc , ne );
	const int e2 = MiscMath::position2rightepoch( tp_stop , elen , einc , ne );
	for (int e=e1; e<=e2; e++)
	  if ( e >= 0 && edf.timeline.masked( e , ch ) ) return true;
      }
    return false;
  }

  struct window_result_t
  {
    bool ok;
    int idx_start_padded;
    int idx_start_report;
    int idx_end;
    window_result_t() : ok(false), idx_start_padded(0), idx_start_report(0), idx_end(0) { }
  };

  // resolves the trailing window ending at 't_tp' for one channel: gap
  // check (over the padded extent, if any), MASK/CHEP check and optional
  // raw-buffer QC gate (both over the reporting sub-range only)
  window_result_t get_window( edf_t & edf , const dpp_trace_t & trace , const std::string & ch ,
			      uint64_t t_tp , double w_sec , double pad_sec ,
			      bool qc , double qc_flat , double qc_clip , double qc_th )
  {
    window_result_t r;

    const int idx_end = nearest_sample( trace.tp , t_tp );
    if ( idx_end < 0 ) return r;

    const int w_samples = (int) std::llround( w_sec * trace.sr );
    const int pad_samples = (int) std::llround( pad_sec * trace.sr );

    const int idx_start_report = idx_end - w_samples + 1;
    const int idx_start_padded = idx_start_report - pad_samples;

    if ( idx_start_padded < 0 ) return r;

    if ( timeline_t::discontinuity( trace.tp , trace.sr , idx_start_padded , idx_end ) ) return r;

    // cover the padded run-in too (not just the reporting sub-range): when
    // pad_sec==0 this is idx_start_report == idx_start_padded, a no-op --
    // but for ENVELOPE/PLV, masked/CHEP-bad samples in the padding would
    // otherwise still feed the filter's settling transient unnoticed
    if ( window_masked( edf , trace.tp[ idx_start_padded ] , trace.tp[ idx_end ] , ch ) ) return r;

    if ( qc )
      {
	std::vector<double> rep( trace.effective.begin() + idx_start_report , trace.effective.begin() + idx_end + 1 );
	if ( qc_bad( rep , qc_flat , qc_clip , qc_th ) ) return r;
      }

    r.ok = true;
    r.idx_start_padded = idx_start_padded;
    r.idx_start_report = idx_start_report;
    r.idx_end = idx_end;
    return r;
  }

  // per-(channel,band,window_sec) cache of the causally-padded, filtered
  // (or raw, if band=="") window buffer and its Hilbert magnitude/phase --
  // cleared at the start of every output time t (nothing is reusable
  // across different t, since the window position itself moves), shared
  // for the rest of that t across every spec needing the identical
  // (channel,band,window) combination. Filtering (a direct-form causal FIR
  // convolution, cost roughly proportional to window-length x filter-taps)
  // and the Hilbert transform are each computed at most once per t per
  // distinct key, however many specs consume them -- confirmed via real
  // timing data (HJORTH/SKEW/KURTOSIS at the identical channel/band/window
  // costing essentially the same, regardless of which cheap statistic is
  // computed afterward) that this redundant recomputation, not filter
  // *design* (comparatively cheap, one Kaiser calculation per call), was
  // the dominant real-world cost once several features share a few bands.
  struct dpp_window_cache_t
  {
    struct entry_t
    {
      bool ok = false;
      // buf: the padded (band!="") or reporting-only (band=="") buffer;
      // drop: samples to trim off the front of buf/hmag/hphase to reach
      // the reporting-only range (0 when band=="", since there's no pad)
      std::vector<double> buf;
      int drop = 0;
      bool have_hilbert = false;
      std::vector<double> hmag , hphase;
    };

    std::map<std::string,entry_t> cache;

    void clear() { cache.clear(); }

    static std::string key( const std::string & ch , const std::string & band , double window_sec )
    {
      return ch + "\x1f" + band + "\x1f" + Helper::dbl2str( window_sec );
    }

    entry_t & resolve( edf_t & edf , const dpp_specs_t & specs , const dpp_trace_t & tr ,
		       const std::string & ch , uint64_t t_tp , double window_sec ,
		       const std::string & band , bool qc , double qc_flat_th ,
		       double qc_clip_th , double qc_th )
    {
      const std::string k = key( ch , band , window_sec );
      auto it = cache.find( k );
      if ( it != cache.end() ) return it->second;

      entry_t e;
      const bool has_band = band != "";
      const dpp_filter_t * filt = has_band ? &specs.filters.at( band ) : NULL;
      const double pad_sec = has_band ? dpp_filters::pad_seconds( *filt , tr.sr ) : 0;

      window_result_t wr = get_window( edf , tr , ch , t_tp , window_sec , pad_sec , qc , qc_flat_th , qc_clip_th , qc_th );
      if ( wr.ok )
	{
	  e.ok = true;
	  e.drop = wr.idx_start_report - wr.idx_start_padded;
	  if ( has_band )
	    {
	      std::vector<double> padded( tr.effective.begin() + wr.idx_start_padded , tr.effective.begin() + wr.idx_end + 1 );
	      e.buf = dpp_filters::apply_band( padded , tr.sr , *filt );
	    }
	  else
	    e.buf.assign( tr.effective.begin() + wr.idx_start_report , tr.effective.begin() + wr.idx_end + 1 );
	}
      return cache.emplace( k , std::move(e) ).first->second;
    }

    // trimmed (reporting-only) time-domain buffer -- PSD/SLOPE (band=="")
    // and HJORTH/SKEW/KURTOSIS/MSE/CATCH22 (band== "" or a named band)
    bool get_window_buf( edf_t & edf , const dpp_specs_t & specs , const dpp_trace_t & tr ,
			 const std::string & ch , uint64_t t_tp , double window_sec ,
			 const std::string & band , bool qc , double qc_flat_th ,
			 double qc_clip_th , double qc_th , std::vector<double> * out )
    {
      entry_t & e = resolve( edf , specs , tr , ch , t_tp , window_sec , band , qc , qc_flat_th , qc_clip_th , qc_th );
      if ( ! e.ok ) return false;
      out->assign( e.buf.begin() + e.drop , e.buf.end() );
      return true;
    }

    // trimmed Hilbert magnitude/phase (front-padding dropped -- asymmetric/
    // causal, NOT ipc_param_t::edge_drop_sec, which trims symmetrically and
    // would cut into the live end of a window that's only padded on the
    // causal/start side; confirmed via dsp/ipc.cpp:592-594) -- ENVELOPE/
    // PLV/PAC (band!="" always, enforced at parse time)
    bool get_hilbert( edf_t & edf , const dpp_specs_t & specs , const dpp_trace_t & tr ,
		      const std::string & ch , uint64_t t_tp , double window_sec ,
		      const std::string & band , bool qc , double qc_flat_th ,
		      double qc_clip_th , double qc_th ,
		      std::vector<double> * mag , std::vector<double> * phase )
    {
      entry_t & e = resolve( edf , specs , tr , ch , t_tp , window_sec , band , qc , qc_flat_th , qc_clip_th , qc_th );
      if ( ! e.ok ) return false;
      if ( ! e.have_hilbert )
	{
	  hilbert_t h( e.buf ); // plain ctor: already filtered
	  e.hmag = *h.magnitude();
	  e.hphase = *h.phase();
	  e.have_hilbert = true;
	}
      if ( mag )   mag->assign( e.hmag.begin() + e.drop , e.hmag.end() );
      if ( phase ) phase->assign( e.hphase.begin() + e.drop , e.hphase.end() );
      return true;
    }
  };

  // segmented-Welch parameters that scale sanely down to short windows:
  // default 4s/2s (matching Luna convention elsewhere, e.g. PWELCH's own
  // default), but capped so at least ~2 segments always fit
  void welch_params( double w_sec , double * segment_sec , int * noverlap_segments , int n_samples , int sr )
  {
    *segment_sec = std::min( 4.0 , w_sec / 2.0 );
    const double overlap_sec = *segment_sec / 2.0;
    const int segment_points = (int) std::llround( *segment_sec * sr );
    const int noverlap_points = (int) std::llround( overlap_sec * sr );
    int n = (int) std::floor( ( n_samples - noverlap_points ) / (double)( segment_points - noverlap_points ) );
    *noverlap_segments = n < 1 ? 1 : n;
  }

}


void dsptools::dpp( edf_t & edf , param_t & param )
{

  // Low-rate/vector observations (e.g. SleepFM embeddings) are already
  // window-level samples and must not enter the scalar DSP/specification
  // path below.  Keep this path separate: it accepts fractional sampling
  // rates and uses actual time-aligned samples rather than integer-Hz
  // window arithmetic.
  if ( dpp_vector::enabled( param ) )
    {
      dpp_vector::run( edf , param );
      return;
    }

  //
  // build the feature spec
  //

  dpp_specs_t specs;

  if ( param.has( "windows" ) )
    {
      std::vector<std::string> tok = Helper::parse( param.value( "windows" ) , ',' );
      double d;
      if ( tok.size() > 0 && Helper::str2dbl( tok[0] , &d ) && d > 0 ) specs.default_window_sec = d;
    }

  if ( param.has( "step" ) ) specs.default_step_sec = param.requires_dbl( "step" );

  const bool has_spec_file = param.has( "spec" );

  if ( has_spec_file )
    {
      specs.read( param.value( "spec" ) );
    }
  else
    {
      signal_list_t sl = edf.header.signal_list( param.requires( "sig" ) , true );
      if ( sl.size() == 0 ) { logger << "  *** no signals selected for DPP\n"; return; }
      std::vector<std::string> channels;
      for (int i=0; i<sl.size(); i++) channels.push_back( sl.label(i) );
      specs.init_default( channels );
      specs.apply_inline_overrides( param.has( "windows" )   ? param.value( "windows" )   : "" ,
				    param.has( "filters" )   ? param.value( "filters" )   : "" ,
				    param.has( "features" )  ? param.value( "features" )  : "" ,
				    param.has( "prefilter" ) ? param.value( "prefilter" ) : "" );
    }

  if ( specs.specs.size() == 0 ) { logger << "  *** no DPP features specified\n"; return; }

  // label_root() (stats/dpp-spec.cpp) does not include the spec's own
  // block name, so two blocks that happen to specify the identical
  // feature/channel(s)/band/window would silently collide on the same
  // output column label (VAR level, and binary-matrix column position is
  // still separate but then ambiguous to any downstream consumer keyed
  // off the label) -- fail loudly instead of merging/overwriting silently
  {
    std::set<std::string> seen;
    for (int i=0; i<(int)specs.specs.size(); i++)
      {
	const std::string full_label = specs.specs[i].label_root() + ".w" + Helper::dbl2str( specs.specs[i].window_sec );
	if ( seen.find( full_label ) != seen.end() )
	  Helper::halt( "duplicate DPP feature specification (two blocks resolve to the same label): " + full_label );
	seen.insert( full_label );
      }
  }

  // total LGBM feature columns (sum of cols() across every spec) -- a
  // materially more useful number than the raw spec count alone, since
  // e.g. one PSD spec is 5 columns, one CATCH22 spec is 22-24; computed
  // here (depends only on specs.specs, not on traces) so it's available
  // for the options summary below, not just later at data=/mat build time
  int n_features_total = 0;
  for (int i=0; i<(int)specs.specs.size(); i++) n_features_total += specs.specs[i].cols();

  const bool any_plv = std::any_of( specs.specs.begin() , specs.specs.end() ,
				    []( const dpp_spec_t & s ) { return s.ftr == DPP_PLV; } );

  const double step_sec = specs.default_step_sec > 0 ? specs.default_step_sec : 30;

  // interactive SEC x VAR output is off by default: a cohort-level
  // data=/hypno= corpus-writing run (the common case, e.g. windows=30
  // step=5 over a full night x many individuals) has no use for it, and
  // without an explicit -o it would otherwise dump the entire feature
  // table as plain text to stdout for every individual (Luna's default,
  // no-'-o' writer behaviour is "print", not "discard"). data=/hypno=
  // themselves are unaffected either way -- they're written directly via
  // dpp_io::save(), independent of the writer/db mechanism entirely.
  const bool show_features = param.has( "show-features" ) ? param.yesno( "show-features" ) : false;

  // instrumentation: verbose=T logs each (t, feature-spec) as it starts
  // computing (useful to see exactly what's running if the recording is
  // slow/appears stuck); a per-spec cumulative timing breakdown and a
  // low-frequency wall-clock progress heartbeat are always printed
  // (negligible overhead -- a couple of std::chrono::now() calls per spec,
  // versus the actual filter/FFT/entropy work being timed) since they're
  // useful for any run, not just when actively debugging
  const bool verbose = param.has( "verbose" ) ? param.yesno( "verbose" ) : false;

  const bool qc = param.has( "qc" ) ? param.yesno( "qc" ) : true;
  const double qc_flat_th = param.has( "qc-flat" ) ? param.requires_dbl( "qc-flat" ) : 0.05;
  const double qc_clip_th = param.has( "qc-clip" ) ? param.requires_dbl( "qc-clip" ) : 0.05;
  const double qc_th      = param.has( "qc-th" )   ? param.requires_dbl( "qc-th" )   : 6;

  // PLV amplitude-weighting/gating (ipc_param_t, dsp/ipc.h): previously
  // always left at struct defaults; now exposed. edge_drop_sec is *not*
  // exposed -- DPP already trims the causal pad manually before calling
  // compute_ipc(), and edge_drop_sec trims symmetrically (see stats/dpp.cpp's
  // filtered_report()/get_window() comments), so it stays 0 here regardless
  ipc_param_t plv_param;
  plv_param.amplitude_weighting = param.has( "plv-weighted" ) ? param.yesno( "plv-weighted" ) : true;
  plv_param.gate_low_amp = param.has( "plv-gate" ) ? param.yesno( "plv-gate" ) : true;
  if ( param.has( "plv-gate-abs" ) )
    {
      plv_param.gate_use_quantile = false;
      plv_param.gate_abs = param.requires_dbl( "plv-gate-abs" );
    }
  else
    {
      plv_param.gate_use_quantile = true;
      plv_param.gate_quantile = param.has( "plv-gate-q" ) ? param.requires_dbl( "plv-gate-q" ) : 0.30;
      if ( plv_param.gate_quantile < 0 || plv_param.gate_quantile >= 1 )
	Helper::halt( "plv-gate-q must be in [0,1)" );
    }

  logger << "  DPP options:\n"
	 << "    spec        = " << ( has_spec_file ? param.value( "spec" ) : "(default)" ) << "\n"
	 << "    step        = " << step_sec << "\n"
	 << "    show-features = " << ( show_features ? "T" : "F" ) << "\n"
	 << "    qc          = " << ( qc ? "T" : "F" ) << "\n";
  if ( any_plv )
    logger << "    plv-weighted= " << ( plv_param.amplitude_weighting ? "T" : "F" ) << "\n"
	   << "    plv-gate    = " << ( plv_param.gate_low_amp ? "T" : "F" ) << "\n";
  logger << "    n feature specs   = " << specs.specs.size() << "\n"
	 << "    n feature columns = " << n_features_total << " (passed to LGBM at --dpp-fit/model= time)\n";

  //
  // pull each declared channel's whole trace once; optionally prefilter
  // (once, whole-trace -- no padding needed, see plan Stage 2 SS1)
  //

  std::map<std::string,dpp_trace_t> traces;

  std::map<std::string,dpp_channel_t>::const_iterator cc = specs.chs.begin();
  while ( cc != specs.chs.end() )
    {
      const std::string & ch = cc->first;

      if ( ! edf.header.has_signal( ch ) )
	Helper::halt( "signal " + ch + " not present in the attached EDF" );

      signal_list_t sl = edf.header.signal_list( ch );
      if ( sl.size() != 1 ) Helper::halt( "could not resolve signal " + ch );

      dpp_trace_t tr;
      tr.slot = sl(0);
      tr.sr = (int) edf.header.sampling_freq( tr.slot );

      // sub-1Hz channels truncate to sr=0 here, which get_window()/
      // timeline_t::discontinuity() then divide by (globals::tp_1sec/sr) --
      // halt with a clear message rather than crash on that division.
      // Sub-1Hz signals (e.g. raw POPS hypnodensity) are explicitly
      // out of scope for DPP's own feature-window engine -- see
      // stats/dpp-hypno.h's separate nearest-timepoint lookup path, used
      // instead for exactly this case via hypno-prefix=
      if ( tr.sr < 1 )
	Helper::halt( "signal " + ch + " has sampling rate < 1 Hz ("
		     + Helper::dbl2str( edf.header.sampling_freq( tr.slot ) ) +
		     " Hz) -- not supported as a DPP feature channel" );

      slice_t wtrace( edf , tr.slot , edf.timeline.wholetrace() );
      tr.raw = *wtrace.pdata();
      tr.tp  = *wtrace.ptimepoints();

      if ( cc->second.has_prefilter )
	tr.effective = dpp_filters::prefilter_trace( tr.raw , tr.tp , tr.sr , cc->second.prefilter_lwr , cc->second.prefilter_upr );
      else
	tr.effective = tr.raw;

      traces[ ch ] = tr;

      ++cc;
    }

  // two-channel features require matching sample rates (same requirement
  // PSI's own existing command already imposes) *and* a shared timepoint
  // grid: get_window() resolves each channel independently by nearest-
  // sample lookup, so two channels must have identical tp[] (same start,
  // same length, same gaps) for a given (t, window) request to select the
  // same absolute samples on both sides and produce equal-length buffers.
  // Within one EDF, two signals of matching Fs necessarily share the same
  // record grid (record count/duration are file-level, not per-signal),
  // so tp.size() equality is both necessary and (given matching Fs, this
  // strongly) sufficient in practice; checking it directly here is more
  // robust than assuming it from Fs alone (e.g. if a channel had been
  // independently resampled/added with a different effective grid)
  for (int i=0; i<(int)specs.specs.size(); i++)
    {
      const dpp_spec_t & s = specs.specs[i];
      if ( s.ch.size() != 2 ) continue;
      const dpp_trace_t & t1 = traces[ s.ch[0] ];
      const dpp_trace_t & t2 = traces[ s.ch[1] ];
      if ( t1.sr != t2.sr )
	Helper::halt( "connectivity features require matching sample rates: " + s.ch[0] + ", " + s.ch[1] );
      if ( t1.tp.size() != t2.tp.size() || ( ! t1.tp.empty() && ( t1.tp.front() != t2.tp.front() || t1.tp.back() != t2.tp.back() ) ) )
	Helper::halt( "connectivity features require channels on the same timepoint grid (matching start/end/length): " + s.ch[0] + ", " + s.ch[1] );
    }

  //
  // output times: elapsed seconds, step_sec apart, based on the longest
  // available channel trace's own last timepoint
  //

  uint64_t last_tp = 0;
  {
    std::map<std::string,dpp_trace_t>::const_iterator ii = traces.begin();
    while ( ii != traces.end() )
      {
	if ( ! ii->second.tp.empty() && ii->second.tp.back() > last_tp ) last_tp = ii->second.tp.back();
	++ii;
      }
  }
  const double duration_sec = last_tp / (double)globals::tp_1sec;
  logger << "    duration    = " << Helper::dbl2str( duration_sec , 1 ) << "s (~"
	 << (int) std::floor( duration_sec / step_sec ) << " output times)\n";

  //
  // optional binary corpus output
  //

  const bool write_binary = param.has( "data" );
  dpp_matrix_t mat;
  mat.id = edf.id;

  //
  // optional hypnodensity lookup (Stage 4: parallel corpus for --dpp-fit's
  // stage-conditioned training via hypno=, and/or per-row blending inputs
  // for a stage-conditioned DPP model= apply) -- either hypno= (write) or
  // hypno-prefix= (apply) being explicitly given triggers this; the
  // underlying per-t lookup is shared by both consumers. See stats/dpp-
  // hypno.h: reads already-attached POPS PP_* signals directly, no
  // hypnogram_t/staging dependency.
  //

  const bool want_hypno = param.has( "hypno" ) || param.has( "hypno-prefix" );
  const bool write_hypno = param.has( "hypno" );
  std::unique_ptr<dpp_hypno::lookup_t> hypno_lookup;
  if ( want_hypno )
    {
      const std::string hypno_prefix = param.has( "hypno-prefix" ) ? param.value( "hypno-prefix" ) : "PP";
      const bool hypno_three_state = param.has( "hypno-three-state" ) ? param.yesno( "hypno-three-state" ) : false;
      hypno_lookup.reset( new dpp_hypno::lookup_t( edf , hypno_prefix , hypno_three_state ) );
    }
  dpp_matrix_t hypno_mat;
  hypno_mat.id = edf.id;

  // instrumentation state (see verbose= above): per-spec cumulative wall-
  // clock time, and a ~10s-interval progress heartbeat
  std::vector<double> spec_time_sec( specs.specs.size() , 0.0 );
  const auto run_start_wall = std::chrono::steady_clock::now();
  auto last_heartbeat_wall = run_start_wall;

  // shared per-(channel,band,window) filtered/Hilbert cache -- cleared at
  // the top of every output time t (see dpp_window_cache_t above)
  dpp_window_cache_t win_cache;

  //
  // main loop
  //

  for (double t = step_sec ; t <= duration_sec ; t += step_sec )
    {

      const uint64_t t_tp = (uint64_t) std::llround( t * globals::tp_1sec );

      win_cache.clear();

      if ( show_features ) writer.level( t , globals::sec_strat );

      {
	const auto now_wall = std::chrono::steady_clock::now();
	const double since_heartbeat = std::chrono::duration<double>( now_wall - last_heartbeat_wall ).count();
	if ( since_heartbeat >= 10.0 )
	  {
	    const double elapsed = std::chrono::duration<double>( now_wall - run_start_wall ).count();
	    const double pct = duration_sec > 0 ? 100.0 * t / duration_sec : 0.0;
	    // ETA extrapolated from the current (elapsed / fraction-done) rate --
	    // only meaningful once some real progress has been made
	    std::string eta_str = "n/a";
	    if ( t > 0 )
	      {
		const double est_total = elapsed * ( duration_sec / t );
		const double remaining = est_total - elapsed;
		eta_str = Helper::dbl2str( std::max( 0.0 , remaining ) , 0 ) + "s (est. total "
		  + Helper::dbl2str( est_total , 0 ) + "s)";
	      }
	    logger << "  ... DPP progress: t=" << t << "/" << duration_sec << "s ("
		   << Helper::dbl2str( pct , 1 ) << "%), elapsed=" << Helper::dbl2str( elapsed , 1 )
		   << "s, ETA remaining=" << eta_str << "\n";
	    last_heartbeat_wall = now_wall;
	  }
      }

      std::vector<double> row_vals;
      row_vals.reserve( n_features_total );

      for (int si=0; si<(int)specs.specs.size(); si++)
	{

	  const dpp_spec_t & spec = specs.specs[si];
	  const int ncol = spec.cols();
	  std::vector<double> vals( ncol , std::numeric_limits<double>::quiet_NaN() );

	  const std::string spec_label = spec.label_root() + ".w" + Helper::dbl2str( spec.window_sec );
	  if ( verbose ) logger << "  [t=" << t << "s] " << spec_label << "\n";
	  const auto spec_start_wall = std::chrono::steady_clock::now();

	  if ( spec.ftr == DPP_PSD || spec.ftr == DPP_SLOPE || spec.ftr == DPP_HJORTH ||
	       spec.ftr == DPP_SKEW || spec.ftr == DPP_KURTOSIS || spec.ftr == DPP_MSE ||
	       spec.ftr == DPP_CATCH22 )
	    {
	      const dpp_trace_t & tr = traces[ spec.ch[0] ];
	      // PSD/SLOPE always pass band="" (rejected at parse time
	      // otherwise); HJORTH/SKEW/KURTOSIS/MSE/CATCH22 optionally carry
	      // one band= value each (dpp-spec.cpp's band_expandable handling
	      // -- one spec per requested band, so spec.band has at most one
	      // element here)
	      const std::string band = ( ! spec.band.empty() ) ? spec.band[0] : "";
	      std::vector<double> win;
	      if ( win_cache.get_window_buf( edf , specs , tr , spec.ch[0] , t_tp , spec.window_sec , band ,
					     qc , qc_flat_th , qc_clip_th , qc_th , &win ) )
		{
		  if ( spec.ftr == DPP_PSD || spec.ftr == DPP_SLOPE )
		    {
		      double segment_sec; int noverlap_segments;
		      welch_params( spec.window_sec , &segment_sec , &noverlap_segments , (int)win.size() , tr.sr );
		      PWELCH pwelch( win , tr.sr , segment_sec , noverlap_segments , WINDOW_TUKEY50 );

		      if ( spec.ftr == DPP_PSD )
			{
			  // natural log, unconditionally -- matches POPS's own
			  // canonical-band feature (pops/indiv.cpp, POPS_BANDS:
			  // X1(en,cols[k]) = log(p_<band>)). Floored (not
			  // POPS's own, unguarded log()) rather than risking
			  // -inf on a degenerate all-zero band
			  const double EPS = 1e-300;
			  vals[0] = log( std::max( pwelch.psdsum( DELTA ) , EPS ) );
			  vals[1] = log( std::max( pwelch.psdsum( THETA ) , EPS ) );
			  vals[2] = log( std::max( pwelch.psdsum( ALPHA ) , EPS ) );
			  vals[3] = log( std::max( pwelch.psdsum( SIGMA ) , EPS ) );
			  vals[4] = log( std::max( pwelch.psdsum( BETA ) , EPS ) );
			}
		      else
			{
			  const double fit_lwr = spec.has( "fit-lwr" ) ? spec.narg( "fit-lwr" ) : 1;
			  const double fit_upr = spec.has( "fit-upr" ) ? spec.narg( "fit-upr" ) : 32;
			  std::vector<double> fr = { fit_lwr , fit_upr };
			  double beta=0, betan=0, intercept=0, rsq=0;
			  const bool sok = spectral_slope_helper( pwelch.psd , pwelch.freq , fr , 0 , false ,
								  &beta , &betan , &intercept , &rsq );
			  if ( sok ) vals[0] = beta;
			}
		    }
		  else if ( spec.ftr == DPP_HJORTH )
		    {
		      double act=0, mob=0, comp=0;
		      MiscMath::hjorth( &win , &act , &mob , &comp );
		      // activity (only) log-transformed, matching POPS's own
		      // convention (pops/indiv.cpp, POPS_HJORTH_LEGACY: log(activity),
		      // floored to avoid log(0) rather than POPS's own unguarded
		      // fixed floor) -- mobility/complexity are already
		      // scale-invariant ratios, left as-is
		      const double EPS = 1e-300;
		      vals[0] = log( std::max( act , EPS ) ); vals[1] = mob; vals[2] = comp;
		    }
		  else if ( spec.ftr == DPP_SKEW ) vals[0] = MiscMath::skewness( win );
		  else if ( spec.ftr == DPP_KURTOSIS ) vals[0] = MiscMath::kurtosis( win );
		  else if ( spec.ftr == DPP_MSE )
		    {
		      mse_t mse( 1 , 1 , 1 , 2 , 0.15 );
		      std::map<int,double> m = mse.calc( win );
		      if ( m.find(1) != m.end() ) vals[0] = m[1];
		    }
		  else if ( spec.ftr == DPP_CATCH22 )
		    {
		      const bool catch24 = spec.has( "catch24" ) && Helper::yesno( spec.arg.at( "catch24" ) );
		      catch22_t c22( catch24 );
		      c22.calc( win.data() , (int)win.size() );
		      const int n = catch24 ? 24 : 22;
		      for (int k=0; k<n; k++)
			vals[k] = c22.valid() ? c22.stat(k) : std::numeric_limits<double>::quiet_NaN();
		    }
		}
	    }
	  else if ( spec.ftr == DPP_ENVELOPE )
	    {
	      const dpp_trace_t & tr = traces[ spec.ch[0] ];
	      std::vector<double> mag_report;
	      if ( win_cache.get_hilbert( edf , specs , tr , spec.ch[0] , t_tp , spec.window_sec , spec.band[0] ,
					  qc , qc_flat_th , qc_clip_th , qc_th , &mag_report , NULL ) )
		{
		  const double m = MiscMath::mean( mag_report );
		  const double sd = mag_report.size() > 1 ? MiscMath::sdev( mag_report , m ) : 0;
		  vals[0] = m; vals[1] = sd; vals[2] = ( m != 0 ) ? sd / m : std::numeric_limits<double>::quiet_NaN();
		}
	    }
	  else if ( spec.ftr == DPP_PLV )
	    {
	      const dpp_trace_t & tr1 = traces[ spec.ch[0] ];
	      const dpp_trace_t & tr2 = traces[ spec.ch[1] ];

	      std::vector<double> mag1, ph1, mag2, ph2;
	      const bool ok1 = win_cache.get_hilbert( edf , specs , tr1 , spec.ch[0] , t_tp , spec.window_sec , spec.band[0] ,
						      qc , qc_flat_th , qc_clip_th , qc_th , &mag1 , &ph1 );
	      const bool ok2 = win_cache.get_hilbert( edf , specs , tr2 , spec.ch[1] , t_tp , spec.window_sec , spec.band[1] ,
						      qc , qc_flat_th , qc_clip_th , qc_th , &mag2 , &ph2 );

	      if ( ok1 && ok2 )
		{
		  if ( ph1.size() == ph2.size() && ph1.size() > 0 )
		    {
		      ipc_phaseamp_t seed; seed.phase = ph1; seed.amp = mag1;
		      ipc_phaseamp_t tgt;  tgt.phase  = ph2; tgt.amp  = mag2;
		      ipc_output_t out = ipc_t::compute_ipc( seed , tgt , tr1.sr , plv_param );
		      vals[0] = out.summary.plv;
		      vals[1] = out.summary.mean_ipc;
		    }
		}
	    }
	  else if ( spec.ftr == DPP_PAC )
	    {
	      // phase-amplitude coupling: band[0]/ch[0] is the phase source,
	      // band[1]/ch[1] the amplitude source (spec.band always has
	      // exactly these 2, distinct, roles -- see dpp-spec.cpp's
	      // band_pac handling). Direct per-window mean-vector-length
	      // (Canolty et al. 2006): z = mean( A(t) * exp(i*phase(t)) ),
	      // reported normalized by mean(A) (dimensionless, comparable
	      // across windows/subjects with different absolute amplitude
	      // scale) plus the coupling's preferred phase. No surrogate/
	      // permutation testing (that's dsp/pac.h's pac_t, a CWT-based,
	      // ~1000-repetition whole-recording significance-testing tool --
	      // wrong shape for a per-window ML feature computed potentially
	      // tens of thousands of times per recording); this is a fast,
	      // deterministic point estimate, same spirit as PLV's own use of
	      // compute_ipc() rather than a permutation-based PLV test.
	      const dpp_trace_t & tr1 = traces[ spec.ch[0] ];
	      const dpp_trace_t & tr2 = traces[ spec.ch[1] ];

	      std::vector<double> ph1;   // phase, from the phase source
	      std::vector<double> mag2;  // magnitude, from the amplitude source
	      const bool ok1 = win_cache.get_hilbert( edf , specs , tr1 , spec.ch[0] , t_tp , spec.window_sec , spec.band[0] ,
						      qc , qc_flat_th , qc_clip_th , qc_th , NULL , &ph1 );
	      const bool ok2 = win_cache.get_hilbert( edf , specs , tr2 , spec.ch[1] , t_tp , spec.window_sec , spec.band[1] ,
						      qc , qc_flat_th , qc_clip_th , qc_th , &mag2 , NULL );

	      if ( ok1 && ok2 )
		{
		  if ( ph1.size() == mag2.size() && ph1.size() > 0 )
		    {
		      std::complex<double> z( 0.0 , 0.0 );
		      double sum_mag = 0;
		      for (size_t i=0; i<ph1.size(); i++)
			{
			  z += mag2[i] * std::polar( 1.0 , ph1[i] );
			  sum_mag += mag2[i];
			}
		      z /= (double) ph1.size();
		      const double mean_mag = sum_mag / (double) ph1.size();
		      if ( mean_mag > 0 )
			{
			  vals[0] = std::abs( z ) / mean_mag;
			  vals[1] = std::arg( z );
			}
		    }
		}
	    }
	  else if ( spec.ftr == DPP_COH )
	    {
	      const dpp_trace_t & tr1 = traces[ spec.ch[0] ];
	      const dpp_trace_t & tr2 = traces[ spec.ch[1] ];
	      window_result_t wr1 = get_window( edf , tr1 , spec.ch[0] , t_tp , spec.window_sec , 0 ,
						qc , qc_flat_th , qc_clip_th , qc_th );
	      window_result_t wr2 = get_window( edf , tr2 , spec.ch[1] , t_tp , spec.window_sec , 0 ,
						qc , qc_flat_th , qc_clip_th , qc_th );
	      if ( wr1.ok && wr2.ok )
		{
		  std::vector<double> x1( tr1.effective.begin() + wr1.idx_start_report , tr1.effective.begin() + wr1.idx_end + 1 );
		  std::vector<double> x2( tr2.effective.begin() + wr2.idx_start_report , tr2.effective.begin() + wr2.idx_end + 1 );

		  // defense in depth: the up-front tp-grid check (near the top of
		  // this function) guarantees this in normal operation, but
		  // coherence_t::process() (fftw/cohfft.cpp) indexes both channels'
		  // per-segment spectra by a shared loop bound with no bounds check
		  // of its own -- never risk that on mismatched buffer sizes
		  if ( x1.size() == x2.size() )
		  {
		  double segment_sec; int noverlap_segments;
		  welch_params( spec.window_sec , &segment_sec , &noverlap_segments , (int)x1.size() , tr1.sr );

		  coherence_t coh( (int)x1.size() , tr1.sr , segment_sec , segment_sec / 2.0 );
		  coh.clear();
		  coh.prepare( 1 , x1 );
		  coh.prepare( 2 , x2 );
		  coh.process( 1 , 2 );
		  coh.res.proc_and_output( coh , false );

		  if ( spec.band.size() == 1 && spec.band[0] != "" )
		    {
		      const dpp_filter_t & filt = specs.filters[ spec.band[0] ];
		      const std::vector<double> & frq = coh.frq();
		      double sum_coh = 0; int n = 0;
		      for (int k=0; k<(int)frq.size(); k++)
			{
			  if ( frq[k] < filt.lwr || frq[k] > filt.upr ) continue;
			  if ( coh.res.bad[k] ) continue;
			  const double phi2 = std::norm( coh.res.sxy[k] );
			  sum_coh += phi2 / ( coh.res.sxx[k] * coh.res.syy[k] );
			  ++n;
			}
		      if ( n > 0 ) vals[0] = sum_coh / n;
		    }
		  else
		    {
		      vals[0] = coh.res.bcoh[ DELTA ];
		      vals[1] = coh.res.bcoh[ THETA ];
		      vals[2] = coh.res.bcoh[ ALPHA ];
		      vals[3] = coh.res.bcoh[ SIGMA ];
		      vals[4] = coh.res.bcoh[ BETA ];
		    }
		  }
		}
	    }
	  else if ( spec.ftr == DPP_PSI )
	    {
	      const dpp_trace_t & tr1 = traces[ spec.ch[0] ];
	      const dpp_trace_t & tr2 = traces[ spec.ch[1] ];
	      window_result_t wr1 = get_window( edf , tr1 , spec.ch[0] , t_tp , spec.window_sec , 0 ,
						qc , qc_flat_th , qc_clip_th , qc_th );
	      window_result_t wr2 = get_window( edf , tr2 , spec.ch[1] , t_tp , spec.window_sec , 0 ,
						qc , qc_flat_th , qc_clip_th , qc_th );
	      if ( wr1.ok && wr2.ok && spec.band.size() == 1 )
		{
		  std::vector<double> x1( tr1.effective.begin() + wr1.idx_start_report , tr1.effective.begin() + wr1.idx_end + 1 );
		  std::vector<double> x2( tr2.effective.begin() + wr2.idx_start_report , tr2.effective.begin() + wr2.idx_end + 1 );

		  // defense in depth (see COH above): psi_t's X(i,1)=x2[i] indexing
		  // below assumes matching lengths with no bounds check of its own
		  if ( x1.size() == x2.size() )
		  {

		  Data::Matrix<double> X( (int)x1.size() , 2 );
		  for (int i=0; i<(int)x1.size(); i++) { X(i,0) = x1[i]; X(i,1) = x2[i]; }

		  const int eplen = std::min( (int)x1.size() , tr1.sr * 4 );
		  const int seglen = std::max( 1 , eplen / 2 );

		  psi_t psi( &X , eplen , seglen , tr1.sr );
		  const dpp_filter_t & filt = specs.filters[ spec.band[0] ];
		  psi.add_freqbin( filt.lwr , filt.upr );
		  psi.calc();

		  if ( psi.n_models > 0 && psi.psi.size() > 0 )
		    {
		      const double EPS = 1e-8;
		      const double raw = psi.psi[0](0,1);
		      const double sd  = psi.std_psi[0](0,1);
		      vals[0] = raw / ( EPS + sd );
		    }
		  }
		}
	    }

	  spec_time_sec[si] += std::chrono::duration<double>( std::chrono::steady_clock::now() - spec_start_wall ).count();

	  for (int k=0; k<ncol; k++) row_vals.push_back( vals[k] );

	  if ( show_features )
	    {
	      writer.level( spec_label , globals::var_strat );
	      for (int k=0; k<ncol; k++)
		if ( Helper::realnum( vals[k] ) )
		  writer.value( ncol == 1 ? "V" : ( "V" + Helper::int2str( k+1 ) ) , vals[k] );
	      writer.unlevel( globals::var_strat );
	    }

	}

      mat.time_sec.push_back( t );
      mat.X.push_back( row_vals );

      if ( hypno_lookup )
	{
	  hypno_mat.time_sec.push_back( t );
	  hypno_mat.X.push_back( hypno_lookup->at( t_tp ) );
	}

      if ( show_features ) writer.unlevel( globals::sec_strat );

    }

  // per-spec cumulative timing breakdown, most expensive first -- run
  // unconditionally (a handful of log lines, negligible relative to the
  // computation just timed), so any run's log shows where the time
  // actually went without needing to re-run with verbose=T
  {
    const double total_wall = std::chrono::duration<double>( std::chrono::steady_clock::now() - run_start_wall ).count();
    std::vector<int> order( specs.specs.size() );
    for (int i=0; i<(int)order.size(); i++) order[i] = i;
    std::sort( order.begin() , order.end() , [&]( int a , int b ) { return spec_time_sec[a] > spec_time_sec[b]; } );
    logger << "  DPP timing (" << Helper::dbl2str( total_wall , 1 ) << "s total, " << specs.specs.size() << " feature specs):\n";
    for (int i=0; i<(int)order.size(); i++)
      {
	const int si = order[i];
	const std::string spec_label = specs.specs[si].label_root() + ".w" + Helper::dbl2str( specs.specs[si].window_sec );
	const double pct = total_wall > 0 ? 100.0 * spec_time_sec[si] / total_wall : 0.0;
	logger << "    " << spec_label << " : " << Helper::dbl2str( spec_time_sec[si] , 2 )
	       << "s (" << Helper::dbl2str( pct , 1 ) << "%)\n";
      }
  }

  if ( write_binary )
    dpp_io::save( param.value( "data" ) , mat , n_features_total , false );

  if ( write_hypno )
    dpp_io::save( param.value( "hypno" ) , hypno_mat , (int)hypno_lookup->labels.size() , false );

  if ( param.has( "model" ) )
    {
#ifdef HAS_LGBM
      dpp_fit::apply( edf , param , specs , mat , hypno_lookup ? &hypno_mat : NULL );
#else
      Helper::halt( "LGBM support not compiled in" );
#endif
    }

}
