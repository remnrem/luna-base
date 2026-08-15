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
#include "fftw/fftwrap.h"
#include "fftw/cohfft.h"
#include "stats/matrix.h"
#include "db/db.h"
#include "defs/defs.h"

#include <vector>
#include <string>
#include <map>
#include <set>
#include <cmath>
#include <limits>
#include <algorithm>

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

  // filtered magnitude/phase for one (channel,band) input, front-padding
  // dropped (asymmetric/causal -- NOT ipc_param_t::edge_drop_sec, which
  // trims symmetrically and would cut into the live end of a window
  // that's only padded on the causal/start side; confirmed via
  // dsp/ipc.cpp:592-594)
  void filtered_report( const dpp_trace_t & trace , const window_result_t & wr ,
			const dpp_filter_t & filt ,
			std::vector<double> * mag_report , std::vector<double> * phase_report )
  {
    std::vector<double> padded( trace.effective.begin() + wr.idx_start_padded ,
				trace.effective.begin() + wr.idx_end + 1 );
    std::vector<double> filtered = dpp_filters::apply_band( padded , trace.sr , filt );
    hilbert_t h( filtered ); // plain ctor: already filtered
    const int drop = wr.idx_start_report - wr.idx_start_padded;
    if ( mag_report )
      mag_report->assign( h.magnitude()->begin() + drop , h.magnitude()->end() );
    if ( phase_report )
      phase_report->assign( h.phase()->begin() + drop , h.phase()->end() );
  }

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

  const double step_sec = specs.default_step_sec > 0 ? specs.default_step_sec : 30;

  const bool qc = param.has( "qc" ) ? param.yesno( "qc" ) : true;
  const double qc_flat_th = param.has( "qc-flat" ) ? param.requires_dbl( "qc-flat" ) : 0.05;
  const double qc_clip_th = param.has( "qc-clip" ) ? param.requires_dbl( "qc-clip" ) : 0.05;
  const double qc_th      = param.has( "qc-th" )   ? param.requires_dbl( "qc-th" )   : 6;

  logger << "  DPP options:\n"
	 << "    spec        = " << ( has_spec_file ? param.value( "spec" ) : "(default)" ) << "\n"
	 << "    step        = " << step_sec << "\n"
	 << "    qc          = " << ( qc ? "T" : "F" ) << "\n"
	 << "    n features  = " << specs.specs.size() << "\n";

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

  //
  // optional binary corpus output
  //

  const bool write_binary = param.has( "data" );
  dpp_matrix_t mat;
  mat.id = edf.id;
  int n_features_total = 0;
  for (int i=0; i<(int)specs.specs.size(); i++) n_features_total += specs.specs[i].cols();

  //
  // main loop
  //

  for (double t = step_sec ; t <= duration_sec ; t += step_sec )
    {

      const uint64_t t_tp = (uint64_t) std::llround( t * globals::tp_1sec );

      writer.level( t , globals::sec_strat );

      std::vector<double> row_vals;
      row_vals.reserve( n_features_total );

      for (int si=0; si<(int)specs.specs.size(); si++)
	{

	  const dpp_spec_t & spec = specs.specs[si];
	  const int ncol = spec.cols();
	  std::vector<double> vals( ncol , std::numeric_limits<double>::quiet_NaN() );

	  if ( spec.ftr == DPP_PSD || spec.ftr == DPP_SLOPE || spec.ftr == DPP_HJORTH ||
	       spec.ftr == DPP_SKEW || spec.ftr == DPP_KURTOSIS || spec.ftr == DPP_MSE )
	    {
	      const dpp_trace_t & tr = traces[ spec.ch[0] ];
	      window_result_t wr = get_window( edf , tr , spec.ch[0] , t_tp , spec.window_sec , 0 ,
					       qc , qc_flat_th , qc_clip_th , qc_th );
	      if ( wr.ok )
		{
		  std::vector<double> win( tr.effective.begin() + wr.idx_start_report , tr.effective.begin() + wr.idx_end + 1 );

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
		      vals[0] = act; vals[1] = mob; vals[2] = comp;
		    }
		  else if ( spec.ftr == DPP_SKEW ) vals[0] = MiscMath::skewness( win );
		  else if ( spec.ftr == DPP_KURTOSIS ) vals[0] = MiscMath::kurtosis( win );
		  else if ( spec.ftr == DPP_MSE )
		    {
		      mse_t mse( 1 , 1 , 1 , 2 , 0.15 );
		      std::map<int,double> m = mse.calc( win );
		      if ( m.find(1) != m.end() ) vals[0] = m[1];
		    }
		}
	    }
	  else if ( spec.ftr == DPP_ENVELOPE )
	    {
	      const dpp_trace_t & tr = traces[ spec.ch[0] ];
	      const dpp_filter_t & filt = specs.filters[ spec.band[0] ];
	      const double pad_sec = dpp_filters::pad_seconds( filt , tr.sr );
	      window_result_t wr = get_window( edf , tr , spec.ch[0] , t_tp , spec.window_sec , pad_sec ,
					       qc , qc_flat_th , qc_clip_th , qc_th );
	      if ( wr.ok )
		{
		  std::vector<double> mag_report;
		  filtered_report( tr , wr , filt , &mag_report , NULL );
		  const double m = MiscMath::mean( mag_report );
		  const double sd = mag_report.size() > 1 ? MiscMath::sdev( mag_report , m ) : 0;
		  vals[0] = m; vals[1] = sd; vals[2] = ( m != 0 ) ? sd / m : std::numeric_limits<double>::quiet_NaN();
		}
	    }
	  else if ( spec.ftr == DPP_PLV )
	    {
	      const dpp_trace_t & tr1 = traces[ spec.ch[0] ];
	      const dpp_trace_t & tr2 = traces[ spec.ch[1] ];
	      const dpp_filter_t & filt1 = specs.filters[ spec.band[0] ];
	      const dpp_filter_t & filt2 = specs.filters[ spec.band[1] ];
	      const double pad1 = dpp_filters::pad_seconds( filt1 , tr1.sr );
	      const double pad2 = dpp_filters::pad_seconds( filt2 , tr2.sr );

	      window_result_t wr1 = get_window( edf , tr1 , spec.ch[0] , t_tp , spec.window_sec , pad1 ,
						qc , qc_flat_th , qc_clip_th , qc_th );
	      window_result_t wr2 = get_window( edf , tr2 , spec.ch[1] , t_tp , spec.window_sec , pad2 ,
						qc , qc_flat_th , qc_clip_th , qc_th );

	      if ( wr1.ok && wr2.ok )
		{
		  std::vector<double> mag1, ph1, mag2, ph2;
		  filtered_report( tr1 , wr1 , filt1 , &mag1 , &ph1 );
		  filtered_report( tr2 , wr2 , filt2 , &mag2 , &ph2 );

		  if ( ph1.size() == ph2.size() && ph1.size() > 0 )
		    {
		      ipc_phaseamp_t seed; seed.phase = ph1; seed.amp = mag1;
		      ipc_phaseamp_t tgt;  tgt.phase  = ph2; tgt.amp  = mag2;
		      ipc_param_t ipar; // default: edge_drop_sec = 0 (already trimmed manually above)
		      ipc_output_t out = ipc_t::compute_ipc( seed , tgt , tr1.sr , ipar );
		      vals[0] = out.summary.plv;
		      vals[1] = out.summary.mean_ipc;
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

	  for (int k=0; k<ncol; k++) row_vals.push_back( vals[k] );

	  const std::string root = spec.label_root() + ".w" + Helper::dbl2str( spec.window_sec );
	  writer.level( root , globals::var_strat );
	  for (int k=0; k<ncol; k++)
	    if ( Helper::realnum( vals[k] ) )
	      writer.value( ncol == 1 ? "V" : ( "V" + Helper::int2str( k+1 ) ) , vals[k] );
	  writer.unlevel( globals::var_strat );

	}

      mat.time_sec.push_back( t );
      mat.X.push_back( row_vals );

      writer.unlevel( globals::sec_strat );

    }

  if ( write_binary )
    dpp_io::save( param.value( "data" ) , mat , n_features_total , false );

}
