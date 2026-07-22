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

#include "dsp/plm.h"

#include "param.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "db/db.h"
#include "edf/edf.h"
#include "edf/slice.h"
#include "annot/annot.h"
#include "defs/defs.h"
#include "miscmath/miscmath.h"
#include "dsp/fir.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

extern writer_t writer;
extern logger_t logger;

namespace {

  const double NA = std::numeric_limits<double>::quiet_NaN();
  inline bool is_na( double x ) { return std::isnan( x ); }

  // 1.4826 * median-absolute-deviation (robust SD estimate)
  double mad_scale( const std::vector<double> & v )
  {
    if ( v.empty() ) return 0.0;
    const double med = MiscMath::median( v );
    std::vector<double> dev( v.size() );
    for ( size_t i = 0 ; i < v.size() ; ++i ) dev[i] = std::fabs( v[i] - med );
    return 1.4826 * MiscMath::median( dev );
  }

  // linear interpolation of grid values (at sample positions gx) onto every
  // sample 0..N-1; flat extrapolation beyond the grid ends
  void interp_grid( const std::vector<int> & gx , const std::vector<double> & gy ,
                    std::vector<double> & out , int N )
  {
    out.assign( N , 0.0 );
    if ( gx.empty() ) return;
    if ( gx.size() == 1 ) { out.assign( N , gy[0] ); return; }
    int g = 0;
    for ( int i = 0 ; i < N ; ++i )
      {
        while ( g + 1 < (int)gx.size() - 1 && gx[g+1] < i ) ++g;
        if ( i <= gx.front() ) out[i] = gy.front();
        else if ( i >= gx.back() ) out[i] = gy.back();
        else
          {
            const int x0 = gx[g] , x1 = gx[g+1];
            const double f = ( x1 > x0 ) ? (double)( i - x0 ) / (double)( x1 - x0 ) : 0.0;
            out[i] = gy[g] + f * ( gy[g+1] - gy[g] );
          }
      }
  }

  std::string mkid( char pfx , int n )
  {
    char buf[16];
    std::snprintf( buf , sizeof(buf) , "%c%06d" , pfx , n );
    return std::string( buf );
  }

  double vmean( const std::vector<double> & v ) { return v.empty() ? NA : MiscMath::mean( v ); }
  double vsd  ( const std::vector<double> & v ) { return v.size() < 2 ? NA : MiscMath::sdev( v ); }
  double vmed ( const std::vector<double> & v ) { return v.empty() ? NA : MiscMath::median( v ); }
  double vmin ( const std::vector<double> & v ) { return v.empty() ? NA : *std::min_element( v.begin() , v.end() ); }
  double vmax ( const std::vector<double> & v ) { return v.empty() ? NA : *std::max_element( v.begin() , v.end() ); }

  //
  // staging: build canonical stage intervals directly from existing stage
  // annotations (same non-enforcing idiom as COMBINE-EMG); codes:
  //   0=W 1=N1 2=N2 3=N3 4=R
  //
  struct stg_iv_t { uint64_t a , b; int code; };

  int label2code( const std::string & lab )
  {
    if ( Helper::iequals( lab , "W" )    || Helper::iequals( lab , "wake" )  ) return 0;
    if ( Helper::iequals( lab , "N1" )   || Helper::iequals( lab , "NREM1" ) ) return 1;
    if ( Helper::iequals( lab , "N2" )   || Helper::iequals( lab , "NREM2" ) ) return 2;
    if ( Helper::iequals( lab , "N3" )   || Helper::iequals( lab , "NREM3" ) ) return 3;
    if ( Helper::iequals( lab , "N4" )   || Helper::iequals( lab , "NREM4" ) ) return 3;
    if ( Helper::iequals( lab , "R" )    || Helper::iequals( lab , "REM" )   ) return 4;
    return -1;
  }

  const char * code2stage( int c )
  {
    switch ( c ) { case 0: return "W"; case 1: return "N1"; case 2: return "N2";
                   case 3: return "N3"; case 4: return "R"; default: return "UNK"; }
  }
  const char * code2state( int c )
  {
    if ( c == 0 ) return "W";
    if ( c >= 1 && c <= 4 ) return "S";
    return "UNK";
  }

  struct staging_t {
    std::vector<stg_iv_t> ivs;   // sorted by start
    bool present;
    double hours[5];             // W,N1,N2,N3,R
    staging_t() : present(false) { for (int i=0;i<5;i++) hours[i]=0; }

    void build( edf_t & edf );
    int code_at( uint64_t tp ) const;
  };

  void staging_t::build( edf_t & edf )
  {
    static const char * labels[] = { "W","wake","N1","NREM1","N2","NREM2",
                                     "N3","NREM3","N4","NREM4","R","REM" };
    for ( int i = 0 ; i < 12 ; ++i )
      {
        annot_t * a = edf.annotations->find( labels[i] );
        if ( a == NULL ) continue;
        const int code = label2code( labels[i] );
        if ( code < 0 ) continue;
        annot_map_t evts = a->extract( edf.timeline.wholetrace() );
        for ( annot_map_t::const_iterator aa = evts.begin() ; aa != evts.end() ; ++aa )
          {
            stg_iv_t s;
            s.a = aa->first.interval.start;
            s.b = aa->first.interval.stop;
            s.code = code;
            ivs.push_back( s );
          }
      }
    std::sort( ivs.begin() , ivs.end() ,
               []( const stg_iv_t & x , const stg_iv_t & y ){ return x.a < y.a; } );
    present = ! ivs.empty();
    for ( size_t i = 0 ; i < ivs.size() ; ++i )
      hours[ ivs[i].code ] += ( ivs[i].b - ivs[i].a ) / (double)globals::tp_1sec / 3600.0;
  }

  int staging_t::code_at( uint64_t tp ) const
  {
    if ( ivs.empty() ) return -1;
    // binary search: last interval with start <= tp
    int lo = 0 , hi = (int)ivs.size() - 1 , best = -1;
    while ( lo <= hi )
      {
        const int mid = ( lo + hi ) / 2;
        if ( ivs[mid].a <= tp ) { best = mid; lo = mid + 1; }
        else hi = mid - 1;
      }
    if ( best >= 0 && tp >= ivs[best].a && tp < ivs[best].b ) return ivs[best].code;
    return -1;
  }

  // ensure a signal is in micro-volts; for method=wasm this is required
  void ensure_units( edf_t & edf , int s , bool require_uV )
  {
    const std::string pdim = Helper::trim( edf.header.phys_dimension[s] );
    if ( Helper::imatch( pdim , "uV" ) ) return;
    const bool is_mV = Helper::imatch( pdim , "mV" );
    const bool is_V  = Helper::imatch( pdim , "V" );
    if ( is_mV || is_V ) { edf.rescale( s , "uV" ); return; }
    if ( require_uV )
      Helper::halt( "LM method=wasm requires micro-volt-interpretable signal; '"
                    + edf.header.label[s] + "' has physical dimension '"
                    + ( pdim.empty() ? "(empty)" : pdim ) + "' -- use method=adaptive or rescale first" );
    else
      logger << "  ** warning: " << edf.header.label[s] << " physical dimension '"
             << ( pdim.empty() ? "(empty)" : pdim ) << "' not a recognised voltage unit\n";
  }

} // anonymous namespace


//
// ---------------------------------------------------------------------------
// parameters
// ---------------------------------------------------------------------------
//

void plm::plm_param_t::set_defaults()
{
  detector = DET_WASM;
  onset = 8; offset = 2; morph = 2;
  k_on = 6; k_off = 2; k_morph = 2;
  offset_dur = 0.5; min_dur = 0.5; morph_win = 0.5; max_clm_dur = 10;
  bilateral_gap = 0.5; bilateral_max_components = 4; bilateral_max_dur = 15;
  imi_min = 10; imi_max = 90; min_seq = 4;
  baseline_window = 60; baseline_step = 1; baseline_iter = 3;
  baseline_exclude_k = 3; baseline_exclude_pad = 0.5;
  qc_base_high = 16;
  prefix = "LM";
  do_filter = false; f_lwr = 0; f_upr = 0;
  verbose = false;
}

void plm::plm_param_t::parse( param_t & param )
{
  if ( param.has( "method" ) )
    {
      const std::string m = Helper::toupper( param.value( "method" ) );
      if ( m == "WASM" ) detector = DET_WASM;
      else if ( m == "ADAPTIVE" ) detector = DET_ADAPTIVE;
      else Helper::halt( "LM: method must be 'wasm' or 'adaptive'" );
    }

  if ( param.has( "onset" ) )      onset = param.requires_dbl( "onset" );
  if ( param.has( "offset" ) )     offset = param.requires_dbl( "offset" );
  morph = offset;
  if ( param.has( "morph" ) )      morph = param.requires_dbl( "morph" );

  if ( param.has( "k-on" ) )       k_on = param.requires_dbl( "k-on" );
  if ( param.has( "k-off" ) )      k_off = param.requires_dbl( "k-off" );
  k_morph = k_off;
  if ( param.has( "k-morph" ) )    k_morph = param.requires_dbl( "k-morph" );

  if ( param.has( "offset-dur" ) ) offset_dur = param.requires_dbl( "offset-dur" );
  if ( param.has( "min-dur" ) )    min_dur = param.requires_dbl( "min-dur" );
  if ( param.has( "morph-win" ) )  morph_win = param.requires_dbl( "morph-win" );
  if ( param.has( "max-clm-dur" ) )max_clm_dur = param.requires_dbl( "max-clm-dur" );

  if ( param.has( "bilateral-gap" ) )            bilateral_gap = param.requires_dbl( "bilateral-gap" );
  if ( param.has( "bilateral-max-components" ) ) bilateral_max_components = param.requires_int( "bilateral-max-components" );
  if ( param.has( "bilateral-max-dur" ) )        bilateral_max_dur = param.requires_dbl( "bilateral-max-dur" );

  if ( param.has( "imi-min" ) )    imi_min = param.requires_dbl( "imi-min" );
  if ( param.has( "imi-max" ) )    imi_max = param.requires_dbl( "imi-max" );
  if ( param.has( "min-seq" ) )    min_seq = param.requires_int( "min-seq" );

  if ( param.has( "baseline-window" ) )      baseline_window = param.requires_dbl( "baseline-window" );
  if ( param.has( "baseline-step" ) )        baseline_step = param.requires_dbl( "baseline-step" );
  if ( param.has( "baseline-iter" ) )        baseline_iter = param.requires_int( "baseline-iter" );
  if ( param.has( "baseline-exclude-k" ) )   baseline_exclude_k = param.requires_dbl( "baseline-exclude-k" );
  if ( param.has( "baseline-exclude-pad" ) ) baseline_exclude_pad = param.requires_dbl( "baseline-exclude-pad" );

  if ( param.has( "qc-base-high" ) ) qc_base_high = param.requires_dbl( "qc-base-high" );

  if ( param.has( "prefix" ) )     prefix = param.value( "prefix" );

  if ( param.has( "f-lwr" ) ) { f_lwr = param.requires_dbl( "f-lwr" ); do_filter = true; }
  if ( param.has( "f-upr" ) ) { f_upr = param.requires_dbl( "f-upr" ); do_filter = true; }

  verbose = param.has( "verbose" );
}

const char * plm::side_str( plm::side_t s )
{
  switch ( s ) { case SIDE_L: return "L"; case SIDE_R: return "R";
                 case SIDE_B: return "B"; default: return "LEG"; }
}


//
// ---------------------------------------------------------------------------
// dynamic iterative baseline / robust-scale estimator (shared by both methods)
// ---------------------------------------------------------------------------
//

namespace {

  void compute_baseline( const std::vector<double> & y , double fs ,
                         const plm::plm_param_t & P ,
                         std::vector<double> & B , std::vector<double> & S ,
                         bool * scale_ok )
  {
    const int N = (int)y.size();
    B.assign( N , 0.0 ); S.assign( N , 0.0 );
    *scale_ok = true;
    if ( N == 0 ) return;

    // recording-wide robust scale (fallback)
    double global_scale = 0.0;
    {
      std::vector<double> fin;
      fin.reserve( N );
      for ( int i = 0 ; i < N ; ++i ) if ( std::isfinite( y[i] ) ) fin.push_back( y[i] );
      global_scale = mad_scale( fin );
    }
    if ( global_scale <= 0 ) *scale_ok = false;

    const int half = std::max( 1 , (int)std::lround( P.baseline_window * fs / 2.0 ) );
    const int step = std::max( 1 , (int)std::lround( P.baseline_step * fs ) );
    const int pad  = std::max( 0 , (int)std::lround( P.baseline_exclude_pad * fs ) );

    // regular grid (+ final sample)
    std::vector<int> gx;
    for ( int g = 0 ; g < N ; g += step ) gx.push_back( g );
    if ( gx.empty() || gx.back() != N - 1 ) gx.push_back( N - 1 );
    const int G = (int)gx.size();

    std::vector<char> excluded( N , 0 );
    std::vector<double> Bg( G , 0.0 ) , Sg( G , 0.0 );
    std::vector<double> win; win.reserve( 2 * half + 1 );

    // one estimation pass over the grid, honouring the current exclusion mask
    auto estimate = [&]()
    {
      for ( int gi = 0 ; gi < G ; ++gi )
        {
          const int c = gx[gi];
          const int lo = std::max( 0 , c - half ) , hi = std::min( N , c + half );
          win.clear();
          for ( int i = lo ; i < hi ; ++i )
            if ( ! excluded[i] && std::isfinite( y[i] ) ) win.push_back( y[i] );
          // if exclusion emptied the window, fall back to all finite samples
          if ( win.empty() )
            for ( int i = lo ; i < hi ; ++i )
              if ( std::isfinite( y[i] ) ) win.push_back( y[i] );
          if ( win.empty() ) { Bg[gi] = 0.0; Sg[gi] = NA; continue; }
          Bg[gi] = MiscMath::median( win );
          const double s = mad_scale( win );
          Sg[gi] = ( s > 0 ) ? s : NA;   // flag negligible scale for repair
        }

      // repair negligible/undefined scales: interpolate from valid grid
      // neighbours, else recording-wide robust scale
      std::vector<int> vx; std::vector<double> vy;
      for ( int gi = 0 ; gi < G ; ++gi )
        if ( ! is_na( Sg[gi] ) ) { vx.push_back( gx[gi] ); vy.push_back( Sg[gi] ); }
      if ( vx.empty() )
        {
          for ( int gi = 0 ; gi < G ; ++gi ) Sg[gi] = ( global_scale > 0 ? global_scale : 1.0 );
        }
      else
        {
          std::vector<double> rep;
          interp_grid( vx , vy , rep , N );
          for ( int gi = 0 ; gi < G ; ++gi )
            if ( is_na( Sg[gi] ) ) Sg[gi] = rep[ gx[gi] ];
        }
    };

    // provisional estimate, then baseline_iter rounds of excursion-exclusion
    estimate();
    interp_grid( gx , Bg , B , N );
    interp_grid( gx , Sg , S , N );

    for ( int it = 0 ; it < P.baseline_iter ; ++it )
      {
        // flag provisional excursions y > B + k*S
        std::vector<char> ex( N , 0 );
        for ( int i = 0 ; i < N ; ++i )
          if ( std::isfinite( y[i] ) && y[i] > B[i] + P.baseline_exclude_k * S[i] )
            ex[i] = 1;
        // dilate by +/- pad
        excluded.assign( N , 0 );
        for ( int i = 0 ; i < N ; ++i )
          if ( ex[i] )
            for ( int j = std::max( 0 , i - pad ) ; j < std::min( N , i + pad + 1 ) ; ++j )
              excluded[j] = 1;

        estimate();
        interp_grid( gx , Bg , B , N );
        interp_grid( gx , Sg , S , N );
      }

    // guarantee strictly positive scale everywhere
    for ( int i = 0 ; i < N ; ++i )
      if ( ! ( S[i] > 0 ) )
        S[i] = ( global_scale > 0 ? global_scale : 1.0 );
  }

  //
  // low-level LM detector on one rectified channel; appends to 'out'
  //
  void detect_components( const std::vector<double> & y ,
                          const std::vector<uint64_t> & tp ,
                          double fs ,
                          const std::vector<double> & B ,
                          const std::vector<double> & S ,
                          plm::side_t side , const std::string & siglabel ,
                          const plm::plm_param_t & P ,
                          int & comp_counter ,
                          std::vector<plm::lm_component_t> & out )
  {
    const int N = (int)y.size();
    if ( N == 0 ) return;
    const bool wasm = ( P.detector == plm::DET_WASM );

    const int off_len   = std::max( 1 , (int)std::lround( P.offset_dur * fs ) );
    const int morph_len = std::max( 1 , (int)std::lround( P.morph_win * fs ) );
    const uint64_t sample_tp = (uint64_t)std::llround( globals::tp_1sec / fs );

    auto on_thr  = [&]( int i ){ return wasm ? B[i] + P.onset  : B[i] + P.k_on  * S[i]; };
    auto off_thr = [&]( int i ){ return wasm ? B[i] + P.offset : B[i] + P.k_off * S[i]; };

    int i = 0;
    while ( i < N )
      {
        if ( ! ( std::isfinite( y[i] ) && y[i] >= on_thr( i ) ) ) { ++i; continue; }

        const int onset = i;

        // offset = first sample beginning a continuous >= offset_dur below-run
        int offset = -1 , run_start = -1 , run_len = 0;
        for ( int j = onset ; j < N ; ++j )
          {
            if ( std::isfinite( y[j] ) && y[j] < off_thr( j ) )
              {
                if ( run_len == 0 ) run_start = j;
                if ( ++run_len >= off_len ) { offset = run_start; break; }
              }
            else { run_len = 0; run_start = -1; }
          }
        if ( offset < 0 ) offset = N;   // never returns to baseline: extends to end

        const uint64_t onset_tp  = tp[ onset ];
        const uint64_t offset_tp = ( offset < N ) ? tp[ offset ] : ( tp[ N - 1 ] + sample_tp );
        const double dur = ( offset_tp - onset_tp ) / (double)globals::tp_1sec;

        // event-window metrics
        double peak = 0 , peak_exc = 0 , peak_z = 0;
        bool first = true;
        for ( int k = onset ; k < offset ; ++k )
          {
            if ( ! std::isfinite( y[k] ) ) continue;
            const double exc = y[k] - B[k];
            const double z   = exc / S[k];
            if ( first ) { peak = y[k]; peak_exc = exc; peak_z = z; first = false; }
            else { peak = std::max(peak,y[k]); peak_exc = std::max(peak_exc,exc); peak_z = std::max(peak_z,z); }
          }

        // morphology: max 0.5s sliding-window median excursion (uV or z)
        double morph_value = 0.0; bool morph_ok = false; bool mfirst = true;
        std::vector<double> mw; mw.reserve( morph_len );
        for ( int w = onset ; w + morph_len <= offset ; ++w )
          {
            mw.clear();
            for ( int k = w ; k < w + morph_len ; ++k )
              {
                if ( ! std::isfinite( y[k] ) ) continue;
                mw.push_back( wasm ? ( y[k] - B[k] ) : ( ( y[k] - B[k] ) / S[k] ) );
              }
            if ( mw.empty() ) continue;
            const double m = MiscMath::median( mw );
            if ( mfirst || m > morph_value ) { morph_value = m; mfirst = false; }
          }
        const double morph_thr = wasm ? P.morph : P.k_morph;
        morph_ok = ( ! mfirst ) && ( morph_value >= morph_thr );

        // advance scan position past the event
        i = ( offset > onset ) ? offset : onset + 1;

        // apply minimum-duration + morphology gates
        if ( dur < P.min_dur || ! morph_ok ) continue;

        plm::lm_component_t c;
        c.comp_id = ++comp_counter;
        c.side = side;
        c.sig = siglabel;
        c.onset_tp = onset_tp;
        c.offset_tp = offset_tp;
        c.onset_sec = onset_tp / (double)globals::tp_1sec;
        c.dur = dur;
        c.base_onset = B[ onset ];
        c.scale_onset = S[ onset ];
        c.on_thr = on_thr( onset );
        c.off_thr = off_thr( onset );
        c.peak = peak;
        c.peak_exc = peak_exc;
        c.peak_z = peak_z;
        c.morph_value = morph_value;
        c.morph_ok = morph_ok;
        c.clm_component = ( dur >= P.min_dur && dur <= P.max_clm_dur );
        c.stage = "UNK"; c.state = "UNK";
        c.qc_base_high = wasm && ( B[ onset ] > P.qc_base_high );
        out.push_back( c );
      }
  }

} // anonymous namespace


//
// ---------------------------------------------------------------------------
// pure WASM 2016 grammar
// ---------------------------------------------------------------------------
//

std::vector<plm::lm_event_t>
plm::combine_single( std::vector<plm::lm_component_t> & comps , const plm::plm_param_t & P )
{
  std::sort( comps.begin() , comps.end() ,
             []( const lm_component_t & a , const lm_component_t & b ){ return a.onset_tp < b.onset_tp; } );

  std::vector<lm_event_t> ev;
  int eid = 0;
  for ( size_t i = 0 ; i < comps.size() ; ++i )
    {
      const lm_component_t & c = comps[i];
      lm_event_t e;
      e.evt_id = ++eid;
      e.side = c.side;                 // SIDE_LEG (or L/R if used monolaterally)
      e.onset_tp = c.onset_tp; e.offset_tp = c.offset_tp;
      e.onset_sec = c.onset_sec; e.dur = c.dur;
      e.comp_ids.push_back( c.comp_id );
      e.n_comp = 1;
      e.n_l = ( c.side == SIDE_L ) ? 1 : 0;
      e.n_r = ( c.side == SIDE_R ) ? 1 : 0;
      e.all_comps_clm = c.clm_component;
      e.clm = false; e.plm = false;
      e.seq = 0; e.seq_pos = 0; e.seq_n = 0;
      e.imi_prev = NA; e.imi_next = NA;
      e.stage = c.stage; e.state = c.state;
      e.peak_exc_max = c.peak_exc; e.peak_z_max = c.peak_z;
      e.qc_base_high = c.qc_base_high;
      ev.push_back( e );
    }
  return ev;
}

std::vector<plm::lm_event_t>
plm::combine_bilateral( std::vector<plm::lm_component_t> & comps , const plm::plm_param_t & P )
{
  const int n = (int)comps.size();
  std::vector<lm_event_t> ev;
  if ( n == 0 ) return ev;

  std::sort( comps.begin() , comps.end() ,
             []( const lm_component_t & a , const lm_component_t & b ){ return a.onset_tp < b.onset_tp; } );

  const uint64_t gap_tp = (uint64_t)std::llround( P.bilateral_gap * globals::tp_1sec );

  // union-find over the cross-leg near-overlap graph (edges only connect
  // opposite sides; strict gap < bilateral-gap)
  std::vector<int> uf( n );
  for ( int i = 0 ; i < n ; ++i ) uf[i] = i;
  std::function<int(int)> find = [&]( int x ){ while ( uf[x] != x ) { uf[x] = uf[uf[x]]; x = uf[x]; } return x; };
  auto unite = [&]( int a , int b ){ uf[ find(a) ] = find(b); };

  for ( int i = 0 ; i < n ; ++i )
    for ( int j = i + 1 ; j < n ; ++j )
      {
        // sorted by onset: once comps[j] starts at/after comps[i] ends + gap,
        // no later j can connect to i either
        if ( comps[j].onset_tp >= comps[i].offset_tp + gap_tp ) break;
        if ( comps[i].side == comps[j].side ) continue;   // same leg: no edge
        // strict gap < bilateral-gap (overlap => gap 0); the break above
        // already enforces onset_tp < offset_tp + gap_tp
        unite( i , j );
      }

  // gather connected components
  std::map<int,std::vector<int> > groups;
  for ( int i = 0 ; i < n ; ++i ) groups[ find(i) ].push_back( i );

  int eid = 0;
  for ( std::map<int,std::vector<int> >::iterator g = groups.begin() ; g != groups.end() ; ++g )
    {
      const std::vector<int> & idx = g->second;
      lm_event_t e;
      e.evt_id = ++eid;
      e.onset_tp = comps[ idx[0] ].onset_tp;
      e.offset_tp = comps[ idx[0] ].offset_tp;
      int n_l = 0 , n_r = 0;
      bool all_clm = true;
      double pemax = 0 , pzmax = 0; bool pfirst = true;
      bool qc = false;
      int earliest = idx[0];
      for ( size_t k = 0 ; k < idx.size() ; ++k )
        {
          const lm_component_t & c = comps[ idx[k] ];
          e.comp_ids.push_back( c.comp_id );
          if ( c.onset_tp < e.onset_tp ) e.onset_tp = c.onset_tp;
          if ( c.offset_tp > e.offset_tp ) e.offset_tp = c.offset_tp;
          if ( c.onset_tp < comps[ earliest ].onset_tp ) earliest = idx[k];
          if ( c.side == SIDE_L ) ++n_l; else if ( c.side == SIDE_R ) ++n_r;
          if ( ! c.clm_component ) all_clm = false;
          if ( pfirst ) { pemax = c.peak_exc; pzmax = c.peak_z; pfirst = false; }
          else { pemax = std::max( pemax , c.peak_exc ); pzmax = std::max( pzmax , c.peak_z ); }
          if ( c.qc_base_high ) qc = true;
        }
      std::sort( e.comp_ids.begin() , e.comp_ids.end() );
      e.n_comp = (int)idx.size();
      e.n_l = n_l; e.n_r = n_r;
      e.side = ( n_l > 0 && n_r > 0 ) ? SIDE_B : ( n_l > 0 ? SIDE_L : SIDE_R );
      e.all_comps_clm = all_clm;
      e.onset_sec = e.onset_tp / (double)globals::tp_1sec;
      e.dur = ( e.offset_tp - e.onset_tp ) / (double)globals::tp_1sec;
      e.clm = false; e.plm = false;
      e.seq = 0; e.seq_pos = 0; e.seq_n = 0;
      e.imi_prev = NA; e.imi_next = NA;
      e.stage = comps[ earliest ].stage; e.state = comps[ earliest ].state;
      e.peak_exc_max = pemax; e.peak_z_max = pzmax;
      e.qc_base_high = qc;
      ev.push_back( e );
    }

  std::sort( ev.begin() , ev.end() ,
             []( const lm_event_t & a , const lm_event_t & b ){ return a.onset_tp < b.onset_tp; } );
  for ( size_t i = 0 ; i < ev.size() ; ++i ) ev[i].evt_id = (int)i + 1;
  return ev;
}

void plm::classify_clm( std::vector<plm::lm_event_t> & events , const plm::plm_param_t & P )
{
  for ( size_t i = 0 ; i < events.size() ; ++i )
    {
      lm_event_t & e = events[i];
      if ( e.n_comp <= 1 )
        e.clm = e.all_comps_clm;    // monolateral: CLM iff its (single) component is CLM
      else
        e.clm = e.all_comps_clm
                && e.n_comp <= P.bilateral_max_components
                && e.dur <= P.bilateral_max_dur;
    }
}

std::vector<plm::plm_sequence_t>
plm::detect_sequences( std::vector<plm::lm_event_t> & events , const plm::plm_param_t & P )
{
  std::sort( events.begin() , events.end() ,
             []( const lm_event_t & a , const lm_event_t & b ){ return a.onset_tp < b.onset_tp; } );
  const int n = (int)events.size();
  for ( int i = 0 ; i < n ; ++i ) { events[i].imi_prev = NA; events[i].imi_next = NA;
                                    events[i].plm = false; events[i].seq = 0;
                                    events[i].seq_pos = 0; events[i].seq_n = 0; }

  // measurable consecutive-CLM IMIs (no intervening non-CLM LM)
  int prev = -1;
  for ( int i = 0 ; i < n ; ++i )
    {
      if ( events[i].clm )
        {
          if ( prev >= 0 )
            {
              const double imi = events[i].onset_sec - events[prev].onset_sec;
              events[i].imi_prev = imi;
              events[prev].imi_next = imi;
            }
          prev = i;
        }
      else prev = -1;  // non-CLM LM breaks measurability
    }

  std::vector<plm_sequence_t> seqs;
  std::vector<int> run;

  auto close_run = [&]()
  {
    if ( (int)run.size() >= P.min_seq )
      {
        plm_sequence_t s;
        s.seq_id = (int)seqs.size() + 1;
        s.n = (int)run.size();
        s.onset_tp = events[ run.front() ].onset_tp;
        s.offset_tp = events[ run.back() ].offset_tp;
        s.dur = ( s.offset_tp - s.onset_tp ) / (double)globals::tp_1sec;
        s.n_l = s.n_r = s.n_b = s.n_leg = 0;
        for ( size_t t = 0 ; t < run.size() ; ++t )
          {
            lm_event_t & e = events[ run[t] ];
            e.plm = true; e.seq = s.seq_id; e.seq_pos = (int)t + 1; e.seq_n = s.n;
            switch ( e.side ) { case SIDE_L: ++s.n_l; break; case SIDE_R: ++s.n_r; break;
                                case SIDE_B: ++s.n_b; break; default: ++s.n_leg; break; }
            if ( t > 0 )
              s.imis.push_back( events[ run[t] ].onset_sec - events[ run[t-1] ].onset_sec );
          }
        seqs.push_back( s );
      }
    run.clear();
  };

  for ( int i = 0 ; i < n ; ++i )
    {
      if ( events[i].clm )
        {
          if ( run.empty() ) run.push_back( i );
          else
            {
              const double imi = events[i].onset_sec - events[ run.back() ].onset_sec;
              if ( imi >= P.imi_min && imi <= P.imi_max ) run.push_back( i );
              else { close_run(); run.push_back( i ); }
            }
        }
      else
        close_run();   // intervening non-CLM LM breaks the sequence
    }
  close_run();

  return seqs;
}


//
// ---------------------------------------------------------------------------
// output helpers
// ---------------------------------------------------------------------------
//

namespace {

  using namespace plm;

  const char * detector_str( const plm_param_t & P ) { return P.detector == DET_WASM ? "WASM" : "ADAPTIVE"; }

  std::string join_ids( const std::vector<int> & v )
  {
    std::string s;
    for ( size_t i = 0 ; i < v.size() ; ++i ) { if ( i ) s += ","; s += mkid( 'C' , v[i] ); }
    return s;
  }

  // attach final-event metadata (shared by LM/CLM/PLM annotations & verbose)
  void set_event_meta( instance_t * inst , const lm_event_t & e ,
                       const plm_param_t & P , const std::string & mode )
  {
    inst->set( "ID" , mkid( 'E' , e.evt_id ) );
    inst->set( "DETECTOR" , std::string( detector_str( P ) ) );
    inst->set( "GRAMMAR" , std::string( "WASM2016" ) );
    inst->set( "MODE" , mode );
    inst->set( "SIDE" , std::string( side_str( e.side ) ) );
    inst->set( "DUR" , e.dur );
    inst->set( "NCOMP" , e.n_comp );
    inst->set( "COMP_IDS" , join_ids( e.comp_ids ) );
    inst->set( "N_L" , e.n_l );
    inst->set( "N_R" , e.n_r );
    inst->set( "CLM" , (int)( e.clm ? 1 : 0 ) );
    inst->set( "PLM" , (int)( e.plm ? 1 : 0 ) );
    if ( e.seq > 0 ) { inst->set( "SEQ" , e.seq ); inst->set( "SEQ_POS" , e.seq_pos ); inst->set( "SEQ_N" , e.seq_n ); }
    if ( ! is_na( e.imi_prev ) ) inst->set( "IMI_PREV" , e.imi_prev );
    if ( ! is_na( e.imi_next ) ) inst->set( "IMI_NEXT" , e.imi_next );
    inst->set( "STAGE" , e.stage );
    inst->set( "STATE" , e.state );
    inst->set( "PEAK_EXC_MAX" , e.peak_exc_max );
    inst->set( "PEAK_Z_MAX" , e.peak_z_max );
    inst->set( "QC_BASE_HIGH" , (int)( e.qc_base_high ? 1 : 0 ) );
  }

  void write_event_vars( const lm_event_t & e , const plm_param_t & P , const std::string & mode )
  {
    writer.value( "DETECTOR" , std::string( detector_str( P ) ) );
    writer.value( "GRAMMAR" , std::string( "WASM2016" ) );
    writer.value( "MODE" , mode );
    writer.value( "SIDE" , std::string( side_str( e.side ) ) );
    writer.value( "START" , e.onset_tp / (double)globals::tp_1sec );
    writer.value( "STOP" , e.offset_tp / (double)globals::tp_1sec );
    writer.value( "DUR" , e.dur );
    writer.value( "NCOMP" , e.n_comp );
    writer.value( "COMP_IDS" , join_ids( e.comp_ids ) );
    writer.value( "N_L" , e.n_l );
    writer.value( "N_R" , e.n_r );
    writer.value( "CLM" , (int)( e.clm ? 1 : 0 ) );
    writer.value( "PLM" , (int)( e.plm ? 1 : 0 ) );
    if ( e.seq > 0 ) { writer.value( "SEQ" , e.seq ); writer.value( "SEQ_POS" , e.seq_pos ); writer.value( "SEQ_N" , e.seq_n ); }
    if ( ! is_na( e.imi_prev ) ) writer.value( "IMI_PREV" , e.imi_prev );
    if ( ! is_na( e.imi_next ) ) writer.value( "IMI_NEXT" , e.imi_next );
    writer.value( "STAGE" , e.stage );
    writer.value( "STATE" , e.state );
    writer.value( "PEAK_EXC_MAX" , e.peak_exc_max );
    writer.value( "PEAK_Z_MAX" , e.peak_z_max );
    writer.value( "QC_BASE_HIGH" , (int)( e.qc_base_high ? 1 : 0 ) );
  }

  // membership predicate for a stage-stratum given an event's stage/state
  bool in_stratum( const lm_event_t & e , const std::string & ss )
  {
    if ( ss == "W" ) return e.stage == "W";
    if ( ss == "S" ) return e.state == "S";
    if ( ss == "NREM" ) return e.stage == "N1" || e.stage == "N2" || e.stage == "N3";
    if ( ss == "REM" ) return e.stage == "R";
    return e.stage == ss;   // N1/N2/N3
  }

  double stratum_hours( const staging_t & stg , const std::string & ss )
  {
    if ( ss == "W" ) return stg.hours[0];
    if ( ss == "N1" ) return stg.hours[1];
    if ( ss == "N2" ) return stg.hours[2];
    if ( ss == "N3" ) return stg.hours[3];
    if ( ss == "REM" ) return stg.hours[4];
    if ( ss == "NREM" ) return stg.hours[1] + stg.hours[2] + stg.hours[3];
    if ( ss == "S" ) return stg.hours[1] + stg.hours[2] + stg.hours[3] + stg.hours[4];
    return 0.0;
  }

} // anonymous namespace


//
// ---------------------------------------------------------------------------
// command driver
// ---------------------------------------------------------------------------
//

leg_movements_t::leg_movements_t( edf_t & edf , param_t & param )
{
  using namespace plm;

  plm_param_t P;
  P.parse( param );

  //
  // signals: exactly 1 (LEG) or 2 (LAT,RAT)
  //

  const std::string sigstr = param.requires( "sig" );
  signal_list_t signals = edf.header.signal_list( sigstr );
  const int ns = signals.size();
  if ( ns < 1 || ns > 2 )
    Helper::halt( "LM requires exactly one signal (LEG) or two signals (LAT,RAT); got " + Helper::int2str( ns ) );

  const bool lr_mode = ( ns == 2 );
  const std::string mode = lr_mode ? "LR" : "LEG";

  logger << "\n  LM : leg-movement / PLM detection (WASM 2016 grammar)\n";
  logger << "    method   : " << ( P.detector == DET_WASM ? "wasm (fixed uV excursions)" : "adaptive (robust local scale)" ) << "\n";
  logger << "    mode     : " << ( lr_mode ? "LR (2 signals: LAT,RAT)" : "LEG (1 signal)" ) << "\n";

  //
  // staging (non-enforcing)
  //

  staging_t stg;
  stg.build( edf );
  logger << "    staging  : " << ( stg.present ? "found" : "not found (whole-recording metrics only)" ) << "\n";

  //
  // per-channel detection
  //

  std::vector<lm_component_t> comps;
  int comp_counter = 0;
  double fs_min = 0 , fs_max = 0;
  int n_comp_side[2] = { 0 , 0 };   // L , R (or LEG in [0])

  for ( int s = 0 ; s < ns ; ++s )
    {
      const int slot = signals( s );
      const double fs = edf.header.sampling_freq( slot );
      if ( fs <= 0 ) Helper::halt( "LM: invalid sample rate for " + signals.label(s) );
      if ( s == 0 ) { fs_min = fs_max = fs; } else { fs_min = std::min(fs_min,fs); fs_max = std::max(fs_max,fs); }

      ensure_units( edf , slot , P.detector == DET_WASM );

      // working copy (never modifies the EDF signal)
      slice_t slice( edf , slot , edf.timeline.wholetrace() );
      const std::vector<double> * pd = slice.pdata();
      const std::vector<uint64_t> * ptp = slice.ptimepoints();
      if ( pd == NULL || ptp == NULL || pd->empty() )
        { logger << "  ** no data for " << signals.label(s) << "\n"; continue; }

      std::vector<double> y = *pd;
      const std::vector<uint64_t> & tp = *ptp;

      // optional pre-filter (reuse Luna FIR); band/low/high per f-lwr,f-upr
      if ( P.do_filter )
        {
          std::vector<double> ripple( 1 , 0.01 ) , tw( 1 , 2.0 );
          const bool hp = P.f_lwr > 0 , lp = P.f_upr > 0;
          if ( hp && lp )
            y = dsptools::apply_fir( y , (int)fs , fir_t::BAND_PASS , 1 , ripple , tw , P.f_lwr , P.f_upr , 0 , fir_t::HAMMING , true );
          else if ( hp )
            y = dsptools::apply_fir( y , (int)fs , fir_t::HIGH_PASS , 1 , ripple , tw , P.f_lwr , 0 , 0 , fir_t::HAMMING , true );
          else if ( lp )
            y = dsptools::apply_fir( y , (int)fs , fir_t::LOW_PASS , 1 , ripple , tw , P.f_upr , 0 , 0 , fir_t::HAMMING , true );
        }

      // full-wave rectify
      for ( size_t i = 0 ; i < y.size() ; ++i ) y[i] = std::fabs( y[i] );

      // shared dynamic baseline / robust scale
      std::vector<double> B , S; bool scale_ok = true;
      compute_baseline( y , fs , P , B , S , &scale_ok );
      if ( P.detector == DET_ADAPTIVE && ! scale_ok )
        Helper::halt( "LM method=adaptive: could not estimate a robust scale for " + signals.label(s)
                      + " (signal appears constant)" );

      side_t side = lr_mode ? ( s == 0 ? SIDE_L : SIDE_R ) : SIDE_LEG;

      const int before = comp_counter;
      detect_components( y , tp , fs , B , S , side , signals.label(s) , P , comp_counter , comps );
      const int nc = comp_counter - before;
      if ( lr_mode ) n_comp_side[ s ] += nc; else n_comp_side[0] += nc;

      logger << "    " << signals.label(s) << " : " << nc << " component event(s)\n";
    }

  // assign stage at each component onset
  for ( size_t i = 0 ; i < comps.size() ; ++i )
    {
      const int code = stg.code_at( comps[i].onset_tp );
      comps[i].stage = code2stage( code );
      comps[i].state = code2state( code );
    }

  //
  // final evaluation-unit LMs
  //

  std::vector<lm_event_t> events = lr_mode ? combine_bilateral( comps , P ) : combine_single( comps , P );

  // re-assign event stage from onset (combine used earliest component's stage;
  // ensure consistency for single-channel too)
  for ( size_t i = 0 ; i < events.size() ; ++i )
    {
      const int code = stg.code_at( events[i].onset_tp );
      events[i].stage = code2stage( code );
      events[i].state = code2state( code );
    }

  classify_clm( events , P );
  std::vector<plm_sequence_t> seqs = detect_sequences( events , P );

  //
  // IMI statistics (measurable consecutive-CLM IMIs)
  //

  std::vector<double> all_imi , periodic_imi , log_periodic_imi;
  int n_short = 0 , n_long = 0;
  for ( size_t i = 0 ; i < events.size() ; ++i )
    {
      if ( is_na( events[i].imi_next ) ) continue;   // count each pair once
      const double imi = events[i].imi_next;
      all_imi.push_back( imi );
      if ( imi < P.imi_min ) ++n_short;
      else if ( imi > P.imi_max ) ++n_long;
      else { periodic_imi.push_back( imi ); if ( imi > 0 ) log_periodic_imi.push_back( std::log( imi ) ); }
    }

  //
  // counts
  //

  int n_lm = (int)events.size() , n_clm = 0 , n_plm = 0;
  int n_lm_side[3] = {0,0,0} , n_clm_side[3] = {0,0,0} , n_plm_side[3] = {0,0,0}; // L,R,B
  for ( size_t i = 0 ; i < events.size() ; ++i )
    {
      const lm_event_t & e = events[i];
      if ( e.clm ) ++n_clm;
      if ( e.plm ) ++n_plm;
      int si = -1;
      if ( e.side == SIDE_L ) si = 0; else if ( e.side == SIDE_R ) si = 1; else if ( e.side == SIDE_B ) si = 2;
      if ( si >= 0 ) { ++n_lm_side[si]; if ( e.clm ) ++n_clm_side[si]; if ( e.plm ) ++n_plm_side[si]; }
    }

  //
  // -----------------------------------------------------------------------
  // annotations
  // -----------------------------------------------------------------------
  //

  const std::string px = P.prefix;
  const std::string comp_L   = px + "_L";
  const std::string comp_R   = px + "_R";
  const std::string comp_LEG = px + "_LEG";
  const std::string lm_class  = ( px == "LM" ) ? "LM"      : px + "_LM";
  const std::string clm_class = ( px == "LM" ) ? "CLM"     : px + "_CLM";
  const std::string plm_class = ( px == "LM" ) ? "PLM"     : px + "_PLM";
  const std::string seq_class = ( px == "LM" ) ? "PLM_SEQ" : px + "_PLM_SEQ";

  annot_t * a_compL   = lr_mode ? edf.annotations->add( comp_L )   : NULL;
  annot_t * a_compR   = lr_mode ? edf.annotations->add( comp_R )   : NULL;
  annot_t * a_compLEG = lr_mode ? NULL : edf.annotations->add( comp_LEG );
  annot_t * a_lm  = edf.annotations->add( lm_class );
  annot_t * a_clm = edf.annotations->add( clm_class );
  annot_t * a_plm = edf.annotations->add( plm_class );
  annot_t * a_seq = edf.annotations->add( seq_class );

  // component annotations (auditability)
  for ( size_t i = 0 ; i < comps.size() ; ++i )
    {
      const lm_component_t & c = comps[i];
      annot_t * a = ( c.side == SIDE_L ) ? a_compL : ( c.side == SIDE_R ) ? a_compR : a_compLEG;
      if ( a == NULL ) continue;
      instance_t * inst = a->add( mkid( 'C' , c.comp_id ) , interval_t( c.onset_tp , c.offset_tp ) , "." );
      inst->set( "ID" , mkid( 'C' , c.comp_id ) );
      inst->set( "SIDE" , std::string( side_str( c.side ) ) );
      inst->set( "SIG" , c.sig );
      inst->set( "DETECTOR" , std::string( detector_str( P ) ) );
      inst->set( "GRAMMAR" , std::string( "WASM2016" ) );
      inst->set( "DUR" , c.dur );
      inst->set( "BASE" , c.base_onset );
      inst->set( "SCALE" , c.scale_onset );
      inst->set( "ON_THR" , c.on_thr );
      inst->set( "OFF_THR" , c.off_thr );
      inst->set( "PEAK" , c.peak );
      inst->set( "PEAK_EXC" , c.peak_exc );
      inst->set( "PEAK_Z" , c.peak_z );
      inst->set( "MORPH" , c.morph_value );
      inst->set( "CLM" , (int)( c.clm_component ? 1 : 0 ) );
      inst->set( "STAGE" , c.stage );
      inst->set( "STATE" , c.state );
      inst->set( "QC_BASE_HIGH" , (int)( c.qc_base_high ? 1 : 0 ) );
    }

  // final LM (and CLM/PLM subsets) — copy full event metadata to each
  for ( size_t i = 0 ; i < events.size() ; ++i )
    {
      const lm_event_t & e = events[i];
      const interval_t iv( e.onset_tp , e.offset_tp );
      const std::string eid = mkid( 'E' , e.evt_id );

      instance_t * i1 = a_lm->add( eid , iv , "." );
      set_event_meta( i1 , e , P , mode );

      if ( e.clm )
        { instance_t * i2 = a_clm->add( eid , iv , "." ); set_event_meta( i2 , e , P , mode ); }
      if ( e.plm )
        { instance_t * i3 = a_plm->add( eid , iv , "." ); set_event_meta( i3 , e , P , mode ); }
    }

  // PLM sequences
  for ( size_t i = 0 ; i < seqs.size() ; ++i )
    {
      const plm_sequence_t & s = seqs[i];
      instance_t * inst = a_seq->add( mkid( 'P' , s.seq_id ) , interval_t( s.onset_tp , s.offset_tp ) , "." );
      inst->set( "SEQ" , s.seq_id );
      inst->set( "N" , s.n );
      inst->set( "DUR" , s.dur );
      inst->set( "IMI_N" , (int)s.imis.size() );
      if ( ! s.imis.empty() )
        {
          inst->set( "IMI_MEAN" , vmean( s.imis ) );
          if ( s.imis.size() >= 2 ) inst->set( "IMI_SD" , vsd( s.imis ) );
          inst->set( "IMI_MIN" , vmin( s.imis ) );
          inst->set( "IMI_MAX" , vmax( s.imis ) );
        }
      inst->set( "N_L" , s.n_l );
      inst->set( "N_R" , s.n_r );
      inst->set( "N_B" , s.n_b );
      inst->set( "N_LEG" , s.n_leg );
      inst->set( "DETECTOR" , std::string( detector_str( P ) ) );
      inst->set( "GRAMMAR" , std::string( "WASM2016" ) );
      inst->set( "MODE" , mode );
    }

  //
  // -----------------------------------------------------------------------
  // normal output: QC / method
  // -----------------------------------------------------------------------
  //

  writer.value( "MODE" , mode );
  writer.value( "METHOD" , std::string( P.detector == DET_WASM ? "wasm" : "adaptive" ) );
  writer.value( "N_SIG" , ns );
  writer.value( "FS_MIN" , fs_min );
  writer.value( "FS_MAX" , fs_max );
  writer.value( "LOW_SR" , (int)( fs_min < 100 ? 1 : 0 ) );
  writer.value( "NONSTANDARD_BW" , (int)( ( P.do_filter && P.f_lwr == 10 && P.f_upr == 100 ) ? 0 : 1 ) );
  writer.value( "NO_STAGING" , (int)( stg.present ? 0 : 1 ) );

  // overall counts
  writer.value( "N_LM" , n_lm );
  writer.value( "N_CLM" , n_clm );
  writer.value( "N_PLM" , n_plm );
  writer.value( "N_PLM_SEQ" , (int)seqs.size() );
  if ( n_clm > 0 ) writer.value( "PERIODICITY" , n_plm / (double)n_clm );

  // IMI statistics
  writer.value( "N_IMI" , (int)all_imi.size() );
  if ( ! all_imi.empty() )
    {
      writer.value( "IMI_MEAN" , vmean( all_imi ) );
      if ( all_imi.size() >= 2 ) writer.value( "IMI_SD" , vsd( all_imi ) );
      writer.value( "IMI_MEDIAN" , vmed( all_imi ) );
      writer.value( "IMI_MIN" , vmin( all_imi ) );
      writer.value( "IMI_MAX" , vmax( all_imi ) );
    }
  writer.value( "N_IMI_PERIODIC" , (int)periodic_imi.size() );
  if ( ! log_periodic_imi.empty() )
    {
      writer.value( "LOG_IMI_MEAN" , vmean( log_periodic_imi ) );
      if ( log_periodic_imi.size() >= 2 ) writer.value( "LOG_IMI_SD" , vsd( log_periodic_imi ) );
    }
  writer.value( "N_IMI_SHORT" , n_short );
  writer.value( "N_IMI_LONG" , n_long );

  // laterality (LR only)
  if ( lr_mode )
    {
      writer.value( "N_LM_L" , n_lm_side[0] );  writer.value( "N_LM_R" , n_lm_side[1] );  writer.value( "N_LM_B" , n_lm_side[2] );
      writer.value( "N_CLM_L" , n_clm_side[0] ); writer.value( "N_CLM_R" , n_clm_side[1] ); writer.value( "N_CLM_B" , n_clm_side[2] );
      writer.value( "N_PLM_L" , n_plm_side[0] ); writer.value( "N_PLM_R" , n_plm_side[1] ); writer.value( "N_PLM_B" , n_plm_side[2] );
      writer.value( "N_COMP_L" , n_comp_side[0] ); writer.value( "N_COMP_R" , n_comp_side[1] );
      const double sh = stratum_hours( stg , "S" );
      if ( stg.present && sh > 0 )
        {
          writer.value( "PLMI_L" , n_plm_side[0] / sh );
          writer.value( "PLMI_R" , n_plm_side[1] / sh );
          writer.value( "PLMI_B" , n_plm_side[2] / sh );
        }
    }

  //
  // stage/state-stratified output
  //

  if ( stg.present )
    {
      static const char * strata[] = { "W" , "S" , "NREM" , "REM" , "N1" , "N2" , "N3" };
      for ( int k = 0 ; k < 7 ; ++k )
        {
          const std::string ss = strata[k];
          const double hrs = stratum_hours( stg , ss );
          int c_lm = 0 , c_clm = 0 , c_plm = 0 , c_iso = 0 , c_short = 0;
          double plm_dur_sum = 0;
          for ( size_t i = 0 ; i < events.size() ; ++i )
            {
              const lm_event_t & e = events[i];
              if ( ! in_stratum( e , ss ) ) continue;
              ++c_lm;
              if ( e.clm ) ++c_clm;
              if ( e.plm ) { ++c_plm; plm_dur_sum += e.dur; }
              if ( e.clm && ! e.plm ) ++c_iso;
              if ( e.clm && ! is_na( e.imi_prev ) && e.imi_prev < P.imi_min ) ++c_short;
            }
          if ( c_lm == 0 && hrs == 0 ) continue;

          writer.level( ss , globals::stage_strat );
          writer.value( "DUR_HR" , hrs );
          writer.value( "N_LM" , c_lm );
          writer.value( "N_CLM" , c_clm );
          writer.value( "N_PLM" , c_plm );
          writer.value( "N_ISOLATED_CLM" , c_iso );
          writer.value( "N_SHORT_IMI_CLM" , c_short );
          if ( hrs > 0 )
            {
              writer.value( "LM_INDEX" , c_lm / hrs );
              writer.value( "CLM_INDEX" , c_clm / hrs );
              writer.value( "PLM_INDEX" , c_plm / hrs );
              writer.value( "ISOLATED_CLM_INDEX" , c_iso / hrs );
              writer.value( "SHORT_IMI_CLM_INDEX" , c_short / hrs );
            }
          if ( c_plm > 0 ) writer.value( "PLM_DUR_MEAN" , plm_dur_sum / c_plm );
        }
      writer.unlevel( globals::stage_strat );

      // convenience top-level indices
      const double h_s = stratum_hours( stg , "S" ) , h_w = stratum_hours( stg , "W" );
      const double h_nrem = stratum_hours( stg , "NREM" ) , h_rem = stratum_hours( stg , "REM" );
      int p_s=0,p_w=0,p_nrem=0,p_rem=0;
      for ( size_t i = 0 ; i < events.size() ; ++i )
        {
          const lm_event_t & e = events[i];
          if ( ! e.plm ) continue;
          if ( in_stratum(e,"S") ) ++p_s;
          if ( in_stratum(e,"W") ) ++p_w;
          if ( in_stratum(e,"NREM") ) ++p_nrem;
          if ( in_stratum(e,"REM") ) ++p_rem;
        }
      if ( h_s > 0 ) writer.value( "PLMI" , p_s / h_s );
      if ( h_w > 0 ) writer.value( "PLMWI" , p_w / h_w );
      if ( h_nrem > 0 ) writer.value( "PLMI_NREM" , p_nrem / h_nrem );
      if ( h_rem > 0 ) writer.value( "PLMI_REM" , p_rem / h_rem );
    }

  //
  // verbose event-level output
  //

  if ( P.verbose )
    {
      // TYPE=COMP
      writer.level( "COMP" , "TYPE" );
      for ( size_t i = 0 ; i < comps.size() ; ++i )
        {
          const lm_component_t & c = comps[i];
          writer.level( mkid( 'C' , c.comp_id ) , "EVENT" );
          writer.value( "SIDE" , std::string( side_str( c.side ) ) );
          writer.value( "SIG" , c.sig );
          writer.value( "DETECTOR" , std::string( detector_str( P ) ) );
          writer.value( "GRAMMAR" , std::string( "WASM2016" ) );
          writer.value( "START" , c.onset_tp / (double)globals::tp_1sec );
          writer.value( "STOP" , c.offset_tp / (double)globals::tp_1sec );
          writer.value( "DUR" , c.dur );
          writer.value( "BASE" , c.base_onset );
          writer.value( "SCALE" , c.scale_onset );
          writer.value( "ON_THR" , c.on_thr );
          writer.value( "OFF_THR" , c.off_thr );
          writer.value( "PEAK" , c.peak );
          writer.value( "PEAK_EXC" , c.peak_exc );
          writer.value( "PEAK_Z" , c.peak_z );
          writer.value( "MORPH" , c.morph_value );
          writer.value( "CLM" , (int)( c.clm_component ? 1 : 0 ) );
          writer.value( "STAGE" , c.stage );
          writer.value( "STATE" , c.state );
          writer.value( "QC_BASE_HIGH" , (int)( c.qc_base_high ? 1 : 0 ) );
        }
      writer.unlevel( "EVENT" );

      // TYPE=LM
      writer.level( "LM" , "TYPE" );
      for ( size_t i = 0 ; i < events.size() ; ++i )
        {
          writer.level( mkid( 'E' , events[i].evt_id ) , "EVENT" );
          write_event_vars( events[i] , P , mode );
        }
      writer.unlevel( "EVENT" );

      // TYPE=SEQ
      writer.level( "SEQ" , "TYPE" );
      for ( size_t i = 0 ; i < seqs.size() ; ++i )
        {
          const plm_sequence_t & s = seqs[i];
          writer.level( mkid( 'P' , s.seq_id ) , "EVENT" );
          writer.value( "N" , s.n );
          writer.value( "START" , s.onset_tp / (double)globals::tp_1sec );
          writer.value( "STOP" , s.offset_tp / (double)globals::tp_1sec );
          writer.value( "DUR" , s.dur );
          writer.value( "IMI_N" , (int)s.imis.size() );
          if ( ! s.imis.empty() )
            {
              writer.value( "IMI_MEAN" , vmean( s.imis ) );
              if ( s.imis.size() >= 2 ) writer.value( "IMI_SD" , vsd( s.imis ) );
              writer.value( "IMI_MIN" , vmin( s.imis ) );
              writer.value( "IMI_MAX" , vmax( s.imis ) );
            }
          writer.value( "N_L" , s.n_l );
          writer.value( "N_R" , s.n_r );
          writer.value( "N_B" , s.n_b );
          writer.value( "N_LEG" , s.n_leg );
        }
      writer.unlevel( "EVENT" );
      writer.unlevel( "TYPE" );
    }

  logger << "  LM complete: " << n_lm << " LM, " << n_clm << " CLM, "
         << n_plm << " PLM in " << seqs.size() << " sequence(s)\n\n";

}
