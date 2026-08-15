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

#include "dynamics/sigdyn.h"

#include "dynamics/qdynam.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "annot/annot.h"
#include "timeline/timeline.h"
#include "param.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "miscmath/miscmath.h"
#include "db/db.h"

#include <vector>
#include <string>
#include <set>
#include <map>
#include <utility>
#include <algorithm>
#include <cctype>
#include <cmath>

extern writer_t writer;
extern logger_t logger;

//
// SIGDYN : summarize the temporal dynamics of an arbitrary, already-existing
// continuous/epoch-level signal, in a way that is comparable across
// individuals.  This command touches neither dynamics/qdynam.{h,cpp},
// dynamics/evtdyn.cpp, nor timeline/hypno.{h,cpp} -- it only calls their
// existing, public interfaces:
//
//   0) simple stratified descriptive statistics (N/MEAN/SD/MIN/MAX/RANGE) of
//      the per-epoch signal mean, by sleep stage, stage x local stability,
//      NREM cycle, and stage x cycle -- new, local glue code only; per-epoch
//      stage/cycle lookups mirror the existing idiom used by
//      hypnogram_t::annotate() (timeline/hypno.cpp ~254-260) and by
//      qdynam_t::dynam_compile_cycles() (dynamics/qdynam.cpp ~42-107), both
//      of which call timeline_t::epoch_annotation( label , epoch ) with the
//      raw 'epoch' value as returned by timeline_t::next_epoch().
//
//   1) whole-recording trend/decile/NREM-cycle-stratified summary, via
//      qdynam_t, exactly the way spectral/welch.cpp already embeds it
//      (qd.init() once, qd.add() per epoch/channel, qd.proc_all() once at
//      the end).
//
//   2) annotation/hypnogram-anchored peri-event averaging, via a native
//      per-offset/bin engine implemented directly in this file (previously
//      this reused dsp/tlock.{h,cpp}'s tlock_t; that had no partial-window,
//      min-n, or bin-width support, and each addition meant bypassing more
//      of it -- so mode 2 is now self-contained here instead). Anchors come
//      from either an explicit annot= annotation class list, or
//      auto-discovered HYPNO-derived landmark/cycle/transition annotations
//      (see hypnogram_t::annotate(), timeline/hypno.cpp ~3241-3787, for the
//      exact naming convention matched below). The one piece still reused
//      is timeline_t::discontinuity() (timeline/timeline.h) -- a general,
//      independent timeline utility, not part of tlock_t -- to exclude any
//      instance whose window spans a genuine gap in the recording. Edge-of-
//      recording truncation is not treated as a discontinuity: each
//      offset/bin simply uses however many instances have valid data there,
//      so partial windows near the start/end of the recording fall out
//      naturally rather than needing special-case handling.
//

namespace {

  // ---- small local string helpers (no changes to helper/helper.h needed) ----

  bool ends_with( const std::string & s , const std::string & suffix )
  {
    if ( suffix.size() > s.size() ) return false;
    return s.compare( s.size() - suffix.size() , suffix.size() , suffix ) == 0;
  }

  // matches "...*_cycle_n<digits>" i.e. hypnogram_t::annotate()'s
  //   prefix + "cycle" + suffix + "n" + <cycle#>
  bool is_cycle_marker( const std::string & s )
  {
    const std::string tag = "_cycle_n";
    size_t p = s.find( tag );
    if ( p == std::string::npos ) return false;
    size_t q = p + tag.size();
    if ( q >= s.size() ) return false;
    for (size_t i=q; i<s.size(); i++)
      if ( ! std::isdigit( (unsigned char)s[i] ) ) return false;
    return true;
  }

  // does 'name' look like one of hypnogram_t::annotate()'s generated
  // annotation classes that make useful peri-event anchors?  i.e.
  //  - stage-transition instants:  <prefix>_tr_<S1>_<S2>
  //  - night landmarks:            <prefix>_t0_start .. <prefix>_t6_stop
  //  - per-cycle onset markers:    <prefix>_cycle_n<N>
  // (annotate() always uses '_' as the default separator between prefix
  //  and these fixed suffixes -- see timeline/hypno.cpp; matching by
  //  suffix therefore works regardless of the user-chosen HYPNO 'annot='
  //  prefix)
  bool is_hypno_anchor_class( const std::string & name )
  {
    static const std::vector<std::string> suffixes = {
      "_tr_NR_R" , "_tr_NR_W" , "_tr_R_NR" , "_tr_R_W" , "_tr_W_NR" , "_tr_W_R" ,
      "_t0_start" , "_t1_lights_out" , "_t2_sleep_onset" , "_t3_sleep_midpoint" ,
      "_t4_final_wake" , "_t5_lights_on" , "_t6_stop"
    };

    for (int i=0; i<(int)suffixes.size(); i++)
      if ( ends_with( name , suffixes[i] ) ) return true;

    return is_cycle_marker( name );
  }

  // nearest-sample index in 'tp' for time-point 't' (tp assumed sorted ascending)
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

  // furthest sample offset from 'anchor', in the direction of 'sign'
  // (+1/-1), reachable without crossing a genuine discontinuity in the
  // recording -- i.e. where timeline_t::discontinuity() first becomes
  // true somewhere between 'anchor' and that point. 'bound' is the
  // desired furthest sample in that direction (already clipped to the
  // recording's own extent). Reachability is monotonic (if sp is
  // reachable, every point between anchor and sp is too), so this is a
  // binary search: one discontinuity() call in the common gap-free case,
  // O(log window-length-in-samples) in the rare case a gap exists
  // somewhere in range. This is what lets one anchor instance contribute
  // to bins on the near side of a gap while correctly excluding bins
  // beyond it, rather than discarding the whole instance (matching the
  // same edge-of-recording philosophy in the file header comment, just
  // applied to internal gaps too).
  int furthest_reachable( const std::vector<uint64_t> & tp , int sr , int anchor , int bound , int sign )
  {
    if ( bound == anchor ) return anchor;

    const int a = sign > 0 ? anchor : bound;
    const int b = sign > 0 ? bound : anchor;
    if ( ! timeline_t::discontinuity( tp , sr , a , b ) ) return bound;

    int lo = 0;                          // offset from anchor, known reachable
    int hi = std::abs( bound - anchor );  // offset from anchor, known not fully reachable
    while ( hi - lo > 1 )
      {
	const int mid = ( lo + hi ) / 2;
	const int sp  = anchor + sign * mid;
	const int aa  = sign > 0 ? anchor : sp;
	const int bb  = sign > 0 ? sp : anchor;
	if ( timeline_t::discontinuity( tp , sr , aa , bb ) ) hi = mid;
	else lo = mid;
      }
    return anchor + sign * lo;
  }

  // start/middle/end alignment, shared by anchor= (where within an
  // annotation instance's interval the t=0 anchor sits) and bin-align=
  // (which part of each bin sits at its labeled offset)
  enum align_t { ALIGN_START , ALIGN_MIDDLE , ALIGN_END };

  align_t parse_align( const std::string & s , const std::string & pname )
  {
    if ( s == "start" ) return ALIGN_START;
    if ( s == "middle" ) return ALIGN_MIDDLE;
    if ( s == "end" ) return ALIGN_END;
    Helper::halt( pname + " must be one of start, middle, end" );
    return ALIGN_START;
  }

  // simple SD-threshold outlier control on a set of values collected at
  // one offset/bin: either drop values beyond 'th' SDs of the mean, or
  // clip ('winsorize') them to 'winsor' SDs -- th/winsor <= 0 disables
  // each respectively.  Deliberately simple (single pass, un-iterated)
  // rather than a general-purpose utility, matching the narrow scope of
  // TLOCK's own th=/win= this replaces.
  void apply_outlier_control( std::vector<double> * v , double th , double winsor )
  {
    if ( v->size() < 2 ) return;
    const double m = MiscMath::mean( *v );
    const double sd = MiscMath::sdev( *v , m );
    if ( sd <= 0 ) return;

    if ( th > 0 )
      {
	std::vector<double> kept;
	for (int i=0; i<(int)v->size(); i++)
	  if ( std::fabs( (*v)[i] - m ) <= th * sd )
	    kept.push_back( (*v)[i] );
	if ( kept.size() > 0 ) *v = kept;
      }
    else if ( winsor > 0 )
      {
	const double lo = m - winsor * sd;
	const double hi = m + winsor * sd;
	for (int i=0; i<(int)v->size(); i++)
	  {
	    if ( (*v)[i] < lo ) (*v)[i] = lo;
	    else if ( (*v)[i] > hi ) (*v)[i] = hi;
	  }
      }
  }

}


void dsptools::sigdyn( edf_t & edf , param_t & param )
{

  signal_list_t signals = edf.header.signal_list( param.requires( "sig" ) , true );

  if ( signals.size() < 1 )
    {
      logger << "  *** no signals selected for SIGDYN\n";
      return;
    }

  const int ns = signals.size();

  // Modes 0 and 1 need an epoched timeline; rather than requiring the
  // caller to have run EPOCH explicitly, epoch with Luna's usual defaults
  // if not already epoched (a no-op, returning the existing count, if
  // epochs already exist) -- same idiom used e.g. by spindles/spindles.cpp,
  // dsp/detrend.cpp, suds/place.cpp.
  edf.timeline.ensure_epoched();

  // SIGDYN's own options (modes 0 and 2, native to this command) --
  // printed separately from qdynam_t::init()'s own "dynam options:" dump
  // (mode 1, inherited wholesale from DYNAM/PSD's dynam= submodule) so the
  // two option sets and their different provenance are both visible.
  {
    logger << "  SIGDYN options:\n"
	   << "    epoch-stats  = " << ( ( param.has("epoch-stats") ? param.yesno("epoch-stats") : true ) ? "T" : "F" ) << "\n"
	   << "    stable-flank = " << ( param.has("stable-flank") ? param.value("stable-flank") : "1" ) << "\n"
	   << "    annot        = " << ( param.has("annot") ? param.value("annot") : "." ) << "\n"
	   << "    hypno-annot  = " << ( ( param.has("hypno-annot") ? param.yesno("hypno-annot") : true ) ? "T" : "F" ) << "\n"
	   << "    w            = " << ( param.has("w") ? param.value("w") : "60" ) << "\n"
	   << "    anchor       = " << ( param.has("anchor") ? param.value("anchor") : "start" ) << "\n"
	   << "    bin-align    = " << ( param.has("bin-align") ? param.value("bin-align") : "middle" ) << "\n"
	   << "    bin          = " << ( param.has("bin") ? param.value("bin") : "." ) << "\n"
	   << "    inc          = " << ( param.has("inc") ? param.value("inc") : ( param.has("bin") ? param.value("bin") : "." ) ) << "\n"
	   << "    require-full = " << ( ( param.has("require-full") ? param.yesno("require-full") : false ) ? "T" : "F" ) << "\n"
	   << "    min-n        = " << ( param.has("min-n") ? param.value("min-n") : "1" ) << "\n";
  }

  //
  // ---------------- Mode 0: simple stratified descriptive statistics ----------------
  //
  //  for each requested signal, tabulate the per-epoch signal mean (same
  //  "epoch interval -> slice_t -> MiscMath::mean()" reduction as Mode 1)
  //  as plain N/MEAN/SD/MIN/MAX/RANGE, overall, and stratified by sleep
  //  stage, by stage x local stability, by NREM cycle, and by stage x
  //  cycle.  Runs whenever the recording is epoched (same precondition as
  //  Mode 1); can be disabled with epoch-stats=F.
  //

  const bool do_epoch_stats = edf.timeline.epoched() &&
    ( param.has( "epoch-stats" ) ? param.yesno( "epoch-stats" ) : true );

  if ( do_epoch_stats )
    {

      const int stable_flank = param.has( "stable-flank" ) ? param.requires_int( "stable-flank" ) : 1;
      if ( stable_flank < 0 ) Helper::halt( "stable-flank must be >= 0" );

      // local factor name for the stability level (does not collide with
      // any globals::*_strat already used in this file or elsewhere: no
      // matches the "REGION" factor convention used by HDSTATS
      // (pops/hdstats.cpp) for its own STABLE/TRANS stratification
      const std::string stab_strat = "REGION";

      //
      // classify every current epoch once (independent of channel): sleep
      // stage (empty string if none of the recognized labels match) and
      // NREM cycle (empty string if none of _NREMC_1.._NREMC_8 match).
      // Mirrors, label-for-label and idiom-for-idiom, hypnogram_t's own
      // per-epoch stage lookup (timeline/hypno.cpp ~254-260) and
      // qdynam_t::dynam_compile_cycles()'s cycle lookup (dynamics/qdynam.cpp
      // ~42-107): both call epoch_annotation( label , epoch ) with the raw
      // 'epoch' value returned by next_epoch() (i.e. the same index used to
      // fetch the epoch's interval_t via timeline.epoch( epoch )).
      //

      std::vector<int> epoch_idx;         // raw epoch index, for timeline.epoch()/epoch_annotation()
      std::vector<std::string> stage;     // "" == no stage
      std::vector<int> cyc;               // 0 == no cycle; globals::cycle_strat is a numeric-only factor (db.cpp)

      bool any_stage = false;
      bool any_cycle = false;

      edf.timeline.first_epoch();

      while ( 1 )
	{

	  int epoch = edf.timeline.next_epoch();
	  if ( epoch == -1 ) break;

	  epoch_idx.push_back( epoch );

	  std::string stg = "";
	  if      ( edf.timeline.epoch_annotation( "W"     , epoch ) ) stg = "W";
	  else if ( edf.timeline.epoch_annotation( "N1"    , epoch ) ) stg = "N1";
	  else if ( edf.timeline.epoch_annotation( "N2"    , epoch ) ) stg = "N2";
	  else if ( edf.timeline.epoch_annotation( "N3"    , epoch ) ) stg = "N3";
	  else if ( edf.timeline.epoch_annotation( "NREM4" , epoch ) ) stg = "NREM4";
	  else if ( edf.timeline.epoch_annotation( "R"     , epoch ) ) stg = "R";

	  if ( stg != "" ) any_stage = true;
	  stage.push_back( stg );

	  int c = 0;
	  if      ( edf.timeline.epoch_annotation( "_NREMC_1" , epoch ) ) c = 1;
	  else if ( edf.timeline.epoch_annotation( "_NREMC_2" , epoch ) ) c = 2;
	  else if ( edf.timeline.epoch_annotation( "_NREMC_3" , epoch ) ) c = 3;
	  else if ( edf.timeline.epoch_annotation( "_NREMC_4" , epoch ) ) c = 4;
	  else if ( edf.timeline.epoch_annotation( "_NREMC_5" , epoch ) ) c = 5;
	  else if ( edf.timeline.epoch_annotation( "_NREMC_6" , epoch ) ) c = 6;
	  else if ( edf.timeline.epoch_annotation( "_NREMC_7" , epoch ) ) c = 7;
	  else if ( edf.timeline.epoch_annotation( "_NREMC_8" , epoch ) ) c = 8;

	  if ( c != 0 ) any_cycle = true;
	  cyc.push_back( c );

	}

      const int ne = (int)epoch_idx.size();

      // stability: epoch 'e' (with a determined stage) is STABLE iff every
      // neighbor within stable_flank epochs on BOTH sides exists and shares
      // the same stage; a missing neighbor (off the edge of the recording)
      // or a differently-staged neighbor makes it TRANS.  Undefined (false)
      // for unstaged epochs, which are simply excluded from the STAGE/REGION
      // tables below.

      std::vector<bool> stable( ne , false );

      for (int e=0; e<ne; e++)
	{
	  if ( stage[e] == "" ) continue;

	  bool ok = true;
	  for (int k=1; k<=stable_flank; k++)
	    {
	      const int lo = e - k;
	      const int hi = e + k;
	      if ( lo < 0 || stage[lo] != stage[e] ) { ok = false; break; }
	      if ( hi >= ne || stage[hi] != stage[e] ) { ok = false; break; }
	    }
	  stable[e] = ok;
	}

      //
      // small local running accumulator for N/MEAN/SD/MIN/MAX/RANGE
      // (stats/running-stats.h's running_stats_t has no MIN/MAX, so this
      // simple by-hand accumulator is used instead, per-stratum)
      //

      struct acc_t {
	long long n;
	double sum, sumsq, mn, mx;
	acc_t() : n(0), sum(0), sumsq(0), mn(0), mx(0) { }
	void add( double x )
	{
	  if ( n == 0 ) { mn = mx = x; }
	  else { if ( x < mn ) mn = x; if ( x > mx ) mx = x; }
	  ++n; sum += x; sumsq += x * x;
	}
	double mean() const { return n > 0 ? sum / (double)n : 0; }
	double sd() const
	{
	  if ( n < 2 ) return 0;
	  const double m = mean();
	  double var = ( sumsq - (double)n * m * m ) / (double)( n - 1 );
	  if ( var < 0 ) var = 0;
	  return std::sqrt( var );
	}
      };

      auto emit = [&]( const acc_t & a )
	{
	  writer.value( "N" , (int)a.n );
	  // n==0 (e.g. every epoch's slice came back empty): MEAN/MIN/MAX/
	  // RANGE would otherwise default to 0, indistinguishable from a
	  // genuine zero-valued signal -- leave them unwritten instead
	  if ( a.n == 0 ) return;
	  writer.value( "MEAN" , a.mean() );
	  writer.value( "SD" , a.sd() );
	  writer.value( "MIN" , a.mn );
	  writer.value( "MAX" , a.mx );
	  writer.value( "RANGE" , a.mx - a.mn );
	};

      for (int s=0; s<ns; s++)
	{

	  if ( edf.header.is_annotation_channel( signals(s) ) ) continue;

	  writer.level( signals.label(s) , globals::signal_strat );

	  acc_t overall;
	  std::map<std::string,acc_t> by_stage;
	  std::map<std::pair<std::string,bool>,acc_t> by_stage_stab;
	  std::map<int,acc_t> by_cycle;
	  std::map<std::pair<std::string,int>,acc_t> by_stage_cycle;

	  for (int e=0; e<ne; e++)
	    {

	      interval_t interval = edf.timeline.epoch( epoch_idx[e] );

	      slice_t slice( edf , signals(s) , interval );
	      const std::vector<double> * d = slice.pdata();

	      if ( d->size() == 0 ) continue;

	      const double value = MiscMath::mean( *d );

	      overall.add( value );

	      if ( stage[e] != "" )
		{
		  by_stage[ stage[e] ].add( value );
		  by_stage_stab[ std::make_pair( stage[e] , stable[e] ) ].add( value );
		}

	      if ( cyc[e] != 0 )
		by_cycle[ cyc[e] ].add( value );

	      if ( stage[e] != "" && cyc[e] != 0 )
		by_stage_cycle[ std::make_pair( stage[e] , cyc[e] ) ].add( value );

	    }

	  // 1) CH : overall, all epochs whether staged or not
	  emit( overall );

	  // 2) CH,STAGE  and  3) CH,STAGE,REGION
	  if ( any_stage )
	    {

	      std::map<std::string,acc_t>::const_iterator ii = by_stage.begin();
	      while ( ii != by_stage.end() )
		{
		  writer.level( ii->first , globals::stage_strat );
		  emit( ii->second );
		  ++ii;
		}
	      writer.unlevel( globals::stage_strat );

	      std::map<std::pair<std::string,bool>,acc_t>::const_iterator jj = by_stage_stab.begin();
	      while ( jj != by_stage_stab.end() )
		{
		  writer.level( jj->first.first , globals::stage_strat );
		  writer.level( jj->first.second ? "STABLE" : "TRANS" , stab_strat );
		  emit( jj->second );
		  ++jj;
		}
	      if ( by_stage_stab.size() > 0 ) writer.unlevel( stab_strat );
	      writer.unlevel( globals::stage_strat );

	    }

	  // 4) CH,CYCLE  and  5) CH,STAGE,CYCLE
	  if ( any_cycle )
	    {

	      std::map<int,acc_t>::const_iterator kk = by_cycle.begin();
	      while ( kk != by_cycle.end() )
		{
		  writer.level( kk->first , globals::cycle_strat );
		  emit( kk->second );
		  ++kk;
		}
	      writer.unlevel( globals::cycle_strat );

	      if ( any_stage )
		{
		  std::map<std::pair<std::string,int>,acc_t>::const_iterator ll = by_stage_cycle.begin();
		  while ( ll != by_stage_cycle.end() )
		    {
		      writer.level( ll->first.first , globals::stage_strat );
		      writer.level( ll->first.second , globals::cycle_strat );
		      emit( ll->second );
		      ++ll;
		    }
		  if ( by_stage_cycle.size() > 0 ) writer.unlevel( globals::cycle_strat );
		  writer.unlevel( globals::stage_strat );
		}

	    }

	  writer.unlevel( globals::signal_strat );

	} // next signal

    }
  else if ( ! edf.timeline.epoched() )
    {
      logger << "  *** data not epoched, skipping SIGDYN stratified descriptive statistics (mode 0)\n";
    }


  //
  // ---------------- Mode 1: whole-recording trend/decile/cycle summary ----------------
  //
  //  reuses qdynam_t exactly as spectral/welch.cpp does (see welch.cpp
  //  ~108-110, ~1210-1230, ~1679-1680): construct once, add() per
  //  epoch/channel under the active CH stratum, proc_all() once at the end.
  //

  if ( edf.timeline.epoched() )
    {

      qdynam_t qd;
      qd.init( edf , param );

      for (int s=0; s<ns; s++)
	{

	  if ( edf.header.is_annotation_channel( signals(s) ) ) continue;

	  writer.level( signals.label(s) , globals::signal_strat );

	  const std::map<std::string,std::string> faclvl = writer.faclvl_notime();

	  edf.timeline.first_epoch();

	  while ( 1 )
	    {

	      int epoch = edf.timeline.next_epoch();
	      if ( epoch == -1 ) break;

	      interval_t interval = edf.timeline.epoch( epoch );

	      slice_t slice( edf , signals(s) , interval );
	      const std::vector<double> * d = slice.pdata();

	      if ( d->size() == 0 ) continue;

	      const double value = MiscMath::mean( *d );

	      // qdynam_t::add() expects the 0-based display-epoch convention
	      // (matches dynam_compile_cycles()'s uepochs / welch.cpp usage)
	      const int qd_epoch = edf.timeline.display_epoch( epoch ) - 1;

	      qd.add( faclvl , "MEAN" , qd_epoch , value );

	    }

	  writer.unlevel( globals::signal_strat );

	}

      qd.proc_all();

    }
  else
    {
      logger << "  *** data not epoched, skipping SIGDYN whole-recording trend/decile analysis (mode 1)\n";
    }


  //
  // ---------------- Mode 2: anchor/peri-event windowed averaging ----------------
  //
  //  native implementation (see file header comment above) -- the only
  //  reused piece is timeline_t::discontinuity(), a general timeline
  //  utility unrelated to any particular command
  //

  std::vector<std::string> anchors;
  std::set<std::string> already;

  // param is 'annot=', matching MEANS/EVTDYN's convention for "which
  // annotation classes to use" (not 'anchor=', which elsewhere in Luna --
  // e.g. EVTDYN's own anchor=MID/START/STOP -- means "which point within
  // an event's interval to use for timing", a different concept)
  if ( param.has( "annot" ) )
    {
      std::vector<std::string> a = param.strvector( "annot" );
      for (int i=0;i<(int)a.size(); i++)
	{
	  if ( already.find( a[i] ) == already.end() )
	    {
	      anchors.push_back( a[i] );
	      already.insert( a[i] );
	    }
	}
    }

  const bool use_hypno_anchors = param.has( "hypno-annot" ) ? param.yesno( "hypno-annot" ) : true;

  if ( use_hypno_anchors )
    {
      std::vector<std::string> names = edf.annotations->names();

      for (int i=0;i<(int)names.size(); i++)
	{
	  if ( is_hypno_anchor_class( names[i] ) && already.find( names[i] ) == already.end() )
	    {
	      anchors.push_back( names[i] );
	      already.insert( names[i] );
	    }
	}
    }

  if ( anchors.size() == 0 )
    {
      logger << "  *** no annot= annotations specified, and no HYPNO-derived anchors found: skipping SIGDYN peri-event analysis (mode 2)\n";
      return;
    }

  logger << "  SIGDYN peri-event anchors: " << Helper::stringize( anchors ) << "\n";

  const double half_window = param.has( "w" ) ? param.requires_dbl( "w" ) : 60;
  if ( half_window <= 0 ) Helper::halt( "w must be a positive number" );

  // anchor= : where within an annotation instance's interval the t=0
  // anchor sits (start/middle/end of the instance; a 0-duration/point
  // instance gives the same anchor regardless). bin-align= : which part of
  // each SEC bin (start/middle/end) sits at that bin's reported offset.
  // SEC offsets are always exactly k*inc (..,-2*inc,-inc,0,inc,2*inc,..),
  // symmetric about the anchor, whatever bin-align is chosen -- bin-align
  // only changes which part of the bin's own sample range that labeled
  // point corresponds to.
  const align_t anchor_point = parse_align( param.has( "anchor" ) ? param.value( "anchor" ) : "start" , "anchor" );
  const align_t bin_align = parse_align( param.has( "bin-align" ) ? param.value( "bin-align" ) : "middle" , "bin-align" );

  const bool take_log = param.has( "tolog" );

  const double outlier_th = param.has( "th" ) ? param.requires_dbl( "th" ) : -1;

  const double outlier_winsor = param.has( "win" ) ? param.requires_dbl( "win" ) : -1;

  // minimum number of contributing instances required before emitting any
  // output for a given offset/bin (and, at the coarser level, for the
  // per-anchor summary row) -- matches the spirit of qdynam_t's own
  // 'dynam-min-ne=' (dynamics/qdynam.cpp); default 1 means no filtering
  const int min_n = param.has( "min-n" ) ? param.requires_int( "min-n" ) : 1;
  if ( min_n < 1 ) Helper::halt( "min-n must be a positive integer" );

  // optional coarser bin width (seconds) for the CH,ANNOT,SEC table;
  // default is native per-sample resolution (bin = 1 sample), same as
  // before
  const double bin_sec = param.has( "bin" ) ? param.requires_dbl( "bin" ) : 0;
  if ( bin_sec < 0 ) Helper::halt( "bin must be a positive number of seconds" );

  // optional step between bins (seconds); < bin gives overlapping/sliding
  // bins, > bin gives gaps, == bin (the default, if bin is set) gives the
  // original non-overlapping/tiled behavior. Only meaningful alongside bin=.
  if ( param.has( "inc" ) && ! param.has( "bin" ) )
    Helper::halt( "inc= requires bin= to also be set" );
  const double inc_sec = param.has( "inc" ) ? param.requires_dbl( "inc" ) : bin_sec;
  if ( param.has( "inc" ) && inc_sec <= 0 ) Helper::halt( "inc must be a positive number of seconds" );

  // if set, an instance contributes only when its *entire* window (bounds
  // and discontinuity) is available -- the original, stricter TLOCK-style
  // all-or-nothing-per-instance behavior; default is the more permissive
  // per-offset/bin partial-window inclusion (see file header comment)
  const bool require_full = param.has( "require-full" ) ? param.yesno( "require-full" ) : false;

  std::vector<double> Fs = edf.header.sampling_freq( signals );

  // tolog=T: log() is undefined for values <= 0; warn once (rather than
  // silently emitting -inf/NaN into every downstream mean/SD/median) and
  // skip only the offending samples
  bool warned_tolog_nonpositive = false;

  for (int s=0; s<ns; s++)
    {

      if ( edf.header.is_annotation_channel( signals(s) ) ) continue;

      writer.level( signals.label(s) , globals::signal_strat );

      // sample time-points for this channel (used to convert anchor
      // instants into nearest-sample indices)
      slice_t wtrace( edf , signals(s) , edf.timeline.wholetrace() );
      const std::vector<uint64_t> * tp = wtrace.ptimepoints();

      for (int a=0; a<(int)anchors.size(); a++)
	{

	  annot_t * annot = edf.annotations->find( anchors[a] );
	  if ( annot == NULL ) continue;

	  // build sample-index offsets for this channel/anchor combination,
	  // and total span (matching MEANS' "S": total annotation span in
	  // seconds, summed here over the same instances used as anchors)
	  std::vector<int> idx;
	  double total_span_sec = 0;
	  annot_map_t::const_iterator ii = annot->interval_events.begin();
	  while ( ii != annot->interval_events.end() )
	    {
	      const uint64_t i_start = ii->first.interval.start;
	      const uint64_t i_stop  = ii->first.interval.stop;

	      uint64_t t_anchor;
	      switch ( anchor_point )
		{
		case ALIGN_START:  t_anchor = i_start; break;
		case ALIGN_END:    t_anchor = i_stop; break;
		case ALIGN_MIDDLE: default: t_anchor = i_start + ( i_stop - i_start ) / 2; break;
		}

	      const int p = nearest_sample( *tp , t_anchor );
	      if ( p >= 0 )
		{
		  idx.push_back( p );
		  total_span_sec += ii->first.interval.duration_sec();
		}
	      ++ii;
	    }

	  if ( idx.size() == 0 ) continue;

	  const int sr_i = (int)Fs[s];
	  const int half_points = (int) std::lround( half_window * sr_i );
	  const int wsize = (int)wtrace.pdata()->size();
	  const int bin_samples = bin_sec > 0 ? std::max( 1 , (int) std::lround( bin_sec * sr_i ) ) : 1;
	  const int inc_samples = bin_sec > 0 ? std::max( 1 , (int) std::lround( inc_sec * sr_i ) ) : 1;
	  const std::vector<double> * wd = wtrace.pdata();

	  // symmetric bin-index grid: k = 0 is always the anchor itself, and
	  // every other bin is exactly k*inc_samples away from it -- so SEC
	  // offsets are always exactly 0, +/-inc, +/-2*inc, .. regardless of
	  // bin-align or bin width. The same grid is used for every instance,
	  // so bins line up across instances regardless of which ones can
	  // complete which.
	  const int k_max = (int) std::floor( half_points / (double)inc_samples );
	  const int global_k_lo = -k_max;
	  const int global_k_hi = k_max;

	  // which part of the bin (start/middle/end) sits at its labeled
	  // k*inc_samples offset -- e.g. bin-align=start means the bin *begins*
	  // at the labeled offset; end means it *ends* there; middle centers it
	  // (biased down by one sample for even widths). No effect in native
	  // (bin_samples==1) mode, where a bin is always a single sample.
	  int align_offset_samples = 0;
	  switch ( bin_align )
	    {
	    case ALIGN_START:  align_offset_samples = 0; break;
	    case ALIGN_END:    align_offset_samples = bin_samples - 1; break;
	    case ALIGN_MIDDLE: default: align_offset_samples = ( bin_samples - 1 ) / 2; break;
	    }

	  // bin index k -> raw values / contributing instances. With
	  // inc_samples == bin_samples (default whenever bin= is set) bins are
	  // non-overlapping/tiled; inc < bin overlaps, inc > bin gaps. Native
	  // (bin= unset) is bin_samples = inc_samples = 1.
	  std::map<int, std::vector<double> > bin_values;
	  std::map<int, std::set<int> >       bin_instances;
	  std::set<int> used; // instances contributing to at least one complete bin

	  for (int i=0; i<(int)idx.size(); i++)
	    {
	      const int lo = idx[i] - half_points;
	      const int hi = idx[i] + half_points;
	      const int lo_c = std::max( 0 , lo );
	      const int hi_c = std::min( wsize - 1 , hi );
	      if ( lo_c > hi_c ) continue;

	      // require-full=T: reject this instance outright unless its
	      // entire desired window fits within the recording AND is
	      // gap-free (original, stricter TLOCK-style all-or-nothing-per-
	      // instance behavior)
	      if ( require_full )
		{
		  if ( lo_c != lo || hi_c != hi ) continue;
		  if ( timeline_t::discontinuity( *tp , sr_i , lo_c , hi_c ) ) continue;
		}

	      // otherwise, find how far this instance reaches on each side
	      // without crossing a genuine gap (edge-of-recording truncation,
	      // already reflected in lo_c/hi_c, is not itself a discontinuity)
	      const int left_reach  = require_full ? lo_c : furthest_reachable( *tp , sr_i , idx[i] , lo_c , -1 );
	      const int right_reach = require_full ? hi_c : furthest_reachable( *tp , sr_i , idx[i] , hi_c , +1 );

	      // only COMPLETE bins contribute: this instance is included in a
	      // given bin only if every sample that bin covers is within its
	      // reachable range -- never a partial/thinned contribution. A
	      // bin too wide to fit in a short contiguous stretch, or that
	      // runs off the recording edge, simply isn't formed for this
	      // instance (it may still be formed by other instances, and this
	      // instance may still complete other, better-covered bins).
	      //
	      // each contributing instance adds exactly ONE value per bin (its
	      // own within-bin mean across however many samples the bin spans)
	      // rather than one push per raw sample -- so N (bin_instances),
	      // and the M/SD/MD/outlier-control computed from bin_values, are
	      // all consistently event-weighted, not sample-weighted (with
	      // bin=unset/native resolution, bin_samples==1 and this is a
	      // no-op: one sample IS one instance's contribution)
	      for (int k=global_k_lo; k<=global_k_hi; k++)
		{
		  const int bin_start = k * inc_samples - align_offset_samples;
		  const int bin_end   = bin_start + bin_samples - 1;

		  if ( idx[i] + bin_start < left_reach || idx[i] + bin_end > right_reach ) continue;

		  double sum = 0;
		  int cnt = 0;
		  for (int sp = idx[i] + bin_start; sp <= idx[i] + bin_end; sp++)
		    {
		      double v = (*wd)[sp];
		      if ( take_log )
			{
			  // log() undefined for v <= 0: skip this sample
			  // rather than contaminating the instance's bin mean
			  // with -inf/NaN
			  if ( v <= 0 )
			    {
			      if ( ! warned_tolog_nonpositive )
				{
				  logger << "  *** tolog=T but signal has values <= 0; skipping those samples\n";
				  warned_tolog_nonpositive = true;
				}
			      continue;
			    }
			  v = log( v );
			}
		      sum += v; ++cnt;
		    }

		  // an instance whose entire bin span was non-positive (and
		  // so contributed zero valid log() samples) does not
		  // contribute to this bin at all
		  if ( cnt == 0 ) continue;

		  used.insert( i );
		  bin_instances[ k ].insert( i );
		  bin_values[ k ].push_back( sum / cnt );
		}
	    }

	  if ( used.size() == 0 ) continue;

	  writer.level( anchors[a] , globals::annot_strat );

	  // per-offset/bin SEC table: only bins with >=min_n contributing
	  // instances are emitted (independent of whether the anchor-level
	  // summary row below meets the same threshold)
	  double sumAll = 0 , sumL = 0 , sumR = 0;
	  int nAll = 0 , nL = 0 , nR = 0;

	  std::map<int, std::vector<double> >::iterator bi = bin_values.begin();
	  while ( bi != bin_values.end() )
	    {
	      const int k = bi->first;
	      std::vector<double> & vals = bi->second;

	      apply_outlier_control( &vals , outlier_th , outlier_winsor );

	      const int bin_n_instances = (int) bin_instances[ k ].size();

	      if ( bin_n_instances >= min_n && vals.size() > 0 )
		{
		  const double bin_m = MiscMath::mean( vals );

		  // label is always exactly k*inc -- i.e. 0, +/-inc, +/-2*inc,
		  // .. -- symmetric about the anchor regardless of bin-align
		  // or bin width; bin-align only changes which part of the
		  // bin's own sample range corresponds to this offset.
		  writer.level( ( k * inc_samples ) / (double)sr_i , globals::sec_strat );
		  writer.value( "N" , bin_n_instances );
		  writer.value( "M" , bin_m );
		  if ( vals.size() > 1 ) writer.value( "SD" , MiscMath::sdev( vals , bin_m ) );
		  writer.value( "MD" , MiscMath::median( vals ) );

		  sumAll += bin_m; ++nAll;
		  if ( k < 0 ) { sumL += bin_m; ++nL; }
		  else if ( k > 0 ) { sumR += bin_m; ++nR; }
		}

	      ++bi;
	    }
	  writer.unlevel( globals::sec_strat );

	  // per-anchor summary (superset of MEANS' M/S/L/R): only emitted if
	  // the anchor as a whole meets min_n on distinct contributing
	  // instances (the SEC rows above already applied min_n per-bin)
	  if ( (int)used.size() >= min_n )
	    {
	      writer.value( "N" , (int)used.size() );
	      if ( nAll > 0 ) writer.value( "M" , sumAll / nAll );
	      if ( nL > 0 ) writer.value( "L" , sumL / nL );
	      if ( nR > 0 ) writer.value( "R" , sumR / nR );
	      writer.value( "S" , total_span_sec );
	    }

	  writer.unlevel( globals::annot_strat );

	} // next anchor

      writer.unlevel( globals::signal_strat );

    } // next signal

}
