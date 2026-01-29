
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

#include "medflt.h"
#include "param.h"

#include "stats/eigen_ops.h"

#include "edf/edf.h"
#include "edf/slice.h"

#include <vector>
#include <queue>
#include <stdexcept>
#include <cstddef>
#include <algorithm>

// Fast sliding median filter (two-heaps + lazy deletion + side tracking)
//   O(N log W), space O(W) plus O(N) bookkeeping
//   - constant-size window at edges via padding (Replicate or Reflect)
//   - even-W median = average of the two middle values
//   - if window > N: either throw (Strict) or shrink to N (Shrink)

enum class EdgePadding { Replicate, Reflect };

enum class WindowTooLarge { Strict, Shrink };

namespace median_filter_detail {

inline int clamp_int(int v, int lo, int hi) {
  return v < lo ? lo : (v > hi ? hi : v);
}

inline int reflect_index(int j, int N) {
  if (N <= 1) return 0;
  const int period = 2 * (N - 1);
  int m = j % period;
  if (m < 0) m += period;
  return (m <= (N - 1)) ? m : (period - m);
}

struct Node {
  double v;
  int id;
};

struct MaxCmp {
  bool operator()(const Node& a, const Node& b) const { return a.v < b.v; } // max-heap by v
};
struct MinCmp {
  bool operator()(const Node& a, const Node& b) const { return a.v > b.v; } // min-heap by v
};

struct DualHeap {
  // lo: max-heap (lower half), hi: min-heap (upper half)
  std::priority_queue<Node, std::vector<Node>, MaxCmp> lo;
  std::priority_queue<Node, std::vector<Node>, MinCmp> hi;

  std::vector<unsigned char> deleted; // deleted[id] = 1 if id is out of window
  std::vector<unsigned char> side;    // side[id] = 0 => lo, 1 => hi  (valid for non-deleted ids)
  int loSize = 0; // active counts (non-deleted)
  int hiSize = 0;

  explicit DualHeap(std::size_t max_ids)
    : deleted(max_ids, 0), side(max_ids, 0) {}

  bool is_deleted(int id) const { return deleted[(std::size_t)id] != 0; }

  void prune_lo() {
    while (!lo.empty() && is_deleted(lo.top().id)) lo.pop();
  }
  void prune_hi() {
    while (!hi.empty() && is_deleted(hi.top().id)) hi.pop();
  }

  // desiredLo = (W+1)/2 for odd W, and W/2 for even W (lo holds "lower middle")
  void rebalance(int desiredLo) {
    prune_lo();
    prune_hi();

    // Move from lo -> hi if lo too big
    while (loSize > desiredLo) {
      prune_lo();
      if (lo.empty()) break;
      Node n = lo.top(); lo.pop();
      if (is_deleted(n.id)) continue;

      hi.push(n);
      side[(std::size_t)n.id] = 1;
      --loSize;
      ++hiSize;
    }

    // Move from hi -> lo if lo too small
    while (loSize < desiredLo) {
      prune_hi();
      if (hi.empty()) break;
      Node n = hi.top(); hi.pop();
      if (is_deleted(n.id)) continue;

      lo.push(n);
      side[(std::size_t)n.id] = 0;
      ++loSize;
      --hiSize;
    }

    prune_lo();
    prune_hi();
  }

  void add(double v, int id, int desiredLo) {
    prune_lo();
    if (lo.empty() || v <= lo.top().v) {
      lo.push(Node{v, id});
      side[(std::size_t)id] = 0;
      ++loSize;
    } else {
      hi.push(Node{v, id});
      side[(std::size_t)id] = 1;
      ++hiSize;
    }
    rebalance(desiredLo);
  }

  void erase(int id, int desiredLo) {
    // Mark deleted and decrement the correct active side count deterministically.
    deleted[(std::size_t)id] = 1;

    if (side[(std::size_t)id] == 0) --loSize;
    else --hiSize;

    // Prune tops if needed and rebalance.
    prune_lo();
    prune_hi();
    rebalance(desiredLo);
  }

  double median_even_avg(int window) {
    prune_lo();
    prune_hi();
    if (window <= 0) Helper::halt( "median_filter: empty window" );

    if ((window & 1) != 0) {
      if (lo.empty()) Helper::halt( "median_filter: invalid state (odd)" );
      return lo.top().v;
    } else {
      if (lo.empty() || hi.empty()) Helper::halt( "median_filter: invalid state (even)");
      return 0.5 * (lo.top().v + hi.top().v);
    }
  }
};

} // namespace median_filter_detail

inline std::vector<double> median_filter_fast(
    const std::vector<double>& x,
    int window,
    double drop_top_frac = 0.0 ,
    EdgePadding padding = EdgePadding::Replicate,
    WindowTooLarge too_large = WindowTooLarge::Shrink )
{
  using namespace median_filter_detail;

  const int N = (int)x.size();
  std::vector<double> y((std::size_t)N);
  if (N == 0) return y;

  if (window <= 0) Helper::halt( "median_filter_fast: window must be >= 1");

  if (window > N) {
    if (too_large == WindowTooLarge::Strict) {
      Helper::halt( "median_filter_fast: window larger than series");
    } else {
      window = N;
    }
  }

  // allow dropping top X% of values, i.e. for baseline not impacted by extremes
  // if "extremes"/events are actually somewhat common - i.e. a trimmed median
  if (drop_top_frac < 0.0 || drop_top_frac >= 1.0)
    Helper::halt( "median_filter_fast: drop_top_frac must be in [0,1)" );
  
  const int r = (int)(drop_top_frac * window); // drop top r values (floor)
  const int M = window - r;                    // effective sample count used for median
  if (M <= 0)
    Helper::halt( "median_filter_fast: drop_top_frac too large for window" );

  
  const bool odd = (window & 1) != 0;
  const int v1 = odd ? (window - 1) / 2 : (window / 2);
  const int v2 = odd ? (window - 1) / 2 : (window / 2 - 1);

  auto map_index = [&](int j) -> int {
    if (padding == EdgePadding::Replicate) {
      return clamp_int(j, 0, N - 1);
    } else {
      return reflect_index(j, N);
    }
  };
  
  
  // adjusted to handle case of trimming
  const int desiredLo = (M % 2) ? ((M + 1) / 2) : (M / 2);
  
  // Total inserted IDs = initial window + (N-1) slides <= N + window
  const std::size_t max_ids = (std::size_t)(N + window + 5);
  DualHeap dh(max_ids);

  int next_id = 0;

  // Initialize window for i=0: indices [-v1 .. +v2]
  for (int k = -v1; k <= v2; ++k) {
    const int src = map_index(k);
    dh.add(x[(std::size_t)src], next_id++, desiredLo);
  }

  y[0] = dh.median_even_avg(M);

  // Slide: remove oldest ID, add one new ID
  for (int i = 1; i < N; ++i) {
    const int out_id = next_id - window; // oldest element currently in window

    const int in_src = map_index(i + v2);
    const double in_v = x[(std::size_t)in_src];
    const int in_id = next_id++;

    dh.erase(out_id, desiredLo);
    dh.add(in_v, in_id, desiredLo);

    y[(std::size_t)i] = dh.median_even_avg(M);
  }

  return y;
}

  

void dsptools::median_filter( edf_t & edf , param_t & param )
{
    
  // signal(s)
  const bool NO_ANNOTS = true; 
  signal_list_t signals = edf.header.signal_list( param.value( "sig" ) , NO_ANNOTS );  
  const int ns = signals.size();
  if ( ns == 0 ) return;

  // always perform this within each contig - so this forces
  // epoching here
  const int ne = edf.timeline.calc_epochs_contig();
  logger << "  iterating over " << ne << " contig-basad epochs\n";

  // window size
  const double hwin_sec = param.requires_dbl( "hwin" );

  // remove median, or return? (default = remove) 
  const bool remove = param.yesno( "remove" , true, true );

  // drop top fraction of values? (e.g. use in snore detector)
  const double trim_frac = param.has( "trim" ) ? param.requires_dbl( "trim" ) : 0 ; 
  
  logger << "  processed:";
  
  // process each signal
  for (int s=0; s<ns; s++)
    {
      
      // get sample rate
      int sr = edf.header.sampling_freq( signals(s) );
      
      int hwin = hwin_sec * sr;
     
      if ( hwin == 0 )
	{
	  logger << "  skipping " << signals.label(s) << ", sample rate too low\n";
	  continue;
	}

      // make a full window
      hwin = 1 + 2 * hwin;
      
      // get whole signal  (although we update epoch-by-epoch)
      // we need to this edit and return back 
      
      slice_t slice0( edf , signals(s) , edf.timeline.wholetrace() );
      std::vector<double> orig = *slice0.pdata();
      const int n = orig.size();
      
      // idx into orig
      int idx = 0;
      
      while ( 1 )
	{
	  
	  // next epoch
	  int epoch = edf.timeline.next_epoch();
	  
	  // all done?
	  if ( epoch == -1 ) break;
	  
	  // get data
	  interval_t interval = edf.timeline.epoch( epoch ) ;	  
	  slice_t slice( edf , signals(s) , interval );	  
	  const std::vector<double> * data = slice.pdata();

	  // filter
	  std::vector<double> flt = median_filter_fast( *data , hwin , trim_frac );

	  // update
	  for (int i=0; i<flt.size(); i++)
	    {
	      if ( idx == n ) Helper::halt( "internal error in median_filter()" );
	      orig[idx] = remove ? orig[idx] - flt[i] : flt[i];
	      ++idx;
	    }
	}
      
      logger << " " << signals.label(s);

      // send back to the EDF
      edf.update_signal( signals(s) , &orig );      
      
    }
  logger << "\n";
  
}


