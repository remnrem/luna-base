
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

#include "edf/insert.h"
#include "param.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "db/db.h"
#include "helper/logger.h"
#include "helper/helper.h"
#include "stats/eigen_ops.h"
#include "dsp/tsync.h"
#include "dsp/xcorr.h"
#include "dsp/spline.h"
#include "dsp/fir.h"

#include <limits>

extern writer_t writer;
extern logger_t logger;

namespace {

int two_pass_trim_outliers( const std::vector<double> & x ,
			    const std::vector<double> & y ,
			    std::vector<bool> * keep )
{
  const int n = x.size();
  if ( y.size() != n || n < 3 ) return 0;

  std::vector<bool> local_keep( n , true );
  std::vector<bool> * k = keep ? keep : &local_keep;
  if ( keep )
    {
      if ( (int)keep->size() != n ) keep->assign( n , true );
      else
	for (int i=0; i<n; i++) (*keep)[i] = true;
    }

  int removed = 0;
  auto fit_ols = [&]( const std::vector<bool> & mask , double * slope , double * intercept ) -> int
    {
      double sx = 0 , sy = 0;
      int nn = 0;
      for (int i=0; i<n; i++) if ( mask[i] ) { sx += x[i]; sy += y[i]; ++nn; }
      if ( nn < 2 ) { *slope = 0; *intercept = 0; return nn; }
      sx /= nn; sy /= nn;
      double sxy = 0 , sxx = 0;
      for (int i=0; i<n; i++) if ( mask[i] )
	{
	  sxy += ( x[i] - sx ) * ( y[i] - sy );
	  sxx += ( x[i] - sx ) * ( x[i] - sx );
	}
      *slope = sxx > 0 ? sxy / sxx : 0;
      *intercept = sy - *slope * sx;
      return nn;
    };

  for (int pass=0; pass<2; pass++)
    {
      double slope = 0 , intercept = 0;
      const int nn = fit_ols( *k , &slope , &intercept );
      if ( nn < 3 ) break;

      double ss = 0;
      for (int i=0; i<n; i++) if ( (*k)[i] )
	ss += std::pow( y[i] - ( intercept + slope * x[i] ) , 2 );
      const double resid_sd = std::sqrt( ss / ( nn - 2 ) );
      if ( resid_sd <= 0 ) break;

      bool any_removed = false;
      for (int i=0; i<n; i++) if ( (*k)[i] )
	{
	  const double resid = y[i] - ( intercept + slope * x[i] );
	  if ( std::fabs( resid ) > 3.0 * resid_sd )
	    {
	      (*k)[i] = false;
	      ++removed;
	      any_removed = true;
	    }
	}
      if ( ! any_removed ) break;
    }

  return removed;
}

}

edf_inserter_t::edf_inserter_t( edf_t & edf , param_t & param )
{
  
  // base EDF, already attached
  // second EDF, with signals to be inserted

  // both EDFs must be continuous (or effectively contiuous, i.e. no gaps)
  
  // two modes
  //   1) given pairs of signals, find the lag for each based on xcorr() or euclidean distance
  //        edf= pairs= w= verbose

  //   2) given a lag value (in seconds, or tp-units) insert the second channel, but with an offset 
  //        edf= sig (from secondary EDF)   align=<tp-units>

  annotation_set_t annotations2;
  edf_t edf2( &annotations2) ;
  
  if ( ! edf2.attach( param.requires( "edf" ) , "." ) )
    Helper::halt( "problem attaching second EDF, edf=" + param.value( "edf" ) );
  
  // currently, both EDFs must be continuous - or at least no gaps
  if ( ! ( edf.header.continuous &&  edf2.header.continuous ) )
    if ( edf.is_actually_discontinuous() || edf2.is_actually_discontinuous() )
      Helper::halt( "neither EDF can be discontinous, with gaps" );

  // compute offset from EDF header timestamps (edf2_start - edf1_start)
  // negative means edf2 started after edf1
  double header_offset = 0;
  bool header_offset_valid = false;
  {
    clocktime_t t1( edf.header.startdate , edf.header.starttime );
    clocktime_t t2( edf2.header.startdate , edf2.header.starttime );
    if ( t1.valid && t2.valid )
      {
	// offset convention: negative means edf2 started after edf1 (matches ystart - minidx)
	header_offset = -clocktime_t::difference_seconds( t1 , t2 );
	header_offset_valid = true;
	logger << "  header-derived offset: " << header_offset << " seconds (negative = edf2 starts after edf1)\n";
      }
    else
      logger << "  warning: could not parse EDF start times, no header-derived offset available\n";
  }

  // insert mode?
  if ( param.has( "offset" ) )
    {
      const std::string signal_label = param.value( "sig" );
      const double offset = param.requires_dbl( "offset" );
      const std::string annot_label = param.has( "annot" ) ? param.value( "annot" ) : "" ;

      // optionally, timestretch secondary signal?   drift sec per 'secs' secs
      const double stretch_denom = param.has( "secs" ) ? param.requires_dbl( "secs" ) : 1 ;
      const double stretch_shift = param.has( "drift" ) ? param.requires_dbl( "drift" ) : 0 ;
      
      // e.g. if offset shifts (constant) by -10 seconds in 8 hours
      //  secs=28800 drift=-10
      // where 28800 = 8 * 60 * 60 
      
      // if just do shift  implies per 1 second... but to avoid floating point issues, probably
      // better to give a reasonable denominator (with secs) 
      const bool timestretch = param.has( "drift" ) && stretch_denom > 0 ;
      const double fac = timestretch ? stretch_shift / stretch_denom : 0 ; 
      
      insert( edf , edf2 , signal_label , offset , timestretch ? &fac : NULL , annot_label );
      // all done
      return;
    }


  // after estimation, optionally splice secondary signals into primary using fitted offset+drift
  const bool do_insert = param.has( "insert" );
  const std::string insert_sig   = param.has( "sig" )   ? param.value( "sig" )   : "*";
  const std::string insert_annot = param.has( "annot" ) ? param.value( "annot" ) : "";

  // method
  const bool euclidean = param.has( "euclidean" );

  // verbose outputs
  const bool verbose = param.has( "verbose" );

  // get signal pairs
  std::vector<std::string> pairs = param.strvector( "pairs" );

  if ( pairs.size() == 0 || pairs.size() % 2 )
    Helper::halt( "expecting an even number of channels pairs=sig1,sig2,sig1,sig2,... (or run in insert mode with offset= arg)" );

  const int np = pairs.size() / 2;

  std::vector<int> slot1, slot2, srs;
  std::vector<std::string> slab1, slab2;

  for (int s=0; s<np*2; s+=2)
    {
      const int s1 = edf.header.signal( pairs[s] );
      const int s2 = edf2.header.signal( pairs[s+1] );
      if ( s1 == -1 ) Helper::halt( "could not find " + pairs[s] + " in primary EDF" );
      if ( s2 == -1 ) Helper::halt( "could not find " + pairs[s+1] + " in secondary EDF" );
      slot1.push_back( s1 );
      slot2.push_back( s2 );
      slab1.push_back( pairs[s] );
      slab2.push_back( pairs[s+1] );
      int sr1 = edf.header.sampling_freq( s1 );
      int sr2 = edf2.header.sampling_freq( s2 );
      if ( sr1 != sr2 )
	Helper::halt( "sample rates must match for " + pairs[s] + " and " + pairs[s+1] );
      srs.push_back( sr1 );
    }

  int sr = srs[0];
  for (int p=1; p<np; p++)
    if ( srs[p] != sr )
      Helper::halt( "sample rates must match for all signal pairs" );

  // windowing parameters (ystart/yend defaults computed after data load to use ny)
  const bool user_start = param.has( "start" );
  const bool user_end   = param.has( "end" );
  double ylen_sec = param.has( "len" ) ? param.requires_dbl( "len" ) : 300;
  double yinc_sec = param.has( "inc" ) ? param.requires_dbl( "inc" ) : 60;
  if ( yinc_sec <= 0 ) Helper::halt( "inc must be greater than 0" );
  const int ysteps = param.has( "steps" ) ? param.requires_int( "steps" ) : INT_MAX;
  const double min_peak = param.has( "min-peak" ) ? param.requires_dbl( "min-peak" ) : 0.3;
  const bool do_warn = ! param.has( "no-warn" );
  const double warn_r2 = param.has( "warn-r2" ) ? param.requires_dbl( "warn-r2" ) : 0.5;
  const double warn_p_ok = param.has( "warn-p-ok" ) ? param.requires_dbl( "warn-p-ok" ) : 0.5;
  const double warn_peak = param.has( "warn-peak" ) ? param.requires_dbl( "warn-peak" ) : std::max( 0.35 , min_peak );
  const bool auto_try_mode = param.has( "auto-try" );
  const bool manual_try_windows = param.has( "try-start" ) || param.has( "try-len" ) || param.has( "try-inc" );
  const bool try_windows = auto_try_mode || manual_try_windows;
  if ( auto_try_mode && manual_try_windows )
    Helper::halt( "auto-try cannot be combined with try-start/try-len/try-inc" );

  // filter params: bandpass applied to segments before comparison (on by default)
  const bool do_filter  = ! param.has( "no-filter" );
  const double filt_low  = param.has( "filt-low" )  ? param.requires_dbl( "filt-low" )  : 0.5;
  const double filt_high = param.has( "filt-high" ) ? param.requires_dbl( "filt-high" ) : 15.0;
  if ( do_filter && filt_low >= filt_high )
    Helper::halt( "filt-low must be less than filt-high" );
  if ( do_filter && filt_high >= sr / 2.0 )
    Helper::halt( "filt-high must be less than Nyquist (" + Helper::dbl2str( sr/2.0 ) + " Hz)" );

  // offset search space:
  //   1) full-search (no offset constraint)
  //   2) explicit absolute interval via offset-range=min,max
  //   3) header-derived center with user-specified half-width via offset-margin=sec
  //   4) fallback header-derived center with default half-width
  const bool full_search = param.has( "full-search" );
  if ( full_search && ( param.has( "offset-range" ) || param.has( "offset-margin" ) ) )
    Helper::halt( "full-search is mutually exclusive with offset-range and offset-margin" );
  if ( param.has( "offset-range" ) && param.has( "offset-margin" ) )
    Helper::halt( "offset-range and offset-margin are mutually exclusive" );

  const bool constrained_offset = !full_search &&
    ( param.has( "offset-range" ) || param.has( "offset-margin" ) || header_offset_valid );
  std::vector<double> offsets;
  if ( full_search )
    {
      logger << "  using full-search over valid primary bounds (no offset constraint)\n";
    }
  else if ( param.has( "offset-range" ) )
    {
      offsets = param.dblvector( "offset-range" );
      if ( offsets.size() != 2 || offsets[1] <= offsets[0] )
	Helper::halt( "expecting offset-range=min,max" );
    }
  else if ( param.has( "offset-margin" ) )
    {
      if ( ! header_offset_valid )
	Helper::halt( "offset-margin requires valid EDF header start times in both files" );
      const double margin = param.requires_dbl( "offset-margin" );
      if ( margin <= 0 )
	Helper::halt( "offset-margin must be > 0" );
      offsets.push_back( header_offset - margin );
      offsets.push_back( header_offset + margin );
      logger << "  using header-derived offset with user margin: "
	     << offsets[0] << " to " << offsets[1] << " seconds\n";
    }
  else if ( header_offset_valid )
    {
      const double margin = 60;
      offsets.push_back( header_offset - margin );
      offsets.push_back( header_offset + margin );
      logger << "  using header-derived offset-range: " << offsets[0] << " to " << offsets[1] << " seconds\n";
    }

  // load all signal data into Eigen matrices
  std::vector<Eigen::VectorXd> dX(np), dY(np);
  for (int p=0; p<np; p++)
    {
      slice_t slice1( edf,  slot1[p], edf.timeline.wholetrace() );
      slice_t slice2( edf2, slot2[p], edf2.timeline.wholetrace() );
      const std::vector<double> * dx = slice1.pdata();
      const std::vector<double> * dy = slice2.pdata();

      // apply bandpass filter in-memory before storing
      if ( do_filter )
	{
	  std::vector<double> ripple_v( 1, 0.02 ), tw_v( 1, 1.0 );
	  std::vector<double> fx = dsptools::apply_fir( *dx , sr , fir_t::BAND_PASS ,
						        1 , ripple_v , tw_v , filt_low , filt_high );
	  std::vector<double> fy = dsptools::apply_fir( *dy , sr , fir_t::BAND_PASS ,
						        1 , ripple_v , tw_v , filt_low , filt_high );
	  dX[p] = Eigen::VectorXd::Zero( fx.size() );
	  dY[p] = Eigen::VectorXd::Zero( fy.size() );
	  for (int i=0; i<(int)fx.size(); i++) dX[p][i] = fx[i];
	  for (int i=0; i<(int)fy.size(); i++) dY[p][i] = fy[i];
	}
      else
	{
	  dX[p] = Eigen::VectorXd::Zero( dx->size() );
	  dY[p] = Eigen::VectorXd::Zero( dy->size() );
	  for (int i=0; i<(int)dx->size(); i++) dX[p][i] = (*dx)[i];
	  for (int i=0; i<(int)dy->size(); i++) dY[p][i] = (*dy)[i];
	}
    }

  const int nx = dX[0].size();
  const int ny = dY[0].size();

  auto zscore = []( Eigen::VectorXd v ) -> Eigen::VectorXd {
    const double m = v.mean();
    const double s = std::sqrt( (v.array() - m).square().mean() );
    Eigen::VectorXd r = v.array() - m;
    if ( s > 0 ) r /= s;
    return r;
  };

  std::vector<std::vector<double>> dX_prefix_sum;
  std::vector<std::vector<double>> dX_prefix_sq;
  if ( full_search && ! euclidean )
    {
      dX_prefix_sum.resize( np );
      dX_prefix_sq.resize( np );
      for (int p=0; p<np; p++)
	{
	  dX_prefix_sum[p].resize( nx + 1 , 0.0 );
	  dX_prefix_sq[p].resize( nx + 1 , 0.0 );
	  for (int i=0; i<nx; i++)
	    {
	      dX_prefix_sum[p][i+1] = dX_prefix_sum[p][i] + dX[p][i];
	      dX_prefix_sq[p][i+1] = dX_prefix_sq[p][i] + dX[p][i] * dX[p][i];
	    }
	}
    }

  auto full_search_best_index = [&]( int p ,
				     const Eigen::VectorXd & y_seg_z ) -> int
    {
      const int xlen = dX[p].size();
      const int yseg_len = y_seg_z.size();
      const int max_start = xlen - yseg_len;
      if ( max_start < 0 )
	return -1;

      std::vector<double> seg_a( xlen ) , seg_b( yseg_len );
      for (int i=0; i<xlen; i++) seg_a[i] = dX[p][i];
      for (int i=0; i<yseg_len; i++) seg_b[i] = y_seg_z[i];

      xcorr_t xc( seg_a , seg_b , 0 , 0 );

      double y_ss = 0;
      for (int i=0; i<yseg_len; i++)
	y_ss += y_seg_z[i] * y_seg_z[i];
      if ( y_ss <= 0 )
	return -1;

      double best_score = -1;
      int best_idx = -1;
      for (int i=0; i<(int)xc.lags.size(); i++)
	{
	  const int lag = xc.lags[i];
	  if ( lag < 0 || lag > max_start )
	    continue;
	  const double sx = dX_prefix_sum[p][lag + yseg_len] - dX_prefix_sum[p][lag];
	  const double sx2 = dX_prefix_sq[p][lag + yseg_len] - dX_prefix_sq[p][lag];
	  const double x_ss = sx2 - ( sx * sx ) / yseg_len;
	  if ( x_ss <= 0 )
	    continue;
	  const double score = std::fabs( xc.C[i] ) / std::sqrt( x_ss * y_ss );
	  if ( score > best_score )
	    {
	      best_score = score;
	      best_idx = lag;
	    }
	}

      return best_idx;
    };

  // default window range: middle 80% of the shorter of the two recordings
  // (avoids noisy leading/trailing segments at device start/stop)
  const int usable = std::min( nx, ny );
  const int ystart_default_sp = user_start ? (int)( param.requires_dbl("start") * sr ) : (int)( 0.10 * usable );
  const int yend_sp            = user_end   ? (int)( param.requires_dbl("end")   * sr ) : (int)( 0.90 * usable );

  double run_start_sec = ystart_default_sp / (double)sr;
  bool auto_tuned = false;
  if ( try_windows )
    {
      if ( euclidean )
	Helper::halt( "try-start/try-len/try-inc are currently only supported for xcorr mode" );

      std::vector<double> try_starts, try_lens, try_incs;
      if ( auto_try_mode )
	{
	  // Local window search around current/default settings:
	  // keep start fixed, vary len, and derive inc deterministically from len.
	  const double min_len = 30.0;
	  try_starts = { run_start_sec };
	  try_lens = {
	    min_len,
	    std::max( min_len, ylen_sec * 0.25 ),
	    std::max( min_len, ylen_sec * 0.5 ),
	    std::max( min_len, ylen_sec )
	  };

	  auto sort_unique = []( std::vector<double> & v ) {
	    std::sort( v.begin(), v.end() );
	    std::vector<double> u;
	    for (int i=0; i<(int)v.size(); i++)
	      if ( u.empty() || std::fabs( v[i] - u.back() ) > 1e-9 ) u.push_back( v[i] );
	    v.swap( u );
	  };
	  sort_unique( try_starts );
	  sort_unique( try_lens );
	  for (int i=0; i<(int)try_lens.size(); i++)
	    try_incs.push_back( std::max( 1.0, try_lens[i] / 5.0 ) );
	  logger << "  auto-try: local window search around defaults "
		 << "(fixed start=" << run_start_sec << "s, len~" << ylen_sec
		 << "s, inc=len/5, min len=" << min_len << "s)\n";
	}
      else
	{
	  try_starts = param.has( "try-start" ) ? param.dblvector( "try-start" ) : std::vector<double>( 1 , run_start_sec );
	  try_lens   = param.has( "try-len" )   ? param.dblvector( "try-len" )   : std::vector<double>( 1 , ylen_sec );
	  try_incs   = param.has( "try-inc" )   ? param.dblvector( "try-inc" )   : std::vector<double>( 1 , yinc_sec );
	}

      struct tune_result_t {
	bool valid = false;
	double score = -1;
	double start_sec = 0;
	double len_sec = 0;
	double inc_sec = 0;
	int total_windows = 0;
	int accepted_windows = 0;
	double p_ok = 0;
	double median_peak = 0;
	double mean_peak = 0;
	double r2 = 0;
	double slope = 0;
	double intercept = 0;
      };

      auto median_of = []( std::vector<double> v ) -> double {
	if ( v.empty() ) return 0;
	std::sort( v.begin() , v.end() );
	return v[ v.size() / 2 ];
      };

      auto try_one = [&]( double start_sec , double len_sec_try , double inc_sec_try ) -> tune_result_t
	{
	  tune_result_t tr;
	  tr.start_sec = start_sec;
	  tr.len_sec = len_sec_try;
	  tr.inc_sec = inc_sec_try;

	  const int ylen_try = len_sec_try * sr;
	  const int yinc_try = inc_sec_try * sr;
	  const int ystart_try0 = start_sec * sr;
	  const int na_try = nx - ylen_try + 1;
	  if ( ylen_try <= 0 || yinc_try <= 0 || na_try <= 0 ) return tr;

	  int ystart_try = ystart_try0;
	  int steps_try = 0;
	  std::vector<double> sec, win, peaks;

	  while ( 1 )
	    {
	      if ( ystart_try + ylen_try > ny || ystart_try > yend_sp ) break;
	      ++steps_try;
	      if ( steps_try > ysteps ) break;
	      ++tr.total_windows;

	      int mina = 0, maxa = na_try;
	      if ( constrained_offset )
		{
		  const int minoff = offsets[0] * sr;
		  const int maxoff = offsets[1] * sr;
		  mina = ystart_try - maxoff;
		  maxa = ystart_try - minoff;
		  if ( mina > maxa ) { int t = mina; mina = maxa; maxa = t; }
		  if ( mina < 0 )    mina = 0;
		  if ( maxa < 0 )    maxa = 0;
		  if ( mina >= na_try )  mina = na_try;
		  if ( maxa >= na_try )  maxa = na_try;
		  if ( mina >= maxa ) break;
		}

	      std::vector<double> pair_peaks( np, 0 );
	      std::vector<int> pair_lags( np, 0 );
	      std::vector<bool> pair_valid( np, false );
	      int n_valid = 0;

	  for (int p=0; p<np; p++)
		{
		  const Eigen::VectorXd eb = zscore( dY[p].segment( ystart_try , ylen_try ) );
		  int p_minidx = -1;

		  if ( full_search )
		    p_minidx = full_search_best_index( p , eb );
		  else
		    {
		      const int pcenter   = ( mina + maxa ) / 2;
		      const int margin_sp = ( maxa - mina ) / 2;
		      const int pstart = std::max( 0, pcenter );
		      const int pend   = std::min( nx, pcenter + ylen_try );
		      if ( pend - pstart < ylen_try / 2 ) continue;

		      Eigen::VectorXd ea = Eigen::VectorXd::Zero( ylen_try );
		      for (int i=0; i<pend-pstart; i++) ea[i] = dX[p][pcenter+i];
		      ea = zscore( ea );

		      std::vector<double> seg_a( ylen_try ) , seg_b( ylen_try );
		      for (int i=0; i<ylen_try; i++) { seg_a[i] = ea[i]; seg_b[i] = eb[i]; }

		      xcorr_t xc( seg_a , seg_b , margin_sp , 0 );
		      pair_lags[p] = xc.lags[xc.mx];
		      p_minidx = std::max( mina, std::min( maxa-1, pcenter + xc.lags[xc.mx] ) );
		    }

		  if ( p_minidx < 0 ) continue;
		  Eigen::VectorXd ax = dX[p].segment( p_minidx , ylen_try );
		  Eigen::VectorXd by = dY[p].segment( ystart_try , ylen_try );
		  const double ax_m = ax.mean();
		  const double ax_s = std::sqrt( (ax.array() - ax_m).square().mean() );
		  ax = ax.array() - ax_m;
		  if ( ax_s > 0 ) ax /= ax_s;
		  const double by_m = by.mean();
		  const double by_s = std::sqrt( (by.array() - by_m).square().mean() );
		  by = by.array() - by_m;
		  if ( by_s > 0 ) by /= by_s;
		  pair_peaks[p] = std::fabs( ( ax.array() * by.array() ).mean() );
		  if ( full_search ) pair_lags[p] = p_minidx;
		  pair_valid[p] = true;
		  ++n_valid;
		}

	      if ( n_valid == 0 ) { ystart_try += yinc_try; continue; }

	      std::vector<double> valid_peaks;
	      std::vector<int> valid_lags;
	      for (int p=0; p<np; p++)
		if ( pair_valid[p] )
		  {
		    valid_peaks.push_back( pair_peaks[p] );
		    valid_lags.push_back( pair_lags[p] );
		  }

	      std::sort( valid_peaks.begin() , valid_peaks.end() );
	      const double median_peak = valid_peaks[ valid_peaks.size() / 2 ];
	      peaks.push_back( median_peak );
	      if ( median_peak >= min_peak )
		{
		  std::sort( valid_lags.begin() , valid_lags.end() );
		  int minidx = full_search ? valid_lags[ valid_lags.size() / 2 ] : 0;
		  if ( ! full_search )
		    {
		      const int pcenter = ( mina + maxa ) / 2;
		      const int median_lag = valid_lags[ valid_lags.size() / 2 ];
		      minidx = pcenter + median_lag;
		      if ( minidx < mina ) minidx = mina;
		      if ( minidx >= maxa ) minidx = maxa - 1;
		    }
		  sec.push_back( ( ystart_try - minidx ) / (double)sr );
		  win.push_back( ystart_try / (double)sr );
		  ++tr.accepted_windows;
		}

	      ystart_try += yinc_try;
	    }

	  tr.p_ok = tr.total_windows > 0 ? tr.accepted_windows / (double)tr.total_windows : 0;
	  if ( ! peaks.empty() )
	    {
	      tr.median_peak = median_of( peaks );
	      for (int i=0; i<(int)peaks.size(); i++) tr.mean_peak += peaks[i];
	      tr.mean_peak /= peaks.size();
	    }
	  if ( sec.size() < 2 ) return tr;

	  const int nw = sec.size();
	  double sx=0, sy=0;
	  for (int i=0;i<nw;i++) { sx += win[i]; sy += sec[i]; }
	  sx /= nw; sy /= nw;
	  double sxy=0, sxx=0;
	  for (int i=0;i<nw;i++) { sxy += ( win[i] - sx ) * ( sec[i] - sy ); sxx += ( win[i] - sx ) * ( win[i] - sx ); }
	  const double slope0 = sxx > 0 ? sxy / sxx : 0;
	  const double intercept0 = sy - slope0 * sx;
	  std::vector<bool> keep( nw, true );
	  two_pass_trim_outliers( win , sec , &keep );

	  double sx2=0, sy2=0; int n2=0;
	  for (int i=0;i<nw;i++) if ( keep[i] ) { sx2 += win[i]; sy2 += sec[i]; ++n2; }
	  if ( n2 < 2 ) return tr;
	  sx2 /= n2; sy2 /= n2;
	  double sxy2=0, sxx2=0;
	  for (int i=0;i<nw;i++) if ( keep[i] ) { sxy2 += ( win[i] - sx2 ) * ( sec[i] - sy2 ); sxx2 += ( win[i] - sx2 ) * ( win[i] - sx2 ); }
	  tr.slope = sxx2 > 0 ? sxy2 / sxx2 : slope0;
	  tr.intercept = sy2 - tr.slope * sx2;

	  double mean_clean=0, ss_res=0, ss_tot=0; int nc=0;
	  for (int i=0;i<nw;i++) if ( keep[i] ) { mean_clean += sec[i]; ++nc; }
	  mean_clean /= nc;
	  for (int i=0;i<nw;i++) if ( keep[i] )
	    {
	      ss_res += std::pow( sec[i] - ( tr.intercept + tr.slope * win[i] ) , 2 );
	      ss_tot += std::pow( sec[i] - mean_clean , 2 );
	    }
	  tr.r2 = ss_tot > 0 ? 1.0 - ss_res / ss_tot : 1.0;
	  tr.score = std::max( 0.0 , tr.r2 ) * tr.p_ok * std::max( 0.0 , tr.median_peak );
	  tr.valid = true;
	  return tr;
	};

      tune_result_t best;
      logger << "  trying window settings:\n";
      for (int a=0; a<(int)try_starts.size(); a++)
	{
	  if ( auto_try_mode )
	    {
	      for (int b=0; b<(int)try_lens.size(); b++)
		{
		  tune_result_t tr = try_one( try_starts[a] , try_lens[b] , try_incs[b] );
		  logger << "    start=" << tr.start_sec << "s len=" << tr.len_sec << "s inc=" << tr.inc_sec
			 << "s  accepted=" << tr.accepted_windows << "/" << tr.total_windows
			 << "  P_OK=" << tr.p_ok
			 << "  peak=" << tr.median_peak
			 << "  R2=" << tr.r2
			 << "  score=" << tr.score
			 << ( tr.valid ? "" : "  [invalid]" ) << "\n";
		  if ( ! tr.valid ) continue;
		  const bool better =
		    ( ! best.valid ) ||
		    ( tr.score > best.score ) ||
		    ( tr.score == best.score && tr.r2 > best.r2 ) ||
		    ( tr.score == best.score && tr.r2 == best.r2 && tr.p_ok > best.p_ok ) ||
		    ( tr.score == best.score && tr.r2 == best.r2 && tr.p_ok == best.p_ok && tr.accepted_windows > best.accepted_windows );
		  if ( better ) best = tr;
		}
	    }
	  else
	    {
	      for (int b=0; b<(int)try_lens.size(); b++)
		for (int c=0; c<(int)try_incs.size(); c++)
		  {
		    tune_result_t tr = try_one( try_starts[a] , try_lens[b] , try_incs[c] );
		    logger << "    start=" << tr.start_sec << "s len=" << tr.len_sec << "s inc=" << tr.inc_sec
			   << "s  accepted=" << tr.accepted_windows << "/" << tr.total_windows
			   << "  P_OK=" << tr.p_ok
			   << "  peak=" << tr.median_peak
			   << "  R2=" << tr.r2
			   << "  score=" << tr.score
			   << ( tr.valid ? "" : "  [invalid]" ) << "\n";
		    if ( ! tr.valid ) continue;
		    const bool better =
		      ( ! best.valid ) ||
		      ( tr.score > best.score ) ||
		      ( tr.score == best.score && tr.r2 > best.r2 ) ||
		      ( tr.score == best.score && tr.r2 == best.r2 && tr.p_ok > best.p_ok ) ||
		      ( tr.score == best.score && tr.r2 == best.r2 && tr.p_ok == best.p_ok && tr.accepted_windows > best.accepted_windows );
		    if ( better ) best = tr;
		  }
	    }
	}

      if ( ! best.valid )
	Helper::halt( "could not find any valid INSERT fit across try-start/try-len/try-inc candidates" );

      run_start_sec = best.start_sec;
      ylen_sec = best.len_sec;
      yinc_sec = best.inc_sec;
      auto_tuned = true;
      logger << "  selected window settings: start=" << run_start_sec
	     << "s len=" << ylen_sec << "s inc=" << yinc_sec
	     << "s  (R2=" << best.r2 << ", P_OK=" << best.p_ok << ", peak=" << best.median_peak << ")\n";
    }

  const int ylen = ylen_sec * sr;
  const int yinc = yinc_sec * sr;
  const int na   = nx - ylen + 1;

  logger << "  method: " << ( euclidean ? "euclidean" : "xcorr" )
	 << ( do_filter ? ", bandpass " + Helper::dbl2str(filt_low) + "-" + Helper::dbl2str(filt_high) + " Hz" : "" )
	 << ";  " << ylen_sec << "s windows every " << yinc_sec << "s"
	 << ",  range " << run_start_sec << "-" << yend_sp/(double)sr << "s\n";

  int ystart = run_start_sec * sr;
  int steps = 0;

  // running offset tracker for full-search xcorr: seed from header, then updated
  // each window so accumulated drift (which can far exceed ±ylen) is tracked
  int running_offset_sp = header_offset_valid ? (int)( header_offset * sr ) : 0;

  std::vector<double> all_sec, all_win;       // combined per-window estimates
  std::vector<std::vector<double>> pair_sec( np ), pair_win( np );  // per-pair per-window
  std::vector<double> peak_win;
  std::vector<double> win_level_sec;
  std::vector<int> win_fit_idx;
  double first_sec = 0;
  bool first_win = true;
  int total_windows = 0;

  // base clock times for both EDFs (for labelling windows)
  clocktime_t base1( edf.header.startdate,  edf.header.starttime );
  clocktime_t base2( edf2.header.startdate, edf2.header.starttime );
  const bool clocks_valid = base1.valid && base2.valid;

  while ( 1 )
    {
      if ( ystart + ylen > ny || ystart > yend_sp )
	{
	  break;
	}

      ++steps;
      if ( steps > ysteps ) break;
      ++total_windows;

      writer.level( ystart / (double)sr , "WIN" );

      // search bounds in primary (offset = ystart - a, so a = ystart - offset)
      int mina = 0, maxa = na;
      if ( constrained_offset )
	{
	  const int minoff = offsets[0] * sr;
	  const int maxoff = offsets[1] * sr;
	  mina = ystart - maxoff;
	  maxa = ystart - minoff;
	  if ( mina > maxa ) { int t = mina; mina = maxa; maxa = t; }
	  if ( mina < 0 )    mina = 0;
	  if ( maxa < 0 )    maxa = 0;
	  if ( mina >= na )  mina = na;
	  if ( maxa >= na )  maxa = na;
	  if ( mina >= maxa )
	    Helper::halt( "offset-range results in an empty search space after clamping to signal bounds" );
	}

      int minidx = mina;
      bool win_ok = true;  // set false by xcorr quality gate if peak too low
      double win_peak = 0;

      if ( euclidean )
	{
	  std::vector<Eigen::VectorXd> sY(np);
	  for (int p=0; p<np; p++)
	    sY[p] = zscore( dY[p].segment( ystart, ylen ) );

	  // combined score across pairs
	  std::vector<double> st( na, 0 );
	  for (int p=0; p<np; p++)
	    for (int a=mina; a<maxa; a++)
	      st[a] += ( zscore( dX[p].segment(a,ylen) ) - sY[p] ).norm();

	  double minst = 0;
	  for (int a=mina; a<maxa; a++)
	    if ( a == mina || st[a] < minst ) { minst = st[a]; minidx = a; }

	  // per-pair best alignment
	  for (int p=0; p<np; p++)
	    {
	      int pidx = mina; double pst = 0;
	      for (int a=mina; a<maxa; a++)
		{
		  const double d = ( zscore( dX[p].segment(a,ylen) ) - sY[p] ).norm();
		  if ( a == mina || d < pst ) { pst = d; pidx = a; }
		}
	      pair_sec[p].push_back( ( ystart - pidx ) / (double)sr );
	      pair_win[p].push_back(   ystart           / (double)sr );
	    }

	  if ( verbose )
	    for (int a=mina; a<maxa; a++)
	      logger << "  " << a << "\t" << st[a] << "\n";
	}
      else
	{
	  // xcorr: if full-search is requested, do an actual global lag search
	  // over the full valid primary trace. Otherwise, search within the
	  // constrained local interval around the expected alignment.
	  const bool xcorr_full = full_search;
	  const int pcenter_full = std::min( std::max( ystart - running_offset_sp, 0 ), na - 1 );
	  const int pcenter   = xcorr_full ? pcenter_full : ( mina + maxa ) / 2;
	  const int margin_sp = xcorr_full ? ylen : ( maxa - mina ) / 2;

	  std::vector<int> pair_lags( np, 0 );
	  std::vector<double> pair_peaks( np, 0 );
	  std::vector<double> pair_offset_sec_tmp( np, 0 );
	  std::vector<bool> pair_valid( np, false );
	  for (int p=0; p<np; p++)
	    {
	      Eigen::VectorXd eb = zscore( dY[p].segment( ystart, ylen ) );
	      int p_minidx = -1;

	      if ( xcorr_full )
		{
		  p_minidx = full_search_best_index( p , eb );
		  if ( p_minidx < 0 )
		    continue;
		  pair_lags[p] = p_minidx;
		}
	      else
		{
		  // clamp primary segment to signal bounds
		  const int pstart = std::max( 0, pcenter );
		  const int pend   = std::min( nx, pcenter + ylen );
		  if ( pend - pstart < ylen / 2 )
		    {
		      logger << "  warning: insufficient primary data at window " << ystart/(double)sr << "s for pair "
			     << slab1[p] << "/" << slab2[p] << ", skipping\n";
		      continue;
		    }

		  Eigen::VectorXd ea = Eigen::VectorXd::Zero( ylen );
		  for (int i=0; i<pend-pstart; i++) ea[i] = dX[p][pcenter+i];
		  ea = zscore( ea );

		  std::vector<double> seg_a( ylen ), seg_b( ylen );
		  for (int i=0; i<ylen; i++) { seg_a[i] = ea[i]; seg_b[i] = eb[i]; }

		  xcorr_t xc( seg_a, seg_b, margin_sp, 0 );
		  pair_lags[p]  = xc.lags[xc.mx];
		  p_minidx = std::max( mina, std::min( maxa-1, pcenter + xc.lags[xc.mx] ) );

		}

	      Eigen::VectorXd ax = dX[p].segment( p_minidx , ylen );
	      Eigen::VectorXd by = dY[p].segment( ystart , ylen );
	      const double ax_m = ax.mean();
	      const double ax_s = std::sqrt( (ax.array() - ax_m).square().mean() );
	      ax = ax.array() - ax_m;
	      if ( ax_s > 0 ) ax /= ax_s;
	      const double by_m = by.mean();
	      const double by_s = std::sqrt( (by.array() - by_m).square().mean() );
	      by = by.array() - by_m;
	      if ( by_s > 0 ) by /= by_s;
	      pair_peaks[p] = std::fabs( ( ax.array() * by.array() ).mean() );
	      pair_offset_sec_tmp[p] = ( ystart - p_minidx ) / (double)sr;
	      pair_valid[p] = true;

	      if ( verbose )
		logger << "  xcorr " << slab1[p] << " x " << slab2[p]
		       << "  waveform_shift=" << pair_offset_sec_tmp[p] << "s  peak=" << pair_peaks[p] << "\n";
	    }

	  // median peak across pairs — gate window quality
	  std::vector<double> sorted_peaks;
	  std::vector<int> sorted_lags;
	  for (int p=0; p<np; p++)
	    if ( pair_valid[p] )
	      {
		sorted_peaks.push_back( pair_peaks[p] );
		sorted_lags.push_back( pair_lags[p] );
	      }
	  if ( sorted_peaks.empty() )
	    {
	      win_ok = false;
	      win_peak = 0;
	      peak_win.push_back( 0 );
	    }
	  else
	    {
	      std::sort( sorted_peaks.begin(), sorted_peaks.end() );
	      const double median_peak = sorted_peaks[ sorted_peaks.size() / 2 ];
	      win_peak = median_peak;
	      peak_win.push_back( median_peak );
	      win_ok = median_peak >= min_peak;
	      if ( ! win_ok && verbose )
		logger << "  window " << ystart/(double)sr << "s: skipped (peak=" << median_peak << " < " << min_peak << ")\n";

	      // median lag across pairs
	      std::sort( sorted_lags.begin(), sorted_lags.end() );
	      if ( xcorr_full )
		minidx = sorted_lags[ sorted_lags.size() / 2 ];
	      else
		{
		  const int median_lag = sorted_lags[ sorted_lags.size() / 2 ];
		  minidx = pcenter + median_lag;
		  if ( minidx < mina ) minidx = mina;
		  if ( minidx >= maxa ) minidx = maxa - 1;
		}

	      // update running offset tracker so next window is centered correctly;
	      // update even for low-quality windows as long as the implied offset shift
	      // is plausible (< half a window length) — prevents tracker from going stale
	      // when few windows pass the quality gate
	      if ( xcorr_full )
		{
		  const int new_offset = ystart - minidx;
		  if ( std::abs( new_offset - running_offset_sp ) <= ylen / 2 )
		    running_offset_sp = new_offset;
		}

	      // Keep per-pair regression on the exact same accepted windows used by
	      // the global summary, so one-pair runs remain internally consistent.
	      if ( win_ok )
		for (int p=0; p<np; p++)
		  if ( pair_valid[p] )
		    {
		      pair_sec[p].push_back( pair_offset_sec_tmp[p] );
		      pair_win[p].push_back( ystart / (double)sr );
		    }
	    }
	}

      if ( verbose && win_ok )
	logger << "  window " << ystart/(double)sr << "s: waveform_shift = "
	       << ( ystart - minidx ) / (double)sr << "s\n";

      const double offset_sec = ( ystart - minidx ) / (double)sr;
      // Report the matched-feature clock difference between secondary and primary.
      // `offset_sec` is the xcorr shift needed to align the secondary waveform back
      // onto the primary, so it must be subtracted from the raw header offset rather
      // than added; otherwise trimmed-leading data get double-counted in both terms.
      const double total_offset_sec = ( header_offset_valid ? header_offset : 0.0 ) - offset_sec;
      const double t1_sec     = minidx / (double)sr;
      const double t2_sec     = ystart / (double)sr;

      // only accumulate high-quality windows for drift regression
      if ( win_ok )
	{
	  if ( first_win ) { first_sec = offset_sec; first_win = false; }
	  all_sec.push_back( offset_sec );
	  all_win.push_back( t2_sec );
	  win_fit_idx.push_back( (int)all_sec.size() - 1 );
	}
      else
	win_fit_idx.push_back( -1 );
      win_level_sec.push_back( t2_sec );

      writer.value( "OK",    (int)win_ok );
      writer.value( "PEAK",  win_peak );
      writer.value( "SP",    ystart - minidx );
      writer.value( "SEC",   offset_sec );
      writer.value( "TOT_SEC", total_offset_sec );
      writer.value( "DSEC",  first_win ? 0.0 : offset_sec - first_sec );
      writer.value( "T1_SEC", t1_sec );
      writer.value( "T2_SEC", t2_sec );

      if ( clocks_valid )
	{
	  clocktime_t ct1 = base1;  ct1.advance_seconds( t1_sec );
	  clocktime_t ct2 = base2;  ct2.advance_seconds( t2_sec );
	  writer.value( "T1_HMS", ct1.as_string(':') );
	  writer.value( "T2_HMS", ct2.as_string(':') );
	}

      ystart += yinc;
    }
  writer.unlevel( "WIN" );

  // Per-window slope-fit status for WIN output:
  //   FIT_USED=1 if this window participated in slope fitting, else 0
  //   FIT_OUTLIER=1 if removed as outlier from fit, 0 if retained, -1 if not used
  std::vector<int> win_fit_used( total_windows , 0 );
  std::vector<int> win_fit_outlier( total_windows , -1 );
  for (int w=0; w<total_windows; w++)
    if ( win_fit_idx[w] >= 0 )
      {
	win_fit_used[w] = 1;
	win_fit_outlier[w] = 0;
      }

  double est_slope = 0, est_intercept = 0;
  bool est_valid = false;
  writer.value( "N_WIN_ALL", total_windows );
  writer.value( "USED_START_SEC", run_start_sec );
  writer.value( "USED_LEN_SEC", ylen_sec );
  writer.value( "USED_INC_SEC", yinc_sec );
  writer.value( "AUTO_TUNED", auto_tuned ? 1 : 0 );
  writer.value( "HDR_OFFSET", header_offset_valid ? header_offset : 0.0 );
  writer.value( "HDR_OFFSET_VALID", header_offset_valid ? 1 : 0 );

  // summary statistics across all windows
  if ( ! all_sec.empty() )
    {
      const int nw = all_sec.size();

      // helper: fit OLS on a subset indicated by 'keep'
      auto ols = [&]( const std::vector<bool> & keep )
	-> std::pair<double,double>  // slope, intercept
	{
	  double sx=0, sy=0; int n=0;
	  for (int i=0;i<nw;i++) if(keep[i]) { sx+=all_win[i]; sy+=all_sec[i]; ++n; }
	  if (n<2) return {0,0};
	  sx/=n; sy/=n;
	  double sxy=0, sxx=0;
	  for (int i=0;i<nw;i++) if(keep[i])
	    { sxy+=(all_win[i]-sx)*(all_sec[i]-sy); sxx+=(all_win[i]-sx)*(all_win[i]-sx); }
	  const double b = sxx>0 ? sxy/sxx : 0;
	  return { b , sy - b*sx };
	};

      // initial fit on all windows
      std::vector<bool> keep( nw, true );
      auto [slope0, intercept0] = ols( keep );

	      // two-pass flagging of outliers: |residual| > 3 * SD
	      const int n_outlier = two_pass_trim_outliers( all_win , all_sec , &keep );

      for (int w=0; w<total_windows; w++)
	{
	  const int idx = win_fit_idx[w];
	  if ( idx >= 0 && idx < nw )
	    win_fit_outlier[w] = keep[idx] ? 0 : 1;
	}

      auto [slope, intercept] = n_outlier > 0 ? ols(keep) : std::make_pair(slope0,intercept0);
      est_slope = slope;  est_intercept = intercept;  est_valid = true;

      // R² on clean set
      double ss_res=0, ss_tot=0, mean_sec_clean=0; int nc=0;
      for (int i=0;i<nw;i++) if(keep[i]) { mean_sec_clean+=all_sec[i]; ++nc; }
      if (nc>0) mean_sec_clean/=nc;
      for (int i=0;i<nw;i++) if(keep[i])
	{
	  ss_res += std::pow( all_sec[i] - (intercept + slope*all_win[i]) , 2 );
	  ss_tot += std::pow( all_sec[i] - mean_sec_clean , 2 );
	}
      const double r2 = ss_tot > 0 ? 1.0 - ss_res/ss_tot : 1.0;

      // raw distribution stats (all windows, before outlier removal)
      std::vector<double> sorted_sec = all_sec;
      std::sort( sorted_sec.begin(), sorted_sec.end() );
      const double min_sec    = sorted_sec.front();
      const double max_sec    = sorted_sec.back();
      const double median_sec = sorted_sec[ nw / 2 ];
      double mean_sec = 0;
      for (int i=0;i<nw;i++) mean_sec += all_sec[i];
      mean_sec /= nw;
      const double p_ok = total_windows > 0 ? nw / (double)total_windows : 0;
      double median_peak = 0, mean_peak = 0, min_peak_obs = 0, max_peak_obs = 0;
      if ( ! peak_win.empty() )
	{
	  std::vector<double> sorted_peak = peak_win;
	  std::sort( sorted_peak.begin(), sorted_peak.end() );
	  min_peak_obs = sorted_peak.front();
	  max_peak_obs = sorted_peak.back();
	  median_peak = sorted_peak[ sorted_peak.size() / 2 ];
	  for (int i=0; i<(int)peak_win.size(); i++) mean_peak += peak_win[i];
	  mean_peak /= peak_win.size();
	}

      // human-readable drift and implied sample rate
      const double slope_per_hr = slope * 3600.0;
      const double implied_sr   = sr / ( 1.0 - slope );

      // Same convention as per-window TOT_SEC above: total matched-feature clock
      // difference equals header offset minus the xcorr-derived intercept.
	      const double total_offset = ( header_offset_valid ? header_offset : 0.0 ) - intercept;
	      logger << "  summary across " << nw << " window(s)"
		     << ( n_outlier>0 ? " (" + Helper::int2str(n_outlier) + " outlier(s) removed from slope fit)" : "" ) << ":\n"
		     << "    quality          accepted=" << nw << "/" << total_windows << " (" << p_ok * 100.0
		     << "%)  peak median=" << median_peak << "  mean=" << mean_peak
		     << "  min=" << min_peak_obs << "  max=" << max_peak_obs << "\n"
		     << "    waveform_shift   median=" << median_sec << "s  mean=" << mean_sec
		     << "s  min=" << min_sec << "s  max=" << max_sec << "s  range=" << max_sec-min_sec << "s\n"
		     << "    offset           " << total_offset << "s"
		     << " (start_shift=" << intercept << "s, header_offset="
		     << ( header_offset_valid ? header_offset : 0.0 ) << "s)\n"
		     << "    drift            slope=" << slope << " s/s  (" << slope_per_hr << " s/hr)"
		     << "  intercept=" << intercept << "s  R2=" << r2 << "\n"
		     << "    implied SR of secondary: " << implied_sr << " Hz  (nominal: " << sr << " Hz)\n"
		     << "    (positive slope = secondary clock running faster than primary)\n";

      writer.value( "N_WIN",       nw );
      writer.value( "N_OUTLIER",   n_outlier );
      writer.value( "P_OK",        p_ok );
      writer.value( "MEDIAN_SEC",  median_sec );
      writer.value( "MEAN_SEC",    mean_sec );
      writer.value( "MIN_SEC",     min_sec );
      writer.value( "MAX_SEC",     max_sec );
      writer.value( "RANGE_SEC",   max_sec - min_sec );
      writer.value( "MEDIAN_PEAK", median_peak );
      writer.value( "MEAN_PEAK",   mean_peak );
      writer.value( "MIN_PEAK",    min_peak_obs );
      writer.value( "MAX_PEAK",    max_peak_obs );
      writer.value( "SLOPE",       slope );
      writer.value( "SLOPE_HR",    slope_per_hr );
      writer.value( "INTERCEPT",     intercept );
      writer.value( "TOTAL_OFFSET", total_offset );
      writer.value( "R2",          r2 );
      writer.value( "IMPLIED_SR",  implied_sr );

      bool quality_ok = true;
      std::vector<std::string> warn_reasons;
      if ( ! euclidean )
	{
	  if ( r2 < warn_r2 ) { quality_ok = false; warn_reasons.push_back( "R2=" + Helper::dbl2str( r2 ) + " < " + Helper::dbl2str( warn_r2 ) ); }
	  if ( p_ok < warn_p_ok ) { quality_ok = false; warn_reasons.push_back( "P_OK=" + Helper::dbl2str( p_ok ) + " < " + Helper::dbl2str( warn_p_ok ) ); }
	  if ( median_peak < warn_peak ) { quality_ok = false; warn_reasons.push_back( "median peak=" + Helper::dbl2str( median_peak ) + " < " + Helper::dbl2str( warn_peak ) ); }
	}
      writer.value( "OKAY", quality_ok ? 1 : 0 );
      if ( do_warn && ! quality_ok )
	{
	  logger << "  warning: alignment quality may be poor: ";
	  for (int i=0; i<(int)warn_reasons.size(); i++)
	    {
	      if ( i ) logger << "; ";
	      logger << warn_reasons[i];
	    }
	  logger << "\n"
		 << "  hint: try a smaller len window; also try a wider offset-range"
		 << " (e.g. offset-range=-360,360) or full-search\n";
	}

      // per-pair slope fitting (same OLS + outlier removal)
      logger << "  per-pair drift:\n";
      for (int p=0; p<np; p++)
	{
	  const std::string chs = slab1[p] + ".." + slab2[p];
	  const int pnw = pair_sec[p].size();
	  if ( pnw < 2 ) continue;

	  // initial OLS
	  std::vector<bool> pk( pnw, true );
	  double psx=0, psy=0;
	  for (int i=0;i<pnw;i++) { psx+=pair_win[p][i]; psy+=pair_sec[p][i]; }
	  psx/=pnw; psy/=pnw;
	  double psxy=0, psxx=0;
	  for (int i=0;i<pnw;i++)
	    { psxy+=(pair_win[p][i]-psx)*(pair_sec[p][i]-psy);
	      psxx+=(pair_win[p][i]-psx)*(pair_win[p][i]-psx); }
	  double pb0 = psxx>0 ? psxy/psxx : 0;
	  double pa0 = psy - pb0*psx;

	      // two-pass outlier removal
	      const int pnout = two_pass_trim_outliers( pair_win[p] , pair_sec[p] , &pk );

	  // refit on clean set
	  double psx2=0,psy2=0; int pnc=0;
	  for (int i=0;i<pnw;i++) if(pk[i]) { psx2+=pair_win[p][i]; psy2+=pair_sec[p][i]; ++pnc; }
	  double pb=pb0, pa=pa0;
	  if ( pnc>=2 )
	    {
	      psx2/=pnc; psy2/=pnc;
	      double pxy2=0,pxx2=0;
	      for (int i=0;i<pnw;i++) if(pk[i])
		{ pxy2+=(pair_win[p][i]-psx2)*(pair_sec[p][i]-psy2);
		  pxx2+=(pair_win[p][i]-psx2)*(pair_win[p][i]-psx2); }
	      pb = pxx2>0 ? pxy2/pxx2 : pb0;
	      pa = psy2 - pb*psx2;
	    }

	  const double p_slope_hr  = pb * 3600.0;
	  const double p_implied_sr = sr / ( 1.0 - pb );

	  logger << "    " << chs << ":  slope=" << pb << " s/s (" << p_slope_hr << " s/hr)"
		 << "  intercept=" << pa << "s"
		 << "  implied SR=" << p_implied_sr << " Hz"
		 << ( pnout>0 ? "  [" + Helper::int2str(pnout) + " outlier(s) removed]" : "" ) << "\n";

	  writer.level( chs , "CHS" );
	  writer.value( "SLOPE",      pb );
	  writer.value( "SLOPE_HR",   p_slope_hr );
	  writer.value( "INTERCEPT",  pa );
	  writer.value( "IMPLIED_SR", p_implied_sr );
	  writer.value( "N_OUTLIER",  pnout );
	  writer.unlevel( "CHS" );
	}
    }
  else
    {
      // all windows failed quality gate — report diagnostics so user knows why
      writer.value( "N_WIN", 0 );
      writer.value( "P_OK",  0.0 );
      writer.value( "OKAY",  0 );
      if ( ! peak_win.empty() )
	{
	  std::vector<double> sp = peak_win;
	  std::sort( sp.begin(), sp.end() );
	  const double med = sp[ sp.size() / 2 ];
	  double mn = 0;
	  for (auto v : sp) mn += v;
	  mn /= sp.size();
	  writer.value( "MEDIAN_PEAK", med );
	  writer.value( "MEAN_PEAK",   mn );
	  writer.value( "MIN_PEAK",    sp.front() );
	  writer.value( "MAX_PEAK",    sp.back() );
	  logger << "  warning: 0/" << total_windows << " windows passed quality gate"
		 << " (min-peak=" << min_peak << ");"
		 << " observed peaks: median=" << med << " mean=" << mn
		 << " min=" << sp.front() << " max=" << sp.back() << "\n"
		 << "  hint: try min-peak=" << Helper::dbl2str_fixed(sp.back()*0.9,2)
		 << ", a smaller len window, a wider offset-range"
		 << " (e.g. offset-range=-360,360), or full-search;"
		 << " add verbose=1 to diagnose per-window results\n";
	}
      else
	logger << "  warning: 0/" << total_windows << " windows passed quality gate (no peak data)\n"
	       << "  hint: try a smaller len window, a wider offset-range"
	       << " (e.g. offset-range=-360,360), or full-search;"
	       << " add verbose=1 to diagnose per-window results\n";
    }

  // Add per-window drift-fit flags after computing final outlier mask.
  for (int w=0; w<total_windows; w++)
    {
      writer.level( win_level_sec[w] , "WIN" );
      writer.value( "FIT_USED",    win_fit_used[w] );
      writer.value( "FIT_OUTLIER", win_fit_outlier[w] );
      writer.unlevel( "WIN" );
    }

  // opt-in: splice secondary into primary using estimated (or manually overridden) offset + drift
  if ( do_insert )
    {
      if ( ! est_valid && ! param.has( "offset" ) )
	Helper::halt( "INSERT insert: no valid estimation results and no offset= provided" );

      const double use_offset = param.has( "offset" ) ? param.requires_dbl( "offset" ) : est_intercept;
      const bool   has_drift  = param.has( "drift" ) || est_valid;
      double use_drift = 0;
      if      ( param.has( "drift" ) ) use_drift = param.requires_dbl( "drift" );
      else if ( est_valid )            use_drift = est_slope;

      logger << "  inserting secondary EDF: offset=" << use_offset << "s"
	     << ", drift=" << use_drift << " s/s"
	     << ( param.has("offset") || param.has("drift") ? " [manual override]" : " [from fit]" ) << "\n";

      insert( edf , edf2 , insert_sig , use_offset , has_drift ? &use_drift : NULL , insert_annot );
    }

}



void edf_inserter_t::insert( edf_t & edf , edf_t & edf2 , const std::string & siglabel ,
			     const double offset , const double * fac, 
			     const std::string annot_label )
{
  
  // insert as much of the signals from edf2 into edf
  // assume both are (effectively) continuous
  // by default, align at 0, i.e. start of edf2 signal equals start of edf, and add as much as we can  

  // if an offset is specified, then we insert after adding an offset
  //   -ve offset implies EDF2 start is AFTER EDF start (i.e. it needs to be shifted backwards) so that starts align
  //   +ve offset implies the reverse:  EDF2 start is BEFORE EDF start, needs to go forward 

  //            S
  // EDF        |-----------------------|

  // EDF2       |-----------------------|      offset = 0
  // EDF2       |-----------------|000000      offset = 0 , pad with zeros and add annotation to indicate missing signal

  // EDF2       000|--------------------|XXX|  offset = -ve (i.e. EDF2 start is after EDF start): pad w/ zeros; truncate at end X as needed
  // EDF2   |XXX|--------------------|000      offset = +ve (i.e. EDF2 start is before EDF start), need to shift forwards 

  
  // optionally, we can add annotations to indicate where the signal is missing from edf2 

  // if fac (secs/shift) is non-null, then apply this timestretch factor to the inserted channel
  // i.e. this is to adjust for linear difference in clock rates.

  const bool timestretch = fac != NULL; 
  
  const bool no_annotations = true;

  signal_list_t signals = edf2.header.signal_list( siglabel , no_annotations );

  const int ns = signals.size();
  
  logger << "  inserting " << ns << " signals from " << edf2.filename << ", ";
  logger << "using an offset of " << offset << " seconds\n";

  if ( timestretch )
    {
      if ( *fac > 0 ) 
	logger << "  shrinking secondary signals by a rate of " <<  (*fac) << " sec per second\n";
      else
	logger << "  stretching secondary signals by a rate of " << -1 * (*fac) << " sec per second\n";
    }
  


  for (int s=0; s<ns; s++)
    {
      // get SR of second signal (EDF2)
      const int Fs = edf2.header.sampling_freq( signals(s) );
      
      // new signal, with length to match EDF  
      const int np = edf.header.nr * edf.header.record_duration * Fs;

      // putative new signal name
      std::string sig = signals.label(s);
      
      // check new signal name will be unique
      if ( edf.header.has_signal( sig ) )
	{
	  int j = 1;
	  while ( 1 )
	    {
	      std::string sig2 = sig + "." + Helper::int2str( j );
	      if ( ! edf.header.has_signal( sig2 ) )
		{
		  sig = sig2;
		  break;
		}
	      ++j;
	    }
	}

      // make a new vector, set to zero-pad
      std::vector<double> d1( np , 0 );
      
      // pull the EDF2 signal
      slice_t slice( edf2 , signals(s) , edf2.timeline.wholetrace() );      
      std::vector<double> d2 = *slice.pdata();
      
      // calculate best offset in sample points
      const int offset_sp = offset * Fs;

      // time-stretch?
      if ( timestretch )
	{

	  // if a time-stretch factor is defined, then
	  // use spline interpolation to stretch or shrink
	  // the secondary signal as needed (based on
	  // the ratio of secs and shift as specified on
	  // the INSERT command line
	  
	  const int n_orig = d2.size();
	  const int n_scaled = n_orig - n_orig * (*fac);
	  
	  if ( n_scaled <= 0 ) Helper::halt( "rescaled signal not defined" );

	  std::vector<double> t( n_orig );
	  for (int i=0; i<n_orig; i++) t[i] = i;
	  
	  tk::spline spline;
	  spline.set_points( t, d2 );
	  
	  d2.clear();
	  d2.resize( n_scaled , 0 );
	  
	  for (int i=0; i<n_scaled; i++)
	    d2[i] = spline( n_orig * ( i / (double)n_scaled ) );	  
	  
	}
            
      // console messages
      logger << "  inserting " << sig << " ( SR = " << Fs << " Hz, offset = " << offset_sp << " samples ) into primary EDF\n";

      // -ve offset : pad w/ new N zeros , skip last N 
      //   ||||||||     d1
      //     ||||||||   d2   <---- need to shift back, aka start early 
      //   00IIIIII
      
      // +ve offset : pad w/ new N zeros , skip first N 
      //   ||||||||
      // ||||||||          -----> need to shift forward, but means start late in 2ndary
      //    IIIIII00
      
      // pointers (sample points to both files)
      int p1 = 0 , p2 = offset_sp ;
      
      // signal lengths
      const int n1 = np;
      const int n2 = d2.size();

      // advance sp in 2ndary signal 
      while ( 1 )
	{

	  // all done?
	  if ( p1 == n1 ) break;
	  
	  // zero-pad if out of input
	  if ( p2 >= n2 ) 
	    d1[ p1 ] = 0 ;
	  else if ( p2 >= 0 ) // else add or 0-pad (i.e. not yet at start of EDF1)
	    d1[ p1 ] = d2[ p2 ];
	  
	  // advance
	  ++p1;
	  ++p2;	  	  
	}

      // add the new signal
      edf.add_signal( sig , Fs , d1 );
      
      // add any annotations 
      // -- todo --
      
    }
  
  
  
}



  
