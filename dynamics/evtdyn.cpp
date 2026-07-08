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

//    --------------------------------------------------------------------
//
//    Event-level dynamics summaries
//
//    --------------------------------------------------------------------

#include "dynamics/evtdyn.h"

#include "edf/edf.h"
#include "annot/annot.h"
#include "param.h"
#include "defs/defs.h"
#include "db/db.h"
#include "helper/helper.h"
#include "helper/logger.h"

#include <algorithm>
#include <iomanip>
#include <cmath>
#include <limits>
#include <numeric>

extern writer_t writer;
extern logger_t logger;

namespace {

bool is_finite( const double x )
{
  return std::isfinite( x );
}

double mean( const std::vector<double> & x )
{
  if ( x.empty() ) return 0;
  double s = 0;
  for (int i=0;i<x.size();i++) s += x[i];
  return s / (double)x.size();
}

double median( std::vector<double> x )
{
  if ( x.empty() ) return 0;
  std::sort( x.begin() , x.end() );
  const int n = x.size();
  return n % 2 ? x[n/2] : 0.5 * ( x[n/2-1] + x[n/2] );
}

double sd( const std::vector<double> & x )
{
  if ( x.size() < 2 ) return 0;
  const double m = mean( x );
  double ss = 0;
  for (int i=0;i<x.size();i++) ss += ( x[i] - m ) * ( x[i] - m );
  return sqrt( ss / (double)( x.size() - 1 ) );
}

double corr( const std::vector<double> & x , const std::vector<double> & y )
{
  if ( x.size() != y.size() || x.size() < 2 ) return 0;
  const double mx = mean( x );
  const double my = mean( y );
  double sxx = 0, syy = 0, sxy = 0;
  for (int i=0;i<x.size();i++)
    {
      const double dx = x[i] - mx;
      const double dy = y[i] - my;
      sxx += dx * dx;
      syy += dy * dy;
      sxy += dx * dy;
    }
  return sxx <= 0 || syy <= 0 ? 0 : sxy / sqrt( sxx * syy );
}

double slope( const std::vector<double> & x , const std::vector<double> & y )
{
  if ( x.size() != y.size() || x.size() < 2 ) return 0;
  const double mx = mean( x );
  const double my = mean( y );
  double sxx = 0, sxy = 0;
  for (int i=0;i<x.size();i++)
    {
      const double dx = x[i] - mx;
      sxx += dx * dx;
      sxy += dx * ( y[i] - my );
    }
  return sxx <= 0 ? 0 : sxy / sxx;
}

std::string fmt_num( const double x )
{
  if ( ! is_finite( x ) ) return "NA";
  const double r = std::round( x );
  if ( std::fabs( x - r ) < 1e-9 )
    {
      std::ostringstream ss;
      ss << (long long)r;
      return ss.str();
    }
  std::ostringstream ss;
  ss.setf( std::ios::fixed );
  ss << std::setprecision( 3 ) << x;
  std::string s = ss.str();
  while ( s.size() > 1 && s.find( '.' ) != std::string::npos && s.back() == '0' ) s.pop_back();
  if ( ! s.empty() && s.back() == '.' ) s.pop_back();
  return s;
}

double log2_ratio( const double num , const double denom )
{
  if ( ! is_finite( num ) || ! is_finite( denom ) || num <= 0 || denom <= 0 ) return 0;
  return log( num / denom ) / log( 2.0 );
}

void write_log2_ratio( const std::string & name , const double num , const double denom )
{
  if ( is_finite( num ) && is_finite( denom ) && num > 0 && denom > 0 )
    writer.value( name , log2_ratio( num , denom ) );
}

std::vector<double> clean_values( const std::vector<double> & x ,
				  const evtdyn_opts_t & opt ,
				  const bool do_log )
{
  std::vector<double> y;
  for (int i=0;i<x.size();i++)
    if ( is_finite( x[i] ) ) y.push_back( x[i] );

  if ( y.empty() ) return y;

  if ( opt.winsor_p > 0 && opt.winsor_p < 0.5 && y.size() > 2 )
    {
      std::vector<double> z = y;
      std::sort( z.begin() , z.end() );
      const int l = std::max( 0 , (int)floor( opt.winsor_p * ( z.size() - 1 ) ) );
      const int u = std::min( (int)z.size()-1 , (int)ceil( ( 1 - opt.winsor_p ) * ( z.size() - 1 ) ) );
      for (int i=0;i<y.size();i++)
	{
	  if ( y[i] < z[l] ) y[i] = z[l];
	  if ( y[i] > z[u] ) y[i] = z[u];
	}
    }

  if ( do_log )
    {
      for (int i=0;i<y.size();i++)
	y[i] = y[i] > 0 ? log( y[i] ) : std::numeric_limits<double>::quiet_NaN();
      std::vector<double> yy;
      for (int i=0;i<y.size();i++)
	if ( is_finite( y[i] ) ) yy.push_back( y[i] );
      y = yy;
    }

  if ( opt.rank )
    {
      std::vector<std::pair<double,int> > r;
      for (int i=0;i<y.size();i++)
	if ( is_finite( y[i] ) ) r.push_back( std::make_pair( y[i] , i ) );
      std::sort( r.begin() , r.end() );
      std::vector<double> yy( y.size() , std::numeric_limits<double>::quiet_NaN() );
      for (int i=0;i<r.size();i++) yy[ r[i].second ] = i + 1;
      y = yy;
    }

  if ( opt.z )
    {
      const double m = mean( y );
      const double s = sd( y );
      if ( s > 0 )
	for (int i=0;i<y.size();i++) y[i] = ( y[i] - m ) / s;
    }

  return y;
}

std::vector<int> clean_keep_idx( const std::vector<double> & x ,
				 const evtdyn_opts_t & opt ,
				 const bool do_log )
{
  std::vector<int> keep;
  for (int i=0;i<x.size();i++)
    {
      if ( ! is_finite( x[i] ) ) continue;
      if ( do_log && x[i] <= 0 ) continue;
      keep.push_back( i );
    }
  return keep;
}

template<typename T>
std::vector<T> apply_keep( const std::vector<T> & x , const std::vector<int> & keep )
{
  std::vector<T> y;
  y.reserve( keep.size() );
  for (int i=0;i<keep.size();i++)
    y.push_back( x[ keep[i] ] );
  return y;
}

struct corr_data_t
{
  std::vector<int> rows;
  std::vector<double> values;
};

}

evtdyn_opts_t::evtdyn_opts_t()
{
  prefix = "";
  bg_none = false;
  hypno = true;
  verbose = false;
  anchor = "MID";
  cluster_sec = 10;
  train_gap_sec = 10;
  train_min = 3;
  refr_lwr = 0;
  refr_upr = 3;
  excit_lwr = 4;
  excit_upr = 9;
  lag_min_sec = 1.0;
  lag_max_sec = 60;
  lag_lin = 0;
  lag_log = 20;
  lag_hist = true;
  winsor_p = 0;
  z = false;
  rank = false;
  min_n = 5;
  min_seg_n = 3;
  min_bg_min = 1.0;
  min_seg_bg_min = 5.0;
}

void evtdyn_t::init( edf_t & edf0 , const param_t & param , const std::string & prefix )
{
  edf = &edf0;
  opt = evtdyn_opts_t();
  opt.prefix = prefix;

  const std::string p = prefix;
  opt.bg_none = param.has( p + "bg-none" ) ? param.yesno( p + "bg-none" ) : false;
  opt.hypno = param.has( p + "hypno" ) ? param.yesno( p + "hypno" ) : true;
  opt.verbose = param.has( p + "verbose" );
  opt.bg = param.has( p + "bg" ) ? param.strvector( p + "bg" ) : std::vector<std::string>();
  opt.vars = param.has( p + "vars" ) ? param.strset( p + "vars" ) : std::set<std::string>();
  if ( param.has( p + "corr" ) ) opt.corr = param.strset( p + "corr" );
  if ( param.has( p + "corr1" ) ) opt.corr1 = param.strset( p + "corr1" );
  if ( param.has( p + "corr2" ) ) opt.corr2 = param.strset( p + "corr2" );
  if ( param.has( p + "log" ) && ! param.empty( p + "log" ) ) opt.log_vars = param.strset( p + "log" );
  opt.anchor = param.has( p + "anchor" ) ? param.value( p + "anchor" , true ) : "MID";
  if ( param.has( p + "cluster" ) ) opt.cluster_sec = param.requires_dbl( p + "cluster" );
  if ( param.has( p + "train-gap" ) ) opt.train_gap_sec = param.requires_dbl( p + "train-gap" );
  if ( param.has( p + "train-min" ) ) opt.train_min = param.requires_int( p + "train-min" );
  if ( param.has( p + "refractory" ) )
    {
      std::vector<double> w = param.dblvector( p + "refractory" );
      if ( w.size() == 2 ) { opt.refr_lwr = w[0]; opt.refr_upr = w[1]; }
    }
  if ( param.has( p + "excitatory" ) )
    {
      std::vector<double> w = param.dblvector( p + "excitatory" );
      if ( w.size() == 2 ) { opt.excit_lwr = w[0]; opt.excit_upr = w[1]; }
    }
  if ( param.has( p + "lag-min" ) ) opt.lag_min_sec = param.requires_dbl( p + "lag-min" );
  if ( param.has( p + "lag-max" ) ) opt.lag_max_sec = param.requires_dbl( p + "lag-max" );
  if ( opt.lag_min_sec >= opt.lag_max_sec )
    Helper::halt( "lag-min (" + Helper::dbl2str( opt.lag_min_sec )
		  + ") must be less than lag-max (" + Helper::dbl2str( opt.lag_max_sec ) + ")" );
  if ( param.has( p + "lag-log" ) ) opt.lag_log = param.requires_int( p + "lag-log" );
  if ( param.has( p + "lag-lin" ) ) opt.lag_lin = param.requires_int( p + "lag-lin" );
  if ( param.has( p + "lag-hist" ) ) opt.lag_hist = param.yesno( p + "lag-hist" );
  if ( param.has( p + "winsor" ) ) opt.winsor_p = param.requires_dbl( p + "winsor" );
  opt.z = param.has( p + "z" );
  opt.rank = param.has( p + "rank" );
  if ( param.has( p + "min-n" ) ) opt.min_n = param.requires_int( p + "min-n" );
  if ( param.has( p + "min-seg-n" ) ) opt.min_seg_n = param.requires_int( p + "min-seg-n" );
  if ( param.has( p + "min-bg-min" ) ) opt.min_bg_min = param.requires_dbl( p + "min-bg-min" );
  if ( param.has( p + "min-seg-bg-min" ) ) opt.min_seg_bg_min = param.requires_dbl( p + "min-seg-bg-min" );

  static bool logged_once = false;
  if ( ! logged_once )
    {
      logged_once = true;
      const std::string pad = "    ";
      const int key_w = 22;
      const int val_w = 11;

      auto show_bool = [&]( const std::string & name , bool value , const std::string & desc )
	{
	  std::ostringstream ss;
	  ss << pad << std::left << std::setw( key_w ) << name << " : "
	     << std::left << std::setw( val_w ) << ( value ? "T" : "F" )
	     << "  " << desc << "\n";
	  logger << ss.str();
	};
      auto show_num = [&]( const std::string & name , const std::string & value , const std::string & desc )
	{
	  std::ostringstream ss;
	  ss << pad << std::left << std::setw( key_w ) << name << " : "
	     << std::left << std::setw( val_w ) << value
	     << "  " << desc << "\n";
	  logger << ss.str();
	};
      auto show_list = [&]( const std::string & name , const std::vector<std::string> & value , const std::string & desc )
	{
	  std::ostringstream joined;
	  for (int i=0;i<value.size();i++) joined << ( i ? "," : "" ) << value[i];
	  std::ostringstream ss;
	  ss << pad << std::left << std::setw( key_w ) << name << " : "
	     << std::left << std::setw( val_w ) << joined.str()
	     << "  " << desc << "\n";
	  logger << ss.str();
	};
      auto show_set = [&]( const std::string & name , const std::set<std::string> & value , const std::string & desc )
	{
	  std::ostringstream joined;
	  for (std::set<std::string>::const_iterator ii = value.begin(); ii != value.end(); ++ii)
	    joined << ( ii == value.begin() ? "" : "," ) << *ii;
	  std::ostringstream ss;
	  ss << pad << std::left << std::setw( key_w ) << name << " : "
	     << std::left << std::setw( val_w ) << joined.str()
	     << "  " << desc << "\n";
	  logger << ss.str();
	};

      logger << "  EVTDYN options";
      if ( ! p.empty() ) logger << " (prefix " << p << ")";
      logger << "\n";
      logger << "    background / scope\n";
      show_bool( "bg-none" , opt.bg_none , "ignore annotation background masks" );
      show_bool( "hypno" , opt.hypno , "use hypnogram for stage-stratified metrics (N2/N3, ASC/DSC, cycles)" );
      show_bool( "verbose" , opt.verbose , "emit extra diagnostics" );
      show_list( "bg" , opt.bg , "background annotations to include" );
      logger << "    event selection / pairwise vars\n";
      show_set( "vars" , opt.vars , "event variables to summarise" );
      show_set( "corr" , opt.corr , "pairwise vars for within-set correlations" );
      show_set( "corr1" , opt.corr1 , "pairwise vars for 1 x 2 correlations: left set" );
      show_set( "corr2" , opt.corr2 , "pairwise vars for 1 x 2 correlations: right set" );
      show_set( "log" , opt.log_vars , "vars to log-transform; if empty, no log transform" );
      logger << "    event structure\n";
      show_num( "anchor" , opt.anchor , "event anchor used for timing summaries" );
      show_num( "cluster" , fmt_num( opt.cluster_sec ) , "cluster window in seconds" );
      show_num( "train-gap" , fmt_num( opt.train_gap_sec ) , "maximum gap within a train" );
      show_num( "train-min" , fmt_num( opt.train_min ) , "minimum events per train" );
      show_num( "refractory" , fmt_num( opt.refr_lwr ) + "," + fmt_num( opt.refr_upr ) , "refractory recurrence window" );
      show_num( "excitatory" , fmt_num( opt.excit_lwr ) + "," + fmt_num( opt.excit_upr ) , "excitatory recurrence window" );
      show_num( "lag-min" , fmt_num( opt.lag_min_sec ) , "lag histogram lower bound (defines window denominator)" );
      show_num( "lag-max" , fmt_num( opt.lag_max_sec ) , "consecutive-pair lag histogram max lag" );
      if ( opt.lag_log > 0 )
	show_num( "lag-log" , fmt_num( opt.lag_log ) , "log-spaced bins [0.5, lag-max]" );
      else
	show_num( "lag-lin" , fmt_num( opt.lag_lin > 0 ? opt.lag_lin : (int)opt.lag_max_sec ) , "linear bins [0, lag-max]" );
      show_bool( "lag-hist" , opt.lag_hist , "output per-bin lag histogram (LAG stratum)" );
      logger << "    transforms\n";
      show_num( "winsor" , fmt_num( opt.winsor_p ) , "winsorization fraction" );
      show_bool( "z" , opt.z , "z-score values before summarizing" );
      show_bool( "rank" , opt.rank , "rank-transform values before summarizing" );
      logger << "    minimum thresholds\n";
      show_num( "min-n" , fmt_num( opt.min_n ) , "min events to compute any metrics" );
      show_num( "min-seg-n" , fmt_num( opt.min_seg_n ) , "min events per subgroup for comparisons" );
      show_num( "min-bg-min" , fmt_num( opt.min_bg_min ) , "min background duration (minutes)" );
      show_num( "min-seg-bg-min" , fmt_num( opt.min_seg_bg_min ) , "min background per subgroup for density comparisons (minutes)" );
    }

  bg = compile_background();
}

void evtdyn_t::add_event( const std::map<std::string,std::string> & faclvl ,
			  const interval_t & interval ,
			  const std::string & ch ,
			  const std::map<std::string,double> & values )
{
  evtdyn_event_t e;
  e.faclvl = faclvl;
  e.interval = interval;
  e.ch = ch;
  e.values = values;
  events.push_back( e );
  if ( seen.insert( faclvl ).second ) order.push_back( faclvl );
}

void evtdyn_t::proc_all()
{
  if ( events.empty() ) return;

  logger << "  running EVTDYN submodule for " << order.size() << " distinct strata\n";

  for (int i=0;i<order.size();i++)
    {
      std::vector<evtdyn_event_t> e;
      for (int j=0;j<events.size();j++)
	if ( events[j].faclvl == order[i] ) e.push_back( events[j] );
      output( order[i] , e );
    }
}

void evtdyn_t::exclude_default_vars( const std::set<std::string> & vars )
{
  opt.default_exclude_vars.insert( vars.begin() , vars.end() );
}

std::vector<interval_t> evtdyn_t::merge( std::vector<interval_t> x )
{
  if ( x.empty() ) return x;
  std::sort( x.begin() , x.end() );
  std::vector<interval_t> y;
  y.push_back( x[0] );
  for (int i=1;i<x.size();i++)
    {
      if ( y.back().stop >= x[i].start )
	{
	  if ( x[i].stop > y.back().stop ) y.back().stop = x[i].stop;
	}
      else y.push_back( x[i] );
    }
  return y;
}

std::vector<interval_t> evtdyn_t::intersect( const std::vector<interval_t> & a ,
					     const std::vector<interval_t> & b )
{
  std::vector<interval_t> y;
  for (int i=0;i<a.size();i++)
    for (int j=0;j<b.size();j++)
      if ( a[i].overlaps( b[j] ) )
	{
	  interval_t z = a[i].intersection_with_overlapping_interval( b[j] );
	  if ( z.stop >= z.start ) y.push_back( z );
	}
  return merge( y );
}

double evtdyn_t::duration_sec( const std::vector<interval_t> & x )
{
  double s = 0;
  for (int i=0;i<x.size();i++) s += x[i].duration_sec();
  return s;
}

std::vector<interval_t> evtdyn_t::epoch_intervals( const bool use_mask ) const
{
  std::vector<interval_t> x;
  edf->timeline.first_epoch();
  while ( 1 )
    {
      const int e = use_mask ? edf->timeline.next_epoch() : edf->timeline.next_epoch_ignoring_mask();
      if ( e == -1 ) break;
      x.push_back( edf->timeline.epoch( e ) );
    }
  return merge( x );
}

std::vector<interval_t> evtdyn_t::annotation_intervals( const std::vector<std::string> & labels ) const
{
  std::vector<interval_t> x;
  for (int i=0;i<labels.size();i++)
    {
      annot_t * a = edf->annotations->find( labels[i] );
      if ( a == NULL ) continue;
      annot_map_t::const_iterator ii = a->interval_events.begin();
      while ( ii != a->interval_events.end() )
	{
	  x.push_back( ii->first.interval );
	  ++ii;
	}
    }
  return merge( x );
}

std::vector<interval_t> evtdyn_t::stage_intervals( const std::set<std::string> & labels ) const
{
  std::vector<interval_t> x;
  if ( ! edf->timeline.epoched() ) return x;
  bool any_stage = false;
  for (std::set<std::string>::const_iterator ss = labels.begin(); ss != labels.end(); ++ss)
    if ( edf->timeline.epoch_annotation( *ss ) ) any_stage = true;
  if ( ! any_stage ) return x;
  edf->timeline.first_epoch();
  int idx = 0;
  while ( 1 )
    {
      const int e = edf->timeline.next_epoch();
      if ( e == -1 ) break;
      bool add = false;
      for (std::set<std::string>::const_iterator ss = labels.begin(); ss != labels.end(); ++ss)
	if ( edf->timeline.epoch_annotation( *ss ) && edf->timeline.epoch_annotation( *ss , idx ) ) add = true;
      if ( add ) x.push_back( edf->timeline.epoch( e ) );
      ++idx;
    }
  return merge( x );
}

std::vector<interval_t> evtdyn_t::sleep_intervals() const
{
  std::set<std::string> sleep;
  sleep.insert( "N1" );
  sleep.insert( "N2" );
  sleep.insert( "N3" );
  sleep.insert( "NREM4" );
  sleep.insert( "R" );
  return stage_intervals( sleep );
}

bool evtdyn_t::sleep_period( interval_t * period , std::vector<interval_t> * sleep ) const
{
  std::vector<interval_t> s = sleep_intervals();
  if ( s.empty() ) return false;
  s = merge( s );
  if ( sleep != NULL ) *sleep = s;
  if ( period != NULL ) *period = interval_t( s.front().start , s.back().stop );
  return true;
}

std::vector<interval_t> evtdyn_t::cycle_intervals( const int cycle ) const
{
  std::vector<interval_t> x;
  if ( ! edf->timeline.epoched() ) return x;
  const std::string c = "_NREMC_" + Helper::int2str( cycle );
  if ( ! edf->timeline.epoch_annotation( c ) ) return x;
  edf->timeline.first_epoch();
  int idx = 0;
  while ( 1 )
    {
      const int e = edf->timeline.next_epoch();
      if ( e == -1 ) break;
      if ( edf->timeline.epoch_annotation( c , idx ) ) x.push_back( edf->timeline.epoch( e ) );
      ++idx;
    }
  return merge( x );
}

evtdyn_bg_t evtdyn_t::compile_background()
{
  evtdyn_bg_t r;

  const std::vector<interval_t> all = epoch_intervals( ! opt.bg_none );

  if ( opt.bg_none )
    r.intervals = all;
  else if ( ! opt.bg.empty() )
    r.intervals = intersect( all , annotation_intervals( opt.bg ) );
  else
    r.intervals = all;

  r.total_sec = duration_sec( r.intervals );
  r.has_hypno = opt.hypno && edf->timeline.epoched()
    && ( edf->timeline.epoch_annotation( "N1" )
	 || edf->timeline.epoch_annotation( "N2" )
	 || edf->timeline.epoch_annotation( "N3" ) );
  r.has_cycles = edf->timeline.epoched()
    && ( edf->timeline.epoch_annotation( "_NREMC_1" )
	 || edf->timeline.epoch_annotation( "_NREMC_2" )
	 || edf->timeline.epoch_annotation( "_NREMC_3" )
	 || edf->timeline.epoch_annotation( "_NREMC_4" )
	 || edf->timeline.epoch_annotation( "_NREMC_5" )
	 || edf->timeline.epoch_annotation( "_NREMC_6" )
	 || edf->timeline.epoch_annotation( "_NREMC_7" )
	 || edf->timeline.epoch_annotation( "_NREMC_8" ) );

  return r;
}

int evtdyn_t::segment( const interval_t & interval , const std::vector<interval_t> & segments ) const
{
  const uint64_t a = anchor( interval );
  for (int i=0;i<segments.size();i++)
    if ( segments[i].contains( a ) ) return i;
  return -1;
}

uint64_t evtdyn_t::anchor( const interval_t & interval ) const
{
  if ( opt.anchor == "START" ) return interval.start;
  if ( opt.anchor == "STOP" ) return interval.stop;
  return interval.mid();
}

double evtdyn_t::elapsed_in_bg( const uint64_t tp , const std::vector<interval_t> & segments ) const
{
  double s = 0;
  for (int i=0;i<segments.size();i++)
    {
      if ( tp >= segments[i].stop )
	s += segments[i].duration_sec();
      else if ( segments[i].contains( tp ) )
	return s + ( tp - segments[i].start ) / (double)globals::tp_1sec;
      else if ( tp < segments[i].start )
	return s;
    }
  return s;
}

double evtdyn_t::elapsed_in_period( const uint64_t tp , const interval_t & period ) const
{
  if ( tp <= period.start ) return 0;
  if ( tp >= period.stop ) return period.duration_sec();
  return ( tp - period.start ) / (double)globals::tp_1sec;
}

void evtdyn_t::output( const std::map<std::string,std::string> & faclvl ,
		       const std::vector<evtdyn_event_t> & e0 )
{
  std::vector<evtdyn_event_t> e;
  for (int i=0;i<e0.size();i++)
    if ( segment( e0[i].interval , bg.intervals ) != -1 ) e.push_back( e0[i] );

  std::sort( e.begin() , e.end() ,
	     [this]( const evtdyn_event_t & a , const evtdyn_event_t & b )
	     { return anchor( a.interval ) < anchor( b.interval ); } );

  writer.level( "EVTDYN" , "DYN" );
  std::map<std::string,std::string>::const_iterator ff = faclvl.begin();
  while ( ff != faclvl.end() )
    {
      writer.level( ff->second , ff->first );
      ++ff;
    }

  const int n = e.size();
  const double mins = bg.total_sec / 60.0;
  writer.value( "N" , n );
  writer.value( "MINS" , mins );
  writer.value( "DENS" , mins > 0 ? n / mins : 0 );

  if ( n < opt.min_n || mins < opt.min_bg_min )
    {
      logger << "  EVTDYN: skipping metrics (n=" << n << " min_n=" << opt.min_n
	     << ", mins=" << mins << " min_bg_min=" << opt.min_bg_min << ")\n";
      for (std::map<std::string,std::string>::const_iterator ff = faclvl.begin(); ff != faclvl.end(); ++ff)
	writer.unlevel( ff->first );
      writer.unlevel( "DYN" );
      return;
    }

  std::vector<double> at;
  std::vector<int> seg;
  for (int i=0;i<n;i++)
    {
      at.push_back( anchor( e[i].interval ) / (double)globals::tp_1sec );
      seg.push_back( segment( e[i].interval , bg.intervals ) );
    }

  std::vector<double> rel;
  for (int i=0;i<n;i++)
    rel.push_back( bg.total_sec > 0 ? elapsed_in_bg( anchor( e[i].interval ) , bg.intervals ) / bg.total_sec : 0 );
  writer.value( "TB" , n ? median( rel ) : 0 );

  interval_t slp_period;
  std::vector<interval_t> slp_intervals;
  if ( n && sleep_period( &slp_period , &slp_intervals ) )
    {
      const double all_sec = slp_period.duration_sec();
      const double sleep_sec = duration_sec( slp_intervals );
      std::vector<double> ta, ts;
      for (int i=0;i<n;i++)
	{
	  const uint64_t a = anchor( e[i].interval );
	  if ( all_sec > 0 )
	    ta.push_back( elapsed_in_period( a , slp_period ) / all_sec );
	  if ( sleep_sec > 0 )
	    ts.push_back( elapsed_in_bg( a , slp_intervals ) / sleep_sec );
	}
      if ( ta.size() ) writer.value( "TA" , median( ta ) );
      if ( ts.size() ) writer.value( "TS" , median( ts ) );
    }

  int early = 0, late = 0;
  for (int i=0;i<rel.size();i++)
    if ( rel[i] < 0.5 ) ++early; else ++late;
  if ( early >= opt.min_seg_n && late >= opt.min_seg_n )
    write_log2_ratio( "EARLY_LATE_LR" , early , late );

  std::vector<double> isi;
  int refr = 0, excit = 0;
  for (int i=1;i<n;i++)
    {
      if ( seg[i] == -1 || seg[i] != seg[i-1] ) continue;
      const double lag = at[i] - at[i-1];
      if ( lag < 0 ) continue;
      isi.push_back( lag );
      if ( lag >= opt.refr_lwr && lag <= opt.refr_upr ) ++refr;
      if ( lag >= opt.excit_lwr && lag <= opt.excit_upr ) ++excit;
	}
  if ( ! isi.empty() ) writer.value( "ISI_MD" , median( isi ) );

  const double rate = bg.total_sec > 0 ? n / bg.total_sec : 0;
  const double refr_exp = rate > 0 ? std::exp( -rate * opt.refr_lwr ) - std::exp( -rate * opt.refr_upr ) : 0;
  const double excit_exp = rate > 0 ? std::exp( -rate * opt.excit_lwr ) - std::exp( -rate * opt.excit_upr ) : 0;
  if ( n > 1 )
    {
      write_log2_ratio( "REFR_OBS_EXP" , refr / (double)( n - 1 ) , refr_exp );
      write_log2_ratio( "EXCIT_OBS_EXP" , excit / (double)( n - 1 ) , excit_exp );
    }

  const bool lag_log_mode = opt.lag_log > 0;
  const int nb = lag_log_mode ? opt.lag_log
                              : ( opt.lag_lin > 0 ? opt.lag_lin : (int)( opt.lag_max_sec - opt.lag_min_sec ) );

  std::vector<double> lag_edges( nb + 1 ), lag_centers( nb );
  if ( lag_log_mode )
    {
      const double lmin = std::log( opt.lag_min_sec );
      const double lmax = std::log( opt.lag_max_sec );
      for (int k=0;k<=nb;k++)
	lag_edges[k] = std::exp( lmin + (double)k / nb * ( lmax - lmin ) );
      for (int k=0;k<nb;k++)
	lag_centers[k] = std::sqrt( lag_edges[k] * lag_edges[k+1] );
    }
  else
    {
      for (int k=0;k<=nb;k++)
	lag_edges[k] = opt.lag_min_sec + (double)k / nb * ( opt.lag_max_sec - opt.lag_min_sec );
      for (int k=0;k<nb;k++)
	lag_centers[k] = ( lag_edges[k] + lag_edges[k+1] ) / 2.0;
    }

  std::vector<int> lh( nb , 0 );
  int n_lag_pairs = 0;
  for (int i=1;i<n;i++)
    {
      if ( seg[i] == -1 || seg[i] != seg[i-1] ) continue;
      const double lag = at[i] - at[i-1];
      if ( lag < lag_edges[0] || lag > lag_edges[nb] ) continue;
      ++n_lag_pairs;
      const int b = (int)( std::upper_bound( lag_edges.begin() , lag_edges.end() , lag ) - lag_edges.begin() ) - 1;
      if ( b >= 0 && b < nb ) ++lh[b];
    }
  // per-bin expected: uniform-in-lag-time null (proportional to bin width in seconds)
  const double window_width = lag_edges[nb] - lag_edges[0];
  std::vector<double> exp_bins( nb , 0 );
  for (int k=0;k<nb;k++)
    exp_bins[k] = window_width > 0 ? (double)n_lag_pairs * ( lag_edges[k+1] - lag_edges[k] ) / window_width : 0;

  int pk = -1;
  for (int i=0;i<nb;i++)
    if ( pk == -1 || lh[i] > lh[pk] ) pk = i;
  if ( pk >= 0 && lh[pk] > 0 && nb > 0 )
    {
      writer.value( "LAG_PEAK_T" , lag_centers[pk] );
      write_log2_ratio( "LAG_PEAK_H" , lh[pk] , exp_bins[pk] );
    }
  if ( opt.lag_hist && nb > 0 )
    {
      for (int b=0;b<nb;b++)
	{
	  writer.level( lag_centers[b] , "LAG" );
	  writer.value( "LAG_N" , lh[b] );
	  if ( exp_bins[b] > 0 )
	    writer.value( "LAG_OE" , log2_ratio( lh[b] + 0.5 , exp_bins[b] ) );
	  writer.unlevel( "LAG" );
	}
    }

  std::vector<int> clustered( n , 0 ), train_id( n , -1 );
  for (int i=0;i<n;i++)
    for (int j=i+1;j<n;j++)
      {
	if ( seg[i] == -1 || seg[i] != seg[j] ) continue;
	const double lag = at[j] - at[i];
	if ( lag > opt.cluster_sec ) break;
	clustered[i] = clustered[j] = 1;
      }

  int cln = std::accumulate( clustered.begin() , clustered.end() , 0 );
  writer.value( "CLST_FRAC" , n ? cln / (double)n : 0 );

  std::vector<int> train_len;
  int i = 0;
  while ( i < n )
    {
      int j = i + 1;
      while ( j < n && seg[j] != -1 && seg[j] == seg[j-1] && at[j] - at[j-1] <= opt.train_gap_sec ) ++j;
      if ( j - i >= opt.train_min )
	{
	  train_len.push_back( j - i );
	  for (int k=i;k<j;k++) train_id[k] = train_len.size() - 1;
	}
      i = j;
    }
  int trn = 0;
  for (int i=0;i<train_id.size();i++) if ( train_id[i] != -1 ) ++trn;
  writer.value( "TRAIN_FRAC" , n ? trn / (double)n : 0 );
  if ( ! train_len.empty() )
    writer.value( "TRAIN_LEN_MN" , mean( std::vector<double>( train_len.begin() , train_len.end() ) ) );
  writer.value( "TRAIN_DENS" , mins > 0 ? train_len.size() / mins : 0 );

  if ( bg.has_hypno )
    {
      std::set<std::string> n2; n2.insert( "N2" );
      std::set<std::string> n3; n3.insert( "N3" );
      std::vector<interval_t> b2 = intersect( bg.intervals , stage_intervals( n2 ) );
      std::vector<interval_t> b3 = intersect( bg.intervals , stage_intervals( n3 ) );
      int n2c = 0, n3c = 0;
      for (int i=0;i<n;i++)
	{
	  if ( segment( e[i].interval , b2 ) != -1 ) ++n2c;
	  if ( segment( e[i].interval , b3 ) != -1 ) ++n3c;
	}
      const double dur2 = duration_sec( b2 ) / 60.0;
      const double dur3 = duration_sec( b3 ) / 60.0;
      const double d2 = dur2 > 0 ? n2c / dur2 : 0;
      const double d3 = dur3 > 0 ? n3c / dur3 : 0;
      if ( n2c >= opt.min_seg_n && n3c >= opt.min_seg_n
	   && dur2 >= opt.min_seg_bg_min && dur3 >= opt.min_seg_bg_min )
	write_log2_ratio( "N2_N3_LR" , d2 , d3 );
    }

  std::vector<interval_t> b_n2asc, b_n2dsc;
  bool has_n2_asc_dsc = false;
  if ( edf->timeline.epoch_annotation( "_N2_ASC" )
       && edf->timeline.epoch_annotation( "_N2_DSC" ) )
    {
      std::set<std::string> n2asc; n2asc.insert( "_N2_ASC" );
      std::set<std::string> n2dsc; n2dsc.insert( "_N2_DSC" );
      b_n2asc = intersect( bg.intervals , stage_intervals( n2asc ) );
      b_n2dsc = intersect( bg.intervals , stage_intervals( n2dsc ) );
      has_n2_asc_dsc = ! b_n2asc.empty() && ! b_n2dsc.empty();
      int na = 0, nd = 0;
      for (int i=0;i<n;i++)
	{
	  if ( segment( e[i].interval , b_n2asc ) != -1 ) ++na;
	  if ( segment( e[i].interval , b_n2dsc ) != -1 ) ++nd;
	}
      const double dura = duration_sec( b_n2asc ) / 60.0;
      const double durd = duration_sec( b_n2dsc ) / 60.0;
      const double da = dura > 0 ? na / dura : 0;
      const double dd = durd > 0 ? nd / durd : 0;
      if ( na >= opt.min_seg_n && nd >= opt.min_seg_n
	   && dura >= opt.min_seg_bg_min && durd >= opt.min_seg_bg_min )
	{
	  writer.value( "N2_ASC_DSC_DIFF" , da - dd );
	  write_log2_ratio( "N2_ASC_DSC_LR" , da , dd );
	}
    }

  // within-cycle relative position for each event (NaN if not in any cycle background)
  std::vector<double> wcyc_rel( n , std::numeric_limits<double>::quiet_NaN() );
  double wcyc_edur = 0 , wcyc_ldur = 0; // total early/late bg-half seconds across all cycles
  if ( bg.has_cycles )
    {
      std::vector<double> cx, cy;
      for (int c=1;c<=4;c++)
	{
	  std::vector<interval_t> bc = intersect( bg.intervals , cycle_intervals( c ) );
	  if ( bc.empty() ) continue;
	  const double dur = duration_sec( bc );
	  if ( dur <= 0 ) continue;
	  wcyc_edur += dur / 2.0;
	  wcyc_ldur += dur / 2.0;
	  int nc = 0;
	  for (int i=0;i<n;i++)
	    if ( segment( e[i].interval , bc ) != -1 )
	      {
		++nc;
		wcyc_rel[i] = elapsed_in_bg( anchor( e[i].interval ) , bc ) / dur;
	      }
	  cx.push_back( c );
	  cy.push_back( nc / ( dur / 60.0 ) );
	}
      writer.value( "CYCLE_SLOPE" , cx.size() > 1 ? slope( cx , cy ) : 0 );

      // within-cycle timing and density split
      std::vector<double> wcyc_valid;
      int wcyc_en = 0 , wcyc_ln = 0;
      for (int i=0;i<n;i++)
	{
	  if ( ! is_finite( wcyc_rel[i] ) ) continue;
	  wcyc_valid.push_back( wcyc_rel[i] );
	  if ( wcyc_rel[i] < 0.5 ) ++wcyc_en; else ++wcyc_ln;
	}
      if ( ! wcyc_valid.empty() )
	writer.value( "WCYC_TB" , median( wcyc_valid ) );
      if ( wcyc_en >= opt.min_seg_n && wcyc_ln >= opt.min_seg_n
	   && ( wcyc_edur / 60.0 ) >= opt.min_seg_bg_min
	   && ( wcyc_ldur / 60.0 ) >= opt.min_seg_bg_min )
	{
	  const double de = wcyc_en / ( wcyc_edur / 60.0 );
	  const double dl = wcyc_ln / ( wcyc_ldur / 60.0 );
	  write_log2_ratio( "WCYC_EARLY_LATE_LR" , de , dl );
	}
    }

  std::set<std::string> vars = opt.vars;
  if ( vars.empty() )
    for (int i=0;i<e.size();i++)
      for (std::map<std::string,double>::const_iterator vv = e[i].values.begin(); vv != e[i].values.end(); ++vv)
	if ( opt.default_exclude_vars.find( vv->first ) == opt.default_exclude_vars.end() )
	  vars.insert( vv->first );

  std::map<std::string,corr_data_t> corr_cache;

  for (std::set<std::string>::const_iterator vv = vars.begin(); vv != vars.end(); ++vv)
    {
      std::vector<double> x, t, relv;
      std::vector<int> rows;
      int raw_n = 0;
      double raw_sum = 0;
      for (int i=0;i<n;i++)
	{
	  std::map<std::string,double>::const_iterator ii = e[i].values.find( *vv );
	  if ( ii == e[i].values.end() || ! is_finite( ii->second ) ) continue;
	  ++raw_n;
	  raw_sum += ii->second;
	  x.push_back( ii->second );
	  t.push_back( rel[i] );
	  relv.push_back( rel[i] );
	  rows.push_back( i );
	}
      if ( x.empty() ) continue;
      const double pre_sum = std::accumulate( x.begin() , x.end() , 0.0 );
      const bool do_log = ! opt.log_vars.empty() && opt.log_vars.find( *vv ) != opt.log_vars.end();
      const std::vector<int> keep = clean_keep_idx( x , opt , do_log );
      x = apply_keep( x , keep );
      t = apply_keep( t , keep );
      relv = apply_keep( relv , keep );
      rows = apply_keep( rows , keep );
      x = clean_values( x , opt , do_log );
      const double post_sum = std::accumulate( x.begin() , x.end() , 0.0 );
      if ( opt.verbose )
	{
	  logger << "  EVTDYN calc " << *vv
		 << " n_raw=" << raw_n
		 << " sum_raw=" << raw_sum
		 << " n_kept=" << x.size()
		 << " sum_kept=" << post_sum
		 << " mean=" << ( x.empty() ? 0 : mean( x ) )
		 << "\n";
	  if ( *vv == "SO_OVERLAP" || *vv == "COUPL_OVERLAP" )
	    {
	      int ones = 0;
	      for (int j=0;j<x.size();j++) if ( x[j] > 0.5 ) ++ones;
	      logger << "    binary detail: ones=" << ones
		     << " zeros=" << (int)x.size() - ones
		     << " pre_clean_sum=" << pre_sum
		     << " post_clean_sum=" << post_sum
		     << "\n";
	    }
	}
      std::vector<double> xc, xs, xa, xd;
      for (int j=0;j<x.size();j++)
	{
	  const int r = rows[j];
	  if ( clustered[r] ) xc.push_back( x[j] ); else xs.push_back( x[j] );
	  if ( ! b_n2asc.empty() && segment( e[r].interval , b_n2asc ) != -1 ) xa.push_back( x[j] );
	  if ( ! b_n2dsc.empty() && segment( e[r].interval , b_n2dsc ) != -1 ) xd.push_back( x[j] );
	}
      writer.level( *vv , globals::var_strat );
      writer.value( "MEAN" , x.empty() ? 0 : mean( x ) );
      writer.value( "MEDIAN" , x.empty() ? 0 : median( x ) );
      writer.value( "SD" , x.empty() ? 0 : sd( x ) );
      if ( x.size() == t.size() && x.size() > 1 )
	{
	  writer.value( "TIME_CORR" , corr( t , x ) );
	  writer.value( "TIME_BETA" , slope( t , x ) );
	}
      std::vector<double> early_x, late_x;
      for (int j=0;j<x.size();j++)
	if ( relv[j] < 0.5 ) early_x.push_back( x[j] ); else late_x.push_back( x[j] );
      if ( (int)early_x.size() >= opt.min_seg_n && (int)late_x.size() >= opt.min_seg_n )
	writer.value( "EARLY_LATE_DIFF" , mean( early_x ) - mean( late_x ) );
      if ( bg.has_hypno )
	{
	  std::set<std::string> n2; n2.insert( "N2" );
	  std::set<std::string> n3; n3.insert( "N3" );
	  std::vector<interval_t> b2 = intersect( bg.intervals , stage_intervals( n2 ) );
	  std::vector<interval_t> b3 = intersect( bg.intervals , stage_intervals( n3 ) );
	  std::vector<double> x2, x3;
	  for (int j=0;j<x.size();j++)
	    {
	      const int r = rows[j];
	      if ( segment( e[r].interval , b2 ) != -1 ) x2.push_back( x[j] );
	      if ( segment( e[r].interval , b3 ) != -1 ) x3.push_back( x[j] );
	    }
	  if ( (int)x2.size() >= opt.min_seg_n && (int)x3.size() >= opt.min_seg_n )
	    writer.value( "N2_N3_DIFF" , mean( x2 ) - mean( x3 ) );
	}
      if ( bg.has_cycles )
	{
	  std::vector<double> cx, cy;
	  for (int c=1;c<=4;c++)
	    {
	      std::vector<interval_t> bc = intersect( bg.intervals , cycle_intervals( c ) );
	      if ( bc.empty() ) continue;
	      std::vector<double> cv;
	      for (int j=0;j<x.size();j++)
		if ( segment( e[ rows[j] ].interval , bc ) != -1 ) cv.push_back( x[j] );
	      if ( cv.empty() ) continue;
	      cx.push_back( c );
	      cy.push_back( mean( cv ) );
	    }
	  if ( cx.size() > 1 )
	    writer.value( "CYCLE_SLOPE" , slope( cx , cy ) );

	  // within-cycle property dynamics (pool across cycles)
	  std::vector<double> wt, wx, wx_e, wx_l;
	  for (int j=0;j<x.size();j++)
	    {
	      const double wr = wcyc_rel[ rows[j] ];
	      if ( ! is_finite( wr ) ) continue;
	      wt.push_back( wr );
	      wx.push_back( x[j] );
	      if ( wr < 0.5 ) wx_e.push_back( x[j] ); else wx_l.push_back( x[j] );
	    }
	  if ( wt.size() > 1 )
	    {
	      writer.value( "WCYC_CORR" , corr( wt , wx ) );
	      writer.value( "WCYC_BETA" , slope( wt , wx ) );
	    }
	  if ( (int)wx_e.size() >= opt.min_seg_n && (int)wx_l.size() >= opt.min_seg_n )
	    writer.value( "WCYC_EARLY_LATE_DIFF" , mean( wx_e ) - mean( wx_l ) );
	}
      if ( (int)xc.size() >= opt.min_seg_n && (int)xs.size() >= opt.min_seg_n )
	writer.value( "CLST_SOL_DIFF" , mean( xc ) - mean( xs ) );
      if ( has_n2_asc_dsc )
	{
	  if ( (int)xa.size() >= opt.min_seg_n && (int)xd.size() >= opt.min_seg_n )
	    writer.value( "N2_ASC_DSC_DIFF" , mean( xa ) - mean( xd ) );
	}
      if ( ! opt.corr.empty() || ! opt.corr1.empty() || ! opt.corr2.empty() )
	{
	  corr_data_t & vc = corr_cache[ *vv ];
	  vc.rows = rows;
	  vc.values = x;
	}
      writer.unlevel( globals::var_strat );
    }

  if ( opt.corr.size() > 1 )
    {
      std::vector<std::string> pv( opt.corr.begin() , opt.corr.end() );
      for (int i=0;i<pv.size();i++)
	for (int j=i+1;j<pv.size();j++)
	  {
	    std::map<std::string,corr_data_t>::const_iterator ia = corr_cache.find( pv[i] );
	    std::map<std::string,corr_data_t>::const_iterator ib = corr_cache.find( pv[j] );
	    if ( ia == corr_cache.end() || ib == corr_cache.end() ) continue;

	    std::map<int,double> ma, mb;
	    for (int k=0;k<ia->second.rows.size();k++)
	      ma[ ia->second.rows[k] ] = ia->second.values[k];
	    for (int k=0;k<ib->second.rows.size();k++)
	      mb[ ib->second.rows[k] ] = ib->second.values[k];

	    std::vector<double> xa, xb;
	    for (std::map<int,double>::const_iterator mm = ma.begin(); mm != ma.end(); ++mm)
	      {
		std::map<int,double>::const_iterator nn = mb.find( mm->first );
		if ( nn == mb.end() ) continue;
		xa.push_back( mm->second );
		xb.push_back( nn->second );
	      }
	    if ( xa.size() < 2 ) continue;
	    writer.level( pv[i] , "VAR1" );
	    writer.level( pv[j] , "VAR2" );
	    writer.value( "CORR" , corr( xa , xb ) );
	    writer.value( "N" , (int)xa.size() );
	    writer.unlevel( "VAR2" );
	    writer.unlevel( "VAR1" );
	  }
    }

  if ( ! opt.corr1.empty() && ! opt.corr2.empty() )
    {
      std::vector<std::string> p1( opt.corr1.begin() , opt.corr1.end() );
      std::vector<std::string> p2( opt.corr2.begin() , opt.corr2.end() );
      for (int i=0;i<p1.size();i++)
	for (int j=0;j<p2.size();j++)
	  {
	    std::map<std::string,corr_data_t>::const_iterator ia = corr_cache.find( p1[i] );
	    std::map<std::string,corr_data_t>::const_iterator ib = corr_cache.find( p2[j] );
	    if ( ia == corr_cache.end() || ib == corr_cache.end() ) continue;

	    std::map<int,double> ma, mb;
	    for (int k=0;k<ia->second.rows.size();k++)
	      ma[ ia->second.rows[k] ] = ia->second.values[k];
	    for (int k=0;k<ib->second.rows.size();k++)
	      mb[ ib->second.rows[k] ] = ib->second.values[k];

	    std::vector<double> xa, xb;
	    for (std::map<int,double>::const_iterator mm = ma.begin(); mm != ma.end(); ++mm)
	      {
		std::map<int,double>::const_iterator nn = mb.find( mm->first );
		if ( nn == mb.end() ) continue;
		xa.push_back( mm->second );
		xb.push_back( nn->second );
	      }
	    if ( xa.size() < 2 ) continue;
	    writer.level( p1[i] , "VAR1" );
	    writer.level( p2[j] , "VAR2" );
	    writer.value( "CORR" , corr( xa , xb ) );
	    writer.value( "N" , (int)xa.size() );
	    writer.unlevel( "VAR2" );
	    writer.unlevel( "VAR1" );
	  }
    }

  for (std::map<std::string,std::string>::const_iterator ff = faclvl.begin(); ff != faclvl.end(); ++ff)
    writer.unlevel( ff->first );
  writer.unlevel( "DYN" );
}

void evtdyn_t::from_annots( edf_t & edf , param_t & param )
{
  evtdyn_t ev;
  ev.init( edf , param );
  if ( ! param.has( "annot" ) ) Helper::halt( "requires annot parameter" );
  std::vector<std::string> anames = param.strvector( "annot" );
  std::set<std::string> vars = param.has( "vars" ) ? param.strset( "vars" ) : std::set<std::string>();
  for (int a=0;a<anames.size();a++)
    {
      annot_t * annot = edf.annotations->find( anames[a] );
      if ( annot == NULL ) continue;
      annot_map_t::const_iterator ii = annot->interval_events.begin();
      while ( ii != annot->interval_events.end() )
	{
	  std::map<std::string,std::string> faclvl;
	  faclvl[ globals::annot_strat ] = anames[a];
	  if ( ii->first.ch_str != "." ) faclvl[ globals::signal_strat ] = ii->first.ch_str;
	  std::map<std::string,double> values;
	  if ( ii->second != NULL )
	    {
	      std::set<std::string> use = vars;
	      if ( use.empty() )
		for (std::map<std::string,avar_t*>::const_iterator dd = ii->second->data.begin();
		     dd != ii->second->data.end(); ++dd) use.insert( dd->first );
	      for (std::set<std::string>::const_iterator vv = use.begin(); vv != use.end(); ++vv)
		{
		  avar_t * av = ii->second->find( *vv );
		  if ( av != NULL && ! av->is_vector() && ! av->is_missing() )
		    values[ *vv ] = av->double_value();
		}
	    }
	  ev.add_event( faclvl , ii->first.interval , ii->first.ch_str , values );
	  ++ii;
	}
    }
  ev.proc_all();
}
