//    --------------------------------------------------------------------
//
//    Event-level dynamics summaries
//
//    --------------------------------------------------------------------

#ifndef __EVTDYN_H__
#define __EVTDYN_H__

#include "intervals/intervals.h"

#include <map>
#include <set>
#include <string>
#include <vector>

struct edf_t;
struct param_t;

struct evtdyn_event_t
{
  std::map<std::string,std::string> faclvl;
  interval_t interval;
  std::string ch;
  std::map<std::string,double> values;
};

struct evtdyn_opts_t
{
  std::string prefix;
  bool bg_none;
  bool hypno;
  bool verbose;
  std::vector<std::string> bg;
  std::set<std::string> vars;
  std::set<std::string> corr;
  std::set<std::string> corr1;
  std::set<std::string> corr2;
  std::set<std::string> log_vars;
  std::set<std::string> default_exclude_vars;
  std::string anchor;
  double cluster_sec;
  double train_gap_sec;
  int train_min;
  double lag_short_lwr;
  double lag_short_upr;
  double refr_lwr;
  double refr_upr;
  double excit_lwr;
  double excit_upr;
  double ac_max_sec;
  double ac_bin_sec;
  double winsor_p;
  bool z;
  bool rank;

  evtdyn_opts_t();
};

struct evtdyn_bg_t
{
  std::vector<interval_t> intervals;
  double total_sec;
  bool has_hypno;
  bool has_cycles;
  evtdyn_bg_t() : total_sec(0), has_hypno(false), has_cycles(false) { }
};

struct evtdyn_t
{
  void init( edf_t & edf , const param_t & param , const std::string & prefix = "" );

  void add_event( const std::map<std::string,std::string> & faclvl ,
		  const interval_t & interval ,
		  const std::string & ch ,
		  const std::map<std::string,double> & values );

  void proc_all();

  void exclude_default_vars( const std::set<std::string> & vars );

  static void from_annots( edf_t & edf , param_t & param );

private:
  edf_t * edf;
  evtdyn_opts_t opt;
  evtdyn_bg_t bg;
  std::vector<evtdyn_event_t> events;
  std::vector<std::map<std::string,std::string> > order;
  std::set<std::map<std::string,std::string> > seen;

  evtdyn_bg_t compile_background();
  std::vector<interval_t> epoch_intervals( const bool use_mask ) const;
  std::vector<interval_t> annotation_intervals( const std::vector<std::string> & labels ) const;
  std::vector<interval_t> stage_intervals( const std::set<std::string> & labels ) const;
  std::vector<interval_t> sleep_intervals() const;
  bool sleep_period( interval_t * period , std::vector<interval_t> * sleep ) const;
  std::vector<interval_t> cycle_intervals( const int cycle ) const;
  static std::vector<interval_t> intersect( const std::vector<interval_t> & a ,
					    const std::vector<interval_t> & b );
  static std::vector<interval_t> merge( std::vector<interval_t> x );
  static double duration_sec( const std::vector<interval_t> & x );
  int segment( const interval_t & interval , const std::vector<interval_t> & segments ) const;
  uint64_t anchor( const interval_t & interval ) const;
  double elapsed_in_bg( const uint64_t tp , const std::vector<interval_t> & segments ) const;
  double elapsed_in_period( const uint64_t tp , const interval_t & period ) const;
  void output( const std::map<std::string,std::string> & faclvl ,
	       const std::vector<evtdyn_event_t> & e );
};

#endif
