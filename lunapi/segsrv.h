\
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

#ifndef __SEGSRV_H__
#define __SEGSRV_H__

#include "luna.h"
#include "stats/Eigen/Dense"

struct lunapi_inst_t;
typedef std::shared_ptr<lunapi_inst_t> lunapi_inst_ptr;

struct evt_t {
  evt_t( const interval_t & interval , const std::string & name )
    : interval(interval) , name(name) { }
  
  interval_t interval;
  std::string name;

  bool operator<( const evt_t & rhs ) const
  {
    if ( interval < rhs.interval ) return true;
    if ( interval > rhs.interval ) return false;
    return name < rhs.name ; 
  }
};

struct fevt_t {
fevt_t( const double start ,
	const double stop ,
	const std::string & name )
  : start(start), stop(stop), name(name) { }
  
  double start, stop;
  std::string name;
  
  bool operator<( const fevt_t & rhs ) const
  {
    if ( start < rhs.start ) return true;
    if ( start > rhs.start ) return false;
    if ( stop < rhs.stop ) return true;
    if ( stop > rhs.stop ) return false;
    return name < rhs.name ; 
  }
};
  
struct segsrv_t {
  
public:

    static Eigen::VectorXf decimate( const Eigen::VectorXf & x0 , const int sr, const int q );
  
  // set up
  segsrv_t( lunapi_inst_ptr ); 

  // initialize (or reset) to current internal state
  int populate( const std::vector<std::string> & chs ,
		const std::vector<std::string> & anns );

  // set physical limits (mins) for a given channel
  void fix_physical_scale( const std::string & ch , const double lwr, const double upr );
  void free_physical_scale( const std::string & ch );
  
  // get overall time-scale
  std::vector<std::pair<double,double> > get_time_scale() const;

  // set window given non-gapped (i.e. elapsed clock time)
  bool set_window( double a , double b );

  // set window given gapped elapsed time (i.e. gapless display)
  bool set_window_ungapped( double a , double b );
  
  // get signals
  Eigen::VectorXf get_signal( const std::string & ch ) const;
  
  // set scaling params
  void set_scaling( const int nchs , const int nanns ,
		    const double yscale , const double ygroup ,
		    const double yheader , const double yfooter ,
		    const double scaling_fixed_annot );


  
  // get scaled signals (0-1 give other annots  
  bool get_yscale_signal( const int n1 , double * lwr, double * upr ) const; // --> can make private
  Eigen::VectorXf get_scaled_signal( const std::string & ch , const int n1 );

  // and times
  Eigen::VectorXf get_timetrack( const std::string & ch ) const;

  // summary stats (mean, min, max, SD)
  Eigen::MatrixXf get_summary_stats( const std::string & ch );
  Eigen::VectorXf get_summary_timetrack( const std::string & ch ) const;

  
  // also, a mask of gapped regions in a given time-window   
  // e.g. if 0 .. 300 is total window and gaps between 30-60 and 270 - 350
  //  --> (30,60) , (270,300) 
  std::set<std::pair<double,double> > get_gaps() const;

  // for internal use when making clock/disp emap
  bool has_gaps( uint64_t atp, uint64_t btp , uint64_t * t = NULL ) const;
  
  // set epoch size (default = 30, at least 4s)
  void set_epoch_size( const double d ) { epoch_sec = d > 4 ? d : 4; } 
  double get_epoch_size() const { return epoch_sec; }
  
  // request bands
  void calc_bands( const std::vector<std::string> & chs );
  void calc_hjorths( const std::vector<std::string> & chs );

  // get summary info back
  int nepochs() const { return epoch_num; }
  int nepochs_clock() const { return clock_epoch_num; }
  
  Eigen::MatrixXf get_bands( const std::string & ch );
  Eigen::MatrixXf get_hjorths( const std::string & ch );
  //Eigen::VectorXf get_epoch_timetrack( ) const { return epoch_sec_starts_clock; } 
  std::map<int,int> clk2sig_emap;
  
  // is a valid window?
  bool is_window_valid() const { return valid_window; }
  double get_window_left() const { return awin; }
  double get_window_right() const { return bwin; }

  /* double get_ewindow_left() const { return aewin; } */
  /* double get_ewindow_right() const { return bewin; } */

  std::string get_window_left_hms() const; 
  std::string get_window_right_hms() const; 
  std::string get_hms( const double s ) const;
  std::map<double,std::string> get_clock_ticks(const int n) const;
  
  clocktime_t edf_start;
  
  //double get_ungapped_total_sec() const { return cumul_sec; }
  double get_total_sec() const; 
  double get_total_sec_original() const; // not impacted by mask - for hypno viz etc

  // initial decimation (on input when populating) if high SR
  int get_input_throttle() const { return max_samples_in; }
  void input_throttle( const int m ) { max_samples_in = m < 0 ? 0 : m; } 

  // on-the-fly decimation on output if high # points (i.e. high SR and/or 
  int get_throttle() const { return max_samples_out; }
  void throttle( const int m ) { max_samples_out = m < 0 ? 0 : m; } 

  // display summary statistics (mean, range, SD) if > minutes
  void summary_threshold_mins( const double s ) { summary_threshold_secs = 60 * s; } 
  
  std::pair<double,double> get_window_phys_range( const std::string & ch ) const
  {
    std::map<std::string,std::pair<double,double> >::const_iterator ss = window_phys_range.find( ch );
    if ( ss == window_phys_range.end() ) return std::pair<double,double>(0,0);
    return ss->second;
  }

  double get_ylabel( const int idx ) const {
    if ( idx < 0 || idx > scaling_upr.size() ) return -1;
    return scaling_upr[idx];
  }
  
  bool serve_raw_signals() const { return bwin - awin > summary_threshold_secs ; } 

  
private:
  
  void init();
  
  // person
  lunapi_inst_ptr p;

  // segments, gaps (clock time)
  std::set<interval_t> segments, gaps;

  // gaps, collapsed/unhapped time (secs)
  //std::set<double> gaps_ungapped;
  
  // to look up indices given two timepoints
  // sr: time->idx
  std::map<int,std::map<double,int> > tidx;

  // current time windows (clock seconds elapsed from EDF starts)
  double awin, bwin;
    
  // track whether a valid window is currently set
  bool valid_window;
  
  // maximum (last time-point in seconds)
  double smax;

  // for throttling  
  int max_samples_in; // decimate at input (populate/add_channel())
  int max_samples_out; // decimate when outputting data (get_signal())
  double summary_threshold_secs; // give summaries at output, instead of raw signals
  
  // current window in idx points (for a given SR)
  std::map<int,int> aidx, bidx; 
  
  // add actual data
  bool add_channel( const std::string & );

  // given two times and a sample rate, get indices
  bool get_tidx( double a, double b , int sr , int * aidx, int *bidx ) const;
  
  // signal data
  std::map<std::string,int> srmap;
  std::map<std::string,Eigen::VectorXf> sigmap;
  std::map<int,Eigen::VectorXf> tmap;

  // post input-decimation, track new implied SR
  std::map<int,double> decimated_srmap;
  
  // scaling
  std::vector<double> scaling_lwr, scaling_upr;
  double scaling_ygroup, scaling_yscale;
  int scaling_nchs, scaling_nanns;
  double scaling_yheader, scaling_yfooter;
  double scaling_fixed_annot;

  // physical scaling (for scaled_signal)
  std::map<std::string, std::pair<double,double> > phys_ranges;

  // store the min/max per signal after a get_scaled_signal()
  std::map<std::string,std::pair<double,double> > window_phys_range;
 
  // cumulative seconds (length of data, w/out gaps)
  //  double cumul_sec;
  
  // epoch size and others
  double epoch_sec; 
  int clock_epoch_num; // simple clk-based viz epochs
  int epoch_num; // defined full epochs w/ signals
  std::vector<double> epoch_sec_starts;
  //  Eigen::VectorXf epoch_sec_starts_clock; // --> get_epoch_timetrack()
  
  // summaries
  void do_summaries( const std::string & ch ,
		     const int sr , const std::vector<double> * data ,
		     const bool do_band, const bool do_hjorth );

  std::map<std::string,Eigen::MatrixXf> bands; // ch -> bands x epochs
  std::map<std::string,Eigen::MatrixXf> hjorth; // ch -> h1/h2/h3 x epochs


  //
  // annotations
  //

public:

  // add annots
  bool add_annot( const std::string & );

  // compile a set of selected events for the current window
  void compile_evts( const std::vector<std::string> & anns );

  // given a compilation (subset of all evts), get evts for a particular class
  std::vector<float> get_evnts_xaxes( const std::string & ann ) const;

  // given a compilation (subset of all evts), get y-axis stacking for a particular class
  std::vector<float> get_evnts_yaxes( const std::string & ann ) const;
  
  // get all events from the current window (not used so much)
  std::map<std::string,std::vector<std::pair<double,double> > > fetch_evts() const;
  
  // for selection window
  std::vector<std::string> fetch_all_evts( const std::vector<std::string> & ) const;
  
private:

  std::set<evt_t> evts;

  // generate by compile_evts(); fetched by get_evnts_xaxes() and get_evnts_yaxes()
  std::map<std::string,std::vector<float> > compiled_annots_times;
  std::map<std::string,std::vector<float> > compiled_annots_stacks;

  
};


#endif
