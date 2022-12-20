
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

#ifndef __HYPNO_H__
#define __HYPNO_H__

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <set>
#include <map>

#include "helper/helper.h"
#include "helper/logger.h"

#include "defs/defs.h"
#include "annot/annot.h"

struct edf_t;

struct timeline_t;

struct interval_t;

extern logger_t logger;

//
// General class to represent a hypnogram and perform calculations on it
// linked to a timeline_t's epochs
//

void dummy_hypno();

struct bout_t
{
  bout_t( const int e1, const int e2, const sleep_stage_t & stg )
    : start(e1), stop(e2), ss(stg) 
  {
    // just track W, S, NR and ?
    // NR --> N2 by default
    if ( ss == NREM1 || ss == NREM3 || ss == NREM4 || ss == NREM2 )
      ss = NREM2;
    else if ( ! ( ss == WAKE || ss == REM ) )
      ss = UNKNOWN;           
  }

  int start, stop;
  sleep_stage_t ss;

  bool operator<( const bout_t & rhs ) const
  {
    if ( start < rhs.start ) return true;
    if ( rhs.start < start ) return false;
    if ( stop < rhs.stop ) return true;
    if ( rhs.stop < stop ) return false;
    return ss < rhs.ss;    
  }
};

struct hypnogram_t
{
  
  hypnogram_t( bool collapse = true ) 
  { 
    collapse_nrem34 = collapse;
  } 
  
  bool collapse_nrem34;
  
  bool construct( timeline_t * t , param_t & param , const bool verbose , const std::vector<std::string> & s );
  bool construct( timeline_t * t , param_t & param , const bool verbose , const std::string sslabel = "SleepStage" );
  void edit(  timeline_t * t , param_t & param  );

  void calc_stats( const bool verbose ); // verbose == STAGES vs HYPNO
  void output( const bool verbose , const bool epoch_lvl_output ,
	       const std::string & eannot = "" ,
	       const std::string & annot_cycle = "" );

  // special case, if analysing a hypnogram with no EDF
  void fudge( double es, int ne );

  timeline_t * timeline;
  
  // sleep stage information
  std::vector<sleep_stage_t> stages;
  std::vector<sleep_stage_t> original_stages;
  std::vector<int> epoch_n;
  std::vector<double> epoch_dur; // to track gap length
  std::vector<double> epoch_start; // elapsed seconds
  std::vector<bool> epoch_gap; // this 'epoch' is actually a gap

  // total epoch count
  int ne;  // total number of observed epochs (i.e. standard, from timeline)
  int ne_gaps; // epochs (ne) plus any extra gaps;
               // epoch_n[], stages[], etc correspond to ne_gaps, not ne;


  // standard epoch sizes (for observed epochs)
  double epoch_sec;
  double epoch_mins;
  double epoch_hrs;
  
  // number of conflicting epochs (set to missing)
  int n_conflicts;

  // number of 'fixed' due to leading/trailing sleep
  int n_fixed;

  // number of 'fixed' epochs due to lights on/off
  int n_lights_fixed;
  int n_lights_fixed_was_sleep;

  // number of W epochs set to ? if they lead/trail sleep
  int n_ignore_wake;
  
  // only calculate stats for first N epochs after X (make rest missing)
  // X is either start of recording, or sleep onset 
  // this is run *after* fixing above 'early sleep' issues
  int n_only_first_mins;
  std::string first_anchor; // SO or EDF


  // times 
  clocktime_t clock_start;
  clocktime_t clock_lights_out;
  clocktime_t clock_lights_on;
  clocktime_t clock_sleep_onset;
  clocktime_t clock_sleep_midpoint;
  clocktime_t clock_wake_time;
  clocktime_t clock_stop;

  // Epoch markers
  
  int first_sleep_epoch;    // (epoch)  
  int first_persistent_sleep_epoch;  // first persistent sleep epoch
  int final_wake_epoch;     // (epoch)

  // how to handle gaps (by default, missing, but could be 'wake')
  static sleep_stage_t gap_treatment;

  // statistics
  
  bool any_sleep;

  double TIB;  // time in bed : length of entire record, ignoring all staging
  double TRT;  // total recording time : from lights out to lights on
  
  double TWT;  // total wake time
  double TST;  // total sleep time
  double TpST;  // total persistent sleep time
  
  double FWT;  // final wake time  ( from final wake to end of test )

  // note - in the Luna output, we exclude FWT from these SPT measures
  double SPT;  // sleep period time [ Amount of time available for sleep after Sleep Onset , SPT = TRT - Sleep Latency ]
  
  double WASO; // wake after sleep onset  [ WASO = Wake epochs/2 - (Sleep Latency {+ Final Wake Time}) ]  (mins)  

  //double LOT; // LightsOn Time --> mins["L"]

  bool   starts_in_sleep; // study starts/ends in sleep
  bool   ends_in_sleep;
  
  double slp_eff_pct;   // sleep efficiency [ TST/TIB X 100 ]  (%)
  double slp_main_pct;  // sleep maintaince [ TST/SPT X 100 ] (%)
  double slp_eff2_pct;  // sleep maintaince 2 

  double slp_lat;   // sleep latency [ Sleep Latency = Lights Out - Sleep Onset ]  (mins)
  double per_slp_lat;  // latency to 10 mins sleep
  double rem_lat_mins;      // REM latency (minutes)
  double rem_lat_nowake_mins;  // REM latency (minutes), excluding W

  std::map<std::string,double> mins;
  std::map<std::string,double> pct;
  std::map<std::string,int> bout_n;
  std::map<std::string,double> bout_mean, bout_med, bout_5, bout_10, bout_max;

  std::set<bout_t> bouts;

  //
  // stage distribution stats (devel=T)
  //
  
  // mean time (i.e. central point, in mins from start) = mt
  // normalized mean time (0..1)                        = zmt
  
  // double mt_n1, mt_n2, mt_n3, mt_r, mt_w;
  // double zmt_n1, zmt_n2, zmt_n3, zmt_r, zmt_w;

  // conditional on stage, what proportion of +/-1, +/-2, +/-3, etc are of the same stage?
  //  -then can summarize those distributions
  //  -i.e. how to get at clumps? 
  
  // int locwin;
  // std::vector<double> loc_n1, loc_n2, loc_n3, loc_r, loc_w;
  
  // runs-test 
  double runs_pv5, runs_pv3;
  
  // sleep cycles (modified Feinberg & Floyd definitions)
  std::vector<bool> in_persistent_sleep;
  std::vector<int> sleep_cycle_number;
  std::vector<int> sleep_code; // 1 NR, 5 REM, 0 other
  
  int num_nremc;  // number of complete (NREM/REM) cycles
  double nremc_mean_duration;  // mean duration (mins)
  
  // duration/timing for each cycle
  std::map<int,double> nremc_duration;   // includes W/? too here
  std::map<int,double> nremc_nrem_duration;
  std::map<int,double> nremc_rem_duration;
  std::map<int,int> nremc_start_epoch;
  std::map<int,int> nremc_epoch_duration;
  
  // Summarize epochs w.r.t cycles  

  std::vector<double> cycle_pos_relative;
  std::vector<double> cycle_pos_absolute;

  // ascending/descending N2
  std::vector<double> n2_ascdesc;

  // other local context
  bool flanking_3class;
  std::vector<int> flanking;      // number of similar flanking epochs as 'x'
                                  // defined as min # of epochs until reach non-'x'
  std::vector<int> flanking_tot;  // total number of similar flanking epochs as 'x'

  std::vector<int> nearest_wake;  // distance to nearest wake

  //
  // Stage transitions (incl. N2-REM and N2-WAKE transitions)
  //

  // X is the number of epochs until next REM/WAKE
  // X_total is the size of that segment

  //                N3 N2  N2  N2  N2  R  R  R
  // NREM2REM        0  4   3   2   1  0  0  0
  // NREM2REM_TOTAL  0  4   4   4   4  0  0  0 

  // to count as a transition, we require to see at least this
  // many epochs pre and post   (the transitions[][] below counts all)

  int req_pre_post_epochs;
  
  std::vector<int> nrem2rem;    
  std::vector<int> nrem2rem_total; 
  std::vector<int> rem2nrem;    
  std::vector<int> rem2nrem_total; 
  
  std::vector<int> nrem2wake;    
  std::vector<int> nrem2wake_total;     
  std::vector<int> wake2nrem;    
  std::vector<int> wake2nrem_total; 

  std::vector<int> rem2wake;    
  std::vector<int> rem2wake_total;     
  std::vector<int> wake2rem;    
  std::vector<int> wake2rem_total; 

  std::map<sleep_stage_t,std::map<sleep_stage_t,int> > transitions;
  std::map<sleep_stage_t,std::map<sleep_stage_t,int> > transitions5;
  
  std::vector<bool> is_waso;      // distinguish WAKE during SLEEP from pre/post
    

};



#endif
