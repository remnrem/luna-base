
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

struct hypnogram_t
{
  
  hypnogram_t() { } 
  void construct( timeline_t * t , const bool verbose , const std::vector<std::string> & s );
  void construct( timeline_t * t , const bool verbose , const std::string sslabel = "SleepStage" );
  void calc_stats( const bool verbose ); // verbose == STAGES vs HYPNO
  void output( const bool verbose );

  // special case, if analysing a hypnogram with no EDF
  void fudge( double es, int ne );

  timeline_t * timeline;
  
  // sleep stage information
  std::vector<sleep_stage_t> stages;
  std::vector<int> epoch_n;
  
  // times 
  clocktime_t clock_lights_out;
  clocktime_t clock_lights_on;

  clocktime_t clock_sleep_onset;
  clocktime_t clock_sleep_midpoint;
  clocktime_t clock_wake_time;

  // Epoch markers
  
  int first_sleep_epoch;    // (epoch)  
  int first_persistent_sleep_epoch;  // first persistent sleep epoch
  int final_wake_epoch;     // (epoch)


  // statistics
  
  double TIB;  // time in bed : length of entire record, ignoring all staging
  double TRT;  // total recording time : from lights out to lights on

  double TWT;  // total wake time
  double TST;  // total sleep time
  double TpST;  // total persistent sleep time
  
  double FWT;  // final wake time  ( from final wake to end of test )
  double SPT;  // sleep period time [ Amount of time available for sleep after Sleep Onset , SPT = TRT - Sleep Latency ]

  double WASO; // wake after sleep onset  [ WASO = Wake epochs/2 - (Sleep Latency {+ Final Wake Time}) ]  (mins)  
  
  
  double slp_eff_pct;   // sleep efficiency [ TST/TIB X 100 ]  (%)
  double slp_main_pct;  // sleep maintaince [ TST/SPT X 100 ] (%)
  double slp_eff2_pct;  // sleep maintaince 2 

  double slp_lat;   // sleep latency [ Sleep Latency = Lights Out - Sleep Onset ]  (mins)
  double per_slp_lat;  // latency to 10 mins sleep
  double rem_lat_mins;      // REM latency (minutes)

  
  double mins_wake;  // minutes awake
  double mins_n1;  // minutes N1
  double mins_n2;  // etc
  double mins_n3;  // 
  double mins_n4;  // 
  double mins_rem;  // 
  double mins_other;

  double pct_n1;   // % of sleep that is N1
  double pct_n2;   // etc
  double pct_n3;   // 
  double pct_n4;   // 
  double pct_rem;  // 
  double pct_other;
  
  
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
  std::vector<int> flanking;      // number of similar flanking epochs as 'x'
  std::vector<int> nearest_wake;  // distance to nearest wake

  // N2-REM and N2-WAKE transitions

  // X is the number of epochs until next REM/WAKE
  // X_total is the size of that segment

  //                N3 N2  N2  N2  N2  R  R  R
  // NREM2REM        0  4   3   2   1  0  0  0
  // NREM2REM_TOTAL  0  4   4   4   4  0  0  0 

  std::vector<int> nrem2rem;    
  std::vector<int> nrem2rem_total; 
  
  std::vector<int> nrem2wake;    
  std::vector<int> nrem2wake_total; 

  std::vector<bool> is_waso;      // distinguish WAKE during SLEEP from pre/post
    

};



#endif
