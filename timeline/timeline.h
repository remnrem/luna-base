
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

#ifndef __TIMELINE_H__
#define __TIMELINE_H__

#include <iostream>
#include <iomanip>
#include "intervals/intervals.h"
#include "annot/annot.h"
#include "helper/helper.h"
#include "defs/defs.h"

struct edf_t;
struct timeline_t;

//
// General class to represent a hypnogram and perform calculations on it
// linked to a timeline_t's epochs
//

void dummy_hypno();

struct hypnogram_t
{
  
  hypnogram_t() { } 
  void construct( timeline_t * t , const std::vector<std::string> & s );
  void construct( timeline_t * t , const std::string sslabel = "SleepStage" );
  void calc_stats();
  void output();

  // special case, if analysing a hypnogram with no EDF
  void fudge( double es, int ne );

  timeline_t * timeline;
  
  std::vector<sleep_stage_t> stages;
  
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


//
// For a single EDF provides a unified timeline
// Annotations and epochs will be defined with respect to this
//

// A 'timeline' coordinates several things for an EDF:
//  1) the primary 'epoch'
//  2) an epoch-based mask
//  3) any annotations
//  4) (?todo: non-epoch based masks)


class timeline_t 
{
  
  friend struct hypnogram_t;
  
 public:
  
  timeline_t( edf_t * p )
    {      
      edf = p;            
      unepoch();      
    } 
  

  //
  // Primary creation/alteration of timeline
  //

  void init_timeline( bool redo = false );

  void re_init_timeline() { init_timeline( true ); }

  void restructure( const std::set<int> & keep );

  
  //
  // Basic duration information for EDF
  //


  // Total file duration (microsecond units)

  uint64_t  total_duration_tp;
  
  // Final timepoint; for EDF+ w/ gaps, then
  // total_duration_tp <= last_time_point_tp + 1 

  uint64_t  last_time_point_tp;
  
  //
  // Iteration over records in the 'retained' timeline
  //

  int first_record() const; // -1 if empty
  int next_record(const int r) const; // -1 if at end
  bool retained(const int r ) const;
  
  
  //
  // Record-level time-point information
  //
  
      
  std::map<uint64_t,int> tp2rec;
  std::map<int,uint64_t>  rec2tp;
  std::map<int,uint64_t>  rec2tp_end;
  
  bool interval2records( const interval_t & interval , 
			 uint64_t srate , 
			 int * start_rec , 
			 int * start_smp , 
			 int * stop_rec , 
			 int * stop_smp ) const;

  // wrapper around the above
  std::set<int> records_in_interval( const interval_t & interval ) const;

  interval_t record2interval( int r ) const;
  
  

  // for a given record/sample pair, give the time-point (in tp)

  uint64_t timepoint( int r , int s = 0 , int nsamples = 0 ) const;

  uint64_t endpoint( int r ) const;


  //
  // Annoations
  //
  
  annotation_set_t annotations;
  
  void list_all_annotations( const param_t & param );


  //
  // Hypnogram
  //

  hypnogram_t hypnogram;
 
 
  //
  // 2) The primary 'epoch'
  //
  
  
  bool epoched() const { return epoch_length_tp != 0; } 
  
  void unepoch() 
  {
    current_epoch = -1;
    
    // No EPOCHs
    epoch_length_tp = 0L;
    epoch_overlap_tp = 0L;
    epochs.clear();
    
    // Masks
    clear_epoch_mask();
    mask_mode = 0; // default (0=mask; 1=unmask; 2=force)

    // Annotations
    clear_epoch_annotations();    
    
    // record/epoch mappings
    rec2epoch.clear();
    epoch2rec.clear();
  }
 
  int reset_epochs()
  {
    // called automatically after RESTRUCTURE only
    // this doesn't wipe any epoch-level annotations, 
    // as they are fixed to the original file-EDF epochs
    // (for this particularly epoch size) in any case

    first_epoch();
    return calc_epochs();
  }
  
  int set_epoch(const double s, const double o ) 
  { 
    clear_epoch_annotations();
    epoch_length_tp = s * globals::tp_1sec;
    epoch_overlap_tp = o * globals::tp_1sec;
    if ( epoch_length_tp == 0 || epoch_overlap_tp == 0 ) 
      Helper::halt( "invalid epoch parameters" );
    first_epoch();
    return calc_epochs();
  }
  
  double epoch_length() const 
  { return (double)epoch_length_tp / globals::tp_1sec; }
  
  double epoch_overlap() const 
  { return (double)epoch_overlap_tp / globals::tp_1sec; }
  
  int calc_epochs();
  
  int first_epoch()  
  { 
    // point to first epoch, and return number of non-masked epochs also
    if ( ! epoched() ) Helper::halt( "no EPOCHs set" );
    current_epoch = -1; 
    return num_epochs();
  } 
  
  int next_epoch()  
  { 
    // return the next unmasked epoch
    while (1)
      {
	++current_epoch;
	if ( current_epoch == epochs.size() ) return -1;
	if ( ! mask_set ) break;
	if ( ! mask[ current_epoch ] ) break; 	
      }
    return current_epoch;
  }
  
  int next_epoch_ignoring_mask()  
  { 
    ++current_epoch;
    if ( current_epoch == epochs.size() ) return -1;
    return current_epoch;
  }

  // all unmasked epochs
  int num_epochs() const 
  { 
    if ( ! mask_set ) return epochs.size();
    int r = 0;
    for (int i=0;i<mask.size();i++)
      if ( ! mask[i] ) ++r;
    return r;
  }

  // all epochs
  int num_total_epochs() const 
  { 
    return epochs.size(); 
  }

  interval_t epoch( const int e ) const
  { 
    return epochs[e]; 
  }
 


  
  //
  // Epoch Masks
  //

  // mask as simply 0/1 in a file, which must match # of epochs
  void load_mask( const std::string & , bool exclude = true );

  // interval list corresponding to original EDF times
  void load_interval_list_mask( const std::string & , bool exclude = true );

  // directly specify epoch-level masks
  void apply_simple_epoch_mask( const std::set<std::string> & , const std::string & , bool );
  
  // convert from full annotation to epoch-level masks
    
  void apply_epoch_include_mask( annot_t * a , std::set<std::string> * values = NULL )
  {
    apply_epoch_mask( a , values , true );
  }

  void apply_epoch_exclude_mask( annot_t * a , std::set<std::string> * values = NULL )
  {
    apply_epoch_mask( a , values , false );
  }
  
  // behavior if annotation missing, i.e. an 'empty' mask
  void apply_empty_epoch_mask( const std::string & , bool include );
  
  // other masks : randomly select up to 'n' epochs from the current set 
  void select_epoch_randomly( int n );

  // mask all after 'n' epochs
  void select_epoch_first( int n );

  // mask range of epochs from a to b inclusive
  void select_epoch_range( int a , int b );
  
  // flip all values of current epoch mask
  void flip_epoch_mask();
  
  // select all EPOCHs until we come across an EPOCH that does /not/ have the 'str' annotation
  void select_epoch_until_isnot( const std::string & str );

  // select EPOCHs within a run of 'n' (either side) other similarly annotated epochs
  void select_epoch_within_run( const std::string & str , int n );
  
  void apply_epoch_mask( annot_t * a , std::set<std::string> * values , bool include );
  
  void clear_epoch_mask( bool b = false ) 
  {
    mask.clear();
    mask_set = b;  // i.e. if b==T, equivalent to masking all entries
    mask.resize( epochs.size() , b );
    if ( epoched() )
      std::cerr << " reset all " << epochs.size() << " epochs to be " << ( b ? "masked" : "included" ) << "\n";
  }
  
  bool is_epoch_mask_set() const { return mask_set; }

  int  set_epoch_mask( const int e , const bool b = true ) 
  {
    mask_set = true;
    if (e < 0 || e >= mask.size() ) Helper::halt( "internal error setting mask" );
    
    bool original = mask[e];

    // implement mask mode
    // only mask
    if      ( mask_mode == 0 ) 
      {
	if ( (!original) && b ) mask[e] = true;  // default
      }
    else if ( mask_mode == 1 ) // 'unmask' --> only unmask
      {
	if ( original && (!b) ) mask[e] = false;
      }
    else if ( mask_mode == 2 ) // 'force' --> set either way
      {
	mask[e] = b ; // force (default)   
      }
    
    // teturn 0 if no change;
    // return +1 if set a mask (N->Y)
    // return -1 if freed a mask (Y->N)

    if ( original == mask[e] ) return 0;
    return mask[e] ? 1 : -1 ;     
  }
  

  void set_epoch_mask_mode( const int m ) 
  {
    mask_mode = m;
  }

  int epoch_mask_mode() const { return mask_mode; } 

  bool masked( const int e ) const
  { 
    return mask[e]; 
  }

  void mask2annot( const std::string & path , const std::string & tag );
  void dumpmask();


  //
  // Generic masks
  //
  
  bool masked_timepoint( uint64_t ) const;

  bool interval_is_completely_masked( const interval_t & interval ) const 
   { return masked_interval( interval , true ); } 
  
  bool interval_is_completely_unmasked( const interval_t & interval ) const 
   { return ! masked_interval( interval , false ); } 
  
  bool interval_overlaps_masked_region( const interval_t & interval ) const 
  { return masked_interval( interval , false ); } 

  bool interval_overlaps_unmasked_region( const interval_t & interval ) const 
  { return ! masked_interval( interval , true ); } 

  bool interval_start_is_masked( const interval_t & interval ) const 
  { 
    bool start_masked = false;
    masked_interval( interval , false , &start_masked );
    return start_masked;
  }

  bool masked_interval( const interval_t & , const bool all_masked , bool * start_masked = NULL ) const;

  bool masked_record( int r ) const;

  bool masked_epoch( int e ) const;


  //
  // EPOCH annotations
  //
  
  // i.e. typically used for 'stage' information

  void annotate_epochs( const std::string & label , 
			const std::string & annot_label , 
			const std::set<std::string> & values );
  
  void annotate_epochs( const std::string & label , bool b )
  {
    const int ne = num_epochs();
    eannots[ label ].clear();
    for (int e = 0 ; e < ne ; e++) eannots[ label ][ e ] = b ;
  }
  
  void annotate_epochs( const std::vector<std::string> & s )
  {
    // if ( has_epoch_mapping() ) Helper::halt( "cannot EPOCH-ANNOT after RESTRUCTURE-ing an EDF" );
    if ( s.size() != num_total_epochs() ) Helper::halt("incorrect number of epochs specified in EPOCH-ANNOT");
    for (int i=0;i<s.size();i++) eannots[ s[i] ][ i ] = true; 
  }
  
  void clear_epoch_annotations() 
  {
    if ( eannots.size() > 0 ) 
      std::cerr << " clearing all epoch-annotations\n";
    eannots.clear();
  }

  // Return all epoch annotations

  std::set<std::string> epoch_annotations() const
    {
      std::set<std::string> r;
      std::map<std::string,std::map<int,bool> >::const_iterator ii = eannots.begin();
      while ( ii != eannots.end() )
	{
	  r.insert( ii->first );
	  ++ii;
	}
      return r;
    }
  

  // does EPOCH 'e' contain annotation 'k'?

  bool epoch_annotation(const std::string & k, int e) const
  {

    std::map<std::string,std::map<int,bool> >::const_iterator ii = eannots.find( k );
    if ( ii == eannots.end() ) return false;
    if ( has_epoch_mapping() ) 
      {
	// off-the-grid
	if ( epoch_curr2orig.find( e ) == epoch_curr2orig.end() ) return false;
	// convert query to the original mapping
	e = epoch_curr2orig.find( e )->second;
      }
    if ( ii->second.find( e ) == ii->second.end() ) return false;
    return ii->second.find( e )->second; 
  }
  

  //
  // Misc
  //

  interval_t wholetrace() const; 
  
  //
  // Epoch mappings
  //

  void clear_epoch_mapping()
  {
    epoch_orig2curr.clear();
    epoch_curr2orig.clear();
  }

  bool has_epoch_mapping() const
  {
    return epoch_orig2curr.size() != 0 ;
  }

  void set_epoch_mapping();
  
  int original_epoch(int e)
  {
    if ( ! has_epoch_mapping() ) return e;
    if ( epoch_curr2orig.find(e) == epoch_curr2orig.end() ) return -1;
    return epoch_curr2orig.find(e)->second;
  }

  // 1-based epoch mapping
  int display_epoch(int e) 
    {
      if ( ! has_epoch_mapping() ) return e+1;
      if ( epoch_curr2orig.find(e) == epoch_curr2orig.end() ) return -1;
      return epoch_curr2orig.find(e)->second + 1 ;
/*       std::stringstream ss; */
/*       ss << std::setfill('0') << std::setw(6) << epoch_curr2orig.find(e)->second + 1 ; */
/*       return ss.str(); */
    }


  static bool discontinuity( const std::vector<uint64_t> & t , int sr, int sp1, int sp2 );

  bool spans_epoch_boundary( const interval_t & ) const ; 


 private:

  
  //
  // Data
  //
  
  edf_t      * edf;
  
  interval_t   window;
  
  uint64_t     epoch_length_tp;  
  
  uint64_t     epoch_overlap_tp;    
  
  std::vector<interval_t> epochs;

  std::vector<bool> mask;
  
  bool mask_set;
  
  int mask_mode;
  
  int current_epoch;

  std::map<int,std::set<int> > epoch2rec;
  std::map<int,std::set<int> > rec2epoch;
  
  // original <--> current epoch mappings
  // i.e. track if masks have been applied
  
  std::map<int,int> epoch_orig2curr;
  std::map<int,int> epoch_curr2orig;
  int               orig_epoch_size;

  // Epoch annotations 
  // type --> [ epoch-2-bool ]
  // where epoch is *always* with regard to the original value

  // boolean epoch-based annotations
  std::map<std::string,std::map<int,bool> > eannots;
  
  // string epoch-based annotations: primary for sleep staging
  //  std::map<std::string,std::vector<std::string> > str_eannots
};



#endif
