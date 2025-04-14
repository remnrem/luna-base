

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
#include <vector>
#include <string>
#include <set>
#include <map>

#include "helper/helper.h"
#include "helper/logger.h"
#include "timeline/hypno.h"
#include "timeline/cache.h"

#include "edf/signal-list.h"
#include "defs/defs.h"
#include "annot/annot.h"

struct edf_t;

struct timeline_t;

struct interval_t;

extern logger_t logger;


//
// For a single EDF provides a unified timeline
// Annotations and epochs will be defined with respect to this
//

// A 'timeline' coordinates several things for an EDF:
//  1) the primary 'epoch' 
//  2) an epoch-based mask
//  3) channel-epoch (CHEP) masks
//  4) contains general annotations
//  5) contains caches

// standard epochs: all of same size and fixed increment
// generic epochs (e.g. typically defined by annotations, etc) can be any size/location
// certain functions (e.g. hypnogram) will require standard (and non-overlapping) epochs
// to make sense



struct timeline_t 
{
  
  friend struct hypnogram_t;

 public:
  
  timeline_t( edf_t * p );
  
  ~timeline_t();
  
  
  //
  // Primary creation/alteration of timeline
  //

  void init_timeline( bool redo = false );

  void re_init_timeline() { init_timeline( true ); }

  void create_discontinuous_timeline( const std::vector<uint64_t> & tps );

  void restructure( const std::set<int> & keep , const bool preserve_cache );

  

  //
  // Basic duration information for EDF
  //
  

  // Total file duration (microsecond units)

  uint64_t  total_duration_tp;
  
  // Final accessible timepoint (0-base) 
  // for EDF+ w/ gaps, total_duration_tp <= last_time_point_tp + 1 
  
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

  // for internal EDF+D, map 'file' rec --> implied rec
  // (if iterating through only retained records)

  std::map<int,int>  rec2orig_rec;
  
  bool interval2records( const interval_t & interval , 
			 uint64_t srate , 
			 int * start_rec , 
			 int * start_smp , 
			 int * stop_rec , 
			 int * stop_smp ) const;

  // wrapper around the above
  std::set<int> records_in_interval( const interval_t & interval ) const;

  interval_t record2interval( int r ) const;
  
  // collapse edf+d time --> edf standard time-point (from EDF start)
  interval_t collapse( const interval_t & interval ) const;
		       
  // for a given record/sample pair, give the time-point (in tp)
  uint64_t timepoint( int r , int s = 0 , int nsamples = 0 ) const;

  // for an internal EDF+D, give current tp (tp1 = in elapsed tp)
  // given an original, external tp
  bool remap_timepoint( const uint64_t & tp , uint64_t * tp1 );



  //
  // Cache (for temporary, inter-command communication)
  //

  caches_t cache;

  
  //
  // Annotations
  //
  
  annotation_set_t * annotations;
  
  void list_all_annotations( const param_t & param );

  void list_spanning_annotations( const param_t & param );

  void annot_crosstabs( const param_t & param );
  
  void annot2signal( const param_t & param );
  
  void signal2annot( const param_t & param );

  void annot2cache( const param_t & param );

  void cache2annot( const param_t & param );
  
  void signal_means_by_annot( const param_t & param );

  void set_annot_metadata( const param_t & param );
  
  int annot2sp( edf_t & edf , const std::string & astr , bool ,
		std::vector<interval_t> * , 
		std::vector<interval_t> * , 
		int * orig_n , 
		std::string ch , int sr = 0 );

  // ------------------------------------------------------------
  //
  // Hypnogram
  //
  // ------------------------------------------------------------
  
  hypnogram_t hypnogram;
 

  
  // ------------------------------------------------------------
  //
  // Epochs
  //
  // ------------------------------------------------------------
  
  int ensure_epoched();
  
  bool epoched() const;

  void unepoch(); 

  // reset_epochs()
  //  called automatically after RESTRUCTURE only
  //  this doesn't wipe any epoch-level annotations, 
  //  as they are fixed to the original file-EDF epochs
  //  (for this particularly epoch size) in any case

  int whole_recording_epoch_dur(); 

  int set_epoch(const double s, const double o , const uint64_t offset = 0LLU , 
		const std::string align_str = "" , 
		const std::vector<std::string> * align_annots = NULL );

  bool generic_epochs() const;

  bool fixed_epoch_length() const;
  
  double epoch_length() const;
  
  double epoch_inc() const;

  double epoch_offset() const;

  bool epoch_any_offset() const;

  std::string align_string() const;
  
  bool exactly_contiguous_epochs() const;
  
  double epoch_len_tp() const;
  
  uint64_t epoch_increment_tp() const;

  uint64_t epoch_len_tp_uint64_t() const;

  // initial setting: standard epochs
  int calc_epochs(); 

  // initial setting: generic (annot/cache-based epochs)
  int calc_epochs_generic_from_annots( param_t & param ); 

  // called after RE to remap epochs
  int reset_epochs(); 
  
  void debug_dump_epochs();

  void output_epoch_info( const bool verbose , const bool show_masked = false );
    
  bool align_epochs( uint64_t * tp , int * rec , const std::set<uint64_t> & annots );
  
  int first_epoch();

  int next_epoch(); 
  
  int next_epoch_ignoring_mask();

  // all unmasked epochs
  int num_epochs() const; 

  // all epochs
  int num_total_epochs() const;

  interval_t epoch( const int e ) const;
 
  bool epoch_records( const int e , int * a , int * b ) const;

  //
  // Epoch mappings
  //

  void clear_epoch_mapping();

  bool has_epoch_mapping() const;

  void set_epoch_mapping();
  
  int original_epoch(int e);

  // 1-based epoch mapping
  int display_epoch(int e) const;

  int display2curr_epoch(int e) const;
  
  std::map<int,bool> spanning_epoch_masks( const int r );


  // Epoch Masks

  // mask as simply 0/1 in a file, which must match # of epochs
  void load_mask( const std::string & , bool exclude = true );

  // interval list corresponding to original EDF times
  void load_interval_list_mask( const std::string & , bool exclude = true );

  // directly specify epoch-level masks - no longer done
  //  void apply_simple_epoch_mask( const std::set<std::string> & , const std::string & , bool );

  // 
  // ----- legacy annotation-based epoch masks -----
  //

  // primary mask function
  // void apply_epoch_mask( annot_t * a , std::set<std::string> * values , bool include );

  // // wrappers to call apply_epoch_mask()
  // // convert from full annotation to epoch-level masks
  // void apply_epoch_include_mask( annot_t * a , std::set<std::string> * values = NULL );

  // void apply_epoch_exclude_mask( annot_t * a , std::set<std::string> * values = NULL );
  
  // // behavior if annotation missing, i.e. an 'empty' mask
  // void apply_empty_epoch_mask( const std::string & , bool include );

  // 
  // ----- revised annotation-based epoch masks -----
  //

  // new generic epoch mask to handle all the above cases
  void apply_epoch_mask2( const std::map<annot_t *,std::set<std::string> > & annots ,
			  const std::set<std::string> & fullspan , 
			  const std::string & alabel , 
			  bool or_mode, bool include );
			  

  //
  // ----- eval-expression epoch masks -----
  //
  
  // eval-based mask
  void apply_eval_mask( const std::string & , int mask_mode , const bool verbose = false );

  //
  // ----- other masks ----- 
  //
  
  // randomly select up to 'n' epochs from the current set 
  void select_epoch_randomly( int n );

  // trim leading and trailing epochs (allow only n)
  void trim_epochs( std::string & , int );
  
  // retain around a core set of annots (kind of inverse of trim, but specify what to keep)
  void retain_epochs( const std::set<std::string> & labels );

  // regional mask
  void regional_mask( int x , int y );

  // unmask interior
  void unmask_interior();

  // mask all after 'n' epochs
  void select_epoch_first( int n );

  // mask range of epochs from a to b inclusive
  void select_epoch_range( int a , int b , bool include );

  // mask set of epochs
  void select_epoch_range( const std::set<int> & e , bool include );
  
  // flip all values of current epoch mask
  void flip_epoch_mask();
  
  // select all EPOCHs until we come across an EPOCH that does /not/ have the 'str' annotation
  void select_epoch_until_isnot( const std::string & str );

  // select EPOCHs within a run of 'n' (either side) other similarly annotated epochs
  void select_epoch_within_run( const std::string & str , int n );

  // clear all epoch masks
  void clear_epoch_mask( bool b = false );

  //
  // ---- epoch-mask helper function -----
  //
  
  bool is_epoch_mask_set() const;

  int  set_epoch_mask( const int e , const bool b = true );

  void set_epoch_mask_mode( const int m );

  int epoch_mask_mode() const;

  bool masked( const int e ) const;

  void dumpmask( const param_t & );

  bool elapsed_seconds_to_spanning_epochs( const double t1, const double t2, int * e1 , int * e2 );

  //
  // ---- legacy -----
  //

  //  MASK2ANNOT function (ANNOT-MASK - currently note used)
  //    see if otherwise implemented?
  
  void add_mask_annot( const std::string & tag );
  

  // ------------------------------------------------------------
  //
  // Channel-specific epoch masks (ch/ep mask)
  //
  // ------------------------------------------------------------
  
  static void proc_chep( edf_t & edf , param_t & param );

  bool is_chep_mask_set() const;
  
  void clear_chep_mask();

  std::map<int,std::set<std::string> > make_chep_copy() const;

  void set_chep_mask( const int e , const std::string & s );

  void merge_chep_mask( const std::map<int,std::set<std::string> > & m ) ;

  bool unset_chep_mask( const int e , const std::string & s ) ;

  void dump_chep_mask( signal_list_t , bool );
  
  bool masked( const int e , const std::string & s ) const;

  void read_chep_file( const std::string & f , bool reset = true );

  void write_chep_file( const std::string & f ) const;

  // sets main epoch mask
  void collapse_chep2epoch( signal_list_t signals , const double pct , const int k );
  
  // doesn't change epoch-mask, returns a channel list
  // optionally, if a channel is designated as 'bad', then set all to bad
  // and/or,     if a channel is designated as 'good', then set all to good
  
  signal_list_t collapse_chep2ch( signal_list_t signals , 
				  const double pct , const int k  , 
				  const bool bad_set_all_bad = true , 
				  const bool good_set_all_good = true );
  
  // query
  std::vector<std::string> masked_channels( const int e , const signal_list_t & ) const;
  
  std::vector<std::string> unmasked_channels( const int e , const signal_list_t & ) const;

  signal_list_t masked_channels_sl( const int e , const signal_list_t & ) const;
  
  signal_list_t unmasked_channels_sl( const int e , const signal_list_t & ) const;


  // ------------------------------------------------------------
  //
  // Generic masks
  //
  // ------------------------------------------------------------
  
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

  uint64_t valid_tps(  const interval_t & );


  // ------------------------------------------------------------
  //
  // OVERLAP command implementation (annots.cpp)
  //
  // ------------------------------------------------------------

  // FIXME... REDUNDANT??? or not?
  //  void annotate( param_t & );
  

  // ------------------------------------------------------------  
  //
  // EPOCH annotations
  //
  // ------------------------------------------------------------

  // ?? redundant?
  
  // i.e. a lightweight annot, separate from annot_t, that is only
  // used internally for store SleepStage information
  
  // *internally*, numbering works with respect to the original, EDF
  // epoch scheme, i.e. and so should follow any remapping
  // all function calls for eannots-related stuff should use the current
  // epoch mapping however (i.e. remapping done internally if needed)
  
  void annotate_epochs( const std::string & label , 
			const std::string & annot_label , 
			const std::set<std::string> & values );
  
  // should be used with current 0..(ne-1) mapping, will
  // be converted to original epoch mapping if needed

  void annotate_epoch( const std::string & label , int e );

  void clear_epoch_annotations();

  // Return all epoch annotations
  
  std::set<std::string> epoch_annotations() const;
  
  // does annotation 'k' exist at all? 
  
  bool epoch_annotation(const std::string & k ) const;
  
  // does EPOCH 'e' contain annotation 'k'?
  // where 'e' is in the current 0..ne epoch form, 
  // and will be remapped if necessary
  
  bool epoch_annotation(const std::string & k, int e) const;

  //
  // Misc helpers
  //

  bool check( const std::string & cmd ) const;

  interval_t wholetrace( const bool silent = false ) const; 
  
  static bool discontinuity( const std::vector<uint64_t> & t , int sr, int sp1, int sp2 );
  
  std::set<interval_t> segments();

  std::set<interval_t> gaps( const std::set<interval_t> & segs );

 private:

  
  //
  // Data
  //
  
  edf_t      * edf;
  
  interval_t   window;
  
  uint64_t     epoch_length_tp;  
  
  uint64_t     epoch_inc_tp;    

  uint64_t     epoch_offset_tp;

  // generic epochs

  std::set<std::string> epoch_generic_param_annots;
  double                epoch_generic_param_w1; // left
  double                epoch_generic_param_w2; // right
  int                   epoch_generic_param_set_point; // 0 none, 1=start, 2=mid, 3=end
  double                epoch_generic_param_min_epoch_size;
  
  double                epoch_generic_param_shift; // secs -/+
  double                epoch_generic_param_trunc; // secs (+ve) from end
  
  // epoch alignment
  
  std::string epoch_align_str; // un-tokenized version of below, for quick check by EVAL when calling proc_epoch()
  std::vector<std::string> epoch_align_annots; // for EDF+D w/ 'align' (i.e. cannot have a single offset number)

  std::vector<interval_t> epochs;

  // for generic-epoch (vs standard epoch) case:
  std::vector<std::string> epoch_labels;

  // T if fully standard epochs
  bool standard_epochs;

  // generic (nonstandard) but still fixed-size epochs?
  bool fixed_size_epochs;
  
  int current_epoch;

  std::vector<bool> mask;
  
  bool mask_set;
  
  int mask_mode;
  
  // epoch -> ch ; presence implies mask set
  std::map<int,std::set<std::string> > chep;
  
  // epoch to record mapping
  std::map<int,std::set<int> > epoch2rec;
  // kludge...
public:
  std::map<int,std::set<int> > rec2epoch;
private:


  // original <--> current epoch mappings
  // i.e. track if masks have been applied
  
  std::map<int,int> epoch_orig2curr;
  std::map<int,int> epoch_curr2orig;

  // Epoch annotations 

  // type --> [ epoch-2-bool ]
  // where epoch is *always* with regard to the original value
  
  // boolean epoch-based annotations
  std::map<std::string,std::map<int,bool> > eannots;
  
};



#endif
