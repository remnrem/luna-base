

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
//  3) any annotations
//  4) (?todo: non-epoch based masks)


struct timeline_t 
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

  void create_discontinuous_timeline( const std::vector<uint64_t> & tps );

  void restructure( const std::set<int> & keep );

  
  //
  // Basic duration information for EDF
  //


  // Total file duration (microsecond units)

  uint64_t  total_duration_tp;
  
  // Final accessible timepoint (0-base) 
  // for EDF+ w/ gaps, then
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

  // for internal EDF+D, map 'file' rec --> implied rec (if iterating
  // through only retained records)
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


  //  uint64_t endpoint( int r ) const;


  //
  // Cache (for temporary, inter-command communication)
  //

  caches_t cache;

  
  //
  // Annoations
  //
  
  annotation_set_t annotations;
  
  void list_all_annotations( const param_t & param );

  void list_spanning_annotations( const param_t & param );

  void annot2signal( const param_t & param );
  
  void signal2annot( const param_t & param );

  void annot2cache( const param_t & param );
  
  void signal_means_by_annot( const param_t & param );

  int annot2sp( edf_t & edf , const std::string & astr , bool ,
		std::vector<interval_t> * , 
		std::vector<interval_t> * , 
		int * orig_n , 
		std::string ch , int sr = 0 );
  
  //
  // Hypnogram
  //

  hypnogram_t hypnogram;
 
 
  //
  // 2) The primary 'epoch'
  //
  
  int ensure_epoched() 
  {

    if ( epoched() ) return num_epochs();

    // otherwise, set to defaults

    int ne = set_epoch( globals::default_epoch_len , globals::default_epoch_len );

    logger << "  set epochs to default " 
	   << globals::default_epoch_len 
	   << " seconds, " << ne << " epochs\n";

    return ne;
  }
  
  bool epoched() const { return epoch_length_tp != 0; } 
  
  void unepoch() 
  {
    current_epoch = -1;
    
    // No EPOCHs
    epoch_length_tp = 0L;
    epoch_inc_tp = 0L;
    epoch_offset_tp = 0L;
    epoch_align_annots.clear();
    epoch_align_str = "";
    epochs.clear();
    
    // Masks
    clear_epoch_mask();
    mask_mode = 0; // default (0=mask; 1=unmask; 2=force)

    // Annotations
    clear_epoch_annotations();    
    
    // old/new epoch mapping
    clear_epoch_mapping();

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


  int whole_recording_epoch_dur(); 

  int set_epoch(const double s, const double o , const double offset = 0 , 
		const std::string align_str = "" , 
		const std::vector<std::string> * align_annots = NULL ) 
  { 
    if ( s <= 0 || o < 0 || offset < 0 ) 
      Helper::halt( "cannot specify negative epoch durations/increments/offsets");
    
    clear_epoch_annotations();
    epoch_length_tp = s * globals::tp_1sec;
    epoch_inc_tp = o * globals::tp_1sec;
    epoch_offset_tp = offset * globals::tp_1sec;
    
    epoch_align_str = align_str;
    if ( align_annots != NULL ) 
      epoch_align_annots = *align_annots;
    
    if ( epoch_length_tp == 0 || epoch_inc_tp == 0 ) 
      Helper::halt( "invalid epoch parameters" );
    first_epoch();
    return calc_epochs();
  }
  
  double epoch_length() const 
  { return (double)epoch_length_tp / globals::tp_1sec; }
  
  double epoch_inc() const 
  { return (double)epoch_inc_tp / globals::tp_1sec; }

  double epoch_offset() const 
  { return (double)epoch_offset_tp / globals::tp_1sec; }

  std::string align_string() const 
  { return epoch_align_str; }
  
  bool exactly_contiguous_epochs() const
  { return epoch_length_tp == epoch_inc_tp; } 
  
  double epoch_len_tp() const 
  { return epoch_length_tp ; }
  
  double epoch_increment_tp() const 
  { return epoch_inc_tp ; } 

  uint64_t epoch_len_tp_uint64_t() const
  { return epoch_length_tp ; }


  int calc_epochs();

  bool align_epochs( uint64_t * tp , int * rec , const std::set<uint64_t> & annots );
  
  int first_epoch()  
  { 

    // point to first epoch, and return number of non-masked epochs also

    if ( ! epoched() ) 
      {

	int ne = set_epoch( globals::default_epoch_len , globals::default_epoch_len );

	logger << "  set epochs to default " 
	       << globals::default_epoch_len 
	       << " seconds, " << ne << " epochs\n";
      }

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
    if ( e < 0 || e >= epochs.size() ) return interval_t(0LLU,0LLU);
    return epochs[e]; 
  }
 
  bool epoch_records( const int e , int * a , int * b ) const 
  {
    *a = *b = 0;
    std::map<int,std::set<int> >::const_iterator rr = epoch2rec.find( e );
    if ( rr == epoch2rec.end() ) return false;
    const std::set<int> & recs = rr->second;
    *a = *recs.begin();
    *b = *recs.rbegin();
    return true;
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

  // eval-based mask
  void apply_eval_mask( const std::string & , int mask_mode , const bool verbose = false );
  
  // other masks : randomly select up to 'n' epochs from the current set 
  void select_epoch_randomly( int n );

  // trim leading and trailing epochs (allow only n)
  void trim_epochs( std::string & , int );

  // regional mask
  void regional_mask( int x , int y );

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
  
  void apply_epoch_mask( annot_t * a , std::set<std::string> * values , bool include );
  
  void clear_epoch_mask( bool b = false );
  
  bool is_epoch_mask_set() const { return mask_set; }

  int  set_epoch_mask( const int e , const bool b = true );

  void set_epoch_mask_mode( const int m ) 
  {
    mask_mode = m;
  }

  int epoch_mask_mode() const { return mask_mode; } 

  bool masked( const int e ) const
  { 
    return mask[e]; 
  }

  //  void mask2annot( const std::string & path , const std::string & tag , const bool with_id = true );

  void add_mask_annot( const std::string & tag );
  
  void dumpmask( const param_t & );

  //
  // Channel-specific epoch masks (ch/ep mask)
  //
  
  static void proc_chep( edf_t & edf , param_t & param );

  bool is_chep_mask_set() const { return chep.size() != 0; } 
  
  void clear_chep_mask() { chep.clear(); } 

  std::map<int,std::set<std::string> > make_chep_copy() const { return chep; }

  void set_chep_mask( const int e , const std::string & s ) { chep[ display_epoch( e ) ].insert(s); } 

  void merge_chep_mask( const std::map<int,std::set<std::string> > & m ) 
  {
    if ( chep.size() == 0 ) { chep = m ; return; } 
    std::map<int,std::set<std::string> >::const_iterator ii = m.begin();
    while ( ii != m.end() )
      {
	std::set<std::string>::const_iterator jj = ii->second.begin();
	while ( jj != ii->second.end() )
	  {
	    chep[ ii->first ].insert( *jj );
	    ++jj;
	  }
	++ii;
      }
  }
  

  bool unset_chep_mask( const int e , const std::string & s ) 
  { 
    // return T if anything removed
    int e1 = display_epoch( e );    
    std::map<int,std::set<std::string> >::iterator ii = chep.find( e1 ) ;
    if ( ii == chep.end() ) return false;
    std::set<std::string>::iterator jj = ii->second.find( s );
    if ( jj == ii->second.end() ) return false;
    ii->second.erase( jj );
    return true; 
  } 

  void dump_chep_mask( signal_list_t , bool );
  
  bool masked( const int e , const std::string & s ) const 
  {
    std::map<int,std::set<std::string> >::const_iterator ee = chep.find( display_epoch( e ) );
    if ( ee == chep.end() ) return false;
    return ee->second.find( s ) != ee->second.end() ;
  }

  // save/load cheps
  void read_chep_file( const std::string & f , bool reset = true );

  void write_chep_file( const std::string & f ) const;

  // query
  std::vector<std::string> masked_channels( const int e , const signal_list_t & ) const;
  
  std::vector<std::string> unmasked_channels( const int e , const signal_list_t & ) const;

  signal_list_t masked_channels_sl( const int e , const signal_list_t & ) const;
  
  signal_list_t unmasked_channels_sl( const int e , const signal_list_t & ) const;

  // sets main epoch mask
  void collapse_chep2epoch( signal_list_t signals , const double pct , const int k );
  
  // doesn't change epoch-mask, returns a channel list
  // optionally, if a channel is designated as 'bad', then set all to bad
  // and/or,     if a channel is designated as 'good', then set all to good
  
  signal_list_t collapse_chep2ch( signal_list_t signals , 
				  const double pct , const int k  , 
				  const bool bad_set_all_bad = true , 
				  const bool good_set_all_good = true );
  
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

  uint64_t valid_tps(  const interval_t & );


  //
  // ANNOTATE command implementation (annots.cpp)
  //

  void annotate( param_t & );
  
  //
  // EPOCH annotations
  //
  
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

  void annotate_epoch( const std::string & label , int e )
  {
    
    // do we need to remap the epoch?
    
    if ( has_epoch_mapping() )
      {
        // off-the-grid  
        if ( epoch_curr2orig.find( e ) == epoch_curr2orig.end() ) 
	  return;
	
        // convert query to the original mapping
        e = epoch_curr2orig.find( e )->second;
      }

    eannots[ label ][ e ] = true;
  }


/*   void annotate_epochs( const std::string & label , bool b ) */
/*   { */
/*     const int ne = num_epochs(); */
/*     eannots[ label ].clear(); */
/*     for (int e = 0 ; e < ne ; e++) eannots[ label ][ e ] = b ; */
/*   } */
  
/*   void annotate_epochs( const std::vector<std::string> & s ) */
/*   { */
/*     if ( s.size() != num_total_epochs() )  */
/*       Helper::halt( "internal error: incorrect number of epochs in annotate_epochs()" ); */
/*     for (int i=0;i<s.size();i++) eannots[ s[i] ][ i ] = true;  */
/*   } */
  
  void clear_epoch_annotations();


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
  
  // does annotation 'k' exist at all? 
  
  bool epoch_annotation(const std::string & k ) const
  {
    return eannots.find( k ) != eannots.end() ;
  }
  
  // does EPOCH 'e' contain annotation 'k'?
  // where 'e' is in the current 0..ne epoch form, 
  // and will be remapped if necessary
  
  bool epoch_annotation(const std::string & k, int e) const
  {

    // look up this annotation 'k'
    std::map<std::string,std::map<int,bool> >::const_iterator ii = eannots.find( k );

    // annotation k does not exist anywhere
    if ( ii == eannots.end() ) return false;
    
    // do we need to remap the epoch? 
    if ( has_epoch_mapping() ) 
      {
	// off-the-grid
	if ( epoch_curr2orig.find( e ) == epoch_curr2orig.end() ) return false;
	// convert query to the original mapping
	e = epoch_curr2orig.find( e )->second;
      }
    
    // now we have the correct original-EDF epoch number, do
    // we have an flag? if not, means FALSE
    if ( ii->second.find( e ) == ii->second.end() ) return false;

    // return the boolean value (i.e. could still be set FALSE explicitly)
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
  int display_epoch(int e) const
    {      
      if ( ! has_epoch_mapping() ) return e+1;
      if ( epoch_curr2orig.find(e) == epoch_curr2orig.end() ) return -1;
      return epoch_curr2orig.find(e)->second + 1 ;
    }


  int display2curr_epoch(int e) const 
  {
    if ( ! has_epoch_mapping() ) return e-1;    
    if ( epoch_orig2curr.find(e-1) == epoch_orig2curr.end() ) return -1;
    return epoch_orig2curr.find(e-1)->second ;
  }
  
  static bool discontinuity( const std::vector<uint64_t> & t , int sr, int sp1, int sp2 );

  std::map<int,bool> spanning_epoch_masks( const int r ) 
  {
    std::map<int,bool> rec;
    std::map<int,std::set<int> >::const_iterator ii = rec2epoch.find( r );
    if ( ii == rec2epoch.end() ) return rec;
    std::set<int>::const_iterator jj = ii->second.begin();
    while ( jj != ii->second.end() )
      {
	rec[ *jj ] = masked_epoch( *jj );
	++jj;
      }
    return rec;
  }
  

 private:

  
  //
  // Data
  //
  
  edf_t      * edf;
  
  interval_t   window;
  
  uint64_t     epoch_length_tp;  
  
  uint64_t     epoch_inc_tp;    

  uint64_t     epoch_offset_tp;
  
  
  std::string epoch_align_str; // un-tokenized version of below, for quick check by EVAL when calling proc_epoch()
  std::vector<std::string> epoch_align_annots; // for EDF+D w/ 'align' (i.e. cannot have a single offset number)

  std::vector<interval_t> epochs;

  int current_epoch;

  std::vector<bool> mask;
  
  bool mask_set;
  
  int mask_mode;
  
  // epoch -> ch ; presence implies mask set
  std::map<int,std::set<std::string> > chep;
  
  // epoch to record mapping
  std::map<int,std::set<int> > epoch2rec;
  std::map<int,std::set<int> > rec2epoch;
  
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
