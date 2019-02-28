

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

#include "timeline.h"

#include "edf/edf.h"
#include "defs/defs.h"
#include "db/db.h"
#include "miscmath/crandom.h"
#include "eval.h"
#include "intervals/intervals.h"
#include "annot/annot.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "helper/token-eval.h"

extern writer_t writer;

extern logger_t logger;


bool is_rem( sleep_stage_t s ) { return s == REM; } 
bool is_nrem( sleep_stage_t s ) { return s == NREM1 || s == NREM2 || s == NREM3 || s == NREM4; } 
bool is_nrem1( sleep_stage_t s ) { return s == NREM1; } 
bool is_nrem23( sleep_stage_t s ) { return s == NREM2 || s == NREM3; } 
bool is_wake( sleep_stage_t s ) { return s == WAKE; }
bool is_sleep( sleep_stage_t s ) { return s == NREM1 || s == NREM2 || s == NREM3 || s == NREM4 || s == REM ; } 


// helper function: check if there is a discontinuity in a timeline
bool timeline_t::discontinuity( const std::vector<uint64_t> & t , int sr , int sp1, int sp2  )
{
  if ( sp2 < sp1 ) return true;
  if ( sp1 < 0 || sp2 >= t.size() ) return true;
  uint64_t x = ( globals::tp_1sec / sr ) * ( sp2-sp1 );
  uint64_t y = t[sp2] - t[sp1];  
  return x != y ; 
}

int timeline_t::first_record() const
{
  if ( rec2tp.size() == 0 ) return -1; //empty
  return rec2tp.begin()->first;
}

int timeline_t::next_record(const int r) const
{
  std::map<int,uint64_t>::const_iterator i = rec2tp.find(r);
  if ( i == rec2tp.end() ) return -1;
  ++i;
  if ( i == rec2tp.end() ) return -1;
  return i->first;
}

bool timeline_t::retained(const int r ) const
{
  std::map<int,uint64_t>::const_iterator i = rec2tp.find(r);
  return i != rec2tp.end();
}


void timeline_t::init_timeline( bool okay_to_reinit ) 
{
  
  if ( rec2tp.size() != 0 && ! okay_to_reinit ) 
    Helper::halt( "internal error: cannot re-init timeline" );
  
  tp2rec.clear();
  rec2tp.clear();
  rec2tp_end.clear();

  
  clear_epoch_mapping();
  orig_epoch_size = -1;
  

  //
  // Continuous timeline?
  //
  
  if ( edf->header.continuous )
    {
      
      total_duration_tp = 
	(uint64_t)edf->header.nr * edf->header.record_duration_tp;
      last_time_point_tp = total_duration_tp - 1LLU;
      
      uint64_t tp = 0;
      
      for (int r = 0;r < edf->header.nr;r++)
	{	  
	  tp2rec[tp] = r;	
       	  rec2tp[r] = tp;
	  rec2tp_end[r] = tp + edf->header.record_duration_tp - 1LLU;
	  tp += edf->header.record_duration_tp;
	}            

    }

  //
  // For a discontinuous EDF (which implies EDF+)
  //
    
  else  
    {
      
      // 1) Does an INDEX edf (.edf.idx) exist for this EDF?
      
      // 2) Otherwise, we need to read the whole EDF...
      // i.e. records are still non-overlapping
      
      total_duration_tp = 
	(uint64_t)edf->header.nr * edf->header.record_duration_tp;
      
      // okay to use header.nr here, as this will only be called
      // once, on first loading the EDF (i.e. so nr==nr_all as
      // no records have yet been removed)
      
      for (int r = 0;r < edf->header.nr;r++)
	{
	  uint64_t tp = edf->timepoint_from_EDF(r);
	  tp2rec[tp] = r;
	  rec2tp[r] = tp;
	  rec2tp_end[r] = last_time_point_tp = tp + edf->header.record_duration_tp - 1LLU;
	  // last_time_point_tp will be updated, 
	  // and end up being thelast (i.e. record nr-1).
	}
    }
}



bool timeline_t::spans_epoch_boundary( const interval_t & interval ) const 
{
  // if the timeline is epoched, return T/F as to whether this interval
  // spans a boundary

  // e.g. use case: for restricting spindle/SO detection only to events that 
  // fall within a single epoch
  
  if ( ! epoched() ) return false; 
  
  // if the interval is discontinuous, it must, by definition span a boundary, 
  // as all epochs must be continuous
  
  return false;
}


void timeline_t::restructure( const std::set<int> & keep )
{

  // the restructured EDF header should be in place at this point
  // here, it should not matter whether the original was continuous or not
  // i.e. it is now discontinuous, if records are being dropped
  // this change will already have been reflected in the header 
  // i.e.  edf_t::restructure

  total_duration_tp = 
    (uint64_t)edf->header.nr * edf->header.record_duration_tp;      
  last_time_point_tp = 0;
  
  std::map<uint64_t,int> copy_tp2rec;    
  std::map<int,uint64_t> copy_rec2tp;    
  std::map<int,uint64_t> copy_rec2tp_end;

  int r = first_record();

  while ( r != -1 ) 
    {
      if ( keep.find(r) != keep.end() )
	{	  
	  uint64_t tp = rec2tp[r];
	  copy_rec2tp[r] = tp;
	  copy_rec2tp_end[r] = rec2tp_end[r];
	  copy_tp2rec[tp] = r;
	  if ( rec2tp_end[r] > last_time_point_tp ) 
	    last_time_point_tp = rec2tp_end[r];
	}
      
      r = next_record(r);
    }
  
  // copy over

  tp2rec     = copy_tp2rec;
  rec2tp     = copy_rec2tp;
  rec2tp_end = copy_rec2tp_end;

  // reset epochs (but retain epoch-level annotations)
  reset_epochs();
  logger << " retaining " << num_epochs() << " epochs\n";
}


interval_t timeline_t::record2interval( int r ) const
{ 
  std::map<int,uint64_t>::const_iterator ll = rec2tp.find(r);
  std::map<int,uint64_t>::const_iterator uu = rec2tp_end.find(r);
  if ( ll == rec2tp.end() ) return interval_t(0,0);
  return interval_t( ll->second , uu->second );
}






bool timeline_t::interval2records( const interval_t & interval , 
				   uint64_t n_samples_per_record , 
				   int * start_rec , 
				   int * start_smp , 
				   int * stop_rec , 
				   int * stop_smp ) const

{
  
  //   std::cout << "i2r: interval = " << interval << "\n";

  //
  // Note: here we want to find records/samples that are inclusive w.r.t. the interval
  // so change coding of stop being 1 unit past the end below
  //

  if ( interval.stop == 0 ) 
    Helper::halt( "internal error in timeline()" );

  uint64_t stop_tp = interval.stop - 1LLU;
  
  if ( interval.start >= stop_tp ) return false;

  //
  // For a continuous timeline, given time-points can
  // straightforwardly calculate record/sample
  //

  if ( edf->header.continuous )
    {
      
      // get initial records/samples

      uint64_t start_record = interval.start / edf->header.record_duration_tp;
      uint64_t start_offset = interval.start % edf->header.record_duration_tp;
      uint64_t start_sample = 
	( start_offset / (double)edf->header.record_duration_tp ) * n_samples_per_record ;
      
      if ( start_sample >= n_samples_per_record ) start_sample = (uint64_t)n_samples_per_record - 1LLU; 
      
      
      // get final records/samples
      
      uint64_t stop_record = stop_tp / edf->header.record_duration_tp;
      uint64_t stop_offset = stop_tp % edf->header.record_duration_tp;
      uint64_t stop_sample = 
	( stop_offset / (double)edf->header.record_duration_tp ) * n_samples_per_record ;

      if ( stop_sample >= n_samples_per_record ) stop_sample = n_samples_per_record - 1LLU;
      
      // pass back to calling function
      
      *start_rec = (int)start_record;
      *start_smp = (int)start_sample;

      *stop_rec = (int)stop_record;
      *stop_smp = (int)stop_sample;

    }
  else
    {

      //
      // For a discontinuous EDF+ we need to search 
      // explicitly across record timepoints
      //
      
      //
      // Get first record that is not less than start search point (i.e. equal to or greater than)
      //
      
      std::map<uint64_t,int>::const_iterator lwr = tp2rec.lower_bound( interval.start ); 
           
      //
      // This will find the first record AFTER the start; thus we should skip one record back;
      // This should never be the first record, but check in case...
      //

      // Does the search point fall outside of a record?
      
      bool in_gap = false;
      
      if ( lwr != tp2rec.begin() ) 
	{

	  // go back one record
	  --lwr;
	  uint64_t previous_rec_start = lwr->first;
	  uint64_t previous_rec_end   = previous_rec_start + edf->header.record_duration_tp - 1LLU;

	  // does the start point fall within this previous record?
	  if ( interval.start >= previous_rec_start && interval.start <= previous_rec_end ) 
	    in_gap = false;
	  else
	    {
	      in_gap = true;
	      ++lwr;
	    }
	}
      else if ( lwr == tp2rec.begin() )
       	{
	  // If the search point occurs before /all/ records, need to indicate that we are in a gap
	  // also
	  if ( interval.start < lwr->first ) 
	    in_gap = true;	      
	}
      
      // problem? return empty record set
      if ( lwr == tp2rec.end() ) 
	{
	  *start_rec = 0;
	  *start_smp = 0;	  
	  *stop_rec = 0;
	  *stop_smp = 0;
	  return false;
	}

      
      *start_rec = lwr->second;
      
      if ( in_gap )
	*start_smp = 0; // i.e. use start of this record, as it is after the 'true' start site
      else
	{	
	  uint64_t start_offset = interval.start % edf->header.record_duration_tp;
	  uint64_t start_sample = 
	    ( start_offset / (double)edf->header.record_duration_tp ) * n_samples_per_record ;
	  if ( start_sample >= n_samples_per_record ) start_sample = n_samples_per_record - 1LLU; 
	  *start_smp = (int)start_sample;
	}
      
      //
      // for upper bound, find the record whose end is equal/greater *greater* 
      // 
      
      std::map<uint64_t,int>::const_iterator upr = tp2rec.upper_bound( stop_tp ); 
      
      //
      // this should have return one past the one we are looking for 
      // i.e. that starts *after* the search point
      //
      
      if ( upr != tp2rec.begin() ) --upr;  
      
      *stop_rec  = upr->second;
      
      // get samples within (as above)      
      uint64_t previous_rec_start = upr->first;
      uint64_t previous_rec_end   = previous_rec_start + edf->header.record_duration_tp - 1;
      in_gap = ! ( stop_tp >= previous_rec_start && stop_tp <= previous_rec_end );
      
      if ( in_gap )
	*stop_smp = n_samples_per_record - 1; // set to last point
      else
	{	  
	  uint64_t stop_offset = stop_tp % edf->header.record_duration_tp;
	  uint64_t stop_sample = 
	    ( stop_offset / (double)edf->header.record_duration_tp ) * n_samples_per_record ;
	  if ( stop_sample >= n_samples_per_record ) stop_sample = n_samples_per_record - 1LLU;
	  *stop_smp = (int)stop_sample;
	}
      
    }

  // std::cout << "recs = " << *start_rec << " " << *stop_rec << "\n";
  // std::cout << "smps = " << *start_smp << " " << *stop_smp << "\n";

  // if the interval is in a gap, we will not get any records here (in fact, stop < start), so fix
  // this

  if ( *start_rec > *stop_rec ) { *stop_rec = *start_rec; } 
  if ( *start_rec == *stop_rec && *start_smp > *stop_smp ) { *stop_smp = *start_smp; } 

  return true;
  
}



int timeline_t::calc_epochs()
{

  // EPOCHS have to be in the oringal time-series units (i.e. that
  // correspond to the CONTINUOUS EDF). In the case of a DISCONTINUOUS
  // EDF, we require that epochs are of the specified time on the
  // /reduced/ time-scale (i.e. nominally, the interval may have >
  // than specified epoch length, i.e. if it contains a gap...)
  
  // we'll also populate the rec2epoch and epoch2rec mappings
  
  epochs.clear();
  
  mask.clear();
  
  rec2epoch.clear();
  
  epoch2rec.clear();

  if ( edf->header.continuous )
    {
      
      uint64_t s = 0;
      
      while ( 1 ) 
	{
	  
	  // get end of interval: for this purpose (of finding records)
	  // we set last point of 
	  uint64_t end = s + epoch_length_tp - 1LLU;
	  
	  // done? [ skip any final epochs that do not fit into the frame ] 	  
	  if ( end >= total_duration_tp ) break;
	  
	  // add to list, but with end as +1 past end
	  interval_t interval( s , end + 1LLU );
	  epochs.push_back( interval );

	  // find matching records (within interval)
	  interval_t search_interval( s , end );
	  
	  // maked records in this epoch
	  int start_record = search_interval.start / edf->header.record_duration_tp;  
	  int stop_record = search_interval.stop / edf->header.record_duration_tp;
	  int e = epochs.size()-1;
	  for (int r=start_record; r<=stop_record; r++) 
	    {
	      epoch2rec[ e ].insert( r );
	      rec2epoch[ r ].insert( e );
	    }
	  
	  // shift to next interval
	  s += epoch_overlap_tp;
	}
    }
  else
    {
      //
      // Epochs for the discontinuous case:
      //

      // Temporary issues, both can be fixed

      // 1) No overlapping epochs allowed (can fix)
      // (nb. 'overlap' better interpreted as 'increment')

      if ( epoch_overlap_tp != epoch_length_tp ) 
	Helper::halt( "cannot have overlapping epochs with EDF+D" );

      // 2) Epoch length must be >= record length

      if ( epoch_length_tp < edf->header.record_duration_tp )
	Helper::halt( "epoch length must be greater or equal to record length" );

      int r = first_record();
      
      if ( r == -1 ) return 0;
      
      uint64_t estart = rec2tp[r];
      uint64_t curr = 0;

      // for epoch2rec, rec2epoch mapping
      int e = 0;
            
      while ( 1 ) 
	{
	  // current EPOCH start 's'
	  // putative EPOCH end -- before or after current record end?

	  uint64_t rec_start = rec2tp[r];
	  uint64_t rec_end   = rec2tp_end[r];

	  uint64_t rec_dur = rec_end - rec_start + 1LLU;
	  
	  // if epoch will end within this record
	  if ( curr + rec_dur >= epoch_length_tp ) 
	    {
	      
	      uint64_t estop = rec_start + ( epoch_length_tp - curr - 1LLU );
	      
	      // add to list of epochs
	      interval_t saved_interval( estart , estop + 1LLU );
	      epochs.push_back( saved_interval );
	      
	      
	      interval_t interval( estart , estop );
	      
	      // record mappings
	      rec2epoch[r].insert(e);
	      epoch2rec[e].insert(r);

	      // move on
	      ++e;
	      
	      // check this is within current record, else get next
	      if ( estop < rec_end )
		{
		  estart = estop + 1LLU;
		  curr = rec_end - estart + 1LLU;
		  
		  // and mark this too
		  rec2epoch[r].insert(e);
		  epoch2rec[e].insert(r);
	      
		  // Note:: this assumes that there will not be another new
		  // epoch within this record...
		  r = next_record(r);
		  if ( r == -1 ) break;		  
		}
	      else
		{
		  // advance to mext record
		  r = next_record(r);
		  if ( r == -1 ) break;
		  curr = 0;
		  estart = rec2tp[r]; 
		}
	    }
	  else 
	    {
	      curr += rec_dur;	      
	      
	      rec2epoch[r].insert(e);
	      epoch2rec[e].insert(r);

	      r = next_record(r);
	      if ( r == -1 ) break;
	    }
	}
    }


  // reset counter
  current_epoch = -1;
  mask.resize( epochs.size() , false );
  mask_set = false;
  mask_mode = 0; 

  // all done
  return epochs.size();    
}


interval_t timeline_t::wholetrace() const
{  
  // end is defined as 1 past the last time point
  return interval_t( 0 , last_time_point_tp + 1LLU );
}


void timeline_t::clear_epoch_mask( bool b )
{
  mask.clear();
  mask_set = b;  // i.e. if b==T, equivalent to masking all entries
  mask.resize( epochs.size() , b );
  if ( epoched() )
    logger << " reset all " << epochs.size() << " epochs to be " << ( b ? "masked" : "included" ) << "\n";
}

int timeline_t::set_epoch_mask( const int e , const bool b ) 
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


void timeline_t::clear_epoch_annotations()
{
  if ( eannots.size() > 0 ) 
    logger << " clearing all epoch-annotations\n";
  eannots.clear();
}

void timeline_t::apply_empty_epoch_mask( const std::string & label , bool include )
{
  
  // this is requested if the annotation is missing
  // i.e. returns match == F for every epoch; treat as specified by include and mask_mode

  // include T/F   means whether a 'match' means having (T) versus not-having (F) the annotation
  // mask_mode will already have been set
  
  mask_set = true;
  
  const int ne = epochs.size();
  
  int cnt_mask_set = 0;
  int cnt_mask_unset = 0;
  int cnt_unchanged = 0;
  int cnt_now_unmasked = 0;
  int cnt_basic_match = 0;  // basic count of matches, whether changes mask or not
  
  for (int e=0;e<ne;e++)
    {
      
      bool matches = false;
      
      // set new potential mask, depending on match_mode
      
      bool new_mask = mask[e];

      if ( include ) 
	{
	  if      ( mask_mode == 0 ) new_mask = matches;   // mask-if
	  else if ( mask_mode == 1 ) new_mask = !matches;  // unmask-if
	  else if ( mask_mode == 2 ) new_mask = matches ;  // if
	}
      else
	{
	  if      ( mask_mode == 0 ) new_mask = !matches;  // mask-ifnot
	  else if ( mask_mode == 1 ) new_mask = matches;   // unmask-ifnot
	  else if ( mask_mode == 2 ) new_mask = ! matches; // ifnot
	}

      int mc = set_epoch_mask( e , new_mask );

      if      ( mc == +1 ) ++cnt_mask_set;
      else if ( mc == -1 ) ++cnt_mask_unset;
      else                 ++cnt_unchanged;
      
      if ( !mask[e] ) ++cnt_now_unmasked;
      
    }

  logger << " based on " << label << " " << cnt_basic_match << " epochs match; ";
  logger << cnt_mask_set << " newly masked, "   
	 << cnt_mask_unset << " unmasked, " 
	 << cnt_unchanged << " unchanged\n";
  logger << " total of " << cnt_now_unmasked << " of " << epochs.size() << " retained\n";
  
  // mask, # epochs masked, # epochs unmasked, # unchanged, # total masked , # total epochs
  
  writer.level( label  , "EPOCH_MASK" );

  writer.var( "N_MATCHES"    , "Number of matching epochs" );
  writer.var( "N_MASK_SET"   , "Number of epochs newly masked" ); 
  writer.var( "N_MASK_UNSET" , "Number of epochs newly unmasked" );
  writer.var( "N_UNCHANGED"  , "Number of epochs unchanged by this mask" );
  writer.var( "N_RETAINED"   , "Number of epochs retained for analysis" );
  writer.var( "N_TOTAL"      , "Total number of epochs" );

  writer.value( "N_MATCHES"    , cnt_basic_match  );
  writer.value( "N_MASK_SET"   , cnt_mask_set     );
  writer.value( "N_MASK_UNSET" , cnt_mask_unset   );
  writer.value( "N_UNCHANGED"  , cnt_unchanged    );
  writer.value( "N_RETAINED"   , cnt_now_unmasked );
  writer.value( "N_TOTAL"      , epochs.size()    );

  writer.unlevel( "EPOCH_MASK" );

}

void timeline_t::apply_epoch_mask( annot_t * a , std::set<std::string> * values , bool include )
{
  
  // include T/F   means whether a 'match' means having (T) versus not-having (F) the annotation
  
  // mask_mode will already have been set
  
  // if 'values' is NULL, then we just use presence of an annotation,
  // rather than looking at the instance ID
  
  bool value_mask = values != NULL;
  
  mask_set = true;

  const int ne = epochs.size();
  
  //
  // We do not clear the mask here, as we want to allow multiple
  // filters to be added on top of oneanther
  //

  int cnt_mask_set = 0;
  int cnt_mask_unset = 0;
  int cnt_unchanged = 0;
  int cnt_now_unmasked = 0;
  int cnt_basic_match = 0;  // basic count of matches, whether changes mask or not
  
  for (int e=0;e<ne;e++)
    {

      interval_t interval = epoch( e );
      
      annot_map_t events = a->extract( interval );
      
      bool matches = false;
      
      if ( value_mask ) 
	{
	  // do any of the instance IDs match any of the values?
	  annot_map_t::const_iterator ii = events.begin();
	  while ( ii != events.end() )
	    {		  
	      const instance_idx_t & instance_idx = ii->first;	      
	      if ( values->find( instance_idx.id ) != values->end() )
		{
		  matches = true;
		  break;
		}
	      ++ii;
	    }
	}
      else 
	{	  
	  matches = events.size() > 0 ;
	}

      // count basic matches

      if ( matches ) ++cnt_basic_match;
      
      // set new potential mask, depending on match_mode
      
      bool new_mask = mask[e];

      if ( include ) 
	{
	  if      ( mask_mode == 0 ) new_mask = matches;   // mask-if
	  else if ( mask_mode == 1 ) new_mask = !matches;  // unmask-if
	  else if ( mask_mode == 2 ) new_mask = matches ;  // if
	}
      else
	{
	  if      ( mask_mode == 0 ) new_mask = !matches;  // mask-ifnot
	  else if ( mask_mode == 1 ) new_mask = matches;   // unmask-ifnot
	  else if ( mask_mode == 2 ) new_mask = ! matches; // ifnot
	}

      int mc = set_epoch_mask( e , new_mask );

      if      ( mc == +1 ) ++cnt_mask_set;
      else if ( mc == -1 ) ++cnt_mask_unset;
      else                 ++cnt_unchanged;
      
      if ( !mask[e] ) ++cnt_now_unmasked;
      
    }
  
  logger << " based on " << a->name << ( value_mask ? "[" + Helper::stringize( *values , "|" ) + "]" : "" )  
	 << " " << cnt_basic_match << " epochs match; ";

  logger << cnt_mask_set << " newly masked, " 
	 << cnt_mask_unset << " unmasked, " 
	 << cnt_unchanged << " unchanged\n";
  logger << " total of " << cnt_now_unmasked << " of " << epochs.size() << " retained\n";

  
  // mask, # epochs masked, # epochs unmasked, # unchanged, # total masked , # total epochs
  
  writer.level( a->name , "EPOCH_MASK" );

  writer.var( "N_MATCHES"    , "Number of matching epochs" );
  writer.var( "N_MASK_SET"   , "Number of epochs newly masked" ); 
  writer.var( "N_MASK_UNSET" , "Number of epochs newly unmasked" );
  writer.var( "N_UNCHANGED"  , "Number of epochs unchanged by this mask" );
  writer.var( "N_RETAINED"   , "Number of epochs retained for analysis" );
  writer.var( "N_TOTAL"      , "Total number of epochs" );

  writer.value( "N_MATCHES"    , cnt_basic_match  );
  writer.value( "N_MASK_SET"   , cnt_mask_set     );
  writer.value( "N_MASK_UNSET" , cnt_mask_unset   );
  writer.value( "N_UNCHANGED"  , cnt_unchanged    );
  writer.value( "N_RETAINED"   , cnt_now_unmasked );
  writer.value( "N_TOTAL"      , epochs.size()    );

  writer.unlevel( "EPOCH_MASK" );
}


bool timeline_t::masked_timepoint( uint64_t a ) const
{

  Helper::halt( "masked_timepoint() not implemented" );

  if ( ! edf->header.continuous ) 
    Helper::halt( "masked_timepoint() not implemented for EDF+D yet" );

    if ( ! mask_set ) return false;
  
  int e1 = MiscMath::position2leftepoch( a , epoch_length_tp, epoch_overlap_tp , mask.size() );
  int e2 = MiscMath::position2rightepoch( a , epoch_length_tp, epoch_overlap_tp , mask.size() );
  
  // above functions return -1 if the tp is off the map
  // (or epochs are not overlapping/contiguous); here it is
  // effectively 'masked'

  if ( e1 == -1 || e2 == -1 ) return true;

  if ( e1 >= mask.size() || e2 >= mask.size() ) 
    Helper::halt( "internal error, timeline : e > mask.size()" 
		  + Helper::int2str(e1) + " " + Helper::int2str(e1) 
		  + " " + Helper::int2str( (int)mask.size() ) );   
  
  


  // do /any/ of these mask epochs that span this position have 
  // a positive mask set? 
  
  bool m = false;

  for (int e=e1;e<=e2;e++) 
    if ( mask[e] ) m = true;
  
  return m;

}


bool timeline_t::masked_interval( const interval_t & interval , bool all_masked , bool * start_masked ) const
{
  
  // if all_masked,   returns T if /all/ of interval falls within masked regions
  // if not,          returns T if interval falls in at least one masked region

  if ( start_masked != NULL ) *start_masked = false;

  if ( ! mask_set ) 
    {
      return false;
    }
  
  if ( edf->header.continuous )
    {

      int eleft = MiscMath::position2leftepoch( interval.start , epoch_length_tp, epoch_overlap_tp , mask.size() );
      int eright = MiscMath::position2rightepoch( interval.stop , epoch_length_tp, epoch_overlap_tp , mask.size() );
      
      //std::cout << "e1e2 = " << eleft << "  " << eright << "\n";

      if ( start_masked != NULL )
	{
	  if ( eleft == -1 || mask[eleft] ) *start_masked = true;
	}

      if ( eleft == -1 || eright == -1 ) return true;
      
      // above functions return -1 if the tp is off the map
      // (or epochs are not overlapping/contiguous); here it is
      // effectively 'masked'
      
      for (int e=eleft;e<=eright;e++) 
	{
	  if ( all_masked && ! mask[e] ) return false;
	  if ( mask[e] && ! all_masked ) return true;
	}

    }

  else // for EDF+D
    {
      std::set<int> records = records_in_interval( interval );
      std::set<int>::const_iterator rr = records.begin();
      while ( rr != records.end() )
	{
	  const std::set<int> & epochs = rec2epoch.find( *rr )->second;

	  std::set<int>::const_iterator ee = epochs.begin();
	  
	  if ( start_masked != NULL )
	    {
	      if ( mask[ *ee ] ) *start_masked = true;
	    }
	  
	  while ( ee != epochs.end() )
	    {
	      if ( all_masked && ! mask[ *ee ] ) return false;
	      if ( mask[ *ee ] && ! all_masked ) return true;
	      ++ee;
	    }
	  ++rr;
	}      
    }
  
  if ( all_masked ) return true;
  else return false;
  
}



std::set<int> timeline_t::records_in_interval( const interval_t & interval ) const
{
  
  int start_rec = 0 , stop_rec = 0;
  int start_smp = 0 , stop_smp = 0;

  const int srate = 100;  // will not matter, as we only consider whole records here

  std::set<int> recs;
  
  bool any = interval2records( interval , srate , &start_rec , &start_smp , &stop_rec , &stop_smp );
  
  if ( ! any ) return recs;
  
  int r = start_rec;
  while ( r != -1 ) 
    {
      recs.insert(r);
      r = next_record(r);
    }
  return recs;
}



bool timeline_t::masked_record( int r ) const
{
  
  if ( ! mask_set ) return false;
  
  std::map<int,std::set<int> >::const_iterator rr = rec2epoch.find(r);
  if ( rr == rec2epoch.end() ) return true; // i.e. out of bounds
  
  const std::set<int> & epochs = rr->second;
  std::set<int>::const_iterator ee = epochs.begin();
  while ( ee != epochs.end() )
    {
      if ( mask[ *ee ] ) return true;
      ++ee;
    }
  return false;   
}


bool timeline_t::masked_epoch( int e ) const
{
  if ( ! mask_set ) return false;
  if ( e < 0 || e >= mask.size() ) return true; // out-of-bounds, so 'masked'
  return mask[e];
}

// flip all values of a mask
// i.e. to /include/ artifactual epochs only
void timeline_t::flip_epoch_mask()
{

  if ( ! mask_set ) return;

  const int ne = epochs.size();
  
  int cnt_mask_set = 0;
  int cnt_mask_unset = 0;
  int cnt_unchanged = 0;
  int cnt_now_unmasked = 0;
  
  // flip all (i.e. every epoch will change)
  for (int e=0;e<ne;e++)
    {
      
      mask[e] = ! mask[e];
      
      if ( mask[e] ) ++cnt_mask_set;
      else ++cnt_mask_unset;

    }
  
  logger << " flipped all epoch masks\n";
  logger << " total of " << cnt_mask_unset << " of " << epochs.size() << " retained\n";

}

// other masks : randomly select up to 'n' epochs from the current set 
void timeline_t::select_epoch_randomly( int n )
{

  mask_set = true;
  
  // from the unmasked set, pick at random 'n' (or as many as possible)
  std::vector<int> unmasked;
  
  const int ne = epochs.size();
  
  for (int e=0;e<ne;e++)
    if ( ! mask[e] ) unmasked.push_back(e);
  
  std::set<int> selected;
  
  int s = 0;
  
  const int num_unmasked = unmasked.size();
  
  const int n_to_select = num_unmasked < n ? num_unmasked : n;

  while ( s < n_to_select )
    {
      int rnd = CRandom::rand( num_unmasked );
      int sel = unmasked[ rnd ];

      if ( selected.find( sel ) == selected.end() )
	{
	  selected.insert( sel );
	  ++s;
	}
    }
  
  
  int cnt_mask_set = 0;
  int cnt_mask_unset = 0;
  int cnt_unchanged = 0;
  int cnt_now_unmasked = 0;


  // mask everything that was /not/ in the selected set
  for (int e=0;e<ne;e++)
    {
      
      if ( selected.find(e) == selected.end() ) 
	{
	  int mc = set_epoch_mask( e , true );
	  if ( mc == +1 ) ++cnt_mask_set;
	  else if ( mc == -1 ) ++cnt_mask_unset;
	  else ++cnt_unchanged;
	}
      
      if ( ! mask[e] ) ++cnt_now_unmasked;
    }

  logger << " randomly selected up to " << n << " epochs; ";

  logger << cnt_mask_set << " newly masked " 
	 << cnt_mask_unset << " unmasked and " 
	 << cnt_unchanged << " unchanged\n";
  logger << " total of " << cnt_now_unmasked << " of " << epochs.size() << " retained\n";

}


// other masks: select epochs from 'a' to 'b' inclusive (include=T)
// otherwise do the opposite
void timeline_t::select_epoch_range( int a , int b , bool include )
{

  if ( a > b ) 
    {
      int tmp = b;
      b = a;
      a = tmp;
    }

  mask_set = true;

  const int ne = epochs.size();
  
  int cnt_mask_set = 0;
  int cnt_mask_unset = 0;
  int cnt_unchanged = 0;
  int cnt_now_unmasked = 0;
  
  // mask everything that was /not/ in the selected set
  for (int e=0;e<ne;e++)
    {
      // use base-1 coding of epochs
      const int epoch = e+1;

      bool match = include ? 
	epoch < a || epoch > b : 
	epoch >= a && epoch <= b ;
      
      if ( match ) 
	{
	  int mc = set_epoch_mask( e , true );
	  if ( mc == +1 ) ++cnt_mask_set;
	  else if ( mc == -1 ) ++cnt_mask_unset;
	  else ++cnt_unchanged;
	}
      
      if ( ! mask[e] ) ++cnt_now_unmasked;
    }

  if ( include )
    logger << " selecting epochs from " << a << " to " << b << "; ";
  else
    logger << " masking epochs from " << a << " to " << b << "; ";

  logger << cnt_mask_set << " newly masked, " 
	 << cnt_mask_unset << " unmasked, " 
	 << cnt_unchanged << " unchanged\n";
  logger << " total of " << cnt_now_unmasked << " of " << epochs.size() << " retained\n";

}




// other masks : select up to 'n' epochs from the start of the record
void timeline_t::select_epoch_first( int n )
{
  
  mask_set = true;
  
  // from the unmasked set, pick at random 'n' (or as many as possible)
  //  std::vector<int> unmasked;
  
  const int ne = epochs.size();
     
  int cnt_mask_set = 0;
  int cnt_mask_unset = 0;
  int cnt_unchanged = 0;
  int cnt_now_unmasked = 0;

  // mask everything that was /not/ in the selected set
  for (int e=0;e<ne;e++)
    {
      if ( e >= n )
	{
	  int mc = set_epoch_mask( e , true );
	  if ( mc == +1 ) ++cnt_mask_set;
	  else if ( mc == -1 ) ++cnt_mask_unset;
	  else ++cnt_unchanged;
	}
      
      if ( ! mask[e] ) ++cnt_now_unmasked;
    }

  logger << " selecting up to " << n << " epochs for start; ";
  logger << cnt_mask_set << " newly masked, "  
	 << cnt_mask_unset << " unmasked, " 
	 << cnt_unchanged << " unchanged\n";
  logger << " total of " << cnt_now_unmasked << " of " << epochs.size() << " retained\n";
}




// select only EPOCHs that are in contiguous runs of EPOCH /str/ (i.e. +1 means one either side)

void timeline_t::select_epoch_within_run( const std::string & str , int b )
{

  if ( b < 1 ) Helper::halt( "epoch border must be 1 or greater" );

  mask_set = true;

  const int ne = epochs.size();  

  int cnt_mask_set = 0;
  int cnt_mask_unset = 0;
  int cnt_unchanged = 0;
  int cnt_now_unmasked = 0;

  for (int e=0;e<ne;e++)
    {  
      
      bool set_mask = false;

      if ( ! epoch_annotation( str , e ) ) 
	set_mask = true;
      
      if ( ! set_mask )
	{	  
	  int cnt = 0;
	  
	  int current = e;
	  for (int bwk=0;bwk<b;bwk++)
	    {
	      --current;
	      if ( epoch_annotation( str , current ) ) ++cnt;
	    }
	  
	  current = e;
	  for (int fwd=0;fwd<b;fwd++)
	    {
	      ++current;
	      if ( epoch_annotation( str , current ) ) ++cnt;
	    }
	  
	  if ( cnt < b * 2 ) set_mask = true;
      
	}
      
      int mc = set_epoch_mask( e , set_mask );
      if ( mc == +1 ) ++cnt_mask_set;
      else if ( mc == -1 ) ++cnt_mask_unset;
      else ++cnt_unchanged;
      
      if ( ! mask[e] ) ++cnt_now_unmasked;
      
    }
  
  logger << " based on " << str << " with " << b << " flanking epochs; ";
  logger << cnt_mask_set << " newly masked, " 
	 << cnt_mask_unset << " unmasked, " 
	 << cnt_unchanged << " unchanged\n";
  logger << " total of " << cnt_now_unmasked << " of " << epochs.size() << " retained\n";

}

// select all EPOCHs until we come across an EPOCH that does /not/ have the 'str' annotation
void timeline_t::select_epoch_until_isnot( const std::string & str )
{

  mask_set = true;
  
  const int ne = epochs.size();

  bool found = false;

  int cnt_mask_set = 0;
  int cnt_mask_unset = 0;
  int cnt_unchanged = 0;
  int cnt_now_unmasked = 0;

  for (int e=0;e<ne;e++)
    {      
      bool a = epoch_annotation( str , e );	  
      if ( ! a ) found = true;	  
      
      int mc = set_epoch_mask( e , found );
      if ( mc == +1 ) ++cnt_mask_set;
      else if ( mc == -1 ) ++cnt_mask_unset;
      else ++cnt_unchanged;
      
      if ( ! mask[e] ) ++cnt_now_unmasked;

    }

  logger << " based on " << str << " leading epochs; ";

  logger << cnt_mask_set << " newly masked, " 
	 << cnt_mask_unset << " unmasked, " 
	 << cnt_unchanged << " unchanged\n";
  logger << " total of " << cnt_now_unmasked << " of " << epochs.size() << " retained\n";

}



void timeline_t::annotate_epochs( const std::string & label , 
				  const std::string & annot_label , 
				  const std::set<std::string> & values )
{


  //
  // Take information from the annot_t class, and make a simple per-epoch annotation  
  //
  
  // this can be performed after a restructure, but (total) # of epochs must match the file exactly
  
  //
  // Get things that we want to display; from an EPOCH-ANNOT-EPOCHS command, 
  // all options are keys;  the key is the display name; the value is the 
  // mask-style option, e.g.   nrem2=SRO::Stage2Sleep[1] 
  //  label=annot_label[values]  
  // where values is a comma-delimited list
  
  //
  // What annotations are present? (i.e. already loaded)
  //
  
  //  std::vector<std::string> annots = annotations.names();
  

  //
  // Point to first epoch, but get the 'total' number of epochs (masked and unmasked), 
  // first_epoch() only returns the unmasked counts
  //
  
  first_epoch();
  
  int ne = num_total_epochs();
  
  //
  // Populate epoch-annotation vectors to the appropriate size
  //
  
  eannots[ label ].clear();
  

  //
  // for each each epoch 
  //
  
  while ( 1 ) 
    {
      
      //
      // Get next epoch
      //
      
      int e = next_epoch_ignoring_mask();      

      if ( e == -1 ) break;
      
      interval_t interval = epoch( e );
      
      //
      // Get annotations
      //
      
      annot_t * annot = annotations( annot_label );
      
      if ( annot == NULL ) continue;

      annot_map_t events = annot->extract( interval );
      
      // search for a matching value (at least one)
      
      annot_map_t::const_iterator ii = events.begin();
      while ( ii != events.end() )
	{	
	  
	  const instance_idx_t & instance_idx = ii->first;
	  const instance_t * instance = ii->second;

	  if ( values.find( instance_idx.id ) != values.end() )
	    {
	      eannots[ label ][ e ] = true;
	      break;
	    }	      
	  
	  ++ii;
	  
	}

    } // next epoch
}


uint64_t timeline_t::timepoint( int r , int s , int nsamples ) const
{

  std::map<int,uint64_t>::const_iterator rr = rec2tp.find(r);
  if ( rr == rec2tp.end() ) return 0;
  
  uint64_t x = s != 0 && nsamples != 0 
    ? edf->header.record_duration_tp * s / nsamples 
    : 0 ;

  return rr->second + x;
}


// uint64_t timeline_t::endpoint( int r ) const
// {
//   // note -- not using rec2tp_end here (which is probably redundant in any case)
//   std::map<int,uint64_t>::const_iterator rr = rec2tp.find(r);
//   if ( rr == rec2tp.end() ) return 0;  
//   return rr->second + edf->header.record_duration_tp - 1;
// }


void timeline_t::mask2annot( const std::string & path , const std::string & tag ) 
{

  if ( ! mask_set ) return;
  
  std::string path2 = path[ path.size() - 1 ] != globals::folder_delimiter 
    ? path + globals::folder_delimiter 
    : path ; 
  
  std::string filename = path2 + tag + "-" + edf->id + ".annot"; 
  
  annot_t * a = annotations.add( tag );
  a->description = tag + "-mask" ;  
  a->types[ "M" ] = globals::A_BOOL_T;
  
  const int ne = mask.size();
  
  for (int e=0;e<ne;e++)
    {
      if ( mask[e] )
	{
	  instance_t * instance = a->add( tag , epoch(e) );
	  instance->set( "M" , true );
	}
    }

  a->save( filename );

  // this will also retain the annotatiton 'tag', so it can be used
  // downstream by explicitly requesting the 'tag' annotation even if
  // the mask changes (i.e. rather than delete the annotation here)
  
}



void timeline_t::dumpmask()
{

  // no mask set: means all clear so display that
  //if ( ! mask_set ) return;

  first_epoch();
  
  logger << " dumping MASK\n";

  while ( 1 ) 
    {
      
      int e = next_epoch_ignoring_mask();      
      
      if ( e == -1 ) break;
      
      interval_t interval = epoch( e );

      // EPOCH_INTERVAL will already have been output by the EPOCH command
      writer.epoch( display_epoch( e ) );
      //      writer.var(   "INTERVAL" , "Epoch time start/stop" );
      writer.var(   "EPOCH_MASK" ,      "Is masked? (1=Y)" );
      //      writer.value( "INTERVAL" , interval.as_string() );
      writer.value( "EPOCH_MASK" , mask_set ? mask[e] : false );

    }

  writer.unepoch();
  
}


void timeline_t::set_epoch_mapping()
{

  bool has_mapping = has_epoch_mapping();

  first_epoch();
  
  // E1  1    1->1         1            1    1->1
  // E2  .                 2            .
  // E3  2    2->3         3            .
  // E4  .                 4            .
  // E5  3    3->5                      2    2->5
  // E6  4    4->6                      3    
  

 
  //
  // First mapping (i.e. not previously set)
  //

  if ( ! has_mapping )
    {

      clear_epoch_mapping();
      
      orig_epoch_size = num_total_epochs();

      int curr = 0; 
      
      while ( 1 ) 
	{      
	  int epoch = next_epoch_ignoring_mask();      
	  
	  if ( epoch == -1 ) break;	  
	  
	  if ( ! masked_epoch( epoch ) )
	    {
	      epoch_orig2curr[ epoch ] = curr;
	      epoch_curr2orig[ curr ] = epoch;
	      ++curr;
	    }
	}
    }
  else // otherwise, already has a mapping
    {
      
      std::map<int,int> copy_curr2orig = epoch_curr2orig;
      clear_epoch_mapping();      
      int curr = 0;       
      while ( 1 ) 
	{      
	  int epoch = next_epoch_ignoring_mask();      
	
	  if ( epoch == -1 ) break;

	  if ( ! masked_epoch( epoch ) )
	    {
	      int orig = copy_curr2orig[ epoch ];
	      epoch_orig2curr[ orig ] = curr;
	      epoch_curr2orig[ curr ] = orig;
	      ++curr;
	    }
	  
	}
    }
  
  
  // output -- note, always a 1-based epoch counting, 
  // not the internal 0-based representation
  
  // Note -- no longer output this, as we now make sure we always
  // output the original epoch (EDF-based) as well as the new 1..n
  // epoch code for any epoch-specific analysis
  
//   if ( false )
//     {
//       for (int e = 0 ; e < orig_epoch_size ; e++ ) 
// 	{
	  
// 	  std::cout << "EPOCH-MAP\t"
// 		    << edf->id << "\t"
// 		    << "[" << globals::current_tag << "]\t"
// 		    << e+1 << "\t";
	  
// 	  if ( epoch_orig2curr.find(e) == epoch_orig2curr.end() )
// 	    std::cout << ".\n";
// 	  else
// 	    std::cout << epoch_orig2curr[e] + 1 << "\n"; 
// 	}
//     }

}


void timeline_t::load_mask( const std::string & f , bool exclude )
{
  
  if ( ! epoched() ) 
    {
      int ne = set_epoch( globals::default_epoch_len , globals::default_epoch_len );
      logger << " set epochs to default " << globals::default_epoch_len << " seconds, " << ne << " epochs\n";
    }
     
  if ( ! Helper::fileExists( f ) ) Helper::halt( "could not find " + f );
  
  logger << " attaching mask file " << f << "\n";
  
  logger << " currently, mask mode set to: ";
  int mm = epoch_mask_mode();
  if ( mm == 0 ) logger << " mask (default)\n";
  else if ( mm == 1 ) logger << " unmask\n";
  else if ( mm == 2 ) logger << " force\n";

  
  // load
  std::ifstream FIN( f.c_str() , std::ios::in );

  int cnt_total = num_total_epochs();
  int cnt_mask0 = 0;
  int cnt_mask1 = 0;
  
  int e = 0;
  
  while ( ! FIN.eof() )
    {
      
      int m = 0;
      
      FIN >> m;

      if ( FIN.eof() ) continue;
      
      if ( ( exclude && m == 1 ) || ( (!exclude) && m == 0 ) )
	{
	  if ( ! masked(e) ) ++cnt_mask1;
	  set_epoch_mask( e );
	  ++cnt_mask0;	  
	}
      
      ++e;

      if ( e > cnt_total )
	{
	  logger << e << " masks read, for " << cnt_total << " existing epochs\n";
	  Helper::halt( "too many epochs specified in " + f );	
	}
    }
  
  
  FIN.close();

  logger << " processed " << e
	    << " lines, with "
	    << cnt_mask0 << " masked epochs\n";

  logger << " changed mask for " << cnt_mask1
	    << " of " << cnt_total << " epochs\n";

  return;

}




void timeline_t::load_interval_list_mask( const std::string & f , bool exclude )
{

  
  Helper::halt( "not supported" );

  // assume format  time1   time2    { meta-data ....  ignored }

  // if +time1 +time2  implies an offset from start of record
  // otherwise, assume implies real clocktime, based on edf.header.starttime
  // need to handle discontinuous EDFs here
  
  
  if ( ! Helper::fileExists( f ) ) Helper::halt( "could not find " + f );
  
  logger << " reading intervals to " << ( exclude ? " exclude" : "retain" ) << " from " << f << "\n";
  
  logger << " currently, mask mode set to: ";
  int mm = epoch_mask_mode();
  if      ( mm == 0 ) logger << " mask (default)\n";
  else if ( mm == 1 ) logger << " unmask\n";
  else if ( mm == 2 ) logger << " force\n";
  
  // load
  std::ifstream FIN( f.c_str() , std::ios::in );
  
  std::vector<interval_t> intervals;
  int cnt = 0;
  while ( ! FIN.eof() )
    {

      std::string line;

      std::getline( FIN , line , '\n' );
      
      std::vector<std::string> tok = Helper::parse( line , "\t" );
      
      if ( tok.size() == 0 ) 
	continue;
      
      if ( tok.size() < 2 ) 
	Helper::halt( "bad format in " + f + ", expecting at least 2 tab-delimited time fields" ) ;
      
      clocktime_t t1( tok[0] );
      clocktime_t t2( tok[1] );
      
      if ( ! t1.valid ) Helper::halt( "invalid HH:MM:SS timestring: " + tok[0] );
      if ( ! t2.valid ) Helper::halt( "invalid HH:MM:SS timestring: " + tok[1] );
      ++cnt;
    }
  
  FIN.close();

  logger << " processed " << cnt << " " << intervals.size() << " intervals\n";


  //
  // figure out start time of EDF... either from header, or from EDF itself, i.e. if it has been editted. 
  //

  
  //
  // Make sure that we have a time-track set
  //
  
  edf->add_continuous_time_track();
  
    
  return;
    
}


void timeline_t::apply_simple_epoch_mask( const std::set<std::string> & labels , const std::string & onelabel , bool include )
{
  
  // if 'ifnot', can only specify a single 
  if ( labels.size() > 1 && ! include ) Helper::halt( "can only specify a single mask for 'ifnot'");

  mask_set = true;
  
  const int ne = epochs.size();
 
  // Note: we do not clear the mask here, as we want to allow multiple
  // filters to be added on top of oneanther
  
  int cnt_mask_set = 0;
  int cnt_mask_unset = 0;
  int cnt_unchanged = 0;
  int cnt_now_unmasked = 0;
  int cnt_basic_match = 0;  // basic count of matches, whether changes mask or not
  

  for (int e=0;e<ne;e++)
    {

      bool matches = false;
      
      std::set<std::string>::const_iterator ii = labels.begin();
      while ( ii != labels.end() )
	{
	  if ( epoch_annotation( *ii , e ) ) { matches = true; break; }
	  ++ii;
	}
      
      // count basic matches
      if ( matches ) ++cnt_basic_match;
      
      // set new potential mask, depending on match_mode      
      bool new_mask = mask[e];
      
      if ( include ) 
	{
	  if      ( mask_mode == 0 ) new_mask = matches;   // mask-if
	  else if ( mask_mode == 1 ) new_mask = !matches;  // unmask-if
	  else if ( mask_mode == 2 ) new_mask = matches ;  // if
	}
      else
	{
	  if      ( mask_mode == 0 ) new_mask = !matches;  // mask-ifnot
	  else if ( mask_mode == 1 ) new_mask = matches;   // unmask-ifnot
	  else if ( mask_mode == 2 ) new_mask = ! matches; // ifnot
	}

      int mc = set_epoch_mask( e , new_mask );

      if      ( mc == +1 ) ++cnt_mask_set;
      else if ( mc == -1 ) ++cnt_mask_unset;
      else                 ++cnt_unchanged;
      
      if ( !mask[e] ) ++cnt_now_unmasked;
      
    }
  
  logger << " based on " << onelabel << " " << cnt_basic_match << " epochs match; ";

  logger << cnt_mask_set << " newly masked, "   
	 << cnt_mask_unset << " unmasked, " 
	 << cnt_unchanged << " unchanged\n";
  logger << " total of " << cnt_now_unmasked << " of " << epochs.size() << " retained\n";

  // mask, # epochs masked, # epochs unmasked, # unchanged, # total masked , # total epochs
  
  writer.level( onelabel , "EPOCH_MASK" );

  writer.var( "N_MATCHES"    , "Number of matching epochs" );
  writer.var( "N_MASK_SET"   , "Number of epochs newly masked" ); 
  writer.var( "N_MASK_UNSET" , "Number of epochs newly unmasked" );
  writer.var( "N_UNCHANGED"  , "Number of epochs unchanged by this mask" );
  writer.var( "N_RETAINED"   , "Number of epochs retained for analysis" );
  writer.var( "N_TOTAL"      , "Total number of epochs" );

  writer.value( "N_MATCHES"    , cnt_basic_match  );
  writer.value( "N_MASK_SET"   , cnt_mask_set     );
  writer.value( "N_MASK_UNSET" , cnt_mask_unset   );
  writer.value( "N_UNCHANGED"  , cnt_unchanged    );
  writer.value( "N_RETAINED"   , cnt_now_unmasked );
  writer.value( "N_TOTAL"      , epochs.size()    );

  writer.unlevel( "EPOCH_MASK" );

}



//
// Hypnogram functions
//

void hypnogram_t::construct( timeline_t * t , const bool verbose , const std::vector<std::string> & s )
{ 
  timeline = t;
  if ( s.size() != timeline->num_total_epochs() ) Helper::halt( "bad number of stages, " + Helper::int2str( (int)s.size() ) + " but expecting " + Helper::int2str( timeline->num_total_epochs() ) );    
  stages.resize( s.size() );
  for (int e=0;e<s.size();e++) stages[e] = globals::stage( s[e] );
  calc_stats( verbose );
} 

void hypnogram_t::construct( timeline_t * t , const bool verbose , const std::string sslabel ) 
{
  
  // point to 'parent' timeline
  timeline = t ;
  
  // get handle
  annot_t * annot = timeline->annotations( sslabel );
  if ( annot == NULL ) Helper::halt( "[" + sslabel + "] not set" );
  
  // set epoch-level annotations
  std::set<std::string> values;
  values.clear(); values.insert( "wake" );
  timeline->annotate_epochs(  globals::stage( WAKE ) , "SleepStage" , values );

  values.clear(); values.insert( "NREM1" );
  timeline->annotate_epochs(  globals::stage( NREM1  )  , "SleepStage" , values );

  values.clear(); values.insert( "NREM2" );
  timeline->annotate_epochs(  globals::stage( NREM2  )  , "SleepStage" , values );

  values.clear(); values.insert( "NREM3" );
  timeline->annotate_epochs(  globals::stage( NREM3  )  , "SleepStage" , values );

  values.clear(); values.insert( "NREM4" );
  timeline->annotate_epochs(  globals::stage( NREM4 )  , "SleepStage" , values );

  values.clear(); values.insert( "REM" );
  timeline->annotate_epochs(  globals::stage( REM ) , "SleepStage" , values );

  //
  // If we've masked the data, epoch count may not start at 0...
  // Although it probably doesn't make sense to use HYPNO then, 
  // we still want to be able to dump STAGE information easily w/ the
  // STAGE command...
  //
  

  // in VERBOSE (HYPNO) mode, we require the FULL epoch set
  
  if ( verbose ) 
    {
      if ( timeline->num_total_epochs() !=  timeline->num_epochs() ) 
	Helper::halt( "cannot run HYPNO on masked data" );
      
      int eprev = -1;
      timeline->first_epoch();
      while ( 1 ) 
	{	  
	  int e = timeline->next_epoch();      
	  if ( e == -1 ) break;
	  if ( eprev >= 0 && timeline->display_epoch( e ) - eprev != 1 ) 
	    Helper::halt( "cannot run HYPNO on masked data" );
	  eprev =  timeline->display_epoch( e ) ;
	}
    }

  const int ne = timeline->num_total_epochs();
  
  timeline->first_epoch();

  stages.clear();
  epoch_n.clear();

  //
  // Need to check how epoch annotations work after a RESTRUCTURE...
  //
  
  while ( 1 ) 
    {

      int e = timeline->next_epoch();
      
      if ( e == -1 ) break;
      
      int e2 = timeline->display_epoch(e) ;
      
      bool wake = timeline->epoch_annotation( "wake"  , e2 );
      bool n1   = timeline->epoch_annotation( "NREM1" , e2 );
      bool n2   = timeline->epoch_annotation( "NREM2" , e2 );
      bool n3   = timeline->epoch_annotation( "NREM3" , e2 );
      bool n4   = timeline->epoch_annotation( "NREM4" , e2 );
      bool rem  = timeline->epoch_annotation( "REM"   , e2 );
      


      bool other = ! ( wake || n1 || n2 || n3 || n4 || rem );
      bool conflict = ( (int)wake + (int)n1 + (int)n2 + (int)n3 + (int)n4 + (int)rem ) > 1;
      if ( conflict ) other = true;

      std::cerr << "ss " << wake << n1 << n2<<n3<<n4<<rem << other << "\n";
      
      if      ( conflict ) stages.push_back( UNSCORED );
      else if ( other ) stages.push_back( UNSCORED );
      else if ( wake ) stages.push_back( WAKE );
      else if ( n1 ) stages.push_back( NREM1 );
      else if ( n2 ) stages.push_back( NREM2 );
      else if ( n3 ) stages.push_back( NREM3 );
      else if ( n4 ) stages.push_back( NREM4 );
      else if ( rem ) stages.push_back( REM );
      else stages.push_back( UNSCORED );
      
      epoch_n.push_back( e2 );
      
    }

   calc_stats( verbose );
}   


void hypnogram_t::calc_stats( const bool verbose )
{

  //
  // epoch size (in minutes) and number
  //

  const double epoch_mins = timeline->epoch_length() / 60.0 ; 
  
  const int ne = stages.size();
  

  //
  // Recode any leading/trailing "?" as "L"
  //
  
  for (int e =0; e < ne ; e++)
    {
      if ( stages[e] == UNSCORED ) stages[e] = LIGHTS_ON;
      if ( stages[e] != UNSCORED && stages[e] != LIGHTS_ON ) break;
    }
  
  for (int e = ne - 1 ; e != 0 ; e--)
    {
      if ( stages[e] == UNSCORED ) stages[e] = LIGHTS_ON;
      if ( stages[e] != UNSCORED && stages[e] != LIGHTS_ON ) break;
    }

  
  //
  // Basic summary statistics per-individual/night
  //

  mins_wake = mins_n1 = mins_n2 = mins_n3 = mins_n4 = mins_rem = mins_other = 0;
  
  for (int e = 0 ; e < ne ; e++ )
    {
      if      ( stages[e] == WAKE  ) mins_wake += epoch_mins;
      else if ( stages[e] == NREM1 ) mins_n1 += epoch_mins;
      else if ( stages[e] == NREM2 ) mins_n2 += epoch_mins;
      else if ( stages[e] == NREM3 ) mins_n3 += epoch_mins;
      else if ( stages[e] == NREM4 ) mins_n4 += epoch_mins;
      else if ( stages[e] == REM   ) mins_rem += epoch_mins;
      else mins_other += epoch_mins;
    }
  
  final_wake_epoch = ne; // i.e. one past end
  for (int e = ne-1 ; e >= 0 ; e-- )
    if ( stages[e] != WAKE ) { final_wake_epoch = e+1; break; }
  
  int first_rem_epoch = ne;
  for (int e = 0 ; e < ne ; e++ )
    if ( stages[e] == REM ) { first_rem_epoch = e; break; }

  // lights out/on
  int lights_out_epoch = 0;
  for (int e=0;e<ne-1;e++) 
    if ( stages[e] == LIGHTS_ON ) 
      { ++lights_out_epoch; break; }
  
  int lights_on_epoch = ne; // by default, one past the end  
  for (int e=ne-1;e>0;e--) 
    if ( stages[e] == LIGHTS_ON ) 
      { --lights_on_epoch; break; }
  
  // requires non-missing SLEEP
  // persistent sleep defined as 10 mins
  int lps_required = 10.0 / (double)epoch_mins;
  int cnt_sleep = 0;

  bool found_first_sleep = false;
  first_sleep_epoch = ne;
  first_persistent_sleep_epoch = ne;

  for (int e = 0 ; e < ne ; e++ )
    {
      if ( is_sleep( stages[e] ) )
	{
	  if ( ! found_first_sleep )
	    {
	      first_sleep_epoch = e; 
	      found_first_sleep = true;
	    }
	  
	  // LPS start?
	  bool lps = true;  
	  for (int e2=e; e2 < e + lps_required; e2++ )
	    {	      
	      if ( e2 >= ne ) { lps=false; break; }
	      if ( ! is_sleep( stages[e2] ) ) 
		{ lps = false; break; }
	    }
	  
	  if ( lps )
	    {
	      first_persistent_sleep_epoch = e;
	      break;
	    }
	}
    }
  
  // last epoch of sleep
  int last_sleep_epoch = 0;
  for (int e = ne - 1; e != 0 ; e-- )
    if ( is_sleep( stages[e] ) ) { last_sleep_epoch = e; break; } 
  
  // total time in bed
  TIB = ne * epoch_mins;
  
  // total recording time (i.e. only from lights out, lights on)
  int TRT_total_epochs = lights_on_epoch - lights_out_epoch + 1;
  TRT =  TRT_total_epochs * epoch_mins;
  
  // total wake time (ignores pre lights out, post lights off)
  TWT = mins_wake;
  
  // final wake time 
  FWT = ( lights_on_epoch - final_wake_epoch ) * epoch_mins; 
  
  // REM latency
  rem_lat_mins = ( first_rem_epoch - first_sleep_epoch ) * epoch_mins;
  
  // Total sleep time (includes 'other')
  TST = TIB - TWT; 
  
  // sleep latency
  slp_lat = ( first_sleep_epoch - lights_out_epoch ) * epoch_mins;
  
  // latency to persistent sleep
  per_slp_lat = ( first_persistent_sleep_epoch - lights_out_epoch ) * epoch_mins;

  // Sleep period time
  SPT = TRT - slp_lat;

  // WASO (ignores leading and also trailing wake)
  WASO = TWT - slp_lat - FWT;

  // sleep efficiency
  slp_eff_pct = ( TST / TRT ) * 100;
  
  // sleep maintainence/efficiency 2 (denom is from initial sleep to final sleep)
  slp_eff2_pct = ( TST / ( epoch_mins * ( last_sleep_epoch - first_sleep_epoch + 1 ) ) ) * 100 ; 
  
  // sleep maintainence
  slp_main_pct = ( TST / SPT ) * 100;

  if ( TST > 0 ) 
    {
      pct_n1  = mins_n1  / TST;
      pct_n2  = mins_n2  / TST;
      pct_n3  = mins_n3  / TST;
      pct_n4  = mins_n4  / TST;
      pct_rem = mins_rem / TST;
    }
  else
    {
      pct_n1  = 0;
      pct_n2  = 0;
      pct_n3  = 0;
      pct_n4  = 0;
      pct_rem = 0;
    }

   


  //
  // Sleep cycles : based on modified Floyd & Feinberg rules
  //
  

  // Thresolds

  // Minimum duration for a NREM period
  const double def_min_nrem_duration_mins = 15;

  // Minimum duration for REM period (cycle 2 and after)
  const double def_min_rem_duration_mins = 5;

  // Maximum duration of NREM/W allowed within a single REM episode
  const double def_rem_period_interuption_mins  = 15;

  // If skipping REM period, minimum W/N1 to terminate a NREM period
  const double def_terminating_waso_duration_mins = 15;

  // Persistent sleep is defined after 10 minutes of sleep
  const double def_persistent_sleep_mins = 10;


  // convert the above to epoch counts
  const int def_persistent_sleep_epochs = def_persistent_sleep_mins / (double)epoch_mins ; 
  const int def_rem_period_interuption_epochs = def_rem_period_interuption_mins / (double)epoch_mins;
  const int def_min_nrem_duration_epochs = def_min_nrem_duration_mins / (double)epoch_mins;
  const int def_terminating_waso_duration_epochs = def_terminating_waso_duration_mins / (double)epoch_mins;
  const int def_min_rem_duration_epochs = def_min_rem_duration_mins / (double)epoch_mins;


  //
  // 0)  Handle movement: 'impute' with the following epoch value
  //     ignored for now... 
  

  //  
  // 1)  Find periods of 'persistent sleep' (default 10mins prior sleep)
  //
  
  TpST = 0;

  std::vector<std::string> persistent_sleep( ne , "" );
  for (int e=0;e<ne;e++)
    {
      
      if ( stages[ e ] == WAKE || stages[ e ] == LIGHTS_ON || stages[e] == UNSCORED )
	{
	  persistent_sleep[e] = "W";
	  continue;
	}
      
      // otherwise, assume all other annotations are consistent with sleep
      bool okay = true;
      int ec = e - def_persistent_sleep_epochs;
      
      while ( okay )
	{
	  if ( ec < 0 ) { okay = false; break; }	  
	  if ( stages[ ec ] == WAKE || stages[ ec ] == LIGHTS_ON ) { okay = false; break; }  	  
	  if ( ++ec == e ) break;
	}

      if ( okay ) 
	{
	  persistent_sleep[e] = "S"; 
	  TpST += epoch_mins;
	}
      else persistent_sleep[e] = "W";
    }
  
  
  
  //
  // 2) Find sleep onset 
  //
  
  std::vector<std::string> sleep_onset( ne , "" );
  bool found_sleep = false;
  for (int e = 0; e < ne ; e++ ) 
    {
      if ( is_sleep( stages[e] ) ) found_sleep = true;
      sleep_onset[ e ] = found_sleep ? "S" : "W";
    }

  // and likewise, same from end of sleep
  found_sleep = false;
  for (int e = ne-1; e >= 0 ; e-- ) 
    {
      if ( is_sleep( stages[e] ) ) break;
      sleep_onset[ e ] = "W";
    }
 


  //
  // 3) Cumulative count of sleep
  //

  // after LightsOn, defined as '-1', if this matters 
  std::vector<int> sleep_count( ne , 0 );
  int cum_sleep = 0;
  for (int e = 0 ; e < ne ; e++)
    {
      if ( persistent_sleep[e] == "S" ) ++cum_sleep;
      if ( stages[e] == LIGHTS_ON && cum_sleep > 0 ) sleep_count[e] = -1;
      else sleep_count[e] = cum_sleep;
    }
  

  //
  // 4) Sleep state
  //
  
  std::vector<std::string> sleep_state( ne , "" );
  for (int e = 0 ; e < ne; e++)
    {
      if ( stages[e] == LIGHTS_ON && sleep_count[e] == 0 ) sleep_state[e] = "Prior";
      else if ( sleep_count[e] == 0 ) sleep_state[e] = "LPS"; // latency before persistent sleep
      else if ( sleep_count[e] == 1 ) sleep_state[e] = "LPO"; // onset of persistent sleep
      else 
	{
	  if ( sleep_count[e] > 1 ) sleep_state[e] = "SPT"; // sleep period time
	  else sleep_state[e] = "After"; 
	}
    }
  

  //
  // 5) Final wake ('WATA')
  //

  std::vector<bool> wata( ne , false );
  for (int e = ne-1 ; e>= 0 ; e-- )
    {
      if ( is_sleep( stages[e] ) ) break;
      if ( is_wake( stages[e] ) ) wata[e] = true;
    }

  
  //
  // 6) Sleep period/cycle 
  //
  
  std::vector<std::string> sleep_period( ne , "" );
  std::vector<bool> cycle_ending_waso( ne , false );
    
  for (int e = 0 ; e < ne ; e++ ) 
    {

      if ( sleep_onset[e] == "W" ) continue;

      bool previous_epoch_defined = e == 0 ? false : sleep_period[ e - 1 ] != "" ;

      if ( is_rem( stages[e] ) && previous_epoch_defined )
	{
	  // continues a new REM stage
	  sleep_period[e] = "REM";	  
	}
      else
	{

	  // check subsequent 15 mins
	  bool has_another_rem = false;
	  const int elimit = ne-1 < e + def_rem_period_interuption_epochs - 1 ? 
	    ne-1 : e + def_rem_period_interuption_epochs - 1 ;
	  for (int e2=e;e2 <= elimit; e2++)
	    {
	      if ( is_rem( stages[e2] ) ) { has_another_rem = true; break; }
	    }


	  // from start (i.e. including this one), next 15 mins has to have at least 1 other REM 
	  if ( ( e>0 && sleep_period[e-1] == "REM" ) && has_another_rem )
	    sleep_period[e] = "REM";
	  else
	    {

	      // else, if previously REM and cycle ended last epoch OR W/N1
	      if ( ( ( e>0 && sleep_period[e-1] == "REM" ) || ( e>0 && cycle_ending_waso[e-1] ) ) // previous REM, or cycle end
		   && ( is_wake( stages[e] ) || is_nrem1( stages[e] ) ) ) // AND W/N1 this epoch
		sleep_period[e] = ""; 
	      else
		{
		  
		  // note, uses potentially different duration from above, but obviosuly all this code 
		  // can be streamlined.
		  
		  bool has_another_rem = false;
		  const int elimit = ne-1 < e + def_min_nrem_duration_epochs - 1 ? 
		    ne-1 : e + def_min_nrem_duration_epochs - 1 ;
		  for (int e2=e;e2 <= elimit; e2++)
		    if ( is_rem( stages[e2] ) ) { has_another_rem = true; break; }
		  
		  // else, if previously in NREM, continue; (Q6!="")
		  if ( ( e==0 || sleep_period[e-1]=="" ) &&
		       ( is_wake( stages[e] ) || is_nrem1( stages[e] )  || has_another_rem ) )
		    sleep_period[e] = ""; 
		  else	// else, initiate a new NREM cycle, only if no REM wihtin the next 15mins
		    {
		      sleep_period[e] = "NREM";		      
		    }
		}
	    }
	}
      

      //
      // Cycle-ending WASO
      //

      //  put a cycle end IF currently NREM and in the next 15mins (including this epoch),
      //  if this is NREM but the next 15 mins has only W/N1, then end cycle here
      //  ELSE, if past epoch was end flag, and is WAKE, also flag (i.e. up to end of WASO W)
      
      
      bool no_near_sleep = true;
      const int elimit = ne-1 < e + def_terminating_waso_duration_epochs - 1 ? ne-1 : e+ def_terminating_waso_duration_epochs - 1 ;
      for (int e2=e;e2<=elimit;e2++)
	if ( is_nrem23( stages[e2] ) || is_rem( stages[e2] ) ) 
	  { no_near_sleep = false; break; } 
      
      if ( sleep_period[e] == "NREM" && no_near_sleep )
	{
	  cycle_ending_waso[e] = true;
	}
      else
	{
	  if ( e > 0 && cycle_ending_waso[ e - 1 ] && is_wake( stages[e] ) ) 
	    cycle_ending_waso[e] = true;
	}
      
    }

  
  //
  // Cycle type, number
  //

  sleep_code.resize( ne , 0 ); // 0, 1, 5 for W, NREM, REM
  sleep_cycle_number.resize( ne , 0 );

  
  // get first REM/cycle-ending epoch
  
  int first_sleep_period_rem = 99999;
  int first_cycle_ending_waso = 99999;
  
  for (int e=0; e<ne; e++)
    if ( sleep_period[e] == "REM" ) { first_sleep_period_rem = e; break; }
    
  for (int e=0; e<ne; e++)
    if ( cycle_ending_waso[e] ) { first_cycle_ending_waso = e; break; }

  for (int e=0; e<ne; e++)
    {

      // skip if a cycle-ending WASO
      if ( cycle_ending_waso[e] ) continue;

      if ( sleep_period[e] == "NREM" )
	{
	  sleep_code[e] = 1;
	}
      else
	{
	  
	  if ( sleep_period[e] == "REM" )
	    {
	      if ( e > 0 && sleep_period[e-1] == "NREM" )
		{
		  
		  // first cycle?
		  if ( e <= first_sleep_period_rem 
		       && e <= first_cycle_ending_waso )
		    sleep_code[e] = 5;
		  else
		    {
		      // check ahead... requires at least 'def_min_rem_duration_epochs' of REM
		      // nb. not sure why this would ever be ">=" the # of epochs... but keep as is for now
		      
		      int count_rem = 0;
		      const int elimit = ne-1 < e + def_min_rem_duration_epochs - 1 ? ne-1 : e+ def_min_rem_duration_epochs -1 ;
		      for (int e2=e; e2<=elimit; e2++)
			if ( sleep_period[e2] == "REM" ) ++count_rem;
		      if ( count_rem >= def_min_rem_duration_epochs ) 
			sleep_code[e] = 5;
		      else
			sleep_code[e] = 1;
		    }
		}
	      else
		{
		  if ( e>0 && sleep_period[e-1] == "REM" && sleep_code[e-1] == 5 ) 
		    sleep_code[e] = 5;
		  else 
		    sleep_code[e] = 1;
		}             
	    }
	  else
	    {
	      if ( e > 0 && sleep_period[e-1] == "REM" && sleep_code[e-1] == 1 )
		sleep_code[e] = 1; 
	      else
		{
		  if ( wata[e] ) 
		    sleep_code[e] = 0;
		  else
		    {
		      if ( sleep_period[e] == "" && e>0 && sleep_code[e-1] == 1 && ! cycle_ending_waso[e] )
			 sleep_code[e] = 1; 
		      else
			sleep_code[e] = 0;
		    }
		  
		}
	      
	    }
	  
	}
      
      // next epoch
    }

        
  //
  // Define cycles
  //

  
  int cnt_cycle = 0;

  for (int e=0; e< ne; e++)
    {
      if ( sleep_code[e] == 0 ) 
	sleep_cycle_number[e] = 0;
      else
	{
	  	  
	  const int previous_code = e == 0 ? 0 : sleep_code[e-1];

	  // change in cycle?
	  // start of a new NREM?
	  if ( sleep_code[e] - previous_code == 1  // into NREM(1) from WASO/N1(0) 
	       || previous_code - sleep_code[e] == 4 ) // from REM(5) to NREM(1)
	    {

	      // requires NREM (15mins) of 
	      
	      // find next REM and WASO epoch [ 'sleep_code' ]
	      // count eppochs (F) [ stages] 
	      
	      int elimit = ne - 1;
	      for (int e2 = e ; e2 < ne ; e2++)
		if ( sleep_code[e2] == 0 || sleep_code[e2] == 5 ) 
		  { elimit = e2; break; } 
	      
	      int cnt_nrem = 0;
	      for (int e2=e; e2 <= elimit; e2++)
		if ( is_nrem( stages[e2] ) ) ++cnt_nrem;
	      
	      // enough NREM for a new cycle?
	      if ( cnt_nrem >= def_min_nrem_duration_epochs ) 
		sleep_cycle_number[e] = ++cnt_cycle;
	      
	    }
	  else
	    sleep_cycle_number[e] = e == 0 ? 0 : sleep_cycle_number[e-1];
	  
	}
    }


  
  //
  // Get cycle/period statistics
  //

  // Count number of cyles
  
  num_nremc = 0;
  nremc_mean_duration = 0;

  std::map<int,int> cmin;
  std::map<int,int> cmax;
  std::map<int,int> counts_rem;
  std::map<int,int> counts_nrem;
  std::map<int,int> counts_other;

  for (int e=0;e<ne;e++)
    {
      const int & sn = sleep_cycle_number[e];
      if ( sn == 0 ) continue;
      if ( sn > num_nremc ) num_nremc = sn;
      if ( cmin.find( sn ) == cmin.end() ) cmin[ sn ] = cmax[sn] = e;
      cmax[sn] = e; // track max
      if ( is_rem( stages[e] ) ) counts_rem[sn]++;
      else if ( is_nrem( stages[e] ) ) counts_nrem[sn]++;
      else counts_other[sn]++;
    }


  
  std::map<int,int>::iterator ii = cmin.begin();
  while ( ii != cmin.end())
    {
      const int & sn = ii->first ;
      
      // total cycle duration
      double dur = cmax[ sn ] - ii->second + 1;
      double dur_mins = dur * epoch_mins ; 
      
      nremc_mean_duration += dur_mins;
 
      nremc_duration[ sn ] = ( counts_rem[sn] + counts_nrem[sn] + counts_other[sn] ) * epoch_mins ; 
      nremc_nrem_duration[ sn ] = counts_nrem[sn] * epoch_mins ; 
      nremc_rem_duration[ sn ] = counts_rem[sn]  * epoch_mins ; 
      
      nremc_start_epoch[ sn ] = ii->second + 1 ;  // output 1-based coding
            
      ++ii;
    }

  if ( num_nremc > 0 ) nremc_mean_duration /= (double)num_nremc;

  // cycle positions
  cycle_pos_relative.resize( ne , -1 );
  cycle_pos_absolute.resize( ne , -1 );
  for (int e=0; e<ne; e++)
    {
      const int & sn = sleep_cycle_number[e];
      if ( sn == 0 ) continue;      
      int cycle_start = cmin[sn];
      
      // position within each cycle.
      cycle_pos_absolute[e] = ( e - cycle_start ) * epoch_mins ; 
      cycle_pos_relative[e] = cycle_pos_absolute[e] / (double)nremc_duration[sn];
    }


  // after the fact, track epoch-level stats
  in_persistent_sleep.resize( ne , false );
  for (int e=0; e<ne; e++)
    if ( persistent_sleep[e] == "S" ) in_persistent_sleep[e] = true;


  
  //
  // Sleep cycle definitions
  //
  
  // Based primarily on Feinberg & Floyd (1979) definitions
  // with some modifications
  
  // Cycle = "NREM phase of at least X mins terminated by end of REM phase OR a preset duration of a bout of wake/N1
  
  // "Wake bout" = periods of wake of any duration after latency to persistent sleep onset
  
  // If REM phase is skipped (i.e. cycle is terminated by a wake bout
  // of a modifiable duration in minutes), bout of wake/N1 is excluded
  // from cycle new NREM phase onset is with first stage of sleep;
  // wake between cycles is excluded from a cycle and is denoted cycle
  // 0, but wake within a cycle is still included


  // REM phase is period of REM sleep of any duration terminated by a
  // preset contiguous duration of any other stage (wake or NREM).
  // REM phase may contain NREM or wake as long as it is less than the
  // preset duration of any combination of NREM or wake.  the user may
  // define whether stage N1 can define the onset of a new NREM phase.
 
   

  //
  // Ascending/descending N2 
  //
  
  // Assign each N2 epoch a value to describe its relative position in
  // the hypnogram: consider 'k' epochs before and after, calculating 
  // the average of non-N2 epochs,  either W/R/N1 or N3/N4  

  // Calculate a single score: 
  // Defined only for N2 epochs
  //   Left epochs    +1  N3              Right epochs    -1  N3
  //                  -1  N1/W/R                          +1  N1/W/R
  
  const int n2_ascdesc_k = 10;  // 5 minutes
  
  // Can select extreme epochs as 'ascending' and 'descending' based on this score 
  // (e.g. >+0.5 and < -0.5)
  
  n2_ascdesc.resize( ne , 0 );
  
  for (int e=0;e<ne;e++)
    {
      
      if ( stages[e] != NREM2 ) continue;
      
      double left_wgt = 0;
      int left_n = 0; 
      int k = e-1;

      while ( k >= 0 )
	{

	  if ( stages[k] == NREM3 || stages[k] == NREM4 )
	    {
	      left_wgt += 1;
	      ++left_n;
	    }
	  
	  if ( stages[k] == NREM1 || stages[k] == REM || stages[k] == WAKE )
	    {
	      left_wgt += -1;
	      ++left_n;
	    }

	  // counted enough?
	  if ( left_n > n2_ascdesc_k ) break;
	  
	  // next left epoch
	  --k;
	}

      //
      // Right-most
      //

      double right_wgt = 0;
      int right_n = 0; 
      k = e+1;

      while ( k < ne )
	{

	  if ( stages[k] == NREM3 || stages[k] == NREM4 )
	    {
	      right_wgt += -1;
	      ++right_n;
	    }
	  
	  if ( stages[k] == NREM1 || stages[k] == REM || stages[k] == WAKE )
	    {
	      right_wgt += +1;
	      ++right_n;
	    }
	  
	  // counted enough?
	  if ( right_n > n2_ascdesc_k ) break;
	  
	  // next right epoch
	  ++k;
	}
      
      if ( left_n  > 0 ) left_wgt /= (double)left_n;
      if ( right_n > 0 ) right_wgt /= (double)right_n;

      // simple average of left/right averages
      // if no data, wgt will be 0, which is fine
      n2_ascdesc[e] = ( left_wgt + right_wgt ) / 2.0;

    }
  

  
  //
  // Flanking epochs
  //

  is_waso.resize( ne , false );
  for (int e = 0 ; e < ne ; e++)
    {
      if ( stages[e] == WAKE 
	   && e > first_sleep_epoch 
	   && e < final_wake_epoch ) is_waso[e] = true;
    }

  flanking.resize( ne , 0 );
  nearest_wake.resize( ne , 0 );

  nrem2rem.resize( ne , 0 ); nrem2rem_total.resize( ne , 0 );
  nrem2wake.resize( ne , 0 ); nrem2wake_total.resize( ne , 0 );

  for (int e = 0 ; e < ne ; e++)
    {
      
      //
      // calculate the number of similar epochs 
      // (FLANKING_SIM)
      //
      
      int sim = 0;  
      
      for (int j=1;j<ne;j++)
	{
	  const int eleft  = e - j;
	  const int eright = e + j;
	  // too much
	  if ( eleft < 0 || eright >= ne ) { sim = j-1; break; }
	  if ( stages[eleft] != stages[e] || stages[eright] != stages[e] ) { sim = j-1; break; }	  
	}

      int nw = 0;

      if ( stages[e] != WAKE )
	{
	  for (int j=1;j<ne;j++)
	    {
	      const int eleft  = e - j;
	      const int eright = e + j;
	      // too much
	      if ( eleft < 0 || eright >= ne ) { nw = j; break; }
	      if ( stages[eleft] == WAKE || stages[eright] == WAKE ) { nw = j; break; }	  
	    }
	}

      flanking[e] = sim;
      nearest_wake[e] = nw;

      //
      // transitions FROM N2?
      //
            
      if ( stages[e] == NREM2 )
	{
	  
	  // n2 to rem
	  int ei = 1;
	  while ( 1 ) 
	    {
	      if ( e + ei == ne ) { ei=0; break; }
	      if ( stages[ e + ei ] == NREM2 ) { ++ei; continue; }
	      if ( stages[ e + ei ] == REM   ) break;
	      ei = 0; break;
	    }
	  nrem2rem[e] = ei;
	  
	  // n2 to wake
	  ei = 1;
	  while ( 1 )
	    {
	      if ( e + ei == ne ) { ei=0; break; }
	      if ( stages[ e + ei ] == NREM2 ) { ++ei; continue; }
	      if ( stages[ e + ei ] == WAKE   ) break;
	      ei = 0; break;
	    }
	  nrem2wake[e] = ei;
	}
    
    } // next epoch
  
  // now figure out the _total values 
  // i.e. move forward and copy largest number until we hit 0      
  int e_rem  = nrem2rem[0];
  int e_wake = nrem2wake[0]; 
  
  for (int e = 1 ; e < ne ; e++ )
    {
      if ( nrem2rem[e] == 0 ) e_rem = 0;
      else if ( nrem2rem[e] > e_rem ) e_rem = nrem2rem[e];
      nrem2rem_total[e] = e_rem;

      if ( nrem2wake[e] == 0 ) e_wake = 0;
      else if ( nrem2wake[e] > e_wake ) e_wake = nrem2wake[e];
      nrem2wake_total[e] = e_wake;
    }
    
  //
  // Clocktime-based measures
  //
  
  clocktime_t starttime( timeline->edf->header.starttime );
  if ( ! starttime.valid ) 
    {
      clock_lights_out.valid = clock_sleep_onset.valid 
	= clock_sleep_midpoint.valid = clock_wake_time.valid 
	= clock_lights_on.valid = false;
    }
  else
    {
      clock_lights_out     = starttime;
        
      double epoch_hrs = epoch_mins / 60.0;

      clock_sleep_onset    = starttime;
      clock_sleep_onset.advance( epoch_hrs * first_sleep_epoch );
      
      clock_wake_time      = starttime;
      clock_wake_time.advance( epoch_hrs * final_wake_epoch );
      
      clock_lights_on      = starttime;
      clock_lights_on.advance( epoch_hrs * ne );
    
      clock_sleep_midpoint.midpoint( clock_sleep_onset , clock_wake_time );      
      
    }


  
   

}

void hypnogram_t::output( const bool verbose )
{

  // currently, this routine is hard-coded to assume 30-second epochs,
  // so for now flag if this is not the case (we can fix downstream)


  if ( ! Helper::similar( timeline->epoch_length() , 30 , 0.001 ) ) 
    Helper::halt( "requires 30-second epochs to be set currently" );

  
  //
  // Per individual level output (VERBOSE MODE ONLY)
  //

  if ( verbose )
    {
      writer.var( "LIGHTS_OUT"  , "Lights out time [0,24)" );
      writer.var( "SLEEP_ONSET" , "Sleep onset time [0,24)" );
      writer.var( "SLEEP_MIDPOINT" , "Sleep mid-point time [0,24)" );
      writer.var( "FINAL_WAKE" , "Final wake time [0,24)" );
      writer.var( "LIGHTS_ON" , "Lights on time [0,24)" );
      
      writer.var( "NREMC" , "Number of NREM cycles" );
      writer.var( "NREMC_MINS" , "Average NREM cycle duration (mins)" );
      
      writer.var( "TIB" , "Time in Bed (hours): LIGHTS_OUT --> LIGHTS_ON" );
      writer.var( "TST" , "Total Sleep Time (hours): SLEEP_ONSET --> FINAL_WAKE" );
      writer.var( "TPST" , "Total persistent Sleep Time (hours): PERSISTENT_SLEEP_ONSET --> FINAL_WAKE" );
      
      writer.var( "TWT" , "Total Wake Time (hours): all WAKE" );
      writer.var( "WASO" , "Wake After Sleep Onset (hours)" );
      
      writer.var( "SLP_LAT" , "Sleep latency" );
      writer.var( "PER_SLP_LAT" , "Persistent sleep latency" );
      
      writer.var( "SLP_EFF" , "Sleep efficiency: LIGHTS_OUT --> LIGHTS_ON" );
      writer.var( "SLP_MAIN_EFF" , "Sleep maintainence efficiency" );
      writer.var( "SLP_EFF2" , "Sleep efficiency: SLEEP_ONSET --> FINAL_WAKE" );

      writer.var( "REM_LAT" , "REM latency (from SLEEP_ONSET)" );
      
      writer.var( "PCT_N1" , "Proportion of sleep that is N1" );
      writer.var( "PCT_N2" , "Proportion of sleep that is N2" );
      writer.var( "PCT_N3" , "Proportion of sleep that is N3" );
      writer.var( "PCT_N4" , "Proportion of sleep that is N4" );
      writer.var( "PCT_REM" , "Proportion of sleep that is REM" );
      
      writer.var( "MINS_N1" , "Proportion of sleep that is N1" );
      writer.var( "MINS_N2" , "Proportion of sleep that is N2" );
      writer.var( "MINS_N3" , "Proportion of sleep that is N3" );
      writer.var( "MINS_N4" , "Proportion of sleep that is N4" );
      writer.var( "MINS_REM" , "Proportion of sleep that is REM" );
      
      // values
      writer.value(  "LIGHTS_OUT" , clock_lights_out.as_numeric_string() );
      writer.value(  "SLEEP_ONSET" , clock_sleep_onset.as_numeric_string() );
      writer.value(  "SLEEP_MIDPOINT" , clock_sleep_midpoint.as_numeric_string() );
      writer.value(  "FINAL_WAKE" , clock_wake_time.as_numeric_string() );
      writer.value(  "LIGHTS_ON" , clock_lights_on.as_numeric_string() );
      
      writer.value(  "NREMC" , num_nremc );
      writer.value(  "NREMC_MINS" , nremc_mean_duration );
      
      writer.value( "TIB" , TIB );
      writer.value( "TST" , TST );
      writer.value( "TPST" , TpST );
      writer.value( "TWT" , TWT );
      writer.value( "WASO" , WASO );
      
      writer.value( "SLP_LAT" , slp_lat );
      writer.value( "PER_SLP_LAT" , per_slp_lat );
      
      writer.value( "SLP_EFF" , slp_eff_pct );
      writer.value( "SLP_MAIN_EFF" , slp_main_pct );
      writer.value( "SLP_EFF2" , slp_eff2_pct );
      
      if ( mins_rem > 0 )
	writer.value( "REM_LAT" , rem_lat_mins );
      
      writer.value( "PCT_N1" , pct_n1 );
      writer.value( "PCT_N2" , pct_n2 );
      writer.value( "PCT_N3" , pct_n3 );
      writer.value( "PCT_N4" , pct_n4 );
      writer.value( "PCT_REM" , pct_rem);
      
      writer.value( "MINS_N1" , mins_n1 );
      writer.value( "MINS_N2" , mins_n2 );
      writer.value( "MINS_N3" , mins_n3 );
      writer.value( "MINS_N4" , mins_n4 );
      writer.value( "MINS_REM" , mins_rem);

    }
  
  //
  // Cycle-specific output (verbose mode only)
  //

  if ( verbose ) 
    {
      
      writer.var( "NREMC_START" , "NREM cycle start epoch" );
      writer.var( "NREMC_NREM_MINS" , "NREM cycle NREM duration (mins)" );
      writer.var( "NREMC_REM_MINS" , "NREM cycle REM duration (mins)" );
      writer.var( "NREMC_OTHER_MINS" , "NREM cycle other duration (mins)" );
      writer.var( "NREMC_MINS" , "NREM cycle total duration (mins)" );
      
      std::map<int,double>::iterator cc = nremc_duration.begin();
      while ( cc != nremc_duration.end() )
	{
	  writer.level( cc->first , globals::cycle_strat );
	  
	  writer.value(  "NREMC_START" , nremc_start_epoch[ cc->first ] );
	  writer.value(  "NREMC_NREM_MINS" , nremc_nrem_duration[ cc->first ] );
	  writer.value(  "NREMC_REM_MINS" , nremc_rem_duration[ cc->first ] );
	  writer.value(  "NREMC_OTHER_MINS" , cc->second - nremc_nrem_duration[ cc->first ] - nremc_rem_duration[ cc->first ] );
	  writer.value( "NREMC_MINS" , cc->second );
	  
	  ++cc;
	}
      
      writer.unlevel( globals::cycle_strat );


    }

  
  //
  // Per epoch level output
  //


  // stage information and time only in non-verbose mode

  std::map<sleep_stage_t,int> stagen;
  stagen[ WAKE ] = 1;
  stagen[ REM ] = 0;
  stagen[ NREM1 ] = -1;
  stagen[ NREM2 ] = -2;
  stagen[ NREM3 ] = -3;
  stagen[ NREM4 ] = -4;
  stagen[ UNSCORED ] = 2;
  stagen[ UNKNOWN ] = 2; // this should not happen
  stagen[ MOVEMENT ] = 2;
  stagen[ ARTIFACT ] = 2;
  stagen[ LIGHTS_ON ] = 2;

  writer.var( "MINS" , "Elapsed time since start of recording (minutes)" );
  writer.var( "CLOCK_TIME" , "Clock time (hh:mm:ss)" );

  writer.var( "STAGE" , "Sleep stage, string label" );
  writer.var( "STAGE_N" , "Sleep stage, numeric encoding" );

  // epoch size (in minutes)
  const double epoch_mins = timeline->epoch_length() / 60.0 ; 
  const int ne = timeline->num_epochs();

  clocktime_t epoch_time( clock_lights_out );

  clocktime_t epoch_duration( "00:00:30" );
  
  std::cerr << "ne2 = " << ne << "\n";

  // output
  for (int e=0;e<ne;e++)
    {
      
      // epoch-level stratification
      writer.epoch( epoch_n[e] );
      
      writer.value( "MINS" ,  epoch_n[e] * epoch_mins );
      writer.value( "CLOCK_TIME" , epoch_time.as_string() );      
      if ( verbose ) 
	writer.value( "CLOCK_HOURS" ,  epoch_time.as_numeric_string() );
      
      // next epoch...
      epoch_time.advance( epoch_duration );           
      
      // stages
      writer.value( "STAGE" , globals::stage( stages[e] ) );
      writer.value( "STAGE_N" , stagen[ stages[e] ] );
    }

  writer.unepoch();

  
  //
  // ... otherwise, the rest of this function is verbose mode only
  //
  
  if ( ! verbose ) return;


  
  // Outputs
  // Per epoch, we have
  //   a) stage (done above)
  //   b) elapsed time
  //   c) elapsed sleep
  //   d) period number
  //   e) N2 measure of direction
  
  
  double elapsed_n1 = 0 , elapsed_n2 = 0 , elapsed_n34 = 0 , elapsed_rem = 0;
  double elapsed_sleep = 0 , elapsed_wake = 0 , elapsed_waso = 0 ;

  
  // header

  writer.var( "CLOCK_HOURS" , "Clock time [0,24) hours" );
        
  writer.var( "E_WAKE" , "Elapsed wake (mins)" );
  writer.var( "E_WASO" , "Elapsed WASO (mins)" );
  writer.var( "E_SLEEP" , "Elapsed sleep (mins)" );

  writer.var( "E_N1" , "Elapsed N1 (mins)" );
  writer.var( "E_N2" , "Elapsed N2 (mins)" );
  writer.var( "E_N3" , "Elapsed N3 (mins)" );
  writer.var( "E_REM" , "Elapsed REM (mins)" );

  writer.var( "PCT_E_SLEEP" , "Elapsed sleep (percent of all sleep)" );
  writer.var( "PCT_E_N1" , "Elapsed N1 (percent of all N1)" );
  writer.var( "PCT_E_N2" , "Elapsed N2 (percent of all N2)" );
  writer.var( "PCT_E_N3" , "Elapsed N3 (percent of all N3)" );
  writer.var( "PCT_E_REM" , "Elapsed REM (percent of all REM)" );

  writer.var( "PERSISTENT_SLEEP" , "Persistent sleep yes/no? (1=Y)" );

  writer.var( "CYCLE" , "NREMC number" );
  writer.var( "PERIOD" , "NREMC period (NREM/REM)" );
  
  writer.var( "CYCLE_POS_REL" , "Position within NREMC, relative" );
  writer.var( "CYCLE_POS_ABS" , "Position within NREMC, absolute (mins)" );

  writer.var( "FLANKING_SIM" , "Number of similar epochs w.r.t. stage" );
  writer.var( "NEAREST_WAKE" , "Number of epochs until the nearest wake" );

  writer.var( "WASO" , "Epoch is WASO (1=Y)" );
    
  
  // these next four are all reported for the  NREM epoch
  writer.var( "NREM2REM" , "If NREM epoch, number of NREM if next non-NREM is REM" );
  writer.var( "NREM2REM_TOTAL" , "If NREM epoch, total number of contiguous NREM if next non-NREM is REM" );

  writer.var( "NREM2WAKE" , "If NREM epoch, number of NREM if next non-NREM is WAKE" );
  writer.var( "NREM2WAKE_TOTAL" , "If NREM epoch, total number of contiguous NREM if next non-NREM is WAKE" );

  writer.var ("N2_WGT" , "Score for descending/ascending N2 epochs (-1 to +1)" );

  
  
  // output
  for (int e=0;e<ne;e++)
    {
      
      // epoch-level stratification
      writer.epoch( timeline->display_epoch( e ) );
      
      // stage stats
      writer.value( "E_WAKE" , elapsed_wake );
      writer.value( "E_WASO" , elapsed_waso );
      writer.value( "E_SLEEP" , elapsed_sleep );
      writer.value( "E_N1" , elapsed_n1 );
      writer.value( "E_N2" , elapsed_n2 );
      writer.value( "E_N3" , elapsed_n34 );
      writer.value( "E_REM" , elapsed_rem );		  

      // and as percentages
      writer.value( "PCT_E_SLEEP" , TST>0 ? elapsed_sleep / TST : 0 );

      writer.value( "PCT_E_N1" , mins_n1 > 0 ? elapsed_n1 / mins_n1 : 0 );
      writer.value( "PCT_E_N2" , mins_n2 > 0 ? elapsed_n2 / mins_n2 : 0 );
      writer.value( "PCT_E_N3" , (mins_n3+mins_n4) > 0 ? elapsed_n34 / (mins_n3+mins_n4) : 0 );
      writer.value( "PCT_E_REM" , mins_rem > 0 ? elapsed_rem / mins_rem : 0 );

      // track elapsed time
      if ( stages[e] == WAKE ) 
	{
	  elapsed_wake += epoch_mins;
	  if ( e > first_sleep_epoch && e < final_wake_epoch ) elapsed_waso += epoch_mins;
	}
      else if ( stages[e] == NREM1 ) { elapsed_sleep += epoch_mins; elapsed_n1  += epoch_mins; }
      else if ( stages[e] == NREM2 ) { elapsed_sleep += epoch_mins; elapsed_n2  += epoch_mins; }
      else if ( stages[e] == NREM3 || stages[e] == NREM4 ) { elapsed_sleep += epoch_mins; elapsed_n34 += epoch_mins; }
      else if ( stages[e] == REM ) { elapsed_sleep += epoch_mins; elapsed_rem += epoch_mins; }
      
      // persistent sleep
      
      writer.value( "PERSISTENT_SLEEP" , in_persistent_sleep[e] );
      
      // cycles
      
      if ( sleep_cycle_number[e] )
	{
	  writer.value( "CYCLE" , sleep_cycle_number[e] );
	  writer.value( "PERIOD" , sleep_code[e] == 5 ? "REMP" : sleep_code[e]==1 ? "NREMP" : "." ) ;
	  writer.value( "CYCLE_POS_REL" , cycle_pos_relative[e] );
	  writer.value( "CYCLE_POS_ABS" , cycle_pos_absolute[e] );
	}
      
      // flanking epochs

      writer.value( "FLANKING_SIM" , flanking[e] );
      writer.value( "NEAREST_WAKE" , nearest_wake[e] );
      writer.value( "WASO" , is_waso[e] );

      writer.value( "NREM2REM" , nrem2rem[e] );
      writer.value( "NREM2REM_TOTAL" , nrem2rem_total[e] );

      writer.value( "NREM2WAKE" , nrem2wake[e] );
      writer.value( "NREM2WAKE_TOTAL" , nrem2wake_total[e] );

      // N2 ascending/descending status

      if ( stages[e] == NREM2 ) 
	writer.value( "N2_WGT" , n2_ascdesc[e] );
           
    } // next epoch

  writer.unepoch();


}

void dummy_hypno()
{
  edf_t edf;
  
  // dummy values
  
  
  hypnogram_t h;
  h.timeline = &edf.timeline;
  
  while( ! std::cin.eof() )
    {
      std::string s;
      std::cin >> s;
      if ( std::cin.eof() ) break;
      if ( s == "W" ) h.stages.push_back( WAKE );
      else if ( s == "N1" ) h.stages.push_back( NREM1 );
      else if ( s == "N2" ) h.stages.push_back( NREM2 );
      else if ( s == "N3" ) h.stages.push_back( NREM3 );
      else if ( s == "N4" ) h.stages.push_back( NREM4 );
      else if ( s == "R"  ) h.stages.push_back( REM );
      else if ( s == "L"  ) h.stages.push_back( LIGHTS_ON );
      else if ( s == "?"  ) h.stages.push_back( UNSCORED );
      else logger << "did not recognize " << s << "\n";
    }

  logger << "read " << h.stages.size() << "\n";

  edf.header.starttime = "10:00:00";
  
  // fudge so it works in this non-standard case...
  edf.id = "_DUMMY_";
  h.fudge( 30 , h.stages.size() );

  h.calc_stats( true );
  h.output( true ); // verbose mode == T 

}

void hypnogram_t::fudge( double es, int ne )
{
  timeline->epoch_length_tp = es * globals::tp_1sec;
  timeline->epochs.resize( ne );
}




void timeline_t::list_all_annotations( const param_t & param )
{

  //
  // Options
  //

  // count annotations per epoch
  bool per_epoch = param.has( "epoch" );
  if ( per_epoch && ! epoched() ) 
    {
      int ne = set_epoch( globals::default_epoch_len , globals::default_epoch_len );
      logger << " set epochs to default " << globals::default_epoch_len << " seconds, " << ne << " epochs\n";
    }

  // how to decide whether an interval overlaps a mask or not?
  //  start  -- keep annotations that start in an unmasked region
  //  any    -- keep annotations that have any overlap in an unmasked region
  //  all    -- only keep annotations that are completely within unmasked regions
  
  int keep_mode = 0; 
  if ( param.has( "any" ) ) keep_mode = 0;
  if ( param.has( "all" ) ) keep_mode = 1;
  if ( param.has( "start" ) ) keep_mode = 2;  
  
  logger << " keeping annotations based on ";
  if ( keep_mode == 0 )      logger << "any overlap with";
  else if ( keep_mode == 1 ) logger << "complete (all) overlap with";
  else if ( keep_mode == 2 ) logger << "starting in";
  logger << " an unmasked region\n";
  
  bool show_masked = param.has("show-masked");

  // annotation names

  std::vector<std::string> names = annotations.names();

  //
  // Per epoch summary of all annotations
  //

  if ( per_epoch ) 
    {
      
      first_epoch();
      
      while ( 1 ) 
	{
	  
	  int e = show_masked ? next_epoch_ignoring_mask() : next_epoch();
	  
	  if ( e == -1 ) break;
	  
	  writer.epoch( display_epoch( e ) );

	  interval_t interval = epoch( e );
	  
	  // get each annotations
	  for (int a=0;a<names.size();a++)
	    {
	      
	      annot_t * annot = annotations.find( names[a] );
	      
	      // get overlapping annotations for this epoch
	      annot_map_t events = annot->extract( interval );

	      // list
	      annot_map_t::const_iterator ii = events.begin();
	      while ( ii != events.end() )
		{	  
		  
		  const instance_idx_t & instance_idx = ii->first;

		  const instance_t * instance = ii->second;

		  const interval_t & interval = instance_idx.interval;
		  
		  bool is_masked = false;
		  
		  // keep if any part of A overlaps any unmasked region
		  if      ( keep_mode == 0 ) is_masked = ! interval_overlaps_unmasked_region( interval );
		  
		  // ...or, only if entire A is in unmasked region
		  else if ( keep_mode == 1 ) is_masked = ! interval_is_completely_unmasked( interval );
		  
		  // ...or, if start of A is in an unmasked region
		  else if ( keep_mode == 2 ) is_masked = interval_start_is_masked( interval ) ;
		  
		  // skip?
		  if ( is_masked && ! show_masked ) {  continue; } 
		  
		  // else display
		  writer.level( instance_idx.id , "INST" );
		  writer.level( interval.as_string() , "INTERVAL" );

		  writer.value( "EPOCH_MASK" , masked( e ) );
		  writer.value( "ANNOT_MASK" , is_masked );
		  
		  
		  ++ii;
		}      

	      writer.unlevel( "INTERVAL" );
	      writer.unlevel( "INST" );

	    }
	}
      
      writer.unepoch();

      // all done now for epoch-stratified listing
      return;
    }
  


  //
  // Get all annotations (i.e. not stratified by epoch)
  //

  
  // sort by time, collapse across events
  std::map<instance_idx_t,const instance_t*> events;
  
  // class
  std::map<std::string,int> counts;
  std::map<std::string,double> dur;
  
  // class x inst
  std::map<std::string,std::map<std::string,int> > counts2;
  std::map<std::string,std::map<std::string,double> > dur2;
  
  
  // iterate over each annotation
  for (int a = 0 ; a < names.size() ; a++ ) 
    {
      
      annot_t * annot = annotations.find( names[a] );
      
      if ( annot == NULL ) Helper::halt( "internal problem in list_all_annotations()" );

      const int num_events = annot->num_interval_events();
      
      if ( 0 ) 
	{
	  std::cout << names[a] << "\n";
	  std::cout << " ne = " << num_events << "\n";
	  std::cout << " file = " << annot->file << "\n";
	  
	  const int nf = annot->types.size();
	  std::cout << " fields = " << nf << "\n";
// 	  for (int f = 0 ; f < nf ; f++ ) 
// 	    std::cout << " " << annot->type_string(f) ;
	  std::cout << "\n";	  

	}

      
      //
      // iterator over interval/event map
      //

      annot_map_t::const_iterator ii = annot->interval_events.begin();
      while ( ii != annot->interval_events.end() )
	{	  

	  const instance_idx_t & instance_idx = ii->first;
	  const instance_t * instance = ii->second;
	  
	  bool keep_this = false;
	  
	  if      ( keep_mode == 0 ) keep_this = interval_overlaps_unmasked_region( instance_idx.interval );
	  else if ( keep_mode == 1 ) keep_this = interval_is_completely_unmasked( instance_idx.interval );
	  else if ( keep_mode == 2 ) keep_this = ! interval_start_is_masked( instance_idx.interval ) ;
	  
	  if ( keep_this )
	    {      
	      events[ instance_idx ] = instance ; 
	      
	      counts[ annot->name ]++;
	      counts2[ annot->name ][ instance_idx.id ] ++;
	      
	      dur[ annot->name ] += instance_idx.interval.duration_sec(); 	      
	      dur2[ annot->name ][ instance_idx.id ] += instance_idx.interval.duration_sec(); 	      
	    }
	      
	  
	  ++ii;
	}
    }


  // now print all by time point
  std::map<instance_idx_t,const instance_t*>::const_iterator aa = events.begin();
  while ( aa != events.end() )
    {
      
      const instance_idx_t & instance_idx = aa->first;
      
      const interval_t & interval = instance_idx.interval;
      
      const instance_t * instance = aa->second;
      
      // stratify output by interval
    
      writer.interval( interval );
      
      writer.level( instance_idx.parent->name , globals::annot_strat );
      
      writer.level( instance_idx.id , globals::annot_instance_strat );	  

      writer.value( "START" , interval.start_sec() );
      
      writer.value( "STOP" , interval.stop_sec() );
      
      if ( ! instance->empty() ) 
	{
	  writer.value(  "VAL" , instance->print() );
	}

      if ( show_masked ) 
	{
	  
	  bool start_masked = interval_start_is_masked( interval ) ; 
	  bool some_masked = interval_overlaps_masked_region( interval );
	  bool all_masked = interval_is_completely_masked( interval );
	  bool some_unmasked = interval_overlaps_unmasked_region( interval );
	  bool all_unmasked = interval_is_completely_unmasked( interval );
	  
	  writer.value( "START_MASKED"  , start_masked );
	  writer.value( "SOME_MASKED"   , some_masked );
	  writer.value( "ALL_MASKED"    , all_masked );
	  writer.value( "SOME_UNMASKED" , some_unmasked );
	  writer.value( "ALL_UNMASKED"  , all_unmasked );
	}
      
      writer.unlevel( globals::annot_instance_strat );

      writer.unlevel( globals::annot_strat );
      
      ++aa;
      
    }
  writer.uninterval();


  //
  // final counts, durations by class
  //
  
  std::map<std::string,int>::const_iterator cc = counts.begin();
  while ( cc != counts.end() ) 
    {

      writer.level( cc->first , globals::annot_strat );
      writer.value( "COUNT" , cc->second );      
      writer.value( "DUR" , dur[ cc->first ] );

      if ( counts2[ cc->first ].size() > 0 )
	{
	  std::map<std::string,int>::const_iterator dd = counts2[ cc->first ].begin();
	  while ( dd != counts2[ cc->first ].end() )
	    {	      
	      writer.level( dd->first , globals::annot_instance_strat );	  
	      writer.value( "COUNT" , dd->second );      
	      writer.value( "DUR" , dur2[ cc->first ][ dd->first ] );
	      ++dd;
	    }
	  writer.unlevel( globals::annot_instance_strat );
	}

      ++cc;
    }
  writer.unlevel( globals::annot_strat );
}


// eval-based mask
void timeline_t::apply_eval_mask( const std::string & str , int mask_mode )
{

  // mask_mode   0   mask
  //             1   unmask
  //             2   force  (mask & unmask) 
  

  //
  // allow both " and # quoting of EVAL expressions
  //

  std::string expression = Helper::trim( Helper::unquote( str , '#' ) );
  

  //
  // Get all existing annotations (overkill...)
  //
  
  std::vector<std::string> names = annotations.names();

  //
  // Keep track of changes
  //
  
  mask_set = true;

  const int ne = epochs.size();
  
  //
  // We do not clear the mask here, as we want to allow multiple
  // filters to be added on top of oneanther
  //

  int cnt_mask_set = 0;
  int cnt_mask_unset = 0;
  int cnt_unchanged = 0;
  int cnt_now_unmasked = 0;
  int cnt_basic_match = 0;  // basic count of matches, whether changes mask or not

  
  //
  // Iterate over epochs
  //

  first_epoch();
  
  int acc_total = 0 , acc_retval = 0 , acc_valid = 0; 
  
  while ( 1 ) 
    {
      
      int e = next_epoch_ignoring_mask() ;

      if ( e == -1 ) break;
      
      interval_t interval = epoch( e );
	  
      std::map<std::string,annot_map_t> inputs;
      
      // get each annotations
      for (int a=0;a<names.size();a++)
	{
	  
	  annot_t * annot = annotations.find( names[a] );
	  
	  // get overlapping annotations for this epoch
	  annot_map_t events = annot->extract( interval );
	  
	  // store
	  inputs[ names[a] ] = events;
	}

      //
      // create a dummy new instance for the output variables (not saved)
      //
      
      instance_t dummy;
      
      //
      // evaluate the expression, but note, this is set to not 
      // allow any assignments.... this makes it cleaner and easier 
      // to spot bad//undefined variables as errors.
      //

      const bool no_assignments = true;
      
      Eval tok( expression , no_assignments );

      tok.bind( inputs , &dummy );

      bool is_valid = tok.evaluate();
      
      bool matches;
       
      if ( ! tok.value( matches ) ) is_valid = false;
      
      //
      // A match must be a valid value
      //
      
      if ( ! is_valid ) matches = false;

      
      //
      // apply mask (or not)
      //
      
      acc_total++;

      acc_valid += is_valid;
      
      if ( acc_valid ) acc_retval += matches;

      
      // count basic matches

      if ( matches ) ++cnt_basic_match;
      
      // set new potential mask, depending on match_mode
      
      bool new_mask = mask[e];

      if      ( mask_mode == 0 ) new_mask = matches;   // mask
      else if ( mask_mode == 1 ) new_mask = !matches;  // unmask
      else if ( mask_mode == 2 ) new_mask = matches ;  // mask/unmask
      
      int mc = set_epoch_mask( e , new_mask );
      
      if      ( mc == +1 ) ++cnt_mask_set;
      else if ( mc == -1 ) ++cnt_mask_unset;
      else                 ++cnt_unchanged;
      
      if ( !mask[e] ) ++cnt_now_unmasked;
      
      
      // next epoch
    } 
  
  
  logger << " based on eval expression [" << expression << "]\n"
	 << "  " << acc_retval << "  true, " << acc_valid - acc_retval << " false and " 
	 << acc_total - acc_valid << " invalid return values\n"
	 << "  " << cnt_basic_match << " epochs match; " 
	 << cnt_mask_set << " newly masked, "
	 << cnt_mask_unset << " unmasked, "
	 << cnt_unchanged << " unchanged\n";
  logger << "  total of " << cnt_now_unmasked << " of " << epochs.size() << " retained\n";

  
  // mask, # epochs masked, # epochs unmasked, # unchanged, # total masked , # total epochs
  
  writer.level( expression , "EPOCH_MASK" );
  
  writer.var( "N_MATCHES"    , "Number of matching epochs" );
  writer.var( "N_MASK_SET"   , "Number of epochs newly masked" ); 
  writer.var( "N_MASK_UNSET" , "Number of epochs newly unmasked" );
  writer.var( "N_UNCHANGED"  , "Number of epochs unchanged by this mask" );
  writer.var( "N_RETAINED"   , "Number of epochs retained for analysis" );
  writer.var( "N_TOTAL"      , "Total number of epochs" );

  writer.value( "N_MATCHES"    , cnt_basic_match  );
  writer.value( "N_MASK_SET"   , cnt_mask_set     );
  writer.value( "N_MASK_UNSET" , cnt_mask_unset   );
  writer.value( "N_UNCHANGED"  , cnt_unchanged    );
  writer.value( "N_RETAINED"   , cnt_now_unmasked );
  writer.value( "N_TOTAL"      , epochs.size()    );

  writer.unlevel( "EPOCH_MASK" );


  // all done 
  return;

}

