
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
#include <cstddef>

extern writer_t writer;

extern logger_t logger;

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

  // clear any cache
  cache.clear();
  
  logger << "  retaining " << num_epochs() << " epochs\n";
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

  if ( interval.stop < interval.start ) 
    Helper::halt( "badly defined interval requested, with stop before start" );

  // empty record set?
  if ( interval.stop == interval.start )
    {
      *start_rec = 0;
      *start_smp = 0;	  
      *stop_rec = 0;
      *stop_smp = 0;
      return false;
    }

  //    std::cout << "i2r: interval = " << interval << "\n";
  
  //
  // Note: here we want to find records/samples that are inclusive w.r.t. the interval
  // so change coding of stop being 1 unit past the end below
  //

  if ( interval.stop == 0 ) 
    Helper::halt( "internal error in timeline()" );

  // subtract 1 off this timepoint

  uint64_t stop_tp = interval.stop - 1LLU;

  // allow 0-duration annotations, but note that these often
  // need fixes if end if 1 TP before the start (as can happen if the
  // annotation start/stop fall between sample-points); this is handled at
  // the end
  
  if ( interval.start > stop_tp ) return false;

  //
  // For a continuous timeline, given time-points can
  // straightforwardly calculate record/sample
  //

  if ( edf->header.continuous )
    {
      //      std::cout << "EDF-C\n";

      // get initial records/samples, nb. use ceil() to get nearest sample *after* start of interval

      uint64_t start_record = interval.start / edf->header.record_duration_tp;
      uint64_t start_offset = interval.start % edf->header.record_duration_tp;
      uint64_t start_sample = 
	ceil( ( start_offset / (double)edf->header.record_duration_tp ) * n_samples_per_record ) ;
      
      //      std::cout << "othr = " << edf->header.record_duration_tp << " " << n_samples_per_record << "\n";
      //      std::cout << "start = " << start_record << " " << start_offset << " " << start_sample << "\n";
      
      // get final records/samples, nb. use floor() to get the nearest sample *prior* to end
      
      uint64_t stop_record = stop_tp / edf->header.record_duration_tp;
      uint64_t stop_offset = stop_tp % edf->header.record_duration_tp;
      uint64_t stop_sample = 
	floor( ( stop_offset / (double)edf->header.record_duration_tp ) * n_samples_per_record  ) ;
      
      // pass back to calling function
      
      *start_rec = (int)start_record;
      *start_smp = (int)start_sample;

      *stop_rec = (int)stop_record;
      *stop_smp = (int)stop_sample;

    }
  else
    {
      
      // std::cout << "in EDF+D\n";
      
      //
      // For a discontinuous EDF+ we need to search 
      // explicitly across record timepoints
      //
      
      //
      // Get first record that is not less than start search point (i.e. equal to or greater than)
      //
      
      std::map<uint64_t,int>::const_iterator lwr = tp2rec.lower_bound( interval.start ); 
           
      //
      // This will find the first record AFTER the start; thus, if the
      // interval is aligned to the sample point, we'll need to skip
      // one back; this should never be the first record, but check in
      // case...
      //
      
      // track whether the start search point falls between sample-points
      
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
	  // If the search point occurs before /all/ records, need to
	  // indicate that we are in a gap also	  
	  
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
	    ceil( ( start_offset / (double)edf->header.record_duration_tp ) * n_samples_per_record ) ;
	  //if ( start_sample >= n_samples_per_record ) start_sample = n_samples_per_record - 1LLU; 
	  *start_smp = (int)start_sample;
	}
      
      
      //
      // for upper bound, find the record whose end is equal/greater *greater* 
      // 
      
      std::map<uint64_t,int>::const_iterator upr = tp2rec.upper_bound( stop_tp ); 
      
      //
      // this should have returned one past the one we are looking for 
      // i.e. that starts *after* the search point
      //
      
      bool ends_before = upr == tp2rec.begin() ;
      
      if ( ! ends_before ) 
	{
	  --upr;  
	  *stop_rec  = upr->second;
	}
      else
	{
	  *stop_rec  = -1; // i.e. flag as bad, so ensure that stop is before the start (which will also be rec 0)
	}

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
	    floor( ( stop_offset / (double)edf->header.record_duration_tp ) * n_samples_per_record ) ;
	  //if ( stop_sample >= n_samples_per_record ) stop_sample = n_samples_per_record - 1LLU;
	  *stop_smp = (int)stop_sample;
	}
      
    }
  
  //
  // If the interval is in a gap, we will not get any records here, and 
  // stop < start;  so check for this and flag if so
  //
  // however, we may have an annotaiton of very small duration (i.e. at the extreme a 0-second annotation
  // which falls between two sample points, and stop_rec may be 1 less than the start rec;  allow this scenario
  // by adding the 1LLU to stop_smp when start_rec and stop_rec are the same
  //

    
    // allow for off-by-one (i.e. if 0-duration point falls between time-points)

  if ( *start_rec == *stop_rec && *start_smp == *stop_smp + 1LLU )
    {
        *stop_smp = *start_smp;
    }
  
  else if ( *start_rec == *stop_rec + 1LLU
   	    && *stop_smp == 0 // first point of new record
   	    && *start_smp == n_samples_per_record - 1 ) // last point of previous record
     {
       *stop_rec = *start_rec;
       *stop_smp = *start_smp;
     }
  
  // otherwise, check for a gap/problem
  else if ( *start_rec > *stop_rec || ( *start_rec == *stop_rec && *start_smp > *stop_smp ) )
    {
      *start_rec = *start_smp = *stop_rec = *stop_smp = 0;
      return false;
    }
  
  //
  // Otherwise, we're all good
  //
  
  // std::cout << "recs = " << *start_rec << " " << *stop_rec << "\n";
  // std::cout << "smps = " << *start_smp << " " << *stop_smp << "\n";
    
  return true;
  
}



int timeline_t::calc_epochs()
{

  // Calculate EPOCH timings (in the original EDF timescale) for both
  // continuous and discontinuous EDFs

  // Also populate rec2epoch and epoch2rec mappings
  

  // **Epoch size has to be a multiple of EDF record size**

  //
  // currently, epoch length and duration must be multiples of the EDF record size
  // we can likely free this constraint easily, as we allow a one-to-many mapping of 
  // epochs to records and records to epochs, i.e. one record could be in > 1 epoch
  //


  //
  // In all cases, epoch length must be >= record length
  //
  
   if ( epoch_length_tp < edf->header.record_duration_tp )
     Helper::halt( "epoch duration must be greater or equal to EDF record size\n         which is currently "
   		  + Helper::dbl2str( edf->header.record_duration_tp * globals::tp_duration ) + " second(s); "
   		  + "see RECORD-SIZE command to change this" );

   if ( epoch_length_tp < epoch_inc_tp ) Helper::halt( "epoch increment cannot be larger than epoch duration" );

  
  epochs.clear();
  
  mask.clear();
  
  rec2epoch.clear();
  
  epoch2rec.clear();


  //
  // Easy case for continuous EDFs
  //

  if ( edf->header.continuous )
    {

      // set first timepoint
      uint64_t s = 0;
      
      while ( 1 ) 
	{
	  
	  // get end of interval: for this purpose (of finding records)
	  // we set last point of 
	  uint64_t end = s + epoch_length_tp - 1LLU;
	  
	  // done? [ skip any final epochs that do not fit into the frame ] 	  
	  if ( end >= total_duration_tp ) break;
	  
	  // otherwise, add to list, but with end as +1 past end
	  interval_t interval( s , end + 1LLU );
	  epochs.push_back( interval );
	  
	  // find matching records (within interval)
	  interval_t search_interval( s , end );
	  
	  // records in this epoch
	  int start_record = search_interval.start / edf->header.record_duration_tp;  
	  int stop_record = search_interval.stop / edf->header.record_duration_tp;
	  int e = epochs.size()-1;
	  for (int r=start_record; r<=stop_record; r++) 
	    {
	      epoch2rec[ e ].insert( r );
	      rec2epoch[ r ].insert( e );
	    }
	  
	  // shift to next interval
	  s += epoch_inc_tp;
	}
    }
  else
    {

      //
      // Epochs for the discontinuous case:
      //
      
      //      std::cout << "HERE\t" << epoch_length_tp << "\t" << epoch_inc_tp << "\n";

      int r = first_record();

      if ( r == -1 ) return 0;
      
      // epochs have to be continuous in clocktime
      uint64_t estart = rec2tp[r];

      // for purpose of searching, skip last point
      // i.e. normally intervals are defined as END if 1 past the last point
      uint64_t estop  = estart + epoch_length_tp - 1LLU;

      // also, track so we can figure out which record the next epoch starts in 
      // for non-overlapping epochs, this will equal estop, otherwise it will come earlier
      
      uint64_t erestart = estart + epoch_inc_tp;
      
      int restart_rec = -1 ;
      
      // for epoch2rec, rec2epoch mapping
      int e = 0;
      
      // for putative epochs, track the epoch-record mapping
      std::map<int,std::set<int> > putative_e2r;
      std::map<int,std::set<int> > putative_r2e;

      // move forward one record at a time
      // to find the record where 
      
      while ( 1 ) 
	{
	  
	  //
	  // Start and end of this current record	  
	  //

	  uint64_t rec_start = rec2tp[r];
	  uint64_t rec_end   = rec2tp_end[r];
	  
// 	  std::cout << "dets " << rec_start << " " << rec_end << "\t"
// 		    << estart << " " << erestart << " " << estop << "\n";

	  //
	  // Will the next epoch potentially come from this record?
	  //
	  
	  if ( erestart >= rec_start && erestart <= rec_end )
	    {
	      // track this is the restarting record
	      restart_rec = r;	      
	    }
	  
	  //
	  // Will epoch end in this record? 
	  //
	  
	  if ( estop <= rec_end )
	    {
	      
	      // add to list of epochs (note: adding +1 when saving interval)
	      interval_t saved_interval( estart , estop + 1LLU );
	      epochs.push_back( saved_interval );
	      
	      // std::cout << "E " << epochs.size() << " adding " << estart << " -- " << estop + 1LLU << "\n";
	      
	      // add this last record to the epoch in any case
	      putative_r2e[ r ].insert( e );
	      putative_e2r[ e ].insert( r );

	      // add record mappings
	      std::map<int,std::set<int> >::const_iterator ii = putative_r2e.begin();
	      while ( ii != putative_r2e.end() )
		{
		  std::set<int>::const_iterator jj = ii->second.begin();
		  while ( jj != ii->second.end() )
		    {
		      rec2epoch[ ii->first ].insert( *jj );
		      ++jj;
		    }
		  ++ii;
		}
	      
	      ii = putative_e2r.begin();
	      while ( ii != putative_e2r.end() )
		{
		  std::set<int>::const_iterator jj = ii->second.begin();
		  while ( jj != ii->second.end() )
		    {
		      epoch2rec[ ii->first ].insert( *jj );
		      ++jj;
		    }
		  ++ii;
		}
	      
	      // clear temporary markers
	      putative_e2r.clear();
	      putative_r2e.clear();
    
	      // move on to the next epoch
	      ++e;
	    
	      
	      
	      //
	      // Now, we're moving to the next epoch and have to figure out
	      // where to restart
	      //
	      
	      //
	      // For non-overlapping, the next epoch may likely start
	      // in the /next/ record (i.e. if this epochs ends right at
	      // the end of the current record).  In this case, we
	      // will not have set restart_rec yet, so it will have its
	      // default value of -1
	      // 

	      if ( restart_rec == -1 ) 
		{		  
		  // we need to jump to the next record
		  // and no need to check for discontinuity
		  // as we are starting a new epoch in any case
		  r = next_record(r);		  		  
		  
		  // all done...
		  if ( r == -1 ) break;		  
		  
		  // set start point here, as this record may skip ahead of
		  // assumed eretsart

		  erestart = rec2tp[r];

		}
	      else
		{
		  // otherwise, we have already encountered the starting record		  
		  // and we have the 
		  r = restart_rec;	      	  		  
		}
	      
	      //
	      // Reset the epoch goals
	      //

	      restart_rec = -1;
	      estart = erestart;
	      estop = estart + epoch_length_tp - 1LLU;
	      erestart = estart + epoch_inc_tp;
	      
	      //
	      // Track this epoch/record pairing 
	      //
	      
	      putative_r2e[ r ].insert( e );
	      putative_e2r[ e ].insert( r );
	  
	      //
	      // all done for this new epoch; do not advance the record counter, 
	      // just let the loop go forwards
	      //

	    }
	  
	  
	  //
	  // if current epoch will *not* end in this epoch
	  // then just add this record to the mapping
	  // and get the next record
	  //
		  
	  else 
	    {

	      // should this not be putative_ cases?

	      // original
// 	      rec2epoch[r].insert(e);
// 	      epoch2rec[e].insert(r);
	      
	      // revised
	      putative_r2e[r].insert(e);
	      putative_e2r[e].insert(r);

	      r = next_record(r);
	      
	      // if out of epochs, then just quit (i.e. we won't be adding this
	      // putative epoch in any case)
	      
	      if ( r == -1 ) break;
	      
	      // check that we did not skip anything
	      // these two values should be EQUAL is
	      // they are contiguous 
	      
	      uint64_t rec2_start = rec2tp[r];
	      
	      //std::cout << "recs " << rec2_start << "\t" << rec_end << "\n";

	      // 
	      // If we've come to a break, we need to give up on the previous
	      // epoch, and start adding one now for this next record
	      //
	      
	      if ( rec2_start - rec_end != 1 ) 
		{		  
		  estart = rec2_start;
                  estop  = estart + epoch_length_tp - 1LLU;
		  erestart = estart + epoch_inc_tp;
		  putative_e2r.clear();
		  putative_r2e.clear();		  
		}
	      
	    }
	  
	  // next record 

	}
    }


  if ( 0 ) 
    {
      std::cout << "REC 2 E\n";
      
      std::map<int,std::set<int> >::const_iterator ii = rec2epoch.begin();
      while ( ii != rec2epoch.end() )
	{
	  std::cout << "r" << ii->first << "->";
	  std::set<int>::const_iterator jj = ii->second.begin();
	  while ( jj != ii->second.end() ) 
	    {
	      std::cout << " " << *jj;
	      ++jj;
	    }
	  std::cout << "\n";
	  ++ii;
	}
      std::cout << "\n";

      std::cout << "E 2 REC\n";
      
      ii = epoch2rec.begin();
      while ( ii != epoch2rec.end() )
	{
	  std::cout << "e" << ii->first << "->";
	  std::set<int>::const_iterator jj = ii->second.begin();
	  while ( jj != ii->second.end() ) 
	    {
	      std::cout << " " << *jj;
	      ++jj;
	    }
	  std::cout << "\n";
	  ++ii;
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
  
  if ( e < 0 || e >= mask.size() ) Helper::halt( "internal error setting mask" );
  
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
  
  writer.level( label  , "EMASK" );

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
  writer.value( "N_TOTAL"      , (int)epochs.size()    );

  writer.unlevel( "EMASK" );

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
  
  writer.level( a->name , "EMASK" );

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
  writer.value( "N_TOTAL"      , (int)epochs.size()    );

  writer.unlevel( "EMASK" );
}


signal_list_t timeline_t::masked_channels_sl( const int e0 , const signal_list_t & signals ) const
{
  const bool silent_mode = true;
  int e = display_epoch( e0 );
  //  std::cout << "e in TLM = " << e0 << " " << e << "\n";
  signal_list_t msigs;
  std::vector<std::string> m = masked_channels( e0 , signals );
  for (int i=0;i<m.size();i++) 
    {
      int chn = edf->header.signal( m[i] , silent_mode );
      if ( chn != -1 ) 
	msigs.add( chn , m[i] );
    }
  return msigs;
}

signal_list_t timeline_t::unmasked_channels_sl( const int e0 , const signal_list_t & signals ) const
{
  const bool silent_mode = true;
  signal_list_t usigs;
  int e = display_epoch( e0 );
  if ( e == -1 ) return usigs;
  std::vector<std::string> u = unmasked_channels( e0 , signals );
  for (int i=0;i<u.size();i++) 
    {
      int chn = edf->header.signal( u[i] , silent_mode );
      if ( chn != -1 ) 
	usigs.add( chn , u[i] );
    }
  return usigs;
}


std::vector<std::string> timeline_t::masked_channels( const int e0 , const signal_list_t & signals ) const
{
  int e = display_epoch( e0 );
  //  std::cerr << " e , e0 = " << e << " " << e0 << "\n";
  std::vector<std::string> m;
  const int ns = signals.size();
  bool any_masked = chep.find( e ) != chep.end() ;
  if ( ! any_masked ) return m; // all good

  const std::set<std::string> & masked = chep.find(e)->second ;
  for (int s=0; s<ns; s++) 
    {
      if ( masked.find( signals.label(s) ) != masked.end() )
	m.push_back( signals.label(s) );
    }
  return m;
}


std::vector<std::string> timeline_t::unmasked_channels( const int e0 , const signal_list_t & signals ) const
{

  int e = display_epoch( e0 );

  std::vector<std::string> u;
  const int ns = signals.size();
  bool any_masked = chep.find( e ) != chep.end() ;
  if ( ! any_masked ) 
    {
      // all good
      for (int s=0; s<ns; s++) u.push_back( signals.label(s) );
      return u; 
    }

  const std::set<std::string> & masked = chep.find(e)->second ;
  for (int s=0; s<ns; s++) 
    {
      if ( masked.find( signals.label(s) ) == masked.end() )
	u.push_back( signals.label(s) );
    }
  return u;
}

void timeline_t::collapse_chep2epoch( signal_list_t signals , const double pct , const int k )
{
  
  // drop any non-data channels (modifies 'signals')
  edf->header.drop_annots_from_signal_list( &signals );

  logger << "  masking epochs";
  if ( k ) logger << " with " << k << " or more masked channels";
  if ( pct < 1 ) logger << ( k ? ", or " : " with >" ) << pct * 100 << "% masked channels: ";

  // if k or more channels are masked --> set epoch mask  [ default 1 ] 
  // if more than pct channels are masked --> set epoch mask [ default 0 ]
  // automatically adjust main 'mask' (using set_mask(), i.e. respecting mask_mode etc)
  
  int epoch = 0; 
  int masked = 0;
  int masks_set = 0;

  std::map<int,std::set<std::string> >::iterator ee = chep.begin();
  while ( ee != chep.end() )
    {
      // **assume** same signals overlap

      int sz = ee->second.size();

      int epoch = ee->first;

      if ( ( k != 0 && sz >= k ) || 
	   ( sz / (double)signals.size() > pct ) ) 
	{
	  // change main epoch mask
	  int epoch0 = display2curr_epoch( epoch );

	  // if this epoch is still present in current file, set mask
	  if ( epoch0 != -1 ) 
	    if ( set_epoch_mask( epoch0 ) ) ++masked;
	  
	  // and also set all CHEP masks (to signals) for this epoch
	  for (int s=0;s<signals.size();s++) ee->second.insert( signals.label(s) );
	  
	}
      
      if ( mask[epoch] ) ++masks_set; 
    
      ++ee;
    }
  
  logger << masked << " epochs\n";
  
}


signal_list_t timeline_t::collapse_chep2ch( signal_list_t signals , 
					    const double pct , const int k  , 
					    bool bad_set_all_bad ,
					    bool good_set_all_good )
{

  // identify which channels (from the set 'signals') have more than 'k' (or more than pct %) of epochs masked
  // return this as a set of 'bad channels'

  // if k or more channels are masked --> set epoch mask  [ default 1 ] 
  // if more than pct channels are masked --> set epoch mask [ default 0 ]
  
  // pct is defined with the denominator as the # of data channels in 'signals'

  // drop any non-data channels (modifies 'signals')
  edf->header.drop_annots_from_signal_list( &signals );

  logger << "  masking channels";
  if ( k ) logger << " with " << k << " or more masked epochs";
  if ( pct < 1 ) logger << (k?", or " : " with > " ) << pct *100 << "% masked epochs:";

  // count of bad epochs per channel
  std::map<std::string,int> c;
  int ns = signals.size();
  int ne = num_epochs();
  for (int i=0; i<ns; i++) c[ signals.label(i) ] = 0;
  
  // get channel slots lookup-table 
  std::map<std::string,int> l2s;
  for (int i=0; i<ns; i++) 
    l2s[ signals.label(i) ] = signals(i);


  std::map<int,std::set<std::string> >::const_iterator ee = chep.begin();
  while ( ee != chep.end() )
    {
      std::set<std::string>::const_iterator ss = ee->second.begin();
      while ( ss != ee->second.end() )
	{
	  if ( c.find( *ss ) != c.end() ) c[ *ss ]++;
	  ++ss;
	}
      ++ee;
    }
  
  signal_list_t good_signals;
  signal_list_t bad_signals;

  const bool silent_mode = true;

  std::map<std::string,int>::const_iterator cc = c.begin(); 
  while ( cc != c.end() ) 
    {
      if ( l2s.find( cc->first ) != l2s.end() ) 
	{
	  if ( ! ( ( k != 0 && cc->second >= k ) || 
		   ( cc->second/(double)ne > pct ) ) )
	    good_signals.add( l2s[ cc->first ] , cc->first );
	  else
	    bad_signals.add( l2s[ cc->first ] , cc->first );
	}
      ++cc;
    }
  
 
  std::set<std::string> good_sigs;
  for (int i=0;i<good_signals.size();i++) 
    good_sigs.insert( good_signals.label(i) );
 
  // set all epochs as masked for a 'bad channel'?
  if ( bad_set_all_bad ) 
    {
      for (int i=0; i<ns; i++) 
	{
	  const std::string label = signals.label(i);
	  if ( good_sigs.find( label ) == good_sigs.end() ) 
	    {
	      logger << " " << label;
	      for (int e=0;e<ne;e++) chep[ display_epoch( e ) ].insert( label );
	    }
	}
    }
      
  // set all eoochs as unmasked for a 'good channel'?
  if ( good_set_all_good )
    {
      for (int i=0; i<ns; i++) 
	{
	  const std::string label = signals.label(i);
	  if ( good_sigs.find( label ) != good_sigs.end() ) 
	    for (int e=0;e<ne;e++) 
	      {
		int e1 = display_epoch( e );
		std::set<std::string>::iterator ii = chep[ e1 ].find( label );
		if ( ii != chep[e].end() ) chep[ e1 ].erase( ii );
	      }
	}
    }

  logger << "\n";
  
  return bad_signals;
}


void timeline_t::dump_chep_mask( signal_list_t signals , bool write_out )
{

  const int ne = first_epoch();

  // for report to console
  int total_masked = 0;
  int total_total = 0;
  std::map<int,int> track_epochs;
  std::map<std::string,int> track_channels;
  
  // use signals list to restrict summary to channels of interest only
  edf->header.drop_annots_from_signal_list( &signals );
  const int ns = signals.size();
  
  std::map<std::string,int> chtots;
  
  while ( 1 ) 
    {
      
      int e = next_epoch_ignoring_mask();      
      
      if ( e == -1 ) break;
      
      int eptot = 0;

      interval_t interval = epoch( e );
      
      int depoch = display_epoch( e );
      
      if ( write_out )
	writer.epoch( depoch );

      if ( chep.find( depoch ) == chep.end() )
	{
	  for (int s=0;s<ns;s++)     
	    {

	      ++total_total;

	      if ( write_out )
		{
		  writer.level( signals.label(s) , globals::signal_strat );
		  writer.value( "CHEP" , false );
		}
	    }
	  
	  if ( write_out )
	    writer.unlevel( globals::signal_strat );
	}
      else
	{

	  // this epoch has 1+ channel masked 
	  
	  track_epochs[ depoch ]++;

	  const std::set<std::string> & ss = chep.find( depoch )->second;

	  for (int s=0;s<ns;s++)     
	    {
	      
	      const std::string label = signals.label(s);

	      // track total
	      ++total_total;
	      
	      bool masked = ss.find( label ) != ss.end() ;
		  
	      if ( write_out )
		{
		  writer.level( label , globals::signal_strat );
		  writer.value( "CHEP" , masked );
		}
	      
	      if ( masked ) 
		{
		  track_channels[ label ]++;
		  ++total_masked;
		  chtots[ label ]++;
		  eptot++;
		}
	      
	    }
	  
	  if ( write_out )
	    writer.unlevel( globals::signal_strat );	  
	}
      
      if ( write_out )
	writer.value( "CHEP" , eptot );
      
    } // next epoch
  
  if ( write_out )
    writer.unepoch();
  
  // ch totals
  if ( write_out )
    {
      for (int s=0;s<ns;s++)     
	{
	  writer.level( signals.label(s) , globals::signal_strat );
	  writer.value( "CHEP" , chtots[ signals.label(s) ] );
	}    
      writer.unlevel( globals::signal_strat );
    }

  //
  // report to console
  //

  int partial_masks = 0; // ignores if epoch or channels complelety masked
  int partial_masks2 = 0; // should match the above..

  int epochs_totally_masked = 0 , channels_totally_masked = 0;
  std::map<int,int>::const_iterator ii = track_epochs.begin();
  while ( ii != track_epochs.end() )
    {
      if ( ii->second == ns ) ++epochs_totally_masked;
      else partial_masks += ii->second;
      ++ii;
    }

  std::map<std::string,int>::const_iterator jj = track_channels.begin();
  while ( jj != track_channels.end() )
    {
      if ( jj->second == ne ) ++channels_totally_masked;
      else partial_masks2 += jj->second ;
      ++jj;
    }

  // ne is the current number of unmasked epochs

  // total epoch count
  const int ne_all = num_total_epochs();

  logger << "  CHEP summary:\n"
	 << "   " << total_masked << " of " << total_total << " channel/epoch pairs masked (" 
	 << round( 100 * ( total_masked / double(total_total ) ) ) << "%)\n"
	 
	 << "   " << track_epochs.size() << " of " << ne_all << " epochs with 1+ masked channel, " 
	 << epochs_totally_masked << " with all channels masked\n"

	 << "   " << track_channels.size() << " of " << ns << " channels with 1+ masked epoch, " 
	 << channels_totally_masked << " with all epochs masked\n";

  // Hmmm... need to revisit this, not clear   
  // ignore totally masked channels/epochs.. how much is left?
  // logger << "   " << partial_masks  << " (" << round( 100 * ( partial_masks / double(total_total ) ) )  << "%) channels/epochs partially masked\n";
  // logger << "   " << partial_masks2 << " (" << round( 100 * ( partial_masks2 / double(total_total ) ) )  << "%) channels/epochs partially masked\n";

}


void timeline_t::read_chep_file( const std::string & f , bool reset )
{

  if ( reset ) clear_chep_mask();
  
  if ( ! Helper::fileExists( f ) ) Helper::halt( f + " does not exist" );

  std::ifstream FIN( f.c_str() , std::ios::in );
  
  // **assumes** same epoch size, channel names
  // user's responsibility to keep that as it should be

  bool silent_mode = true;

  while ( 1 ) 
    {
      std::string ch;
      int e;      
      FIN >> e >> ch ; 
      if ( FIN.eof() ) break;
      if ( ch == "" ) break;      
      int chn = edf->header.signal( ch , silent_mode );      
      if ( chn != -1 ) chep[ e ].insert( ch );  // i.e. expecting display epoch encoding (1-based)
    }
  
  FIN.close();
}

void timeline_t::write_chep_file( const std::string & f ) const
{
  std::ofstream FOUT( f.c_str() , std::ios::out );
  if ( FOUT.bad() ) Helper::halt( "could not open " + f );
  std::map<int,std::set<std::string> >::const_iterator ee = chep.begin();
  while ( ee != chep.end() )
    {
      const std::set<std::string> & chs = ee->second;
      std::set<std::string>::const_iterator cc = chs.begin();
      while ( cc != chs.end() )
	{
	  FOUT << ee->first << "\t" 
	       << *cc << "\n";
	  //	       << edf->header.label[ *cc ] << "\n";
	  ++cc;
	}
      ++ee;
    }
  FOUT.close();
}




bool timeline_t::masked_timepoint( uint64_t a ) const
{

  Helper::halt( "masked_timepoint() not implemented" );

  if ( ! edf->header.continuous ) 
    Helper::halt( "masked_timepoint() not implemented for EDF+D yet" );

    if ( ! mask_set ) return false;
  
  int e1 = MiscMath::position2leftepoch( a , epoch_length_tp, epoch_inc_tp , mask.size() );
  int e2 = MiscMath::position2rightepoch( a , epoch_length_tp, epoch_inc_tp , mask.size() );
  
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


uint64_t timeline_t::valid_tps(  const interval_t & interval )
{

  // given an interval, how many of its time-points fall within the EDF?
  
  // nb:  last_time_point_tp is the actual last time point 
  //      interval.stop is one-past-the-end 

  if ( edf->header.continuous )
    {
      // all out? (last time point tp is  dur - 1 , so GT rather than GE ) 
      if ( interval.start > last_time_point_tp ) return 0;

      // all in?  stop is 1 past last interval 
      if ( interval.stop <= last_time_point_tp + 1LLU ) return interval.duration();
      
      // otherwise, we have partial overlap
      return last_time_point_tp - interval.start + 1LLU;
    }
  else  // for the discontinuous case
    {      
      std::set<int> records = records_in_interval( interval );
      std::set<int>::const_iterator rr = records.begin();
      int tpin = 0;//, tpout = 0;
      while ( rr != records.end() )
	{
	  // start/stop for this record	(as above, does not use 1-past-end encoding here)
	  interval_t rec = record2interval( *rr );
	  
	  // all out? (last time point tp is  dur - 1 , so GT rather than GE )
	  //if ( interval.start > rec.stop ) tpout += rec.duration() + 1LLU;
	  
	  // all in?  stop is 1 past last interval
	  if ( interval.stop <= rec.stop + 1LLU ) tpin += rec.duration() + 1LLU;
	  
	  // otherwise, we have partial overlap
	  tpin += rec.stop - interval.start + 1LLU;
	  //tpout += interval.stop - ( rec.stop + 1LLU );

	  ++rr;
	}
      
      return tpin;
      
    }
  
  return 0;

}


bool timeline_t::masked_interval( const interval_t & interval , bool all_masked , bool * start_masked ) const
{
  
  // if all_masked,   returns T if /all/ of interval falls within masked regions
  // if not,          returns T if interval falls in at least one masked region

  if ( start_masked != NULL ) *start_masked = false;  

  
  if ( edf->header.continuous )
    {

      if ( ! mask_set ) 
	{
	  return false;
	}
      
      int eleft = MiscMath::position2leftepoch( interval.start , epoch_length_tp, epoch_inc_tp , mask.size() );

      // end of interval is one past end of region:
      int eright = MiscMath::position2rightepoch( interval.stop-1LLU , epoch_length_tp, epoch_inc_tp , mask.size() );
      
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

      // if ( ! mask_set )
      // 	{
      //     return false;
      // 	}

          
      std::set<int> records = records_in_interval( interval );

      // falls off edge of the map
      if ( records.size() == 0 ) return true;
      
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

  //  std::cerr << "searching " << interval.as_string() << "\n";
  
  bool any = interval2records( interval , srate , &start_rec , &start_smp , &stop_rec , &stop_smp );

  if ( ! any ) return recs;
  
  int r = start_rec;
  while ( r != -1 ) 
    {
      recs.insert(r);
      r = next_record(r);
      if ( r > stop_rec ) break;
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
  std::set<int> e;
  if ( a > b ) 
    {
      int tmp = b;
      b = a;
      a = tmp;
    }

  for (int i=a; i<=b; i++) e.insert( i );

  if ( include )
    logger << " selecting epochs from " << a << " to " << b << "; ";
  else
    logger << " masking epochs from " << a << " to " << b << "; ";

  return select_epoch_range( e , include );

}

void timeline_t::select_epoch_range( const std::set<int> & specified_epochs , bool include )
{

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

      bool inset = specified_epochs.find( epoch ) != specified_epochs.end() ;

      bool match = include ? 
	! inset : inset ;
      
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
    logger << " selecting";
  else
    logger << " masking";

  logger << " from set of " << specified_epochs.size() << " epochs; ";
  
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
  // Take information from the annot_t class, and make a simple
  // per-epoch annotation this can be performed after a restructure,
  // but (total) # of epochs must match the file exactly
  //


  //
  // Point to first epoch, and get the 'total' number of epochs
  // (i.e. both masked and unmasked), as first_epoch() only returns
  // the unmasked counts; nb. first_epoch() actually sets the pointer
  // *before* the first epoch, so fine to use this whether actual
  // first epoch is masked or not
  //
  
  first_epoch();
  
  int ne = num_total_epochs();
  
  //
  // Populate epoch-annotation vectors to the appropriate size
  //
  
  eannots[ label ].clear();
  

  //
  // Get annotations
  //
  
  annot_t * annot = annotations( annot_label );

  // if not found, then all eannots are effectively false
  // (i.e. missing)

  if ( annot == NULL ) return;


  //
  // for each epoch 
  //
  
  while ( 1 ) 
    {
      
      //
      // Get next epoch
      //
      
      int e = next_epoch_ignoring_mask();      

      if ( e == -1 ) break;


      // use 'e' for look-up here; but need to use orginal EDF
      // encoding for eannots, internally i.e. the zero-based version
      
      int e0 = original_epoch( e ) ;

      if ( e0 == -1 ) 
	Helper::halt( "internal error in annotate_epochs()" );

      interval_t interval = epoch( e );
      
      annot_map_t events = annot->extract( interval );
      
      // search for a matching value (at least one)
      
      annot_map_t::const_iterator ii = events.begin();

      while ( ii != events.end() )
	{	
	  
	  const instance_idx_t & instance_idx = ii->first;
	  const instance_t * instance = ii->second;

	  if ( values.find( instance_idx.id ) != values.end() )
	    {	      
	      // nb. store w.r.t. original epoch encoding e0
	      eannots[ label ][ e0 ] = true;
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


void timeline_t::mask2annot( const std::string & path , const std::string & tag , bool with_id ) 
{

  if ( ! mask_set ) return;
  
  std::string path2 = path[ path.size() - 1 ] != globals::folder_delimiter 
    ? path + globals::folder_delimiter 
    : path ; 
  
  std::string filename = with_id ? ( path2 + tag + "-" + edf->id + ".annot" ) : ( path2 + tag + ".annot" ) ;
  
  annot_t * a = annotations.add( tag );
  a->description = "Mask based on " + tag ;
  //a->types[ "M" ] = globals::A_BOOL_T;
  
  const int ne = mask.size();
  
  for (int e=0;e<ne;e++)
    {
      if ( mask[e] )
	{
	  instance_t * instance = a->add( tag , epoch(e) , "." );
	  //instance->set( "M" , true );
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
      writer.var(   "EMASK" ,      "Is masked? (1=Y)" );
      //      writer.value( "INTERVAL" , interval.as_string() );
      writer.value( "EMASK" , mask_set ? mask[e] : false );

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
  
}


void timeline_t::load_mask( const std::string & f , bool exclude )
{
  
  if ( ! epoched() ) 
    {
      int ne = set_epoch( globals::default_epoch_len , globals::default_epoch_len );
      logger << "  set epochs to default " << globals::default_epoch_len << " seconds, " << ne << " epochs\n";
    }
     
  if ( ! Helper::fileExists( f ) ) Helper::halt( "could not find " + f );
  
  logger << "  attaching mask file " << f << "\n";
  
  logger << "  currently, mask mode set to: ";
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

      Helper::safe_getline( FIN , line );
      
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
  
  writer.level( onelabel , "EMASK" );

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
  writer.value( "N_TOTAL"      , (int)epochs.size()    );

  writer.unlevel( "EMASK" );

}


void timeline_t::annot2signal( const param_t & param )
{
  // create a new signal based on one or more annotations
  if ( ! param.has( "annot" ) ) Helper::halt( "no annotations specified: e.g. annot=A1,A2" );
  std::vector<std::string> anames = param.strvector( "annot" );

  // SR of new signals
  int sr = param.requires_int( "sr" );

  // if not otherwise specified, use annot names as new channel labels
  std::vector<std::string> labels = param.has( "label" ) ? param.strvector( "label" ) : anames;

  if ( anames.size() != labels.size() )
    Helper::halt( "label size does not match annot size" );
  
  // get whole signal size for this SR

  const int np = sr * edf->header.record_duration * edf->header.nr;  

  const uint64_t srtp = (1.0/sr) * globals::tp_1sec;

  // create synthetic signal, 0/1 for presence/absence of the annotation
  for (int a=0; a<anames.size(); a++)
    {
      // does annot exist?
      annot_t * annot = edf->timeline.annotations( anames[a] );
      if ( annot == NULL ) continue;

      // get all events
      const annot_map_t & events = annot->interval_events;

      // new channel to be populated/added to EDF
      std::vector<double> adat( np , 0 );

      annot_map_t::const_iterator aa = events.begin();
      while ( aa != events.end() )
	{

	  const interval_t & interval = aa->first.interval;
	  
	  // convert from time-points to (nearest) sample-points
	  // (after removing the N+1 end point in annotations)
	  
	  int start = interval.start / srtp;
	  int stop  = (interval.stop-1LLU) / srtp;
	  
	  if ( start < 0 || stop >= np )
	    Helper::halt( "internal error in timeline_t::annot2signal()" );
	  
	  // populate (up to and including the start/stop, as we removed the final +1 TP above)
	  for (int p=start; p<=stop; p++) adat[p] = true;

	  // next annotation
	  ++aa;
	}

      //
      // track total time implicated
      //
      
      int points = 0;
      for (int i=0;i<adat.size();i++) if ( adat[i] > 0 ) ++points;
      double seconds = points / sr;
      int minutes = seconds / 60.0;
      if ( minutes > 0 )
	{
	  seconds -= minutes * 60.0;
	}

      
      //
      // write as a new signal
      //
      
      logger << "  adding " << events.size() << " "
	     << anames[a] << " annotations (spanning ";
      if ( minutes > 0 ) logger << minutes << " min " << seconds << " sec)";
      else logger << seconds << " sec)";
      logger << " as 0/1 signal " << labels[a] << "\n";
      
      edf->add_signal( labels[a]  , sr , adat );

    }
      
}


void timeline_t::list_spanning_annotations( const param_t & param )
{
    
  if ( mask_set ) 
    Helper::halt( "cannot run SPANNING with a MASK set... use RE" );
  
  // currently, SPANNING only for continuous EDFs
  if ( ! edf->header.continuous )
    Helper::halt( "currently, can only run SPANNING on continuous EDF" );
  

  // given a /set/ of annotations, determine 
  //   - seconds outside of EDF
  //   - total duration of signal covered by these (seconds)
  //   - coverage as a proportion of EDF file
  //   - coverage as a proportion of in-memory representation
  // etc
  

  //
  // which signals: either look at all, or the requested set
  //


  std::vector<std::string> requested = param.has( "annot" ) 
    ? param.strvector( "annot" ) 
    : annotations.names() ;
  

  //
  // Get all annotations (i.e. not stratified by epoch), sort by time and collapse
  //
  

  std::set<instance_idx_t> events;


  //
  // iterate over each annotation
  //

  for (int a = 0 ; a < requested.size() ; a++ ) 
    {
      
      annot_t * annot = annotations.find( requested[a] );
      
      if ( annot == NULL ) continue;
      
      const int num_events = annot->num_interval_events();
      
      //
      // iterator over interval/event map
      //
      
      annot_map_t::const_iterator ii = annot->interval_events.begin();
      while ( ii != annot->interval_events.end() )
	{	  
	  const instance_idx_t & instance_idx = ii->first;	  
	  events.insert( instance_idx );	  
	  ++ii;
	}
      
      
    }


  //
  // track total coverage, etc
  //

  uint64_t total = 0;
  
  uint64_t total_collapsed = 0;
  
  uint64_t invalid_tps = 0;

  int over_extended = 0;


  // keep track of where longest spanning annot reaches to 
  // or 0 if past the previous spanning annot

  uint64_t earliest = 0;

  uint64_t furthest = 0;


  std::set<instance_idx_t>::const_iterator aa = events.begin();
  while ( aa != events.end() )
    {
      
      const interval_t & interval = aa->interval;
      
      //
      // track total (uncollapsed) duration across all ANNOTs
      //

      total += interval.duration();
      
      // 
      // what overlap, if any?
      //
      
      uint64_t vtp = valid_tps( interval );
      
      bool is_valid = interval.duration() == vtp;

      if ( ! is_valid ) 
	{

	  // duration of annots that do not map to a EDF region
	  invalid_tps +=  interval.duration() - vtp; 
	  
	  // count of intervals that do not perfectly match valid regions
	  ++over_extended;
	  
	  // report
	  writer.level( over_extended , globals::count_strat );

	  writer.value( "ANNOT" , aa->parent->name );
	  writer.value( "INST" , aa->id );
	  writer.value( "START" , interval.start_sec() );
	  writer.value( "STOP" , interval.stop_sec() );
	  writer.unlevel( globals::count_strat );
	}
      
          
      //
      // track collapsed duration, but here only consider completely 'valid' intervals
      //

      if ( is_valid ) 
	{

	  // start of a 'new' region?
	  
	  if ( furthest == 0 ) 
	    {
	      earliest = interval.start;
	      furthest = interval.stop;
	    }
	  else // we already have at least one region counted
	    {
	      
	      // is the old region finished?  if so, add
	      if ( interval.start > furthest ) 
		{
		  total_collapsed += furthest - earliest ;
		  earliest = interval.start;
		  furthest = interval.stop;
		}
	      else // add to current region 
		{
		  if ( interval.stop > furthest )
		    furthest = interval.stop;
		}
	      
	    }
	 	  
	}
	  
      // next segment
      ++aa;
      
    }
  
  // add final interval(s)
  
  total_collapsed +=  furthest - earliest ;
  

  //
  // Report
  //
  

  writer.value( "REC_SEC" , Helper::tp2sec( total_duration_tp )  );
  writer.value( "REC_HMS" , Helper::timestring( total_duration_tp , ':' )  );

  writer.value( "ANNOT_N" , (int)events.size() );
  writer.value( "ANNOT_SEC" , Helper::tp2sec( total )  );
  writer.value( "ANNOT_HMS" , Helper::timestring( total , ':' )  );

  // do any (valid) annots overlap each other?
  writer.value( "ANNOT_OVERLAP" , total_collapsed < total );

  // how many annots over-extended beyond range of EDF?
  writer.value( "INVALID_N" , over_extended );
  writer.value( "VALID_N" , (int)events.size() - over_extended );

  // extent of this over-extension
  writer.value( "INVALID_SEC" , Helper::tp2sec( invalid_tps ) );
  
  writer.value( "SPANNED_PCT" , 100 * ( Helper::tp2sec( total_collapsed ) / Helper::tp2sec( total_duration_tp ) )  );
  writer.value( "SPANNED_SEC" , Helper::tp2sec( total_collapsed )  );
  writer.value( "SPANNED_HMS" , Helper::timestring( total_collapsed , ':' ) );

  writer.value( "UNSPANNED_SEC" , Helper::tp2sec( total_duration_tp - total_collapsed ) );
  writer.value( "UNSPANNED_PCT" , 100 * ( 1 - Helper::tp2sec( total_collapsed ) / Helper::tp2sec( total_duration_tp ) )  );
  writer.value( "UNSPANNED_HMS" , Helper::timestring( total_duration_tp - total_collapsed , ':' ) ); 
    
}


void timeline_t::list_all_annotations( const param_t & param )
{

  //
  // Options
  //

  // count annotations per epoch
  bool per_epoch = param.has( "epoch" );

  // do this either way, as EDF+ mode requires epochs to locate annots
  if ( ! epoched() ) 
    {
      int ne = set_epoch( globals::default_epoch_len , globals::default_epoch_len );
      logger << "  set epochs to default " << globals::default_epoch_len << " seconds, " << ne << " epochs\n";
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
		  if ( is_masked && ! show_masked ) { ++ii; continue; } 
		  
		  // else display
		  writer.level( instance_idx.id , "INST" );
		  writer.level( interval.as_string() , "INTERVAL" );
		  writer.level( instance_idx.ch_str , globals::signal_strat );

		  writer.value( "EMASK" , masked( e ) );
		  writer.value( "AMASK" , is_masked );
		  
		  
		  ++ii;
		}      

	      writer.unlevel( "INTERVAL" );
	      writer.unlevel( "INST" );
	      writer.unlevel( globals::signal_strat );

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

	  std::map<std::string, globals::atype_t>::const_iterator tt = annot->types.begin();
	  while ( tt != annot->types.end() )
	    {
	      std::cout << "  " << tt->first << ", is " << globals::type_name[ tt->second ] << "\n";
	      ++tt;
	    }
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

  //
  // report HMS?
  //
  
  clocktime_t starttime( edf->header.starttime );
  bool hms = true;  
  if ( ! starttime.valid )
    {
      logger << " ** could not find valid start-time in EDF header **\n";
      hms = false;
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

      // NOTE: not sure why we previously did this... for output only, keep consistent form
      // do not +1 time-unit
      //writer.value( "STOP" , interval.stop_sec_exact() );

      writer.value( "STOP" , interval.stop_sec() );

      // channel label
      writer.value( "CH" , instance_idx.ch_str );

      // HMS : elapsed
      // HMS : clock

      if ( hms )
	{
	  
	  double tp1_sec = interval.start_sec();
	  clocktime_t present1 = starttime;
	  present1.advance( tp1_sec / 3600.0 );
	  // add down to 1/100th of a second
	  double tp1_extra = tp1_sec - (long)tp1_sec;

	  // Not sure why we used this form previously... to be consistent, stick with STOP being +1 end
	  
	  // stop_sec_exact() return last time point (rather than usual 1-past-the-end)
	  //double tp2_sec = interval.stop_sec_exact();

	  double tp2_sec = interval.stop_sec();

	  clocktime_t present2 = starttime;
	  present2.advance( tp2_sec / 3600.0 );
	  double tp2_extra = tp2_sec - (long)tp2_sec;
	   
	  writer.value( "START_HMS"  , present1.as_string(':') +  Helper::dbl2str_fixed( tp1_extra , 4  ).substr(1) );
	  writer.value( "STOP_HMS"   , present2.as_string(':') +  Helper::dbl2str_fixed( tp2_extra , 4  ).substr(1) );

	  // elapsed time (00:00:00 is start of EDF)
	  clocktime_t present3;
	  present3.advance( tp1_sec / 3600.0 );
	  // add down to 1/100th of a second
	  tp1_extra = tp1_sec - (long)tp1_sec;

	  clocktime_t present4;
	  present4.advance( tp2_sec / 3600.0 );
	  tp2_extra = tp2_sec - (long)tp2_sec;
	   
	  writer.value( "START_ELAPSED_HMS"  , present3.as_string(':') +  Helper::dbl2str_fixed( tp1_extra , 4  ).substr(1) );
	  writer.value( "STOP_ELAPSED_HMS"   , present4.as_string(':') +  Helper::dbl2str_fixed( tp2_extra , 4  ).substr(1) );


	}
	

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
void timeline_t::apply_eval_mask( const std::string & str , int mask_mode , const bool verbose )
{

  // mask_mode   0   mask
  //             1   unmask
  //             2   force  (T = mask    &  F = unmask)   [ drop mode ] 
  //            -2   force  (T = unmask  &  F = mask)     [ keep mode ]
  

  // is this 'KEEP' mode? 

  bool flip = false;
  
  if ( mask_mode == - 2 ) { mask_mode = 2 ; flip = true; } 


  //
  // Set mask mode
  //

  if ( mask_mode > -1 ) 
    {
      set_epoch_mask_mode( mask_mode );  
      logger << " set masking mode to " << ( mask_mode == 2 ? "'force'" : mask_mode == 1 ? "'unmask'" : "'mask' (default)" ) << "\n";
    }


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
  int cnt_basic_match = 0;  

  
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

      bool is_valid = tok.evaluate( verbose );
      
      bool matches;
      
      if ( ! tok.value( matches ) ) is_valid = false;



      //
      // Flip?
      //

      if ( flip ) matches = ! matches;
            
      
      //
      // A match must be a valid value
      //
      
      if ( ! is_valid ) matches = false;


      //
      // apply mask (or not)
      //
      
      acc_total++;

      acc_valid += is_valid;
      
      //
      // Only for valud results
      //

      if ( is_valid ) 	
	{
	  acc_retval += matches;
      
	  
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
	}
      else ++cnt_unchanged;

      // count current mask state

      if ( !mask[e] ) ++cnt_now_unmasked;
      
      // next epoch
    } 
  
  
  logger << " based on eval expression [" << expression << "]\n"
	 << "  " << acc_retval << " true, " << acc_valid - acc_retval << " false and " 
	 << acc_total - acc_valid << " invalid return values\n"
	 << "  " << cnt_basic_match << " epochs match; " 
	 << cnt_mask_set << " newly masked, "
	 << cnt_mask_unset << " unmasked, "
	 << cnt_unchanged << " unchanged\n";
  logger << "  total of " << cnt_now_unmasked << " of " << epochs.size() << " retained\n";

  
  // mask, # epochs masked, # epochs unmasked, # unchanged, # total masked , # total epochs
  
  writer.level( expression , "EMASK" );
  
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
  writer.value( "N_TOTAL"      , (int)epochs.size()    );

  writer.unlevel( "EMASK" );


  // all done 
  return;

}


int timeline_t::whole_recording_epoch_dur() {
  // only allow this for a continuous EDF 
  if ( ! edf->header.continuous ) return 0;    
  return floor( edf->header.nr * edf->header.record_duration_tp * globals::tp_duration );
}

