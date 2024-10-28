
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
#include "edf/slice.h"
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

  // observed difference between sample points (in tp-units)
  const uint64_t observed = t[sp2] - t[sp1];
  
  // expected difference between sample points, given sample rate and nominal number
  // of sample points between t1 and t2
  
  const uint64_t one_sample_tp = globals::tp_1sec / sr;
  
  const uint64_t expected = one_sample_tp * ( sp2-sp1 );
  
  // are these ~equal?
  const uint64_t diff = observed > expected ? observed - expected : expected - observed ;

  // test to within resolution of half 1/SR 
  return diff > one_sample_tp / 2 ;

}


void timeline_t::init_timeline( bool okay_to_reinit ) 
{
  
  if ( rec2tp.size() != 0 && ! okay_to_reinit ) 
    Helper::halt( "internal error: cannot re-init timeline" );
  
  tp2rec.clear();
  rec2tp.clear();
  rec2tp_end.clear();
  
  rec2orig_rec.clear();
  
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
	  rec2orig_rec[ r ] = r; // 1-to-1 mapping in the continuous case
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
	  
	  // this will be cached if we've done a 'preload'
	  uint64_t tp = edf->timepoint_from_EDF(r);
	  
	  tp2rec[tp] = r;
	  rec2tp[r] = tp;
	  rec2orig_rec[ r ] = r; // 1-to-1 mapping before any RE
	  rec2tp_end[r] = last_time_point_tp = tp + edf->header.record_duration_tp - 1LLU;
	  // last_time_point_tp will be updated, 
	  // and end up being thelast (i.e. record nr-1).
	}
    }

}



void timeline_t::restructure( const std::set<int> & keep , const bool preserve_cache )
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

  // set rec2orig_rec
  int cnt = 0;
  std::map<int,uint64_t>::const_iterator rr = rec2tp.begin();
  while ( rr != rec2tp.end() )
    {
      rec2orig_rec[ rr->first ] = cnt++;
      ++rr;
    }

  // reset epochs (but retain epoch-level annotations)
  reset_epochs();
  
  // clear any cache
  if ( ! preserve_cache ) 
    {
      logger << "  clearing any cached values and recording options\n";
      cache.clear();
      writer.no_cache();
    }
  else
    logger << "  preserving any cached values and recording options\n";
  
  logger << "  retaining " << num_epochs() << " epochs\n";
}


void timeline_t::create_discontinuous_timeline( const std::vector<uint64_t> & tps )
{
  
  // this is only used when making a new merged EDF+D from multiple standard EDFs
  // that contain gaps
  
  // we can assume that header.nr and header.record_duration_tp will have been set 
  
  total_duration_tp = 
    (uint64_t)edf->header.nr * edf->header.record_duration_tp;      
  last_time_point_tp = 0;

  // check
  if ( edf->header.nr != tps.size() )
    Helper::halt( "internal error in timeline_t::create_discontinuous_timeline()" );
  
  // okay to use header.nr here, as this will only be called
  // once, on first creating the new in-memory EDF+D
        
  for (int r = 0;r < edf->header.nr; r++)
	{	  
	  uint64_t tp = tps[r] ;
	  tp2rec[ tp ] = r;
	  rec2tp[ r ] = tp;
	  rec2orig_rec[ r ] = r; // 1-to-1 mapping before any RE
	  rec2tp_end[r] = last_time_point_tp = tp + edf->header.record_duration_tp - 1LLU;
	  // last_time_point_tp will be updated, 
	  // and end up being thelast (i.e. record nr-1).
	}
  
  logger << "  set EDF+D timeline for " << edf->header.nr << " records\n";
}


interval_t timeline_t::wholetrace( const bool silent ) const
{  
  //std::cout << "LTP = " << last_time_point_tp + 1LLU << "\n";
  // end is defined as 1 past the last time point
  
  // check that we don't have a mask set::: if we do, give a warning to the console

  if ( mask_set && ! silent )
    logger << "\n"
	   << "  *** warning - running a command that pulls the whole trace\n"
	   << "  ***           but currently an epoch mask set has been set;\n"
	   << "  ***           for this operation to skip masked epochs,\n"
	   << "  ***           you need to run RE (RESTRUCTURE) beforehand\n";
  
  return interval_t( 0 , last_time_point_tp + 1LLU );
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

      uint64_t tpin = 0;

      while ( rr != records.end() )
	{
	  
	  // start/stop for this record	(as above, does not use 1-past-end encoding here)
	  interval_t rec = record2interval( *rr );
	  
	  // make +1 encoding for record (same as interval)
	  ++rec.stop;
	  
	  // std::cout << " rec " << *rr << "  --> " << rec.start << " " << rec.stop
	  //  	    << "  int " << interval.start << " " << interval.stop << "\n";

	  // REC     |-----------------------|
	  // INT   |-------------------------------|

	  // REC     |-----------------------|
	  // INT          |--------------|

	  // all in?  stop is 1 past last interval in both cases,so use equality test
	  if ( rec.start >= interval.start &&  rec.stop <= interval.stop )
	    {
	      // allow for interval < record size
	      uint64_t tt = interval.duration() < rec.duration() ? interval.duration() : rec.duration() ;
	      tpin += tt;	      
	    }
	  else // otherwise, we have partial overlap 
	    {
	      uint64_t partial = 
		( interval.stop < rec.stop ? interval.stop : rec.stop ) - 
		( interval.start > rec.start ? interval.start : rec.start ) ;	      
	      tpin += partial;
	    }

	  ++rr;
	}
      
      return tpin;
      
    }
  
  return 0;

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




    

int timeline_t::whole_recording_epoch_dur() {
  // only allow this for a continuous EDF 
  if ( ! edf->header.continuous ) return 0;    
  return floor( edf->header.nr * edf->header.record_duration_tp * globals::tp_duration );
}


bool timeline_t::align_epochs( uint64_t * tp , int * rec , const std::set<uint64_t> & annots )
{

  // std::cout << " tp = " << *tp << "\t"
  // 	    << " rec = " << *rec << "\t"
  // 	    << " annots.s() = " << annots.size() << "\n";
  
  std::set<uint64_t>::const_iterator ii = annots.begin();
  while ( ii != annots.end() )
    {
      //      std::cout << "considering annot = " << *ii << "\n";

      // if this annot starts before the record
      if ( *ii < *tp ) { ++ii; continue; }
      //std::cout << " setting tp = " << *tp << "\n";
      *tp = *ii;
      break;
    }

  // is this tp in the same record?  
  // if not, we also need to advance rec
  
  while ( 1 ) 
    {
      interval_t reci = record2interval( *rec );
      // unable to find this record
      if ( reci.start == 0 && reci.stop == 0 )  return false; 
      if ( *tp >= reci.start && *tp <= reci.stop ) return true;
      // advance to the next record and check
      //std::cout << " advancing rec = " << *rec << "\n";
      (*rec)++;
    }

  //  std::cout << " DONE\n";
  
  // all done
  return true;
  
}




