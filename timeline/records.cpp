
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
#include "db/db.h"
#include "helper/logger.h"
#include <iterator>     // std::distance

extern writer_t writer;
extern logger_t logger;


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

interval_t timeline_t::record2interval( int r ) const
{ 
  std::map<int,uint64_t>::const_iterator ll = rec2tp.find(r);
  std::map<int,uint64_t>::const_iterator uu = rec2tp_end.find(r);
  if ( ll == rec2tp.end() ) return interval_t(0,0);
  return interval_t( ll->second , uu->second );
}

bool timeline_t::remap_timepoint( const uint64_t & tp , uint64_t * tp1 ) 
{

  
  // note - this should be called w/ all stop timepoints tp = tp--;

  // std::cout << "\n\n IN remap_timepoint() - to map tp = " << tp << "\n";

  // std::map<uint64_t,int>::const_iterator qq = tp2rec.begin();
  // while ( qq != tp2rec.end() )
  //   {
  //     std::cout << "dets: " << qq->first << "\t"
  //   		<< qq->second << "\t"
  //   		<< rec2orig_rec[ qq->second ] << "\n";
  //     ++qq;
  //   }
  
  
  //
  // find an exact mapping of a single timepoint
  // using tp2rec to count records
  //

  // note - this should be called w/ all stop timepoints tp = tp--;
    
  

  // degenerate case: nothing to do w/ a continuous recording
  if ( edf->header.continuous )
    {
      *tp1 = tp;
      return true;
    }

  // for the EDF+D case, get first record that is equal to greater than
  
  if ( tp2rec.size() == 0 ) return false;
  
  std::map<uint64_t,int>::const_iterator pp = tp2rec.lower_bound( tp ); 


  // if ( pp != tp2rec.end() )
  //    std::cout << " initial = " << pp->first << " " << pp->second << "\n";
  
  int recnum = -1;
  uint64_t offset = 0LLU;

  // point after last point?
  if ( pp == tp2rec.end() )
    {
      // get the last record
      --pp;
      
      //std::cout << " after start of last record...\n";
      
      // but is tp within the last record?
      if ( tp <= pp->first + edf->header.record_duration_tp )
	{
	  //std::cout << " but wihin rec...\n";
	  if ( tp < pp->first ) Helper::halt( "internal logic error in remap_timepoint()" );	  
	  offset = tp - pp->first ;
	  //std::cout << " offset into this rec = " << offset << "\n";
	  recnum = rec2orig_rec[ pp->second ];
	  //std::cout << "recnum = " << recnum << " given " << pp->second << "\n";
	  *tp1 = recnum * edf->header.record_duration_tp + offset ;
	  return true;
	}
      else
	return false;
    }


  //     recs     |------|    |----|----|----|   
  //  tp     X       Y      X              Y   X
  
  bool in_gap = false;

  //   std::cout << " tp = " << tp << " " << pp->first << "\n";

  // exact match?
  if ( tp == pp->first )
    {
      offset = 0LLU;
      recnum = rec2orig_rec[ pp->second ];

      // 	std::cout << " ex match recnum " << recnum << " " << offset << "\n";
      
      // set tp1  
      *tp1 = recnum * edf->header.record_duration_tp ; 
      return true;
    }
  
  if ( pp != tp2rec.begin() ) 
    {
      // go back one record
      --pp;
      uint64_t previous_rec_start = pp->first;
      uint64_t previous_rec_end   = previous_rec_start + edf->header.record_duration_tp - 1LLU ;
      
      // 	std::cout << " previous_rec_start = " << previous_rec_start << " " << previous_rec_end << "\n";
      
      // does the start point fall within this previous record?
      if ( tp >= previous_rec_start && tp <= previous_rec_end ) 
	{
	  //   std::cout << " not in GAP\n";
	  in_gap = false;
	}
      else
	{
	  //	  if ( verb ) std::cout << "gapper\n";
	  // is in a gap
	  return false;
	}
    }
  else if ( pp == tp2rec.begin() )
    {
      
      // If the search point occurs before /all/ records, need to
      // indicate that we are in a gap also	  
      
      if ( tp  < pp->first )
	return false;
      
    }

  // if here, we are not in a gap, so set tp1

  offset = tp - pp->first ;
  recnum = rec2orig_rec[ pp->second ];	  
  //     std::cout << " recnum " << recnum << " " << offset << "\n";
  // set tp1
  *tp1 = recnum * edf->header.record_duration_tp + offset ;
  
  return true;
}

bool timeline_t::interval2records( const interval_t & interval , 
				   uint64_t n_samples_per_record , 
				   int * start_rec , 
				   int * start_smp , 
				   int * stop_rec , 
				   int * stop_smp ) const

{

  // if n_samples_per_record == 0 implies the SR is actually time-point resoluation
  //  i.e. SR = 1e9
  // almost certainly not necessary, but below we avoid some divisions etc in this special case
  
  // if interval starts inbetween gaps, effectively shift it back to the one just before
  //  i.e. and similarly shift the end point
  //  for continunous cases, this means that all segments with same interval size will return
  //  the same number of samples;  otherwise, prior span rules could lead to +1 additional sample
  //  e.g. if SR=100Hz, pulling interval 0.015 seconds:

  //   S   0       1       2        ...
  //       0.000   0.010   0.020 

  //  as intervals are defined [a,b), this means that under old encoding: 

  //  T                  S
  //  0.0000 0.0150  ->  [ 0 , 1 )    L=1
  //  0.0025 0.0175  ->  [ 0 , 1 )    L=1
  //  0.0050 0.0200  ->  [ 0 , 2 )    L=2
  //  0.0075 0.0225  ->  [ 0 , 2 )    L=2
  //  0.0100 0.0250  ->  [ 1 , 2 )    L=1


  // now: 
  //  T                                     S
  //  0.0000 0.0150                     ->  [ 0 , 1 )    L=1
  //  0.0025 0.0175  ->  0.0000 0.0150  ->  [ 0 , 1 )    L=1
  //  0.0050 0.0200  ->  0.0000 0.0150  ->  [ 0 , 1 )    L=1
  //  0.0075 0.0225  ->  0.0000 0.0150  ->  [ 0 , 1 )    L=1
  //  0.0100 0.0250                     ->  [ 1 , 2 )    L=1

  // for EDF+D cases (i.e. w/ gaps which might be of any length) it is not guaranteed that a similarly-sized
  // pair of intervals will return the same number of samples, if the intervals span gaps, and so we don't
  // case.    But at least for epochs, every epoch is defined to be a contiguous set of samples, and so even
  // for EDF+D this ensures that we won't get slices returned w/ +1 extra sample point ever in epoch-based
  // analyses

  // for most cases (i.e. where intervals are aligned with records/samples, as is typical case, this was never
  // and issue;  but this provides a nicer handling of fractional epochs, etc.   The primary impact was for TLOCK
  // where we'd get one additional sample at the end depending on alignment as above; although not really a substantive
  // problem, this led to data-handling issues.
  
   // std::cout << " search: " << interval.as_string() << "\n";
   // std::cout << " search (tp): " << interval.start << " " << interval.stop << "\n";
   // std::cout << " n-samples per rec: " << n_samples_per_record << "\n";
   // std::cout << " EDF contin: " <<  edf->header.continuous << "\n";

  //
  // badly-defined interval?
  //
  
  if ( interval.stop < interval.start ) 
    Helper::halt( "internal error: badly defined interval requested, with stop before start" );
  
  //
  // 0-duration interval: i.e. an empty record set
  //
  
  if ( interval.stop == interval.start )
    {
      *start_rec = 0;
      *start_smp = 0;	  
      *stop_rec = 0;
      *stop_smp = 0;
      return false;
    }
  
  //
  // Otherwise, requires that stop is >0; should be the case given the
  // above checks, but whatever...
  //
  
  if ( interval.stop == 0 ) 
    Helper::halt( "internal error in timeline()" );
    
  // sample rate, in tp units: special case, SR=1e-9
  uint64_t sample_tp = n_samples_per_record == 0 ? 1LLU : edf->header.record_duration_tp / n_samples_per_record ; 
  
  // pull out stop as well may need to adjust this below
  uint64_t stop_tp = interval.stop;
  
  // ensure not past end  (last_time_point is the last accessible time-point, so go +1 on this)
  if ( stop_tp > edf->timeline.last_time_point_tp )
    stop_tp = edf->timeline.last_time_point_tp + 1LLU ;
    
  //
  // For a continuous timeline, given time-points can
  // straightforwardly calculate record/sample
  //

  if ( edf->header.continuous )
    {

      // get start record 
      uint64_t start_record = interval.start / edf->header.record_duration_tp;
      
      // tp-offset into this record
      uint64_t start_offset = interval.start % edf->header.record_duration_tp;

      // get sample point at or just prior (floor). nb special case SR = 1e9 Hz -> value = 0
      uint64_t start_sample = n_samples_per_record == 0 ? start_offset :
       	floor( ( start_offset / (double)edf->header.record_duration_tp ) * n_samples_per_record ) ;

      // if interval started after this sample (i.e. started between samples), then shift whole
      // interval (including stop) back too)
      uint64_t shift = start_offset - start_sample * sample_tp ; 
      
      // for intervals starting between exact samples, shift whole interval back to start at the
      // sample just before

      // this should not happen, just in case
      if ( shift > stop_tp )
	Helper::halt( "internal error in interval2records(), with unaligned interval" );

      stop_tp -= shift; 
      
      // get final records/samples, nb. use floor() to get the sample at or just prior to the end
      // record/sample selection returned is *inclusive of end* so base this on stop-1 rather than stop
      
      uint64_t stop_record = stop_tp / edf->header.record_duration_tp;
      uint64_t stop_offset = stop_tp % edf->header.record_duration_tp;
      uint64_t stop_sample = n_samples_per_record == 0 ? stop_offset : 
	floor( ( stop_offset / (double)edf->header.record_duration_tp ) * n_samples_per_record  ) ;


      //      std::cout <<" provisional = " << stop_record << "\t" << stop_sample << " offset = " << stop_offset << "\n";
      //
      // Shift one sample backwards, i.e.
      //   i.e.  1) only include "whole" sample-sample intervals spanned by the search interval
      //         2) search interval stop is defined as 1-past end, and so we don't want to include the
      //            final sample point in any case, if it lands exactly on a sample-point
      //
      
      if ( stop_sample == 0 )
	{
	  --stop_record;

	  // huh, special case in reverse...
	  
	  stop_sample = n_samples_per_record == 0 ? globals::tp_1sec - 1LLU : 
	    n_samples_per_record - 1LLU ;
	}
      else
	--stop_sample;
      
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

	  // interval start is before /all/ records?
	  //   -> flag as a gap
	  
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

      // figure out start sample, and also whether we need to shift
      // for fractional cases (starts between sample points)

      uint64_t shift = 0LLU;
      
      if ( in_gap )
	{
	  // i.e. use start of this record, as it is after the 'true' start site
	  *start_smp = 0; 
	  // shift stays at 0
	}
	else
	{	
	  
	  uint64_t start_offset = interval.start - lwr->first;
		  
	  uint64_t start_sample = n_samples_per_record == 0 ? start_offset :
	    floor( ( start_offset / (double)edf->header.record_duration_tp ) * n_samples_per_record ) ;
	  
	  *start_smp = (int)start_sample;
	  
	  // any shift required?
	  shift = start_offset - start_sample * sample_tp ;
	  
	}
      
      //
      // Is start misaligned with a sample point?   Even though this will not matter necessarily in a
      // gapped case, for when we are working w/ contiguous intervals (the most common) scenario, we
      // want to employ the same logic as above, i.e. to shift the interval implicitly back (including the
      // stop) so that we always get a fixed number of samples (at a given SR) for a given search interval
      // length. 
      //
      
      if ( shift > stop_tp )
        Helper::halt( "internal error in interval2records(), with unaligned interval" );
      
      stop_tp -= shift;


      //
      // Now find stop record/sample
      // For upper bound, find the record whose end is *greater* than the interval stop
      //

      // nb: easier to base this on the last time-point in the interval, rather than +1
      uint64_t stop_tp_m1 = stop_tp > 0 ? stop_tp - 1LLU : 0 ; 
      
      std::map<uint64_t,int>::const_iterator upr = tp2rec.upper_bound( stop_tp_m1 ); 
      

      
      //
      // this should have returned one past the one we are looking for 
      // i.e. that starts *after* the search point, as we used stop rather than stop+1 as usual
      //
      // this handles the case of the search being past the last
      // record, as we skip back one i.e. to last - 1, below
      //

      // special case; if the interval ends before the records even start
      bool ends_before = upr == tp2rec.begin() ;

      if ( ! ends_before ) 
	{
	  --upr;  
	  *stop_rec  = upr->second;
	}
      else
	{
	  // i.e. flag as bad, to ensure that stop is before the start (which will also be rec 0)
	  //*stop_rec  = -1;

	  // just return now, empty set:
	  *start_rec = *start_smp = *stop_rec = *stop_smp = 0;
	  return false;	  
	}

    
      //
      // get samples within (as above)      
      //

      uint64_t previous_rec_start = upr->first;
      uint64_t previous_rec_end   = previous_rec_start + edf->header.record_duration_tp - 1; // inclusive defn.


      //
      // Does this end point fall in a gap? i.e. after the previous record ends
      // Note: stop_tp_m1 is exact end (not +1 )
      //
      
      // assuming stop_tp_m1 == inclusive
      in_gap = ! ( stop_tp_m1 >= previous_rec_start && stop_tp_m1 <= previous_rec_end );
      
      // assuming stop_tp is end+1
      //in_gap = ! ( stop_tp > previous_rec_start && stop_tp < previous_rec_end );
      

      if ( in_gap )
	{
	  // stop falls in gap, so set to last sample of the previous record
	  // handling the special SR=0 case
	  *stop_smp = n_samples_per_record == 0 ? globals::tp_1sec - 1LLU :
	    n_samples_per_record - 1LLU ;	  
	}
      else
	{
	  
	  // else, determine the offset into the prior record for stop_tp
	  
	  // how far into this record (TP)? nb. here using end+1 time-point as above for continuous case
	  uint64_t stop_offset = stop_tp - upr->first ; 
	  
	  // convert to sample points
	  uint64_t stop_sample = n_samples_per_record == 0 ? stop_offset : 
	    floor( ( stop_offset / (double)edf->header.record_duration_tp ) * n_samples_per_record ) ;
	  
	  *stop_smp = (int)stop_sample;
	}

      //
      // Shift one sample back, following the same logic as above for
      // the continuous case
      //

      if ( *stop_smp == 0 )
        {
	  if ( *stop_rec == 0 )
	    {
	      *start_rec = *start_smp = *stop_rec = *stop_smp = 0;
	      return false;
	    }
	  
          *stop_rec = *stop_rec - 1;  

          *stop_smp = n_samples_per_record ? globals::tp_1sec - 1LLU :
	    n_samples_per_record - 1LLU;
        }
      else
	*stop_smp = *stop_smp - 1 ; 
      
    }

  
  // std::cout << "recs = " << *start_rec << " " << *stop_rec << "\n";
  // std::cout << "smps = " << *start_smp << " " << *stop_smp << "\n";

  //
  // old code handling scenarios of short intervals (e.g. annotations)
  // that might be 0-duration and/or fall between gaps.   Given new
  // align-to-sample logic, not sure we need these any longer.  Thus,
  // remove for now.
  //

  // if ( *start_rec == *stop_rec && *start_smp == *stop_smp + 1 )
  //   {
  //     *stop_smp = *start_smp;
  //   }

  // NOTE - if revisiting this below, would need to change n_samples_per_record to
  // handle the special case of 0
  
  // else if ( *start_rec == *stop_rec + 1LLU
  //  	    && *stop_smp == 0 // first point of new record
  //  	    && *start_smp == n_samples_per_record - 1 ) // last point of previous record
  //    {
  //      *stop_rec = *start_rec;
  //      *stop_smp = *start_smp;
  //    }
  

  
  //
  // check for a gap/problem: the 'stop' might come before 'start' in this case,
  // so just return as an empty set;   this can happen is the interval starts and
  // stops entirely within a gap between two sample points
  //
  
  if ( *start_rec > *stop_rec || ( *start_rec == *stop_rec && *start_smp > *stop_smp ) )
    {
      *start_rec = *start_smp = *stop_rec = *stop_smp = 0;
      return false;
    }
  
  //
  // All done
  //

    
  return true;
  
}

interval_t timeline_t::collapse( const interval_t & interval ) const
{
  int start_rec = 0 , stop_rec = 0;
  int start_smp = 0 , stop_smp = 0;
  
  // to get 1/100,000 second resolution
  const int srate = 100000;  
  
  bool any = interval2records( interval , srate , &start_rec , &start_smp , &stop_rec , &stop_smp );

  //  std::cout << " start rec smp = " << start_rec << " " << start_smp << "\n";
  //std::cout << " stop rec smp = " << stop_rec << " " << stop_smp << "\n";
  // std::cout << " any = " << any << "\n";
  
  // interval has to fall completely in a valid area
  if ( ! any ) return interval_t( 1LLU , 0LLU );

  if ( rec2orig_rec.find( start_rec ) == rec2orig_rec.end() ) return interval_t( 1LLU , 0LLU );
  if ( rec2orig_rec.find( stop_rec ) == rec2orig_rec.end() ) return interval_t( 1LLU , 0LLU );

  start_rec = rec2orig_rec.find( start_rec )->second;
  stop_rec = rec2orig_rec.find( stop_rec )->second;
  
  uint64_t start = start_rec * edf->header.record_duration_tp + ( start_smp / (double)srate ) * globals::tp_1sec ;

  // nb + add one extra STOP smp here, +1 end
  uint64_t stop = stop_rec * edf->header.record_duration_tp + ( (stop_smp+1) / (double)srate ) * globals::tp_1sec ;

  // std::cout << " start = " << interval.start / globals::tp_1sec  << " --> " << interval.stop / globals::tp_1sec << "\n";
  // std::cout << " diff = " << ( interval.start - start ) / globals::tp_1sec << "\t" <<  ( interval.stop - stop ) /	globals::tp_1sec << "\n";
  // std::cout << " XXart = " << start / globals::tp_1sec  << " --> " << stop / globals::tp_1sec << "\n";
  
  return interval_t( start , stop );
    
}


std::set<int> timeline_t::records_in_interval( const interval_t & interval ) const
{
  
  int start_rec = 0 , stop_rec = 0;
  int start_smp = 0 , stop_smp = 0;

  //  const int srate = 100;  // will not matter, as we only consider whole records here

  // set to 1tp-unit for accurate lookups
  // (we do only select whole records, but can be edge cases when using certain epoch
  // offsets (only function that calls records_in_interval() ) if SR == 100 or
  // some such, i.e. as the interval will align to sample points and be shifted, meaning
  // that it was possible to miss the last record that would otherwise have been spanned by
  // the interval.   Setting 'infinite' sample rate (i.e. 1e=9 = time-point resolution)
  // fixes this.
  
  // const uint64_t srate = edf->header.record_duration_tp;
  const uint64_t srate = 0LLU;  // setting SR implies the above for interval2records()

  std::set<int> recs;

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

  // an EDF record is 'masked' if it is out-of-bounds
  // of if *any* of the spanning epochs are masked
  //     REC   |-------| 
  //    E1111|E22222|E33333|E44444
  //  e.g. 

  if ( ! mask_set ) return false;
  
  std::map<int,std::set<int> >::const_iterator rr = rec2epoch.find(r);
  if ( rr == rec2epoch.end() ) return true; // i.e. out of bounds
  

  // OLD behavior : MASKED if spanned by ANY masked epoch
  // const std::set<int> & epochs = rr->second;
  // std::set<int>::const_iterator ee = epochs.begin();
  // while ( ee != epochs.end() )
  //   {
  //     if ( mask[ *ee ] ) return true;
  //     ++ee;
  //   }
  // return false;

  // NEW behaviour : MASKED if not spanned by ANY unmasked epoch
  const std::set<int> & epochs = rr->second;
  std::set<int>::const_iterator ee = epochs.begin();
  while ( ee != epochs.end() )
    {
      if ( ! mask[ *ee ] ) return false;
      ++ee;
    }
  return true;


  
}
