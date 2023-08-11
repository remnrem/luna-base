
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

  //  std::cout << "\n\n IN remap_timepoint() - to map tp = " << tp << "\n";

  // std::map<uint64_t,int>::const_iterator qq = tp2rec.begin();
  // while ( qq != tp2rec.end() )
  //   {
  //     std::cout << "dets: " << qq->first << "\t"
  //  		<< qq->second << "\t"
  //  		<< rec2orig_rec[ qq->second ] << "\n";
  //     ++qq;
  //   }

  
  // find an exact mapping of a single timepoint
  // using tp2rec to count records

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
  //   std::cout << " initial = " << pp->first << " " << pp->second << "\n";

  int recnum = -1;
  uint64_t offset = 0LLU;

  // point after last point?
  if ( pp == tp2rec.end() )
    {
      // get the last record
      --pp;
      
      // but is tp within the last record?
      if ( tp <= pp->first + edf->header.record_duration_tp )
	{
	  if ( tp < pp->first ) Helper::halt( "internal logic error in remap_timepoint()" );	  
	  offset = tp - pp->first ;
	  recnum = rec2orig_rec[ pp->second ];	  
	  *tp1 = recnum * edf->header.record_duration_tp + offset ;
	  return true;
	}
      else
	return false;
    }


  //     recs     |------|    |----|----|----|   
  //  tp     X       Y      X              Y   X
  
  bool in_gap = false;

  //  std::cout << " tp = " << tp << " " << pp->first << "\n";
  // exact match?
  if ( tp == pp->first )
    {
      offset = 0LLU;
      recnum = rec2orig_rec[ pp->second ];
      //      std::cout << " ex match recnum " << recnum << " " << offset << "\n";
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
      
      //      std::cout << " previous_rec_start = " << previous_rec_start << " " << previous_rec_end << "\n";
      
      // does the start point fall within this previous record?
      if ( tp >= previous_rec_start && tp <= previous_rec_end ) 
	{
	  //std::cout << " not in GAP\n";
	  in_gap = false;
	}
      else
	{
	  //std::cout << "gapper\n";
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
  //  std::cout << " recnum " << recnum << " " << offset << "\n";
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

  // std::cout << " search: " << interval.as_string() << "\n";
  // std::cout << " search (tp): " << interval.start << " " << interval.stop << "\n";
  // std::cout << " n-samples per rec: " << n_samples_per_record << "\n";
  // std::cout << " EDF contin: " <<  edf->header.continuous << "\n";
  
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

  //
  // Note: here we want to find records/samples that are inclusive w.r.t. the interval
  // so change coding of stop being 1 unit past the end below
  //

  if ( interval.stop == 0 ) 
    Helper::halt( "internal error in timeline()" );

  //
  // Interval defines one-past-the end;   we therefore subtract 1-tp unit from it;
  //  Hpwever, if matching fractional intervals (i.e. interval.start and .stop do not
  //  fall exactly on a sample point, then we will get 1 extra sample point;  we therefore
  //  track both stop+1 and stop, and if they give the same implied sample point, we subtract one
  //
  //  e.g.  if 100 Hz sample
  //   okay     0.02    0.92
  //                    = 9200000000 is one past  (of whatever, 1e9 tp = 1sec)
  //                    = 9199999999 is exact end
  //                    = means we get samples #3 (0.02) to #92 (0.91)

  //  not okay: 0.025   0.925
  //                   = 9255000000
  //                   = 9254999999
  //                   = implies #3 to #93 (0.92) .. and so one extra; therefore, here subtract 1
  //   

  uint64_t stop_plus1_tp = interval.stop;
  
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
      //std::cout << "EDF-C\n";
      
      // old version:  get initial records/samples, nb. use ceil() to get nearest sample *after* start of interval
      
      // now (v0.25+) for fractionally split things (i.e. if interval boundary does not align w/ a sample point exactly)
      // for both start and end, take the point prior to the fractional point;  this should preserve the number of samples
      // selected, i.e. in ALIGN
      
      uint64_t start_record = interval.start / edf->header.record_duration_tp;
      uint64_t start_offset = interval.start % edf->header.record_duration_tp;

      // prior to v0.25
      // uint64_t start_sample = 
      // 	ceil( ( start_offset / (double)edf->header.record_duration_tp ) * n_samples_per_record ) ;
      
      // change in v0.25 

      uint64_t start_sample = 
       	floor( ( start_offset / (double)edf->header.record_duration_tp ) * n_samples_per_record ) ;
      
      // std::cout << "othr = " << edf->header.record_duration_tp << " " << n_samples_per_record << "\n";
      // std::cout << "start = " << start_record << " " << start_offset << " " << start_sample << "\n";
      
      // get final records/samples, nb. use floor() to get the nearest sample *prior* to end
      
      uint64_t stop_record = stop_tp / edf->header.record_duration_tp;
      uint64_t stop_offset = stop_tp % edf->header.record_duration_tp;
      uint64_t stop_sample = 
	floor( ( stop_offset / (double)edf->header.record_duration_tp ) * n_samples_per_record  ) ;

      uint64_t stop_plus1_record = stop_plus1_tp / edf->header.record_duration_tp;
      uint64_t stop_plus1_offset = stop_plus1_tp % edf->header.record_duration_tp;
      uint64_t stop_plus1_sample = 
	floor( ( stop_plus1_offset / (double)edf->header.record_duration_tp ) * n_samples_per_record  ) ;

      // move one sample back (i.e. if the spanning interval did not exactly land on a sample point)
      // fixed: added 'record' check in instance of SR == 1 Hz
      if ( stop_sample == stop_plus1_sample && stop_record == stop_plus1_record )  
	{	
	  if ( stop_sample == 0 )
	    {
	      --stop_record;
	      stop_sample = n_samples_per_record - 1;
	    }
	  else
	    --stop_sample;	  
	}
      
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
      
      //      std::cout << " -- EDF-D \n";
      
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
	  //	  std::cout << "lwr != begin\n";
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
	  //std::cout << "lwr = begin\n";
	  // If the search point occurs before /all/ records, need to
	  // indicate that we are in a gap also	  
	  
	  if ( interval.start < lwr->first ) 
	    {
	      in_gap = true;	      
	      //std::cout << "start in gap\n";
	    }
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
	{
	  //std::cout << "  getting start ... after gap\n";
	  *start_smp = 0; // i.e. use start of this record, as it is after the 'true' start site
	}
	else
	{	
	  // OLD code: broken for EDF+D w/ gaps that are not multiples of the record size
	  //uint64_t start_offset = interval.start % edf->header.record_duration_tp;


	  // NEW version for EDF+D
	  uint64_t start_offset = interval.start - lwr->first;
	  
	  //std::cout << "tp start = " << lwr->first << " " << interval.start << "\n";
	  
	  // nb: old code prior to v0.25
	  // uint64_t start_sample = 
	  //   ceil( ( start_offset / (double)edf->header.record_duration_tp ) * n_samples_per_record ) ;

	  // nb: new version: floor() not ceil() , same behavior as for the end-points
	  uint64_t start_sample = 
	    floor( ( start_offset / (double)edf->header.record_duration_tp ) * n_samples_per_record ) ;
	  
	  //if ( start_sample >= n_samples_per_record ) start_sample = n_samples_per_record - 1LLU; 
	  *start_smp = (int)start_sample;
	}
      
      

      //
      // for upper bound, find the record whose end is equal/greater *greater* 
      // 
      
      std::map<uint64_t,int>::const_iterator upr = tp2rec.upper_bound( stop_tp ); 

      //if ( upr == tp2rec.end() ) std::cout << " AT END\n";
	
      //
      // this should have returned one past the one we are looking for 
      // i.e. that starts *after* the search point
      //
      
      bool ends_before = upr == tp2rec.begin() ;

      //std::cout << " ends before: " << ends_before << "\n";
      
      if ( ! ends_before ) 
	{
	  --upr;  
	  *stop_rec  = upr->second;
	}
      else
	{
	  *stop_rec  = -1; // i.e. flag as bad, to ensure that stop is before the start (which will also be rec 0)
	}

      //std::cout << "stop_rec = " << upr->second << "\n";

      //
      // get samples within (as above)      
      //

      uint64_t previous_rec_start = upr->first;
      uint64_t previous_rec_end   = previous_rec_start + edf->header.record_duration_tp - 1;


      //
      // Does this end point fall in a gap?
      //

      in_gap = ! ( stop_tp >= previous_rec_start && stop_tp <= previous_rec_end );

      // std::cout << " previous_rec_start: " << previous_rec_start << "\n";
      // std::cout << " previous_rec_end: " << previous_rec_end << "\n";

      // if so, set to the last point
      
      if ( in_gap )
	*stop_smp = n_samples_per_record - 1; 
      else
	{

	  // else, determine the offset into this record

	  // std::cout << "stop tp: " << stop_tp << "\n";
	  // std::cout << "edf->header.record_duration_tp: " << edf->header.record_duration_tp << "\n";
	  
	  // OLD code::: broken when EDF+D has gaps that are not multiples of record size...
	  //uint64_t stop_offset = stop_tp % edf->header.record_duration_tp;

	  // how far into this record (TP)?
	  uint64_t stop_offset = stop_tp - upr->first ; 

	  // convert to sample points
	  uint64_t stop_sample = 
	    floor( ( stop_offset / (double)edf->header.record_duration_tp ) * n_samples_per_record ) ;
	  
	  //if ( stop_sample >= n_samples_per_record ) stop_sample = n_samples_per_record - 1LLU;
	  *stop_smp = (int)stop_sample;
	}
      
    }


  //
  // ?TODO: potentially, also require fix for the case of fractional spanning intervals in the case of an EDF+C?
  //
  
  
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

interval_t timeline_t::collapse( const interval_t & interval ) const
{
  int start_rec = 0 , stop_rec = 0;
  int start_smp = 0 , stop_smp = 0;
  
  // to get 1/100,000 second resolution
  const int srate = 100000;  
  
  bool any = interval2records( interval , srate , &start_rec , &start_smp , &stop_rec , &stop_smp );

  // std::cout << " start rec smp = " << start_rec << " " << start_smp << "\n";
  // std::cout << " stop rec smp = " << stop_rec << " " << stop_smp << "\n";
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

  return interval_t( start , stop );
    
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
