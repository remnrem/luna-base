
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

#include "edf/edf.h"
#include "edf/slice.h"

bool edf_t::align( const std::vector<std::string> & annots )
{

  // 0 ) will be keeping EDF record size as is
  // 1 ) get every time point, and figure out whether it is spanned by a completely included annotation as listed in 'annots'
  // 2 ) this will involves changing all records (i.e. reading the entire record) and then forcing a WRITE, i.e. in the same way 
  //     that RECORD-SIZE does
  // 3 ) we allow records must be whole to be included

  // As different signals may have different sample rates, a potential
  // issue that if splitting an EDF record, there may be different
  // number of EDF samples in a partial record... how to count /
  // ensure all match if allowing an arbitrary start point?  This may be common if 
  //
 
  // conditions;
  // 1) there must not be any overlap (at the sample / time-point level) between annotations in 'annots'
  // 2) instances in 'annots' must be whole multiple of EDF record size: i.e. 30-second or 20 annotataion blocks
  
  // check for overlap:: put all annots in a set
  std::set<interval_t> aset;

  // EDF record size: header.record_duration_tp;
  
  // track implied number of records in the new EDF
  int new_nr = 0;
  
  // track # of skipped annots
  int skipped_disc = 0 , skipped_dur = 0;

  for (int a=0; a<annots.size(); a++)
    {
      annot_t * annot = timeline.annotations.find( annots[ a ] );
      annot_map_t::const_iterator ii = annot->interval_events.begin();
      while ( ii != annot->interval_events.end() )
	{
	  
	  // only consider this annotation is it is completely within a valid time-point
	  // range (i.e. no discontinuities, not out-of-range)
	  uint64_t adur =  ii->first.interval.duration();
	  
	  if ( timeline.valid_tps( ii->first.interval ) != adur ) 
	    {
	      //std::cerr << "  skipping " << ii->first.interval << "\n";
	      ++skipped_disc;
	      ++ii;
	      continue;
	    }
	  
	  // check duration is a multiple of the EDF record duration
	  if ( adur % header.record_duration_tp ) 
	    {
	      //std::cerr << "  skipping, not an exact multiple of EDF recdur: " << ii->first.interval << "\n";
	      ++skipped_dur;
	      ++ii;
	      continue;
	    }

	  // track number of implied new records
	  new_nr += adur / header.record_duration_tp; 

	  aset.insert( ii->first.interval );
	  ++ii;
	}
    }

  logger << "  expecting " << new_nr 
	 << " (of " << header.nr << ") "
	 << " records in the new EDF\n";

  if ( skipped_disc ) logger << "  skipped " << skipped_disc << " annotations that span discontinuities\n";
  if ( skipped_dur ) logger << "  skipped " << skipped_dur << " annotations that do not align with EDF record size\n";

  if ( new_nr == 0 ) 
    {
      logger << "  leaving ALIGN, nothing to do (leaving dataset as is)\n";
      return false;
    }

  // check all annotations for a) overlap and b) an acceptable duration 
  //  (i.e. a multiple of the EDF record size duration)
    
  std::set<interval_t>::const_iterator ii = aset.begin();
  interval_t last = *ii;
  ++ii;
  while ( ii != aset.end() )
    {
      if ( ii->start < last.stop ) 
	Helper::halt( "cannot specify overlapping annotations to EDF re-ALIGN-ment" );
      last = *ii;
      ++ii;
    }
  
  // make new records; and fill them up.
  // with different sample rates, possible the implied times will shift slightly??
  // e.g. if fractional start
  //   EDF record time t  :  sig1 100 samples
  //                         sig2 5   samples
  //   if we start a new record at t/0.25, this will copy 75 samples from sig1, but only 3 from sig2
  //   as we insist the annotations must be multiples of the EDF record size, the final record will have 25, or 2.
  //   i.e. these all adds up to make complete records;  but if the start offset isn't completely landing on an
  //    sample point, then there will be some (minor) implied shift in the timing,  i.e. sig2 sample 3 will start at (e.g.) 0.25 seconds rather than 0.6 seconds;

  // i.e. this is a consequence of us tracking only the offset of each EDF record, which may contain signals at different rates

  //  0 0.0   --> (start at 0.25 seconds)    --->  first record is now    2       0.25  (EDF record offset from EDF start of 0.0)
  //  1 0.2                                                               3       0.45  (implied SR )
  //  2 0.4                                                               4       0.65  
  //  3 0.6                                                               0(b)    0.85  
  //  4 0.8                                                               1(b)    1.05  
  
  //  0(b) 1.0                               --->                         2(b)    1.25   Record POS = 1.25
  //  1(b) 1.2                                                            3(b)    1.45
  //  2(b) 1.4                                                            4(b)    1.65
  //  3(b) 1.6                                                            0(c)    1.85
  //  4(b) 1.8                                                            1(c)    2.05 



  //
  // Drop any other EDF annotations also, as it doesn't make sense
  // to try to keep those in the re-aligned EDF
  //

  drop_annots();
  
  //
  // Create a buffer for the new data
  //

  // write all new records here
  std::map<int,edf_record_t> new_records;
  
  // buffer for new records
  edf_record_t new_record( this );

  // track time-point of each new record (versus EDF header start)
  std::vector<uint64_t> tps;
  
  // allocate space
  for (int r=0;r<new_nr;r++) 
    new_records.insert( std::map<int,edf_record_t>::value_type( r , new_record ) );
  
  //
  // For each sample separaetly, copy Iterate over each annotation; add to new
  //
  
  bool got_tps = false;

  for (int s=0; s<header.ns; s++)
    {
      
      // skip the time-track
      if ( s == header.time_track() ) continue;

      // process this signal
      int fs = header.n_samples[s];
      int curr_rec = 0; // this should count 0 .. (new_nr-1)
      int curr_smp = 0; // this should count 0 .. (fs-1) 
      
      std::set<interval_t>::const_iterator aa = aset.begin();
      
      while ( aa != aset.end() )
	{
	  
	  const int downsample = 1; // i.e. no downsampling
	  const bool return_ddata = true ;

	  slice_t slice( *this , s , *aa , downsample , return_ddata );
	  
	  // point to digital values rather than physical; i.e. keep bitvalue/offset 
	  // as is in the EDF header for each signal, just change the number of records
	  // so no need to do digital->physical->digital conversion
	  
	  std::vector<int16_t> * d = slice.nonconst_ddata();
	  
	  const int n = d->size();
	  
	  // add to records
	  
	  // fetch record:
	  std::map<int,edf_record_t>::iterator rr = new_records.find( curr_rec );

	  if ( rr == new_records.end() ) Helper::halt( "internal error" );

	  edf_record_t & new_record = rr->second;
	  
	  for ( int i=0; i<n; i++)
	    {
	      //	      std::cout << " rec_ " << curr_rec << " " << curr_smp << " " << (*d)[i] << "\n";
	      // add data point
	      new_record.data[ s ][ curr_smp ] = (*d)[i]; 
	      
	      // add time-point for EDF record?
	      if ( curr_smp == 0 && ! got_tps ) 
		{
		  uint64_t rec_offset = aa->start + (uint64_t)curr_rec * header.record_duration_tp ;
		  
		  //		  std::cerr << " adding tps " << i << " " << rec_offset << "\t" << globals::tp_duration * rec_offset << "\n";
		  
		  tps.push_back( rec_offset );
		}
	      
	      // move on
	      ++curr_smp;

	      // pointing to a new record?
	      if ( curr_smp == fs )
		{
		  curr_smp = 0;
		  ++curr_rec;
		  // fetch record:
		  if ( curr_rec < new_nr )
		    {
		      rr = new_records.find( curr_rec );
		      if ( rr == new_records.end() ) Helper::halt( "internal error" );
		      new_record = rr->second;
		    }
		}	      
	    }
	  
	  // next annotation instance
	  ++aa;
	}

      // check that we correctly loaded all the annotations as expected
      
      //      std::cout << "curr_rec curr_smp = " << curr_rec << " " << curr_smp << " " << new_nr << "\n";

      if ( ! ( curr_rec == new_nr && curr_smp == 0 ) )
	Helper::halt( "problem loading up newly ALIGN'ed records" );
      
      got_tps = true;
      
      // next signal
    }
  
  if ( tps.size() != new_nr ) 
    Helper::halt( "internal error in ALIGN: time-point and record counts do not align" );


  //
  // Update EDF headers; only changing number of records here
  // i.e. record size, etc, remaining constant
  //
  
  header.nr = new_nr;
  

  //
  // Copy over new records
  //

  records = new_records;
  new_records.clear();


  //
  // Add time-track back?  We'll do this as a special case, rather
  // than use existing functions (i.e. as the core rec->timeline
  // mapping will change, which is an unusal sitation
  //
  
  timeline.tp2rec.clear();
  timeline.rec2tp.clear();
  timeline.rec2tp_end.clear();
  timeline.clear_epoch_mapping();

  // 
  // Iterate over all new records
  //
  
  //  std::cerr << "new_nr = " << new_nr << " " << tps.size() << "\n";

  for ( int r = 0;r < new_nr; r++)
    {
      // get the previously stored time-point 
      uint64_t tp = tps[r];
      timeline.tp2rec[tp] = r;
      timeline.rec2tp[r] = tp;
      timeline.rec2tp_end[r] = timeline.last_time_point_tp = tp + header.record_duration_tp - 1LLU;
    }

  timeline.total_duration_tp = (uint64_t)header.nr * header.record_duration_tp;

  //
  // Update time track  
  //
  
  if ( ! header.edfplus )
    {
      logger << "  restructuring as an EDF+ : ";
      set_edfplus();
    }
  
  set_discontinuous();
  

  int tt = header.time_track();  
  if ( tt == -1 ) Helper::halt( "internal error: could not find time-track" );

  std::cerr << "  updating header time-track: " << tt << "\n";

  // need to set a record size -- this should be enough?
  // by default 15-> up to 30 characters
  const int n_samples = globals::edf_timetrack_size;
  const int chars = 9;

  // for each new record
  int r = timeline.first_record();
  
  while ( r != -1 ) 
    {
      
      if ( r >= tps.size() ) 
	Helper::halt( "internal error in ALIGN: counting tps/new records" );

      double onset = globals::tp_duration * tps[ r ]; 
      
      //std::string ts = "+" + Helper::dbl2str_fixed( onset , chars) + "\x14\x14\x00";

      std::string ts = "+" + Helper::dbl2str( onset , chars) + "\x14\x14\x00";

      if ( r < 100 ) 
	std::cerr << "tps [" << tps[r] << "]\n";
      
	// std::cerr << "adding rec " << r << "\t" << "+" << Helper::dbl2str( onset , chars ) << "\n";
      //x      std::cout << "header = " << header.t_track << "  .. " << onset << "\n";

      records.find(r)->second.add_annot( ts , header.t_track );
      
      r = timeline.next_record(r);
    }

  
  //
  // All done
  //

  return true;
}
