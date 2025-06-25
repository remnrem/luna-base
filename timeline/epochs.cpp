
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


#include "timeline/timeline.h"

#include "edf/edf.h"
#include "edf/slice.h"
#include "db/db.h"
#include "helper/logger.h"
#include "helper/helper.h"
#include "annot/annotate.h"

extern writer_t writer;
extern logger_t logger;

bool timeline_t::generic_epochs() const
{
  return ! standard_epochs ; 
}

bool timeline_t::fixed_epoch_length() const
{
  return fixed_size_epochs;
}

bool timeline_t::check( const std::string & cmd ) const
{

  // fine if no epochs are set - this is dealt w/ by the command
  if ( ! epoched() ) return true;
  
  // track list of commands that cannot be used w/ nonstandard epochs
  //   not present in xlist  -- no constraints
  //       present    xlist == 1   requires standard and non-overlapping epochs
  //       present    xlist == 2   requires fixed duration epochs (i.e. no can overlapping or from epochs)

  // may want to update this w/ some commands that assume uniform AND contiguous epochs
  // e.g. IRASE, MOVING-AVERAGE, etc
  std::map<std::string,int> xlist;
    
  xlist[ "HYPNO" ] = 1; 
  xlist[ "ARTIFACTS" ] = 1; 
  xlist[ "EVAL-STAGES" ] = 1; 
  xlist[ "PLACE" ] = 1; 
  xlist[ "POPS" ] = 1; 
  xlist[ "REBASE" ] = 1; 
  xlist[ "SOAP" ] = 1; 
  xlist[ "STAGE" ] = 1; 
  
  xlist[ "CC" ] = 1; 
  xlist[ "COH" ] = 2; 
  xlist[ "IRASA" ] = 1; 
  xlist[ "LINE-DENOISE" ] = 1; 
  xlist[ "MOVING-AVERAGE" ] = 1; 
  xlist[ "PEAKS" ] = 1; 
  xlist[ "ROBUST-NORM" ] = 1; 
  xlist[ "SUPPRESS-ECG" ] = 1; 
  xlist[ "ZC" ] = 1; 

  std::map<std::string,int>::const_iterator cc = xlist.find( cmd );
  if ( cc == xlist.end() ) return true;

  const int etype = cc->second;

  // requires only fixed epoch durations
  if ( etype == 2 ) return fixed_epoch_length();
  
  // more stringent, not allowing generic epochs of any kind
  if ( generic_epochs() ) return false;
  
  return true;
}


int timeline_t::calc_epochs()
{

  // wipe any generic epoch params
  standard_epochs = true;
  fixed_size_epochs = true;
  epoch_generic_param_annots.clear();
  epoch_generic_param_w1 = epoch_generic_param_w2 = 0;
  epoch_generic_param_set_point = 0;
  epoch_generic_param_min_epoch_size = 0.1;
  
  // Calculates EPOCH timings in the original EDF timescale
  // for both continuous and discontinuous EDFs

  // Also populates rec2epoch and epoch2rec mappings

  
  //
  // DISABLED: In all cases, epoch length must be >= record length
  //   --> the above constrain is no longer necessary (and does not make sense in context of
  //       generic epochs in any case)
  // if ( epoch_length_tp < edf->header.record_duration_tp )
  //   Helper::halt( "epoch duration must be greater or equal to EDF record size\n         which is currently "
  // 		  + Helper::dbl2str( edf->header.record_duration_tp * globals::tp_duration ) + " second(s); "
  // 		  + "see RECORD-SIZE command to change this" );
  
  if ( epoch_length_tp < epoch_inc_tp ) 
    Helper::halt( "epoch increment cannot be larger than epoch duration" );
  
  epochs.clear();

  epoch_labels.clear();
  
  mask.clear();
  
  rec2epoch.clear();
  
  epoch2rec.clear();

  //
  // Easy case for continuous EDFs ( == cumul_ungapped)
  //
  
  if ( edf->header.continuous )
    {
      
      // set first timepoint (if any offset)

      uint64_t s = epoch_offset_tp;
      
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

	  // label is E1, E2, ... 1-based encoding
	  epoch_labels.push_back( "E" + Helper::int2str( (int)epochs.size() ) );
	  
	  // find matching records (within epoch interval)
	  interval_t search_interval( s , end );
	  
	  int start_record = search_interval.start / edf->header.record_duration_tp;
	  int stop_record = search_interval.stop / edf->header.record_duration_tp;
	  
	  int e = epochs.size()-1;
	  for (int r=start_record; r<=stop_record; r++) 
	    {
	      epoch2rec[ e ].insert( r );
	      rec2epoch[ r ].insert( e );
	    }
	  
	  // shift to next interval
	  // i.e. keeping the same offsets (if any) across the entire recording
	  //      in the continous (EDF/EDF+C) case
	  s += epoch_inc_tp;
	}
    }
  else if ( gap_spanning_epochs )
    {
      //
      // Epochs for the discontinuous case: but ignoring gaps
      //  i.e. 30sec means 30 seconds of cumulative signal only
      //       as if the study was being read as a standard EDF
      //       but this preserves time-stamps / we don't have to
      //       change the signals
      //

      //  Segments    111111111111111      2222222222222    3333333
      //  Epochs
      //   Gapped     111112222233333      4444455555       66666   
      //   Ungapped   111112222233333      4444455555666    6677777


      // initial record
      int r = first_record();

      // track epoch number
      int e = 0;
      
      // initial epoch
      uint64_t estart = rec2tp[r];

      // cumulative epoch (whole length; also inc step)
      uint64_t ecumm = 0LLU;

     
      // increment / restarts
      uint64_t estart2 = 0LLU; 
      uint64_t ecumm2 = 0LLU;
      int rstart2 = 0;
      
      // under this mode, epoch length/inc must be >= record length
      if ( epoch_length_tp < edf->header.record_duration_tp )
	Helper::halt( "epoch length cannot be less than record length" );
      if ( epoch_inc_tp < edf->header.record_duration_tp )
	Helper::halt( "epoch inc cannot be less than record length" );
      
      while ( 1 )
	{

	  // std::cout << " +++++ epoch " << e << " consider rec " << r << " at " << rec2tp[r] << "\n";
	  // std::cout << "   cumul = " << ecumm 
	  // 	    << " " << ecumm2 << "\n";
	  
	  // tag this record for this epoch if it spans it
	  if ( ecumm + edf->header.record_duration_tp <= epoch_length_tp )
	    {
	      epoch2rec[ e ].insert( r );
	      rec2epoch[ r ].insert( e );
	    }
	  
	  // is the increment (start of next epoch) in this record?
	  if ( estart2 == 0LLU && ecumm2 + edf->header.record_duration_tp > epoch_inc_tp )
	    {
	      uint64_t diff = epoch_inc_tp - ecumm2 ;
	      estart2 = rec2tp[r] + diff;
	      rstart2 = r;
	    }
	  else
	    {
	      // move this along
	      ecumm2 += edf->header.record_duration_tp;
	    }

	  
	  //
	  // can we add this record and continue building the epoch? 
	  //

	  if ( ecumm + edf->header.record_duration_tp <= epoch_length_tp )
	    {
	      // accumulate and skip to the next record
	      ecumm += edf->header.record_duration_tp;
	      r = next_record(r);
	      if ( r == -1 ) break;
	      continue;
	    }


	  // otherwise, implies the epoch ends in this record
	  
	  // start if this record
	  const uint64_t tp_start = rec2tp[r];
	  
	  // amount to add
	  const uint64_t diff = epoch_length_tp - ecumm ;

	  // ?? don't need +1 end-point here
	  interval_t saved_interval( estart , tp_start + diff );
	  epochs.push_back( saved_interval );
	  //	  std::cout << " adding " << e << " " << estart << " " << tp_start + diff << "\n";
	  
	  // label is E1, E2, ... 1-based encoding
	  epoch_labels.push_back( "E" + Helper::int2str( (int)epochs.size() ) );

	  // track epoch count
	  ++e;

	  //
	  // restart start
	  //

	  // tp of this next epoch
	  estart = estart2;

	  // reset this
	  estart2 = 0LLU;
	  
	  // record that this next epoch is starting in
	  r = rstart2;
	  
	  // offset into this record
	  uint64_t offset = estart - rec2tp[r];
	  
	  // add rest of this record (both inc. and 
	  ecumm = ecumm2 = edf->header.record_duration_tp - offset ; 

	  if ( ecumm != 0LLU ) {
              epoch2rec[ e ].insert( r );
              rec2epoch[ r ].insert( e );	    
	      // std::cout << " restarting:: " << estart << " " << r << "\n";
	      // std::cout << " restarting w/ " << ecumm << " on the clock\n";
	  }
	  
	  // advance to next record
	  r = next_record(r);
	  
	  // all done...
	  if ( r == -1 ) break;
	  
	}

      
    }  
  else
    {

      //
      // Epochs for the discontinuous case:
      //

      // annot_alignment: optionally, align all start epochs
      // (i.e. each first epoch within each contiguous region) to one
      // of a set of annotations
      // i.e. this is the more complex version of the single
      // 'epoch_offset_tp' trick, which only works for the continuous
      // timeline case.  This makes use the first_in_interval() function
      // which returns the time-point offset given an interval and a set of
      // annot labels
      
      const bool annot_alignment = epoch_align_annots.size() != 0 ; 

      // overview of the algorithm below:
      //   rec2tp[] start of each record
      //    - by default, epochs are aligned with records
      //      alternatively, they can be aligned with annotation sets in epoch_align_annots[] (i.e. N1, N2, ...) 
      
      //   epoch_length_tp length of epochs (standard case) -- only full epochs added
      //   epoch_inc_tp                                     -- increment
      //     (if < epoch_length_tp, overlapping epochs)
      //   erestart - where the /next/ epoch will start (might be over-lapping w/ the current epoch)
      
      // epoch-to-record mappings:
      //  putative_e2r 
      //  putative_r2e

      // r = current record
      // e = current epoch
      
      int r = first_record();
      
      if ( r == -1 ) return 0;
      
      // epochs have to be continuous in clocktime
      // putative start (i.e. start of record)
      uint64_t estart = rec2tp[r];

      // but, if we are allowing annot offset alignment, we may need to skip ahead a bit? (i.e. will go the 
      // the start of the next valid annot
      // moves start (and record) forward to the next in the annot set
      //  align_annots( &estart , &record, start points... set<uint64_t> ) 
      
      std::set<uint64_t> astarts;
      
      if ( annot_alignment )
	{
	  
	  astarts = annotations->starts( epoch_align_annots , epoch_length_tp );
	  
	  // this function allows for case where annotations are longer than the epoch size... cut into epoch duration
	  // starts
	  //     e.g.
	  //      0 -- 30 
	  //     60 -- 180
	  
	  // if epoch size = 30, implies possible starts of: 0 , 60 , 90 , 120 and 150 

  	  logger << "  within each segment, aligning epochs to " << astarts.size()
		 << " possible starting points from (" << epoch_align_str << ")\n";

	  
	  // next align_epochs() will advance the estart (tp) and spanning record (r) as needed
	  
	  // if we don't find any instances of an alignment annot in this segment
	  // then fall back on the default (i.e. segment start);  as the f()
	  // modifies estart, we should save here, and check retcode below

	  const uint64_t estart0 = estart;
	  const int r0 = r;
	  const bool retcode = align_epochs( &estart , &r , astarts );
	  
	  // std::cout << " DONE align_epochs() :::  retcode = " << retcode << "\t";
	  // std::cout << " estart = " << estart << "\n";
	  
	  // could not find any of the specified annotations to start alignment
	  // within the first segment; fall back to start of segment

	  if ( ! retcode )
	    {
	      estart = estart0;
	      r = r0;
	    }
	  
	}
      
      // we've now found a start point:: see where it would end (i.e. will it fit in contiguous 
      // time, and figure out which records are spanned)
      
      // for purpose of searching, skip last point
      // i.e. normally intervals are defined as END if 1 past the last point
      uint64_t estop  = estart + epoch_length_tp - 1LLU;
      
      // also, track so we can figure out which record the next epoch starts in 
      // for non-overlapping epochs, this will equal estop, otherwise it will come earlier
      // we assume that /within/ a contiguous segment, alignment-annotations are uniform
      // i.e. only try to re-align at the start of a new segment

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

	  // std::cout << "dets " << rec_start << " " << rec_end << "\t"
 	  //  	    << estart << " " << erestart << " " << estop << "\n";

	  //
	  // Will the next epoch potentially come from this record?
	  //
	  //std::cout << " erestart " << erestart << "\n";
	  if ( erestart >= rec_start && erestart <= rec_end )
	    {
	      //std::cout << "  track this is the restarting record\n";
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

	      // label is E1, E2, ... 1-based encoding
	      epoch_labels.push_back( "E" + Helper::int2str( (int)epochs.size() ) );
	      
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

	      //
	      // Any annot-alignment?
	      //
	      
	      if ( annot_alignment )
		{
		  // skip ahead?
		  const uint64_t estart0 = estart;
		  const int r0 = r;
		  const bool retcode = align_epochs( &estart , &r , astarts );
		  if ( ! retcode )
		    {
		      // fall-back 
		      estart = estart0;
		      r = r0;
		      //break;
		    }
		}
	      

	      //
	      // Complete the rest of the epoch
	      //

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

	      putative_r2e[r].insert(e);
	      putative_e2r[e].insert(r);

	      r = next_record(r);
	      
	      // if out of epochs, then just quit (i.e. we won't be adding this
	      // putative epoch in any case)
	      
	      if ( r == -1 ) break;
	      
	      // check that we did not skip anything
	      // these two values should be EQUAL if
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
		  
		  //
		  // Any annot-alignment?
		  //
		  
		  if ( annot_alignment )
		    {
		      // skip ahead?
		      const uint64_t estart0 = estart;
		      const int r0 = r;
		      const bool retcode = align_epochs( &estart , &r , astarts );
		      if ( ! retcode )
			{
			  estart = estart0;
			  r = r0;
			  //break;
			}
		    }
		  
		  // complete the rest of the epoch definition
                  estop  = estart + epoch_length_tp - 1LLU;
		  erestart = estart + epoch_inc_tp;
		  putative_e2r.clear();
		  putative_r2e.clear();
		}
	      
	    }
	  
	  // next record 

	}
    }

  // reset counter
  current_epoch = -1;
  mask.clear();
  mask.resize( epochs.size() , false );
  mask_set = false;
  mask_mode = 0; 

  if ( 0 )
    debug_dump_epochs();
  
  // all done
  return epochs.size();    
}


int timeline_t::ensure_epoched() 
{
  
  if ( epoched() ) return num_epochs();
  
  // otherwise, set to defaults
  
  int ne = set_epoch( globals::default_epoch_len , globals::default_epoch_len );
  
  logger << "  set epochs to default " 
	 << globals::default_epoch_len 
	 << " seconds, " << ne << " epochs\n";
  
  return ne;
}

bool timeline_t::epoched() const
{
  //  std::cout << " standard = " << standard_epochs << "\t" << epoch_length_tp << "\n";
  return (!standard_epochs) || epoch_length_tp != 0;
} 


void timeline_t::unepoch() 
{

  current_epoch = -1;
  
  // No EPOCHs
  epoch_length_tp = 0L;
  epoch_inc_tp = 0L;
  epoch_offset_tp = 0L;
  epoch_align_annots.clear();
  epoch_align_str = "";
  epochs.clear();
  epoch_labels.clear();

  // assume will be standard epochs next set
  standard_epochs = true;
  fixed_size_epochs = true;
  gap_spanning_epochs = false;

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
 

int timeline_t::reset_epochs()
{

  // called automatically after RESTRUCTURE only
  // this doesn't wipe any epoch-level annotations, 
  // as they are fixed to the original file-EDF epochs
  // (for this particular epoch size in any case)

  //  directly transfer epochs boundaries of all unmasked epochs;
  //  after a RE, all retained epochs should have a proper mapping

  std::vector<interval_t> new_epochs;
  std::vector<std::string> new_epoch_labels;
  
  if ( mask.size() != epochs.size() )
    Helper::halt( "internal error in timeline_t::reset_epochs() - mask size does not match epoch size" );
  
  for (int e=0; e<epochs.size(); e++)
    {
  
      // only copy unmasked epochs
      if ( ! mask[e] )
	{
	  interval_t e1 = epochs[e] ;
	  interval_t e2( 0LLU , 0LLU );

	  // nb, step back 1 from stop point, as will do a record-inclusive search
	  //     but then need to add it back on...

	  bool okay1 = remap_timepoint( e1.start , &(e2.start) );
	  
	  bool okay2 = remap_timepoint( e1.stop - 1LLU , &(e2.stop) );
	  ++e2.stop;


	  
	  // logger << "e" << e << "\tM"
	  //  	 << mask[e] << "\t start/end mapping = "
	  //  	 << okay1 << okay2
	  //  	 << "\t" << e1.as_string() << "\t" << e2.as_string() << "\t"
	  //  	 << e2.start << " " << e2.stop << "\n";
	  
	  if ( okay1 && okay2 )
	    {
	      // push old epoch back
	      new_epochs.push_back( e1 );
	      new_epoch_labels.push_back( epoch_labels[e] );
	    }
	  else
	    {
	      logger << "e" << e << "\tM"
		     << mask[e] << "\t start/end mapping = "
		     << okay1 << okay2 
		     << "\t" << e1.as_string() << "\t" << e2.as_string() << "\n";
	      Helper::halt( "internal error in timeline_t::reset_epochs()" );
	    }
	}      
    }

  // copy retained, remapped epochs back
  
  epochs = new_epochs;

  epoch_labels = new_epoch_labels;

  
  // populate record-epoch mappings

  rec2epoch.clear();
  epoch2rec.clear();

  for (int e=0; e<epochs.size(); e++)
    {
      const interval_t interval = epochs[e];
      
      // get the records in this epoch      
      std::set<int> records = records_in_interval( interval );

      //std::cout << " got " << records.size() << "  in epoch " << e << "\n";
      
      // get records
      std::set<int>::const_iterator rr = records.begin();
      while ( rr != records.end() )
	{
	  epoch2rec[ e ].insert( *rr );
	  rec2epoch[ *rr ].insert( e );	  
	  ++rr;
	}      
    }

  
  // reset counter & mask
  current_epoch = -1;
  mask.clear();
  mask.resize( epochs.size() , false );
  mask_set = false;
  mask_mode = 0; 
  
  first_epoch();
    
  return epochs.size();
    
}


int timeline_t::set_epoch(const double s, const double o , const uint64_t offset , 
			  const std::string align_str  , 
			  const std::vector<std::string> * align_annots ) 
{ 
  if ( s <= 0 || o < 0 )
      Helper::halt( "cannot specify negative epoch durations/increments");
  
  clear_epoch_annotations();
  epoch_length_tp = s * globals::tp_1sec;
  epoch_inc_tp = o * globals::tp_1sec;
  standard_epochs = true;
  fixed_size_epochs = true;
  
  // pass offset as uint64_t
  epoch_offset_tp = offset ; // * globals::tp_1sec;
  
  epoch_align_str = align_str;
  if ( align_annots != NULL ) 
    epoch_align_annots = *align_annots;
  
  if ( epoch_length_tp == 0 || epoch_inc_tp == 0 ) 
    Helper::halt( "invalid epoch parameters" );
  first_epoch();
  return calc_epochs();
}

double timeline_t::epoch_length() const 
{
  if ( standard_epochs ) 
    return (double)epoch_length_tp / globals::tp_1sec;
  if ( current_epoch != -1 && epochs.size() > current_epoch )
    return epochs[ current_epoch ].duration_sec();
  return 0;
}
  
double timeline_t::epoch_inc() const 
{
  return (double)epoch_inc_tp / globals::tp_1sec;
}

double timeline_t::epoch_offset() const 
{
  return (double)epoch_offset_tp / globals::tp_1sec;
}

bool timeline_t::epoch_any_offset() const 
{
  return epoch_offset_tp != 0L;
}

std::string timeline_t::align_string() const 
{
  return epoch_align_str;
}
  
bool timeline_t::exactly_contiguous_epochs() const
{
  return epoch_length_tp == epoch_inc_tp;
} 
  
double timeline_t::epoch_len_tp() const 
{
  if ( standard_epochs )
    return (double)epoch_length_tp;
  if ( current_epoch != -1 )
    return epochs[ current_epoch ].duration();
  return 0LLU;

}
  
uint64_t timeline_t::epoch_increment_tp() const 
{
  return epoch_inc_tp ;
} 

uint64_t timeline_t::epoch_len_tp_uint64_t() const
{
  return epoch_length_tp ;
}

  
int timeline_t::first_epoch()  
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


int timeline_t::next_epoch()  
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

int timeline_t::next_epoch_ignoring_mask()  
{ 
  ++current_epoch;
  if ( current_epoch == epochs.size() ) return -1;
  return current_epoch;
}

// all unmasked epochs
int timeline_t::num_epochs() const 
{
  if ( ! mask_set ) return epochs.size();
  int r = 0;
  for (int i=0;i<mask.size();i++)
    if ( ! mask[i] ) ++r;
  return r;
}

// all epochs
int timeline_t::num_total_epochs() const 
{ 
  return epochs.size(); 
}

interval_t timeline_t::epoch( const int e ) const
{
  if ( e < 0 || e >= epochs.size() ) return interval_t(0LLU,0LLU);
  return epochs[e]; 
}

bool timeline_t::epoch_records( const int e , int * a , int * b ) const 
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
// Epoch mappings
//

void timeline_t::clear_epoch_mapping()
{
  epoch_orig2curr.clear();
  epoch_curr2orig.clear();
}

bool timeline_t::has_epoch_mapping() const
{
  return epoch_orig2curr.size() != 0 ;
}

int timeline_t::original_epoch(int e)
{
  if ( ! has_epoch_mapping() ) return e;
  if ( epoch_curr2orig.find(e) == epoch_curr2orig.end() ) return -1;
  return epoch_curr2orig.find(e)->second;
}

// 1-based epoch mapping
int timeline_t::display_epoch(int e) const
{      
  if ( ! has_epoch_mapping() ) return e+1;
  if ( epoch_curr2orig.find(e) == epoch_curr2orig.end() ) return -1;
  return epoch_curr2orig.find(e)->second + 1 ;
}


int timeline_t::display2curr_epoch(int e) const 
{
  if ( ! has_epoch_mapping() ) return e-1;    
  if ( epoch_orig2curr.find(e-1) == epoch_orig2curr.end() ) return -1;
  return epoch_orig2curr.find(e-1)->second ;
}



//
// Helpers
//

std::map<int,bool> timeline_t::spanning_epoch_masks( const int r ) 
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

      //logger << "  adding epoch mapping table\n";
      
      clear_epoch_mapping();
      
      int curr = 0; 
      
      while ( 1 ) 
	{      
	  int epoch = next_epoch_ignoring_mask();      
	  
	  if ( epoch == -1 ) break;	  
	  
	  if ( ! masked_epoch( epoch ) )
	    {
	      //	      std::cout << " adding epoch-mapping: " << epoch << " --> " << curr << "\n";
	      epoch_orig2curr[ epoch ] = curr;
	      epoch_curr2orig[ curr ] = epoch;
	      ++curr;
	    }
	}
    }
  else // otherwise, already has a mapping
    {
      //logger << "  using existing epoch-mapping\n";
      
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
	      //std::cout << " remapping " << orig << " --> " << curr << "\n";
	      epoch_orig2curr[ orig ] = curr;
	      epoch_curr2orig[ curr ] = orig;
	      ++curr;
	    }
	  
	}
    }
  
}
  


int timeline_t::calc_epochs_generic_from_annots( param_t & param )
{
    
  // clear everything first
  unepoch();

  // first time this is called, from EPOCH command
  standard_epochs = false;
  gap_spanning_epochs = false;
  
  // option must specify the fixed duration of these (from annot start)
  fixed_size_epochs = param.has( "fixed" );
  if ( fixed_size_epochs ) 
    {
      const double f = param.requires_dbl( "fixed" );
      if ( f <= 0.001 ) Helper::halt( "fixed duration must be positive (secs)" );
      epoch_length_tp = f * globals::tp_1sec;
      epoch_inc_tp = 0LLU; // not defined/used
    }

  // if running with 'fixed', then 'only-one' means do not add multiple
  // eppchs from the same annot, i.e. if they would fit
  const bool add_all_fixed = ! param.has( "only-one" );
  if ( (!add_all_fixed) && ! fixed_size_epochs ) Helper::halt( "can only add 'only-one' with 'fixed' "); 

  //  ANNOT [-----------------------------]

  // default
  //  FIXED [12345][12345][12345][12345]--]
  // with 'only-one'
  //  FIXED [12345]-----------------------]
  
  
  // main annotations to use, allowing wildcards
  epoch_generic_param_annots = annotate_t::root_match( param.strset( "annot" ) , annotations->names() );
  
  if ( (!param.has( "annot" )) || epoch_generic_param_annots.size() == 0 )
    Helper::halt( "no 'annot' specified to define epochs" );

  // 'else' annotations - will define a new annot and corresponding epochs
  // based on being 'not' listed in 'annot' - but this respects any EDF+D
  // segments/gaps

  const bool else_epochs = param.has( "else" ) && ! param.empty( "else" );  
  const std::string else_epoch_label = else_epochs ? param.value( "else" ) : "";

  // can't already exist
  if ( else_epochs )
    {
      annot_t * annot = annotations->find( else_epoch_label );      
      if ( annot != NULL )
	Helper::halt( "'else' cannot specify an existing annotation label" );
    }

  // add this as a new annotation - i.e. so it can be used in masks, etc
  // downstream
  annot_t * else_annots = else_epochs ? annotations->add( else_epoch_label ) : NULL ; 
  
  // add flanking / set mid-points or start/end points?
  epoch_generic_param_set_point = 0;
  if      ( param.has( "midpoint" ) ) epoch_generic_param_set_point = 2;
  else if ( param.has( "start" ) ) epoch_generic_param_set_point = 1;
  else if ( param.has( "stop" ) ) epoch_generic_param_set_point = 3;

  // window/flanking around points
  const bool has_w = param.has( "w" );
  const bool has_w_before = param.has( "w-before" );
  const bool has_w_after = param.has( "w-after" );
  const bool some_w = has_w || has_w_before || has_w_after;
  if ( has_w && ( has_w_before || has_w_after ) )
    Helper::halt( "can only specify w or ( w-before and/or w-after )" );
    
  epoch_generic_param_w1 = epoch_generic_param_w2 = 0;
  if ( has_w ) 
    epoch_generic_param_w1 = epoch_generic_param_w2 = param.requires_dbl( "w" );
  else
    {
      if ( has_w_before )
	epoch_generic_param_w1 = param.requires_dbl( "w-before" );
      
      if ( has_w_after )
	epoch_generic_param_w2 = param.requires_dbl( "w-after" );
    }
  
  if ( epoch_generic_param_w1 < 0 ) Helper::halt( "'w' (or w-before/w-after) cannot be negative" );
  if ( epoch_generic_param_w2 < 0 ) Helper::halt( "'w' (or w-before/w-after) cannot be negative" );
  if ( epoch_generic_param_set_point && ( (!some_w) || fabs( epoch_generic_param_w1 ) < 0.001 ) ) 
    Helper::halt( "epochs too small: need larger 'w' (or w-before/w-after) if using 'midpoint/start/stop'" ); 
  if ( epoch_generic_param_set_point && ( (!some_w) || fabs( epoch_generic_param_w2 ) < 0.001 ) ) 
    Helper::halt( "epochs too small: need larger 'w' (or w-before/w-after) if using 'midpoint/start/stop'" ); 

  // window shift
  const bool has_shift = param.has( "shift" );
  epoch_generic_param_shift = 0;
  if ( has_shift ) 
    epoch_generic_param_shift = param.requires_dbl( "shift" );

  // truncate last N seconds 
  const bool has_trunc = param.has( "trunc" );
  epoch_generic_param_trunc = 0;
  if ( has_trunc )
    epoch_generic_param_trunc = param.requires_dbl( "trunc" );
  if ( epoch_generic_param_trunc < 0 )
    Helper::halt( "trunc must be positive" );


  // verbose debug output
  const bool debug = param.has( "debug" );
  
  // handling merging of contiguous/overlapping epochs
  const bool flatten = param.has( "flatten" ) ? param.yesno( "flatten" ) : true; 


  //
  // for now, if using fixed size epochs
  //
  
  if ( fixed_size_epochs ) 
    {
      if ( else_epochs )
	Helper::halt( "cannot use else with fixed" );

      // trunc and windows okay, as the fixed size extraction/check happens AFTER those
    }


  // require a minimum epoch size: set to 1/10th of a second by default
  epoch_generic_param_min_epoch_size = param.has( "min" ) ? param.requires_dbl( "min" ) : 0.1;
  if ( epoch_generic_param_min_epoch_size < 0.001 ) Helper::halt( "'min' must be 0.001 or greater" );

  // below, pull all annots and populate epochs[] and epoch_labels[]
  // also need to populate rec2epoch[] and epoch2rec[]


  // step 1, segment based on the background (i.e.
  //   EDF    -------------- [ GAP ] ------------
  //   ANNOT            ----------------
  //    --> becomes
  //                    ----         ---
  //

  std::set<interval_t> background = segments();

  if ( debug )
    logger << "found " << background.size() << " background segments\n";
  
  // compile all intervals from annots first
  std::set<interval_t> intervals0;
  
  // label for generic annots
  std::string egen_label = Helper::stringize( epoch_generic_param_annots );
  
  // pull each annot class
  std::set<std::string>::const_iterator aa = epoch_generic_param_annots.begin();
  while ( aa != epoch_generic_param_annots.end() )
    {
      
      annot_t * annot = annotations->find( *aa );
      
      if ( annot == NULL )
	{
	  logger << "  *** could not find annotation " << *aa << "\n";
	  ++aa;
	  continue;
	}
      
      annot_map_t::const_iterator ii = annot->interval_events.begin();
      while ( ii != annot->interval_events.end() )
        {
	  intervals0.insert( ii->first.interval ) ;
          ++ii;
        }

      ++aa;
    }

  if ( debug )
    logger << "  considering " << intervals0.size() << " initial intervals\n";

  // split by background?
  if ( flatten ) 
    intervals0 = annotate_t::apairs( intervals0 , background , "intersection" );
  else
    {
      // do not use that as we don't want to flatten; means we also don't have a background splicing
      // .. (nothing currently) ..
    }

  if ( debug )
    logger << "  given " << background.size()
	   << " background segments, split intervals to "
	   << intervals0.size() << " intervals\n";
  
  // make an ordered list
  std::map<interval_t,std::string> intervals;

  int aidx = 0;
  
  // pull each annot class
  std::set<interval_t>::const_iterator ii2 = intervals0.begin();
  while ( ii2 != intervals0.end() )
    {
      
      interval_t interval = *ii2;

      if ( debug )
	logger << "\n  considering interval " << aidx++
	       << "\t" << interval.as_string() << "\n";
	  
      // adjustments?
      if ( epoch_generic_param_set_point )
	{
	  // zero-duration midpoint marker
	  if ( epoch_generic_param_set_point == 1 )
	    {
	      interval.stop = interval.start;
	    }
	  else if ( epoch_generic_param_set_point == 3 )
	    {
	      interval.start = interval.stop;
	    }
	  else // 2 == midpoint
	    {
	      uint64_t m = interval.mid();
	      interval.start = interval.stop = m;
	    }
	}
      
      // exapand?
      if ( some_w && ( epoch_generic_param_w1 > 0 || epoch_generic_param_w2 > 0 ) )
	{

	  if ( has_w ) 
	    {
	      if ( epoch_generic_param_w1 == epoch_generic_param_w2 )
		interval.expand( epoch_generic_param_w1 * globals::tp_1sec );
	      else
		{
		  interval.expand_left( epoch_generic_param_w1 * globals::tp_1sec );
		  interval.expand_right( epoch_generic_param_w2 * globals::tp_1sec ); 
		}
	    }

	  if ( has_w_before )
	    interval.expand_left( epoch_generic_param_w1 * globals::tp_1sec ); 

	  if ( has_w_after )
	    interval.expand_right( epoch_generic_param_w2 * globals::tp_1sec );	      
	}

      // shift?
      if ( has_shift ) 
	{
	  if ( epoch_generic_param_shift < 0 ) interval.shift_left( -epoch_generic_param_shift * globals::tp_1sec );
	  else if ( epoch_generic_param_shift > 0 ) interval.shift_right( epoch_generic_param_shift * globals::tp_1sec );
	}
      
      // truncate last N seconds?
      if ( has_trunc )
	{
	  // zero-dur event if truncating too much 
	  if ( epoch_generic_param_trunc >= interval.duration_sec() )
	    interval.stop = interval.start;
	  else
	    interval.stop -= epoch_generic_param_trunc * globals::tp_1sec;
	}

      // handle fixed-size as special case, as we may add 1+ epochs per annot
      
      // general case:
      if ( ! fixed_size_epochs )
	{
	  // add as an epoch, if large enough
	  if ( interval.duration_sec() >= epoch_generic_param_min_epoch_size )
	    {
	      if ( debug )
		logger << "  ++ adding (non-fixed size) as " << interval.as_string()
		       << "\t" << interval.duration_sec() << " secs\n";
	      intervals[ interval ] = egen_label;
	    }
	  else
	    {
	      if ( debug )
		logger << "  -- rejected (non-fixed size), not above "
		       << epoch_generic_param_min_epoch_size
		       << "\t (is only " << interval.duration_sec() << " secs)\n";
	      
	    }
	}
      else // i.e. fixed_size_epochs == T  
	{
	  // here, adding 1+ fixed size epochs from this single annot
	  
	  interval_t original = interval;
	  
	  while ( 1 )
	    {
	      
	      // adjust
	      if ( interval.duration() >= epoch_length_tp )
		interval.stop = interval.start + epoch_length_tp; // enforce end
	      else
		break;
	      
	      // add as an epoch, if large enough
	      if ( interval.duration_sec() >= epoch_generic_param_min_epoch_size ) 
		intervals[ interval ] = egen_label;
	      
	      // all done, or going back to try to add more?
	      if ( ! add_all_fixed ) break;
	      
	      if ( debug )
		logger << "  ++ adding, fixed size epoch "
		       << interval_t( interval.stop , original.stop ).as_string() << "\n";
	      
	      // update remaining interval
	      interval = interval_t( interval.stop , original.stop );
	    }	      
	}
      
      // next instance
      ++ii2;
    }

  if ( debug )
    logger << "\n---------------------------\n"
	   << "*** found " << intervals.size() << " intervals that meet size criteria\n";
  
  //
  // Optionally, if adding 'else' annot, need to know
  // breakpoints/gaps; add the new annot as well as the epochs
  // based on those annotations
  //
  
  std::set<interval_t> segs;
  
  if ( else_epochs )
    {
      
      if ( edf->header.continuous )	
	{
	  // start/end
	  uint64_t duration_tp = globals::tp_1sec * (uint64_t)(edf->header.nr) * edf->header.record_duration ;
	  segs.insert( interval_t( 0LLU , duration_tp ) );
	}
      else
	{
	  
	  int r = first_record();	  
	  uint64_t tp0 = rec2tp[r];	  
	  uint64_t tp_start = tp0;
	  
	  while ( r != -1 )
	    {
	      r = next_record( r );
	      uint64_t tp; // start of this next record
	      bool segend = false;

	      if ( r == -1 ) // end of recording
		{
		  // make this the 'previous'
		  tp0 = tp;
		  segend = true;
		}
	      else // or a discontinuity?
		{
		  tp = rec2tp[r] ;		  
		  // discontinuity / end of segment?
		  segend = tp - tp0 != edf->header.record_duration_tp ;
		}
	      
	      // record this segment 
	      if ( segend )
		{
		  uint64_t tp_stop = tp0 + edf->header.record_duration_tp;
		  interval_t interval( tp_start , tp_stop );
		  segs.insert( interval );
		  // current point becomes start of the next segment
		  tp_start = tp;		  
		}	      
	      // current point becomes the last one, for next lookup
	      tp0 = tp;	      	      
	    } // next record	  
	}


      //
      // We've now populated segs[] and intervals[] with
      // the background, and also with the intervals of interest
      //  -- go through and make new 'else' annot which is the
      //     complement of the annotations we have 
      //

      // first, will be helpful to flatten the intervals
      
      std::set<interval_t> eanns;

      std::map<interval_t,std::string>::const_iterator ii = intervals.begin();
      while ( ii != intervals.end() )
	{
	  eanns.insert( ii->first );
	  ++ii;
	}

      // remove epoch-defining annotations from the background set:
      segs = annotate_t::excise( segs , eanns );

      // what is left are the 'else' annots      
      logger << "  adding " << segs.size() << " else annotations, with label " << else_epoch_label << "\n";
      std::set<interval_t>::const_iterator ss = segs.begin();
      while ( ss != segs.end() )
	{
	  // add as annotation
	  else_annots->add( else_epoch_label , *ss , "." );
	  
	  // add as a to-be-created epoch
	  intervals[ *ss ] = else_epoch_label;
	  
	  ++ss;
	}
	     
    }
  
  
  //
  // add epochs, now they are in interval sort order
  //  -- optionally, splice in 'else' epochs too
  //
  
  std::map<interval_t,std::string>::const_iterator ii = intervals.begin();
  while ( ii != intervals.end() )
    {
      
      const interval_t & interval = ii->first ; 
      
      // check that the epoch is valid - i.e. does not span gaps
      uint64_t vtp = valid_tps( interval );

      if ( debug )	
	logger << "\n  checking interval gap spanning " << vtp << " vs " << interval.duration() << "\n";
      
      if ( vtp != interval.duration() )
	{
	  logger << "  skipping interval that falls in a gap " << interval.as_string() << "\n";
	  ++ii;
	  continue;
	}

      //
      // add epoch
      //
      
      // 0-based epoch index
      const int e = epochs.size();
      
      epochs.push_back( interval );
      
      epoch_labels.push_back( ii->second );

      // get the records in this epoch      
      std::set<int> records = records_in_interval( interval );

      std::set<int>::const_iterator rr = records.begin();
      while ( rr != records.end() )
	{
	  epoch2rec[ e ].insert( *rr );
	  rec2epoch[ *rr ].insert( e );	  
	  ++rr;
	}
      
      ++ii;
    }

  // reset counter
  current_epoch = -1;
  mask.clear();
  mask.resize( epochs.size() , false );
  mask_set = false;
  mask_mode = 0; 

  if ( 0 ) 
    debug_dump_epochs();

  return epochs.size();
}



void timeline_t::output_epoch_info( const bool verbose , const bool show_masked )
{
  
  int n_masked = 0, n_unmasked = 0; 
  
  // track clock time
  
  clocktime_t starttime( edf->header.starttime );
  
  bool hms = starttime.valid;
  
  first_epoch();

  uint64_t total_epoched = 0LLU;
  std::set<interval_t> fepochs; // --> flattened epochs

  while ( 1 ) 
    {

      // get next epoch (either masked or unmasked)
      int epoch0 = show_masked ? next_epoch_ignoring_mask() : next_epoch();      
      
      if ( epoch0 == -1 ) break;
      
      interval_t interval = epoch( epoch0 );
      // std::cout << " E int " << interval.as_string() << "\n";
      
      // original encoding (i.e. to allows epochs to be connected after the fact)
      if ( verbose ) 
	writer.epoch( display_epoch( epoch0 ) );
      
      // if present, original encoding
      if ( verbose )
	writer.value( "E1" , epoch0+1 );

      // is this epoch masked?
      const bool is_masked = mask_set && mask[epoch0] ; 

      if ( verbose ) 
	if ( show_masked )
	  writer.value( "EMASK" , is_masked );
      
      if ( is_masked )
	++n_masked;
      else
	++n_unmasked;

      total_epoched += interval.duration();
      fepochs.insert( interval );
      
      if ( verbose )
	{
	  writer.value( "LABEL" , epoch_labels[ epoch0 ] );
	  writer.value( "INTERVAL" , interval.as_string() );      
	  writer.value( "START"    , interval.start_sec() );
	  writer.value( "MID"      , interval.mid_sec() );
	  writer.value( "STOP"     , interval.stop_sec() );
	  writer.value( "TP"       , interval.as_tp_string() );
	  writer.value( "DUR"      , interval.duration_sec() );

	  // original time-points
	  
	  if ( hms )
	    {
	      const double sec0 = interval.start * globals::tp_duration;
	      clocktime_t present = starttime;
	      present.advance_seconds( sec0 );
	      std::string clocktime = present.as_string( ':' , true ); // disp. fractional seconds
	      writer.value( "HMS" , clocktime );
	    }
	}

      // next epoch
    }		  

  if ( verbose ) 
    writer.unepoch();
    
  writer.value( "NE" , n_unmasked );

  if ( show_masked ) 
    writer.value( "NE_MASKED" , n_masked );
  
  if ( standard_epochs ) {
    writer.value( "DUR" , epoch_length() );
    writer.value( "INC" , epoch_inc() );
    writer.value( "OFFSET" , epoch_offset() );
  }
  else if ( fixed_size_epochs )
    {
      writer.value( "DUR" , epoch_length() );
    }
  writer.value( "GENERIC" , (int)(!standard_epochs) );
  writer.value( "FIXED_DUR" , (int)(epoch_length() ) );

  // total duration of recording
  //  a) spanned by an epoch  
  //  b) not spanned

  fepochs = annotate_t::flatten( fepochs );
  uint64_t total_fepoched = 0LLU;
  std::set<interval_t>::const_iterator ff = fepochs.begin();
  while ( ff != fepochs.end() )
    {
      total_fepoched += ff->duration();
      ++ff;
    }

  double p_spanned = (double)total_fepoched / (double)total_duration_tp;
  
  writer.value( "TOT_DUR" , (double)total_epoched / (double)globals::tp_1sec );
  writer.value( "TOT_SPANNED" , (double)total_fepoched / (double)globals::tp_1sec );
  writer.value( "TOT_UNSPANNED" ,
		(double)( total_duration_tp - total_fepoched ) / (double)globals::tp_1sec );
  writer.value( "TOT_REC" , (double)total_duration_tp / (double)globals::tp_1sec );
  writer.value( "TOT_PCT" , p_spanned );
  
 
}


void timeline_t::debug_dump_epochs()
{
  std::cout << "records2epochs:\n";      

  std::map<int,std::set<int> >::const_iterator ii = rec2epoch.begin();
  while ( ii != rec2epoch.end() )
    {
      std::cout << "r" << ii->first << " -> ";
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
  
  std::cout << "\nepochs2records:\n";
  
  ii = epoch2rec.begin();
  while ( ii != epoch2rec.end() )
    {
      std::cout << "e" << ii->first << " " 
		<< epoch_labels[ ii->first ] << " " 
		<< " -> ";
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


