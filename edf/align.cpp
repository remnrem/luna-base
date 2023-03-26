
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


void edf_t::set_timestamps( param_t & param )
{

  if ( header.nr == 0 ) return;
  
  const std::string filename = Helper::expand( param.requires( "file" ) );

  if ( ! Helper::fileExists( filename ) )
    Helper::halt( "could not find " + filename );

  std::vector<uint64_t> tps;
  std::ifstream IN( filename.c_str(), std::ios::in );
  while ( ! IN.eof() )
    {
      std::string x;
      Helper::safe_getline( IN , x );
      if ( IN.eof() ) break;
      if ( x == "" ) continue;

      double secs;
      if ( ! Helper::str2dbl( x , &secs ) )
	Helper::halt( "bad numeric value: " + x );

      tps.push_back( (uint64_t)( secs * globals::tp_1sec ) );
          
    }
  IN.close();

  logger << "  read " << tps.size() << " timestamps\n";

  // check this lines up
  if ( header.nr != tps.size() )
    Helper::halt( "expecting " + Helper::int2str( header.nr ) + " timestamps (i.e. to match number of EDF records" );

  // check all ascending
  for (int i=1; i<header.nr; i++)    
    if ( tps[i] <= tps[i-1] )
      Helper::halt( "found non-increasing consecutive time-points" );

  // make EDF+ (adds a time-track and at)  
  set_edfplus(); 

  // set as EDF+D explicitly
  set_discontinuous();
  
  // now update in-memory time-track
  timeline.create_discontinuous_timeline( tps );
  
  // now add EDF annotations w/ explcitly calculated tps
  add_time_track( &tps );

  logger << "  updated EDF+D time-track\n";
  
}


bool edf_t::edf_minus()
{

  // takes an EDF+D and makes a new standard EDF with gaps added (zeros)
  // it will also create an annot_offset_table[] for use w/ a subsequen WRITE-ANNOTS
  
  if ( ! header.edfplus )
    {
      logger << "  not already a standard EDF -- nothing for EDF-MINUS to do\n";
      return false;
    }

  // this function can work with original EDF+D but also in-memory 'EDF+D'
  //  i.e. after restructuring

  // special case where EDF (in-memory) is already continuous -- in this case, nothing to do,
  // so just force EDF in the standard manner
  
  if ( header.continuous )
    {
      logger << "  no discontinuities found -- peforming simple 'EDF' operation instead to force EDF\n";
      set_edf();
      //reset_start_time();
      return false;
    }

  //
  // otherwise, we have some form of a discontinuous EDF 
  //

  // get number of signals
  int nsigs = 0;
  for (int s=0; s<header.ns; s++)
    if ( header.is_data_channel(s) ) ++nsigs;
  
  logger << "  making a standard EDF with " << nsigs << " data channels\n";


  // record size          = header.record_duration   (sec)
  // implied new EDF size =  

  // actual signal in seconds
  double observed_sec = header.nr * header.record_duration ;

  // implied standaard EDF durtion
  double implied_sec = timeline.last_time_point_tp * globals::tp_duration ;

  double spliced_sec = implied_sec - observed_sec;


  //
  // keeping EDF record size constant, 
  // i.e. as per SEGMENTS command
  //
  
  std::set<interval_t> segs, gaps;

  int r = timeline.first_record();  
  uint64_t tp0 = timeline.rec2tp[r];
  uint64_t tp_start = tp0;
  
  uint64_t gap_start = 0; // i.e. always start at EDF starttime
  
  while ( r != -1 )
    {      
      // next record
      r = timeline.next_record( r );

      // start of this next record
      uint64_t tp;
      bool segend = false;
      
      // end?
      if ( r == -1 )
	{
	  // make this the 'previous'
	  tp0 = tp;
	  segend = true;
	}
      else
	{
	  tp = timeline.rec2tp[r] ;	  
	  // discontinuity / end of segment?
	  segend = tp - tp0 != header.record_duration_tp ;
	}
      
      // record this segment 
      if ( segend )
	{	  

	  uint64_t tp_stop = tp0 + header.record_duration_tp;

	  segs.insert( interval_t( tp_start , tp_stop ) );
	  
	  // did we observe a gap prior to this?
	  if ( tp_start > gap_start )
	    {
	      gaps.insert( interval_t( gap_start , tp_start ) );	      
	    }
	  
	  // reset start of next gap to end to this segment	      
	  gap_start = tp0 + header.record_duration_tp;
	  
	  // current point becomes start of the next segment
	  tp_start = tp;
	  
	}
      
      // current point becomes the last one, for next lookup
      tp0 = tp;

    }

  int num_segments = segs.size();
  int num_gaps = gaps.size();

  // the total segment duration should be an exact multiple of EDF record size
  // the gap durations need not be, however;
  // in the new standard EDF, each gap may need to be expanded to be an even
  // EDF record size

  // if we extend gaps, then we need to track the time added and appropriately
  // build an annot offset table for subsequent WRITE-ANNOTS

  //  std::map<
  
  

  return true;
}

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

      if ( annot != NULL ) 
	{
	  annot_map_t::const_iterator ii = annot->interval_events.begin();
	  while ( ii != annot->interval_events.end() )
	    {
	      //std::cout << "ANNOT = " << ii->first.interval << "\t" << ii->first.interval.duration() << "\n";
	      
	      // only consider this annotation is it is completely within a valid time-point
	      // range (i.e. no discontinuities, not out-of-range)
	      uint64_t adur =  ii->first.interval.duration();
	      
	      if ( timeline.valid_tps( ii->first.interval ) != adur ) 
		{
		  std::cerr << "  skipping " << ii->first.interval << "\n";
		  ++skipped_disc;
		  ++ii;
		  continue;
		}
	      
	      // check duration is a multiple of the EDF record duration
	      if ( adur % header.record_duration_tp ) 
		{
		  std::cerr << "  skipping, not an exact multiple of EDF recdur: " << ii->first.interval << "\n";
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
	{
	  logger << "  *** annotation overlapped prior:\n"
		 << last << "\n"
		 << *ii << "\n";

	  Helper::halt( "cannot specify overlapping annotations to EDF re-ALIGN-ment" );
	}
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

	  // as we enforce that ANNOTs should be multiples of EDF record size, when starting a new
	  // annot block, we should always also be starting a new record;   use in_annot_rec to track the
	  // offset (relative to the start of this annotation) for each added EDF record
	  
	  int in_annot_rec = 0;

	  if ( curr_smp != 0 )
	    Helper::halt( "internal logic error in ALIGN: new annotation not aligning with new EDF record;"
			  " potential floating-point inconistencies in interval specifications?\n" );
	  
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

	  if ( rr == new_records.end() ) Helper::halt( "internal error1" );

	  for ( int i=0; i<n; i++)
	    {

	      // add data point
	      rr->second.data[ s ][ curr_smp ] = (*d)[i]; 
	      
	      // add time-point for EDF record?
	      if ( curr_smp == 0 && ! got_tps ) 
		{
		  uint64_t rec_offset = aa->start + (uint64_t)in_annot_rec * header.record_duration_tp ;
		  tps.push_back( rec_offset );
		}
	      
	      // move on
	      ++curr_smp;

	      // pointing to a new record?
	      if ( curr_smp == fs )
		{
		  curr_smp = 0;
		  ++curr_rec;
		  ++in_annot_rec;
		  // fetch record:
		  if ( curr_rec < new_nr )
		    {
		      rr = new_records.find( curr_rec );
		      if ( rr == new_records.end() ) Helper::halt( "internal error2" );
		      new_record = rr->second;
		    }
		}	      
	    }
	  
	  // next annotation instance
	  ++aa;
	}

      // check that we correctly loaded all the annotations as expected
      
      if ( ! ( curr_rec == new_nr && curr_smp == 0 ) )
	Helper::halt( "problem loading up newly ALIGN'ed records" );
      
      got_tps = true;
      
      // next signal
    }

  
  // no signal channels? --> all done here
  if ( ! got_tps ) return false;
  
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
  // Note that EDF start times can only be in exact seconds::: so, if the new implied EDF start time is not a fraction of a second,
  // remove this extra part;   but then also remove for a) all other records (tps[] prior to writing) , and b) all annotations
  //

  uint64_t edf_start_tp = tps[0];
  //  uint64_t edf_start_fraction = edf_start_tp % globals::tp_1sec;

  logger << "  to obtain an EDF starting on an exact second, adjusting this and any annotations by -"
	 << edf_start_tp * globals::tp_duration << " seconds\n";
  
  // adjust all records to start at zero
  for ( int r = 0;r < new_nr; r++)
    tps[r] -= edf_start_tp;
  
  // adjust all loaded annotations for subseuent WRITE-ANNOTS
  // that will be called; this ensures they will start at elapsed time = 0 
  
  timeline.annotations.set_annot_offset( edf_start_tp );

  // change EDF header starttime also

  clocktime_t et( header.starttime );
  if ( et.valid )
    {
      double time_sec = ( edf_start_tp * globals::tp_duration ) ;
      et.advance_seconds( time_sec );
      header.starttime = et.as_string();
      logger << "  resetting EDF header starttime to " << header.starttime << "\n";
    }
  else
    {
      logger << "  no valid EDF header startime: setting to null 00.00.00\n";
      header.starttime = "00.00.00";
    }

  // and fix any special variables in the annotation_set_t
  
  timeline.annotations.duration_sec = new_nr * header.record_duration ;

  // no fractional seconds                                                     
  timeline.annotations.duration_hms = Helper::timestring( globals::tp_1sec * timeline.annotations.duration_sec , '.' , false ); 

  timeline.annotations.start_hms = header.starttime ;

  timeline.annotations.epoch_sec = timeline.epoched() ? timeline.epoch_length() : globals::default_epoch_len ;
  
  
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
      logger << "  restructuring as an EDF+\n";
      set_edfplus();
    }
  
  set_discontinuous();
  
  int tt = header.time_track();  
  if ( tt == -1 ) Helper::halt( "internal error: could not find time-track" );

  // global variable stores size of time-track to use
  //  by default = 15 (up to 30 chars, so plenty)
  const int n_samples = globals::edf_timetrack_size;
  const int chars = 5; // number of DP
  
  // for each new record
  int r = timeline.first_record();
  
  while ( r != -1 ) 
    {
      
      if ( r >= tps.size() ) 
	Helper::halt( "internal error in ALIGN: counting tps/new records" );

      double onset = globals::tp_duration * tps[ r ]; 
      
      std::string ts = "+" + Helper::dbl2str( onset , chars) + "\x14\x14\x00";
      
      records.find(r)->second.add_annot( ts , header.t_track );
      
      r = timeline.next_record(r);
    }

  
  //
  // All done
  //

  return true;
}
