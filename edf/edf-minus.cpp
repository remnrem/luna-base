
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

#include "helper/helper.h"
#include "db/db.h"
#include "helper/logger.h"

extern writer_t writer;

extern logger_t logger;

// helpers
bool edf_minus_helper_align( const std::set<instance_idx_t> & e ,
			     const interval_t & s ,
			     const std::set<std::string> & a ,
			     const uint64_t d, 
			     interval_t * t );

bool edf_minus_helper_has_annot( const std::set<instance_idx_t> & e ,
				 const interval_t & s ,
				 const bool whole , 
				 const uint64_t d,				 
				 const std::set<std::string> & a );

std::map<std::string,int> edf_minus_helper_count_annots( const std::set<instance_idx_t> & e ,
							 const interval_t & s ,
							 const bool whole,
							 const uint64_t d,							 
							 const std::set<std::string> & a ,
							 std::map<std::string,int> * tots );

uint64_t edf_minus_helper_whole_sec( uint64_t tp , uint64_t * f );

bool edf_t::edf_minus( param_t & param )
{
  
  // Given an EDF(+D) and optionally annotations:
  //   - create a new edf_t with
  //     - record size of 1 second
  //       ZOH for signals with SRs < 1 Hz
  //     - always standard EDF
  //       any gaps from an EDF+D are either:
  //          1) zero-padded
  //          2) segments concatenated (but w/ annotations made/changed)
  //          3) single (largest) segment selected
  //     - if alignment annotations, ensure we align (chop signals) to match

  //   - when times are generated, shifted so that always start at 0 seconds

  // this function can work with original EDF+D but also in-memory 'EDF+D'
  //  i.e. after restructuring


  // In practice, to handle edge cases: 

  // If any signals have SR < 1 Hz, needs to be handled (i.e. given 1-sec record size). THis will be detected
  // and a message returned, to apply ZOH prior:
  //
  //     ZOH osr=1 sr=1

  // If any signals have SR over max value (max-sr), will be flagged: can run this first
  //
  //     ENFORCE-SR dur=1 range=${f1},${f2}

  
  const double max_sr = param.has( "max-sr" ) ? param.requires_dbl( "max-sr" ) : 1024 ; 
  
  
  //
  // annotations to align to: i.e. start each segment at first observation
  //
  
  std::set<std::string> alignments = { "N1", "N2", "N3", "R", "W", "?" };
  
  if ( param.has( "align" ) )
    alignments = param.strset( "align" );
  else if ( param.has( "unaligned" ) )
    alignments.clear();
  
  
  //
  // annot/epoch duration: e.g. assuming (multiples of) fixed 30s epochs
  //   - sensible values: 30 (epochs) 1 (records)

  const uint64_t alignment_unit = param.has( "dur" )
    ? param.requires_dbl( "dur" ) * globals::tp_1sec 
    : globals::tp_1sec * 30;

  if ( alignment_unit % globals::tp_1sec )
    logger << "  *** warning, advised that 'dur' is an integer numbr of records/seconds\n";
  

  
  //
  // only include segments w/ at least some of these annotations (e.g. required=N1,N2,N3,R for at least some sleep scored)
  //  this can be combined w/ any policy (i.e. if largest doesn't have any staging, then nothing will be emitted)
  //  this is also distinct from alignment (i.e. adjusting segment times to match), which can also require that these must exist
  //  i.e. this means
  //
  
  std::set<std::string> requirements;
  if ( param.has( "require" ) )
    requirements = param.strset( "require" );
  else if ( param.has( "require-whole" ) )
    requirements = param.strset( "require-whole " );
  
  const bool requirement_whole = param.has( "require-whole" );
  const bool requirement_min_dur = param.has( "require-dur" );
  const uint64_t requirement_unit = requirement_min_dur ? (uint64_t)(param.requires_dbl( "require-dur" ) * globals::tp_1sec) : 0LLU;  
  if ( requirement_unit % globals::tp_1sec )	
    logger << "  *** warning, advised that 'require-dur' is an integer numbr of records/seconds\n";


  
  //
  // output filename root
  //  --> <out>.edf
  //      <out>.annot

  if ( ( ! param.has( "out" ) ) || param.empty( "out" ) )
    Helper::halt( "requires 'out' argument to specify filename root" );
  
  const std::string out_root = Helper::expand( param.requires( "out" ) );

  //
  // segment annotation prefix (for new, added annots) 
  //

  const std::string aprefix = param.has( "prefix" ) ? param.value( "prefix" ) : "";

      

  // --------------------------------------------------------------------------------
  
  //
  // gap policy
  //

  // segments=all                                    0
  //          largest                               -1
  //          1,2,3,4 for 1st, 2nd, 3rd, etc         2  ( save in keepsegs )

  // policy=splice (default)                       1
  //        zero-pad                               0

  
  int join_policy = 1; // splice
  
  if ( param.has( "policy" ) && ! param.empty( "policy" ) )
    {
      const std::string p =  param.value( "policy" ) ;
      if ( p == "0" || p == "zero-pad" || p == "zero" || p == "pad" ) join_policy = 0;
      else if ( p == "splice" ) join_policy = 1; // default
    }


  // --------------------------------------------------------------------------------
  
  //
  // which segments to retain?
  //

  int segment_policy = 0; // all
  
  std::set<int> keeps;
  
  if ( param.has( "segments" ) && ! param.empty( "segments" ) )
    {

      const std::string p =  param.value( "segments" ) ;

      if      ( p == "largest" ) segment_policy = -1;
      else if ( p == "all" ) segment_policy = 0;
      else
	{
	  std::vector<std::string> tok = Helper::parse( p , "," );
	  for (int i=0;i<tok.size();i++)
	    {
	      
	      int k;
	      if ( Helper::str2int( tok[i] , &k ) )
		if ( k >= 1 )
		  keeps.insert( k );
	    }
	}
      if ( keeps.size() != 0 ) segment_policy = 2; // look-up in keeps
    }


  // --------------------------------------------------------------------------------  
  //
  // Report parameters
  //

  logger << "\n  settings:\n"
	 << "     join-policy (policy)                   = " << ( join_policy == 0 ? "zero-pad" : "splice" ) << "\n"
	 << "     retained segments (segments)           = " << ( segment_policy == 0 ? "all" : param.value( "segments" ) ) << "\n"
	 << "     maximum sample rate allowed (max-sr)   = " << max_sr << " Hz\n"
	 << "     segment alignment annotations (align)  = " << Helper::stringize( alignments ) << "\n"
	 << "       alignment duration unit (dur)        = " << alignment_unit / globals::tp_1sec << "s\n"
	 << "     required annotations (require)         = " << Helper::stringize( requirements ) << "\n"
	 << "       require whole annots (require-whole) = " << ( requirement_whole ? "T" : "F" ) << "\n"
	 << "       require at least (require-dur)       = " << requirement_unit / globals::tp_1sec << "s\n"      	 
	 << "     annotation prefix (prefix)             = " << aprefix << "\n"
	 << "     output file-root (out)                 = " << out_root << "\n";
    
    
  // --------------------------------------------------------------------------------  
  //
  // information on the current EDF
  //
  
  std::set<interval_t> segments1 = timeline.segments();
  std::map<interval_t,int> seg2num;
  
  const bool gapped = segments1.size() != 1;

  if ( gapped && header.continuous )
    Helper::halt( "internal inconsistency: gapped EDF is marked as continuous" );

  const uint64_t rec_size_tp = header.record_duration_tp;

  // --------------------------------------------------------------------------------  
  //
  // forcing selection of segments (i.e. if not all) - if so, do that first, so that
  // gaps are properly defined (i.e. non-selected segments become part of gaps)
  //

  std::set<interval_t> segments;

  if ( keeps.size() )
    {
      logger << "\n  initial segment retention:\n";
      int segn = 1;
      std::set<interval_t>::const_iterator kk = segments1.begin();
      while ( kk != segments1.end() )
	{
	  const bool okay = keeps.find( segn ) != keeps.end() ;
	  logger << "  seg #" << segn << ": " << kk->as_string(2,"-") << " " << ( okay ? "[retained]" : "[skipped]" ) << "\n";
	  seg2num[ *kk ] = segn; // for output below
	  if ( okay )
	    segments.insert( *kk );
	  ++segn;
	  ++kk;
	}
    }
  else
    segments = segments1; // copy all
  
  
  // --------------------------------------------------------------------------------  
  //
  // store all signals
  //

  std::map<std::string, std::vector<double> > sdat;
  std::map<std::string, int> sr;
  std::map<int, std::vector<uint64_t> > tdat;

  signal_list_t signals = header.signal_list( param.value( "sig" ) );
  
  const int ns = signals.size();

  for (int s=0; s<ns; s++)
    {
      if ( header.is_annotation_channel( signals(s) ) ) continue;
      
      slice_t slice( *this , signals(s) , timeline.wholetrace() );

      const std::string slab = signals.label(s);
	
      sdat[ slab ] = *slice.pdata();

      // enforce sample-rate rules

      int    nsamples = header.n_samples[ signals(s) ];
      double fs = (double)nsamples / header.record_duration;

      // too low?
      if ( fs < 1 )
	{
	  writer.value( "MSG" , "sample rate <1Hz for " + slab + ": run 'ZOH osr=1 sr=1' first" );
	  writer.value( "OKAY" , (int)0 );
	  return false;
	}

      // too high?
      if ( fs > max_sr )
	{
	  writer.value( "MSG" , "sample rate >" + Helper::dbl2str(max_sr) + ": set max-sr higher ");
	  writer.value( "OKAY" , (int)0 );
	  return false;
	}

      // non-integer value of samples (for a 1-second record size)
      int new_nsamples1 = fs; 
      if ( fabs( (double)new_nsamples1 - fs ) > 0.000001 )
	{
	  writer.value( "MSG" , "non-integer # of samples per record (given "
			+ Helper::dbl2str(fs) + "Hz for " + slab + "): run RESAMPLE first" );
	  writer.value( "OKAY" , (int)0 );
          return false;
	}
      
      // otherwise, okay
      sr[ slab ] = fs;


      //
      // track implied timepoints (uniquify for a given sample rate)
      //
      
      if ( tdat.find( new_nsamples1 ) == tdat.end() )
	tdat[ new_nsamples1 ] = *slice.ptimepoints();
            
    }


  
  // --------------------------------------------------------------------------------  
  //
  // Get all annotations
  //

  std::set<instance_idx_t> events0;

  std::map<std::string,int> ecounts;
  int alignment_counts = 0;
  int req_counts = 0;
  int nonstandard_alignment_annot_counts = 0;
  
  std::vector<std::string> anames = timeline.annotations.names();  
  
  for (int a = 0 ; a < anames.size() ; a++ )
    {
      
      annot_t * annot = timeline.annotations.find( anames[a] );
      
      if ( annot == NULL ) continue;
      
      const bool align_annot = alignments.find( anames[a] ) != alignments.end();
      const bool req_annot = requirements.find( anames[a] ) != requirements.end();
      
      annot_map_t::const_iterator ii = annot->interval_events.begin();
      while ( ii != annot->interval_events.end() )
        {
          const instance_idx_t & instance_idx = ii->first;
          events0.insert( instance_idx );
	  ecounts[ instance_idx.parent->name ]++;
	  if ( align_annot )
	    {
	      ++alignment_counts;
	      // check whether exact multiple of alignment unit (i.e. 30s)
	      if ( instance_idx.interval.duration() % alignment_unit )
		++nonstandard_alignment_annot_counts;
	    }
	  if ( req_annot ) ++req_counts;
	  ++ii;
        }
    }
  
  logger << "\n  dataset contains " << sr.size() << " signals and "
	 << anames.size() << " annotation classes ("
	 << events0.size() << " instances)\n";

  
  if ( alignments.size() )
    {
      logger << "  specified " << alignments.size()
	     << " annotation classes ("
	     << Helper::stringize( alignments ) 
	     << ") for alignment ("
	     << alignment_counts << " instances found)\n";
      
      if ( nonstandard_alignment_annot_counts ) 
	logger << "  *** warning, " << nonstandard_alignment_annot_counts
	       << " alignment annotations not multiples of standard 'dur'\n";     
      
    }
  
  
  if ( requirements.size() )
    {

      logger << "  specified " << requirements.size()
             << " required annotation classes ("
	     << Helper::stringize( requirements )
             << " with " << req_counts << " instances found)\n";
      
      if ( nonstandard_alignment_annot_counts )
        logger << "  *** warning, " << nonstandard_alignment_annot_counts
               << " alignment annotations not multiples of standard 'dur'\n";
      
    }
  
  
  if ( req_counts == 0 && requirements.size() != 0  )
    {
      writer.value( "MSG" , "no 'require' annotations found" );
      writer.value( "OKAY" , (int)0 );
      return false;
    }
  
  
 
  // --------------------------------------------------------------------------------  
  //
  // select & edit segments to retain (based on gap policies, and/or annotation alignment)
  //
  
  std::set<interval_t> retained;
  std::map<interval_t,interval_t> orig2edit;
  std::map<interval_t,int> included;
  
  // consider each segment
  int sidx = 0;
  
  // track largest (segment_policy -1)
  interval_t largest(0LLU,0LLU);
  
  std::set<interval_t>::const_iterator ii = segments.begin();
  while ( ii != segments.end() )
    {
      
      interval_t seg = *ii; 
      interval_t seg0 = *ii; // orig
      
      //
      // if requiring 1+ from a set of annotations, check we have them here:
      //

      if ( requirements.size() )
	{
	  if ( ! edf_minus_helper_has_annot( events0 , seg , requirement_whole, requirement_unit, requirements ) )
	    {
	      logger << "  skipping segment " << sidx + 1 << " as it does not meet annotation requirements\n";
	      ++ii;
	      ++sidx;
	      orig2edit[ seg0 ] = seg;
	      included[ seg0 ] = 0;
	      continue;
	    }
	}
      
      //
      // add this segment?
      //
      
      if ( segment_policy == 0 || segment_policy == 2 ) 
	{
	  // all (0) or requested (2 - already done filtering)
	  
	  interval_t original_seg = seg;
	  
	  // nudge start forward for any alignment?
	  edf_minus_helper_align( events0 , seg, alignments, alignment_unit, &seg );

	  if ( seg != original_seg ) 
	    logger << "  aligned segment " << sidx + 1 << " : "
		   << original_seg.as_string( 2, "-" ) << " --> "
		   << seg.as_string( 2, "-" ) << "\n";
	      
	  // add to list
	  retained.insert( seg );

	  // track
	  orig2edit[ seg0 ] = seg ;
	  included[ seg0 ] = 1;
	}
      
      // retain only the largest
      if ( segment_policy == -1 )
	{

	  interval_t original_seg = seg;
	  
	  // edit first
          edf_minus_helper_align( events0 , seg, alignments, alignment_unit, &seg );

	  if ( seg != original_seg ) 
	    logger << "  aligned segment " << sidx + 1 << " : "
		   << original_seg.as_string( 2, "-" ) << " --> "
		   << seg.as_string( 2, "-" ) << "\n";

	  // is this largest
	  if ( seg.duration() > largest.duration() )
	    largest = seg;
	  
	  // track anyway
	  included[ seg0 ] = 0; // we'll update the actual largest below (-->1)
	  orig2edit[ seg0 ] = seg;
	}
      
      // next segment
      ++sidx;
      ++ii;
    }
  
  
  //
  // add largest single segment
  //
  
  if ( segment_policy == -1 && largest.duration() != 0LLU )
    {
      included[ largest ] = 1; // redundant prob - but keep, can use for other stats (e.g. mins of annot)
      retained.insert( largest );
    }
    

  
  //
  // No retained (valid) segments?
  //

  if ( retained.size() == 0 )
    {
      writer.value( "MSG" , "no valid retained segments found" );
      writer.value( "OKAY" , (int)0 );
      return false;
    }
  

   
  // --------------------------------------------------------------------------------    
  //
  // Checks
  //
    
  int nsigs = sdat.size();
  
  

  
  // --------------------------------------------------------------------------------  
  //
  // Find start time & date for the new EDF
  //  a) i.e. find first segment
  //  b) these will already be annot-aligned, if that was requested
  //  b) shift to exact # of seconds (for EDF start time)
  //  c) update start time & date values
  
  // get start of first segment
  
  interval_t first = *retained.begin();
  
  uint64_t start_tp = first.start;
  
  // EDF start times are required to be an integer number of seconds
  
  uint64_t frac_tp;

  start_tp = edf_minus_helper_whole_sec( start_tp , &frac_tp );
  
   //
  // EDF start date: note, 1.1.85 is null date, so is not advanced even if
  //   aligned annots imply a new start next day.
  // i.e. use 2.1.85 as the null if you want to track days
  //      SET-HEADERS start-date=2.1.85 
  //
  
  clocktime_t startdatetime( header.startdate, header.starttime );
  
  clocktime_t starttime( header.starttime );
  
  if ( ! startdatetime.valid )
    {
      
      if ( ! starttime.valid )
	{
	  logger << "  *** invalid EDF start time: setting start to [01-01-85-00:00:00]\n";
	  startdatetime.parse_string( "01-01-85-00:00:00" );
	}
      else
	{	  
	  startdatetime.parse_string( "01-01-85-" + header.starttime ); // valid time, enter dummy date
	  logger << "  *** invalid EDF start date: setting start to [" << startdatetime.as_string() << "\n";
	}
    }

  std::string stime1 = startdatetime.as_string();
  std::string sdate1 = startdatetime.as_date_string();
  
  startdatetime.advance_tp( start_tp );
  
  std::string stime2 = startdatetime.as_string();
  std::string sdate2 = startdatetime.as_date_string();

  
  
  // --------------------------------------------------------------------------------
  //
  // Create new annotations and signals
  //
  
  std::set<instance_idx_t> events1; // --> all annots in new EDF/annot
  std::map<std::string, std::vector<double> > sdat1; // --> all signals in new EDF

  
  //
  // splice mode
  //

  std::map<uint64_t,int64_t> splice_deltas;  // define here, as used below in annots
  
  if ( join_policy == 1 ) 
    {

      //   - signals stay 'as is'
      //   - we alter all annotations
      //   - but trim annotations so they fall in a segment
      //     (i.e. in-gap annotations are left removed)
        
      // at or after point 'key' and shift 'value' to annots
      // this is constructed based on the gaps
      
      interval_t prior;
      std::set<interval_t>::const_iterator qq = retained.begin();
      while ( qq != retained.end() )
	{
	  // seed on start of this segment, store gap from prior 
	  if ( qq == retained.begin() )
	    splice_deltas[ qq->start ] = qq->start;
	  else
	    splice_deltas[ qq->start ] = qq->start - prior.stop ;

	  //	  std::cout << " segment offset = " << splice_deltas[ qq->start ] << "\n";

	  prior = *qq;
	  ++qq;
	}
      
      //
      // annotations
      //

      // iterate over segments
      std::set<interval_t>::const_iterator ii = retained.begin();
      
      std::set<instance_idx_t>::const_iterator evt = events0.begin();      

      int64_t offset = splice_deltas[ ii->start ]; // any initial offset?
      
      while ( ii != retained.end() )
	{

	  // iterate over (sorted) annotations (evt)
	  
	  // if no events (e.g. empty), then all done
	  if ( evt == events0.end() )
	    break;
	  
	  // current segment
	  interval_t seg = *ii;

	  // need to head to the next segment
	  if ( evt->interval.start >= seg.stop )
	    {
	      // next segment
	      ++ii;
	      
	      // all done?
	      if ( ii == retained.end() ) break;
	      
	      // update offset (cumulative) 
	      offset += splice_deltas[ ii->start ];
	      
	      // loop back
	      continue;
	    }
	  
	  // does this event belong here?
	  const bool overlap_complete = evt->interval.is_completely_spanned_by( seg );
	  const bool overlap_any      = overlap_complete || evt->interval.overlaps( seg );
	  
	  if ( overlap_complete || overlap_any )
	    {
	      
	      instance_idx_t e1 = *evt;
	      
	      // truncate?
	      if ( ! overlap_complete )
		{
		  if ( e1.interval.start < seg.start ) e1.interval.start = seg.start;
		  if ( e1.interval.stop > seg.stop ) e1.interval.stop = seg.stop;
		}
	      
	      // adjust
	      e1.interval.shift_left( offset );
	      
	      // std::cout << " mapping " << seg.as_string() << "\t"
	      // 		<< evt->parent->name << "\t" << evt->interval.as_string() << " --> "
	      // 		<< e1.interval.as_string() << " " << offset << "\n";
	      
	      // store
	      events1.insert( e1 );
	      
	    }
	  
	  // advance to next annotation
	  ++evt;
	  
	  // loop back 
	}
      

      //
      // Signals
      //
      
      std::map<std::string, std::vector<double> >::const_iterator ss = sdat.begin();
      while ( ss != sdat.end() )
	{

	  // tp-map
	  int Fs = sr[ ss->first ];
	  const std::vector<uint64_t> & tp = tdat[ Fs ];
	  const std::vector<double> & x  = ss->second;
	  const int n = x.size();
	  
	  // create a new signal
	  std::vector<double> x1;

	  // same logic as for annotations, i.e. iterate over
	  // segments, then samples

	  int p = 0;

	  // iterate over segments
	  std::set<interval_t>::const_iterator ii = retained.begin();
	  
	  while ( ii != retained.end() )
	    {
	      // no more signal	      
	      if ( p == n ) break;
	  
	      uint64_t pos = tp[p];
	      
	      // need to head to the next segment
	      if ( pos >= ii->stop )
		{
		  // next segment
		  ++ii;
		  
		  // all done?
		  if ( ii == retained.end() ) break;
		  
		  // loop back
		  continue;
		}
	      
	      // in segment?
	      if ( pos >= ii->start && pos < ii->stop )
		x1.push_back( x[p] );
	      
	      // next sample
	      ++p;
	      
	      // loop back 
	    }

	  // store
	  sdat1[ ss->first ] = x1;

	  // std::cout << " sig " << ss->first << "\t"
	  // 	    << ss->second.size() << " to " << x1.size() << "\n";

	  // next signal
	  ++ss;
	}
      
    }
  

  // --------------------------------------------------------------------------------
  //
  // zero-pad mode
  //
  
  std::map<uint64_t,int64_t> zpad_deltas;
  std::map<uint64_t,int> zpad_recs; 

  if ( join_policy == 0 ) 
    {
      
      //   - same as splice mode, except fill gaps w/ zeros
      //     after expanding/contracting the gap to the nearest
      //     record (1s) unit

      //   - all 'alignment' annotations stay as is (i.e. typically staging, and that
      //     we are doing this because staging is constant, i.e. so we don't want to change
      //     it, rather we'll stretch out the signals
      
      //   - however, all non-alignment annotations will be shifted by a small delta, to reflect the
      //     increase/decrease of the samples, i.e. so that any micro-events stay w/ the recording;
      //
      //   - so, naturally, this induces small changes which may impact some edge cases - e.g. a small
      //     annotation may now fall into a different stage, if the <1 sec shift changes boundaries...

      //   - given that typically we'll just be dealing w/ stage annotations at this point, should not be
      //     a big deal... but perhaps we can think about alternative policies to change this, etc, downstream

      
      interval_t prior;
      int segn=1;
      std::set<interval_t>::const_iterator qq = retained.begin();
      while ( qq != retained.end() )
	{
	  
	  int64_t gap;
	  
	  // seed on start of this segment, store gap from prior 
	  if ( qq == retained.begin() )
	    gap = qq->start;
	  else
	    gap = qq->start - prior.stop ;

	  // fixed 1-second record (rounds down)
	  int nrecs = gap / globals::tp_1sec ; 
	  
	  // do we need to shrink or stretch the gap to the
	  // nearest number of whole (1-sec) records?
	  int64_t diff = gap - ( (int64_t)nrecs * (int64_t)globals::tp_1sec );
	  	  
	  // e.g   2.22 seconds
	  //    --> diff = 0.22; as < 0.5, we should shrink
	  //       2.99
	  //    --> diff = 0.99; as > 0.5 we should expand (by 3-2.99 = 0.01 seconds) 
	  
	  const bool stretch = diff > (int64_t)globals::tp_1sec * 0.5 ;
	  
	  // add an extra record if we are closer to a complete record
	  // otherwise, we'll end up shrinking (removing) the partial record
	  // if less than 0.5 seconds has been observed
	  
	  if ( stretch )
	    {
	      ++nrecs;
	      diff = (int64_t)nrecs * (int64_t)globals::tp_1sec - gap ; 
	    }
	  else
	    diff = -diff; // for shrinking, store as a negative
	   
	  // if shrinking/stretching a segment, flag to console
	  if ( diff )
	    {
	      logger << "\n";
	      if ( segn == 1 )
		logger << "  aligned gap before first segment";
	      else
		logger << "  aligned gap between segments " << segn-1 << " and " << segn;
	      
	      logger << " is not a multiple of 1s (EDF record size)";
	      
	      if ( stretch ) 
		logger << ", so stretching by " << diff / (double)globals::tp_1sec << "s\n";
	      else
		logger << ", so shrinking by " << -diff / (double)globals::tp_1sec << "s\n";
	      
	      logger << "  subsequent annots will be shifted to align w/ signals (except alignment-annots)\n";
	      logger << "\n"; 
	    }

	  
	  // store: delta is the difference between the original gap
	  // and the whole-record sized gap (i.e. which might be negative)
	  zpad_recs[ qq->start ] = nrecs;
	  zpad_deltas[ qq->start ] = diff;
	  
	  // next segment
	  prior = *qq;
	  ++qq;
	  ++segn;
	}
      
      
      //
      // annotations
      //

      // iterate over segments
      std::set<interval_t>::const_iterator ii = retained.begin();
      
      std::set<instance_idx_t>::const_iterator evt = events0.begin();      
    
      int64_t offset = zpad_deltas[ ii->start ]; // any initial offset?
      
      while ( ii != retained.end() )
	{
	  
	  // iterate over (sorted) annotations (evt)
	  
	  // if no events (e.g. empty), then all done
	  if ( evt == events0.end() )
	    break;
	  
	  // current segment
	  interval_t seg = *ii;

	  // need to head to the next segment
	  if ( evt->interval.start >= seg.stop )
	    {
	      // next segment
	      ++ii;
	      
	      // all done?
	      if ( ii == retained.end() ) break;
	      
	      // update offset (cumulative) 
	      offset += zpad_deltas[ ii->start ];
	      
	      // loop back
	      continue;
	    }
	  
	  // now we don't care whether the annotation overlapped a segment or not
	  // i.e. as we are keeping the entire timeline;  however, we do need to
	  // decide whether we should shift non-alignment annotations
	  
	  const bool aligment_annot = alignments.find( evt->parent->name ) != alignments.end();
	  
	  instance_idx_t e1 = *evt;
	  
	  // adjust non-alignment annots
	  // here offset might be negative, so we'll select to either shift
	  // left or right (as these functions expected an unsigned value)
	  
	  if ( ! aligment_annot )
	    {
	      //std::cout << " shifting " << evt->parent->name << " by " << offset << "\n";
	      if ( offset < 0 )
		e1.interval.shift_left( -offset ) ;
	      else if ( offset > 0 ) 
		e1.interval.shift_right( offset );
	    }
	  
	  // store
	  events1.insert( e1 );
	  
	  
	  // advance to next annotation
	  ++evt;
	  
	  // loop back 
	}
      
      
      //
      // Signals
      //      
      
      std::map<std::string, std::vector<double> >::const_iterator ss = sdat.begin();
      while ( ss != sdat.end() )
	{

	  // tp-map
	  int Fs = sr[ ss->first ];
	  const std::vector<uint64_t> & tp = tdat[ Fs ];
	  const std::vector<double> & x  = ss->second;
	  const int n = x.size();
	  
	  // create a new signal
	  std::vector<double> x1;

	  // same logic as for annotations, i.e. iterate over
	  // segments, then samples
	  
	  int p = 0;

	  // iterate over segments
	  std::set<interval_t>::const_iterator ii = retained.begin();

	  // zero-pad before first segment?
	  // i.e. if we are keeping the timeline constant,
	  // this is unlike the splice-join case

	  int nrecs = zpad_recs[ ii->start ] ;
	  for (int i=0; i< nrecs * Fs; i++) // i.e. assumes rec = 1s                                                
	    x1.push_back( 0.0 );
	  
	  
	  while ( ii != retained.end() )
	    {
	      // no more signal	      
	      if ( p == n ) break;
	  
	      uint64_t pos = tp[p];
	      
	      // need to head to the next segment
	      if ( pos >= ii->stop )
		{
	  
		  // next segment
		  ++ii;
		  
		  // all done?
		  if ( ii == retained.end() ) break;

		  // zero-pad?
		  int nrecs = zpad_recs[ ii->start ] ;
		  for (int i=0; i< nrecs * Fs; i++) // i.e. assumes rec = 1s
		    x1.push_back( 0.0 );
		  
		  // loop back
		  continue;
		}
	      
	      // in segment?
	      if ( pos >= ii->start && pos < ii->stop )
		x1.push_back( x[p] );
	      
	      // next sample
	      ++p;
	      
	      // loop back 
	    }

	  // store
	  sdat1[ ss->first ] = x1;

	  // std::cout << " sig " << ss->first << "\t"
	  // 	    << ss->second.size() << " to " << x1.size() << "\n";

	  // next signal
	  ++ss;
	}
      
    }

  
  
  // --------------------------------------------------------------------------------
  //
  // a few checks
  //

  int nr = -1;
  
  std::map<std::string,std::vector<double> >::const_iterator ss1 = sdat1.begin();
  while ( ss1 != sdat1.end() )
    {
      int Fs = sr[ ss1->first ];
      int implied_nr = ss1->second.size() / Fs;
      if ( ss1->second.size() % Fs )
	{
	  logger << "  *** problem - incomplete record found for " << ss1->first << "\n";
	  Helper::halt( "internal error" );
	}

      if ( nr == -1 ) nr = implied_nr;
      else if ( nr != implied_nr )
	{
	  logger << "  *** problem - varying record count across signals\n";
	  Helper::halt( "internal error" );	      
	}
    
      ++ss1;
    }

  
  // --------------------------------------------------------------------------------
  //
  // Console outputs (also determine the actual final placements and show that)
  //
  
  logger << "\n  found " << segments.size() << " segment(s)\n";
  if ( join_policy == 1 )
    logger << "    [ original segments ] -> [ aligned, editted ] --> [ final segments ]\n";
  else 
    logger << "    [ original segments ] --> [ aligned, editted final segments ]\n";
  
  double duration_secs0 = 0 , duration_secs1 = 0;

  std::map<interval_t,interval_t> gaps; // gap before this segment (w.r.t. original)
  std::map<interval_t,interval_t> gaps_edit; // gap before this segment (in final edit)
  std::map<interval_t,interval_t> spliced; // placed segment (in new EDF, if spliced)

  // any initial zpad offsets?
  int64_t offset = zpad_deltas[ ii->start ];
    
  uint64_t last = 0LLU;
  uint64_t last_edit = 0LLU;
  uint64_t running = 0LLU;
  
  sidx = 1;
  std::set<interval_t>::const_iterator ss = segments.begin();
  while( ss != segments.end() )
    {
     
      interval_t orig = *ss;
      interval_t edit = orig2edit[ orig ];
      bool inc = included[ orig ];

      // any gap prior to this? [ note this is accurate for splice-mode only
      // for zero-padding, we may have tweaked gap dur to nearest record ]
      gaps[ orig ] = interval_t( last , orig.start );
      gaps_edit[ orig ] = interval_t( last_edit , edit.start );
      
      // new placed values (duration from the 'editted' version) 
      spliced[ orig ] = interval_t( running , running + edit.duration() );

      // track total duration
      duration_secs0 += orig.duration_sec();
      
      if ( included[ orig ] )
	{
	  duration_secs1 += edit.duration_sec() ;

	  // if zero-padding, then the final also includes the gaps
	  if ( join_policy == 0 ) duration_secs1 += gaps[ orig ].duration_sec();
	}
      

      // outputs

      // gap first
      interval_t g = gaps[ orig ];

      //offset += zpad_deltas[ ii->start ];

      if ( g.duration() ) { // don't call 0-dur gap before first seg a 'gap'                                                         
	
	logger << "    - gap #" << sidx << " : " << g.as_string( 2, "-" ) << " (" << g.duration_sec() << "s)";
	
	if ( join_policy == 1 ) logger << " [spliced]";
	else
	  {
	    logger << " [zero-padded";
	    if ( zpad_deltas[ edit.start ] != 0 ) logger << ", w/ shift " << zpad_deltas[ edit.start ] / (double)globals::tp_1sec << "s";
	    logger << "]";
	  }

	if ( join_policy == 0 )
	  {
	    logger << " --> " << gaps_edit[ orig ].as_string( 2 , "-" );
	    double diff = gaps_edit[ orig ].duration_sec() - gaps[ orig ].duration_sec();
	    if ( diff < -0.001 ) 
	      logger << " (" << -diff << "s shorter)";
	    else if ( diff > 0.001 )
	      logger << " (" << diff << "s longer)";
	  }
	logger << "\n";
      }
	
      
      // seg
      logger << "   " << ( included[ orig ] ? "+" : " " )
	     << "+ seg #" << sidx << " : " << ss->as_string( 2 , "-") << " (" << orig.duration_sec() << "s)";
      
      if ( included[ orig ] )
	{
	  double reduction = orig.duration_sec() - spliced[ orig ].duration_sec();
	  
	  logger << " [included] --> ";
	  
	  if ( join_policy == 1 ) 
	    logger << edit.as_string( 2, "-" ) << " --> " << spliced[ orig ].as_string( 2 , "-" );
	  else
	    logger << edit.as_string( 2 , "-" ); // add offset??

	  if ( fabs( reduction ) > 0.001 )
	    logger << " (" << reduction << "s shorter)";
	  
	  logger << "\n";

	}
      else
	logger << " [excluded]\n";
      
      
      // update
      running += edit.duration();
      last = orig.stop;
      last_edit = edit.stop;
	
      // placement?
      ++ss; ++sidx;
    }
  
  logger << "  original total duration = " << duration_secs0 << "s\n"
	 << "  retained total duration = " << duration_secs1 << "s";
  if ( duration_secs0 > duration_secs1 )
    logger << " (" << duration_secs0 - duration_secs1 << "s shorter)\n" ;
  else if ( duration_secs1 > duration_secs0 )
    logger << " (" << duration_secs1 - duration_secs0 << "s longer)\n" ;



  
  // --------------------------------------------------------------------------------
  //
  // Set up the new EDF
  //
  
  logger << "\n  creating a new EDF " << out_root << ".edf with " << nsigs << " channels\n";
   
  // new EDF outputs
  
  if ( stime1 != stime2 )
    logger << "  updating EDF start-time from " << stime1 << " to " << stime2 << "\n";
  else
    logger << "  retaining original EDF start-time of " << stime2 << "\n";

  if ( sdate1 != sdate2 )
    logger << "  updating EDF start-date from " << sdate1 << " to " << sdate2 << "\n";
  else
    logger << "  retaining original EDF start-date of " << sdate2 << "\n";

  
  edf_t e;
  
  e.init_empty( header.patient_id ,
   		nr ,        // estimated just above
   		1 ,         // record size fixed at 1 second
		startdatetime.as_date_string('.',2) , // date only, YY format
		startdatetime.as_string('.') ); // time only
  
  //
  // add signals
  //

  std::map<std::string, std::vector<double> >::const_iterator xx = sdat1.begin();
  while ( xx != sdat1.end() )
    {
      int Fs = sr[ xx->first ];
      e.add_signal( xx->first , Fs , xx->second );		   	
      ++xx;
    }

  //
  // copy over transducer info and other header information
  //

  e.header.recording_info = header.recording_info;
    
  xx = sdat1.begin();
  while ( xx != sdat1.end() )
    {      
      const int slot0 = header.signal( xx->first );
      const int slot1 = e.header.signal( xx->first );      
      e.header.transducer_type[ slot1 ] = header.transducer_type[ slot0 ];
      e.header.phys_dimension[ slot1 ] = header.phys_dimension[ slot0 ];
      e.header.prefiltering[ slot1 ] = header.prefiltering[ slot0 ];
      ++xx;
    }


  //
  // add annotations
  //
  
  std::set<std::string> anns;

  std::set<instance_idx_t>::const_iterator ee = events1.begin();
  while ( ee != events1.end() )
    {
      anns.insert( ee->parent->name );
      ++ee;
    }
  
  std::set<std::string>::const_iterator aa = anns.begin();
  while ( aa != anns.end() )
    {
      annot_t * annot = e.timeline.annotations.add( *aa );
      std::set<instance_idx_t>::const_iterator ee = events1.begin();
      while ( ee != events1.end() )
	{	  
	  if ( ee->parent->name == *aa )
	    annot->add( ee->id , ee->interval , ee->ch_str );
	  ++ee;
	}
      ++aa;
    }

  
  if ( events1.size() )
    logger << "  creating annotation file " << out_root << ".annot with "
	   << events1.size() << " annotations from " << anns.size() << " classes\n";


  // --------------------------------------------------------------------------------
  //  
  // Outputs
  //

  //
  // Core segments
  //
  
  int segn = 1;
  std::map<interval_t,interval_t>::const_iterator oo = orig2edit.begin();
  while ( oo != orig2edit.end() )
    {
      writer.level( segn, globals::segment_strat );
      writer.value( "ORIG" , oo->first.as_string() );
      writer.value( "EDIT" , oo->second.as_string() );
      writer.value( "INCLUDED" , included[ oo->first ] );
      writer.value( "DUR_ORIG" , oo->first.duration_sec() );
      writer.value( "DUR_EDIT" , oo->second.duration_sec() );
           
      ++oo;
      ++segn;
    }
  writer.unlevel( globals::segment_strat );



  // --------------------------------------------------------------------------------
  //
  // Annotation counts per segment
  //

  sidx = 0;
  ii = segments.begin();
  while ( ii != segments.end() )
    {

      writer.level( sidx+1 , globals::segment_strat );
      
      std::map<std::string,int> rr_tot;
      std::map<std::string,int> rr = edf_minus_helper_count_annots( events0 , *ii ,
								    requirement_whole, requirement_unit,
								    requirements , &rr_tot ) ;
      
      std::map<std::string,int> ra_tot;
      std::map<std::string,int> ra = edf_minus_helper_count_annots( events0 , *ii ,
								    false, alignment_unit,
								    alignments , &ra_tot) ;
      
      std::map<std::string,int>::const_iterator qq = rr.begin();
      while ( qq != rr.end() )
	{
	  writer.level( qq->first , globals::annot_strat );
	  writer.value( "N_ALIGN" , qq->second );
	  writer.value( "N_ALL" , rr_tot[ qq->first ] );
	  ++qq;
	}
      
      qq = ra.begin();
      while ( qq != ra.end() )
	{
	  writer.level( qq->first , globals::annot_strat );
          writer.value( "N_REQ" , qq->second );
	  // this may overwrite rr_tot[] if same annot given, but fine, is the same value
	  // defined twiced as we may have annots in ra_tot not in rr_tot and vice versa
          writer.value( "N_ALL" , ra_tot[ qq->first ] ); 	  
	  ++qq;
	}
      writer.unlevel( globals::annot_strat );
      ++sidx;
      ++ii;
    }
  writer.unlevel( globals::segment_strat );
  
  

  // --------------------------------------------------------------------------------
  //
  // Create new book-keeping annotations 
  //
  
  
  // segments + breakpoints
  
  annot_t * annot_segs = e.timeline.annotations.add( aprefix + "segment" );
  annot_t * annot_gaps = e.timeline.annotations.add( aprefix + "gap" );

  //  interval_t seg1 = 

  sidx = 1;
  int gidx = 1;
  oo = orig2edit.begin();
  while ( oo != orig2edit.end() )
    {
      
      interval_t orig = oo->first;
      interval_t edit = oo->second;

      if ( ! included[ orig ] )
	{
	  ++oo;
	  ++sidx;
	  continue;
	}
      
      // zpad-delta
      if ( join_policy == 0 )
	{
	  // segment (editted) and store original in meta
	  instance_t * s1 = annot_segs->add( Helper::int2str( sidx ) , edit , "." ); 	  
	  s1->set( "orig" , orig.as_string() );

	  // zero-padded gap
	  int grecs = zpad_recs[ edit.start ];
	  int64_t g = zpad_deltas[ edit.start ];	  
	  if ( grecs || g ) { // don't call 0-dur gap before first seg a 'gap'
	    int64_t put_gap = edit.start - grecs * globals::tp_1sec + g;
	    if ( put_gap < 0 ) Helper::halt( "internal error in writing annots, pls contact luna.remnrem@gmail.com" );
	    // gap goes up until the start of this editted segment
	    interval_t zgap( put_gap , edit.start );
	    instance_t * g1 = annot_gaps->add( Helper::int2str( gidx ) , zgap , "." );            
            g1->set( "orig_dur" , grecs );
	    g1->set( "adj" , (double)g / (double)globals::tp_1sec );
	    ++gidx;
          }
	  
	}
      
      // splice-mode
      if ( join_policy == 1 )
	{
	  interval_t s = spliced[ orig ];
	  // segment (editted) and store original in meta
	  instance_t * s1 = annot_segs->add( Helper::int2str( sidx ) , s , "." ); 	  
	  s1->set( "orig" , orig.as_string() );
	  s1->set( "edit" , edit.as_string() );
	  
	  // gap show breakpoint and where any sub-second padding was expanded or stretched
	  // (0-duration time-point at start)
	  interval_t g = gaps[ orig ];
	  if ( g.duration() ) { // don't call 0-dur gap before first seg a 'gap' 
	    instance_t * g1 = annot_gaps->add( Helper::int2str( gidx ) , interval_t( s.start, s.start ), "." );
	    g1->set( "gap" , g.as_string() );
	    g1->set( "dur" , g.duration_sec() );
	    ++gidx;
	  }
	}	  

      ++oo;
      ++sidx;
    }
    



  
  // --------------------------------------------------------------------------------
  //
  // save as a standard EDF 
  //

  bool saved = e.write( out_root + ".edf" ,
			false ,  // standard EDF, not EDFZ
			1 ,      // should be standard EDF so no need to force
			false ,  // always EDF+D = false
			NULL );  // no need to set channel order here


  
  // --------------------------------------------------------------------------------
  //
  // save annotations
  //

  param_t param_write_annots;

  // pass on options for dhms or hms timing options (vs elapsed secs)
  if ( param.has( "dhms" ) )
    param_write_annots.add( "dhms" );  
  else if ( param.has( "hms" ) )
    param_write_annots.add( "hms" );  

  e.timeline.annotations.write( out_root + ".annot" , param_write_annots , e );
    


  // --------------------------------------------------------------------------------
  //
  // all done
  //

  writer.value( "OKAY" , (int)1 );

  return true;
}



// ------------------------------------------------------------
// helper functions

// find start of first instance of an annot 'a' from list 'e' that
// falls within segment 's'
// optionally, impose units-ize (e.g. integer # of 30s epochs from start)
// to prune segment

bool edf_minus_helper_align( const std::set<instance_idx_t> & e ,
			     const interval_t & s ,
			     const std::set<std::string> & a ,
			     const uint64_t d , 
			     interval_t * t )
{

  if ( e.size() == 0 || a.size() == 0 ) return false;
  
  std::set<instance_idx_t>::const_iterator ee = e.begin();
  while ( ee != e.end() )
    {      
      if ( ee->interval.start >= s.stop )
	return false;
      if ( ee->interval.start >= s.start )
	{
	  if ( ee->parent != NULL && a.find( ee->parent->name ) != a.end() )
	    {
	      // set start
	      t->start = ee->interval.start;
	      // also impose a fixed/whole numbr of epochs?
	      if ( d != 0 )
		{
		  // new implied segment duration (now we've aligned w/ start of first annot)
		  uint64_t td = t->duration();
		  int ne = td / d;
		  t->stop = t->start + ne * d;
		}
	      return true;
	    }
	}
      ++ee;
    }
  return false;
}

uint64_t edf_minus_helper_whole_sec( uint64_t tp , uint64_t * f )
{
  *f = tp % globals::tp_1sec;
  return tp - *f;
}

bool edf_minus_helper_has_annot( const std::set<instance_idx_t> & e ,
                                 const interval_t & s ,
				 const bool whole , // requires whole
				 const uint64_t d ,  // requires at least 
                                 const std::set<std::string> & a )
{

  // can require the whole annot and/or at least 'd' tps of the fit
  // for at least one annot in 'a' in this segment 's'
  
  if ( e.size() == 0 || a.size() == 0 ) return false;
  
  std::set<instance_idx_t>::const_iterator ee = e.begin();
  while ( ee != e.end() )
    {
      // gone too far
      if ( ee->interval.start >= s.stop )
	return false;

      // not in set
      if ( ee->parent != NULL && a.find( ee->parent->name ) == a.end() )
	{
	  ++ee;
	  continue;
	}
      
      // require at least some overlap
      if ( s.overlaps( ee->interval ) )
	{
	  
	  // whole annot completely contained in segment
	  const bool fit_whole = ee->interval.is_completely_spanned_by( s );

	  // duration requirement (may have this as well as whole req.)
	  const uint64_t overlap = s.overlap( ee->interval );
	  const bool fit_dur     = overlap >= d;

	  if ( whole )
	    {
	      if ( fit_whole && fit_dur )
		return true;
	    }
	  else
	    {
	      if ( fit_dur )
		return true;
	    }
	}

      // keep searching...
      ++ee;
    }
  // nothing found matching criteria
  return false;
}



std::map<std::string,int> edf_minus_helper_count_annots( const std::set<instance_idx_t> & e ,
							 const interval_t & s ,
							 const bool whole , // requires whole
							 const uint64_t d ,  // requires at least 
							 const std::set<std::string> & a ,
							 std::map<std::string,int> * tots )
{
  
  std::map<std::string,int> r;
  
  //init
  std::set<std::string>::const_iterator aa = a.begin();
  while ( aa != a.end() ) { r[*aa] = 0; ++aa; } 
  
  std::set<instance_idx_t>::const_iterator ee = e.begin();
  while ( ee != e.end() )
    {
      // gone too far
      if ( ee->interval.start >= s.stop )
	break;;
      
      // not in set                                                                                                                           
      if ( ee->parent != NULL && a.find( ee->parent->name ) == a.end() )
        {
          ++ee;
          continue;
        }
      
      // require at least some overlap
      if ( s.overlaps( ee->interval ) )
	{

	  // total count
	  (*tots)[ ee->parent->name ]++;

	  // whole annot completely contained in segment
	  const bool fit_whole = ee->interval.is_completely_spanned_by( s );
	  
	  // duration requirement (may have this as well as whole req.)
	  const uint64_t overlap = s.overlap( ee->interval );
	  const bool fit_dur     = overlap >= d;

	  
	  if ( whole )
	    {
	      if ( fit_whole && fit_dur )
		r[ ee->parent->name ]++;
	    }
	  else
	    {
	      if ( fit_dur )
		r[ ee->parent->name ]++;
	    }
	}

      ++ee;
    }

  return r;

}
