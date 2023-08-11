
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

extern writer_t writer;
extern logger_t logger;

// implementations of various timeline-based functions - or those
// that work with combinations of annots/signals/epochs/records
// and so do not obviously fit somewhere else
//

//  S2A
//  A2S
//  SPANNING
//  ANNOTS
//  MEANS

// also:
//  internal annot2sp() function (used by spindles analysis)

void timeline_t::annot2signal( const param_t & param )
{
  // create a new signal based on one or more annotations
  if ( ! param.has( "annot" ) ) Helper::halt( "no annotations specified: e.g. annot=A1,A2" );
  std::vector<std::string> anames = param.strvector( "annot" );

  // SR of new signals
  const int sr = param.requires_int( "sr" );

  // use instance ID as a numeric value e.g. for NREMC 1 , 2 , 3 
  const bool numeric_instance = param.has( "numeric-inst" );
  
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

	  double value = 1;
	  if ( numeric_instance )
	    {
	      if ( aa->first.id == "" || aa->first.id == "." )
		value = 0;
	      else
		{
		  // set value to numeric
		  if ( ! Helper::str2dbl( aa->first.id , &value ) )
		    Helper::halt( "requires numeric instance IDs" ); 
		}
	    }
	  
	  // populate (up to and including the start/stop, as we removed the final +1 TP above)
	  for (int p=start; p<=stop; p++)
	    adat[p] = value;

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
	seconds -= minutes * 60.0;

      
      //
      // write as a new signal
      //
      
      logger << "  adding " << events.size() << " "
	     << anames[a] << " annotations (spanning ";
      if ( minutes > 0 ) logger << minutes << " min " << seconds << " sec)";
      else logger << seconds << " sec)";

      if ( numeric_instance )
	logger << " as numeric instance-ID signal " << labels[a] << "\n";
      else
	logger << " as 0/1 signal " << labels[a] << "\n";
      
      edf->add_signal( labels[a]  , sr , adat );

    }
      
}



void timeline_t::signal2annot( const param_t & param )
{

  //
  // signal to use
  //

  std::string signal_label = param.requires( "sig" );
  
  signal_list_t signals = edf->header.signal_list( signal_label );
  
  if ( signals.size() == 0 ) Helper::halt( "could not find any signals: " + signal_label );

  const int ns = signals.size();
  
  
  //
  // S2A encoding
  //
  
  // encoding=LABEL,lwr,upr
  // encoding=label,val,+win
  // bins=min,max,n

  //  VALUE :  X    // --> X+EPS
  //           X-Y
  //           X+Y  // eps

  if ( ! ( param.has( "encoding" ) || param.has( "encoding2" ) || param.has( "bins" ) ) )
    Helper::halt( "no encoding=label,value,... or encoding2=label,value1,value2,... or bins=min,max,n" );

  bool e2 = param.has( "encoding" );
  bool e3 = param.has( "encoding2" );
  bool eb = param.has( "bins" );
  if ( e2 + e3 + eb > 1 ) Helper::halt( "must either specify encoding or encoding2 or bins");
  const std::string bin_label = param.has( "bin-label" ) ? param.value( "bin-label" ) : "B" ; 

  std::vector<std::string> enc; 
  int nxy = -1;

  if ( e2 ) 
    {
      enc = param.strvector( "encoding" );
      nxy = 2;
    }
  else if ( e3 ) 
    {
      enc = param.strvector( "encoding2" );
      nxy = 3;
    }
  else
    {
      // make 'encoding2' style string
      std::vector<double> b = param.dblvector("bins");
      if ( b.size() != 3 ) Helper::halt( "expecting bins=min,max,n" );
      const int n = b[2];
      const double bmin = b[0];
      const double bmax = b[1];
      if ( bmin >= bmax || n == 0 ) Helper::halt( "expecting bins=min,max,m" );
      const double binc = (b[1] - b[0] ) / (double)n;
      
      nxy = 3;
      enc.clear();
      for (int i=0; i<n; i++)
	{
	  enc.push_back( bin_label + Helper::int2str( i+1 ) ); // annotate label
	  enc.push_back( Helper::dbl2str( bmin + i * binc ) );
	  enc.push_back( Helper::dbl2str( bmin + (i+1) * binc ) );
	}
    }

  if ( enc.size() % nxy != 0 )
    Helper::halt( "requires " + Helper::int2str( nxy) + " args per encoding value" );


  //
  // Either make one annot class (and labels are instances)
  //  or each label --> a distinct class
  //

  bool use_class = param.has( "class" );
  
  std::string class_name = use_class ? param.value( "class" ) : "" ; 

  //
  // Span EDF discontinuities or no?
  //

  bool span_disc = param.has( "span-gaps" );

  
  //
  // Parse encodings
  //

  std::map<std::string,std::pair<double,double> > e;

  for (int i=0; i<enc.size(); i += nxy )
    {
      std::string label = enc[i];
      
      double ex = 0 ;
      if ( ! Helper::str2dbl( enc[i+1] , &ex ) )
	Helper::halt( "bad numeric value for encoding" + enc[i+1] );

      // default window
      bool window = true;
      double ey = 0.05;
      
      if ( e3 )
	{
	  if ( ! Helper::str2dbl( enc[i+2] , &ey ) )
	    Helper::halt( "bad numeric value for encoding" + enc[i+2] );
	  
	  window = enc[i+2].substr(0,1) == "+" ;
      
	  if ( ! window ) 
	    {
	      if ( ey < ex )
		{
		  double t = ex;
		  ex = ey;
		  ey = t;
		}
	    }
	}

      if ( window )
	{
	  double w = ey;
	  ey = ex + w;
	  ex -= w;	  
	}           
      
      // record
      e[ label ] = std::make_pair( ex , ey );

    }

  logger << "  encoding " << e.size() << " annotation instances\n";

  //
  // For each signal
  //

  for (int s=0; s<ns; s++)
    {

      if ( edf->header.is_annotation_channel( signals(s) ) )
	Helper::halt( "can only use S2A for data channels" );
      
      //
      // get signal data
      //

      slice_t slice( *edf , signals(s) , wholetrace() );  

      std::vector<double> * d = slice.nonconst_pdata();  
      
      const std::vector<uint64_t> * tp = slice.ptimepoints();
      
      //
      // Add annot class?
      //
      
      if ( use_class )
	annotations.add( class_name );
      
      int sr = edf->header.sampling_freq( signals(s) );
      
      std::map<std::string,std::pair<double,double> >::const_iterator ee = e.begin();
      
      while ( ee != e.end() )
	{
	  //std::cout << " label = " << ee->first << "\n";
	  
	  const std::string & label = ee->first; 
	  double ex = ee->second.first;
	  double ey = ee->second.second;
	  
	  // get annot_t to add to
	  annot_t * a = use_class ? annotations.find( class_name ) : annotations.add( label );
	  
	  if ( a == NULL ) Helper::halt( "internal error in signal2annot()" );
	  
	  // iterate over signal points, find in-range intervals
	  
	  const int n = d->size();
	  if ( n == 0 ) {++ee; continue; }
	  
	  bool in = (*d)[0] >= ex && (*d)[0] <= ey;
	  uint64_t start = (*tp)[0];
	  
	  int cnt = 0;
	  
	  for (int i=0; i<n; i++)
	    {
	      // did we just cross a gap, or is this the last data-point?
	      bool gap = span_disc ? false : ( i != 0 ? discontinuity( *tp , sr , i-1 , i ) : false ) ; 
	      
	      // last observed sample?
	      bool end = i == n - 1;
	      
	      // still in region?
	      bool in1 = (*d)[i] >= ex && (*d)[i] <= ey; 
	      
	      // end of an interval? 
	      if ( in && ( gap || end || ! in1 ) ) 
		{	      
		  // 1-past-end encoding
		  uint64_t stop = end ? last_time_point_tp + 1LLU : (*tp)[i] ;
		  a->add( use_class ? label : "." , interval_t( start , stop ) , signals.label(s) );
		  
		  // update status (i.e. may still be a new interval after a gap)
		  in = in1;
		  
		  if ( gap && in1 ) 
		    {
		  start = (*tp)[i];
		  // unlikely, but could be gap and then last single sample
		  if ( end )
		    a->add( use_class ? label : "." , interval_t( start , last_time_point_tp + 1LLU ) , signals.label(s) );		  
		}	      
	      ++cnt;
	    }
	  else if ( in1 && ! in ) // ... or start a new interval
	    {
	      start = (*tp)[i];
	      in = true;
	      if ( i == n - 1 ) // single point interval?
		a->add( use_class ? label : "." , interval_t( start , last_time_point_tp + 1LLU ) , signals.label(s) );
	    }
	}
      
      logger << "  added " << cnt << " intervals for " << label << " based on " << ex << " <= " << signals.label(s) << " <= " << ey << "\n";
      
      // next label
      ++ee;
	}
      
      // next signal
    }

}


void timeline_t::list_spanning_annotations( const param_t & param )
{
    
  if ( mask_set ) 
    Helper::halt( "cannot run SPANNING with a MASK set... use RE" );

  
  // currently, SPANNING only for continuous EDFs
  // if ( ! edf->header.continuous )
  //   Helper::halt( "currently, can only run SPANNING on continuous EDF" );
  
  // given a /set/ of annotations, determine 
  //   - seconds outside of EDF
  //   - total duration of signal covered by these (seconds)
  //   - coverage as a proportion of EDF file
  //   - coverage as a proportion of in-memory representation
  //   - number of contiguous blocks of the requested annotations
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
  
  uint64_t total_all = 0;

  uint64_t total_collapsed = 0;
  
  uint64_t invalid_tps = 0;

  int over_extended = 0;

  int annot_blocks = 0;

  // keep track of where longest spanning annot reaches to 
  // or 0 if past the previous spanning annot
  
  uint64_t earliest = 0;

  uint64_t furthest = 0;
  
  // TODO: track parts of recording that are not spanned
  //std::set<interval_t> unspanned;

  std::set<instance_idx_t>::const_iterator aa = events.begin();
  while ( aa != events.end() )
    {
      
      const interval_t & interval = aa->interval;
      
      //
      // track total (uncollapsed) duration across all ANNOTs
      // i.e. whether valid or not
      //

      total_all += interval.duration();
      
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
	  
	  // only count whole annotations for this total 
	  // (i.e. entire annot must be contained in a contiguous segment
	  //   of the record )

	  total += interval.duration();

	  // start of a 'new' region?
	  
	  if ( furthest == 0 ) 
	    {
	      earliest = interval.start;
	      furthest = interval.stop;
	      ++annot_blocks;
	    }
	  else // we already have at least one region counted
	    {
	      
	      // is the old region finished?  if so, add
	      if ( interval.start > furthest ) 
		{
		  total_collapsed += furthest - earliest ;
		  earliest = interval.start;
		  furthest = interval.stop;
		  ++annot_blocks; // track that this starts a new block
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
  writer.value( "ANNOT_OVERLAP" , ( total_collapsed < total ? "YES" : "NO" )  );

  // how many annots over-extended beyond range of EDF?
  writer.value( "INVALID_N" , over_extended );
  writer.value( "VALID_N" , (int)events.size() - over_extended );
  
  // number of annotation segments, i.e. annotation-based analog of 
  // the SEGMENTS command  
  writer.value( "NSEGS" , annot_blocks );

  // // total annotation duration, whether overlaping or not
  // writer.value( "ALL_ANNOT_SEC" , Helper::tp2sec( total_all )  );
  
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
  
  logger << "  keeping annotations based on ";
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
      logger << "  *** could not find valid start-time in EDF header ***\n";
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
	  //present1.advance( tp1_sec / 3600.0 );
	  present1.advance_seconds( tp1_sec );
	  
	  // add down to 1/100th of a second
	  double tp1_extra = tp1_sec - (long)tp1_sec;

	  // Not sure why we used this form previously... to be consistent, stick with STOP being +1 end
	  
	  // stop_sec_exact() return last time point (rather than usual 1-past-the-end)
	  //double tp2_sec = interval.stop_sec_exact();

	  double tp2_sec = interval.stop_sec();

	  clocktime_t present2 = starttime;
	  //present2.advance( tp2_sec / 3600.0 );
	  present2.advance_seconds( tp2_sec );
	  
	  double tp2_extra = tp2_sec - (long)tp2_sec;
	   
	  writer.value( "START_HMS"  , present1.as_string(':') +  Helper::dbl2str_fixed( tp1_extra , globals::time_format_dp  ).substr(1) );
	  writer.value( "STOP_HMS"   , present2.as_string(':') +  Helper::dbl2str_fixed( tp2_extra , globals::time_format_dp  ).substr(1) );

	  // elapsed time (00:00:00 is start of EDF)
	  clocktime_t present3;
	  present3.advance_seconds( tp1_sec );	  	  
	  // add down to 1/100th of a second
	  tp1_extra = tp1_sec - (long)tp1_sec;

	  clocktime_t present4;	  
	  present4.advance_seconds( tp2_sec );	  
	  tp2_extra = tp2_sec - (long)tp2_sec;
	   	  
	  // std::cout << "xx\t" << present3.as_string(':') 
	  // 	    << "\t" 
	  // 	    << tp1_sec << "\t" 
	  // 	    << (long)tp1_sec << "\t"
	  // 	    << tp1_extra << "\t" 
	  // 	    << Helper::dbl2str_fixed( tp1_extra , globals::time_format_dp ) << "\n";

	  writer.value( "START_ELAPSED_HMS"  , present3.as_string(':') +  Helper::dbl2str_fixed( tp1_extra , globals::time_format_dp ).substr(1) );
	  writer.value( "STOP_ELAPSED_HMS"   , present4.as_string(':') +  Helper::dbl2str_fixed( tp2_extra , globals::time_format_dp ).substr(1) );


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


int timeline_t::annot2sp( edf_t & edf , const std::string & astr ,
			  bool only_this_channel ,
			  std::vector<interval_t> * sample_points , 
			  std::vector<interval_t> * time_points , 
			  int * orig_n , 
			  std::string ch , int sr )
{
  
  sample_points->clear();
  time_points->clear();
  
  // use "CH" to get SR (unless it is otherwise specified)
  //  but read all annots (irrespective of channel) unless only_this_channel == 1
  //  (in which case, we require a specified CH ratehr than a SR (where we simply
  //   find the first match) 

  if ( only_this_channel && ( ch == "" || ch == "." ) )
    Helper::halt( "require a specified channel for annot2sp() " );
  
  // either, find the SR of the given channel:
  if ( sr == 0 )
    {
      signal_list_t signals = edf.header.signal_list( ch );
      if ( signals.size() == 0 ) return 0;
      if ( signals.size() != 1 ) Helper::halt( "problem matching a single channel" );
      std::vector<double> Fs = edf.header.sampling_freq( signals );      
      sr = Fs[0];
    }
  else
    {      
      signal_list_t signals = edf.header.signal_list( "*" );
      std::vector<double> Fs = edf.header.sampling_freq( signals );
      for (int s=0; s<Fs.size(); s++)
	{
	  if ( (int)(Fs[s]) == sr )
	    {
	      ch = signals.label( s );
	      break;
	    }
	}
    }

  
  if ( sr == 0 || ch == "" || ch == "." )
    Helper::halt( "problem finding a channel w/ SR matching" );


  signal_list_t signals = edf.header.signal_list( ch );
  if ( signals.size() != 1 ) Helper::halt( "problem matching a single channel" );
    
  logger << "  using " << ch << " (SR = " << sr << ") to align annotations to sample-points\n";

  
  // must map to within 1 sample (n.b. if at edge, ignored)
  const double max_diff = 1/(double)sr;
  logger << "  mapping to closest sample-point within " << max_diff << " seconds\n";
  
  
  //
  // get the annotation
  //
  
  annot_t * annot = annotations( astr );
  
  if ( annot == NULL )
    Helper::halt( "could not find annotation class " + astr );

  
  
  //
  // get time-points for this SR (via pull of a dummy channel) 
  //

  slice_t slice( edf , signals(0) , edf.timeline.wholetrace() );
  
  const std::vector<uint64_t> * tp = slice.ptimepoints();

  const int np = tp->size();
  
  //
  // Iterate over elements, and build a single ordered table of all times
  //
  
  std::map<uint64_t,int> times;

  *orig_n = 0;
  
  const annot_map_t & events = annot->interval_events;  
  annot_map_t::const_iterator aa = events.begin();
  while ( aa != events.end() )
    {
      // get this annot interval
      const instance_idx_t instance = aa->first;
      
      bool add = ( ! only_this_channel ) || instance.ch_str == ch ; 
      if ( add )
	{
	  // track count
	  (*orig_n)++;

	  times[ instance.interval.start ] = -1;
	  times[ instance.interval.stop ] = -1;
	}
      ++aa;      
    }
  
  //
  // now map all *unique & sorted* times 
  //
  
  // index of tp-map (starts at 1)
  int idx = 1;

  std::map<uint64_t,int>::iterator tt = times.begin();
  while ( tt != times.end() )
    {
      
      const uint64_t & curr = tt->first;
      const uint64_t & prior = (*tp)[idx-1];
      const uint64_t & next = (*tp)[idx];
      
      // shift sample-point window up
      if ( next < curr )
	{	  
	  ++idx;
	  if ( idx == np ) break;
	  continue; // i.e. bounce back but do not update ++tt
	}

      //      std::cout << "\n considering = " << curr << "\n";      
      //std::cout << " idx now = " << idx << " ; curr , prior next = " << curr << " " << prior << " " << next << "\n";
      
      // is in-between these two points?
      if ( curr >= prior && curr <= next )
	{
	  const uint64_t d1 = curr - prior;
	  const uint64_t d2 = next - curr ;
	  
	  const bool first = d1 < d2 ; 
	  
	  //	  std::cout << " closest first = " << first  << "\n";
	  
	  double df = ( first ? d1 : d2 ) * globals::tp_duration ; 
	  
	  //std::cout << "  in nbetween = " << df << "\n";
	  
	  // close enough?
	  if ( df <= max_diff )
	    {
	      int sp = first ? idx-1 : idx ;

	      //std::cout << " storing " << sp << "\n";

	      // store
	      tt->second = sp ;	      
	    }

	}
	  

      // advance to next point 
      
      ++tt;
    }
  
  
  //
  // map back to starts and stops
  //
      
  aa = events.begin();
  while ( aa != events.end() )
    {
      
      int start = -1 , stop = -1;

      const instance_idx_t & instance = aa->first;
      
      bool add = ( ! only_this_channel ) || instance.ch_str == ch ; 
      
      if ( add )
	{
	  if ( times.find( instance.interval.start ) != times.end() )
	    start = times[ instance.interval.start ];
	  
	  if ( times.find( instance.interval.stop ) != times.end() )
	    stop = times[ instance.interval.stop ];
	  
	  // std::cout << " start , stop = " << start << "\t" << stop << " <---- "
	  // 	    << instance.interval.start << "  " << instance.interval.stop << "\n";
	  
	  // add this event (both SP and TP to ensure these are aligned
	  // in the returned value (i.e. to go to makeing a spindle_t ) 
	  if ( start != -1 && stop != -1 )
	    {
	      sample_points->push_back( interval_t( (uint64_t)start , (uint64_t)stop ) );
	      time_points->push_back( instance.interval );
	    }
	}
      
      ++aa;

    }

  
  return sample_points->size();
  
}


void timeline_t::signal_means_by_annot( const param_t & param )
{

  //
  // annots
  //

  if ( ! param.has( "annot" ) ) Helper::halt( "no annotations specified: e.g. annot=A1,A2" );
  std::vector<std::string> anames = param.strvector( "annot" );

  //
  // ignore annotation instannce IDs?
  //

  const bool ignore_instance_ids = ! param.has( "by-instance" );

  //
  // flanking windows
  //

  const bool flanking = param.has( "w" );

  const double flanking_tp = flanking ? param.requires_dbl( "w" ) * globals::tp_1sec : 0 ;
  
  //
  // signals
  //
  
  std::string signal_label = param.requires( "sig" );
  const bool no_annotations = true;
  signal_list_t signals = edf->header.signal_list( signal_label , no_annotations );  
  const int ns = signals.size();

  if ( ns == 0 ) return;
  
  const int Fs = edf->header.sampling_freq( signals(0) );
  for (int s=1; s<ns; s++)
    if ( edf->header.sampling_freq( signals(s) ) != Fs )
      Helper::halt( "signals must have similar sampling rates" );
  
    
  //
  // stores
  //

  // class -> [instance] -> N   [ assumes same SR across channel ]
  // class -> [instance] -> channel -> sum 
  
  std::map<std::string,std::map<std::string,int> > an;
  std::map<std::string,std::map<std::string,std::map<int,double> > > ax;

  // flanking : prior
  std::map<std::string,std::map<std::string,int> > left_an;
  std::map<std::string,std::map<std::string,int> > right_an;

  std::map<std::string,std::map<std::string,std::map<int,double> > > left_ax;
  std::map<std::string,std::map<std::string,std::map<int,double> > > right_ax;

  //
  // iterate over annots
  //


  for (int a=0; a<anames.size(); a++)
    {

      // does annot exist?
      annot_t * annot = edf->timeline.annotations( anames[a] );
      if ( annot == NULL ) continue;
      const std::string & class_name = anames[a];
	  
      // get all events
      const annot_map_t & events = annot->interval_events;

      annot_map_t::const_iterator aa = events.begin();
      while ( aa != events.end() )
	{
	  // instance ID (or not)
	  const std::string inst_id = ignore_instance_ids ? "." : aa->first.id ;
	  
	  // get main interval
	  const interval_t & interval = aa->first.interval;
	  eigen_matslice_t mslice( *edf , signals , interval );	  
          const Eigen::MatrixXd & X = mslice.data_ref();
          const int rows = X.rows();
          const int cols = X.cols();

	  // add to count, accumulate mean
	  an[ class_name ][ inst_id ] += rows;	  
	  Eigen::ArrayXd sum = X.array().colwise().sum();	  
	  for (int s=0; s<ns; s++)
	    ax[ class_name ][ inst_id ][ s ] += sum[ s ];

	  // repeat for flanking regions?
	  if ( flanking )
	    {
	      // left
	      interval_t left = interval;
	      left.shift_left( flanking_tp );	      
	      eigen_matslice_t left_mslice( *edf , signals , left );
	      const Eigen::MatrixXd & left_X = left_mslice.data_ref();
	      left_an[ class_name ][ inst_id ] += left_X.rows();
	      Eigen::ArrayXd left_sum = left_X.array().colwise().sum();
	      for (int s=0; s<ns; s++)
		left_ax[ class_name ][ inst_id ][ s ] += left_sum[ s ];

	      // right
	      interval_t right = interval;
	      right.shift_right( flanking_tp );	      
	      eigen_matslice_t right_mslice( *edf , signals , right );
	      const Eigen::MatrixXd & right_X = right_mslice.data_ref();
	      right_an[ class_name ][ inst_id ] += right_X.rows();
	      Eigen::ArrayXd right_sum = right_X.array().colwise().sum();
	      for (int s=0; s<ns; s++)
		right_ax[ class_name ][ inst_id ][ s ] += right_sum[ s ];	      	      
	    }
	  
	  // next annotation
	  ++aa;
	}

    } // next annotation


  //
  // Report means
  //

  // by annotation class
  std::map<std::string,std::map<std::string,int> >::const_iterator ii = an.begin();
  while ( ii != an.end() )
    {
      writer.level( ii->first , globals::annot_strat );

      // by instance ID 
      
      std::map<std::string,int>::const_iterator jj = ii->second.begin();
      while ( jj != ii->second.end() )
	{

	  // if ignoring instance IDs, then only a single '.' here, so skip
	  // adding as a factor
	  if ( ! ignore_instance_ids ) 
	    writer.level( ii->first , globals::annot_instance_strat );

	  // by channel
	  const std::map<int,double> & chs = ax[ ii->first ][ jj->first ];
	  std::map<int,double>::const_iterator kk = chs.begin();
	  while ( kk != chs.end() )
	    {	      

	      writer.level( signals.label( kk->first ) , globals::signal_strat );

	      // main mean
	      writer.value( "M" , kk->second / (double)jj->second );

	      // flanking regions?
	      if ( flanking )
		{
		  writer.value( "L" , left_ax[ ii->first ][ jj->first ][ kk->first ] / (double)left_an[ ii->first ][ jj->first ] ) ;
		  writer.value( "R" , right_ax[ ii->first ][ jj->first ][ kk->first ] / (double)right_an[ ii->first ][ jj->first ] ) ;
		}

	      // next signal
	      ++kk;
	    }

	  writer.unlevel( globals::signal_strat );
	  ++jj;
	}
      
      if ( ! ignore_instance_ids )
	writer.unlevel( globals::annot_instance_strat );
      
      // next class
      ++ii;	  
    }
  
  writer.unlevel( globals::annot_strat );
  
  // all done
}
