
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

#include "annot/annotate.h"

#include "timeline/timeline.h"
#include "db/db.h"
#include "helper/logger.h"
#include "annot/annot.h"
#include "eval.h"
#include "edf/edf.h"

#include "miscmath/crandom.h"

extern logger_t logger;
extern writer_t writer;


// initiate from a single attached EDF/timeline ( 'OVERLAP' command )
annotate_t::annotate_t( edf_t & edf1 , param_t & param )
{
  edf = &edf1;
  single_indiv_mode = true;
  set_options( param );
  prep();
  loop();
  output();
}

  
annotate_t::annotate_t( param_t & param )
{

  // command-line, multi-sample invocation
  
  // create a 'super' individual, where we simply concatenate segments
  // across individuals; we pull the same annotations from all sets 
  
  // this means we cannot write out annotations as matched/unmatched 
  single_indiv_mode = false;


  // nb. we only allow this is 'bg' mode is specified -- this will
  // implicitly ensure that annotations are only shuffled within
  // individuals

  if ( ! param.has( "bg" ) )
    Helper::halt( "bg specification is required in multi-sample mode" );
  
  // expect a file: indiv -- annot
  
  std::string alist = param.requires( "a-list" );

  std::map<std::string,std::set<std::string> > annots;

  if ( ! Helper::fileExists( alist ) )
    Helper::halt( "could not open " + alist );

  // read each 
  std::ifstream IN1( alist.c_str() , std::ios::in );
  int acnt = 0;
  while ( 1 )
    {
      std::string line;
      Helper::safe_getline( IN1 , line );
      if ( IN1.eof() || IN1.bad() ) break;
      if ( line == "" ) continue;
      std::vector<std::string> tok = Helper::quoted_parse( line , "\t " );
      if ( tok.size() != 2 )
	Helper::halt( "expecting two tab/space delimited fields: ID  annot-file" );
      annots[ tok[0] ].insert( tok[1] );
      ++acnt;
    }
  IN1.close();
  
  logger << "  expecting " << acnt << " annotation files from " << annots.size() << " individuals\n";

  //
  // read in these annotations and make a super-individual annotation file
  //

  // build up a new annotation file, and re-read
  const std::string aggregated = param.requires( "merged" ) + ".annot";

  std::ofstream OUT1( aggregated.c_str() , std::ios::out );
  
  // get size of each 
  std::map<std::string,double> ind2dur;
  
  double offset = 0 ; 
  
  std::map<std::string,std::set<std::string> >::const_iterator aa = annots.begin();
  while ( aa != annots.end() )
    {

      const std::string & indiv = aa->first ;

      logger << "\n  processing " << indiv;
      
      const std::set<std::string> & afiles = aa->second;
      std::set<std::string>::const_iterator bb = afiles.begin();
      while ( bb != afiles.end() )
	{
	  
	  if ( ! Helper::fileExists( *bb ) )
	    Helper::halt( "could not open " + *bb );
	  
	  // notes:
	  //  expects a 'duration_sec' filed in the .annot
	  //  ignores any meta-data
	  //  drops headers

	  logger << " " << * bb ;
	  
	  const bool seen_indiv = ind2dur.find( indiv ) != ind2dur.end();
	  
	  std::ifstream IN1( bb->c_str() , std::ios::in );
	  bool seen_dur = false;
	  while ( 1 )
	    {
	      std::string line;
	      Helper::safe_getline( IN1 , line );
	      if ( IN1.eof() || IN1.bad() ) break;
	      if ( line == "" ) continue;
	      if ( line[0] == '#' ) continue;
	      std::vector<std::string> tok = Helper::quoted_parse( line , "\t " );
	      if ( tok.size() != 6 ) Helper::halt( "expecting standard 6-field annotations:" + line );

	      if ( tok[0] == "duration_sec" )
		{
		  seen_dur = true;
		  
		  double sec;

		  if ( ! Helper::str2dbl( tok[1] , &sec ) )
		    Helper::halt( "problem reading duration_sec field: " + tok[1] );
		  
		  if ( seen_indiv )
		    {
		      if ( fabs( ind2dur[ indiv ] - sec ) > 0.1 )
			Helper::halt( "different duration_sec observed for individual: " + indiv );
		    }
		  else
		    {
		      ind2dur[ indiv ] = sec ; 		      
		    }		  
		}
	      
	      
	      // otherwise skip headers, etc
	      if ( tok[0] == "class" || tok[3] == "." || tok[4] == "." ) continue;

	      // requires we have duration before reading annots
	      if ( ! seen_dur ) continue;

	      // otherwise, we can expect to parse this line
	      //  - assumption, these will be elapsed seconds (i.e. standard WRITE-ANNOTS output form)
	      double start, stop;

	      if ( ! Helper::str2dbl( tok[3] , &start ) )
		Helper::halt( "invalid start field (secs): " + tok[3] );

	      if ( ! Helper::str2dbl( tok[4] , &stop ) )
		Helper::halt( "invalid start field (secs): " + tok[4] );

	      // move forward
	      start += offset;
	      stop += offset;

	      // write out to mega-file

	      OUT1 << tok[0] << "\t"
		   << tok[1] << "\t"
		   << tok[2] << "\t"
		   << start << "\t"
		   << stop << "\t"
		   << tok[5] << "\n";
	      
	    }
	  IN1.close();

	  // next annotation file
	  ++bb;
	}

      logger << "\n"
	     << "  annotations aligned from " << offset << " to " << offset + ind2dur[ indiv ] + 10.0 << " seconds\n";
	
      // shift offset along, with a spacer (arbitrary, 10 seconds)
      // to avoid any flattening of contiguous regions (stop of ind A -- start of ind B)
      offset += ind2dur[ indiv ] + 10.0; 
      
      // next individual
      ++aa;
    }
  
  // done writing
  
  OUT1.close();
  
  //
  // done reading all inputs
  //  -- now create the dummy EDF and attach the mega-file
  //

  //
  // make an empty dummy EDF
  //

  // record size does not matter - set to 100 seconds
  // and allow an extra record at end to make sure all fits

  int nr = (int)( offset / 100 ) + 1 ;
  int rs = 100;
  std::string startdate = "01.01.85";
  std::string starttime = "00.00.00";

  logger << "\n";

  edf_t edfm;
  
  bool okay = edfm.init_empty( "_aggregate_" , nr , rs , startdate , starttime );
  
  if ( ! okay ) Helper::halt( "problem creating the new aggrgete EDF" );
  
  //  
  // load the aggregated file
  //

  edfm.timeline.annotations.set( &edfm );
  
  edfm.load_annotations( aggregated );
  
  //
  // now run as for single-individual mode
  //

  edf = &edfm;

  set_options( param );
  prep();
  loop();
  output();
  
}

void annotate_t::set_options( param_t & param )
{
  
  midpoint = param.has( "midpoint" );

  // for *seeds* only, add flanking values
  flanking_sec = param.has( "f" ) ? param.requires_dbl( "f" ) : 0 ;
  
  window_sec = param.has( "w" ) ? param.requires_dbl( "w" ) : 10 ; 
  
  include_overlap_in_dist = param.has( "dist-includes-overlapping" );
  
  overlap_th = param.has( "overlap" ) ? param.requires_dbl( "overlap" ) : 0 ;

  // flattened all channels to a single event
  pool_channels = param.has( "pool-channels" ) || param.has( "pool-specific-channels" ) ;
  if ( param.has( "pool-specific-channels" ) )
    pool_channel_sets = param.strset( "pool-specific-channels" );
  
  if ( pool_channels ) logger << "  pooling annotations across channels\n";
  else logger << "  retaining channel-level information\n";
  
  // keep channels separate, but permute similarly
  aligned_permutes = param.strset( "align" );
  
  ordered_groups = param.has( "ordered" );
  
  nreps = param.has( "nreps" ) ? param.requires_int( "nreps" ) : 1000 ; 

  fixed.clear();
  if ( param.has( "fixed" ) ) fixed = param.strset( "fixed" );

  chs_inc.clear();
  if ( param.has( "chs-inc" ) )
    proc_chlist( param.value( "chs-inc" ) , true );

  chs_exc.clear();
  if ( param.has( "chs-exc" ) )
    proc_chlist( param.value( "chs-exc" ) , false );

  if ( param.has( "chs-inc" ) && param.has( "chs-exc" ) )
    Helper::halt( "cannot specify by chs-inc and chs-exc lists" );
  
  
  // meta-data filters
  // flt=wgt,10,. 
  
  if ( param.has( "flt" ) )
    {
      std::vector<std::string> w = param.strvector( "flt" );

      // must be divisible by 3
      if ( w.size() % 3 != 0 ) Helper::halt( "expecting flt=label,lwr,upr{,label2,lwr2,upr2}" );
      
      const int nw = w.size() / 3 ;
      int p = 0;
      for (int i=0; i<nw; i++)
	{
	  std::string label = w[p];
	  // lwr?
	  if ( w[p+1] != "." )
	    if ( ! Helper::str2dbl( w[p+1] , &flt_lwr[ label ] ) )
	      Helper::halt( "bad numeric format for " + label + w[p+1] );
	  // upr?
	  if ( w[p+2] != "." )
	    if ( ! Helper::str2dbl( w[p+2] , &flt_upr[ label ] ) )
	      Helper::halt( "bad numeric format for " + label + w[p+2] );

	  p += 3;
	}

      logger << "  set " << flt_lwr.size() << " lower-limits, and " << flt_upr.size() << " upper-limits w/ 'flt'\n";
    }
  
  // checks
  if ( flanking_sec < 0 || window_sec < 0 ) Helper::halt( "invalid negative values for 'f' and/or 'w'" );

  if ( overlap_th < 0 || overlap_th > 1 ) Helper::halt( "invalid value for 'overlap' (0 - 1)" );

  if ( midpoint ) logger << "  reducing all annotations to midpoints\n";
  if ( flanking_sec ) logger << "  adding f=" << flanking_sec << " seconds to each annotation\n";
  if ( window_sec ) logger << "  using window w=" << window_sec << " seconds to search for local intervals\n";
  logger << "  " << ( include_overlap_in_dist ? "" : "not" )
	 << " including overlapping events in nearest-neighbor distances\n";
  
  if ( ordered_groups ) logger << "  ordered=T, so preserving order of seed-seed overlap groups (A,B != B,A)\n";
  else logger << "  ordered=F, so pooling seed-seed permutations, i.e. A,B == B,A (default)\n"; 
  
  // annotations    
  if ( ! param.has( "seed" ) )
    Helper::halt( "require seed argument" );
  
  // requires 1+ seed: look at enrichment of ALL combinations of seeds
  sseeds = param.strset( "seed" );
  
  // from each seed, look at all enrichment w/ all other annots
  //  non-seed annotations are not permuted
  if ( param.has( "annot" ) )
    sannots = param.strset( "annot" );

  // background (i.e. defines the space; only select/permute within contiguous blocks of these regions)
  if ( param.has( "bg" ) ) 
    sbgs = param.strset( "bg" );
    
  // edge for background -- i.e. if elements could not be placed within X seconds of background segment edge,
  // then denote here
  if ( param.has( "edges" ) )
    edge_sec = param.requires_dbl( "edges" );
  
  // outputs - i.e. seed annotations that are w/ or w/out a 'matched' annot
  make_anew = false;
  if ( param.has( "matched" ) )
    {
      if ( param.has( "unmatched" ) )
	Helper::halt( "cannot specify both 'matched' and 'unmatched'" );
      make_anew = true;
      out_include = true;
      out_tag = param.value( "matched" );      
    }
  else if ( param.has( "unmatched" ) )
    {
      make_anew = true;
      out_include = false;
      out_tag = param.value( "unmatched" );
    }

  // 1+ matching... or need more?
  mcount = param.has( "m-count" ) ? param.requires_int( "m-count" ) : 1 ;

  // seed-annot(non-seed) matching only (versus seed-seed)
  seed_nonseed = ! param.has( "seed-seed" ) ;

  
  // unless explicitly specified, do not do perms if
  // getting these outputs
  if ( make_anew && ! param.has( "nreps" ) ) nreps = 0;
  
}



void annotate_t::prep()
{
  
  //
  // Indiv ID
  //

  // for now, assume a single indiv only
      
  iid = edf->id;
  

  //
  // Get annotations from attached timeline
  //

  //
  // backgrounds
  //

  bgs.clear();

  uint64_t edge_tp = edge_sec * globals::tp_1sec;

  std::set<std::string>::const_iterator aa = sbgs.begin();
  while ( aa != sbgs.end() )
    {
      annot_t * a = edf->timeline.annotations.find( *aa );
      if ( a == NULL ) logger << "  ** warning, could not find " << *aa << "\n";
      bgs.insert( a );
      ++aa;
    }
  
  
  //
  // seeds
  //
  
  seeds.clear();
  
  aa = sseeds.begin();
  while ( aa != sseeds.end() )
    {
      annot_t * a = edf->timeline.annotations.find( *aa );
      if ( a == NULL ) logger << "  ** warning, could not find " << *aa << "\n";
      else seeds.insert( a );
      ++aa;
    }


  //
  // other annots
  //

  annots.clear();

  aa = sannots.begin();
  while ( aa != sannots.end() )
    {
      annot_t * a = edf->timeline.annotations.find( *aa );
      if ( a == NULL ) logger << "  ** warning, could not find " << *aa	<< "\n";
      else annots.insert( a );
      ++aa;
    }
  
  if ( seeds.size() == 0 )
    Helper::halt( "no matching seed annotations found" );

  
  //
  // contruct the background set
  //

  std::set<interval_t> abg; // all backgrouns
  
  std::set<annot_t*>::const_iterator bb = bgs.begin();
  while ( bb != bgs.end() )
    {
      annot_t * annot = *bb;
      annot_map_t::const_iterator ii = annot->interval_events.begin();
      while ( ii != annot->interval_events.end() )
        {
	  instance_idx_t instance_idx = ii->first;	  	  	  
	  abg.insert( instance_idx.interval );
	  ++ii;
	}      
      ++bb;
    }

  std::set<interval_t> mbg; // merged backgrounds (merge contiguous regions)

  //
  // combine and flatten all backgrounds
  // 

  if ( abg.size() != 0 )
    {
      mbg = flatten( abg );      
    }

  
  //
  // reduce edges?
  //

  if ( edge_sec )
    {
      std::set<interval_t> b = mbg;
      mbg.clear();
      std::set<interval_t>::const_iterator bb = b.begin();
      while ( bb != b.end() )
	{
	  interval_t b2 = *bb;
	  b2.start += edge_tp;
	  if ( b2.stop - edge_tp < b2.start )
	    b2.stop = b2.start; // null interval
	  else
	    b2.stop -= edge_tp;
	  
	  mbg.insert( b2 );
	  ++bb;
	}
    }

  
  //
  // final summary
  //
  
  if ( mbg.size() != 0 )    
    {

      tottp = total_duration( mbg );
      
      logger << "  background intervals reduced to " << mbg.size()
	     << " contiguous segments, spanning " << tottp * globals::tp_duration << " seconds\n";

      if ( edge_tp != 0 )
	logger << "  background intervals reduced by " << edge_sec << " seconds at edges\n";
    }
  else
    logger << "  no background intervals ('bg'), will assume a single region from 0 to last annotation end-point\n";
  

  
  
  //
  // make a combined set
  //

  std::set<annot_t*> all_annots = seeds;
  
  std::set<annot_t*>::const_iterator pp = annots.begin();
  while ( pp != annots.end() )
    {
      // check not already present
      if ( all_annots.find( *pp ) != all_annots.end() )
	Helper::halt( "cannot specifiy an annotation via both 'seed' and 'annot'" );
      all_annots.insert( *pp );
      ++pp;
    }

  
  //
  // make a list of 'break-points' (brk) - in permutation, we are not
  // allowed to throw down any annotation that will span a break-point
  //

  brk.clear();
  
  if ( mbg.size() == 0 )
    {      
      // find end of last annot
      tottp = 0;
      pp = all_annots.begin();
      while ( pp != all_annots.end() )
	{
	  annot_t * annot = *pp;
	  annot_map_t::const_iterator ii = annot->interval_events.begin();
	  while ( ii != annot->interval_events.end() )
	    {
	      const instance_idx_t & instance_idx = ii->first;
	      if ( instance_idx.interval.stop >= tottp ) tottp = instance_idx.interval.stop;
	      ++ii;
	    }
	  ++pp;
	}

      // 0 1 2 3 4 5 6 7 8 9 10
      // sz = 10,  i.e. random from 0 to 9	  

      seg[ iid ][ 0LLU ] = tottp;
      
      brk.insert( 0LLU );
      brk.insert( tottp );
    }
  else
    {
      std::set<interval_t>::const_iterator bb = mbg.begin();
      while ( bb != mbg.end() )
	{
	  brk.insert( bb->start );
	  brk.insert( bb->stop );	  
	  seg[ iid ][ bb->start ] = bb->duration();
	  ++bb;	  
	}
    }
  
  
  if ( 0 )
    {
      std::cout << "breaks\n";
      std::set<uint64_t>::const_iterator ff = brk.begin();
      while ( ff != brk.end() )
	{
	  std::cout << *ff << "\n";
	  ++ff;
	}
    }


  //
  // any filters to apply?
  //

  const bool filters = flt_lwr.size() || flt_upr.size();
  int filtered_out = 0;
  
  //
  // Pull all intervals and construct the primary datastore
  //

  // indiv -> annot -> set of intervals
  events.clear();
  
  int cnt = 0;
  
  pp = all_annots.begin();
  while ( pp != all_annots.end() )
    {
      
      annot_t * annot = *pp;
      
      const bool is_seed = seeds.find( *pp ) != seeds.end();
      
      annot_map_t::const_iterator ii = annot->interval_events.begin();
      while ( ii != annot->interval_events.end() )
        {

	  const instance_idx_t & instance_idx = ii->first;
	  
	  // annot class ID
	  //  plus/minus the channel
	  
	  bool pool = pool_channels &&
	    ( pool_channel_sets.size() == 0 || pool_channel_sets.find( instance_idx.parent->name ) != pool_channel_sets.end() ) ; 

	  
	  
	  const std::string aid = pool ?
	    instance_idx.parent->name :
	    instance_idx.parent->name + "_" + instance_idx.ch_str ;

	  // skip if not in chs-inc list
	  // or skip if in chs-exc list
	  if ( ! process_channel( instance_idx.parent->name , instance_idx.ch_str ) )
	    {
	      ++ii;
	      continue;
	    }
	  
	  
	  // need to add channel-specific version to seed fix-list?
	  if ( ( ! pool ) && fixed.find( instance_idx.parent->name ) != fixed.end() )
	    {
	      logger << "  adding "
	       	     << instance_idx.parent->name + "_" + instance_idx.ch_str
	       	     << " to fixed list\n";
	      fixed.insert( instance_idx.parent->name + "_" + instance_idx.ch_str );
	    }
	  
	  // track actual AIDs for analysis
	  
	  if ( is_seed ) sachs.insert( aid );
	  achs.insert( aid );

	  // track, if we need to split back (writing new seed annots only)
	  // if explicitly told to pool channels, then we drop this info in
	  // the new output
	  achs_name_ch[ aid ] =
	    std::make_pair( instance_idx.parent->name ,
			    pool ? "." : instance_idx.ch_str );
	  
	  
	  // actual interval

	  interval_t interval = instance_idx.interval;

	  // manipulations?

	  // set to midpoint?
	  
	  if ( midpoint )
	    {
	      uint64_t m = interval.mid();
	      // zero-duration midpoint marker
	      interval.start = interval.stop = m;
	    }
	  
	  // add to the map?
	  
	  uint64_t offset;	  

	  bool okay = segment( iid , interval , &offset );

	  // did not map
	  if ( ! okay ) { ++ii; continue; } 

	  // weights to inc/exc?
	  if ( filters ) 
	    {
	      const instance_t * instance = ii->second;

	      bool include = true;
	      
	      std::map<std::string,double>::const_iterator ll = flt_lwr.begin();
	      while ( ll != flt_lwr.end() )
		{
		  avar_t * mw = instance->find( ll->first ) ;
		  if ( mw != NULL ) 
		    {
		      double w = mw->double_value();
		      if ( w < ll->second ) { include = false ; break; }
		    }
		  ++ll;
		}

	      std::map<std::string,double>::const_iterator uu =	flt_upr.begin();
	      if ( include ) 
		while ( uu != flt_upr.end() )
		  {
		    avar_t * mw = instance->find( uu->first );
		    if ( mw != NULL )
		      {
			double w = mw->double_value();
			if ( w > uu->second ) { include = false ; break; }
		      }   
		    ++uu;
		  }

	      // skip this annotation?
	      if ( ! include ) { ++ii; ++filtered_out; continue; }
	      
	    }
	  
	  // adjust interval by offset
	  // i.e. all times relative to the bounding segment	  

	  if ( offset > interval.start ) Helper::halt( "logic error (1)" );
	  if ( offset > interval.stop ) Helper::halt( "logic error (2)" );

	  // re-register interval relative to bounding segment
	  interval.start -= offset;
	  interval.stop -= offset;
	  
	  // add flanking regions to start/stop to seeds only
	  //  - but do not overstep bounding edges
	  
	  if ( is_seed && flanking_sec > 0 )
	    {
	      uint64_t f = globals::tp_1sec * flanking_sec ;

	      // expand left edge: 
	      interval.start = f < interval.start ? interval.start - f : 0;

	      // expand right edge:
	      // if segment() return T, as above, we are guaranteed to
	      // have duration of bounding segment in seg[ offset ]

	      if ( seg[iid].find( offset ) == seg[iid].end() ) Helper::halt( "logic error 3" );
	      uint64_t dur = seg[ iid ][ offset ];
	      interval.stop = interval.stop + f > dur ? dur : interval.stop + f ;
	    }
	  
	  // add this segment to the primary list
	  events[ iid ][ offset ][ aid ].insert( interval );
	  
	  ++cnt;
	  
	  // next instance
          ++ii;
        }
    
      // next annotation
      ++pp;
    }
  
  
  logger << "  registered " << cnt
	 << " intervals across " << achs.size()
	 << " annotation classes, including " << sachs.size() << " seed(s)\n";

  if ( filters )
    logger << "  excluded " << filtered_out << " of "
	   << filtered_out + cnt
	   << " annotations based on filters, leaving " << cnt << "\n";
  
  
  // review

  if ( 0 )
    {
      std::map<std::string,std::map<uint64_t,std::map<std::string,std::set<interval_t> > > >::const_iterator ee = events.begin();
      
      while ( ee != events.end() )
	{

	  const std::map<uint64_t,std::map<std::string,std::set<interval_t> > > & region = ee->second;
	  std::map<uint64_t,std::map<std::string,std::set<interval_t> > >::const_iterator rr = region.begin();
	  while ( rr != region.end() )
	    {
	      const std::map<std::string,std::set<interval_t> > & annots = rr->second;
	      std::map<std::string,std::set<interval_t> >::const_iterator qq = annots.begin();
	      while ( qq != annots.end() )
		{
		  const std::set<interval_t> & ints = qq->second;
		  std::set<interval_t>::const_iterator ii = ints.begin();
		  while ( ii != ints.end() )
		    {
		      std::cout << " ee-> " << ee->first << "\t" 
				<< rr->first << "\t"
				<< qq->first << "\t"
				<< ii->as_string() << "\n";

		      ++ii;
		    }
		    
		  ++qq;
		}
	      ++rr;
	    }
	  ++ee;	  
	}
    }

  
}



void annotate_t::loop()
{
  
  // evaluate the original dataset
  annotate_stats_t s = eval();
  observed( s );

  if ( make_anew )
    {
      // write out, then clear/turn off this mechanism for null evals
      new_seeds(); 
      hits.clear();
      make_anew = false;
    }
  
  // evaluate 'nreps' shuffled datasets
  for (int r=0; r<nreps; r++)
    {
      // console reports
      if ( r == 0 ) logger << "  ";
      logger << ".";
      if ( r % 50 == 49 ) logger << " " << r+1 << " of " << nreps << " replicates done\n  ";
      else if ( r % 10 == 9 ) logger << " ";
      
      // permute
      shuffle();

      // calc statistics for null data
      annotate_stats_t s = eval();

      // track null distribution
      build_null( s );
    }
  
}



void annotate_t::shuffle()
{

  //
  // shuffle each indiv/seed independently
  //
  
  std::map<std::string,std::map<uint64_t,std::map<std::string,std::set<interval_t> > > >::const_iterator ee = events.begin();

  // individuals
  while ( ee != events.end() )
    {
      
      // contiguous regions
      const std::map<uint64_t,std::map<std::string,std::set<interval_t> > > & region = ee->second;
      
      std::map<uint64_t,std::map<std::string,std::set<interval_t> > >::const_iterator rr = region.begin();
      while ( rr != region.end() )
	{
	  
	  // shuffle (w/ wrapping) each seed annotation independently
	  const uint64_t maxshuffle = seg[ iid ][ rr->first ];
	  
	  // each seed for this region
	  std::set<std::string>::const_iterator ss = sachs.begin();
	  while ( ss != sachs.end() )
	    {
	      
	      // skipping this annotation? (fixed), then nothing to do here
	      if ( fixed.find( *ss ) != fixed.end() )
		{
		  ++ss;		  
		  continue;
		}


	      // get the factor by which to shuffle these annot in this region

	      uint64_t pp = 0;
	      
	      
	      // if using aligned permutations, have we already
	      // generated a shuffle for an aligned annot? if so use
	      // that


	      if ( 0 ) // aligned already 
		{
		  // if ( aligned )
		  // 	{
		  /// TODO...
		}
	      
	      
	      // get a random offset
	      //  - which results in no annots that span the end of this segment
	      //  - keep going unitl we get one... ouch 
	      
	      else 
		{

		  int iter = 0;
		  
		  while ( 1 )
		    {
		      // we having to try too hard?   tells us the data are
		      // not appropriate for this
		      
		      ++iter;
		      if ( iter > 1000 )
			Helper::halt( "cannot find shuffle sets for " + *ss );
		      
		      // putative shuffle 
		      pp = CRandom::rand( maxshuffle );
		      
		      // check all events
		      bool okay = true;
		      const std::set<interval_t> & original = events[ ee->first ][ rr->first ][ *ss ];
		      std::set<interval_t>::const_iterator ii = original.begin();
		      while ( ii != original.end() )
			{
			  interval_t i = *ii;
			  i.start += pp;
			  i.stop += pp;
			  // check - spans segment break?
			  if ( i.start < maxshuffle && i.stop >= maxshuffle )
			    {
			      okay = false;
			      break;
			    }
			  ++ii;
			}
		      
		      // need to try again?
		      if ( okay ) break;
		    }
		  

		  // save this aligned shuffle?
		  
		  
		}

	      
	      //
	      // now we have a valid shuffle value... do the shuffle 
	      //
	      
	      // copy over events
	      std::set<interval_t> original = events[ ee->first ][ rr->first ][ *ss ];
	      
	      std::set<interval_t> shuffled;
	      
	      std::set<interval_t>::const_iterator ii = original.begin();
	      while ( ii != original.end() )
		{
		  
		  interval_t i = *ii;
		  
		  i.start += pp;
		  i.stop += pp;
		  
		  // need to wrap?
		  if ( i.start >= maxshuffle )
		    {
		      i.start -= maxshuffle;
		      i.stop -= maxshuffle;
		    }
		  
		  //std::cout << " chng " << ii->as_string() << " --> " << i.as_string() << "\n";
		  // add this new one to the list
		  shuffled.insert( i );
		  	      
		  ++ii;
		}
	    
	      //std::cout << " resizing " << events[ ee->first ][ rr->first ][ *ss ].size() << " to " << shuffled.size() << "\n";

	      // update
	      events[ ee->first ][ rr->first ][ *ss ] = shuffled;
	      
	      // next seed
	      ++ss;
	    }

	  // next region
	  ++rr;
	}

      // next person
      ++ee;
    }
  
}


annotate_stats_t annotate_t::eval()
{
  // main task: aggregate overlap stats in 'r'
  annotate_stats_t r;

  // secondary: optionally (only w/ true data) output new
  // seed annotations (if 'make_anew' set)
  
  
  std::map<std::string,std::map<uint64_t,std::map<std::string,std::set<interval_t> > > >::const_iterator ii = events.begin();
  
  while ( ii != events.end() )
    {
      
      // each region
      const std::map<uint64_t,std::map<std::string,std::set<interval_t> > > & region = ii->second;
      std::map<uint64_t,std::map<std::string,std::set<interval_t> > >::const_iterator rr = region.begin();
      while ( rr != region.end() )
	{
	  
	  // for seed-pileup
	  std::set<named_interval_t> puints;
	  
	  // each seed for this region
	  std::set<std::string>::const_iterator aa = sachs.begin();
	  while ( aa != sachs.end() )
	    {
	      // does this interval have any seeds?
	      if ( rr->second.find( *aa ) == rr->second.end() ) { ++aa; continue; }
	      
	      // get all seed events
	      const std::set<interval_t> & a = rr->second.find( *aa )->second;
	      
	      // track for pile-up
	      std::set<interval_t>::const_iterator qq = a.begin();
	      while ( qq != a.end() )
		{
		  puints.insert( named_interval_t( *qq , *aa ) );
		  ++qq;
		}
	      
	      // track # of (flattened) annots
	      std::set<interval_t> flata = flatten( a );	      
	      r.ns[ *aa ] += flata.size();
 	      
	      // consider all other annots
	      std::set<std::string>::const_iterator bb = achs.begin();
	      while ( bb != achs.end() )
		{		  
		  // skip self comparison
		  if ( *aa == *bb ) { ++bb; continue; }
		  
		  // does this interval have any seeds?
		  if ( rr->second.find( *bb ) == rr->second.end() ) { ++bb; continue; }
		  
		  // get all other annots
		  const std::set<interval_t> & b = rr->second.find( *bb )->second;
		  		  
		  // calc and record stats on flattened lists 
		  seed_annot_stats( flata , *aa , flatten( b ) , *bb , &r );
		  		  
		  // next annot
		  ++bb;
		}
	      
	      // next seed
	      ++aa;
	    }
	  
	  // seed-seed pileup 	  
	  std::map<std::string,double> pu = pileup( puints );
	  std::map<std::string,double>::const_iterator pp = pu.begin();
	  while ( pp != pu.end() )
	    {
	      if ( 1 || pp->first != "_O1" )
		{
		  r.nss[ pp->first ] += pp->second;
		}
	      ++pp;
	    }
	  
	  // region
	  ++rr;
	}
      // person
      ++ii;
    }


  // final
  // std::map<std::string,double>::const_iterator pp = r.nss.begin();
  // while ( pp != r.nss.end() )
  //   {
  //     std::cout << pp->first << " = " << pp->second << "\n";
  //     ++pp;
  //   }
  // std::cout << "\n\n.......\n\n";

  
  // all done
  return r;
  
}


void annotate_t::output()
{

  //
  // seed-seed group overlap
  //

  std::map<std::string,double>::const_iterator ss = obs.begin();
  while ( ss != obs.end() )
    {
      writer.level( ss->first , "SEEDS" );
      writer.value( "OBS" , obs[ ss->first ] );
      if ( nreps )
	{
	  double mean = exp[ ss->first ] / (double)nreps;
	  double var = expsq[ ss->first ] / (double)nreps - mean * mean;	  
	  writer.value( "EXP" , mean );
	  writer.value( "P" , ( pv[ ss->first ] + 1 ) / (double)( nreps + 1 ) );
	  writer.value( "Z" , ( obs[ ss->first ] - mean ) / sqrt( var ) );
	}
      ++ss;
    }
  writer.unlevel( "SEEDS" );
  
  //
  // seed-* prop overlap
  //

  std::map<std::string,double>::const_iterator pp = prop_obs.begin();
  while ( pp != prop_obs.end() )
    {
      writer.level( pp->first , "SEED" );
      writer.value( "PROP" , prop_obs[ pp->first ] );
      if ( nreps )
        {
          double mean = prop_exp[ pp->first ] / (double)nreps;
          double var = prop_expsq[ pp->first ] / (double)nreps - mean * mean;	  
          writer.value( "PROP_EXP" , mean );
          writer.value( "PROP_P" , ( prop_pv[ pp->first ] + 1 ) / (double)( nreps + 1 ) );
          writer.value( "PROP_Z" , ( prop_obs[ pp->first ] - mean ) / sqrt( var ) ); 
        }
      ++pp;
    }
  writer.unlevel( "SEED" );

  //
  // seed-annot overlap & distances; seed on distance measures, as that will be a
  // superset of the overlap measures
  //
  
  std::map<std::string,std::map<std::string,double> >::const_iterator sa = absd_obs.begin();
  while ( sa != absd_obs.end() )
    {
      writer.level( sa->first , "SEED" );
      const std::map<std::string,double> & p = sa->second;
      std::map<std::string,double>::const_iterator pp = p.begin();
      while ( pp != p.end() )
	{
	  writer.level( pp->first , globals::annot_strat );

	  //
	  // seed-annot overlap: count
	  //
	  
	  if ( p_obs[ sa->first ][ pp->first ] != 0 )
	    {
	      writer.value( "N_OBS" , p_obs[ sa->first ][ pp->first ]  );
	      if ( nreps )
		{
		  double mean = p_exp[ sa->first ][ pp->first ] / (double)nreps;
		  double var = p_expsq[ sa->first ][ pp->first ] / (double)nreps - mean * mean;	  
		  writer.value( "N_EXP" , mean );
		  writer.value( "N_P" , ( p_pv[ sa->first ][ pp->first ]  + 1 ) / (double)( nreps + 1 ) );
		  writer.value( "N_Z" , ( p_obs[ sa->first ][ pp->first ] - mean ) / sqrt( var ) );
		}
	    }
	  	  
	  // seed-annot distances
	  
	  writer.value( "D1_OBS" , absd_obs[ sa->first ][ pp->first ]  );
	  writer.value( "D_N" , dn_obs[  sa->first ][ pp->first ]  );
	  if ( nreps )
	    {
	      double mean = absd_exp[ sa->first ][ pp->first ] / (double)nreps;
	      double var = absd_expsq[ sa->first ][ pp->first ] / (double)nreps - mean * mean;
	      writer.value( "D1_EXP" , mean );
	      writer.value( "D1_P" , ( absd_pv[ sa->first ][ pp->first ] + 1 ) / (double)( nreps + 1 ) );
	      writer.value( "D1_Z" , ( absd_obs[ sa->first ][ pp->first ] - mean ) / sqrt( var ) );

	      writer.value( "D_N_EXP" , dn_exp[  sa->first ][ pp->first ] / (double)nreps  );	      
	      
	    }

	  writer.value( "D2_OBS" , sgnd_obs[ sa->first ][ pp->first ]  );
	  if ( nreps )
	    {
	      double mean = sgnd_exp[ sa->first ][ pp->first ] / (double)nreps;
	      double var = sgnd_expsq[ sa->first ][ pp->first ] / (double)nreps - mean * mean;
	      writer.value( "D2_EXP" , mean );
	      writer.value( "D2_P" , ( sgnd_pv[ sa->first ][ pp->first ] + 1 ) / (double)( nreps + 1 ) );
	      writer.value( "D2_Z" , ( sgnd_obs[ sa->first ][ pp->first ] - mean ) / sqrt( var ) );
	    }
	  
	  ++pp;
	}	  
      writer.unlevel( globals::annot_strat );
      ++sa;
    }
  writer.unlevel( "SEED" );
  
}

 

bool annotate_t::place_interval( const interval_t & i ,  uint64_t * offset ) const 
{
  // to test if interval i spans any point in 'brk' (where 'brk' contains 0 and end also)
  // test whether the start and stop have the same iterator from upper_bound search
  
  std::set<uint64_t>::const_iterator u1 = brk.upper_bound( i.start );

  // for end, check on END - 1, i.e. if end of annot is last of all, we need the end of the GREATER than the last event
  std::set<uint64_t>::const_iterator u2 = brk.upper_bound( i.stop == 0 ? 0 : i.stop - 1LLU ); 

  // if we span a break-point, we know this segment is no good
  if ( u1 != u2 ) return false;

  // but what if we are in a gap? (including starting before or after all segments)
  if ( u1 == brk.begin() || u1 == brk.end() ) return false;

  // track back, and check whether value is a key in seg[]
  --u1;

  // is in gap between two breaks? (i.e. offset for start not tracked) 
  const std::map<uint64_t,uint64_t> & seg2 = seg.find( iid )->second;
  if ( seg2.find( *u1 ) == seg2.end() ) return false;
  
  // seems okay, return offset for bounding segment
  *offset = *u1;
  
  return true;
}


bool annotate_t::segment( const std::string & id , const interval_t & i , uint64_t * segoff ) const
{
  uint64_t offset = 0;
  
  // if spans a break point, or falls in a gap, no good
  if ( ! place_interval( i , &offset ) ) return false;

  // return start to of this containing segment
  *segoff = offset;
  
  return true;
}




std::map<std::string,double> annotate_t::pileup( const std::set<named_interval_t> & allints ) const
{

  std::map<std::string,double> r;
  if ( allints.size() == 0 ) return r;

  std::set<named_interval_t> basket;
  
  std::set<named_interval_t>::const_iterator ii = allints.begin();
  
  // add first
  named_interval_t last = *ii;
  basket.insert( last );
  
  ++ii; 

  while ( ii != allints.end() )
    {
      // std::cout << ii->n << "\t" << ii->i.as_string() << "\n";

      // gap?
      if ( ii->i.start >= last.i.stop ) 
	{
	  // simple count
	  ++r[ "_O" + Helper::int2str( (int)basket.size() ) ];
	  // actual combo
	  ++r[  Helper::int2str( (int)basket.size() ) + ":" + stringize( basket ) ];

	  // start a new basket
	  basket.clear();
	  basket.insert( *ii );

	  // track last added
	  last = *ii;
	  
	}
      else // add / expand last 
	{
	  //std::cout << " expanding...\n";
	  basket.insert( *ii );
	  last.i.stop = last.i.stop > ii->i.stop ? last.i.stop : ii->i.stop;
	}
      
      ++ii;
    }

  // last one
  ++r[ "_O" + Helper::int2str( (int)basket.size() ) ];
  ++r[ Helper::int2str( (int)basket.size() ) + ":" + stringize( basket ) ];
   
  return r;
}

std::string annotate_t::stringize( const std::set<named_interval_t> & t ) const
{

  // is A,B different from B,A ? 
  //  i.e. can test for differential enrichment
  
  if ( ordered_groups )
    {
      std::stringstream ss;
      std::set<named_interval_t>::const_iterator tt = t.begin();
      while ( tt != t.end() )
	{
	  if ( tt != t.begin() ) ss << ",";
	  ss << tt->n;
	  ++tt;
	}
      return ss.str();
    }

  // otherwise, default is to pool all perms
  
  // first reduce to nams, i.e. so B,A --> A,B 
  std::set<std::string> names;
  std::set<named_interval_t>::const_iterator tt = t.begin();
  while ( tt != t.end() )
    {
      names.insert( tt->n );
      ++tt;
    }
  return Helper::stringize( names );
}


void annotate_t::seed_annot_stats( const std::set<interval_t> & a , const std::string & astr , 
				   const std::set<interval_t> & b , const std::string & bstr , 
				   annotate_stats_t * r )
{

  // if no b annots in this segment, then nothing to do
  // stats remain as they are
  
  if ( b.size() == 0 ) return;

  // is 'b' also a seed?
  const bool bseed = sachs.find( bstr ) != sachs.end();
  
  // tmp debug verbose output 
  
  if ( 0 )
    {
      std::cout << " \nshowing seed = " << astr << " ann = " << bstr << "\n";
      std::set<interval_t>::const_iterator iaa = a.begin();
      while ( iaa != a.end() )
	{
	  std::cout << "SEED = " << astr << "\t" << iaa->as_string() << " " << iaa->duration_sec() << "\n";
	  ++iaa;
	}
      
      std::set<interval_t>::const_iterator ibb = b.begin();
      while ( ibb != b.end() )
	{
	  std::cout << "ANNOT = " << bstr << "\t" << ibb->as_string() << "\n";
	  ++ibb;
	}
    }

  
  // consider each seed
  std::set<interval_t>::const_iterator aa = a.begin();
  while ( aa != a.end() )    
    {
      
      double dist = 0;
      bool overlap = false;
      
      // find the first annot not before (at or after) the seed
      std::set<interval_t>::const_iterator bb = b.lower_bound( *aa );
      
      // edge cases:
      // no annot at or past seed? : bb == b.end() 
      // no annot before seed?     : bb == b.begin()

      
      //              |-----| SEED           LB?
      //    |-----|   |     |
      //            |-|-|   |
      //              |--|  |                 Y
      //              |-----|                 Y
      //              |   |-|                 Y
      //              |   |-|--|              Y 
      //              |     |----|            Y
      //              |     |        |---|    Y
      
      // does first take overlap seed?
      
      if ( bb != b.end() && bb->overlaps( *aa ) )
	{
	  // we're done, found complete overlap
	  dist = 0;
	  overlap = true;
	}
      else // it must come afterwards
	{
	  
	  uint64_t seed_start = aa->start;
	  uint64_t seed_stop  = aa->stop - 1LLU;	  
	  
	  // nb. here -1 means that no right distance is defined
	  dist = bb != b.end() ?
	    ( bb->start - seed_stop ) * globals::tp_duration : 
	    -1; 
	  
	  // step back, if we can - is there a closer annot /before/ the seed?
	  if ( bb != b.begin() )
	    {
	      --bb;
	      
	      // nb - this may overlap seed
	      // i.e. starts before, but ends after seed-start, and so
	      //  was not captured by the lower_bound()
	      
	      if ( bb->stop >= aa->start )
		{
		  dist = 0;
		  overlap = true;
		}
	      else
		{

		  // annot stops before start of seed
		  double left_dist = ( aa->start - bb->stop ) * globals::tp_duration; 
		  
		  // closer? (or equal to)
		  // nb. store as signed here (neg -> before)

		  // nb. use dist < 0 condition to track that there was no
		  // lower_bound (i.e. the final annot occurs before this seed)
		  // in which case, this left annot will be the closest 
		  if ( dist < 0 || left_dist <= dist )
		    dist = - left_dist; 
		}	  
	    }
	}

      // track: overlap = dist == 0,
      // but use bool overlap to avoid floating-point equality test
      if ( overlap ) r->nsa[ astr ][ bstr ] += 1 ;
      
      // to track proprtion of seeds w/ at least one (non-seed) annot overlap
      if ( overlap && ! bseed )
	r->psa[ astr ].insert ( *aa );
      
      // for mean distance -- do we meet the window criterion?
      const double adist = fabs( dist );
      if ( adist <= window_sec )
	{
	  // do we include complete overlap as "nearest"?
	  if ( include_overlap_in_dist || ! overlap )
	    {
	      r->adist[ astr ][ bstr ] += fabs( dist );
	      r->sdist[ astr ][ bstr ] += dist;
	      r->ndist[ astr ][ bstr ] += 1; // denom for both the above
	    }
	}
      
      // are we tracking hits
      if ( make_anew )
	{
	  if ( overlap || adist <= window_sec )
	    {
	      // only tracjing seed-nonseed matches? or all?
	      const bool okay = seed_nonseed ? ! bseed : true ; 
	      //	      std::cout << " okay = " << astr << " " << bstr << " = " << okay << " " << overlap << " " << window_sec << " " << adist << "\n";
	      if ( okay )
		{
		  named_interval_t named( *aa , astr );
		  hits[ named ]++;
		}
	    }
	}
      
      // next seed annot
      ++aa;
    }
  
  
}


std::set<interval_t> annotate_t::flatten( const std::set<interval_t> & x )
{

  std::set<interval_t> m;
  
  if ( x.size() == 0 ) return m;
  
  interval_t curr = *x.begin();
  
  std::set<interval_t>::const_iterator xx = x.begin();
  while ( xx != x.end() )
    {
      // because of +1 end encoding, contiguous regions will (a,b) ... (b,c) 
      const interval_t & pro = *xx;
      if ( pro.start > curr.stop )
	{	      
	  m.insert( curr );	      
	  	      
	  // and update the current
	  curr = pro;
	}
      else // expand background as needed
	{
	  if ( pro.stop > curr.stop ) curr.stop = pro.stop;
	}

      // consider next element
      ++xx;
    }
  
  // add final element
  m.insert( curr );

  return m;
}

uint64_t annotate_t::total_duration( const std::set<interval_t> & x )
{
  uint64_t d = 0;
  std::set<interval_t>::const_iterator xx = x.begin();
  while ( xx != x.end() )
    {
      d += xx->duration();
      ++xx;
    }
  return d;
}


void annotate_t::observed( const annotate_stats_t & s )
{
  
  // seed-seed group overlap ( std::map<std::string,double> )
  obs = s.nss;
  
  // seed-annot pairwise overlap (std::map<std::string,std::map<std::string,double> > )
  p_obs = s.nsa;
  
  // seed-annot proportion spanned
  std::map<std::string,std::set<interval_t> >::const_iterator pp = s.psa.begin();
  while ( pp != s.psa.end() )
    {
      prop_obs[ pp->first ] = pp->second.size() / s.ns.find(  pp->first )->second;
      //std::cout <<"  prop_obs = " << prop_obs[ pp->first ] << "\n";
      ++pp;
    }
  
  // absolute distance (from each seed to nearest annot) std::map<std::string,std::map<std::string,double> >
  absd_obs = s.adist;
  
  // signed distance (from each seed to nearest annot) std::map<std::string,std::map<std::string,double> >
  sgnd_obs = s.sdist;
  
  // get average distances (these may be < S-A count, because of window_sec threshold)
  std::map<std::string,std::map<std::string,double> >::const_iterator dd = s.ndist.begin();
  while ( dd != s.ndist.end() )
    {
      const std::map<std::string,double> & e = dd->second;
      std::map<std::string,double>::const_iterator ee = e.begin();
      while ( ee != e.end() )
	{
	  absd_obs[ dd->first ][ ee->first ] /= (double)ee->second;
	  sgnd_obs[ dd->first ][ ee->first ] /= (double)ee->second;
	  dn_obs[ dd->first ][ ee->first ] = (double)ee->second;
	  ++ee;
	}
      ++dd;
    }
  
}


void annotate_t::build_null( const annotate_stats_t & s )
{
  // consider only the observed configurations
  
  // seed-seed group overlap
  std::map<std::string,double>::const_iterator ss = obs.begin();
  while ( ss != obs.end() )
    {
      const bool is_seen = s.nss.find( ss->first ) != s.nss.end();      
      if ( is_seen )
	{
	  double val = s.nss.find( ss->first )->second;
	  exp[ ss->first ] += val;
	  expsq[ ss->first ] += val * val;	  
	  if ( val >= obs[ ss->first ] ) ++pv[ ss->first ];
	}      
      ++ss;
    }

  //
  // seed-annot overlap
  //
  
  std::map<std::string,std::map<std::string,double> >::const_iterator sa = p_obs.begin();
  while ( sa != p_obs.end() )
    {
      // should always be okay, but check just in case some weirdness
      const bool is_seen = s.nsa.find( sa->first ) != s.nsa.end();
      if ( ! is_seen ) { ++sa; continue; }
      
      const std::map<std::string,double> & p = sa->second;
      const std::map<std::string,double> & pe = s.nsa.find( sa->first )->second;
      
      std::map<std::string,double>::const_iterator pp = p.begin();
      while ( pp != p.end() )
	{	  
	  const bool is_seen = pe.find( pp->first ) != pe.end();
	  if ( is_seen )
	    {
	      double val = pe.find( pp->first )->second;
	      p_exp[ sa->first ][ pp->first ] += val;
	      p_expsq[ sa->first ][ pp->first ] += val * val;
	      if ( val >= p_obs[ sa->first ][ pp->first ] ) ++p_pv[ sa->first ][ pp->first ];
	    }
	  ++pp;
	}
      ++sa;
    }

  //
  // prop-seed overlap
  //

  std::map<std::string,double>::const_iterator pp = prop_obs.begin();
  while ( pp != prop_obs.end() )
    {
      const bool is_seen = s.psa.find( pp->first ) != s.psa.end();
      if ( is_seen )
        {
          double val = s.psa.find( pp->first )->second.size() / s.ns.find( pp->first )->second;
	  // std::cout << " prop_exp " << pp->first << " = "
	  // 	    << val << " = " << s.psa.find( pp->first )->second.size() << " / " << s.ns.find( pp->first )->second <<  "\n";
	  //	  std::cout << " val = " << val << "\n";

	  prop_exp[ pp->first ] += val;
          prop_expsq[ pp->first ] += val * val;
          if ( val >= prop_obs[ pp->first ] ) ++prop_pv[ pp->first ];
        }
      ++pp;
    }

  
  //
  // seed-annot distances : sgnd_obs and absd_obs will always have the same keys, so do just once
  //
  
  sa = absd_obs.begin();
  while ( sa != absd_obs.end() )
    {

      const bool is_seen = s.ndist.find( sa->first ) != s.ndist.end();
      if ( ! is_seen ) { ++sa; continue; }
      
      const std::map<std::string,double> & p = sa->second;
      const std::map<std::string,double> & pe_abs = s.adist.find( sa->first )->second;
      const std::map<std::string,double> & pe_sgn = s.sdist.find( sa->first )->second;
      const std::map<std::string,double> & pe_n   = s.ndist.find( sa->first )->second;
      
      std::map<std::string,double>::const_iterator pp = p.begin();
      while ( pp != p.end() )
        {
          const bool is_seen = pe_n.find( pp->first ) != pe_n.end();
          if ( is_seen )
            {
	      const double n = pe_n.find( pp->first )->second;
	      // normalize here:
	      const double a = pe_abs.find( pp->first )->second / n;
	      const double s = pe_sgn.find( pp->first )->second / n;

	      // expecteds: sums
              absd_exp[ sa->first ][ pp->first ] += a;
	      sgnd_exp[ sa->first ][ pp->first ] += s;

	      // sum of sqs
	      absd_expsq[ sa->first ][ pp->first ] += a * a;
	      sgnd_expsq[ sa->first ][ pp->first ] += s * s;

	      // track counts
	      dn_exp[ sa->first ][ pp->first ] += n;
	      
	      // pvals : testing whether *closer* so stat is LE rather than GE
	      if ( a <= absd_obs[ sa->first ][ pp->first ] ) ++absd_pv[ sa->first ][ pp->first ];
	      
	      // nb. although we've calculated mean signed-dist, here the test is 2-sided, so take abs(x)
	      // nb. test if closer (smaller dist) is more significant, thus reversed sign 
	      if ( fabs(s) <= fabs( sgnd_obs[ sa->first ][ pp->first ] ) ) ++sgnd_pv[ sa->first ][ pp->first ];
	    }
          ++pp;
        }
      ++sa;
    }
  
}


void annotate_t::new_seeds()
{

  if ( ! single_indiv_mode )
    {
      logger << "  *** cannot add a new seed annotation when running in multi-individual mode ***\n";      
      return;
    }
  
  // hits[ named_interval_t ] -> # of annots that 'matched'
  // depending on out_include T/F, include of exclude annots w/ 1+ match
  
  std::set<std::string>::const_iterator ss = sachs.begin();
  while ( ss != sachs.end() )
    {
      // get original name/ channel 
      //  nb - here we any instance ID info
      const std::string aname = achs_name_ch[ *ss ].first;
      const std::string chname = achs_name_ch[ *ss ].second;
      
      logger << "  creating new annotation " << aname << out_tag << " ( channel = " << chname << " )\n";
      
      annot_t * a = edf->timeline.annotations.add( aname + out_tag );
      int acnt = 0 , tcnt =0 ;
      
      // iterate over all iids/regions/etc
     interval_map_t::const_iterator ee = events.begin();
      while ( ee != events.end() )
	{
	  // get regions
	  const std::map<uint64_t,std::map<std::string,std::set<interval_t> > > & regions = ee->second;
	  std::map<uint64_t,std::map<std::string,std::set<interval_t> > >::const_iterator rr = regions.begin();
	  while ( rr != regions.end() )
	    {
	      // track offset (i.e. need to add back in for the output)
	      uint64_t offset = rr->first;
	      
	      // pick this seed only
	      std::map<std::string,std::set<interval_t> >::const_iterator this_seed = rr->second.find( *ss );
	      
	      // may not be any seeds in this region - skip to next region
	      if ( this_seed == rr->second.end() ) { ++rr; continue; }

	      // get intervals
	      const std::set<interval_t> & intervals = this_seed->second;
	      
	      std::set<interval_t>::const_iterator ii = intervals.begin();
	      while ( ii != intervals.end() )
		{
		  named_interval_t named( *ii , *ss );

		  // requisite # of hits
		  const bool write_this = out_include ?
		    hits[ named ] >= mcount : // include
		    hits[ named ] <  mcount ; // exclude
		  
		  // nb - we drop any instance ID information
		  // for now... can fix this up later?
		  //  but... given we've a) flattened, and b) perhaps
		  //  pooled channels, this is probably the 'right'
		  //  thing to do...
		  
		  if ( write_this )
		    {
		      interval_t mapped( ii->start + offset , ii->stop + offset );
		      a->add( "." , mapped , chname );
		      ++acnt;
		    }
		  ++tcnt;
		  
		  // next interval
		  ++ii;
		}
	      ++rr;
	    }
	  ++ee;
	}
      
      logger << "   - wrote " << acnt << " (of " << tcnt << ") seed events, based on ";
      if ( ! out_include ) logger << "not ";
      logger << "matching " << mcount << " or more other annots\n";
      
      // all done, next seed
      ++ss;
    }
}

  
void annotate_t::proc_chlist( const std::string & s , const bool inc )
{
  if ( inc ) chs_inc.clear();
  else chs_exc.clear();

  //  annot:ch,annot:ch
  // format chs-inc=SP11:C3,SP11:C4,RIP:LHH1

  std::vector<std::string> tok = Helper::parse( s , "," );
  for (int i=0; i<tok.size(); i++)
    {
      // expect annot:ch
      std::vector<std::string> tok2 = Helper::parse( tok[i] , ":" );
      if ( tok2.size() != 2 ) Helper::halt( "expecting annot:ch format for chs-inc and chs-exc" );

      if ( inc ) chs_inc[ tok2[0] ].insert( tok2[1] );
      else chs_exc[ tok2[0] ].insert( tok2[1] );      
    }
  
}

bool annotate_t::process_channel( const std::string & a , const std::string & ch )
{

  std::map<std::string,std::set<std::string> >::const_iterator aa = chs_inc.find( a );  
  if ( aa != chs_inc.end() )
    {
      const std::set<std::string> & chs = aa->second;
      if ( chs.find( ch ) == chs.end() ) return false;      
    }

  aa = chs_exc.find( a );
  if ( aa != chs_exc.end() )
    {
      const std::set<std::string> & chs	= aa->second;
      if ( chs.find( ch ) != chs.end() ) return	false;
    }
  
  return true;

}
