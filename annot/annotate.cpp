
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
		   << Helper::dbl2str( start , globals::time_format_dp ) << "\t"
		   << Helper::dbl2str( stop , globals::time_format_dp ) << "\t"
		   << tok[5] << "\n";
	      
	    }
	  IN1.close();

	  // next annotation file
	  ++bb;
	}

      logger << "\n"
	     << "   annotations aligned from " << offset << " to " << offset + ind2dur[ indiv ] + 10.0 << " seconds\n";
	
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

  std::vector<std::string> names = edfm.timeline.annotations.names() ; 

  logger << "  have " << names.size() << " annotations in " << aggregated << " :";
  for (int i=0; i<names.size(); i++) logger << " " << names[i] ;
  logger << "\n";
  
  set_options( param );
  prep();
  loop();
  output();
  
}

void annotate_t::set_options( param_t & param )
{

  // number of permutations
  nreps = param.has( "nreps" ) ? param.requires_int( "nreps" ) : 1000 ; 

  // verbose/debug mode
  debug_mode = param.has( "verbose" ) || param.has( "debug" );
  
  // reduce annotations to midpoints
  // midpoint     do all
  // midpoint=SP,SO   do only these (plus w/ channel specifiers)  
  midpoint = false;  
  if ( param.has( "midpoint" ) )
    {
      if ( param.empty( "midpoint" ) )
	midpoint = true; // do all 
      else
	midpoint_annot = param.strset( "midpoint" ); // do some
    }
    
  // for *seeds* only, add flanking values
  flanking_sec = 0;
  if ( param.has( "f" ) )
    {
      // could be f=0.5   means all
      // or       f=A1:0,A2:1,A3:0
      
      std::vector<std::string> tok = param.strvector( "f" );
      
      // single, numeric value
      double fval = 0;
      if ( tok.size() == 1 && Helper::str2dbl( tok[0] , &fval ) )
	flanking_sec = fval;
      else
	{
	  // assume annot:f,annot:f format	  
	  for (int t=0; t<tok.size(); t++)
	    {
	      std::vector<std::string> tok2 = Helper::parse( tok[t] , ":" );
	      if ( tok2.size() != 2 ) Helper::halt( "expecting annot:second format for 'f'" );
	      if ( ! Helper::str2dbl( tok2[1] , &fval ) ) Helper::halt( "expecting annot:second format for 'f'" );
	      flanking_sec_annot[ tok2[0] ] = fval ; 
	    }
	}
    }
  
  // distance to neighbour stats
  window_sec = param.has( "w" ) ? param.requires_dbl( "w" ) : 10 ; 
  
  include_overlap_in_dist = ! param.has( "dist-excludes-overlapping" );
  
  overlap_th = param.has( "overlap" ) ? param.requires_dbl( "overlap" ) : 0 ;
  
  
  //
  // flatten all channels to a single annotation class
  //
  
  pool_channels = param.has( "pool-channels" );
  
  if ( ! param.empty( "pool-channels" ) )
    pool_channel_sets = param.strset( "pool-channels" );
  
  if ( pool_channels ) logger << "  pooling annotations across channels\n";
  else logger << "  retaining channel-level information\n";
  
  
  //
  // only compare within similar channels?
  //  (for unpooled comparisons only)
  //
  
  only_within_channel = param.has( "within-channel" );

  if ( only_within_channel && pool_channels )
    Helper::halt( "cannot specify within-channel and pool-channel together" );

  //
  // Include/exclude particular channels
  //

  chs_inc.clear();
  if ( param.has( "chs-inc" ) )
    proc_chlist( param.value( "chs-inc" ) , true );

  chs_exc.clear();
  if ( param.has( "chs-exc" ) )
    proc_chlist( param.value( "chs-exc" ) , false );
  
  if ( param.has( "chs-inc" ) && param.has( "chs-exc" ) )
    Helper::halt( "cannot specify by chs-inc and chs-exc lists" );
    
  
  //
  // keep channels separate, but permute similarly
  // align=A1,A2|R1|Z1,Z2
  //  here | delimits different sets (R1 may imply all channels within R1)
  //
  
  if ( param.has( "align" ) )
    {
      std::string ap = param.value( "align" );
      std::vector<std::string> ap_tok = Helper::parse( ap , "|" );
      
      for (int aidx=0; aidx<ap_tok.size(); aidx++)
	{
	  // nb. insert self --> self so that channel-expansion works
	  std::vector<std::string> tok = Helper::parse( ap_tok[ aidx ] , "," );
	  const int na = tok.size();
	  for (int i=0; i<na; i++)
	    for (int j=0; j<na; j++)
	      aligned_permutes[ tok[i] ].insert( tok[j] );
	}
      logger << "  using " << ap_tok.size() << " alignment groups\n";
    }


  //
  // also permute annots (as well as seeds)
  //
  
  shuffle_annots = param.has( "shuffle-others" ) || param.has( "shuffle-others" );
  
  //
  // do not permute a particular seed/annotation
  //

  fixed.clear();
  if ( param.has( "fixed" ) ) fixed = param.strset( "fixed" );

  //
  // do seed-seed pileup?
  //
  
  do_pileup = param.has( "pileup" ) ? param.yesno( "pileup" ) : true ; 
  
  //
  // constrained shuffle duration (i.e. only up to X seconds in either direction)?
  //
  
  constrained_shuffle_dur = param.has( "max-shuffle" );
  if ( constrained_shuffle_dur )
    {
      max_shuffle_sec = param.requires_dbl( "max-shuffle" );
      if ( max_shuffle_sec < 0 ) Helper::halt( "max-shuffle must be positive" );
    }
  else
    max_shuffle_sec = 0;


  //
  // Misc
  //
      
  ordered_groups = param.has( "ordered" );
  

  //
  // meta-data filters
  // flt=wgt,10,. 
  //
  
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
  if ( window_sec ) logger << "  truncating distance at w=" << window_sec << " seconds for nearest neighbours\n";
  logger << "  " << ( include_overlap_in_dist ? "" : "not " )
	 << "including overlapping events in nearest-neighbor distances\n";
  
  if ( ordered_groups ) logger << "  ordered=T, so preserving order of seed-seed overlap groups (A,B != B,A)\n";
  else logger << "  ordered=F, so pooling seed-seed permutations, i.e. A,B == B,A (default)\n"; 
  
  // annotations    
  if ( ! param.has( "seed" ) )
    Helper::halt( "require seed argument" );
  
  // requires 1+ seed: look at enrichment of ALL combinations of seeds
  sseeds = param.strset( "seed" );
  
  // from each seed, look at all enrichment w/ all other annots
  //  non-seed annotations are not permuted
  if ( param.has( "other" ) )
    sannots = param.strset( "other" );

  // background (i.e. defines the space; only select/permute within contiguous blocks of these regions)
  if ( param.has( "bg" ) ) 
    sbgs = param.strset( "bg" );
    
  // edge for background -- i.e. if elements could not be placed within X seconds of background segment edge,
  // then denote here
  if ( param.has( "edges" ) )
    edge_sec = param.requires_dbl( "edges" );

  // exclusionary background (xbg) -- i.e. xbg is the converse of bg, and thus
  // specifies gaps rather than allowed intervals
  if ( param.has( "xbg" ) )
    sxbgs = param.strset( "xbg" );

  if ( param.has( "xbg" ) && ! param.has( "bg" ) )
    Helper::halt( "xbg requires bg to be explicitly specified" );
  
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
      else bgs.insert( a );
      ++aa;
    }

  //
  // exclusionary backgrounds
  //
  
  xbgs.clear();
  
  aa = sxbgs.begin();
  while ( aa != sxbgs.end() )
    {
      annot_t * a = edf->timeline.annotations.find( *aa );
      if ( a == NULL ) logger << "  ** warning, could not find " << *aa << "\n";
      else xbgs.insert( a );
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

  std::set<interval_t> abg; // all backgrounds
  
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
  //  ( true --> join adjacent/contiguous intervals for the BG)
  //
  
  if ( abg.size() != 0 )
    {
      mbg = flatten( abg , true );
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
	    {
	      b2.stop = b2.start; // null interval
	      //std::cout << " making a NULL interval...\n";
	    }
	  else
	    b2.stop -= edge_tp;
	  
	  if ( b2.stop > b2.start ) 
	    mbg.insert( b2 );
	  
	  ++bb;
	}

      if ( edge_tp != 0 )
	logger << "  background intervals reduced by " << edge_sec << " seconds at edges\n";
      
    }

  //
  // remove exclusionary backgrounds (i.e. make holes)
  //

  if ( xbgs.size() )
    {
      
      std::set<interval_t> xs; // all excisions
  
      std::set<annot_t*>::const_iterator bb = xbgs.begin();
      while ( bb != xbgs.end() )
	{
	  annot_t * annot = *bb;
	  annot_map_t::const_iterator ii = annot->interval_events.begin();
	  while ( ii != annot->interval_events.end() )
	    {
	      instance_idx_t instance_idx = ii->first;	  	  	  
	      xs.insert( instance_idx.interval );
	      ++ii;
	    }      
	  ++bb;
	}

      // flatten 
      //  ( true --> join adjacent/contiguous intervals for the BG)      
      xs = flatten( xs , true );

      logger << "  excising " << xs.size() << " unique xbg intervals\n";
      
      // and now remove these from the primary background
      mbg = excise( mbg , xs );
      
      if ( mbg.size() == 0 )
	Helper::halt( "no valid background intervals left after exclusions" );
      
    }

  // std::set<interval_t>::const_iterator mm = mbg.begin();
  // while ( mm != mbg.end() )
  //   {
  //     std::cout << "  BB = " << mm->start << " .. " << mm->stop << "\n";
  //     ++mm;
  //   }
    
  
  //
  // final summary
  //
  
  if ( mbg.size() != 0 )    
    {
      
      tottp = total_duration( mbg );
      
      logger << "  background intervals reduced to " << mbg.size()
	     << " contiguous segments, spanning " << tottp * globals::tp_duration << " seconds\n";

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

      seg[ 0LLU ] = tottp;
      
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
	  seg[ bb->start ] = bb->duration();
	  ++bb;	  
	}
    }
  
  
  if ( debug_mode )
    {
      
      std::cout << "background: # of discontinuities = " << brk.size() << "\n";
      std::set<uint64_t>::const_iterator ff = brk.begin();
      while ( ff != brk.end() )
	{
	  std::cout << " background discontinuity tp = " << *ff << "\tsec = " << *ff * globals::tp_duration << "\n";
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
	  
	  // if no channel specified, then always 'pool' (i.e. in simple case, do not add "_." to end
	  if ( instance_idx.ch_str == "." ) pool = true;
	  
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
	  
	  // need to add channel-specific version for 'midpoint'?
	  if ( ( ! pool ) && midpoint_annot.find( instance_idx.parent->name ) != midpoint_annot.end() )
	    {
	      midpoint_annot.insert( instance_idx.parent->name + "_" + instance_idx.ch_str );
	    }
	  
	  // need to add channel-specific version for 'f'?
	  if ( ( ! pool ) && flanking_sec_annot.find( instance_idx.parent->name ) != flanking_sec_annot.end() )
	    {
	      flanking_sec_annot[ instance_idx.parent->name + "_" + instance_idx.ch_str ]
		= flanking_sec_annot[ instance_idx.parent->name ];
	    }
	  
	  // need to add channel-specific version to seed fix-list?
	  if ( ( ! pool ) && fixed.find( instance_idx.parent->name ) != fixed.end() )
	    {
	      // logger << "  adding "
	      //  	     << instance_idx.parent->name + "_" + instance_idx.ch_str
	      //  	     << " to fixed list\n";
	      fixed.insert( instance_idx.parent->name + "_" + instance_idx.ch_str );
	    }

	  // need to add channel-specific version to alignment groups?
	  if ( ( ! pool ) && aligned_permutes.find( instance_idx.parent->name ) != aligned_permutes.end() )
            {

	      const std::string ch_name = instance_idx.parent->name + "_" + instance_idx.ch_str;

	      if ( aligned_permutes.find( ch_name ) == aligned_permutes.end() )
		{
		  // logger << "  adding "
		  // 	 << instance_idx.parent->name + "_" + instance_idx.ch_str
		  // 	 << " to alignment group\n";
		  
		  const std::string root_name = instance_idx.parent->name;		  
		  
		  // update key:
		  aligned_permutes[ ch_name ] = aligned_permutes[ root_name ];
		  // and all other members
		  std::map<std::string,std::set<std::string> >::iterator pp = aligned_permutes.begin();
		  while ( pp != aligned_permutes.end() )
		    {
		      if ( pp->second.find( root_name ) != pp->second.end() )
			pp->second.insert( ch_name );
		      ++pp;
		    }
		}
            }

	  // track with channel this annotaton belongs to 
	       
	  if ( only_within_channel )
	    {
	      label2channel[ aid ] = instance_idx.ch_str ; 
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

	  //
	  // manipulations?
	  //

	  // track original 

	  interval_t original = interval;

	  //
	  // set to midpoint?
	  //
	  
	  if ( midpoint || midpoint_annot.find( aid ) != midpoint_annot.end() )
	    {
	      uint64_t m = interval.mid();
	      // zero-duration midpoint marker
	      interval.start = interval.stop = m;
	    }

	  //
	  // add to the map?
	  //
	  
	  uint64_t offset;	  

	  bool okay = segment( interval , &offset );

	  if ( ! okay ) { ++ii; continue; } 

	  //
	  // weights to inc/exc?
	  //
	  
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

	  bool generic_seed_flank = is_seed && flanking_sec > 0 ;
	  bool annot_specific_flank = flanking_sec_annot.find( aid ) != flanking_sec_annot.end();
	  
	  if ( generic_seed_flank || annot_specific_flank )
	    {
	      double fsec = generic_seed_flank ? flanking_sec : flanking_sec_annot[ aid ];
	      uint64_t f = globals::tp_1sec * fsec ;

	      // expand left edge: 
	      interval.start = f < interval.start ? interval.start - f : 0;

	      // expand right edge:
	      // if segment() return T, as above, we are guaranteed to
	      // have duration of bounding segment in seg[ offset ]

	      if ( seg.find( offset ) == seg.end() ) Helper::halt( "logic error 3" );
	      uint64_t dur = seg[ offset ];
	      interval.stop = interval.stop + f > dur ? dur : interval.stop + f ;
	    }
	  
	  // track the original (used only if outputting matched annots)
	  if ( make_anew ) 
	    unmanipulated[ named_interval_t( offset , interval , aid ) ] = original;
	  
	  // add this segment to the primary list
	  events[ offset ][ aid ].insert( interval );
	  
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


  //
  // some more information on which annotations
  //

  std::map<std::string,int> annot_n;
  std::map<std::string,double> annot_s;
  
  std::map<uint64_t,std::map<std::string,std::set<interval_t> > >::const_iterator rr = events.begin();
  while ( rr != events.end() )
    {
      const std::map<std::string,std::set<interval_t> > & annots = rr->second;
      std::map<std::string,std::set<interval_t> >::const_iterator qq = annots.begin();
      while ( qq != annots.end() )
	{
	  const std::set<interval_t> & ints = qq->second;
	  std::set<interval_t>::const_iterator ii = ints.begin();
	  while ( ii != ints.end() )
	    {
	      annot_n[ qq->first ]++;
	      annot_s[ qq->first ] += ii->duration_sec();
	      ++ii;
	    }	  
	  ++qq;
	}
      ++rr;
    }
  
  std::map<std::string,int>::const_iterator qq = annot_n.begin();
  while ( qq != annot_n.end() )
    {
      logger << "  " << qq->first;
      if ( sachs.find( qq->first ) != sachs.end() ) logger << " [seed]"; else logger << " [other]";
      logger << " : n = " << qq->second
	     << " , mins = " << annot_s[ qq->first ] / 60.0
	     << " , avg. dur (s) = " << annot_s[ qq->first ] / (double)qq->second;
      
      if ( aligned_permutes.find( qq->first ) != aligned_permutes.end() )
	logger << " [aligned shuffle across channels]";

      if ( fixed.find( qq->first ) != fixed.end() || ( ! shuffle_annots && sachs.find( qq->first ) == sachs.end() ) )
	logger << " [fixed]";

      if ( midpoint_annot.find( qq->first ) != midpoint_annot.end() ) logger << " [midpoint]";

      if ( flanking_sec > 0 && sachs.find( qq->first ) != sachs.end() )
	logger << " [f=" << flanking_sec << "]";
      else if ( flanking_sec_annot.find( qq->first ) != flanking_sec_annot.end() )
	logger << " [f=" << flanking_sec_annot[ qq->first ] << "]";
      
      logger << "\n";
      
      ++qq;
    }
  
 
  // info about aligned permutes
  
  if ( aligned_permutes.size() )
    {
      std::map<std::string,std::set<std::string> >::const_iterator aa = aligned_permutes.begin();
      while ( aa != aligned_permutes.end() )
	{
	  logger << " aligned permute : " << aa->first ;
	  std::set<std::string>::const_iterator bb = aa->second.begin();
	  while ( bb != aa->second.end() )
	    {
	      logger << " " << *bb ; 
	      ++bb;
	    }
	  logger << "\n";
	  ++aa;
	}
    }
  

  if ( constrained_shuffle_dur )
    logger << "  shuffling constrained to +/- " << max_shuffle_sec << "s within each contiguous background interval\n";  
  else
    logger << "  unconstrained shuffling within each contiguous background interval\n";
  
}



void annotate_t::loop()
{

  if ( debug_mode )
    {
      std::cout << "--- observed data ---\n";
      view();
    }
  
  // evaluate the original dataset
  annotate_stats_t s = eval();

  // record
  observed( s );

  // do we need to track the original, observed events?
  if ( constrained_shuffle_dur )
    observed_events = events;
  
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

      // do we need to replace originals?
      if ( constrained_shuffle_dur )
	events = observed_events;
      
      // permute
      shuffle();

      // verbose output?
      if ( debug_mode )
	{
	  std::cout << "--- shuffled data, replicate " << r + 1 << " ---\n";
	  view();
	}
      
      // calc statistics for null data
      annotate_stats_t s = eval();

      // track null distribution
      build_null( s );
    }
  
}



void annotate_t::shuffle()
{

  // constrained shuffle (duration-wise)?
  const uint64_t maxshuffle_constrained_tp = constrained_shuffle_dur ? max_shuffle_sec * globals::tp_1sec : 0 ; 
  
  //
  // shuffle each seed/region independently
  //
  
  // contiguous regions
  std::map<uint64_t,std::map<std::string,std::set<interval_t> > >::const_iterator rr = events.begin();
  while ( rr != events.end() )
    {
      
      //
      // shuffle (w/ wrapping) each seed class independently
      // ( potentially aligning across channels ('align') or doing channels independently)
      //   
      // by default, can shuffle fully across the entire region;  optionally ('shuffle-sec')
      // we can constrain the shuffle to within an X second window (in either direction)
      //

      const uint64_t region_size = seg[ rr->first ] ; 
      
      const uint64_t maxshuffle = constrained_shuffle_dur ?
	( maxshuffle_constrained_tp * 2 > region_size ? region_size : maxshuffle_constrained_tp * 2 )
	: region_size ;
      
      // if permuting blocks of seed together
      std::map<std::string,uint64_t> aligned_shuffle;
      
      // permuting seeds only, or seeds + annots too
      const std::set<std::string> & permset = shuffle_annots ? achs : sachs;
      
      // each seed for this region
      std::set<std::string>::const_iterator ss = permset.begin();
      while ( ss != permset.end() )
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
	  
	  if ( aligned_shuffle.find( *ss ) != aligned_shuffle.end() ) // aligned already 
	    {
	      pp = aligned_shuffle[ *ss ];	      
	    }
	  
	  // *otherwise*, get a new random offset
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
		  if ( iter > 500 )
		    Helper::halt( "cannot find any valid shuffle sets for " + *ss
				  + "\n please sanity-check the number/size of background/event intervals" );
		  
		  // putative shuffle 
		  
		  pp = maxshuffle != 0 ? CRandom::rand( maxshuffle ) : 0 ; 
		  
		  // in constrained shuffle
		  //  REGION 0 .....   1000
		  //  constrined = +/- 50
		  //      ORIG = 200
		  //        range  -->  150 .. 250 
		  //      
		  //  perm max = 2 * 50 = 100
		  //   if <  50 , then just shuffle forward
		  //   if >= 50 , then shuffle backward
		  //    *but* to keep things positive, 0  <= p < 50  -> p = 0 - 50
		  //                                   50 <= p < 100 -> p = REG-SIZE - ( p - 50 ) 
		  //     i.e. we do a full wrap to do the 'backwards' shuffle

		  if ( pp >= maxshuffle_constrained_tp )
		    pp = region_size - ( pp - maxshuffle_constrained_tp ) ; 

		  // check all events
		  bool okay = true;
		  const std::set<interval_t> & original = events[ rr->first ][ *ss ];
		  std::set<interval_t>::const_iterator ii = original.begin();
		  while ( ii != original.end() )
		    {
		      interval_t i = *ii;
		      i.start += pp;
		      i.stop += pp;
		      // check - spans segment break?
		      if ( i.start < region_size && i.stop >= region_size )
			{
			  okay = false;
			  break;
			}
		      ++ii;
		    }
		  
		  // need to try again for this one?
		  if ( ! okay ) continue;

		  // will this be an aligned-group shuffle? if so, we need to check that
		  // all the members are good too before accepting.
		  
		  if ( aligned_permutes.find( *ss ) == aligned_permutes.end() ) break;

		  const std::set<std::string> & apset = aligned_permutes.find( *ss )->second;
		  std::set<std::string>::const_iterator kk = apset.begin();
		  while ( kk != apset.end() )		    
		    {
		      if ( *kk == *ss ) { ++kk; continue; }

		      // test this seed for break-point overlap
		      
		      bool okay1 = true;
		      const std::set<interval_t> & original = events[ rr->first ][ *kk ];
		      std::set<interval_t>::const_iterator ii = original.begin();
		      while ( ii != original.end() )
			{
			  interval_t i = *ii;
			  i.start += pp;
			  i.stop += pp;
			  // check - spans segment break?
			  if ( i.start < region_size && i.stop >= region_size )
			    {
			      okay1 = false;
			      break;
			    }
			  ++ii;
			}

		      // failed on an aligned annot?
		      if ( ! okay1 )
			{
			  okay = false;
			  break;
			}
		      ++kk;
		    }

		  // we good overall?
		  if ( okay ) break;
		  
		  // otherwise, loop back and try again...		  		  
		}	      
	      
	      // save this aligned shuffle?
	      
	      if ( aligned_permutes.find( *ss ) != aligned_permutes.end() )
		{
		  const std::set<std::string> & apset = aligned_permutes.find( *ss )->second;
		  
		  std::set<std::string>::const_iterator kk = apset.begin();
		  while ( kk != apset.end() )
		    {
		      aligned_shuffle[ *kk ] = pp;
		      ++kk;
		    }
		}
	      
	    }
	  
	      
	  //
	  // now we have a valid shuffle value... do the shuffle 
	  //
	
	  // copy over events
	  std::set<interval_t> original = events[ rr->first ][ *ss ];
	  
	  std::set<interval_t> shuffled;
	  
	  std::set<interval_t>::const_iterator ii = original.begin();
	  while ( ii != original.end() )
	    {
	      
	      interval_t i = *ii;
	      
	      i.start += pp;
	      i.stop += pp;
	      
	      // need to wrap?
	      if ( i.start >= region_size )
		{
		  i.start -= region_size ;
		  i.stop -= region_size ;
		}
	      
	      // add this new one to the list
	      shuffled.insert( i );
	      
	      ++ii;
	    }

	  // update
	  events[ rr->first ][ *ss ] = shuffled;
	  
	  // next seed
	  ++ss;
	}
      
      // next region
      ++rr;
    }
  
}


annotate_stats_t annotate_t::eval()
{
  // main task: aggregate overlap stats in 'r'
  annotate_stats_t r;

  // secondary: optionally (only w/ true data) output new
  // seed annotations (if 'make_anew' set)

  int scnt = 0;
  
  // each region
  std::map<uint64_t,std::map<std::string,std::set<interval_t> > >::const_iterator rr = events.begin();
  while ( rr != events.end() )
    {

      // track offset			       
      uint64_t offset = rr->first;
      
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
	  
	  //
	  // track for pile-up
	  //

	  if ( do_pileup )
	    {
	      std::set<interval_t>::const_iterator qq = a.begin();
	      while ( qq != a.end() )
		{
		  puints.insert( named_interval_t( offset , *qq , *aa ) );
		  ++qq;
		}
	    }
	  
	  //
	  // track # of (flattened) annots
	  //

	  // by default, only flatten annotations that actually overlap here
	  // i.e. keep contiguous annots 'as is' ;  this should not really matter,
	  // but it will avoid a possible edge case where annots start and end at
	  // the edge of the region -- when wrapped, these two would be merged
	  // and so this would change the overall number of annots.
	  
	  const bool flatten_mode = false;
	  
	  //std::set<interval_t> flata = flatten( a , flatten_mode );
	  r.ns[ *aa ] += a.size();
	  
	  // ensure s2a mapping is initialized for each A
	  std::set<interval_t>::const_iterator ff = a.begin();
	  while ( ff != a.end() )
	    {	  
	      // if already exists, leave as is; otherwise need to insert an empty set
	      // so that the first key exists for each seed event	      
	      
	      named_interval_t na( offset, *ff , *aa );
	      
	      if ( r.s2a_mappings.find( na ) == r.s2a_mappings.end() )
		{
		  std::set<std::string> dummy;
		  r.s2a_mappings[ na ] = dummy;
		}

	      // next seed annot
	      ++ff;
	    }

	  //
	  // consider all other annots
	  //
	  
	  std::set<std::string>::const_iterator bb = achs.begin();
	  while ( bb != achs.end() )
	    {		  
	      // skip self comparison
	      if ( *aa == *bb ) { ++bb; continue; }
	      
	      // requiring only intra-channel comparisons?
	      if ( only_within_channel && ! same_channel( *aa , *bb ) ) { ++bb; continue; }
	      
	      // does this interval have any annots?
	      if ( rr->second.find( *bb ) == rr->second.end() ) { ++bb; continue; }
	      
	      // get all other annots
	      const std::set<interval_t> & b = rr->second.find( *bb )->second;
	      
	      // calc and record stats on flattened b lists (so that nearest neighbour search
	      // only needs to go to lower-bound (at/after) and possibly one step before)
	      seed_annot_stats( a , *aa , flatten( b , flatten_mode ) , *bb , offset , &r );

	      // next annot
	      ++bb;
	    }
	  
	  // next seed
	  ++aa;
	}

      //
      // seed-seed pileup 	  
      //

      if ( do_pileup )
	{
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
	}


      //
      // next region
      //

      ++rr;
    }
 
   // all done
  return r;
  
}


void annotate_t::output()
{

  //
  // seed-seed group overlap
  //

  if ( do_pileup )
    {
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
    }

  
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
	  if ( var > 0 ) 
	    writer.value( "PROP_Z" , ( prop_obs[ pp->first ] - mean ) / sqrt( var ) ); 
        }
      ++pp;
    }
  writer.unlevel( "SEED" );


  //
  // one-to-many seed/annot overlap
  //
  
  std::map<std::string,std::map<std::string,uint64_t> >::const_iterator saa = s2a_obs.begin();
  while ( saa != s2a_obs.end() )
    {
      writer.level( saa->first , "SEED" );

      const std::map<std::string,uint64_t> & p = saa->second;
      std::map<std::string,uint64_t>::const_iterator pp = p.begin();
      while ( pp != p.end() )
        {
          writer.level( pp->first , "OTHERS" );
	  
	  //writer.value( "N_OBS" , s2a_obs[ saa->first ][ pp->first ]  );
	  writer.value( "N_OBS" , (int)pp->second );
	  
	  if ( nreps )
	    {
	      double mean = s2a_exp[ saa->first ][ pp->first ] / (double)nreps;
	      double var = s2a_expsq[ saa->first ][ pp->first ] / (double)nreps - mean * mean;
	      writer.value( "N_EXP" , mean );
	      //std::cout << " var = " << var << " " << s2a_expsq[ saa->first ][ pp->first ]  << "\n";
	      if ( var > 0 )
		writer.value( "N_Z" , ( (double)s2a_obs[ saa->first ][ pp->first ] - mean ) / sqrt( var ) );
	    }
	  
	  ++pp;
	}
      writer.unlevel( "OTHERS" );
      
      ++saa;
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
	  writer.level( pp->first , "OTHER" );

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
		  if ( var > 0 ) 
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
	      if ( var > 0 ) 
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

	      if ( var > 0 ) 
		writer.value( "D2_Z" , ( sgnd_obs[ sa->first ][ pp->first ] - mean ) / sqrt( var ) );
	    }
	  
	  ++pp;
	}	  
      writer.unlevel( "OTHER" );
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
  if ( seg.find( *u1 ) == seg.end() ) return false;
  
  // seems okay, return offset for bounding segment
  *offset = *u1;
  
  return true;
}


bool annotate_t::segment( const interval_t & i , uint64_t * segoff ) const
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
	  basket.insert( *ii );
	  last.i.stop = last.i.stop > ii->i.stop ? last.i.stop : ii->i.stop;
	}
      
      ++ii;
    }

  // last one
  ++r[ "_O" + Helper::int2str( (int)basket.size() ) ];
  ++r[ Helper::int2str( (int)basket.size() ) + ":" + stringize( basket ) ];

  // print
  // std::cout << "\nbasket size = " << basket.size() << "\n";
  // std::set<named_interval_t>::const_iterator bb =  basket.begin();
  // while ( bb != basket.end() )
  //   {
  //     std::cout << bb->n << "\t"
  // 		<< bb->i.start << "\t"
  // 		<< bb->i.stop << "\n";
      
  //     ++bb;
  //   }

  
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
  
  // first reduce to names, i.e. so B,A --> A,B 
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
				   uint64_t offset , 
				   annotate_stats_t * r )
{
  //  debug_mode = true;
  // if ( debug_mode ) std::cout << "\nseed_annot_stats( "
  //  			      << astr << " n = " << a.size() << " -- "
  //  			      << bstr << " n = " << b.size() << " )\n";
  
  // if no b annots in this segment, then nothing to do
  // stats remain as they are  
  if ( b.size() == 0 ) return;
  
  // is 'b' also a seed?
  const bool bseed = sachs.find( bstr ) != sachs.end();
  
  // consider each seed
  std::set<interval_t>::const_iterator aa = a.begin();
  while ( aa != a.end() )    
    {
      
      double dist = 0;
      bool overlap = false;
      
      // find the first annot not before (at or after) the seed
      std::set<interval_t>::const_iterator bb = b.lower_bound( *aa );

      // track matched B (verbose output only)
      interval_t bmatch = *bb;
      std::string mtype = ".";
      
      // verbose output
      // if ( debug_mode )
      //   std::cout << "a = " << aa->as_string() << "\n";

      // edge cases:
      // no annot at or past seed? : bb == b.end() 
      // no annot before seed?     : bb == b.begin()
      
      
      //              |-----| SEED           LB?  MTYPE
      //    |-----|   |     |                     B2
      //            |-|-|   |                     B1
      //              |--|  |                 Y   A1 
      //              |-----|                 Y   A1
      //              |   |-|                 Y   A1
      //              |   |-|--|              Y   A1
      //              |     |----|            Y   A2
      //              |     |        |---|    Y   A2

      // if ( debug_mode )
      // 	{
      // 	  std::cout << " initial b = ";
      // 	  if ( bb == b.end() ) std::cout << " -END-";
      // 	  else
      // 	    {
      // 	      std::cout << bb->as_string() ;
      // 	      if ( bb == b.begin() ) std::cout << " (begin)";
      // 	    }	  
      // 	  std::cout << "\n";
      // 	}

      
      // does first take overlap seed?
      
      if ( bb != b.end() && bb->overlaps( *aa ) )
	{
	  // we're done, found complete overlap
	  dist = 0;
	  overlap = true;

	  // if ( debug_mode )
	  //   std::cout << "  found overlap, dist = 0 \n";
	  
	}
      else // it must come afterwards
	{
	  
	  uint64_t seed_start = aa->start;
	  uint64_t seed_stop  = aa->stop - 1LLU;
	  
	  // nb. here -1 means that no right distance is defined
	  dist = bb != b.end() ?
	    ( bb->start - seed_stop ) * globals::tp_duration : 
	    -1; 
	  
	  // if ( debug_mode )
	  //   std::cout << " prov dist = " << dist << "\n";
	  
	  // step back, if we can - is there a closer annot /before/ the seed?
	  if ( bb != b.begin() )
	    {
	      --bb;

	      // if ( debug_mode )
	      // 	std::cout << " stepping back, b -> " << bb->as_string() << "\n";
	      
	      // nb - this may overlap seed
	      // i.e. starts before, but ends after seed-start, and so
	      //  was not captured by the lower_bound()
	      
	      if ( bb->stop > aa->start )
		{
		  dist = 0;
		  overlap = true;

		  // if ( debug_mode )
		  //   std::cout << "  back b overlaps a, done\n"; 
		}
	      else
		{
		  
		  // annot stops before start of seed
		  double left_dist = ( aa->start - ( bb->stop - 1LLU ) ) * globals::tp_duration; 
		  
		  // closer? (or equal to)
		  // nb. store as signed here (neg -> before)
		  
		  // nb. use dist < 0 condition to track that there was no
		  // lower_bound (i.e. the final annot occurs before this seed)
		  // in which case, this left annot will be the closest 
		  if ( dist < 0 || left_dist <= dist )
		    dist = - left_dist;

		  // if ( debug_mode )
		  //   {
		  //     std::cout << " back b is before, so dist = " << left_dist << "\n";
		  //     std::cout << " final dist = " << dist << "\n";
		  //   }
		}
	    }
	}
      
      // track: overlap = dist == 0,
      // but use bool overlap to avoid floating-point equality test
      if ( overlap ) r->nsa[ astr ][ bstr ] += 1 ;
      
      // to track proprtion of seeds w/ at least one (non-seed) annot overlap
      if ( overlap && ! bseed )
	{
	  //	  std::cout << " found overlap : " << astr << " " << bstr << " = " << aa->as_string() << " by " << bb->as_string() << "\n";
 	  r->psa[ astr ].insert ( named_interval_t( offset, *aa , astr ) );
	  //	  std::cout << "   size = " << r->psa[ astr ].size() << "\n";
	}
      
      // save original distance
      const double adist_orig = fabs( dist );

      // truncate at window length
      if ( dist > window_sec ) dist = window_sec;
      else if ( dist < -window_sec ) dist = -window_sec;
      
      // for mean distance -- do we meet the window criterion?
      const double adist = fabs( dist );
      
      // do we include complete overlap as "nearest"?
      if ( include_overlap_in_dist || ! overlap )
	{
	  r->adist[ astr ][ bstr ] += adist ; 
	  
	  // track only -1 or +1 for 'before' or 'after'
	  // overlap == 0 here, so add that qualifier
	  if ( ! overlap ) 
	    r->sdist[ astr ][ bstr ] += dist > 0 ? +1 : -1 ;
	  
	  r->ndist[ astr ][ bstr ] += 1; // denom for both the above
	}
      
      // are we tracking hits
      if ( make_anew )
	{
	  
	  // use adist_orig as we truncated the other adist
	  
	  //if ( overlap || adist_orig <= window_sec )
	  if ( overlap ) // only report based on absolute overlap
	    {
	      // only tracking seed-nonseed matches? or all?
	      
	      const bool okay = seed_nonseed ? ! bseed : true ; 
	      
	      if ( okay )
		{
		  named_interval_t named( offset, *aa , astr );
		  hits[ named ]++;
		}
	    }
	}

      // build up 1-to-many seed-annot mappings
      named_interval_t na( offset, *aa , astr ) ;
      if ( ! bseed )
	{
	  if ( overlap )
	    r->s2a_mappings[ na ].insert( bstr );
	  else
	    {
	      // if already exists, leave as is; otherwise need to insert an empty set
	      // so that the first key exists	  
	      if ( r->s2a_mappings.find( na ) == r->s2a_mappings.end() )
		{
		  std::set<std::string> empty;
		  r->s2a_mappings[ na ] = empty;
		}
	    }
	}


      // next seed annot
      ++aa;
    }
  
  
}



std::set<interval_t> annotate_t::flatten( const std::set<interval_t> & x , const bool join_neighbours )
{

  std::set<interval_t> m;
  
  if ( x.size() == 0 ) return m;
  
  interval_t curr = *x.begin();
  
  std::set<interval_t>::const_iterator xx = x.begin();
  while ( xx != x.end() )
    {
      // because of +1 end encoding, contiguous regions will (a,b) ... (b,c)
      // truly overlapping events will be (a,b) .. (c,d) where c < b

      const interval_t & pro = *xx;

      // ORIGINAL: if ( pro.start > curr.stop )
      // NEW: option to not join contiguous events
      if ( join_neighbours ? pro.start > curr.stop : pro.start >= curr.stop )
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


std::set<interval_t> annotate_t::excise( const std::set<interval_t> & y , const std::set<interval_t> & x )
{

  if ( x.size() == 0 || y.size() == 0 ) return y;
  
  // ensure that exclusions are flattened
  std::set<interval_t> fx = flatten( x , true );
  
  std::set<interval_t> z;

  std::set<interval_t>::const_iterator yy = y.begin();
  while ( yy != y.end() )
    {
      const interval_t & interval = *yy;
      
      // find the first annot not before (at or after) this 
      std::set<interval_t>::const_iterator xx = fx.lower_bound( interval );

      // as flattened, slide back one, i.e. to see whether
      // x starts before y but overlaps
      if ( xx != fx.begin() )
	{
	  // slide back one
	  --xx; 
	  
	  // too early? revert to first 
	  if ( xx->stop <= interval.start ) ++xx;
	}
      
      // no overlap?
      if ( xx == fx.end() ) { z.insert( interval ); ++yy; continue; }
      
      // this is after target?
      if ( xx->start >= interval.stop ) { z.insert( interval ); ++yy; continue; }

      
      // if here, means we have at least one region to exclude 
      
      uint64_t curr = interval.start;
      
      while ( 1 )
	{
	  // putative new interval, up until the start of the exclusion
	  // (as we know that bb starts within this spanning interval)
	  if ( curr < xx->start )
	    z.insert( interval_t( curr , xx->start ) );
	  
	  // update curr to after this hole (i.e. +1 = start of new)
	  curr = xx->stop;
	  
	  // but check this is not past end; if so, all done
	  if ( curr >= interval.stop ) { break; }

	  // search for another hole?
	  ++xx;
	  
	  // nothing found? then jump out 
	  if ( xx == fx.end() ) { break;}
	  if ( xx->start >= interval.stop ) { break; }
	  
	}

      // do we have any remaining span we want to add?
      if ( curr < interval.stop )
	z.insert( interval_t( curr , interval.stop ) );
      
      // all done, find next interval to prune
      ++yy;
    }
  
  return z;
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
  std::map<std::string,std::set<named_interval_t> >::const_iterator pp = s.psa.begin();
  while ( pp != s.psa.end() )
    {
      // std::cout << " prop obs " << pp->first << " " << pp->second.size() << " " << s.ns.find(  pp->first )->second << " = "
      // 		<< pp->second.size() / s.ns.find(  pp->first )->second << "\n";
      prop_obs[ pp->first ] = pp->second.size() / s.ns.find(  pp->first )->second;
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

  // 1-to-many seed-to-annots mappings
  s2a_obs = s2a_proc( s.s2a_mappings );

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
      //      std::cout << " TESTING " << pp->first << " " << is_seen << "\n";
      if ( is_seen )
        {
          double val = s.psa.find( pp->first )->second.size() / s.ns.find( pp->first )->second;
	  //  std::cout << " val (exp) = " << val << " " << s.psa.find( pp->first )->second.size() << " " << s.ns.find( pp->first )->second << "\n";
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
	      
	      // nb. we've calculated mean pre/prior, but here the test is 2-sided, so take abs(x)
	      if ( fabs(s) > fabs( sgnd_obs[ sa->first ][ pp->first ] ) ) ++sgnd_pv[ sa->first ][ pp->first ];
	    }
          ++pp;
        }
      ++sa;
    }
  
  //
  // 1-to-many seed-to-annots mappings
  //

  // summarize permuted values
  std::map<std::string,std::map<std::string,uint64_t> > s2a = s2a_proc( s.s2a_mappings );
  
  std::map<std::string,std::map<std::string,uint64_t> >::const_iterator qq = s2a_obs.begin();
  while ( qq != s2a_obs.end() )
    {
      const std::map<std::string,uint64_t> & k = qq->second;
      std::map<std::string,uint64_t>::const_iterator kk = k.begin();
      while ( kk != k.end() )
	{
	  // aggregate
	  uint64_t perm = s2a[ qq->first ][ kk->first ];
	  s2a_exp[  qq->first ][ kk->first ] += perm;
	  s2a_expsq[  qq->first ][ kk->first ] += perm * perm;
	  ++kk;
	}
      ++qq;
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

      std::map<uint64_t,std::map<std::string,std::set<interval_t> > >::const_iterator rr = events.begin();
      while ( rr != events.end() )
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
	      named_interval_t named( offset, *ii , *ss );
	      
	      // requisite # of hits
	      const bool write_this = out_include ?
		hits[ named ] >= mcount : // include
		hits[ named ] <  mcount ; // exclude
	      
	      // n.b. - we drop any instance ID information given we've
	      //  a) flattened, and b) perhaps pooled channels, this
	      //  is probably the right thing to do...
	      
	      if ( write_this )
		{
		  // get original, mapped/unmanipulated version
		  // i.e. will not have +/- 'f' or mid-point reduction
		  // and will be in normal time space (not mapped relative
		  // to start of this background interval)
		  
		  if ( unmanipulated.find( named ) == unmanipulated.end() )
		    Helper::halt( "internal problem tracking named_interval_t in making a new annotatipn" );
		  
		  interval_t mapped = unmanipulated[ named ];

		  a->add( "." , mapped , chname );
		  ++acnt;
		}
	      ++tcnt;
	      
	      // next interval
	      ++ii;
	    }
	  ++rr;
	}
      
      logger << "   - wrote " << acnt << " (of " << tcnt << ") seed events, based on ";
      if ( ! out_include ) logger << "not ";
      logger << "matching " << mcount << " or more other annots, f=" << flanking_sec << "\n";
      
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


std::map<std::string,std::map<std::string,uint64_t> > annotate_t::s2a_proc( const std::map<named_interval_t,std::set<std::string> > & s )
{

  std::map<std::string,std::map<std::string,uint64_t> > r;
  
  //  std::cout << " s2a_proc S.size() " << s.size() << "\n";
  
  int cnt = 0;
  
  std::map<named_interval_t,std::set<std::string> >::const_iterator ss = s.begin();
  while ( ss != s.end() )
    {
      // seed name
      const std::string & seed = ss->first.n;

      //      std::cout << " seed = " << ++cnt << " " << ss->first.i.as_string() << "\n";

      // stringize annots -- default (none --> .)
      std::string mapped = ".";
      
      const std::set<std::string> & as = ss->second;

      if ( as.size() != 0 )
	{
	  std::stringstream sstr;
	  std::set<std::string>::const_iterator aa = as.begin();
	  while ( aa != as.end() )
	    {
	      if ( aa != as.begin() )
		sstr << ",";
	      sstr << *aa;
	      ++aa;
	    }
	  mapped = sstr.str();
	}

      // count 
      r[ seed ][ mapped ]++;
      
      // next seed event
      ++ss;
    }

    
  return r;
    
}

void annotate_t::view()
{
  
  std::map<uint64_t,std::map<std::string,std::set<interval_t> > >::const_iterator rr = events.begin();
  while ( rr != events.end() )
    {
      const std::map<std::string,std::set<interval_t> > & annots = rr->second;
      std::map<std::string,std::set<interval_t> >::const_iterator qq = annots.begin();
      while ( qq != annots.end() )
        {
          const std::set<interval_t> & ints = qq->second;
          std::set<interval_t>::const_iterator ii = ints.begin();
          while ( ii != ints.end() )
            {
	      std::cout << "region = " << rr->first << "\t"
			<< "annot = " << qq->first << "\t"
			<< "event = " << ii->as_string() << "\n";
              ++ii;
            }
          ++qq;
        }
      ++rr;
    }
  std::cout << "\n";
}
