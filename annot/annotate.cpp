
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

std::map<std::string,int> annotate_t::global_alignment_group;
std::map<int,std::set<std::string> > annotate_t::group_global_alignment;


// initiate from a single attached EDF/timeline ( 'OVERLAP' command )
annotate_t::annotate_t( edf_t & edf1 , param_t & param )
{
  edf = &edf1;
  single_indiv_mode = true;
  set_options( param );
  prep();
  if ( event_perm ) init_event_permutation();
  annotate_stats_t s = loop();
  output(s);
}

  
annotate_t::annotate_t( param_t & param )
{

  // command-line, multi-sample invocation
  
  // create a 'super' individual, where we simply concatenate segments
  // across individuals; we pull the same annotations from all sets 
  
  // this means we cannot write out annotations as matched/unmatched 
  single_indiv_mode = false;

  // in event perm mode, need to build some structures to ensure only within-person
  // event permutations

  event_perm = param.has( "event-perm" );

  if ( event_perm )
    {
      if ( ! param.empty( "event-perm") )
	event_neighbourhood_sec = param.requires_dbl( "event-perm" );
      else
	event_neighbourhood_sec = 5; 
    }
  
  // nb. we only allow this if 'bg' mode is specified, or if using block permutation -- this will
  // implicitly ensure that annotations are only shuffled within
  // individuals

  if ( ! ( param.has( "bg" ) || event_perm ) )
    Helper::halt( "bg specification, or event permutation, is required in multi-sample mode" );

  
  // expect a file: indiv -- annot
  
  std::string alist = Helper::expand( param.requires( "a-list" ) );

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
      annots[ tok[0] ].insert( Helper::expand( tok[1] ) );
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
	      if ( tok.size() != 6 ) Helper::halt( "expecting standard 6-field annotations:" + line +
						   "\nuse WRITE-ANNOTS to generate" );

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

      
      // for block-based permutation, special code for individuals
      // to allow for intra-indiv perms only

      if ( event_perm )
	OUT1 << "indiv_int_mrkr" << "\t"
	     << indiv << "\t"
	     << "." << "\t"
	     << offset << "\t"
	     << offset + ind2dur[ indiv ] << "\t"
	     << "." << "\n";
            
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
  if ( event_perm ) init_event_permutation();
  annotate_stats_t s = loop();
  output(s);
  
}

void annotate_t::set_options( param_t & param )
{

  // number of permutations
  nreps = param.has( "nreps" ) ? param.requires_int( "nreps" ) : 1000 ; 

  // verbose/debug mode
  debug_mode = param.has( "verbose" ) || param.has( "debug" );

  // add shuffled annots (only in INDIV mode)
  add_shuffled_annots = param.has( "add-shuffled-annots" );
  if ( add_shuffled_annots && ! single_indiv_mode )
    Helper::halt( "cannot add-shuffled-annots in multi-individual mode" );
  shuffled_annots_names = param.strset( "add-shuffled-annots" );
  shuffled_annots_tag = param.has( "shuffled-annots-tag" ) ? param.value( "shuffled-annots-tag" ) : "s_" ; 
  
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

  // reduce annotations to a relative position (rp_) specified meta
  // field (this will be 0..1 as a relative position w/in the interval)
  // e.g. SO peaks, SP midpoints, etc
  use_rps = false;
  if ( param.has( "rp" ) )
    {
      // expecting rp=<annot>|<tag>,<annot>|<tag>
      //   to get rp_<tag> for <annot>
      std::vector<std::string> tok = param.strvector( "rp" );
      for (int i=0; i<tok.size(); i++)
	{
	  std::vector<std::string> tok2 = Helper::parse( tok[i] , "|" );
	  if ( tok2.size() != 2 ) Helper::halt( "expecting rp=<annot>|<rp-tag>" );
	  rp_annot[ tok2[0] ] = tok2[1] ;
	}
      use_rps = rp_annot.size() != 0;
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
	  // assume annot~f,annot~f format, or | delim
	  for (int t=0; t<tok.size(); t++)
	    {
	      std::vector<std::string> tok2 = Helper::parse( tok[t] , "~|" );
	      if ( tok2.size() != 2 ) Helper::halt( "expecting annot|second format for 'f'" );
	      if ( ! Helper::str2dbl( tok2[1] , &fval ) ) Helper::halt( "expecting annot|second format for 'f'" );
	      flanking_sec_annot[ tok2[0] ] = fval ; 
	    }
	}
    }
  

  // specify flanking offset windows for 'nosa' overlap?
  flanking_overlap_intervals.clear();
  flanking_overlap_mx = 0LLU;
  n_flanking_offsets = 0;
  
  // nb. we only extend around seeds
  if ( param.has( "offset-seed" ) ) flanking_overlap_seeds = param.strset( "offset-seed" );
  if ( param.has( "offset-other" ) ) flanking_overlap_others = param.strset( "offset-other" );
  
  if ( param.has( "offset" ) )
    {
      std::vector<double> s = param.dblvector("offset");
      if ( s.size() != 3 ) Helper::halt( "expecting offset=size,inc,max in seconds" );
      if ( s[0] < 0 || s[1] <= 0 || s[2] <= 0 ) Helper::halt( "offset values must be positive" );
          
      flanking_overlap_intervals.clear();
      uint64_t dur = s[0] * globals::tp_1sec;
      for (double ss = 0 ; ss <= s[2] ; ss += s[1] )
	{
	  uint64_t start = ss * globals::tp_1sec; 
	  flanking_overlap_intervals.push_back( interval_t( start , start + dur ) );
	  flanking_overlap_mid.push_back( ss + s[0]/2.0 );
	  flanking_overlap_desc.push_back( Helper::dbl2str( ss + s[0]/2.0 ) + "+/-" + Helper::dbl2str( s[0]/2.0 ) ) ;
	  flanking_overlap_mx = start + dur; // sets to max: specifies search range around each seed
	}
      
      n_flanking_offsets = flanking_overlap_intervals.size(); 
    }

  if ( n_flanking_offsets )
    logger << "  evaluating " << n_flanking_offsets << " flanking intervals for overlap (offset arg)\n";
  
  
  // distance to neighbour stats
  window_sec = param.has( "w" ) ? param.requires_dbl( "w" ) : 10 ;
  
  // distance to marker stats: default 12 hours = 12*60^2 effectively means no limit
  marker_window_sec = param.has( "mw" ) ? param.requires_dbl( "mw" ) : 12 * 60 * 60;  ;

  // is D2 based on 'before' / 'after' ratio only, or the actual distance ( by default, have the counts +1/-1)
  d2_signed = param.has( "d2-quant" ) ? ( ! param.yesno( "d2-quant" ) ) : true ; 
  
  include_overlap_in_dist = ! param.has( "dist-excludes-overlapping" );
  
  overlap_th = param.has( "overlap" ) ? param.requires_dbl( "overlap" ) : 0 ;
  
  // contrasts 

  contrasts.clear();
  if ( param.has( "contrasts" ) )
    {
      // contrasts=seed1|annot1-seed2|annot2,
      std::vector<std::string> tok = param.strvector( "contrasts" );
      for (int i=0; i<tok.size(); i++)
	{
	  // parse by | and -
	  std::vector<std::string> tok2 = Helper::parse( tok[i] , "|-" );
	  if ( tok2.size() != 4 ) Helper::halt( "expecting contrasts=seed|annot-seed|annot" );
	  // a1 b1 a2 b2 
	  contrasts.push_back( annot_contrast_t( tok2[0], tok2[1], tok2[2], tok2[3] ) );  
	}
    }
  
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
  // non-segment markers?
  //

  mannots.clear();
  has_markers = false;
  if ( param.has( "marker" ) )
    {
      if ( ! param.has( "event-perm" ) )
	Helper::halt( "marker requires event-perm" );
      has_markers = true;
      // root_match() allows expanding, ie.  hyp_clock_*
      mannots = annotate_t::root_match( param.strset( "marker" ) , edf->timeline.annotations.names() );
    }
  
  //
  // do seed-seed pileup?
  //
  
  do_pileup = param.has( "pileup" ) ? param.yesno( "pileup" ) : false ; 

  //
  // include seed-seed pairwise comparisons in SEED-ANNOT pairs? default = no
  //

  do_seed_seed = param.has( "seed-seed" ) ? param.yesno( "seed-seed" ) : false; 
  
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
  // block permutation across intervals (instead of circular shuffle w/ in intervals)
  //

  event_perm = param.has( "event-perm" );

  if ( event_perm )
    {
      if ( ! param.empty( "event-perm" ) )
	event_neighbourhood_sec = param.requires_dbl( "event-perm" );
      else
	event_neighbourhood_sec = 5; 
    }

    
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
  if ( use_rps ) logger <<"  reducing some annotations to rp-tags\n";
  if ( flanking_sec ) logger << "  adding f=" << flanking_sec << " seconds to each annotation\n";
  if ( window_sec ) logger << "  truncating distance at w=" << window_sec << " seconds for nearest neighbours\n";
  if ( marker_window_sec ) logger << "  truncating seed-marker distance at mw=" << marker_window_sec << " seconds\n";

  logger << "  " << ( include_overlap_in_dist ? "" : "not " )
	 << "including overlapping events in nearest-neighbor distances\n";
  
  if ( ordered_groups ) logger << "  ordered=T, so preserving order of seed-seed overlap groups (A,B != B,A)\n";
  else logger << "  ordered=F, so pooling seed-seed permutations, i.e. A,B == B,A (default)\n"; 
  
  // annotations    
  if ( ! param.has( "seed" ) )
    Helper::halt( "require seed argument" );
  
  // requires 1+ seed: look at enrichment of ALL combinations of seeds
  // allow wildcards:
  //sseeds = param.strset( "seed" );
  sseeds = annotate_t::root_match( param.strset( "seed" ) , edf->timeline.annotations.names() );

  
  // from each seed, look at all enrichment w/ all other annots
  //  non-seed annotations are not permuted
  if ( param.has( "other" ) )
    {
      //sannots = param.strset( "other" );
      sannots = annotate_t::root_match( param.strset( "other" ) , edf->timeline.annotations.names() );

    }
  
  // background (i.e. defines the space; only select/permute within contiguous blocks of these regions)
  if ( param.has( "bg" ) ) 
    {
      ///sbgs = param.strset( "bg" );
      sbgs = annotate_t::root_match( param.strset( "bg" ) , edf->timeline.annotations.names() );      
    }
  
  // edge for background -- i.e. if elements could not be placed within X seconds of background segment edge,
  // then denote here
  if ( param.has( "edges" ) )
    edge_sec = param.requires_dbl( "edges" );
  else
    edge_sec = 0;
  
  // exclusionary background (xbg) -- i.e. xbg is the converse of bg, and thus
  // specifies gaps rather than allowed intervals
  if ( param.has( "xbg" ) )
    {
      //sxbgs = param.strset( "xbg" );
      sxbgs = annotate_t::root_match( param.strset( "xbg" ) , edf->timeline.annotations.names() );
    }
  
  if ( param.has( "xbg" ) && ! param.has( "bg" ) )
    Helper::halt( "xbg requires bg to be explicitly specified" );

  // Note -- these are redundant now... can use MAKE-ANNOTS
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
  // other markers; initially put all these in a single bin
  // starting (0); in the multi_indiv case, these will be
  // split out into different bins based on indiv (in init_event_perm())
  //
  
  markers.clear();
  
  aa = mannots.begin();
  while ( aa != mannots.end() )
    {
      annot_t * a = edf->timeline.annotations.find( *aa );
      if ( a == NULL ) logger << "  ** warning, could not find " << *aa << "\n";
      else
	{
	  annot_map_t::const_iterator ii = a->interval_events.begin();
	  while ( ii != a->interval_events.end() )
	    {
	      instance_idx_t instance_idx = ii->first;
	      markers[ 0 ][ *aa ].insert( instance_idx.interval );
	      ++ii;
	    }

	  // flatten, joining adjacent regions
	  if ( markers[0][ *aa ].size() != 0 ) 
	    markers[0][ *aa ] = flatten( markers[0][ *aa ] , true );
	}
      
      ++aa;
    }

  
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

  if ( markers[0].size() != 0 )
    logger << "  " << markers[ 0 ].size() << " marker annotations also included\n";

  
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

      // 
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
	  std::cout << " background discontinuity tp = " << *ff << "\t@sec = " << *ff * globals::tp_duration << "\n";
	  ++ff;
	}
    }

  //
  // track orig/final counts
  //

  std::map<std::string,int> orig_counts;
  
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
  int off_target = 0;
  
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

	  // need to add channel-specific version for 'rp'?
	  if ( ( ! pool ) && rp_annot.find( instance_idx.parent->name ) != rp_annot.end() )
	    {
	      rp_annot[ instance_idx.parent->name + "_" + instance_idx.ch_str ]
		= rp_annot[ instance_idx.parent->name ];
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
	  

	  // count as an original

	  orig_counts[ aid ]++;
	  
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
	  // set to rp-tag time-point?
	  //

	  if ( use_rps && rp_annot.find( aid ) != rp_annot.end() )
	    {

	      // if this is already a 0-dur tp, then skip (but print warning message)
	      uint64_t dur = interval.stop - interval.start;
	      
	      if ( dur == 0LLU )
		{
		  logger << "  warning - requested a rp-tag recoding of a 0-duration time-stamp, skipping...\n";
		}
	      else
		{
		  const instance_t * instance = ii->second;
		  
		  // expecting special encoding: rp_<tag>
		  const std::string & rp_tag = "rp_" + rp_annot[ aid ];
		  
		  avar_t * rp = instance->find( rp_tag );
		  
		  if ( rp == NULL )
		    Helper::halt( "annotation " + aid + " found that does not have required rp-tag " + rp_tag );
		  
		  double rpval = rp->double_value();
		  
		  // can be negative, i.e. implies before
		  
		  // rp      -1      0   0.5   1       2
		  //          .      |---------|       .  
		  
		  if ( rpval >= 0 )
		    {
		      uint64_t p = interval.start + dur * rpval;
		      interval.start = interval.stop = p;
		    }
		  else
		    {
		      uint64_t x = dur * -rpval;
		      if ( interval.start > x )
			interval.start -= x;
		      else
			interval.start = 0LLU;
		      interval.stop = interval.start ; 
		    }		
		  
		}
	    }
	  
	  //
	  // add to the map?
	  //
	  
	  uint64_t offset;	  
	  
	  bool okay = segment( interval , &offset );
	  
	  //std::cout << " ii->f " << aid << " " << okay << "\n";
	  
	  if ( ! okay ) { ++ii; off_target++; continue; } 

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

	  // for event-based perm (unlike circular perm), retain the original (absolute)
	  // tp-encoding i.e. rather than reduced to 0..N for each segment
	  // we do it this way (lazy) just so we don't have to change the right/left edge expansion code
	  // above
	  
	  if ( event_perm )
	    {
	      interval.start += offset;
	      interval.stop += offset;
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

  logger << "  " << off_target << " events fell outside of the background and were rejected\n";
  
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
	     << " of " << orig_counts[ qq->first ] 
	     << " , mins = " << annot_s[ qq->first ] / 60.0
	     << " , avg. dur (s) = " << annot_s[ qq->first ] / (double)qq->second;
      
      if ( aligned_permutes.find( qq->first ) != aligned_permutes.end() )
	logger << " [aligned shuffle across channels]";

      if ( fixed.find( qq->first ) != fixed.end() || ( ! shuffle_annots && sachs.find( qq->first ) == sachs.end() ) )
	logger << " [fixed]";

      if ( use_rps && rp_annot.find( qq->first ) != rp_annot.end() ) logger << " [rp=" << rp_annot[qq->first ] << "]";
      
      if ( midpoint || midpoint_annot.find( qq->first ) != midpoint_annot.end() ) logger << " [midpoint]";

      if ( flanking_sec > 0 && sachs.find( qq->first ) != sachs.end() )
	logger << " [f=" << flanking_sec << "]";
      else if ( flanking_sec_annot.find( qq->first ) != flanking_sec_annot.end() )
	logger << " [f=" << flanking_sec_annot[ qq->first ] << "]";
      
      logger << "\n";
      
      ++qq;
    }

  if ( mannots.size() )
    {
      std::set<std::string>::const_iterator mm = mannots.begin();
      while ( mm != mannots.end() )
	{
	  if ( markers[ 0 ].find( *mm ) != markers[ 0 ].end() )
	    {
	      const std::set<interval_t> & mrk = markers[ 0 ][ *mm ];

	      double secs = 0;
	      std::set<interval_t>::const_iterator qq = mrk.begin();
	      while ( qq != mrk.end() )
		{
		  secs += qq->duration_sec();
		  ++qq;
		}
	      
	      logger << "  " << *mm
		     << " [marker]";
	      logger << " : n = " << mrk.size() 
		     << " , mins = " << secs / 60.0 
		     << " , avg. dur (s) = " << secs / (double)mrk.size() << "\n";	      
	      
	    }
	  ++mm;
	}
    }
  
 
  // info about aligned permutes
  // and make 'global_alignment_group' object (used in pinstance_t/event-perm) 
  
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
  
  // construct global_alignment_group & group_alignment_alignment

  //  seed-name --> map 
  const std::set<std::string> & permset = shuffle_annots ? achs : sachs;
  global_alignment_group.clear();
  
  std::set<std::string>::const_iterator ppp = permset.begin();
  while ( ppp != permset.end() )
    {
      // already done?
      if ( global_alignment_group.find( *ppp ) != global_alignment_group.end() )
	{
	  ++ppp;
	  continue;
	}
      
      // is this in an alignment set? add as is (unique group)? 
      std::map<std::string,std::set<std::string> >::const_iterator aa = aligned_permutes.find( *ppp );
      
      // if not, add as is (unique group)
      if ( aa == aligned_permutes.end() )
	{
	  int gnum = global_alignment_group.size();
	  global_alignment_group[ *ppp ] = gnum;
	}
      else
	{
	  // get all alignments;
	  const std::set<std::string> & g = aa->second;
	  std::set<std::string>::const_iterator gg = g.begin();
	  while ( gg != g.end() )
	    {
	      // has this been assigned to a group already?
	      if ( global_alignment_group.find( *gg ) != global_alignment_group.end() )
		{
		  global_alignment_group[ *ppp ] = global_alignment_group[ *gg ];
		  break;
		}
	      ++gg;
	    }
	  
	  // if still not added, add as new group
	  if ( global_alignment_group.find( *ppp ) == global_alignment_group.end() )
	    {
	      int gnum = global_alignment_group.size();
	      global_alignment_group[ *ppp ] = gnum;
	    }	      
	}
      ++ppp;
    }

  
  // make reverse lookup
  group_global_alignment.clear();
  std::map<std::string,int>::const_iterator mm = global_alignment_group.begin();
  while ( mm != global_alignment_group.end() )
    {
      group_global_alignment[ mm->second ].insert( mm->first );      
      //std::cout << " alignment group " << mm->first << "\t" << mm->second << "\n";
      ++mm;
    }

  

  if ( event_perm )
    logger << "  event-based permutation with " << event_neighbourhood_sec << "s spans\n";
  else if ( constrained_shuffle_dur )
    logger << "  shuffling constrained to +/- " << max_shuffle_sec << "s within each contiguous background interval\n";  
  else
    logger << "  unconstrained shuffling within each contiguous background interval\n";
  
}



annotate_stats_t annotate_t::loop()
{

  if ( debug_mode )
    {
      std::cout << "--- observed data ---\n";
      view();
    }
  
  // evaluate the original dataset
  annotate_stats_t s = eval();

  // any contrasts
  add_contrasts( &s );
  
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
      if ( event_perm )
	event_permutation();
      else
	shuffle();

      // save new annots? (first perm only)
      if ( add_shuffled_annots && r == 0 )
	add_permuted_annots();
      
      // verbose output?
      if ( debug_mode )
	{
	  std::cout << "--- shuffled data, replicate " << r + 1 << " ---\n";
	  view();
	}
      
      // calc statistics for null data
      annotate_stats_t p = eval();

      // any contrasts
      add_contrasts( &p );

      // track null distribution
      build_null( p );
    }

  // return observed stats for output keys
  return s;
}



void annotate_t::add_permuted_annots()
{

  // std::set<std::string> shuffled_annots_names;
  // std::string shuffled_annots_tag;
  
  std::map<uint64_t,std::map<std::string,std::set<interval_t> > >::const_iterator rr = events.begin();
  while ( rr != events.end() )
    {
      // each annot
      
      const uint64_t offset = rr->first;
      
      std::set<std::string>::const_iterator aa = shuffled_annots_names.begin();
      while ( aa !=  shuffled_annots_names.end() )
	{
	  // valid annot
	  if ( rr->second.find( *aa ) != rr->second.end() )
	    {
	      
	      annot_t * a = edf->timeline.annotations.add( shuffled_annots_tag + *aa );
	      
	      const std::set<interval_t> & ints = rr->second.find( *aa )->second;
	      
	      logger << "  adding shuffled/permutation annotation class " << (shuffled_annots_tag + *aa)  << " (" << ints.size() << " events)\n";

	      std::set<interval_t>::const_iterator ii = ints.begin();
	      while ( ii != ints.end() )
		{

		  // nb. event-perm mode stores exact versions, 
		  //     otherwise, we need to add the offset back in
		  if ( event_perm )
		    {		      
		      a->add( "." , *ii , "." );
		    }
		  else // this will likely be the default anyway
		    {
		      interval_t pint( ii->start + offset , ii->stop + offset );
                      a->add( "." , pint , "." );		      
		    }

		  // std::cout << "region = " << rr->first << "\t"
		  // 	    << "annot = " << *aa << "\t"
		  // 	    << "interval = " << ii->as_string() << "\t"
		  // 	    << "dur = " << ii->duration_sec() << "\t"
		  // 	    << ii->start << "\n";
		  
		  ++ii;
		}
	    }
	  ++aa;
	}
      ++rr;
    }

  //  std::cout << "\n";
  
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


void annotate_t::init_event_permutation()
{

  //  view();
  
  //
  // first handle the scenario where we have >1 individual
  // i.e. indiv_int_mrkr annotation will be present in merged.annot
  //
  
  annot_t * aindivs = edf->timeline.annotations.find( "indiv_int_mrkr" );

  multi_indiv = aindivs != NULL; 

  if ( multi_indiv )
    {
      indiv_segs.clear();
      
      annot_map_t::const_iterator aa = aindivs->interval_events.begin();
      while ( aa != aindivs->interval_events.end() )
        {
          const instance_idx_t & instance_idx = aa->first;
          indiv_segs.push_back( instance_idx.interval );
	  //std::cout << " pushing indiv_segs " << instance_idx.interval.as_string() << "\n";
          ++aa;
        }
    }

  
  // neighbourhood size around each seed
  uint64_t neighbourhood_tp = event_neighbourhood_sec * globals::tp_1sec;
  
  // if we need to track which indivs these segments come from (multi_indiv==T)
  std::vector<interval_t>::const_iterator ii = indiv_segs.begin();
  const int nindiv = multi_indiv ? indiv_segs.size() : 1 ;   
  int indiv = 0;
  
  if ( multi_indiv )
    logger << "  running in multi-individual mode, for N=" << nindiv << "\n";
  else
    logger << "  running in single-individual mode\n";

  // input:
  //  std::map<uint64_t,uint64_t> seg      : all segments for all indivs (start/length)
  //  std::vector<interval_t> indiv_segs   : for each indiv, their total span     
  //  
  
  // outputs:
  //  
  // indiv-> start of each segment for that indiv  
  // std::map<int,std::set<uint64_t> > indiv2segs
  //  
  // by indiv, under aligned permutations, identify local friends (will be permuted together)
  // std::map<int,std::map<pinstance_t,std::set<pinstance_t> > > event2friends;

  // by indiv, for each event, the implied neighbourhood span (for permuting)
  // std:map<int, std::map<pinstance_t,interval_t> > event2neighbourhood;;

  //
  // pool all events that will be permuted (seeds, not fixed, etc)
  //

  std::set<pinstance_t> pevents;
  
  
  //
  // consider all segments (starst/offset -> segment size)
  //
  
  std::map<uint64_t,uint64_t>::const_iterator ss = seg.begin();
  
  while ( ss != seg.end() )
    {
      
      uint64_t segstart = ss->first;
      uint64_t seglen = ss->second;

      //      std::cout <<" looking at segstart " << segstart << "\t" << seglen << "\n";
      
      // which 'indiv'? (else fixed to '0')
      //      std::cout << "is multi indiv " << multi_indiv << "\n";
      
      if ( multi_indiv )
	{
	  while (1)
	    {
	      // current indiv okay
	      if ( segstart >= ii->start && segstart < ii->stop )
		break;
	      // advance to next
	      ++ii; ++indiv;
	      if ( indiv > nindiv )
		Helper::halt( "mismatch between indivs and blocks in annotate_t::init_block_permutation()" );
	    }
	}
      
      //      std::cout << " indiv = " << indiv << "\n";
      
      // track segments by indiv
      indiv2segs[ indiv ].insert( segstart );
      seg2indiv[ segstart ] = indiv;
      
      // get actual events, to populate pevents[]
      // std::cout << " events.size() " << events.size() << "\t" << segstart << "\n";
      
      // std::map<uint64_t,std::map<std::string,std::set<interval_t> > >::const_iterator jj1 = events.begin();
      // while ( jj1 != events.end() ) { std::cout << "  events key " << jj1->first << "\n"; ++jj1; } 

      // might not be any events in this segment? -> skip to the next segment
      std::map<uint64_t,std::map<std::string,std::set<interval_t> > >::const_iterator jj = events.find( segstart );      
      if ( jj == events.end() )
	{	  
	  ++ss;
	  continue;
	}
        
      const std::map<std::string,std::set<interval_t> > & annot2event = jj->second;

      // now consider all events that will be permuetd only
      const std::set<std::string> & permset = shuffle_annots ? achs : sachs;
      
      std::set<std::string>::const_iterator pp = permset.begin();
      while ( pp != permset.end() )
	{
	  //	  std::cout << " considering permset " << *pp << "\n";
	  
	  // skip?
          if ( fixed.find( *pp ) != fixed.end() )
            {
              ++pp;
              continue;
            }
	  
	  // pull the events
	  
	  std::map<std::string,std::set<interval_t> >::const_iterator aa = annot2event.find( *pp );

	  if ( aa == annot2event.end() )
	    {
	      ++pp;
	      continue;
	    }

	  // iterate over events
	  const std::set<interval_t> & evts = aa->second;
	  std::set<interval_t>::const_iterator ee = evts.begin();
	  while ( ee != evts.end() )
	    {
	      pevents.insert( pinstance_t( aa->first ,
					   *ee ,					  
					   segstart ) ); 	      
	      ++ee;
	    }
	  
	  // next annotation
	  ++pp;
	}


      //
      // Populate fixed_events structure too
      //

      std::map<std::string,std::set<interval_t> >::const_iterator ff = annot2event.begin();
      while ( ff != annot2event.end() )
	{

	  // skip if not fixed / not permuting

	  bool is_fixed = permset.find( ff->first ) == permset.end();
	  if ( ( ! is_fixed ) && fixed.find( ff->first ) != fixed.end() ) is_fixed = true;

	  if ( ! is_fixed )
	    {
	      ++ff;
              continue;
            }

	  // iterate over events                                                                                                                       
          const std::set<interval_t> & evts = ff->second;
          std::set<interval_t>::const_iterator ee = evts.begin();
          while ( ee != evts.end() )
            {	      
	      fixed_events[ segstart ][ ff->first ].insert( *ee );
	      ++ee;
	    }
	  
	  ++ff;
	}
      
      //      std::cout << "\n\nskipping to next segment....\n";
      
      // next segment
      ++ss;
    }

  //  std::cout << " finished step 1\n";

  //
  // Align any markers to create segments == indivs
  //

  // initially, all markers are in a single bin, starting 0
  // which will always be the start of the first/only indiv
  // in the non-multi_indiv case

  // otherwise, we should create some segments that map to the indiv
  // i.e. to make sure that seeds are not compared to other people's
  // markers
  
  if ( multi_indiv )
    {
      
      interval_map_t pooled_markers = markers;
      markers.clear();
      
      std::set<std::string>::const_iterator mm = mannots.begin();
      while ( mm != mannots.end() )
	{
	  
	  // no markers exist?
	  if ( pooled_markers[ 0 ].find( *mm ) == pooled_markers[ 0 ].end() ) { ++mm; continue; }
	  
	  // get single marker pool
	  const std::set<interval_t> & mrk = pooled_markers[ 0 ][ *mm ];
	  
	  // deal out into indiv pools (if here, will always have at least two people)
	  int indiv = 0;
	  uint64_t indiv_curr = indiv_segs[0].start;
	  uint64_t indiv_next = indiv_segs[1].start;
	  
	  std::set<interval_t>::const_iterator qq = mrk.begin();
	  while ( qq != mrk.end() )
	    {
	      
	      while ( 1 )
		{
		  // is before the next person (or this is the last person)
		  // and so always assign to the last person in that case
		  if ( qq->start < indiv_next || indiv == nindiv - 1 ) break;
		  
		  // next indiv?
		  if ( indiv < nindiv - 1 ) 
		    {
		      ++indiv;
		      
		      indiv_curr = indiv_segs[ indiv ].start;
		      
		      // only need the next marker if there is a person afterwards
		      if ( indiv < nindiv - 1 )
			indiv_next = indiv_segs[ indiv + 1 ].start;
		    }
		}
	      
	      // now 'indiv_curr' should point to segment start for this person
	      // add to the new list
	      markers[ indiv_curr ][ *mm ].insert( *qq );

	      // next marker event
	      ++qq;
	    }
	  
	  // next marker class
	  ++mm;
	}

      logger << "  split markers up into " << markers.size() << " regions/individuals\n";
    }
  
  
  //  std::cout << "\n\n-------------- step2\n\n";
  
  //
  // Now we've populated pevents with all permutable events, define neighbourhoods
  // (w.r.t. other permutable events only)
  //
  
  // under aligned permutations, identify local friends (will be permuted together)
  // std::map<instance_t,std::set<instance_t> > event2friends;
  event2friends.clear();
  
  // for each event, the implied neighbourhood span (for permuting)
  // std::map<instance_t,interval_t> event2neighbourhood;;
  event2neighbourhood.clear();
  
  std::set<pinstance_t>::const_iterator ee = pevents.begin();
  while ( ee != pevents.end() )
    {
      
      // 'index' is ee

      // search for neighbours, using ff
      std::set<pinstance_t>::const_iterator ff = ee;

      //      std::cout << " searching " << ee->name << " " << ee->interval << "\n";
	
      // make friends
      std::set<pinstance_t> friends;
      
      // back
      while ( 1 ) 
	{
	  
	  if ( ff == pevents.begin() ) break;	  
	  
	  // move back
	  --ff;

	  // different segment/indiv?
	  if ( ff->seg != ee->seg ) break;
	  
	  // too far away?
	  if ( ee->interval.start - ff->interval.stop > neighbourhood_tp ) break;
	  
	  // else add, if aligned
	  if ( aligned( ee->name , ff->name ) )
	    friends.insert( *ff );
	}
      
      // now reset ff and move forwards
      ff = ee;
      while ( 1 )
        {

	  // move forward
	  ++ff;

	  if ( ff == pevents.end() ) break;

          // different segment/indiv?
	  if ( ff->seg != ee->seg ) break;
	  
          // too far away?
	  if ( ff->interval.start - ee->interval.stop > neighbourhood_tp ) break;
	  
          // else add, if aligned
	  if ( aligned( ee->name , ff->name ) )
	    friends.insert( *ff );
        }
      
      // now define interval, starting w/ self
      interval_t n = ee->interval ;
      
      std::set<pinstance_t>::const_iterator gg = friends.begin();
      while ( gg != friends.end() )
	{
	  if ( gg->interval.start < n.start )
	    n.start = gg->interval.start;

	  if ( gg->interval.stop > n.stop )
	    n.stop = gg->interval.stop;

	  ++gg;
	}

      // get indiv
      int indiv = seg2indiv[ ee->seg ];

      // std::cout << " adding span "
      // 		<< indiv << " " << ee->name << " " << ee->interval.as_string() << " " << friends.size() <<  "  span " << n << "\n";
      
      // track friends
      event2friends[indiv][ *ee ] = friends;

      //      std::cout << " event2neighbourhood (span) = " << n << "\n";

      // track total neighbourhood span
      event2neighbourhood[indiv][ *ee ] = n;

      // next permutable event
      ++ee;
    }


  if ( debug_mode )
    {
      std::cout << "\n\n----- Verbose display\n\n";
      
      // verbose display
      std::map<int,std::map<pinstance_t,std::set<pinstance_t> > >::const_iterator kk = event2friends.begin();
      while ( kk != event2friends.end() )
	{
	  
	  int indiv = kk->first;
	  std::cout << "\n\n----- indiv " << indiv << "\n";
	  
	  std::map<pinstance_t,std::set<pinstance_t> >::const_iterator gg = kk->second.begin();
	  
	  while ( gg != kk->second.end() )
	    {
	      
	      std::cout << " pseed "
			<< " " << gg->first.interval.as_string() << " "
			<< " " << gg->first.interval << " " 
			<< " " << gg->first.name << " "
			<< " " << gg->first.seg << "\n";
	      std::cout << "   " << gg->second.size() << " friends, spanning " << event2neighbourhood[ indiv ][ gg->first ] << "\n";
	      std::set<pinstance_t>::const_iterator ff = gg->second.begin();
	      while ( ff != gg->second.end() )
		{
		  std::cout << "   " << ff->interval.as_string() << " " << ff->name << " " << ff->seg << "\n";
		  ++ff;
		}
	      ++gg;
	    }
	  
	  std::cout << " nexxt indiv...\n";
	  ++kk;
	}
      std::cout << "leaving init_event_perm()\n";
    }
  

}




void annotate_t::event_permutation()
{

  // initiate event boundaries with fixed events only
  events = fixed_events;

  // std::cout << "\n\n\n---------------------------------------\nevent_permutation()  -- showing fixed events only\n";
  // view();
  // std::cout << "\n";
  
  // indiv-wise, permutation of seed, eventwise (rather than circular shuffle) 
  std::map<int,std::map<pinstance_t,std::set<pinstance_t> > >::const_iterator kk = event2friends.begin();
  while ( kk != event2friends.end() )
    {
      
    // which individual?
    const int indiv = kk->first;

    // segment starts : indiv2segs[ indiv ] -> set of segments (start tps)
    // segment info: seg[ start_tp ]        -> duration of that segment (tps)

    // algorithm: for each to-be-placed interval (of fixed duration), create a unique map
    //   of possible placements, taking into account the currently dropped events as well
    //   as the core boundaries


    // 0) initialize blackout map w/ gaps between segments for this individual
    // [ in loop [
    // 1) pick event to permute; find friends; determine total span = t
    // 2) create a list of possible placements and track overall span, given current blackout map and 't'
    //   0                ... blackout[0].start  = put1
    //   blackout[0].stop ... blackout[1].start  = put2
    //   blackout[1].stop ... blackout[2].start  = put3 ...
    //   etc.
    // 3) for each putative region (put1, put2, ...) only consider is length of put is >= 't' (or perhaps t+2 to ensure gap each side)
    // 4)  record  put1.start -> put1.length - t 
    //             put2.start -> put2.length - t
    //             ...
    //     and aggregate lengths of each = R
    // 5) pick random seed r = 0 ... R
    // 6) map to which put1 region, and place
    // 7) add this placed interval to the blackout map, i.e. as if a new gap is added
    
    // breaks for this indiv
    const std::set<uint64_t> & isegs = indiv2segs.find( indiv )->second;

    // make original segment intervals (convenience used below)
    std::set<interval_t> iints;
    std::set<uint64_t>::const_iterator aa = isegs.begin();
    while ( aa != isegs.end() )
      {
	interval_t i( *aa , *aa + seg[ *aa ] );
	iints.insert( i );
	++aa;
      }
    
    // initial offset for this indiv
    uint64_t offset = multi_indiv ? indiv_segs[ indiv ].start : 0LLU ; 
    
    // calc total duration from 0 to last segment end below
    uint64_t last_offset = 0;

    
    // create initial blackout map (0-based)
    // only defined if >1 segment
    // for N segments, we should now have N-1 blackout regions
    //  these will a) be added to by any new placed events,
    //             b) used to construct a list of valid placements give span of 't'
    
    std::set<interval_t> blackouts0;
    
    if ( isegs.size() > 1 )
      {
	std::set<uint64_t>::const_iterator ii = isegs.begin();
	// S111   S2222  S3333 
	//    |   |
	//  prior p
	
	// start of first seg + length of first seg = start of first black out
	uint64_t prior = *ii + seg[ *ii ] ;
	
	// go to start of next segment
	++ii;
	
	while ( ii != isegs.end() )
	  {
	    // 0-based start of next segment / end of this blackout
	    uint64_t p = *ii ; 
	    
	    // add to the list
	    //std::cout << " initializing blackouts with " << prior << " " << p << "\n";
	    blackouts0.insert( interval_t( prior , p ) ) ;

	    // update for next time - end of this segment = start of next blackout
	    prior = *ii + seg[ *ii ];

	    // update last time point 
	    last_offset = prior;
	    
	    ++ii;
	  }
      }
    else
      last_offset = *isegs.begin() + seg[ *isegs.begin() ];

    //    std::cout << " last_offset -> " << last_offset << "\n";
    

    // iterate over events in order of alignment groups;
    //  contruct blackouts for each group;  when encountering a new group, wipe
    //  blackout map to original (based only on indiv/background = blackouts0 )
    // this way, only evens w/in the same alignment group cannot be placed on
    // top of each other

    std::map<int,std::set<std::string> >::const_iterator g2sets = group_global_alignment.begin();
    while ( g2sets != group_global_alignment.end() )
      {
	
	
	// events for this indiv
	const std::map<pinstance_t,std::set<pinstance_t> > & ievents_all = kk->second;

	// construct subset of these ievents[] which only includes events from the current group
	std::map<pinstance_t,std::set<pinstance_t> > ievents;
	std::map<pinstance_t,std::set<pinstance_t> >::const_iterator vv = ievents_all.begin();
	while ( vv != ievents_all.end() )
	  {
	    if ( g2sets->second.find( vv->first.name ) != g2sets->second.end() )
	      ievents[ vv->first ] = vv->second;
	    ++vv;
	  }

	// reset the blackouts to default, i.e. fresh for this alignment group
	std::set<interval_t> blackouts = blackouts0;
	
	// randomize order in which we pick events to permute
	const int nev = ievents.size();

	// list events for easy access
	std::vector<pinstance_t> pidx;
	std::map<pinstance_t,std::set<pinstance_t> >::const_iterator ff = ievents.begin();
	while ( ff != ievents.end() ) { pidx.push_back( ff->first ); ++ff; } 
    
	// also idx as an int
	std::vector<int> idx( nev );

	// randomize event picks
	CRandom::random_draw( idx );
	
	// track which events have already been moved
	std::set<pinstance_t> moved;    
	
	// randomly place each event (+friends) in the above (random) order
	for (int e=0; e<nev; e++)
	  {
	    
	    // locate the event to permute
	    pinstance_t & fidx = pidx[ idx[e] ];
	    
	    // all ready done?	
	    if ( moved.find( fidx ) != moved.end() )
	      continue;
	    
	    // make friends - some original friends may have already been
	    // moved (by another index, as this is not a fully transitive
	    // process)
	    
	    const std::set<pinstance_t> & oldfriends = ievents.at( fidx );
	    
	    std::set<pinstance_t> friends;
	    std::set<pinstance_t>::const_iterator ff = oldfriends.begin();
	    while ( ff != oldfriends.end() )
	      {
		if ( moved.find( *ff ) == moved.end() )
		  friends.insert( *ff );
		++ff;
	      }
	    
	    
	    // total interval spanned by index and all true friends
	    
	    interval_t span = event2neighbourhood.at( indiv ).at( fidx );
	    uint64_t span_tp = span.duration();
	    
	    // given the current blackout map and span_tp, create a list of possible
	    // potential placements = places --> place-segement start
	    // but also need to track 'original' segment start (i.e. place-segments will be
	    // created by adding new events -> new gaps -> new place-segments
	    
	    //  SEGMENT1---------------   SEGMENT2-----------------
	    //            |placed|
	    
	    //  a)
	    //  P1--------       P2----   P3-----------------------
	    
	    //  b) becomes:
	    //  P1-------- P2----P3-----------------------
	    
	    //  P1,P2,P3 --> point back to start positions in a)
	    //  but also we want the original SEGMENT the event belongs to
	    //   i.e P1 == P2 --> S1
	    //       P3       --> S2
	    //  but we'll calculate this when needed below, rather than try to track here
	    
	    std::map<interval_t,uint64_t> places; 
	    
	    // start at first segment
	    uint64_t place0 = *isegs.begin();
	    
	    // get aggregate length
	    uint64_t totr = 0LLU;
	    
	    uint64_t max_dur = 0LLU;
	    
	    std::set<interval_t>::const_iterator bb = blackouts.begin();
	    while ( bb != blackouts.end() )
	      {
		
		uint64_t place1 = bb->start;
		
		if ( place0 >= place1 ) Helper::halt( "internal error in creating the blackout list" );
		
		uint64_t dur = place1 - place0 ;
		
		// duration sufficient?
		// this will also ensure that if the indiv interval starts in
		// a blackout, then this segment will not be added
		// nb. contiguous regions should not be merged during eval(), i.e.
		// so should be okay to allow the (unlikely) event that two indpendent
		// regions are permuted to be completely contiguous - they will still
		// count as two separate events
		
		//if ( ( dur - span_tp ) > max_dur ) max_dur = dur - span_tp ;
		
		if ( dur >= span_tp )
		  {
		    // add interval, adter reducing by span_tp, and track the start point
		    places[ interval_t( totr , totr + dur - span_tp ) ] = place0; 
		    // update latest
		    totr += dur - span_tp ;				
		  }
		
		// update for the next potential
		place0 = bb->stop;
		
		++bb;
	      }
	    
	    // need to add final segment?	
	    uint64_t place1 = last_offset;
	    //if ( place0 > place1 ) std::cout <<	"HAH2 MAJOR BLUNDER!\n";
	    
	    uint64_t dur = place1 - place0;
	    //if ( ( dur - span_tp ) > max_dur ) max_dur = dur - span_tp ;
	    if ( dur >= span_tp )
	      {
		places[ interval_t( totr , totr + dur - span_tp ) ] = place0;
		totr += dur - span_tp;	    
	      }
	    
	    //	std::cout << " MAX dur = " << max_dur * globals::tp_duration << "\n";
	    
	    if ( debug_mode )
	      {
		std::cout << " verbose----- event " << e << "\n";
		std::cout << "isegs\n";
		
		std::set<uint64_t>::const_iterator ii = isegs.begin();
		while ( ii != isegs.end() )
		  {
		    std::cout << " " << *ii << " (dur " << seg[ *ii ] << ")\n";
		    ++ii;
		  }
		
		std::cout << " blackouts\n";
		std::set<interval_t>::const_iterator bb = blackouts.begin();
		while ( bb != blackouts.end() )
		  {
		    std::cout << "  " << bb->as_string() << "\n";
		    ++bb;
		  }
		
		std::cout << "span " << span.as_string() << " " << span.duration_sec() << "\n";
		
		std::cout << "places\n";
		std::map<interval_t,uint64_t>::const_iterator pp = places.begin();
		while ( pp != places.end() )
		  {
		    std::cout << " " << pp->first.as_string() << " --> " << pp->second << "\n";
		    ++pp;
		  }
		std::cout << "\n";
	      }
	    

	    
	    //
	    // now find a placement 0 .. totr
	    //
	    
	    // get rr which is 0 .. totr in places[] space	
	    uint64_t rr = CRandom::rand( totr ) ;
	    if ( rr == totr ) --rr;
	    //std::cout << " picked rr = " << rr << " , with totr = " << totr << "\n";
	    
	    // get actual time offset (in tp from 'offset', i.e. t==0 for this person)
	    uint64_t pp = 0LLU;
	    
	    // place-segment start
	    uint64_t place_segstart = 0LLU;
	    
	    // nb. could speed this up...
	    std::map<interval_t,uint64_t>::const_iterator qq = places.begin();
	    //std::cout << "  considering " << places.size() << " places\n";
	    while ( qq != places.end() )
	      {
		//std::cout << " qq->first.start / stop = " << qq->first.start << " " << qq->first.stop << "\n";
		
		if ( rr >= qq->first.start  && rr < qq->first.stop )
		  {
		    // 0 1 2 3 4 5 6 7 8 9 [10]		 
		    uint64_t df = rr - qq->first.start ;
		    
		    // std::cout << "located places- place interval " << qq->first.as_string() << " points to "
		    // 	  << qq->second << " as true segment start, with offset df = " << df << "\n";
		    
		    // set the point: this should /always/  work...
		    pp = qq->second + df; 
		    
		    //		std::cout << " setting pp = " << pp << "\n";
		    
		    // get place-segment start 
		    place_segstart = qq->second;
		    
		    break;
		  }
		++qq;
	      }
	    
	    // get the original segment start that this place-segment starts in
	    //  i.e. this gets us back to the original events[] keying
	    
	    uint64_t segstart;
	    
	    if ( ! annotate_t::get_segment_start( iints , place_segstart , &segstart ) )
	      Helper::halt( "internal error from annotate_t::get_segment_start()" );
	    
	    // std::cout << "  placing event at pp = " << fidx.interval << " " << pp << "; "
	    // 	  << " segstart = " << segstart << "  from " << place_segstart << "\n";
	    
	    //
	    // track this these events have been moved
	    //
	    
	    // track index
	    // and build up the map
	    moved.insert( fidx );
	    
	    // new position : pp is start of new span
	    // thus also add the offset from span start to the event of interest
	    //   start = pp + ( fidx.interval.start - span.start )
	    //   end   = start + fidx.interval.duration
	    
	    uint64_t new_start = pp + ( fidx.interval.start - span.start );
	    uint64_t new_stop  = new_start + fidx.interval.duration();
	    interval_t new_interval = interval_t( new_start , new_stop );
	    events[ segstart ][ fidx.name ].insert( new_interval );
	    
	    //	std::cout << " ADDING A INDEX... " << fidx.interval.as_string() << " --> " << new_interval.as_string() << "\n";
	
	    // & any friends that come along for the ride
	    std::set<pinstance_t>::const_iterator gg = friends.begin();
	    while ( gg != friends.end() )
	      {
		moved.insert( *gg );

		// as above, need to consider local offset relative to start of span
		// (which pp now places)
		//std::cout << "  ++++ A FRIEND " << gg->interval.as_string() << "\n";
		uint64_t new_start = pp + ( gg->interval.start - span.start );
		uint64_t new_stop  = new_start + gg->interval.duration();
		events[ segstart ][ gg->name ].insert( interval_t( new_start , new_stop ) );
		
		++gg;
	      }
	    

	    //
	    // Update the blackouts list to include placed event as a new gap
	    //
	    
	    // new span
	    uint64_t new_span_start = pp ;
	    uint64_t new_span_stop  = pp + span.duration();
	    blackouts.insert( interval_t( new_span_start , new_span_stop ) );
	    
	    // next index event
	  }
	
	// next alignent-group of annots
	++g2sets;
      }
    
    // next indiv
    ++kk;
    }
  

  // std::cout << "\n\n---------- FINAL VIEW ---------------\n\n";
  // view();
       
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
	  // and so this would change the overall number of annots.  Also, will work
	  // better with event-permutation, where in theory two events could be placed
	  // right next to each other
	  
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

	      // skip seed-seed comparison?
	      if ( (!do_seed_seed ) && sachs.find( *bb ) != sachs.end() ) { ++bb; continue; }
		   
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

	  //
	  // consider any non-segment markers
	  //

	  std::set<std::string>::const_iterator mm = mannots.begin();
	  while ( mm != mannots.end() )
	    {
	      
	      // should not happen, but skip any self comparison
	      if ( *aa == *mm ) { ++mm; continue; }
	      
	      // markers do not have channel information, (perhaps can change in future)
	      // but for now no need to check the 'same_channel()' test
	      
	      // the offset for this seed segment can be used to figure out the
	      // indiv/bin for
	      
	      int indiv = multi_indiv ? indiv_segs[ seg2indiv[ offset ] ].start : 0 ; 
	      
	      // get the markers
	      const std::set<interval_t> & mrk = markers[ indiv ][ *mm ];
	      
	      // calculate stats - can use the same seed_annot_stats()
	      // but swap out window_sec first (then replace)
	      double orig_window_sec = window_sec;
	      window_sec = marker_window_sec;
	      seed_annot_stats( a , *aa ,  mrk , *mm , offset , &r );
	      window_sec = orig_window_sec;
	      
	      ++mm;
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


void annotate_t::add_contrasts( annotate_stats_t * r )
{

  // contrast=SP11-SO|SP15,SO
  // initially, only modify nsa[ a ][ b ]
  
  //          seed,annot|new-seed,annot
  //          seed,annot|seed,new-annot
  //          seed,seed2|seed,new-seed2

  //  logger << "  attempting to add " << contrasts.size() << " contrasts\n";
  for (int i=0; i<contrasts.size(); i++)
    {

      annot_contrast_t c = contrasts[i];

      // search on adist[][]
      // then make contrasts for adist, sdist and nsa
      
      //std::cout << " c " << c.a1 << " " << c.b1 << " ----> " << c.a2 << " " << c.b2 << "\n";
      //std::cout << " sz " << r->nsa.size() << "\n";
      // std::map<std::string,std::map<std::string,double> >::const_iterator qq = r->nsa.begin();
      // while ( qq != r->nsa.end() ) { std::cout << " ---> " << qq->first << "\n"; ++qq ; } 
	    
      // do we find a1 and b1?
      std::map<std::string,std::map<std::string,double> >::const_iterator nn = r->adist.find( c.a1 );
      
      if ( nn == r->adist.end() )
	continue;

      // std::map<std::string,double>::const_iterator qq2 = nn->second.begin();
      // while ( qq2 != nn->second.end() ) { std::cout << " ---> " << qq2->first << "\n"; ++qq2 ; }
      
      if ( nn->second.find( c.b1 ) == nn->second.end() )
	continue;

      std::map<std::string,std::map<std::string,double> >::const_iterator mm = r->adist.find( c.a2 );

      if ( mm == r->adist.end() )
	continue;
      
      // qq2 = nn->second.begin();
      // while ( qq2 != nn->second.end() ) { std::cout << " ---> " << qq2->first << "\n"; ++qq2 ; }

      if ( mm->second.find( c.b2 ) == mm->second.end() )
	continue;

      // make new scores - need denom added here to get diff in average scores

      double sc_adist = r->adist[ c.a1 ][ c.b1 ] / ( r->nadist[ c.a1 ][ c.b1 ] == 0 ? 1 : r->nadist[ c.a1 ][ c.b1 ] )
	- r->adist[ c.a2 ][ c.b2 ] / ( r->nadist[ c.a2 ][ c.b2 ] == 0 ? 1 : r->nadist[ c.a2 ][ c.b2 ] ) ;
      
      double sc_sdist = r->sdist[ c.a1 ][ c.b1 ] / ( r->nsdist[ c.a1 ][ c.b1 ] == 0 ? 1 : r->nsdist[ c.a1 ][ c.b1 ] )
	- r->sdist[ c.a2 ][ c.b2 ] / ( r->nsdist[ c.a2 ][ c.b2 ] == 0 ? 1 : r->nsdist[ c.a2 ][ c.b2 ] ) ;
      
      //double sc_sdist = r->sdist[ c.a1 ][ c.b1 ] - r->sdist[ c.a2 ][ c.b2 ];

      double sc_nsa   = r->nsa[ c.a1 ][ c.b1 ]   - r->nsa[ c.a2 ][ c.b2 ];
           
      const std::string l1 = c.a1 == c.a2 ? c.a1 : c.a1 + "-" + c.a2;
      const std::string l2 = c.b1 == c.b2 ? c.b1 : c.b1 + "-" + c.b2;

      //      std::cout << " CONTRASTS = " << l1 << " " << l2 << " -- " << sc_adist << " " << sc_sdist << " " << sc_nsa << "\n";
      
      // make a new label and store
      r->adist[ l1 ][ l2 ] = sc_adist;
      r->sdist[ l1 ][ l2 ] = sc_sdist;
      r->nsa[ l1 ][ l2 ] = sc_nsa;
      
      // need to set ndist to 1.0 - i.e. has no meaning/effect, but means that distance measures are seen:
      r->nadist[ l1 ][ l2 ] = 1.0;
      r->nsdist[ l1 ][ l2 ] = 1.0;
      
    }

}



void annotate_t::output( const annotate_stats_t & s )
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

      //      std::cout << " sa obs2 " << sa->first << "\n";
      
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


  //
  // seed-annot flanking overlap window
  //

  if ( n_flanking_offsets )
    {

      std::vector<int> fidx;
      std::vector<double> fmid;
      std::vector<std::string> fstr;
      for (int fi=n_flanking_offsets; fi!=0; fi--)
	{
	  fidx.push_back( -fi );
	  fmid.push_back( -1 * flanking_overlap_mid[ fi-1 ] );
	  fstr.push_back( "-" + flanking_overlap_desc[ fi-1 ] );
	}
      for (int fi=1; fi<=n_flanking_offsets; fi++)
	{
	  fidx.push_back( fi );
	  fmid.push_back( flanking_overlap_mid[ fi-1 ] );
	  fstr.push_back( flanking_overlap_desc[ fi-1 ] );
	}
      
      std::map<std::string,std::map<std::string,std::map<int,double> > >::const_iterator sf = pf_obs.begin();
      while ( sf != pf_obs.end() )
	{
	  
	  if ( ! ( flanking_overlap_seeds.size() == 0
		   || flanking_overlap_seeds.find( sf->first ) != flanking_overlap_seeds.end() ) )
	    {
	      ++sf; continue;
	    }    
	    
	  writer.level( sf->first , "SEED" );

	  const std::map<std::string,std::map<int,double> > & p = sf->second;
	  std::map<std::string,std::map<int, double> >::const_iterator pp = p.begin();
	  while ( pp != p.end() )
	    {

	      if ( ! ( flanking_overlap_others.size() == 0
		       || flanking_overlap_others.find( pp->first ) != flanking_overlap_others.end() ) )
		{
		  ++pp; continue;
		}
	      
	      writer.level( pp->first , "OTHER" );

	      // always do all rather than just observed
	      

	      for (int fi=0; fi<fidx.size(); fi++)
		{
		  
		  int pos = fidx[fi];		  
		  
		  // pos
		  writer.level( fmid[fi] , "OFFSET" );
		  writer.value( "INT" , fstr[fi] );
		  writer.value( "N_OBS" , pf_obs[ sf->first ][ pp->first ][ pos ]  );
		  if ( nreps )
		    {
		      double mean = pf_exp[ sf->first ][ pp->first ][ pos ] / (double)nreps;
		      double var = pf_expsq[ sf->first ][ pp->first ][ pos ] / (double)nreps - mean * mean;	  
		      writer.value( "N_EXP" , mean );
		      writer.value( "N_P" , ( pf_pv[ sf->first ][ pp->first ][ pos ]  + 1 ) / (double)( nreps + 1 ) );
		      if ( var > 0 ) 
			writer.value( "N_Z" , ( pf_obs[ sf->first ][ pp->first ][ pos ] - mean ) / sqrt( var ) );
		    }

		}
	      writer.unlevel( "OFFSET" );

	      // next annot
	      ++pp;
	    }	  
	  writer.unlevel( "OTHER" );

	  // next seed
	  ++sf;
	}
      writer.unlevel( "SEED" );
    }

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

  
  // read from r->nsa[][] to ensure that nsa will
  // always contain all necessary keys (whether or not combos are seen
  // in the observed data)
  
  int dummy = r->nsa[ astr ][ bstr ] ; 
  
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

      // track closet match (used in offsets[] calcs below)
      // as we tweak *bb
      std::set<interval_t>::const_iterator closestb = bb;
      
      // verbose output
      //if ( debug_mode )
      //      std::cout << "a = " << aa->as_string() << "\n";

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
      // std::cout << " initial b = ";
      // if ( bb == b.end() ) std::cout << " -END-";
      // else
      //   {
      //     std::cout << bb->as_string() ;
      //     if ( bb == b.begin() ) std::cout << " (begin)";
      //   }	  
      // std::cout << "\n";
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
	  //  std::cout << " prov dist = " << dist << "\n";
	  
	  // step back, if we can - is there a closer annot /before/ the seed?
	  if ( bb != b.begin() )
	    {
	      --bb;

	      // if ( debug_mode )
	      // std::cout << " stepping back, b -> " << bb->as_string() << "\n";
	      
	      // nb - this may overlap seed
	      // i.e. starts before, but ends after seed-start, and so
	      //  was not captured by the lower_bound()
	      
	      if ( bb->stop > aa->start )
		{
		  dist = 0;
		  overlap = true;

		  // track
		  closestb = bb;
		  
		  // if ( debug_mode )
		  // std::cout << "  back b overlaps a, done\n"; 
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
		    {
		      dist = - left_dist;
		      closestb = bb;
		    }
		  // if ( debug_mode )
		  //   {
		  // std::cout << " back b is before, so dist = " << left_dist << "\n";
		  // std::cout << " final dist = " << dist << "\n";
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
	  //std::cout << "   size = " << r->psa[ astr ].size() << "\n";
	}
      
      // truncate at window length?
      const bool truncate = dist > window_sec || dist < -window_sec ;

      // sets abs-dist to max;
      // ignores in calc of signed-dist
      double adist = truncate ? window_sec : fabs( dist ) ;
      
      //      std::cout << "a, closest = " << aa->as_string() << "\t" << closestb->as_string() << "\t" << overlap << "\t" << dist << "\n";
            
      // do we include complete overlap as "nearest"?
      if ( include_overlap_in_dist || ! overlap )
	{

	  r->adist[ astr ][ bstr ] += adist ; 
	  r->nadist[ astr ][ bstr ] += 1; 
	  
	  // track only -1 or +1 for 'before' or 'after'
	  // overlap == 0 here, so add that qualifier

	  // nb. above:
	  //   +1 implies ANNOT is AFTER SEED
	  //   -1         ANNOT is BEFORE SEED
	  //  we want the output to be SEED-CENTRIC,
	  //   i..e. -1 means SEED comes BEFORE
	  //       and so change sign of sdist below

	  const double seed_centric_dist = - dist ; 
	  
	  if ( ! overlap )
	    {
	      if ( ! truncate )
		{
		  if ( d2_signed )  // reduce to -1/+1
		    r->sdist[ astr ][ bstr ] += seed_centric_dist > 0 ? +1 : -1 ;
		  else
		    r->sdist[ astr ][ bstr ] += seed_centric_dist ;

		  // denom
		  r->nsdist[ astr ][ bstr ] += 1;

		}
	    }
	  
	}


      //
      // evaluate any offsets windows? 
      //

      if ( n_flanking_offsets 
	   && ( flanking_overlap_seeds.size() == 0 || flanking_overlap_seeds.find( astr ) != flanking_overlap_seeds.end() )
	   && ( flanking_overlap_others.size() == 0 || flanking_overlap_others.find( bstr ) != flanking_overlap_others.end() ) 
	   )
	{

	  // seed          is   *aa
	  // closest match is   *closestb
	  // max_range     is   flanking_overlap_mx
	  // windows       are  flanking_overlap_intervals[]

	  //                  |---aa---|
	  //  |   |   |   |   |        |   |   |   |   |
	  //    -4  -3  -2  -1           +1  +2  +3  +4
	  // -mx                                      +mx

	  std::set<interval_t>::const_iterator cc = closestb;
	  
	  // std::cout << "\n\n CHECKING " << astr << " " << bstr << "\n";
	  // std::cout << " aa " << aa->as_string() << "\n";
	  // std::cout << " cc " << cc->as_string() << "\n";
	  
	  // forwards:: events must span after
	  while ( 1 )
	    {
	      // nothing left?
	      if ( cc == b.end() ) break;	      
	      
	      //	      std::cout << "  checking -> cc " << cc->as_string() << "\n";
	      
	      // comes before, i.e. this cc does not extend after, then advance 
	      if ( cc->stop <= aa->stop )
		{
		  //		  std::cout << " advancing...\n";
		  ++cc;
		  continue;
		}

	      // if here, cc must extend from aa forwards; check we have
	      // not gone too far : from end of seed to start of cc; but
	      // need to check that cc does actually start after end of aa
	      
	      if ( cc->start >= aa->stop && cc->start - aa->stop > flanking_overlap_mx )
		break;
	      
	      // otherwise, we have at least some of the range between end of *aa and
	      // final search region spanned by this *cc.   Count which bins have overlap
	      
	      // zero-bounded start of overlap
	      uint64_t s1 = cc->start < aa->stop ? 0LLU : cc->start - aa->stop ;
	      
	      // end of overlap: we know cc ends after 
	      uint64_t s2 = cc->stop - aa->stop;

	      //  std::cout << " s12 " << s1 << " " << s2 << "\n";
	      
	      for (int fi=0; fi<n_flanking_offsets; fi++)
		{
		  const interval_t & win = flanking_overlap_intervals[fi];
		  
		  // gone past (s2 is end+1 still)
		  if ( win.start >= s2 ) break;
		  
		  //     |aaaaaa|  [s1]------[s2]
		  //            |   |   |   |   |   |   |   win[] interval_t
		  // overlaps?
		  //		  std::cout << " win " << fi << " of " << n_flanking_offsets << " = " << win.as_string() << "\n";
		  
		  if ( s1 < win.stop && s2 > win.start ) 
		    {
		      //std::cout << "overlaps!\n";
		      r->nosa[ astr ][ bstr ][ fi + 1 ] += 1 ; // +1 based offset counts:
		    }
		}
	      
	      // loop back to consider next putative *cc in forwards direction
	      ++cc;				
	    }

	  //std::cout << " now going back\n";
	  
	  // now consider going backwards: reset to closest
	  cc = closestb;

	  while ( 1 )
	    {

	      //std::cout << " chking neg " << cc->as_string() << "\n";

	      // starts after start?, i.e. this cc does not extend before?
	      // then roll back (unless we are already at the start of the list) 
	      if ( cc->start >= aa->start )
		{
		  //std::cout << "  cc starts after aa, roll back...\n";
		  if ( cc == b.begin() ) break; 
		  --cc;
		  continue;
		}
	      
	      // if here, cc must start before aa start; check we have
	      // not gone too far : from start of seed to end of cc;
	      //  but need to check that end of seed to start of cc; but
	      // need to check that cc does actually start after end of aa
	      
	      if ( aa->start >= cc->stop && aa->start - cc->stop > flanking_overlap_mx )
		{
		  //std::cout << " gone too far::: " << cc->stop << " " << aa->start << " " << aa->start - cc->stop << " " << flanking_overlap_mx << "\n";
		  break;
		}

	      // otherwise, we have at least some of the range between start of *aa and
	      // prior final search region spanned by this *cc.   Count which bins have overlap
	      
	      // start of offset (we know cc starts before aa starts if here)
	      uint64_t s2 = aa->start - cc->start;

	      // zero-bounded end of overlap (going back in time)
	      uint64_t s1 = cc->stop < aa->start ? aa->start - cc->stop : 0LLU ;
	      
	      //std::cout << " s1 s2 " << s1 << " " << s2 << "\n";
	      
	      for (int fi=0; fi<n_flanking_offsets; fi++)
		{
		  const interval_t & win = flanking_overlap_intervals[fi];
		  
		  //std::cout << " win " << win.as_string() << "\n";
		  
		  // gone past , based on end of win: no point in looking here.
		  if ( win.start >= s2 ) 
		    {
		      // std::cout << " huh, " << win.stop << "  " << s2 << "\n";
		      // std::cout << "gone past?\n";
		      break;		      
		    }

		  //     |aaaaaa|  [s1]------[s2]
		  //            |   |   |   |   |   |   |   win[] interval_t
		  // overlaps?
		  //std::cout << " checking neg pos  " << -( fi + 1 ) << "\n";
		  if ( s1 < win.stop && s2 > win.start ) 
		    {
		      r->nosa[ astr ][ bstr ][ -( fi + 1 ) ] += 1 ; // -ve +1 based offset counts:
		      //std::cout << " FOUND ONE OH MY!\n";
		    }
		  
		}
	      
	      // loop back to consider next putative *cc in forwards direction
	      if ( cc == b.begin() ) break;
	      --cc;

	    }

	  //std::cout << " done1\n";
	}
      
      // are we tracking hits
      if ( make_anew )
	{
	  
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

bool annotate_t::overlaps_flattened_set(  const interval_t & a , const std::set<interval_t> & b )
{

  // std::cout << " qry = " << a << "\n";
  // std::set<interval_t>::const_iterator bb1 = b.begin();
  // while ( bb1 != b.end() )
  //   {
  //     std::cout << " --> " << *bb1 << "\n";
  //     ++bb1;
  //   }

	    
  
  // this assumes that 'b' is flattened already

  // nothing to overlap?
  if ( b.size() == 0 ) return false;
  
  // find the first annot not before (at or after) the seed
  std::set<interval_t>::const_iterator bb = b.lower_bound( a );


  // does first take overlap seed?
      
  if ( bb != b.end() && bb->overlaps( a ) )
    {
      // we're done, found complete overlap
      return true;      
    }
  else // it must come afterwards
    {
      
      // step back, if we can - is there a closer annot /before/ the seed?
      if ( bb != b.begin() )
	{
	  
	  --bb;
	  
	  // nb - this may overlap seed
	  // i.e. starts before, but ends after seed-start, and so
	  //  was not captured by the lower_bound(); as b is flattened,
	  // only need to step back once
	  
	  if ( bb->stop > a.start )
	    return true;
	  else
	    return false;
	}
    } 
  return false;
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
  
  // for pairwise combinations (nsa and nosa, nadist, nsdist, etc), only nsa is guaranteed to 
  // contain all keys;  therefore, here we should ensure that zero's are entered for nosa
  // and nxdist[]

  // seed-seed group overlap ( std::map<std::string,double> )
  obs = s.nss;
  
  // seed-annot pairwise overlap (std::map<std::string,std::map<std::string,double> > )
  p_obs = s.nsa;
  
  // flanking offset seed-annot pairwise overlaps
  pf_obs = s.nosa;
  
  // seed-annot proportion spanned
  std::map<std::string,std::set<named_interval_t> >::const_iterator pp = s.psa.begin();
  while ( pp != s.psa.end() )
    {
      // std::cout << " prop obs " << pp->first << " " << pp->second.size() << " " << s.ns.find(  pp->first )->second << " = "
      //  		<< pp->second.size() / s.ns.find(  pp->first )->second << "\n";
      prop_obs[ pp->first ] = pp->second.size() / s.ns.find(  pp->first )->second;
      ++pp;
    }
  
  // absolute distance (from each seed to nearest annot) std::map<std::string,std::map<std::string,double> >
  absd_obs = s.adist;
  
  // signed distance (from each seed to nearest annot) std::map<std::string,std::map<std::string,double> >
  sgnd_obs = s.sdist;
  
  // get average distances (these may be < S-A count, because of window_sec threshold)
  // isa nsa[] keys, but check that nxdist is >0 

  std::map<std::string,std::map<std::string,double> >::const_iterator dd = s.nsa.begin();
  while ( dd != s.nsa.end() )
    {
      
      const std::map<std::string,double> & e = dd->second;
      std::map<std::string,double>::const_iterator ee = e.begin();
      while ( ee != e.end() )
	{
	  
	  // handle nosa here, to ensure all keys are populated (based on nsa)
          // i.e. ensure that pf_obs[] can be used as fully populated
	  for (int fi=0; fi<n_flanking_offsets; fi++)
	    {
	      int dummy1 = pf_obs[ dd->first ][ ee->first ][ fi + 1 ]; 
	      int dummy2 = pf_obs[ dd->first ][ ee->first ][ - ( fi + 1 )  ]; 
	    }
	  
	  // handle absd_obs, sgnd_obs and dn_obs to ensure all keys are populated
	  
	  double nad = -1, nsd = -1;
	  
	  std::map<std::string,std::map<std::string,double> >::const_iterator nn = s.nadist.find( dd->first );
          if ( nn != s.nadist.end() )
	    {
	      const std::map<std::string,double> & n2 = nn->second;
	      std::map<std::string,double>::const_iterator nnn = n2.find( ee->first );
	      
	      if ( nnn != n2.end() )  
		nad = nnn->second; // set a (non-zero) value
	    }

	  nn = s.nsdist.find( dd->first );
	  if ( nn != s.nsdist.end() )
	    {
	      const std::map<std::string,double> & n2 = nn->second;
	      std::map<std::string,double>::const_iterator nnn = n2.find( ee->first );
	      
	      if ( nnn != n2.end() )  
		nsd = nnn->second; // set a (non-zero) value
	    }
	  
	  if ( nad > 0 ) 
	    {
	      absd_obs[ dd->first ][ ee->first ] /= nad;
	      dn_obs[ dd->first ][ ee->first ] = nad; // hmm? use abs (nad not nsd) - we prob don't need this...
	    }
	  else
	    {
	      absd_obs[ dd->first ][ ee->first ] = window_sec; // set to max (for an annot: Q markers?)	      
	      dn_obs[ dd->first ][ ee->first ] = 0;
	    }

	  if ( nsd > 0 )	    
	    sgnd_obs[ dd->first ][ ee->first ] /= nsd;
	  else
	    sgnd_obs[ dd->first ][ ee->first ] = 0; // null value 'not-defined'
	  
	  ++ee;
	}
      ++dd;
    }

  // 1-to-many seed-to-annots mappings
  s2a_obs = s2a_proc( s.s2a_mappings );

}


void annotate_t::build_null( annotate_stats_t & s )
{

  //
  // seed-seed group overlap
  //
  
  // consider only observed values (observed combos for 1+ seed pileup)
  
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
      const std::map<std::string,double> & p = sa->second;      
      std::map<std::string,double>::const_iterator pp = p.begin();
      while ( pp != p.end() )
	{	  
	  // nsa[] should always exist, but even if not, would get '0' here, so okay
	  double val = s.nsa[ sa->first ][ pp->first ];
	  
	  p_exp[ sa->first ][ pp->first ] += val;
	  p_expsq[ sa->first ][ pp->first ] += val * val;
	  if ( val >= p_obs[ sa->first ][ pp->first ] ) ++p_pv[ sa->first ][ pp->first ];
	  
	  ++pp;
	}
      ++sa;
    }


  //
  // seed-annot flanking window overlap
  //
  
  std::map<std::string,std::map<std::string,std::map<int,double> > >::const_iterator sf = pf_obs.begin();
  while ( sf != pf_obs.end() )
    {
      
      const std::map<std::string,std::map<int,double> > & p = sf->second;
      std::map<std::string,std::map<int,double> >::const_iterator pp = p.begin();
      while ( pp != p.end() )
	{	  
	  // for each flanking val: all these will exist in the obs
	  std::map<int,double>::const_iterator ff = pp->second.begin();
	  while ( ff != pp->second.end() )
	    {
	      // get the permuted value; it might not exist, but then this gets set to 0 so okay		  
	      double val = s.nosa[ sf->first ][ pp->first ][ ff->first ] ; 
	      pf_exp[ sf->first ][ pp->first ][ ff->first ] += val;
	      pf_expsq[ sf->first ][ pp->first ][ ff->first ] += val * val;
	      if ( val >= pf_obs[ sf->first ][ pp->first ][ ff->first ] ) ++pf_pv[ sf->first ][ pp->first ][ ff->first ];
	      ++ff;	    
	    }
	  ++pp;
	}
      ++sf;
    }

  //  std::cout << "in BN(2)\n";

  //
  // prop-seed overlap
  //
  
  std::map<std::string,double>::const_iterator pp = prop_obs.begin();
  while ( pp != prop_obs.end() )
    {

      const bool is_seen = s.psa.find( pp->first ) != s.psa.end();

      // set to 0.0 if not any observed in the permuted
      double val = is_seen ? s.psa.find( pp->first )->second.size() / s.ns.find( pp->first )->second : 0 ;       
      prop_exp[ pp->first ] += val;
      prop_expsq[ pp->first ] += val * val;
      if ( val >= prop_obs[ pp->first ] ) ++prop_pv[ pp->first ];	  
    
      ++pp;
    }

  
  //
  // seed-annot distances : sgnd_obs and absd_obs may now have different keys/denoms
  //    but absd should always be the superset, and sgnd = 0 if not, so just use absd
  //
  
  sa = absd_obs.begin();
  while ( sa != absd_obs.end() )
    {
     
      const std::map<std::string,double> & p = sa->second;
      
      std::map<std::string,double>::const_iterator pp = p.begin();
      while ( pp != p.end() )
        {
          
	  // keys: sa->first  pp->first 
	  //  access s.ndist[][]   s.sdist[][]   s.adist[][] 
	  
	  // if not defined, ndist and sdist should be zero
	  // otherwise, adist should be max window size (i.e. 'undefined')
	  // Q/TODO: annot vs marker max value issue
	  	  
	  const double na = s.nadist[ sa->first ][ pp->first ];
	  const double ns = s.nsdist[ sa->first ][ pp->first ];
	  
	  double a1 = window_sec;
	  double s1 = 0; 
	  
	  // normalize here if defined

	  if ( na > 0 ) 
	    a1 = s.adist[ sa->first ][ pp->first ] / na ; 
	  
	  if ( ns > 0 )
	    s1 = s.sdist[ sa->first ][ pp->first ] / ns ; 		  
	  
	  //	  std::cout << "sgnd\t" << s1 << "\t" << a1 << "\t" << na << "\t" << ns  << "\t" << s.nsa[ sa->first ][ pp->first ]  << "\n";
	  
	  // expecteds: sums
	  absd_exp[ sa->first ][ pp->first ] += a1;
	  sgnd_exp[ sa->first ][ pp->first ] += s1;
	  
	  // sum of sqs
	  absd_expsq[ sa->first ][ pp->first ] += a1 * a1;
	  sgnd_expsq[ sa->first ][ pp->first ] += s1 * s1;
	      
	  // track counts (prob. don't need??)
	  // just use abs() value here
	  dn_exp[ sa->first ][ pp->first ] += na;
	  
	  // pvals : one-sided p-values
	  if ( a1 <= absd_obs[ sa->first ][ pp->first ] ) ++absd_pv[ sa->first ][ pp->first ];
	  
	  // nb. we've calculated mean pre/prior, but here the test is 2-sided, so take abs(x)
	  if ( fabs(s1) >= fabs( sgnd_obs[ sa->first ][ pp->first ] ) ) ++sgnd_pv[ sa->first ][ pp->first ];
	  
          ++pp;
        }
      ++sa;
      //      std::cout << "s7\n";
    }
  //  std::cout << "s-done\n";


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

  //  std::cout << "done BN\n";  
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
      // expect annot~ch or annot|ch
      std::vector<std::string> tok2 = Helper::parse( tok[i] , "~|" );
      if ( tok2.size() != 2 ) Helper::halt( "expecting annot|ch or annot~ch format for chs-inc and chs-exc" );

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
			<< "interval = " << ii->as_string() << "\t"
			<< "dur = " << ii->duration_sec() << "\t"
			<< ii->start << "\n";
              ++ii;
            }
          ++qq;
        }
      ++rr;
    }

  
  rr = markers.begin();
  while ( rr != markers.end() )
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
			<< "marker = " << qq->first << "\t"
			<< "interval = " << ii->as_string() << "\t"
			<< "dur = " << ii->duration_sec() << "\n";
              ++ii;
            }
          ++qq;
        }
      ++rr;
    }

  std::cout << "\n";
  
}



std::set<std::string> annotate_t::root_match( const std::string & s , const std::vector<std::string> & names )
{
  std::vector<std::string> tok = Helper::parse( s , "," );
  std::set<std::string> ss;
  for (int i=0; i<tok.size(); i++)
    ss.insert( tok[i] );
  return root_match( ss , names );
}

std::set<std::string> annotate_t::root_match( const std::set<std::string> & s , const std::vector<std::string> & names )
{
  std::set<std::string> r;

  std::set<std::string>::const_iterator ss = s.begin();
  while ( ss != s.end() )
    {
      // empty
      if ( ss->size() == 0 ) { ++ss; continue; }

      // vanilla? add as is (whether in the annot list or not) 
      if ( (*ss)[ ss->size() - 1 ] != '*' )
	{
	  r.insert( *ss );
	}
      else
	{
	  // ignore if only a *
	  if ( ss->size() == 1 ) { ++ss; continue; }

	  // expand if ends in wildcar: 'annot*'
	  const std::string root = ss->substr( 0 , ss->size() - 1 );
	  const int rootsize = root.size();
	  
	  for (int i=0; i<names.size(); i++)
	    {
	      if ( names[i].size() < rootsize ) continue;
	      if ( names[i].substr( 0 , rootsize ) == root )
		r.insert( names[i] );
	    }	  
	}
      
      ++ss;
    }
  return r;
}


bool annotate_t::get_segment_start( const std::set<interval_t> & y , uint64_t x , uint64_t * start )
{
  if ( y.size() == 0 ) return false;
  
  // assumption: y is non-overlapping segments;
  // return start of segment in which x falls
  interval_t xx( x , x );

  //An iterator to the the first element in the container which is considered to go after val, or set::end if no elements are considered to go after val.
  std::set<interval_t>::const_iterator ll = y.upper_bound( xx );
  
  if ( ll != y.end() )
    {
      if ( x >= ll->start && x < ll->stop )
	{
	  *start = ll->start;
	  return true;
	}
    }

  if ( ll == y.begin() ) return false;
  
  // else wind back
  --ll;

  if ( x >= ll->start && x < ll->stop )
    {
      *start = ll->start;
      return true;
    }

  return false;
    
}
