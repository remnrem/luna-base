
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
#include "timeline/cache.h"

#include "edf/edf.h"
#include "edf/slice.h"

#include "db/db.h"
#include "helper/logger.h"

extern writer_t writer;
extern logger_t logger;

void ctest()
{
  
  writer.level( "L1" , "F1" );
  
  writer.level( 123 , "FFE" );
  writer.epoch( 222 );
  
  cache_t<double> cache( "my1" );

  ckey_t ckey1( "y" , writer.faclvl() ) ;
  ckey_t ckey2( "z" , writer.faclvl() ) ;

  std::vector<double> y( 10 , 22 );
  std::vector<double> z( 10 , 23 );
  
  cache.add( ckey1 , y );
  cache.add( ckey2 , z );
  
  writer.unlevel();

  std::cout << cache.print();
  
}

void ctest2( edf_t & edf )
{
  std::string cmd = "PSD";
  std::string var = "RELPSD";

  std::map<std::string,std::string> fac;
  fac[ "CH" ] = "C3";
  fac[ "B" ] = "SIGMA";

  cache_t<double> * c = edf.timeline.cache.find_num( "c1" );

  std::vector<double> d = c->fetch( cmd , var , fac );

  std::cout << " d size = " << d.size() << "\n";
  for (int i=0; i<d.size(); i++)
    std::cout << " d = " << d[i] << "\n";

  std::cout<< "fetch1()...\n";
  double d1 = 0;
  if ( c->fetch1( cmd , var , fac , &d1 ) )
    std::cout << " d1 = " << d1 << "\n";
    
}


void caches_t::load( const std::string & filename )
{
  
  std::ifstream IN1( filename.c_str() , std::ios::in );
  

  // 
  // format : 
  //
  // cache: peaks[int]
  // strata: fac=lvl
  // strata: clear
  // value: var=123.456
  

  std::map<std::string,std::string> curr_strata;
  
  cache_t<std::string> * str_cache = NULL ;
  cache_t<double> * num_cache = NULL ;
  cache_t<int> * int_cache = NULL ;
  cache_t<uint64_t> * tp_cache = NULL ;
  
  int cnt = 0 ;

  while ( ! IN1.eof() )
    {
      
      std::string line;

      Helper::safe_getline( IN1 , line );
      
      if ( IN1.eof() ) break;
      if ( line == "" ) continue;
      
      std::vector<std::string> tok = Helper::parse( line , "\t " );
      if ( tok.size() != 2 ) Helper::halt( "problem with cache format: " + line );
      
      if ( tok[0] == "cache:" ) 
	{
	  std::vector<std::string> tok2 = Helper::parse( line , "[]" );

	  str_cache = NULL;
	  int_cache = NULL;
	  num_cache = NULL;
	  tp_cache = NULL;

	  if ( tok2.size() != 2 ) Helper::halt( "problem with cache format: " + line );
	  if ( tok2[1] == "int" )
	    int_cache = find_int( tok2[0] );
	  else if ( tok2[1] == "num" )
	    num_cache = find_num( tok2[0] );
	  else if ( tok2[1] == "str" )
	    str_cache = find_str( tok2[0] );
	  else if ( tok2[1] == "tp" )
	    tp_cache = find_tp( tok2[0] );
	  else
	    Helper::halt( "problem with cache format: " + line );
	  logger << "reading into " << tok2[0] << "\n";
	}
      else if ( tok[0] == "strata:" )
	{
	  if ( tok[1] == "clear" ) curr_strata.clear();
	  else {
	    std::vector<std::string> tok2 = Helper::parse( tok[1], "=" );
	    if ( tok2.size() != 2 )  Helper::halt( "problem with cache format: " + line);
	    curr_strata[ tok2[0] ] = tok2[1];
	  }
	}
      else if ( tok[0] == "value:" )
	{
	  std::vector<std::string> tok2 = Helper::parse( tok[1], "=" );
	  if ( tok2.size() !=2 )  Helper::halt( "problem with cache format: " + line);
	  
	  if ( num_cache )
	    {
	      double num_value;
	      if ( ! Helper::str2dbl( tok2[1] , &num_value ) ) Helper::halt( "problem with cache format: " + line);
	      num_cache->add( ckey_t( tok2[0] , curr_strata ) , num_value ) ;
	      std::cout << " adding " << tok2[0] << " --> " << tok2[1] << "\n";
	    }
	  else if ( int_cache )
	    {
	      int int_value;
	      if ( ! Helper::str2int( tok2[1] , &int_value ) ) Helper::halt( "problem with cache format: " + line);
	      int_cache->add( ckey_t( tok2[0] , curr_strata ) , int_value ) ;
	    }
	  else if ( str_cache )
	    {
	      str_cache->add( ckey_t( tok2[0] , curr_strata ) , tok2[1] );
	    }
	  else if ( tp_cache )
	    {
	      uint64_t tp_value;
	      if ( ! Helper::str2int64( tok2[1] , &tp_value ) ) Helper::halt( "problem with cache format: " + line);
	      tp_cache->add( ckey_t( tok2[0] , curr_strata ) , tp_value ) ; 
	    }
	  else
	    Helper::halt( "problem with cache format: " + line);
	  
	  ++cnt;

	}
      else
	Helper::halt( "problem with cache format: " + line);

      // next row
    }
  
  IN1.close();

  logger << "  read " << cnt << " values from " << filename << "\n";

  std::cout << " print \n\n" << num_cache->print() << "\n\n---\n";
  
}



void caches_t::import( const std::string & filename , 
		       const std::string & cache_name , 
		       const std::string & id , 
		       const std::set<std::string> & factors , 
		       const std::set<std::string> * variables )
{
  
  // we assume all values are NUMERIC   
  cache_t<double> * num_cache = find_num( cache_name );
  
  std::ifstream IN1( filename.c_str() , std::ios::in );
  
  // process headers
  std::map<std::string,int> fslot, vslot;
  
  std::string line;
  
  Helper::safe_getline( IN1 , line );

  if ( IN1.eof() ) Helper::halt( "problem reading " + filename );
  if ( line == "" ) Helper::halt( "problem reading " + filename );
  std::vector<std::string> tok = Helper::parse( line , "\t " );
  if ( tok.size() <= 2 ) Helper::halt( "problem with imported format: need at least two cols:\n" + line );
  if ( tok[0] != "ID" ) Helper::halt( "bad header row: first col should be ID" );
  
  for (int i=1;i<tok.size();i++)
    {
      if ( factors.find( tok[i] ) != factors.end() )
	fslot[ tok[i] ] = i;
      else
	{
	  if ( variables == NULL || variables->find( tok[i] ) != variables->end() ) 
	    vslot[ tok[i] ] = i;
	}
    }
  
  // all factors found?
  if ( fslot.size() != factors.size() ) 
    Helper::halt( "problem finding all factors in " + filename );
  
  if ( vslot.size() == 0 ) 
    Helper::halt( "no variables to import in " + filename );
  
  // extract individual with ID == 'id' only
  
  int cnt = 0 , cnt2 = 0 ;
  
  while ( ! IN1.eof() )
    {      
      std::string line;
      Helper::safe_getline( IN1 , line );      
      if ( IN1.eof() ) break;
      if ( line == "" ) continue;
      std::vector<std::string> tok = Helper::parse( line , "\t " );
      if ( tok.size() <= 2 ) Helper::halt( "problem with imported format: need at least two cols:\n" + line );      
      
      // only read for this individual ; do not assume sorted, so will have to parse all lines
      // (repeatedly).   This should not be too bad for most purposes, but if needed we can 
      // store a static version of the file in memory

      if ( tok[0] != id ) continue;
      
      //
      // build strata
      //

      std::map<std::string,std::string> curr_strata;      
      std::map<std::string,int>::const_iterator ff = fslot.begin();
      while ( ff != fslot.end() ) 
	{
	  curr_strata[ ff->first ] = tok[ ff->second ];
	  ++ff;
	}

      //
      // insert variables
      //

      std::map<std::string,int>::const_iterator vv = vslot.begin();
      while ( vv != vslot.end() )
        {
	  double x;
	  if ( Helper::str2dbl( tok[ vv->second ] , &x ) )
	    {
	      num_cache->add( ckey_t( vv->first , curr_strata ) , x );
	      ++cnt2;
	    }
          ++vv;
        }
      ++cnt;

      // next row
    }
  
  IN1.close();
  
  logger << "  read " << cnt << " strata (" << cnt2 << " distinct values) for " << id << " from " << filename << "\n";
  
}



void timeline_t::cache2annot( const param_t & param )
{
  
  // create a set of annotations based on a cache value
  
  // expecting cache value to be in seconds 'seconds'
  
  // add a window 'w'

  // window/midpoint options
  const bool add_window = param.has( "w" )  ; 
  const double w = add_window ? param.requires_dbl( "w" ) : 0 ;
  const uint64_t w_tp = w * globals::tp_1sec;


  //
  // set the annotation class
  //

  const std::string aname = param.requires( "annot" );

  annot_t * a = edf->annotations->add( aname );

  //
  // select the cache
  //
  
  std::string cache_name = param.requires( "cache" );
  
  if ( ! edf->timeline.cache.has_num( cache_name ) )
    Helper::halt( "cache not found for this individual: " + cache_name );
  
  cache_t<double> * cache = edf->timeline.cache.find_num( cache_name );
  
  //  std::cout <<" cache: " << cache->print() << "\n\n";
  
  // always expecting 'seconds'
  std::set<ckey_t> ckeys = cache->keys( "seconds" );
  
  std::set<ckey_t>::const_iterator cc = ckeys.begin();
  while ( cc != ckeys.end() )
    {

      // pull times
      std::vector<double> cx = cache->fetch( *cc );

      // channel label
      std::string ch_label = ".";
      if ( cc->stratum.find( globals::signal_strat ) != cc->stratum.end() )
	ch_label = cc->stratum.find( globals::signal_strat )->second ;

      // instance ID
      std::string inst_id = "";      
      std::map<std::string,std::string>::const_iterator ss = cc->stratum.begin();
      while ( ss != cc->stratum.end() )
	{
	  if ( ss->first != globals::signal_strat )
	    {
	      if ( inst_id != "" ) inst_id += ";"; 
	      inst_id += ss->first + "=" + ss->second ;
	    }
	  ++ss;
	}      
      if ( inst_id == "" ) inst_id = ".";

      // add annots
      for (int i=0; i<cx.size(); i++)
	{

	  uint64_t mid = cx[i] * globals::tp_1sec ;

	  interval_t interval( mid , mid );

	  // add 'w' each side
	  if ( add_window )
	    interval.expand( w_tp );

	  // add as annot
	  a->add( inst_id , interval , ch_label );
	}

      // next set of strata
      ++cc;
    }
  
}

void timeline_t::annot2cache( const param_t & param )
{
  
  // create a set of cache<int> sample "points" based on an annotation meta-data,

  // where we expect the time of each point in seconds elapsed from EDF start
  // format e.g. p=7737.4886
  
  // expected by TLOCK, etc.   need to make these globals more explicit
  const std::string cache_points_label = "points"; 
  
  // get annotations
  if ( ! param.has( "annot" ) ) Helper::halt( "no annotations specified: e.g. annot=A1,A2" );
  std::vector<std::string> anames = param.strvector( "annot" );
  
  // add channel stratifier
  const bool ignore_channel = param.yesno( "ignore-channel" );

  if ( ignore_channel ) logger << "  ignoring channel from annotations\n";
  else logger << "  tracking channel from annotations (add 'ignore-channel' to ignore)\n";
  
  // which meta-data field contains the seconds time stamp?
  std::string meta = param.requires( "meta" ); 
  logger << "  looking for tp:XXX data in annotation meta-field " << meta << "\n";
    
  // if not otherwise specified, use annot names as new channel labels, otherwise map to a single label  
  const bool single_cache = param.has( "cache" );
  const std::string single_cache_name = single_cache ? param.value( "cache" ) : "" ; 
  std::vector<std::string> cnames = anames;
  if ( single_cache ) logger << "  mapping to a single cache " << single_cache_name << "\n";
  else logger << "  mapping to caches based on annotation names\n";    
  
  // requires a sample rate to be specified
  //  - this done by 'attaching' a channel (sig) OR by specifying a SR
  //  - but we then need to find a channel with that SR (got get TP from slice)

  int slot = -1;
  int sr = -1;

  if ( param.has( "sr" ) )
    {
      sr = param.requires_int( "sr" );

      // pick first channel w/ this SR (allow for floating point wobble)
      signal_list_t signals = edf->header.signal_list( "*" );
      std::vector<double> srs = edf->header.sampling_freq( signals );
      for (int s=0;s<signals.size();s++)
	{
	  if ( fabs( srs[s] - sr ) < 0.0001 )
	    {
	      slot = s;
	      break;
	    }
	}
    }
  else
    {
      const std::string attached_sig = param.requires( "sig" );
      signal_list_t signals = edf->header.signal_list( attached_sig );
      if ( signals.size() != 1 ) Helper::halt( "expecting a single channel (present in EDF) for 'sig' " );
      
      sr = edf->header.sampling_freq( signals )[0];
      logger << "  using " << attached_sig << " (Fs = " << sr
	     << ") to anchor annotations to sample-points\n";

      slot = signals(0);
    }

  if ( sr == -1 || slot == -1 )
    Helper::halt( "need to specify sig or sr to get a sample rate (and EDF needs channel w/ that sr)");
  
  // get current time-point channel (for all)
  slice_t slice( *edf , slot , edf->timeline.wholetrace() );
  
  const std::vector<uint64_t> * tp = slice.ptimepoints();
  
  // must map to within 1 sample (i.e. if at edge?)
  const double max_diff = param.has( "diff" ) ? param.requires_dbl( "diff" ) : 1/(double)sr;
  logger << "  mapping to closest sample-point within " << max_diff << " seconds\n";
  
  // store int peaks (derived from second-level times)
  // in a channel-specific map
  std::map<std::string,std::vector<int> > d;
  
  
  struct chpt_t {
    chpt_t( const std::string & ch , uint64_t tp )
      : ch(ch) , tp(tp) { }
    
    std::string ch;
    uint64_t tp;
    
    bool operator< (const chpt_t & rhs ) const
    {
      if ( tp < rhs.tp ) return true;
      if ( tp > rhs.tp ) return false;
      return ch < rhs.ch;
    }
  };
  
  // for each annotation
  for (int a=0; a<anames.size(); a++)
    {
      
      // does annot exist?
      annot_t * annot = (*annotations)( anames[a] );
      if ( annot == NULL ) continue;
      
      int cnt = 0 ;
      
      // get all events (which are sorted)
      const annot_map_t & events = annot->interval_events;

      // tp index
      std::set<chpt_t> tps;
      
      // look at each annotation event
      annot_map_t::const_iterator aa = events.begin();
      while ( aa != events.end() )
	{
	  
	  // get instance
	  const instance_t * instance = aa->second;
	  
	  // does it have the requisite field?
	  avar_t * m = instance->find( meta );
	  
	  if ( m != NULL )
	    {
	      // get as 'tp:XXXXXXX' string
	      std::string t = m->text_value() ;

	      if ( t.size() >= 4 && t.substr(0,3) == "tp:" )
		{
		  
		  std::string t2 = t.substr(3);
		  uint64_t tpval;
		  if ( ! Helper::str2int64( t2 , &tpval ) )
		    Helper::halt( "invalid tp: format in annotation file:" + t );
		  
		  // add as channel-specific mapping or no?
		  tps.insert( chpt_t( ignore_channel ? "." : aa->first.ch_str , tpval ) );
				      
		}
	    }
	  
	  // next instance
	  ++aa;
	  
	}

      
      //
      // now we have a sorted set of time-points
      //
      
      int idx = 1;
      const int np = tp->size();
      
      // std::cout << " np = " << np << "\n";
      // std::cout << " num annots = " << tps.size() << "\n";
      
      std::set<chpt_t>::const_iterator tt = tps.begin();
      while ( tt != tps.end() )
	{
	  
	  const uint64_t & curr = tt->tp;
	  const uint64_t & prior = (*tp)[idx-1];
	  const uint64_t & next = (*tp)[idx];
  
	  // shift sample-point window up
	  if ( next < curr )
	    {
	      ++idx;
	      if ( idx == np ) break;
	      continue; // i.e. bounce back but do not update ++tt
	    }
	  
	  // is in-between these two points?
	  if ( curr >= prior && curr <= next )
	    {
	      const uint64_t d1 = curr - prior;
	      const uint64_t d2 = next - curr ;
	      
	      const bool first = d1 < d2 ; 

	      double df = ( first ? d1 : d2 ) * globals::tp_duration ; 

	      // close enough?
	      if ( df <= max_diff )
		{
		  int sp = first ? idx-1 : idx ;
		  
		  // store
		  d[ tt->ch ].push_back( sp );
		  
		  //std::cout << " adding " << anames[a] <<" " << sp << " " << ( first ? prior : next ) << "\n";
		  
		  ++cnt;
		}
	    }
	  
	  // advance to next point 
	  ++tt;
	}
      
           
      //
      // done mapping
      //
      
      logger << "  added " << cnt << " (of " << tps.size() << ") cache points from " << anames[a] << "\n";

      //
      // add to this cache, or save to add to a combined cache?
      //
      
      if ( ! single_cache )
	{
	  cache_t<int> * c = cache.find_int( cnames[a] );

	  std::map<std::string,std::vector<int> >::const_iterator dd = d.begin();
	  while ( dd != d.end() )
	    {
	      writer.level( dd->first , globals::signal_strat );
	      c->add( ckey_t( "points" , writer.faclvl() ) , dd->second );
	      ++dd;
	    }
	  writer.unlevel( globals::signal_strat );
	  d.clear();
	}
      
      // next annotation class
    }
    
  // add all points to a single cache
  if ( single_cache )
    {
      cache_t<int> * c = cache.find_int( single_cache_name );
      
      std::map<std::string,std::vector<int> >::const_iterator dd = d.begin();
      while ( dd != d.end() )
	{
	  writer.level( dd->first , globals::signal_strat );
	  c->add( ckey_t( "points" , writer.faclvl() ) , dd->second );
	  ++dd;
	}   
      writer.unlevel( globals::signal_strat );
      
    }
}
