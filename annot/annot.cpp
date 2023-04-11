
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

#include "annot.h"
#include "eval.h"
#include "edf/edf.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "defs/defs.h"
#include "tinyxml/xmlreader.h"
#include "db/db.h"
#include "nsrr-remap.h"
#include "helper/token-eval.h"

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>

extern writer_t writer;

extern logger_t logger;

extern globals global;


void annot_t::wipe()
{
  std::set<instance_t *>::iterator ii = all_instances.begin();
  while ( ii != all_instances.end() )
    {	 
      if ( *ii != NULL ) 
	delete *ii;
      ++ii;
    }    
  all_instances.clear();
}


instance_t * annot_t::add( const std::string & id , const interval_t & interval , const std::string & ch )
{

  // swap in hh:mm:ss for null instance ID?
  
  bool id2hms = 
    globals::set_annot_inst2hms_force ||
    ( globals::set_annot_inst2hms &&
      ( id == "." || id == "" || id == name ) );

  std::string id2 = id;
  
  if ( id2hms )
    {
      clocktime_t t = parent->start_ct ;
      double start_secs = interval.start_sec(); 
      t.advance_seconds( start_secs );
      id2 = t.as_string( ':' , true ); // T = include any fractional seconds
    }
  
  // get index
  instance_idx_t idx( this , interval , id2 , ch );

  // does this already exist?
  if ( interval_events.find( idx ) != interval_events.end() )
    return interval_events[ idx ];

  // else create a new instance
  instance_t * instance = new instance_t ;  

  interval_events[ instance_idx_t( this , interval , id2 , ch ) ] = instance;
  
  // track (for clean-up)
  all_instances.insert( instance );
    
  return instance; 
  
}

void annot_t::remove( const std::string & id , const interval_t & interval , const std::string & ch )
{

  instance_idx_t key = instance_idx_t( this , interval , id , ch );

  std::map<instance_idx_t,instance_t*>::iterator ii = interval_events.find( key );

  if ( ii == interval_events.end() ) return;

  // clean up instance
  if ( ii->second != NULL ) {
    
    // remove pointer from global instance tracker
    std::set<instance_t*>::iterator kk = all_instances.find( ii->second );
    if ( kk != all_instances.end() )
      all_instances.erase( kk );

    // release actual instance
    delete ii->second;
  }
  
  // clean up idx
  interval_events.erase( key );

}


std::string instance_t::print( const std::string & delim , const std::string & prelim ) const
{
  std::stringstream ss;

  std::map<std::string,avar_t*>::const_iterator dd = data.begin();
  while ( dd != data.end() )
    {
      
      if ( dd != data.begin() ) ss << delim;

      ss << prelim;
      
      if ( dd->second == NULL )
	ss << dd->first;
      else if ( dd->second->atype() == globals::A_BOOLVEC_T )
	ss << dd->first << "=" << Helper::stringize( dd->second->text_vector() , "," ) ;
      else if ( dd->second->atype() == globals::A_INTVEC_T )
	ss << dd->first << "=" << Helper::stringize( dd->second->int_vector() , "," ) ;
      else if ( dd->second->atype() == globals::A_DBLVEC_T )
	ss << dd->first << "=" << Helper::stringize( dd->second->double_vector() , "," ) ;
      else if ( dd->second->atype() == globals::A_TXTVEC_T )
	ss << dd->first << "=" << Helper::stringize( dd->second->text_vector() , "," ) ;

      else
 	ss << dd->first << "=" << dd->second->text_value();
      ++dd;

    }

  return ss.str();
}

globals::atype_t instance_t::type( const std::string & s ) const 
{
  std::map<std::string,avar_t*>::const_iterator ii = data.find( s );
  if ( ii == data.end() ) return globals::A_NULL_T;  
  return ii->second->atype();
}

void instance_t::check( const std::string & name )
{

  std::map<std::string,avar_t*>::iterator dd = data.find( name );

  if ( dd == data.end() ) return;

  if ( dd->second == NULL ) return;  // flag, so no storage set
  
  // erase actual storage...
  delete dd->second; 

  // erase place in the tracker
  std::set<avar_t*>::iterator ff = tracker.find( dd->second );
  if ( ff != tracker.end() )
    tracker.erase( tracker.find( dd->second ) ); 
  else 
    Helper::halt( "internal error in instance_t::check()... avar_t not tracked" );
  
  // and erase from this data map instance
  data.erase( dd );  
  
  return;
}

void instance_t::set( const std::string & name ) 
{
  check( name );
  avar_t * a = new flag_avar_t ;
  tracker.insert( a );
  data[ name ] = a;    
}

void instance_t::set( const std::string & name , const int i ) 
{
  check( name );
  avar_t * a = new int_avar_t( i ) ;
  tracker.insert( a );
  data[ name ] = a;    
}

void instance_t::set( const std::string & name , const std::string & s ) 
{
  check( name );
  avar_t * a = new text_avar_t( s ) ;
  tracker.insert( a );
  data[ name ] = a;    
}

void instance_t::set( const std::string & name , const bool b )
{
  check( name );
  avar_t * a = new bool_avar_t( b ) ;
  tracker.insert( a );
  data[ name ] = a;    
}

void instance_t::set_mask( const std::string & name , const bool b )
{
  check( name );
  avar_t * a = new mask_avar_t( b ) ;
  tracker.insert( a );
  data[ name ] = a;    
}

void instance_t::set( const std::string & name , const double d )
{
  check( name );
  avar_t * a = new double_avar_t( d ) ;
  tracker.insert( a );
  data[ name ] = a;    
}


// vectors
void instance_t::set( const std::string & name , const std::vector<int> &  i ) 
{
  check( name );
  avar_t * a = new intvec_avar_t( i ) ;
  tracker.insert( a );
  data[ name ] = a;    
}

void instance_t::set( const std::string & name , const std::vector<std::string> & s ) 
{
  check( name );
  avar_t * a = new textvec_avar_t( s ) ;
  tracker.insert( a );
  data[ name ] = a;    
}

void instance_t::set( const std::string & name , const std::vector<bool> & b )
{
  check( name );
  avar_t * a = new boolvec_avar_t( b ) ;
  tracker.insert( a );
  data[ name ] = a;    
}

void instance_t::set( const std::string & name , const std::vector<double> & d )
{
  check( name );
  avar_t * a = new doublevec_avar_t( d ) ;
  tracker.insert( a );
  data[ name ] = a;    
}



instance_t::~instance_t()
{
  std::set<avar_t *>::iterator ii = tracker.begin();
  while ( ii != tracker.end() )
    {	 
      delete *ii;
      ++ii;
    }
}  

std::ostream & operator<<( std::ostream & out , const avar_t & a )
{
  out << a.text_value();
  return out;
}


void summarize_annotations( edf_t & edf , param_t & param )
{
  
  writer.var( "ANNOT_N" , "Number of occurrences of an annotation" );
  
  std::map<std::string,int>::const_iterator ii = edf.aoccur.begin();
  while ( ii != edf.aoccur.end() ) 
    {
      // annot as 'level'
      writer.level( ii->first , globals::annot_strat );
      writer.value( "ANNOT_N" , ii->second );
      ++ii;
    }
}


bool annot_t::map_epoch_annotations(   edf_t & parent_edf , 
				       const std::vector<std::string> & ann , 
				       const std::string & filename , 
				       uint64_t elen , 
				       uint64_t einc )
{


  // static function that will create multiple annotations (one per class label)
  // no associated channel labels

  bool unepoched = elen == 0 ;
  
  if ( unepoched )  
    {
      elen = Helper::sec2tp( globals::default_epoch_len );
      einc = Helper::sec2tp( globals::default_epoch_len );
    }

  // get implied number of epochs
  double seconds = (uint64_t)parent_edf.header.nr * parent_edf.header.record_duration ;      
  const int ne = seconds / ( unepoched ? globals::default_epoch_len : elen / globals::tp_1sec );
  
  const int delta = abs( (int)( ne - ann.size() ) );

  if ( delta > globals::enforce_epoch_check ) 
    Helper::halt( "expecting " + Helper::int2str(ne) + " epoch annotations, but found " + Helper::int2str( (int)ann.size() ) );
  
  if ( delta != 0 )
    logger << "  ** warning: expecting " << ne << " epochs but found " << ann.size()
	   << "; will allow given epoch-check=" 
	   << globals::enforce_epoch_check << "\n";
  
  //
  // because we otherwise might have a discontinuous EDF, we need to look up the proper epoch 
  // intervals below , if epoched 
  //
  
  //
  // map of all labels/annot classes to be added
  //

  std::map<std::string,annot_t*> amap;

  for (int e=0;e<ann.size();e++) 
    {
     
      //
      // skip this annotation?
      //
      
      if ( globals::specified_annots.size() > 0 && 
	   globals::specified_annots.find( ann[e] ) == globals::specified_annots.end() ) 
	continue;

      //
      // we may have a null annot here, e.g if a whitelist was
      // specified
      //

      if ( ann[e] == "" )
	continue;
      
      //
      // ignore this annotation if past the end?
      //
      
      if ( e >= ne ) continue;
      
      //
      // otherwise, create the new annotation class
      //
      
      annot_t * a = parent_edf.timeline.annotations.add( ann[e] );

      amap[ ann[e] ] = a;
      
      a->description = ann[e];
      
      a->file = filename ;
      
      a->type = globals::A_FLAG_T;  // no meta-data from .eannot 
      
      a->types.clear();
  
    }
  
  
  //
  // Populate intervals
  //

  if ( unepoched ) 
    {

      for ( int e = 0 ; e < ann.size() ; e++ )
	{
	  
	  if ( amap.find( ann[e] ) != amap.end() )
	    {
	      
	      interval_t interval( e * elen , e * elen + einc );
	      
	      annot_t * a = amap[ ann[e] ];
	      
	      // . indicates no specifically assigned channel
	      instance_t * instance = a->add( ann[e] , interval , "." ); 
	      
	      // track how many annotations we add
	      parent_edf.aoccur[ a->name ]++;
	      
	    }
	  
	} // next epoch 
      
    }
  else
    {

      // but if we do already have an in-memory EDF, which might be
      // discontinuous, we need to use the timeline to get the 
      // proper interval for the e'th epoch
      
      parent_edf.timeline.first_epoch();
      
      std::vector<int> epoch_counts;

      int e = 0;

      while ( 1 ) 
	{

	  int epoch = parent_edf.timeline.next_epoch_ignoring_mask();
	  
	  if ( epoch == -1 ) break;
          
	  if ( e >= ann.size() ) Helper::halt( "internal error map_epoch_annot()" );

	  interval_t interval = parent_edf.timeline.epoch( epoch );
	  
	  annot_t * a = amap[ ann[e] ];
	  
	  // . indicates no assigned channel
	  instance_t * instance = a->add( ann[e] , interval , "." );
	  
	  // track how many annotations we add
	  parent_edf.aoccur[ a->name ]++;
	  
	  // next row in file
	  ++e;

	}

    }

  
  //
  // all done
  //
  
  return true;
}


bool annot_t::special() const
{
  if ( name == "duration_hms" ) return true;
  if ( name == "duration_sec" ) return true;
  if ( name == "epoch_sec" ) return true; 
  if ( name == "start_hms" ) return true; 
  return false;
}


bool annot_t::process_special( const std::string & , const std::string & )
{
  // look for special flags, and add to 
  return true;
}

bool annot_t::load( const std::string & f , edf_t & parent_edf )
{

  
  //
  // static annot_t function, which will create multiple annot_t for each new 'name'
  // encountered
  //
  
  // 
  // Check file exists and is of the correct type
  //

  if ( ! Helper::fileExists(f) ) return -1;

  if ( Helper::file_extension( f , "xml" ) ) 
    {
      Helper::halt( f + " is an XML file... should already have been loaded (internal error)" );
      return -1;
    }
  
  if ( Helper::file_extension( f , "ftr" ) ) 
    {
      Helper::halt( f + " is an FTR file... should already have been loaded (internal error)" );
      return -1;
    }
  

  //
  // A simple epoch annotation file? (based on extension)
  // These are not allowed for EDF+ files.
  //
  
  // nb. by default, file_extension() matches w/ a period
  // before the extension, so "annot" is different from "eannot"
  
  bool is_eannot = Helper::file_extension( f , "eannot" ) 
    || Helper::file_extension( f , "stages" )
    || Helper::file_extension( f , "eannot.txt" );
  
  // this also matches file.annot.txt , file.txt, file.tsv, etc
  bool is_annot = Helper::file_extension( f , "annot" )    
    || Helper::file_extension( f , "txt" ) 
    || Helper::file_extension( f , "tsv" ); 
  
  if ( is_eannot && ! parent_edf.header.continuous ) 
    Helper::halt( "cannot use .eannot files with discontinuous (EDF+) files" );


  // otherwise, need to figure this out by looking at the file?
  
  if ( ! ( is_eannot || is_annot ) )
    {
      std::ifstream IN1( f.c_str() , std::ios::in );

      while ( 1 )
	{
	  
	  std::string x;
	  Helper::safe_getline( IN1 , x );

	  if ( IN1.eof() ) break;

	  // no blank lines allowed for .eannot, so assume .annot if so
	  // which is more flexible

	  if ( x == "" )	    
	    {
	      is_eannot = false;
	      is_annot = true;
	      break;
	    }
	  
	  // start w/ a # header?
	  if ( x[0] != '#' )
	    {
	      is_eannot = false;
	      is_annot = true;
	      break;
	    }
	  
	  // number of columns? (based on defined delimiters: either tab or tab/space
	  std::vector<std::string> tok = Helper::parse( x , globals::allow_space_delim ? " \t" : "\t" );

	  if ( tok.size() > 1 )
	    {
	      is_eannot = false;
	      is_annot = true;	      
	      break;
	    }
	}

      IN1.close();

      // hmm, not sure what this file is...
      if ( is_annot == is_eannot )
	Helper::halt( "unable to determine whether " + f + " is .annot or .eannot format" );
      
    }
  


  //
  // parse as a simple .eannot file
  // 
  
  if ( is_eannot )
    {
      std::vector<std::string> a;
      
      std::ifstream IN1( f.c_str() , std::ios::in );
      while ( ! IN1.eof() )
	{
	  std::string x;
	  Helper::safe_getline( IN1 , x );
	  if ( IN1.eof() ) break;
	  if ( x == "" ) continue;
	  
	  x = Helper::unquote( x );
	  
	  // sanitize? [ done in remap now ] 
	  //	  if ( globals::sanitize_everything )
	  //  x = Helper::sanitize( x );
	  
	  // remap? (and if so, track)
	  std::string y = nsrr_t::remap( x ) ;

	  // empty annots (e.g. from white-list) will be skipped in map_epoch_annotations()
	  //if ( y == "" ) continue;
	  
	  if ( y != x && y != "" ) parent_edf.timeline.annotations.aliasing[ y ] = x ;
	  
	  // store
	  a.push_back( y );
	}
      IN1.close();

      annot_t::map_epoch_annotations( parent_edf , 
				      a , 
				      f , 
				      parent_edf.timeline.epoch_len_tp() , 
				      parent_edf.timeline.epoch_increment_tp() );
      
      
      return true;
      
    }


  
  
  //
  // Otherwise, this is an .annot file   
  //

  std::ifstream FIN( f.c_str() , std::ios::in );

  // header with # character
  
  // types: int, str, dbl, bool
  // [txt] [str] 
  // [dbl] [num]
  // [int] 
  // [yn] [bool]
  // [.] none
  
  // # name1 | description | col1(int) col2(str) col3(dbl) 
  // # name2 | description | col1(str)
  // # name3 | description                              [ just means a bool ] 
  
  // then rows are either interval or epoch-based, but must **always** have 6 tab-delim columns 
  // exception: if we find 4 (old format, make channel and meta missing)

  // name  id1  ch  sec1  sec2  { vars }
  // name  id   .   e:1   {e:2} { vars } 
  // name  id   ch  hh:mm:ss  hh:mm:ss { vars }
  // special .  .   . . .  // i.e. no time specified... store as 0/0 

  // assume e:1   means 30-second epoch, no overlap  [ hard-code this ] 
  // e:1:20 is first of a 20-second epoch
  // e:1:20:10 is similar, but w/ epoch overlap of 10 seconds (10=increment)

  // if second time-point is '...' this means go until the next annot (or end of EDF)
  // (note: a single '.' can be a missing variable)

  // align certain epochs to the start of each second (on the assumption that
  // EDF record size == 1 seconds;   i.e. this is to shift staging annotations forward by <1 sec
  // to get an aligned EDF;  we can then trim unstaged/unepoched data at the end of the recordig

  const bool align_annots = globals::annot_alignment.size() > 0 ;
    
  // check EDF starttime, which might be needed
  //  -- but add EDF start date here too, to allow dhms printing
  
  clocktime_t startdatetime( parent_edf.header.startdate, parent_edf.header.starttime );
  clocktime_t starttime( parent_edf.header.starttime );
    
  // read header then data
  
  bool epoched = false;

  bool has_class_class = false;

  int line_count = 0;

  std::map<std::string,annot_t*> annot_map;
  
  std::map<annot_t*,std::vector<std::string> > cols;

  // to allow '...' in the second line, we need to read ahead
  // and so store the read line here
  std::string buffer = "";
  
  while ( ! FIN.eof() )
    {
      
      if ( FIN.bad() ) continue;
      std::string line;      

      // read from buffer or disk?
      if ( buffer != "" )
	{
	  line = buffer;
	  // now clear buffer
	  buffer = "";
	}
      else // get fresh line from disk
	Helper::safe_getline( FIN , line );      

      if ( FIN.eof() || line == "" ) continue;

      
      //
      // header or data row? , or type header (optionally, this is skipped)  
      //

      if ( line[0] == '#' ) 
	{
	  
	  // skip initial # here
	  // quoted parse, in case annotation class has a `|` character
	  std::vector<std::string> tok = Helper::quoted_parse( line.substr(1) , "|" );
	  
	  if ( tok.size() < 1 || tok.size() > 3 )
	    Helper::halt( "bad header for format\n" + line );

	  //
	  // Get the name and ID 
	  //

	  std::string orig_name = Helper::unquote( Helper::trim( tok[0] ) ) ;
	  
	  // by default, '.' is the class.inst delimiter and this is
	  // not sanitized;  so we can skip this here, as remap() handles
	  // sanitization (but respects keeping spaces, if requested)
	  // so, comment out this step below
	  // want to keep '.' symbol here though... 
	  // if ( globals::sanitize_everything )
	  //   orig_name = Helper::sanitize( orig_name , globals::class_inst_delimiter );
	  
	  std::string name = nsrr_t::remap( orig_name );
	  
	  if ( name == "" ) continue;
	  
	  //
	  // track any aliasing (pre-split class/inst)
	  //
	  
	  if ( name != orig_name )
	    {
	      //std::cout << "alias mapping " << name << " --> " << orig_name << "\n";
	      parent_edf.timeline.annotations.aliasing[ name ] = orig_name ;
	    }

	  //
	  // If, post-remapping, this is a clas/inst split, only record the class annot
	  //

	  if ( name.find( globals::class_inst_delimiter ) != std::string::npos )
	    name = name.substr( 0 , name.find( globals::class_inst_delimiter ) );
	  
	  
	  //
	  //
	  // skip this annotation
	  //

	  if ( globals::specified_annots.size() > 0 && 
	       globals::specified_annots.find( name ) == globals::specified_annots.end() ) 
	    continue;

	  //std::cout << " loading " << name << "\n";
	  
	  //
	  // otherwise, create the annotation if it doesn't exist
	  //

	  annot_t * a = parent_edf.timeline.annotations.add( name );


	  //
	  // store a temporary lookup table
	  //

	  annot_map[ name ] = a;
	  
	  //
	  // special case: if an explicit class label is 'class', then do not
	  // allow a header row starting 'class' also
	  //

	  if ( name == "class" ) has_class_class = true;
	  
	  //
	  // Other details
	  //
		  
	  a->description = tok.size() >= 2 ? Helper::trim( tok[1] ) : name ; 
	  
	  a->file = f;
	  
	  a->type = globals::A_FLAG_T; // unless we learn otherwise, i.e. just below when parsing header line

	  a->types.clear();
	  
	  // columns specified
	  if ( tok.size() == 3 ) 
	    {
	      std::vector<std::string> type_tok = Helper::parse( tok[2] , " \t" );
	      
	      for (int j=0;j<type_tok.size();j++)
		{
		  std::vector<std::string> type_tok2 = Helper::parse( type_tok[j] , "[(" );
		  if ( type_tok2.size() > 2 ) Helper::halt( "bad type '" + type_tok[j] + "'" );
		  
		  // column name
		  const std::string var_name = type_tok2[0] ;
		  
		  // track order of columns for this annotation 
		  cols[a].push_back( var_name );

		  // type		  
		  globals::atype_t t = globals::A_NULL_T;
		  
		  if ( type_tok2.size() == 1 ) 
		    {
		      // this is commented out to mean that we now force a type 
		      // specification, i.e. what does 'FLAG' mean as meta-data?
		      // t = globals::A_FLAG_T;
		    }
		  else
		    {
		      char c = type_tok2[1][ type_tok2[1].size() - 1 ] ;
		      if ( c != ']' && c != ')' ) Helper::halt( "bad type '" + type_tok[j] + "' -> " + c );
		      std::string tstr  = type_tok2[1].substr( 0 , type_tok2[1].size() - 1 );
		      
		      if ( globals::name_type.find( tstr ) != globals::name_type.end() )
			t = globals::name_type[ tstr ];

		    }
		  
		  if ( t == globals::A_NULL_T )
		    Helper::halt( "unsupported annotation type from\n" + line );

		  a->types[ var_name ] = t ; 
		
		  // if only a single TYPE has been specified, assign this to the 
		  // annotation class, otherwise set it as undefined
		  
		  if ( type_tok.size() == 1 ) 
		    a->type = t ; 
		  else // i.e. instead of 'FLAG', this means that we have multiple types
		    a->type = globals::A_NULL_T ; 
		  
		} // next column for this annotation
	      
	    }
	  
	}
      
      

      //
      // optional case-insensitive header (e.g. for stats package, that is skipped); if an annotation 
      // named 'class' was defined above, we do not allow this header row
      //

      else if ( line_count == 0 && (! has_class_class) && line.size() > 5 && Helper::iequals( line.substr(0,5) , "class" ) )
	{
	  std::vector<std::string> tok = Helper::parse( line , globals::allow_space_delim ? " \t" : "\t" );
	  if ( tok.size() ==6 ) 
	    {
	      if ( ! Helper::iequals( tok[0] , "class" ) ) Helper::halt( "expecting column 1 to be 'class':\n" + line );
	      if ( ! Helper::iequals( tok[1] , "instance" ) ) Helper::halt( "expecting column 2 to be 'instance':\n" + line );
	      if ( ! Helper::iequals( tok[2] , "channel" ) ) Helper::halt( "expecting column 3 to be 'channel':\n" + line );
	      if ( ! Helper::iequals( tok[3] , "start" ) ) Helper::halt( "expecting column 4 to be 'start':\n" + line );
	      if ( ! Helper::iequals( tok[4] , "stop" ) ) Helper::halt( "expecting column 5 to be 'stop':\n" + line );
	      if ( ! Helper::iequals( tok[5] , "meta" ) ) Helper::halt( "expecting column 6 to be 'meta':\n" + line );
	    } 
	  else if ( tok.size() == 4 ) // exception: still allowing old .annot format
	    {
	      if ( ! Helper::iequals( tok[0] , "class" ) ) Helper::halt( "expecting column 1 to be 'class':\n" + line );
	      if ( ! Helper::iequals( tok[1] , "instance" ) ) Helper::halt( "expecting column 2 to be 'instance':\n" + line );
	      if ( ! Helper::iequals( tok[2] , "start" ) ) Helper::halt( "expecting column 3 to be 'start':\n" + line );
	      if ( ! Helper::iequals( tok[3] , "stop" ) ) Helper::halt( "expecting column 4 to be 'stop':\n" + line );	      
	    }
	  else if ( tok.size() == 3 ) // exception: still allowing old .annot format minus instance
	    {
	      if ( ! Helper::iequals( tok[0] , "class" ) ) Helper::halt( "expecting column 1 to be 'class':\n" + line );
	      if ( ! Helper::iequals( tok[1] , "start" ) ) Helper::halt( "expecting column 2 to be 'start':\n" + line );
	      if ( ! Helper::iequals( tok[2] , "stop" ) ) Helper::halt( "expecting column 3 to be 'stop':\n" + line );	      
	    }
	  else
	    Helper::halt( "invalid header line:\n" + line );
	  
	}

      //
      // otherwise, assume this is a data row
      //

      else 
	{
	  
	  // data-rows are tab-delimited typically, but optionally allow spaces; also allow these to be quoted
	  
	  std::vector<std::string> tok = globals::allow_space_delim ? 
	    Helper::quoted_parse( line , " \t" ) :
	    Helper::parse( line , "\t" );
	  
	  if ( tok.size() == 0 ) continue; 
	  
	  if ( tok.size() == 1 ) 
	    Helper::halt( "invalid data line:\n" + line + "\n (hint: use the 'tab-only' option to ignore space delimiters)" );
	  

	  //
	  // Sanitize class name, but keep '/' and '.' symbol 
	  //
	  
	  // std::string aname = globals::sanitize_everything
	  //   ? Helper::sanitize( Helper::unquote( tok[0] ) , globals::class_inst_delimiter )
	  //   : Helper::unquote( tok[0] );

	  std::string aname = Helper::unquote( tok[0] );
	  
	  //
	  // Remap (and sanitize) term?
	  //
	  
	  std::string tname = nsrr_t::remap( aname );
	  
	  if ( tname == "" ) continue;
	  
	  if ( tname != aname )
	    parent_edf.timeline.annotations.aliasing[ tname ] = aname;

	  aname = tname;
	  

	  //
	  // save original class name (prior to any combining)
	  // as this is what any in-file header information is based on for meta-data
	  //

	  const std::string cls_root = aname;

	  //
	  // Instance label
	  //

	  std::string iname = tok[1];


	  //
	  // Sanitize instance ID?
	  //

	  if ( globals::sanitize_everything )
	    iname = Helper::sanitize( iname );
	  
	  //
	  // Combine class & instance ID? (unless missing, or same as class)
	  //
	  
	  if ( globals::combine_annot_class_inst && iname != "." && iname != aname )
	    aname += globals::annot_class_inst_combiner + iname ;
	  
	  //
	  // Is this an aggregate class/inst form?
	  //

	  const bool split_annot = aname.find( globals::class_inst_delimiter ) != std::string::npos;
	  
	  std::string new_inst_id = ".";
	  
	  if ( split_annot ) 
	    {
	      // old : class=A/B inst=X
	      // new : class=A   inst=B    meta:inst=X
	      // if original inst is null, ignore
	      
	      std::vector<std::string> toks =
		Helper::parse( aname , std::string( 1 , globals::class_inst_delimiter ) );

	      if ( toks.size() != 2 ) Helper::halt( "bad format for class/inst pairing: " + aname ); 

	      // update class ID now; update meta-data (if needed) below
	      // any exclusions are based on that / also any meta-data look up 
	      aname = toks[0];
	      
	      // have to save inst ID, as might not be slot[1], e.g. if 3-col format
	      // so switch that in below also
	      new_inst_id = toks[1];
	    }


	  //
	  // are we skipping this annotation anyway?
	  //
	  
	  if ( globals::specified_annots.size() > 0 && 
	       globals::specified_annots.find( aname ) == globals::specified_annots.end() ) 
	    continue;

	  //
	  // was this annotation specified in the header?  If not, create on-the-fly as
	  // a vanilla annot (i.e. same as '# annot' in header row, so implies no 
	  // meta-data
	  //

	  std::map<std::string,annot_t*>::iterator aa = annot_map.find( aname );
	  
	  if ( aa == annot_map.end() ) 
	    {
	      // instead of complaining, we now create a new annot, copying over
	      // any existing header info if availanle from the original class name

	      std::map<std::string,annot_t*>::iterator oo = annot_map.find( cls_root );
	      const bool has_original = oo != annot_map.end();

	      // make a new annot_t class
	      annot_t * a = parent_edf.timeline.annotations.add( aname );
	      
	      // add to the current map
	      annot_map[ aname ] = a;

	      // copy description, or set to new class name
	      a->description = has_original ? oo->second->description : aname;

	      // set file code
	      a->file = f;
	      a->type = globals::A_FLAG_T; 

	      // copy meta-data over
	      if ( has_original ) 
		a->types = oo->second->types;
	      else
		a->types.clear();

	      // update local tracker of meta-col order
	      cols[ a ] = cols[ oo->second ] ;
	      
	      // and now update this iterator
	      aa = annot_map.find( aname );

	      
	    }

	  //
	  // track any aliasing and/or sanitizing
	  // nb. this will also flag split/combining instances.. that's okay
	  //
	  
	  if ( cls_root  != tok[0] )
	    parent_edf.timeline.annotations.aliasing[ cls_root ] = tok[0];
	  

	  //
	  // Check size
	  //

	  if ( tok.size() != 6 ) 
	    {
	      
	      // exception #1 : allow old 4-col formatting
	      //   class inst start stop
	      if ( tok.size() == 4 ) 
		{
		  tok.resize( 6 );
		  // 0  1  2  3  4  5
		  // cl in ch bg ed mt
		  // cl in    bg ed    	  
		  tok[5] = ".";
		  tok[4] = tok[3];
		  tok[3] = tok[2];
		  tok[2] = ".";		  
		}
	      else if (  tok.size() == 3 )
		{
		  // exception #2: 3-col format
		  //   class start stop		  
		  tok.resize( 6 );
		  // 0  1  2  3  4  5
		  // cl in ch bg ed mt
		  // cl       bg ed		  
                  tok[5] = ".";
                  tok[4] = tok[2];
                  tok[3] = tok[1];
                  tok[2] = ".";
		  tok[1] = ".";
		}
	      else
		Helper::halt ( "expecting 6/4/3 columns, but found " 
			       + Helper::int2str( (int) tok.size() ) 
			       + "\n  (hint: use the 'tab-only' option to ignore space delimiters)\n"
			       + "line [ " + line + "]" );
	    }

	  //
	  // Update instance ID if needed
	  //

	  std::string original_inst_id = ".";
		  
	  if ( split_annot && new_inst_id != "." ) 
	    {
	      // save any existing instance ID (--> meta data, below)
	      original_inst_id = iname;
	      
	      // update actual instance ID
	      iname = new_inst_id;
	    }

	  //
	  // Get annot_t pointer 
	  //

	  annot_t * a = aa->second; 			


	  //
	  // Instance ID (unless already combined w/ class)
	  //

	  const std::string id = globals::combine_annot_class_inst ? "." : iname ;
	  
	  
	  //
	  // Get interval implied 
	  //

	  // if the second time-point is '...', i.e. in which case we
	  // need to read the next line (and store in buffer for the
	  // subsequent row)
	  
	  bool readon = false;
	  
	  std::string ch; 

	  interval_t interval = get_interval( line , tok ,
					      &ch , 
					      &readon , 
					      parent_edf , a ,
					      starttime , startdatetime, 
					      f ,
					      align_annots );
	  
	  //
	  // Check this isn't prior to EDF start, i.e. can happen
	  // if hh:mm:ss format is given, but EDF start if after 
	  // that time;  in this case, get_interval() returns a special code
	  // of 
	  //
	  
	  if ( interval.start == 1 && interval.stop == 0 ) 
	    {
	      //logger << "  *** warning, skipping annot\n";
	      continue;
	      
	    }
	  
	  	  
	  //
	  // a single time-point
	  //

	  if ( readon )
	    {
	      
	      if ( FIN.bad() )
		Helper::halt( "problem reading from " + f );

	      // read into the read-ahead buffer
	      Helper::safe_getline( FIN , buffer );
	      
	      // was this the last line?
	      if ( FIN.eof() )
		{		  
		  // setting stop to 1 time-point past the last accessible
		  // time-point in the EDF/EDF+
		  interval.stop = parent_edf.timeline.last_time_point_tp + 1LLU;		  
		}
	      else
		{
		  std::vector<std::string> ntok = Helper::parse( buffer , globals::allow_space_delim ? " \t" : "\t" );

		  if ( ntok.size() == 0 ) 
		    Helper::halt( "invalid line following '...' end timepoint" );

		  // allow diff formats
		  if ( ntok.size() != 6 ) 
		    {
		      if ( ntok.size() == 4 ) 
			{
			  ntok.resize( 6 );
			  ntok[5] = ".";
			  ntok[4] = ntok[3];
			  ntok[3] = ntok[2];
			  ntok[2] = ".";			  			  
			}
		      else if (  ntok.size() == 3 )
			{
			  ntok.resize( 6 );
			  ntok[5] = ".";
			  ntok[4] = ntok[2];
			  ntok[3] = ntok[1];
			  ntok[2] = ".";
			  ntok[1] = ".";
			}
		      else
			Helper::halt ( "expecting 6/4/3 columns, but found " 
				       + Helper::int2str( (int) ntok.size() ) 
				       + "\n  (hint: use the 'tab-only' option to ignore space delimiters)\n"
				       + "line [ " + buffer + "]" );
		      
		    }
		  
		  std::string nch;
		  bool dummy;
		  
		  // only need to get the start of the interval here:
		  interval_t ninterval = get_interval( line , ntok ,
						       &nch, 
						       &dummy, 
						       parent_edf , NULL ,
						       starttime , startdatetime, 
						       f ,
						       align_annots );
		  
		  
		  // Check for valid ordering: i.e. start of next cannot be before start of prior
		  
		  if ( interval.start >= ninterval.start )
		    Helper::halt( "invalid '...' interval, next line starts too soon: \n" 
				  + line + "\n"
				  + buffer + "\n"
				  + Helper::int2str( interval.start ) + " >= " 
				  + Helper::int2str( ninterval.start ) );

		  // update with start of next annotation
		  // (i.e. as stop is encoded as +1 end, this will
		  //  imply we go right up until the point just before the
		  //  next epoch, which is what we want
		  
		  interval.stop = ninterval.start;		  

		}

	    }

	  
	  //
	  // for a FLAG for this annotation (with same name as primary annotaton)
	  //

	  instance_t * instance = a->add( id , interval , ch );

	  
	  //
	  // track how many annotations we add
	  //
	  
	  parent_edf.aoccur[ a->name ]++;

	  
	  //
	  // Also add any other columns (as separate events under this same annotation)
	  //

	  const int n = tok.size();
	  

	  //
	  // We've now fixed to 6 columns
	  //

	  if ( n != 6 ) 
	    Helper::halt( f + " has a non-blank row with other than 6 tab-delimited fields:\n" + line );

	  //
	  // Special case: if we split a class/inst ID, put any old
	  // instance ID info into a new meta-field 'inst'
	  //

	  if ( split_annot && original_inst_id != "." )
	    {
	      // ensure this special type is added
	      a->types[ "_inst" ] = globals::A_TXT_T ;
	      instance->set( "_inst" , original_inst_id );
	    }


	  //
	  // If var column is '.' we can skip ahead to next line now (or if
	  // we are not expecting any variables)
	  //
	  
	  if ( tok[5] == "." ) continue;

	  // relax this assumption for now, i.e. if we let key=value
	  // pairs be defined on the command line
	  //if ( cols[a].size() == 0 ) continue;
	  
	  //
	  // Otherwise, parse |-delimited values, that should match the header
	  //   - nb. still allow for quoted `|` characters in meta-data
	  //
	  
	  std::vector<std::string> vartok = Helper::quoted_parse( tok[5] , "|" );
	  
	  const int nobs = vartok.size();
	  const int nexp = cols[a].size();
	  
	  //
	  // Are meta-data in key=value pair mode? Look at the first |- delimited element, do it contain key=value
	  // ( Can escape this by quoting the element ) 
	  //
	  
	  bool key_value = Helper::quoted_parse( vartok[0] , std::string(1,globals::annot_keyval_delim ) ).size() == 2 ;

	  // vartok[0][0] != '"' && vartok[0].find( globals::annot_keyval_delim ) != std::string::npos; 

	  //
	  // need at least this many fields, i.e. if we've pre-specifed above	  
	  //
	  
	  if ( nobs > nexp && ! key_value )
	    Helper::halt( "expecting at most "
			  + Helper::int2str( nexp )
			  + " |-delimited fields for " + aname + "\n" + line );
	  
	  // 
	  // Read expected fields, with specified types
	  //
	  
	  for (int j=0; j<nobs; j++)
	    {
	      
	      // skip missing values
	      if ( vartok[j] == "." ) continue;
	      
	      // key=value pair?
	      std::vector<std::string> kv;
	      if ( key_value ) 
		{
		  kv = Helper::quoted_parse( vartok[j] , std::string(1,globals::annot_keyval_delim ) );
		  
		  if ( kv.size() != 2 ) 
		    Helper::halt( "expecting key"
				  + std::string(1,globals::annot_keyval_delim)
				  + "value pair: " + vartok[j] );
		}

	      // get label
	      const std::string & label = key_value ? kv[0] : cols[a][j];	      
	      
	      // if this key not declare or previously seen, add now as TXT
	      if ( key_value && a->types.find( label ) == a->types.end() )
		{
		  a->types[ label ] = globals::A_TXT_T;
		  //Helper::halt( "could not read undefined type from annotation file for " + label + "\n" + line );
		}
	      
	      // get type
	      globals::atype_t t = a->types[label];
	      
	      if ( t == globals::A_MASK_T )
		{
		  // accepts F and T as well as long forms (false, true)
		  bool value = Helper::yesno( key_value ? kv[1] : vartok[j] );
		  instance->set_mask( label , value );		
		}
	      
	      else if ( t == globals::A_BOOL_T )
		{
		  // accepts F and T as well as long forms (false, true)
		  bool value = Helper::yesno( key_value ? kv[1] : vartok[j] );
		  instance->set( label , value );
		}

	      else if ( t == globals::A_INT_T )
		{
		  int value = 0;
		  if ( ! Helper::str2int( key_value ? kv[1] : vartok[j] , &value ) )
		    Helper::halt( "invalid E line, bad numeric value" );
		  instance->set( label , value );
		}

	      else if ( t == globals::A_DBL_T )
		{
		  double value = 0;
		  
		  if ( Helper::str2dbl( key_value ? kv[1] : vartok[j] , &value ) )
		    instance->set( label , value );
		  else
		    if ( vartok[j] != "NA" )
		      Helper::halt( "invalid line, bad numeric value:\n" + line );
		}

	      else if ( t == globals::A_TXT_T )
		{		  
		  instance->set( label , key_value ? kv[1] : Helper::unquote( vartok[j] ) );
		}
	      
	      //
	      // TODO.. add vector readers; for now, they can be encoded as, e.g. comma-delimited strings
	      //
	      
	      else
		logger << "could not read undefined type from annotation file for " << label << "\n";
	      
	      //
	      // next value
	      //
	    }


	  //
	  // Done processing this line
	  //

	  ++line_count;
	}
      
    } // next line

  FIN.close();
  
  return line_count;
}



interval_t annot_t::get_interval( const std::string & line ,
				  const std::vector<std::string> & tok ,
				  std::string * ch , 
				  bool * readon , 
				  const edf_t & parent_edf , 
				  annot_t * a ,
				  const clocktime_t & starttime ,
				  const clocktime_t & startdatetime , 
				  const std::string & f ,
				  const bool align_annots 
				  )
{

  // std::cout << "[" << line << "]\n";
  
  // 0 class
  // 1 instance
  // 2 channel
  // 3 start
  // 4 stop
  // 5 meta

  // epoch (single or range) or an interval (range)? 
  
  if ( tok.size() != 6 ) 
    Helper::halt( "bad line format, need exactly 6 columns:\n" + line );
  
  // specification in terms of one (or two) epochs?  
  bool eline  = tok[3][0] == 'e' ;
  bool eline2 = tok[4][0] == 'e' ;

  if ( eline2 && ! eline ) 
    Helper::halt( "not a valid epoch row if only second field has e:N encoding");
  
  // if the second timepoint is '...' or '-' then this will be set to true
  // i.e. indicating that we need to look at the next record in order to
  // figure out the end of this annotation
  
  *readon = tok[4] == "..." || tok[4] == "-";
  
  //
  // Get channel label
  //
  
  *ch = globals::sanitize_everything ? Helper::sanitize( tok[2] ) : tok[2] ;


  //
  // Get interval
  //
  
  interval_t interval;

  
  // if hh:mm:ss is given that is before the EDF start::
  //  assign a special code to the returned interval, so that is can be skipped
  //   interval_t(1,0) 
          
  bool before_edf_start = false;

  //
  // for special annotations, can use thisL just insert at start of recording
  //

  if ( tok[3] == "." && tok[4] == "." )
    {  
      interval.start = 0LLU;
      interval.stop = 0LLU;
    } 
  else if ( eline )
    {
      
      if ( parent_edf.header.edfplus || ! parent_edf.header.continuous ) 
	Helper::halt( "cannot use e:1 notation in .annot files with (discontinuous) EDF+ files" );
            
      //  e:1        assumes 30
      //  e:30:1     assumes no overlap
      //  e:15:5:4   specifies all
      
      std::vector<std::string> tok2 = Helper::parse( tok[3] , ":" );
      
      if ( tok2.size() < 2 || tok2.size() > 4 ) 
	Helper::halt( "bad epoch specification, expecting e:1, e:30:1, e:30:30:1, etc" );
      
      if ( tok2[0] != "e" ) 
	Helper::halt( "bad epoch specification, expecting e:1, e:30:1, e:30:30:1, etc" );	    
      
      int epoch_length = globals::default_epoch_len;
      int epoch_increment = globals::default_epoch_len; // i.e. non-overlapping
      int epoch;
      
      if ( ! Helper::str2int( tok2[ tok2.size() - 1 ] , &epoch ) ) 
	Helper::halt( "invalid epoch: " + tok[2] );
      
      if ( epoch == 0 ) 
	Helper::halt( "invalid E value of '0' (first epoch should be '1')" );
      
      if ( tok2.size() >= 3 ) 
	if ( ! Helper::str2int( tok2[ 1 ] , &epoch_length ) ) 
	  Helper::halt( "invalid epoch length:  " + tok[2] );
      
      if ( tok2.size() == 4  ) 
	if ( ! Helper::str2int( tok2[ 1 ] , &epoch_increment ) ) 
	  Helper::halt( "invalid epoch increment:  " + tok[2] );
      
      uint64_t epoch_length_tp = Helper::sec2tp( epoch_length );
      
      uint64_t epoch_increment_tp = Helper::sec2tp( epoch_increment );
      
      // set interval from current line
      // last point is defined as point *past* end of interval
      
      interval.start = epoch_increment_tp * (epoch-1);
      interval.stop  = interval.start + epoch_length_tp;
      
      //
      // A second epoch?   in this case, overwrite interval.stop from above
      //
      
      if ( eline2 )
	{
	  std::vector<std::string> tok2 = Helper::parse( tok[4] , ":" );
	  
	  if ( tok2.size() < 2 || tok2.size() > 4 ) 
	    Helper::halt( "bad epoch specification, expecting e:1, e:30:1, e:30:30:1, etc" );
	  
	  if ( tok2[0] != "e" ) 
	    Helper::halt( "bad epoch specification, expecting e:1, e:30:1, e:30:30:1, etc" );	    
	  
	  int epoch;
	  
	  if ( ! Helper::str2int( tok2[ tok2.size() - 1 ] , &epoch ) ) 
	    Helper::halt( "invalid epoch: " + tok[2] );
	  
	  if ( epoch == 0 ) 
	    Helper::halt( "invalid E value of '0' (first epoch should be '1')" );
	  
	  if ( tok2.size() >= 3 ) 
	    if ( ! Helper::str2int( tok2[ 1 ] , &epoch_length ) ) 
	      Helper::halt( "invalid epoch length:  " + tok[2] );
	  
	  if ( tok2.size() == 4  ) 
	    if ( ! Helper::str2int( tok2[ 1 ] , &epoch_increment ) ) 
	      Helper::halt( "invalid epoch increment:  " + tok[2] );
	  
	  uint64_t epoch_length_tp = Helper::sec2tp( epoch_length );
	  
	  uint64_t epoch_increment_tp = Helper::sec2tp( epoch_increment );
	  	  
	  // set interval from current line
	  // last point is defined as point *past* end of interval
	      
	  uint64_t start_of_last_epoch = epoch_increment_tp * (epoch-1);
	  interval.stop  = start_of_last_epoch + epoch_length_tp;
	}
      
    }
  else // an INTERVAL (i.e. starting with other than an e:N value)
    {

      // assume this is either: 
      //    single numeric value (in seconds) which is an offset past the EDF start
      // OR in clock-time, in hh:mm:ss (24-hour) format or dd-mm-yy-hh:mm:ss
      // OR in elapsed clock-time in format 0+hh:mm:ss
      
      // with either format, the second column can be a read-on  ('...')
      // with either format, the second column can be a duration (in secs) if second col starts +
      
      bool col2dur = tok[4][0] == '+';

      //
      // Check: are these times? hh:mm:ss or dd:hh:mm:ss 
      //  1) allow microseconds to be specified
      //  2) allow optional dd prefix 
      //  3) assume 24-hour clock time by default ...      
      //  4) ... unless it starts with 0+00:01:20, in which case it means elapsed
      //  5) ... or unless it has a AM/PM modifier
      
      bool is_elapsed_hhmmss_start = tok[3].size() > 2 && tok[3][0] == '0' && tok[3][1] == '+';
      bool is_elapsed_hhmmss_stop  = tok[4].size() > 2 && tok[4][0] == '0' && tok[4][1] == '+';

      std::string start_str = is_elapsed_hhmmss_start ? 
	tok[3].substr(2) : tok[3] ;

      std::string stop_str = is_elapsed_hhmmss_stop ? 
	tok[4].substr(2) : tok[4] ;
      
      // get :-delimited hh:mm:ss values ( potentially stripping [ and ] from start/end
      std::vector<std::string> tok_start_hms = Helper::parse( start_str , ":" );

      std::vector<std::string> tok_stop_hms;
      if ( ! ( *readon || col2dur ) ) 
	tok_stop_hms = Helper::parse( stop_str  , ":" );
      
      // does this look like hh:mm:ss or dd:hh:mm:ss?   (nb can be hh:mm:ss.ssss) 
      bool is_hms1 = tok_start_hms.size() == 3 || tok_start_hms.size() == 4; 
      bool is_hms2 = ( *readon || col2dur ) ? false : ( tok_stop_hms.size() == 3 || tok_stop_hms.size() == 4 );

      // if so, check that there is a valid EDF header starttime
      if ( is_hms1 && (!is_elapsed_hhmmss_start) && ( !starttime.valid) )
	Helper::halt( "specifying hh:mm:ss clocktime start, but no valid EDF header starttime" );
      
      if ( is_hms2 && (!is_elapsed_hhmmss_stop) && ( !starttime.valid) )
	Helper::halt( "specifying hh:mm:ss clocktime stop, but no valid EDF header starttime" );

    
      // convert to / read as seconds 
      
      double dbl_start = 0 , dbl_stop = 0;
      
      // start time
      
      if ( is_hms1 )
	{

	  clocktime_t atime( start_str );

	  // 0+hh:mm:ss format
	  if ( is_elapsed_hhmmss_start )
	    {

	      if ( atime.d != 0 )
		Helper::halt( "elapsed clock-times cannot contain dates: format = 0+hh:mm:ss" );

	      // i.e. seconds past 'midnight == start of EDF' 
	      dbl_start = atime.seconds();
	      
	    }
	  else
	    {

	      // if dates are specified, check that annot does not start before the EDF start
	      // otherwose, *assume* that it starts afterwards
	      //  i.e. if start = 10pm,  then 9pm --> 23 hours later, assumed the next day
	      
	      // if start is before EDF start, flag that ( special flag interval interval_t(1,0)   

	      // day information specified?
	      
	      if ( startdatetime.d != 0 && atime.d != 0 ) 
		{		  
		  int earlier = clocktime_t::earlier( startdatetime , atime );		  
		  if ( earlier == 2 )
		    before_edf_start = true;
		  else
		    dbl_start = clocktime_t::ordered_difference_seconds( startdatetime , atime ) ;
		}
	      else if ( startdatetime.d == 0 && atime.d != 0 )
		{
		  // do not allow date info in annot if EDF start is null
		  Helper::halt( "cannot specify annotations with date-times if the EDF start date is null (1.1.85)" );
		}
	      else
		{
		  
		  // otherwise, no date information for the annotation, so
		  //  a) ignore date of EDF start and
		  //  b) assume that the time is the next to occur		   
		  
		  dbl_start = clocktime_t::ordered_difference_seconds( starttime , atime ) ;
		  
		}
	      	      
	    }

	}
      else 
	{
	  // if here, we are assuming this is not a (dd-mm-yy-)hh:mm:ss format time, 
	  // so assume this is seconds 

	  if ( ! Helper::str2dbl( start_str , &dbl_start ) )
	    Helper::halt( "invalid interval (start) : " + line );
	  
	}

      
      // stop time:
      if ( is_hms2 )
	{
	  
	  clocktime_t btime( stop_str );
	  
	  if ( is_elapsed_hhmmss_stop )
	    {
	      // was elapsed 0+hh:mm:ss
	      
	      if ( btime.d != 0 )
                Helper::halt( "elapsed clock-times cannot contain dates: format = 0+hh:mm:ss" );
	      
	      dbl_stop = btime.seconds(); 
	    }
	  else
	    {
	      // date-time available for stop and EDF start?

              if ( startdatetime.d != 0 && btime.d != 0 )
                {
		  
                  int earlier = clocktime_t::earlier( startdatetime , btime );
		  
		  if ( earlier == 2 )
	            before_edf_start = true;
                  else
                    dbl_stop = clocktime_t::ordered_difference_seconds( startdatetime , btime ) ;		  
		}
	      else if ( startdatetime.d == 0 && btime.d != 0 )
		{
		  Helper::halt( "cannot specify annotations with date-times if the EDF start date is null (1.1.85)" );
		}
	      else
		{
		  // otherwise, no date information for the annotation, so
		  //  a) ignore date of EDF start and
		  //  b) assume that the time is the next to occur		   
		  
		  dbl_stop = clocktime_t::ordered_difference_seconds( starttime , btime ) ;
		  
		}
	      
	    }
	  
	}
      else if ( col2dur ) // expecting ""
	{
	  // if a + if stop column, ALWAYS has to be in seconds 
	  double dur = 0;
	  if ( ! Helper::str2dbl( tok[4].substr(1) , &dur ) ) // skip '+' not that it matters
	    Helper::halt( "could not parse stop time for line:\n" + line );
	  dbl_stop = dbl_start + dur;
	}
      else if ( ! *readon )
	{	  
	  if ( ! Helper::str2dbl( tok[4] , &dbl_stop ) )
	    Helper::halt( "invalid interval (stop): " + line );
	}
      
      if ( dbl_start < 0 )
	{
	  //std::cout << " S1 " << dbl_start << "\n";
	  Helper::halt( f + " contains row(s) with negative time points" ) ;
	}

      if ( ( !*readon ) && dbl_stop < 0 )
	{
	  //std::cout << " S2 " << dbl_start << "\n";
	  Helper::halt( f + " contains row(s) with negative time points" ) ;
	}
      
      // annot(epoch)/record alignment (to the leftmost second)
      // i.e. to handle staging annots that start at a fractional second onset

      if ( align_annots )
	{
	  // shift backwards to start of second
	  if ( globals::annot_alignment.find( a->name ) !=  globals::annot_alignment.end() )
	    {
	      dbl_start = floor( dbl_start );
	      dbl_stop = floor( dbl_stop );
	    }	       
	}
      
      // convert to uint64_t time-point units
      
      interval.start = Helper::sec2tp( dbl_start );
      
      // assume stop is already specified as 1 past the end, e.g. 30 60
      // a zero-duration interval will therefore have start == stop (i.e. duration = 0)
      // which should be fine

      if ( ! *readon )
	interval.stop  = Helper::sec2tp( dbl_stop );

    }

  
  
  if ( ! *readon )
    {

  
      // TMP remove this
      if ( 0 )  // TMP CHANGE
	{
	  if ( interval.stop == interval.start )
	    ++interval.stop;
	}

      // Check for valid ordering
      
      if ( interval.start > interval.stop )
	Helper::halt( "invalid interval: stop is before start\n" + line );

    }

  // special code for interval that starts before EDF start (have to ignore)
  if ( before_edf_start ) 
    {
      interval.start = 1LL;
      interval.stop  = 0LL;
    }

  return interval;
  
}



int annot_t::load_features( const std::string & f )
{
  
  // set basic values for this annotation type, then add events/features  
  
  // logger << " attaching feature-list file " << f << "\n";
  
  std::ifstream FIN( f.c_str() , std::ios::in );
  
  int line_count = 0;
  
  // tab-delimited

  // tp1 tp2 label key=value key=value
  // with special values  _rgb=255,255,255
  //                      _value={float}

  while ( ! FIN.eof() )
    {
      
      if ( FIN.bad() ) continue;
      std::string line;      
      Helper::safe_getline( FIN , line );      
      if ( FIN.eof() || line == "" ) continue;
      
      std::vector<std::string> tok = Helper::parse( line , "\t" );
      const int n = tok.size();
      if ( n < 3 ) continue;
      
      feature_t feature;
      
      // features work directly in interval-TP coding, so no need to change the end-point
      if ( ! Helper::str2int64( tok[0] , &feature.feature.start ) ) Helper::halt( "bad format " + line + "\n" );
      if ( ! Helper::str2int64( tok[1] , &feature.feature.stop  ) ) Helper::halt( "bad format " + line + "\n" );
      feature.label = tok[2];

      if ( feature.feature.start > feature.feature.stop ) Helper::halt( "bad format, start > stop : " + line + "\n" );
      
      for (int t=3;t<tok.size();t++)
	{
	  std::vector<std::string> tok2 = Helper::parse( tok[t] , "=" );
	  if ( tok2.size() == 1 ) feature.data[ tok2[0] ] = "";
	  else 
	    {
	      
	      feature.data[ tok2[0] ] = tok2[1];
	      
	      if ( tok2[0] == "_rgb" ) 
		{
		  feature.has_colour = true;
		  feature.colour = tok2[1];
		}
	      else if ( tok2[0] == "_val" ) 
		{
		  feature.has_value = Helper::str2dbl( tok2[1] , &feature.value ) ;
		}
		
	    }
	}

      //
      // Add this interval
      //
	  
      instance_t * instance = add( feature.label , feature.feature , "." );
      
      //
      // and append meta-data
      //
      
      instance->add( feature.data );
      
      //
      // and add the type information in too  (even though not every instance may have all types)
      //
      
      std::map<std::string,std::string>::const_iterator ss = feature.data.begin();
      while ( ss != feature.data.end() ) 
	{
	  types[ ss->first ] = globals::A_TXT_T;
	  ++ss;
	}
      
      //
      // row count
      //
      
      ++line_count;

    } // next line
  
  
  //  logger << "  processed " << line_count << " lines\n";
  
  FIN.close();
  
  return line_count;

}


bool annot_t::save( const std::string & t)
{

  //
  // Note: would prefer 'write()' option from annotation_set_t::write()
  //  as this better handles all new formats, etc.
  //

  std::ofstream O1( t.c_str() , std::ios::out );
  
  bool has_vars = types.size() > 0 ;

  //
  // Headers
  //

  O1 << "# " << name;
  
  if ( description != "" )
    O1 << " | " << description;
  else if ( has_vars ) // need a dummy description here 
    O1 << " | " << description ;
  if ( has_vars ) 
    O1 << " |";
  
  std::map<std::string, globals::atype_t>::const_iterator aa = types.begin();
  while ( aa != types.end() )
    {
      O1 << " " << aa->first 
	 << "[" 
	 << globals::type_name[ aa->second ] 
	 << "]";	      
      ++aa;
    }
  
  O1 << "\n";
  
  
  //
  // Interval-based annotation
  //


  annot_map_t::const_iterator ii = interval_events.begin();
  while ( ii != interval_events.end() )
    {
      
      const instance_idx_t & instance_idx = ii->first;

      const instance_t * inst = ii->second;

      O1 << name << "\t";
      
      if ( instance_idx.id != "." && instance_idx.id != "" ) 
	O1 << instance_idx.id << "\t";
      else 
	O1 << ".\t";
      
      if ( instance_idx.ch_str != "." && instance_idx.ch_str != "" )
	O1 << instance_idx.ch_str << "\t";
      else
	O1 << ".\t";

      // start/stop in seconds, with 4 d.p.

      O1 << Helper::dbl2str( instance_idx.interval.start_sec() , globals::time_format_dp ) << "\t"
	 << Helper::dbl2str( instance_idx.interval.stop_sec() , globals::time_format_dp );
      
      if ( inst->data.size() == 0 ) 
	O1 << "\t.";
      else
	{
	  O1 << "\t";
	  
	  std::map<std::string,avar_t*>::const_iterator dd = inst->data.begin();
	  
	  while ( dd != inst->data.end() )
	    {
	      // pipe-delimiter
	      if ( dd != inst->data.begin() ) O1 << "|";
	      // meta-data value
	      O1 << *dd->second;
	      ++dd;
	    }
	}
      
      O1 << "\n";

      //
      // next annotation
      //
      
      ++ii;
    }
	  
  //
  // All done
  //
    
  O1.close();

  return true;
}  


void annot_t::dumpxml( const std::string & filename , bool basic_dumper )
{

  std::map<interval_t,std::vector<std::string> > res;
  
  XML xml( filename );
  if ( ! xml.valid() ) Helper::halt( "invalid annotation file: " + filename );
  
  if ( basic_dumper )
    {
      xml.dump();
      return;
    }
  
  // automatically determine format
  
  std::vector<element_t*> nsrr_format_test = xml.children( "PSGAnnotation" );
  bool profusion_format = nsrr_format_test.size() == 0 ;
  
  const std::string EventConcept = profusion_format ? "Name"           : "EventConcept" ;
  const std::string epoch_parent = profusion_format ? "CMPStudyConfig" : "PSGAnnotation" ;
  
  //
  // Epoch Length
  //

  int epoch_sec = -1;

  // Document --> CMPStudyConfig --> EpochLength 
  // PSGAnnotation --> EpochLength
    
  std::vector<element_t*> elements = xml.children( epoch_parent );
  for (int e=0;e<elements.size();e++)
    {
      if ( elements[e]->name == "EpochLength" ) 
	{
	  if ( ! Helper::str2int(  elements[e]->value , &epoch_sec ) ) 
	    Helper::halt( "bad EpochLength" ) ;
	  std::stringstream ss ;
	  ss << ".\t.\tEpochLength\t" << epoch_sec << "\n";
	  res[ interval_t(0,0) ].push_back( ss.str() );
	  break;
	}
    }
  
  if ( epoch_sec == -1 ) 
    {
      Helper::warn( "did not find EpochLength in XML, defaulting to " 
		    + Helper::int2str( globals::default_epoch_len ) + " seconds" );
      epoch_sec = globals::default_epoch_len;
    }


  //
  // Scored Events
  //

  std::vector<element_t*> scored = xml.children( "ScoredEvents" );
  
  // assume all annotations will then be under 'ScoredEvent'
  // Profusion: with children: 'EventConcept' , 'Duration' , 'Start' , and optionally 'Notes'
  // NSRR     : with children: 'Name'         , 'Duration' , 'Start' , and optionally 'Notes'
  
  for (int i=0;i<scored.size();i++)
    {

      element_t * e = scored[i];
      
      if ( ! Helper::iequals( e->name , "ScoredEvent" ) ) continue;
      
      element_t * concept  = (*e)( EventConcept );
      if ( concept == NULL ) concept = (*e)( "name" );
      
      element_t * start    = (*e)( "Start" );
      if ( start == NULL ) start = (*e)( "time" );

      element_t * duration = (*e)("Duration" );

      element_t * notes    = (*e)("Notes" );
     
      element_t * type    = (*e)( "EventType" );
      

      //      if ( concept == NULL || start == NULL || duration == NULL ) continue;
      if ( concept == NULL ) continue;

      double start_sec = 0, stop_sec = 0 , duration_sec = 0;
      uint64_t start_tp = 0 , stop_tp = 0;

      if ( duration != NULL )
	{
	  if ( ! Helper::str2dbl( duration->value , &duration_sec ) ) 
	    Helper::halt( "bad value in annotation" );	  		  
	}
      
      if ( start != NULL ) 
	{
	  if ( ! Helper::str2dbl( start->value , &start_sec ) ) 
	    Helper::halt( "bad value in annotation" );
	  stop_sec = start_sec + duration_sec;
	  start_tp = Helper::sec2tp( start_sec );
	  stop_tp = start_tp + Helper::sec2tp( duration_sec ) ; 
	  
	  // EDIT OUT	  
// 	  // in case duration was 0, make this a 1-time-unit event
// 	  if ( start_tp == stop_tp ) ++stop_tp;
	  
	  // MAKE ALL points one past the end
	  ++stop_tp;

	  //stop_tp = start_tp + (uint64_t)( globals::tp_1sec * duration_sec ) - 1LLU ; 
	  
	}

      interval_t interval( start_tp , stop_tp );
      
      std::stringstream ss;      

      if ( start != NULL ) 
	{
	  ss << start_sec ;
	  if ( duration != NULL ) ss << " - " << stop_sec << "\t"
				     << "(" << duration_sec << " secs)\t";
	  else ss << ".\t";
	}
      else ss << ".\t.\t";
      
      if ( type != NULL ) 
	ss << type->value << "\t";
      else 
	ss << ".\t";
      ss << concept->value << "\t";
      if ( notes != NULL ) ss << "\t" << notes->value ;      
      ss << "\n";
      res[ interval ].push_back( ss.str() );

    }
  

  //
  // Sleep Stages (Profusion only: in NSRR format, staging is incorporated as ScoredEvent)
  //
    
  // Profusion: under 'SleepStages', with children 'SleepStage' which comprises an integer value
  // 0  wake
  // 1  NREM1
  // 2  NREM2
  // 3  NREM3
  // 4  NREM4
  // 5  REM
  // Otherwise 'Unscored'
  

  if ( profusion_format )
    {

      std::vector<element_t*> scored = xml.children( "SleepStages" );
      
      int seconds = 0;

      for (int i=0;i<scored.size();i++)
	{

	  element_t * e = scored[i];

	  if ( e->name != "SleepStage" ) continue;
	  
	  std::string stg = "?";
	  if      ( e->value == "0" ) stg = "W";
	  else if ( e->value == "1" ) stg = "N1";
	  else if ( e->value == "2" ) stg = "N2";
	  else if ( e->value == "3" ) stg = "N3";
	  else if ( e->value == "4" ) stg = "NREM4";
	  else if ( e->value == "5" ) stg = "R";	 

	  // remap?
	  stg = nsrr_t::remap( stg );

	  if ( stg == "" ) continue;
	  
	  interval_t interval( Helper::sec2tp( seconds ) , 
			       Helper::sec2tp( seconds + epoch_sec ) );

// 	  interval_t interval( (uint64_t)(seconds * globals::tp_1sec ) , 
// 			       (uint64_t)(( seconds + epoch_sec ) * globals::tp_1sec ) );

	  std::stringstream ss;      
	  ss << seconds << " - " << seconds + epoch_sec << "\t"
	     << "(" << epoch_sec << " secs)\t"
	     << "SleepStage" << "\t"	     
	     << stg << "\n";
	  res[ interval ].push_back( ss.str() );
	  
	  // advance to the next epoch
	  seconds += epoch_sec;
	  
	}
           
    }

  //
  // Report
  //
  
  std::map<interval_t,std::vector<std::string> >::const_iterator ii = res.begin();
  while ( ii != res.end() )
    {
      std::vector<std::string>::const_iterator jj = ii->second.begin();
      while ( jj != ii->second.end() )
	{
	  std::cout << *jj;
	  ++jj;
	}
      ++ii;
    }
 
}

bool annot_t::loadxml( const std::string & filename , edf_t * edf )
{

  XML xml( filename );

  if ( ! xml.valid() ) Helper::halt( "invalid annotation file: " + filename );

  //
  // Determine format: Profusion or NSRR or Luna ? 
  //
  
  std::vector<element_t*> nsrr_format_test = xml.children( "PSGAnnotation" );
  std::vector<element_t*> luna_format_test = xml.children( "Annotations" );
  
  bool profusion_format = nsrr_format_test.size() == 0 ;
  bool luna_format = luna_format_test.size() > 0 ;

  if ( globals::param.has( "profusion" ) ) profusion_format = true;
  
  if ( luna_format ) return loadxml_luna( filename , edf );

  const std::string EventConcept = profusion_format ? "Name"           : "EventConcept" ;
  const std::string epoch_parent = profusion_format ? "CMPStudyConfig" : "PSGAnnotation" ;
  
  
  std::vector<element_t*> scored = xml.children( "ScoredEvents" );
  

  //
  // NSRR format:
  //
  
  // assume all annotations will then be under 'ScoredEvent'
  // with children: 'EventConcept' , 'Duration' , 'Start' , and optionally 'SignalLocation' and 'Notes'
  // for all other children, add as meta-data (type = 'str') 

  //
  // Profusion format
  //
  
  // assume all annotations will then be under 'ScoredEvent'
  // with children: 'Name' , 'Duration' , 'Start' , and optionally 'Notes'
  
  // SleepStages: under separate 'SleepStages' parent
  // children elements 'SleepStage' == integer
  // ASSUME these are 30-s epochs, starting at 0
  

  //
  // First pass through all 'ScoredEvent's,, creating each annotation
  //

  std::set<std::string> added;
  
  for (int i=0;i<scored.size();i++)
    {
      
      element_t * e = scored[i];
      
      if ( ! Helper::iequals( e->name , "ScoredEvent" ) ) continue;
      
      element_t * concept  = (*e)( EventConcept );
      if ( concept == NULL ) concept = (*e)( "name" );
      
      if ( concept == NULL ) continue;
      
      // skip this..
      if ( concept->value == "Recording Start Time" ) continue;
      
      
      // annotation remap?
      std::string original_label = concept->value;
      concept->value = nsrr_t::remap( concept->value );
      if ( concept->value == "" ) continue;
      
      // are we checking whether to add this file or no? 
      if ( globals::specified_annots.size() > 0 && 
	   globals::specified_annots.find( concept->value ) == globals::specified_annots.end() ) continue;
      
      // already found?
      if ( added.find( concept->value ) != added.end() ) continue;

      // otherwise, add

      if ( original_label != concept->value )
	  edf->timeline.annotations.aliasing[ concept->value ] = original_label ;

      annot_t * a = edf->timeline.annotations.add( concept->value );
      a->description = "XML-derived";
      a->file = filename;
      a->type = globals::A_FLAG_T; // not expecting any meta-data
      added.insert( concept->value );
    }

  //
  // Profusion-formatted sleep-stages?
  //
  
  if ( profusion_format )
    {
      
      std::vector<element_t*> scored = xml.children( "SleepStages" );
      
      for (int i=0;i<scored.size();i++)
	{
	  element_t * e = scored[i];

	  if ( e->name != "SleepStage" ) continue;
	  
	  std::string ss = "Unscored";
	  if      ( e->value == "0" ) ss = "wake";
	  else if ( e->value == "1" ) ss = "NREM1";
	  else if ( e->value == "2" ) ss = "NREM2";
	  else if ( e->value == "3" ) ss = "NREM3";
	  else if ( e->value == "4" ) ss = "NREM4";
	  else if ( e->value == "5" ) ss = "REM";	 

	  // annotation remap?
	  ss = nsrr_t::remap( ss );
	  if ( ss == "" ) continue;
	  
	  // are we checking whether to add this file or no? 
	  
	  if ( globals::specified_annots.size() > 0 && 
	       globals::specified_annots.find( ss ) == globals::specified_annots.end() ) continue;
	  
	  // already found?
	  if ( added.find( ss ) != added.end() ) continue;
	  
	  // otherwise, add
	  annot_t * a = edf->timeline.annotations.add( ss );
	  a->description = "XML-derived";
	  a->file = filename;
	  a->type = globals::A_FLAG_T; // not expecting any meta-data from XML
	  added.insert( ss );
	  
	}
    }


  
  //
  // Back through, adding instances now we've added all annotations
  //

  for (int i=0;i<scored.size();i++)
    {
      
      element_t * e = scored[i];
      
      if ( ! Helper::iequals( e->name , "ScoredEvent" ) ) continue;
      
      element_t * concept  = (*e)( EventConcept );
      if ( concept == NULL ) concept = (*e)( "name" );

      // skip if we are not interested in this element
      if ( added.find( concept->value ) == added.end() ) continue;
      
      element_t * start    = (*e)( "Start" );
      if ( start == NULL ) start = (*e)( "time" );

      element_t * duration = (*e)( "Duration" );
      element_t * notes    = (*e)( "Notes" );
      element_t * signal   = (*e)( "SignalLocation" );

      if ( concept == NULL || start == NULL || duration == NULL ) continue;
      

      // Luna format XML can also specify the channel
      element_t * channel = (*e)( "Channel" );
      
      // otherwise, add 
      double start_sec, duration_sec;
      if ( ! Helper::str2dbl( start->value , &start_sec ) ) Helper::halt( "bad value in annotation" );
      if ( ! Helper::str2dbl( duration->value , &duration_sec ) ) Helper::halt( "bad value in annotation" );

      uint64_t start_tp = Helper::sec2tp( start_sec );

//       uint64_t stop_tp  = duration_sec > 0 
// 	? start_tp + (uint64_t)( duration_sec * globals::tp_1sec ) - 1LLU 
// 	: start_tp;

      // std::cout << "xxx " <<  concept->value << "\t" << start_sec << "\t" << duration_sec << "\t"
      // 		<< duration_sec * globals::tp_1sec  << "\t" <<  (uint64_t)( duration_sec ) * globals::tp_1sec  << "\n";
      
      // stop is defined as 1 unit past the end of the interval
      uint64_t stop_tp  = duration_sec > 0 
	? start_tp + Helper::sec2tp( duration_sec )
	: start_tp ; // for zero-point interval (starts and stops at same place)
      
      interval_t interval( start_tp , stop_tp );
      
      annot_t * a = edf->timeline.annotations.add( concept->value );
      
      if ( a == NULL ) Helper::halt( "internal error in loadxml()");

      std::string sigstr = signal != NULL ? signal->value : ( channel != NULL ? channel->value : "." ) ; 

      // swap spaces from sigstr (channel label)?
      if ( globals::replace_channel_spaces )
	sigstr = Helper::search_replace( sigstr , ' ' , globals::space_replacement );
      
      //instance_t * instance = a->add( concept->value , interval , sigstr );
      // class name is <ConceptValue> tag, so make instance ID null
      instance_t * instance = a->add( "."  , interval , sigstr );
      
      // any notes?  set as TXT, otherwise it will be listed as a FLAG
      if ( notes ) 
	{
	  instance->set( concept->value , notes->value );  
	}

      //
      // any other children of ScoredEvent?  add as string key/value meta-data
      //
      
      const std::vector<element_t*> & kids = e->child;
      
      for (int i=0;i<kids.size();i++)
	{
          element_t * ee = kids[i];
	  if ( ee->name == "EventConcept" ) continue;
	  if ( ee->name == "EventType" ) continue;
	  if ( ee->name == "Notes" ) continue;
	  if ( ee->name == "Channel" ) continue;
	  if ( ee->name == "SignalLocation" ) continue;
	  if ( ee->name == "Start" ) continue;
	  if ( ee->name == "Duration" ) continue;
	  if ( ee->name == "name" ) continue;
	  if ( ee->name == "time" ) continue;

	  // add as meta-data to this instance
	  instance->set( ee->name , ee->value );

	}

	       
      
    }
  
  
  //
  // Profusion-formatted sleep-stages?
  //
  
  if ( profusion_format )
    {
  
      std::vector<element_t*> scored = xml.children( "SleepStages" );
      
      int start_sec = 0;
      int epoch_sec = 30;

      // assume 30-second epochs, starting from 0...

      for (int i=0;i<scored.size();i++)
	{
	  element_t * e = scored[i];

	  if ( e->name != "SleepStage" ) continue;
	  
	  std::string ss = "Unscored";
	  if      ( e->value == "0" ) ss = "wake";
	  else if ( e->value == "1" ) ss = "NREM1";
	  else if ( e->value == "2" ) ss = "NREM2";
	  else if ( e->value == "3" ) ss = "NREM3";
	  else if ( e->value == "4" ) ss = "NREM4";
	  else if ( e->value == "5" ) ss = "REM";	 

	  // annotation remap?
	  ss = nsrr_t::remap( ss );
	  if ( ss == "" ) continue;
	  
	  // skip if we are not interested in this element
	  
	  if ( added.find( ss ) == added.end() ) continue;
	  
	  // otherwise, add
	  
	  uint64_t start_tp = Helper::sec2tp( start_sec );
	  uint64_t stop_tp  = start_tp + Helper::sec2tp( epoch_sec ) ; // 1-past-end encoding
	  
	  // advance to the next epoch
	  start_sec += epoch_sec;
	  
	  interval_t interval( start_tp , stop_tp );	  
	  
	  annot_t * a = edf->timeline.annotations.add( ss );
		  
	  // . indicates no associated channel
	  instance_t * instance = a->add( ss , interval , "." );      
      
	  instance->set( ss );
	  
	}
            
    }



  //
  // Misc. test code: have signal descriptions in XMLs
  //

  if ( false )
    {
      std::vector<element_t*> signals = xml.children( "Signals" );

      //std::cout << "signals size = " << signals.size() << "\n";
      
      for (int i=0; i<signals.size(); i++)
	{
	  
	  element_t * e = signals[i];
	  
	  //       cmd_t::signal_alias( str )
	  // 	"canonical|alias1|alias2"
	  
	  element_t * label = (*e)("Label");
	  element_t * canonical_label = (*e)("CanonicalLabel");
	  element_t * desc = (*e)("Description");
	  
	  if ( label ) std::cout << "label = " << label->value << "\n";
	  if ( canonical_label ) std::cout << "canon " << canonical_label->value << "\n";
	  if ( desc ) std::cout << "Desc " << desc->value << "\n";
	  
	  if ( label->value != "" && 
	       canonical_label->value != "" && 
	       label->value != canonical_label->value )
	    {
	      logger << "  changing " << label->value << " to canonical label " << canonical_label->value << "\n";
	      edf->header.rename_channel( label->value , canonical_label->value );
	      cmd_t::signal_alias( canonical_label->value + "|" + label->value ); // include this????
	    }
	  
	  std::vector<element_t*> attr = e->children( "Attributes" );
	  
	  for (int j=0;j<attr.size();j++)
	    {

	      element_t * ee = attr[j];
	      
	      if ( ee->name != "Attribute" ) continue;
	      
	      element_t * aname = (*ee)( "AttributeKey" );
	      element_t * aval  = (*ee)( "AttributeValue" );
	      if ( aname != NULL && aval != NULL ) 
		std::cout << aname->value << " = " << aval->value << "\n";
	      
	    }
	
	  // Signal
	  //  -Label
	  //  -CanonicalLabel
	  //  -Description
	  //  -Attributes
	  //    Attribute
	  //     -AttrinuteKey 
	  //     -AttributeLabel      
	  
	}
    }



  return true;
}






bool annot_t::savexml( const std::string & f )
{
  Helper::halt( "not yet implemented" );
  return true;
}




annot_map_t annot_t::extract( const interval_t & window ) 
{
  
  //
  // Fetch all annotations that overlap this window
  // where overlap is defined as region A to B-1 for interval_t(A,B)
  //

  annot_map_t r; 
  
  // urghhh... need to implement a much better search... 
  // but for now just use brute force... :-(
  
  annot_map_t::const_iterator ii = interval_events.begin();
  while ( ii != interval_events.end() )
    {
      const interval_t & a = ii->first.interval;
      if ( a.overlaps( window ) ) r[ ii->first ] = ii->second;
      else if ( a.is_after( window ) ) break;
      ++ii;
    }
  
  return r;

  
  // // find first interval just /before/ then step over deciding which
  // // should fit in

  // interval_evt_map_t::const_iterator ii 
  //   = interval_events.lower_bound( window );

  // interval_evt_map_t::const_iterator jj = ii; 
  
  // // if returns iterator to end(), could still overlap the 
  // // one or more events at the end, so need to check for this
  // // special case of looking backwards either way; 
  
  // // back
  // while ( 1 ) 
  //   {

  //     if ( jj == interval_events.begin() ) break;
  //     --jj;
  //     const interval_t & a = jj->first;
      
  //     //      std::cout << "extract: considering " << a.start << " - " << a.stop << "\n";
      
  //     if ( a.overlaps( window ) ) 
  // 	{
  // 	  r[ a ] = jj->second;
  // 	  //std::cout << " found overlap\n";
  // 	}
  //     else if ( a.is_before( window ) )
  // 	{
  // 	  // also need to consider intervals that start before 'a' but still span window
  // 	  // i.e. go back until test interval doesn't overlap 'a'
  // 	  interval_evt_map_t::const_iterator kk = jj; 
  // 	  interval_t prev = a;
  // 	  while ( 1 ) 
  // 	    {
  // 	      if ( kk == interval_events.begin() ) break;
  // 	      --kk;
  // 	      interval_t b = kk->first;
  // 	      if ( ! prev.overlaps( b ) ) break; // really all done now
  // 	      if ( b.overlaps( window ) ) 
  // 		{ 
  // 		  r[ b ] = kk->second; 
  // 		  std::cout << " XXX added an extra!\n"; 		  
  // 		}
  // 	      prev = b;
  // 	    }
  // 	  break;      
  // 	}
  //   }

  // std::cout << "now forward...\n";
  // // forward
  // while ( 1 ) 
  //   {      
  //     if ( ii == interval_events.end() ) break;
  //     const interval_t & a = ii->first;
  //     std::cout << "extract: considering " << a.start << " - " << a.stop << "\n";
  //     if ( a.overlaps( window ) ) 
  // 	{
  // 	  r[ a ] = ii->second;
  // 	  std::cout << " found overlap\n";
  // 	}
      
  //     else if ( a.is_after( window ) ) break;      
  //     ++ii;
  //   }
  // std::cout << "  done...\n";
  // return r;
  
}


bool globals::is_stage_annotation( const std::string & s )
{
  // try to guess if this is a stage duration;
  // this will be called on loading annotations, but after any remapping
  // thus, we should be safe to use canonical forms  

  sleep_stage_t ss = globals::stage( s );
  return ss != UNKNOWN ;

}


bool annotation_set_t::make_sleep_stage( const timeline_t & tl ,
					 const bool force_remake , 
					 const std::string & a_wake , 
					 const std::string & a_n1 , 
					 const std::string & a_n2 , 
					 const std::string & a_n3 , 
					 const std::string & a_n4 , 
					 const std::string & a_rem ,
					 const std::string & a_light , 
					 const std::string & a_other )
{


  //
  // force a remake?
  //
  
  if ( force_remake )
    clear( "SleepStage" );
      
  //
  // already made?
  //
  
  if ( find( "SleepStage" ) != NULL )
    return false; 
  
   
  //
  // Use default annotation labels, if not otherwise specified
  // 

  std::string dwake, dn1, dn2, dn3, dn4, drem, dlight, dother, dunknown;

  std::map<std::string,annot_t*>::const_iterator ii = annots.begin();
  while ( ii != annots.end() )
    {
      
      const std::string & s = ii->first;
      
      // this function takes care of any prefix specified via ss-prefix
      // i.e.   if prefix is 'p' than 'pN1' will match to 'N1' etc
      
      sleep_stage_t ss = globals::stage( s );
      
      if      ( ss == WAKE )     dwake = s;
      else if ( ss == NREM1 )    dn1 = s;
      else if ( ss == NREM2 )    dn2 = s;
      else if ( ss == NREM3 )    dn3 = s;
      else if ( ss == NREM4 )    dn4 = s;
      else if ( ss == REM )      drem = s;
      else if ( ss == LIGHTS_ON ) dlight = s;
      else if ( ss == UNSCORED ) dother = s;
      else if ( ss == MOVEMENT ) dother = s;
      else if ( ss == ARTIFACT ) dother = s;
      ++ii;
    }
  
  std::vector<std::string> v_wake  = Helper::parse( a_wake , "," );
  std::vector<std::string> v_n1    = Helper::parse( a_n1 , "," );
  std::vector<std::string> v_n2    = Helper::parse( a_n2 , "," );
  std::vector<std::string> v_n3    = Helper::parse( a_n3 , "," );
  std::vector<std::string> v_n4    = Helper::parse( a_n4 , "," );
  std::vector<std::string> v_rem   = Helper::parse( a_rem , "," );
  std::vector<std::string> v_light = Helper::parse( a_light , "," );
  std::vector<std::string> v_other = Helper::parse( a_other , "," );

  // add defaults
  if ( v_wake.size() == 0 ) v_wake.push_back( dwake );
  if ( v_n1.size() == 0 ) v_n1.push_back( dn1 );
  if ( v_n2.size() == 0 ) v_n2.push_back( dn2 );
  if ( v_n3.size() == 0 ) v_n3.push_back( dn3 );
  if ( v_n4.size() == 0 ) v_n4.push_back( dn4 );
  if ( v_rem.size() == 0 ) v_rem.push_back( drem );
  if ( v_light.size() == 0 ) v_light.push_back( dlight );
  if ( v_other.size() == 0 ) v_other.push_back( dother );

  //
  // find annotations, allowing a comma-delimited list
  //

  std::vector<annot_t *> wakes, n1s, n2s, n3s, n4s, rems, lights, others;
  
  for (int a=0;a<v_wake.size();a++)
    wakes.push_back( find( v_wake[a] ) );
  
  for (int a=0;a<v_n1.size();a++)
    n1s.push_back( find( v_n1[a] ) );
  
  for (int a=0;a<v_n2.size();a++)
    n2s.push_back( find( v_n2[a] ) );
  
  for (int a=0;a<v_n3.size();a++)
    n3s.push_back( find( v_n3[a] ) );

  for (int a=0;a<v_n4.size();a++)
    n4s.push_back( find( v_n4[a] ) );

  for (int a=0;a<v_rem.size();a++)
    rems.push_back( find( v_rem[a] ) );
  
  for (int a=0;a<v_light.size();a++)
    lights.push_back( find( v_light[a] ) );

  for (int a=0;a<v_other.size();a++)
    others.push_back( find( v_other[a] ) );



  //
  // Check we had sensible annotations
  //

  int assigned = 0;
  for (int a=0;a<n1s.size();a++) if ( n1s[a] != NULL ) ++assigned;
  for (int a=0;a<n2s.size();a++) if ( n2s[a] != NULL ) ++assigned;
  for (int a=0;a<n3s.size();a++) if ( n3s[a] != NULL ) ++assigned;
  for (int a=0;a<rems.size();a++) if ( rems[a] != NULL ) ++assigned;
  for (int a=0;a<wakes.size();a++) if ( wakes[a] != NULL ) ++assigned;
  for (int a=0;a<lights.size();a++) if ( lights[a] != NULL ) ++assigned;
  if ( assigned == 0 ) return false;
  
  //
  // Align all putative stages, and if we see point markers, extend to the next (end) annot
  //
  
  std::map<interval_t,sleep_stage_t> stages;

  for ( int i=0; i<wakes.size(); i++ )
    {
      annot_t * wake = wakes[i];
      if ( wake )
	{
	  annot_map_t & events = wake->interval_events;
          annot_map_t::const_iterator ee = events.begin();
          while ( ee != events.end() )
            {
              stages[ ee->first.interval ] = WAKE;
              ++ee;
            }
	}
    }

  
  for ( int i=0; i<n1s.size(); i++ )
    {
      annot_t * stg = n1s[i];
      if ( stg )
        {
          annot_map_t & events = stg->interval_events;
          annot_map_t::const_iterator ee = events.begin();
          while ( ee != events.end() )
            {
              stages[ ee->first.interval ] = NREM1;
              ++ee;
            }
	}
    }

  for ( int i=0; i<n2s.size(); i++ )
    {
      annot_t * stg = n2s[i];
      if ( stg )
        {
          annot_map_t & events = stg->interval_events;
          annot_map_t::const_iterator ee = events.begin();
          while ( ee != events.end() )
            {
              stages[ ee->first.interval ] = NREM2;
              ++ee;
            }
	}
    }

  for ( int i=0; i<n3s.size(); i++ )
    {
      annot_t * stg = n3s[i];
      if ( stg )
        {
          annot_map_t & events = stg->interval_events;
          annot_map_t::const_iterator ee = events.begin();
          while ( ee != events.end() )
            {
              stages[ ee->first.interval ] = NREM3;
              ++ee;
            }
        }
    }

  for ( int i=0; i<n4s.size(); i++ )
    {
      annot_t * stg = n4s[i];
      if ( stg )
        {
          annot_map_t & events = stg->interval_events;
          annot_map_t::const_iterator ee = events.begin();
          while ( ee != events.end() )
            {
              stages[ ee->first.interval ] = NREM4;
              ++ee;
            }
        }
    }

  for ( int i=0; i<rems.size(); i++ )
    {
      annot_t * stg = rems[i];
      if ( stg )
	{
          annot_map_t & events = stg->interval_events;
          annot_map_t::const_iterator ee = events.begin();
          while ( ee != events.end() )
            {
              stages[ ee->first.interval ] = REM;
              ++ee;
            }
        }
    }
  
  for ( int i=0; i<lights.size(); i++ )
    {
      annot_t * stg = lights[i];
      if ( stg )
	{
          annot_map_t & events = stg->interval_events;
          annot_map_t::const_iterator ee = events.begin();
          while ( ee != events.end() )
            {
              stages[ ee->first.interval ] = LIGHTS_ON;
              ++ee;
            }
        }
    }

  for ( int i=0; i<others.size(); i++ )
    {
      annot_t * stg = others[i];
      if ( stg )
        {
          annot_map_t & events = stg->interval_events;
          annot_map_t::const_iterator ee = events.begin();
          while ( ee != events.end() )
            {
              stages[ ee->first.interval ] = UNSCORED;
              ++ee;
            }
        }
    }


  //
  // Now, look through stages[] and extend any zero-point stages; also, flag conflicts 
  //

  interval_t prior;
  std::vector<interval_t> vec_intervals;
  std::vector<sleep_stage_t> vec_stages;
  
  std::map<interval_t,sleep_stage_t>::const_iterator jj = stages.begin();
  while ( jj != stages.end() )
    {
      
      interval_t curr = jj->first;

      // ensure no overlaps
      if ( jj != stages.begin() )
	{
	  if ( curr.start < prior.stop )
	    {
	      std::cerr << " current interval : " << curr.start << " .. " << curr.stop << "  S = " << jj->second << "\n"
			<< " prior            : " << prior.start << " .. " << prior.stop << "  S = " << "\n";
	      
	      Helper::halt( "bad, overlapping sleep stages" );
	    }
	}
      
      vec_intervals.push_back( curr );
      vec_stages.push_back( jj->second );
      
      // save prior interval
      prior = curr;
      
      ++jj;
    }

  
  //  std::cout	<<  "stg " << vec_intervals[0].start <<	" -- " << vec_intervals[0].stop << "\t" << vec_stages[0] << "\n";

  // start from the second entry
  for (int j=1; j<vec_stages.size(); j++)
    {
      //std::cout <<  "stg " << vec_intervals[j].start <<  " -- " << vec_intervals[j].stop << "\t" << vec_stages[j] << "\n";
      if ( vec_intervals[j-1].duration() == 0 )
	{
	  //std::cout << " adjusting...\n";
	  vec_intervals[j-1].stop = vec_intervals[j].start; // i.e. 1-tp past end
	}
    }
  
  // also handle 'all wake' condition
  if ( vec_stages.size() == 1 )
    {
      if ( vec_intervals[0].duration() == 0 )
	{
	  // set to end (1 past the last TP)
	  vec_intervals[0].stop = tl.last_time_point_tp + 1LLU;
	}
    }
  
  
  
  //
  // Create the 'SleepStage' unified annotation (used by HYPNO, STAGE, and POPS)
  //

  // ensure cleared if already exists; if it doesn't this command won't
  // do anything, so okay
  clear( "SleepStage" );
  
  annot_t * ss = add( "SleepStage" );
  
  ss->description = "SleepStage";
  
  for ( int i=0; i<vec_stages.size(); i++ )
    ss->add( globals::stage( vec_stages[i] ) , vec_intervals[i] , "." );
  

  return true;
  
}


uint64_t annot_t::minimum_tp() const
{
  if ( interval_events.size() == 0 ) return 0;
  return interval_events.begin()->first.interval.start;
}


uint64_t annot_t::maximum_tp() const
{
  if ( interval_events.size() == 0 ) return 0;
  return (--interval_events.end())->first.interval.stop;
}


std::set<std::string> annot_t::instance_ids() const
{
  std::set<std::string> r;
  annot_map_t::const_iterator ii = interval_events.begin();
  while ( ii != interval_events.end() )
    {
      r.insert( ii->first.id );
      ++ii;
    }
  return r;
}



//
// Silly helper functions...
//

std::vector<bool> annot_t::as_bool_vec( const std::vector<int> & x ) { 
  std::vector<bool> y( x.size() );
  for (int i=0;i<x.size();i++) y[i] = (bool)x[i];
  return y;
} 

std::vector<bool> annot_t::as_bool_vec( const std::vector<double> & x ) { 
  std::vector<bool> y( x.size() );
  for (int i=0;i<x.size();i++) y[i] = (bool)x[i];
  return y;
} 

std::vector<bool> annot_t::as_bool_vec( const std::vector<std::string> & x ) { 
  std::vector<bool> y( x.size() );
  for (int i=0;i<x.size();i++) y[i] = Helper::yesno( x[i] );
  return y;
} 


// to int
std::vector<int> annot_t::as_int_vec( const std::vector<bool> & x ) { 
  std::vector<int> y( x.size() );
  for (int i=0;i<x.size();i++) y[i] = (int)x[i];
  return y;
} 

std::vector<int> annot_t::as_int_vec( const std::vector<double> & x ) { 
  std::vector<int> y( x.size() );
  for (int i=0;i<x.size();i++) y[i] = (int)round(x[i]);
  return y;
} 

std::vector<int> annot_t::as_int_vec( const std::vector<std::string> & x ) { 
  std::vector<int> y( x.size() );
  for (int i=0;i<x.size();i++) y[i] = (int)Helper::yesno( x[i] );
  return y;
} 

// to dbl
std::vector<double> annot_t::as_dbl_vec( const std::vector<bool> & x ) { 
  std::vector<double> y( x.size() );
  for (int i=0;i<x.size();i++) y[i] = x[i];
  return y;
} 

std::vector<double> annot_t::as_dbl_vec( const std::vector<int> & x ) { 
  std::vector<double> y( x.size() );
  for (int i=0;i<x.size();i++) y[i] = x[i];
  return y;
} 

std::vector<double> annot_t::as_dbl_vec( const std::vector<std::string> & x ) { 
  std::vector<double> y( x.size() );
  for (int i=0;i<x.size();i++) y[i] = Helper::yesno( x[i] );
  return y;
} 

// to txt
std::vector<std::string> annot_t::as_txt_vec( const std::vector<bool> & x ) { 
  std::vector<std::string> y( x.size() );
  for (int i=0;i<x.size();i++) y[i] = x[i] ? "true" : "false" ;
  return y;
} 

std::vector<std::string> annot_t::as_txt_vec( const std::vector<int> & x ) { 
  std::vector<std::string> y( x.size() );
  for (int i=0;i<x.size();i++) y[i] = x[i] == 0 ? "false" : "true" ;
  return y;
} 

std::vector<std::string> annot_t::as_txt_vec( const std::vector<double> & x ) { 
  std::vector<std::string> y( x.size() );
  for (int i=0;i<x.size();i++) y[i] = x[i] == 0 ? "false" : "true" ;
  return y;
} 




void annotation_set_t::write( const std::string & filename , param_t & param , edf_t & edf )
{
  
  // write all annotations here as a single file; 
  // either XML or .annot file format
  // default is as .annot
  // order by time

  const bool xml_format = param.has( "xml" ) || Helper::file_extension( filename , "xml" ) ;
  
  //
  // use hh:mm:ss, if possible, instead of elapsed seconds (for .annot only)
  //

  bool hms = param.has( "hms" ) || param.has( "dhms" );
  const bool dhms = param.has( "dhms" );

  //
  // If from internal EDF+D, write w/ time-stamps for standard EDF
  // i.e. collapse times
  //

  const bool collapse_disc = param.has( "collapse" );
  
  //
  // Min duration (e.g. to ensure we have 30 sec epochs)
  //

  const bool has_min_dur = param.has( "min-dur" ) && param.requires_dbl( "min-dur" ) > 0 ;
  
  const double min_dur = has_min_dur ?  param.requires_dbl("min-dur" )  : 0 ;
    
    
  //
  // for complete XML compatibility
  //

  // no outputs other than data rows (no class line)  
  const bool minimal = param.has( "minimal" ) || param.has( "min" ) ;
  
  const bool add_specials = param.has( "specials" );
  
  // in .annot mode only, skip # headers
  const bool add_headers = param.has( "headers" );
  
  // ensure date is here too (to allow 'dhms' mode printing)
  clocktime_t starttime( edf.header.startdate , edf.header.starttime );
  
  if ( hms && ! starttime.valid ) 
    {
      logger << " ** could not find valid start-time in EDF header **\n";
      hms = false;
    }

  
  

  
  //
  // Any offsets specified to annotations for output? (i.e. via ALIGN)
  //

  // nb. this could have been set by the ALIGN command -- although, we
  //    are taking that as redundant / too complicated now... so
  //    here also allow direct specificaiton (in secs) for WRITE-ANNOT
  //    but in that case it is interpeted here as +ve (i.e. to add, not subjtract)
  
  if ( param.has( "offset" ) )
    {
      double s1 = param.requires_dbl( "offset" );      
      annot_offset = s1 * globals::tp_1sec;
      annot_offset_dir = +1;
    }
  
  if ( annot_offset )
    logger << "  applying a offset of " << ( annot_offset_dir == 1 ? "+" : "-" ) 
	   << annot_offset * globals::tp_duration
	   << " seconds to all annotations when writing out\n";

  //
  // either for all annots, or just a subset
  //

  std::set<std::string> annots2write = param.strset( "annot" );
  if ( annots2write.size() > 0 )
    logger << "  writing a subset of all annotations, based on " << annots2write.size() << " specified\n";
  
  if ( filename == "" ) Helper::halt( "bad filename for WRITE-ANNOTS" );

  logger << "  writing annotations (" 
	 << ( xml_format ? ".xml" : ".annot" ) 
	 << " format) to " 
	 << filename << "\n";

  std::ofstream O1( filename.c_str() , std::ios::out );
  
  if ( xml_format ) 
    {
      
      // XML header
      
      O1 << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\n";
      O1 << "<Annotations>\n\n";
      O1 << "<SoftwareVersion>luna-" << globals::version << "</SoftwareVersion>\n\n";
      
      O1 << "<StartTime>" << start_hms << "</StartTime>\n";
      O1 << "<Duration>" << duration_hms << "</Duration>\n";
      O1 << "<DurationSeconds>" << duration_sec << "</DurationSeconds>\n";      
      O1 << "<EpochLength>" << epoch_sec << "</EpochLength>\n";
      
      O1 << "\n";
      
      //
      // Loop over each annotation
      //
      
      std::vector<std::string> anames = names();

      //
      // Track all instances (in order)      
      //

      std::set<instance_idx_t> events;
      
      //
      // Annotation header
      //
      
      O1 << "<Classes>\n";
      
      for (int a=0;a<anames.size();a++)
	{
	  
	  //   <Class>
	  //    <Name>Annotation Name</Name>
	  //    <Description>Annotation Description</Description>
	  //     <Variable name="label" type="type">Numeric variable name</Variable>
	  //     <Variable name="label" type="type">Numeric variable name</Variable>
	  //     <Variable name="label" type="type">Numeric variable name</Variable>
	  //   </Class>


	  if ( annots2write.size()
	       && annots2write.find( anames[a] ) == annots2write.end() )
	    continue;
	  	  
	  annot_t * annot = find( anames[a] );
	  
	  if ( annot == NULL ) continue;
	  
	  O1 << "<Class name=\"" << annot->name << "\">\n"
	     << " <Description>" << annot->description << "</Description>\n";
	  
	  std::map<std::string, globals::atype_t>::const_iterator aa = annot->types.begin();
	  while ( aa != annot->types.end() )
	    {
	      O1 << "  <Variable type=\"" 
		 << globals::type_name[ aa->second ] 
		 << "\">" 
		 << aa->first 
		 << "</Variable>\n";
	      ++aa;
	    }
	  
	  O1 << "</Class>\n\n";
	  

	  //
	  // Enter all events in to a single table
	  //

	  annot_map_t::const_iterator ii = annot->interval_events.begin();
          while ( ii != annot->interval_events.end() )
            {
              const instance_idx_t & instance_idx = ii->first;
	      events.insert( instance_idx );
	      ++ii;
	    }
	  

	  //
	  // Next annotation/class header
	  // 
	}
      
      O1 << "</Classes>\n\n";  
      
      
      //
      // Loop over all annotation instances
      //
      
      O1 << "<Instances>\n\n";
  
  
      //   <Instance>   
      //      <Class>Recording Start Time</Class>
      //      <Name>Recording Start Time</Name>
      //      <Start>0</Start>
      //      <Duration>32820.0</Duration>
      //      <Channel>Optional channel label(s)</Channel>
      //      <Value var="name">0.123</Value>
      //      <Value var="name">0.123</Value>
      //      <Value var="name">0.123</Value>
      //    </Instance>


      std::set<instance_idx_t>::const_iterator ee = events.begin();
      while ( ee != events.end() )
	{

	  const instance_idx_t & instance_idx = *ee;
	  
	  const annot_t * annot = instance_idx.parent;

	  annot_map_t::const_iterator ii = annot->interval_events.find( instance_idx );
          if ( ii == annot->interval_events.end() ) {++ee; continue; }
	  instance_t * inst = ii->second;
	  
	  // skip if too short?
	  
	  interval_t interval = instance_idx.interval;

	  if ( has_min_dur )
	    {	      
	      const bool too_short = interval.duration_sec() < min_dur ; 	      
	      if ( too_short ) { ++ee; continue;}
	    }

	  // output
	  
	  O1 << "<Instance class=\"" << annot->name << "\">\n";
	  
	  if ( instance_idx.id != "." && instance_idx.id != "" ) 
	    O1 << " <Name>" << instance_idx.id << "</Name>\n";
	  
	  // adjuts by offset, if needed (ALIGN)
	  
	  if ( annot_offset )
	    {
	      if ( annot_offset_dir == +1 )
		{
                  interval.start += annot_offset;
                  interval.stop += annot_offset;
		}
	      else // from ALIGN
		{
		  if ( interval.start < annot_offset ) interval.start = 0;
		  else interval.start -= annot_offset;
		  if ( interval.stop < annot_offset ) interval.stop = 0;
		  else interval.stop -= annot_offset;
		}
	    }
	  
	  O1 << " <Start>" << interval.start_sec() << "</Start>\n"
	     << " <Duration>" << interval.duration_sec() << "</Duration>\n";
	  	  
	  if ( instance_idx.ch_str != "." && instance_idx.ch_str != "" )
	    O1 << " <Channel>" << instance_idx.ch_str << "</Channel>\n";

	  // name : instance_idx.id
	  // start : instance_idx.interval.start_sec()
	  // duration : instance_idx.interval.duration_sec()
	  
	  std::map<std::string,avar_t*>::const_iterator dd = inst->data.begin();
	  
	  while ( dd != inst->data.end() )
	    {
	      // var-name : dd->first
	      // value : 
	      
	      O1 << " <Value name=\"" << dd->first << "\">" 
		 << *dd->second 
		 << "</Value>\n"; 
	      ++dd;
	    }
	  
	  O1 << "</Instance>\n\n";


	  //
	  // Next instance (from events)
	  //

	  ++ee;

	}
      
      //
      // End of all annotation instatances
      //
      
      O1 << "</Instances>\n\n";

      //
      // Root node, close out the XML
      //
      
      O1 << "</Annotations>\n";
    }


  if ( ! xml_format ) 
    {

      //
      // Track all instances (in order)      
      //
      
      std::set<instance_idx_t> events;
      
      //
      // Annotation header
      //
      
      std::vector<std::string> anames = names();      
      
      for (int a=0;a<anames.size();a++)
	{
	  
	  if ( annots2write.size()
	       && annots2write.find( anames[a] ) == annots2write.end() )
	    continue;
	  
	  annot_t * annot = find( anames[a] );
	  
	  if ( annot == NULL ) continue;

	  // skip 'specials'
	  if ( annot->name == "start_hms" ) continue; 
	  if ( annot->name == "duration_hms" ) continue;
	  if ( annot->name == "duration_sec" ) continue;
	  if ( annot->name == "epoch_sec" ) continue; 

	  
	  // skip if empty
	  if ( annot->interval_events.size() == 0 ) continue;
	  
	  bool has_vars = annot->types.size() > 0 ;
	  
	  
	  if ( add_headers )
	    {
	  
	      // nb. ensure class name is quoted if contains `|` delimiter here 
	      O1 << "# " << Helper::quote_if( annot->name, '|' ) ; 
	      
	      if ( annot->description != "" )
		O1 << " | " << Helper::quote_if( annot->description, '|' );
	      else if ( has_vars ) // need a dummy description here 
		O1 << " | " << Helper::quote_if( annot->description, '|' );
	      
	      if ( has_vars ) 
		O1 << " |";
	      
	      std::map<std::string, globals::atype_t>::const_iterator aa = annot->types.begin();
	      while ( aa != annot->types.end() )
		{
		  O1 << " " << aa->first 
		     << "[" 
		     << globals::type_name[ aa->second ] 
		     << "]";	      
		  ++aa;
		}
	      
	      O1 << "\n";
	      
	    }


	  //
	  // Enter all events in to a single table
	  //
	  
	  annot_map_t::const_iterator ii = annot->interval_events.begin();
          while ( ii != annot->interval_events.end() )
            {
              const instance_idx_t & instance_idx = ii->first;
	      events.insert( instance_idx );
	      ++ii;
	    }
	  
	  //
	  // Next annotation class
	  //

	}

      //
      // dummy markers first 
      //

      if ( add_specials && add_headers )
	{
	  if ( start_hms != "." ) 
	    O1 << "# start_hms | EDF start time\n";
	  if ( duration_hms != "." )
	    O1 << "# duration_hms | EDF duration (hh:mm:ss)\n";
	  if ( duration_sec != 0 )
	    O1 << "# duration_sec | EDF duration (seconds)\n";
	  if ( epoch_sec != 0 )
	    O1 << "# epoch_sec | Default epoch duration (seconds)\n";
	}

      //
      // (optional, but nice to have) header row for data 
      //

      if ( ! minimal )
	O1 << "class" << "\t"
	   << "instance" << "\t"
	   << "channel" << "\t"
	   << "start" << "\t"
	   << "stop" << "\t"
	   << "meta" << "\n";
      
      //
      // Now, the data rows
      //

      // ensure 6-col format for .annot output

      if ( add_specials )
	{
	  if ( start_hms != "." )  
	    O1 << "start_hms\t" << start_hms << "\t.\t.\t.\t.\n";
	  if ( duration_hms != "." )
	    O1 << "duration_hms\t" << duration_hms << "\t.\t.\t.\t.\n";
	  if ( duration_sec != 0 )
	    O1 << "duration_sec\t" << duration_sec << "\t.\t.\t.\t.\n";
	  if ( epoch_sec != 0 )
	    O1 << "epoch_sec\t" << epoch_sec << "\t.\t.\t.\t.\n";
	}

      //
      // Loop over all annotation instances
      //
      
      std::set<instance_idx_t>::const_iterator ee = events.begin();
      while ( ee != events.end() )
	{

	  const instance_idx_t & instance_idx = *ee;
	  
	  const annot_t * annot = instance_idx.parent;

	  if ( annot == NULL ) { ++ee; continue; } 

	  if ( annots2write.size()
	       && annots2write.find( annot->name ) == annots2write.end() )
	    continue;
	  
	  
	  annot_map_t::const_iterator ii = annot->interval_events.find( instance_idx );
          if ( ii == annot->interval_events.end() )  { ++ee; continue; } 
	  instance_t * inst = ii->second;
	  
	  // skip special variables
	  
	  if ( annot->name == "start_hms" ) { ++ee; continue; } 
	  if ( annot->name == "duration_hms" ) { ++ee; continue; } 
	  if ( annot->name == "duration_sec" ) { ++ee; continue; } 
	  if ( annot->name == "epoch_sec" ) { ++ee; continue; } 
	  


	  //
	  // Get interval
	  //

	  interval_t interval = instance_idx.interval;
	  
	  // any duration reqs?

	  if ( has_min_dur )
	    {	      
	      const bool too_short = interval.duration_sec() < min_dur ; 	      
	      if ( too_short ) { ++ee; continue; } 
	    }
	  

	  //
	  // start/stop in seconds, with 4 d.p.
	  //
	  
	  // any re-ALIGNment ? 


          if ( annot_offset )
            {
	      if ( annot_offset_dir == 1 )
		{
		  interval.start += annot_offset;
		  interval.stop += annot_offset;
		}
	      else // from ALIGN
		{
		  if ( interval.start < annot_offset ) interval.start = 0;
		  else interval.start -= annot_offset;
		  if ( interval.stop < annot_offset ) interval.stop	= 0;
		  else interval.stop -= annot_offset;	      
		}
	    }
	  
	  // collapse from EDF+D to elapsed time in standard EDF ?

	  if ( collapse_disc && ! edf.header.continuous ) 
	    {
	      //	      std::cout << " pre  " << interval.start << " -- " << interval.stop << "\n";
	      interval = edf.timeline.collapse( interval );
	      //std::cout << " post " << interval.start << " -- " << interval.stop << "\n";
	      
	      // if the annotation doesn't completely fit in a region, skip it
	      if ( interval.start == 1LLU && interval.stop == 0LLU )
		{ ++ee; continue; }
	      // note -- if start and stop have a gap in the middle, this will still output
	      // (although duration will be shorter)
	      //   +++++ START ++++   |---GAP---|   +++++ STOP +++++
	      //  becomes
	      //   +++++ START +++++++++ STOP +++++
	    }
	  
	  // write an ... instead of the second time-poimt (for 0-duration intervals)
	  const bool add_ellipsis = globals::set_0dur_as_ellipsis && interval.start == interval.stop ;

	  //std::cout << " xx " << interval.start <<  " " << interval.stop << " " << interval.stop - interval.start << "   .... " << add_ellipsis << "\n";

	  
	  // output row (nb. no need to quote class, `|` allowed here
	  
	  O1 << annot->name << "\t";
	  
	  if ( instance_idx.id != "." && instance_idx.id != "" ) 
	    O1 << instance_idx.id << "\t";
	  else 
	    O1 << ".\t";

	  if ( instance_idx.ch_str != "." && instance_idx.ch_str != "" )
            O1 << instance_idx.ch_str << "\t";
          else
            O1 << ".\t";

	  
	  // write in hh:mm:ss format
	  if ( hms ) 
	    {

	      double tp1_sec =  interval.start / (double)globals::tp_1sec;
	      // add down to 1/1000th of a second
	      double tp1_extra = tp1_sec - (long)tp1_sec;
	      // but if we be round up to 1.000 then we need to add +1 to tp1_sec
	      if ( tp1_extra >= 0.9995 )
		{
		  ++tp1_sec;
		  tp1_extra = 0;
		}

	      // and get clock time
	      clocktime_t present1 = starttime;
	      present1.advance_seconds( tp1_sec );

	      // stop time
	      double tp2_sec =  interval.stop / (double)globals::tp_1sec;
	      double tp2_extra = tp2_sec - (long)tp2_sec;
	      if ( tp2_extra >= 0.9995 )
		{
		  ++tp2_sec;
                  tp2_extra = 0;
		}
	      clocktime_t present2 = starttime;
	      present2.advance_seconds( tp2_sec );
		  
	      // std::cout << " tp = " << interval.as_string() << "\n";
	      // std::cout << " times = " << interval.start_sec() << " " << interval.stop_sec() << "\t"
	      // 		<< present1.as_string(':')  << " -- " << present2.as_string(':') << "\n";

	      // std::cout << " fl = " << tp1_sec << " " << tp1_extra << " ---- " << tp2_sec << " " << tp2_extra << "\n";
	      
	      // add dd-mm-yy-hh:mm:ss
	      
	      // hh:mm:ss.ssss
	      if ( globals::time_format_dp ) 
		O1 << ( dhms ? present1.as_datetime_string(':') : present1.as_string(':') )
		   << Helper::dbl2str_fixed( tp1_extra , globals::time_format_dp  ).substr(1)
		   << "\t"
		   << ( add_ellipsis ? "..." : ( dhms ? present2.as_datetime_string(':') : present2.as_string(':') ) + Helper::dbl2str_fixed( tp2_extra , globals::time_format_dp  ).substr(1) ) ;
	      else // or truncate to hh:mm:ss
		O1 << ( dhms ? present1.as_datetime_string(':') : present1.as_string(':') ) << "\t"
		   << ( add_ellipsis ? "..." : ( dhms ? present2.as_datetime_string(':') : present2.as_string(':') ) ) ;
	      
	    }
	  else // write as elapsed seconds
	    {
	      O1 << Helper::dbl2str( interval.start_sec() , globals::time_format_dp ) << "\t"
		 << ( add_ellipsis ? "..." : Helper::dbl2str( interval.stop_sec() , globals::time_format_dp ) ) ;
	    }

	  if ( inst->data.size() == 0 ) 
	    O1 << "\t.";
	  else
	    {
	      O1 << "\t";

	      std::map<std::string,avar_t*>::const_iterator dd = inst->data.begin();
	      
	      while ( dd != inst->data.end() )
		{
		  // pipe-delimiter
		  if ( dd != inst->data.begin() ) O1 << "|";
		  
		  // meta-data value, always key/value pairing
		  // as there may be missing data
		  // if the meta-data is string and contains a pipe, then
		  // we need to quote this
		  
		  std::stringstream ss;
		  ss << *dd->second;

		  O1 << dd->first << "=" << Helper::quote_spaced( Helper::quote_if( ss.str() , '|' , '=' ) );
		  ++dd;
		}
	    }

	  O1 << "\n";

	  //
	  // Next instance/event
	  //

	  ++ee;
	  
	}
      
    }
        
  //
  // All done
  //
  
  O1.close();
    

}




bool annot_t::loadxml_luna( const std::string & filename , edf_t * edf )
{
  
  XML xml( filename );
  
  if ( ! xml.valid() ) Helper::halt( "invalid annotation file: " + filename );

  //
  // Annotation classes
  //
  
  std::vector<element_t*> classes = xml.children( "Classes" );

  std::set<std::string> added;
  
  for (int i=0;i<classes.size();i++)
    {
      
      element_t * cls = classes[i];
      
      if ( ! Helper::iequals( cls->name , "Class" ) ) continue;
      
      std::string cls_name = cls->attr.value( "name" );
      
      //
      // alias remapping?
      //

      std::string original_label = cls_name;
      cls_name = nsrr_t::remap( cls_name );
      if ( cls_name == "" ) continue;
      
      //
      // ignore this annotation?
      //
      
      if ( globals::specified_annots.size() > 0 && 
	   globals::specified_annots.find( cls_name )
	   == globals::specified_annots.end() ) continue;

      //
      // track aliasing
      //
      
      if ( cls_name != original_label )
	edf->timeline.annotations.aliasing[ cls_name ] = original_label ;
      
      std::string desc = "";
      std::map<std::string,std::string> atypes;

      std::vector<element_t*> kids = cls->child;
      
      for (int j=0; j<kids.size(); j++)
        {
	  
          const std::string & key = kids[j]->name;

	  if ( key == "Description" ) 
	    {
	      desc = kids[j]->value;
	    }
	  else if ( key == "Variable" ) 
	    {
	      atypes[ kids[j]->value ] = kids[j]->attr.value( "type" );
	    }
	  
        }

      
      //       <Class name="a3">
      // 	  <Name>a3</Name>
      // 	  <Description>This annotation also specifies meta-data types</Description>
      // 	  <Variable type="txt">val1</Variable>
      // 	  <Variable type="num">val2</Variable>
      // 	  <Variable type="bool">val3</Variable>
      //       </Class>
      
      //
      // add this annotation
      //
      
      annot_t * a = edf->timeline.annotations.add( cls_name );
      
      a->description = desc;
      a->file = filename;
      a->type = globals::A_FLAG_T; // not expecting any meta-data (unless changed below)

      std::map<std::string,std::string>::const_iterator aa = atypes.begin();
      while ( aa != atypes.end() )
	{
	  // if a recognizable type, add
	  if ( globals::name_type.find( aa->second ) != globals::name_type.end() )
	    a->types[ aa->first ] = globals::name_type[ aa->second ];
	  ++aa;
	}
      
      // as with .annot files; if only one variable, set annot_t equal to the one instance type
      // otherwise, set as A_NULL_T ; in practice, don't think we'll ever use annot_t::type 
      // i.e. will always use annot_t::atypes[]

      if ( a->types.size() == 1 ) a->type = a->types.begin()->second;
      else if ( a->type > 1 ) a->type = globals::A_NULL_T; 
      // i.e. multiple variables/types set, so set overall one to null

    }


  //
  // Annotation Instances
  //
  
  std::vector<element_t*> instances = xml.children( "Instances" );
  

  //
  // First pass through all instances
  //
  
  for (int i=0; i<instances.size(); i++) 
    {
      
      element_t * ii = instances[i];
      
      std::string cls_name = ii->attr.value( "class" );

      //
      // alias remapping?
      //

      std::string original_label = cls_name;
      cls_name = nsrr_t::remap( cls_name );
      if ( cls_name == "" ) continue;

      
      //
      // ignore this annotation?
      //
           
      if ( globals::specified_annots.size() > 0 && 
	   globals::specified_annots.find( cls_name ) 
	   == globals::specified_annots.end() ) continue;

      if ( cls_name != original_label )
	edf->timeline.annotations.aliasing[ cls_name ] = original_label ;

    
      //
      // get a pointer to this class
      //

      annot_t * a = edf->timeline.annotations.find( cls_name );
      
      if ( a == NULL ) continue;
      
      // pull information for this instance:
      
      element_t * name     = (*ii)( "Name" );
      
      element_t * start    = (*ii)( "Start" );
      
      element_t * duration = (*ii)( "Duration" );

      element_t * channel  = (*ii)( "Channel" );
      
      //
      // Get time interval
      //

      double dbl_start = 0 , dbl_dur = 0 , dbl_stop = 0;
    
      if ( ! Helper::str2dbl( start->value , &dbl_start ) )
	Helper::halt( "invalid interval: " + start->value );
    
      if ( ! Helper::str2dbl( duration->value , &dbl_dur ) ) 
	Helper::halt( "invalid interval: " +  duration->value );
      
      dbl_stop = dbl_start + dbl_dur; 
	      
      if ( dbl_start < 0 )
	{
	  //std::cout << " S3 " << dbl_start << "\n";
	  Helper::halt( filename + " contains row(s) with negative time points" ) ;
	}
      
      if ( dbl_dur < 0 )
	{
	  //std::cout << " S4 " << dbl_dur << "\n";
	  Helper::halt( filename + " contains row(s) with negative durations" );
	}
      
      // convert to uint64_t time-point units
      
      interval_t interval;

      interval.start = Helper::sec2tp( dbl_start );
      
      // assume stop is already specified as 1 past the end, e.g. 30 60
      // *unless* it is a single point, e.g. 5 5 
      // which is handled below
      
      interval.stop  = Helper::sec2tp( dbl_stop );

      
      // given interval encoding, we always want one past the end
      // if a single time-point given (0 duration)
      //
      // otherwise, assume 30 second duration means up to 
      // but not including 30 .. i..e  0-30   30-60   60-90 
      // in each case, start + duration is the correct value

      if ( interval.start == interval.stop ) ++interval.stop;

      //
      // Create the instance; only add the instance ID/Name if it is different from the class ID
      //
    
      instance_t * instance = a->add( name ? ( name->value != cls_name ? name->value : "." ) : "." , 
				      interval ,
				      channel ? channel->value : "." );

      //
      // Add any additional data members
      //

      std::vector<element_t*> kids = ii->child;
      
      for (int j=0; j<kids.size(); j++) 
	{
	  
	  const std::string & key = kids[j]->name;
	  
	  if ( key == "Value" ) 
	    {
	      std::string var = kids[j]->attr.value( "name" );
	      std::string val = kids[j]->value;
	      
	      if ( a->types.find( var ) != a->types.end() ) 
		{
		  
		  globals::atype_t t = a->types[ var ];
	      
		  if ( t == globals::A_FLAG_T ) 
		    {
		      instance->set( var );
		    }
		  
		  else if ( t == globals::A_MASK_T )
		    {
		      if ( var != "." )
			{
			  // accepts F and T as well as long forms (false, true)
			  instance->set_mask( var , Helper::yesno( val ) );
			}
		    }
		  
		  else if ( t == globals::A_BOOL_T )
		    {
		      if ( val != "." )
			{
			  // accepts F and T as well as long forms (false, true)
			  instance->set( var , Helper::yesno( val ) );
			}
		    }
		  
		  else if ( t == globals::A_INT_T )
		    {
		      int value = 0;
		      if ( ! Helper::str2int( val , &value ) )
			Helper::halt( "bad numeric value in " + filename );
		      instance->set( var , value );
		    }

		  else if ( t == globals::A_DBL_T )
		    {
		      double value = 0;
		      
		      if ( Helper::str2dbl( val , &value ) )		    
			instance->set( var , value );
		      else
			if ( var != "." && var != "NA" ) 
			  Helper::halt( "bad numeric value in " + filename );		  
		    }
		  
		  else if ( t == globals::A_TXT_T )
		    {
		      instance->set( var , val );
		    }
		  
		}

	    } // added this data member
	  
	}
     
      //
      // Next instance
      //
    }
  
  return true;
}



void annotation_set_t::clear() 
{ 
  std::map<std::string,annot_t*>::iterator ii = annots.begin();
  while ( ii != annots.end() ) 
    {
      delete ii->second;
      ++ii;
    }
  
  annots.clear(); 
  
  start_ct.reset();
  
  start_hms = ".";
  
  duration_hms = ".";
  
  duration_sec = 0 ;
  
  epoch_sec = 0 ; 

  annot_offset = 0LLU;

  annot_offset_dir = -1;
}


//
// Initiate annotation set from EDF, to seed with a few key values
//

void annotation_set_t::set( edf_t * edf ) 
{
  // populate start_hms
  // duration_hms,
  // duration_sec
  // and epoch_sec

  
  if ( edf != NULL )
    {
      
      duration_sec = edf->header.nr_all * edf->header.record_duration ;
      
      duration_hms = Helper::timestring( globals::tp_1sec * duration_sec , '.' , false ); // no fractional seconds
      
      clocktime_t etime( edf->header.starttime );
      
      if ( etime.valid )
	{

	  start_ct = etime;
	  start_hms = edf->header.starttime ;
	  
	  // double time_hrs = ( edf->timeline.last_time_point_tp * globals::tp_duration ) / 3600.0 ;
	  // etime.advance( time_hrs );	  
	}
      
      epoch_sec = edf->timeline.epoched() ?
	edf->timeline.epoch_length() :
	globals::default_epoch_len ; 
      
    }
  
}


//
// Convert from EDF Annotations track(s) to Luna-format annotations
//

annot_t * annotation_set_t::from_EDF( edf_t & edf , edfz_t * edfz )
{
  
  if ( ! globals::skip_edf_annots ) 
    logger << "  extracting 'EDF Annotations' track "
	   << ( edfz == NULL ? "from EDF+" : "from EDFZ .idx" )
	   << "\n";  
  else
    logger << "  extracting only EDF+D time-track 'EDF Annotations' track\n";  
  
  // create a single annotation (or bind to it, if it already exists)
  // by default, this is edf_annot_t and the entries here are added as the instance ID

  // however, we can specify that certain EDF annotations are entered as a class
  //  these will be remapped etc as above;  if edf-annot-class-all=T, then /all/
  //  EDF+ annotations are added at the class level
  
  annot_t * a =  NULL;

  // only need edf_annot_t if edf-annot-class-all=F
  if ( ! nsrr_t::all_edf_class )
    {
      a = edf.timeline.annotations.add( globals::edf_annot_label );
      a->name = globals::edf_annot_label;
      a->description = "EDF Annotations";
      a->file = edf.filename;
      a->type = globals::A_FLAG_T; 
    }
  
  // if we need to expand 0-duration stages
  uint64_t epoch_len = globals::tp_1sec *
    ( edf.timeline.epoch_len_tp_uint64_t() == 0 ?
      globals::default_epoch_len :
      edf.timeline.epoch_len_tp_uint64_t() );


  //
  // when reading a typical EDF+ we need to read from disk
  // but when reading a compressed EDF+, the annotations will
  // already be duplicated in the .idx, and so we can pull
  // directly from that, which is much quicker
  //

  //
  // Parse EDFZ (from index)
  //

  if ( edfz != NULL )
    {
      int r = edf.timeline.first_record();
      while ( r != -1 )
	{
	  std::string s = edfz->get_annots( r );

	  if ( s == "." )
	    {
	      // skip to next record
	      r = edf.timeline.next_record( r );
	      continue;
	    }

	  // quoted, comma-delimited
	  // "onset|dur|text","onset|dur|text"
	  std::vector<std::string> tok = Helper::quoted_parse( s , "," );
	  
	  for (int j=0; j<tok.size(); j++)
	    {
	      // track that this has actual EDF Annotations 
	      edf.has_edf_annots = true;
	      
	      std::vector<std::string> tok2 = Helper::parse( Helper::unquote( tok[j] ) , "|" );
	      if ( tok2.size() < 3 ) Helper::halt( "bad format for EDF .idx annots (vec-len):\n" + tok[j] );
	      
	      double onset = 0, dur = 0;
	      if ( ! Helper::str2dbl( tok2[0] , &onset ) )
		Helper::halt( "bad format for EDF .idx annots (onset):\n" + tok[j]  );

	      if ( ! Helper::str2dbl( tok2[1] , &dur ) )
		Helper::halt( "bad format for EDF .idx annots (dur):\n" + tok[j]  );
	      
	      std::string txt = tok2[2];
	      
	      // add this annotation (clunky, but keep this code in sync w/
	      // what we do below when reading from the EDF+ directly)
	      
	      uint64_t start_tp = Helper::sec2tp( onset );
	      
	      uint64_t dur_tp = Helper::sec2tp( dur );
		      
	      // stop is one past the end 
	      // NOTE: zero-lengh annot is [a,a),
	      uint64_t stop_tp  = start_tp + dur_tp ;
	      
	      // get the annotation label
	      std::string aname = Helper::trim( txt );
	      
	      // // sanitize?
	      // if ( globals::sanitize_everything )
	      // 	aname = Helper::sanitize( aname );
	      
	      // do any remapping
	      const std::string tname = nsrr_t::remap( aname );

	      // track aliasing?
	      if ( tname != aname )
		edf.timeline.annotations.aliasing[ tname ] = aname;	      

	      aname = tname;
	      
	      
	      // fix stage duration (if 0-dur point)?  (unless
	      // adding ellipsis, i.e. here change points
	      // might not map to even epochs... 30, 90, 30,
	      // 180, etc... )
	      
	      if ( globals::sleep_stage_assume_epoch_duration
		   && globals::is_stage_annotation( aname )
		   && ( ! globals::set_0dur_as_ellipsis )
		   && start_tp == stop_tp )
		stop_tp += epoch_len;
	      
	      // make interval
	      interval_t interval( start_tp , stop_tp );
	      
	      // is this a class?
	      bool edf_class =  nsrr_t::as_edf_class( aname );
	      
	      if ( aname != "" )
		{			  
		  // add as standard edf_annot_t ?
		  if ( ! edf_class )
		    {
		      if ( ! nsrr_t::whitelist )
			{
			  instance_t * instance = a->add( aname , interval , "." );
			  // track how many annotations we add
			  edf.aoccur[ globals::edf_annot_label ]++;
			}
		    }
		  else // ... else add as a class
		    {
		      // add as a new class
		      // (no meta-info)
		      annot_t * a = edf.timeline.annotations.add( aname );
		      instance_t * instance = a->add( "." , interval , "." );
		      edf.aoccur[ aname ]++;
		    }
		}
	    }

	  
	  // next record
	  r = edf.timeline.next_record( r );		
	}

      // all done
      return a;
    }
  
  
  //
  // read from main file
  //
  
  int r = edf.timeline.first_record();
  
  while ( r != -1 )
    {
      
      for ( int s = 0 ; s < edf.header.ns; s ++ )
	{
	  
	  if ( edf.header.is_annotation_channel( s ) )
	    {	      
	      
	      tal_t t = edf.tal( s , r );
	      
	      // store (for use if WRITE edfz is later called,
	      //  i.e. to populate the .idx)
	      
	      edf.edf_annots[ r ] = t.export_annots();
	      
	      //std::cout << " edf-annot s,r = " << s << " " << r << "\n" << t << "\n";
	      
	      const int na = t.size();
	      
	      for (int i=0; i<na; i++)
		{
		  
		  tal_element_t & te = t.d[i];
		  
		  if ( te.name != globals::edf_timetrack_label )
		    {
		      
		      // track that this has actual EDF Annotations 
		      edf.has_edf_annots = true;
		      
		      uint64_t start_tp = Helper::sec2tp( te.onset );
		      
		      uint64_t dur_tp = Helper::sec2tp( te.duration );
		      
		      // stop is one past the end 
		      // NOTE: zero-lengh annot is [a,a),
		      uint64_t stop_tp  = start_tp + dur_tp ;
		      
		      // get the annotation label
		      std::string aname = Helper::trim( te.name );

		      // sanitize? (done in remap() )
		      //	      if ( globals::sanitize_everything )
		      //   aname = Helper::sanitize( aname );
		      
		      // do any remapping                                                                                                                         
		      const std::string tname = nsrr_t::remap( aname );
		      
		      // track aliasing?                                                                                                                          
		      if ( tname != aname )
			edf.timeline.annotations.aliasing[ tname ] = aname;
		      
		      aname = tname;
		      
		     	      
		      // fix stage duration (if 0-dur point)?  (unless
		      // adding ellipsis, i.e. here change points
		      // might not map to even epochs... 30, 90, 30,
		      // 180, etc... )
		      
		      if ( globals::sleep_stage_assume_epoch_duration
			   && globals::is_stage_annotation( aname )
			   && ( ! globals::set_0dur_as_ellipsis )
			   && start_tp == stop_tp )
			stop_tp += epoch_len;
				      
		      // make interval
		      interval_t interval( start_tp , stop_tp );
		      
		      // is this a class?
		      bool edf_class =  nsrr_t::as_edf_class( aname );
		      // if not, do we ignore? : nsrr_t::only_add_named_EDF_annots 
		      
		      if ( aname != "" )
			{			  
			  // add as standard edf_annot_t ?
			  if ( ! edf_class )
			    {
			      if ( ! nsrr_t::whitelist )
				{
				  instance_t * instance = a->add( aname , interval , "." );
				  // track how many annotations we add
				  edf.aoccur[ globals::edf_annot_label ]++;
				}
			    }
			  else // ... else add as a class
			    {
			      // add as a new class
			      // (no meta-info)
			      annot_t * a = edf.timeline.annotations.add( aname );
			      instance_t * instance = a->add( "." , interval , "." );
			      edf.aoccur[ aname ]++;
			    }
			}
		    }
		  
		}
	      
	    } 
	  
	} // next signal

      r = edf.timeline.next_record( r );
      
    } // next record
  
  return a;
}



      
// of a list of annotations, find the first start time 
// i.e. for aligning staging with a start that is not 0-seconds

double annotation_set_t::first(const std::vector<std::string> & requested ) const
{
  std::set<double> starts;
  
  for (int a=0; a < requested.size(); a++)
    {
      // find this annotation
      annot_t * annot = find( requested[a] );
      if ( annot == NULL ) continue;
      // sorted first by time
      annot_map_t::const_iterator ii = annot->interval_events.begin();
      if (  ii == annot->interval_events.end() ) continue;
      starts.insert( ii->first.interval.start_sec() );
    }
  
  // none found: return 0
  if ( starts.size() == 0 ) return 0;
  
  // return smallest
  return *starts.begin();
  
}

uint64_t annotation_set_t::first_in_interval( const std::vector<std::string> & requested , 
					      const interval_t & range ) const
{
    
  // within 'range', get first to of the annots in 'requested'
  //  - this is initially for making epochs align well for EDF+D
  //  - i.e. if we have staging and annots
  
  //  annot_map_t extract( const interval_t & window );
  
  std::set<uint64_t> starts;
  
  for (int a=0; a < requested.size(); a++)
    {
      // find this annotation
      annot_t * annot = find( requested[a] );
      if ( annot == NULL ) continue;
      
      // get annots in this window only
      annot_map_t amap = annot->extract( range );
      
      // get the first 
      annot_map_t::const_iterator ii = annot->interval_events.begin();
      if (  ii == annot->interval_events.end() ) continue;
      starts.insert( ii->first.interval.start );
    }
  
  // none found: return 0 (i.e. start of range)
  if ( starts.size() == 0 ) return range.start;
  
  // return smallest
  return *starts.begin();

}


std::set<uint64_t> annotation_set_t::starts( const std::vector<std::string> & requested , uint64_t dur ) const
{

  // get start points from these requested epochs;
  // but add in extra start points for each 'dur' period within that annotation
  // i.e. for epoch-alignment, this handles the case of annotations that are >1 multiples
  // of the epoch size (e.g. 90s REM)  and adds in possible starts at 0, 30, 60 s
  
  std::set<uint64_t> sts;
  
  for (int a=0; a < requested.size(); a++)
    {
      // find this annotation                                                                                                                
      annot_t * annot = find( requested[a] );
      if ( annot == NULL ) continue;
      
      annot_map_t::const_iterator ii = annot->interval_events.begin();
      while ( ii != annot->interval_events.end() )
	{

	  if ( dur == 0 )
	    sts.insert( ii->first.interval.start );
	  else
	    {
	      	  
	      uint64_t pos = ii->first.interval.start ;
	      uint64_t end = ii->first.interval.stop ;
	      while ( 1 )
		{
		  if ( pos + dur <= end )
		    {
		      sts.insert( pos );
		      pos += dur;
		    }
		  else break;
		}
	    }
	  ++ii;
	}
    }
  return sts;
}


void annotation_set_t::extend( param_t & param )
{
  // for a set of annotations, if 0-duration, extend it until the
  // next set instance: e.g. for staging

  if ( ! param.has( "annots" ) )
    Helper::halt( "requires annots argument" ) ;
  std::set<std::string> a = param.strset( "annots" );

  
}


int annotation_set_t::remap( const std::vector<std::string> & files , int remap_field , bool remap_spaces , bool verbose )
{

  if ( verbose )
    {
      logger << "  REMAP annotations:\n";
      if ( remap_spaces )  logger << "   - allowing space-delimited & tab-delimited fields\n";
      else logger << "   - only allowing tab-delimited fields\n";
      
      if      ( remap_field == 0 ) logger << "   - assuming no 'remap' column 1 fields\n";
      else if ( remap_field == 1 ) logger << "   - assuming 'remap' column 1 fields present\n";
      else if ( remap_field == 2 ) logger << "   - optionally allowing but not requiring 'remap' columns\n";
    }	
  
  // clear all prior 'aliasing' info
  aliasing.clear();
  
  int mapped = 0;

  // remap_field
  //   0  
  //  primary  second|third                 [ remap_field == F ]

  //   1 
  //  remap     primary|second|third        [ remap_field == T ]

  //   2 -- moonlight mode
  //  could be either -- scan first field to determine
  //  ignore 'remap' and also 'nsrr-remap'
  

  // if remap_spaces == T , then allow space-delimiters (i.e. assumes that
  // annots w/ spaces are quoted)
  
  std::map<std::string,std::string> old2new;
  
  for (int fi=0; fi<files.size(); fi++)
    {
      const std::string fname = Helper::expand( files[fi] );
      
      if ( ! Helper::fileExists( fname ) )
	Helper::halt( "could not find " + fname );

      std::ifstream IN1( fname.c_str() , std::ios::in );
      while ( ! IN1.eof() )
	{
	  std::string x;
	  Helper::safe_getline( IN1 , x );
          if ( IN1.eof() ) break;
	  x = Helper::trim( x , ' ' , '\t' );
	  
	  if ( x == "" ) continue;
	  if ( x[0] == '%' ) continue;
	  
	  std::vector<std::string> tok = remap_spaces ? Helper::quoted_parse( x, " \t" ) : Helper::parse( x , "\t" );

	  // requires 'remap' field?
	  if ( remap_field == 1 && ! Helper::iequals( tok[0] , "remap" ) )
	    continue;

	  // skip special term in NSRR annot files
	  if ( Helper::iequals( tok[0] , "nsrr-remap" ) )
	    continue;
	  
	  // allow
	  //   remap    pri|sec|third    [ remap_field == T ] 
	  //   pri      sec|third
	  //   pri|sec|third

	  bool has_remap = remap_field == 1 || ( remap_field == 2 && Helper::iequals( tok[0] , "remap" )  ) ;
	  
	  if ( tok.size() > 2 ) Helper::halt( "bad format: " + x );
	  
	  if ( has_remap && tok.size() != 2 ) Helper::halt( "bad format: " + x );
	  
	  std::string tok1 = has_remap ? tok[1] : ( tok.size() == 1 ? tok[0] : tok[0] + "|" + tok[1] ) ; 

	  // swap out spaces?
	  if ( globals::replace_annot_spaces )
	    tok1 = Helper::search_replace( tok1 , ' ' , globals::space_replacement );

	  // sanitize?
	  if ( globals::sanitize_everything )
	    {
	      if ( globals::replace_annot_spaces )
		tok1 = Helper::sanitize( tok1 );
	      else // allow spaces in a sanitized version still, and keeps | and "
		tok1 = Helper::sanitize( tok1 , ' ' );
	    }
	  
	  std::vector<std::string> tok2 = Helper::quoted_parse( tok1 , "|" );
	  
	  if ( tok2.size() < 2 ) Helper::halt( "problem with line: " + x );

	  // trims spaces and underscores
	  std::string snew = Helper::trim( Helper::unquote( tok2[0] ) , '_' );
	  for (int j=1; j<tok2.size(); j++)
	    {
	      std::string sorig = Helper::trim( Helper::unquote( tok2[j] ) , '_' );
	      sorig = Helper::squash( Helper::squash( sorig , ' ' ) , '_' );
	      old2new[ sorig ] = snew ;
	      if ( verbose )
		logger << "  adding mapping [" << sorig << "] --> [" << snew << "]\n";
	    }
	  
	}
    }

  //
  // requires a one-to-one mapping, i.e. cannot merge;  check this quickly here
  //
  
  std::map<std::string,std::string> target2orig;
  
  std::map<std::string,std::string>::const_iterator aa = old2new.begin();
  while ( aa != old2new.end() )
    {
      if ( aa->first == "start_hms" || aa->second == "start_hms" )
	Helper::halt( "cannot remap to a special annotation term: start_hms" );

      if ( aa->first == "duration_hms" || aa->second == "duration_hms" )
	Helper::halt( "cannot remap to a special annotation term: duration_hms" );

      if ( aa->first == "duration_sec" || aa->second == "duration_sec" )
	Helper::halt( "cannot remap to a special annotation term: duration_sec" );
      
      if ( aa->first == "epoch_sec" || aa->second == "epoch_sec" )
	Helper::halt( "cannot remap to a special annotation term: epoch_sec" );

      if ( aa->first == "annot_offset" || aa->second == "annot_offset" )
	Helper::halt( "cannot remap to a special annotation term: annot_offset" );

      // does the original term actually exist?
      if ( annots.find( aa->first ) != annots.end() )
	{
	  // check that the new term does *not* already exist
	  if ( annots.find( aa->second ) != annots.end() )
	    Helper::halt( "cannot map to an existing term: " + aa->first + " " + aa->second );

	  // check that another original has not already pointed to the same new term, ( and exists in the data)
	  if ( target2orig.find( aa->second ) != target2orig.end() )
	    Helper::halt( "cannot map multiple existing terms to the same target: "
			  + aa->first + " and " + target2orig[ aa->second ] + " --> " + aa->second );
	  
	  // otherwise, this will be okay to map
	  target2orig[ aa->second ] = aa->first ; 
	  
	}	
      
      ++aa;
    }


  //
  // Now do the actual remapping
  //

  // MAIN:: 
  //  std::map<std::string,annot_t*> annots;
  // update :  track alias swaps for this person
  //    std::map<std::string,std::string> aliasing;

  std::map<std::string,std::string>::const_iterator tt = target2orig.begin();
  while ( tt != target2orig.end() )
    {
      // original : tt->second
      // new      : tt->first
      
      logger << "  remapping " << tt->second << " to " << tt->first << "\n";
      
      // 1) copy index in pointer map
      annots[ tt->first ] = annots[  tt->second ];

      // 2) erase old version
      annots.erase( annots.find(  tt->second ) );

      // 3) update annot_t name
      annot_t * a = annots[ tt->first ];
      a->name = tt->first;

      // 4) track in aliasing map [ new -> old ] 
      //    i.e. so results will show in an ALIASES command
      aliasing[ tt->first ] = tt->second ;
      
      ++mapped;
      
      ++tt;
    }
  
  return mapped;
}

