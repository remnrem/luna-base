
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
#include "main.h"
#include "edf/edf.h"
#include "helper/helper.h"
#include "defs/defs.h"
#include "tinyxml/xmlreader.h"
#include "db/db.h"

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>

extern writer_t writer;

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


instance_t * annot_t::add( const std::string & id , const interval_t & interval )
{
  
  instance_t * instance = new instance_t ;
  
  // track (for clean-up)
  all_instances.insert( instance );
  
  interval_events[ instance_idx_t( interval , id ) ] = instance; 
  
  return instance; 
  
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
  if ( dd->second == NULL ) return; // flag, so no storage set
  tracker.erase( dd->second ); // otherwise, need to clear old storage
  data.erase( dd );  // and erase from list
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

void instance_t::set( const std::string & name , const double d )
{
  check( name );
  avar_t * a = new double_avar_t( d ) ;
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
  // Original .annot file format
  //
      
  // logger << " attaching annotation file " << f << "\n";

  //
  // record filename
  //

  
    
  std::ifstream FIN( f.c_str() , std::ios::in );

  
  // New annotation format
  
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
  
  // then rows are either interval or epoch-based 
  
  // name  id1  sec1  sec2  { vars }
  // name  id   e30:1       { vars } 
  
  // for now, assume e:1   means 30-second epoch, no overlap  [ hard-code this ] 
  
  bool epoched = false;

  int line_count = 0;

  std::map<std::string,annot_t*> annot_map;
  
  std::map<annot_t*,std::vector<std::string> > cols;
    
  while ( ! FIN.eof() )
    {
      
      if ( FIN.bad() ) continue;
      std::string line;      
      std::getline( FIN , line );      
      if ( FIN.eof() || line == "" ) continue;
      
      //
      // header or data row?
      //

   
      if ( line[0] == '#' ) 
	{
	  
	  // skip initial # here
	  std::vector<std::string> tok = Helper::parse( line.substr(1) , "|" );
	  
	  if ( tok.size() < 1 || tok.size() > 3 ) Helper::halt( "bad header for format\n" + line );

	  //
	  // Get the name and ID 
	  //

	  std::string name = Helper::trim( tok[0] );
	  
	  //
	  // skip this annotation
	  //

	  if ( globals::specified_annots.size() > 0 && 
	       globals::specified_annots.find( name ) == globals::specified_annots.end() ) 
	    continue;
	    
	  
	  //
	  // otherwise, create the annotation if it doesn't exist
	  //

	  annot_t * a = parent_edf.timeline.annotations.add( name );


	  //
	  // store a temporary lookup table
	  //

	  annot_map[ name ] = a;
	  
	  //
	  // Other details
	  //
	  
	  a->description = tok.size() >= 2 ? tok[1] : name ; 
	  
	  a->file = f;
	  
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
		    t = globals::A_FLAG_T;
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
		  
		} // next column for this annotation

	    }
	  
	}


      //
      // otherwise, assume this is a data row
      //

      else 
	{

	  // data-rows are tab-delimited
	  
	  std::vector<std::string> tok = Helper::parse( line , "\t" );
	  
	  if ( tok.size() == 0 ) continue; 
	  
	  // was this annotation specified? 
	  
	  std::map<std::string,annot_t*>::iterator aa = annot_map.find( tok[0] );
	  
	  if ( aa == annot_map.end() ) 
	    Helper::halt( "annotation " + tok[0] + " not in header of " + f );
	  
	  annot_t * a = aa->second; 			
	  
	  std::string id = tok[1];

	  // epoch or interval? 

	  if ( tok.size() < 3 ) Helper::halt( "bad line format, need at least 3 cols:\n" + line );
	  
	  bool eline = tok[2][0] == 'e' ;
	  
	  // expected # of cols
	  
	  const int expected = 2 + ( eline ? 1 : 2 ) + a->types.size() ; 
	  
	  if ( tok.size() != expected ) 
	    Helper::halt( "bad line format: saw " + Helper::int2str( (int)tok.size() ) + " cols but expecting " 
			  + Helper::int2str( expected) + " for " + a->name + "\n" + line );
	  
	  
	  interval_t interval;
	  
	  if ( eline )
	    {
	      
	      // 2 e:1        assumes 30
	      // 3 e:30:1     assumes no overlap
	      // 4 e:15:5:4   specifies all
	      
	      std::vector<std::string> tok2 = Helper::parse( tok[2] , ":" );
	      
	      if ( tok2.size() < 2 || tok2.size() > 4 ) 
		Helper::halt( "bad epoch specification, expecting e:1, e:30:1, e:30:30:1, etc" );
	      
	      if ( tok2[0] != "e" ) 
		Helper::halt( "bad epoch specification, expecting e:1, e:30:1, e:30:30:1, etc" );
	      
	      int epoch_length = 30 ;
	      int epoch_increment = 30;
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

	      uint64_t epoch_length_tp = epoch_length * globals::tp_1sec;
	      
	      uint64_t epoch_increment_tp = epoch_increment * globals::tp_1sec;
	      
	      // set interval from current line
	      // last point is defined as point *past* end of interval
	      
	      interval.start = epoch_increment_tp * (epoch-1);
	      interval.stop  = interval.start + epoch_length_tp;

	      //interval.stop  = interval.start + epoch_length_tp - 1;
	      
	    }
	  else // an INTERVAL
	    {
	      
	      // read as seconds and convert to uint64_t 
	      double dbl_start = 0 , dbl_stop = 0;
	      
	      if ( ! Helper::str2dbl( tok[2] , &dbl_start ) )
		Helper::halt( "invalid interval: " + line );
	      
	      if ( ! Helper::str2dbl( tok[3] , &dbl_stop ) ) 
		Helper::halt( "invalid interval: " + line );
	      
	      interval.start = globals::tp_1sec * dbl_start;
	      interval.stop  = globals::tp_1sec * dbl_stop;	      
	      
	      //	      std::cout << "reading " << interval.start << " - " << interval.stop << "\n";
	      
	      if ( interval.start > interval.stop )
		Helper::halt( "invalid interval: " + line );
	      
	    }
	  
	  
	  //
	  // for a FLAG for this annotation (with same name as primary annotaton)
	  //
	  
	  instance_t * instance = a->add( id , interval );
	  
	  
	  // track how many annotations we add
	  parent_edf.aoccur[ a->name ]++;
	  
	  
	  const int n = tok.size();
	          

	  //
	  // Also add any other columns (as separate events under this same annotation)
	  //
	  
	  for (int j = (eline ? 3 : 4 ) ;j<n; j++)
	    {

	      const int idx = j - (eline ? 3 : 4 );
	      
	      const std::string & label = cols[a][idx];	      

	      globals::atype_t t = a->types[label];
	      
	      if ( t == globals::A_FLAG_T ) 
		{
		  instance->set( label );
		}
	      
	      else if ( t == globals::A_BOOL_T )
		{
		  bool value = ! ( tok[j] == "0" || tok[j] == "F" || tok[j] == "." );
		  instance->set( label , value );
		}
	      
	      else if ( t == globals::A_INT_T )
		{
		  int value = 0;
		  if ( ! Helper::str2int( tok[j] , &value ) )
		    Helper::halt( "invalid E line, bad numeric value" );
		  instance->set( label , value );
		}

	      else if ( t == globals::A_DBL_T )
		{
		  double value = 0;
		  
		  if ( Helper::str2dbl( tok[j] , &value ) )		    
		    instance->set( label , value );
		  else
		    if ( tok[j] != "." && tok[j] != "NA" ) 
		      Helper::halt( "invalid E line, bad numeric value" );		  
		}

	      else if ( t == globals::A_TXT_T )
		{
		  instance->set( label , tok[j] );
		}

	      else
		logger << "could not read undefined type from annotation file for " << label << "\n";

	      
// 	      // track how many annotations we add
// 	      parent_edf.aoccur[ a->name + ":" + label ]++;
	      
	    }
	  
	  ++line_count;
	}
      
    } // next line

  //  logger << "  processed " << line_count << " lines\n";
  
  FIN.close();
  
  return line_count;
}


int annot_t::load_features( const std::string & f )
{
  
  // set basic values for this annotation type, then add events/features  

  //logger << " attaching feature-list file " << f << "\n";
  
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
      std::getline( FIN , line );      
      if ( FIN.eof() || line == "" ) continue;
      
      std::vector<std::string> tok = Helper::parse( line , "\t" );
      const int n = tok.size();
      if ( n < 3 ) continue;
      
      feature_t feature;
      
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
	  
      instance_t * instance = add( feature.label , feature.feature );
      
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

  std::ofstream FOUT( t.c_str() , std::ios::out );

  FOUT << "# "
       << name << " | "
       << description ;

  std::map<std::string,globals::atype_t>::const_iterator aa = types.begin();
  while ( aa != types.end() )
    {
      if ( aa == types.begin() ) FOUT << " |";
      FOUT << " " << aa->first << "[" << globals::type_name[ aa->second ] << "]";
      ++aa;
    }
  
  FOUT << std::fixed << std::setprecision(4);
  
  //
  // Interval-based annotation
  //
  
  
  annot_map_t::const_iterator ii = interval_events.begin();
  while ( ii != interval_events.end() )
    {
      
      const instance_idx_t & instance_idx = ii->first;
      const instance_t * instance = ii->second;
      
      FOUT << name << "\t"
	   << instance_idx.id << "\t"
	   << instance_idx.interval.start/(double)globals::tp_1sec << "\t" 
	   << instance_idx.interval.stop/(double)globals::tp_1sec; 
      
      std::map<std::string,avar_t*>::const_iterator ti = instance->data.begin();
      while ( ti != instance->data.end() )
	{
	  FOUT << "\t" << ti->second->text_value();
	  ++ti;
	}

      FOUT << "\n";
      ++ii;
    }
  
  FOUT.close();
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
      Helper::warn( "did not find EpochLength in XML, defaulting to 30 seconds" );
      epoch_sec = 30;
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
	  start_tp = globals::tp_1sec * start_sec; 
	  stop_tp = start_tp + (uint64_t)( globals::tp_1sec * duration_sec ) - 1LLU ; 
	  
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
  // 0  Wake
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
	  
	  std::string stg = "Unscored";
	  if      ( e->value == "0" ) stg = "Wake";
	  else if ( e->value == "1" ) stg = "NREM1";
	  else if ( e->value == "2" ) stg = "NREM2";
	  else if ( e->value == "3" ) stg = "NREM3";
	  else if ( e->value == "4" ) stg = "NREM4";
	  else if ( e->value == "5" ) stg = "REM";	 
	 
	  interval_t interval( (uint64_t)(seconds * globals::tp_1sec ) , 
			       (uint64_t)(( seconds + epoch_sec ) * globals::tp_1sec - 1LLU ) );
 
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

  //  logger << "  reading XML annotations from " << filename << "\n";
  
  XML xml( filename );

  if ( ! xml.valid() ) Helper::halt( "invalid annotaiton file: " + filename );
  
  // Determine format: Profusion or NSRR ? 

  std::vector<element_t*> nsrr_format_test = xml.children( "PSGAnnotation" );
  bool profusion_format = nsrr_format_test.size() == 0 ;
  if ( globals::param.has( "profusion" ) ) profusion_format = true;

  const std::string EventConcept = profusion_format ? "Name"           : "EventConcept" ;
  const std::string epoch_parent = profusion_format ? "CMPStudyConfig" : "PSGAnnotation" ;
  
  
  std::vector<element_t*> scored = xml.children( "ScoredEvents" );
  

  //
  // NSRR format:
  //
  
  // assume all annotations will then be under 'ScoredEvent'
  // with children: 'EventConcept' , 'Duration' , 'Start' , and optionally 'Notes'
  
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
      
      // are we checking whether to add this file or no? 
      if ( globals::specified_annots.size() > 0 && 
	   globals::specified_annots.find( concept->value ) == globals::specified_annots.end() ) continue;
      
      // already found?
      if ( added.find( concept->value ) != added.end() ) continue;

      // otherwise, add
      annot_t * a = edf->timeline.annotations.add( concept->value );
      a->description = "XML-derived";
      a->file = filename;
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
	  if      ( e->value == "0" ) ss = "Wake";
	  else if ( e->value == "1" ) ss = "NREM1";
	  else if ( e->value == "2" ) ss = "NREM2";
	  else if ( e->value == "3" ) ss = "NREM3";
	  else if ( e->value == "4" ) ss = "NREM4";
	  else if ( e->value == "5" ) ss = "REM";	 
	  
	  // are we checking whether to add this file or no? 
	  
	  if ( globals::specified_annots.size() > 0 && 
	       globals::specified_annots.find( ss ) == globals::specified_annots.end() ) continue;
	  
	  // already found?
	  if ( added.find( ss ) != added.end() ) continue;
	  
	  // otherwise, add
	  annot_t * a = edf->timeline.annotations.add( ss );
	  a->description = "XML-derived";
	  a->file = filename;
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
      
      element_t * start    = (*e)("Start" );
      if ( start == NULL ) start = (*e)( "time" );

      element_t * duration = (*e)("Duration" );
      element_t * notes    = (*e)("Notes" );

//       if ( concept  == NULL ) std::cout << "concept NULL\n";
//       if ( start    == NULL ) std::cout << "start   NULL\n";
//       if ( duration == NULL ) std::cout << "duration NULL\n";
      
      if ( concept == NULL || start == NULL || duration == NULL ) continue;
      
      // otherwise, add 

      double start_sec, duration_sec;
      if ( ! Helper::str2dbl( start->value , &start_sec ) ) Helper::halt( "bad value in annotation" );
      if ( ! Helper::str2dbl( duration->value , &duration_sec ) ) Helper::halt( "bad value in annotation" );

      uint64_t start_tp = start_sec * globals::tp_1sec;
      uint64_t stop_tp  = duration_sec > 0 
	? start_tp + (uint64_t)( duration_sec * globals::tp_1sec ) - 1LLU 
	: start_tp;
      
      interval_t interval( start_tp , stop_tp );
      
      annot_t * a = edf->timeline.annotations.add( concept->value );
      
      instance_t * instance = a->add( concept->value , interval );      
      
      // any notes?  set as TXT, otherwise it will be listed as a FLAG
      if ( notes ) 
	instance->set( concept->value , notes->value ); 
      else
	instance->set( concept->value );
      
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
	  if      ( e->value == "0" ) ss = "Wake";
	  else if ( e->value == "1" ) ss = "NREM1";
	  else if ( e->value == "2" ) ss = "NREM2";
	  else if ( e->value == "3" ) ss = "NREM3";
	  else if ( e->value == "4" ) ss = "NREM4";
	  else if ( e->value == "5" ) ss = "REM";	 


	  // skip if we are not interested in this element
	  
	  if ( added.find( ss ) == added.end() ) continue;
	  
	  // otherwise, add

	  uint64_t start_tp = start_sec * globals::tp_1sec;
	  uint64_t stop_tp  = start_tp + (uint64_t)( epoch_sec * globals::tp_1sec ) - 1LLU ;
	  
	  // advance to the next epoch
	  start_sec += epoch_sec;
	  
	  interval_t interval( start_tp , stop_tp );	  


	  
	  annot_t * a = edf->timeline.annotations.add( ss );
	  
	  instance_t * instance = a->add( ss , interval );      
      
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
  
  // need to implement a better search... but for now just use brute force... :-(
  
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


bool annotation_set_t::make_sleep_stage( const std::string & a_wake , 
					 const std::string & a_n1 , 
					 const std::string & a_n2 , 
					 const std::string & a_n3 , 
					 const std::string & a_n4 , 
					 const std::string & a_rem ,
					 const std::string & a_other )
{

  //
  // already made?
  //
  
  if ( find( "SleepStage" ) != NULL ) return false; 
  
  std::string dwake, dn1, dn2, dn3, dn4, drem, dother;
  
  std::map<std::string,annot_t*>::const_iterator ii = annots.begin();
  while ( ii != annots.end() )
    {
      const std::string & s = ii->first;
      
      sleep_stage_t ss = globals::stage( s );
      
      if      ( ss == WAKE ) dwake = s;
      else if ( ss == NREM1 ) dn1 = s;
      else if ( ss == NREM2 ) dn2 = s;
      else if ( ss == NREM3 ) dn3 = s;
      else if ( ss == NREM4 ) dn4 = s;
      else if ( ss == REM ) drem = s;
      else if ( ss == UNKNOWN ) dother = s;
      else if ( ss == MOVEMENT ) dother = s;

      ++ii;
    }


  //
  // find annotations
  //

  annot_t * wake  = find( a_wake != ""    ? a_wake : dwake );
  annot_t * n1    = find( a_n1 != ""      ? a_n1   : dn1 );
  annot_t * n2    = find( a_n2 != ""      ? a_n2   : dn2 );
  annot_t * n3    = find( a_n3 != ""      ? a_n3   : dn3 );
  annot_t * n4    = find( a_n4  != ""     ? a_n4   : dn4 );
  annot_t * rem   = find( a_rem  != ""    ? a_rem  : drem );
  annot_t * other = find( a_other != ""   ? a_other : dother );
  
  if ( wake == NULL 
       && n1 == NULL 
       && n2 == NULL 
       && n3 == NULL 
       && rem == NULL ) 
    return false;
  
  annot_t * ss = add( "SleepStage" );
  ss->description = "SleepStage";
  
  if ( wake ) 
    {
      annot_map_t & events = wake->interval_events;
      annot_map_t::const_iterator ee = events.begin();
      while ( ee != events.end() )
	{	  
	  instance_t * instance = ss->add( globals::stage( WAKE ) , ee->first.interval );
	  ++ee;
	}
    }


  if ( n1 ) 
    {
      annot_map_t & events = n1->interval_events;
      annot_map_t::const_iterator ee = events.begin();
      while ( ee != events.end() )
	{	  
	  instance_t * instance = ss->add( globals::stage( NREM1 ) , ee->first.interval );
	  ++ee;
	}
    }

  if ( n2 ) 
    {
      annot_map_t & events = n2->interval_events;
      annot_map_t::const_iterator ee = events.begin();
      while ( ee != events.end() )
	{	  
	  instance_t * instance = ss->add( globals::stage( NREM2 ) , ee->first.interval );
	  ++ee;
	}
    }

  if ( n3 ) 
    {
      annot_map_t & events = n3->interval_events;
      annot_map_t::const_iterator ee = events.begin();
      while ( ee != events.end() )
	{	  
	  instance_t * instance = ss->add( globals::stage( NREM3 ) , ee->first.interval );
	  ++ee;
	}
    }


  if ( n4 ) 
    {
      annot_map_t & events = n4->interval_events;
      annot_map_t::const_iterator ee = events.begin();
      while ( ee != events.end() )
	{	  
	  instance_t * instance = ss->add( globals::stage( NREM4 ) , ee->first.interval );
	  ++ee;
	}
    }
  

  if ( rem ) 
    {
      annot_map_t & events = rem->interval_events;
      annot_map_t::const_iterator ee = events.begin();
      while ( ee != events.end() )
	{	  
	  instance_t * instance = ss->add( globals::stage( REM ) , ee->first.interval );
	  ++ee;
	}
    }


  if ( other ) 
    {
      annot_map_t & events = other->interval_events;
      annot_map_t::const_iterator ee = events.begin();
      while ( ee != events.end() )
	{	  
	  instance_t * instance = ss->add( globals::stage( UNKNOWN ) , ee->first.interval );
	  ++ee;
	}
    }
  
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



