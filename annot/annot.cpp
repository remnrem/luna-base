
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

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>

std::ostream & operator<<( std::ostream & out , const event_t & e )
{
  out << e.label;
  if ( e.has_value ) out << "=" << e.double_value();  
  return out;
}


// annotation format
// sets of intervals in MSEC 
// tab separated, one event / epoch per line

// Need to define EPOCH if E is used below
// EPOCHS are always specified as base-1

// NAME
// DESC
// EPOCH   LENGTH (sec)    OVERLAP (INCREMENT)

// col 1 :  I / E    ( 
// start \t


bool annot_t::load( const std::string & f )
{

  reset();

  // 
  // Check file exists
  //
  
  if ( ! Helper::fileExists(f) ) return false;
  
  
  //
  // XML files -- these should already have been loaded separately
  //
  
  if ( Helper::file_extension( f , "xml" ) ) 
    {
      std::cerr << "Problem :: " << f << " is an XML file... should already have been loaded...\n";
      return false;
    }
  
  
  //
  // .ftr file
  //

  if ( Helper::file_extension( f , "ftr" ) ) 
    {
      std::cerr << "Problem :: " << f << " is an FTR file... should already have been loaded...\n";
      return false;
    }
  
  
  //
  // Original .annot file format
  //
      
  std::cerr << " attaching annotation file " << f << "\n";

  //
  // record filename
  //

  file = f ;
  
  std::ifstream FIN( f.c_str() , std::ios::in );

  // if EPOCHS are specified in the command
  uint64_t epoch_length_tp;
  uint64_t epoch_overlap_tp;
  bool epoched = false;

  int line_count = 0;

  while ( ! FIN.eof() )
    {
      
      if ( FIN.bad() ) continue;
      std::string line;      
      std::getline( FIN , line );      
      if ( FIN.eof() || line == "" ) continue;

      std::vector<std::string> tok = Helper::parse( line , "\t" );
      const int n = tok.size();
      if ( n == 0 ) continue;
      
      uint64_t epoch = 0;
      interval_t interval;
      
      if ( tok[0] == "NAME" )
	{
	  //name = tok[1];+ "["+f+"]";
	  name = tok[1];
	}
      else if ( tok[0] == "TYPE" )
	{
	  type.clear();
	  for (int j=1;j<tok.size();j++) 
	    {
	      ANNOTATION t;
	      if      ( tok[j].substr(0,1) == "S" ) t = ATYPE_SLEEP_STAGE;
	      else if ( tok[j].substr(0,1) == "Q" ) t = ATYPE_QUANTITATIVE;
	      else if ( tok[j].substr(0,1) == "B" ) t = ATYPE_BINARY;
	      else if ( tok[j].substr(0,1) == "T" ) t = ATYPE_TEXTUAL;
	      else if ( tok[j].substr(0,1) == "C" ) t = ATYPE_COLCAT;
	      else if ( tok[j].substr(0,1) == "I" ) t = ATYPE_INTEGER;
	      else if ( tok[j].substr(0,1) == "M" ) t = ATYPE_MASK;
	      else t = ATYPE_UNKNOWN;
	      
	      if ( t == ATYPE_UNKNOWN )		   
		Helper::halt( "unsupported annotation type" );
	      
	      type.push_back(t);
	    }
	}
      else if ( tok[0] == "DESC" )
	{
	  description = tok[1];
	}
      else if ( tok[0] == "RANGES" )
	{
	  if ( cols.size() == 0 ) 
	    Helper::halt( "no COLS specified before RANGES" );
	  if ( tok.size() != 3 ) 
	    Helper::halt( "incorrect # of RANGES values, expecting RANGES min max" );
	  has_range = true;
	  if ( ! Helper::str2dbl( tok[1] , &min ) ) 
	    Helper::halt( "bad numeric value in RANGES" );
	  if ( ! Helper::str2dbl( tok[2] , &max ) ) 
	    Helper::halt( "bad numeric value in RANGES" );	      
	}
      else if ( tok[0] == "COLS" )
	{	  
	  for (int j=1;j<tok.size();j++) cols.push_back( tok[j] );
	}
      else if ( tok[0] == "EPOCH" )
	{
	  // EPOCH LENGTH OVERLAP
	  if ( n != 3 ) Helper::halt("invalid EPOCH line");
	  epoched = true;
	  double a,b;
	  Helper::str2dbl( tok[1] , &a );
	  Helper::str2dbl( tok[2] , &b );
	  if ( a <= 0 || b <= 0 ) Helper::halt( "bad EPOCH specs" );
	  epoch_length_tp = a * globals::tp_1sec;
	  epoch_overlap_tp = b * globals::tp_1sec;
	}

      
      //
      // Data row ('E' or 'I')
      //
      
      if ( tok[0] == "E" || tok[0] == "I" ) 
	{
	  
	  bool eline = tok[0] == "E";
	  
	  if ( cols.size() == 0 ) Helper::halt( "no COLS specified" );

	  interval_t interval;
	  
	  if ( eline )
	    {
	      
	      if ( ! epoched ) Helper::halt( "no EPOCH defined" );

	      // E  epoch-number  #cols
	      if ( n != cols.size() + 2 ) Helper::halt( "invalid E line" );
	      
	      if ( ! Helper::str2int64( tok[1] , &epoch ) ) 
		Helper::halt( "invalid E line, based epoch number" );
	      
	      if ( epoch == 0 ) 
		Helper::halt( "invalid E value of '0' (first epoch should be '1')" );
	      
	      // set interval from current line
	      interval.start = epoch_overlap_tp * (epoch-1);
	      interval.stop  = interval.start + epoch_length_tp - 1;
	      
	    }
	  else // an INTERVAL
	    {
	      
	      // I  start stop   #cols
	      if ( n != cols.size() + 3 ) Helper::halt( "invalid I line" );
	      
	      // Read as seconds -- convert internally
	      double dbl_start = 0 , dbl_stop = 0;

	      if ( ! Helper::str2dbl( tok[1] , &dbl_start ) )
		Helper::halt( "invalid interval: " + line );
	      
	      if ( ! Helper::str2dbl( tok[2] , &dbl_stop ) ) 
		Helper::halt( "invalid interval: " + line );
	      
	      interval.start = globals::tp_1sec * dbl_start;
	      interval.stop  = globals::tp_1sec * dbl_stop;	      
	      
	      //	      std::cout << "reading " << interval.start << " - " << interval.stop << "\n";
	      
	      if ( interval.start > interval.stop )
		Helper::halt( "invalid interval: " + line );
	      
	    }
	  
	  
	  //
	  // Now process values 
	  //

	  for (int j = (eline ? 2 : 3 ) ;j<n; j++)
	    {
	      
	      const int idx = j - (eline ? 2 : 3 );
	      
	      if ( idx >= type.size() || idx >= cols.size() )
		Helper::halt( "problem with annotation file" );
	      
	      ANNOTATION t = type[idx];

	      const std::string & label = cols[idx];
	      
	      event_t * e = NULL;
	      
	      if ( t == ATYPE_BINARY )
		{
		  bool value = ! ( tok[j] == "0" || tok[j] == "F" || tok[j] == "." );
		  e = new bool_event_t( label , value );
		}

	      else if ( t == ATYPE_INTEGER )
		{
		  int value = 0;
		  if ( ! Helper::str2int( tok[j] , &value ) )
		    Helper::halt( "invalid E line, bad numeric value" );
		  e = new int_event_t( label ,  value );
		}

	      else if ( t == ATYPE_SLEEP_STAGE )
		{
		  int value = globals::stage( tok[j] );
// 		  if ( ! Helper::str2int( tok[j] , &value ) )
// 		    Helper::halt( "invalid E line, bad numeric value" );
		  std::cout << "adding " << value << " for " << tok[j] << "\n";
		  e = new int_event_t( label , value );
		}
	      
	      else if ( t == ATYPE_QUANTITATIVE )
		{
		  double value = 0;

		  if ( Helper::str2dbl( tok[j] , &value ) )		    
		    e = new double_event_t( label , value );
		  else
		    if ( tok[j] != "." && tok[j] != "NA" ) 
		      Helper::halt( "invalid E line, bad numeric value" );
		}
	      else if ( t == ATYPE_MASK )
		{
		  bool value = tok[j] == "1";
		  e = new bool_event_t ( label , value );
		}
	      else if ( t == ATYPE_TEXTUAL || t == ATYPE_COLCAT )
		{
		  e = new text_event_t ( label , tok[j] );
		}
	      
	      else if ( t == ATYPE_BINARY ) 
		{
		  bool value = tok[j] == "1";
		  e = new bool_event_t ( label , value );
		}
	      else
		std::cerr << "could not read undefined type from annotation file\n";

	      if ( e ) 
		{
		  add( interval , e );
		  delete e; // remove temporary event
		}

	    }

	  ++line_count;
	}
      
    } // next line

  std::cerr << " processed " << line_count << " lines\n";

  FIN.close();
  
  if ( type.size() != cols.size() ) 
    Helper::halt( "incorrectly specified annot file: COLS not equal TYPE" );

  return true;
}


int annot_t::load_features( const std::string & f )
{

  // set basic values for this annotation type, then add events/features  

  std::cerr << " attaching feature-list file " << f << "\n";
  
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
      // Create annotation, with 'features' as metadata
      //
      
      event_t * e = new bool_event_t( feature.label , true );

      //
      // copy metadata
      //
      
      e->metadata = feature.data;
      
      //
      // Add this interval
      //
	  
      add( feature.feature , e );
      
      delete e; // remove temporary event
    
      ++line_count;
    } // next line
  
  
  std::cerr << " processed " << line_count << " lines\n";
  
  FIN.close();
  
  return line_count;

}


bool annot_t::save( const std::string & t)
{

  std::ofstream FOUT( t.c_str() , std::ios::out );
  
  FOUT << "NAME" << "\t" << name << "\n"
       << "DESC" << "\t" << description << "\n";
  
  FOUT << "COLS";
  for (int i=0;i<cols.size();i++) FOUT << "\t" << cols[i];
  FOUT << "\n";
  
  FOUT << "TYPE";
  for (int i=0;i<type.size();i++) FOUT << "\t" << type_name[type[i]];
  FOUT << "\n";
  
  if ( has_range )
    FOUT << "RANGES" << "\t" << min << "\t" << max << "\n";
  

  FOUT << std::fixed << std::setprecision(4);
  
  //
  // Interval-based annotation
  //
  
  

  interval_evt_map_t::const_iterator ii = interval_events.begin();
  while ( ii != interval_events.end() )
    {

      const interval_t & interval = ii->first;
      
      FOUT << "I" << "\t"
	   << interval.start/(double)globals::tp_1sec << "\t" 
	   << interval.stop/(double)globals::tp_1sec;
      
      evt_table_t::const_iterator ti = ii->second.begin();
      while ( ti != ii->second.end() )
	{
	  const event_t * evt = *ti;
	  FOUT << "\t" << evt->text_value();
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
      std::cerr << "**warning, did not find EpochLength in XML, defaulting to 30 seconds**\n";
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

  std::cerr << "  reading annotations from " << filename << "\n";
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
  
  // Create a separate annot_t for each type of EventConcept
  std::map<std::string,std::vector<event_t*> > evts;
  std::map<std::string,std::vector<interval_t> > ints;
  
  for (int i=0;i<scored.size();i++)
    {
      element_t * e = scored[i];
      
      if ( ! Helper::iequals( e->name , "ScoredEvent" ) ) continue;
      
      element_t * concept  = (*e)( EventConcept );
      if ( concept == NULL ) concept = (*e)( "name" );

      element_t * start    = (*e)("Start" );
      if ( start == NULL ) start = (*e)( "time" );

      element_t * duration = (*e)("Duration" );
      element_t * notes    = (*e)("Notes" );

      if ( concept  == NULL ) std::cout << "concept NULL\n";
      if ( start    == NULL ) std::cout << "start   NULL\n";
      if ( duration == NULL ) std::cout << "duration NULL\n";

      if ( concept == NULL || start == NULL || duration == NULL ) continue;
      
      // skip this..
      if ( concept->value == "Recording Start Time" ) continue;

      event_t * evt = new text_event_t( concept->value , notes ? notes->value : "1" ) ;
      double start_sec, duration_sec;
      if ( ! Helper::str2dbl( start->value , &start_sec ) ) Helper::halt( "bad value in annotation" );
      if ( ! Helper::str2dbl( duration->value , &duration_sec ) ) Helper::halt( "bad value in annotation" );

      uint64_t start_tp = start_sec * globals::tp_1sec;
      uint64_t stop_tp  = duration_sec > 0 
	? start_tp + (uint64_t)( duration_sec * globals::tp_1sec ) - 1LLU 
	: start_tp;
      
      interval_t interval( start_tp , stop_tp );
      
      evts[ concept->value ].push_back( evt );
      ints[ concept->value ].push_back( interval );
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
	  
	  event_t * evt = new text_event_t( ss , "1" );

	  uint64_t start_tp = start_sec * globals::tp_1sec;
	  uint64_t stop_tp  = start_tp + (uint64_t)( epoch_sec * globals::tp_1sec ) - 1LLU ;
	  
	  // advance to the next epoch
	  start_sec += epoch_sec;
	  
	  interval_t interval( start_tp , stop_tp );	  
	  evts[ ss ].push_back( evt );
	  ints[ ss ].push_back( interval );
	}
      
      
    }

  //
  // Add to annotation
  //

  std::map<std::string,std::vector<interval_t> >::iterator ii = ints.begin();
  while ( ii != ints.end() ) 
    {
      
      annot_t * a = edf->timeline.annotations.add( ii->first );
      
      // always qualifies as single TEXTUAL event
      a->description = "XML-derived";
      a->type.resize(1);
      a->type[0] = ATYPE_TEXTUAL;
      a->cols.resize(1);
      a->cols[0] = ".";
      
      const std::vector<interval_t> & aints = ii->second;
      const std::vector<event_t*> & aevts  = evts[ ii->first ];
      
      //      std::cerr << "Scored Event Concept [" << ii->first << "], with " << aints.size() << " entries\n";
      
      // map annotations to this file
      edf->alist[ ii->first ] = filename;
      edf->aoccur[ ii->first ] += aints.size();
      
      for (int j=0;j<aints.size();j++)
	{	  
	  a->add( aints[j] , aevts[j] );
	  delete aevts[j];
	}

      ++ii;
    }


  //
  // Misc. text code: have signal descriptions
  //

  if ( true )
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
	      std::cerr << "  changing " << label->value << " to canonical label " << canonical_label->value << "\n";
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




interval_evt_map_t annot_t::extract( const interval_t & window ) 
{

  //
  // Fetch all annotations that overlap this window
  //
  
  interval_evt_map_t r; 
  
  // need to implement a better search... but for now just use brute force... :-(
  
  interval_evt_map_t::const_iterator ii = interval_events.begin();
  while ( ii != interval_events.end() )
    {
      const interval_t & a = ii->first;
      if ( a.overlaps( window ) ) r[ a ] = ii->second;
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

  // already made?
  if ( find( "SleepStage" ) != NULL ) return false; 

  std::string dwake, dn1, dn2, dn3, dn4, drem, dother;
  
  std::map<std::string,annot_t*>::const_iterator ii = annots.begin();
  while ( ii != annots.end() )
    {
      const std::string & s = ii->first;
      
      //      std::cout << "s = " << s << "\n";

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
  ss->type.push_back( ATYPE_SLEEP_STAGE );
  ss->cols.push_back( "SleepStage" );
  ss->description = "SleepStage";
  
  if ( wake ) 
    {
      interval_evt_map_t & events = wake->interval_events;
      interval_evt_map_t::const_iterator ee = events.begin();
      while ( ee != events.end() )
	{	  
	  event_t * e = new text_event_t( "W" , globals::stage( WAKE ) );
	  ss->add( ee->first , e );
	  ++ee;
	}
    }

  if ( n1 ) 
    {
      interval_evt_map_t & events = n1->interval_events;
      interval_evt_map_t::const_iterator ee = events.begin();
      while ( ee != events.end() )
	{	  
	  event_t * e = new text_event_t( "N1" , globals::stage( NREM1 ) );
	  ss->add( ee->first , e );
	  ++ee;
	}
    }

  if ( n2 ) 
    {
      interval_evt_map_t & events = n2->interval_events;
      interval_evt_map_t::const_iterator ee = events.begin();
      while ( ee != events.end() )
	{	  
	  event_t * e = new text_event_t( "N2" , globals::stage( NREM2 ) );
	  ss->add( ee->first , e );
	  ++ee;
	}
    }
  
  if ( n3 ) 
    {
      interval_evt_map_t & events = n3->interval_events;
      interval_evt_map_t::const_iterator ee = events.begin();
      while ( ee != events.end() )
	{	  
	  event_t * e = new text_event_t( "N3" , globals::stage( NREM3 ) );
	  ss->add( ee->first , e );
	  ++ee;
	}
    }
  
  if ( n4 ) 
    {
      interval_evt_map_t & events = n4->interval_events;
      interval_evt_map_t::const_iterator ee = events.begin();
      while ( ee != events.end() )
	{	  
	  event_t * e = new text_event_t( "N4" , globals::stage( NREM4 ) );
	  ss->add( ee->first , e );
	  ++ee;
	}
    }

  if ( rem ) 
    {
      interval_evt_map_t & events = rem->interval_events;
      interval_evt_map_t::const_iterator ee = events.begin();
      while ( ee != events.end() )
	{	  
	  event_t * e = new text_event_t( "R" , globals::stage( REM ) );
	  ss->add( ee->first , e );
	  ++ee;
	}
    }

  if ( other ) 
    {
      interval_evt_map_t & events = other->interval_events;
      interval_evt_map_t::const_iterator ee = events.begin();
      while ( ee != events.end() )
	{	  
	  event_t * e = new text_event_t( "?" , globals::stage( UNKNOWN ) );
	  ss->add( ee->first , e );
	  ++ee;
	}
    }

  return true;
  
}


uint64_t annot_t::minimum_tp() const
{
  if ( interval_events.size() == 0 ) return 0;
  return interval_events.begin()->first.start;
}


uint64_t annot_t::maximum_tp() const
{
  if ( interval_events.size() == 0 ) return 0;
  return (--interval_events.end())->first.stop;
}



