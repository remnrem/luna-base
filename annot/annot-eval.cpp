
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


//
// Implement the EVAL command
//


void proc_eval( edf_t & edf , param_t & param )
{
  
  // expects:
  //   annot=name
  //   expr=# expression #
  //   globals=J,K,L

  // create a new annotation 'name' 
  // EITHER epoch-by-epoch (i.e. evaluate the expression once per epoch) 
  // OR for each unqiue set of annotations, i.e. whereby the new annotation may span several sets)
  
  // Epoch mode: annot (and meta-data) is either present or not
  // E1   Y 
  // E2   N
  // E3   N
  // E4   Y

  // Interval mode:  (if expression is A3 = A1 && A2)
  // A1
  // A1  A2  --> A3
  //     A2
  // A1  A2  --> A3
  // A1  A2      A3
  // A1
  // 
  // i.e. here two new A3 annotations, that start stop at the respective junctions
  // n.b. currently, this will consider intervals based on ALL existing annotations
  //      this may not be desirable, e.g. may cause redundant evaluations, but should
  //      be okay for now
  
  // epoch versus interval mode?
  bool interval_mode = param.has( "interval" );
  
  std::string new_annot_class = interval_mode ? param.value( "interval" ) : param.requires( "annot" );
  
  std::string expression = Helper::unquote( param.requires( "expr" ) , '#' );

  
  std::set<std::string> acc_vars;
  
  bool use_globals = param.has( "globals" );
  if ( use_globals ) 
    acc_vars = param.strset( "globals" );
  
  logger << "  evaluating expression           : " << expression << "\n";
  logger << "  derived values annotation class : " << new_annot_class ;
  if ( use_globals ) logger << " (and " << new_annot_class << "_global)";
  logger << "\n";
  
  //
  // Get all existing annotations
  //
  
  std::vector<std::string> names = edf.timeline.annotations.names();


  //
  // Create/attach new annotation class, which will have multiple
  // epoch-level instances 
  //
  
  annot_t * new_annot = edf.timeline.annotations.add( new_annot_class );
  

  // 
  // Make global annotation an entirely separate class of annotation
  //
  
  annot_t * global_annot = use_globals ? edf.timeline.annotations.add( new_annot_class + "_global" ) : NULL ;
 
  instance_t dummy;
  
  instance_t * accumulator = use_globals ? global_annot->add( "." , edf.timeline.wholetrace() , "." ) : &dummy ;
  
  //
  // We need to initialize any global variables that will appear in the main expression
  // Assume these are all floats for now, and will have the form _var 
  //
  
  if ( use_globals ) 
    {
      std::set<std::string>::const_iterator ii = acc_vars.begin();
      while ( ii != acc_vars.end() ) 
	{
	  accumulator->set( *ii , 0 );
	  ++ii;
	}
    }


  
  //
  // Iterate over epochs? Or Intervals?
  //

  int acc_total = 0 , acc_retval = 0 , acc_valid = 0; 
  int added_intervals = 0; // post concatenation

  if ( ! interval_mode ) 
    {
  
      edf.timeline.first_epoch();
      
      while ( 1 ) 
	{
	  
	  // consider _ALL_ epochs
	  
	  int e = edf.timeline.next_epoch_ignoring_mask() ;
	  
	  if ( e == -1 ) break;
	  
	  interval_t interval = edf.timeline.epoch( e );
	  
	  std::map<std::string,annot_map_t> inputs;
	  
	  // get each annotations
	  for (int a=0;a<names.size();a++)
	    {
	      
	      annot_t * annot = edf.timeline.annotations.find( names[a] );
	      
	      // get overlapping annotations for this epoch
	      annot_map_t events = annot->extract( interval );
	      
	      // store
	      inputs[ names[a] ] = events;
	    }
	  
	  
	  //
	  // create new annotation
	  //
	  
	  instance_t * new_instance = new_annot->add( "e:" + Helper::int2str( edf.timeline.display_epoch(e) ) , interval , "." );
	  
	  //
	  // evaluate the expression
	  //
	  
	  Eval tok( expression );
	  
	  tok.bind( inputs , new_instance , accumulator , &acc_vars );
	  
	  bool is_valid = tok.evaluate();
	  
	  bool retval;
	  
	  if ( ! tok.value( retval ) ) is_valid = false;
	  
	  //
	  // Output
	  //
	  
	  acc_total++;
	  
	  acc_valid += is_valid;
	  
	  if ( is_valid ) 
	    acc_retval += retval;
	  
	  // remove instance if expression was F or invalid
	  if ( ( ! is_valid ) || ( ! retval ) ) 
	    {
	      new_annot->remove( "e:" + Helper::int2str( edf.timeline.display_epoch(e) ) , interval , "." );
	    }
	  
	  // next epoch
	} 
      
    }


  //
  // Interval-level evaluation
  //

  if ( interval_mode )
    {

      //
      // step 1: get all sets of annotation change-points
      //
      
      // sort by time, collapse across events

      std::set<uint64_t> changepoints;
            
      for (int a = 0 ; a < names.size() ; a++ ) 
	{	  
	  annot_t * annot = edf.timeline.annotations.find( names[a] );	  
	  if ( annot == NULL ) Helper::halt( "internal problem in eval, cannot map annot " + names[a] );	  
	  annot_map_t::const_iterator ii = annot->interval_events.begin();
	  while ( ii != annot->interval_events.end() )
	    {	  
	      changepoints.insert( ii->first.interval.start );
	      changepoints.insert( ii->first.interval.stop ); // 1-past end, == start of new segment
	      //std::cout << " adding start, stop = " << ii->first.interval.start << "  " << ii->first.interval.stop << "\n";
	      ++ii;
	    }
	}
      
      //
      // start and end of recording also
      //
      
      changepoints.insert( 0LLU );
      changepoints.insert( edf.timeline.last_time_point_tp );

      //
      // Get all unique intervals
      //
      
      std::set<interval_t> uniq;
      uint64_t prior_start = 0LLU;
      std::set<uint64_t>::const_iterator pp = changepoints.begin();
      ++pp; // skip start
      while ( pp != changepoints.end() )
	{
	  uniq.insert( interval_t( prior_start , *pp ) );
	  prior_start = *pp;
	  ++pp;
	}
      
      //
      // iterate over each unique interval (i.e. will be spanned by the same annotations)
      // note: how we are doing this will not automatically merge contiguous regions of 
      // the new annotation; that should be fine, & if needed, we can add a 'MERGE' function
      // to apply to any generic set of annotations (i.e. if the same annot class/instance 
      // are found to be contiguous)
      //
      
      
      std::set<interval_t> new_annots;

      std::set<interval_t>::const_iterator uu = uniq.begin();
      while ( uu != uniq.end() )
	{
	  //	  std::cout << " uniq = " << uu->start << " " << uu->stop << "\n";

	  interval_t interval = *uu;
	  
	  std::map<std::string,annot_map_t> inputs;
	  
	  // get each annotations for this interval

	  for (int a=0;a<names.size();a++)
	    {
	      
	      annot_t * annot = edf.timeline.annotations.find( names[a] );
	      
	      // get overlapping annotations for this interval
	      annot_map_t events = annot->extract( interval );
	      
	      // store
	      inputs[ names[a] ] = events;
	    }

	  //
	  // create new annotation (here we don't allow assignments to the new VAR.. so not needed)
	  //
	  
	  // instance_t * new_instance = new_annot->add( "." , interval , "." );
	  instance_t dummy_instance;

	  //
	  // evaluate the expression
	  //
	  
	  Eval tok( expression );
	  
	  //tok.bind( inputs , new_instance , accumulator , &acc_vars );
	  tok.bind( inputs , &dummy_instance , accumulator , &acc_vars );
	  
	  bool is_valid = tok.evaluate();
	  
	  bool retval;
	  
	  if ( ! tok.value( retval ) ) is_valid = false;
	  
	  //
	  // Output
	  //
	  
	  acc_total++;
	  
	  acc_valid += is_valid;
	  
	  if ( is_valid ) 
	    acc_retval += retval;
	  
	  // add if valid and T
	  if ( is_valid && retval )
	    {
	      //	      std::cout << " adding NEW ANNOT " << interval.start << " " << interval.stop << "\n";
	      new_annots.insert( interval );
	    }
	  
	  // next unique interval
	  
	  ++uu;
	}

      //
      // now concatenate and add new annotations
      //
      
      if ( new_annots.size() != 0 ) 
	{
	  std::set<interval_t>::const_iterator nn = new_annots.begin();
	  interval_t curr = *nn;
	  ++nn;
	  while ( nn != new_annots.end() )
	    {
	       // std::cout << "\n\n";
	       // std::cout << "CURR = " << curr.start << " - " << curr.stop  << "\n";
	       // std::cout << "THIS = " << nn->start << " - " << nn->stop  << "\n";

	      // if this is contiguous w/ the prior one, extend
	      // current interval
	      if ( curr.stop == nn->start ) 
		{
		  //std::cout << " extending " << curr.start << " - " << curr.stop << "    to    " << nn->stop << "\n";
		  curr.stop = nn->stop;
		}
	      else // add current interval
		{
		  //std::cout << "adding " << curr.start << " - " << curr.stop << "\n";
		  new_annot->add( "." , curr , "." );
		  ++added_intervals;
		  curr = *nn;
		}
	      ++nn;
	    }

	  // and add final one
	  //std::cout << "adding final " << curr.start << " - " << curr.stop <<"\n";	  
	  new_annot->add( "." , curr , "." );	  
	  ++added_intervals;
	}
    }
  

  //
  // show accumulator output in log
  //

  logger << "  evaluated expressions/epochs  " 
	 << acc_total << " ("
	 << acc_valid << " valid, " 
	 << acc_retval << " true)\n";
  if ( interval_mode ) 
    logger << "  added " << added_intervals << " distinct " << new_annot_class << " interval-annotations\n";
  logger << "  global variables (if any):\n" << accumulator->print( "\n" , "\t" ) ;
  logger << "\n";

  // all done 
  return;
  
}

