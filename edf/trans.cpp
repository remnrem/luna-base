
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
#include "helper/logger.h"
#include "defs/defs.h"

#include "helper/token-eval.h"
#include "annot/annot.h"

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>

extern logger_t logger;

extern globals global;


//
// Implement the TRANS command
//


void proc_trans( edf_t & edf , param_t & param )
{

  // operates only on channels
  // (annotations can be mapped via prior A2S command)
  
  // expects:
  //   expr=# expression #

  // 
  // can run in one of two modes::
  //
  //    sig=XX     either create or update an existing signal; expecting assignment in the expression
  //               (return value is ignored)
  //    annot=XX   create a new annotation based on the expression return value
  //      
  
  std::string siglab = param.requires( "sig" );

  //
  // create/update a channel?
  //
  
  const bool channel_mode = siglab != "*" ; 
  
  const bool update_existing_channel = channel_mode && edf.header.has_signal( siglab ); 

  //
  // get the expression to evaluate
  //
  
  std::string expression = Helper::unquote( param.requires( "expr" ) , '#' );

  //
  // if evaluating a channel, ensure that the final return value is for that channel
  //

  if ( channel_mode ) expression += " ; " + siglab ;
  
  //
  // options:
  //

  const bool verbose = param.has( "verbose" );


  //
  // set epoch flags based on this?
  //
  
  logger << "  evaluating expression  : " << expression << "\n";

  //
  // Get SR
  //

  int sr = 0;

  if ( update_existing_channel )
    sr = edf.header.sampling_freq( edf.header.signal( siglab ) ) ;
  
  //
  // output
  //
  
  instance_t out;

  //
  // expression
  //
  
  Eval tok( expression );

  //
  // inputs: get & bind any symbols (i.e. channel vevtors) required by the expression
  //
  
  std::map<std::string,std::vector<double> > inputs;
  
  std::set<std::string> symbols = tok.symbols();

  std::set<std::string>::const_iterator ss = symbols.begin();
  while ( ss != symbols.end() )
    {
      
      if ( edf.header.has_signal( *ss ) )
	{
	  
	  int slot = edf.header.signal( *ss );

	  // this should not happen... but just in case... skip here,
	  // the expr will return the error
	  if ( edf.header.is_annotation_channel( slot ) ) continue;

	  // check SR
	  int sr1 = edf.header.sampling_freq( slot );
	  
	  if ( sr != 0 && sr != sr1 )
	    Helper::halt( "all channels need to have similar sampling rates" );
	  else
	    sr = sr1;
	  
	  logger << "  attaching " << *ss << "...\n";
	  
          slice_t slice( edf , slot , edf.timeline.wholetrace() );

          const std::vector<double> * d = slice.pdata();

	  logger << "  bind " << *ss << " " << d->size() << " sample-points\n";
	  
	  inputs[ *ss ] = *d;

	}
      ++ss;
    }

  //
  // perhaps no channels were used?
  //

  if ( sr == 0 )
    Helper::halt( "no channels attached: i.e. no sample rate value attached" );
  
  //
  // bind input/output data to token evaluator
  //
  
  tok.bind( inputs , &out );
  
  //
  // evaluate expression
  //
  
  bool is_valid = tok.evaluate( verbose );

  bool retval;

  bool is_boolean_retval = true;

  if ( ! tok.value( retval ) ) is_boolean_retval = false;

  //
  // update channel  
  //

  if ( channel_mode )
    {

      std::vector<double> rdat = tok.value().as_float_vector();

      logger << "  returned " << rdat.size() << " sample-points\n";
      
      if ( update_existing_channel )
	{
	  logger << "  updating " << siglab << "...\n";
	  edf.update_signal( edf.header.signal( siglab ) , &rdat );
	}
      else
	{
	  logger << "  creating new channel " << siglab << "...\n";
          edf.add_signal( siglab , sr , rdat );
	}
      
    }
  
  //
  // Final output to log
  //

  if ( verbose )
    {
      logger << "parsed as a valid expression : " << ( is_valid ? "yes" : "no" ) << "\n";
      logger << "return value                 : " << tok.result() << "\n";
      if ( is_boolean_retval ) 
	logger << "return value (as T/F)        : " << ( retval ? "true" : "false" ) << "\n";
      logger << "assigned meta-data           : " << out.print() << "\n";  
    }
  
  // all done 
  return;
  
}

