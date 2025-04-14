
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
  // annotation mode: create an annotation track based on true values
  //

  const std::string annot = ! channel_mode ? param.requires( "annot" ) : "" ; 
  
  
  //
  // get the expression to evaluate
  //
  
  std::string expression = Helper::unquote( param.requires( "expr" ) , '#' );

  //
  // if evaluating a channel, ensure that the final return value is for that channel
  // (and ensuree it is sanitized, just in case this is needed)
  //
  
  if ( channel_mode ) expression += " ; " + Helper::sanitize( siglab ) ;
  
  //
  // options:
  //
  
  const bool verbose = param.has( "verbose" );

  
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

  // i.e. incase C3-M2, etc, Helper::sanitize all labels .. i.e. C3_M2
  
  std::map<std::string,std::string> clean2dirty, dirty2clean;
  for (int s=0; s<edf.header.ns; s++)
    {
      if ( edf.header.is_annotation_channel( s ) ) continue;
      dirty2clean[ edf.header.label[s] ] = Helper::sanitize( edf.header.label[s] );
      clean2dirty[ Helper::sanitize( edf.header.label[s] ) ] = edf.header.label[s] ;
    }
    
  std::map<std::string,std::vector<double> > inputs;
  
  std::set<std::string> symbols = tok.symbols();

  std::vector<uint64_t> tp;
    
  std::set<std::string>::const_iterator ss = symbols.begin();
  while ( ss != symbols.end() )
    {
      
      // assume expression should use the clean form

      if ( clean2dirty.find( *ss ) != clean2dirty.end() )
	{

	  // in case the original label is different
	  const std::string ch_label = clean2dirty[ *ss ] ; 
	  
	  int slot = edf.header.signal( ch_label );

	  // this should not happen... but just in case... skip here,
	  // the expr will return the error
	  if ( edf.header.is_annotation_channel( slot ) ) continue;

	  // check SR
	  int sr1 = edf.header.sampling_freq( slot );
	  
	  if ( sr != 0 && sr != sr1 )
	    Helper::halt( "all channels need to have similar sampling rates" );
	  else
	    sr = sr1;
	  
          slice_t slice( edf , slot , edf.timeline.wholetrace() );
	  
          const std::vector<double> * d = slice.pdata();
	  
	  if ( ( ! channel_mode ) && tp.size() == 0 )
	    tp = *slice.ptimepoints();
	  
	  if ( ch_label != *ss ) 
	    logger << "  attaching " << ch_label << " (mapped to " << *ss << ") for " << d->size() << " sample-points...\n";
	  else
	    logger << "  attaching " << ch_label << " for " << d->size() << " sample-points...\n";
	  
	  // here, original (clean) encoding expected by Eval()
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
  // else, add as an annotation
  //

  if ( ! channel_mode )
    {
      std::vector<bool> b;      
      
      if ( ! tok.value().is_bool_vector(&b) )  
	Helper::halt( "expression does not evaluate to a boolean vector" );
      
      if ( b.size() != tp.size() )
	Helper::halt( "problem aligning time-points and length of return value" );

      uint64_t start = 0;

      bool inseg = false;

      const int np = b.size();

      annot_t * a = edf.annotations->add( annot );
      
      for (int i=0; i<np; i++)
	{
	  if ( b[i] && ! inseg )
	    {
	      inseg = true;
	      start = tp[i];
	    }
	  else if ( inseg && ! b[i] )
	    {
	      // i.e. +1 encoding for ends already
	      a->add( annot , interval_t( start , tp[i] ) , "." );
	      inseg = false;
	    }
	}

      if ( inseg ) 
	{
	  // i.e. one-past end of last point
	  a->add( annot , interval_t( start , tp[ tp.size() - 1 ] + 1LLU  ) , "." );
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

