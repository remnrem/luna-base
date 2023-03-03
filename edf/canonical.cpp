
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

#include "canonical.h"

#include "edf.h"
#include "defs/defs.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "db/db.h"
#include "dsp/resample.h"
#include "eval.h"

#include <iostream>
#include <fstream>

extern writer_t writer;
extern logger_t logger;

std::vector<canon_rule_t> canonical_t::rules;
std::map<std::string,std::string> canonical_t::aliases;
std::map<int,std::string> canonical_t::scale_codes;

// special case: single label means 'close out'
canon_rule_t::canon_rule_t( const std::string & label  )
{
  canonical_label = label;
  closed = true;
}

canon_rule_t::canon_rule_t( const std::vector<std::string> & lines )
{
  const int l = lines.size();
  
  // logger << "rule:\n";
  // for (int i=0; i<l; i++) logger << i << "\t[" << lines[i] << "]\n";
  
  // default req values (i.e. do not require)
  req_sr_min = req_sr_max = 0;
  // 0 none, -1 all NEG, +1 all POS, +2 = NEG & POS
  req_scale = 0;

  // default set values (i.e. do not change)
  set_sr = 0; // do not change
  set_unit = ".";

  // <<- rules
  relabel_canonical = false;

  // special closed rule [ignore_generics]
  closed = false;
    
  std::string current_rule = "";

  bool rules_set = false;
  
  for (int i=0; i<l; i++)
    {
      
      const std::string & line = lines[i];
      
      if ( Helper::trim(line) == "" ) continue;
      if ( line[0] == '%' ) continue;
      
      // no indentation: new canonical label
      //  -- which might also use the canon <- alias1 alias2 ... syntax
      // will only have one of these
      if ( line[0] != ' ' )
	{
	  
	  // special case 1: (handled below)
	  //  canon 

	  // equals

	  //  canon
	  //   req:
	  //    sig = canon

	  // special case 2: (handled here)
	  //  canon <- alias1 alias2 

	  // equals
	  
	  //  canon
	  //   req:
	  //    sig = canon,alias1,alias2


	  // special case 3: (handled here)
	  //  canon <<- canon1 canon2 canon3 

	  // this relabels 'canon1', 'canon2' etc all to 'canon'
	  //  i.e. to handle case where same EDF has > 1 version
	  //  of a canonical signal;   hopefully, if handled properly,
	  //  the transducer field, or other, will disambiguate downstream
	  //  When loading back into Luna, they will automatically be uniquified
	  
	  std::vector<std::string> t1 = Helper::quoted_parse( line , ", \t" );

	  if ( t1.size() > 2 && t1[1] == "<-" )
	    {	      
	      canonical_label = t1[0];
	      std::vector<std::string> tfin;
	      for (int j=2; j<t1.size(); j++) // skips canon <-
		{
		  std::vector<std::string> t2 =
		    Helper::quoted_parse( canonical_t::swap_in_alias( t1[j] ) , ", \t" );
		  for (int k=0; k < t2.size(); k++)
		    tfin.push_back( t2[k] );
		}

	      // add as req_sigs...
	      for (int j=0; j<tfin.size(); j++)
		{
		  req_sig.push_back( Helper::unquote( Helper::toupper( tfin[j] ) ) );		  
		}
	      
	      // also add canonical label itself
	      req_sig.push_back( Helper::unquote( Helper::toupper( canonical_label ) ) );
	      
	    }

	  // otherwise, handled "<<-" special rule
	  else if ( t1.size() > 2 && t1[1] == "<<-" )
            {	      
	      canonical_label = t1[0];
	      relabel_canonical = true;
	      original_canonical_label.clear();
	      for (int j=2; j<t1.size(); j++)
		original_canonical_label.push_back( t1[j] );
	    }	  

	  // otherwise, not a special case.. this is the start of a rule
	  else 
	    {	      
	      canonical_label = Helper::sanitize( Helper::trim( line ) );
	    }
	  
	  // add self name to unless list
	  // (i.e. only do each rule once)
	  unless.insert( canonical_label );
	  current_rule = "";
	}
      else
	{

	  if ( canonical_label == "" )
	    Helper::halt( "canonical signal not identified yet:\n" + line );

	  if ( line.size() <= 2 )
	    Helper::halt( "invalid line:\n" + line );

	  rules_set = true;
	  
	  const bool rtype = line[1] != ' ';
	  
	  // the rule type? 
	  if ( rtype )
	    {
	      std::string r = Helper::trim( Helper::toupper( line ) , ' ' , '\t' ) ;
	      if ( r == "GROUP:" )
		current_rule = "group";
	      else if ( r == "REQ:" || r == "REQUIRES:" )
		current_rule = "req";
	      else if ( r == "UNLESS:" )
		current_rule = "unless";
	      else if ( r == "SET:" || r == "SETS:" )
		current_rule = "set";
	      else
		Helper::halt( "unrecogized type of rule:\n" + r );
	    }
	  else
	    {
	      if ( current_rule == "" )
		Helper::halt( "no current rule type (group:, req:, unless: or set:) specified:\n" + line );

	      // values here -- allow aliases to be swapped in for all cases

	      if ( current_rule == "group" )
		{
		  std::vector<std::string> tokb = Helper::quoted_parse( line , "," );
		  for (int j=0; j<tokb.size(); j++)
		    {
		      std::vector<std::string> tokc = Helper::quoted_parse( canonical_t::swap_in_alias( Helper::trim( tokb[j] ) ) , "," );
		      for (int k=0; k<tokc.size(); k++)
			group.insert( Helper::trim( tokc[k] )) ;
		    }
		}
	      else if ( current_rule == "unless" )
		{
		  std::vector<std::string> tokb = Helper::quoted_parse( line , "," );
		  for (int j=0; j<tokb.size(); j++)
		    {
		      std::vector<std::string> tokc = Helper::quoted_parse( canonical_t::swap_in_alias( Helper::trim( tokb[j] ) ) , "," );
		      for (int k=0; k<tokc.size(); k++)
                        unless.insert( Helper::trim( tokc[k] ) );		      
		    }
		  }
	      else if ( current_rule == "req" )
		{
		  std::vector<std::string> tok = Helper::quoted_parse( line , "=" );
		  if ( tok.size() != 2 )
		    Helper::halt( "expecting key = value format:\n" + line );
		  
		  std::string key = Helper::trim( Helper::toupper( tok[0] ) , ' ', '\t' );
		  std::string value = Helper::trim( tok[1] , ' ', '\t' );
		  std::vector<std::string> tokb = Helper::quoted_parse( value , "," );
		  std::vector<std::string> tok2;

		  // swap in aliases?
		  if ( key == "SIG" || key == "REF" || key == "TRANS" || key == "UNIT" )
		    for (int j=0; j<tokb.size(); j++)
		      {
			std::vector<std::string> tokc =
			  Helper::quoted_parse( canonical_t::swap_in_alias( tokb[j] ) , "," );
			for (int k=0; k<tokc.size(); k++)
			  {
			    tok2.push_back( tokc[k] );
			    //std::cout << " addding [" << tokc[k] << "]\n";
			  }
		      }
		    
		  // SIG, REF, TRANS, SR-MIN, SR-MAX, UNIT, SCALE
		  // n.b. signals/references are vectors, as order is important
		  
		  if ( key == "SIG" )
		    {
		      for (int j=0; j<tok2.size(); j++)
			req_sig.push_back( Helper::sanitize( Helper::trim( Helper::unquote( Helper::toupper( tok2[j] ) ) ) ) );
		    }
		  else if ( key == "REF" )
		    {
		      for (int j=0; j<tok2.size(); j++)
                        req_ref.push_back( Helper::sanitize( Helper::trim( Helper::unquote( Helper::toupper( tok2[j] ) ) ) ));
		    }
		  else if ( key == "TRANS" )
		    {
		      // n.b. first version cannot be an alias
		      // n.b. . means missing (empty)
		      // check that a field hasn't already been specified as something else i.e. cannot have one-to-many mapping
		      
		      std::string pref_str = Helper::sanitize( Helper::trim( Helper::unquote( Helper::toupper( tok2[0] ) ) ) , '*' );
		      
		      for (int j=0; j<tok2.size(); j++)
                        {
			  // allow transducer fields to retain wildcard '*' in sanitization
			  std::string str = Helper::sanitize( Helper::trim( Helper::unquote( Helper::toupper( tok2[j] ) ) ) , '*' );
			  //if ( str == "*" && j == 0 ) Helper::halt( "first field cannot be '*' for " + canonical_label );
			  if ( canonical_t::empty_field( str ) )
			    {
			      if ( j == 0 ) Helper::halt( "first field cannot be '.' for " + canonical_label ); 
			      str = ".";
			    }

			  if ( req_transducer.find( str ) != req_transducer.end() && req_transducer[ str ] != pref_str )
			    Helper::halt( "cannot specify a transducer type multiple times: " + str + " for " + canonical_label );
			  req_transducer[ str ] = pref_str;
			}
		    }
		  else if ( key == "UNIT" )
		    {
		      // n.b. first version cannot be an alias
		      // n.b. '.' means missing field 
		      // as above, cannot specify a value multiple times
		      // allow '*' as a pref str too 

		      std::string pref_str = Helper::sanitize( Helper::trim( Helper::unquote( Helper::toupper( tok2[0] ) ) ) , '*' );
		      
		      for (int j=0; j<tok2.size(); j++)
			{
			  // allow unit fields to retain wildcard '*' in sanitization
			  std::string str = Helper::sanitize( Helper::trim( Helper::unquote( Helper::toupper( tok2[j] ) ) ) , '*' );
			  
			  //if ( str == "*" && j == 0 ) Helper::halt( "first field cannot be '*' for " + canonical_label );
			  if ( canonical_t::empty_field( str ) )
                            {
                              if ( j == 0 ) Helper::halt( "first field cannot be '.' for " + canonical_label );
                              str = ".";
                            }

			  if ( req_unit.find( str ) != req_unit.end() && req_unit[ str ]  != pref_str )
			    Helper::halt( "cannot specify different unit types : " + str + " and " + req_unit[ str ] + " for " + canonical_label );
			  req_unit[ str ] = pref_str;
			}
		    }
		  else if ( key == "SR-MIN" || key == "MIN-SR" )
		    {
		      if ( ! Helper::str2int( value , &req_sr_min ) )
			Helper::halt( "invalid integer minimum sample rate requirement:\n" + line );
		    }
		  else if ( key == "SR-MAX" || key == "MAX-SR" )
		    {
		      if ( ! Helper::str2int( value , &req_sr_max ) )
			Helper::halt( "invalid integer maximum sample rate requirement:\n" + line );
		    }
		  else if ( key == "SCALE" )
		    {
		      std::string value2 = Helper::toupper( value );
		      
		      if ( value2 == "POSNEG" || value2 == "AC" )
			req_scale = 2;
		      else if ( value2.substr(0,3) == "POS" )
			req_scale = 1;
		      else if ( value2.substr(0,3) == "NEG" )
			req_scale = -1;
		      else if ( value2 == "NONE" ) // the default
			req_scale = 0;	      
		      else
			Helper::halt( "bad scale requirement code:\n" + line );
		      
		    }
		  else
		    Helper::halt( "did not recognized required value:\n " + line );
		}		  
	      else if ( current_rule == "set" )
		{
		  
		  std::vector<std::string> tok = Helper::quoted_parse( line , "=" );
		  if ( tok.size() != 2 )
                    Helper::halt( "expecting key = value format:\n" + line );
		  
                  std::string key = Helper::trim( Helper::toupper( tok[0] ) );
                  std::string value = Helper::trim( tok[1] );
                  std::vector<std::string> tok2 = Helper::quoted_parse( value , "," );

		  //
                  // SR, UNIT, TRANS
		  //
		  
		  if ( key == "UNIT" )
                    {
		      set_unit = value;
		      if ( set_unit != "uV" && set_unit != "mV" && set_unit != "V" )
			Helper::halt( "currently can only set units to uV, mV or V : " + canonical_label );
		    }
		  else if ( key == "SR" )
		    {
		      if ( ! Helper::str2int( value , &set_sr ) )
			Helper::halt( "invalid integer sample rate value:\n" + line );
		      if ( set_sr <= 0 || set_sr > 10000 )
			Helper::halt( "invalid value for setting the sample rate for " + canonical_label );
		    }
		  else
		    Helper::halt( "did not recognized set value:\n " + line );
		  
		}	      
	    } 
	}
    
    } // next line

  //
  // check for special case, where only a canonical label is specified...
  //  i.e. add canonical label as its own signal requirement:
  //

  // also add canonical label itself 
  if ( ! rules_set ) 
    req_sig.push_back( Helper::unquote( Helper::toupper( canonical_label ) ) );

  //
  // all done
  //
}
  
    

canon_edf_signal_t::canon_edf_signal_t( edf_header_t & hdr , const int slot )
{
  if ( slot < 0 || slot >= hdr.ns )
    Helper::halt( "bad EDF header slot" );
  
  label = Helper::sanitize( Helper::trim( Helper::toupper( hdr.label[ slot ] ) ) );
  sr = hdr.sampling_freq( slot );
  unit = Helper::sanitize( Helper::trim( Helper::toupper( hdr.phys_dimension[ slot ] ) ) );
  transducer = Helper::sanitize( Helper::trim( Helper::toupper( hdr.transducer_type[ slot ] ) ) );

  if ( canonical_t::empty_field( unit ) ) unit = ".";
  if ( canonical_t::empty_field( transducer ) ) transducer = ".";
  
  scale = 0; // -1, 0, +1
  double phys_min = hdr.physical_min[slot] < hdr.physical_max[slot] ? hdr.physical_min[slot] : hdr.physical_max[slot];
  double phys_max = hdr.physical_min[slot] < hdr.physical_max[slot] ? hdr.physical_max[slot] : hdr.physical_min[slot];

  // TODO: ?? need to allow floating-point variation here?
  // (i.e. non-zero zero constant)
  const double zero_constant = 0;
  if ( phys_max < zero_constant ) scale = -1; // all negative
  else if ( phys_min >= -zero_constant ) scale = 1; // all positive  
  if ( phys_min < zero_constant && phys_max > zero_constant ) scale = 2;
};


// alternatively, in 'dry_run' mode, allow 'new' canonical signals (that are not in
// EDF) to be entered here - i.e. if they feature as for subsequent CS

canon_edf_signal_t::canon_edf_signal_t( const std::string & label , 
					const int sr ,
					const std::string & unit , 
					const std::string & transducer ,
					const int scale )
  : label(label) , sr(sr) , unit( unit ) , transducer( transducer ) , scale( scale ) 
{

}


int canonical_t::read( const std::string & filename )
{
  if ( ! Helper::fileExists( filename ) )
    Helper::halt( "could not open " + filename );
  
  std::ifstream IN1( filename.c_str() , std::ios::in );

  std::vector<std::string> lines;

  while ( 1 )
    {
      std::string line;

      Helper::safe_getline( IN1 , line );
      if ( IN1.eof() ) break;
      if ( line == "" ) continue;
      if ( line == "_quit" ) break;
      if ( line[0] == '%' ) continue;
      if ( line[0] == '#' ) continue;

      // process any aliases separately
      std::vector<std::string> tok = Helper::quoted_parse( line , "\t " );

      //
      // Variable assignment:
      //   let XX=A,B,C,D
      //
      
      if ( line.size() >= 4
	   && tok.size() >= 2
	   && Helper::toupper( line.substr(0,4) ) == "LET " 
	   && line.find( "=" ) != std::string::npos )
	{
	  // comma or space delimited for the second set
	  std::vector<std::string> str = Helper::quoted_parse( line.substr(4) , " ,=" );

	  if ( str.size() < 2 ) Helper::halt( "requires A=X or A=X,Y,Z or X = X Y Z" ); 
	  
	  // build up a single command-delimited list
	  // (i.e. do not replace if the same key used)
	  // let s = x,y,z
	  // let s = a b c
	  //  s --> is swapped in for 'x,y,z,a,b,c'

	  const std::string & key = str[0];
	  
	  for (int i=1; i<str.size(); i++)
	    {
	      if ( aliases[ key ] != "" )
		aliases[ key ] = aliases[ key ] + "," + str[i] ;
	      else
		aliases[ key ] = str[i] ;
	    }
	  
	  // this is not needed to be part of a rule now -- aliases are generic across all rules
	  continue;
	}


      //
      // one-line close rule (i.e. for group-specific rules, to stop generics from being run)
      //

      // 'closed: '
      if ( line.size() >= 8
	   && tok.size() >= 2
	   && Helper::toupper( line.substr(0,8) ) == "CLOSED: " ) 
	{
	  // comma or space delimited for the second set
	  std::vector<std::string> str = Helper::quoted_parse( line.substr(8) , " ," );
	  
	  // add special 'close-out' rule (i.e. single string constructer)
	  for (int i=0; i<str.size(); i++)
	    rules.push_back( canon_rule_t( str[i] ) );
	  
	  // nothing else to do here
	  continue;
	}

      
      //
      // process line as part of a rule
      //

      const bool is_name = line[0] != ' ';
      
      if ( is_name )
	{
	  if ( lines.size() == 0 )
	    lines.push_back( line );
	  else
	    {	      
	      rules.push_back( canon_rule_t( lines ) );
	      lines.clear();
	      lines.push_back( line );		      
	    }
	}
      else
	lines.push_back( line );
    }
  
  // flush buffer and add last rule
  if ( lines.size() != 0 )
    rules.push_back( canon_rule_t( lines ) );
  
  IN1.close();
  
  return rules.size();
}


bool canonical_t::apply_this( const std::string & rule )
{
  // if inclusions specified, it needs to be in here
  const bool inc_rule = canins.size() == 0 || canins.find( rule ) != canins.end();

  // if exclusions specified, it needs to not be in here
  const bool exc_rule = canouts.size() == 0 || canouts.find( rule ) == canouts.end();

  return inc_rule && exc_rule;
}

canonical_t::canonical_t( edf_t & edf , param_t & param )
  : edf(edf)
{

  //
  // some set up
  //

  scale_codes[ 0 ] = "NONE";
  scale_codes[ 1 ] = "POS";
  scale_codes[ -1 ] = "NEG";
  scale_codes[ 2 ] = "AC";

  //
  // spcify a subset of canonical rules to apply? 
  //

  if ( param.has( "inc" ) ) canins = param.strset( "inc" );
  if ( param.has( "exc" ) ) canouts = param.strset( "exc" );
  
  
  //
  // read in rules, if not already done, from 1 or more files
  //

  const bool has_prefix = param.has( "prefix" );
  const std::string prefix = has_prefix ? Helper::expand( param.value( "prefix" ) ) : "" ;

  if ( rules.size() == 0 )
    {
      
      if ( ! param.has( "file" ) )
	Helper::halt( "CANONICAL requires a 'file' argument" );
      
      std::vector<std::string> filenames = param.strvector( "file" ) ;
      
      for (int i=0; i<filenames.size(); i++)
	{
	  std::string filename = Helper::expand( filenames[i] );

	  // add a prefix to any relative paths, if one was specified
	  if ( filename.size() > 1 && filename[0] != globals::folder_delimiter )
	    filename = prefix + filename;

	  // read rules
	  int nrules = read( filename );
	  
	  logger << "  read " << nrules << " rules from " << filename << "\n";
	  
	}
      

      
      logger << "  in total, read " << rules.size()
	     << " rules and " << aliases.size() << " variables\n";
      
      // number that will actually be applied
      if ( canins.size() != 0 || canouts.size() != 0 )
	{
	  int napp = 0;
	  for (int r=0; r<rules.size(); r++)
	    {
	      const canon_rule_t rule = rules[r];
	      if ( apply_this( rule.canonical_label ) ) ++napp;
	    }
	  logger << "  of these, " << napp << " rules will be applied to the dataset, based on inc/exc options\n";
	}

      logger << "\n";
    }
  
  
  //
  // Other options
  //

  if ( param.has( "group" ) )
    group = param.strset( "group" );
  
  drop_originals = param.yesno( "drop-originals" );
  
  dry_run = param.yesno( "check" );

  only_check_labels = param.has( "mapper-util-mode" );

  verbose = param.has( "verbose" );

  retain_prefiltering = param.yesno( "prefiltering" );
  
  //
  // Get info on the available channels
  //

  const int ns = edf.header.ns;

  for (int s=0; s<ns; s++)
    {
      if ( edf.header.is_annotation_channel( s ) ) continue;
      
      canon_edf_signal_t sig( edf.header , s );
      
      signals.insert( sig );
      
    }

  logger << "  " << signals.size() << " signals from EDF\n";

  // 
  // do the processing
  //
  
  proc();

  return;
}


void canonical_t::add_alias( const std::string & primary , const std::string & terms )
{
  std::vector<std::string> tok = Helper::quoted_parse( terms , "," );
  for (int i=0; i<tok.size(); i++) aliases[ Helper::toupper( tok[i] ) ] = primary;
  return;
}


void canonical_t::proc( )
{
  
  //
  // track attempted & completed canonical signals
  //

  std::set<std::string> attempted, completed;
  
  //
  // track which original signals were used / should not be used
  //

  std::set<std::string> used;
  std::set<std::string> do_not_drop;

  //
  // track that, if a group has been specified, but did not match, and
  // we encountered a closed command, then we do not allow any
  // subsequent generic rules to match for that target canonical
  // signal
  //

  std::set<std::string> ignore_generics;  

  //
  // in case we are dropping originals later, get a list of those now
  //

  const bool only_data_signals = true;

  signal_list_t osignals = edf.header.signal_list( "*" , only_data_signals );

  
  //
  // sequentially track through each rule
  //
  
  for (int r=0; r<rules.size(); r++)
    {

      const canon_rule_t rule = rules[r];


      //
      // is this a special close-out rule?
      //

      if ( rule.closed )
	{
	  if ( verbose ) logger << "\n  - attempting rule " << r+1 << " of " << rules.size() << " : target = " << rule.canonical_label << "\n";
 	  if ( verbose ) logger << "   closing out all generic rules for " << rule.canonical_label << "\n";
 	  ignore_generics.insert( rule.canonical_label );
	  continue;
	}
      
      //
      // are we skipping this?
      //
      
      if ( ! apply_this( rule.canonical_label ) )
	continue;
      
      //
      // start attempting to do this rule
      //
      
      if ( verbose ) logger << "\n  - attempting rule " << r+1 << " of " << rules.size() << " : target = " << rule.canonical_label << "\n";
      
      attempted.insert( rule.canonical_label );
      
      //
      // has this rule already been satisfied?
      //
      
      const bool already_processed = is_in( rule.canonical_label , completed );
      
      if ( already_processed )
	{
	  if ( verbose ) logger << "   already processed " << rule.canonical_label << "\n";
	  continue;
	}

      //
      // handle special case: <<- canonical relabelling
      //

      if ( rule.relabel_canonical )
	{
	  int v = 0;
	  for (int j = 0; j < rule.original_canonical_label.size(); j++)
	    {
	      const std::string & olab = rule.original_canonical_label[j];
	      
	      if ( edf.header.has_signal( olab ) )
		{
		  std::string nlab = rule.canonical_label + ( j ? "." + Helper::int2str(j) : "" ) ; 
		  
		  if ( ! dry_run )
		    logger << "   renaming canonical " << olab
			   << " as " << nlab << "\n";
		  edf.header.rename_channel( olab , nlab );		  
		}
	    }

	  // now advance to the next rule
	  continue;
	}
      
      //
      // group specifier?
      //
      
      if ( rule.group.size() != 0 )
       	{

	  if ( verbose )
	    logger << "   rule group(s) [ " << print( rule.group ) << " ]\n";
	  
       	  // if a group has been specified for this rule, then
       	  // we need to match it (with at least one of the rules)

	  if ( ! is_in( group , rule.group ) )
       	    {
	      if ( verbose )
		logger << "   bailing: EDF group(s) did not match [ " << print( group ) << " ]\n";
	      continue;
	    }
	  
          // track that we have encountered a group-specific rule for
	  // a specified group: this means that any generic rules will
	  // be ignored for this canonical signal

	  //  NO ... this is now done explicitly via closed: special rule
	  // ignore_generics.insert( rule.canonical_label );

	}
      else
	{
	  // for all generic rules, ensure that we have not already encountered
	  // an unmatched, non-generic version
	  
	  if ( ignore_generics.find( rule.canonical_label ) != ignore_generics.end() )
	    {
	      if ( verbose ) logger << "   bailing: this rule has been previously closed out\n";
	      continue;
	    }
	}
      

      //
      // can we find a matching primary signal?
      //

      std::string matched_sig = "";

      // logger << " req [ " << print( rule.req_sig )  << " ]\n";
      // //      logger << " sig [ " << print( signals )  << " ]\n";
      // std::set<canon_edf_signal_t>::const_iterator ss = signals.begin();
      // while ( ss != signals.end() )
      // 	{
      // 	  logger << " sig [" << ss->label << "]\n";
      // 	  ++ss;
      // 	}

      bool matched = match( rule.req_sig , signals , &matched_sig );
      
      if ( ! matched )
	{
	  if ( verbose ) logger << "   bailing: no EDF channel matched required sig [ " << print( rule.req_sig ) << " ]\n";
	  continue;
	}
      
      if ( verbose ) logger << "   matched " << matched_sig << " from sig [ " << print( rule.req_sig ) << " ]\n";
			   
      std::set<canon_edf_signal_t>::const_iterator sig = signals.find( matched_sig );

      if ( sig == signals.end() ) Helper::halt( "internal error in canonical match sig" );
      
      //
      // a reference signal, if required
      //
      
      std::string matched_ref = "";
 
      if ( rule.req_ref.size() != 0 )
	{
	  
	  matched = ref_match( rule.req_ref , signals , &matched_ref );
	  
	  if ( ! matched )
	    {	  
	      if ( verbose ) logger << "   bailing: no EDF channel matched required ref [ " << print( rule.req_ref ) << " ]\n";
	      continue;
	    }
	  
	  if ( verbose ) logger << "   matched " << matched_ref << " from ref [ " << print( rule.req_ref ) << " ]\n";	  
	  
	}
      
      //
      // transducer rules?
      //

      bool wild_trans = false;
      bool wild_trans_pref = false;
      
      if ( rule.req_transducer.size() != 0 )
	{
	  if ( rule.req_transducer.find( sig->transducer ) == rule.req_transducer.end() )
	    {
	      
	      // wildcard?  (note: sanitization spared * -> _ conversion in trans special case)
	      if ( rule.req_transducer.find( "*" ) != rule.req_transducer.end() && sig->transducer != "." )
		{
		  
		  wild_trans_pref = rule.req_transducer.find( "*" )->second == "*";
		  
		  if ( verbose )
                    {
                      if ( wild_trans_pref )
                        logger << "   allowing wildcard '*' match for " << sig->transducer << ", will keep as is\n";
                      else
                        logger << "   allowing wildcard '*' match for " << sig->transducer << ", will set to "
                               << rule.req_transducer.find( "*" )->second << "\n";
                    }


		  if ( verbose )
		    logger << "   allowing wildcard '*' match for " << sig->transducer << ", will set to "
			   << rule.req_transducer.find( "*" )->second << "\n";

		  wild_trans = true;		  
		}
	      else
		{
		  if ( verbose ) logger << "   bailing: did not match " << sig->transducer
					<< " from transducer [ " << print( rule.req_transducer ) << " ]\n";
		  continue;
		}
	      
	    }
	  if ( verbose ) logger << "   matched " << sig->transducer << " from transducer [ " << print( rule.req_transducer ) << " ]\n";

	}


      //
      // unit rules?
      //
      
      bool wild_unit = false; // match
      bool wild_unit_pref = false; // pref is wild (i.e. keep as is) 
      
      if ( rule.req_unit.size() != 0 )
	{
          if ( rule.req_unit.find( sig->unit ) == rule.req_unit.end() )
	    {
	      // wildcard? (note: sanitization spared * -> _ conversion in unit special case)
	      // only match non-missing values
	      if ( rule.req_unit.find( "*" ) != rule.req_unit.end() && sig->unit != "." )
                {
		  wild_unit_pref = rule.req_unit.find( "*" )->second == "*";

                  if ( verbose )
                    {
		      if ( wild_unit_pref )
			logger << "   allowing wildcard '*' match for " << sig->unit << ", will keep as "
			       << sig->unit << "\n";
		      else
			logger << "   allowing wildcard '*' match for " << sig->unit << ", will set to "
			       << rule.req_unit.find( "*" )->second << "\n";
		    }
		  
                  wild_unit = true;
                }
              else
                {
		  if ( verbose ) logger << "   bailing: did not match " << sig->unit << " from unit [ " << print( rule.req_unit ) << " ]\n";
		  continue;
		}
	    }
	  if ( verbose ) logger << "   matched " << sig->unit << " from unit [ " << print( rule.req_unit ) << " ]\n";
	}

      //
      // sample rate rules?
      //

      if ( rule.req_sr_min != 0 )
	{
	  if ( sig->sr < rule.req_sr_min )
	    {
	      if ( verbose ) logger << "   bailing: did not satisfy min sr " <<  sig->sr << " < " <<  rule.req_sr_min << "\n";
	      continue;
	    }
	  if ( verbose ) logger << "   sample rate satisfies min sr " << sig->sr << " >= " <<  rule.req_sr_min << "\n";
	}

      if ( rule.req_sr_max != 0 )
	{

	  if ( sig->sr > rule.req_sr_max )
	    {
	      if ( verbose ) logger << "   bailing: did not satisfy max sr " <<  sig->sr << " > " <<  rule.req_sr_max << "\n";
	      continue;
	    }
	  if ( verbose ) logger << "   sample rate satisfies max sr " << sig->sr << " <= " <<  rule.req_sr_max << "\n";

	}

      //
      // Scale?
      //

      if ( rule.req_scale != 0 )
	{
	  if ( rule.req_scale != sig->scale )
	    {
	      if ( verbose ) logger << "   bailing: did not satisfy scale "
				    << scale_codes[ sig->scale ]
				    << " != " << scale_codes[ rule.req_scale ] << "\n";
	      continue;
	    }
	  if ( verbose ) logger << "   satisfies scale "
				<< scale_codes[ sig->scale ]
				<< " == " << scale_codes[ rule.req_scale ] << "\n";

	}


      //
      // If here, we have a match
      //

      if ( dry_run ) 
	logger << "  matched rule for " << rule.canonical_label << "\n";

      //
      // Construct the CS
      //

      if ( ! dry_run )
	{
	  logger << "  + generating canonical signal " << rule.canonical_label
		 << " from existing signal(s) " << matched_sig;
	  if ( matched_ref != "" ) logger << " / " << matched_ref ;
	  logger << "\n";
	}
      
      //
      // Does the canonical already exist? 
      //

      const bool already_present = edf.header.has_signal( rule.canonical_label );
      
      
      //
      // track that we are using these channels
      //
      
      used.insert( Helper::toupper( matched_sig ) );

      if ( matched_ref != "" )
	{
	  // in case we have linked references
	  //std::cout << " matched ref = " << matched_ref << "\n";
	  std::vector<std::string> tok = Helper::parse( matched_ref , "," );
	  for (int i=0; i<tok.size(); i++)
	    used.insert( Helper::toupper( tok[i] ) );
	}
      
      //
      // Track that the original should not be dropped, as it features
      //
      
      if ( drop_originals && already_present )
	do_not_drop.insert( Helper::toupper( rule.canonical_label ) );


      //
      // Copy signal(s)
      //
      
      signal_list_t siglst = edf.header.signal_list( matched_sig );
      
      signal_list_t reflst;
      if ( matched_ref != "" ) 
	reflst = edf.header.signal_list( matched_ref );
      
      //
      // Generate the new signal (w/ re-referencing optionally)
      //   false , false = not dereference , not verbose
      // 
      
      if ( ! dry_run )
	{

	  //
	  // create a new channel?
	  //
	  
	  if ( ! already_present )
	    {

	      if ( verbose )
		logger << "   creating a new EDF signal "
		       << rule.canonical_label << "\n";
	      
	      edf.reference( siglst ,    // original channel
			     reflst ,    // reference(s) [ or none ] 
			     true ,      // create a new channel? 
			     rule.canonical_label,  // new channel name
			     rule.set_sr , // new channel SR (0 = do not chang)
			     false ,  // do not de-reference		       
			     false ); // not in verbose mode

	      
	    }
	  
	}
	  

      //
      // Get the canonical signal (either newly created, or already existing)
      //
      
      signal_list_t canonical_signal = edf.header.signal_list( rule.canonical_label );

      const int canonical_slot = dry_run ? -1 : canonical_signal(0);
      
      
      //
      // re-reference existing channel?
      //
      
      if ( matched_ref != "" && already_present && ! dry_run )
	{
	  // i.e. do not make a new channel, but if we need to re-reference existing one
	  if ( verbose ) logger << "   re-referencing "
				<< rule.canonical_label
				<< " against " << matched_ref << "\n";
	  
	  edf.reference( siglst , // original channel
			 reflst , // reference(s)
			 false ,  // do not create a new channel
			 "" , // new channel name (ignored)
			 0 , // new channel SR (ignored)
			 false , // do not de-ref
			 false ); // not in verbose mode

	}
      
      //
      // resample an existing channel? (i.e. if not already done above)
      //
      
      if ( rule.set_sr != 0 && already_present && ! dry_run )
	{
	  if ( verbose ) logger << "   re-sampling "
                                << rule.canonical_label
                                << " to SR = " << rule.set_sr << " Hz\n";

	  dsptools::resample_channel( edf , canonical_slot , rule.set_sr );
	}

      
      //
      // Units
      //

      if ( ! dry_run )
	{
	  // copy existing unit (or if empty, set to '.'
	  std::string ustr = ".";
	  
	  // if unit was a requirement, get the preferred label
	  if ( wild_unit )
	    {
	      if ( wild_unit_pref )
		ustr = Helper::sanitize( Helper::trim( edf.header.phys_dimension[ canonical_slot ] ) );
	      else
		ustr = rule.req_unit.find( "*" )->second;
	    }	    
	  else if ( rule.req_unit.find( sig->unit ) != rule.req_unit.end() )
	    ustr = rule.req_unit.find( sig->unit )->second;
	  else // copy and clean original, if not requirement
	    ustr = Helper::sanitize( Helper::trim( edf.header.phys_dimension[ canonical_slot ] ) );
	  
	  if ( empty_field( ustr ) )
	    ustr = ".";
	  
	  // update EDF header
	  if ( verbose && edf.header.phys_dimension[ canonical_slot ] != ustr )
	    logger << "   changing physical unit from ["
		   << edf.header.phys_dimension[ canonical_slot ]
		   << "] to " << ustr << "\n";
	  edf.header.phys_dimension[ canonical_slot ] = ustr;
	      
	  // secondarily, for voltages, convert to either volts, millivolts, microvolts (set-unit)
	  if ( rule.set_unit != "." ) 
	    if ( ustr == "V" || ustr == "uV" || ustr == "mV" )
	      {
		if ( verbose )
		  logger << "   setting voltage scale to " << rule.set_unit << "\n";
		edf.rescale(  canonical_signal(0) , rule.set_unit , true ); // T -> quietly
	      }
	}

      
      //
      // Transducer field
      //

      if ( ! dry_run )
	{
	  
	  std::string transducer = ".";
	  // if transducer was a requirement, get the preferred label
	  if ( wild_trans )
	    {
	      if ( wild_trans_pref )
		transducer = Helper::sanitize( Helper::trim( edf.header.transducer_type[ canonical_slot ] ) );
	      else
		transducer = rule.req_transducer.find( "*" )->second;
	    }
	  else if ( rule.req_transducer.find( sig->transducer ) != rule.req_transducer.end() )
	    transducer = rule.req_transducer.find( sig->transducer )->second; 
	  else // copy and clean original, if not requirement
	    transducer = Helper::sanitize( Helper::trim( edf.header.transducer_type[ canonical_slot ] ) );
	  
	  if ( empty_field( transducer ) )
	    transducer = ".";
	  
	  if ( verbose && edf.header.transducer_type[ canonical_slot ]  != transducer )
	    logger << "   changing transducer field from ["
		   << edf.header.transducer_type[ canonical_slot ]
		   << "] to " << transducer << "\n";
	  
	  edf.header.transducer_type[ canonical_slot ] = transducer;	  
	}
      
      //
      // Clear pre-filtering field?
      //

      if ( ! dry_run )
	{
	  if ( ! retain_prefiltering )
	    {
	      if ( verbose ) logger << "   clearing prefiltering field\n"; 
	      edf.header.prefiltering[ canonical_slot ] = ".";
	    }
	}
      
      //
      // If keeping existing channel, update label?
      //
      
      if ( already_present && ! dry_run )
	{
	  // nb. not updating header.label_all[] , but this is now in memory
	  // and effectively a new, derived channel, so this is not a problem.
	  // i.e. *should* never be reading this from disk again in any case.
	  
	  edf.header.label[ canonical_slot ] = rule.canonical_label;
	  
	  edf.header.label2header[ rule.canonical_label ] = canonical_slot ;
	}
      
      
      //
      // output
      //
      
      if ( ! only_check_labels )
	{
	  writer.level( rule.canonical_label , "CS" );
	  writer.value( "DEFINED" , 1 );
	  writer.value( "SIG" , matched_sig );
	  if ( matched_ref != "" )
	    writer.value( "REF" , matched_ref );
	}
      
      if ( only_check_labels )
	{
	  retval.okay[ rule.canonical_label ] = true;
	  retval.sig[ rule.canonical_label ] = matched_sig;
	  retval.ref[ rule.canonical_label ] = matched_ref == "" ? "." : matched_ref;
	}
      
      
      //
      // track that we have completed this rule
      //

      completed.insert( rule.canonical_label );
      

      //
      // and update the EDF signal list...
      //
      
      if ( ! already_present ) 
	{

	  if ( dry_run ) 
	    {
	      // TODO: need to pass SR, scale, unit, transducer info about
	      //  dry-run 'dummy' new CS (for subsequent CS rule matching as a sig)
	      canon_edf_signal_t new_sig( rule.canonical_label );					  
	      signals.insert( new_sig );	      
	    }
	  else
	    {
	      canon_edf_signal_t new_sig( edf.header , canonical_slot );	  
	      signals.insert( new_sig );
	    }
	}


      //
      // Verbose output
      //
      
      // logger << "  applied rule: " << rule.canonical_label << "\n"
      // 	     << " unless : " << print( rule.unless ) << "\n"
      // 	     << " group  : " << print( rule.group ) << "\n"
      // 	     << " req sig: " << print( rule.req_sig ) << "\n"
      // 	     << " req ref: " << print( rule.req_ref ) << "\n"
      // 	     << " req trs: " << print( rule.req_transducer ) << "\n"
      // 	     << " req unt: " << print( rule.req_unit ) << "\n"
      // 	     << " req scl: " << rule.req_scale << "\n"
      // 	     << " req sr : " << rule.req_sr_min << " - " << rule.req_sr_max << "\n"
      // 	     << " set sr : " << rule.set_sr << "\n"
      // 	     << " set unt: " << rule.set_unit 
      // 	     << "\n\n";

      
    } // next rule 
  
  
  // tidy if any output
  if ( ! only_check_labels )
    {
      if ( completed.size() != 0 )
	writer.unlevel( "CS" );

      // some additional summaries
      writer.value( "CS_SET" , (int)completed.size() );
      writer.value( "CS_NOT" , (int)(attempted.size() - completed.size()) );
    }
  
  

  if ( verbose )
    logger << "\n  finished processing all rules\n";
      

  
  //
  // report on any failed canonical signals
  //
  
  if ( ! only_check_labels ) 
    {
      bool any_incomplete = false;

      // report on which canonical signals we did not complete
      std::set<std::string>::const_iterator ii = attempted.begin();
      while ( ii != attempted.end() )
	{
	  if ( completed.find( *ii ) == completed.end() )
	    {	      
	      any_incomplete = true;
	      writer.level( *ii , "CS" );
	      writer.value( "DEFINED" , 0 );
	    }
	  ++ii;
	}

      if ( any_incomplete )
	writer.unlevel( "CS" );

      
    }


  

   //
   // Drop original signals?
   //

  if ( ! only_check_labels )
    {

      if ( drop_originals && ! dry_run )
 	logger << "  now dropping all (non-canonical) original signals\n";

      const int ns = osignals.size();
      
      int sigs_used = 0;
      int sigs_unused = 0;
      
      for (int s=0; s<ns; s++)
 	{
 	  const std::string label = osignals.label(s) ;
	  
 	  if ( do_not_drop.find( Helper::toupper( label ) ) == do_not_drop.end() )
 	    {

 	      int slot = edf.header.signal( label );

 	      if ( slot == -1 )
 		Helper::halt( "internal error in edf_t::canonical()" );

 	      if ( drop_originals && ! dry_run )
		edf.drop_signal( slot );
	      
	      bool was_used = used.find( Helper::toupper( label ) ) != used.end() ? 1 : 0 ;
	      
	      // report output
	      writer.level( label , globals::signal_strat );	      
	      if ( drop_originals ) 
		writer.value( "DROPPED" , 1 );	      	      
	      writer.value( "USED" , was_used );
	      
	      if ( was_used ) ++sigs_used;
	      else ++sigs_unused;	      
	      
 	    }
	  else
	    {
	      // presumably always 'used' here.. but do this way just in case...
	      bool was_used = used.find( Helper::toupper( label ) ) != used.end() ? 1 : 0 ;
	      
	      writer.level( label , globals::signal_strat );
              writer.value( "DROPPED" , 0 );
              writer.value( "USED" , was_used );
	      
	      if ( was_used ) ++sigs_used;
	      else ++sigs_unused;	      
	      
	    }
	  writer.unlevel( globals::signal_strat ); 
 	}      

      writer.value( "USED_CH" , sigs_used );
      writer.value( "UNUSED_CH" , sigs_unused );
      
    }
  
}




bool canonical_t::ref_match( const std::vector<std::string> & a , std::set<canon_edf_signal_t> & b , std::string * match )
{
  
  // find first instance in 'a' (req sig list) that is in b (EDF) 
  for (int i=0; i<a.size(); i++)
    {
      // allow for "M1,M2" linked mastoid matching
      // (as prior search was a quoted one
      //std::cout << " checking ref [" << a[i] << "]\n";

      std::vector<std::string> tok = Helper::parse( a[i] , "," );

      // all must match
      bool all_match = true;
      for (int j=0; j<tok.size(); j++)
	{
	  canon_edf_signal_t s1( tok[j] );
	  if ( b.find( s1 ) == b.end() )
	    all_match = false;	  
	  if ( ! all_match ) break;
	}

      if ( all_match )
	{
	  // i.e. set the string as the original M1,M2 version, if a linked ref.
	  *match = a[i];
	  //std::cout << "MATYCHED!\n";
	  return true;
	}
    }
  //  std::cout << "NOT MATYCHED!\n";
  return false;    
}


