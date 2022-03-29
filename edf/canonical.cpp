
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


canon_rule_t::canon_rule_t( const std::vector<std::string> & lines )
{
  const int l = lines.size();

  req_sr_min = req_sr_max = 0;

  // 0 none, -1 all NEG, +1 all POS, +2 = NEG & POS
  req_scale = 0;
  
  std::string current_rule = "";
  
  for (int i=0; i<l; i++)
    {
      const std::string & line = lines[i];
      
      if ( Helper::trim(line) == "" ) continue;
      if ( line[0] == '%' ) continue;
      
      // no indentation: new canonical label
      // will only have one of these
      if ( line[0] != ' ' )
	{
	  canonical_label = line;
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
	  
	  const bool rtype = line[1] != ' ';
	  
	  // the rule type? 
	  if ( rtype )
	    {
	      std::string r = Helper::trim( Helper::toupper( line ) ) ;
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
		    group.insert( canonical_t::swap_in_alias( Helper::trim( tokb[j] ) ) );
		}
	      else if ( current_rule == "unless" )
		{
		  std::vector<std::string> tokb = Helper::quoted_parse( line , "," );
		  for (int j=0; j<tokb.size(); j++)
		    unless.insert( canonical_t::swap_in_alias( Helper::trim( tokb[j] ) ) );
		}
	      else if ( current_rule == "req" )
		{
		  std::vector<std::string> tok = Helper::quoted_parse( line , "=" );
		  if ( tok.size() != 2 )
		    Helper::halt( "expecting key = value format:\n" + line );
		  
		  std::string key = Helper::trim( Helper::toupper( tok[0] ) );
		  std::string value = Helper::trim( tok[1] );
		  std::vector<std::string> tokb = Helper::quoted_parse( value , "," );
		  std::vector<std::string> tok2;

		  // swap in aliases?
		  if ( key == "SIG" || key == "REF" || key == "TRANS" || key == "UNIT" )
		    for (int j=0; j<tokb.size(); j++)
		      {
			std::vector<std::string> tokc =
			  Helper::quoted_parse( canonical_t::swap_in_alias( tokb[j] ) , "," );
			for (int k=0; k<tokc.size(); k++)
			  tok2.push_back( tokc[k] );
		      }
		    
		  // SIG, REF, TRANS, SR-MIN, SR-MAX, UNIT, SCALE
		  // n.b. signals/references are vectors, as order is important
		  
		  if ( key == "SIG" )
		    {
		      for (int j=0; j<tok2.size(); j++)
			req_sig.push_back( Helper::toupper( canonical_t::swap_in_alias( tok2[j] ) ) );
		    }
		  else if ( key == "REF" )
		    {
		      for (int j=0; j<tok2.size(); j++)
                        req_ref.push_back( Helper::toupper( canonical_t::swap_in_alias( tok2[j] ) ) );
		    }
		  else if ( key == "TRANS" )
		    {
		      // n.b. first version cannot be an alias
		      for (int j=0; j<tok2.size(); j++)
                        req_transducer[ Helper::toupper( canonical_t::swap_in_alias( tok2[j] ) ) ] = tok2[0];
		    }
		  else if ( key == "UNIT" )
		    {
		      // n.b. first version cannot be an alias
		      for (int j=0; j<tok2.size(); j++)
			req_unit[ canonical_t::swap_in_alias( Helper::toupper( tok2[j] ) ) ] = tok2[0];
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
		  
                  // SR, UNIT
		  if ( key == "UNIT" )
                    {
		      set_unit = value;
		    }
		  else if ( key == "SR" )
		    {
		      if ( ! Helper::str2int( value , &set_sr ) )
			Helper::halt( "invalid integer sample rate value:\n" + line );		      
		    }
		  else
		    Helper::halt( "did not recognized set value:\n " + line );
		  
		}	      
	    }	  
	}
      
    } // next line

  // all done
}
  
    

canon_edf_signal_t::canon_edf_signal_t( edf_header_t & hdr , const int slot )
{
  if ( slot < 0 || slot > hdr.ns )
    Helper::halt( "bad EDF header slot" );
  
  label = Helper::trim( Helper::sanitize( Helper::toupper( hdr.label[ slot ] ) ) );
  sr = hdr.sampling_freq( slot );
  unit = Helper::trim( Helper::sanitize( Helper::toupper( hdr.phys_dimension[ slot ] ) ) );
  transducer = Helper::trim( Helper::sanitize( Helper::toupper( hdr.transducer_type[ slot ] ) ) );

  scale = 0; // -1, 0, +1
  double phys_min = hdr.physical_min[slot] < hdr.physical_max[slot] ? hdr.physical_min[slot] : hdr.physical_max[slot];
  double phys_max = hdr.physical_min[slot] < hdr.physical_max[slot] ? hdr.physical_max[slot] : hdr.physical_min[slot];
  if ( phys_max < 0 ) scale = -1; // all negative
  else if ( phys_min >= 0 ) scale = 1; // all positive
  
 };


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
      if ( line[0] == '%' ) continue;
      if ( line[0] == '#' ) continue;

      // process any aliases separately
      std::vector<std::string> tok = Helper::quoted_parse( line , "\t " );

      // let XX=A,B,C,D
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
      // process line as part of a rule
      //

      const bool is_name = line[0] != ' ';
      
      if ( is_name )
	{
	  if ( lines.size() == 0 )
	    lines.push_back( line );
	  else
	    rules.push_back( canon_rule_t( lines ) );
	}
      else
	lines.push_back( line );
    }
  
  // flush buffer and add last rule
  if ( lines.size() != 0 )
    rules.push_back( canon_rule_t( lines ) );
  
  // logger << "  read " << rules.size() << " rule(s), and "
  // 	 << aliases.size() << " variables\n";
  
  IN1.close();
  
  return rules.size();
}


canonical_t::canonical_t( edf_t & edf , param_t & param )
{

  //
  // read in rules, if not already done, from 1 or more files
  //
  
  if ( rules.size() == 0 )
    {

      if ( ! param.has( "file" ) )
	Helper::halt( "CANONICAL requires a 'file' argument" );
      
      std::vector<std::string> filenames = param.strvector( "file" ) ;
      
      for (int i=0; i<filenames.size(); i++)
	{
	  std::string filename = Helper::expand( filenames[i] );
	  int nrules = read( filename );
	  logger << "  read " << nrules << " rules from " << filename << "\n";
	}

      logger << "  in total, read " << rules.size()
	     << " rules and " << aliases.size() << " variables\n\n";
    }
  

  //
  // Other options
  //

  if ( param.has( "group" ) )
    group = param.strset( "group" );
  
  drop_originals = param.yesno( "drop-originals" );
  
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
  // track completed canonical signals
  //

  std::set<std::string> completed;
  
  for (int r=0; r<rules.size(); r++)
    {
            
      const canon_rule_t rule = rules[r];
      
      //
      // has this rule already been satisfied?
      //
      
      const bool already_processed = is_in( rule.canonical_label , completed );
      
      if ( already_processed ) continue;
      
      //
      // group specifier?
      //
      
      if ( rule.group.size() != 0 )
       	{
       	  // if a group has been specified for this rule, then
       	  // we need to match it (with at least one of the rules)
       	  if ( ! is_in( group , rule.group ) )
       	    continue;
       	}
      

      //
      // can we find a matching primary signal?
      //

      std::string matched_sig = "";
      
      bool matched = match( rule.req_sig , signals , &matched_sig );

      if ( ! matched ) continue;
      
      std::set<canon_edf_signal_t>::const_iterator sig = signals.find( matched_sig );
      if ( sig == signals.end() ) Helper::halt( "internal error in canonical match sig" );
      
      //
      // a reference signal, if required
      //
      
      std::string matched_ref = "";
 
      if ( rule.req_ref.size() != 0 )
	{
	  
	  matched = match( rule.req_ref , signals , &matched_ref );
	  
	  if ( ! matched ) continue;
	  
	}
      
      //
      // transducer rules?
      //

      if ( rule.req_transducer.size() != 0 )
	{
	  if ( rule.req_transducer.find( sig->transducer ) == rule.req_transducer.end() )
	    continue;
	}


      //
      // unit rules?
      //

      if ( rule.req_unit.size() != 0 )
	{
          if ( rule.req_unit.find( sig->unit ) == rule.req_unit.end() )
            continue;
	}

      //
      // sample rate rules?
      //

      if ( rule.req_sr_min != 0 )
	{
	  if ( sig->sr < rule.req_sr_min )
	    continue;
	}

      if ( rule.req_sr_max != 0 )
	{
	  if ( sig->sr > rule.req_sr_max )
	    continue;
	}

      //
      // Scale?
      //

      if ( rule.req_scale != 0 )
	{
	  if ( rule.req_scale != sig->scale ) continue;
	}

      //
      // If here, we have a match
      //
      
      logger << "  matced rule for " << rule.canonical_label << "\n";

      //
      // Construct the CS
      //

      
      std::cout << rule.canonical_label << "\n"
		<< " unless : " << print( rule.unless ) << "\n"
		<< " group  : " << print( rule.group ) << "\n"
		<< " req sig: " << print( rule.req_sig ) << "\n"
		<< " req ref: " << print( rule.req_ref ) << "\n"
		<< " req trs: " << print( rule.req_transducer ) << "\n"
		<< " req unt: " << print( rule.req_unit ) << "\n"
		<< " req scl: " << rule.req_scale << "\n"
		<< " req sr : " << rule.req_sr_min << " - " << rule.req_sr_max << "\n"
		<< " set sr : " << rule.set_sr << "\n"
		<< " set unt: " << rule.set_unit 
		<< "\n\n";
      
    }

  
  return;
  
}






// END OF NEW SECTION
// --------------------------------------------------------------------------------

//
// Canonical file format
//
// whitespace delimited; one rule per line; spaces can be put in quotes
//
// primary rules:
//
//  group canonical primary (reference) (set-SR) (set-unit) (opt-notes)
// 
// the 'canonical' label can be in the form label|transducer,
//  meaning that the channel name is set to 'label' and the EDF header transducer field is set to 'transducer'
//  no white-space allowed in canonical channel names or transducer types (i.e. use underscores)
//
// comment lines: start with % in first character
//
// aliasis lines: on-the-fly variables, based on a '=' sign in the single field
// additional lines append rather than overwrite
//
// a=v1,v2
// a=v3
//
// means a is 'v1,v2,v3' (i.e. this will be swapped in to subsequent rules below)
//
// requirements:
//  - specifictions for a canonical label|transducer pair, or label implying for all label|(transducers)
//  - i.e. these requirements must be met /for the primary signal/
//
// + label|transducer SR 10 100
// + label|transducer  SR 10 .
// + label : SR 10 .

//  --> ? match anything, but must be non-null
//      * match anything, including null (space/empty)
// + label|transducer1 : UNIT mV,V
// + label|transducer1 : UNIT uV,?
// + label|transducer2 : UNIT mmH2O

// + label|transducer1: POSITIVE
// + label|transducer1: RANGE  0 .
// + label: RANGE  0 100

// make_canonical() arguments:

// if group is != '.' then only attach rows matching
// if the fist column has . then it matches any value of group;
// but if a group is specified, then only match for that
// i.e. group-specific rules should be put before generic rules

 


cansigs_t edf_t::make_canonicals( const std::vector<std::string> & files,
				  const std::string &  group , 
				  const bool make_signals , 
				  const bool drop_originals , 
				  const std::string & prefix , 
				  const std::set<std::string> * cs ,
				  const bool only_check_labels )
{  

  cansigs_t retval;

  if ( files.size() == 0 ) return retval;
    
  // GROUP   CANONICAL|trans   CH   REF   SR  UNITS    NOTES
  // alias=v1,v2,v3
  
  // if cs is non-null, only make the CS in that set ('EEG')

   
  // can have multiple versions of a rule, will pick the first match
  // nb. use of first group column:: if a group is specified on the command line,
  //     then will be preferentially picked
  
  //  STUDY_A  csEEG   EEG3         .       100
  //  STUDY_B  csEEG   C4-REF       M1-REF  100  
  //  .        csEEG   C4,EEG1      M1,A1   100  
  //  .        csEEG   C3,EEG2      M2,A2   100  Using C3, not C4
  //  .        csEEG   C4_A1,C4_M1  .       100
  //  .        csEEG   C3_A2,C3_M2  .       100  Using C3, not C4
  
  // the to-be-created CS go here:
  std::vector<std::string> canons;
  std::set<std::string> canons_uniq;
  
  // if we are dropping originals later, get a list of those now
  const bool only_data_signals = true;
  signal_list_t osignals = header.signal_list( "*" , only_data_signals );
  
  // but, if existing signal is canonical already, flag not to drop
  // (this is allowed if no other transformations, i.e. keep as is)

  std::set<std::string> do_not_drop;
  
  // also, track whether a signal was used or not (i.e. might
  // be dropped, 'C3' but 'C3' --> 'C3-M1' in which was it was
  // still 'used' in the canonical set (versus original channels
  // that basically do not feature at all)  

  std::set<std::string> used;
  
  // if group is non-null, do not process any generic rules after coming
  // across some group-specific rules (whether these worked or not)
  // i.e. if over-riding, must specify the full complement

  std::set<std::string> ignore_generics;

  // on-the-fly defined aliases (w/in the canonical sigs file)
  // n.b. these append, so second line adds to first, not overwrites
  //  M1=M1-ref,M1_ref
  //  M1=A1,A1_ref,A1-ref

  std::map<std::string,std::vector<std::string> > aliases;

  // requirements for canonical signals (i.e. these must be
  // met by a putative signal; they can be specific to a
  // particular label / transducer type pair)

  std::map<std::pair<std::string,std::string>,double> req_sr_min;
  std::map<std::pair<std::string,std::string>,double> req_sr_max; 

  // also includes unit-label --> harm-unit mapping at end
  std::map<std::pair<std::string,std::string>,std::map<std::string,std::string> > req_unit;
  std::map<std::pair<std::string,std::string>,std::set<std::string> > req_trans;
  
  std::map<std::pair<std::string,std::string>,int> req_scale;  // -1, 0 , +1 


  //
  // read in definitions
  //
  
  std::map<std::string,std::vector< std::vector<std::string> > > sigs, refs;
  std::map<std::string,std::vector<std::string> > srs, notes, units, trns;
  
  //
  // iterate over files
  //
  
  for (int f = 0 ; f < files.size() ; f++ )
    {

      std::string file = Helper::expand( files[f] );
      
      if ( ! Helper::fileExists( file ) )
	Helper::halt( "could not find " + file );
      
      std::ifstream IN1( file.c_str() , std::ios::in );
      while ( ! IN1.eof() )
	{
	  std::string line = "";
	  Helper::safe_getline( IN1 , line );
	  
	  if ( line == "" ) continue;
	  if ( Helper::toupper( line ) == "STOP" ) break;
	  if ( IN1.eof() ) break;
	  if ( line.size() >= 1 && line[0] == '%' ) continue;
	  
	  std::cout << " line[" << line << "]\n";

	  // white-space quoted parsing of this line
	  std::vector<std::string> tok = Helper::quoted_parse( line , "\t " );
	  if ( tok.size() == 0 ) continue;
	  
	  // is this an aliases line? if so, add then skip
	  if ( tok.size() == 1 )
	    {
	      std::vector<std::string> tok2 = Helper::quoted_parse( tok[0] , "=" );
	      if ( tok2.size() != 2 )		
		Helper::halt( "expecting 1 column for ch=x,y,z aliases fields:\n" + line );
	      std::vector<std::string> tok3 = Helper::quoted_parse( tok2[1] , "," );
	      for (int a=0; a<tok3.size(); a++) aliases[ tok2[0] ].push_back( tok3[a] ) ;
	      continue;
	    }

	  // is this a requirement line? (starts with '+') if so, add then skip
	  if ( tok[0] == "+" )
	    {
	      if ( tok.size() < 4 ) Helper::halt( "bad format for: " + line );
	      std::vector<std::string> label_trans = Helper::quoted_parse( tok[1] , "|" );

	      std::string label = label_trans[0];
	      std::string trans = label_trans.size() == 2 ? label_trans[1] : "." ; 
	      std::cout << " labe, trans = " << label << " " << trans << "\n";
	      if ( label_trans.size() > 2 )
		Helper::halt( "bad format label|transducer" );

	      // rule type
	      std::string req = Helper::toupper( tok[2] );

	      // . means null
	      // ? means non-null

	      // + VAL|TRNS TRANS .

	      // + VAL    UNIT  canonical,alias1,alais2,...
	      
	      if ( req == "UNIT" )
		{
		  // slot 3 (i==3) will be the preferred unit label
		  if ( tok[3] == "." && tok.size() > 4 )
		    Helper::halt( "cannot specify '.' as the primary unit if others specified" );
						     
		  for (int i=3; i<tok.size(); i++)
		    {
		      //std::cout << " [" << tok[i]  << "] [" << tok[3] << "]...\n";
		      req_unit[ std::make_pair( label , trans ) ][ Helper::toupper( tok[i] ) ] = tok[3] ;  
		    }
		}
	      
	      else if ( req == "SR" )
		{
		  if ( tok.size() != 5 )
		    Helper::halt( "bad format for: " + line );

		  // min SR?
		  if ( tok[3] != "." )
		    {
		      
		      int sr;
		      if ( ! Helper::str2int( tok[3] , &sr ) )
			Helper::halt( "bad format for: " + line );
		      std::cout << "adding rule = " << label << " | " << trans << "\n";
		      req_sr_min[ std::make_pair( label , trans ) ] = sr;
		    }

		  // max SR?
		  if ( tok[4] != "." )
		    {
		      int sr;
		      if ( ! Helper::str2int( tok[4] , &sr ) )
			Helper::halt( "bad format for: " + line );
		      req_sr_max[ std::make_pair( label , trans ) ] = sr;
		    }
	
		}

	      else if ( req == "SCALE" )
		{
		  if ( Helper::toupper( tok[3] ) == "POSNEG" )
		    req_scale[ std::make_pair( label , trans ) ] = 0; // requires POS /and/ NEG range
		  else if ( Helper::toupper( tok[3] ).substr(0,3) == "POS" )
		    req_scale[ std::make_pair( label , trans ) ] = +1; // requires only POS
		  else if ( Helper::toupper( tok[3] ).substr(0,3) == "NEG" )
		    req_scale[ std::make_pair( label , trans ) ] = -1; // requires only NEG
		  else if  ( Helper::toupper( tok[3] ) == "AC" ) // == POSNEG
		    req_scale[ std::make_pair( label , trans ) ] = 0;
		  else
		    Helper::halt( "bad value for scale:\n" + line );
		}
	      
	      else if ( req == "TRANS" )
		{
		  for (int i=3; i<tok.size(); i++)
		    req_trans[ std::make_pair( label , trans ) ].insert( Helper::toupper( tok[i] ) );
		}
	      else
		Helper::halt( "did not recognize requirement: " + req );
	      
	      // skip to next line
	      continue;
	    }
	  

	  // otherwise, a rule requires normal columns
	  if ( tok.size() < 3 && tok.size() > 7 ) 
	    Helper::halt( "expecting 3 to 7 white-delimited columns (use quotes for values with spaces)\nfile: " 
			  + file + "\nline: [" + line + "]\n" );
	  
	  
	  // ensure we have full set - populating w/ '.' as needed
	  if ( tok.size() < 6 ) tok.resize( 6 , "." );
	  
	  // ignore group-specific rules, unless that group has been specifically requested
	  if ( tok[0] != "." && tok[0] != group ) continue;

	  // process transducer type:  canonical|trans_text
	  std::string transducer = ".";
	  std::vector<std::string> st = Helper::parse( tok[1] , "|" );
	  if ( st.size() > 2 ) Helper::halt( "expecting canonical|transducer" );
	  if ( st.size() == 2 )
	    {
	      tok[1] = st[0]; // extract label only
	      transducer = st[1];
	    }
	  
	  // track that we are seeing a matching group-specific rule
	  if ( group != "." && tok[0] == group ) ignore_generics.insert( tok[1] );
	  
	  // ignore any generic rules for a canonical signal if we have
	  // already encountered a matching group-specific rule for that
	  // canonical signal
	  if ( group != "." && tok[0] == "."
	       && ignore_generics.find( tok[1] ) != ignore_generics.end() )
	    continue;
	  
	  // if cs not specified, take all canonical signals as given in 
	  // the file      
	  if ( cs == NULL || cs->find( tok[1] ) != cs->end() )
	    {
	      // only add once to canons[], but preserved file order (first seen)
	      if ( canons_uniq.find( prefix + tok[1] ) == canons_uniq.end() )
		{
		  canons.push_back( prefix + tok[1] );
		  canons_uniq.insert( prefix + tok[1] ); 
		}
	      
	    }
	  
	  // skip if a specific list requested?
	  if ( cs != NULL && cs->find( tok[1] ) == cs->end() ) continue;
	  
	  // otherwise, add to the set of things to be calculated
	  sigs[ prefix + tok[1] ].push_back( Helper::quoted_parse( tok[2] , "," ) );
	  refs[ prefix + tok[1] ].push_back( Helper::quoted_parse( tok[3] , "," ) );
	  trns[ prefix + tok[1] ].push_back( transducer );
	  srs[ prefix + tok[1] ].push_back( tok[4] ) ;
	  units[ prefix + tok[1] ].push_back( tok[5] );
	  notes[ prefix + tok[1] ].push_back( tok.size() == 7 ? tok[6] : "." );
	  
	}
      
      IN1.close();
      
      // read in the next file
    }
  
  //
  // no valid rules found?
  //
  
  if ( sigs.size() == 0 ) 
    Helper::halt( "no valid rules (given group " + group + ")" );
  
    
  //
  // For each canonical signal (ie the order in which we encountered them)
  //
  
  std::vector<std::string>::const_iterator cc = canons.begin();
  while ( cc != canons.end() )
    {
      
      std::string canon = *cc;
      
      if ( ! only_check_labels ) 
	writer.level( canon , "CS" );
      else
	retval.okay[ canon ] = false;
      
      if ( sigs.find( canon ) == sigs.end() )
	{
	  if ( ! only_check_labels ) 
	    writer.value( "DEFINED" , 0 );
	  ++cc;
	  continue;
	}
      
      //
      // check whether canonical form already exists
      //  - we allow this, but ONLY if no transformations
      //    are requested (i.e. no resampling etc)
      //    meaning we basically just leave as is
      //
      
      bool already_present = header.has_signal( canon );      
	  
      // as soon as we find a matching rule, we stop
      bool done = false;
      
      const int n_rules = sigs.find( canon )->second.size() ; 
      
      for (int j=0; j<n_rules; j++ )
	{
	  
	  // find best choice of signals
	  std::string sigstr = "";
	  std::vector<std::string> v = sigs.find( canon )->second[j];
	  for (int k=0; k<v.size(); k++)
	    {
	      
	      if ( header.has_signal( v[k] )  ) // case-insensitive match 
		{
		  sigstr = v[k];
		  break;
		}

	      // or, is this an alias, in which case check (in order)
	      // the equivalent terms
	      if ( aliases.find( v[k] ) != aliases.end() )
		{
		  std::vector<std::string> av =  aliases.find( v[k] )->second;
		  for (int ak=0; ak<av.size(); ak++)
		    {
		      if ( header.has_signal( av[ak] )  ) 
			{
			  sigstr = av[ak];
			  break;
			}
		    }
		}
	    }
	  
	  if ( sigstr == "" ) 
	    continue;
	  
	  //
	  // Reference: nb. if we get quoted "A1,A2", then check both and flag to take average of those
	  //
	  
	  std::string refstr = "";
	  v = refs.find( canon )->second[j];
	  if ( v.size() == 1 && v[0] == "." ) 
	    refstr = ".";
	  else
	    {
	      for (int k=0; k<v.size(); k++)
		{
		  //std::cout << " RESTING [" << v[k] << "]\n";
		  if ( v[k] == "." ) 
		    Helper::halt( "cannot mix '.' and non-'.' references" );

		  // if this quoted (averagig?)
		  if ( v[k][0] == '"' )
		    {
		      std::string reftarget = Helper::unquote( v[k] );
		      // we need to find /all/ members of this "ref1,ref2,..." term
		      // typically, for PSG, expectring this will just be two (linked mastoids)
		      std::vector<std::string> rtok = Helper::parse( reftarget , "," );

		      std::vector<std::string> refstr_vec;
			
		      bool okay = true;
		      for (int r=0; r<rtok.size();r++)
			{
			  bool match1 = false;
			  // self-match?
			  if ( header.has_signal( rtok[r] ) )
			    {
			      refstr_vec.push_back( rtok[r] );
			      match1 = true;
			    }
			  else
			    {
			      if ( aliases.find( rtok[r] ) != aliases.end() )
				{
				  //std::cout << "found aliases for " << rtok[r]  << "\n";
				  std::vector<std::string> av =  aliases.find( rtok[r] )->second;
				  for (int ak=0; ak<av.size(); ak++)
				    {
				      //std::cout << " --> " << av[ak] << "\n";
				      if ( header.has_signal( av[ak] )  )
					{
					  match1 = true;
					  refstr_vec.push_back( av[ak] );
					  break;
					}
				    }
				}
			    }
			  
			  // if this one term does not match, then the whole match is bad
			  if ( ! match1 ) okay = false;
			}
		      
		      // found a match here -- pass in the whole comma-delimited string (i.e. this will work w/ reference() )
		      if ( okay )
			{
			  refstr = "";
			  for (int i=0;i<refstr_vec.size(); i++)
			    {
			      if ( i ) refstr += ",";
			      refstr += refstr_vec[i];
			    }
			  break;
			}
		      
		    }
		  else // not averaged - same checks as above
		    {
		      
		      if ( header.has_signal( v[k] ) ) // case-insensitve, alias-aware match
			{
			  refstr = v[k];
			  break;
			}
		  
		      // or, is this an alias, in which case check (in order)
		      // the equivalent terms
		      
		      if ( aliases.find( v[k] ) != aliases.end() )
			{
			  std::vector<std::string> av =  aliases.find( v[k] )->second;
			  for (int ak=0; ak<av.size(); ak++)
			    {
			      if ( header.has_signal( av[ak] )  )
				{
				  refstr = av[ak];
				  break;
				}
			    }
			}
		    }
		}      
	    }

	  //
	  // matches?
	  //
	  
	  if ( sigstr == "" || refstr == "" )
	    continue;



	  // --------------------------------------------------------------------------------
	  //
	  // does this confirm to any other requirements?
	  //
	  // --------------------------------------------------------------------------------
	
	  // get info on the putative primary signal (not reference)
	  signal_list_t putative = header.signal_list( sigstr );	  
	  std::string canon_transducer = trns.find( canon )->second[j];
	  std::pair<std::string,std::string> key = std::make_pair( canon , canon_transducer );

	  logger << "  putative match for " << canon << " (transducer type " << canon_transducer << ")\n";

	  // transducer?
	  if ( req_trans.find( key ) != req_trans.end() )
	    {
	      const std::string edf_transducer = Helper::toupper( header.transducer_type[ putative(0) ] );
	      const bool empty_field = Helper::trim( edf_transducer ) == "" ;

	      const std::set<std::string> & req_transducers = req_trans[ key ];
	      bool matches = false;
	      std::set<std::string>::const_iterator tt = req_transducers.begin();
	      while ( tt != req_transducers.end() )
		{
		  if ( edf_transducer == *tt ) { matches = true; break; }
		  if ( empty_field && *tt == "." ) { matches = true; break; }
 		  ++tt;
		}
	      if ( ! matches )
		{
		  logger << "  did not meet transducer type requirement: EDF = " << edf_transducer << "\n";
		  continue;
		}
	    }

          // physical dimension
          if ( req_unit.find( key ) != req_unit.end() )
            {
              const std::string edf_phys_dim = Helper::toupper( header.phys_dimension[ putative(0) ] );
              const bool empty_field = Helper::trim( edf_phys_dim ) == "" ;
	      
              const std::map<std::string,std::string> & req_phys_dim = req_unit[ key ];
              bool matches = false;
              std::map<std::string,std::string>::const_iterator tt = req_phys_dim.begin();
              while ( tt != req_phys_dim.end() )
                {
                  if ( edf_phys_dim == tt->first ) { matches = true; break; }
                  if ( empty_field && tt->first == "." ) { matches = true; break; }
                  ++tt;
                }
              if ( ! matches )
		{
		  logger << "  did not meet physical dimension requirement: EDF = " << edf_phys_dim << "\n"; 
		  continue;
		}
	      
            }

	  	 
	  // min sample rate?
	  if ( req_sr_min.find( key ) != req_sr_min.end() )
            {
              const double sr = header.sampling_freq( putative(0) );
	      if ( sr < req_sr_min[ key ] )
		{
		  logger << "  did not meet minimum SR requirement: EDF SR = " << sr << "\n"; 
		  continue;
		}
	    }

	  // max sample rate?
	  if ( req_sr_max.find( key ) != req_sr_max.end() )
            {
              const double sr = header.sampling_freq( putative(0) );
	      if ( sr > req_sr_max[ key ] )
		{
		  logger << "  did not meet maximum SR requirement: EDF SR = " << sr << "\n"; 
		  continue;
		}
	    }

	  // scale 
	  if ( req_scale.find( key ) != req_scale.end() )
            {
	      const double edf_phys_min = header.physical_min[ putative(0) ];
	      const double edf_phys_max = header.physical_max[ putative(0) ];
	      // -1 all neg
	      // +1 only pos
	      // 0  requires both NEG and POS
	      const int scale = req_scale[ key ];
	      const bool has_neg = edf_phys_min < 0 ;
	      const bool has_pos = edf_phys_max > 0 ;
	      
	      if ( scale == 0 && ! ( has_neg && has_pos ) )
		{
		  logger << "  did not meet SCALE POSNEG/AC requirement: EDF min/max = "
			 << edf_phys_min << " - " << edf_phys_max << "\n";
		  continue;
		}

	      if ( scale == 1 && has_neg )
		{
                  logger << "  did not meet SCALE POS requirement: EDF min/max = "
                         << edf_phys_min << " - " << edf_phys_max << "\n";
                  continue;
                }

	      if ( scale == -1 && has_pos )
		{
                  logger << "  did not meet SCALE NEG requirement: EDF min/max = "
                         << edf_phys_min << " - " << edf_phys_max << "\n";
                  continue;
                }

	    }


	  //
	  // ---- this rule has passed ... now add -----
	  //

	  // if ( already_present && refstr != "." )
	  //   Helper::halt( "cannot specify existing canonical name "
	  // 		  + canon + " and a re-reference" );
	  
	  if ( ! only_check_labels ) 
	    logger << "  generating canonical signal " << canon 
		   << " from " << sigstr << "/" << refstr << "\n";
	  
	  //
	  // track that we are using these channels
	  //
	  
	  used.insert( Helper::toupper( sigstr ) );
	  used.insert( Helper::toupper( refstr ) );

	  retval.used.insert( Helper::toupper( sigstr ) );
	  retval.used.insert( Helper::toupper( refstr ) );
	  
	  
	  //
	  // Track that the original should not be dropped, as it features
	  //
	  
	  if ( drop_originals && already_present )
	    do_not_drop.insert( Helper::toupper( canon ) );

	  //
	  // Sample rate changes?
	  //

	  std::string srstr = srs.find( canon )->second[j];
	  
	  // if SR == '.' means do not change SR
	  
	  int sr = 0;
	  if ( srstr != "." )
	    if ( ! Helper::str2int( srstr , &sr ) )
	      Helper::halt( "non-integer SR for " + canon );
	  
          // if ( already_present && srstr != "." )
          //   Helper::halt( "cannot specify existing canonical name "
          //                 + canon + " and re-sample" );

	  
	  // copy signal --> canonical form
	  signal_list_t ref;
	  if ( refstr != "." ) ref = header.signal_list( refstr );
	  
	  signal_list_t sig = header.signal_list( sigstr );


	  //
	  // Generate the new signal (w/ re-referencing optionally)
	  //   false , false = not dereference , not verbose
	  // 

	  //std::cout << "MS = " << make_signals << "\n";
	  
	  if ( make_signals )
	    {
	      
	      // create a new channel?

	      if ( ! already_present )
		{
		  reference( sig ,    // original channel
			     ref ,    // reference (or "." for none)
			     true ,   // create a new channel? 
			     canon ,  // new channel name
			     sr ,     // new channel SR
			     false ,  // do not de-reference		       
			     false ); // not in verbose mode
		}
	      
	    }
	  
	  //
	  // Get the canonical signal
	  //
	  
	  signal_list_t canonical_signal = header.signal_list( canon );


	  //
	  // re-reference existing channel?
	  //
	  
	  // or modify an existing one?
	  if ( refstr != "." && already_present )
	    {
	      // do not make a new channel
	      reference( sig , ref , false , "" , 0 );
	    }

	  //
	  // resample an existing channel? (i.e. if not already done above)
	  //

	  if ( sr != 0 && already_present )
	    {
	      dsptools::resample_channel( *this , canonical_signal(0) , sr );
	    }
	  
	  //
	  // rescale units? ignore if other than V/uV/mV
	  //
	  
	  std::string ustr = units.find( canon )->second[j];
	  
	  // if ( already_present && ustr != "." )
          //   Helper::halt( "cannot specify existing canonical name "
          //                 + canon + " and transform units" );
	  
	  if ( make_signals && req_unit.find( key ) != req_unit.end() )
            {
	      // get implied canonical unit
	      std::string can_unit = Helper::toupper( header.phys_dimension[ canonical_signal(0) ] );
	      // get the preferred value	      
	      std::string pref_unit = req_unit[ key ][ can_unit ];
	      logger << "  swapping " << can_unit << " to " << pref_unit << "\n";

	      // set the new EDF header
	      header.phys_dimension[ canonical_signal(0) ] = pref_unit;
	    }
	  
	  // convert from volts, millivolts, microvolts
	  if ( make_signals ) 
	    if ( ustr == "V" || ustr == "uV" || ustr == "mV" ) 
	      rescale(  canonical_signal(0) , ustr );

	  //
	  // add transducer label?
	  //

	  std::string transducer = trns.find( canon )->second[j];
	  if ( make_signals && transducer != "." )
	    header.transducer_type[ canonical_signal(0) ] = transducer;	  
	  
	  //
	  // If keeping existing channel, update label?
	  //

	  if ( already_present )
	    {
	      // nb. not updating header.label_all[] , but this is now in memory
	      // and effectively a new, derived channel, so this is not a problem.
	      // i.e. *should* never be reading this from disk again in any case.

	      header.label[ canonical_signal(0) ] = canon;

	      header.label2header[ canon ] = canonical_signal(0) ;
	    }
	  
	  //
	  // output
	  //

	  if ( ! only_check_labels )
	    {
	      writer.value( "DEFINED" , 1 );
	      writer.value( "SIG" , sigstr );
	      writer.value( "REF" , refstr );
	      writer.value( "SR" , srstr );
	      writer.value( "UNITS" , ustr );
	      
	      std::string notesstr = notes.find( canon )->second[j];
	      
	      if ( notesstr != "" )
		writer.value( "NOTES" , notesstr );
	    }


	  if ( only_check_labels )
	    {
	      retval.okay[ canon ] = true;
	      retval.sig[ canon ] = sigstr;
	      retval.ref[ canon ] = refstr;
	    }
	  
	  //
	  // at this point, rule was found, so quit
	  //

	  done = true;
	  
	  break;
	
	} // next rule for this CS


      // if we failed to get a match, flag here

      if ( ! only_check_labels ) 
	if ( ! done ) 
	  writer.value( "DEFINED" , 0 );

      
      // next CS
      
      ++cc;
      
    } 
  
  if ( ! only_check_labels )
    writer.unlevel( "CS" );


  //
  // Drop original signals?
  //

  if ( drop_originals && ! only_check_labels )
    {
      if ( make_signals ) 
	logger << "  now dropping the original signals\n";

      const int ns = osignals.size();
      for (int s=0; s<ns; s++)
	{
	  const std::string label = osignals.label(s) ;
	  
	  if ( do_not_drop.find( Helper::toupper( label ) ) == do_not_drop.end() )
	    {

	      int slot = header.signal( label );

	      if ( slot == -1 )
		Helper::halt( "internal error in edf_t::canonical()" );

	      if ( make_signals ) 
		drop_signal( slot );

	      // report output
	      writer.level( label , globals::signal_strat );	      
	      writer.value( "DROPPED" , 1 );	      	      
	      writer.value( "USED" , used.find( Helper::toupper( label ) ) != used.end() ? 1 : 0 ) ; 
	      
	    }
	  writer.unlevel( globals::signal_strat ); 
	}
    }
  
  
  return retval;
}


void edf_t::guess_canonicals( param_t & param , bool make_signals )
{

  // get available signals
  std::string signal_label = param.requires( "sig" );    
  signal_list_t signals = header.signal_list( signal_label );    
  const int ns = signals.size();

  int eeg = -1;
  int eeg_ref = -1 , eeg_ref2 = -1;   // possible two refs, for CZ / (M1+M2)/2
  
  // specified EEG: quick fix, e.g. for EEG3 
  bool has_label = param.has( "eeg" );
  if ( has_label )
    {
      // INCORRECT: need to look up w.r.t. signals rather than EDF header
      //eeg = header.signal( param.value( "eeg" ) );

      std::string ts = param.value( "eeg" );

      eeg = -1;

      for (int s=0; s<ns; s++)
	{
	  if ( signals.label(s) == ts ) { eeg = s ; break; } 
	}

      if ( eeg == -1 ) logger << "  could not find " <<  param.value( "eeg" ) << " -- will try to guess csEEG\n";
    }

  // still not found? then guess
  if ( eeg == -1 ) 
    {
      // EEG: central electrode (and reference)
      int c3 = -1 , m1 = -1 , a1 = -1 , c4 = -1 , a2 = -1 , m2 = -1 , cz = -1 ; 
      
      for (int s=0; s<ns; s++)
	{
	  // C3
	  if ( Helper::contains( signals.label(s) , "C3" ) )
	    { if ( c3 != -1 ) Helper::halt( "matched more than one C3" ); else c3 = s; }

	  // C4
	  if ( Helper::contains( signals.label(s) , "C4" ) )
	    { if ( c4 != -1 ) Helper::halt( "matched more than one C4" ); else c4 = s; }

	  // A1
	  if ( Helper::contains( signals.label(s) , "A1" ) )
	    { if ( a1 != -1 ) Helper::halt( "matched more than one A1" ); else a1 = s; }
	  
	  // A2
	  if ( Helper::contains( signals.label(s) , "A2" ) )
	    { if ( a2 != -1 ) Helper::halt( "matched more than one A2" ); else a2 = s; }

	  // M1
	  if ( Helper::contains( signals.label(s) , "M1" ) )
	    { if ( m1 != -1 ) Helper::halt( "matched more than one M1" ); else m1 = s; }
	  
	  // M2
	  if ( Helper::contains( signals.label(s) , "M2" ) )
	    { if ( m2 != -1 ) Helper::halt( "matched more than one M2" ); else m2 = s; }

	  // Cz
	  if ( Helper::contains( signals.label(s) , "CZ" ) )
	    { if ( cz != -1 ) Helper::halt( "matched more than one CZ" ); else cz = s; }

	}

      // collapse mastoid channels
      if ( a1 != -1 && m1 != -1 ) Helper::halt( "both A1 and M1 present" );
      if ( a2 != -1 && m2 != -1 ) Helper::halt( "both A2 and M2 present" );
      if ( m1 == -1 ) m1 = a1;
      if ( m2 == -1 ) m2 = a2;
      	   
      // C4/M1
      if      ( c4 == m1 && c4 != -1 ) eeg = c4;
      // C4 and M1 (separate ref.)
      else if ( c4 != -1 && m1 != -1 ) { eeg = c4; eeg_ref = m1; }
      // C4 as a single channel
      else if ( c4 != -1 ) { eeg = c4; }
      // C3/M2
      else if ( c3 == m2 && c3 != -1 ) eeg = c3;
      // C3 and M2 (separate ref.)
      else if ( c3 != -1 && m2 != -1 ) { eeg = c3; eeg_ref = m2; }
      // C3 as a single channel
      else if ( c3 != -1 ) { eeg = c3; }

      // else CZ, linked mastoid
      else if ( cz != -1 && m1 != -1 && m2 != -1 ) { eeg = cz ; eeg_ref = m1 ; eeg_ref2 = m2; }
      // else just CZ
      else if ( cz != -1 ) eeg = cz;

    }


  //
  // At this point, we should have eeg, eeg_ref and eeg_ref2 completed
  //

  if ( eeg == -1 )
    {
      logger << "  could not guess the canonical csEEG from any matching signals\n";

      std::string canon = "csEEG";
      writer.level( canon , "CS" );      
      writer.value( "DEFINED" , 0 );
      writer.unlevel( "CS" );

      return;
    }
      

  //
  // Sample rate (default 100 Hz)
  //

  const int sr = param.has( "sr" ) ? param.requires_int( "sr" ) : 100;

  //
  // Note -- this duplicates code below, but keep for now...
  //

  std::string sigstr = signals.label( eeg );

  signal_list_t sig = header.signal_list( sigstr );
  
  std::string refstr = ".";
  if ( eeg_ref2 != -1 ) refstr = signals.label( eeg_ref ) + "," + signals.label( eeg_ref2 );
  else if ( eeg_ref != -1 ) refstr = signals.label( eeg_ref );

  signal_list_t ref;
  if ( refstr != "." ) ref = header.signal_list( refstr );

  logger << "  creating csEEG using signal [" << sigstr << "] and reference [" << refstr << "]\n";

  
  //
  // Rerefence and make canonical signal
  //
  
  std::string canon = "csEEG";

  if ( make_signals ) 
    reference( sig , ref , true , canon , sr );
  
  signal_list_t canonical_signal = header.signal_list( canon );

  
  //
  // rescale units?
  //
  
  // EEG in uV
	  
  std::string units = "uV";
  
  if ( make_signals ) 
    if ( units == "uV" || units == "mV" ) 
      rescale( canonical_signal(0) , units );
  
  

  //
  // output
  //
  
  writer.level( canon , "CS" );
  
  writer.value( "DEFINED" , 1 );
  writer.value( "SIG" , sigstr );
  writer.value( "REF" , refstr );
  writer.value( "SR" , sr );
  writer.value( "UNITS" , units );
  
  writer.unlevel( "CS" );
  
}

