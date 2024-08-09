
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

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <sstream>
#include <cstring>

#include "luna.h"

// #include "defs/defs.h"
// #include "helper/helper.h"
// #include "db/db.h"

extern writer_t writer;
extern globals global;

struct options_t {
  bool print_empty_rows; // -e 
  bool print_cmd_name;   // -n  
  bool long_format;      // -l
  double prec; // -p 
  std::string prepend; // -a
  bool full;   // -f
  bool show_progress; // -
  bool cmd_hash;
  
  char strata_delim;
  char faclvl_delim;
  
  options_t() 
    :
    print_empty_rows( false ) , 
    print_cmd_name( false ) ,
    long_format( false ) ,
    prec(3) ,
    prepend( "" ), 
    full( true ) ,
    cmd_hash( false ) ,
    strata_delim( '.' ) ,
    faclvl_delim( '_' )
  { }     
};

options_t options;


typedef std::map<std::string, // indiv
	std::map<std::string, // r-strata
	std::map<std::string, // var
	std::map<std::string, // c-strata
		 value_t> > > > indexed_value_t;

//
// Global values to be populated across all DB, a unique list of all
// individuals (rows) and variables (cols), and col/row-stratifers
//

indexed_value_t val;  

//
// global ordering information for levels (within individual)
// track level order (i.e. as per packet order, which preserves numeric indices)
// rather than the mapped value (1, 10, 11, ..., 2, 20, ... ) 
//

// only used for rows, as cols have to preserve the same ordering
// across all individuals....  so, for now at least, all cols
// will just be sorted alpha-numeric

std::map<std::string, std::map<std::string,int> > rlvl_keys;//, clvl_keys;
std::map<std::string, std::map<int,std::string> > rlvl_order; //, clvl_order;


std::string print( const std::string & i, 
		   const std::string & r, 
		   const std::string & v, 
		   const std::string & c )
{

  if ( val.find( i ) == val.end() ) return "NA";
  
  std::map<std::string, 
    std::map<std::string, 
    std::map<std::string,value_t> > > & val1 = val[ i ];
  
  if ( val1.find( r ) == val1.end() ) return "NA";
  
  std::map<std::string, 
    std::map<std::string,value_t> > & val2 = val1[ r ];
  
  if ( val2.find( v ) == val2.end() ) return "NA";
  
  std::map<std::string,value_t> & val3 = val2[ v ];
  if ( val3.find( c ) == val3.end() ) return "NA";

  std::string rval = val3[ c ].str();
  if ( rval == "" ) return "NA";
  
  // no formating for numerics
  if ( options.full ) return rval;

  // is numeric?
  double dval;
  if ( Helper::str2dbl( rval , &dval ) )
    {
      std::stringstream ss;
      ss << std::fixed  << std::setprecision( options.prec ) << dval ;
      rval = ss.str();
    }
  
  return rval;

}

		   

std::set<std::string> o_ind, o_var, o_col, o_row;

std::map<std::string,std::string> strata2col_label; // i.e. as indexed in val[][]...
std::map<std::string,std::string> strata2row_label;
std::map<std::string,std::map<std::string,std::string> > row2fac2level;   // o_row ID --> fac / lvl

//
// Functions and structs
//

void dictionary();
void extract();
void display();
void summary();
void pre_summary();
void get_matching_strata( bool show_table = true );

struct request_t;
struct reqvar_t;

bool req_epoch, req_interval, req_timepoints;
bool rvar_timepoint, cvar_timepoint;

strata_t merge_strata( const strata_t & s1 , const strata_t & s2 );

std::vector<std::string> databases;

std::set<request_t> rvars ;
std::set<request_t> cvars ;
std::set<reqvar_t> vars ;
std::set<std::string> cfacs, rfacs;

std::set<int> cmds_id;
std::set<int> vars_id;
std::set<int> inds_id;

// populated for each dataset
std::set<int> match_strata_ids;

bool run_summary;
bool run_dictionary;

struct request_t 
{
  request_t( const std::string & r )
  {
    // split on first '/' only 
    std::vector<std::string> tok = Helper::parse( r , "/" );
    if ( tok.size() == 1 ) { fac = r; return; }
    // merge back, e.g. ANNOT/apnea/obstructive,apnea/central
    if ( tok.size() > 2 )
      {
	for (int i=2;i<tok.size();i++)
	  tok[1] += "/" + tok[i];	
      }
    std::vector<std::string> tok2 = Helper::parse( tok[1] , "," );
    fac = tok[0];
    for (int i=0;i<tok2.size();i++) levels.insert(tok2[i]);
  }
  // a factor -> 0 or more levels
  std::string fac;  
  std::set<std::string> levels;
  bool is_level_specific() const { return levels.size(); } 
  bool includes( const std::string & l ) const { return levels.find( l) != levels.end(); } 

  bool operator< ( const request_t & rhs ) const
  {
    if ( fac < rhs.fac ) return true;
    if ( fac > rhs.fac ) return false;
    if ( levels.size() < rhs.levels.size() ) return true;
    if ( levels.size() > rhs.levels.size() ) return false;

    std::set<std::string>::const_iterator ii = levels.begin();
    std::set<std::string>::const_iterator jj = rhs.levels.begin();
    while ( ii != levels.end() )
      {
	if ( *ii < *jj ) return true;
	if ( *jj < *ii ) return false;
	++ii; ++jj;
      }
    return false;
  }

};

struct reqvar_t {

  // can be command/var
  // or command/*   -- > command only
  // or var         -- > var only (of all commands)
  
  reqvar_t( const std::string & t )
  {
    var = cmd = "";
    std::vector<std::string> tok = Helper::parse( t , "/" );
    if ( tok.size() == 1 ) 
      { 
	const int n = tok[0].size()-1;
	if ( tok[0][n] == '/' ) { cmd = tok[0].substr(0,n); std::cout << "added " << cmd << "\n" ; } 
	else var = tok[0];
      } 
    else if ( tok.size() == 2 )
      {
	var=tok[0]; 
	cmd=tok[1]; 
      }
  }
  
  std::string cmd;
  std::string var;
  
  std::string str() const
  {
    return cmd == "" ? var : cmd + "/" + var;
  }

  bool operator< ( const reqvar_t & rhs ) const 
  {
    if ( cmd == rhs.cmd ) return var < rhs.var;
    return cmd < rhs.cmd;
  }
};





int main(int argc , char ** argv )
{
  
  // turn off logging
  global.api();


  if ( argc < 2 ) 
    Helper::halt( "usage: destrat stout.db {-f|-d|-s|-v|-i|-r|-c|-n|-e}" );
  
  //
  // Get command line options
  //

  char mode = 'D';

  run_summary = false; // dump by default, unless '-x'

  run_dictionary = false;

  std::string cmd_spec = "."; // which luna command, if taking output from -r/-c , i.e. we need a -s too ('statement') 

  bool any_opt = false; // either -x, -l, -d or -r/-c : otherwise do summary

  std::set<std::string> args_rvar, args_cvar, args_ind, args_var;
  
  for (int i=1;i<argc;i++)
    {

      if      ( strcmp( argv[i] , "-x" ) == 0 ) { run_summary = true; any_opt = true; mode = '0'; }
      else if ( strcmp( argv[i] , "-l" ) == 0 ) { run_summary = false; any_opt = true; options.long_format = true; mode = '0'; }
      else if ( strcmp( argv[i] , "-d" ) == 0 ) { run_dictionary = true; any_opt = true; mode = '0'; }
 
      else if ( strcmp( argv[i] , "-n" ) == 0 ) { options.print_cmd_name = true; mode = '0'; }
      else if ( strcmp( argv[i] , "-e" ) == 0 ) { options.print_empty_rows = true; mode = '0'; }

      else if ( strcmp( argv[i] , "-f" ) == 0 ) mode = 'D'; // database
      else if ( strcmp( argv[i] , "-a" ) == 0 ) mode = 'A'; // add prepend
      else if ( strcmp( argv[i] , "-s" ) == 0 ) { any_opt = true; mode = 'S'; } // luna statement
      else if ( strcmp( argv[i] , "-r" ) == 0 ) { any_opt = true; mode = 'R'; } // row-stratifier
      else if ( strcmp( argv[i] , "-c" ) == 0 ) { any_opt = true; mode = 'C'; } // col-stratifier
      
      else if ( strcmp( argv[i] , "-v" ) == 0 ) mode = 'V'; // variable name
      else if ( strcmp( argv[i] , "-i" ) == 0 ) mode = 'I'; // individual name
      else if ( strcmp( argv[i] , "-p" ) == 0 ) mode = 'P'; // set precision
            
      else // assume a variable name
	{
	  
	  // check it is not a command, e.g. [STATS] or #STATS
	  std::string s = argv[i] ;

	  if ( s[0] == '[' && s[ s.size()-1 ] == ']' )
	    {
	      if ( cmd_spec != "." ) Helper::halt( "cannot specify more than one [command] or -s command" );

	      std::string cmd_factor = "_" + s.substr(1,s.size()-2);
	      
	      if ( args_rvar.find( cmd_factor ) != args_rvar.end() ) 
		Helper::halt( "cannot have factor as both row and col stratifier " + std::string( cmd_factor ) );
	      args_rvar.insert( cmd_factor );
	      cmd_spec = cmd_factor;
	      
	      any_opt = true;

	      mode = '0';
	      
	    }
	  
	  // used +STATS instead of [STATS]
	  if ( s[0] == '+' )
	    {

	      
	      if ( cmd_spec != "." ) Helper::halt( "cannot specify more than one #command or [command]" );
	      
	      options.cmd_hash = true;
	      
	      std::string cmd_factor = "_" + s.substr(1,s.size()-1); // ignore lead #

	      if ( args_rvar.find( cmd_factor ) != args_rvar.end() ) 
		Helper::halt( "cannot have factor as both row and col stratifier " + std::string( cmd_factor ) );
	      args_rvar.insert( cmd_factor );
	      cmd_spec = cmd_factor;
	      
	      any_opt = true;

	      mode = '0';
	      
	    }

	  if ( mode == 'A' )
	    {
	      options.prepend = argv[i];
	    }
	  
	  if ( mode == 'D' ) 
	    {
	      databases.push_back( argv[i] );	      
	    }
	  
	  else if ( mode == 'R' ) 
	    { 
	      if ( args_cvar.find( argv[i] ) != args_cvar.end() ) 
		Helper::halt( "cannot have factor as both row and col stratifier " + std::string( argv[i] ) );
	      args_rvar.insert( argv[i] );
	    }
	  
	  else if ( mode == 'C' ) 
	    {
	      if ( args_rvar.find( argv[i] ) != args_rvar.end() ) 
		Helper::halt( "cannot have factor as both row and col stratifier " + std::string( argv[i] ) );
	      args_cvar.insert( argv[i] );
	    }
	  
	  else if ( mode == 'S' ) 
	    {
	      if ( cmd_spec != "." ) Helper::halt( "cannot specify more than one [command] or -s command" );

	      // preprend w/ underscore to indicate this 'factor' is 
	      // in fact a command string
	      std::string cmd_factor = "_" + std::string(argv[i]); 
	      
	      if ( args_rvar.find( cmd_factor ) != args_rvar.end() ) 
		Helper::halt( "cannot have factor as both row and col stratifier " + std::string( cmd_factor ) );
	      args_rvar.insert( cmd_factor );
	      cmd_spec = cmd_factor;

	    }

	  else if ( mode == 'V' ) 
	    {
	      vars.insert( reqvar_t( argv[i] ) );
	    }
	  
	  else if ( mode == 'I' ) 
	    {
	      args_ind.insert( argv[i] );
	    }

	  else if ( mode == 'P' )
	    {
	      int p = 2;
	      if ( ! Helper::str2int( argv[i] , &p ) ) Helper::halt( "expecting integer after -p for precision" );
	      options.prec = p;
	      options.full = false;
	    }
	}      
    }
  
  
  //
  // No options is same as '-x'
  //
  
  if ( ! any_opt )
    {
      run_summary = true; 
      mode = '0';
    }

  //
  // Check that DB were connected
  //

  if ( databases.size() == 0 ) 
    Helper::halt( "no STOUT databases specified" );

  //
  // For now, if >1 databse, no not allow col-stratifiers
  // (i.e. no check in place to get uniform cols across db yet)
  //

  if ( databases.size()>1 && args_cvar.size() > 0 ) 
    Helper::halt(" cannot specify -c with multiple attached databases currently") ;

  const bool IS_READONLY = true;


  //
  // Check variables
  //   
  

  std::set<reqvar_t> all_vars;

  bool verbose = databases.size() > 1; 

  if ( verbose )
    std::cerr << "attaching databases";
  
  for (int d=0;d<databases.size();d++)
    {
      
      if ( verbose )
	std::cerr << ".";
      
      if ( ! Helper::fileExists( databases[d] ) ) 
	Helper::halt( "could not find stout file " + databases[d] );
      
      if ( ! writer.attach( databases[d] , IS_READONLY ) )
	Helper::halt( "could not attach stout-file " + databases[d] );      
      
      // get all variables
      std::set<std::string> these_vars = writer.variable_names();
      std::set<std::string>::const_iterator vv = these_vars.begin();
      while ( vv != these_vars.end() ) { all_vars.insert( reqvar_t( *vv ) ); ++vv; }
      writer.close();
    }
  
  // if no variables explicitly, specified, then just add all
  if ( vars.size() == 0 ) 
    vars = all_vars;
  else
    {
      // otherwise, if some variables specified, each one has to exist in at least 1 DB
      std::set<reqvar_t>::const_iterator vv = vars.begin();
      while ( vv != vars.end() )
	{
	  if ( all_vars.find( *vv ) == all_vars.end() ) 
	    Helper::halt("could not find variable " + vv->str() + " in any databases" );
	  ++vv;
	}
    }
  
  if ( verbose ) 
    std::cerr << "\n";
  
  //
  // Map of strata-labels to row/col specifics: cfacs and rfacs
  //
    
  std::set<std::string>::const_iterator cc = args_cvar.begin();
  while ( cc != args_cvar.end() )
    {
      std::vector<std::string> tok = Helper::parse( *cc , "/" );
      cfacs.insert( tok[0] );   
      ++cc;                                                                                                                                                                                                           
    }                                                                                                                                                                                                                 
  
  std::set<std::string>::const_iterator rr = args_rvar.begin();
  while ( rr != args_rvar.end() )
    {
      std::vector<std::string> tok = Helper::parse( *rr , "/" );
      rfacs.insert( tok[0] );   
      ++rr;                                                                                                                                                                                                           
    }                                                                                                                                                                                                                 

  //
  // Iterate over each database
  //

  for (int d = 0 ; d < databases.size(); d++ )
    {
      
      if ( databases.size() > 1 ) 
	std::cerr << "scanning " << d+1 << " of " << databases.size() << ": " << databases[d] << "\n";

      //
      // ensure all tracking variables here are cleared
      //

      match_strata_ids.clear();


      //
      // Attach and read all information except value-store
      //
      
      writer.attach( databases[d] , IS_READONLY );

      // set index, i.e. for reading mode
      writer.index();

      // get all current information
      // NOTE -- attach() above will already do this
      //      writer.read_all();

      //
      // Check that factors are present
      //

      std::set<std::string>::const_iterator rr = args_rvar.begin();
      while ( rr != args_rvar.end() )
	{
	  std::vector<std::string> tok = Helper::parse( *rr , "/" );
	  std::string s = tok[0];
	  if ( writer.factors_idmap.find( s ) == writer.factors_idmap.end() && s != "E" && s != "T" ) 
	    {
	      if ( s[0] == '_' ) s = "[" + s.substr(1) + "] (command)"; 
	      Helper::halt( "could not find factor " + s );
	    }
	  rvars.insert( request_t( *rr  ) );	  
	  ++rr;
	}
      
      std::set<std::string>::const_iterator cc = args_cvar.begin();
      while ( cc != args_cvar.end() )
	{
	  std::vector<std::string> tok = Helper::parse( *cc , "/" );
	  std::string s = tok[0];
	  if ( writer.factors_idmap.find( s ) == writer.factors_idmap.end() && s != "E" && s != "T" ) 
	    {
	      if ( s[0] == '_' ) s = "[" + s.substr(1) + "] (command)"; 
	      Helper::halt( "could not find factor " + s );
	    }
	  cvars.insert( request_t( *cc ) );	  
	  ++cc;
	}


      //
      // Requested individuals
      //

      inds_id.clear();
      
      std::set<std::string>::const_iterator ii = args_ind.begin();
      while ( ii != args_ind.end() )
	{
	  if ( writer.individuals_idmap.find( *ii ) != writer.individuals_idmap.end() )
	    inds_id.insert( writer.individuals_idmap[ *ii ] );	      	  
	  ++ii;
	}


      //
      // Perform actions
      //
      
      if ( run_dictionary ) 
	{
	  dictionary();	  
	  // and skip to next database
	  writer.close();
	  continue;
	}

      //
      // Summary mode?
      //

      if ( run_summary ) 
	{
	  // only show the main table if we have 
	  // not specified *any* arguments, 
	  if ( ! any_opt ) 
	    pre_summary();
	}

      //
      // Check a command has been specified, if one is needed
      //
      
      if ( (!run_summary) && cmd_spec == "." ) 
	std::cerr << "*** did you forget to type the [COMMAND]?\n"
		  << "\n"
		  << "*** if not, this may be an old-format DB\n"
		  << "*** it should still be processed correctly\n"
		  << "*** but please update Luna and destrat\n";

      //
      // identify which rows we are interested in;  this function also will 
      // print the general table, but only if not any other options have
      // been given
      //
      
      get_matching_strata( !any_opt );


      //
      // Specific commands/variables?
      //
            
      std::set<std::string> req_vars, req_cmds;
      std::set<reqvar_t>::const_iterator vv = vars.begin();
      while ( vv != vars.end() ) 
	{
	  if ( vv->var != "" ) req_vars.insert( vv->var );
	  if ( vv->cmd != "" ) req_cmds.insert( vv->cmd );
	  ++vv;
	}
      
      vars_id = writer.all_matching_vars( req_vars );
      cmds_id = writer.all_matching_cmds( req_cmds );
      
      
      //
      // generate output
      //

      if ( run_summary ) 
	summary();
      else 
	extract(); // continue populating 'val' 
      
      //
      // Done, move to the next DB      
      //

      // close DB connection and wipe all caches

      writer.close();

    }
  
  // all done?

  if ( options.long_format ) std::exit(0);
  if ( run_summary ) std::exit(0);
  if ( run_dictionary ) std::exit(0);

  //
  // Display all (unless already done via long-format)
  //
  
  display();
  
  //
  // All done
  //

  std::exit(0);

}
  

struct fstrata_t 
{ 

  // second value denotes whether this factor is a
  // timepoint
  
  std::set<factor_t> factors; 
  
  fstrata_t() { } 
  
  std::string print() const 
  {
    std::stringstream ss;
    std::set<factor_t>::const_iterator ff = factors.begin();
    while ( ff != factors.end() )
      {
	if ( ff != factors.begin() ) ss << "x";
	ss << ff->factor_name;
	++ff;
      }
    return ss.str();
  }
  
  bool operator<( const fstrata_t & rhs ) const 
  {
    // time-stratified strata go last
    //    if ( epoch && ! rhs.epoch ) return false;
    //    if ( (!epoch) && rhs.epoch ) return true;

    //    if ( interval && ! rhs.interval ) return false;
    //    if ( (!interval) && rhs.interval ) return true;
    
    // find commands first
    std::string lcmd = "", rcmd = "";

    std::set<factor_t>::const_iterator ff = factors.begin();
    while ( ff != factors.end() )
      {
	if ( ff->factor_name.substr(0,1) == "_" ) 
	  {
	    lcmd = ff->factor_name.substr(1);
	    break;
	  }
	++ff;
      }

    std::set<factor_t>::const_iterator gg = rhs.factors.begin();
    while ( gg != rhs.factors.end() )
      {
	if ( gg->factor_name.substr(0,1) == "_" ) 
	  {
	    rcmd = gg->factor_name.substr(1);
	    break;
	  }
	++gg;
      }
    
    if ( lcmd < rcmd ) return true;
    if ( rcmd < lcmd ) return false;

    if ( factors.size() == rhs.factors.size() )
      {
	std::set<factor_t>::const_iterator ff = factors.begin();
	std::set<factor_t>::const_iterator gg = rhs.factors.begin();
	while ( ff != factors.end() )
	  {
	    if ( ff->factor_id < gg->factor_id ) return true;
	    if ( gg->factor_id < ff->factor_id ) return false;
	    ++ff; ++gg;
	  }
	return false;
      }
    return factors.size() < rhs.factors.size();
  }


  int matches( const std::set<factor_t> & fac , bool req_epoch , bool req_interval ) const
  {
    // +1 exact match
    //  0 does not contain all 'fac'
    // -1 contains additional factors beyond 'fac'
    
    // if ( req_epoch && ! epoch ) return 0;
    // if ( epoch && ! req_epoch ) return -1;
    
    // if ( req_interval && ! interval ) return 0;
    // if ( interval && ! req_interval ) return -1;
    
    bool additional = false;
    int match = 0 ;    
    std::set<factor_t>::const_iterator ff = factors.begin();
    while ( ff != factors.end() )
      { 
	if ( fac.find( *ff ) == fac.end() ) additional = true;
	else ++match;            
	++ff;
      }  
    if ( match < fac.size() ) return 0;    
    return additional ? -1 : +1 ;
  }
  
};

fstrata_t fmatch;


void pre_summary()
{

  std::cerr << "--------------------------------------------------------------------------------\n";
  std::cerr << writer.name() << ": ";
  std::cerr << writer.num_commands() << " command(s), ";
  std::cerr << writer.num_individuals() << " individual(s), ";
  std::cerr << writer.num_variables() << " variable(s), ";
  std::cerr << writer.num_values() << " values\n";
  std::cerr << "--------------------------------------------------------------------------------\n";

  //
  // Commands
  //
  
  std::map<int,command_t>::const_iterator cc = writer.commands.begin();
  while ( cc != writer.commands.end() )
    {
      std::cerr << "  command #" << cc->first << ":\t"
		<< "c" << cc->second.cmd_number << "\t"
		<< cc->second.timestamp << "\t"
		<< cc->second.cmd_name << "\t"
		<< cc->second.cmd_parameters << "\n";
      ++cc;
    }

  std::cerr << "--------------------------------------------------------------------------------\n";

}


void summary()
{


  //
  // Get variables/levels for the fstrata of interest
  //
    
  std::set<int>::const_iterator kk = match_strata_ids.begin();

  std::map<std::string,int> e_inds, e_cmds, e_vars;
  std::map<std::string,std::map<std::string,int> > e_faclvl;
  
  while ( kk != match_strata_ids.end() )
    {
      
      // strata details
      const strata_t & strata = writer.strata[ *kk ];
      
      // enumerate variables, for this strata
      packets_t packets = writer.enumerate( *kk );

      packets_t::const_iterator pp = packets.begin();
      while ( pp != packets.end() )
	{

	  const std::string & indiv_name = writer.individuals[ pp->indiv_id ].indiv_name ; 
	  const std::string & var_name = writer.variables[ pp->var_id ].var_name;
	  const std::string & cmd_name = writer.commands[ pp->cmd_id ].cmd_name;
	  
	  double count = pp->value.d;
	  
	  e_inds[ indiv_name ] += count;
	  e_vars[ cmd_name + "/" + var_name ] += count;
	  e_cmds[ cmd_name ] += count;

	  // Also count level/factor instances
	  std::map<factor_t,level_t>::const_iterator ll = strata.levels.begin();
	  while ( ll != strata.levels.end() )
	    {
	      e_faclvl[ ll->first.factor_name ][ ll->second.level_name ] += count;
	      ++ll;
	    }
	  
	  // next packet
	  ++pp;
	}


      // next strata            
      ++kk;
    }
  
  //
  // Report for fstrata of interest
  //
  
  bool empty = e_inds.size() == 0 || e_vars.size() == 0 ;
  
  if ( empty ) return;
  
  bool baseline_level = fmatch.factors.size() <= 1 ; 


  // -1 from factors size to exclude [COMMAND]
  std::cerr << "Factors: " << ( baseline_level ? "NA" : Helper::int2str( (int)fmatch.factors.size() -1 ) ) << "\n";
  
  if ( baseline_level ) 
    std::cerr << "     [ default/baseline ]\n\n";
  else
    {
      
      
      std::set<factor_t>::const_iterator gg = fmatch.factors.begin();
      while ( gg != fmatch.factors.end() )
	{
	  // command
	  if ( gg->factor_name[0] == '_' ) { ++gg; continue; } 

	  // timepoint?
	  bool is_tp = gg->factor_name == "E" || gg->factor_name == "T";
	  
	  if ( ! is_tp )
	    {
	      
	      std::cerr << "     [" << gg->factor_name << "] " << e_faclvl[ gg->factor_name ].size() << " levels\n     ->";

	      std::map<std::string,int> & ss = e_faclvl[ gg->factor_name ];
	      int cnt = 0;
	      std::map<std::string,int>::const_iterator ii = ss.begin();
	      while ( ii != ss.end() )
		{
		  if ( ii == ss.begin() )
		    std::cerr << " ";
		  else
		    std::cerr << ", ";

		  std::cerr << ii->first ;
		  ++cnt;
		  if( cnt > 12 ) { std::cerr << " ..." ; break; } 
		  ++ii;
		}
	      std::cerr << "\n";
	    }
	  else
	    std::cerr << "     [" << gg->factor_name << "] (time/epoch marker)\n";
	  
	  std::cerr << "\n";

	  ++gg;

	}
      
    }

  
  std::cerr << "Individuals: " << e_inds.size() << "\n";

  int c = 0;
  std::cerr << "    ";
  std::map<std::string,int>::const_iterator ii = e_inds.begin();
  while ( ii != e_inds.end() )
    {
      std::cerr << " " << ii->first ;
      //if ( ++c > 8 ) { std::cerr << "..."; break; } 
      if ( ++c > 8 ) { std::cerr << "\n     "; c = 0; } 
      ++ii;
    }
  
  std::cerr << "\n\n";
      
  std::cerr << "Commands: " << e_cmds.size() << "\n    ";
  c = 0;
  ii = e_cmds.begin();
  while ( ii != e_cmds.end() )
    {
      std::cerr << " " << ii->first ; // << "(" << ii->second << ")";
      //      if ( ++c > 6 ) { std::cerr << " ...";  break; } 
      if ( ++c > 6 ) { std::cerr << "\n     "; c = 0; } 
      ++ii;
    }
  std::cerr << "\n\n";

  std::cerr << "Variables: " << e_vars.size() << "\n    ";
  c = 0;
  ii = e_vars.begin();
  while ( ii != e_vars.end() )
    {
      std::cerr << " " << ii->first ; // << "(" << ii->second << ")";
      //      if ( ++c > 6 ) { std::cerr << " ...";  break; } 
      if ( ++c > 6 ) { std::cerr << "\n     "; c = 0; } 
      ++ii;
    }
  std::cerr << "\n";  
  
}



void get_matching_strata( bool show_table)
{

  
  //
  // factor-only strata (i.e. collapse across levels)
  //

  std::map<fstrata_t,int> fstrata;
  std::map<fstrata_t,std::set<int> > fstrata2strata_id;

  //
  // Always insert the baseline strata
  //

  fstrata_t baseline;
  fstrata[ baseline ] = 1;

  // and then actual strata
    
  std::map<int,strata_t>::const_iterator ss = writer.strata.begin();
  while ( ss != writer.strata.end() )
    {
      std::set<factor_t> factors;
      std::map<factor_t,level_t>::const_iterator ff = ss->second.levels.begin();
      while ( ff != ss->second.levels.end() )
	{
	  factors.insert( ff->first );
	  ++ff;
	}

      // store/increase count
      fstrata_t f;
      f.factors = factors;

      // tmp fix,i.e.. only add non-baseline factors
      if ( factors.size() > 0 ) 
	{
	  fstrata[ f ]++;	  
	  fstrata2strata_id[ f ].insert( ss->first );
	}

      ++ss;
    }
  

  //
  // Count number of actually observed strata 
  //
  
  //  std::map<int,int> StratOutDBase::count_strata()

  if ( run_summary && show_table ) 
    {
      
      std::map<int,std::set<int> > vars_by_strata = writer.dump_vars_by_strata();

      // ignore baseline
      std::cerr << "distinct strata group(s):\n";
      
      std::cerr << "  commands      : factors           : levels        : variables \n";
      std::cerr << "----------------:-------------------:---------------:---------------------------\n";

      std::map<fstrata_t,int>::const_iterator ff = fstrata.begin();
      while ( ff != fstrata.end() ) 
	{	  
	  
	  
	  // get variable names
	  std::set<std::string> vars;
	  std::set<int> strata_ids = fstrata2strata_id[ ff->first ];

	  std::set<int>::const_iterator ii = strata_ids.begin();
	  while ( ii != strata_ids.end() )
	    {
	      if ( vars_by_strata.find( *ii ) == vars_by_strata.end() ) 
		{
		  ++ii; continue;
		}
	      const std::set<int> & vs = vars_by_strata[ *ii ];
	      std::set<int>::const_iterator vv = vs.begin();
	      while ( vv != vs.end() )
		{
		  vars.insert( writer.variables[ *vv ].var_name );
		  ++vv;
		}
	      ++ii;
	    }
	  
	  // strata groups -- and expecting one "_COMMAND" strata

	  // baseline, handle separaetly
	  if ( ff->first.factors.size() == 0 )
	    {
	      // special ID == '1'
	      if ( vars_by_strata.find( 1 ) == vars_by_strata.end() ) 
		{
		  ++ff; continue;
		}
	      else
		{
		  const std::set<int> & vs = vars_by_strata[ 1 ];
		  std::set<int>::const_iterator vv = vs.begin();
		  while ( vv != vs.end() )
		    {
		      vars.insert( writer.variables[ *vv ].var_name );
		      ++vv;
		    }
		}

	      if ( vars.size() == 0 )
		{ 
		  ++ff; continue;
		}

	      std::cerr << "  " 
			<< std::left << std::setw( 14 ) << "[ NA ]" 
			<< std::left << std::setw( 20 ) << ": ." 
			<< std::left << std::setw( 16 ) << ": ."
			<< ":";
	    }

	  // display strata group
	  
	  if ( ff->first.factors.size() > 0 ) 
	    {
	      bool has_tp = false;
	      std::cerr << "  " ;	      
	      std::string msg = "[ NA ]";
	      
	      // first show _COMMAND strata if present
	      std::set<factor_t>::const_iterator gg = ff->first.factors.begin();
	      while ( gg != ff->first.factors.end() )
		{
		  if ( gg->factor_name[0] == '_' ) 
		    { 
		      msg = "[" + gg->factor_name.substr(1) + "]";		      
		    } 
		  ++gg;
		}
	      
	      std::cerr << std::left << std::setw(14) << msg;
	      
	      // then show normal strata
	      msg = ":";
	      gg = ff->first.factors.begin();
	      while ( gg != ff->first.factors.end() )
		{
		  if ( gg->factor_name[0] == '_' ) { ++gg; continue; }
		  msg += " " + gg->factor_name ;
		  if ( gg->factor_name == "E" || gg->factor_name == "T" ) has_tp = true;
		  ++gg;
		}
	      
	      if ( msg == ":" ) msg += " .";

	      std::cerr << std::left << std::setw(20) << msg;
	      
	      if ( has_tp ) 
		std::cerr << std::left << std::setw(16) << ": (...)" << ":";
	      else
		{
		  std::string msg = ": " + Helper::int2str(ff->second) + " level(s)";
		  std::cerr << std::left << std::setw(16) << msg << ":";
		}
	      
	    }
	  
	  int w = 0;
	  std::set<std::string>::const_iterator vv = vars.begin();
	  while ( vv != vars.end() )
	    {
	      std::cerr << " " << *vv;
	      w += 1 + vv->size();
	      if ( w > 30 ) 
		{
		  std::cerr << "\n" 
			    << "                :                   :               :";
		  w = 0;
		}
	      ++vv;
	    }
	  
	  // if ( ff->first.epoch ) { if ( ff->first.factors.size() > 0 ) { std::cerr << " x "; } std::cerr << "E"; } 
	  // if ( ff->first.interval ) { if ( ff->first.factors.size() > 0 ) { std::cerr << " x "; } std::cerr << "T"; } 
	  
	  std::cerr << "\n"
		    << "                :                   :               : \n";
	  ++ff;
	}

      std::cerr << "----------------:-------------------:---------------:---------------------------\n";

    }
  
  
  //
  // Specifics factors requested?
  //

  std::set<factor_t> requested;

  std::set<request_t>::const_iterator ii = rvars.begin();
  while ( ii != rvars.end() )
    {
      if ( writer.factors_idmap.find( ii->fac ) == writer.factors_idmap.end() ) 
	{ ++ii; continue; } 
      requested.insert( writer.factors[ writer.factors_idmap[ ii->fac ] ] );
      ++ii;
    }

  ii = cvars.begin();
  while ( ii != cvars.end() )
    {
      if ( writer.factors_idmap.find( ii->fac ) == writer.factors_idmap.end() ) 
	{ ++ii; continue; } 
      requested.insert( writer.factors[ writer.factors_idmap[ ii->fac ] ] );
      ++ii;
    }

  req_epoch = rvars.find( request_t("E") ) != rvars.end() || cvars.find( request_t("E") ) != cvars.end();
  req_interval = rvars.find( request_t("T") ) != rvars.end() || cvars.find( request_t("T") ) != cvars.end();
  req_timepoints = req_epoch || req_interval;
  
  rvar_timepoint = rvars.find( request_t("E") ) != rvars.end() || rvars.find( request_t("T") ) != rvars.end();
  cvar_timepoint = cvars.find( request_t("E") ) != cvars.end() || cvars.find( request_t("T") ) != cvars.end();

  //
  // Nothing to do?  ... should show base level here instead of returning..
  //
  
  if ( requested.size() == 0 && ! req_timepoints ) return;


  //
  // Consider each fstrata -- does this match with input variables?
  // It can, by definition, only match one fstrata (or none)
  //
    
  bool match_found = false;

  std::map<fstrata_t,int>::const_iterator mm = fstrata.begin();
  while ( mm != fstrata.end() )
    {
      if ( mm->first.matches( requested , req_epoch , req_interval ) == 1 ) 
	{
	  if ( match_found ) Helper::halt( "internal error") ;
	  match_found = true; 
	  fmatch = mm->first;
	}
      ++mm;
    }

  if ( ! match_found ) 
    {
      std::cerr << "No matching strata found.\n"; 
      std::exit(0);
    }


  //
  // Now report variable summary for matching strata
  // 
  
  // Find matching stata_t instances
  
  std::set<int> m0 = fstrata2strata_id[ fmatch ];

  //
  // And now prune out based on any level-specific criteria
  //
  
  std::set<int>::const_iterator qq = m0.begin();
  while ( qq != m0.end() )
    {
      const int strata_id = *qq;
      const strata_t & strata = writer.strata[ strata_id ];
      bool okay = true;
      
      std::set<request_t>::const_iterator rr = rvars.begin();
      while ( rr != rvars.end() )
	{	  
	  if ( rr->is_level_specific() ) 
	    {	      
	      const factor_t & factor = writer.factors[ writer.factors_idmap[ rr->fac ] ];
	      if ( strata.levels.find( factor ) == strata.levels.end() ) Helper::halt( "internal error" );
	      const level_t & level = strata.levels.find( factor )->second;
	      if ( ! rr->includes( level.level_name ) ) okay = false; 
	    }
	  ++rr;
	}

      std::set<request_t>::const_iterator cc = cvars.begin();
      while ( cc != cvars.end() )
	{	  
	  if ( cc->is_level_specific() ) 
	    {	      
	      const factor_t & factor = writer.factors[ writer.factors_idmap[ cc->fac ] ];
	      if ( strata.levels.find( factor ) == strata.levels.end() ) Helper::halt( "internal error" );
	      const level_t & level = strata.levels.find( factor )->second;
	      if ( ! cc->includes( level.level_name ) ) okay = false; 
	    }
	  ++cc;
	}
      
      if ( okay ) match_strata_ids.insert( strata_id );

      ++qq;
    }
 

}



void dictionary()
{  
  std::map<int,var_t>::const_iterator vv = writer.variables.begin();
  while ( vv != writer.variables.end() )
    {
      const var_t & var = vv->second;
      std::cout << writer.name() << "\t" 
		<< var.var_name << "\t"
		<< writer.commands[ var.cmd_id ].cmd_name << "\t"
		<< var.var_label << "\n";
      ++vv;
    }  
}


void extract()
{

  // note: vars_id will always be populated.
  
  std::set<int> * p_inds = inds_id.size() > 0 ? &inds_id : NULL ;
  std::set<int> * p_vars = vars_id.size() > 0 ? &vars_id : NULL ;
  std::set<int> * p_cmds = cmds_id.size() > 0 ? &cmds_id : NULL ;
    
  
  //
  // Fetch packets
  //

  packets_t packets;
  
  if ( match_strata_ids.size() > 0 ) 
    {
      std::set<int>::const_iterator kk = match_strata_ids.begin();
      while ( kk != match_strata_ids.end() )
	{	  
	  writer.fetch( *kk , req_timepoints, &packets, p_inds , p_cmds , p_vars );
	  ++kk;
	}
    }
  else
    {
      // strata_id of -1 implies root level (i.e. no stratifying variables)
      writer.fetch( -1 , req_timepoints, &packets, p_inds , p_cmds , p_vars );
    }
  


  
//   //
//   // Get time-points
//   //
  
//   std::set<int> timepoints;
//   std::set<int> empty_set;
//   if ( req_timepoints )
//     {
//       // get /all/ timepoints
//       std::map<int,timepoint_t>::const_iterator ii = writer.timepoints.begin();
//       while ( ii != writer.timepoints.end() )
// 	{
// 	  timepoints.insert( ii->first );
// 	  ++ii;
// 	}      
//     }
  

//   //
//   // Create row- and col-specific strata
//   //
  
//   // strata -> timepoints
//   std::map<strata_t,std::set<int> > rstrata, cstrata;
  
//   std::set<int> rvars_id, cvars_id;
  
//   std::set<request_t>::const_iterator rr = rvars.begin();
//   while ( rr != rvars.end() )
//     {
//       rvars_id.insert( writer.factors_idmap[ rr->fac ] );      
//       ++rr;
//     }
//   std::set<request_t>::const_iterator cc = cvars.begin();
//   while ( cc != cvars.end() )
//     {
//       cvars_id.insert( writer.factors_idmap[ cc->fac ] );      
//       ++cc;
//     }
  
   
//   //
//   // consider each matching strata
//   //

//   std::set<int>::const_iterator kk = match_strata_id.begin();
//   while ( kk != match_strata_ids.end() )
//     {
//       strata_t & s = writer.strata[ *kk ];
      
//       // row strata
//       strata_t rs;
//       rs.strata_id = 1;
//       bool radd = false;
      
//       std::map<factor_t,level_t>::const_iterator jj = s.levels.begin();
//       while ( jj != s.levels.end() )
// 	{
// 	  if ( rvars_id.find( jj->first.factor_id ) != rvars_id.end() ) { rs.levels[ jj->first ] = jj->second; radd = true; }
// 	  ++jj;
// 	}
//       if ( radd ) 
// 	{
// 	  if ( rvar_timepoint )
// 	    rstrata[ rs ] = timepoints;
// 	  else
// 	    rstrata[ rs ] = empty_set;
// 	}
      
//       // col strata
//       strata_t cs;
//       cs.strata_id = 2;
//       bool cadd = false;
//       jj = s.levels.begin();
//       while ( jj != s.levels.end() )
// 	{
// 	  if ( cvars_id.find( jj->first.factor_id ) != cvars_id.end() ) { cs.levels[ jj->first ] = jj->second; cadd = true; } 
// 	  ++jj;
// 	}

//       if ( cadd ) 
// 	{
// 	  if ( cvar_timepoint )
// 	    cstrata[ cs ] = timepoints;
// 	  else
// 	    cstrata[ cs ] = empty_set;
// 	}

//       ++kk;
//     }
  
  

  //
  // Convert packets_t to indexed_value_t
  //  
  
  // track strata-string (strata/timepoint) -> label for printing
  std::map<std::pair<int,int> ,std::string> strata_labels;
  
  // also build a unique list of all col-stratifiers 
  std::map<std::pair<int,int> ,std::string> col_strata_labels;
  
  packets_t::const_iterator pp = packets.begin();
  while ( pp != packets.end() )
    {

      std::pair<int,int> pr( pp->strata_id , pp->timepoint_id) ;
      
      // have we built this strata label already?
      if ( strata_labels.find( pr ) == strata_labels.end() )
	{
	  
	  strata_t strata = writer.strata[ pp->strata_id ];	  
	  
	  timepoint_t timepoint;
	  if ( pp->timepoint_id != -1 ) timepoint = writer.timepoints[ pp->timepoint_id ];
	  
	  std::stringstream ss;
	  if ( strata.levels.size() == 0 ) 
	    {
	      ss << ".";
	    }
	  else
	    {
	      
	      std::map<factor_t,level_t>::const_iterator ll = strata.levels.begin();
	      while ( ll != strata.levels.end() )
		{
		  
		  const std::string & fac = ll->first.factor_name ;
		  
		  // overall strata label
		  if      ( ll != strata.levels.begin() ) ss << ".";
		  if      ( timepoint.is_epoch() && fac == "E" ) ss << "E." << timepoint.epoch;
		  else if ( timepoint.is_interval() && fac == "T" ) ss << "E." << timepoint.start << "_" << timepoint.stop;
		  else    ss << fac << "." << ll->second.level_name;
		  
		  ++ll;
		}
	    }
	  
	  // save label
	  std::string strata_label = ss.str();
	  strata_labels[ pr ] = strata_label;
	  
	  //	  std::cerr << "made strata level [" << strata_label << "]\n";
	  // map onto col- and row specific labels
	  std::string clab, rlab;
	  std::map<std::string,std::string> rlabs;
	  std::map<factor_t,level_t>::const_iterator ll = strata.levels.begin();
	
	  //	  std::cerr << "considering " << strata.levels.size() << " levels\n";

	  while ( ll != strata.levels.end() )
	    {
	      const std::string & fac = ll->first.factor_name ;

	      std::string lvl = ll->second.level_name;

	      if      ( fac == "E" && timepoint.is_epoch() ) 
		lvl = Helper::int2str( timepoint.epoch );
	      else if ( fac == "T" && timepoint.is_interval() ) 
		lvl = Helper::int2str( timepoint.start ) + "_" + Helper::int2str( timepoint.stop );
	      
	      // col or row specific factors?
	      if ( cfacs.find( fac ) != cfacs.end() )
		{
		  if ( fac[0] != '_' ) 
		    { 
		      if ( clab.size() > 0 ) clab += options.strata_delim ;
		      clab += fac + options.faclvl_delim + lvl;
		    }
		}
	      else
		{
		  
		  if ( rlab.size() > 0 ) rlab += ".";
		  rlab += fac + options.faclvl_delim + lvl;
		  rlabs[ fac ] = lvl;

		}
	      
	      ++ll;
	    }
	
	  // no col/row strata?
	  if ( clab == "" ) clab = ".";
	  if ( rlab == "" ) rlab = ".";
	  
	  // track in global vars
	  
	  strata2col_label[ strata_label ] = clab;
	  strata2row_label[ strata_label ] = rlab;
	  row2fac2level[ rlab ] = rlabs;

	}

      std::string db_name     = writer.name();
      std::string indiv_name  = writer.individuals[ pp->indiv_id ].indiv_name;
      std::string cmd_name    = writer.commands[ pp->cmd_id ].cmd_name;
      std::string var_name    = writer.variables[ pp->var_id ].var_name;
      std::string strata_name = strata_labels[ pr ];
      std::string value_name  = pp->value.str();

      // build col- and row- specific stratifier labels
      std::string rstrata_name = strata2row_label[ strata_name ];
      std::string cstrata_name = strata2col_label[ strata_name ];

      //      std::cout << "[" << rstrata_name << "/" << cstrata_name << "]\n";

      // track observed individuals, etc
      o_ind.insert( indiv_name );
      o_var.insert( var_name );

      if ( cstrata_name != "." )
	o_col.insert( cstrata_name );
      
      if ( rstrata_name != "." )
	o_row.insert( rstrata_name );
      
      // long-format output?
      if ( options.long_format )
	{
	  
	  std::cout << db_name << "\t"
		    << indiv_name << "\t"
		    << cmd_name << "\t"
		    << strata_name << "\t"
	    // 		    << rstrata_name << "\t"
	    // 		    << cstrata_name << "\t"
		    << var_name << "\t"
		    << value_name << "\n";	  
	  
	}
      
      
      //
      // save in merged 'indexed_value_t' space
      //
      
      val[ indiv_name ][ rstrata_name ][ var_name ][ cstrata_name ] = value_name ;
      

      //
      // save ordering
      //
      
      if ( rlvl_keys[ indiv_name ].find( rstrata_name ) ==  rlvl_keys[ indiv_name ].end() )
	{
	  int rn = rlvl_keys[ indiv_name ].size();
	  rlvl_keys[ indiv_name ][ rstrata_name ] = rn ;
	  rlvl_order[ indiv_name ][ rn ] = rstrata_name ;
	  //	  std::cout << "rs = " <<indiv_name << " " << rstrata_name << " " << rn << "\t" << rlvl_keys[ indiv_name ].size() << " " << rlvl_order[ indiv_name ].size() << "\n"; 
	}

//       int cn = clvl_keys[ indiv_name ].size();
//       clvl_keys[ indiv_name ][ cstrata_name ] = cn ;
//       clvl_order[ indiv_name ][ cn ] = cstrata_name ;
      
      //
      // next packet
      //

      ++pp;
      
    }
  
  //
  // If in long-format mode, all done here.
  //
  
  if ( options.long_format ) return;


}



void display()
{

  // nothing to display?
  if ( o_ind.size() == 0 ) return;
  

  //
  // We now assume that all DB have been read through, so that val is populated
  //
  
  //
  // Header row, ID
  //

  std::cout << "ID";
  
  //
  // Row-stratifiers
  //
  
  std::set<std::string>::const_iterator ff = rfacs.begin();
  while ( ff != rfacs.end() )
    {
      if ( (*ff)[0] != '_' )  // skip command strata
	std::cout << "\t" << *ff ;
      ++ff;
    }
  

  //
  // Variables, looped by any column-stratifiers (cstrata)
  //
  
  std::set<std::string>::const_iterator vv = o_var.begin();
  while ( vv != o_var.end() )
    {
      
      const std::string & var_name = options.prepend + *vv;
      
      if ( o_col.size() == 0 ) 
	std::cout << "\t" << var_name;
      else
	{
	  std::set<std::string>::const_iterator cc = o_col.begin();
	  while ( cc != o_col.end() )
	    {
	      std::cout << "\t" << var_name << "." << *cc;
	      ++cc;
	    }
	}
      
      // 	      // any col-factor modifiers?
      // 	      std::map<strata_t,std::set<int> >::const_iterator h2 = cstrata.begin();
      // 	      while ( h2 != cstrata.end() )
      // 		{
      
      // 		  std::stringstream ss;
      // 		  std::map<factor_t,level_t>::const_iterator h3 = h2->first.levels.begin();
      // 		  while ( h3 != h2->first.levels.end() )
      // 		    {	      
      // 		      if ( (!cvar_timepoint) || ! ( h3->first.factor_name == "E" || h3->first.factor_name == "T" ) )
      // 			ss << "." << h3->first.factor_name << "." << h3->second.level_name;
      // 		      ++h3;
      // 		    } 
		  
// 		  if ( cvar_timepoint )
// 		    {
// 		      std::set<int>::const_iterator tt = h2->second.begin();
// 		      while ( tt != h2->second.end() )
// 			{
// 			  std::stringstream ss;
// 			  std::map<factor_t,level_t>::const_iterator h3 = h2->first.levels.begin();
// 			  while ( h3 != h2->first.levels.end() )
// 			    {	      
// 			      ss << "." << h3->first.factor_name << "." << h3->second.level_name;
// 			      ++h3;
// 			    } 
			  
// 			  // display
// 			  std::cout << "\t" << var_name << ss.str() << ".E." << writer.timepoints[ *tt ].print();
		      
// 			  ++tt;
// 			}
// 		    }
// 		  else
// 		    {		  
// 		      // display w/out time-point information
// 		      std::cout << "\t" << var_name << ss.str();
// 		    }
      
		  
// 		  ++h2;
// 		}
//	    }
      
      ++vv;
    }
  
  std::cout << "\n";
  
  

  //
  // Now process each row
  //
  
  std::set<std::string>::const_iterator oo = o_ind.begin();
  while ( oo != o_ind.end() )
    {
      
      const std::string indiv_name = *oo;
      
      //
      // no row-strata
      //
      
      if ( o_row.size() == 0 ) 
	{
	  
	  // start row
	  std::cout << indiv_name;
	  
	  // Variables, looped by any column-stratifiers (cstrata)
	  std::set<std::string>::const_iterator vv = o_var.begin();
	  while ( vv != o_var.end() )
	    {
	      
	      const std::string & var_name = *vv;
	      
	      if ( cfacs.size() == 0 ) 
		{
		  std::cout << "\t" << print( indiv_name , "." , var_name , "." );		  
		}
	      else
		{		  
		  std::set<std::string>::const_iterator cc = o_col.begin();
		  while ( cc != o_col.end() )
		    {
		      std::cout << "\t" << print( indiv_name , "." , var_name , *cc );
		      ++cc;
		    }
		}
	      
	      ++vv;
	    }
	  
	  // end of row
	  std::cout << "\n";
	  
	} 
      
      //
      // otherwise, loop over row-strata
      //
      
      else // loop over each row-strata 
	{

	  // doing this some that we preserve a better (original) ordering of rows (i.e. not 
	  // sorted as if strings)
	  
	  const std::map<int,std::string> & rorder = rlvl_order.find( indiv_name )->second;

	  std::map<int,std::string>::const_iterator rord = rorder.begin();
	  while ( rord != rorder.end() )
	    // 	  std::map<std::string,std::map<std::string,std::string> >::iterator rr = row2fac2level.begin();
	    // 	  while ( rr != row2fac2level.end() )
	    {
	      
	      std::map<std::string,std::map<std::string,std::string> >::iterator rr = row2fac2level.find( rord->second );	 

	      //
	      // Check whether this individual has the requiste row-variables
	      //
	      
	      bool has_this_rstrata = val[ indiv_name ].find( rr->first ) != val[ indiv_name ].end();
	      
	      if ( ( ! has_this_rstrata ) && ( ! options.print_empty_rows ) )
		{		  
		  ++rord;
		  continue;
		}
	      	      
	      //
	      // Start row
	      //
	      
	      std::cout << indiv_name;
	      
	      //
	      // row stratifiers
	      //
	      
	      std::map<std::string,std::string> & fac2lvl = rr->second;
	      
	      std::map<std::string,std::string>::iterator ff = fac2lvl.begin();
	      while ( ff != fac2lvl.end() )
		{
		  if ( ff->first[0] != '_' ) 
		    {
		      std::cout << "\t" << ff->second; // display level-value for this row-strata
		    }
		  ++ff;
		}
	      
	      //
	      // variables, w/ or w/out col-stratifiers
	      //
	      
	      std::set<std::string>::const_iterator vv = o_var.begin();
	      while ( vv != o_var.end() )
		{
		  
		  const std::string & var_name = *vv;
		  
		  if ( cfacs.size() == 0 ) 
		    {
		      std::cout << "\t" << print( indiv_name , rr->first , var_name , "." );		  
		    }
		  else
		    {		  
		      std::set<std::string>::const_iterator cc = o_col.begin();
		      while ( cc != o_col.end() )
			{
			  std::cout << "\t" << print( indiv_name , rr->first , var_name , *cc );
			  ++cc;
			}
		    }
		  
		  ++vv;
		}
	      	      
	      // end of line
	      std::cout << "\n";
	      
	      // next row strata	      
	      ++rord;
	    }
	}
      
      // next individual
      ++oo; 
    }

            
}




strata_t merge_strata( const strata_t & s1 , const strata_t & s2 )
{

  strata_t s;  
  std::map<factor_t,level_t>::const_iterator ii = s1.levels.begin();
  while ( ii != s1.levels.end() )
    {
      s.levels[ ii->first ] = ii->second;
      ++ii;
    }
  
  ii = s2.levels.begin();
  while ( ii != s2.levels.end() )
    {
      s.levels[ ii->first ] = ii->second;
      ++ii;
    }
  
  return s;
}
