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
#include <algorithm>
#include <unordered_map>
#include <utility>
#include <filesystem>
#include <cstdio>
#ifdef WINDOWS
#include <io.h>
#include <fcntl.h>
#else
#include <unistd.h>
#endif

// Default destrat reader: the original implementation is preserved in reader.cpp.

#include "defs/defs.h"
#include "helper/helper.h"
#include "db/db.h"

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
  bool compressed = false; // -z
  std::string dump_directory; // -t
  
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

// Data output is independent of diagnostics, which always stay on stderr.
std::ostream * data_output = &std::cout;
std::ostream & output() { return *data_output; }

struct output_file_t
{
  std::ofstream plain;
  gzofstream gzip;
  bool is_gzip = false;
  std::string label = "stdout";

  void open( const std::string & filename, bool compressed )
  {
    label = filename.empty() ? "stdout" : filename;
    is_gzip = compressed;
    if ( compressed )
      {
        if ( filename.empty() )
          {
#ifdef WINDOWS
            const int fd = _dup( _fileno(stdout) );
            if ( fd >= 0 ) _setmode( fd, _O_BINARY );
#else
            const int fd = dup( fileno(stdout) );
#endif
            if ( fd < 0 ) Helper::halt( "could not duplicate stdout for gzip output" );
            gzip.attach( fd, std::ios::out | std::ios::binary );
          }
        else
          gzip.open( filename, std::ios::out | std::ios::binary );
        if ( ! gzip ) Helper::halt( "could not open gzip output: " + label );
        data_output = &gzip;
      }
    else if ( ! filename.empty() )
      {
        plain.open( filename, std::ios::out | std::ios::binary );
        if ( ! plain ) Helper::halt( "could not open output: " + label );
        data_output = &plain;
      }
    else
      data_output = &std::cout;
  }

  void finish()
  {
    output().flush();
    const bool failed = ! output();
    data_output = &std::cout;
    if ( is_gzip )
      {
        gzip.close(); // Finalize the gzip trailer, including for empty output.
        if ( failed || ! gzip ) Helper::halt( "error writing gzip output: " + label );
      }
    else if ( plain.is_open() )
      {
        plain.close();
        if ( failed || ! plain ) Helper::halt( "error writing output: " + label );
      }
    else if ( failed )
      Helper::halt( "error writing stdout" );
  }
};


//
// Functions and structs
//

void dictionary();
void extract();
void display();
void summary();
void pre_summary();
bool get_matching_strata( bool show_table = true );

struct request_t;

bool req_epoch, req_interval, req_timepoints;
bool rvar_timepoint, cvar_timepoint;

bool baseline_request = false;
std::string current_database;
SQL reader_sql;
void open_reader_database( const std::string &, bool create_index = true );
void close_reader_database();
void dump_all( const std::string &, const std::set<std::string> &, bool );
bool dump_include_root = false;

std::vector<std::string> databases;

std::set<request_t> rvars ;
std::set<request_t> cvars ;
std::set<std::string> vars; // literal stored variable names; -s selects a command
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





int main(int argc , char ** argv )
{
  
  // Output uses C++ streams exclusively.
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);

  // This is a CLI: Helper::halt must terminate on invalid input/output errors.
  globals::bail_on_fail = true;
  // turn off logging
  global.api();


  if ( argc < 2 ) 
    Helper::halt( "usage: destrat file.db [more.db ...] [-s CMD] [-r FAC ...] [-c FAC ...] [-v VAR ...] [-i ID ...] [-e] [-z] [-t DIRECTORY]" );
  
  //
  // Get command line options
  //

  char mode = 'D';

  run_summary = false; // dump by default, unless '-x'

  run_dictionary = false;

  std::string cmd_spec = "."; // which luna command, if taking output from -r/-c , i.e. we need a -s too ('statement') 

  bool any_opt = false; // either -x, -l, -d or -r/-c : otherwise do summary

  std::set<std::string> args_rvar, args_cvar, args_ind, args_var;
  std::set<std::string> warned_ambiguous_specs;

  // A slash introduces a level-qualified factor request, whose suffix is
  // already comma-delimited (e.g. STG/N2,N3).  Consequently, only expand a
  // comma-separated command-line list when the token has no slash; otherwise
  // preserving it as one request avoids changing that established syntax.
  const auto option_list = []( const std::string & arg )
    {
      if ( arg.find( '/' ) != std::string::npos )
        return std::vector<std::string>{ arg };
      return Helper::parse( arg , "," );
    };
  
  for (int i=1;i<argc;i++)
    {

      if      ( strcmp( argv[i] , "-x" ) == 0 ) { run_summary = true; any_opt = true; mode = '0'; }
      else if ( strcmp( argv[i] , "-l" ) == 0 ) { run_summary = false; any_opt = true; options.long_format = true; mode = '0'; }
      else if ( strcmp( argv[i] , "-d" ) == 0 ) { run_dictionary = true; any_opt = true; mode = '0'; }
 
      else if ( strcmp( argv[i] , "-n" ) == 0 ) { options.print_cmd_name = true; mode = '0'; }
      else if ( strcmp( argv[i] , "-e" ) == 0 ) { options.print_empty_rows = true; mode = '0'; }
      else if ( strcmp( argv[i] , "-z" ) == 0 ) options.compressed = true;
      else if ( strcmp( argv[i] , "-t" ) == 0 )
        {
          if ( i + 1 >= argc || argv[i+1][0] == '\0' || argv[i+1][0] == '-' )
            Helper::halt( "expecting an output directory after -t" );
          if ( ! options.dump_directory.empty() ) Helper::halt( "specify -t only once" );
          options.dump_directory = argv[++i];
          any_opt = true;
          mode = 'D';
        }

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
	      for ( const auto & item : option_list( argv[i] ) )
	        {
	          if ( args_cvar.find( item ) != args_cvar.end() )
	            Helper::halt( "cannot have factor as both row and col stratifier " + item );
	          args_rvar.insert( item );
	        }
	    }
	  
	  else if ( mode == 'C' ) 
	    {
	      for ( const auto & item : option_list( argv[i] ) )
	        {
	          if ( args_rvar.find( item ) != args_rvar.end() )
	            Helper::halt( "cannot have factor as both row and col stratifier " + item );
	          args_cvar.insert( item );
	        }
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
	      for ( const auto & item : option_list( argv[i] ) ) vars.insert( item );
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

  if ( ! options.dump_directory.empty() )
    {
      bool row_selection = false;
      for ( const auto & factor : args_rvar )
        if ( factor.empty() || factor[0] != '_' ) row_selection = true;
      if ( run_summary || run_dictionary || options.long_format || row_selection || ! args_cvar.empty() )
        Helper::halt( "-t writes all factor sets as rows; it cannot be combined with -r, -c, -l, -x or -d" );
    }

  //
  // Multiple input files may have different individuals, variables and levels.
  //

  // Column identities are merged by factor names and levels across files.

  const bool IS_READONLY = true;


  //
  // Check variables
  //   
  

  std::set<std::string> all_vars;
  const bool requested_variables = ! vars.empty();

  bool verbose = databases.size() > 1; 

  if ( verbose )
    std::cerr << "attaching databases";
  
  for (int d=0;d<databases.size();d++)
    {
      
      if ( verbose )
	std::cerr << ".";
      
      if ( ! Helper::fileExists( databases[d] ) ) 
	Helper::halt( "could not find stout file " + databases[d] );
      
      // Discovery needs only names, not the writer's full metadata caches.
      SQL discovery;
      discovery.open( databases[d] );
      sqlite3_stmt * names = discovery.prepare(
        "SELECT DISTINCT variable_name FROM variables;" );
      if ( ! names ) Helper::halt( "could not read variables from " + databases[d] );
      while ( discovery.step( names ) )
        all_vars.insert( discovery.get_text( names, 0 ) );
      discovery.finalise( names );
      discovery.close();
    }
  
  // if no variables explicitly, specified, then just add all
  if ( vars.size() == 0 ) 
    vars = all_vars;
  else
    {
      // otherwise, if some variables specified, each one has to exist in at least 1 DB
      std::set<std::string>::const_iterator vv = vars.begin();
      while ( vv != vars.end() )
	{
	  if ( all_vars.find( *vv ) == all_vars.end() ) 
	    Helper::halt("could not find variable " + *vv + " in any databases; -v takes literal names, use -s CMD to select a command" );
	  ++vv;
	}
    }
  
  if ( verbose ) 
    std::cerr << "\n";

  if ( ! options.dump_directory.empty() )
    {
      dump_all( cmd_spec, args_ind, requested_variables );
      return 0;
    }

  output_file_t stdout_output;
  stdout_output.open( "", options.compressed );
  
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
      rvars.clear();
      cvars.clear();


      //
      // Attach and read all information except value-store
      //
      
      current_database = databases[d];
      if ( run_summary || run_dictionary )
        {
          writer.attach( databases[d], IS_READONLY );
          writer.index();
        }
      else
	open_reader_database( databases[d] );

      // A factor/level request can legally have slash-containing levels (for
      // example, ANNOT/apnea/obstructive,apnea/central).  Warn only when a
      // comma-delimited component also begins with a factor that exists in
      // this database: that is the common accidental form
      // STG/N2,N3,CH/Cz,Fz, which is otherwise interpreted as STG levels.
      const auto warn_ambiguous_comma_list = [&]( const std::string & spec , const char option )
	{
	  const size_t slash = spec.find( '/' );
	  if ( slash == std::string::npos ) return;
	  const std::vector<std::string> levels = Helper::parse( spec.substr( slash + 1 ) , "," );
	  for ( const auto & level : levels )
	    {
	      const size_t nested_slash = level.find( '/' );
	      if ( nested_slash == std::string::npos ) continue;
	      const std::string possible_factor = level.substr( 0 , nested_slash );
	      if ( writer.factors_idmap.find( possible_factor ) == writer.factors_idmap.end() ) continue;
	      const std::string warning_key = std::string( 1 , option ) + ":" + spec;
	      if ( warned_ambiguous_specs.insert( warning_key ).second )
		std::cerr << "warning : " << option << " " << spec
			  << " is interpreted as levels of " << spec.substr( 0 , slash )
			  << "; separate factor requests with spaces (e.g. -" << option
			  << " STG/N2,N3 CH/Cz,Fz)\n";
	      return;
	    }
	};

      bool unavailable_factors = false;

      //
      // Check that factors are present
      //

      std::set<std::string>::const_iterator rr = args_rvar.begin();
      while ( rr != args_rvar.end() )
	{
	  warn_ambiguous_comma_list( *rr , 'r' );
	  std::vector<std::string> tok = Helper::parse( *rr , "/" );
	  std::string s = tok[0];
	  if ( writer.factors_idmap.find( s ) == writer.factors_idmap.end() && s != "E" && s != "T" ) 
	    {
	      if ( s[0] == '_' ) s = "[" + s.substr(1) + "] (command)"; 
	      if ( databases.size() == 1 ) Helper::halt( "could not find factor " + s );
              unavailable_factors = true;
	    }
	  rvars.insert( request_t( *rr  ) );	  
	  ++rr;
	}
      
      std::set<std::string>::const_iterator cc = args_cvar.begin();
      while ( cc != args_cvar.end() )
	{
	  warn_ambiguous_comma_list( *cc , 'c' );
	  std::vector<std::string> tok = Helper::parse( *cc , "/" );
	  std::string s = tok[0];
	  if ( writer.factors_idmap.find( s ) == writer.factors_idmap.end() && s != "E" && s != "T" ) 
	    {
	      if ( s[0] == '_' ) s = "[" + s.substr(1) + "] (command)"; 
	      if ( databases.size() == 1 ) Helper::halt( "could not find factor " + s );
              unavailable_factors = true;
	    }
	  cvars.insert( request_t( *cc ) );	  
	  ++cc;
	}


      if ( unavailable_factors )
        {
          std::cerr << "skipping " << current_database << ": requested factors not present\n";
          close_reader_database();
          continue;
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
	  close_reader_database();
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
      
      if ( ! get_matching_strata( !any_opt ) )
        {
          close_reader_database();
          continue;
        }


      //
      // Literal variable selection; -s has already selected the command strata.
      //
            
      vars_id.clear();
      cmds_id.clear();
      for ( const auto & v : writer.variables )
        if ( vars.count( v.second.var_name ) ) vars_id.insert( v.first );

      // An explicit filter with no local matches must not select everything.
      if ( ( ! args_ind.empty() && inds_id.empty() )
           || ( ! vars.empty() && vars_id.empty() ) )
        {
          close_reader_database();
          continue;
        }
      
      
      //
      // generate output
      //

      if ( run_summary ) 
	summary();
      else 
	extract(); // stream long output, or merge cells for wide output
      
      //
      // Done, move to the next DB      
      //

      // close DB connection and wipe all caches

      close_reader_database();

    }
  
  // all done?

  if ( options.long_format || run_summary || run_dictionary )
    {
      stdout_output.finish();
      return 0;
    }

  //
  // Display all (unless already done via long-format)
  //
  
  display();
  
  //
  // All done
  //

  stdout_output.finish();
  return 0;

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

  std::map<std::string,int64_t> e_inds, e_cmds, e_vars;
  std::map<std::string,std::map<std::string,int64_t> > e_faclvl;
  
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
	  
	  const int count = pp->value.i;
	  
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

	      std::map<std::string,int64_t> & ss = e_faclvl[ gg->factor_name ];
	      int cnt = 0;
	      std::map<std::string,int64_t>::const_iterator ii = ss.begin();
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
  std::map<std::string,int64_t>::const_iterator ii = e_inds.begin();
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



bool get_matching_strata( bool show_table)
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
  
  baseline_request = requested.empty() && ! req_timepoints;
  if ( baseline_request ) return true;


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
      std::string requested_cmd = ".";
      std::stringstream reqss;
      bool first = true;

      std::set<request_t>::const_iterator ri = rvars.begin();
      while ( ri != rvars.end() )
        {
          if ( ri->fac.size() > 0 && ri->fac[0] == '_' )
            {
              requested_cmd = ri->fac;
              break;
            }
          ++ri;
        }

      if ( requested_cmd == "." )
        {
          ri = cvars.begin();
          while ( ri != cvars.end() )
            {
              if ( ri->fac.size() > 0 && ri->fac[0] == '_' )
                {
                  requested_cmd = ri->fac;
                  break;
                }
              ++ri;
            }
        }

      if ( requested_cmd != "." )
        {
          reqss << "[" << requested_cmd.substr(1) << "]";
          first = false;
        }

      ri = rvars.begin();
      while ( ri != rvars.end() )
        {
          if ( ri->fac == "E" || ri->fac == "T" ) { ++ri; continue; }
          if ( ri->fac.size() > 0 && ri->fac[0] == '_' ) { ++ri; continue; }
          if ( ! first ) reqss << " ";
          reqss << ri->fac;
          first = false;
          ++ri;
        }

      ri = cvars.begin();
      while ( ri != cvars.end() )
        {
          if ( ri->fac == "E" || ri->fac == "T" ) { ++ri; continue; }
          if ( ri->fac.size() > 0 && ri->fac[0] == '_' ) { ++ri; continue; }
          if ( ! first ) reqss << " ";
          reqss << ri->fac;
          first = false;
          ++ri;
        }

      std::stringstream availss;
      bool first_avail = true;
      std::map<fstrata_t,int>::const_iterator fm = fstrata.begin();
      while ( fm != fstrata.end() )
        {
          bool has_requested_cmd = requested_cmd == ".";
          std::stringstream tblss;
          bool first_tbl = true;

          std::set<factor_t>::const_iterator ff = fm->first.factors.begin();
          while ( ff != fm->first.factors.end() )
            {
              const std::string & fac = ff->factor_name;

              if ( fac.size() > 0 && fac[0] == '_' )
                {
                  if ( requested_cmd != "." && fac == requested_cmd ) has_requested_cmd = true;
                  ++ff;
                  continue;
                }

              if ( ! first_tbl ) tblss << ",";
              tblss << fac;
              first_tbl = false;
              ++ff;
            }

          if ( ! has_requested_cmd ) { ++fm; continue; }

          const std::string label = tblss.str().size() == 0 ? "." : tblss.str();
          if ( ! first_avail ) availss << " | ";
          availss << label;
          first_avail = false;
          ++fm;
        }

      std::string msg = "no matching strata found for request "
        + std::string( reqss.str().size() ? reqss.str() : "[baseline]" );

      if ( availss.str().size() )
        msg += "; available strata: " + availss.str();

      if ( databases.size() > 1 )
        {
          std::cerr << "skipping " << current_database << ": " << msg << "\n";
          return false;
        }
      Helper::halt( msg );
    }

  std::set<request_t>::const_iterator ri = rvars.begin();
  while ( ri != rvars.end() )
    {
      if ( ri->fac == "E" || ri->fac == "T" ) { ++ri; continue; }
      if ( ri->fac.size() > 0 && ri->fac[0] == '_' ) { ++ri; continue; }
      if ( writer.factors_idmap.find( ri->fac ) == writer.factors_idmap.end() ) { ++ri; continue; }

      const factor_t & factor = writer.factors[ writer.factors_idmap[ ri->fac ] ];
      if ( fmatch.factors.find( factor ) == fmatch.factors.end() )
        Helper::halt( "requested row factor " + ri->fac + " is not present in the matched strata" );

      ++ri;
    }

  ri = cvars.begin();
  while ( ri != cvars.end() )
    {
      if ( ri->fac == "E" || ri->fac == "T" ) { ++ri; continue; }
      if ( ri->fac.size() > 0 && ri->fac[0] == '_' ) { ++ri; continue; }
      if ( writer.factors_idmap.find( ri->fac ) == writer.factors_idmap.end() ) { ++ri; continue; }

      const factor_t & factor = writer.factors[ writer.factors_idmap[ ri->fac ] ];
      if ( fmatch.factors.find( factor ) == fmatch.factors.end() )
        Helper::halt( "requested column factor " + ri->fac + " is not present in the matched strata" );

      ++ri;
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
 

  return true;
}



void dictionary()
{  
  std::map<int,var_t>::const_iterator vv = writer.variables.begin();
  while ( vv != writer.variables.end() )
    {
      const var_t & var = vv->second;
      output() << writer.name() << "\t" 
		<< var.var_name << "\t"
		<< writer.commands[ var.cmd_id ].cmd_name << "\t"
		<< var.var_label << "\n";
      ++vv;
    }  
}


// Only extraction uses this connection. Summary/dictionary retain the shared
// writer interface; no changes to the original reader or database library.
template <typename F>
void read_metadata( const std::string & query, F consume )
{
  sqlite3_stmt * stmt = reader_sql.prepare( query );
  if ( ! stmt ) Helper::halt( "could not read metadata from " + current_database );
  while ( reader_sql.step( stmt ) ) consume( stmt );
  reader_sql.finalise( stmt );
}

void open_reader_database( const std::string & filename, bool create_index )
{
  writer.clear();
  reader_sql.open( filename );
  reader_sql.synchronous( false );
  // Keep the original on-demand index creation, including its warning if the
  // filesystem does not permit creation. An existing index needs no write.
  if ( create_index )
    reader_sql.query( "CREATE INDEX IF NOT EXISTS vIndex ON datapoints(strata_id);" );

  read_metadata( "SELECT indiv_id,indiv_name FROM individuals;", [](sqlite3_stmt * s) {
    indiv_t v;
    v.clear();
    v.indiv_id = reader_sql.get_int( s, 0 );
    v.indiv_name = reader_sql.get_text( s, 1 );
    writer.individuals_idmap[v.indiv_name] = v.indiv_id;
    writer.individuals.emplace( v.indiv_id, std::move(v) );
  } );
  read_metadata( "SELECT cmd_id,cmd_name,cmd_number,cmd_timestamp,cmd_parameters FROM commands;", [](sqlite3_stmt * s) {
    command_t v;
    v.cmd_id = reader_sql.get_int( s, 0 );
    v.cmd_name = reader_sql.get_text( s, 1 );
    v.cmd_number = reader_sql.get_int( s, 2 );
    v.timestamp = reader_sql.get_text( s, 3 );
    v.cmd_parameters = reader_sql.get_text( s, 4 );
    writer.commands_idmap[v.cmd_name] = v.cmd_id;
    writer.commands.emplace( v.cmd_id, std::move(v) );
  } );
  read_metadata( "SELECT factor_id,factor_name,is_numeric FROM factors;", [](sqlite3_stmt * s) {
    factor_t v;
    v.factor_id = reader_sql.get_int( s, 0 );
    v.factor_name = reader_sql.get_text( s, 1 );
    v.is_numeric = reader_sql.get_int( s, 2 ) == 1;
    writer.factors_idmap[v.factor_name] = v.factor_id;
    writer.factors.emplace( v.factor_id, std::move(v) );
  } );
  read_metadata( "SELECT level_id,factor_id,level_name FROM levels;", [](sqlite3_stmt * s) {
    level_t v;
    v.level_id = reader_sql.get_int( s, 0 );
    v.factor_id = reader_sql.get_int( s, 1 );
    v.level_name = reader_sql.get_text( s, 2 );
    writer.levels.emplace( v.level_id, std::move(v) );
  } );
  read_metadata( "SELECT strata_id,level_id FROM strata;", [](sqlite3_stmt * s) {
    const int id = reader_sql.get_int( s, 0 );
    const int level = reader_sql.get_int( s, 1 );
    auto & stratum = writer.strata[id];
    stratum.strata_id = id;
    if ( level )
      {
        const auto ll = writer.levels.find( level );
        if ( ll == writer.levels.end() ) Helper::halt( "undefined stratum level in " + current_database );
        const auto ff = writer.factors.find( ll->second.factor_id );
        if ( ff == writer.factors.end() ) Helper::halt( "undefined factor in " + current_database );
        stratum.insert( ll->second, ff->second );
      }
  } );
  read_metadata( "SELECT variable_id,variable_name,command_name,variable_label FROM variables;", [](sqlite3_stmt * s) {
    var_t v;
    v.var_id = reader_sql.get_int( s, 0 );
    v.var_name = reader_sql.get_text( s, 1 );
    const auto cmd = writer.commands_idmap.find( reader_sql.get_text( s, 2 ) );
    v.cmd_id = cmd == writer.commands_idmap.end() ? -1 : cmd->second;
    v.var_label = reader_sql.get_text( s, 3 );
    writer.variables.emplace( v.var_id, std::move(v) );
  } );
  // Timepoints are joined only for selected datapoints. Reverse maps used by
  // the writer to allocate IDs are unnecessary for extraction.
}

void close_reader_database()
{
  if ( run_summary || run_dictionary ) writer.close();
  else
    {
      reader_sql.close();
      writer.clear();
    }
}

using factor_levels_t = std::map<std::string,std::string>;
using cell_key_t = std::pair<size_t,size_t>; // variable ID, column-stratum ID

struct cell_hash_t
{
  size_t operator()( const cell_key_t & k ) const
  {
    return std::hash<size_t>()( k.first )
      ^ ( std::hash<size_t>()( k.second ) + 0x9e3779b9
          + ( k.first << 6 ) + ( k.first >> 2 ) );
  }
};

struct output_row_t
{
  std::unordered_map<cell_key_t,std::string,cell_hash_t> cells;
};

struct output_individual_t
{
  std::unordered_map<size_t,output_row_t> rows;
  std::vector<size_t> row_order;
};

std::map<std::string,output_individual_t> output_individuals;
std::map<std::string,size_t> output_variables;
std::map<factor_levels_t,size_t> output_row_ids, output_column_ids;
// std::map keys have stable addresses, so row metadata is stored only once.
std::vector<const factor_levels_t *> output_row_levels;
std::vector<std::string> output_column_labels;
std::vector<std::string> column_factor_order;
bool have_column_factor_order = false;

std::string sql_ids( const std::set<int> & ids )
{
  std::string result;
  for ( const int id : ids )
    {
      if ( ! result.empty() ) result += ',';
      result += std::to_string( id );
    }
  return result;
}

struct labels_t
{
  size_t row = 0, column = 0;
  std::string long_label;
};

labels_t make_labels( int strata_id, const timepoint_t & timepoint )
{
  labels_t labels;
  factor_levels_t rows, columns;
  const auto stratum = writer.strata.find( strata_id );
  if ( stratum != writer.strata.end() )
    for ( const auto & entry : stratum->second.levels )
      {
        const std::string & factor = entry.first.factor_name;
        if ( ! options.dump_directory.empty() && ! factor.empty() && factor[0] == '_' )
          continue;
        std::string level = entry.second.level_name;
        if ( factor == "E" && timepoint.is_epoch() )
          level = Helper::int2str( timepoint.epoch );
        else if ( factor == "T" && timepoint.is_interval() )
          level = Helper::int2str( timepoint.start )
            + ( options.dump_directory.empty() ? "_" : "-" )
            + Helper::int2str( timepoint.stop );

        if ( options.long_format )
          {
            if ( ! labels.long_label.empty() ) labels.long_label += '.';
            // Preserve legacy long-output labels, including its E. prefix for T.
            const std::string prefix = factor == "T" && timepoint.is_interval() ? "E" : factor;
            labels.long_label += prefix + "." + level;
          }
        else if ( cfacs.count( factor ) )
          {
            if ( ! factor.empty() && factor[0] != '_' ) columns.emplace( factor, level );
          }
        else
          rows.emplace( factor, level );
      }

  if ( options.long_format )
    {
      if ( labels.long_label.empty() ) labels.long_label = ".";
      return labels;
    }

  auto row = output_row_ids.emplace( std::move(rows), output_row_ids.size() );
  labels.row = row.first->second;
  if ( row.second ) output_row_levels.push_back( &row.first->first );

  auto column = output_column_ids.emplace( std::move(columns), output_column_ids.size() );
  labels.column = column.first->second;
  if ( column.second )
    {
      std::string label;
      for ( const auto & factor : column_factor_order )
        {
          const auto level = column.first->first.find( factor );
          if ( level == column.first->first.end() ) continue;
          if ( ! label.empty() ) label += options.strata_delim;
          label += factor + options.faclvl_delim + level->second;
        }
      output_column_labels.push_back( label.empty() ? "." : std::move(label) );
    }
  return labels;
}

void extract()
{
  if ( ! baseline_request && match_strata_ids.empty() && ! dump_include_root ) return;

  std::string query = "SELECT d.indiv_id,d.cmd_id,d.variable_id,d.strata_id,d.timepoint_id,d.value";
  if ( req_timepoints ) query += ",t.epoch,t.start,t.stop";
  query += " FROM datapoints d";
  if ( req_timepoints ) query += " LEFT JOIN timepoints t ON t.timepoint_id=d.timepoint_id";
  if ( ! options.dump_directory.empty() )
    {
      query += " WHERE (";
      if ( dump_include_root ) query += "d.strata_id IS NULL";
      if ( ! match_strata_ids.empty() )
        {
          if ( dump_include_root ) query += " OR ";
          query += "d.strata_id IN (" + sql_ids( match_strata_ids ) + ")";
        }
      query += ")";
    }
  else
    {
      query += baseline_request ? " WHERE d.strata_id IS NULL" :
        " WHERE d.strata_id IN (" + sql_ids( match_strata_ids ) + ")";
      query += req_timepoints ? " AND d.timepoint_id IS NOT NULL" : " AND d.timepoint_id IS NULL";
    }
  if ( ! inds_id.empty() ) query += " AND d.indiv_id IN (" + sql_ids( inds_id ) + ")";
  if ( ! vars_id.empty() ) query += " AND d.variable_id IN (" + sql_ids( vars_id ) + ")";
  if ( ! cmds_id.empty() ) query += " AND d.cmd_id IN (" + sql_ids( cmds_id ) + ")";
  // Legacy extraction visits ascending strata IDs, then insertion order within
  // each stratum. Do not alphabetically sort E/N labels or reorder duplicate cells.
  query += baseline_request ? " ORDER BY d.rowid;" : " ORDER BY d.strata_id,d.rowid;";

  sqlite3_stmt * stmt = reader_sql.prepare( query );
  if ( ! stmt ) Helper::halt( "could not prepare extraction from " + current_database );

  int previous_stratum = -2, previous_timepoint = -2;
  labels_t labels;
  // Reuse label IDs within a stratum; long output retains just its last label.
  std::unordered_map<int,labels_t> label_cache;
  std::unordered_map<int,output_individual_t *> individual_cache;
  std::unordered_map<int,size_t> variable_cache;

  while ( reader_sql.step( stmt ) )
    {
      if ( ! options.long_format && ! have_column_factor_order )
        {
          // Preserve a single file's existing column labels; use that same
          // factor order for subsequent files, irrespective of their local IDs.
          for ( const auto & factor : writer.factors )
            if ( cfacs.count( factor.second.factor_name )
                 && ! factor.second.factor_name.empty()
                 && factor.second.factor_name[0] != '_' )
              column_factor_order.push_back( factor.second.factor_name );
          have_column_factor_order = true;
        }

      const int stratum = reader_sql.is_null( stmt, 3 ) ? -1 : reader_sql.get_int( stmt, 3 );
      const int tp = req_timepoints && ! reader_sql.is_null( stmt, 4 ) ? reader_sql.get_int( stmt, 4 ) : -1;
      if ( stratum != previous_stratum || tp != previous_timepoint )
        {
          if ( stratum != previous_stratum ) label_cache.clear();
          auto cached = label_cache.find( tp );
          if ( ! options.long_format && cached != label_cache.end() )
            labels = cached->second;
          else
            {
              timepoint_t timepoint;
              if ( req_timepoints )
                {
                  timepoint.epoch = reader_sql.is_null( stmt, 6 ) ? -1 : reader_sql.get_int( stmt, 6 );
                  if ( ! reader_sql.is_null( stmt, 7 ) )
                    {
                      timepoint.start = reader_sql.get_uint64( stmt, 7 );
                      timepoint.stop = reader_sql.get_uint64( stmt, 8 );
                    }
                }
              labels = make_labels( stratum, timepoint );
              if ( ! options.long_format ) label_cache.emplace( tp, labels );
            }
          previous_stratum = stratum;
          previous_timepoint = tp;
        }

      const int indiv_id = reader_sql.get_int( stmt, 0 );
      const int var_id = reader_sql.get_int( stmt, 2 );
      std::string value = reader_sql.get_text( stmt, 5 );
      if ( options.long_format )
        {
          output() << current_database << '\t'
                    << writer.individuals.at( indiv_id ).indiv_name << '\t'
                    << writer.commands.at( reader_sql.get_int( stmt, 1 ) ).cmd_name << '\t'
                    << labels.long_label << '\t'
                    << writer.variables.at( var_id ).var_name << '\t'
                    << value << '\n';
          continue;
        }

      auto person = individual_cache.find( indiv_id );
      if ( person == individual_cache.end() )
        person = individual_cache.emplace( indiv_id,
          &output_individuals[writer.individuals.at( indiv_id ).indiv_name] ).first;
      auto variable = variable_cache.find( var_id );
      if ( variable == variable_cache.end() )
        {
          const auto global_var = output_variables.emplace(
            writer.variables.at( var_id ).var_name, output_variables.size() );
          variable = variable_cache.emplace( var_id, global_var.first->second ).first;
        }

      auto & individual = *person->second;
      auto row = individual.rows.try_emplace( labels.row );
      if ( row.second ) individual.row_order.push_back( labels.row );
      // Later values replace earlier ones, including across input files.
      row.first->second.cells[cell_key_t( variable->second, labels.column )] = std::move(value);
    }
  reader_sql.finalise( stmt );
}

std::string format_value( const std::string & value )
{
  if ( value.empty() ) return "NA";
  if ( options.full ) return value;
  double number;
  if ( ! Helper::str2dbl( value, &number ) ) return value;
  std::ostringstream formatted;
  formatted << std::fixed << std::setprecision( options.prec ) << number;
  return formatted.str();
}

void clear_output_table()
{
  output_individuals.clear();
  output_variables.clear();
  output_row_levels.clear();
  output_row_ids.clear();
  output_column_ids.clear();
  output_column_labels.clear();
  column_factor_order.clear();
  have_column_factor_order = false;
}

void display()
{
  if ( output_individuals.empty() ) return;

  std::vector<size_t> columns;
  for ( size_t i = 0; i < output_column_labels.size(); ++i ) columns.push_back( i );
  std::stable_sort( columns.begin(), columns.end(), [](size_t a, size_t b) {
    return output_column_labels[a] < output_column_labels[b];
  } );

  output() << "ID";
  for ( const auto & factor : rfacs )
    if ( ! factor.empty() && factor[0] != '_' ) output() << '\t' << factor;

  std::vector<cell_key_t> cells;
  for ( const auto & variable : output_variables )
    for ( const size_t column : columns )
      {
        output() << '\t' << options.prepend << variable.first;
        if ( output_column_labels[column] != "." )
          output() << '.' << output_column_labels[column];
        cells.emplace_back( variable.second, column );
      }
  output() << '\n';

  // Padding traverses shared row IDs in first-encounter order. Missing cells
  // are printed directly, never inserted into the sparse value store.
  const output_row_t empty_row;
  for ( const auto & person : output_individuals )
    for ( size_t r = 0, n = options.print_empty_rows ? output_row_levels.size()
            : person.second.row_order.size(); r < n; ++r )
      {
        const size_t row_id = options.print_empty_rows ? r : person.second.row_order[r];
        const auto found = person.second.rows.find( row_id );
        const auto & row = found == person.second.rows.end() ? empty_row : found->second;
        const auto & levels = *output_row_levels[row_id];
        output() << person.first;
        for ( const auto & factor : rfacs )
          if ( ! factor.empty() && factor[0] != '_' )
            {
              const auto level = levels.find( factor );
              output() << '\t' << ( level == levels.end() ? "NA" : level->second );
            }
        for ( const auto & key : cells )
          {
            const auto cell = row.cells.find( key );
            output() << '\t';
            if ( cell == row.cells.end() || cell->second.empty() ) output() << "NA";
            else if ( options.full ) output() << cell->second;
            else output() << format_value( cell->second );
          }
        output() << '\n';
      }
}

// The catalog stores only table identities and source strata IDs. Values are
// merged one table at a time; input indexes are created during discovery.
struct dump_source_t
{
  std::set<int> strata, commands;
  bool root = false;
};

struct dump_table_t
{
  std::vector<std::string> factor_order;
  std::map<size_t,dump_source_t> sources;
  std::filesystem::path filename;
};

using dump_key_t = std::pair<std::string,std::set<std::string>>;

// Preserve ordinary Luna names while keeping database labels inside -t's
// directory. Encode '%' too so literal escape sequences remain distinguishable.
std::string filename_component( const std::string & name )
{
  static const char hex[] = "0123456789ABCDEF";
  std::string encoded;
  for ( const unsigned char c : name )
    if ( ( c >= 'a' && c <= 'z' ) || ( c >= 'A' && c <= 'Z' )
         || ( c >= '0' && c <= '9' ) || c == '_' || c == '-' || c == '.' )
      encoded += c;
    else
      {
        encoded += '%';
        encoded += hex[c >> 4];
        encoded += hex[c & 15];
      }
  return encoded;
}

bool dump_selection( const std::set<std::string> & people, bool filter_variables )
{
  inds_id.clear();
  vars_id.clear();
  for ( const auto & person : writer.individuals )
    if ( people.count( person.second.indiv_name ) ) inds_id.insert( person.first );
  for ( const auto & variable : writer.variables )
    if ( ! filter_variables || vars.count( variable.second.var_name ) ) vars_id.insert( variable.first );
  return ( people.empty() || ! inds_id.empty() ) && ! vars_id.empty();
}

void dump_all( const std::string & command, const std::set<std::string> & people, bool filter_variables )
{
  namespace fs = std::filesystem;
  std::map<dump_key_t,dump_table_t> tables;
  for ( size_t d = 0; d < databases.size(); ++d )
    {
      current_database = databases[d];
      open_reader_database( current_database );
      if ( ! dump_selection( people, filter_variables ) )
        {
          close_reader_database();
          continue;
        }
      std::string query = "SELECT DISTINCT cmd_id,strata_id FROM datapoints WHERE variable_id IN ("
        + sql_ids( vars_id ) + ")";
      if ( ! inds_id.empty() ) query += " AND indiv_id IN (" + sql_ids( inds_id ) + ")";
      query += " ORDER BY cmd_id,strata_id;";
      read_metadata( query, [&](sqlite3_stmt * stmt) {
        const int cmd_id = reader_sql.get_int( stmt, 0 );
        const auto cmd = writer.commands.find( cmd_id );
        if ( cmd == writer.commands.end() ) Helper::halt( "undefined command in " + current_database );
        if ( command != "." && command != "_" + cmd->second.cmd_name ) return;
        const bool root = reader_sql.is_null( stmt, 1 );
        const int sid = root ? -1 : reader_sql.get_int( stmt, 1 );
        std::vector<std::string> order;
        if ( ! root )
          {
            const auto stratum = writer.strata.find( sid );
            if ( stratum == writer.strata.end() ) Helper::halt( "undefined stratum in " + current_database );
            for ( const auto & level : stratum->second.levels )
              {
                const auto & name = level.first.factor_name;
                if ( ! name.empty() && name[0] != '_' ) order.push_back( name );
              }
          }
        dump_key_t key( cmd->second.cmd_name, std::set<std::string>( order.begin(), order.end() ) );
        auto added = tables.try_emplace( key );
        auto & table = added.first->second;
        if ( added.second ) table.factor_order = std::move(order);
        auto & source = table.sources[d];
        source.commands.insert( cmd_id );
        if ( root ) source.root = true;
        else source.strata.insert( sid );
      } );
      close_reader_database();
    }

  const fs::path directory( Helper::expand( options.dump_directory ) );
  std::set<fs::path> destinations;
  // Check all names before writing any table: separate factor sets must never
  // overwrite each other, and existing exports are left intact.
  for ( auto & entry : tables )
    {
      std::string name = filename_component( entry.first.first );
      if ( name.empty() ) name = "NA";
      for ( const auto & factor : entry.second.factor_order )
        name += "_" + filename_component( factor );
      entry.second.filename = directory / ( name + ( options.compressed ? ".txt.gz" : ".txt" ) );
      if ( ! destinations.insert( entry.second.filename ).second )
        Helper::halt( "different tables produce the same output filename: " + entry.second.filename.string() );
      std::error_code error;
      const auto status = fs::symlink_status( entry.second.filename, error );
      if ( error && error != std::errc::no_such_file_or_directory )
        Helper::halt( "cannot inspect output path: " + entry.second.filename.string() + ": " + error.message() );
      if ( fs::exists( status ) )
        Helper::halt( "output already exists: " + entry.second.filename.string() );
    }
  std::error_code error;
  fs::create_directories( directory, error );
  if ( error ) Helper::halt( "cannot create output directory: " + directory.string() + ": " + error.message() );

  size_t active = databases.size();
  size_t written = 0;
  for ( const auto & entry : tables )
    {
      clear_output_table();
      rfacs = entry.first.second;
      cfacs.clear();
      baseline_request = false;
      // Include NULL as well as non-NULL timepoints when dumping stored data.
      req_timepoints = true;
      for ( const auto & part : entry.second.sources )
        {
          if ( active != part.first )
            {
              if ( active != databases.size() ) close_reader_database();
              active = part.first;
              current_database = databases[active];
              open_reader_database( current_database, false );
            }
          if ( ! dump_selection( people, filter_variables ) ) continue;
          match_strata_ids = part.second.strata;
          cmds_id = part.second.commands;
          dump_include_root = part.second.root;
          extract();
        }
      if ( output_individuals.empty() ) continue;
      output_file_t file;
      file.open( entry.second.filename.string(), options.compressed );
      display();
      file.finish();
      ++written;
    }
  if ( active != databases.size() ) close_reader_database();
  clear_output_table();
  std::cerr << "wrote " << written << " combined table(s) to " << directory.string() << "\n";
}
